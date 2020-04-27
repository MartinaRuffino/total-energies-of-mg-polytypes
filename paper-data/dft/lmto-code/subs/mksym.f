      subroutine mksym(mode,slabl,ssymgr,s_ctrl,s_lat,s_spec,s_site)
C- Setup for symmetry group
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nspec nsite npadl npadr ldlm lpgf ipc idcc
Co     Stored:     idcc nccomp ndlmcl ipcp nclass nclasp
Co     Allocated:  *
Cio    Elts passed:ips ipc ics nl idcc ncomp ipcp dclabl nrc clabl
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat ldist dist
Co     Stored:     s_sym npgrp nsgrp afmt
Co     Allocated:  ag symgr istab
Cio    Elts passed:pos symgr ag istab lsym
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  nang ncpa nthet iscpa ncomp
Co     Stored:     nthet ncomp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  upd_nthet dlmidx dlmclbl
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  class
Co     Stored:     dlmcl pos class clabel
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvsym2
Ci Inputs
Ci   mode  : 1s digit
Ci           0  make space group only
Ci           1  add inversion to point group (see aginv below)
Ci           2  Same as 1, but make additionally ag,istab for extra
Ci              operations, using -g for rotation part; see Remarks
Ci           10s digit
Ci           0  do nothing about classes
Ci           1  split species into classes
Ci           2  Also assign class labels to each class
Ci           4  Assume class and species are the same.
Ci           100s digit
Ci           0  Leave s_lat%s_sym untouched
Ci           1  Allocate s_lat%s_sym and fill with g,ag,inv
Ci           2  Also fill s_sym%mib, s_sym%tib
Ci   slabl : species labels
Ci   ssymgr: string containing symmetry group generators.
Ci           if ssymgr contains 'find', mksym will add basis atoms as
Ci           needed to guarantee generators are valid, and generate
Ci           internally any additonal group operations needed to
Ci           complete the space group.
Co Outputs
Co  ssymgr:  may be enlarged if it contains 'find' initially
Cs Command-line switches
Cs   --fixlat : Adjust lattice vectors and group operations to render
Cs            : them internally consistent
Cs   --fixpos : Make small adjustments to given positions if
Cs            : they are compatiable with, or increase symmetry
Cr Remarks
Cr   In certain cases the inversion operation may be added to the space
Cr   group, for purposes of k integration.  This is permissible when the
Cr   hamiltonian has the form h(-k) = h*(k).  In that case, the
Cr   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
Cr
Cr   Also, the Green's functions are related G(-k) = Gtranspose(k).
Cr   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
Cr   then G(-k) = Gtranspose(k), and the same (g,ag) information is
Cr   needed for either rotation.
Cr
Cr   In addition there is a special "antiferromagnetic" group operation.
Cr   It acts in the same way as a regular space group operation except that
Cr   the spins are exchanged under the operation.  This special operation
Cr   has the restriction that each site gets mapped into
Cr   When this special operation is used for symmetrization,
Cr   it is performed after the regular symmetrization, which combines
Cr   space group op 1 (unit operator) and the special afm-op
Cr   A check is made that under this op sites (ib,jb) combined in pairs,
Cr   so that for symmetry purposes every site belongs to a class with 2 elements.
Cr   The two sites are symmetrized with the second atom density spin-flipped.
Cr   The special afm-op is appended to s_lat%istab, s_lat%symgr, s_lat%ag
Cr   in column s_lat%nsafm.
Cu Updates
Cu   07 Aug 17 Revise AFM symmetrization
Cu   22 Nov 14 Synchronize pgf with Kirill's dlm revisions
Cu   19 Aug 13 (Liqin Ke) Reduction in symmetry from spin orbit coupling
Cu             parting of RHOPOS and GRP2 removed from this routine
Cu   01 Aug 13 Parse SYMOPS_RHOPOS and SYMOPS_GRP2
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   11 Mar 13 species which have no sites are assigned no class
Cu   01 Feb 13 s_lat%s_sym is optionally allocated and filled
Cu   04 Jan 13 s_site%dlmcl matches s_site%class in absence of DLM
Cu   01 Dec 12 Completed migration to f90 structures
Cu   01 Sep 11 Begin migration to f90 structures
Cu   04 Jan 10 Enabled keyword 'AFM:#,#,#' in ssymgr
Cu   04 Jan 06 Enabled keyword 'point' in ssymgr
Cu             Generates extra group ops when generators create more
Cu             than the maximum expected
Cu   07 Oct 05 Put altered positions into site structure
Cu   03 Nov 01 Adapted to altered gensym
Cu   26 Jan 01 Added new '2' option to 1s digit of mode
Cu   31 May 00 revised meaning of mode
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character*8 slabl(1),ssymgr*(*)
      integer mode,nsgrp,npgrp
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_lat):: s_lat
      type(str_site):: s_site(*)
      type(str_spec):: s_spec(*)
C ... Dynamically allocated arrays
      integer,allocatable :: nrspec(:),nrcs(:),ips2(:)
      real(8),allocatable :: pos2(:)
C     character(len=8), allocatable :: clabel(:)
C ... Local parameters
      logical T,F,ltmp,lafms
      integer idest,ig,incli,j,j1,j2,ldist,lpgf,
     .  mks_sym,nbas,nbas0,nbasp,nclasp,nclass,nclspp,ngen,nggen,ngmx,
     .  ngnmx,npadl,npadr,nsite,nspec,nspec0,stdo,usegen,ngenso,nafmg
      integer isym(10)
      procedure(logical) :: cmdopt,a2bin
      procedure(integer) :: a2vec,iprint,nglob,aginv
      procedure(real(8)) :: dsum
      parameter (T=.true., F=.false., ngnmx=10)
      character*120 gens,strn*72,genso,afsg*72
      double precision gen(9,ngnmx),plat(3,3),dist(3,3),qlat(3,3),xx,
     .  pltol,fptol,agen(3,ngnmx),wk(3,3)
C ... For DLM
      integer ldlm,nccomp,ndlmcl,ic

C ... Debugging check ... to compare against mptauof
C      real(8),allocatable:: tiat(:,:,:),shtvg(:,:),dlmm(:,:,:,:)
C      integer,allocatable:: invg(:),miat(:,:)

      nbas = s_ctrl%nbas; nspec = s_ctrl%nspec
      nsite = s_ctrl%nsite  ! nsite includes padding elements (lmpg)
      plat = s_lat%plat; ldist = s_lat%ldist; dist = s_lat%dist
      ngmx = 48
      stdo = nglob('stdo')
      if (iprint() >= 30) write(stdo,'(1x)')

C ... Re-entry when ngmx has been increased
    5 continue

C ... Arrays to be allocated and initialized
      call ptr_ctrl(s_ctrl,1,'ics',max(nsite,nspec),0,0,xx)
      call ptr_ctrl(s_ctrl,1,'ipc',nsite,0,0,xx)
      call ptr_lat(s_lat,1,'ag',3,ngmx+1,0,xx)
      call ptr_lat(s_lat,1,'symgr',9,ngmx+1,0,xx)
      call ptr_lat(s_lat,1,'istab',10*nbas,ngmx+1,0,xx)

C ... Extract keywords from input string; copy rest to gens
C     keyword 'point' is appended to gens
      call words(ssymgr,ngen)
      j1 = 1
      idest = 1
      usegen = 2
      gens = ' '; afsg = ' '
      ltmp = .false.; lafms = .false.
      do  ig = 1, ngen
        call word(ssymgr,ig,j1,j2)
        if (ssymgr(j1:j2) == 'find') then
          usegen = 0
        else if (ssymgr(j1:j2) == 'point') then
          ltmp = .true.
        else if (ssymgr(j1:j1+3) == 'afm:' .or.
     .           ssymgr(j1:j1+3) == 'AFM:') then
          afsg = ssymgr(j1+4:j2)
        else if (ssymgr(j1:j1+4) == 'afmx:' .or.
     .           ssymgr(j1:j1+4) == 'AFMX:') then
          afsg = ssymgr(j1+5:j2)
          lafms = .true.
        else
          call strncp(gens,ssymgr,idest,j1,j2-j1+2)
          idest = idest+j2-j1+2
        endif
      enddo
      if (ltmp) call strncp(gens,'point',idest,1,5)

C --- Generate space group ---
      nbas0 = nbas
      if (cmdopt('--fixpos',8,0,strn)) then
        j1 = 8+1
        if (strn(9:13) == ':tol=') then
          j1 = 13
        endif
        if (strn(9:9) /= ':' .or.
     .    .not. a2bin(strn,fptol,4,0,' ',j1,len(strn))) fptol = 1d-5
      else
        fptol = 0
      endif
      if (cmdopt('--fixlat',8,0,strn)) then
        j1 = 8+1
        if (strn(9:13) == ':tol=') then
          j1 = 13
        endif
        if (strn(9:9) /= ':' .or.
     .    .not. a2bin(strn,pltol,4,0,' ',j1,len(strn))) pltol = 1d-6
      else
        pltol = 0
      endif

C ... When generating the group the basis may become enlarged ...
C     copy larger files relevant arrays to larger workspace
      allocate(pos2(10*3*ngmx*nbas))
      allocate(ips2(10*ngmx*nbas))
      call icopy(nbas,s_ctrl%ips,1,ips2,1)
      call dcopy(3*nbas,s_lat%pos,1,pos2,1)
      nspec0 = nspec
      allocate(nrspec(nspec+1))
      call gensym(slabl,gens,usegen,T,F,pltol,fptol,F,nbas,nspec0,ngmx,plat,
     .  plat,ldist,dist,pos2,ips2,nrspec,nsgrp,
     .  s_lat%symgr,s_lat%ag,ngen,gen,ssymgr,nggen,isym,s_lat%istab)
      if (pltol /= 0) s_lat%plat = plat
      deallocate(nrspec)

      if (nbas > nbas0) call rxs('gensym: the basis was enlarged.',
     .  ' Check group operations.')
      if (nggen > nsgrp) then
        call info2(10,0,0,' MKSYM (warning): generators create more '
     .    //'than ngmx=%i group ops ...',ngmx,0)
        ngmx = ngmx*16
        deallocate(pos2,ips2)
        goto 5
      endif
C     Update pos.  For now pos is stored in all of these places:
C     s_site%pos, s_lat%pos, s_ctrl%pos => s_lat%pos
      if (fptol /= 0) then
        call dcopy(3*nbas,pos2,1,s_lat%pos,1)
        call sitepack(s_site,1,nbas,'-pos',3,xx,s_lat%pos)
C       call prmx('s_lat%pos',s_lat%pos,3,3,nbas)
      endif
      deallocate(pos2,ips2)

C --- Add inversion to point group ---
      incli = -1
      npgrp = nsgrp
      if (mod(mode,10) /= 0) then
        ngen = ngen+1
        call dpzero(gen(1,ngen),9)
        gen(1,ngen) = -1
        gen(5,ngen) = -1
        gen(9,ngen) = -1
        call pshpr(max(iprint()-40,1))
        call grpgen(gen(1,ngen),1,s_lat%symgr,npgrp,ngmx)
        call poppr
        incli = npgrp-nsgrp
      endif

C --- Reduce operations not conforming to magnetic symmetry ---
      if (mod(s_lat%lsym,4) /= 0) then
        j = 1
        if (mod(s_lat%lsym,4) > 1) j = -1
        j1 = npgrp ! hang on to original number of ops
        call sogrp(s_lat%symgr,s_lat%ag,npgrp,j)
C       Assume corresponding reduction in nsgrp
        if ((j1/npgrp)*npgrp /= j1) call rx('bug in sogrp')
        nsgrp = (nsgrp * npgrp) / j1

        if (iprint() >= 30) then
          ngenso = 0
          call groupg(0,nsgrp,s_lat%symgr,s_lat%ag,
     .      plat,ngenso,gen,agen,genso,nggen)

C         Printout for original coordinate system
          if (ldist == 2) then
            call dinv33(dist,0,wk,xx)
            call pshpr(1)
            call lattdf(-ldist,wk,plat,0,xx,ngenso,gen,agen)
            gens = ' '; genso = ' ';
            call asymopn(ngenso,gen,agen,plat,gens(9:),genso(9:))
            call poppr
            call info0(2,0,0,'%9f... in unrotated coordinates:')
            call awrit0('%a',gens,len(gens),-stdo)
            call awrit0('%a',genso,len(genso),-stdo)
          endif
        endif
      endif

C --- Printout of symmetry operations ---
      call info5(20,0,0,' MKSYM:  found %i space group operations'//
     .  '%?#n#; adding inversion generated %i ops#%j#'//
     .  '%?#(n==0)# ... includes inversion##',nsgrp,npgrp-nsgrp,npgrp,
     .  incli,0)
      if (incli == -1) then
        call info0(30,0,0,'         no attempt to add inversion symmetry')
      endif
      if (mod(mode/10,10) == 0) goto 100

C --- Split species into classes, or copy species into classes ---
      call icopy(nsite,s_ctrl%ips,1,s_ctrl%ipc,1)
      allocate(nrcs(nspec)) ! nrcs(is) = number of classes in species is
      call splcls(mod(mode/10,10) == 4,s_lat%pos,nbas,nsgrp,s_lat%istab,
     .  nspec,slabl,nclass,s_ctrl%ipc,s_ctrl%ics,nrcs)
      deallocate(nrcs)  ! Not kept for anything for now

C ... Generate s_lat%s_sym
      if (mod(mode/100,10) /= 0) then
        allocate(p_sym(nsgrp))
        j = 11; if (mod(mode/100,10) > 1) j = 1011
C       j = j + 10000
        j = mks_sym(j,plat,s_lat%symgr,s_lat%ag,nsgrp,nbas,s_lat%pos,
     .    s_ctrl%nl,p_sym)
        s_lat%s_sym => p_sym
      endif

C ... Debugging check ... to compare against mptauof
C      allocate(miat(nbas,nsgrp),tiat(3,nbas,nsgrp),invg(nsgrp),
C     .  shtvg(3,nsgrp))
C      call pshpr(51)
C      call mptauof(s_lat%symgr,nsgrp,plat,nbas,s_lat%pos,s_ctrl%ipc,
C     .  miat,tiat,invg,shtvg)
C      call shoist(1,miat,nbas,s_lat%ag,s_lat%symgr,nsgrp,0)
C      stop

      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      nclasp = nclass
      nclspp = nclass
      nbasp  = nbas + npadl + npadr

C ... (pgf) Extra classes for padding layers
      lpgf = s_ctrl%lpgf(1)
      if (lpgf > 0 .and. lpgf < 10) then
C   ... Remake ipc,ics for padding and double-padding layers
        call pgfpp(nbas,npadl,npadr,1,1,.true.,.false.,s_ctrl%ipc,s_ctrl%ipc)
        call pvsym1(nbas,npadl+npadr,nclasp,s_ctrl%ipc,s_ctrl%ics)
        nclspp = nclasp
        call pvsym1(nbasp,npadl+npadr,nclspp,s_ctrl%ipc,s_ctrl%ics)
        call info2(30,0,0,'%9f%i classes for embedding region + %i for padding layers',
     .    nclass,nclasp-nclass)

C   ... Separate ipc, ipcd (dble padded ic eqiv =sngl pad vs not equiv)
        call ptr_ctrl(s_ctrl,4+1,'ipcp',nsite,0,0,s_ctrl%ipc)
        call icopy(npadl+npadr,s_ctrl%ipc(nbas+1),1,s_ctrl%ipc(nbasp+1),1)
      else
C   ... Reallocate arrays as permanent arrays
        s_ctrl%ipcp => s_ctrl%ipc
      endif

C ... Find the number of DLM classes
      ldlm = s_ctrl%ldlm
      nccomp = 0
      ndlmcl = 0
      call ptr_ctrl(s_ctrl,1,'idcc',nclspp,0,0,xx)
      do  ic  = 1, nclspp
        s_ctrl%idcc(ic) = ic
      enddo
      if (ldlm > 0) then
        call ptr_ctrl(s_ctrl,8+1,'ncomp',nclspp,0,0,xx)
        call upd_nthet(nspec,s_spec)
        call dlmidx(nclspp,s_spec,s_ctrl%ics,nccomp,ndlmcl,s_ctrl%idcc,s_ctrl%ncomp,xx)
      endif
      call ptr_ctrl(s_ctrl,2,'ics',nclspp+nccomp,0,0,xx) ! Reallocate ics
      s_ctrl%nccomp = nccomp
      s_ctrl%ndlmcl = ndlmcl

C ... Remake istab
!       if (mod(mode,10) == 2) then
!         call ptr_lat(s_lat,1,'istab',nbas,npgrp+1,0,xx)
!       else
!         call ptr_lat(s_lat,1,'istab',nbas,nsgrp+1,0,xx)
!       endif
      if (associated(s_lat % istab)) deallocate(s_lat % istab)
      if (mod(mode,10) == 2) then
        allocate(s_lat % istab(nbas,npgrp+1))
      else
        allocate(s_lat % istab(nbas,nsgrp+1))
      endif
      call dinv33(plat,1,qlat,xx)
      call symtbl(1,fptol,nbas,xx,s_lat%pos,s_lat%symgr,s_lat%ag,nsgrp,qlat,s_lat%istab,xx)
C     Append AFM grp op
      nafmg = 0
      if (afsg /= ' ') then
        call info2(30,0,0,'%9fspecial AFM symmetry: '//trim(afsg)//
     .    '%?#n# (map spin 2 from spin 1)##',lafms,2)
        call psymop(afsg,plat,s_lat%symgr(1,npgrp+1),s_lat%ag(1,npgrp+1),ic)
C       call asymop(s_lat%symgr(1,npgrp+1),s_lat%ag(1,npgrp+1),':',afsg)
        call symtbl(201,fptol,nbas,xx,s_lat%pos,s_lat%symgr(1,npgrp+1),s_lat%ag(1,npgrp+1),1,qlat,
     .    s_lat%istab(1,npgrp+1),[xx])
        nafmg = 1
C       if (abs(dlength(9,s_lat%symgr(1,npgrp+1),1)**2-3) > 1d-10)
        if (lafms .and. abs(dsum(3,s_lat%symgr(1,npgrp+1),4)-3) > 1d-10)
     .    call rx('AFMX available only with unit rotations for now, sorry')
      endif

C ... For each i*grpop fix translation ag to be that of grpop
      if (mod(mode,10) == 2) then
        if (aginv(s_lat%symgr,nsgrp,npgrp,nbas,s_lat%ag,s_lat%istab) /= 0) call rx('bug in mksym')
      endif

C ... Symops summary printout
      if (iprint() >= 50) then
         if (mod(mode/100,10) > 0) then
           call shosym(1,nbas,nsgrp,0,plat,p_sym)
           call shoist(0,s_lat%istab,nbas,s_lat%ag,s_lat%symgr,nsgrp,nafmg)
         else
          j = 0; if (iprint() >= 55) j = 1
          call shoist(j,s_lat%istab,nbas,s_lat%ag,s_lat%symgr,npgrp,nafmg)
        endif
      endif

C ... Poke s_ctrl%ipc into ssite
      call sitepack(s_site,1,nbasp+npadl+npadr,'-class',1,s_ctrl%ipcp,xx)

C --- Create s_ctrl%nrc and class labels ---
C     If this branch skipped, s_ctrl%nrc, s_ctrl%dclabl not allocated
      if (mod(mode/10,10) == 2 .or. mod(mode/10,10) == 4) then
        call ptr_ctrl(s_ctrl,1,'dclabl',nclspp+nccomp,0,0,xx)
        call ptr_ctrl(s_ctrl,1,'nrc',nclspp+nccomp,0,0,xx)
        call pvsym2(7,2*nbasp-nbas,nclspp,s_ctrl%ics,s_ctrl%ipcp,
     .    nspec,slabl,s_site,s_ctrl%dclabl,s_ctrl%nrc)
      endif

C --- Set up indexing arrays for DLM
      if (ldlm > 0) then
        call iinit(s_ctrl%ncomp,nclspp)
        call ptr_ctrl(s_ctrl,2,'idcc',nclspp+nccomp,0,0,xx)
        call dlmidx(nclspp,s_spec,s_ctrl%ics,nccomp,ndlmcl,s_ctrl%idcc,s_ctrl%ncomp,s_ctrl%nrc)

C  ...  Make DLM class labels
        call dlmclbl(nclass,nccomp,s_spec,s_ctrl%ics,s_ctrl%dclabl)

C  ...  Record offsets to first DLM class for each site (ic if not DLM)
        do  j = 1, nsite
          ic = s_ctrl%ipc(j)
          s_site(j)%dlmcl = s_ctrl%idcc(ic)
        enddo
      else
        forall (j = 1:nsite) s_site(j)%dlmcl = s_site(j)%class
      endif

C     If this branch skipped, s_ctrl%clabl not allocated
      if (mod(mode/10,10) == 2 .or. mod(mode/10,10) == 4) then
        if (associated(s_ctrl%clabl)) deallocate(s_ctrl%clabl)
        allocate(s_ctrl%clabl(nclspp+nccomp))
        do  ic = 1, (nclspp+nccomp)
          call r8tos8(s_ctrl%dclabl(ic),s_ctrl%clabl(ic))
        enddo
      endif

C      call awrit2('%n:1i',' ',100,6,nsite,s_ctrl%ipc)
C      call awrit2('%n:1i',' ',100,6,nsite,s_ctrl%ipcp)
C      call awrit2('%n:1i',' ',100,6,nclspp+nccomp,s_ctrl%nrc)
C      stop

C --- Cleanup: poke class and symmetry info into structures ---
      s_ctrl%nclass = nclass
      s_ctrl%nclasp = nclasp

C --- Re-entry point when skipping class assignment block
  100 continue
      s_lat%npgrp = npgrp
      s_lat%nsgrp = nsgrp
      j = 0; if (nafmg /= 0) j = npgrp+nafmg
      s_lat%nsafm = j; if (lafms) s_lat%nsafm = -j
      call ptr_lat(s_lat,2,'ag',3,npgrp+nafmg,0,xx)      ! Rellocate with exact size
      call ptr_lat(s_lat,2,'symgr',9,npgrp+nafmg,0,xx)   ! Rellocate with exact size
      end
      subroutine pvsym1(nbas,npad,nclass,ipc,ics)
C- Fix classes for padding layers
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas   number of basis atoms without padding layers
Ci   nclass number of classes without padding layers
Ci Outputs
Co   nclass updated
Co   ipc(nbas..nbas+npad) updated
Co   ics updated
C ----------------------------------------------------------------------
      implicit none
      integer nbas,npad,nclass,ipc(nbas+2*npad),ics(*)
      integer i,nc,is,ic,j
      if (npad == 0) return

C ... ipc->-ipc to flag which classes have been renumbered
      do  i = nbas+1, nbas+npad
        ipc(i) = -ipc(i)
      enddo

C ... For each site in padding layer, assign class
      nc = nclass
C ... Look for next occurrence of class not assigned
      do  i = nbas+1, nbas+npad

        if (ipc(i) > 0) cycle
C   ... Class ipc(i) need to be assigned:  do for each site in class
C       ic is (-)original class, is is species of original class
        ic = ipc(i)
        is = ics(-ic)

        nc = nc+1
          do  j = i, nbas+npad
            if (ipc(j) == ic) ipc(j) = nc
          enddo

C   ... Species index is same as that of original class
        ics(nc) = is

      enddo

      nclass = nc
C      call awrit2('%n:1i',' ',100,6,nbas+npad,ipc)
C      call awrit2('%n:1i',' ',100,6,nclass,ics)
C      pause

      end
      subroutine pvsym2(mode,nbas,nclass,ics,ipc,nspec,slabl,s_site,dclabl,nrc)
C- Create class labels from species labels (double precision format)
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: *
Co     Stored:    clabel
Cio    Passed to: *
Ci Inputs:
Ci   mode  : 0 do nothing
Ci           1 create class labels, store in s_site%clabel
Ci           2 create class labels, store labels in dclabl
Ci           4 create number of sites in each class nrc
Ci   nbas  : size of basis
Ci   nclass: number of classes
Ci   ics   : species table: class ic belongs to species ics(ic)
Ci   ipc   : class index: site ib belongs to class ipc(ib) (mksym.f)
Ci         : Not used unless 1s bit mode is set
Ci   nspec : number of species
Ci   slabl : vector of species labels
Co Outputs:
Co   dclabl: class labels in double precision format, if 2s bit mode is set
Co      nrc: number of sites in each class, created if 4s bit mode is set
Cu Updates
Cu   05 Apr 18 Fixed bug and changed convention for mode a little
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Dec 01 Packs class label into ssite->clabel
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nclass,nspec,ics(nclass),ipc(nbas),nrc(nclass)
      double precision dclabl(nclass)
      character*8 slabl(nspec)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Dynamically allocated arrays
      character(len=8), allocatable :: clabl(:)
C ... Local parameters
      integer ic,iclbsj,idx,is,ib
      character*8 clabll

C --- Make class labels from species labels ---
      if (mod(mode,4) > 0) then
        allocate(clabl(nclass))

C   ... Initial class labels for all sites
        do  is = 1, nspec
C         Loop should be big enough to encompass all occurences of species
C         idx = idx'th repetition of this spec
C         ic =  index to class pointing to idx'th repetition of spec

          do  idx = 1, max(nbas,nspec)
            ic = iclbsj(is,ics,-nclass,idx)
            if (ic <= 0) exit ! no such occurrence
            call clabel0(slabl,is,idx,clabl(ic))

C           Check for collisions of new class with other species
C            do  js = 1, nspec
C              if (js == is) cycle
C              if (clabl(ic) == slabl(js)) then
C                call rxi('mksym: class label '//trim(clabl(ic))//
C     .            ' collides with species',js)
C              endif
C            enddo
          enddo
        enddo

C   ... Render class labels unique
        do  ic = 2, nclass
          clabll = clabl(ic)
          call uniqclabl(0,ic-1,0,[0],1,clabl,clabll)
C         Debugging
C         if (clabl(ic) /= clabll) then
C           print *, 'class',ic,' replace label',clabl(ic),' with',clabll
C         endif
          clabl(ic) = clabll
        enddo

C   ... Copy to dclabl
        if (mod(mode/2,2) == 1) then
          do  ic = 1, nclass
            call s8tor8(clabl(ic),dclabl(ic))
          enddo
        endif

C   ... Poke into s_site
        if (mod(mode,2) == 1) then
          do  ib = 1, nbas
            ic = ipc(ib)
            s_site(ib)%clabel = clabl(ic)
          enddo
        endif

        deallocate(clabl)
      endif

C --- Create nrc ---
      if (mod(mode/4,2) == 1) then
        call iinit(nrc,nclass)
        do  ib = 1, nbas
        ic = ipc(ib)
        nrc(ic) = nrc(ic)+1
        enddo
      endif

C ... debugging : check class labels
C      call awrit2('%n:1i',' ',100,6,nclass,nrc)
C      do  ic = 1, nclass
C        call r8tos8(dclabl(ic),clabll)
C        print 333, ic, nrc(ic), clabll
C  333   format(i4,i4,1x,a)
C   50 continue
C      enddo
C      stop

      end

      integer function aginv(g,nsg,npg,nbas,ag,istab)
C- Finds translation part for extra inversion-generated group ops
C ----------------------------------------------------------------------
Ci Inputs
Ci   g     :point group operations
Ci   nsg   :number of space group operations
Ci   npg   :number of point group operations
ci   nbas  :dimensions istab.  If zero, istab is not altered.
Co Outputs
Co   ag    :for each group op ig in range npg+1..nsg,
Co         :the corresponding group op jg is found such that
Co         :g(jg) = -g(ig), and ag(jg) is copied to ag(ig); see Remarks
Co   istab :istab(jg) is copied into istab(ig) just in the same way as
Co         :ag(jg) is copied into ag(ig).
Co         :The sign of istab(1) is set negative to indicate that this
Co         :operation is derived from a true space group operation.
Cr Remarks
Cr   In certain cases the inversion operation may be added to the space
Cr   group, for purposes of k integration.  This is permissible when the
Cr   hamiltonian has the form h(-k) = h*(k).  In that case, the
Cr   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
Cr
Cr   Also, the Green's functions are related G(-k) = Gtranspose(k).
Cr   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
Cr   then G(-k) = Gtranspose(k), and the same (g,ag) information is
Cr   needed for either rotation.  Let g~ = -g:
Cr     G0(g~^-1 -k) into G(k),
Cr
Cr
Cr   Routine aborts with aginv=ig if for group ig it cannot find
Cr   corresponding jg such that g(jg) = -g(ig)
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nsg,npg,nbas,istab(nbas,*)
      double precision g(3,3,npg),ag(3,npg)
C ... Local parameters
      integer ig,jg
      double precision ginv(3,3)
      logical grpeql

      do  ig = nsg+1, npg
        call dpcopy(g(1,1,ig),ginv,1,9,-1d0)
        do  jg = 1, npg
          if (grpeql(g(1,1,jg),ginv)) then
            call dcopy(3,ag(1,jg),1,ag(1,ig),1)
            if (nbas > 0) then
              call icopy(nbas,istab(1,jg),1,istab(1,ig),1)
              istab(1,ig) = -istab(1,jg)
            endif
C           call awrit2(' aginv mapping %,2i ->%,2i',' ',80,6,ig,jg)
            goto 10
          endif
        enddo
        aginv = ig
        return
   10   continue
      enddo
      aginv = 0
      end

