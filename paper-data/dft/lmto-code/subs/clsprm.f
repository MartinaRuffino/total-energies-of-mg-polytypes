C#define unix
      subroutine clsprm(mode,s_ctrl,s_ham,s_pot,s_spec,s_lat,s_bz,s_strn)
C- Allocate memory for, and try to read, class-related data
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: nbas nclass nl nspin nspec ldlm nccomp lncol loptc
Ci                lscr nmap lves lstonr zbak
Co     Stored:    zbak
Co     Allocated: *
Cio    Elts Passed:lsx lasa
Cio    Passed to: ioqpp aiocls clsprp
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read: nlibu neula nbf
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: aiocls clsprp
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read: rhrmx
Co     Stored:    nrhos rhrmx
Co     Allocated: pnu qnu pp pprel qc qt ves vdif rhrmx vrmax dlmwt
Co                thet gibbs mxy bxc sop rhos grrme vintr pmpol qpp GFr
Cio    Elts Passed:qpp pnu qnu pp pprel rhrmx vrmax ves vdif vintr pmpol
Cio               sop grrme GFr
Cio    Passed to: ioqpp aiocls clsprp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: idmod z lmxa rmt p q
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: aiocls clsprp
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: plat vol
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:pos
Cio    Passed to: aiocls clsprp
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: *
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: clsprp
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read: strn
Co     Stored:    strn
Co     Allocated: *
Cio    Elts Passed:strn
Cio    Passed to: str_pack
Ci Inputs
Ci    mode :1s digit
Ci         : 1 allocates ASA arrays (see spot,sham,sarray)
Cu Updates
Cu   07 Aug 16 Dirac qnur is always allocated
Cu   27 Jun 14 Some adjustments for relativistic (Alena)
Cu   14 Nov 13 Some adjustments in preparation for relativistic CPA
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   01 Sep 11 Begin migration to f90 structures
Cu   05 May 08 Changed the allocation of DLM arrays (Kirill):
Cu             Added storage for each DLM angle at the end of all arrays
Cu             nccomp is the total number of DLM angles, stored in sctrl
Cu   08 Nov 07 (J. Xu) LDA+U implementation; qpp is complex
Cu   02 Aug 07 (Derek Stewart) Added declaration of variables
Cu             ppd, pnud, qnud, for DLM sites
Cu   30 Sep 04 Reads/writes relativistic ppar's
Cu   02 Apr 04 SO array redimensioned; make and pack pot->nrhos
Cu   18 Mar 03 Allocate space for relativistic potential parameters
Cu   07 Feb 03 SO array redimensioned
Cu   19 Feb 02 entry clsprp copies class data to doubly padded layers
Cu   28 Apr 98 allocate space for radial matrix elements of grad
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode
C ... For structures
!      include 'structures.h'
      type(str_bz):: s_bz
      type(str_ctrl):: s_ctrl
      type(str_lat):: s_lat
      type(str_ham):: s_ham
      type(str_pot):: s_pot
      type(str_spec):: s_spec(*)
C     real(8),pointer :: p_rmax(:)
      type(str_strn) :: s_strn(*)
C ... Local parameters
      integer,parameter :: n0=10, NULLI=-99999
      logical cmdopt
      character*80 strn,clabl*8,maprul*(2048)
      integer i1,ib,ic,iclbsj,idmod(n0),iprint,is,lmx,lnsph,nbas,nbf,
     .  nclasp,nclass,nclspd,nclspp,nclsppd,neula,nl,nlibu,nlspc,
     .  nelrel,nmap,nqpp,nrhos,nsp,nspc,nspec,str_pack
      double precision plat(3,3),z,xx
C     DLM variables
      integer ldlm,nccomp
      procedure(integer) :: nglob

      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      nspc = nglob('nspc')      ! 2 for noncollinear case
      nspec = s_ctrl%nspec
      nlibu = s_ham%nlibu
      nclasp = s_ctrl%nclasp
      nclspp = 2*nclasp-nclass
      ldlm = s_ctrl%ldlm

      nccomp = s_ctrl%nccomp
C     Reserve space for additional DLM classes
      if (ldlm == 0) then
        nccomp = 0
        nclspd  = max(nclasp,nspec)
        nclsppd = max(nclspp,nspec)
      else
        nclspd  = nclasp+nccomp
        nclsppd = nclspp+nccomp
      endif
      nlspc = nl*nsp*nclsppd
      nelrel = 8*nl*nl          !*nclsppd

C ... Allocate ASA arrays
      if (mod(mode,10) == 1) then

        call ptr_pot(s_pot,8+1,'pnu',nl*nsp,max(nclsppd,nspec),xx)
        call ptr_pot(s_pot,8+1,'qnu',3*nl*nsp,max(nclsppd,nspec),xx)
        call ptr_pot(s_pot,8+1,'qnur',4*nelrel,max(nclsppd,nspec),xx)
        s_pot%qnur(1,1) = NULLI ! => array has no content yet
        call ptr_pot(s_pot,8+1,'pp',6*nl*nsp,max(nclsppd,nspec),xx)
        call ptr_pot(s_pot,8+1,'pprel',5*nelrel,max(nclsppd,nspec),xx)
        s_pot%pprel = NULLI
        call ptr_pot(s_pot,8+1,'qc',nclspd,0,xx)
        call ptr_pot(s_pot,8+1,'qt',nclspd,0,xx)
        call ptr_pot(s_pot,8+1,'ves',nclsppd,0,xx)
        call ptr_pot(s_pot,8+1,'vdif',nclsppd,0,xx)
        call ptr_pot(s_pot,8+1,'rhrmx',nclsppd,0,xx)
        call ptr_pot(s_pot,8+1,'vrmax',2,nclsppd,xx)

        call ptr_ctrl(s_ctrl,8+1,'initc',nclspd,0,0,xx)
        call ptr_ctrl(s_ctrl,8+1,'group',2*nclspd,0,0,xx)
        call ptr_pot(s_pot,8+1,'shfac',3,nclsppd,xx)
        call ptr_pot(s_pot,8+1,'gfr',4*4*4*nelrel,nclsppd,xx)

C   ... Allocate and set up DLM angles and weights
C       if (ldlm > 0) then
          call ptr_pot(s_pot,8+1,'dlmwt',nclspd,0,xx)
          call ptr_pot(s_pot,8+1,'thet',max(1,nccomp),2,xx)
          call ptr_pot(s_pot,8+1,'gibbs',nclspd,0,xx)
          call ptr_pot(s_pot,8+1,'mxy',nl*2*nclsppd,0,xx)
C       endif
        call ptr_ctrl(s_ctrl,8+1,'rmax',nclsppd,0,0,xx)
        call spec2class(s_spec,nclspp+nccomp,s_ctrl%ics,'rmt',1,xx,
     .    s_ctrl%rmax)

C       Even in collinear case, bxc might be written to
        call ptr_pot(s_pot,8+1,'bxc',3*nclspd,0,xx)
        s_pot%nrhos = 0
        if (s_ctrl%lncol /= 0) then
          neula = s_ham%neula
          nbf = s_ham%nbf
          nrhos = max(nl,neula,nbf)
!         if (mod(s_ctrl%lrel,10) == 2) nrhos = nl*nl
          call ptr_pot(s_pot,8+1,'sop',nl*nsp*nsp*9,nclsppd,xx)
          call ptr_pot(s_pot,8+1,'rhos',2*3*4*nrhos*nclspd,0,xx)
          s_pot%nrhos = nrhos
        else
          allocate(s_pot%sop(1,1))
          allocate(s_pot%rhos(1))
        endif
        if (s_ctrl%loptc > 0)
     .    call ptr_pot(s_pot,8+1,'grrme',4*2*nl*nsp*nclasp,0,xx)
        if (IAND(s_ctrl%lsx,3) /= 0 .or. s_ctrl%lscr /= 0) then
          call ptr_pot(s_pot,8+1,'vintr',nclasp*(nl*nsp)**2,0,xx)
        endif
        lnsph = (IAND(s_ctrl%lasa,32))

C   ... Allocate memory for, attempt to read qpp from file
        if (lnsph /= 0 .or. nlibu /= 0) then
          i1 = nl**2
          nqpp = (i1*(i1+1))/2
          call ptr_pot(s_pot,8+1,'pmpol',
     .      (2*nl-1)*nl**2*3*nsp*nclasp,0,xx)
!           call ptr_pot(s_pot,1,'qpp',2*nqpp*4*nsp*nbas,0,xx)
!           call dvset(s_pot%qpp,1,1,-1d0)
          if (associated(s_pot%qpp)) deallocate(s_pot%qpp)
          allocate(s_pot%qpp(2*nqpp*4*nsp*nbas))
          s_pot%qpp = 0
          s_pot%qpp(1) = -1
          call ioqpp(.false.,s_ctrl,s_pot)
        else
!           call ptr_pot(s_pot,1,'qpp',1,0,xx)
          if (.not. associated(s_pot%qpp)) allocate(s_pot%qpp(1))
        endif

        plat = s_lat%plat

C   --- Create maps ---
        nmap = s_ctrl%nmap
        if (nmap > 0) then
C          call rx('clsprm strn map')
          i1 = str_pack('map',-2,s_strn,maprul)
          do  ic = 1, nclass
            xx = s_ctrl%dclabl(ic)
            call r8tos8(xx,clabl)
            is = s_ctrl%ics(ic)
            ib = iclbsj(ic,s_ctrl%ipc,-nbas,1)
            idmod(1:n0) = s_spec(is)%idmod(1:n0)
            if (ib < 0) then
              if (iprint() >= 20) print *,
     .          'CLSPRM: no basis corresponding to class ' //clabl
              cycle
            endif
            z = s_spec(is)%z
            lmx = s_spec(is)%lmxa
            call wrmap(maprul,nmap,clabl,ic,ib,z,plat,s_lat%pos)
            call mapq(clabl,ic-1,nl,nsp,lmx,idmod,s_pot%pnu,s_pot%qnu,
     .        s_pot%pp,s_pot%pprel,s_pot%rhrmx,s_pot%vrmax,
     .        IAND(s_ctrl%lves,1)==1,s_pot%ves,s_ctrl%initc)
          enddo
          if (cmdopt('-spin1',6,0,strn))
     .      call spinav(0,nclass,nl,nsp,s_pot%pnu,s_pot%qnu)
          if (cmdopt('-map0',5,0,strn)) then
            call aiocls(.true.,0,s_ctrl,s_ham,s_pot,s_spec,
     .        s_lat,1,nclass)
            call rx0('mapped moments written')
          endif
        endif
      endif

      if (s_ctrl%lstonr(1) /= 0 .or. .true.) then
      endif

      return

      entry clsprp(mode,s_ctrl,s_ham,s_pot,s_spec,s_lat,s_bz)
C- Copy class data to doubly padded case

      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      nspec = s_ctrl%nspec
      nclasp = s_ctrl%nclasp
      if (nclasp == nclass) return
      nclspp = max(2*nclasp-nclass,nspec)
C ... Allocate ASA arrays
      if (mod(mode,10) == 1) then

        i1 = nclasp-nclass
        call dpscop(s_pot%pnu,s_pot%pnu,nl*nsp*i1,
     .    nl*nsp*nclass+1,nl*nsp*nclasp+1,1d0)
        call dpscop(s_pot%qnu,s_pot%qnu,3*nl*nsp*i1,
     .    3*nl*nsp*nclass+1,3*nl*nsp*nclasp+1,1d0)
        call dpscop(s_pot%qnur,s_pot%qnur,32*nl*nl*i1,
     .    32*nl*nl*nclass+1,32*nl*nl*nclasp+1,1d0)

        call dpscop(s_pot%ves,s_pot%ves,i1,nclass+1,nclasp+1,1d0)
        call dpscop(s_pot%vdif,s_pot%vdif,i1,nclass+1,nclasp+1,1d0)
        call dpscop(s_pot%vrmax,s_pot%vrmax,i1,nclass+1,nclasp+1,1d0)
        call dpscop(s_pot%rhrmx,s_pot%rhrmx,i1,nclass+1,nclasp+1,1d0)

      else
        call rxi('clsprp: bad mode',mode)
      endif

      end
      subroutine mapq(clabl,ic,nl,nsp,lmx,idmod,pnu,qnu,pp,pprel,rhrmx,vrmax,swmad,ves,initc)
C- Maps moments and pp's from MAP
C ----------------------------------------------------------------------
Ci Inputs:
Co Outputs:
Co   pnu,qnu,pp
Co   initc is updated to reflect changed status
Cr Remarks
Cr   Map has syntax
Cr     name.ext or path/name.ext (one file)
Cr     nfiles: c1 c2 ... cn name1 name1 ... namen (multiple files)
Cr   The 10's digit of c* is interpreted as a switch, to swap
Cr     the up and down spins
Cr   Note: '$' for path is translated into string LIBPTH, below
Cr         '.' uses value of current data
Cr   idmod,rhrmx,vrmax,ves set by 1st mapped file
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      character*8 clabl
      integer ic,nl,nsp,lmx,idmod(nl),initc(0:*)
      double precision pnu(nl,nsp,0:*),qnu(3,nl,nsp,0:*),ves(0:*),
     .  pp(6,nl,nsp,0:*),pprel(5,nl,2*nl,2,2,ic),rhrmx(0:*),vrmax(2,0:*)
      logical swmad
C Local variables
      real(8),parameter :: NULLR=-99999d0
      integer a2vec,idwk(6),fext,fextg,fopna,fpath,i,i1mach,lrel,
     .  ia,ia2,iend,ifi,imap,ip,iterm,ix(10),j,nmap,ipr,initx(2)
      logical scat,a2bin,aiomom,aiopar,last,map2me
      double precision fcof(20),pwk(6*2),qwk(3*6*2),ppwk(6*6*2),
     .  pprwk(32*36*2),rhwk,pout(6*2),qout(3*6*2),ppout(6*6*2),
     .  vrwk(2),pprout(32*36*2),veswk,cof,scof,bxc(3)
      character*1 s(0:71),ch(2),chp(4),ss*72,fnam*40,alabl*8
      equivalence (s,ss)
      procedure(integer) :: nglob
C#ifdef unix
      character*13 LIBPTH
      data LIBPTH /'/usr/lmto/lib'/
C#endif
      data ch /' ','.'/, chp /'$','~','/','.'/


      call getpr(ipr)
      lrel = mod(nglob('lrel'),10)

C --- Open file, exit if no map ---
      ifi = fopna(clabl,30,0)
      if (.not. scat(ifi,'MAP:',':',.true.)) then
        call fclose(ifi)
        return
      endif

C --- Zero out accumulated moms, pp's ---
      initx(1) = 1
      initx(2) = 1
      call dpzero(pwk,  6*2)
      call dpzero(qwk, 3*6*2)
      call dpzero(ppwk,6*6*2)
      call dpzero(pprwk,32*36*2)
      call dpzero(pout, nl*nsp)
      call dpzero(qout, nl*nsp*3)
      call dpzero(ppout,nl*nsp*6)
      call dpzero(pprout,32*36*2)
      rhwk = 0
      veswk = 0
      vrwk(1) = 0
      vrwk(2) = 0
      scof = 0


C --- Parse the file map and admix coefficients ---
      read (ifi,1) s
    1 format(72A1)
      call fclose(ifi)
      ia = 0
      call skipbl(ss,len(ss),ia)
      call skpblb(ss,len(ss),i)
      if (ipr >= 30) print 2,clabl,ss(1:i+1)
    2 format(' MAPQ: mapping class ',a4,' as ',a)
      ia2 = ia
      fcof(1) = 1
      if (a2bin(ss,nmap,2,0,':',ia,len(ss))) then
        if (a2vec(ss,len(ss),ia,4,' ',1,2,nmap,ix,fcof) /= nmap)
     .    call fexit(-1,1,' Exit -1 mapq, class ' // clabl //
     .   'syntax error seeking %i coffs in map '//ss(1:ia)//' ...',nmap)
      else
        nmap = 1
        ia = ia2
      endif

C --- Loop over all mapping files ---
      do  imap = 1, nmap

C ---   Extract file names: full name in s(ia:iend) ---
        call skipbl(s,72,ia)
        iend = ia
        call chrpos(s,' ',72,iend)
        map2me = iend == ia+1 .and. s(ia) == '.'

C ---   Case map is '.': mapped data is current data ---
        if (map2me) then
          if (mod(initc(ic),2) == 1) then
            call dcopy(nl*nsp,pnu(1,1,ic),1,pwk,1)
            call dcopy(3*nl*nsp,qnu(1,1,1,ic),1,qwk,1)
          else
            if (ipr >= 20) print *, 'mapq: file missing moments'
            initx(1) = 0
          endif
          if (mod(initc(ic)/2,2) == 1) then
            call dcopy(6*nl*nsp,pp(1,1,1,ic),1,ppwk,1)
            call dcopy(40*nl*nl,pprel(1,1,1,1,1,ic),1,pprwk,1)
          else
            if (ipr >= 20) print *, 'mapq: file missing ppars'
            initx(2) = 0
          endif
          do  i = 1, nl
            idwk(i) = idmod(i)
          enddo
          vrwk(1)   = vrmax(1,ic)
          vrwk(nsp) = vrmax(nsp,ic)
          rhwk      = rhrmx(ic)
          veswk     = ves(ic)
C --- Case map to another file : extract path name ---
        else
C ...     Path name ...
          fnam = ' '
          ia2 = ia
          i = ia
          call chrps2(s,chp,4,iend,ia2,iterm)
          if (iterm == 1) then
            fnam = LIBPTH
          ia = ia+1
          elseif (iterm == 2 .or. iterm == 3) then
          elseif (iterm == 4) then
            if (ia /= ia2) goto 6
          else
            goto 6
          endif
C ...     Extract path name, looping until last '/' found
    4     continue
          ip = i
          call chrpos(s,'/',iend,i)
          i = i+1
C print *, i, ip, iend, ss(ia+1:i), ' ', ss(ia+1:ip)
          if (i < iend) goto 4
          i = 0
          call strcat(fnam,32,' ',ss(ia+1:ip),ip-ia,' ',i)
          if (iterm == 2) call rxx(ia /= ia2,'mapq: bad path name')
          i = fpath(fnam)
          ia = ip

C --- Open file with local name, ext; then restore them ---
    6     continue
          ia2 = ia
          call chrps2(s,ch,2,iend,ia2,iterm)
C  print *, ia, ia2, iend, iterm, ss(ia:ia2), ' .. ', ss(ia2:iend)
          i = fextg(fnam)
          if (iterm == 2) i = fext(ss(ia2+1:iend))
C  print *, 'Opening file ',ss(ia:ia2),', coff=',sngl(fcof(imap))
          ifi = fopna(ss(ia+1:ia2),31,1)
          i = fext(fnam(1:20))
          i = fpath(' ')

C --- Read in moments to wk ---
          alabl = clabl
          rewind ifi
          if ( .not. (aiomom(alabl,pwk,qwk,idwk,nl,lmx,nsp,rhwk,vrwk,bxc
     .        ,ifi))) then
            initx(1) = 0
            if (ipr >= 20) print *,
     .        'mapq: mapped file '//ss(ia+1:ia2)//' missing moments'
          endif

          if ( .not. (aiopar(alabl,lrel,ppwk,pprwk,veswk,nl,lmx,nsp,ifi)
     .        )) then
            initx(2) = 0
            if (ipr >= 20) print *,
     .        'mapq: mapped file '//ss(ia+1:ia2)//' missing ppars'
          endif
          call fclose(ifi)
        endif

C --- Now wk arrays hold current mapped data ---
        if (imap == 1) then
          do  i = 1, nl
            idmod(i) = idwk(i)
          enddo
          vrmax(1,ic) = vrwk(1)
          vrmax(nsp,ic) = vrwk(nsp)
          if (dabs(fcof(1)) >= 10) then
            vrmax(1,ic) = vrwk(nsp)
            vrmax(nsp,ic) = vrwk(1)
          endif
          rhrmx(ic) = rhwk
          if (.not. swmad) ves(ic)   = veswk
        endif
C ...   Swap up and down spins
        cof = dsign(dmod(dabs(fcof(imap)),10d0),fcof(imap))
        scof = scof + cof
        if (dabs(fcof(imap)) >= 10 .and. nsp == 2) then
          call daxpy(nl,cof,pwk(1+nl),1,pout,1)
          call daxpy(nl,cof,pwk,1,pout(1+nl),1)
          call daxpy(3*nl,cof,qwk(1+3*nl),1,qout,1)
          call daxpy(3*nl,cof,qwk,1,qout(1+3*nl),1)
          call daxpy(6*nl,cof,ppwk(1+6*nl),1,ppout,1)
          call daxpy(6*nl,cof,ppwk,1,ppout(1+6*nl),1)
          if (lrel == 2) call rx('mapq: spin swap not implemented')
        else
          call daxpy(nl*nsp,cof,pwk,1,pout,1)
          call daxpy(3*nl*nsp,cof,qwk,1,qout,1)
          call daxpy(6*nl*nsp,cof,ppwk,1,ppout,1)
          call daxpy(32*nl**2,cof,pprwk,1,pprout,1)
        endif

        ia = iend+1
      enddo

C --- Copy mapped data back to p's, q's and ppar's ---
      call dcopy(nl*nsp,pout,1,pnu(1,1,ic),1)
      call dcopy(nl*nsp*3,qout,1,qnu(1,1,1,ic),1)
      call dcopy(nl*nsp*6,ppout,1,pp(1,1,1,ic),1)
      call dcopy(5*nl*2*nl*2*2,pprout,1,pprel(1,1,1,1,1,ic),1)
      initc(ic) = 4*(initc(ic)/4) + initx(1) + 2*initx(2)

C --- Printout ---
      cof = 0
      do  i = 1, nsp
      vrwk(i) = 0
        do  j = 1, nl
          vrwk(i) = vrwk(i)+qnu(1,j,i,ic)
        enddo
      enddo
      ss =
     . ' mapq: sum of coffs= %d   valence charge= %1;6d   moment= %1;6d'
      if (ipr >= 30 .and. nsp == 1)
     .  call awrit2(ss(1:47),' ',80,i1mach(2),scof,vrwk)
      if (ipr >= 30 .and. nsp == 2) call awrit3(ss(1:63),
     .  ' ',80,i1mach(2),scof,vrwk(1)+vrwk(2),vrwk(2)-vrwk(1))
      if (ipr <= 40) return
      if (mod(initc(ic),2) == 1)
     .  last = aiomom(clabl,pnu(1,1,ic),qnu(1,1,1,ic),[NULLR],
     .  idmod(1),nl,lmx,nsp,0d0,rhrmx(ic),vrmax(1,ic),(/0d0,0d0,1d0/),
     .  -i1mach(2))
      if (mod(initc(ic)/2,2) == 1)
     .  last = aiopar(clabl,lrel,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),
     .                ves(ic),nl,lmx,nsp,-i1mach(2))


      end
      subroutine wrmap(maprul,nmap,clabl,ic,ib,z,planvc,bas)
C- Appends map to atom file according to maprul
C ----------------------------------------------------------------------
Ci Inputs:
Ci   maprul,ic, ific:file handle
Ci   ib:     any basis atom corresponding to class
Co Outputs:
Co   Map appended to atom files
Cr Remarks
Cr   Each file is tested against each map until a match is found.
Cr   If no match found, nothing is done.
Cr   Map syntax:  key:  map-pattern
Cr   where mapnam is wild-card expanded by the shell;
Cr   there is a match when the present atom file is one of the files.
Cr   When there is a match, a map is appended to the current file,
Cr   using the map pattern.  All occurences of strings '\h' and '\t'
Cr   in the map pattern are substituted by the current class name and
Cr   extension, respectively.  A backslash followed by an expresssion
Cr   in parenthesis, eg '\{job}', is evaluated as an expression and
Cr   whose value is inserted as an ascii string.  Example: atom
Cr   name=A, extension=x, variable job=31.  the map pattern
C    \h\t:g\{job+1}z is translated into A.x:g32z .
C    Example of a map:
C    MAP     F  [ac]:     \h.1:g\{job}
C               b:         c.2:g\{job}
C               *:         a.3:g\{job}
C    If extension=x and job=3, this maps classes a and c into
C    a.1:g3 and c.1:g3, class b into b.2:g3 and all other classes
C    into file a.3:g3
C    Argument list and convention differ slightly from FP version.
C    For each class, class-specific variables ic,ib,x,y,z,n are loaded
C    into the variables table before coefficients are parsed.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer ib,ic,ific,nmap,reclen
      double precision planvc(3,3),bas(3,1),z
      parameter (reclen=120)
      character maprul(1)*(reclen),clabl*8
C Local Variables
      integer imap,ich,ln,iout,i,ich0,fextg,iprint,ncof,i2,j,fopna,
     .  ix(10),ival,iv0,a2vec,i1mach
      double precision fcof(10),ddot,xx
      logical foundm,a2bin,scat
      character tstr*40,chx(2),ext*40
      character, parameter :: sch = '\'
! I don't think this excessive escaping is necessary but let's not leave loose ends.
      character(len=2), parameter :: hch = '\'//'h', tch = '\'//'t', pch = '\'//'{'
      character*1 s(0:71),ss*72,fmt*40
      equivalence (s,ss)

      data chx /':',' '/
C#ifdefC CRAY
C      data hch /'\h'/ tch /'\t'/ sch /'\'/ pch /'\{'/
! C#else
!       data hch /'\\h'/ tch /'\\t'/ sch /'\\'/ pch /'\\{'/
C#endif

C --- Until a key matches atom file, do ---
      if (iprint() >= 30) print *
      do  imap = 1, nmap
        if (maprul(imap) == ' ') cycle
        ich = 0
        call skipbl(maprul(imap),reclen,ich)
        ich0 = ich
C   ... Quick check to avoid system call
        if (maprul(imap)(ich+1:ich+3) == '*: ') then
          ich = ich+2
          foundm = .true.
          goto 2
        endif
        call chrps2(maprul(imap),chx,2,72,ich,i)
        call rxx(i /= 1,'WRMAP: bad key in map:'//maprul(imap)(1:20))
        tstr = ' '
        call strcop(tstr,clabl,8,' ',i)
        if (tstr(i:i) == ' ') i = i-1
        j = fextg(ext)
        call strcop(tstr(i+1:40),ext,40,' ',j)
        call locase(tstr)
C       call chrpos(tstr,' ',40,ich)
        ss = 'ls '// maprul(imap)(ich0+1:ich) // ext(1:j)
     .    // '|grep -s '// tstr
C       print *, ss
        call fsystm(ss,j)
        foundm = j == 0
        ich = ich+1
C ...   Regather, with foundm and ich set
    2   continue
        if (iprint() > 40) print *, 'WRMAP: sought map for class ',
     .    clabl,' in key ', maprul(imap)(ich0+1:ich),
     .    '  ... found=',foundm
        if (foundm) then
        call skipbl(maprul(imap),reclen,ich)
        if (ich >= reclen) then
          if (iprint() > 40) print *, 'WRMAP: empty map'
          return
        endif

C ---   Map found:  parse coffs and maprul substituting \h,\t,\{ ---
        iout = 0
        ss = ' '
        ich0 = ich
C ...   first look for coffs, calling a2bin to sub
        if (a2bin(maprul(imap),ncof,2,0,':',ich0,reclen)) then
          call skipbl(maprul(imap),reclen,ich0)
          ich = ich0
          write(ss(2:3),'(i1,'':'')') ncof
          fmt = ' '
          write(fmt(1:9),'('':%'',i1,'':1;6d'')') ncof
          iout = 3

          call numsyv(iv0)
          call lodsyv('ic',1,dble(ic),ival)
          call lodsyv('zc',1,z,ival)
          call lodsyv('ib',1,dble(ib),ival)
          call lodsyv('x', 1,bas(1,ib),ival)
          call lodsyv('y', 1,bas(2,ib),ival)
          call lodsyv('z', 1,bas(3,ib),ival)
          xx = ddot(3,bas(1,ib),1,planvc(1,3),1) /
     .         dsqrt(ddot(3,planvc(1,3),1,planvc(1,3),1))
          call lodsyv('n', 1,bas(3,ib),ival)
C         call shosyv(0,0,0,6)
          i2 = a2vec(maprul(imap),reclen,ich,4,',; ',3,3,-10,ix,fcof)
          j = 0
          do  i = 1, i2
            if (ix(i) /= 2) then
              j = j+1
              fcof(j) = fcof(i)
            endif
          enddo
          if (j /= ncof) then
            call awrit3(' WRMAP map %i:  sought %i coefficients but'//
     .        ' found only %i',' ',72,i1mach(2),imap,ncof,j)
            call rx('Failed to parse map'//maprul(imap)(1:ich)//' ...')
          endif
          call awrit1(fmt(1:9),ss(iout:72),72-iout+1,0,fcof)
          call skpblb(ss,72,iout)
          iout = iout+2
          call clrsyv(iv0)
        endif

    3   continue
        ich0 = ich
        call chrpos(maprul(imap),sch,reclen,ich)
        ln = ich-ich0
        call strcop(s(iout),maprul(imap)(ich0+1:reclen),ln,sch,i)
        iout = iout+i
        if (ich >= reclen) exit
        if (maprul(imap)(ich+1:ich+2) == hch) then
          call strcat(s(iout),72-iout,' ' ,clabl,8,' ',i)
          iout = iout+i
          ich = ich+2
          if (s(iout-1) == ' ') iout = iout-1
        endif
        if (maprul(imap)(ich+1:ich+2) == tch) then
          i = fextg(tstr)
          call strcat(s(iout),72-iout,' ' ,tstr,20,' ',i)
          iout = iout+i
          ich = ich+2
          if (s(iout-1) == ' ') iout = iout-1
        endif
        if (maprul(imap)(ich+1:ich+2) == pch) then
          ich = ich+2
          if (.not. a2bin(maprul(imap),fcof,4,0,'}',ich,reclen))
     .      call rx('WRMAP: bad expr in map'//maprul(imap)(1:ich))
          call bin2a(' ',0,0,fcof,4,0,72,ss,iout)
        endif
C        if (maprul(imap)(ich+1:ich+1) == sch) then
C          s(iout) = sch
C          ich = ich+1
C          iout = iout+1
C        endif
c       print *, ss
        goto 3
      endif
      enddo

      if (.not. foundm) then
        if (iprint() >= 20)
     .    print *, 'WRMAP:  no map found for class ', clabl
        return
      endif

C --- Supersede existing map, or append map to atom file ---
      ific = fopna(clabl,30,0)
      if (.not. scat(ific,'MAP:',':',.true.)) then
        call poseof(ific)
        write(ific,'(''MAP:'')')
      endif
      write(ific,'(a72)') ss
      if (iprint() >= 30) print *, 'WRMAP:  wrote map for class ',
     .    clabl,' ', ss
      call fclose(ific)

      end
