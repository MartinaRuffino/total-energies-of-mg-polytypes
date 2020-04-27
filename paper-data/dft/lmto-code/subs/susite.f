      subroutine susite(s_ctrl,s_ham,s_pot,s_lat,s_spec,s_site)
C- Sets up permanent arrays related to basis
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nsite nbas nl lordn lpgf pgfsl pgfvl
Co     Stored:     ncl pgfsl pgfvl ips npl nsite nbasp pgplp sdxsi
Co                 npadl npadr
Co     Allocated:  *
Cio    Elts passed:ips clssl cllst clp pgfsl pgfvl pgplp lncol
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     seref neula nbf
Co     Allocated:  eula magf
Cio    Elts passed:eula magf
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  vshft aamom
Cio    Elts passed:vshft
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat platl platr pos plate
Co     Stored:     pos
Co     Allocated:  pos
Cio    Elts passed:pos
Cio    Passed to:  clsset
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name eref
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pgfset
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos pl plv bfield eula
Co     Stored:     clabel pos pl plv bfield eula
Co     Allocated:  *
Cio    Elts passed:pl plv
Cio    Passed to:  siteperm pvsub2
Ci Inputs
Co Outputs
Cl Local variables
Cl   nbasmem : On input, s_site and s_lat%pos are allocated
Cl             with space for additional site data.  This routine
Cl             may create extra sites; it checks that the total
Cl             number of sites created does not exceed nbasmem,
Cl             which is inferred from the size of s_lat%pos.
Cr Remarks
Cr  *Arrays which susite sets:
Cr      s_lat%pos,s_ctrl%ips,s_ctrl%pgplp,s_ham%eula
Cr      s_ctrl%pgfsl,s_ctrl%pgfvl,pgord for PGF
Cr
Cr  *pgplp holds information about crystal subblocks, and its meaning
Cr   depends on the context.  pgplp is dimensioned pgplp(6,-1:*)
Cr   pgplp(1,-1) contains information about the context-dependence,
Cr   to permit certain routines that depend on the context to be
Cr   applied across different ones.
Cr
Cr    *For the usual crystal case, there are no subblocks.
Cr     In this case, pgplp(1,-1) = -1.  susite creates and sets pgplp.
Cr
Cr    *For the planar Green's function, subblocks are principal layers.
Cr     In this case, it two additional subblocks are added to treat
Cr     the boundaries: one to the "left" of the first PL, and one to the
Cr     "right" of the last.  Here, pgplp(1,-1) = 0.  Some rows
Cr     of pgplp are initially set by pgfset; others are set in pgfasa.
Cr     Information is kept about npl+2 layers.
Cr
Cr    *For the order-N Green's function, subblocks are the crystal
Cr     subblocks.  Here, pgplp(1,-1) = -2.
Cr
Cu Updates
Cu   13 Aug 13 s_ham%eula made into a 2D array
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   01 Sep 11 Begin migration to f90 structures
Cu   24 May 08 eula and magf read from master only, MPI mode
Cu   07 Feb 03 looks for, and reads magnetic field from file
Cu   23 Jan 02 susite packs seref into sham.  Altered argument list.
Cu   18 Dec 01 susite packs species labels into ssite
Cu   13 Oct 98 convention for pgplp extended to non-layer case
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_ham):: s_ham
      type(str_lat):: s_lat
      type(str_pot):: s_pot
      type(str_spec):: s_spec(*)
      type(str_site):: s_site(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: clord(:)
      real(8), allocatable :: wk(:)
      integer, allocatable :: pgord(:)
      type(str_site),pointer:: s_site0(:)
C ... Local parameters
      integer procid,master,mpipid,nbasmem,err
      logical pass1,lorder,mlog,cmdopt
      integer lpgf,nsite,nbas,nl,nbasp,npadl,npadr,ncl,npl,i,ib,
     .  jb,ix(3),neul,nbf,nclp,iprint,ifi,fopn
      parameter (nclp=9)
      double precision plat(3,3),platl(3,3),platr(3,3),pgfn(3),xsi(5),xx
      double precision dglob,seref
      character outs*80
      integer, parameter :: lBextz=256

C ... Setup
      nsite = s_ctrl%nsite
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      pass1 = .false.
      nbasp = nbas
      npl   = nbas
      master = 0
      procid = mpipid(1)
      mlog = .false.
      nbasmem = size(s_lat%pos)/3

C ... Pack species labels; make sumeref
      seref = 0
      do  ib = 1, nsite
        i = s_site(ib)%spec
C       Skip over sites with no species (e.g. multipoles)
        if (i > 0) then
          s_site(ib)%clabel = s_spec(i)%name
          seref = seref + s_spec(i)%eref
        endif
      enddo
      s_ham%seref = seref

C ... Start again if ssite was reallocated
   10 pass1 = .not. pass1
      plat = s_lat%plat
      platl = s_lat%platl
      platr = s_lat%platr
      call ptr_lat(s_lat,1,'pos',3,nsite,0,xx)
C     NB: if 2nd pass, pos only given to nbas so far
      call sitepack(s_site,1,nsite,'pos',3,xx,s_lat%pos)
C     call prmx('pos',s_lat%pos,3,3,nsite)
      call ptr_ctrl(s_ctrl,1,'ips',nsite,0,0,xx)
      call sitepack(s_site,1,nbas,'spec',1,s_ctrl%ips,xx)
C     p_pos => s_lat%pos
C     p_ips => s_ctrl%ips
C     Structure s_site should be allocated with
C     sufficient space.  If not, exit with error
      if (nsite > nbasmem) call rxi('increase size of '//
     .  's_site in main program to at least',nsite)
C      do  ib = 1, nsite
C        p_pos(1:3,ib) = s_site(ib)%pos(1:3)
CC       p_ips(ib)   = s_site(ib)%spec
C      enddo

C --- Order-N: group all sites in a crystal into clusters ---
      if (s_ctrl%lordn /= 0) then
        call ptr_ctrl(s_ctrl,1,'cllst',nsite+1,0,0,xx)
        call ptr_ctrl(s_ctrl,1,'clp',(nsite+1)*nclp,0,0,xx)
        call ptr_ctrl(s_ctrl,1,'clssl',nsite,0,0,xx)
        allocate(clord(nsite))
        call clsset(11,s_lat,nsite,s_lat%pos,ncl,s_ctrl%clssl,
     .    s_ctrl%cllst,clord,s_ctrl%clp)

C   ... Reorder the site structure according to clord
C   ... Unpack arrays, now permuted; shorten pos
        ix(1) = 1
        ix(2) = 1
        ix(3) = 1
        call rx('susite: fix old ssite')
        call ptr_ctrl(s_ctrl,2,'clp',nclp*(ncl+1),0,0,xx)
        s_ctrl%ncl = ncl
        deallocate(clord)
      endif

C --- Check and order principal layers, and sites by PL ---
      lpgf = s_ctrl%lpgf(1)
      call ptr_pot(s_pot,8+1,'vshft',8+2*nbasp,0,xx)
      if (lpgf /= 0) then
        call ptr_ctrl(s_ctrl,8+1,'pgfsl',nsite,0,0,xx)
        call ptr_ctrl(s_ctrl,8+1,'pgfvl',nsite,0,0,xx)
        call ptr_ctrl(s_ctrl,8+1,'pgplp',6*(npl+2),0,0,xx)
        call icopy(nsite,int(s_site(1:nsite)%pl),1,s_ctrl%pgfsl,1)
        call icopy(nsite,int(s_site(1:nsite)%plv),1,s_ctrl%pgfvl,1)
        allocate(pgord(2*nbas))
        if (pass1) call pshpr(min(iprint(),1))
C       if (pass1) call pshpr(1)
C       call prmx('pos before pgfset',s_lat%pos,3,3,nbas)
        lorder = .not. pass1
        call pgfset(s_spec,nbas,s_lat%pos,plat,lorder,.true.,
     .    s_ctrl%ips,s_pot%vshft,s_ctrl%pgfsl,s_ctrl%pgfvl,pgord,
     .    pgfn,npl,npadl,npadr,s_ctrl%pgplp)
        nbasp = nbas + (npadl + npadr)
        ix(1) = dglob('nbasp',dble(nbasp),1)
        ix(2) = dglob('npl',dble(npl),1)
        nsite = nbas + 2*(npadl + npadr)
C       No corresponding modification needed for s_site
        if (pass1) then
          deallocate(pgord)
          call poppr
          goto 10
        endif
C   ... Reset pos,pfgsl, since they may be shifted by pgfset
        do  i = 1, nsite
          s_site(i)%pos(1:3) = s_lat%pos(1:3,i)
          s_site(i)%pl = s_ctrl%pgfsl(i)
        enddo
C   ... Reorder the site structure according to pgord
        call siteperm(nbas,pgord,s_site)
        allocate(s_site0(nbas))
        s_site0(1:nbas) = s_site(1:nbas)
        do  ib = 1, nsite
          jb = ib
          if (jb > nbasp+npadl) then
            jb = jb - 2*npadl - 2*npadr
          elseif (jb > nbasp) then
            jb = jb - nbasp
          elseif (jb > nbas+npadl) then
            jb = jb - npadl - npadr
          elseif (jb > nbas) then
            jb = jb - nbas
          endif
          s_site(ib) = s_site0(jb)
C          print *, jb,ib,s_site(ib)%pl
C          print *, 'pgfp',ib,jb
        enddo
        deallocate(s_site0)

C   ... Unpack arrays, now permuted and padded
        do  i = 1, nsite
          s_lat%pos(1:3,i) = s_site(i)%pos(1:3)
          s_ctrl%pgfsl(i) = s_site(i)%pl
          s_ctrl%pgfvl(i) = s_site(i)%plv
          s_ctrl%ips(i) = s_site(i)%spec
        enddo

C   ... Shift doubly padded bas; repack
        if (nbasp > nbas) then
          allocate(wk(3*nsite))
          call pgbasp(nbas,npadl,npadr,s_lat%pos,plat,platl,platr,wk)
          call dscal(9,2d0,platl,1)
          call pgbasp(nbasp,npadl,npadr,wk,plat,platl,platr,s_lat%pos)
          call dscal(9,.5d0,platl,1)
          deallocate(wk)

C     ... Shift pgfsl in padding layers
          call pvsub1(s_ctrl%pgfsl,nbas,npadl,npadr)

C     ... Allocate extra space for pgplp
          call ptr_ctrl(s_ctrl,2,'pgplp',6*(npl+2),0,0,xx)

          do  i = 1, nsite
            s_site(i)%pos(1:3) = s_lat%pos(1:3,i)
            s_site(i)%pl  = s_ctrl%pgfsl(i)
            s_site(i)%plv = s_ctrl%pgfvl(i)
C           s_site(i)%spec = s_ctrl%ips(i)
          enddo
          s_ctrl%npl = npl

C        call prmx('pos, double pad',s_lat%pos,3,3,nsite)
C        call awrit2('%n:1i',' ',100,6,nsite,s_ctrl%pgfsl)
C        call awrit2('%n:1i',' ',100,6,nsite,s_ctrl%pgfvl)
C        call awrit2('%n:1i',' ',100,6,nsite,s_ctrl%ips)

        endif
        deallocate(pgord)
        s_ctrl%nsite = nsite
        s_ctrl%nbasp = nbasp

      else
C   ... No partitioning of hamiltonian: make global pgplp
        npadl = 0
        npadr = 0
        call ptr_ctrl(s_ctrl,8+1,'pgplp',12,0,0,xx)
        call ptr_ctrl(s_ctrl,8+1,'pgfsl',1,0,0,xx)
        s_ctrl%pgplp(1) = -1
        s_ctrl%pgplp(6+1) = nbas
        s_ctrl%npl = 1
      endif

C --- Euler angles and external magnetic fields ---
      if (IAND(s_ctrl%lncol,1+2) /= 0) then
        neul = 1
        call ptr_ham(s_ham,8+1,'eula',nbasp,nl**2*3+3,xx)
        call dpzero(xsi,5)
        if (procid == master) then
          call pvsub2(0,s_site,nbas,nbasp,nl,s_ham%eula,neul,xsi)
          if (cmdopt('--weula',6,0,outs)) then
            ifi = fopn('EULA')
            rewind ifi
            call ioeula(nbasp,-nl,s_ham%eula,neul,0d0,-ifi)
            call fclose(ifi)
            call rx0('wrote Euler angles to file')
          endif
        endif
        call mpibc1(neul,1,2,mlog,'susite','neul')
        call mpibc1(xsi,5,4,mlog,'susite','xsi')
        call ptr_ham(s_ham,2,'eula',nbasp,neul*3+3,xx)
!         call mpibc1(s_ham%eula,nbas*neul*3,4,mlog,'susite','eula')
        call mpi_bcast(s_ham%eula,nbasp*neul*3,mpi_real8,0,mpi_comm_world,err)

C   ... Pack into structures
        s_ctrl%sdxsi = xsi(1:4)
        s_ham%neula = neul
        call ptr_pot(s_pot,8+1,'aamom',nbasp,0,xx)
      else
        neul = 0
        s_ham%neula = 0
        call ptr_ham(s_ham,8+1,'eula',1,1,xx)
      endif
      s_ham%nbf = 0
      if (IAND(s_ctrl%lncol,8+lBextz) /= 0) then
        nbf = 0
        call ptr_ham(s_ham,8+1,'magf',nbasp*nl**2*3,0,xx)
        if (procid == master) then
          call pvsub2(1,s_site,nbas,nbasp,nl,s_ham%magf,nbf,xsi)
        endif
        call mpibc1(nbf,1,2,mlog,'susite','nbf')
        if (nbf == 0) call rx('susite:  external B field sought '//
     .    'but none supplied')
C       call redfrr(omagf, nbasp*nbf*3)
        call ptr_ham(s_ham,2,'magf',nbasp*nbf*3,0,xx)
        call mpibc1(s_ham%magf,nbasp*nbf*3,4,mlog,'susite','magf')
C        call prmx('b field',s_ham%magf,nbasp*nbf,nbasp*nbf,3)
        if (nbf > 0 .and. neul > 0 .and. nbf /= neul) then
          call rx('  input bfield and Euler'//
     .      ' angles must have consistent lm decomposition')
        endif
        s_ham%nbf = nbf
      else
        call ptr_ham(s_ham,8+1,'magf',1,0,xx)
      endif
      s_ctrl%npadl = npadl
      s_ctrl%npadr = npadr

      end
      subroutine pvsub1(pgfsl,nbas,npadl,npadr)
C- Sets PL indices for padding layers
      implicit none
      integer nbas,npadl,npadr,pgfsl(nbas+2*npadl+2*npadr)
      integer nbasp,i

      nbasp = nbas + (npadl + npadr)
      do  i = 1, npadl
        pgfsl(nbas+i) = pgfsl(nbas+i) - 1
        pgfsl(nbasp+i) = pgfsl(nbasp+i) - 2
      enddo
      do  i = 1, npadr
        pgfsl(nbas+npadl+i) = pgfsl(nbas+npadl+i) + 1
        pgfsl(nbasp+npadl+i) = pgfsl(nbasp+npadl+i) + 2
      enddo
      end
      subroutine pvsub2(mode,s_site,nbas,nbasp,nl,eula,neul,xsi)
C- Unpacks Euler angles or magnetic fields; reads from disk
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 for euler angles, 1 for bfield
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read:
Ci     Stored:    eula
Ci     Passed to: spackv
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   nl    :(global maximum l) + 1
Co Outputs
Co   eula  :(mode 0) Euler angles for noncollinear spins
Co         :(mode 1) Magnetic field read
Co   neul  :1 if Euler angles (Bfield) are l-independent, nl otherwise
Co         :unchanged if nothing read
Co   xsi   :global deamon parameters for spin dynamics (not read mode 1)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   29 Jan 03 Added mode to also read b-field
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nbasp,nl,neul
      double precision eula(nbasp,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ifi,fopn,fxst
      double precision xsi(5),xx
      character name*6
      real(8), allocatable :: eulal(:)

      call sanrg(.true.,mode,0,1,'pvsub2','mode')
      if (mode == 0) name = 'eula  '
      if (mode == 1) name = 'bfield'

      allocate(eulal(3*nbasp))
      if (mode == 1) then
        call sitepack(s_site,1,nbasp,'bfield',3,xx,eulal)
      else
        call sitepack(s_site,1,nbasp,'eula',3,xx,eulal)
      endif
      call dmcpy(eulal,1,3,eula,nbasp,1,nbasp,3)
C     call prmx(name//' transposed',eula,nbasp,nbasp,3)

C ... Read angles from disk, if available
      if (fxst(name) == 1) then
        neul = 1
        ifi = fopn(name)
        rewind ifi
        if (mode == 0) call ioeula(nbasp,nl,eula,neul,xsi,ifi)
C       if (mode == 1) call iomagf(0,nbasp,nl,eula,eula,neul,ifi)
        if (mode == 1) call ioextf(0,trim(name),nbasp,nl,eula,eula,3,neul,ifi)
        call fclose(ifi)
      else
        call info0(10,0,0,' susite:   (warning) no file '//trim(name)//
     .    ' ... nothing read')
      endif
C     call prmx(name//' yet again',eula,nbasp*neul,nbasp*neul,3)

C     Repack into structure if neula is 1
      if (neul == 1) then
        call dmcpy(eula,nbasp,1,eulal,1,3,nbasp,3)
        if (mode == 1) then
          call sitepack(s_site,1,nbasp,'-bfield',3,xx,eulal)
        else
          call sitepack(s_site,1,nbasp,'-eula',3,xx,eulal)
        endif
      endif

      deallocate(eulal)

      end

