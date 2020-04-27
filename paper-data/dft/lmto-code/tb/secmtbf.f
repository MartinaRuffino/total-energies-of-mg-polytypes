      subroutine secmtbf(s_ctrl,
     .  plat,nbas,nl,nsp,nspu,ispu,lmx,ipc,indxsh,
     .  ldim,nevmx,efmax,ikp,nkp,bk,nsite,iax,npr,hrs,vso,hso,srs,pot0,
     .  rl,ipr,leig,nev,z,eb,rhrs,rors)
C- Set up tight-binding Hamiltonian; diagonalize secular matrix
C ----------------------------------------------------------------------
Ci Inputs:
Ci   plat,nbas,nl,nsp,lmx,ipc,indxsh
Ci   nsp =  2 for coupled spins (empirical spin-orbit, lso=T)
Ci   nspu = 2 for TB+U, ispu is current spin
Ci   ldim: dimension of l - wave block of hamiltonian matrix
Ci   nevmx, max no. of e-vec's; efmax, highest eigenvalue to seek
Ci   ikp,nkp,bk: current k-point, total number of k-points, k-points
Ci   nsite, total number of neighbors in all clusters;
Ci   iax, neighbor lists; npr, see tbham;
Ci   hrs,srs: real-space hamiltonian and overlap matrices
Ci   vso,hso: table of spin-orbit parameters and the hamiltonian
Ci   ipr:  (verbosity) 0, print nothing; 1, print 1st 9 e'vals and
Ci         timings; 2, print all e'vals and timings
Ci   leig: true: calculate eigenvectors; false: do not
Ci   pot0: monopole potentials, from tbesel
Co Outputs:
Co   nev: number of eigenvectors found from diagno
Co   eigenvalues and eigenvectors are returned in eb, z (leig true)
Co   (with command line --invbl)
Co   rhrs, rsrs: H and O in real space formatted like strux.
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cr   31 Oct 11 (ATP) If switch --wvecs is set then overlap, evecs and
Cr                  evals are written to file 'EV'
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nl,nsp,nspu,ispu,ldim,nevmx,ikp,nkp,nsite,ipr,nev
      integer lmx(0:*),ipc(nbas),indxsh(nl**2*nbas),iax(*),npr(*)
      double precision efmax,plat(3,3),bk(3),hrs(*),vso(*),
     .  hso(*),srs(*),z(*),eb(*),pot0(nbas),rhrs(*),rors(*)
      logical leig,rl
C ... Dynamically allocated arrays
      integer, allocatable :: isup(:),iwk(:),iwk2(:),ifail(:)
      real(8), allocatable :: rwk(:),diawk(:),hk(:),sk(:),zwk(:)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
C ... Local parameters
      integer bitand,nbl,lwork,lrwork,liwork,info,ILAENV
      integer i,j,iprint,lgunit,l2,nevec,ipr0,i1mach,linv,ltb,
     .  lncol,lov,fopn,ifi
      logical bittst,lx,lso,lgamma,cmdopt,lapack,RRR,DC,X,allv
      double precision abstol
      real(8), parameter :: dlamch_s  = 2.2250738585075e-308_8
      character*80 outs
      integer :: real_t

! C heap:
!       integer w(1)
!       common / w /   w

      lapack = .true.
      if (cmdopt('--eispack',9,0,outs)) then
        lapack = .false.
      endif
      RRR    = cmdopt('--RRR',5,0,outs)
      DC     = cmdopt('--DC',4,0,outs)
      X      = cmdopt('--X',3,0,outs)
      allv   = cmdopt('--allvecs',8,0,outs)

      call tcn('secmtb')

C --- Set up dimensions and allocate arrays ---
      ltb = s_ctrl%ltb
      lncol = s_ctrl%lncol
      lov = bitand(ltb,1)
      lso = bittst(lncol,4)
      lgamma = bittst(ltb,2**17) .or. bittst(ltb,2**18)
      real_t = 2
      if (lgamma) real_t = 1

      nevec = 0
      if (leig) then
        if (nevmx == 0 .or. RRR .or. allv) then
          nevec = ldim
        else
          nevec = min(nevmx,ldim)
        endif
      endif
C     osk = 1
      l2 = ldim**2
      allocate(iwk2(nbas))
      allocate(hk(l2*real_t))

      if (lov /= 0) then
        allocate(sk(l2*real_t))
      else
        allocate(sk(1))
      endif

C --- Bloch-transform Hamiltonian and overlap ---
      call tcn('assemble H')
      call tbloch(lgamma,bk,nl,nsp,nspu,ispu,nbas,plat,lmx,ipc,indxsh,
     .  nsite,iax,npr,hrs,vso,hso,lso,ldim,hk,iwk2)
      if (lov /= 0) then
C --- ispu = 1 here 'cos S is spin-independent
        call tbloch(lgamma,bk,nl,nsp,1,1,nbas,plat,lmx,ipc,indxsh,
     .    nsite,iax,npr,srs,vso,hso,.false.,ldim,sk,iwk2)
      endif
      call tcx('assemble H')

C --- add off-site self consistent shifts ---
      if (lov /= 0 .and. rl) then
        call addoff(nl,nbas,lmx,ipc,indxsh,ldim,pot0,sk,hk,iwk2)
      endif

C --- Printout ---
      if (iprint() > 100) then
        if (lgamma) then
          call yprm('ham',1,hk,ldim*ldim,ldim,ldim,ldim)
          if (lov /= 0)
     .    call yprm('ovlp',1,sk,ldim*ldim,ldim,ldim,ldim)
        else
          call yprm('ham',12,hk,ldim*ldim,ldim,ldim,ldim)
          if (lov /= 0)
     .    call yprm('ovlp',12,sk,ldim*ldim,ldim,ldim,ldim)
        endif
      endif

      if (cmdopt('--wvecs',7,0,outs)) then
        ifi = fopn('EV')
        if (ispu == 1) then
          rewind ifi
          if (lov /= 0) then
            call awrit0('Overlap',' ',7,ifi)
            call ywrm(0,'',1,ifi,'(9f15.10)',sk,ldim*ldim,
     .                ldim,ldim,ldim)
          endif
        endif
      endif

C --- Diagonalize Hamiltonian ---
      lx = .true.
      linv = 0
C     if (nevec > 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1
      if (nevec > 0 .and. IAND(s_ctrl%lqp,2) /= 0) linv = 1
      if (linv /= 0) then
        allocate(diawk(ldim*11))
      else
C ...   Use 5*ldim for parallel implementations ...
        allocate(diawk(ldim*5))
      endif
      ipr0 = 0
      if (ipr >= 1 .and. ikp == 1) ipr0 = ipr
      if (lgamma) then
        if (lapack) then
!           abstol = 2*DLAMCH('S')
          abstol = 2*dlamch_s
          liwork = 10*ldim+3
          allocate(iwk(liwork))
          allocate(ifail(ldim))
          if (lov /= 0) then
            if (RRR) then
              call rx('SECMTB: No RRR routine for OVLP')
            elseif (DC) then
              lwork = 1 + 6*ldim + 2*ldim**2
            else
              nbl = ILAENV(1,'DSYTRD','U',ldim,-1,-1,-1)
              lwork = max (max(1,3*ldim-1), (nbl+3)*ldim)
            endif
            allocate(rwk(lwork))
            call tcn('DSYGV')
            if (DC) then
              call DSYGVD(1,'V','U',ldim,hk,ldim,sk,ldim,eb,
     .                    rwk,lwork,iwk,liwork,info)
              call dcopy(ldim**2,hk,1,z,1)
              nev = ldim
            else
              if (X) then
                call DSYGVX(1,'V','I','U',ldim,hk,ldim,sk,ldim,
     .               0d0,0d0,1,nevec,abstol,nev,eb,z,ldim,rwk,
     .               lwork,iwk,ifail,info)
              else
                call DSYGV(1,'V','U',ldim,hk,ldim,sk,ldim,
     .               eb,rwk,lwork,info)
                call dcopy(ldim**2,hk,1,z,1)
                nev = ldim
              endif
            endif
            call tcx('DSYGV')
          else
            if (RRR) then
              nbl = ILAENV(1,'DSYTRD','U',ldim,-1,-1,-1)
              nbl = max( nbl, ILAENV(1,'DORMTR','U',ldim,-1,-1,-1) )
              lwork = max (ldim*(nbl+6) , 26*ldim )
            elseif (DC) then
              lwork = 1 + 6*ldim + 2*ldim**2
            else
              nbl = ILAENV(1,'DSYTRD','U',ldim,-1,-1,-1)
              nbl = max( nbl, ILAENV(1,'DORMTR','U',ldim,-1,-1,-1) )
              lwork = max( ldim*(nbl+3) , 3*ldim )
            endif
            allocate(rwk(lwork))
            call tcn('DSYEV')
            if (RRR) then
              allocate(isup(2*ldim))
              call DSYEVR('V','I','U',ldim,hk,ldim,0d0,0d0,1,nevec,
     .                  abstol,nev,eb,z,ldim,isup,rwk,lwork,
     .                  iwk,liwork,info)
              deallocate(isup)
            elseif (DC) then
              call DSYEVD('V','U',ldim,hk,ldim,eb,
     .                    rwk,lwork,iwk,liwork,info)
              call dcopy(ldim**2,hk,1,z,1)
              nev = ldim
            else
              if (X) then
                call DSYEVX('V','I','U',ldim,hk,ldim,0d0,0d0,1,
     .               nevec,abstol,nev,eb,z,ldim,rwk,lwork,iwk,
     .               ifail,info)
                nev = ldim
              else
                call DSYEV('V','U',ldim,hk,ldim,eb,rwk,lwork,
     .               info)
                call dcopy(ldim**2,hk,1,z,1)
                nev = ldim
              endif
            endif
            deallocate(rwk)
            call tcx('DSYEV')
          endif
          deallocate(iwk,ifail)
        else
          call dsev1(ldim,hk,sk,diawk,ipr0,lx,
     .               lov /= 0,linv,nevec,efmax,nev,z,eb)
        endif
      else
        if (lapack) then
          call tcn('ztoyy')
          call ztoyy(hk,ldim,ldim,ldim,ldim,0,1)
          if (lov /= 0) then
            call ztoyy(sk,ldim,ldim,ldim,ldim,0,1)
          endif
          call tcx('ztoyy')
!           abstol = 2*DLAMCH('S')
          abstol = 2*dlamch_s
          liwork = 10*ldim+3
          allocate(ifail(ldim))
          allocate(iwk(liwork))
          if (lov /= 0) then
            if (RRR) then
              call rx('SECMTB: No RRR routine for OVLP')
            endif
            if (DC) then
              lwork = 2*ldim + ldim**2
              lrwork = 1 + 5*ldim + 2*ldim**2
            else
              nbl = ILAENV(1,'ZHETRD','U',ldim,-1,-1,-1)
              lwork = max (ldim*(nbl+1), 2*ldim)
              lrwork = 7*ldim
            endif
            allocate(rwk(lrwork))
            allocate(zwk(lwork))
            call tcn('ZHEGV')
            if (DC) then
              call ZHEGVD(1,'V','U',ldim,hk,ldim,sk,ldim,eb,
     .             zwk,lwork,rwk,lrwork,iwk,liwork,info)
              call dcopy(2*ldim**2,hk,1,z,1)
              nev = ldim
            else
              if (X) then
                call ZHEGVX(1,'V','I','U',ldim,hk,ldim,sk,
     .               ldim,0d0,0d0,1,nevec,abstol,nev,eb,z,ldim,
     .               zwk,lwork,rwk,iwk,ifail,info)
              else
                call ZHEGV(1,'V','U',ldim,hk,ldim,sk,ldim,eb,
     .               zwk,lwork,rwk,info)
                call dcopy(2*ldim**2,hk,1,z,1)
                nev = ldim
              endif
            endif
            deallocate(rwk,zwk)
            call tcx('ZHEGV')
          else
            if (DC) then
              lwork = 2*ldim + ldim**2
              lrwork = 1 + 5*ldim + 2*ldim**2
            elseif (RRR) then
              nbl = ILAENV(1,'ZHETRD','U',ldim,-1,-1,-1)
              nbl = max( nbl, ILAENV(1,'ZUNMTR','U',ldim,-1,-1,-1) )
              lwork = ldim*(nbl+1)
              lrwork = 24*ldim
            else
              nbl = ILAENV(1,'ZHETRD','U',ldim,-1,-1,-1)
              nbl = max( nbl, ILAENV(1,'ZUNMTR','U',ldim,-1,-1,-1) )
              lwork = ldim*(nbl+1)
              lrwork = 7*ldim
            endif
            allocate(rwk(lrwork))
            allocate(zwk(lwork))
            call tcn('ZHEEV')
            if (DC) then
              call ZHEEVD('V','U',ldim,hk,ldim,eb,zwk,lwork,
     .                    rwk,lrwork,iwk,liwork,info)
              call dcopy(2*ldim**2,hk,1,z,1)
              nev = ldim
            elseif (RRR) then
              allocate(isup(2*ldim))
              call ZHEEVR('V','I','U',ldim,hk,ldim,0d0,0d0,1,
     .             nevec,abstol,nev,eb,z,ldim,isup,zwk,
     .             lwork,rwk,lrwork,iwk,liwork,info)
              deallocate(isup)
            else
              if (X) then
                call ZHEEVX('V','I','U',ldim,hk,ldim,0d0,0d0,1,
     .               nevec,abstol,nev,eb,z,ldim,zwk,lwork,
     .               rwk,iwk,ifail,info)
              else
                call ZHEEV('V','U',ldim,hk,ldim,eb,zwk,lwork,
     .               rwk,info)
                call dcopy(2*ldim**2,hk,1,z,1)
                nev = ldim
              endif
            endif
            deallocate(rwk,zwk)
            call tcx('ZHEEV')
          endif
          deallocate(iwk,ifail)
          call tcn('ztoyy')
          call ztoyy(z,ldim,ldim,ldim,ldim,1,0)
          call tcx('ztoyy')
        else
          call diagno(ldim,hk,sk,diawk,lx,lov,linv,nevec,
     .                efmax,nev,z,eb)
        endif
      endif
      deallocate(diawk)

      if (cmdopt('--wvecs',7,0,outs)) then
        if (nspu == 1) then
          call awrit0('Eigenvectors:',' ',56,ifi)
        else
          call awrit1('Eigenvectors, spin %i:',' ',56,ifi,ispu)
        endif
        call ywrm(0,'',1,ifi,'(9f15.10)',z,ldim*ldim,
     .            ldim,ldim,ldim)
        if (nspu == 1) then
          call awrit1(' %i Eigenvalues:',' ',56,ifi,nev)
        else
          call awrit2(' %i Eigenvalues, spin %i:',' ',56,ifi,nev,ispu)
        endif
        write(ifi,'(9f15.10)') (eb(i), i=1,nev)
      endif
      if (iprint() > 100) then
        if (lgamma) then
          call yprm('z',1,z,ldim*ldim,ldim,ldim,ldim)
        else
          call yprm('z',12,z,ldim*ldim,ldim,ldim,ldim)
        endif
      endif

C --- Printout ---
      if (ipr >= 1) then
        if (lgamma) then
          print *,' SECMTB: hamiltonian real '
        endif
        j = min(9*nsp,ldim)
        if (ipr >= 2) j = nev
        write(lgunit(1),'()')
        call awrit3(' SECMTB:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(1),ikp,nkp,bk)
        write(lgunit(1),'(255(9f8.4:/))') (eb(i), i=1,j)
        if (ipr >= 2) call awrit5
     .    (' nev, nevmx, nevec, ldim=  %i  %i  %i  %i  efmax= %1;5d',
     .    ' ',80,i1mach(2),nev,nevmx,nevec,ldim,efmax)
        call ftflsh(lgunit(1))
      endif
      deallocate(iwk2)
      deallocate(hk,sk)

C --- exit ---
      call tcx('secmtb')
      end
