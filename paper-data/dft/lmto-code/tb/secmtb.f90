      subroutine secmtb(s_ctrl,tbc,plat,nbas,nl,nsp,nspu,ispu,lmx,ipc,indxsh, &
     &  ldim,nevmx,efmax,ikp,nkp,bk,nsite,iax,npr,hrs,vso,hso,srs,pot0,rl,ipr,leig,nev,z,eb)
!- Set up tight-binding Hamiltonian; diagonalize secular matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   plat,nbas,nl,nsp,lmx,ipc,indxsh
!i   nsp =  2 for coupled spins (empirical spin-orbit, lso=T)
!i   nspu = 2 for TB+U, ispu is current spin
!i   ldim: dimension of l - wave block of hamiltonian matrix
!i   nevmx, max no. of e-vec's; efmax, highest eigenvalue to seek
!i   ikp,nkp,bk: current k-point, total number of k-points, k-points
!i   nsite, total number of neighbors in all clusters;
!i   iax, neighbor lists; npr, see tbham;
!i   hrs,srs: real-space hamiltonian and overlap matrices
!i   vso,hso: table of spin-orbit parameters and the hamiltonian
!i   ipr:  (verbosity) 0, print nothing; 1, print 1st 9 e'vals and
!i         timings; 2, print all e'vals and timings
!i   leig: true: calculate eigenvectors; false: do not
!i   pot0: monopole potentials, from tbesel
!o Outputs:
!o   nev: number of eigenvectors found from diagno
!o   eigenvalues and eigenvectors are returned in eb, z (leig true)
!o   (with command line --invbl)
!o   rhrs, rsrs: H and O in real space formatted like strux.
!r Remarks
!u Updates
!u   10 Nov 11 Begin migration to f90 structures
!r   31 Oct 11 (ATP) If switch --wvecs is set then overlap, evecs and
!r                  evals are written to file 'EV'
! ----------------------------------------------------------------------
      use fcuda
      use mpi
      use structures
      use tbprl
      use mod_ctx
!       use mod_crs

      implicit none
! ... Passed parameters
      type(str_ctrl)::  s_ctrl
      type(tbc_t), intent(in), target :: tbc
      integer nbas,nl,nsp,nspu,ispu,ldim,nevmx,ikp,nkp,nsite,ipr,nev
      integer lmx(0:*),ipc(nbas),indxsh(nl**2*nbas),iax(*),npr(*)
      real(8) efmax,plat(3,3),bk(3),hrs(*),vso(*),hso(*),srs(*),eb(*),pot0(nbas)
      real(8), target, intent(inout) :: z(*)
      logical leig,rl

! ... Local parameters
      integer bitand,nbl,lwork,lrwork,liwork,info,ILAENV
!       integer owk,odiawk,ozwk,orwk,oiwk,oifail,oisup
      integer i,j,iprint,lgunit,ipr0,i1mach,linv,ltb,lncol,fopn,ifi
      logical bittst,lx,lso, lov, lgamma,cmdopt,lapack,RRR,DC,X,allv,herm
      double precision abstol,DLAMCH
      character*80 outs
      real(8), allocatable, target, dimension(:) :: hk,sk
!       real(8), pointer, dimension(:) :: lochk,locsk,locz
      integer, allocatable :: wk(:)
      integer, pointer :: desc(:)
      integer :: ctx, typesz, err, nevecs_found, nevals_found, lsz, mtype, ldh, u
      type(cart_ctx_t), pointer :: c2d
      procedure(real(8)) :: ddot
!       real(8), allocatable :: hread(:,:)
!       real(8) :: tmp


      call tcn('secmtb')



! --- Set up dimensions and allocate arrays ---
!     call upack('ctrl ltb lncol',sctrl,ltb,lncol,0,0,0)
      ltb = s_ctrl%ltb
      lncol = s_ctrl%lncol
      lov = bitand(ltb,1) /= 0
      lso = bittst(lncol,4)
      lgamma = bittst(ltb,2**17) .or. bittst(ltb,2**18)

      c2d => tbc % c2d
      typesz = 1
      mtype = mpi_real8
      if (.not. lgamma) typesz = 2
      if (.not. lgamma) mtype = mpi_complex16

      lsz = product(tbc%loclsz)

      ldh = tbc%loclsz(1)


      allocate(hk(typesz*lsz))
      if (lov) then
         allocate(sk(typesz*lsz))
      else
         allocate(sk(typesz))
      end if

!       print *, tbc % c3d % id, ' in secmtb'


! --- Bloch-transform Hamiltonian and overlap ---
      call tcn('assemble H')
      call ztblochd(tbc,lgamma,bk,nl,nsp,nspu,ispu,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,hrs,vso,hso,lso,ldim,hk,ldh,typesz)
      if (lov) then
! --- ispu = 1 here 'cos S is spin-independent
        call ztblochd(tbc,lgamma,bk,nl,nsp,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,srs,vso,hso,.false.,ldim,sk,ldh,typesz)
      endif
      call tcx('assemble H')

!
! !       print *, 'Tr[H], ||H||_F: ', ddot(ldh, 1.0_8, 0, hk, ldh+1), sqrt(ddot(lsz, hk, 1, hk, 1))
!       write(673, "('ikp:',i0,'ispu:',i0)") ikp, ispu
!       call r8_printu(ldh,ldh,hk,ldh,673)

!       write(674, "('ikp:',i0,'ispu:',i0)") ikp, ispu
!       call r8_printu(ldh,ldh,sk,ldh,674)


!       call r8_print(ldh,ldh,hk,'h_dense')
!       call tdense_print_as_crs(hk,ldh,ldh,ldh, 'h_spcrs')
!
!
!       call read_crs_as_dense('h_spcrs', hread)
!
!       tmp = 0.0_8
!       do j = 0, ldh-1
!           do i = 0, ldh-1
!               tmp = tmp + abs(hread(i,j) - hk(i+j*ldh+1))
!           end do
!       end do
!       print *, 'hdiff: ', tmp
!         call rx('stop')
!       u = 2300+tbc%c3d%id
!       do i = 1, ldh
!          do j = 1, ldh
!             write(u,'("(",es10.2,x,es10.2,")")', advance='no') hk((i-1)*ldh+2*j-1),hk((i-1)*ldh+2*j)
!          end do
!          write(u,'(" ")')
!       end do
!       write(u,'(" ")')
!       flush(u)

! --- add off-site self consistent shifts ---

      if (lov .and. rl) call addoff(tbc,nl,nbas,lmx,ipc,indxsh,ldh,pot0,sk,hk,typesz)


!
! ! --- Printout ---
!       if (iprint() > 100) then
!         if (lgamma) then
!           call yprm('ham',1,hk,ldim*ldim,ldim,ldim,ldim)
!           if (lov) call yprm('ovlp',1,sk,ldim*ldim,ldim,ldim,ldim)
!         else
!           call yprm('ham',12,hk,ldim*ldim,ldim,ldim,ldim)
!           if (lov) call yprm('ovlp',12,sk,ldim*ldim,ldim,ldim,ldim)
!         endif
!       endif

!       if (cmdopt('--wvecs',7,0,outs)) then
!         ifi = fopn('EV')
!         if (ispu == 1) then
!           rewind ifi
!           if (lov) then
!             call awrit0('Overlap',' ',7,ifi)
!             call ywrm(0,'',1,ifi,'(9f15.10)',sk,ldim*ldim,ldim,ldim,ldim)
!           endif
!         endif
!       endif

      desc => tbc%desc
!       ctx = desc(2)


      herm = .not. lgamma

!       nevmx = ldh
!       if (.not. elpa) then
      call unievp(lov, herm, tbc%c2d, 1, 'v', 'i', 'u', ldim, hk, 1, 1, desc, sk, 1, 1, desc, &
                     & 0.0_8, 0.0_8, 1, nevmx, nev, nevecs_found, eb, z, 1, 1, desc)
!                      & 0.0_8, 0.0_8, 1, nevmx, nev, nevecs_found, eb, locz, 1, 1, desc)

!         print *, 'evals', nev, eb(1:nev)
!         call rx('evals')

!       write(c2d%id+600,*) eb(1:ldh);  flush(c2d%id+600)
!         call rx('evals')
!       call elpa_dsyev2(desc(3), nevmx, lochk, desc(9), &
!                                  & eb, z, desc(9), desc(5), tbc%c2dr%comm, tbc%c2dc%comm, c2d%comm)
!       if
!       call elpa_dsyev2(desc(3), nevmx, lochk, desc(9), &
!                                  & eb, z, desc(9), desc(5), tbc%c2dr%comm, tbc%c2dc%comm, c2d%comm)

!       nev = nevmx

!       if (tbc%sl) then
!          call tcn('writeback & bcast')
!          call blacs_barrier(ctx, 'a' )
!          call darray_gather(mtype, z, locz, desc, c2d)
!          call mpi_bcast(z, ldh2*typesz, mpi_real8, c2d%rt, c2d%comm, err)
!          deallocate(locz)
!          deallocate(lochk)
!          if (lov) deallocate(locsk)
!          call tcx('writeback & bcast')
!       else
         deallocate(hk)
!          if (lov) deallocate(sk)
         deallocate(sk)
!       end if
      call tcx('secmtb')

      end subroutine secmtb








!
! ! --- Diagonalize Hamiltonian ---
! !#ifdef BLAS3
!       lx = .true.
! !#else
!       lx = .false.
! !#endif
!       linv = 0
! !     if (nevmx > 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1
!       if (nevmx > 0 .and. IAND(s_ctrl%lqp,2_8) /= 0) linv = 1
!       if (linv /= 0) then
!         call defdr(odiawk,ldim*11)
!       else
! ! ...   Use 5*ldim for parallel implementations ...
!         call defdr(odiawk,ldim*5)
!       endif
!       ipr0 = 0
!       if (ipr >= 1 .and. ikp == 1) ipr0 = ipr
!       if (lgamma) then
!         if (lapack) then
!           abstol = 2*DLAMCH('S')
!           liwork = 10*ldim+3
!           call defi(oiwk, liwork)
!           call defi(oifail, ldim)
!           if (lov /= 0) then
!             if (RRR) then
!               call rx('SECMTB: No RRR routine for OVLP')
!             elseif (DC) then
!               lwork = 1 + 6*ldim + 2*ldim**2
!             else
!               nbl = ILAENV(1,'DSYTRD','U',ldim,-1,-1,-1)
!               lwork = max (max(1,3*ldim-1), (nbl+3)*ldim)
!             endif
!             call defrr(orwk, lwork)
!             call tcn('DSYGV')
!             if (DC) then
!               call DSYGVD(1,'V','U',ldim,hk,ldim,sk,ldim,eb,w(orwk),lwork,w(oiwk),liwork,info)
!               call dcopy(ldim**2,hk,1,z,1)
!               nev = ldim
!             else
!               if (X) then
!                 call DSYGVX(1,'V','I','U',ldim,hk,ldim,sk,ldim,0d0,0d0,1, &
!                   & nevmx,abstol,nev,eb,z,ldim,w(orwk),lwork,w(oiwk),w(oifail),info)
!               else
!                 call DSYGV(1,'V','U',ldim,hk,ldim,sk,ldim,eb,w(orwk),lwork,info)
!                 call dcopy(ldim**2,hk,1,z,1)
!                 nev = ldim
!               endif
!             endif
!             call tcx('DSYGV')
!           else
! !#ifdef MPI & NOCUDAD
!            call dsolevp('pdsyevd','v','a','u',ldim,hk,eb,z)
!            nev = ldim
! !#else
!             if (RRR) then
!               nbl = ILAENV(1,'DSYTRD','U',ldim,-1,-1,-1)
!               nbl = max( nbl, ILAENV(1,'DORMTR','U',ldim,-1,-1,-1) )
!               lwork = max (ldim*(nbl+6) , 26*ldim )
!             elseif (DC) then
!               lwork = 1 + 6*ldim + 2*ldim**2
!             else
!               nbl = ILAENV(1,'DSYTRD','U',ldim,-1,-1,-1)
!               nbl = max( nbl, ILAENV(1,'DORMTR','U',ldim,-1,-1,-1) )
!               lwork = max( ldim*(nbl+3) , 3*ldim )
!             endif
!             call defrr(orwk, lwork)
!             call tcn('DSYEV')
!             if (RRR) then
!               call defi(oisup, 2*ldim)
!               call DSYEVR('V','I','U',ldim,hk,ldim,0d0,0d0,1,nevmx,abstol, &
!                   & nev,eb,z,ldim,w(oisup),w(orwk),lwork, w(oiwk),liwork,info)
!             elseif (DC) then
!               call DSYEVD('V','U',ldim,hk,ldim,eb,w(orwk),lwork,w(oiwk),liwork,info)
!               call dcopy(ldim**2,hk,1,z,1)
!               nev = ldim
!             else
!               if (X) then
!                 call DSYEVX('V','I','U',ldim,hk,ldim,0d0,0d0,1,nevmx,abstol, &
!                         & nev,eb,z,ldim,w(orwk),lwork,w(oiwk),w(oifail),info)
!                 nev = ldim
!               else
! ! enter perturbation branch for a test
! !                 call eigapprox(ldim, hk, eb)
! !                 call carpar(ldim, hk, eb)
!                 call DSYEV('V','U',ldim,hk,ldim,eb,w(orwk),lwork,info)
!                 call dcopy(ldim**2,hk,1,z,1)
!                 nev = ldim
!               endif
!             endif
!             call tcx('DSYEV')
! !#endif
!           endif
!           call rlse(oiwk)
!         else
!           call dsev1(ldim,hk,sk,w(odiawk),ipr0,lx,lov /= 0,linv,nevmx,efmax,nev,z,eb)
!         endif
!       else
!         if (lapack) then
!           call tcn('ztoyy')
!           call ztoyy(hk,ldim,ldim,ldim,ldim,0,1)
!           if (lov /= 0) then
!             call ztoyy(sk,ldim,ldim,ldim,ldim,0,1)
!           endif
!           call tcx('ztoyy')
!           abstol = 2*DLAMCH('S')
!           liwork = 10*ldim+3
!           call defi(oiwk, liwork)
!           call defi(oifail, ldim)
!           if (lov /= 0) then
!             if (RRR) then
!               call rx('SECMTB: No RRR routine for OVLP')
!             endif
!             if (DC) then
!               lwork = 2*ldim + ldim**2
!               lrwork = 1 + 5*ldim + 2*ldim**2
!             else
!               nbl = ILAENV(1,'ZHETRD','U',ldim,-1,-1,-1)
!               lwork = max (ldim*(nbl+1), 2*ldim)
!               lrwork = 7*ldim
!             endif
!             call defcc(ozwk, lwork)
!             call defrr(orwk, lrwork)
!             call tcn('ZHEGV')
!             if (DC) then
!               call ZHEGVD(1,'V','U',ldim,hk,ldim,sk,ldim,eb, &
!                      & w(ozwk),lwork,w(orwk),lrwork,w(oiwk),liwork,info)
!               call dcopy(2*ldim**2,hk,1,z,1)
!               nev = ldim
!             else
!               if (X) then
!                 call ZHEGVX(1,'V','I','U',ldim,hk,ldim,sk,ldim,0d0,0d0,1, &
!                   & nevmx,abstol,nev,eb,z,ldim,w(ozwk),lwork,w(orwk),w(oiwk),w(oifail),info)
!               else
!                 call ZHEGV(1,'V','U',ldim,hk,ldim,sk,ldim,eb,w(ozwk),lwork,w(orwk),info)
!                 call dcopy(2*ldim**2,hk,1,z,1)
!                 nev = ldim
!               endif
!             endif
!             call tcx('ZHEGV')
!           else
! ! ! C#ifdefC MPI
! ! ! C            call dsolevp('pzheevd','v','a','u',ldim,hk,eb,z)
! ! ! C#else
!             if (DC) then
!               lwork = 2*ldim + ldim**2
!               lrwork = 1 + 5*ldim + 2*ldim**2
!             elseif (RRR) then
!               nbl = ILAENV(1,'ZHETRD','U',ldim,-1,-1,-1)
!               nbl = max( nbl, ILAENV(1,'ZUNMTR','U',ldim,-1,-1,-1) )
!               lwork = ldim*(nbl+1)
!               lrwork = 24*ldim
!             else
!               nbl = ILAENV(1,'ZHETRD','U',ldim,-1,-1,-1)
!               nbl = max( nbl, ILAENV(1,'ZUNMTR','U',ldim,-1,-1,-1) )
!               lwork = ldim*(nbl+1)
!               lrwork = 7*ldim
!             endif
!             call defcc(ozwk, lwork)
!             call defrr(orwk, lrwork)
!             call tcn('ZHEEV')
!             if (DC) then
!               call ZHEEVD('V','U',ldim,hk,ldim,eb,w(ozwk),lwork,w(orwk),lrwork,w(oiwk),liwork,info)
!               call dcopy(2*ldim**2,hk,1,z,1)
!               nev = ldim
!             elseif (RRR) then
!               call defi(oisup, 2*ldim)
!               call ZHEEVR('V','I','U',ldim,hk,ldim,0d0,0d0,1,nevmx,abstol,&
!                   & nev,eb,z,ldim,w(oisup),w(ozwk),lwork,w(orwk),lrwork,w(oiwk),liwork,info)
!             else
!               if (X) then
!                 call ZHEEVX('V','I','U',ldim,hk,ldim,0d0,0d0,1, &
!                    & nevmx,abstol,nev,eb,z,ldim,w(ozwk),lwork,w(orwk),w(oiwk),w(oifail),info)
!               else
!                 call ZHEEV('V','U',ldim,hk,ldim,eb,w(ozwk),lwork,w(orwk),info)
!                 call dcopy(2*ldim**2,hk,1,z,1)
!                 nev = ldim
!               endif
!             endif
!             call tcx('ZHEEV')
! ! ! C#endif
!           endif
! !           call tcn('ztoyy')
! !           call ztoyy(z,ldim,ldim,ldim,ldim,1,0)
! !           call tcx('ztoyy')
!           call rlse(oiwk)
!         else
!           call diagno(ldim,hk,sk,w(odiawk),lx,lov,linv,nevmx,efmax,nev,z,eb)
!         endif
!       endif
!
!       if (cmdopt('--wvecs',7,0,outs)) then
!         if (nspu == 1) then
!           call awrit0('Eigenvectors:',' ',56,ifi)
!         else
!           call awrit1('Eigenvectors, spin %i:',' ',56,ifi,ispu)
!         endif
!         call ywrm(0,'',1,ifi,'(9f15.10)',z,ldim*ldim,ldim,ldim,ldim)
!         if (nspu == 1) then
!           call awrit1(' %i Eigenvalues:',' ',56,ifi,nev)
!         else
!           call awrit2(' %i Eigenvalues, spin %i:',' ',56,ifi,nev,ispu)
!         endif
!         write(ifi,'(9f15.10)') (eb(i), i=1,nev)
!       endif
!       if (iprint() > 100) then
!         if (lgamma) then
!           call yprm('z',1,z,ldim*ldim,ldim,ldim,ldim)
!         else
!           call yprm('z',12,z,ldim*ldim,ldim,ldim,ldim)
!         endif
!       endif
!
! ! --- Printout ---
!       if (ipr >= 1) then
!         if (lgamma) then
!           print *,' SECMTB: hamiltonian real '
!         endif
!         j = min(9*nsp,ldim)
!         if (ipr >= 2) j = nev
!         write(lgunit(1),'()')
!         call awrit3(' SECMTB:  kpt %i of %i, k=%3:2,5;5d',' ',80,lgunit(1),ikp,nkp,bk)
!         write(lgunit(1),'(255(9f8.4:/))') (eb(i), i=1,j)
!         if (ipr >= 2) call awrit5(' nev, nevmx, nevmx, ldim=  %i  %i  %i  %i  efmax= %1;5d', &
!                                                    &  ' ',80,i1mach(2),nev,nevmx,nevec,ldim,efmax)
!         call ftflsh(lgunit(1))
!       endif
!       call rlse(owk)
!
! ! --- exit ---
