   subroutine unievp(ovlp, herm, c2d, itype, jobz, rng, uplo, n, a, ia, ja, desca, b, ib, jb, descb, &
                                                      & vl, vu, il, iu, m, nz, w, z, iz, jz, descz)

!- Unified (SCA)LAPACK eigenvalue problem solution
! Significant portion of the following is sourced from the SCALAPACK manpages and source files.
! For further details please refer to the specific routines documentation there.
!------------------------------------------------------------------------------------
!i Arguments:
!i ovlp : Overlap available, uses a generalised routine
!i
!i herm : Notes the the matrices are of complex hermitian type
!i
!i itype: Specifies the problem type to be solved:
!i       =1  A*x   = (lambda)*B*x (mostly used here)
!i       =2  A*B*x = (lambda)*x
!i       =3  B*A*x = (lambda)*x
!i
!i jobz  :='N'  Compute eigenvalues only;
!i       :='V'  Compute eigenvalues and eigenvectors.
!i
!i rng   := 'A' all eigenvalues will be found.
!i        = 'V' all eigenvalues in the half-open interval (VL,VU] will be found.
!i        = 'I' the IL-th through IU-th eigenvalues will be found.
!i
!i uplo  := 'U'  Upper triangle of A and B are stored;
!i        = 'L'  Lower triangle of A and B are stored.
!i
!i n     : The order of the matrices A and B.  N >= 0.
!i
!i a     : Contains the local pieces of the N-by-N Hermitian distributed A to be diagonalised.
!i
!i ia    : The row index in the global array A indicating the first row of sub( A ).
!i         Values other than 1 have not been tested.
!i
!i ja    : The row index in the global array A indicating the first column of sub( A ).
!i         Values other than 1 have not been tested.
!i
!i desca : BLACS/SCALAPACK descriptor for 'a'. Here is the detailed explanation from the source itself
!i
!i          Each global data object is described by an associated description
!i          vector.  This vector stores the information required to establish
!i          the mapping between an object element and its corresponding process
!i          and memory location.
!i
!i          Let A be a generic term for any 2D block cyclicly distributed array.
!i          Such a global array has an associated description vector DESCA.
!i          In the following comments, the character _ should be read as
!i          "of the global array".
!i
!i          NOTATION        STORED IN      EXPLANATION
!i          --------------- -------------- --------------------------------------
!i          DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!i                                         DTYPE_A = 1.
!i          CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!i                                         the BLACS process grid A is distribu-
!i                                         ted over. The context itself is glo-
!i                                         bal, but the handle (the integer
!i                                         value) may vary.
!i          M_A    (global) DESCA( M_ )    The number of rows in the global
!i                                         array A.
!i          N_A    (global) DESCA( N_ )    The number of columns in the global
!i                                         array A.
!i          MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!i                                         the rows of the array.
!i          NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!i                                         the columns of the array.
!i          RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!i                                         row of the array A is distributed.
!i          CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!i                                         first column of the array A is
!i                                         distributed.
!i          LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!i                                        array.  LLD_A >= MAX(1,LOCr(M_A)).
!i
!i          Let K be the number of rows or columns of a distributed matrix,
!i          and assume that its process grid has dimension p x q.
!i          LOCr( K ) denotes the number of elements of K that a process
!i          would receive if K were distributed over the p processes of its
!i          process column.
!i          Similarly, LOCc( K ) denotes the number of elements of K that a
!i          process would receive if K were distributed over the q processes of
!i          its process row.
!i          The values of LOCr() and LOCc() may be determined via a call to the
!i          ScaLAPACK tool function, NUMROC:
!i                  LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!i                  LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!i          An upper bound for these quantities may be computed by:
!i                  LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!i                  LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!i
!i b     : The local pieces of the overlap matrix B.
!i
!i ib    : See ia in context of B
!i
!i jb    : See ja in context of B
!i
!i descb : Descriptor for B. See desca for details
!i
!i vl    : If RANGE='V', the lower bound of the interval to be searched for eigenvalues.
!i           Not referenced if RANGE = 'A' or 'I'.
!i
!i vu    : If RANGE='V', the upper bound of the interval to be searched for eigenvalues.
!i           Not referenced if RANGE = 'A' or 'I'.
!i
!i il    : If RANGE='I', the index (from smallest to largest) of the
!i           smallest eigenvalue to be returned.  IL >= 1.
!i           Not referenced if RANGE = 'A' or 'V'.
!i
!i iu    : If RANGE='I', the index (from smallest to largest) of the
!i         largest eigenvalue to be returned.  min(IL,N) <= IU <= N.
!i         Not referenced if RANGE = 'A' or 'V'.
!i
!i iz, jz, descz    : see ia, ja, desca
!i
!o Outputs:
!o m     : Total number of eigenvalues found.  0 <= M <= N.
!o
!o nz    : Total number of eigenvectors computed.  0 <= NZ <= M.
!o         The number of columns of Z that are filled.
!o         If JOBZ /= 'V', NZ is not referenced.
!o         If JOBZ == 'V', NZ = M unless the user supplies
!o         insufficient space and PZHEGVX is not able to detect this
!o         before beginning computation.  To get all the eigenvectors
!o         requested, the user must supply both sufficient
!o         space to hold the eigenvectors in Z (M .LE. DESCZ(N_))
!o         and sufficient workspace to compute them.  (See LWORK below.)
!o         PZHEGVX is always able to detect insufficient space without
!o         computation unless RANGE .EQ. 'V'.
!o
!o w     : On normal exit, the first M entries contain the selected eigenvalues in ascending order.
!o
!o
!o z     : If JOBZ = 'V', then on normal exit the first M columns of Z
!o         contain the orthonormal eigenvectors of the matrix
!o         corresponding to the selected eigenvalues.  If an eigenvector
!o         fails to converge, then that column of Z contains the latest
!o         approximation to the eigenvector, and the index of the
!o         eigenvector is returned in IFAIL. If JOBZ = 'N', then Z is not referenced.
!o
!o
!s Command-line switches:
!s
!s '--X'   : use the 'advanced' routines, '-x' suffixes (default)
!s
!s '--RRR' : use the Multiple Relatively Robust Representations (MRRR a.k.a. MR^3)
!s            routines, limited availability to some orthogonal cases only,
!s           If the full spectrum of eigenvalues is needed these
!s           will likely be the fastest routines with O(N^2) scalling.
!s
!s '--DC'  : Full spectrum 'divide and conquer'
!s
!r Remarks:
!r The following routines are covered
!r
!r  dsyevx:  dsyevx( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info )
!r  dsyevr:  dsyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info )
!r  dsyevd:  dsyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info )
!r  dsygvx:  dsygvx( itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info )
!r  dsygvr:
!r  dsygvd:  dsygvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info )
!r  zheevx:  zheevx( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info )
!r  zheevr:  zheevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info )
!r  zheevd:  zheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info )
!r  zhegvx:  zhegvx( itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info )
!r  zhegvr:
!r  zhegvd:  zhegvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info )
!r pdsyevx: pdsyevx( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info )
!r pdsyevr: pdsyevr( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work, lwork, iwork, liwork, info )
!r pdsyevd: pdsyevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, iwork, liwork, info )
!r pdsygvx: pdsygvx( ibtype, jobz, range, uplo, n, a, ia, ja, desca, b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info )
!r pdsygvr:
!r pdsygvd:
!r pzheevx: pzheevx( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info )
!r pzheevr: pzheevr( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, info )
!r pzheevd: pzheevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, info )
!r pzhegvx: pzhegvx( ibtype, jobz, range, uplo, n, a, ia, ja, desca, b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info )
!r pzhegvr:
!r pzhegvd:
!r
!r So far only the error messages from pdsyevx are processed. Nonzero info from all others will lead to exit.
!r If SCALAPACK is enabled routines expect block cyclicly distributed local pieces of the virtual global array.
!r The input arrays A and B are NOT destroyed or overwritten on exit!
!r Worspace inquiry is made on each call and the approriate amounts and types are allocated.
!r The 'advanced' routines are used by default if no switch is given.
!r In the special case of (ovlp=f and rng='a' and jobz='n') the RRR routines are used
!r In case the comminicator/context size is 1 the LAPACK routines are used instead of SCALAPACS's
!
!  2013 DP

      use mpi
      use mod_ctx
      use fmagma

      implicit none
      integer, parameter :: dp = 8

      logical, intent(in) :: ovlp, herm
      character, intent(in) :: jobz, rng, uplo
      integer, intent(in)  :: ia, ib, itype, il, iu, iz, ja, jb, jz, n, desca(9), descb(9), descz(9)
      integer, intent(inout) :: m, nz
      real(dp), intent(in) :: vl, vu
      real(dp), intent(inout) :: w(*)
      real(dp), intent(in), target :: a(*), b(*)
      real(dp), intent(inout), target :: z(*)

      type(cart_ctx_t), intent(in) :: c2d


      real(dp), allocatable :: gap(:), work(:), rwork(:)
      integer, allocatable :: iclustr(:), ifail(:), isuppz(:), iwork(:)
      integer :: typesz
      character(len=7) :: task
      character(len=16) :: outs
      real(dp), pointer :: loca(:), locb(:)
      logical :: alg_adv, alg_dac, alg_rrr, sl, alg_el2, alg_m2s, mg
      integer :: alg_idx, ctx, psizes(2), pcrds(2), lsizes(2), bsizes(2), rsizes(2), bl_prc(2), gsizes(2)  &
               , info, lda, ldb, ldz, nb, mb, lwork, liwork, lrwork, iproc, nproc, llinsz, blocks(2), i,j,u, err, ng &
               , descal(9), descbl(9), desczl(9)
      real(dp) :: abstol, orfac
!       procedure(real(dp)) :: dlamch, pdlamch
      real(dp), parameter :: dlamch_s = 2.2250738585075e-308_dp
      procedure(logical) :: cmdopt
      type(cart_ctx_t) :: rcomm, ccomm

      info = 0
      ctx    = desca(2)
      rsizes = desca(3:4)
!       call blacs_gridinfo(ctx, psizes(1), psizes(2), pcrds(1), pcrds(2))
!       iproc = pcrds(1)*psizes(2) + pcrds(2)
!       nproc = product(psizes)
      psizes = c2d % szs(0:1)
      pcrds  = c2d % crd(0:1)
      nproc = c2d % sz
      sl = nproc > 1

      mg = cmdopt('--MAGMA',7,0,outs)

      if (sl .and. mg) call rx('MAGMA and ScaLAPACK are not simultaneously supported.')
!       sl = .true.

!       print *, 'unievp: ', iu
      alg_idx = 1
! Looks like the rrr routines suffer more fp issues so the default will now be the dcx for a while.
!       if ( (.not. ovlp) .and. (rng=='a' .or. rng=='A' .or. (il == 1 .and. iu == n))) alg_idx = 2

      if      (cmdopt('--X' ,3,0,outs))  then
         alg_idx = 1
      else if (cmdopt('--RRR',5,0,outs)) then
         alg_idx = 2
      else if (cmdopt('--DC',4,0,outs)) then
         alg_idx = 3
      else if (cmdopt('--ELPA2',7,0,outs) .and. sl) then
         alg_idx = 4
      else if (cmdopt('--MAGMA2',8,0,outs)) then
         alg_idx = 5
      end if

      alg_adv = alg_idx == 1
      alg_rrr = alg_idx == 2
      alg_dac = alg_idx == 3
      alg_el2 = alg_idx == 4
      alg_m2s = alg_idx == 5

      if (sl) then
         bsizes = desca(5:6)
         blocks = (rsizes - 1) / bsizes + 1
         bl_prc = (blocks - 1) / psizes + 1
         lsizes = bl_prc * bsizes
         gsizes = lsizes * psizes
      else
         lsizes = rsizes
         gsizes = rsizes
      end if

      llinsz = product(lsizes)

      lda = desca(9)
      ldb = descb(9)
      ldz = descz(9)

      if (lsizes(1) > lda .or. (ovlp .and. lsizes(1) > ldb) .or. lsizes(1) > ldz) &
                  & call rx('UNIEVP: Calculated and passes leading matrix dimensions do not match')

!                  1234567
      task(1:7) = ' dsyevx'

      if (sl) task(1:1) = 'p'
      if (mg) task(1:1) = 'm'
      if (herm) task(2:4) = 'zhe'
      if (ovlp) task(5:5) = 'g'
      if (alg_m2s) task(6:6) = '2'
      if (alg_rrr) task(7:7) = 'r'
      if (alg_dac) task(7:7) = 'd'
      if (alg_el2) task(5:7) = 'el2'


      call tcn(task)

      if (alg_adv) then
         !orfac = 0.001_dp
         orfac = -1
         allocate(ifail(n))
         ifail = 0
         if (sl) then
            allocate(iclustr(2*nproc),gap(nproc))
            iclustr = 0; gap = 0
         end if
      end if

      if (alg_rrr) then
         allocate(isuppz(2*n))
         isuppz = 0
      end if

      if (alg_adv .or. alg_rrr) then
!          if (sl) then
!             abstol = 2*pdlamch(ctx, 's')
!          else
!             abstol = 2*dlamch('s')
!          endif
!          abstol = 2*dlamch_s ! this gains little but forces zheevx to near serial mode.
         abstol = 0
      end if

      if (alg_el2) then
         m = iu-il+1
         rcomm = sub_cart_ctx(c2d, [.true. , .false.])
         ccomm = sub_cart_ctx(c2d, [.false., .true. ])
      end if


      typesz = 1
      if (herm) typesz = 2

      locb => null()
      if (ovlp) then
         allocate(locb(typesz*llinsz))
!          call dcopy(typesz*llinsz, b, 1, locb, 1)
         call dcopy2d(typesz*lsizes(1), lsizes(2), b, typesz*ldb, locb, typesz*lsizes(1))
         ldb = lsizes(1)
         descbl = descb
         descbl(9) = ldb
      end if

      if ((alg_dac .and. (.not. sl)) .or. mg) then
         loca => z(1:typesz*ldz*lsizes(2))
         lda = ldz
      else
         allocate(loca(typesz*llinsz))
         lda = lsizes(1)
      end if

!       call dcopy(typesz*llinsz, a, 1, loca, 1)
      call dcopy2d(typesz*lsizes(1), lsizes(2), a, typesz*desca(9), loca, typesz*lda)

      descal = desca
      descal(9) = lda

      if (.not. alg_el2) then
         lwork = -1
         allocate(work(4*typesz))
         work = 1.0_dp;


         if (alg_adv .and. (.not. sl) .and. (.not. mg)) then
            liwork = 5*n
            allocate(iwork(liwork))
            if (herm) then
               lrwork = 7*n;
               allocate(rwork(lrwork))
            end if
         else
            liwork = -1
            allocate(iwork(4))
            iwork = 1.0_dp
            if (herm) then
               lrwork = -1;
               allocate(rwork(4))
               rwork = 1.0_dp
            end if
         end if
!
!       u = 3000+iproc
!       write(u, '(2(x,i0))') lsizes
!       do i = 1, lsizes(1)
!          do j = 1, lsizes(2)
!             write(u, '(x,es10.2)' ,advance='no') loca((j-1)*lsizes(1)+i)
!          end do
!          write(u,'("")')
!       end do
!       flush(u)
!       call mpi_barrier(mpi_comm_world, err)
!
!       stop


      if      (task == ' dsyevx') then
         call  dsyevx(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      else if (task == ' dsyevr') then
         call  dsyevr(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, &
            & info)
      else if (task == ' dsyevd') then
         call  dsyevd(jobz, uplo, n, loca, lda, w, work, lwork, iwork, liwork, info)
      else if (task == ' dsygvx') then
         call  dsygvx( itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, &
            & iwork, ifail, info)
      else if (task == ' dsygvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == ' dsygvd') then
         call  dsygvd(itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, iwork, liwork, info)
      else if (task == ' zheevx') then
         call  zheevx(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, &
            & info)
      else if (task == ' zheevr') then
         call  zheevr(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, &
            & iwork, liwork, info)
      else if (task == ' zheevd') then
         call  zheevd(jobz, uplo, n, loca, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == ' zhegvx') then
         call  zhegvx(itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork,&
            & iwork, ifail, info)
      else if (task == ' zhegvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == ' zhegvd') then
         call  zhegvd(itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'pdsyevx') then
         call pdsyevx(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work,&
            & lwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pdsyevr') then
         call pdsyevr(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work, lwork, iwork, &
            & liwork, info)
      else if (task == 'pdsyevd') then
         call pdsyevd(jobz, uplo, n, loca, ia, ja, descal, w, z, iz, jz, descz, work, lwork, iwork, liwork, info)
      else if (task == 'pdsygvx') then
         call pdsygvx(itype, jobz, rng, uplo, n, loca, ia, ja, descal, locb, ib, jb, descbl, vl, vu, il, iu, abstol, m, nz, w, &
            & orfac, z, iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pdsygvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == 'pdsygvd') then
         call rx('PDSYGVD not available')
      else if (task == 'pzheevx') then
         call pzheevx(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work,&
            & lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pzheevr') then
         call pzheevr(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work, lwork, rwork, &
            & lrwork, iwork, liwork, info)
      else if (task == 'pzheevd') then
         call pzheevd(jobz, uplo, n, loca, ia, ja, descal, w, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'pzhegvx') then
         call pzhegvx(itype, jobz, rng, uplo, n, loca, ia, ja, descal, locb, ib, jb, descbl, vl, vu, il, iu, abstol, m, nz, &
            &  w, orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pzhegvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == 'pzhegvd') then
         call rx('PZHEGVD not available (2013)')
!-- The new magma-1.4.0 is decent
      else if (task == 'mdsyevd') then
         call mdsyevd(ng, jobz, uplo, n, loca, lda, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsyevx') then
         call mdsyevx(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsye2x') then
         call mdsye2x(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsygvd') then
         call mdsygvd(ng, itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsygvx') then
         call mdsygvx(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il,iu,m,w,work,lwork,iwork,liwork,info)
      else if (task == 'mdsyg2x') then
         call mdsyg2x(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il,iu,m,w,work,lwork,iwork,liwork,info)
      else if (task == 'mzheevd') then
         call mzheevd(ng, jobz, uplo, n, loca, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'mzheevx') then
         call mzheevx(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, rwork,lrwork,iwork,liwork,info)
      else if (task == 'mzhee2x') then
         call mzhee2x(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, rwork,lrwork,iwork,liwork,info)
      else if (task == 'mzhegvd') then
         call mzhegvd(ng, itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'mzhegvx') then
         call mzhegvx(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, m, w, work, lwork, rwork, &
                                                                                              lrwork, iwork, liwork, info)
      else if (task == 'mzheg2x') then
         call mzheg2x(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork,&
                                                                                                         iwork, liwork, info)
      else
         call rx(task//' unknown. Check the supported routines in unievp.f90')
      end if

         lwork  = max(4,int(work(1)))
! allow reorthogonalisation of up to 9 vectors (lwork += 8*n) to avoid having to repeat pdsy?vx due to clustered evals
         if (sl .and. alg_adv .and. (.not. herm)) lwork = lwork + nint(n/sqrt(real(nproc,dp))-1)*n

         deallocate(work)
         allocate(work(lwork*typesz))

         if ((.not. alg_adv) .or. sl .or. mg) then
            liwork = max(4,int(iwork(1)))
            deallocate(iwork)
            allocate(iwork(liwork))
            if (herm) then
               lrwork = max(4,int(rwork(1)))
! to avoid having to repeat pzhe?vx due to clustered evals allow reorthogonalisation of up to 9 vectors (lrwork += 8*n)
               if (alg_adv) lrwork = lrwork + 8*n
               deallocate(rwork)
               allocate(rwork(lrwork))
            end if
         end if
!       print *, 'lwork, liwork, lrwork', lwork, liwork, lrwork
      end if



      if      (task == ' dsyevx') then
         call  dsyevx(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      else if (task == ' dsyevr') then
         call  dsyevr(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, &
             info)
      else if (task == ' dsyevd') then
         call  dsyevd(jobz, uplo, n, loca, lda, w, work, lwork, iwork, liwork, info)
      else if (task == ' dsygvx') then
         call  dsygvx( itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, &
             iwork, ifail, info)
      else if (task == ' dsygvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == ' dsygvd') then
         call  dsygvd(itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, iwork, liwork, info)
      else if (task == ' zheevx') then
         call  zheevx(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, &
             info)
      else if (task == ' zheevr') then
         call  zheevr(jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, &
             iwork, liwork, info)
      else if (task == ' zheevd') then
         call  zheevd(jobz, uplo, n, loca, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == ' zhegvx') then
         call  zhegvx(itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork,&
             iwork, ifail, info)
      else if (task == ' zhegvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == ' zhegvd') then
         call  zhegvd(itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'pdsyevx') then
         call pdsyevx(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work,&
             lwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pdsyevr') then
         call pdsyevr(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work, lwork, iwork, &
             liwork, info)
      else if (task == 'pdsyevd') then
         call pdsyevd(jobz, uplo, n, loca, ia, ja, descal, w, z, iz, jz, descz, work, lwork, iwork, liwork, info)
      else if (task == 'pdsygvx') then
         call pdsygvx(itype, jobz, rng, uplo, n, loca, ia, ja, descal, locb, ib, jb, descbl, vl, vu, il, iu, abstol, m, nz, w, &
             orfac, z, iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pdsygvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == 'pdsygvd') then
         call rx('PDSYGVD not available')
      else if (task == 'pzheevx') then
         call pzheevx(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz, work,&
             lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pzheevr') then
         call pzheevr(jobz, rng, uplo, n, loca, ia, ja, descal, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work, lwork, rwork, &
             lrwork, iwork, liwork, info)
      else if (task == 'pzheevd') then
         call pzheevd(jobz, uplo, n, loca, ia, ja, descal, w, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'pzhegvx') then
         call pzhegvx(itype, jobz, rng, uplo, n, loca, ia, ja, descal, locb, ib, jb, descbl, vl, vu, il, iu, abstol, m, nz, &
              w, orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
      else if (task == 'pzhegvr') then
         call rx('No RRR routine for the generalised eigenvalue problem (2013)')
      else if (task == 'pzhegvd') then
         call rx('PZHEGVD not available (2013)')
      else if (task == 'pdsyel2') then
         call pdsyel2(n, m, loca, lda, w, z, ldz, descal(5), rcomm%comm, ccomm%comm, c2d%comm)
      else if (task == 'pzheel2') then
         call pzheel2(n, m, loca, lda, w, z, ldz, descal(5), rcomm%comm, ccomm%comm, c2d%comm)
!-- The new magma-1.4.0 is decent
      else if (task == 'mdsyevd') then
         call mdsyevd(ng, jobz, uplo, n, loca, lda, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsyevx') then
         call mdsyevx(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsye2x') then
         call mdsye2x(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsygvd') then
         call mdsygvd(ng, itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, iwork, liwork, info)
      else if (task == 'mdsygvx') then
         call mdsygvx(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il,iu,m,w,work,lwork,iwork,liwork,info)
      else if (task == 'mdsyg2x') then
         call mdsyg2x(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il,iu,m,w,work,lwork,iwork,liwork,info)
      else if (task == 'mzheevd') then
         call mzheevd(ng, jobz, uplo, n, loca, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'mzheevx') then
         call mzheevx(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, rwork,lrwork,iwork,liwork,info)
      else if (task == 'mzhee2x') then
         call mzhee2x(ng, jobz, rng, uplo, n, loca, lda, vl, vu, il, iu, m, w, work, lwork, rwork,lrwork,iwork,liwork,info)
      else if (task == 'mzhegvd') then
         call mzhegvd(ng, itype, jobz, uplo, n, loca, lda, locb, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
      else if (task == 'mzhegvx') then
         call mzhegvx(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, m, w, work, lwork, rwork, &
                                                                                              lrwork, iwork, liwork, info)
      else if (task == 'mzheg2x') then
         call mzheg2x(ng, itype, jobz, rng, uplo, n, loca, lda, locb, ldb, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork,&
                                                                                                         iwork, liwork, info)
      else
         call rx(task//' unknown')
      end if

      if (alg_dac) then
         m = n
         if (jobz == 'v') nz = n
      end if

      if (.not. sl .and. jobz == 'v') nz = n

      if (info /= 0 .and. all(pcrds == 0)) then
      if (alg_adv .and. sl) then
         if (info < 0) then
            print '(a)', ' unievp: '//trim(adjustl(task))//' returned error: ', info
         else
            if (mod(info,2) /= 0) then
               print '(a)', ' WARNING: the following eigenvectors failed to converge: '
               i=1; do while (ifail(i) /= 0); write (*,'(x,i0)',advance='no') ifail(i); i=i+1; end do
            else if (mod(info/2,2) /= 0) then
               print '(a)', ' WARNING: Eigenvectors corresponding to the following eigenvalue cluster' &
                   //' could not be reorthogonalised due to insufficient workspace: '
               i=1; do while (iclustr(i) /= 0); write (*,'(x,i0)',advance='no') iclustr(i); i=i+1; end do
               print *, ''
            else if (mod(info/4,2) /= 0) then
               print '(a,i0,a,i0,a,i0)', ' WARNING: Space limit prevented '//trim(adjustl(task)) &
                   //' from computing all of the eigenvectors between ',VL,' and ',VU,           &
                   '.  The number of eigenvectors computed is returned is', nz
            else if (mod(info/8,2) /= 0) then
               print '(a)', ' WARNING: pdstebz failed to compute eigenvalues.'
            end if
         end if
      end if
      end if

      if (.not. alg_el2) then
         deallocate(work, iwork)
         if (herm) deallocate(rwork)
      else
         call del_cart_ctx(rcomm)
         call del_cart_ctx(ccomm)
      end if
      if (ovlp) deallocate(locb)
      if (((.not. alg_dac) .or. sl) .and. (.not. mg)) deallocate(loca)
      if (alg_adv) deallocate(ifail)
      if (alg_rrr) deallocate(isuppz)
      if (sl .and. alg_adv) deallocate(iclustr,gap)


      if (info /= 0 .and. (.not. (alg_adv .and. mod(info/2,2) /= 0)) .and. (.not. alg_el2)) &
               & call rxi('unievp: '//trim(adjustl(task))//' returned error: ', info)



      call tcx(task)

   end subroutine unievp






      subroutine unievp1(ovlp, herm, itype, jobz, rng, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, nz, w, z, ldz)
! 'universal' diagonalisation routine, serial specialisation on top of the more general unievp.
! This should have been done in the other direction. In stead of unievp1 calling unievp, unievp should be calling unievp1 in the special case of sl == .false.

      use mod_ctx
      use structures, only : unievp

      implicit none
      integer, parameter :: dp = 8

      logical, intent(in) :: ovlp, herm
      character, intent(in) :: jobz, rng, uplo
      integer, intent(in)  :: itype, n, il, iu, lda, ldb, ldz
      integer, intent(inout) :: m, nz
      real(dp), intent(in) :: vl, vu
      real(dp), intent(inout) :: w(*)
      real(dp), intent(in), target :: a(*), b(*)
      real(dp), intent(inout), target :: z(*)

      type(cart_ctx_t) :: c2d

      integer :: desca(9), descb(9), descz(9)

      call init_cart_ctx_0(c2d, 2)

      desca = [1, c2d % comm, n, n, n, n, 0, 0, lda]
      descb = [1, c2d % comm, n, n, n, n, 0, 0, ldb]
      descz = [1, c2d % comm, n, n, n, n, 0, 0, ldz]

!       desc(1) = 1
!       desc(3:4) = n
!       desc(5:6) = n
!       desc(7:8) = 0
!       desc(9) = n

      call unievp(ovlp, herm, c2d, itype, jobz, rng, uplo, n, a, 1, 1, desca, b, 1, 1, descb, &
                                                   vl, vu, il, iu, m, nz, w, z, 1, 1, descz)


      end subroutine unievp1


      subroutine dcopy2d(n,m,a,lda,b,ldb)
        implicit none

        integer, intent(in)    :: n,m,lda,ldb
        real(8), intent(in)    :: a(lda,*)
        real(8), intent(inout) :: b(ldb,*)

        integer :: i,j

        do j = 1, m
            do i = 1, n
                b(i,j) = a(i,j)
            end do
        end do

      end subroutine dcopy2d
