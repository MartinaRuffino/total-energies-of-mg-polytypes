      subroutine pzhev(linv,lc,lov,n,oh,os,nb,nprow,npcol,emx,nmx,
     .                 nev,e,ot)
C- MPI parallel diagonaliser (also inverts a real symmetric matrix)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   linv: invert a real symmetric matrix
Ci   lc:  true if complex
Ci   lov: true if non orthogonal
Ci   n :  dimension
Ci   oh:  pointer to h allocated from heap in calling program
Ci   os:  pointer to s allocated from heap in calling program
Ci        not referenced if lov=F
Ci   nb,nprow,npcol: BLACS process configuration
Ci                   (created here if nprow=-1)
Ci   emx,nmx,nev: as usual, see zhev
Co Outputs:
Co   e : eigenvalues
Co   ot: pointer to eigenvectors allocated here
Cr Remarks
Cr   pzhev needs to allocate local arrays from the heap which are
Cr   passed to PZHEGVX in place of h, o, and t. This can be done without
Cr   additional memory as follows. On entry oh and os are allocated but
Cr   not ot. An array oa is allocated and assigned after which it is
Cr   copied back into the heap at address oh. oa is then released and
Cr   ob allocated which is assigned and then copied back at the address
Cr   os. ob is released. Then a local array is allocated at oz; and oh,
Cr   os and oz passed to PZHEGVX. On exit the local arrays at oz have to
Cr   be assembled and returned at the address ot. However we don't need
Cr   oh or os anymore; so the local arrays at oz are copied back to oh
Cr   and oz is released and then ot is allocated. Finally the local
Cr   eigenvector arrays now at oh are distributed into the global
Cr   eigenvector array at ot.
Cr
Cr   Process configuration: nb is a blocking factor; the processes can
Cr   be built into an nprow X npcol array (as long
Cr   as nprow*npcol=numprocs). This can speed up PZHEGVX on some
Cr   architectures (see http://www.netlib.org/blacs/BLACS/QRef.html).
Cr   By default, if nprow=-1 on entry pzhev calls MPI_DIMS_CREATE to
Cr   find determine the configuration.
Cu Updates
Cu   25 Apr 04 (Kirill Belashchenko) workaround scalapack bug
C ----------------------------------------------------------------------
      use mpi
      implicit none
      integer procid, master, numprocs, ierr
C     integer status(MPI_STATUS_SIZE)
      integer MAX_PROCS,dims(2)
      parameter (MAX_PROCS = 1028)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision sdtime, edtime, sctime, ectime
      logical mlog,cmdopt
      character*120 strn
C Passed
      logical linv,lov,lc
      integer nmx,nev,n
      double precision emx
C BLACS process configuration
      integer nb,nprow,npcol
C E-vals (output)
      double precision e(n)
C Pointers to H, S (input) and Z (output)
      integer oh,os,ot
C Pointers to local distributed arrays
      integer oa,ob,oz
C Work arrays
      integer lrwork,lwork,liwork,ifail
      integer owork,orwork,oiwork,oifail
C Work array sizes
      double complex swork
      double precision srwork(3),work(3)

C Local
      double precision zero,VL,VU
      parameter (zero = 0d0)
      integer BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     .        MB_, NB_, RSRC_, CSRC_, LLD_
      parameter ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     .          CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     .          RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      integer context, i, iam, ibtype, info, m, mycol, myrow,
     .        nprocs, nz,
     .        NB_A, NB_B, NB_Z, CSRC_A, CSRC_B, CSRC_Z,
     .        lda, ldb, ldz, mda, mdb, mdz
      character jobz, range
      double precision abstol, d1mach, PDLAMCH
      integer desca(DLEN_), descb(DLEN_), descz(DLEN_),
     .        iclustr( MAX_PROCS*2 )
      double precision gap( MAX_PROCS )
      integer IU
      integer lgunit, numroc
      external blacs_exit, blacs_get, blacs_gridexit,
     .         blacs_gridinfo, blacs_gridinit, blacs_pinfo,
     .         blacs_setup, descinit, pzhegvx, pzlaprnt
C ... Heap
      integer w(1)
      common /w/ w
C length of real / complex number
      integer lrc
      if (lc) then
        lrc = 2
      else
        lrc = 1
      endif

      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs == 1) return
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
C MPI process configuration
      if (nprow == -1) then
        dims(1) = 0
        dims(2) = 0
        call MPI_DIMS_CREATE(numprocs,2,dims,ierr)
        npcol = dims(1)
        nprow = dims(2)
        if (mlog) then
          call awrit2(
     .      ' MPI creating process configuration .. nprow=%i npcol=%i',
     .      ' ',256,lgunit(3),nprow,npcol)
        endif
      endif
C     Initialize the BLACS
      call blacs_pinfo( iam, nprocs )
      if ( nprocs < 1 ) then
         call blacs_setup( iam, nprow*npcol )
      endif
      if (mlog) then
        call gettime(datim)
        call awrit6(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' initialising BLACS; nprow=%i npcol=%i iam=%i nprocs=%i',
     .    ' ',256,lgunit(3),procid,numprocs,nprow,npcol,iam,nprocs)
        call ftflsh(-1)
      endif
C     Initialize a single BLACS context
      call blacs_get( -1, 0, context )
      call blacs_gridinit( context, 'r', nprow, npcol )
      call blacs_gridinfo( context, nprow, npcol, myrow, mycol )
C     Bail out if this process is not a part of this context.
      if ( myrow == -1 ) then
        if (mlog) then
          call gettime(datim)
          call awrit2(' pzhev '//datim//' Process %i of %i on '
     .       //shortname(procid)(1:namelen(procid))//
     .       ' is not in context, aborting ..',' ',256,lgunit(3),
     .        procid,numprocs)
          call ftflsh(-1)
        endif
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' is not in context, aborting ..',' ',256,lgunit(1),
     .    procid,numprocs)
        call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
        call fexit(0,0,' ',0)
      endif
C     These are basic array descriptors
      call descinit( desca, n, n, nb, nb, 0, 0, context, n, info )
      call descinit( descz, n, n, nb, nb, 0, 0, context, n, info )
      call descinit( descb, n, n, nb, nb, 0, 0, context, n, info )
C Get dimensions of local matrices a, b and z
      lda = desca( LLD_ )
      ldb = descb( LLD_ )
      ldz = descz( LLD_ )
      NB_A = desca ( NB_ )
      NB_B = descb ( NB_ )
      NB_Z = descz ( NB_ )
      CSRC_A = desca ( CSRC_ )
      CSRC_B = descb ( CSRC_ )
      CSRC_Z = descz ( CSRC_ )
      mda = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL )
      mdb = NUMROC( N, NB_B, MYCOL, CSRC_B, NPCOL )
      mdz = NUMROC( N, NB_Z, MYCOL, CSRC_Z, NPCOL )
      if (mlog) then
        call gettime(datim)
        if (lov) then
          call awrit7(' pzhev '//datim//
     .                ' getting local matrix dimensions:%N'//
     .                '   n=%i '//
     .                ' a:(%ix%i)'//
     .                ' b:(%ix%i)'//
     .                ' z:(%ix%i)'//
     .                ' allocating from heap ..',' ',256,lgunit(3),
     .                n,lda,mda,ldb,mdb,ldz,mdz)
        else
          call awrit5(' pzhev '//datim//
     .                ' getting local matrix dimensions:%N'//
     .                '   n=%i '//
     .                ' a:(%ix%i)'//
     .                ' z:(%ix%i)'//
     .                ' allocating from heap ..',' ',256,lgunit(3),
     .                n,lda,mda,ldz,mdz)
        endif
        call ftflsh(-1)
      endif
C Distribute h and s into local arrays
      if (lc) then
        call defcc(oa, lda*mda)
        call dstmtz(desca,n,w(oh),w(oa))
      else
        call defrr(oa, lda*mda)
        call dstmtd(desca,n,w(oh),w(oa))
      endif
      call dcopy(lrc*lda*mda,w(oa),1,w(oh),1)
      call rlse(oa)
      if (lov) then
        if (lc) then
          call defcc(ob, ldb*mdb)
          call dstmtz(descb,n,w(os),w(ob))
        else
          call defrr(ob, ldb*mdb)
          call dstmtd(descb,n,w(os),w(ob))
        endif
        call dcopy(lrc*ldb*mdb,w(ob),1,w(os),1)
        call rlse(ob)
      endif
      if (.not. linv) then
        if (lc) then
          call defcc(oz, ldz*mdz)
        else
          call defrr(oz, ldz*mdz)
        endif
      endif
      ibtype = 1
      if (nmx <= 0) then
        jobz = 'N'
        range = 'V'
        VL = -1d12
        VU = emx
        IU = n
      else
        jobz = 'V'
        range = 'I'
        VL = zero
        VU = zero
        IU = min(n,nmx)
      endif
      abstol = 2.0*PDLAMCH(context,'u')
      call defi (oifail,  n)
C Workspace query (use liwork as temporary pointer to iwork ..)
      if (.not. linv) then
        if (lov) then
          if (lc) then
            call PZHEGVX(ibtype,jobz,range,'U',n,w(oh),1,1,desca,w(os),
     .                   1,1,descb,VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,swork,-1,srwork,-1,
     .                   liwork,-1,w(oifail),iclustr,gap,info)
          else
            call PDSYGVX(ibtype,jobz,range,'U',n,w(oh),1,1,desca,w(os),
     .                   1,1,descb,VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,work,-1,liwork,
     .                   -1,w(oifail),iclustr,gap,info)
          endif
        else
          if (lc) then
            call PZHEEVX(jobz,range,'U',n,w(oh),1,1,desca,
     .                   VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,swork,-1,srwork,-1,
     .                   liwork,-1,w(oifail),iclustr,gap,info)
          else
            call PDSYEVX(jobz,range,'U',n,w(oh),1,1,desca,
     .                   VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,work,-1,liwork,
     .                   -1,w(oifail),iclustr,gap,info)
          endif
        endif
        if (lc) then
          lwork = int(swork)
          lrwork = int(srwork(1))
          call defcc(owork,   lwork)
          call defrr(orwork,  lrwork)
        else
          lwork = int(work(1))
          call defrr(owork,   lwork)
        endif
        call defi (oiwork,  liwork)
        if (mlog) then
          if (lc) then
            call gettime(datim)
            call awrit3(' pzhev '//datim/
     .        /' Optimal scalapack worksizes:'/
     .        /'%N   lwork=%i lrwork=%i liwork=%i. '/
     .        /' Allocated from heap.',' ',256,lgunit(3),lwork,lrwork,
     .        liwork)
          else
            call gettime(datim)
            call awrit2(' pzhev '//datim/
     .        /' Optimal scalapack worksizes:'/
     .        /'%N   lwork=%i liwork=%i. '//' Allocated from heap.',' ',
     .        256,lgunit(3),lwork,liwork)
          endif
          call ftflsh(-1)
        endif
      endif
      if (procid == master .and. mlog) then
        sdtime = MPI_WTIME()
      endif
      if (linv) then
C Invert
        call PDPOTRF('U',n,w(oh),1,1,desca,info)
        if (info /= 0 .and. procid == master) then
          call awrit1(' **** in pzhev, PDPOTRF returned info=%i',
     .                ' ',128,lgunit(1),info)
        endif
        call PDPOTRI('U',n,w(oh),1,1,desca,info)
        if (info /= 0 .and. procid == master) then
          call awrit1(' **** in pzhev, PDPOTRI returned info=%i',
     .                ' ',128,lgunit(1),info)
        endif
      else
C Diagonalise
        if (lov) then
          if (lc) then
            call PZHEGVX(ibtype,jobz,range,'U',n,w(oh),1,1,desca,w(os),
     .                   1,1,descb,VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,w(owork),lwork,w(orwork),
     .                   lrwork,w(oiwork),liwork,w(oifail),iclustr,gap,
     .                   info)
          else
            call PDSYGVX(ibtype,jobz,range,'U',n,w(oh),1,1,desca,w(os),
     .                   1,1,descb,VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,w(owork),lwork,w(oiwork),
     .                   liwork,w(oifail),iclustr,gap,info)
          endif
        else
          if (lc) then
            call PZHEEVX(jobz,range,'U',n,w(oh),1,1,desca,
     .                   VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,w(owork),lwork,w(orwork),
     .                   lrwork,w(oiwork),liwork,w(oifail),iclustr,gap,
     .                   info)
          else
            call PDSYEVX(jobz,range,'U',n,w(oh),1,1,desca,
     .                   VL,VU,1,IU,abstol,m,nz,e,
     .                   zero,w(oz),1,1,descz,w(owork),lwork,w(oiwork),
     .                   liwork,w(oifail),iclustr,gap,info)
          endif
        endif
        if (info /= 0 .and. procid == master) then
          if (linv) then
          else
            if (lov) then
              if (lc) then
                call awrit1(' **** in pzhev, PZHEGVX returned info=%i',
     .                      ' ',128,lgunit(1),info)
              else
                call awrit1(' **** in pzhev, PDSYGVX returned info=%i',
     .                      ' ',128,lgunit(1),info)
              endif
            else
              if (lc) then
                call awrit1(' **** in pzhev, PZHEEVX returned info=%i',
     .                      ' ',128,lgunit(1),info)
              else
                call awrit1(' **** in pzhev, PDSYEVX returned info=%i',
     .                      ' ',128,lgunit(1),info)
              endif
            endif
          endif
        endif
      endif
      if (procid == master .and. mlog) then
        edtime = MPI_WTIME()
      endif
      nev = nz
      call rlse(oifail)
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' is at the barrier',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if (procid == master .and. mlog) then
        sctime = MPI_WTIME()
      endif
C Poke distributed array into t (use heap location oh for temp)
      if (.not. linv) then
        call dcopy(lrc*ldz*mdz,w(oz),1,w(oh),1)
        if (mlog) then
          call gettime(datim)
          call awrit2(' pzhev '//datim//' Process %i of %i on '
     .                //shortname(procid)(1:namelen(procid))//
     .                ' poked oz to oh',' ',256,lgunit(3),
     .                procid,numprocs)
          call ftflsh(-1)
        endif
        call rlse(oz)
      endif
      if (lc) then
        call defcc(ot, n*n)
      else
        call defrr(ot, n*n)
      endif
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' ready to distribute ot',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
      if (lc) then
        call udstmtz(descz,n,w(oh),w(ot))
      else
        call udstmtd(descz,n,w(oh),w(ot))
      endif
      if (mlog) then
        call gettime(datim)
        call awrit2(' pzhev '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' done distribute ot',' ',256,lgunit(3),
     .    procid,numprocs)
        call ftflsh(-1)
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if (linv) then
        call symm(n,w(ot))
      endif
C Don't do this!
C      call blacs_gridexit(context)
C      call blacs_exit(1)
C#ifdefC SUN
C      call ieee_flags( 'clear', 'exception', 'underflow', '')
C#endif
      if (procid == master .and. mlog) then
        ectime = MPI_WTIME()
        call gettime(datim)
        call awrit2(' pzhev '//datim//' walltimes: compute %gs,'//
     .    ' distribute %gs.',' ',256,lgunit(3),ectime-sctime,
     .    edtime-sdtime)
      endif
      end

      subroutine dstmtz(desc,n,ag,al)
C Distribute global matrix ag into local matrix al
      implicit none
      integer desc(1),n
      double complex ag(n,n),al(n,n)
      integer i,j
      do  i = 1, n
        do  j = 1, n
           CALL PZELSET( al, i, j, desc, ag(i,j))
         enddo
       enddo
       end
      subroutine udstmtz(desc,n,al,ag)
C Undistribute local matrix al into global matrix ag
      implicit none
      integer desc(1),n
      double complex al(n,n),ag(n,n)
      integer i,j
      double complex alpha
      do  i = 1, n
        do j = 1, n
          call PZELGET( 'A', ' ', alpha, al, i, j, desc)
          ag(i,j) = alpha
        enddo
      enddo
      end
      subroutine dstmtd(desc,n,ag,al)
C Distribute global matrix ag into local matrix al
      implicit none
      integer desc(1),n
      double precision ag(n,n),al(n,n)
      integer i,j
      do  i = 1, n
        do  j = 1, n
          CALL PDELSET( al, i, j, desc, ag(i,j))
        enddo
      enddo
      end
      subroutine udstmtd(desc,n,al,ag)
C Undistribute local matrix al into global matrix ag
      implicit none
      integer desc(1),n
      double precision al(n,n),ag(n,n)

      integer i,j
      double precision alpha
      do  i = 1, n
        do j = 1, n
          call PDELGET( 'A', ' ', alpha, al, i, j, desc)
          ag(i,j) = alpha
        enddo
      enddo
      end
      subroutine symm(n,a)
C Copy upper triangle into lower triangle
      implicit none
      integer n
      double precision a(n,n)
      integer i, j
      do  i = 2, n
        do j = 1, i-1
          a(i,j) = a(j,i)
        enddo
      enddo
      end
