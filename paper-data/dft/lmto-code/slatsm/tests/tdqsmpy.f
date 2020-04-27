      subroutine fmain
      implicit none
      logical cmdopt,a2bin
      INTEGER   lda, ldb, ldc, ldv, ic, count, count0, isw, nmin
      INTEGER   lda2, ldb2, ldc2, ldv2
      double precision mflops,s_time,s_mflops,p_time,p_mflops,cpusec
      PARAMETER (lda = 5, ldb = 7, ldv = ldb, ldc = lda)
C     PARAMETER (lda2 = 1501, ldb2 = 2501, ldv2 = ldb2, ldc2 = lda2)
*     On MAC OS, gfortran, i7@2.7 GHZ, the following PARAMETER yields (GFlop)
*     Parallel threads phpdl1, link with thread lib and set MKL_NUM_THREADS
*              MAC(1)       henry(1)  (4) (12)  samuel(1)  (4)
*     dgemm:  23 (isw=12)     11      44   111    12       48
*     dqsmpy: 28 (isw=12)     15      55   121    16       62
*     dqsmpy: 38 (isw=11)
      PARAMETER (lda2 = 1501, ldb2 = 2501, ldv2 = ldb2, ldc2 = lda2)
*     On MAC OS, gfortran, i7@2.7 GHZ, the following PARAMETER yields
*     dgemm:  15 GFlop isw=12
*     dqsmpy: 24 GFlop isw=12
*     PARAMETER (lda2 = 201, ldb2 = 2501, ldv2 = ldb2, ldc2 = lda2)
*     On MAC OS, gfortran, i7@2.7 GHZ, the following PARAMETER yields (GFlop)
*              MAC           henry
*     dgemm:   9 (isw=12)     9 scalar  26 (4 thr)   37 (12 thr)
*     dqsmpy: 12 (isw=12)    11 scalar  27 (4 thr)   38 (12 thr)
*     PARAMETER (lda2 = 51,  ldb2 = 2501, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 51,  ldb2 = 1101, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 50, ldb2 = 1501, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 50, ldb2 = 1501, ldv2 = ldb2, ldc2 = lda2)
*     On MAC OS, gfortran, i7@2.7 GHZ, the following PARAMETER yields (GFlop)
*              MAC           henry
*     dgemm:  17 (isw=12)     9 scalar  22 (4 thr)   66 (12 thr)
*     dqsmpy: 28 (isw=12)     9 scalar  32 (4 thr)   63 (12 thr)
*     PARAMETER (lda2 = 1501, ldb2 = 50, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 5, ldb2 = 7, ldv2 = ldb2, ldc2 = lda2)
      character*(120) outs
      double precision a(lda,ldb), b(ldb,lda), v(ldv,ldv), avb(lda,lda)
      double precision av(lda,ldb),c(ldc,ldc), avt(ldb2,lda2), beta
      double precision a2(lda2,ldb2), b2(ldb2,lda2), v2(ldv2,ldv2)
      double precision av2(lda2,ldb2),c2(ldc2,ldc2),ava2(lda2,lda2)
      common /static/ a2,b2,v2,av2,c2,ava2,avt

C      data a / 1d0,  7d0, -1d0,  4d0,  5d0,
C     .         2d0, -3d0,  8d0,  0d0,  4d0,
C     .         5d0, 11d0, -3d0, -7d0,  3d0,
C     .        -6d0, -2d0,  1d0,  4d0, -2d0,
C     .         0d0,  9d0, -1d0, -2d0,  1d0,
C     .        -6d0,  5d0,  2d0,  8d0, -7d0,
C     .        -2d0, -5d0, -4d0,  3d0, -9d0/
C      data v / 17d0,  4d0, 11d0, 16d0, 12d0,  3d0,  6d0,
C     .          4d0,  4d0,  6d0,  7d0,  9d0, 10d0, 13d0,
C     .         11d0,  6d0, 11d0,  8d0,  5d0,  4d0,  3d0,
C     .         16d0,  7d0,  8d0, 19d0, 10d0, 11d0,  7d0,
C     .         12d0,  9d0,  5d0, 10d0, 17d0, 10d0, 13d0,
C     .          3d0, 10d0,  4d0, 11d0, 10d0,  7d0,  9d0,
C     .          6d0, 13d0,  3d0,  7d0, 13d0,  9d0,  8d0/
C      data avb / 1154d0,  -725d0,  -581d0, -1357d0,   749d0,
C     .           -725d0,  5288d0,  -228d0,  1378d0,   404d0,
C     .           -581d0,  -228d0,  -344d0,   474d0, -1208d0,
C     .          -1357d0,  1378d0,   474d0,  1529d0,  -887d0,
C     .            749d0,   404d0, -1208d0,  -887d0,   704d0/
C
C --- Setup for simple matrix, debugging ---
CC     b = transpose of a
C      call dmcpy(a,lda,1,b,1,ldb,lda,ldb)
CC     Make a*v
C      call dgemm('N','N',lda,ldb,ldb,1d0,a,lda,v,ldv,0d0,av,lda)
CC     call prmx('a',a,lda,lda,ldb)
CC     call prmx('av',av,lda,lda,ldb)
CC     call prmx('b',b,ldb,ldb,lda)
CC     call prmx('v',v,ldb,ldb,ldb)
CC     call prmx('avb',avb,lda,lda,lda)
C      call dqsmpy(11,'N','N',lda,ldb,av,lda,b,ldb,c,ldc)
C      CALL diff('compare to reference', c, avb, ldc, ldc, ldc, ldb )
C      stop

      call dpzero(ava2,lda2*lda2)
      call dpzero(c2,lda2*lda2)
      beta = 0.1d0
      nmin = 128
      isw = 11
      isw = 12

C --- Setup standard test ---
      call dqsmpy0(nmin,1)
C     call dqsmpy0(64,1)

      WRITE (*,100) lda2, lda2, lda2, ldb2, ldb2, lda2, isw, nmin
  100 format( 1X,'tdqsmpy Matrix Multiply: (',I4,',',I4,')',
     .        ' <- (',I4,',',I4,') x (',I4,',',I4,')    isw=',i2,'  nmin=',i3)

      call ran1in(1)
      call init(a2, lda2, ldb2, 'n')
C     call prmx('a',a2,lda2,lda2,ldb2)
      call init(v2, ldb2, ldb2, 'h')
C     call prmx('v',v2,ldb2,ldb2,ldb2)
C     Make a*v
      call dgemm('N','N',lda2,ldb2,ldb2,1d0,a2,lda2,v2,ldv2,0d0,av2,lda2)
C     call prmx('av2',av2,lda2,lda2,ldb2)

C     avt = conjugate of av
      call dmcpy(av2,lda2,1,avt,1,ldb2,lda2,ldb2)
C     call prmx('av+',avt,ldb2,ldb2,lda2)

C     b = transpose of a
      call dmcpy(a2,lda2,1,b2,1,ldb2,lda2,ldb2)
C     call prmx('a+',b2,ldb2,ldb2,lda2)

C     Make a*va
      count0 = 2
      ic = 7
      if (cmdopt('-count=',7,0,outs)) then
        if (.not. a2bin(outs,count0,2,0,' ',ic,-1)) stop 'bad -count='
      endif
      s_time = cpusec()
      if (.not. cmdopt('-zqonly',7,0,outs)) then
      print *, 'start dgemm multiply ...'
      do  count = 1, count0
        call dgemm('N','N',lda2,lda2,ldb2,1d0,av2,lda2,b2,ldb2,beta,ava2,lda2)
C       call dgemm('T','N',lda2,lda2,ldb2,1d0,avt,ldb2,b2,ldb2,beta,ava2,lda2)
C       call dgemm('N','T',lda2,lda2,ldb2,1d0,av2,lda2,a2,lda2,beta,ava2,lda2)
C       call dgemm('T','T',lda2,lda2,ldb2,1d0,avt,ldb2,a2,lda2,beta,ava2,lda2)
C       if (cpusec() - s_time < 1d0 .and. count < 60) goto 5
      enddo
      endif
      s_time = (cpusec() - s_time)/count0
C     call zprm('ava+',2,ava2,lda2,lda2,lda2)
      p_time = cpusec()
      print 127, count0
  127 format(' executed multiply',i5,' times')

C First call std, second slowest, third fastest, fourth like std
      if (.not. cmdopt('-dgonly',7,0,outs)) then
      print *, 'start dqsmpy multiply ...'
      do  ic = 1, count0
        call dqsmpy(isw,'N','N',lda2,ldb2,av2,lda2,b2,ldb2,beta,c2,ldc2)
C        call dqsmpy(isw,'T','N',lda2,ldb2,avt,ldb2,b2,ldb2,beta,c2,ldc2)
C       call dqsmpy(isw,'N','T',lda2,ldb2,av2,lda2,a2,lda2,beta,c2,ldc2)
C       call dqsmpy(isw,'T','T',lda2,ldb2,avt,ldb2,a2,lda2,beta,c2,ldc2)
C       call prmx('c',c2,ldc2,ldc2,ldc2)
C       call zprm('c',2,c2,ldc2,ldc2,ldc2)
C       stop
      enddo
      endif
      p_time = (cpusec() - p_time)/count0
      if (.not. cmdopt('-dgonly',7,0,outs) .and.
     .    .not. cmdopt('-dqonly',7,0,outs)) then
        CALL diff('compare to reference', c2, ava2, ldc2,ldc2,ldc2,ldb2)
      endif

      mflops = 2.0d-6 * lda2 * lda2 * ldb2
      s_mflops = mflops / s_time
      p_mflops = mflops / p_time
      print 110, s_time, s_mflops, p_time, p_mflops
  110 format(/1X,' dgemm  Time: ',F7.3,'   dgemm  MFlops: ',F7.0,
     .       /1X,' dqsmpy Time: ',F7.3,'   dqsmpy MFlops: ',F7.0)

      print 111, s_time*count0, s_time, p_time*count0, p_time
  111 format(/' time elapsed, dgemm  multiply',F7.2,
     .    ' sec (',F9.4,' sec/product)'/
     .        ' time elapsed, dqsmpy multiply',F7.2,
     .    ' sec (',F9.4,' sec/product)')

      end

      SUBROUTINE init ( a, lda, lda2, cs)
C- Initialize arrays
      implicit none
      integer lda, lda2
      double precision a( lda, lda2)
      integer i,j
      character *1 cs
      real ran1

      do 10 i = 1, lda
      do 10 j = 1, lda2
   10 a(i,j) = ran1()

      if (cs == 'h' .and. lda == lda2) then
        do 20 i = 1, lda
        do 20 j = 1, lda
   20   a(i,j) =  a(j,i)
      endif

      end

      SUBROUTINE diff (strn, sc, pc, ldc, n, m, l )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER ldc, n, m, l, i, j
      double precision sc( ldc, m )
      double precision pc( ldc, m )
      double precision tol,maxdiff

      tol = 1e-7
      if (dble(n)*dble(m)*dble(l) > 1e8) tol = 1d-6

      maxdiff = 0
      DO i = 1, n
        DO j = 1, m
          maxdiff = max(maxdiff,dABS( sc(i,j) - pc(i,j)))
          IF ( dABS( sc(i,j) - pc(i,j) ) .GT. tol ) THEN
            WRITE (6,100) strn, tol
            RETURN
          END IF
        END DO
      END DO
      WRITE (6,110) strn, tol, maxdiff
      RETURN
  100 FORMAT(1X,a,'*** ERROR ***   Arrays Have Different Results',
     .  ' tol=',1pe8.1)
  110 FORMAT(1X,a,'... arrays same within tol=',1pe8.1,'  max diff=',1pe8.1)
      END
