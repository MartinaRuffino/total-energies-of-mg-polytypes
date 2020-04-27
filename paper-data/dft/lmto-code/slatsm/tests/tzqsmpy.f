C /usr/bin/time -f "\t%E real,\t%U user,\t%S sys" a.out -count=# -setup -zqonly -zgonly
      subroutine fmain
      implicit none
      logical cmdopt,a2bin
      INTEGER   lda, ldb, ldc, ldv, ic, count, count0, isw, nmin
      INTEGER   lda2, ldb2, ldc2, ldv2
      double precision mflops,s_time,s_mflops,p_time,p_mflops,cpusec
      PARAMETER (lda = 5, ldb = 7, ldv = ldb, ldc = lda)
*     PARAMETER (lda2 = 1501, ldb2 = 5501, ldv2 = ldb2, ldc2 = lda2)
*     On phpdl1.ph.kcl.ac.uk, the following PARAMETER yields:
*     zgemm:  3.4 sec, OMP_NUM_THREADS=1, 0.95, OMP_NUM_THREADS=4
*     zqsmpy: 1.9 sec, OMP_NUM_THREADS=1, 0.62, OMP_NUM_THREADS=4
      PARAMETER (lda2 = 1501, ldb2 = 2501, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 201, ldb2 = 2501, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 51,  ldb2 = 2501, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 51,  ldb2 = 1101, ldv2 = ldb2, ldc2 = lda2)
*     On phpdl1.ph.kcl.ac.uk, the following PARAMETER yields:
*     zgemm:  0.6 msec, OMP_NUM_THREADS=1, 0.22 msec, OMP_NUM_THREADS=4
*     zqsmpy: 0.5 msec, OMP_NUM_THREADS=1, 0.30 msec, OMP_NUM_THREADS=4
C     PARAMETER (lda2 = 41,  ldb2 = 51, ldv2 = ldb2, ldc2 = lda2)
*     On phpdl1.ph.kcl.ac.uk, the following PARAMETER yields:
*     zgemm:  4.4 msec, OMP_NUM_THREADS=1, 1.4 msec, OMP_NUM_THREADS=4
*     zqsmpy: 2.1 msec, OMP_NUM_THREADS=1, 0.6 msec, OMP_NUM_THREADS=4
*     PARAMETER (lda2 = 50, ldb2 = 1501, ldv2 = ldb2, ldc2 = lda2)
*     PARAMETER (lda2 = 5, ldb2 = 7, ldv2 = ldb2, ldc2 = lda2)
C     double complex a(lda,ldb), b(ldb,lda), v(ldv,ldv), avb(lda,lda)
      character*(120) outs
      double complex av(lda,ldb),c(ldc,ldc), avt(ldb2,lda2), beta
      double complex a2(lda2,ldb2), b2(ldb2,lda2), v2(ldv2,ldv2)
      double complex av2(lda2,ldb2),c2(ldc2,ldc2),ava2(lda2,lda2)
      common /static/ a2,b2,v2,av2,c2,ava2,avt

      call dpzero(ava2,lda2*lda2*2)
      call dpzero(c2,lda2*lda2*2)
      beta = .1d0
      nmin = 128
      isw = 11
C     isw = 12

C --- Setup standard test ---
      call zqsmpy0(nmin,1)
C     call zqsmpy0(64,1)

      WRITE (*,100) lda2, lda2, lda2, ldb2, ldb2, lda2, isw, nmin
  100 format( 1X,'tzqsmpy Matrix Multiply: (',I4,',',I4,')',
     .        ' <- (',I4,',',I4,') x (',I4,',',I4,')    isw=',i2,'  nmin=',i3)

      call ran1in(1)
      call init(a2, lda2, ldb2, 'n')
C     call zprm('a',2,a2,lda2,lda2,ldb2)

      call init(v2, ldb2, ldb2, 'h')
C     call zprm('v',2,v2,ldb2,ldb2,ldb2)
C     Make a*v
      call zgemm('N','N',lda2,ldb2,ldb2,(1d0,0d0),a2,lda2,v2,ldv2,
     .  0d0,av2,lda2)
C     call zprm('av',2,av2,lda2,lda2,ldb2)

C     avt = conjugate of av
      call zmcpy('C',av2,lda2,1,avt,1,ldb2,lda2,ldb2)
C     call zprm('av+',2,avt,ldb2,ldb2,lda2)

C     b = transpose of a
      call zmcpy('C',a2,lda2,1,b2,1,ldb2,lda2,ldb2)
C     call zprm('a+',2,b2,ldb2,ldb2,lda2)
      if (cmdopt('-setup',6,0,outs)) stop 'setup'

C     Make a*va
      count0 = 2
      ic = 7
      if (cmdopt('-count=',7,0,outs)) then
        if (.not. a2bin(outs,count0,2,0,' ',ic,-1)) stop 'bad -count='
      endif
      s_time = cpusec()
      if (.not. cmdopt('-zqonly',7,0,outs)) then
      print *, 'start zgemm multiply ...'
      do  count = 1, count0
C   5 continue
C      call zgemm('N','N',lda2,lda2,ldb2,(1d0,0d0),av2,lda2,b2,ldb2,
C     .  beta,ava2,lda2)
      call zgemm('C','N',lda2,lda2,ldb2,(1d0,0d0),avt,ldb2,b2,ldb2,
     .  beta,ava2,lda2)
C      call zgemm('N','C',lda2,lda2,ldb2,(1d0,0d0),av2,lda2,a2,lda2,
C     .  beta,ava2,lda2)
C     if (cpusec() - s_time < 1d0 .and. count < 60) goto 5
      enddo
      endif
      s_time = (cpusec() - s_time)/count0
C     call zprm('ava+',2,ava2,lda2,lda2,lda2)
      p_time = cpusec()
      print 127, count0
  127 format(' executed multiply',i5,' times')

C First call std, second slowest, third fastest, fourth like std
      if (.not. cmdopt('-zgonly',7,0,outs)) then
      print *, 'start zqsmpy multiply ...'
      do  ic = 1, count0
C       print *, 'hi',ic
C       call zqsmpy(isw,'N','N',lda2,ldb2,av2,lda2,b2,ldb2,beta,c2,ldc2)
C       call zqsmpy(isw,'N','C',lda2,ldb2,av2,lda2,a2,lda2,beta,c2,ldc2)
        call zqsmpy(isw,'C','N',lda2,ldb2,avt,ldb2,b2,ldb2,beta,c2,ldc2)
C       call zqsmpy(isw,'C','C',lda2,ldb2,avt,ldb2,a2,lda2,beta,c2,ldc2)
C       call zprm('c',2,c2,ldc2,ldc2,ldc2)
C       stop
      enddo
      endif
      p_time = (cpusec() - p_time)/count0
      if (.not. cmdopt('-zgonly',7,0,outs) .and.
     .  .not. cmdopt('-zqonly',7,0,outs)) then
        CALL diff('compare to reference', c2, ava2, ldc2,ldc2,ldc2,ldb2)
      endif

      mflops = 2.0d-6 * lda2 * lda2 * ldb2 * 4
      s_mflops = mflops / s_time
      p_mflops = mflops / p_time
      print 110, s_time, s_mflops, p_time, p_mflops
  110 format(/1X,' zgemm  Time: ',F7.3,'   zgemm  MFlops: ',F7.0,
     .       /1X,' zqsmpy Time: ',F7.3,'   zqsmpy MFlops: ',F7.0)

      print 111, s_time*count0, s_time, p_time*count0, p_time
  111 format(/' time elapsed, zgemm  multiply',F7.2,
     .    ' sec (',F9.4,' sec/product)'/
     .        ' time elapsed, zqsmpy multiply',F7.2,
     .    ' sec (',F9.4,' sec/product)')

      end

      SUBROUTINE init ( a, lda, lda2, cs)
C- Initialize arrays
      implicit none
      integer lda, lda2
      double precision a( 2, lda, lda2)
      integer i,j
      character *1 cs
      real ran1
C     call ran1in(1)

      do 10 i = 1, lda
      do 10 j = 1, lda2
      a(1,i,j) = ran1()
   10 a(2,i,j) = ran1()

      if (cs == 'h') then
        do 20 i = 1, min(lda,lda2)
        do 20 j = 1, i
        a(1,i,j) =  a(1,j,i)
   20   a(2,i,j) = -a(2,j,i)
        do  22 i = 1, min(lda,lda2)
   22   a(2,i,i) = 0
      endif

      end
      SUBROUTINE diff (strn, sc, pc, ldc, n, m, l )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER ldc, n, m, l, i, j
      double complex sc( ldc, * )
      double complex pc( ldc, * )
      double precision tol,maxerr

      tol = 1e-10
      if (dble(n)*dble(m)*dble(l) > 1e8) tol = 1d-8

      DO i = 1, n
         DO j = 1, m
           maxerr = max(maxerr,ABS(sc(i,j)-pc(i,j)))
         END DO
       END DO
       IF ( maxerr .GT. tol ) THEN
         WRITE (6,100) strn, maxerr
         RETURN
       else
         WRITE (6,110) strn, tol
       endif
  100 FORMAT(1X,a,' *** ERROR ***   Arrays Have Different Results,',
     .   ' max err=',1pe8.1)
  110 FORMAT(1X,a,'... arrays agree to within tol=',1pe8.1)

C      DO i = 1, n
C         DO j = 1, m
C            IF ( ABS(sc(i,j)-pc(i,j)) .GT. tol ) THEN
C               WRITE (6,100) strn, tol
C               RETURN
C               END IF
C            END DO
C        END DO
C      WRITE (6,110) strn, tol
C      RETURN
C  100 FORMAT(1X,a,'*** ERROR ***   Arrays Have Different Results',
C     .  ' tol=',1pe8.1)
C  110 FORMAT(1X,a,'... arrays have the same results',
C     .  ' tol=',1pe8.1)
      END
