      subroutine fmain
C dgemm timings for 1500 1500 1500 matrix
C               E5-2640,v3,@2.60GHz
C MKL_NUM_THREADS 1   32623
      implicit none
      logical cmdopt,a2bin
      character outs*80
      INTEGER   lda, ldb, ldc, ldw, n, m, k, l, i, ic, count
C     integer(8) nl,ml,ll,ldal,ldbl,ldcl
      double precision mflops
*     PARAMETER (lda = 1501, ldb = 1501, ldc = 1501, ldw=1501)
      PARAMETER (lda = 5501, ldb = 5501, ldc = 5501, ldw=5501)
C     PARAMETER (lda = 401, ldb = 401, ldc = 401, ldw=401)
C      PARAMETER (n = 600, m = 600, l = 600)
*
      double precision sa( lda, lda ), sb( ldb, ldb ), sc( ldc, ldc )
      double precision pa( lda, lda ), pb( ldb, ldb ), pc( ldc, ldc ),
     .                 pw( ldw, ldw)
      double precision s_time, p_time, temp, s_mflops, p_mflops,
     .  q_time, q_mflops, walltime0
      double precision cpusec,swalltime,pwalltime,qwalltime,dwtime,delwc
      common /static/ sa, sb, sc, pa, pb, pc, pw

C --- Setup ---
      n = 300
      m = 300
      l = 300
      print *, 'n,m,l?'
      read(*,*) n,m,l
      mflops = 2.0d-6 * n * m * l
      call finits(2,0,0,i)
      call nada
      WRITE (*,100) n, m, n, l, l, m
      walltime0 = delwc()
      walltime0 = dwtime()

      CALL init ( sc,sa,sb,ldc,lda,ldb,n,m,l )
C     call ywrm(0,'a',1,6,'(4f16.10)',sa,1,lda,lda,n)
C     call ywrm(0,'b',1,6,'(4f16.10)',sb,1,ldb,l,m)
      swalltime = dwtime()
      temp = cpusec()
      count=0
      i = 7
      if (cmdopt('-count=',i,0,outs)) then
        if (.not. a2bin(outs,count,2,0,' ',i,72))
     .    call rxs('cannot parse',outs)
      endif
      if (count == 0) then
    5 CALL scalar_matrix_multiply( sc, sa, sb, ldc, lda, ldb, n, m, l )
      count = count+1
      if (cpusec() - temp < 1d0 .and. count < 10) goto 5
C     if (cpusec() - temp < 3d0 .and. count < 50) goto 5
      qwalltime = dwtime() - swalltime
      if (qwalltime /= 0 .and. qwalltime < 1) goto 5
      endif

      print 127, count
  127 format(' execute multiply',i4,' times')

C --- In-line matrix multiplication ---
      swalltime = delwc()
      print 125, swalltime
  125 format(' starting in-line call ... elapsed wall time',F9.4)
      temp = cpusec()
      do  ic = 1, count
      CALL scalar_matrix_multiply( sc, sa, sb, ldc, lda, ldb, n, m, l )
      enddo
      s_time = (cpusec() - temp)/count
      swalltime = delwc()
      if (s_time == 0d0) s_time = 1
      s_mflops = mflops / s_time
*     call ywrm(0,'c',1,6,'(4f16.10)',sc,1,ldc,n,m)
      print 128, swalltime
  128 format(' finished in-line call ... elapsed wall time',F9.4)

C --- dgemm matrix multiplication ---
      CALL init ( pc,pa,pb,ldc,lda,ldb,n,m,l )
      pwalltime = delwc()
      pwalltime = dwtime()
      temp = cpusec()
      print 129, dwtime() - walltime0
  129 format(' starting dgemm call ... total elapsed wall time',F9.4)
      do  ic = 1, count
C       call dgemm ('N','N',n,m,l,1d0,pa,lda,pb,ldb,0d0,pc,ldc)
C       call dgemm ('N','N',nl,ml,ll,1d0,pa,ldal,pb,ldbl,0d0,pc,ldcl)
       CALL fast_matrix_multiply (pc,pa,pb,ldc,lda,ldb,n,m,l)
c      call mmul(n, m, l, pa, lda, pb, ldb, pc, ldc)
      enddo
      p_time = (cpusec() - temp)/count
      pwalltime = delwc()
      print 130, pwalltime
  130 format(' finished dgemm call ... elapsed wall time',F9.4)
      if (p_time == 0d0) p_time = 1
      p_mflops = mflops / p_time
      CALL diff('compare in-line to dgemm', sc, pc, ldc, n, m, l)

C --- qmpy matrix multiplication ---
      if (cmdopt('-qmpy',5,0,outs)) then
      CALL init ( pc,pa,pb,ldc,lda,ldb,n,m,l )
      qwalltime = dwtime()
      temp = cpusec()
      do  9  ic = 1, count
    9 call qmpy('N','N',n,m,l,pw,ldw,pa,lda,pb,ldb,pc,ldc)
      q_time = (cpusec() - temp)/count
      CALL diff('compare in-line to qmpy ', sc, pc, ldc, n, m, l )
      if (q_time == 0d0) q_time = 1
      q_mflops = mflops / q_time
      qwalltime = dwtime() - qwalltime
      print 131, dwtime() - walltime0
  131 format(' finished qmpy call ... total elapsed wall time',F9.4)
      print 110, s_time, s_mflops, swalltime, swalltime/count,
     .           p_time, p_mflops, pwalltime, pwalltime/count,
     .           q_time, q_mflops, qwalltime
      else
      print 110, s_time, s_mflops, swalltime, swalltime/count,
     .           p_time, p_mflops, pwalltime, pwalltime/count
      endif

  100 format( 1X,'Matrix Multiply: ',I4,'x',I4,
     .        ' <- ',I4,'x',I4,' x ',I4,'x',I4)
  110 format(/1X,'Serial CPU time: ',F7.3,'  Serial MFlops: ',F7.0,
     .                 '   Wall time: ',F7.3,' (',F9.4,' /product)'
     .       /1X,' dgemm CPU time: ',F7.3,'   dgemm MFlops: ',F7.0,
     .                 '   Wall time: ',F7.3,' (',F9.4,' /product)':
     .       /1X,'  qmpy CPU time: ',F7.3,'    qmpy MFlops: ',F7.0
     .                 '   Wall time: ',F7.3)

      end

      SUBROUTINE init ( c, a, b, ldc, lda, ldb, n, m, l )
C- Initialize arrays
      integer lda, ldb, ldc, n, m, l
      double precision a( lda, * ), b( ldb, * ), c( ldc, * )
*
      do i = 1, n
         do j = 1, m
            c(i,j) = dble( i+j ) * 0d0
            end do
         end do
*
      do i = 1, n
         do j = 1, l
            a(i,j) = dble( i-j )
            end do
         end do
*
      do i = 1, l
         do j = 1, m
            b(i,j) = dble( i )/j
            end do
         end do
      end

      SUBROUTINE diff (strn, sc, pc, ldc, n, m, l )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER ldc, n, m, l, i, j
      double precision sc( ldc, * )
      double precision pc( ldc, * )
      double precision tol,maxerr

      tol = 1e-6
      if (dble(n)*dble(m)*dble(l) > 1e8) tol = 1d-5
      maxerr = 0

      DO i = 1, n
         DO j = 1, m
           maxerr = max(maxerr,dABS(sc(i,j)-pc(i,j)))
         END DO
       END DO
       IF ( maxerr .GT. tol ) THEN
         WRITE (6,100) strn, maxerr
         RETURN
       else
         WRITE (6,110) strn, tol
       endif
  100 FORMAT(1X,a,'*** ERROR ***   Arrays Have Different Results,',
     .   ' max err=',1pe8.1)
  110 FORMAT(1X,a,'... arrays agree to within tol=',1pe8.1)
      END
      SUBROUTINE fast_matrix_multiply(c, a, b, ldc, lda, ldb, n, m, l)
C- Matrix Multiplication with dgemm
      implicit none
      integer lda, ldb, ldc, n, m, l
      double precision a (lda,*), b(ldb,*), c(ldc,*)

C     Use this call for 32 bit integers
      call dgemm ('N','N',n,m,l,1d0,a,lda,b,ldb,0d0,c,ldc)

C     Use this call for 64 bit integers
C     integer(8) nl,ml,ll,ldal,ldbl,ldcl
C      nl = n
C      ml = m
C      ll = l
C      ldal = lda
C      ldbl = ldb
C      ldcl = ldc
C      call dgemm ('N','N',nl,ml,ll,1d0,a,ldal,b,ldbl,0d0,c,ldcl)
      end

      SUBROUTINE scalar_matrix_multiply(c, a, b, ldc, lda, ldb, n, m, l)
C- In-line matrix multipliation
      integer lda, ldb, ldc, n, m, l
      double precision a (lda,*), b(ldb,*), c(ldc,*)
      do j = 1, m
        do i = 1, n
          c(i,j) = 0
        end do
      end do
      do j = 1, m
        do k = 1, l
          do i = 1, n
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do
      end
      subroutine qmpy(transa,transb,m,n,k,w,ldw,a,lda,b,ldb,c,ldc)
C- Quick matrix multiply
C  See Numerical Recipes, 2.11.
C  In practice, at best 10% faster than dgemm for SGI, HP, AIX
C  The same result on 4x4 matrices a,b
C  mc -f4f16.10 a -split a 1,3,5 1,3,5 b -split b 1,3,5 1,3,5 -x\
C  a12 a22 -- b21 b22 -+ -x -a q7  a11 a22 -+ b11 b22 -+ -x -a q1\
C  a21 a22 -+ b11 -x -a q2 a11 b12 b22 -- -x  -a q3\
C  a22 b21 b11 -- -x -a q4  a11 a12 -+ b22 -x -a q5\
C  a21 a11 -- b11 b12 -+ -x -a q6  q1 q4 -+ q5 -- q7 -+ -a c11\
C  q2 q4 -+ -a c21 q3 q5 -+ -a c12 q1 q3 -+ q2 -- q6 -+ -a c22\
C  c11 c12 -ccat c21 c22 -ccat -rcat
      implicit none
      character*1 transa, transb
      integer m,n,k,lda,ldb,ldc,ldw
      double precision a(lda,1),b(ldb,1),c(ldc,1),w(ldw,1)
      integer mmin,npass,m1,m2,i,j
      parameter(mmin=4)
C ... debugging
      character*10 fmt
      integer ip
      data ip /0/
      data fmt /'(4f16.10)'/
      save ip

      if (m /= n .or. m /= k)  then
        if (ip == 0) print *, 'qmpy m must be n,k for now'
        ip = 1
        return
      endif

      if ((transa /= 'n' .and. transa /= 'N') .or.
     .    (transb /= 'n' .and. transb /= 'N'))
     .  stop 'qmpy: not implemented trans'

      npass = dlog(dble(m)/mmin)/dlog(2d0) + 1d-9
      if (npass == 0) then
        call dgemm(transa,transb,m,n,k,1d0,a,lda,b,ldb,0d0,c,ldc)
        return
      endif

      m1 = m/2
      m2 = 2*m1

C --- Q6 in c22 ---
      do  10  j = 1, m1
      do  10  i = 1, m1
      w(i,j)    = a(m1+i,j) - a(i,j)
   10 w(i+m1,j) = b(i,j) + b(i,j+m1)
*     call ywrm(0,'a21-a11',1,6,fmt,w,1,ldw,m1,m1)
*     call ywrm(0,'b11+b12',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
     .  c(1+m1,1+m1),ldc)
*     call ywrm(0,'Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)

C --- Q7 in c11 ---
      do  20  j = 1, m1
      do  20  i = 1, m1
      w(i,j)    = a(i,m1+j) - a(m1+i,m1+j)
   20 w(i+m1,j) = b(m1+i,j) + b(m1+i,j+m1)
*     call ywrm(0,'a12-a22',1,6,fmt,w,1,ldw,m1,m1)
*     call ywrm(0,'b21+b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
     .  c,ldc)
*     call ywrm(0,'Q7',1,6,fmt,c,1,ldc,m1,m1)

C --- Q1 in c12 ---
      do  30  j = 1, m1
      do  30  i = 1, m1
      w(i,j)    = a(i,j) + a(i+m1,j+m1)
   30 w(i+m1,j) = b(i,j) + b(i+m1,j+m1)
*     call ywrm(0,'a11+a22',1,6,fmt,w,1,ldw,m1,m1)
*     call ywrm(0,'b11+b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
     .  c(1,1+m1),ldc)
*     call ywrm(0,'Q1',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)

C ... c11 += Q1 and c22 += Q1
      do  40  j = 1, m1
      do  40  i = 1, m1
      c(i,j)       = c(i,j)       + c(i,j+m1)
   40 c(i+m1,j+m1) = c(i+m1,j+m1) + c(i,j+m1)
*     call ywrm(0,'Q1+Q7',1,6,fmt,c(1,1),1,ldc,m1,m1)
*     call ywrm(0,'Q1+Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)

C --- Q2 in c21, Q4 in c12 ---
      do  50  j = 1, m1
      do  50  i = 1, m1
      w(i,j)    = a(i+m1,j) + a(i+m1,j+m1)
   50 w(i+m1,j) = b(i+m1,j) - b(i,j)
*     call ywrm(0,'a21+a22',1,6,fmt,w,1,ldw,m1,m1)
*     call ywrm(0,'b21-b11',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,b,ldb,0d0,
     .  c(1+m1,1),ldc)
      call dgemm(transa,transb,m1,m1,m1,1d0,a(1+m1,1+m1),lda,
     .  w(1+m1,1),ldw,0d0,c(1,1+m1),ldc)
*     call ywrm(0,'Q2',1,6,fmt,c(1+m1,1),1,ldc,m1,m1)
*     call ywrm(0,'Q4',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)

C ... c11 += Q4, c22 -= Q2, c21 = Q2+Q4 (watch optimizer here)
      do  60  j = 1, m1
      do  60  i = 1, m1
      c(i,j)       = c(i,j)       + c(i,j+m1)
      c(i+m1,j+m1) = c(i+m1,j+m1) - c(i+m1,j)
      c(i+m1,j) = c(i+m1,j) + c(i,j+m1)
   60 continue
*     call ywrm(0,'Q1+Q7+Q4',1,6,fmt,c(1,1),1,ldc,m1,m1)
*     call ywrm(0,'Q1+Q6-Q2',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)
*     call ywrm(0,'c21=Q2+Q4',1,6,fmt,c(1+m1,1),1,ldc,m1,m1)

C --- Q5 in c12, Q3 in w ---
      do  70  j = 1, m1
      do  70  i = 1, m1
      w(i,j)    = a(i,j) + a(i,j+m1)
   70 w(i+m1,j) = b(i,j+m1) - b(i+m1,j+m1)
*     call ywrm(0,'a11+a12',1,6,fmt,w,1,ldw,m1,m1)
*     call ywrm(0,'b12-b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,b(1+m1,1+m1),ldb,0d0,
     .  c(1,1+m1),ldc)
      call dgemm(transa,transb,m1,m1,m1,1d0,a,lda,
     .  w(1+m1,1),ldw,0d0,w,ldw)
*     call ywrm(0,'Q5',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)
*     call ywrm(0,'Q3',1,6,fmt,w,1,ldw,m1,m1)

C ... c11 -= Q5, c22 += Q3, c12 = Q3+Q5
      do  80  j = 1, m1
      do  80  i = 1, m1
      c(i,j)       = c(i,j)       - c(i,j+m1)
      c(i+m1,j+m1) = c(i+m1,j+m1) + w(i,j)
   80 c(i,j+m1) = c(i,j+m1) + w(i,j)

*     call ywrm(0,'c',1,6,fmt,c,1,ldc,m2,m2)
      end
