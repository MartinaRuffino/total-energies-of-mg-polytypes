      subroutine fmain
      implicit none
      logical cmdopt,a2bin
      character outs*80
      INTEGER   lda, ldb, ldc, ldw, n, m, i, ic, count, omppid
      double precision mflops
C      PARAMETER (lda = 1501, ldb = 1501, ldc = 1501, ldw=1501)
*      PARAMETER (lda = 5501, ldb = 5501, ldc = 5501, ldw=5501)
      PARAMETER (lda = 401, ldb = 401, ldc = 401, ldw=401)
C      PARAMETER (n = 600, m = 600, l = 600)

      complex(8):: sa( lda*lda ), sb( ldb*ldb ), sc( ldc* ldc )
      complex(8):: pa( lda* lda ), pb( ldb* ldb ), pc( ldc* ldc ),
     .             pw( ldw* ldw)
      complex(8):: zdotu
      double precision s_time, p_time, temp, s_mflops, p_mflops,
     .  walltime0
      double precision cpusec,swalltime,pwalltime,qwalltime,dwtime,delwc
      common /static/ sa, sb, sc, pa, pb, pc, pw

C --- Setup ---
      n = 300
      m = 300
      print *, 'n,m?'
      read(*,*) n,m
      mflops = 2.0d-6 * n * m * 4
      call finits(2,0,0,i)
      call nada
      WRITE (*,100) n, m
  100 format( 1X,'zdotu: n,m',2I6)

      walltime0 = delwc()
      walltime0 = dwtime()

C --- In-line matrix multiplication ---
      CALL init ( sc,sa,sb,ldc,lda,ldb,n,m )
C     call ywrm(0,'a',1,6,'(4f16.10)',sa,1,lda,lda,n)
C     call ywrm(0,'b',1,6,'(4f16.10)',sb,1,ldb,l,m)
      swalltime = dwtime()
      temp = cpusec()
      count = 0
      i = 7
      if (cmdopt('-count=',i,0,outs)) then
        if (.not. a2bin(outs,count,2,0,' ',i,72))
     .    call rxs('cannot parse',outs)
      endif
      if (count == 0) then
    5 CALL scalar_zdotu( sc, sa, sb, ldc*ldc, lda*lda, ldb*ldb, n, m )
      count = count+1
      if (cpusec() - temp < 1d0 .and. count < 100000) goto 5
C     if (cpusec() - temp < 3d0 .and. count < 50) goto 5
      qwalltime = dwtime() - swalltime
      if (qwalltime /= 0 .and. qwalltime < 1) goto 5
      endif

      print 127, count
  127 format(' execute multiply',i7,' times')

C --- In-line matrix multiplication ---
      CALL init ( sc,sa,sb,ldc,lda,ldb,n,m )
      swalltime = delwc()
      print 125, swalltime
  125 format(' starting in-line call ... elapsed wall time',F9.4)
      temp = cpusec()
      do  ic = 1, count
      CALL scalar_zdotu(sc, sa, sb, ldc, lda, ldb, n, m )
      enddo
      s_time = (cpusec() - temp)/count
      swalltime = delwc()
      if (s_time == 0d0) s_time = 1
      s_mflops = mflops / s_time
*     call ywrm(0,'c',1,6,'(4f16.10)',sc,1,ldc,n,m)
      print 128, swalltime
  128 format(' finished in-line call ... elapsed wall time',F9.4)

C --- zgemm matrix multiplication ---
      CALL init ( pc,pa,pb,ldc,lda,ldb,n,m )
      pwalltime = delwc()
      pwalltime = dwtime()
      temp = cpusec()
      print 129, dwtime() - walltime0
  129 format(' starting zdotu call ... total elapsed wall time',F9.4)
      do  ic = 1, count
        pb(1) = zdotu(n*m,pa,1,pc,1)
      enddo
c      call mmul(n, m, l, pa, lda, pb, ldb, pc, ldc)
      p_time = (cpusec() - temp)/count
      pwalltime = delwc()
      print 130, pwalltime
  130 format(' finished zdotu call ... elapsed wall time',F9.4)
      if (p_time == 0d0) p_time = 1
      p_mflops = mflops / p_time
      CALL diff('compare in-line to zdotu', sb, pb, ldc, n, m)

      print 110, s_time, s_mflops, swalltime, swalltime/count,
     .           p_time, p_mflops, pwalltime, pwalltime/count

  110 format(/1X,'Serial CPU time: ',F7.3,'  Serial MFlops: ',F7.0,
     .                 '   Wall time: ',F7.3,' (',f14.8,' /product)'
     .       /1X,' zdotu CPU time: ',F7.3,'   zdotu MFlops: ',F7.0,
     .                 '   Wall time: ',F7.3,' (',f14.8,' /product)':
     .       /1X,'  qmpy CPU time: ',F7.3,'    qmpy MFlops: ',F7.0
     .                 '   Wall time: ',F7.3)

      i = omppid(0)
      if (i /= 0) print *, 'OMP_NUM_THREADS =',i


      end

      SUBROUTINE init ( c, a, b, ldc, lda, ldb, n, m )
C- Initialize arrays
      integer lda, ldb, ldc, n, m
      double complex a(lda), b(ldb), c(ldc)
*
      do i = 1, n*m
        c(i) = dble( i+j )
      end do
*
      do i = 1, n*m
        a(i) = dble( i-j ) + (0d0,1d0)
      end do

      b(1) = 2

      end

      SUBROUTINE diff (strn, sc, pc, ldc, n, m )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER ldc, n, m, i
      double complex sc( ldc)
      double complex pc( ldc)
      double precision tol,maxerr

      tol = 1d-6
C      if (dble(n)*dble(m) > 1e8) tol = 1d-5
C      maxerr = 0

      DO i = 1, 1
        maxerr = max(maxerr,ABS(sc(i)-pc(i)))
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

      subroutine scalar_zdotu(c, a, b, ldc, lda, ldb, n, m)
C- In-line matrix zdotu
      integer lda, ldb, ldc, n, m, nm
      double complex a (lda), b(ldb), c(ldc), sum

      nm = n*m

      sum = 0
C$OMP parallel do reduction(+:sum)
      do i = 1, nm
        sum = sum + (a(i) * c(i))
      enddo
C$OMP end parallel do
      b(1) = sum

      end
C      subroutine qmpy(transa,transb,m,n,k,w,ldw,a,lda,b,ldb,c,ldc)
CC- Quick matrix multiply
CC  See Numerical Recipes, 2.11.
CC  In practice, at best 10% faster than dgemm for SGI, HP, AIX
CC  The same result on 4x4 matrices a,b
CC  mc -f4f16.10 a -split a 1,3,5 1,3,5 b -split b 1,3,5 1,3,5 -x\
CC  a12 a22 -- b21 b22 -+ -x -a q7  a11 a22 -+ b11 b22 -+ -x -a q1\
CC  a21 a22 -+ b11 -x -a q2 a11 b12 b22 -- -x  -a q3\
CC  a22 b21 b11 -- -x -a q4  a11 a12 -+ b22 -x -a q5\
CC  a21 a11 -- b11 b12 -+ -x -a q6  q1 q4 -+ q5 -- q7 -+ -a c11\
CC  q2 q4 -+ -a c21 q3 q5 -+ -a c12 q1 q3 -+ q2 -- q6 -+ -a c22\
CC  c11 c12 -ccat c21 c22 -ccat -rcat
C      implicit none
C      character*1 transa, transb
C      integer m,n,k,lda,ldb,ldc,ldw
C      double complex a(lda,1),b(ldb,1),c(ldc,1),w(ldw,1)
C      integer mmin,npass,m1,m2,i,j
C      parameter(mmin=4)
CC ... debugging
C      character*10 fmt
C      integer ip
C      data ip /0/
C      data fmt /'(4f16.10)'/
C      save ip
C
C      if (m /= n .or. m /= k)  then
C        if (ip == 0) print *, 'qmpy m must be n,k for now'
C        ip = 1
C        return
C      endif
C
C      if ((transa /= 'n' .and. transa /= 'N') .or.
C     .    (transb /= 'n' .and. transb /= 'N'))
C     .  stop 'qmpy: not implemented trans'
C
C      npass = dlog(dble(m)/mmin)/dlog(2d0) + 1d-9
C      if (npass == 0) then
C        call dgemm(transa,transb,m,n,k,1d0,a,lda,b,ldb,0d0,c,ldc)
C        return
C      endif
C
C      m1 = m/2
C      m2 = 2*m1
C
CC --- Q6 in c22 ---
C      do  10  j = 1, m1
C      do  10  i = 1, m1
C      w(i,j)    = a(m1+i,j) - a(i,j)
C   10 w(i+m1,j) = b(i,j) + b(i,j+m1)
C*     call ywrm(0,'a21-a11',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b11+b12',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
C     .  c(1+m1,1+m1),ldc)
C*     call ywrm(0,'Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)
C
CC --- Q7 in c11 ---
C      do  20  j = 1, m1
C      do  20  i = 1, m1
C      w(i,j)    = a(i,m1+j) - a(m1+i,m1+j)
C   20 w(i+m1,j) = b(m1+i,j) + b(m1+i,j+m1)
C*     call ywrm(0,'a12-a22',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b21+b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
C     .  c,ldc)
C*     call ywrm(0,'Q7',1,6,fmt,c,1,ldc,m1,m1)
C
CC --- Q1 in c12 ---
C      do  30  j = 1, m1
C      do  30  i = 1, m1
C      w(i,j)    = a(i,j) + a(i+m1,j+m1)
C   30 w(i+m1,j) = b(i,j) + b(i+m1,j+m1)
C*     call ywrm(0,'a11+a22',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b11+b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
C     .  c(1,1+m1),ldc)
C*     call ywrm(0,'Q1',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)
C
CC ... c11 += Q1 and c22 += Q1
C      do  40  j = 1, m1
C      do  40  i = 1, m1
C      c(i,j)       = c(i,j)       + c(i,j+m1)
C   40 c(i+m1,j+m1) = c(i+m1,j+m1) + c(i,j+m1)
C*     call ywrm(0,'Q1+Q7',1,6,fmt,c(1,1),1,ldc,m1,m1)
C*     call ywrm(0,'Q1+Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)
C
CC --- Q2 in c21, Q4 in c12 ---
C      do  50  j = 1, m1
C      do  50  i = 1, m1
C      w(i,j)    = a(i+m1,j) + a(i+m1,j+m1)
C   50 w(i+m1,j) = b(i+m1,j) - b(i,j)
C*     call ywrm(0,'a21+a22',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b21-b11',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,b,ldb,0d0,
C     .  c(1+m1,1),ldc)
C      call dgemm(transa,transb,m1,m1,m1,1d0,a(1+m1,1+m1),lda,
C     .  w(1+m1,1),ldw,0d0,c(1,1+m1),ldc)
C*     call ywrm(0,'Q2',1,6,fmt,c(1+m1,1),1,ldc,m1,m1)
C*     call ywrm(0,'Q4',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)
C
CC ... c11 += Q4, c22 -= Q2, c21 = Q2+Q4 (watch optimizer here)
C      do  60  j = 1, m1
C      do  60  i = 1, m1
C      c(i,j)       = c(i,j)       + c(i,j+m1)
C      c(i+m1,j+m1) = c(i+m1,j+m1) - c(i+m1,j)
C      c(i+m1,j) = c(i+m1,j) + c(i,j+m1)
C   60 continue
C*     call ywrm(0,'Q1+Q7+Q4',1,6,fmt,c(1,1),1,ldc,m1,m1)
C*     call ywrm(0,'Q1+Q6-Q2',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)
C*     call ywrm(0,'c21=Q2+Q4',1,6,fmt,c(1+m1,1),1,ldc,m1,m1)
C
CC --- Q5 in c12, Q3 in w ---
C      do  70  j = 1, m1
C      do  70  i = 1, m1
C      w(i,j)    = a(i,j) + a(i,j+m1)
C   70 w(i+m1,j) = b(i,j+m1) - b(i+m1,j+m1)
C*     call ywrm(0,'a11+a12',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b12-b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,b(1+m1,1+m1),ldb,0d0,
C     .  c(1,1+m1),ldc)
C      call dgemm(transa,transb,m1,m1,m1,1d0,a,lda,
C     .  w(1+m1,1),ldw,0d0,w,ldw)
C*     call ywrm(0,'Q5',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)
C*     call ywrm(0,'Q3',1,6,fmt,w,1,ldw,m1,m1)
C
CC ... c11 -= Q5, c22 += Q3, c12 = Q3+Q5
C      do  80  j = 1, m1
C      do  80  i = 1, m1
C      c(i,j)       = c(i,j)       - c(i,j+m1)
C      c(i+m1,j+m1) = c(i+m1,j+m1) + w(i,j)
C   80 c(i,j+m1) = c(i,j+m1) + w(i,j)
C
C*     call ywrm(0,'c',1,6,fmt,c,1,ldc,m2,m2)
C      end