      subroutine fmain
      implicit none
      logical cmdopt,a2bin
      character outs*80
      INTEGER   lda, ldb, ldc, ldw, n, m, l, i, ic, count, omppid
      double precision mflops
*      PARAMETER (lda = 1501, ldb = 1501, ldc = 1501, ldw=1501)
*      PARAMETER (lda = 5501, ldb = 5501, ldc = 5501, ldw=5501)
       PARAMETER (lda = 901, ldb = 901, ldc = 901, ldw=901)
C     PARAMETER (n = 600, m = 600, l = 600)

      complex(8):: sa( lda, lda ), sb( ldb, 1 ), sc( ldc, 1 )
      complex(8):: pa( lda, lda ), pb( ldb, 1 ), pc( ldc, 1 ),
     .             pw( ldw, ldw)
      double precision s_time, p_time, temp, s_mflops, p_mflops,
     .  walltime0
      double precision cpusec,swalltime,pwalltime,qwalltime,dwtime,delwc
      common /static/ sa, sb, sc, pa, pb, pc, pw

C --- Setup ---
      n = 300
      m = 1
      l = 300
      print *, 'n,l?'
      read(*,*) n,l
      mflops = 2.0d-6 * n * m * l * 4
      call finits(2,0,0,i)
      call nada
      WRITE (*,100) n, m, n, l, l, m
      walltime0 = delwc()
      walltime0 = dwtime()

C --- In-line matrix multiplication ---
      CALL init ( sc,sa,sb,ldc,lda,ldb,n,m,l )
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
    5 CALL scalar_matrix_multiply( sc, sa, sb, ldc, lda, ldb, n, m, l )
      count = count+1
      if (cpusec() - temp < 1d0 .and. count < 1000) goto 5
C     if (cpusec() - temp < 3d0 .and. count < 50) goto 5
      qwalltime = dwtime() - swalltime
      if (qwalltime /= 0 .and. qwalltime < 1) goto 5
      endif

      print 127, count
  127 format(' execute multiply',i6,' times')

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

C --- zgemv matrix multiplication ---
      CALL init(pc,pa,pb,ldc,lda,ldb,n,m,l)
      pwalltime = delwc()
      pwalltime = dwtime()
      temp = cpusec()
      print 129, dwtime() - walltime0
  129 format(' starting zgemv call ... total elapsed wall time',F9.4)
      do  ic = 1, count
      call zgemv('N',n,l,(1d0,0d0),pa,lda,pb,1,(0d0,0d0),pc,1)
      enddo
      p_time = (cpusec() - temp)/count
      pwalltime = delwc()
      print 130, pwalltime
  130 format(' finished zgemv call ... elapsed wall time',F9.4)
      if (p_time == 0d0) p_time = 1
      p_mflops = mflops / p_time
      CALL diff('compare in-line to zgemv', sc, pc, ldc, n, m, l)

      print 110, s_time, s_mflops, swalltime, swalltime/count,
     .           p_time, p_mflops, pwalltime, pwalltime/count

  100 format( 1X,'Matrix Multiply: ',I4,'x',I4,
     .        ' <- ',I4,'x',I4,' x ',I4,'x',I4)
  110 format(/1X,'Serial CPU time: ',F7.3,'  Serial MFlops: ',F7.0,
     .                 '   Wall time: ',F7.3,' (',F12.6,' /product)'
     .       /1X,' dgemv CPU time: ',F7.3,'   dgemv MFlops: ',F7.0,
     .                 '   Wall time: ',F7.3,' (',F12.6,' /product)':
     .       /1X,'  qmpy CPU time: ',F7.3,'    qmpy MFlops: ',F7.0
     .                 '   Wall time: ',F7.3)

C$OMP parallel private(i)
      i = omppid(0)
      if (i > 1) print *, 'OMP_NUM_THREADS =',i
C$OMP end parallel

      end

      SUBROUTINE init ( c, a, b, ldc, lda, ldb, n, m, l )
C- Initialize arrays
      integer lda, ldb, ldc, n, m, l
      double complex a( lda, * ), b( ldb, * ), c( ldc, * )
*
      do i = 1, n
         do j = 1, m
            c(i,j) = dble( i+j )
            end do
         end do
*
      do i = 1, n
         do j = 1, l
            a(i,j) = dcmplx(dble(i-j),dble(i+j))
            end do
         end do
*
      do i = 1, l
         do j = 1, m
            b(i,j) = dcmplx(dble(i*i),dble(i*i))
            end do
         end do
      end

      SUBROUTINE diff (strn, sc, pc, ldc, n, m, l )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER ldc, n, m, l, i, j
      double complex sc( ldc, * )
      double complex pc( ldc, * )
      double precision tol,maxerr

      tol = 1e-6
      if (dble(n)*dble(m)*dble(l) > 1e8) tol = 1d-5
      maxerr = 0

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
  100 FORMAT(1X,a,'*** ERROR ***   Arrays Have Different Results,',
     .   ' max err=',1pe8.1)
  110 FORMAT(1X,a,'... arrays agree to within tol=',1pe8.1)
      END

      SUBROUTINE scalar_matrix_multiply(c, a, b, ldc, lda, ldb, n, m, l)
C- In-line matrix multipliation
      integer lda, ldb, ldc, n, m, l
      double complex a (lda,*), b(ldb,*), c(ldc,*)
      do j = 1, m
        do i = 1, n
          c(i,j) = 0
        end do
      end do
      do j = 1, m
C$OMP   parallel do  firstprivate(j) private(k,i) shared(a,b,c)
        do k = 1, l
          do i = 1, n
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
C$OMP   end parallel do
      end do
      end
