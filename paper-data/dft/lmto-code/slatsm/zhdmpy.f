      recursive subroutine zhdmpy(n,k,a,lda,diag,beta,c,ldc)
C- Multiplies  a(+) d a , where d is a real and diagonal matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of rows in a
Ci   k     :number of columns in a
Ci   a     :complex matrix to be multiplied
Ci   lda   :leading dimension of a
Ci   diag  :real diagonal matrix of same dimension as a
Ci   beta  :add beta*c to final result
Ci   ldc   :leading dimension of c
Co Outputs
Co   c     :c is overwritten by beta*c + (a+ d a)
Cl Local variables
Cl   b     :internal array holding (da)
Cr Remarks
Cr   As written, and internal array of dimension (n,n) is allocated
Cr   Some memory could be saved by alocating smaller internal
Cr   arrays b and repeating steps in zqsmpy.f.
Cu Updates
Cu   12 Apr 13 Added argument k
Cu   18 Apr 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,k,lda,ldc
      double complex a(lda,n),c(ldc,n),beta
      double precision diag(n)
C ... Local parameters
      integer i,j
      complex(8), allocatable :: b(:,:)

      allocate(b(k,n))

      do  i = 1, k
        do  j = 1, n
          b(i,j) = diag(i)*a(i,j)
        enddo
      enddo
C     call zprm('b',2,b,n,n,k)
      call zqsmpy(0,'C','N',n,k,a,lda,b,k,beta,c,ldc)
C     call zprm('c',2,c,ldc,n,n)

      deallocate(b)
      end

C      subroutine fmain
CC-    Test zhdmpy
C      implicit none
C      INTEGER   lda, ldb, ldc, n, ic, count
C      PARAMETER (lda = 1501, ldb = lda, ldc = lda, n=lda-1)
C      integer i,j,k
C      double precision mflops,s_time,s_mflops,p_time,p_mflops,cpusec
C      double complex a(lda,lda),c(ldc,ldc),cref(ldc,ldc),b(lda,lda),beta
C      double precision diag(lda*2)
C      common /static/ a,diag
C
C      call dpzero(a,lda*lda*2)
C      call dpzero(c,ldc*ldc*2)
C      beta = 0d0
C
C      k = lda/2 + 1
C
CC --- Setup standard test ---
C      WRITE (*,100) lda,k
C  100 format( 1X,
C     .'tzhdmpy matrix multiply  s+ d s  w/ d real and diagonal: lda =',
C     .  I4,'  k = ',I4)
C
C      call ran1in(1)
C      call init(a, lda, lda, 'n')
C      call z2herm('U',lda,lda,a)
CC     call zprm('a',2,a,lda,lda,lda)
C      call init(c, ldc, ldc, 'n')
C      call z2herm('U',ldc,ldc,c)
CC     call zprm('c',2,c,ldc,ldc,ldc)
C      call init(diag, lda, 1, 'n')
CC     call prmx('diag',diag,lda,lda,1)
C
C      do  i = 1, lda
C        do  j = 1, lda
C          b(i,j) = diag(i)*a(i,j)
C        enddo
C      enddo
CC     call zprm('b',2,b,ldb,ldb,ldb)
C
C      s_time = cpusec()
C      count = 0
C    5 continue
C      call zgemm('C','N',n,n,k,(1d0,0d0),a,lda,b,ldb,(0d0,0d0),cref,ldc)
C      count = count+1
C      if (cpusec() - s_time < 1d0 .and. count < 60) goto 5
C      s_time = (cpusec() - s_time)/count
C      print 127, count
C  127 format(' execute multiply',i4,' times')
C
CC     Make reference
C      call zcopy(ldc*ldc,c,1,cref,1)
C      call zgemm('C','N',n,n,k,(1d0,0d0),a,lda,b,ldb,beta,cref,ldc)
C      call zhdmpy(n,k,a,lda,diag,beta,c,ldc)
C      CALL diff('compare to reference', c, cref, ldc,ldc,ldc,ldc)
C
CC     Run last time for counting
C      p_time = cpusec()
C      do  ic = 1, count
C        call zhdmpy(n,k,a,lda,diag,(0d0,0d0),c,ldc)
C      enddo
C      p_time = (cpusec() - p_time)/count
C
C      mflops = 2.0d-6 * n * n * k * 4
C      s_mflops = mflops / s_time
C      p_mflops = mflops / p_time
C      print 110, s_time, s_mflops, p_time, p_mflops
C  110 format(/1X,' zgemm  Time: ',F7.3,'   zgemm  MFlops: ',F7.0,
C     .       /1X,' zqsmpy Time: ',F7.3,'   zqsmpy MFlops: ',F7.0)
C
C      print 111, s_time*count, s_time, p_time*count, p_time
C  111 format(/' time elapsed, zgemm  multiply',F7.2,
C     .    ' sec (',F9.4,' sec/product)'/
C     .        ' time elapsed, zqsmpy multiply',F7.2,
C     .    ' sec (',F9.4,' sec/product)')
C
C      end
C
C      SUBROUTINE init ( a, lda, lda2, cs)
CC- Initialize arrays
C      implicit none
C      integer lda, lda2
C      double precision a( 2, lda, lda2)
C      integer i,j
C      character *1 cs
C      real ran1
CC     call ran1in(1)
C
C      do 10 i = 1, lda
C      do 10 j = 1, lda2
C      a(1,i,j) = ran1()
C   10 a(2,i,j) = ran1()
C
C      if (cs == 'h') then
C        do 20 i = 1, min(lda,lda2)
C        do 20 j = 1, i
C        a(1,i,j) =  a(1,j,i)
C   20   a(2,i,j) = -a(2,j,i)
C        do  22 i = 1, min(lda,lda2)
C   22   a(2,i,i) = 0
C      endif
C
C      end
C      SUBROUTINE diff (strn, sc, pc, ldc, n, m, l )
CC- Compare the two arrays for differences
C      implicit none
C      character*(*) strn
C      INTEGER ldc, n, m, l, i, j
C      double complex sc( ldc, * )
C      double complex pc( ldc, * )
C      double precision tol
C
C      tol = 1e-12
C      if (dble(n)*dble(m)*dble(l) > 1e8) tol = 1d-10
C
C      DO i = 1, n
C         DO j = 1, m
C            IF ( ABS( sc(i,j) - pc(i,j) ) .GT. 1.0d-6 ) THEN
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
C      END
