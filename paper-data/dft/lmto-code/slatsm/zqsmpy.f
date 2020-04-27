      recursive subroutine zqsmpy(isw,transa,transb,m,k,a,lda,b,ldb,
     .  beta,c,ldc)
C- Multiplication of two matrices whose result is hermitian
C ----------------------------------------------------------------
Ci Inputs:
Ci   isw   :1s digit:
Ci         :0 straight dgemm call: full c generated (no partitioning)
Ci         :1 Recursive call for upper triangle in diagonal subblocks,
Ci            while n>nmin and nmin>0
Ci         :  else:: straight dgemm for (full) diagonal subblocks
Ci         :  Default value of nmin may be changed by calling zqsmpy0
Ci         :2 Calculate 11, 12, 22 subblocks of c; no recursion for c11, c22
Ci         :  This is a good choice for threaded MPI zgemm.
Ci   isw   :10s digit:
Ci         :0 only upper triangle is guaranteed to be returned
Ci         :1 entire matrix is filled in
Ci   transa: specifies the form of op( A ) to be used in multiplication
Ci           'N' or 'n',  op( A ) = A.
Ci           'T' or 't',  op( A ) = A', with A' = transpose A
Ci           'C' or 'c',  op( A ) = conjg( A' ).
Ci         : NB: current implementation requires transa='N' if
Ci         : in-line or Strassen algorithm is used
Ci   transb: specifies the form of op( B ) to be used in multiplication
Ci           'N' or 'n',  op( B ) = B.
Ci           'T' or 't',  op( B ) = B', with B' = transpose B
Ci           'C' or 'c',  op( B ) = conjg( B' ).
Ci         : NB: current implementation requires transb='N' if
Ci         : in-line or Strassen algorithm is used
Ci   m     :2nd dimension of b; dimension of result matrix
Ci         :It is an error for m>lda.
Ci   k     :2nd dimension of a; length of inner product
Ci   a     :left matrix
Ci   lda   :leading dimension of a
Ci   b     :right matrix
Ci   ldb   :leading dimension of b
Ci   beta  :add beta*c to final result
Ci   ldc   :leading dimension of c
Co Outputs:
Co   c     :c is overwritten by = a*b + beta*c (an mxm matrix)
Cl Local variables
Cl   nmin  :If >0 smallest matrix which on diagonal block calculated recursively
Cr Remarks:
Cr   This is a complex analog of dqsmpy
Cu Updates
Cu   09 Jul 10 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*1 transa, transb
      integer lda,ldb,ldc,k,m
      double complex a(lda,k),b(ldb,m),c(ldc,m),beta
C ... Local parameters
      logical nota,notb,lsame
      integer nmin,i,j,isw,m1,m1p,nmin0,iswl
      save nmin
C#ifdefC DEBUG
C      character*10 fmt
C      data fmt /'(8f16.9)'/
C#endif
      data nmin /256/

C     print *, 'enter zqsmpy',isw,m

      if (m == 0 .or. k == 0) return
      if (m > lda) call rx('zqsmpy: m>lda')

C --- Straight dgemm if m<nmin ---
      if (k < nmin .or. mod(isw,10) == 0) then
        call zgemm(transa,transb,m,m,k,(1d0,0d0),a,lda,b,ldb,
     .    beta,c,ldc)
C#ifdefC DEBUG
C        call ywrm(0,'c',3,6,fmt,c,1,ldc,m,m)
C#endif
C       print *, 'exit zqsmpy',isw,m
        return
      endif

      m1 = (m+1)/2
      m1p = m/2
      nota  = lsame( transa, 'N' )
      notb  = lsame( transb, 'N' )
      iswl = mod(isw,10); if (m1 < nmin .or. nmin == 0 .or. iswl >= 2) iswl = 0

C --- c12 = a11 b12 + a12 b22 ---
      if (notb) then
        call zgemm(transa,transb,m1,m1p,k,(1d0,0d0),a,lda,
     .    b(1,1+m1),ldb,beta,c(1,1+m1),ldc)
      else
        call zgemm(transa,transb,m1,m1p,k,(1d0,0d0),a,lda,
     .    b(1+m1,1),ldb,beta,c(1,1+m1),ldc)
      endif
C     call ywrm(0,'c12',3,6,fmt,c(1,1+m1),1,ldc,m1,m1p)

C --- Upper triangle of c11 = a11 b11 + a12 b21 with recursive call ---
      call zqsmpy(iswl,transa,transb,m1,k,a,lda,b,ldb,beta,c,ldc)
*     call ywrm(0,'c11',3,6,fmt,c,1,ldc,m1,m1)

C --- Upper triangle of c22 = a21 b12 + a22 b22 ---
      if (m1p == 0) then
        return
C ... Upper triangle with recursive call
      elseif (m1p >= nmin) then
        if (nota) then
        if (notb) then
          call zqsmpy(iswl,transa,transb,m1p,k,a(1+m1,1),lda,
     .      b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call zqsmpy(iswl,transa,transb,m1p,k,a(1+m1,1),lda,
     .      b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        else
        if (notb) then
          call zqsmpy(iswl,transa,transb,m1p,k,a(1,1+m1),lda,
     .      b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call zqsmpy(iswl,transa,transb,m1p,k,a(1,1+m1),lda,
     .      b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        endif
C ... Entire square of c22, with dgemm call
      else
        if (nota) then
        if (notb) then
          call zgemm(transa,transb,m1p,m1p,k,(1d0,0d0),
     .      a(1+m1,1),lda,b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call zgemm(transa,transb,m1p,m1p,k,(1d0,0d0),
     .      a(1+m1,1),lda,b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        else
        if (notb) then
          call zgemm(transa,transb,m1p,m1p,k,(1d0,0d0),
     .      a(1,1+m1),lda,b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call zgemm(transa,transb,m1p,m1p,k,(1d0,0d0),
     .      a(1,1+m1),lda,b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        endif
      endif
*     call ywrm(0,'c22',3,6,fmt,c(1+m1,1+m1),1,ldc,m1p,m1p)

      if (isw >= 10) then
        do  i = 1, m
        do  j = i+1, m
          c(j,i) = dconjg(c(i,j))
        enddo
        enddo
      endif

C#ifdefC DEBUG
C      print '(1x,a,i4)', 'dimension of c =',m
C      call ywrm(0,'c',3,6,fmt,c,1,ldc,m,m)
C#endif

      return

      entry zqsmpy0(nmin0,isw)
      if (isw > 0) then
        nmin = nmin0
      else
        nmin0 = nmin
      endif

      end
