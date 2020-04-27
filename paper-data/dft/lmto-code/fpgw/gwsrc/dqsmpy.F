      recursive subroutine dqsmpy(isw,transa,transb,m,k,a,lda,b,ldb,
     .  beta,c,ldc)
C- Multiplication of two matrices whose result is symmetric
C ----------------------------------------------------------------
Ci Inputs:
Ci   isw   :1s digit:
Ci         :0 straight dgemm call: full c generated (no partitioning)
Ci         :1 Recursive call for upper triangle in diagonal subblocks,
Ci            while n>nmin and nmin>0
Ci         :  else:: straight dgemm for (full) diagonal subblocks
Ci         :  Default value of nmin may be changed by calling dqsmpy0
Ci         :2 Calculate 11, 12, 22 subblocks of c; no recursion for c11, c22
Ci         :  This is a good choice for threaded MPI dgemm.
Ci         :2 (old) Strassen algorithm (commented out --- too slow)
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
Cr   For Strassen algorithm, see Numerical Recipes, 2.11.
Cr   (For debugging:)
Cr   Start w/ matrices a(5,7) and symmetric v(7,7)  Then do:
Cr   set su = (" a -t -a b a v -x -a a a " \
Cr            "-split a 1,4,7 1,5,9 -pop b -split b 1,5,9 1,4,7 -pop")
Cr   set q1 = "a11 a22 -+ b11 b22 -+ -x -a q1"
Cr   set q2 = "a21 a22 -+  b11 -x -a q2"
Cr   set q3 = "a11 b12 b22 -- -x -a q3"
Cr   set q4 = "a22 b21 b11 -- -x -a q4"
Cr   set q5 = "a11 a12 -+ b22 -x -a q5"
Cr   set q6 = "a21 a11 -- b11 b12 -+ -x -a q6"
Cr   set q7 = "a12 a22 -- b21 b22 -+ -x -a q7"
Cr   set c11 = "q1 q4 -+ q5 -- q7 -+"
Cr   set c21 = "q2 q4 -+"
Cr   set c12 = "q3 q5 -+"
Cr   set c22 = "q1 q3 -+ q2 -- q6 -+"
Cr   set red = "-sub 1,5,1,5"
Cr   mc -f9f7.0 $su $q1 $q2 $q3 $q4 $q5 $q6 $q7 \
Cr              $c11 $c12 -ccat $c21 $c22 -ccat -rcat $red -w . a b -x
Cu Updates
Cu   09 Jul 10 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*1 transa, transb
      integer lda,ldb,ldc,k,m
      double precision a(lda,k),b(ldb,m),c(ldc,m),beta
C ... Local parameters
      logical nota,notb,lsame
      integer nmin,i,j,isw,m1,m1p,nmin0,iswl
C     integer k1,k1p
C     real(8),allocatable:: w1(:,:),w2(:,:)
      save nmin
C#ifdefC DEBUG
C      character*10 fmt
C      data fmt /'(8f16.9)'/
C#endif
      data nmin /256/

C     print *, 'enter dqsmpy',isw,m

      if (m == 0 .or. k == 0) return
      if (m > lda) call rx('dqsmpy: m>lda')

C --- Straight dgemm if m<nmin ---
      if (k < nmin .or. mod(isw,10) == 0) then
        call dgemm(transa,transb,m,m,k,1d0,a,lda,b,ldb,beta,c,ldc)
C#ifdefC DEBUG
C      call ywrm(0,'c',1,6,fmt,c,1,ldc,m,m)
C#endif
C       print *, 'exit dqsmpy',isw,m
        return
      endif

      m1 = (m+1)/2
      m1p = m/2
      nota  = lsame( transa, 'N' )
      notb  = lsame( transb, 'N' )
      iswl = mod(isw,10); if (m1 < nmin .or. nmin == 0 .or. iswl >= 2) iswl = 0

C --- c12 = a11 b12 + a12 b22 ---
      if (notb) then
        call dgemm(transa,transb,m1,m1p,k,1d0,a,lda,
     .    b(1,1+m1),ldb,beta,c(1,1+m1),ldc)
      else
        call dgemm(transa,transb,m1,m1p,k,1d0,a,lda,
     .    b(1+m1,1),ldb,beta,c(1,1+m1),ldc)
      endif
C     call ywrm(0,'c12',1,6,fmt,c(1,1+m1),1,ldc,m1,m1p)

C --- Upper triangle of c11 = a11 b11 + a12 b21 with recursive call ---
      call dqsmpy(iswl,transa,transb,m1,k,a,lda,b,ldb,beta,c,ldc)
C     call ywrm(0,'c11',1,6,fmt,c,1,ldc,m1,m1)

C --- Upper triangle of c22 = a21 b12 + a22 b22 ---
      if (m1p == 0) then
        return
C ... Upper triangle with recursive call
      elseif (m1p >= nmin) then
        if (nota) then
        if (notb) then
          call dqsmpy(iswl,transa,transb,m1p,k,a(1+m1,1),lda,
     .      b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call dqsmpy(iswl,transa,transb,m1p,k,a(1+m1,1),lda,
     .      b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        else
        if (notb) then
          call dqsmpy(iswl,transa,transb,m1p,k,a(1,1+m1),lda,
     .      b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call dqsmpy(iswl,transa,transb,m1p,k,a(1,1+m1),lda,
     .      b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        endif
C ... Entire square of c22, with dgemm call
      else
        if (nota) then
        if (notb) then
          call dgemm(transa,transb,m1p,m1p,k,1d0,
     .      a(1+m1,1),lda,b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call dgemm(transa,transb,m1p,m1p,k,1d0,
     .      a(1+m1,1),lda,b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        else
        if (notb) then
          call dgemm(transa,transb,m1p,m1p,k,1d0,
     .      a(1,1+m1),lda,b(1,1+m1),ldb,beta,c(1+m1,1+m1),ldc)
        else
          call dgemm(transa,transb,m1p,m1p,k,1d0,
     .      a(1,1+m1),lda,b(1+m1,1),ldb,beta,c(1+m1,1+m1),ldc)
        endif
        endif
      endif
C     call ywrm(0,'c22',1,6,fmt,c(1+m1,1+m1),1,ldc,m1p,m1p)

      if (isw >= 10) then
        do  i = 1, m
        do  j = i+1, m
          c(j,i) = c(i,j)
        enddo
        enddo
      endif

C#ifdefC DEBUG
C      print '(1x,a,i4)', 'dimension of c =',m
C      call ywrm(0,'c',1,6,fmt,c,1,ldc,m,m)
C#endif

      return

C ----------------------- Strassen Algorithm -------------------------
C   20 continue
C      k1 = (k+1)/2
C      k1p = k/2
C      allocate(w1(m1,k1),w2(k1,m1))
C
CC --- c11 += Q1 = (a11+a22)*(b11+b22) ---
C      do  j = 1, k1p
C        do  i = 1, m1p
C          w1(i,j) = a(i,j) + a(i+m1,j+k1)
C          w2(j,i) = b(j,i) + b(j+k1,i+m1)
C        enddo
C        if (m1p < m1) then
C          w1(m1,j) = a(m1,j)
C          w2(j,m1) = b(j,m1)
C        endif
C      enddo
C      if (k1p < k1) then
C        do  i = 1, m1
C          w1(i,k1)= a(i,k1)
C          w2(k1,i) = b(k1,i)
C        enddo
C      endif
C*     call ywrm(0,'a11+a22',1,6,fmt,w1,1,m1,m1,k1)
C*     call ywrm(0,'b11+b22',1,6,fmt,w2,1,k1,k1,m1)
C      call dgemm(transa,transb,m1,m1,k1,1d0,w1,m1,w2,k1,0d0,c,ldc)
C*     call ywrm(0,'Q1',1,6,fmt,c,1,ldc,m1,m1)
C
CC --- c11 -= Q2 = Q1-(a21+a22)*b11.  Result dimensioned (m1p,k1) ---
C      do  j = 1, k1p
C        do  i = 1, m1p
C          w1(i,j) = a(i+m1,j) + a(i+m1,j+k1)
C        enddo
C      enddo
C      if (k1p < k1) then
C        do  i = 1, m1
C          w1(i,k1)= a(i+m1,k1)
C        enddo
C      endif
C*     call ywrm(0,'a21+a22',1,6,fmt,w1,1,m1,m1p,k1)
C      call dgemm(transa,transb,m1p,m1,k1,-1d0,w1,m1,b,ldb,1d0,c,ldc)
C*     call ywrm(0,'Q1-Q2',1,6,fmt,c,1,ldc,m1,m1)
C
CC ... Copy Q1-Q2 to c22
C      call dmcpy(c,ldc,1,c(1+m1,1+m1),ldc,1,m1p,m1p)
C
CC --- c22 += Q6 = (a21-a11)*(b11+b12) ---
C      do  j = 1, k1
C        do  i = 1, m1p
C          w1(i,j) = a(i+m1,j) - a(i,j)
C          w2(j,i) = b(j,i) + b(j,i+m1)
C        enddo
C      enddo
C*     call ywrm(0,'a21-a11',1,6,fmt,w1,1,m1,m1p,k1)
C*     call ywrm(0,'b11+b12',1,6,fmt,w2,1,k1,k1,m1p)
C      call dgemm(transa,transb,m1p,m1p,k1,1d0,w1,m1,w2,k1,1d0,
C     .  c(1+m1,1+m1),ldc)
C*     call ywrm(0,'Q1-Q2+Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1p,m1p)
C
CC --- c11 += Q7 = (a12-a22)*(b21+b22) => c11=Q1-Q2+Q7 ---
C      do  j = 1, k1p
C        do  i = 1, m1p
C          w1(i,j) = a(i,j+k1) - a(i+m1,j+k1)
C          w2(j,i) = b(j+k1,i) + b(j+k1,i+m1)
C        enddo
C        if (m1p < m1) then
C          w1(m1,j) = a(i,k1+j)
C          w2(j,m1) = b(j+k1,m1)
C        endif
C      enddo
C*     call ywrm(0,'a12-a22',1,6,fmt,w1,1,m1,m1,k1p)
C*     call ywrm(0,'b21+b22',1,6,fmt,w2,1,k1,k1p,m1)
C      call dgemm(transa,transb,m1,m1,k1p,1d0,w1,m1,w2,k1,1d0,c,ldc)
C*     call ywrm(0,'Q1-Q2+Q7',1,6,fmt,c,1,ldc,m1,m1)
C
CC --- c12 <- Q5=(a11+a12)*b22, (c21)+ <- Q3=a11*(b12-b22) ---
CC     c22 += Q3 => c22 = Q1 - Q2 + Q3 + Q6
CC     c11 += (Q3)+
C      do  j = 1, k1p
C        do  i = 1, m1p
C          w1(i,j) = a(i,j) + a(i,j+k1)
C          w2(j,i) = b(j,i+m1) - b(j+k1,i+m1)
C        enddo
C        if (m1p < m1) then
C          w1(m1,j) = a(m1,j) + a(m1,j+k1)
C          w2(j,m1) = b(j+k1,m1)
C        endif
C      enddo
C      if (k1p < k1) then
C        do  i = 1, m1p
C          w2(k1,i) = b(k1,i+m1)
C        enddo
C      endif
C*     call ywrm(0,'a11+a12',1,6,fmt,w1,1,m1,m1,k1p)
C*     call ywrm(0,'b12-b22',1,6,fmt,w2,1,k1,k1,m1p)
C
C      call dgemm(transa,transb,m1,m1p,k1,1d0,a,lda,w2,k1,0d0,
C     .  c(1,1+m1),ldc)
C*     call ywrm(0,'Q3',1,6,fmt,c(1,1+m1),1,ldc,m1,m1p)
CC     Copy transpose of Q3 to c(2,1)
C      call dmcpy(c(1,1+m1),ldc,1,c(1+m1,1),1,ldc,m1,m1p)
C*     call ywrm(0,'Q3+',1,6,fmt,c(1+m1,1),1,ldc,m1p,m1)
CC     c22 += Q3
C      call dmadd(c(1,1+m1),ldc,1,1d0,c(1+m1,1+m1),ldc,1,1d0,
C     .  c(1+m1,1+m1),ldc,1,m1p,m1p)
C*     call ywrm(0,'c22=Q1-Q2+Q3+Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1p,m1p)
CC     c12 <- Q5
C      call dgemm(transa,transb,m1,m1p,k1p,1d0,w1,m1,b(1+k1,1+m1),ldb,
C     .  0d0,c(1,1+m1),ldc)
C*     call ywrm(0,'Q5',1,6,fmt,c(1,1+m1),1,ldc,m1,m1p)
C
CC ... c11 -= Q5
C      call dmadd(c(1,1+m1),ldc,1,-1d0,c,ldc,1,1d0,c,ldc,1,m1,m1p)
C*     call ywrm(0,'Q1-Q2-Q5+Q7',1,6,fmt,c,1,ldc,m1,m1)
C
CC ... c21 += (Q5)+ => c21 += (Q3+Q5)+
C      do  j = 1, m1
C        do  i = 1, m1p
C          c(i+m1,j) = c(i+m1,j) + c(j,i+m1)
C        enddo
C      enddo
C*     call ywrm(0,'c21=(Q3+Q5)+',1,6,fmt,c(1+m1,1),1,ldc,m1p,m1)
C
CC ... c11 += (Q3+Q5)+ => c11 = Q1 - Q2 + (Q3+Q5)+ - Q5 + Q7
C      call dmadd(c(1+m1,1),ldc,1,1d0,c,ldc,1,1d0,c,ldc,1,m1p,m1)
C*     call ywrm(0,'c11=Q1-Q2+(Q3+Q5)+ -Q5+Q7',1,6,fmt,c,1,ldc,m1,m1)
C
CC ... c12 = (c21)+ => c12 = Q3+Q5
C      call dmcpy(c(1+m1,1),1,ldc,c(1,1+m1),ldc,1,m1,m1p)
C*     call ywrm(0,'c21=Q3+Q5',1,6,fmt,c(1+m1,1),1,ldc,m1p,m1)
C
C      deallocate(w1,w2)
C
CC#ifdefC DEBUG
CC      call ywrm(0,'c',1,6,fmt,c,1,ldc,m,m)
CC#endif

      entry dqsmpy0(nmin0,isw)
      if (isw > 0) then
        nmin = nmin0
      else
        nmin0 = nmin
      endif

      end
