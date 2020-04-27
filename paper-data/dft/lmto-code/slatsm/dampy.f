C#define BLAS
C#define BLAS3
      subroutine dampy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
C- matrix multiplication, adding result into destination
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   c,ncc,nrc is the product matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   Adapted from dmpy.
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of (n,l) subblock of a into (l,m) subblock of b
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nca,nra,ncb,nrb,ncc,nrc,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
C Local parameters
      double precision ar
      integer i,j,k,nccj,ncbj
C#ifdefC APOLLO | HP
C      double precision vec_$ddot_i
C#endif
C#ifdef BLAS3
      integer lda,ldb
      character*1 transa,transb
C#endif

C#ifdefC TIMING
C      call dmpytm(0)
C#endif

C      if (min(nra,nca,nrb,ncb,nrc,ncc,n,m,l) < 0)
C     .  call rx('dmpy: bad integer argument')

C#ifdef BLAS3
      if (nra == 1) then
        lda = nca
        transa = 'n'
      elseif (nca == 1) then
        lda = nra
        transa = 't'
      else
        lda = -1
      endif
      if (nrb == 1) then
        ldb = ncb
        transb = 'n'
      elseif (ncb == 1) then
        ldb = nrb
        transb = 't'
      else
        ldb = -1
      endif
      if (min(lda,ldb) < 0 .or. nrc /= 1) goto 10
C#ifdefC PARALLEL
C      call pp_$dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,1d0,c,ncc)
C#else
      call dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,1d0,c,ncc)
C#endif
C#ifdefC TIMING
C      call dmpytm(1)
C#endif
      return
   10 continue
C#endif

C#ifdefC APOLLO | HP
C      do  200  i = n-1, 0, -1
C        do  200  j = m-1, 0, -1
C        c(i*nrc+j*ncc) = c(i*nrc+j*ncc)
C     .      + vec_$ddot_i(a(nra*i),nca,b(ncb*j),nrb,l)
C  200 continue
C#else

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = n-1, 0, -1
        ar = a(nra*i + nca*k)
        if (ar == 0) goto 20
C#ifdef BLAS
       call daxpy(m,ar,b(nrb*k),ncb,c(nrc*i),ncc)
C#elseC
C        nccj = -ncc + nrc*i
C        ncbj = -ncb + nrb*k
C        do  25  j = m-1, 0, -1
C          nccj = nccj + ncc
C          ncbj = ncbj + ncb
C          c(nccj) = c(nccj) + ar*b(ncbj)
C   25   continue
C#endif
   20 continue
C#endif
C#ifdefC TIMING
C      call dmpytm(1)
C#endif
      end
