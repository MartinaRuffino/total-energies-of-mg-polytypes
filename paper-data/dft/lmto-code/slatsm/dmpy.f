C#define BLAS
C#define BLAS3
      subroutine dmpy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
C- matrix multiplication
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
      integer i,j,k,nccj,nrci,nrcicj,ncbj
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
      if (min(lda,ldb) < 0 .or. nrc /= 1) goto 11
C#ifdefC PARALLEL
C      call pp_$dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
C#else
      call dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
C#endif
C#ifdefC TIMING
C      call dmpytm(1)
C#endif
      return
   11 continue
C#endif

C#ifdefC CRAY
C      call mxma(a,nra,nca,b,nrb,ncb,c,nrc,ncc,n,l,m)
C#elseifC APOLLO | HP
C      do  200  i = n-1, 0, -1
C        do  200  j = m-1, 0, -1
C        c(i*nrc+j*ncc) = vec_$ddot_i(a(nra*i),nca,b(ncb*j),nrb,l)
CC The equivalent in-line code ...
Cc        double precision sum
Cc        sum = 0
Cc        nakpi = nra*i
Cc        nbjpk = ncb*j
Cc        do  210  k = l-1, 0, -1
Cc          sum = sum + a(nakpi)*b(nbjpk)
Cc          nakpi = nakpi + nca
Cc          nbjpk = nbjpk + nrb
Cc  210   continue
C  200 continue
C#else

C --- Initialize array to zero ---
      do  10  i = n-1, 0, -1
        nrci = nrc*i
        nccj = -ncc
        do  10  j = m-1, 0, -1
        nccj = nccj + ncc
        nrcicj = nrci + nccj
        c(nrcicj) = 0
   10 continue

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = n-1, 0, -1
        ar = a(nra*i + nca*k)
        if (ar == 0) goto 20
C#ifdef BLAS
        call daxpy(m,ar,b(nrb*k),ncb,c(nrc*i),ncc)
C#elseC
C        nrci = nrc*i
C        nccj = -ncc
C        ncbj = -ncb + nrb*k
C        do  15  j = m-1, 0, -1
C          nccj = nccj + ncc
C          ncbj = ncbj + ncb
C          nrcicj = nrci + nccj
C          c(nrcicj) = c(nrcicj) + ar*b(ncbj)
C   15   continue
C#endif
   20 continue
C#endif
C#ifdefC TIMING
C      call dmpytm(1)
C#endif
      end
C#ifdefC TIMING
C      subroutine dmpytm(i)
C      implicit none
C      integer i
C      integer i1mach,iprint
C      double precision xx,nettim
C      save nettim
C      data nettim /0d0/
C
C      if (iprint() < 52) return
C#ifdefC DEBUG
C      if (i == 0) then
C        call cpudel(6,'Enter dmpytm',xx)
C      else
C        call cpudel(-1,'Enter dmpytm',xx)
C      endif
C#elseC
C      call cpudel(-1,' ',xx)
C#endifC
C      if (i == 0) return
C      nettim = nettim+xx
C      print 333, xx, nettim
C  333 format(22x,'Done dmpy:  time:',g10.3,' dmpy t:',g10.3)
C      end
C#endif
