      subroutine yqinv(cs,a,ofai,lda,nlev,n,w,ldw,ierr)
C- Inversion of a complex matrix using Strassen's algorithm
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   :if 'h', a is assumed hermitian
Ci         :1st character is an integer => call LAPACK
Ci         :Integer is kcmplx (see ztoyy)
Ci         :0: a has real, imaginary separated and
Ci         :   a = a(lda,lda,1..2) with Re a = a(lda,nc)
Ci         :1: a is in complex*16 format, i.e. a is dimensioned
Ci         :   a = a(1..2,lda,nc)
Ci         :2: Real, imaginary parts separated by columns; a dimensioned
Ci         :   a = a(lda,2,nc), with a(*,1..2,*) = real..imag
Ci         :Thus cs = '0h' calls LAPACK for hermitian matrix
Ci         :               with real,imag parts of a separated
Ci   a     :matrix to be inverted
Ci   ofai  :number of elements separating real, imaginary parts of a
Ci   lda   :leading dimension of a
Ci   n     :rank of the matrix to be inverted
Ci   nlev  :the maximum number of recursion levels allowed.
Ci          To avoid roundoff errors, nlev=2 is suggested.
Ci   w     :double precision work array of dimension ldw*(n+1)
Ci   ldw   :leading dimension of w
Co Outputs:
Co   a     :is overwritten by inverse of input a
Co   ierr  :returned nonzero if matrix was not fully inverted.
Cb Bugs
Cb   The algorithm fails if a22 is singular, even if a is not.
Cb   Similarly, if smaller subblocks are singular, yyqinv may fail
Cb   when called recursively.
Cr Remarks:
Cr   This is a front end for yyqinv, which which is passed the
Cr   real and imaginary parts of arrays.
Cu Updates
Cu   12 Nov 10 Add option for LAPACK call
C ----------------------------------------------------------------
      implicit none
      character*(*) cs
      integer n,lda,ldw,ierr,nlev,ofai
      double precision a(lda,n),w(ldw,1)
      integer kcmplx,ldb
      character*(10) css

C ... LAPACK call
      css = cs(1:1)
      if (css >= '0' .and. css <= '2') then
        read(css,'(I1)') kcmplx
        ldb = ofai/lda
        if (lda*ldb /= ofai) call rx('yqinv: bad ofai')
C       call yprm('before',2,a,ofai,lda,lda,lda)
C       Fast, but destroys subblock a(n+1..,*) and a(*,n+1..)
C       call ztoyy(a,lda,ldb,n,n,kcmplx,1)
        call ztoyy(a,lda,ldb,lda,lda,kcmplx,1)
C       call ztoyy(a,lda,ldb,lda,lda,1,kcmplx)
C       call yprm('after',2,a,ofai,lda,lda,lda)
        css = cs(2:)
        call zqinv(css,a,lda,-ldw*(n+1),n,w,ldw,ierr)
        call ztoyy(a,lda,ldb,lda,lda,1,kcmplx)

C ... Strassen call
      else
        call yyqinv(cs,a,a(1+ofai,1),lda,nlev,n,w,ldw,ierr)
      endif

      end

