      subroutine ymscpk(mode,kcplxi,kcplxf,nr,nc,s,lsa,lsb,d,lda,ldb)
C- Copies a matrix in one complex storage mode to one in another mode
C ----------------------------------------------------------------------
Ci   mode : 0 copies s into d
Ci          1 adds s into d
Ci          2 copies -s into d
Ci          3 adds -s into d
Ci          Add 4 if s is real
Ci   nr    :number of rows to copy
Ci   nc    :number of columns to copy
Ci   s     :source matrix
Ci   lsa   :formal leading dimension of s; see kcplxi
Ci   lsb   :formal second dimension of s;  see kcplxi
Ci   d     :destination matrix
Ci   lda   :formal leading dimension of d
Ci   ldb   :formal second dimension of d
Ci   kcplxi:source complex storage mode.  Matrix s is formally
Ci          dimensioned s(lsa,nc).  The true dimensions of s are:
Ci          0: s has real, imaginary separated and
Ci             Re s = s(lsa,nc) and  Im s offset from Re s by lsa*lsb
Ci          1: s is in complex*16 format, i.e. s is dimensioned
Ci             s = s(1..2,lsa,nc)
Ci          2: Real, imaginary parts separated by columns; s dimensioned
Ci             s = s(lsa,2,nc), with s(*,1..2,*) = real..imag
Ci   kcplxf:final complex storage mode, with conventions as in kcplxi
Cu Updates
Cu   23 Oct 13 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,kcplxi,kcplxf,nr,nc,lsa,lsb,lda,ldb
      double precision s(lsa,lsb),d(lda,ldb)
C ... Local parameters

C ... Copy to destination array, kcplxi format
      if (kcplxi == 0) then
        call ymscop(mode,nr,nc,lsa,lda,0,0,0,0,s,lsa*lsb,d,lda*ldb)
      elseif (kcplxi == 1) then
        call zmscop(mode,nr,nc,lsa,lda,0,0,0,0,s,d)
      elseif (kcplxi == 2) then
        call ymscop(mode,nr,nc,lsa*2,lda*2,0,0,0,0,s,lsa,d,lda)
      endif

      if (kcplxf == kcplxi .or. mode >= 4) return
      call ztoyy(d,lda,ldb,nr,nc,kcplxi,kcplxf)

      end

C     Test ymscpk ... to run, create a complex file 'sin'
C      subroutine fmain
C      implicit none
C      integer i,rdm,nr,nc,ifi,fopng,lsa,lsb,kcplx,kcplxf,lda,ldb
C      parameter (lsa=81,lsb=77,lda=55,ldb=99)
C      real(8) :: sin(lsa,lsb,2),sout(lda,ldb,2)
C      real(8), allocatable :: srd(:,:,:)
C
C      nr = 0; nc = 0
C      ifi = fopng('sin',-1,1)
C      i = rdm(ifi,0,0,' ',srd,nr,nc)
C      if (nr*nc > 55*99) stop 'oops'
C      allocate(srd(nr,nc,2)); srd = 0
C      rewind ifi
C      i = rdm(ifi,30,nr*nc*2,' ',srd,nr,nc)
C      sin(1:nr,1:nc,1:2) = srd(:,:,1:2)
C
C      call yprm('s before copy',2,sin,lsa*lsb,lsa,nr,nc)
C      kcplx = 2
C      call ztoyy(sin,lsa,lsb,nr,nc,0,kcplx)
C      print *, kcplx
C      call yprm('s in kcplx format',2+kcplx,sin,lsa*lsb,lsa,nr,nc)
C
C      kcplxf = 2
C      call ymscpk(0,kcplx,kcplxf,nr,nc,sin,lsa,lsb,sout,lda,ldb)
C      print *, kcplxf
C      call yprm('d in kcplx format',2+kcplxf,sout,lda*ldb,lda,nr,nc)
C
C      end
