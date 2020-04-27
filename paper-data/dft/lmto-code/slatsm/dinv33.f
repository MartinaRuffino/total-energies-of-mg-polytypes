      subroutine dinv33(matrix,iopt,invrse,det)
C- Inverts 3x3 matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   matrix:  matrix to be inverted
Ci   iopt:  if 0, usual inverse
Ci             1, transpose of inverse
Ci             2, 2*pi*inverse
Ci             3, 2*pi*transpose of inverse
Co Outputs
Co   invrse    see iopt
Co   det:      determinant, or det/2*pi (sign ok ??)
Cr Remarks
Cr   To generate reciprocal lattice vectors, call dinv33(plat,3,rlat)
C ----------------------------------------------------------------
      implicit none
      integer iopt,i,j
      double precision matrix(3,3),invrse(3,3)

      double precision det,xx,invm(3,3)
      procedure(real(8)) :: ddot

      call cross(matrix(1,2),matrix(1,3),invm     )
      call cross(matrix(1,3),matrix     ,invm(1,2))
      call cross(matrix     ,matrix(1,2),invm(1,3))
      det = ddot(3,matrix,1,invm,1)
      if (det == 0d0) then
        write(*,350) ((matrix(i,j),i=1,3),j=1,3)
  350   format(3f10.5)
        call rx('INV33: 3x3 matrix has vanishing determinant')
      endif
      call dcopy(9,invm,1,invrse,1)

      if (iopt >= 2) det = det/(8*datan(1d0))
      if (mod(iopt,2) == 0) then
        do  10  i = 1, 3
          do  10  j = i+1, 3
          xx = invrse(i,j)
          invrse(i,j) = invrse(j,i)
          invrse(j,i) = xx
   10   continue
      endif
      call dscal(9,1/det,invrse,1)
      end

       subroutine dmpy31(mode,a,b,c)
C- Multiplies 3x3 matrix a into 3x1 matrix b.  Result stored in c
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0  c = a*b
Ci         :1  c = c + a*b
Ci         :10 c = trans(a)*b
Ci         :11 c = c + trans(a)*b
Ci   a     :left matrix
Ci   b     :right matrix
Co Outputs
Co   c     :Product a*b
Cr Remarks
Cr   It is permissible for any of a,b,c to use the same address space
Cu Updates
Cu   30 Jan 13 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode
      double precision a(3,3), b(3), c(3)
C ... Local parameters
      double precision cloc(3)

      if (mode/10 == 0) then
        cloc(1) = a(1,1)*b(1) + a(1,2)*b(2) + a(1,3)*b(3)
        cloc(2) = a(2,1)*b(1) + a(2,2)*b(2) + a(2,3)*b(3)
        cloc(3) = a(3,1)*b(1) + a(3,2)*b(2) + a(3,3)*b(3)
      else
        cloc(1) = a(1,1)*b(1) + a(2,1)*b(2) + a(3,1)*b(3)
        cloc(2) = a(1,2)*b(1) + a(2,2)*b(2) + a(3,2)*b(3)
        cloc(3) = a(1,3)*b(1) + a(2,3)*b(2) + a(3,3)*b(3)
      endif
      if (mod(mode,10) == 0) then
        c(1) = 0; c(2) = 0; c(3) = 0
      endif
      c(1) = c(1) + cloc(1)
      c(2) = c(2) + cloc(2)
      c(3) = c(3) + cloc(3)

      end

C     Test
C      subroutine fmain
C
C      double precision a(3,3), b(3), c(3)
C      logical lerr
C
C      a(1,1) = 1d0
C      a(1,2) =-2d0
C      a(1,3) = 4d0
C      a(2,1) = 2d0
C      a(2,2) = 3d0
C      a(2,3) = 0d0
C      a(3,1) = 1d0
C      a(3,2) = 4d0
C      a(3,3) = 1d0
C
C      b(1) = -5
C      b(2) =  7
C      b(3) =  4
C
C      c(1) = 1
C      c(2) = 2
C      c(3) = 3
C
CC      call prmx('a',a,3,3,3)
CC      call prmx('b',b,3,3,1)
C
C      call dmpy31(1,a,b,c)
C      print 333, c
C
C      call dmpy31(0,a,b,c)
C      print 333, c
C
C      lerr = .false.
C      if (c(1) == -3 .and. c(2) == 11 .and. c(3) == 27) then
C        print *, 'product ok'
C      else
C        lerr = .true.
C        print *, 'product FAILED'
C      endif
C
C      call dmpy31(11,a,b,c)
C      print 333, c
C      if (c(1) == 13-3 .and. c(2) == 47+11 .and. c(3) == -16+27) then
C        print *, 'product ok'
C      else
C        lerr = .true.
C        print *, 'product FAILED'
C      endif
C
C      if (lerr) call fexit(-1,119,'test failed',0d0)
C
C  333 format(3f10.2)
C      end
