      subroutine gintz(g1,g2,a,b,nr,z,e,l,v,rofi,sum)
C- Integrate inner product of two wave equations
C ----------------------------------------------------------------
Ci Inputs:
Ci   g1,g2 :First and second radial wave functions (complex)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   z     :nuclear charge
Ci   e     :energy (complex)
Ci   l     :l quantum number of g1,g2
Ci   v     :spherical potential
Ci   rofi  :radial mesh points
Co Outputs:
Co   sum   :inner product
Cr Remarks:
Cr   This is a complex analog of gintsr
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,l
      double precision a,b,z,v(nr),rofi(nr)
      double complex e,g1(nr,2),g2(nr,2),sum
C ... Local parameters
      integer i,ir
      double precision fllp1,c,r
      double complex tmc,fi
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c
      tmc(i,r) = c - (v(i) - 2d0*z/r - e)/c
      fi(i,r) = (r+b)*(
     .  dconjg(g1(i,1))*g2(i,1)*(1+fllp1/(tmc(i,r)*r)**2) +
     .  dconjg(g1(i,2))*g2(i,2))

      fllp1 = l*(l+1)
      sum = 0d0
      do  10  ir = 2, nr-1, 2
   10 sum = sum + fi(ir,rofi(ir))
      sum = 2*sum
      do  11  ir = 3, nr-2, 2
   11 sum = sum + fi(ir,rofi(ir))
      sum = (2*sum + fi(nr,rofi(nr)))*a/3

      end
