      subroutine ropj(r,e,lmax,y,h,b,n)
C- Vector of bessel functions
      implicit none
      integer n,lmax
      double precision r(1),e(1),y(1),h(1),b(n,0:1)
      integer ir,l
      call ropbes(r,e,lmax,y,h,b,n,1)
      h(1) = 0d0
      if (r(1) /= 0) h(1) = 1/r(1)
      do  1  ir = 2, n
    1 h(ir) = 1/r(ir)
      do  12  l = 0, lmax
      do  2  ir = 1, n
      h(ir) = h(ir)*r(ir)
    2 b(ir,l) = b(ir,l)*h(ir)
   12 continue
      end
