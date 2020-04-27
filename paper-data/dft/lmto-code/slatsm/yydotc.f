      subroutine yydotc(n,dx,dxi,incx,dy,dyi,incy,dotr,doti)
c
c     forms the dot product of two complex vectors (first complex conj)
c
      implicit none
      double precision dx(*),dxi(*),dy(*),dyi(*),dotr,doti
      integer i,incx,incy,ix,iy,n

      dotr = 0d0
      doti = 0d0
      ix = 1
      iy = 1
      if (incx < 0) ix = (1-n)*incx + 1
      if (incy < 0) iy = (1-n)*incy + 1
      do  i = 1, n
        dotr = dotr + dx(ix)*dy(iy) + dxi(ix)*dyi(iy)
        doti = doti + dx(ix)*dyi(iy) - dxi(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      enddo
      end
