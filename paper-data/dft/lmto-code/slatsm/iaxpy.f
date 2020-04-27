      subroutine iaxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      integer dx(*),dy(*),da
      integer i,incx,incy,ix,iy,n
c
      if (n <= 0  .or.  da == 0) return
      ix = 1
      iy = 1
      if (incx < 0) ix = (1-n)*incx + 1
      if (incy < 0) iy = (1-n)*incy + 1
      if (incx == 1 .and. incy == 1) then
        do  10  ix = 1, n
   10   dy(ix) = dy(ix) + da*dx(ix)
        return
      endif
      do  20  i = 1, n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   20 continue
      end
