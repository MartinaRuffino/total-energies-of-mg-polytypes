      subroutine yydot(n,dx,dxi,incx,dy,dyi,incy,dotr,doti)
c
c     forms the dot product of two complex vectors.
c     Adapted from: jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dxi(1),dy(1),dyi(1),dotr,doti
      integer i,incx,incy,ix,iy,n
c
C#ifdefC APOLLO | HP
C      double precision vec_$ddot,vec_$ddot_i
C      dotr = 0d0
C      doti = 0d0
C      if (n <= 0) return
C      if (incx == 1 .and. incy == 1) then
C        dotr = vec_$ddot(dx,dy,n)  - vec_$ddot(dxi,dyi,n)
C        doti = vec_$ddot(dx,dyi,n) + vec_$ddot(dxi,dy,n)
C      else
C        dotr = vec_$ddot_i(dx,incx,dy,incy,n)
C     .        -vec_$ddot_i(dxi,incx,dyi,incy,n)
C        doti = vec_$ddot_i(dx,incx,dyi,incy,n)
C     .        +vec_$ddot_i(dxi,incx,dy,incy,n)
C      endif
C#else
      dotr = 0d0
      doti = 0d0
      ix = 1
      iy = 1
      if (incx < 0) ix = (1-n)*incx + 1
      if (incy < 0) iy = (1-n)*incy + 1
      do  10  i = 1, n
        dotr = dotr + dx(ix)*dy(iy) - dxi(ix)*dyi(iy)
        doti = doti + dx(ix)*dyi(iy) + dxi(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
C#endif
      end
