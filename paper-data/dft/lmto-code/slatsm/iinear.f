      integer function iinear(n,da,dx,incx)
C- Finds the index of element closest to specified value
      integer dx(*),dmax,da
      integer i,incx,ix,n
c
      iinear = 0
      if ( n < 1 ) return
      iinear = 1
      if ( n == 1 ) return
      ix = 1
      dmax = iabs(dx(1)-da)
      if (dmax == 0) return
      ix = ix + incx
      do  10  i = 2, n
        if (iabs(dx(ix)-da) < dmax) then
          iinear = i
          dmax = iabs(dx(ix)-da)
          if (dmax == 0) return
        endif
        ix = ix + incx
   10 continue
      end
