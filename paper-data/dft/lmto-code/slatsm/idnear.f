      integer function idnear(n,da,dx,incx)
C- Finds the index of element closest to specified value
      integer i,incx,ix,n
      double precision dx(n),dmax,da

c
      idnear = 0
      if ( n < 1 ) return
      idnear = 1
      if ( n == 1 ) return
      ix = 1
      dmax = dabs(dx(1)-da)
      if (dmax == 0) return
      ix = ix + incx
      do  10  i = 2, n
        if (dabs(dx(ix)-da) < dmax) then
          idnear = i
          dmax = dabs(dx(ix)-da)
          if (dmax == 0) return
        endif
        ix = ix + incx
   10 continue
      end
