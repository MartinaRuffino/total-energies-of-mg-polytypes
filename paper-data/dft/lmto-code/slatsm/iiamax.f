      integer function iiamax(n,idx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx <= 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      integer idx(*),imax
      integer i,incx,ix,n
c
      iiamax = 0
      if( n < 1 .or. incx <= 0 ) return
      iiamax = 1
      if(n == 1)return
      if(incx == 1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      imax = iabs(idx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(iabs(idx(ix)) <= imax) go to 5
         iiamax = i
         imax = iabs(idx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 imax = iabs(idx(1))
      do 30 i = 2,n
         if(iabs(idx(i)) <= imax) go to 30
         iiamax = i
         imax = iabs(idx(i))
   30 continue
      return
      end

