      subroutine ropech(c,nc,n,u2,f,wk)
c  Evaluates a chebyshev polynomial for a list of points.
      implicit real*8 (a-h,p-z)
      dimension c(nc),u2(n),wk(n,3),f(n)
c ------ case nc=1 -----------
      if(nc == 1) then
      do 1 i=1,n
  1   f(i)=0.5d0*c(1)
      endif
c ------ case nc=2 -----------
      if(nc == 2) then
      do 2 i=1,n
  2   f(i)=0.5d0*c(1)+(0.5d0*c(2))*u2(i)
      endif
c ------ case nc=3 -----------
      if(nc == 3) then
      do 3 i=1,n
      dd=u2(i)*c(3)+c(2)
  3   f(i)=0.5d0*u2(i)*dd+(0.5d0*c(1)-c(3))
      endif
      if(nc <= 3) return
c ------ unrolled for first two steps ---------
      do 5 i=1,n
      wk(i,1)=u2(i)*c(nc)+c(nc-1)
  5   wk(i,2)=u2(i)*wk(i,1)+(c(nc-2)-c(nc))
      ib=1
      ia=2
c ------ remaining steps ---------
      do 21 j=nc-3,2,-1
        ix=ia+1
        if(ix == 4) ix=1
        do 6 i=1,n
  6     wk(i,ix)=u2(i)*wk(i,ia)-wk(i,ib)+c(j)
        ib=ia
  21    ia=ix
      do 8 i=1,n
  8   f(i)=0.5d0*u2(i)*wk(i,ia)-wk(i,ib)+0.5d0*c(1)
      end
