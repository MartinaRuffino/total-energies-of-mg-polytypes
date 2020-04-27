      subroutine ropecs(y,m,n,wk,c,f)
C- evaluate vector of Chebyshev fits, all at argument y
      implicit real*8 (a-h,p-z), integer (o)
      dimension c(n,m),f(n),wk(n,2)
      y2=2d0*y
c ------ case m=1 -----------
      if(m == 1) then
      do 1 i=1,n
  1   f(i)=0.5d0*c(i,1)
      endif
c ------ case m=2 -----------
      if(m == 2) then
      do 2 i=1,n
  2   f(i)=0.5d0*c(i,1)+y*c(i,2)
      endif
c ------ case m=3 -----------
      if(m == 3) then
      xx=y*y2-1d0
      do 3 i=1,n
  3   f(i)=xx*c(i,3)+y*c(i,2)+0.5d0*c(i,1)
      endif
      if(m <= 3) return
c ------ unrolled for first three steps ---------
      do 4 i=1,n
  4   wk(i,1)=y2*c(i,m)+c(i,m-1)
      do 5 i=1,n
  5   wk(i,2)=y2*wk(i,1)+c(i,m-2)-c(i,m)
      ib=1
      ia=2
c ------ remaining steps ---------
      do 21 j=m-3,2,-1
        do 6 i=1,n
  6     wk(i,ib)=y2*wk(i,ia)-wk(i,ib)+c(i,j)
        ix=ia
        ia=ib
  21    ib=ix
      do 8 i=1,n
  8   f(i)=y*wk(i,ia)-wk(i,ib)+0.5d0*c(i,1)
      end
