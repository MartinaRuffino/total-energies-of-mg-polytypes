      subroutine ropecd(y,m,n,wk,dc,c,f)
C- Evaluate vector of derivatives of Cheb fits, all at argument y
C  Assumes interval is (0,1); otherwise multiply by 1/(b-a) (?)
C  dc and wk are work spaces.
      implicit real*8 (a-h,p-z), integer (o)
      dimension c(n,m),f(n),wk(n,2),dc(n,2)
      y2=2d0*y
c ------ case m=1 -----------
      if(m == 1) then
      do 1 i=1,n
  1   f(i)=0.d0
      endif
c ------ case m=2 -----------
      if(m == 2) then
      do 2 i=1,n
  2   f(i)=c(i,2)
      endif
c ------ case m=3 -----------
      if(m == 3) then
      do 3 i=1,n
  3   f(i)=(4d0*y)*c(i,3)+c(i,2)
      endif
      if(m <= 3) return
c ------ unrolled for first three steps ---------
      do 4 i=1,n
      dc(i,1)=(2d0*(m-1)) * c(i,m)
  4   wk(i,1)=dc(i,1)
      do 5 i=1,n
      dc(i,2)=(2d0*(m-2)) * c(i,m-1)
  5   wk(i,2)=y2*wk(i,1) + dc(i,2)
      ib=1
      ia=2
c ------ remaining steps ---------
      do 21 j=m-3,2,-1
        do 6 i=1,n
        dc(i,ib)=dc(i,ib)+(2d0*j)*c(i,j+1)
  6     wk(i,ib)=y2*wk(i,ia)-wk(i,ib)+dc(i,ib)
        ix=ia
        ia=ib
  21    ib=ix
      do 8 i=1,n
      dc(i,ib)=dc(i,ib)+(2d0*j)*c(i,j+1)
  8   f(i)=y*wk(i,ia)-wk(i,ib)+0.5d0*dc(i,ib)
      end
