      subroutine ftcdis(d0,ndist,adec,dist)
C- Set the ndist distances for making chebyshev coeffs.
c  Distances run from d0 to infinity; x runs from 0 to 1.
c  Transformation:  x=dexp(adec*(d0-d))
c  special branch: for ndist=1, sets distance to d0.
      implicit real*8 (a-h,p-z), integer (o)
      dimension dist(1)
      call getpr(ipr)
      if(ipr >= 20) write(6,200) adec,d0,ndist
  200 format(/' ftcdis:  adec=',f7.3,'   d0=',f7.3,'   ndist=',i3)
      dist(1)=d0
      if(ndist == 1) return
      if(ipr >= 30) write(6,782)
      pi=4d0*datan(1d0)
      do 10 k=1,ndist
      y=dcos(pi*(k-0.5d0)/ndist)
      x=0.5d0*(y+1d0)
      dist(k)=d0-dlog(x)/adec
      lp=0
      if(ipr >= 30.and.(k == 1.or.k == ndist)) lp=1
      if(ipr >= 40) lp=1
      if(lp == 1) write(6,781) k,y,x,dist(k)
  781 format(7x,i5,3f12.6)
  782 format(7x,'   id',6x,'y',11x,'x',11x,'dist')
  10  continue
      end
