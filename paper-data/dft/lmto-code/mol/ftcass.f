      subroutine ftcass(ncof,ndist,nalf,dist,adec,d0,che,cof)
C- re-assemble coefficients from cheby fit
      implicit real*8 (a-h,p-z), integer (o)
      dimension dist(120),che(ncof,nalf),cof(ncof,ndist)
      real w(1)
      common /w/ w
      call defrr(owk,  2*ncof)
C|    write(6,330) ndist,dist(1),dist(ndist)
C|330 format(' ftcass:  ndist=',i4,'   d:',2f10.4)
      do 10 id=1,ndist
      d=dist(id)
      x=dexp(adec*(d0-d))
      y=2d0*x-1d0
      call ropecs(y,nalf,ncof,w(owk),che,cof(1,id))
  10  continue
      call rlse(owk)
      end
