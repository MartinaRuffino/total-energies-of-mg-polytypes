      subroutine ftcmch(ncof,ndist,nalf,cof,che)
C- Make chebyshev fit to cof
      implicit real*8 (a-h,p-z), integer (o)
      dimension cof(ncof,ndist),che(ncof,nalf)
      call getpr(ipr)
      pi=4d0*datan(1d0)
      if(ipr >= 20) write(6,221) nalf
  221 format(' ftcmch:  make chebyshev coeffs,   nalf=',i5)
      fac=2d0/ndist
      call dpzero(che,ncof*nalf)
      do 12 j=1,nalf
      do 12 k=1,ndist
      wkj=fac*dcos((pi*(j-1))*((k-0.5d0)/ndist))
      do 15 i=1,ncof
  15  che(i,j)=che(i,j)+wkj*cof(i,k)
  12  continue
      end
