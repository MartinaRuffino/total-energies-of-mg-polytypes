      subroutine hyfipd(d,d0,a1,cofa,cof,dof,ncof,nalf)
c  Interpolate to get coeffs at distance d, and derivative
      implicit real*8 (a-h,p-z), integer (o)
      dimension cofa(ncof,nalf),cof(ncof),dof(ncof)
      real w(1)
      common /w/ w
      call getpr(ipr)
      x=dexp(a1*(d0-d))
      if(ipr >= 70) write(6,440) d,d0,a1,x,nalf
  440 format(/' hyfipd:  d=',f9.6,'   d0=',f9.6,'   a1=',f7.4,
     .   '   x=',f9.5,'   nalf=',i2)
c ------ here for sum of Chebyshev polynomials -----
      y=2d0*x-1d0
      call defrr(owk,  2*ncof)
      call defrr(odc,  2*ncof)
      call ropecs(y,nalf,ncof,w(owk),cofa,cof)
      call ropecd(y,nalf,ncof,w(owk),w(odc),cofa,dof)
      xx=-2d0*a1*x
      do 1 i=1,ncof
  1   dof(i)=dof(i)*xx
      call rlse(owk)
c ------ here for sum of monomials ---------
C|    dxdd=-a1*x
C|    do 10 i=1,ncof
C|    dof(i)=0d0
C|10  cof(i)=0.d0
C|    xm=1.d0
C|    do 11 ialf=1,nalf
C|    ym=dxdd*ialf*xm
C|    xm=xm*x
C|    do 11 i=1,ncof
C|    dof(i)=dof(i)+ym*cofa(i,ialf)
C|11  cof(i)=cof(i)+xm*cofa(i,ialf)

      end
