      subroutine hyfipl(d,d0,a1,cofa,cof,ncof,nalf)
C- Evaluates 2cf expansion to get 2cf at distance d
      implicit real*8 (a-h,p-z), integer (o)
      dimension cofa(ncof,nalf),cof(ncof)
      real w(1)
      common /w/ w
      call getpr(ipr)
      x=dexp(a1*(d0-d))
      if(ipr >= 60) write(6,440) d,d0,a1,x,nalf
  440 format(/' hyfipl:  d=',f9.6,'   d0=',f9.6,'   a1=',f7.4,
     .   '   x=',f9.5,'   nalf=',i2)
c ------ here for sum of Chebyshev polynomials -----
      y=2d0*x-1d0
      call defrr(owk,  2*ncof)
      call ropecs(y,nalf,ncof,w(owk),cofa,cof)
      call rlse(owk)
c ------ here for sum of monomials ---------
C|    do 10 icof=1,ncof
C|10  cof(icof)=0.d0
C|    xm=1.d0
C|    do 11 ialf=1,nalf
C|    xm=xm*x
C|    do 11 icof=1,ncof
C|11  cof(icof)=cof(icof)+xm*cofa(icof,ialf)
      end
