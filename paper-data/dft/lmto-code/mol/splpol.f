      subroutine splpol(nn,nrad,nlml,rmax,rhal,nr1,rofi,rhol,job)
c  Interpolate tail-rho to fine radial mesh by expanding
c  in Legendre polynomials up to order nn.
c  Set job=1 to multiply result by rofi**2.
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhal(nrad,nlml),rofi(nr1),rhol(nr1,nlml)
      real w(1)
      common /w/ w
      err0=1d-6
      call dpzero(rhol,     nr1*nlml)
c ------- make legendre coefficients -------------
      call defrr(ocof,    (nn+1)*nlml)
      call defrr(oxx,     nrad)
      call defrr(oww,     nrad)
      call mklegw(nrad,w(oxx),w(oww),0)
      call legprj(nrad,w(oxx),w(oww),rhal,nn,nlml,w(ocof),err)
      if(err > err0) write(6,331) err,err0
  331 format(' ---- splpol warning: err=',1p,d8.1,'  gt',d8.1)
      if(iprint() >= 40.or.err > err0) write(6,330) nrad,nn,err
  330 format(' splpol:  nrad=',i3,'   nn=',i3,1p,'   max err=',1p,d9.1)

c ------- evaluate legendre expansion to get rhol -----------
      call defrr(op,     nr1*3)
      call defrr(ox,     nr1)
      call xxlpol(nn,nlml,w(ocof),nr1,rofi,rhol,rho0,w(op),w(ox),rmax)
      call rlse(ocof)

c ------- multiply result by rofi**2 ---------------
      if(job /= 1) return
      do 30 ilm=1,nlml
      do 30 i=1,nr1
  30  rhol(i,ilm)=rhol(i,ilm)*rofi(i)**2
      end

c -------- sub xxlpol -------------------
      subroutine xxlpol(nn,nlml,cof,nr1,rofi,rhol,rho0,p,x,rmax)
      implicit real*8 (a-h,p-z)
      dimension cof(0:nn,nlml),rhol(nr1,nlml),p(nr1,3),x(1),rofi(1)
      do 5 i=1,nr1
  5   x(i)=(2d0/rmax)*rofi(i)-1d0
c ------- add terms p0,p1 ----------------------
      if(nn < 2) call rx('splpol: expect nn >= 2')
      do 1 i=1,nr1
      p(i,1)=1d0
  1   p(i,2)=x(i)
      do 10 ilm=1,nlml
      do 2 i=1,nr1
  2   rhol(i,ilm)=rhol(i,ilm)+cof(0,ilm)*p(i,1)+cof(1,ilm)*p(i,2)
  10  continue
      i1=1
      i2=2
c ------- add terms p2 to pnn ------------------
      do 20 m=2,nn
      i3=i2+1
      if(i3 == 4) i3=1
      do 3 i=1,nr1
  3   p(i,i3)=((2*m-1d0)/m)*x(i)*p(i,i2) - ((m-1d0)/m)*p(i,i1)
      do 21 ilm=1,nlml
      do 21 i=1,nr1
  21  rhol(i,ilm)=rhol(i,ilm)+cof(m,ilm)*p(i,i3)
      i1=i2
  20  i2=i3
      rho0=rhol(1,1)
      end
