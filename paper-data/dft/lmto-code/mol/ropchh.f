      subroutine ropchh(e,rsm,lmax,r1,r2,atr,n,c,tol,n1,kpr)
c  Sets up chebyshev polys for smooth hankel functions, xi=h*/r**l.
c  Fit is good from r1 to r2. Transformation is defined by atr.
c  Returns n1 so that all higher coeffs are smaller than tol.
      implicit real*8 (a-h,p-z), integer(o)
      dimension xi(20),c(n,lmax+1)
      asm=1d0/rsm
      pi=4d0*datan(1d0)
      call ropg0(atr,r1,y1)
      call ropg0(atr,r2,y2)
      bma=0.5d0*(y2-y1)
      bpa=0.5d0*(y2+y1)
c ------- make coeffs -----------------
      do 1 m=1,lmax+1
      do 1 k=1,n
  1   c(k,m)=0d0
      fac=2d0/n
      do 12 k=1,n
        u=dcos(pi*(k-0.5d0)/n)
        y=u*bma+bpa
        call ropgi(atr,y,r)
        call hansmr(r,e,asm,xi,lmax)
        do 12 j=1,n
          xxx=fac*dcos((pi*(j-1))*((k-0.5d0)/n))
          do 12 m=1,lmax+1
  12      c(j,m)=c(j,m)+xxx*xi(m)
c ------- find n1 ---------------------
      do 30 k=n,1,-1
      top=0d0
      do 31 m=1,lmax+1
  31  top=dmax1(top,dabs(c(k,m)))
      n1=k
      if(top > tol) goto 32
  30  continue
  32  continue
      n1=max0(n1,2)
c ------- print coeffs --------------
      if(kpr <= 0) return
      write(6,717) r1,r2,tol,n,n1
  717 format(' ropchh:  r1,r2=',2f8.3,'   tol=',1p,d10.3,
     .  '    n0=',i5,'    n1=',i5)
      do 20 k=1,n1
  20  write(6,220) k,(c(k,m),m=1,lmax+1)
  220 format(i4,7f10.4/(4x,7f10.4))
      end
