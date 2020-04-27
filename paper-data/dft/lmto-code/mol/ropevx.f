      subroutine ropevx(rho,exc,vxc,n)
c  Ceperly-Alder xc energy density and potential, vectorizes.
c  Sets exc,vxc to zero if rho is negative.
      implicit real*8 (a-h,p-z)
      dimension rho(n),exc(n),vxc(n)
      parameter( pi=3.141592653589793238d0 )
      parameter( au=0.0311d0, bu=-0.048d0, cu=0.002d0, du=-0.0116d0,
     .   aph=-0.1423d0, bt1=1.0529d0, bt2=0.3334d0 )
      parameter( aa=-3d0/(4d0*pi), bu1=bu-au/3d0,
     .   cu1=2d0*cu/3d0, du1=(2d0*du-cu)/3d0,
     .   p13=1d0/3d0, pi4=4d0*pi, c76=7d0/6d0,
     .   c43=4d0/3d0 )
      fac=2d0
      ax=aa*(9d0*pi/4d0)**p13
      do 11 i=1,n
        rr=dmax1(rho(i),0d0)
        xs=((pi4*rr)/3d0)**p13
        rs=1d0/dmax1(xs,1d0)
        xxs=dsqrt(xs)
        dpr=1d0/(xs+bt1*xxs+bt2)
        e1=fac*aph*dpr*xs
        u1=(c76*bt1*xxs+xs+c43*bt2)*e1*dpr
        rsl=dlog(rs)
        e2=fac*(du*rs+bu+au*rsl+cu*rs*rsl)
        u2=fac*(du1*rs+bu1+au*rsl+cu1*rs*rsl)
        xx=(1d0-xs)/dabs(1d0-xs)
        xx=0.5d0*(xx+1d0)
        ex=xs*(ax*fac)
        vx=xs*(c43*ax*fac)
        exc(i)=  xx*e1+(1d0-xx)*e2   +   ex
        vxc(i)=  xx*u1+(1d0-xx)*u2   +   vx
  11  continue
      end
