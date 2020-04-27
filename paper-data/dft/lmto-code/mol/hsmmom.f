      subroutine hsmmom(r,rsm,e,lmax,qg,qh)
c  Multipole moments of smooth Hankels and gaussians inside sphere.
c  The gaussian is that for energy zero. To get the moment of the
c  energy-dependent gaussian, multiply qg by dexp(e/(4*a*a)).
      implicit real*8 (a-h,p-z), integer (o)
      dimension qg(0:lmax),qh(0:lmax),xi(0:20),x0(0:20)
      pi=4d0*datan(1d0)
      a=1d0/rsm
      gam=1d0/(4d0*a*a)
      call hansmr(r,0d0,a,x0,lmax+1)
      call hansmr(r,e,a,xi,lmax+1)
      do 10 l=0,lmax
      qg(l)=r**(2*l+3)*x0(l+1)/(4d0*pi)
  10  qh(l)=r**(2*l+3)*( xi(l+1)-dexp(gam*e)*x0(l+1) ) / e
      end
