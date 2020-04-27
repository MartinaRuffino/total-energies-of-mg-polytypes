      subroutine qsmuz3(rsm,nxi,lxi,exi,n0,nbas,ips,rhoi,rho0,qsm)
c  integral over all space of smooth density
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(n0,1),nxi(1),lxi(n0,1),ips(1),rhoi(1),rsm(1),
     .   rho0(1)
      pi=4d0*datan(1d0)
      srfpi=dsqrt(4d0*pi)
      qsm=0d0
c ------ add together smooth charge by atoms -----------
      i=1
      j=1
      do 50 ib=1,nbas
      is=ips(ib)
      gam=0.25d0*rsm(is)*rsm(is)
      qh=0d0
      lmax=0
      do 51 ie=1,nxi(is)
      lmax=max0(lmax,lxi(ie,is))
      e=exi(ie,is)
      qh=qh+rhoi(i)*srfpi*dexp(gam*e)/(-e)
  51  i=i+(lxi(ie,is)+1)**2
      qg=rho0(j)/srfpi
      j=j+(lmax+1)**2
      write(6,821) ib,qh,qg
  821 format(' qsmuz3:   ib=',i5,'   qh,qg=',2f12.6)
  50  qsm=qsm+qh+qg
      write(6,763) qsm
  763 format(' qsm=',f12.6)
      end
