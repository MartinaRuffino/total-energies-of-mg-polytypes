      subroutine qsmuz(rsm,nxi,lxi,exi,n0,nbas,ips,rhoi,qsm,asmag)
c  integral over all space of smooth density
C  28 Dec 94 made spin-pol (MvS)
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(n0,1),nxi(1),lxi(n0,1),ips(1),rhoi(1),rsm(1),q1(2)
      pi=4d0*datan(1d0)
      srfpi=dsqrt(4d0*pi)
      qsm=0d0
      asmag=0d0
c ------ add together smooth charge by atoms -----------
      nsp = lsp()+1
      i=1
      do 50 isp=1,nsp
      do 50 ib=1,nbas
      is=ips(ib)
      gam=0.25d0*rsm(is)*rsm(is)
      q1(isp)=0d0
      do 51 ie=1,nxi(is)
      e=exi(ie,is)
      q1(isp)=q1(isp)+rhoi(i)*srfpi*dexp(gam*e)/(-e)
  51  i=i+(lxi(ie,is)+1)**2
      if(iprint() >= 50 .and. nsp == 1) print 821, ib,q1(1)
      if(iprint() >= 50 .and. isp == 2)
     .  print 821, ib,q1(1)+q1(2),q1(2)-q1(1)
  821 format(' ib=',i5,'   qsm(ib)=',f12.6:'   smooth mag. mom=',f12.6)
      asmag=asmag+q1(isp)*(2*isp-3)
  50  qsm=qsm+q1(isp)
      if (nsp == 1) asmag=0d0
      if(iprint() >= 30 .and. nsp == 1) write(6,763) qsm
      if(iprint() >= 30 .and. nsp == 2) write(6,763) qsm,asmag
  763 format(' qsmuz:  qsm=',f12.6:'   smoothed mag. mom=',f12.6)
      end
