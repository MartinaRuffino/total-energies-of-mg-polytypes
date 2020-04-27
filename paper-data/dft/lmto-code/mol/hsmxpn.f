      subroutine hsmxpn(taub,taua,e,rsm,rsmg,nlmb,nlma,npwr,
     .   cg,indxcg,jcg,cy,c)
c  Makes the coefficients to expand the smooth Hankel functions at taub
c  into polynomials P_kL at point taua.
      parameter( lx=16, lpx=52 )
      implicit real*8 (a-h,p-z), integer(o)
      dimension c(0:npwr,nlma,nlmb),cy(1),cg(1),indxcg(1),jcg(1),dr(3),
     .  yl((lx+1)**2),phi(0:lx),taua(3),taub(3),f(0:lpx,0:lx),
     .  g(0:lpx,0:lx)

      fpi=16.d0*datan(1.d0)
      gamh=0.25d0*rsm*rsm
      gamg=0.25d0*rsmg*rsmg
      gam=gamh+gamg
      a=0.5d0/dsqrt(gam)
      do 1 m=1,3
    1 dr(m)=taub(m)-taua(m)
      lmaxa=ll(nlma)
      lmaxb=ll(nlmb)
      lmmax=lmaxa+lmaxb
      kmax=max0(lmaxa,lmaxb)+npwr
      if(lmmax > lx) call rx('hsmxpn: increase lx')
      if(kmax > lpx) call rx('hsmxpn: increase lpx')
      call sylm(dr,yl,lmmax,r2)
      do 2 ilm=1,(lmmax+1)**2
    2 yl(ilm)=yl(ilm)*cy(ilm)
      r=dsqrt(r2)
C --- Set up f(k,l) --------------
      call hansmr(r,e,a,phi,lmmax)
      call gauskl(r,e,a,g,lmmax,kmax,lpx)
      do 10 l=0,lmmax
      f(0,l)=phi(l)
      do 10 k=1,kmax
   10 f(k,l)=-e*f(k-1,l)-fpi*g(k-1,l)
C --- Combine with Clebsch-Gordan coefficients ---
      ee=dexp(-gamg*e)
      call dpzero(c,(npwr+1)*nlmb*nlma)
      do 11 ila=1,nlma
      la=ll(ila)
      do 11 ilb=1,nlmb
      lb=ll(ilb)
      ii=max0(ila,ilb)
      indx=(ii*(ii-1))/2+min0(ila,ilb)
      icg1=indxcg(indx)
      icg2=indxcg(indx+1)-1
      do 11 icg=icg1,icg2
      ilm=jcg(icg)
      lm=ll(ilm)
      k=(la+lb-lm)/2
      fac=ee*(-1d0)**lb*cg(icg)*yl(ilm)
      do 11 ip=0,npwr
   11 c(ip,ila,ilb)=c(ip,ila,ilb)+fac*f(k+ip,lm)

      end
