      subroutine hsg2ci(tauh,taug,e,rsm,rsmg,nlmh,nlmg,
     .  cg,indxcg,jcg,cy,s)
C- Integral of smoothed Hankels times Gaussians, centered on separate
C  sites.
      parameter( lx=16, lpx=12 )
      implicit real*8 (a-h,p-z), integer(o)
      dimension s(nlmh,nlmg),cy(1),cg(1),indxcg(1),jcg(1),dr(3),
     .  yl((lx+1)**2),phi(0:lx),taug(3),tauh(3),f(0:lpx,0:lx),
     .  g(0:lpx,0:lx)
      fpi=16.d0*datan(1.d0)
      gamh=0.25d0*rsm*rsm
      gamg=0.25d0*rsmg*rsmg
      gam=gamh+gamg
      a=0.5d0/dsqrt(gam)
      do 1 m=1,3
    1 dr(m)=tauh(m)-taug(m)
      lmaxg=ll(nlmg)
      lmaxh=ll(nlmh)
      lmmax=lmaxg+lmaxh
      kmax=max0(lmaxg,lmaxh)
      if(lmmax > lx) call rx('hsg2ci: increase lx')
      if(kmax > lpx) call rx('hsg2ci: increase lpx')
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
      call dpzero(s,nlmh*nlmg)
      do 11 ilg=1,nlmg
      lg=ll(ilg)
      do 11 ilh=1,nlmh
      lh=ll(ilh)
      ii=max0(ilg,ilh)
      indx=(ii*(ii-1))/2+min0(ilg,ilh)
      icg1=indxcg(indx)
      icg2=indxcg(indx+1)-1
      do 11 icg=icg1,icg2
      ilm=jcg(icg)
      lm=ll(ilm)
      k=(lg+lh-lm)/2
      fac=ee*(-1d0)**lh*cg(icg)*yl(ilm)
   11 s(ilh,ilg)=s(ilh,ilg)+fac*f(k,lm)

      end
