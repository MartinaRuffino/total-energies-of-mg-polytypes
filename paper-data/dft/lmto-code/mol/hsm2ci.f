      subroutine hsm2ci(tau1,tau2,e1,e2,rsm1,rsm2,s,nlm1,nlm2,ndim,
     .   cg,indxcg,jcg,cy)
c  Two-center integrals of smoothed hankel fcts, all cases
      parameter( lx=16 )
      implicit real*8 (a-h,p-z), integer(o)
      dimension s(ndim,1),cy(1),cg(1),indxcg(1),jcg(1),dr(3),
     .   psi1(0:lx+1),psi2(0:lx),phi(0:lx),tau1(3),tau2(3),
     .   f(0:lx,0:lx),g(0:lx,0:lx),xip(0:lx),yl((lx+1)**2)
      fpi=16.d0*datan(1.d0)
      gam1=0.25d0*rsm1*rsm1
      gam2=0.25d0*rsm2*rsm2
      gam=gam1+gam2
      a=0.5d0/dsqrt(gam)
      leq=0
      if(dabs(e1-e2) < 1.d-5) leq=1
      do 1 m=1,3
    1 dr(m)=tau2(m)-tau1(m)
      l1=ll(nlm1)
      l2=ll(nlm2)
      lmax=l1+l2
      kmax=max0(l1,l2)
      if(lmax > lx) call rx('hsm2ci: change dimensions''')
      call sylm(dr,yl,lmax,r2)
      do 2 ilm=1,(lmax+1)**2
  2   yl(ilm)=yl(ilm)*cy(ilm)
      r=dsqrt(r2)
c --------- set up f(k,l) for case e1 /= e2 --------------
      if(leq == 0) then
      fac1=dexp(gam2*(e2-e1))/(e1-e2)
      fac2=dexp(gam1*(e1-e2))/(e2-e1)
      call hansmr(r,e1,a,phi,lmax)
      call gauskl(r,e1,a,g,lmax,kmax,lx)
      do 10 l=0,lmax
      hkl=phi(l)
      f(0,l)=hkl*fac1
      do 10 k=1,kmax
      hkl=-e1*hkl-fpi*g(k-1,l)
  10  f(k,l)=hkl*fac1
      call hansmr(r,e2,a,phi,lmax)
      call gauskl(r,e2,a,g,lmax,kmax,lx)
      do 14 l=0,lmax
      hkl=phi(l)
      f(0,l)=f(0,l)+hkl*fac2
      do 14 k=1,kmax
      hkl=-e2*hkl-fpi*g(k-1,l)
  14  f(k,l)=f(k,l)+hkl*fac2
      endif
c --------- set up f(k,l) for case e1 == e2 --------------
      if(leq == 1) then
      e=0.5d0*(e1+e2)
      akap=dsqrt(-e)
      uminus=derfc(0.5d0*akap/a-r*a)*dexp(-akap*r)
      uplus =derfc(0.5d0*akap/a+r*a)*dexp(+akap*r)
      call hansmr(r,e,a,phi,lmax)
      call gauskl(r,e,a,g,lmax,kmax,lx)
      hklp=(uminus+uplus)/(4d0*akap)
      do 20 l=0,lmax
      hkl= phi(l)
      f(0,l)=hklp-gam*hkl
      do 19 k=1,kmax
      hklp=-e*hklp-fpi*gam*g(k-1,l)-hkl
      hkl= -e*hkl-fpi*g(k-1,l)
  19  f(k,l)=hklp-gam*hkl
  20  hklp=0.5d0*phi(l)
      endif
c -------- combine with clebsch-gordan coefficients ----------
      do 11 mlm=1,nlm1
      lm=ll(mlm)
      do 11 klm=1,nlm2
      lk=ll(klm)
      sum=0.d0
      ii=max0(mlm,klm)
      indx=(ii*(ii-1))/2+min0(mlm,klm)
      icg1=indxcg(indx)
      icg2=indxcg(indx+1)-1
      do 21 icg=icg1,icg2
      jlm=jcg(icg)
      lj=ll(jlm)
      k=(lm+lk-lj)/2
  21  sum=sum+cg(icg)*yl(jlm)*f(k,lj)
  11  s(mlm,klm)=fpi*sum* (-1d0)**lk
      return
      end
