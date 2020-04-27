      subroutine mstrxp(e,dr,s,sp,nlma,nlmb,ndim,cg,indxcg,jcg,cy)
C- Mol strux and energy derivatives, standard definition iv-43.
C  Exp. theorem is:  k(k,r-dr) = sum(m) s(m,k,dr)*j(m,r)
      implicit real*8 (a-h,p-z), integer(o)
      parameter (lx=14)
      dimension s(ndim,nlmb),cy(1),cg(1),indxcg(1),jcg(1),
     .  hl((lx+2)**2),hlp((lx+2)**2),dr(3),
     .  psi(0:lx+1),phi(0:lx+1),efac(0:lx+1),sig(0:lx+1),sp(ndim,nlmb)
      if(nlma > ndim) call rx('mstrux: nlma gt ndim')
      if(nlma <= 0.or.nlmb <= 0) return
      fpi=16.d0*datan(1.d0)
      lmax=ll(nlma)+ll(nlmb)
      if(lmax > lx) call rx('mstrux: lmax gt lx')
      call sylm(dr,hl,lmax,r2)
      call bessl(e*r2,lmax+1,phi,psi)
      ilm=0
      rfac=dsqrt(r2)
      xxx=1.d0/r2
      do 10 l=0,lmax
      psidot = ((l+l+1)*psi(l)-psi(l+1))/(e+e)
      rfac=rfac*xxx
      do 10 m=1,2*l+1
      ilm=ilm+1
      hlp(ilm)=rfac*psidot*cy(ilm)*hl(ilm)
  10  hl(ilm)=rfac*psi(l)*cy(ilm)*hl(ilm)
      efac(0)=1.d0
      sig(0)=1.d0
      do 1 l=1,lmax
      efac(l)=-e*efac(l-1)
  1   sig(l)=-sig(l-1)
c ---------------------------------------
      ube=1.d0/e
      do 11 mlm=1,nlma
      lm=ll(mlm)
      do 11 klm=1,nlmb
      lk=ll(klm)
      sup=0.d0
      sum=0.d0
      ii=max0(mlm,klm)
      indx=(ii*(ii-1))/2+min0(mlm,klm)
      icg1=indxcg(indx)
      icg2=indxcg(indx+1)-1
      do 21 icg=icg1,icg2
      llm=jcg(icg)
      ipow=(lm+lk-ll(llm))/2
      sum=sum+cg(icg)*efac(ipow)*hl(llm)
  21  sup=sup+cg(icg)*efac(ipow)*(hlp(llm)+ipow*ube*hl(llm))
      s(mlm,klm)=fpi*sum*sig(lk)
  11  sp(mlm,klm)=fpi*sup*sig(lk)
      return
      end
