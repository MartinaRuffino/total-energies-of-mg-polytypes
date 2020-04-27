      subroutine mstrux(e,dr,s,nlma,nlmb,ndim,cg,indxcg,jcg,cy)
C- Mol strux, standard definition iv-43.
C  Exp. theorem is:  k(k,r-dr) = sum(m) s(m,k,dr)*j(m,r)
      implicit real*8 (a-h,p-z), integer(o)
      parameter (lx=14)
      dimension s(ndim,nlmb),cy(1),cg(1),indxcg(1),jcg(1),
     .  hl((lx+1)**2),psi(0:lx),phi(0:lx),dr(3),efac(0:lx),sig(0:lx)
      if(nlma > ndim) call rx('mstrux: nlma gt ndim')
      fourpi=16.d0*datan(1.d0)
      lmax=ll(nlma)+ll(nlmb)
      if(lmax > lx) call rx('mstrux: lmax gt lx')
      call sylm(dr,hl,lmax,r2)
      if(r2 < 1.d-10) then
        do 5 jlm=1,nlmb
        do 5 ilm=1,nlma
  5     s(ilm,jlm)=0.d0
        return
        endif
      call bessl(e*r2,lmax,phi,psi)
      ilm=0
      rfac=dsqrt(r2)
      xxx=1.d0/r2
      do 10 l=0,lmax
      rfac=rfac*xxx
      do 10 m=1,2*l+1
      ilm=ilm+1
  10  hl(ilm)=rfac*psi(l)*cy(ilm)*hl(ilm)
      efac(0)=1.d0
      sig(0)=1.d0
      do 1 l=1,lmax
      efac(l)=-e*efac(l-1)
  1   sig(l)=-sig(l-1)
c ---------------------------------------
      do 11 mlm=1,nlma
      lm=ll(mlm)
      do 11 klm=1,nlmb
      lk=ll(klm)
      sum=0.d0
      ii=max0(mlm,klm)
      indx=(ii*(ii-1))/2+min0(mlm,klm)
      icg1=indxcg(indx)
      icg2=indxcg(indx+1)-1
      do 21 icg=icg1,icg2
      llm=jcg(icg)
      ipow=(lm+lk-ll(llm))/2
  21  sum=sum+cg(icg)*efac(ipow)*hl(llm)
  11  s(mlm,klm)=fourpi*sum*sig(lk)
      return
      end
