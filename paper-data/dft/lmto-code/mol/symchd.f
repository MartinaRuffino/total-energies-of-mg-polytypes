      subroutine symchd(lmaxa,nlml,nr,u1,u2,rhol,q,
     .   ipbas,ioffp,sym,nrc,cg,indxcg,jcg,thet,wk)
c  Adds a contribution to symmetrized charge density.
c  thet and wk are work spaces.
      implicit real*8 (a-h,p-z), integer(o)
      dimension thet(nlml,nrc,0:lmaxa,0:lmaxa),wk(nr),ioffp(1),
     .   u1(nr,0:lmaxa),u2(nr,0:lmaxa),rhol(nr,nlml),ipbas(1),
     .   cg(1),indxcg(1),jcg(1),q(1),sym(nlml,nlml,1)
      nlma=(lmaxa+1)**2
c ------ make help coefficients in thet ------------
      call dpzero(thet,   nlml*nrc*(lmaxa+1)**2)
      iq=0
      do 30 jlm=1,nlma
      l2=ll(jlm)
      do 30 ilm=1,nlma
      l1=ll(ilm)
      iq=iq+1
      ind=max0(ilm,jlm)
      ind=(ind*(ind-1))/2+min0(ilm,jlm)
      do 32 icg=indxcg(ind),indxcg(ind+1)-1
      mlm=jcg(icg)
      if(mlm <= nlml) then
       do 5 ja=1,nrc
       jb=ipbas(ja)
  5    thet(mlm,ja,l1,l2)=thet(mlm,ja,l1,l2)+cg(icg)*q(iq+ioffp(jb))
       endif
  32  continue
  30  continue
c ------ combine thet with projectors, add to rhol ----
      do 36 l1=0,lmaxa
      do 36 l2=0,lmaxa
      do 33 ir=1,nr
  33  wk(ir)=u1(ir,l1)*u2(ir,l2)
      do 36 klm=1,nlml
      sum=0d0
      do 38 ja=1,nrc
      do 38 mlm=1,nlml
  38  sum=sum+sym(klm,mlm,ja)*thet(mlm,ja,l1,l2)
      if(dabs(sum) > 1d-9) then
        do 39 ir=1,nr
  39    rhol(ir,klm)=rhol(ir,klm)+sum*wk(ir)
        endif
  36  continue
      end
