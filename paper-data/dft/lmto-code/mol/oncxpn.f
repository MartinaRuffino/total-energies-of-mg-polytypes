      subroutine oncxpn(lx,lmax,nx,nf,lp1,lp2,nlm1,nlm2,cof,x,
     .   cg,jcg,indxcg)
c  on-site expansion of products of smoothed hankels
      implicit real*8 (a-h,p-z), integer(o)
      dimension cof(0:lmax,nx,0:lp1,0:lp2),x(nf,nlm1,nlm2),lx(8),
     .   cg(1),jcg(1),indxcg(1),ix0(20)
      real w(1)
      common /w/ w
      call getpr(ipr)
c ------- offsets to help find the xi's ----------
      ix0(1)=0
      do 17 k=1,nx-1
  17  ix0(k+1)=ix0(k)+(lx(k)+1)**2
c ------ start loop over the phi1,phi2 ----------
      call dpzero(x,   nf*nlm1*nlm2)
      do 10 ilm1=1,nlm1
      l1=ll(ilm1)
      do 10 ilm2=1,nlm2
      l2=ll(ilm2)
      iii=max0(ilm1,ilm2)
      iii=(iii*(iii-1))/2+min0(ilm1,ilm2)
      icg1=indxcg(iii)
      icg2=indxcg(iii+1)-1
      do 11 icg=icg1,icg2
      mlm=jcg(icg)
      lm=ll(mlm)
      do 12 ix=1,nx
      if(lm <= lx(ix)) then
        jf=mlm+ix0(ix)
        x(jf,ilm1,ilm2)=cg(icg)*cof(lm,ix,l1,l2)
        endif
  12  continue
  11  continue
  10  continue
      return
      end
