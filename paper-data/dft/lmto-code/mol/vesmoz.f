      subroutine vesmoz(nxi,lxi,exi,n0,rmt,rint,lmxl,nbas,alat,pos,
     .  ips,rhoi,ioff,cg,jcg,indxcg,cy,ioffv0,nvi,ixp,dxp,qmoi)
c  Multipole moments of interstitial tail density inside spheres
c  28 Dec 94 spin pol (MvS)
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(n0,1),nxi(1),ips(1),lxi(n0,1),ioff(1),ioffv0(1),
     .  rhoi(1),pos(3,1),rint(1),cg(1),jcg(1),indxcg(1),cy(1),
     .  rmt(1),lmxl(1),qmoi(nvi),dr(3),ojkl(0:20),ixp(1),dxp(1)
      real w(1)
      common /w/ w
      y0=1d0/dsqrt(16d0*datan(1d0))
      call getpr(ipr)
      call tcn('vesmoz')
      call dpzero(qmoi,nvi)
      nsp = lsp()+1
      nri = ioff(nbas+1)
      if (nsp == 2) then
        call defrr(owk, nri)
        call dpcopy(rhoi,w(owk),1,nri,1d0)
        call dpadd(rhoi,rhoi(1+nri),1,nri,1d0)
      endif
      if(ipr >= 50) write(6,221)
c -------- start loop over atoms ---------------------
      do 10 ib=1,nbas
      is=ips(ib)
      nlml=(lmxl(is)+1)**2
      r1=rmt(is)
      i1=ioffv0(ib)+1
      call defrr(odum,   10)
      call strxsu(nlml,nxi,lxi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      call defrr(os,    nlml*nlmbx)
c -------- start loop over other atoms ---------------
      do 20 jb=1,nbas
      if(jb == ib .and. ixp(1) == 0) goto 20
      js=ips(jb)
      call dpdist(pos(1,ib),pos(1,jb),3,dd)
      if(alat*dd > rint(js)+r1) goto 20
      dr(1)=alat*(pos(1,jb)-pos(1,ib))
      dr(2)=alat*(pos(2,jb)-pos(2,ib))
      dr(3)=alat*(pos(3,jb)-pos(3,ib))
c -------- loop over energies, add to qmoi --------
      j1=ioff(jb)+1
      call defrr(ohl,nlmp*nxi(js))
      call rstr0(nxi(js),lxi(1,js),exi(1,js),nlmp,1,dr(1),dr(2),dr(3),
     .  lmxl(is),1,w(ohl),w(ohl))
      do 11 je=1,nxi(js)
      e=exi(je,js)
      lb=lxi(je,js)
      nlmb=(lb+1)**2
      ojj=ojkl(lb)
C|    call mstrux(e,dr,w(os),nlml,nlmb,nlml,cg,indxcg,jcg,cy)
C     call nstrux(e,dr,nlml,nlmb,nlmp,npow,w(oikl),w(ojj),w(oip),
C    .   w(ocf),cy,w(os))
C      ih = nlmp*(je-1)
C      call hstrux(e,nlml,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
C     .  w(ocf),w(ohl),w(os))
      call hstrux(e,nlml,nlmb,nlmp,npow,je,je,w(oikl),w(ojj),w(oip),
     .  w(ocf),w(ohl),w(os))
      call vesmi1(e,r1,nlml,nlmb,w(os),rhoi(j1),qmoi(i1))
  11  j1=j1+nlmb
      call rlse(ohl)
  20  continue
      call rlse(odum)
  10  continue
      if (nsp == 2) call dpcopy(w(owk),rhoi,1,nri,1d0)
      if (nsp == 2) call rlse(owk)
      call tcx('vesmoz')

c -------- printout -------------
      if(ipr < 50) return
      do 30 ib=1,nbas
      i0=ioffv0(ib)
      write(6,220) i0+1,ib,1,qmoi(i0+1),qmoi(i0+1)/y0
      do 31 ilm=2,nlml
      i=i0+ilm
  31  if(dabs(qmoi(i)) > 1d-6) write(6,220) i,ib,ilm,qmoi(i)
  220 format(i12,2i6,2f12.6,f8.2)
  221 format(/' vesmoi:   i    ib   ilm      qmom',8x,'q')
  30  continue
      end

c ---------- sub vesmi1 ------------------------
      subroutine vesmi1(e,r1,nlml,nlmb,s,rhoi,qmoi)
Cr Use integral_0_R (J_L Y_L r**l) d^3r = j_l+1 R^(l+2) ??
      implicit real*8 (a-h,p-z), integer (o)
      dimension qmoi(1),phi(0:20),psi(0:20),s(nlml,nlmb),rhoi(1)
      lmaxl=ll(nlml)
      call bessl(e*r1*r1,lmaxl+1,phi,psi)
      do 14 ilml=1,nlml
      l=ll(ilml)
      bmom=phi(l+1)*r1**(2*l+3)
      do 14 ilmb=1,nlmb
  14  qmoi(ilml)=qmoi(ilml)+bmom*s(ilml,ilmb)*rhoi(ilmb)
      end
