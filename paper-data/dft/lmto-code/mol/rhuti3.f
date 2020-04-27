      subroutine rhuti3(wgt,el,nel,nphi,lphi,nxi,lxi,exi,pos,r,
     .  n0,nlmx,nri,nvi,nbas,ips,alat,tspec,tdata,
     .  cg,jcg,indxcg,ioff,ioffv0,nhs,qval,t,rhoi,rhoix,rho0,rho0x)
c  Adds together output interstitial density in rhoi.
c  Tails of rhoix(..,ib) are to be added inside sphere ib.
      parameter( nrot=81, nfx=300 )
      implicit real*8 (a-h,p-z), integer (o)
      dimension nxi(1),lxi(n0,1),exi(n0,1),pos(3,1),dr(3),r(1),
     .  nphi(1),lphi(n0,1),cg(1),jcg(1),indxcg(1),
     .  tspec(1),tdata(1),ips(1),t(nhs,nhs),el(1),ioff(1),rhoi(1),
     .  rhoix(nri,nbas),rho0(nvi),rho0x(nvi,nbas),ioffv0(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('rhuti3')
      nqval=qval+0.01
      nstate=nqval/2
      istate = 1
      jstate = nstate
      if (2*nstate /= nqval) call rx('odd electron number')

      call defrr(ormat,   nrot*nrot)
      call defrr(owk,     nfx)
      call defrr(owr,     nfx)
      call defrr(oro1,    nfx)
      call defrr(oro2,    nfx)
      call defrr(org1,    nfx)

c --- start loop over hamiltonian blocks -----
      iham=1
      do 10 ie=1,nel
      e1=el(ie)
      do 10 ib=1,nbas
      is=ips(ib)
      r1=r(is)
      jham=1
      do 12 je=1,nel
      e2=el(je)
      do 12 jb=1,nbas
      js=ips(jb)
      r2=r(js)
      do 14 i=1,3
   14 dr(i)=alat*(pos(i,jb)-pos(i,ib))
      l1=lphi(ie,is)
      l2=lphi(je,js)
      nx1=nxi(is)
      nx2=nxi(js)
      nlm1=(l1+1)**2
      nlm2=(l2+1)**2
      nfi=ioff(ib+1)-ioff(ib)
      nfj=ioff(jb+1)-ioff(jb)
      nfig=nfi+(lxi(nx1+1,is)+1)**2
      nfjg=nfj+(lxi(nx2+1,js)+1)**2
      i1=ioff(ib)+1
      j1=ioff(jb)+1
      if (l1 < 0 .or. l2 < 0) goto 12
      if(jham > iham) goto 12
      if(nfig > nfx.or.nfjg > nfx) call rx('rhuti3 increase nfx')
      wgt1=wgt
      if(iham /= jham) wgt1=2d0*wgt1
c --- case ib ne jb ---
      if (ib /= jb) then
        call defrr(oc1,     nfi*nlm1*nlm2)
        call defrr(oc2,     nfj*nlm1*nlm2)
        call hyfevl(dr,r1,e1,nlm1,r2,e2,nlm2,w(oc1),w(oc2),
     .    w(ormat),nrot,nfi,nfj,nxi(is),lxi(1,is),exi(1,is),
     .    nxi(js),lxi(1,js),exi(1,js),tspec,tdata)
        call hhrui2(nlm1,nlm2,nfi,nfj,w(oc1),w(oc2),nhs,t(iham,1),
     .    t(jham,1),istate,jstate,wgt1,w(owk),w(owr),
     .    w(oro1),w(oro2),w(ormat),nrot,nxi(is),lxi(1,is),
     .    nxi(js),lxi(1,js))
        call dpadd(rhoi(i1),w(oro1),1,nfi,1d0)
        call dpadd(rhoi(j1),w(oro2),1,nfj,1d0)
        do 20 kb=1,nbas
        if(kb /= ib.and.kb /= jb) then
          call dpadd(rhoix(i1,kb),w(oro1),1,nfi,1d0)
          call dpadd(rhoix(j1,kb),w(oro2),1,nfj,1d0)
        endif
  20    continue

c --- case ib eq jb ---
      else
        call defrr(oc1,     nfig*nlm1*nlm2)
        i1g=ioffv0(ib)+1
        mfi=nfig-nfi
        call oncget(r1,e1,nlm1,e2,nlm2,w(oc1),nfig,nxi(is)+1,
     .    lxi(1,is),exi(1,is),tspec,tdata,cg,jcg,indxcg)
        call hhrui1(nlm1,nlm2,nfi,nfig,w(oc1),nhs,t(iham,1),t(jham,1),
     .    istate,jstate,wgt1,w(owk),w(oro1),w(org1))
        call dpadd(rhoi(i1),w(oro1),1,nfi,1d0)
        call dpadd(rho0(i1g),w(org1),1,mfi,1d0)
        do 22 kb=1,nbas
        if(kb /= ib) call dpadd(rhoix(i1,kb),w(oro1),1,nfi,1d0)
  22    if(kb /= ib) call dpadd(rho0x(i1g,kb),w(org1),1,mfi,1d0)
      endif
   12 jham=jham+nlm2
   10 iham=iham+nlm1
      call rlse(ormat)
      call tcx('rhuti3')
      end
c -------- sub hhrui1 ---------------
      subroutine hhrui1(nlm1,nlm2,nfi,nfig,c1,nhs,t1,t2,istate,jstate,
     .  wgt1,wk,rhoi1,rho01)
      implicit real*8 (a-h,p-z), integer (o)
      dimension c1(nfig,nlm1,nlm2),rho01(1),
     .  t1(nhs,1),t2(nhs,1),rhoi1(1),wk(1)
      call dpzero(rhoi1,nfi)
      call dpzero(rho01,nfig-nfi)
      do 30 is=istate,jstate
      call dpzero(wk,nfig)
      do 10  i1=1,nlm1
      do 10  i2=1,nlm2
      xx = wgt1*t1(i1,is)*t2(i2,is)
      do 10 m=1,nfig
   10 wk(m) = wk(m) + xx*c1(m,i1,i2)
      do 11 m=1,nfi
   11 rhoi1(m)=rhoi1(m) + wk(m)
      do 12 m=1,nfig-nfi
   12 rho01(m) = rho01(m) + wk(m+nfi)
   30 continue
      end
c -------- sub hhrui2 ---------------
      subroutine hhrui2(nlm1,nlm2,nfi,nfj,c1,c2,nhs,t1,t2,
     .  istate,jstate,wgt1,wk,wr,rhoi1,rhoi2,
     .  rmat,nrot,nx1,lx1,nx2,lx2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension c1(nfi,nlm1,nlm2),c2(nfj,nlm2,nlm1),
     .  t1(nhs,1),t2(nhs,1),rhoi1(1),rhoi2(1),wk(1),
     .  b1(50),b2(50),wr(1),lx1(1),lx2(1),rmat(nrot,nrot)
      call dpzero(rhoi1,nfi)
      call dpzero(rhoi2,nfj)
c  start loop over states
      do 30 is=istate,jstate
      call prerot(t1(1,is),b1,nlm1,rmat,nrot)
      call prerot(t2(1,is),b2,nlm2,rmat,nrot)
c  add together unrotated density on site 1
      call dpzero(wk,nfi)
      do 10  i1=1,nlm1
      do 10  i2=1,nlm2
      xx = wgt1*b1(i1)*b2(i2)
      do 10 m=1,nfi
   10 wk(m)=wk(m) + xx*c1(m,i1,i2)
c  rotate density
      i1=1
      do 15 ie=1,nx1
      nlm=(lx1(ie)+1)**2
      call posrot(wk(i1),wr(i1),nlm,rmat,nrot)
  15  i1=i1+nlm
c  add into rhoi1
      do 11 m=1,nfi
   11 rhoi1(m)=rhoi1(m) + wr(m)
c  add together unrotated density on site 2
      call dpzero(wk,nfj)
      do 20  i1=1,nlm1
      do 20  i2=1,nlm2
      xx = wgt1*b1(i1)*b2(i2)
      do 20 m=1,nfj
   20 wk(m)=wk(m) + xx*c2(m,i2,i1)
c  rotate density
      i1=1
      do 16 ie=1,nx2
      nlm=(lx2(ie)+1)**2
      call posrot(wk(i1),wr(i1),nlm,rmat,nrot)
  16  i1=i1+nlm
c  add into rhoi2
      do 22 m=1,nfj
   22 rhoi2(m)=rhoi2(m) + wr(m)
   30 continue

      end
