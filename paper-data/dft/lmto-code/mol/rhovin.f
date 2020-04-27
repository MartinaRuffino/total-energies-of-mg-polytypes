      subroutine rhovin(nri,nvi,nbas,ips,lmxa,njj,nv,ioffp,
     .  rhout,rhoux,rhotg,rhoxg,zeta,zets,zet0,zet0s,qjj,vab0,
     .  pjj,quu,qus,qss,vab,puu,pus,pss,zetai,zet0i,sumev,rhov,tv)
C- Makes rhout * vin
      implicit none
      integer nri,nvi,nbas,njj,nv,ioffp(1),ips(1),lmxa(1)
      double precision quu(njj),qus(njj),qss(njj),rhov,qjj0(1),
     .  puu(njj),pus(njj),pss(njj),qjj(njj,1),pjj(njj,1),vab0(9,nbas,1),
     .  rhout(1),rhoux(nri,1),rhotg(nvi),rhoxg(nvi,1),vab(4,9,1),
     .  zeta(1),zets(nri,1),zet0(1),zet0s(nvi,1),zetai(1),zet0i(1)
      integer i,ib,ilen,ibp,iq,ilm,jlm,nlma,lp1,ll,is,iv,iprint
      double precision rzeta,rzet0,rvaug,rvaugn,
     .  rhtv,rhtvns,rxiv,rgv,xx,sumev,tv

C --- rho*v for the interstitial ---
      do  12  i = 1, nri
   12 zetai(i) = zeta(i)
      do  14  i = 1, nvi
   14 zet0i(i) = zet0(i)
      do  10  ib = 1, nbas
        do  16  i = 1, nri
   16   zetai(i) = zetai(i) - zets(i,ib)
        do  18  i = 1, nvi
   18   zet0i(i) = zet0i(i) - zet0s(i,ib)
   10 continue
      call dpdot(rhout,zetai,nri,rzeta)
      call dpdot(rhotg,zet0i,nvi,rzet0)

C --- rho*v for the augmentation and tail densities (l=0 pot) ---
      iq = 0
      rvaug = 0
      rhtv = 0
      do  20  ib = 1, nbas
        is = ips(ib)
        nlma = (lmxa(is)+1)**2
        do  24  jlm = 1, nlma
        do  24  ilm = 1, nlma
        iq = iq+1
        if (ilm == jlm) then
          lp1 = ll(ilm)+1
          rvaug = rvaug + quu(iq)*vab(1,lp1,ib) +
     .      qus(iq)*vab(2,lp1,ib) + qss(iq)*vab(4,lp1,ib)
          do   26  iv = 1, nv
   26     rhtv = rhtv + qjj(iq,iv)*vab0(lp1,ib,iv)
        endif
   24 continue
   20 continue

C --- rho*v for the augmentation and xi-tail densities (nonsph pot) ---
C ... augmentation density * nonsph pot for all spheres
      call dpdot(quu,puu,njj,xx)
      rvaugn = xx
      call dpdot(qus,pus,njj,xx)
      rvaugn = rvaugn + xx
      call dpdot(qss,pss,njj,xx)
      rvaugn = rvaugn + xx
      rvaug = rvaug + rvaugn
C ... tail density * nonsph pot for all spheres
      call dpdot(qjj,pjj,njj*nv,rhtvns)
C ... xi-tail density * smoothed pot for all spheres
      call dpdot(rhoux,zets,nri*nbas,rxiv)
      call dpdot(rhoxg,zet0s,nvi*nbas,rgv)
      rxiv = rxiv+rgv

C --- Kinetic energy and printout ---
      rhov = rzeta+rzet0 + rvaug + rxiv+rhtv+rhtvns
      tv = sumev - rhov
      if (iprint() <= 10) return
      print 333,
     .  rvaug,  rzeta+rzet0, rxiv, rhtv+rhtvns, rxiv+rhtv+rhtvns,
     .  rvaugn, rxiv+rhtv, rhtvns
  333 format(/'rhovin:',5x,
     .'rho-aug*v      rhoi*v     xi-tail*v    -qtail*v     diff'/
     .' total',f15.6,f14.6,3f12.6/'  nsph',f15.6,14x,3f12.6)
      print 334, rhov, sumev, tv
  334 format(' rho*v',f15.6,'   sumev=',f14.6,'   ekin=',f15.6)

      end
