      subroutine ppot3(nri,nxi,lxi,nvi,n0,nbas,ips,zeta,zets,rhoi,
     .  zeta0,zet0s,rho0,modep)
C- Prints out Hartree and exchange interstitial potential.
C  For printout, hankels renormalized by dividing by (2l+1)!!
C  For now, take out renormalizing ...
      implicit real*8 (a-h,p-z), integer (o)
      dimension zeta(nri,1),zets(nri,nbas,1),ips(1),nxi(1),lxi(n0,1),
     .  rhoi(nri,2),rho0(nvi,1),zeta0(nvi,1),zet0s(nvi,nbas,1)

      if (iprint() < 40) return
      nsp = lsp()+1
      print 334
  334 format(/'    n ib ie ilm    Rhoi      Zeta       Zets ...')
      do 12 isp=1,nsp
      if (nsp == 2) print *, 'spin 2:'
      i=0
      iv=0
      do 10 ib=1,nbas
      is=ips(ib)
      lmax = 0
C --- Printout for rhoi and associated quantities ---
      do 20 ie=1,nxi(is)
      lmax = max(lmax,lxi(ie,is))
      df=1d0
      ilm=0
      do 20 l=0,lxi(ie,is)
      df=df*(2*l+1)
      f=1d0/df
      f = 1
      do 20 m=1,2*l+1
      i=i+1
      ilm=ilm+1
      if (modep == 0 .and. dabs(zeta(i,isp)) > 1d-5)
     .  print 333,i,ib,ie,ilm,rhoi(i,isp),f*zeta(i,isp)
      if (modep == 1 .and. dabs(zeta(i,isp)) > 1d-5) print 333,i,ib,
     .  ie,ilm,rhoi(i,isp),f*zeta(i,isp),(f*zets(i,jb,isp),jb=1,nbas)
  333 format(i5,3i3,6f11.6:/100(11x,'...',4f11.6:/))
   20 continue
C --- Printout for gaussians and associated quantities ---
      ilm=0
      do 30 l=0,lmax
      do 30 m=1,2*l+1
      iv=iv+1
      ilm=ilm+1
      if (dabs(zeta0(iv,isp)) > 1d-5) print 332,iv,ib,ilm,
C     if (dabs(zeta0(iv,isp))+dabs(zet0s(iv,1,isp))+
C     .  dabs(zet0s(iv,2,isp)) > 1d-5) print 332,iv,ib,ilm,
     .  rho0(iv,isp),f*zeta0(iv,isp),(f*zet0s(iv,jb,isp), jb=1, nbas)
  332 format(i5,i3,'  G',i3,6f11.6:/100(11x,'...',4f11.6:/))
   30 continue

   10 continue
   12 continue
      end
