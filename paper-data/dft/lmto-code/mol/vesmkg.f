      subroutine vesmkg(nxi,lxi,exi,rmt,rsm,lmxl,n0,nri,nvi,rhoi,rho0,
     .   nbas,ips,qmom,poti,pot0,pot00)
c  Solves poisson equation and makes interstitial e-static potential.
c  The density is the smooth density as given by rhoi,rho0 plus point
c  multipoles to correct the moments inside the spheres.
c  Output: poti  = coeffs to smooth Hankels, e /= 0.
c          pot0  = coeffs to smooth Hankels, e == 0.
c          pot00 = coeffs to normal Hankels, e == 0.
      implicit real*8 (a-h,p-z), integer (o)
      dimension rmt(1),nxi(1),lxi(n0,1),exi(n0,1),qmom(1),lmxl(1),
     .   rsm(1),poti(nri),pot0(nvi),pot00(nvi),qg(0:20),qh(0:20),
     .   rhoi(nri),rho0(1),ips(1),iii(10),pout(10)
      call tcn('vesmkg')
      call getpr(ipr)
      call dpscop(rho0,pot0,nvi,1,1,2d0)
      call dpzero(pot00,nvi)
      pi=4d0*datan(1d0)
      atepi=8d0*pi
c -------- start loop over atoms ----------------
      i=0
      j0=0
      do 10 ib=1,nbas
      is=ips(ib)
      r=rmt(is)
      lmax=0
c ---------- start loop over xi-energies -----------------
      do 20 ie=1,nxi(is)
        e=exi(ie,is)
        lx=lxi(ie,is)
        lmaxl=lmxl(is)
        lmax=max(lmax,lx)
        call hsmmom(r,rsm(is),e,lx,qg,qh)
        beta=dexp(e*rsm(is)**2/4d0)
c ---------- start loop over ilm -------------------------
        ilm=0
        do 20 l=0,lx
        do 20 m=-l,l
        ilm=ilm+1
        i=i+1
        j=j0+ilm
        poti(i)=(atepi/e)*rhoi(i)
        pot0(j)=pot0(j)-(atepi/e)*beta*rhoi(i)
   20   if(l <= lmaxl) pot00(j)=pot00(j)+rhoi(i)*qh(l)
c ---------- contribution from gaussians; e arbitrary ------------
        call hsmmom(r,rsm(is),e,lmaxl,qg,qh)
        j=j0
        do 21 l=0,lmaxl
        do 21 m=-l,l
        j=j+1
   21   pot00(j) = pot00(j) + rho0(j)*qg(l)
c ---------- contribution from point multipoles ----------
        ilm=0
        df=1d0
        do 25 l=0,lmaxl
        df=df*(2*l+1)
        do 25 m=-l,l
        ilm=ilm+1
        j=j0+ilm
   25   pot00(j)=(qmom(j)-pot00(j))*atepi/df
c -------- end loop over atoms --------------
  10  j0=j0+(lmax+1)**2

c -------- printout ---
      if(ipr < 50) return
      write(6,220)
  220 format(/' vesmak:'/'   ib  ilm',5x,'V00+V0',6x,'V0',8x,'V(EI)')
      i0=0
      j0=0
      do 30 ib=1,nbas
      is=ips(ib)
      lmax=0
      jjj=0
      do 31 ie=1,nxi(is)
      iii(ie)=jjj
      lx=lxi(ie,is)
      lmax=max(lmax,lx)
   31 jjj=jjj+(lx+1)**2
      nlm=(lmax+1)**2
      do 32 m=1,nxi(is)
   32 pout(m)=0
      do 33 ilm=1,nlm
        pmax=0
        do 34 ie=1,nxi(is)
        nlm1=(lxi(ie,is)+1)**2
        if (ilm <= nlm1) pout(ie)=poti(iii(ie)+ilm+i0)
   34   pmax=max(pmax,dabs(pout(ie)))
        j=j0+ilm
        if (pmax >= 1d-5.or.ilm == 1) write(6,335) ib,ilm,
     .    pot00(j)+pot0(j),pot0(j),(pout(m),m=1,nxi(is))
  335   format(2i5,2f11.4,5f11.4)
   33 continue
      j0=j0+nlm
   30 i0=i0+jjj

      call tcx('vesmkg')
      end
