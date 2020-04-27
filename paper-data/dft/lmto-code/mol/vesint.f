      subroutine vesint(nxi,lxi,exi,rmt,lmxl,n0,nri,nvi,rhoi,
     .   nbas,ips,qmom,qmoi,poti,pot0)
c  Solves poisson equation and makes interstitial e-static potential.
c  The density is the smooth density as given by rhoi, plus point
c  multipoles to correct the moments inside the spheres.
c  Output: poti  = coeffs to Hankels, e /= 0.
c          pot0  = coeffs to Hankels, e == 0.
c  28 Dec 94 made spin pol (MvS)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rmt(1),nxi(1),lxi(n0,1),exi(n0,1),qmom(1),lmxl(1),
     .   poti(nri),pot0(nvi),phi(0:20),psi(0:20),
     .   rhoi(nri),ips(1),iii(10),pout(10),qmoi(1)
      real w(1)
      common /w/ w
      call tcn('vesint')
      call getpr(ipr)
      call dpzero(pot0,nvi)
      pi=4d0*datan(1d0)
      atepi=8d0*pi
C -------- copy rho+ to work and add rho- into rho+ (nsp=2) -----------
      nsp = lsp()+1
      if (nsp == 2) then
        call defrr(owk, nri)
        call dpcopy(rhoi,w(owk),1,nri,1d0)
        call dpadd(rhoi,rhoi(1+nri),1,nri,1d0)
      endif

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
        call bessl(e*r*r,lx+1,phi,psi)
c ---------- start loop over ilm -------------------------
        ilm=0
        df=1d0
        do 20 l=0,lx
        df=df*(2*l+1)
        qh=(psi(l+1)-df)/e
        do 20 m=-l,l
        ilm=ilm+1
        i=i+1
        j=j0+ilm
        poti(i)=(atepi/e)*rhoi(i)
        pot0(j)=pot0(j)-(atepi/e)*rhoi(i)
  20    pot0(j)=pot0(j)-rhoi(i)*qh*atepi/df
c ---------- add to pot0 from tail and true multipoles ------
      ilm=0
      df=1d0
      do 21 l=0,lmax
      df=df*(2*l+1)
      do 21 m=-l,l
      ilm=ilm+1
      j=j0+ilm
  21  pot0(j)=pot0(j)+(qmom(j)-qmoi(j))*atepi/df
c -------- end loop over atoms --------------
  10  j0=j0+(lmax+1)**2

C -------- restore rho+, copy poti into poti- -----------
      if (nsp == 2) then
        call dpcopy(w(owk),rhoi,1,nri,1d0)
        call rlse(owk)
        call dpcopy(poti,poti(1+nri),1,nri,1d0)
        call dpcopy(pot0,pot0(1+nvi),1,nvi,1d0)
      endif
      call tcx('vesint')

c -------- printout ---
      if(ipr < 50) return
      write(6,220)
  220 format(/' vesint:'/'   ib  ilm',7x,'V0',7x,'V(EI)')
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
     .     pot0(j),(pout(m),m=1,nxi(is))
  335   format(2i5,2f11.4,5f11.4)
   33 continue
      j0=j0+nlm
   30 i0=i0+jjj

      end
