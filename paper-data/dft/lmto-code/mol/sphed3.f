      subroutine sphed3(nrad,rad,wrad,lmaxl,nxi,lxi,exi,
     .   rsm,xi,rhoi,rho0,poti,pot0,pot00,rhal,vhal)
C  adds interstitial head density to l-decomposed density in rhal,
c  and interstitial head e-static potential to vhal.
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(1),lxi(1),rhoi(1),poti(1),pot0(1),pot00(1),rho0(1),
     .  xi(nrad,0:lmaxl),rad(1),wrad(1),rhal(nrad,1),vhal(nrad,1)
      call getpr(ipr)
      nlml=(lmaxl+1)**2
C ------ contribution to rho,v from smooth Hankels, e /= 0 ------
      i0=0
      do 10 ie=1,nxi
      lx=min0(lxi(ie),lmaxl)
      nlmx=(lx+1)**2
      call rophs0(exi(ie),rsm,lx,nrad,rad,xi,001)
      do 11 ilm=1,nlmx
      l=ll(ilm)
      do 6 ir=1,nrad
      rhal(ir,ilm)=rhal(ir,ilm)+rhoi(ilm+i0)*xi(ir,l)
    6 vhal(ir,ilm)=vhal(ir,ilm)+poti(ilm+i0)*xi(ir,l)
   11 continue
   10 i0=i0+(lxi(ie)+1)**2
c ------- contribution to rho from gaussians --------
      call ropgau(rsm,lmaxl,nrad,rad,xi,001)
      do 9 ilm=1,nlml
      l=ll(ilm)
      do 9 ir=1,nrad
    9 rhal(ir,ilm)=rhal(ir,ilm)+rho0(ilm)*xi(ir,l)
c ------- contribution to v from smooth hankels with e=0 --------
      call rophs0(0d0,rsm,lmaxl,nrad,rad,xi,001)
      do 8 ilm=1,nlml
      l=ll(ilm)
      do 8 ir=1,nrad
    8 vhal(ir,ilm)=vhal(ir,ilm)+pot0(ilm)*xi(ir,l)
c ------- contribution to v from normal hankels with e=0 --------
      do 22 i=1,nrad
   22 xi(i,0)=1d0/rad(i)
      do 23 l=1,lmaxl
      do 23 i=1,nrad
   23 xi(i,l)=(2*l-1d0)*xi(i,l-1)/rad(i)
      do 24 ilm=1,nlml
      l=ll(ilm)
      do 24 ir=1,nrad
   24 vhal(ir,ilm)=vhal(ir,ilm)+pot00(ilm)*xi(ir,l)

c -------- printout -----------------
      if(ipr < 40) return
      qtot=0d0
      vtot=0d0
      do 42 ir=1,nrad
      qtot=qtot+rad(ir)**2*wrad(ir)*rhal(ir,1)
   42 vtot=vtot+rad(ir)**2*wrad(ir)*vhal(ir,1)
      srfpi=dsqrt(16d0*datan(1d0))
      write(6,992) qtot*srfpi,vtot*srfpi
  992 format(' integrals over y-xpn:  qtot=',f11.6,'    vtot=',f11.6)
      if(ipr <= 50) return
      write(6,991)
      do 40 j=1,nlml
      sum=0d0
      do 41 ir=1,nrad
   41 sum=sum+rad(ir)**2*wrad(ir)*rhal(ir,j)*vhal(ir,j)
   40 if(dabs(sum) > 1d-5) write(6,990) j,sum
  990 format(9x,i5,3x,f14.6)
  991 format(/' sphead:    ilm',5x,'int(vhal*rhal)')
      end
