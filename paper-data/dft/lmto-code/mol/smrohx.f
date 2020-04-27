      subroutine smrohx(nr,nr1,nr2,rofi,rwgt,nxi,lxi,exi,rsm,xi,
     .   rhoi,nri,nlml,nsp,rl1,rl2,drl2)
C- Adds interstitial head density to l-decomposed density in rl1,rl2.
C  Also makes drl2: diff of on-site unsmoothed and smoothed density.
C  30 Dec 94 made spin pol (MvS)
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(1),lxi(1),rhoi(nri,2),xi(nr,0:1),rofi(1),rwgt(1),
     .   rl1(nr1,nlml,nsp),drl2(nr2,nlml,nsp),rl2(nr2,nlml,nsp),
     .   q1(2),q2(2)
      call getpr(ipr)
      lmaxl=ll(nlml)
      srfpi=dsqrt(16d0*datan(1d0))
      do 41 isp = 1,nsp
      q1(isp)=0d0
      do 41 ir=1,nr1
  41  q1(isp)=q1(isp)+rofi(ir)**2*rwgt(ir)*rl1(ir,1,isp)

C ------ loop over Hankel energies ------------
      call dpzero(drl2,   nr2*nlml*nsp)
      i0=0
      do 10 ie=1,nxi
      lx=min0(lxi(ie),lmaxl)
      nlmx=(lx+1)**2
      call rophs0(exi(ie),rsm,lx,nr,rofi,xi,001)
      do 11 isp = 1,nsp
      do 11 ilm=1,nlmx
      l=ll(ilm)
      do 5 ir=1,nr1
  5   rl1(ir,ilm,isp)=rl1(ir,ilm,isp)+rhoi(ilm+i0,isp)*xi(ir,l)
      do 6 ir=1,nr2
  6   rl2(ir,ilm,isp)=rl2(ir,ilm,isp)+rhoi(ilm+i0,isp)*xi(ir+nr1,l)
      do 7 ir=1,nr2
  7   drl2(ir,ilm,isp)=drl2(ir,ilm,isp)-rhoi(ilm+i0,isp)*xi(ir+nr1,l)
  11  continue
      call rophs0(exi(ie),0d0,lx,nr,rofi,xi,001)
      do 14 isp = 1,nsp
      do 14 ilm=1,nlmx
      l=ll(ilm)
      do 8 ir=1,nr2
  8   drl2(ir,ilm,isp)=drl2(ir,ilm,isp)+rhoi(ilm+i0,isp)*xi(ir+nr1,l)
  14  continue
  10  i0=i0+(lxi(ie)+1)**2
c ------ integrals for info --------
      q3 = 0d0
      do 40 isp=1,nsp
      q2(isp)=0d0
      do 42 ir=1,nr1
   42 q2(isp)=q2(isp)+rofi(ir)**2*rwgt(ir)*rl1(ir,1,isp)
      do 43 ir=1,nr2
   43 q3=q3+rofi(ir+nr1)**2*rwgt(ir+nr1)*drl2(ir,1,isp)
      q1(isp)=q1(isp)*srfpi
      q2(isp)=q2(isp)*srfpi
      if(ipr > 40.and.nsp == 2) print 992, q1(isp),q2(isp),isp
  992 format(' qsm inside rmt:  tails',f10.6,
     .  '  tot',f10.6:' (spin ',i1,')')
      amom = 0
      if (isp == 2) then
        amom =q2(2)-q2(1)
        q1(1)=q1(1)+q1(2)
        q2(1)=q2(1)+q2(2)
      endif
   40 continue
      if(ipr >= 30) print 993, q1(1), q2(1), amom
  993 format(' qsm inside rmt:  tails',f10.6,
     .  '  tot',f10.6,'  smoothed mag mom=',f8.5)

      end
