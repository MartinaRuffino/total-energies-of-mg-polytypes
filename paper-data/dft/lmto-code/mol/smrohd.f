      subroutine smrohd(nr1,rofi,rwgt,nxi,lxi,exi,rsm,xi,
     .   rhoi,nlml,rl)
C  adds interstitial head density to l-decomposed density in rl
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(1),lxi(1),rhoi(1),xi(nr1,0:1),rofi(1),rwgt(1),
     .   rl(nr1,1)
      call getpr(ipr)
      lmaxl=ll(nlml)
      srfpi=dsqrt(16d0*datan(1d0))
      q1=0d0
      do 41 ir=1,nr1
  41  q1=q1+rofi(ir)**2*rwgt(ir)*rl(ir,1)
C ------ sum over smooth Hankel energies ------------
      i0=0
      do 10 ie=1,nxi
      lx=min0(lxi(ie),lmaxl)
      nlmx=(lx+1)**2
      call rophs0(exi(ie),rsm,lx,nr1,rofi,xi,001)
      do 11 ilm=1,nlmx
      l=ll(ilm)
      do 6 ir=1,nr1
  6   rl(ir,ilm)=rl(ir,ilm)+rhoi(ilm+i0)*xi(ir,l)
  11  continue
  10  i0=i0+(lxi(ie)+1)**2
      q2=0d0
      do 42 ir=1,nr1
  42  q2=q2+rofi(ir)**2*rwgt(ir)*rl(ir,1)
      if(ipr >= 40) write(6,992) q1*srfpi,q2*srfpi
  992 format(' smooth charge in sphere:     tails=',f11.6,
     .   '     tot=',f11.6)
      end
