      subroutine zthead(nrad,rad,nlml,nxi,lxi,exi,rmt,rsm,xi,np,
     .   ves,yl,vhal,zets,sxi0)
c  On-site head terms of zets (integrals of xi_m times i-pot).
c  ves must be: potential on mesh with weights multiplied in.
c  On output, vhal is l-decomposed pot times weights and rofi**2.
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(1),lxi(1),vhal(nrad,nlml),zets(1),sxi0(1),rad(1),
     .   xi(nrad,0:1),yl(nlml,np),ves(nrad,np),xi1(0:2),xi0(0:2)
      lmaxl=ll(nlml)
c ------- in work array vhal, make yl's times ves -----------
      call dpzero(vhal,   nrad*nlml)
      do 3 ir=1,nrad
      do 3 ilm=1,nlml
      do 3 ip=1,np
  3   vhal(ir,ilm)=vhal(ir,ilm)+ves(ir,ip)*yl(ilm,ip)
C ------ start loop over Hankels ---------------
      i0=0
      do 10 ie=1,nxi
      lx=min0(lxi(ie),lmaxl)
      nlmx=(lx+1)**2
      call rophs0(exi(ie),rsm,lx,nrad,rad,xi,001)
      do 11 ilm=1,nlmx
      l=ll(ilm)
      do 6 ir=1,nrad
  6   zets(ilm+i0)=zets(ilm+i0)+xi(ir,l)*vhal(ir,ilm)
  11  continue
  10  i0=i0+(lxi(ie)+1)**2
C ------ this part adds to sxi0 ------------------
      srfpi=dsqrt(16d0*datan(1d0))
      asm=1d0/rsm
      gam=1d0/(4d0*asm*asm)
      i1=1
      do 14 ie=1,nxi
      nlm=(lxi(ie)+1)**2
      call hansmr(rmt,exi(ie),asm,xi1,1)
      call hansmr(rmt,0d0,asm,xi0,1)
      ee=dexp(gam*exi(ie))
      s00=srfpi*(-ee-rmt**3*(xi1(1)-ee*xi0(1)))/exi(ie)
      sxi0(i1)=sxi0(i1)+s00
  14  i1=i1+nlm

      end
