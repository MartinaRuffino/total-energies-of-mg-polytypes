      subroutine xcpmsx(nlml,np,nr2,rofi,rwgt,yl,wp,rl,
     .   rp,exc,vxc,vxcl,dvxcl,drl,drp,q2,rep,rmu)
c  Makes yl-expansion of vxcl and drl * deriv of mu
      implicit real*8 (a-h,p-z), integer (o)
      dimension rofi(1),rwgt(1),yl(nlml,np),wp(np),rl(nr2,nlml),
     .   rp(nr2,np),exc(nr2,np),vxc(nr2,np),vxcl(nr2,nlml),
     .   dvxcl(nr2,nlml),drl(nr2,nlml),drp(nr2,np)
      call getpr(ipr)
      drho=1d-6
      srfpi=dsqrt(16d0*datan(1d0))

c ------- evaluate density,exc,vxc on point-wise mesh --------
      call dpzero(rp,   nr2*np)
      call dpzero(drp,  nr2*np)
      do 10 ip=1,np
      do 10 ilm=1,nlml
      do 10 ir=1,nr2
      rp(ir,ip)=rp(ir,ip)+rl(ir,ilm)*yl(ilm,ip)
  10  drp(ir,ip)=drp(ir,ip)+drl(ir,ilm)*yl(ilm,ip)
      call ropevx(rp,exc,vxc,nr2*np)

c ------- put yl-projection of vxc into vxcl -------
      call dpzero(vxcl,nr2*nlml)
      do 20 ilm=1,nlml
      do 20 ip=1,np
      do 20 ir=1,nr2
  20  vxcl(ir,ilm)=vxcl(ir,ilm)+vxc(ir,ip)*wp(ip)*yl(ilm,ip)

c ------- numerical derivative of vxc ----------------
      call dpzero(dvxcl,nr2*nlml)
      xxx=0.5d0/drho
      do 30 ip=1,np
      do 30 ir=1,nr2
  30  rp(ir,ip)=rp(ir,ip)+drho
      call ropevx(rp,exc,vxc,nr2*np)
      do 31 ilm=1,nlml
      do 31 ip=1,np
      do 31 ir=1,nr2
      dvxcl(ir,ilm)=dvxcl(ir,ilm)
     .    +vxc(ir,ip)*drp(ir,ip)*wp(ip)*yl(ilm,ip)*xxx
  31  continue
      do 32 ip=1,np
      do 32 ir=1,nr2
  32  rp(ir,ip)=rp(ir,ip)-2d0*drho
      call ropevx(rp,exc,vxc,nr2*np)
      do 33 ilm=1,nlml
      do 33 ip=1,np
      do 33 ir=1,nr2
      dvxcl(ir,ilm)=dvxcl(ir,ilm)
     .    -vxc(ir,ip)*drp(ir,ip)*wp(ip)*yl(ilm,ip)*xxx
  33  continue

c ------- integrals -----------
      q2=0d0
      do 23 ir=1,nr2
 23   q2=q2+rofi(ir)**2*rwgt(ir)*drl(ir,1)
      rep=0d0
      rmu=0d0
      do 28 ilm=1,nlml
      do 28 ir=1,nr2
      ww=rofi(ir)**2*rwgt(ir)
      rep=rep+ww*drl(ir,ilm)*vxcl(ir,ilm)
  28  rmu=rmu+ww*rl(ir,ilm)*dvxcl(ir,ilm)
      rmu=rmu+rep
      q2=q2*srfpi

      end
