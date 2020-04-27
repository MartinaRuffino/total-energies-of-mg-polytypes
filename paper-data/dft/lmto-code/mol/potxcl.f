      subroutine potxcl(nlml,np,nr1,rofi,rwgt,yl,wp,rl,
     .   rp,exc,vxc,qsp,rep,rmu,vxcl,excl)
c  Makes yl-expansion of vxcl, given yl-expansion of rho in rl,
c  by integrating over point mesh.
      implicit real*8 (a-h,p-z), integer (o)
      dimension rofi(1),rwgt(1),yl(nlml,np),wp(np),rl(nr1,nlml),
     .   rp(nr1,np),exc(nr1,np),vxc(nr1,np),vxcl(nr1,1),sum(0:20),
     .   excl(nr1,1)
      call getpr(ipr)
      srfpi=dsqrt(16d0*datan(1d0))
      qsp=0d0
      do 1 ir=2,nr1
  1   qsp=qsp+rwgt(ir)*rl(ir,1)
      qsp=qsp*srfpi
c ------- evaluate density on point-wise mesh through sphere ---
      call dpzero(rp,   nr1*np)
      do 10 ip=1,np
      do 10 ilm=1,nlml
      do 10 ir=2,nr1
  10  rp(ir,ip)=rp(ir,ip)+rl(ir,ilm)*yl(ilm,ip)/rofi(ir)**2

c ------- make eps and mu on mesh, integrals --------
      call ropevx(rp,exc,vxc,nr1*np)
      rep=0d0
      rmu=0d0
      do 12 ip=1,np
      do 12 ir=1,nr1
      weight=rofi(ir)**2*rwgt(ir)*wp(ip)
      rep=rep+exc(ir,ip)*rp(ir,ip)*weight
  12  rmu=rmu+vxc(ir,ip)*rp(ir,ip)*weight

c ------- put yl-projection of vxc into vxcl -------
      call dpzero(excl,nr1*nlml)
      call dpzero(vxcl,nr1*nlml)
      do 20 ilm=1,nlml
      do 20 ip=1,np
      do 20 ir=2,nr1
      vxcl(ir,ilm)=vxcl(ir,ilm)+vxc(ir,ip)*wp(ip)*yl(ilm,ip)
  20  excl(ir,ilm)=excl(ir,ilm)+exc(ir,ip)*wp(ip)*yl(ilm,ip)
      end
