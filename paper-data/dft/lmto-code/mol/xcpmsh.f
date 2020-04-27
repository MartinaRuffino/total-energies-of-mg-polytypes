      subroutine xcpmsh(nlml,np,nr1,rofi,rwgt,yl,wp,rl,
     .   rp,exc,vxc,qsm,rep,rmu,vxcl)
c  Makes yl-expansion of vxcl, given yl-expansion of rho in rl,
c  by integrating over point mesh.
      implicit real*8 (a-h,p-z), integer (o)
      dimension rofi(1),rwgt(1),yl(nlml,np),wp(np),rl(nr1,nlml),
     .   rp(nr1,np),exc(nr1,np),vxc(nr1,np),vxcl(nr1,nlml),sum(0:20)
      call getpr(ipr)
c ------- evaluate density on point-wise mesh through sphere ---
      call dpzero(rp,   nr1*np)
      do 10 ip=1,np
      do 10 ilm=1,nlml
      do 10 ir=1,nr1
  10  rp(ir,ip)=rp(ir,ip)+rl(ir,ilm)*yl(ilm,ip)

c ------- make eps and mu on mesh, integals --------
      call ropevx(rp,exc,vxc,nr1*np)
      qsm=0d0
      rep=0d0
      rmu=0d0
      do 12 ip=1,np
      do 12 ir=1,nr1
      weight=(rofi(ir)**2*rwgt(ir))*wp(ip)
      qsm=qsm+rp(ir,ip)*weight
      rep=rep+exc(ir,ip)*rp(ir,ip)*weight
  12  rmu=rmu+vxc(ir,ip)*rp(ir,ip)*weight
      if(ipr >= 30) write(6,725) qsm,rep,rmu
  725 format(' smooth rho:   qsm=',f11.6,'    rep=',f11.6,
     .   '    rmu=',f10.6)

c ------- put yl-projection of vxc into vxcl -------
      call dpzero(vxcl,nr1*nlml)
      do 20 ilm=1,nlml
      do 20 ip=1,np
      do 20 ir=1,nr1
  20  vxcl(ir,ilm)=vxcl(ir,ilm)+vxc(ir,ip)*wp(ip)*yl(ilm,ip)

c ------- print out rmu by angular momentum ------
      if(ipr < 30) return
      lmaxl=ll(nlml)
      do 33 l=0,lmaxl
  33  sum(l)=0d0
      do 30 ilm=1,nlml
      l=ll(ilm)
      do 30 ir=1,nr1
  30  sum(l)=sum(l)+ rl(ir,ilm)*vxcl(ir,ilm) * rofi(ir)**2*rwgt(ir)
      write(6,340) (sum(l),l=0,lmaxl)
  340 format(' rmu l-decomposed: ',4f12.6:/(19x,4f12.6))

      end
