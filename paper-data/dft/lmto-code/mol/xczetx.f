      subroutine xczetx(nr,nr1,nr2,rofi,rwgt,nxi,lxi,exi,rsm,xi,xi0,
     .   nlml,vxc1,vxc2,dvxc2,zetxc,nri)
C  Subtracts head of xi_m times smooth xc pot from zetxc.
C  4 Jan 95 spin polarized (MvS)
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(1),lxi(1),rofi(1),rwgt(1),zetxc(nri,1),
     .  vxc1(nr1,nlml,1),vxc2(nr2,nlml,1),dvxc2(nr2,nlml,1),
     .  xi(nr,0:1),xi0(nr,0:1)
      call getpr(ipr)
      lmaxl=ll(nlml)
      nsp=lsp()+1
C ------ loop over Hankel energies ------------
      do 10 isp=1,nsp
      i0=0
      do 10 ie=1,nxi
      lx=min0(lxi(ie),lmaxl)
      nlmx=(lx+1)**2
      call rophs0(exi(ie),rsm,lx,nr,rofi,xi,001)
      call rophs0(exi(ie),0d0,lx,nr,rofi,xi0,001)
      do 11 ilm=1,nlmx
      l=ll(ilm)
      sum=0d0
C ... sum: smoothed xi * smoothed vxc inside rmt
C     sam: smoothed xi * (smoothed - true vxc) outside rmt
C     sim: (true - smoothed) xi * true vxc outside rmt
      do 6 ir=1,nr1
  6   sum=sum + vxc1(ir,ilm,isp)*xi(ir,l) * rwgt(ir)*rofi(ir)**2
      sam=0d0
      sim=0d0
      do 7 ir=1,nr2
      jr=ir+nr1
      ww=rwgt(jr)*rofi(jr)**2
      sam=sam + ww*dvxc2(ir,ilm,isp) * xi(jr,l)
  7   sim=sim + ww*vxc2(ir,ilm,isp)  * (xi0(jr,l)-xi(jr,l))
  11  zetxc(ilm+i0,isp)=zetxc(ilm+i0,isp)-sum+sam+sim
C|11  zetxc(ilm+i0,isp)=zetxc(ilm+i0,isp)-sum
  10  i0=i0+(lxi(ie)+1)**2
      end
