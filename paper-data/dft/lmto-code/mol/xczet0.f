      subroutine xczet0(nr1,rofi,rwgt,nxi,lxi,exi,rsm,xi,
     .   nlml,vxcl,zetxc)
C  Subtracts head of xi_m times smooth xc pot from zetxc.
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(1),lxi(1),xi(nr1,0:1),rofi(1),rwgt(1),
     .   vxcl(nr1,1),zetxc(1)
      call getpr(ipr)
      lmaxl=ll(nlml)
C ------ loop over smooth Hankel energies ------------
      i0=0
      do 10 ie=1,nxi
      lx=min0(lxi(ie),lmaxl)
      nlmx=(lx+1)**2
      call rophs0(exi(ie),rsm,lx,nr1,rofi,xi,001)
      do 11 ilm=1,nlmx
      l=ll(ilm)
      sum=0d0
      do 6 ir=1,nr1
  6   sum=sum+vxcl(ir,ilm)*xi(ir,l) * rwgt(ir)*rofi(ir)**2
  11  zetxc(ilm+i0)=zetxc(ilm+i0)-sum
  10  i0=i0+(lxi(ie)+1)**2
      end
