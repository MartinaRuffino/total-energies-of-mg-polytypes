      subroutine prpot(nri,nxi,lxi,n0,nbas,ips,zeta,zetxc,rhoi)
C- Prints out Hartree and exchange interstitial potential.
C  For printout, hankels renormalized by dividing by (2l+1)!!
      implicit real*8 (a-h,p-z), integer (o)
      dimension zeta(nri,1),zetxc(nri,1),rhoi(nri,1),
     .  ips(1),nxi(1),lxi(n0,1)
      integer lgunit,stdo
      stdo = lgunit(1)
      nsp = lsp()+1
      do 1 isp=1,nsp
      call dpdot(rhoi(1,isp),zeta(1,isp),nri,sum1)
      call dpdot(rhoi(1,isp),zetxc(1,isp),nri,sum2)
    1 if(iprint() >= 30) write(stdo,331) sum1,sum2
  331 format(/' prpot:  dot products rhoi*zeta, rhoi*zetxc:',2f13.6)
      if (iprint() < 60) return
      write(stdo,334)
  334 format(/'    n ib ie ilm    Rhoi      Zeta',7x,'Zetxc',6x,'Sum')
      do 10 isp=1,nsp
      if (nsp == 2) write (stdo,'(''spin '',i1)') isp
      i=0
      do 10 ib=1,nbas
      is=ips(ib)
      do 10 ie=1,nxi(is)
      df=1d0
      ilm=0
      do 10 l=0,lxi(ie,is)
      df=df*(2*l+1)
      f=1d0/df
C|    f = 1
      do 10 m=1,2*l+1
      i=i+1
      ilm=ilm+1
      top=dmax1(dabs(rhoi(i,isp)),dabs(zeta(i,isp)),dabs(zetxc(i,isp)))
      if (top > 1d-5) then
        write (stdo,333) i,ib,ie,ilm,rhoi(i,isp)/f,
     .  f*zeta(i,isp),f*zetxc(i,isp),f*(zeta(i,isp)+zetxc(i,isp))
      endif
  333 format(i5,3i3,6f11.6)
   10 continue
      end
