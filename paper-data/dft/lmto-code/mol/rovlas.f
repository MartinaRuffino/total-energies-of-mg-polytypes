      subroutine rovlas(nxi,lxi,exi,n0,rmt,rint,a,nr,lmxl,
     .  nbas,alat,pos,ips,rhoi,ioff,cg,jcg,indxcg,ixp,dxp,orho,iov)
C  Overlap or decompose interstitial tail density for all spheres.
C  iov 1 for overlapping, -1 for decomposition
C  23 Dec 94 Mvs made spin pol.
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(n0,1),nxi(1),ips(1),lxi(n0,1),ioff(1),
     .  rhoi(1),pos(3,1),rint(1),cg(1),jcg(1),indxcg(1),
     .  rmt(1),lmxl(1),a(1),nr(1),orho(1),dr(3),ojkl(0:20),
     .  ixp(1),dxp(1),q1(2),qtot(2)
      integer intopt,nglob
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('rovlas')
      pi=4d0*datan(1d0)
      srfpi=dsqrt(4d0*pi)
      nri=ioff(nbas+1)
      nsp=lsp()+1
      if (ipr >= 35) print
     .  '('' rovlas: ib      q1          q2          qtot'')'
c ------ start loop over atoms ---------------------
      do 10 ib=1,nbas
      is=ips(ib)
      lmaxl=lmxl(is)
      nlml=(lmaxl+1)**2
      r1=rmt(is)
      a1=a(is)
      nr1=nr(is)
      oro=orho(ib)
      call defrr(orofi,       nr1)
      call defrr(orwgt,       nr1)
      call defrr(oxi,         nr1*(lmaxl+1))
      call defrr(oh,          nr1)
      call defrr(oy,          nr1)
      call radmsh(r1,a1,nr1,w(orofi))
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,r1,a1,nr1,w(orwgt))
C
C      print *, 'nr1 = ',nr1
C      call prmx('w(oro):',w(oro),nr1,nr1,1)
C      call prmx('w(orwgt):',w(orwgt),nr1,nr1,1)
C
      call dpdot(w(oro),w(orwgt),nr1,q1)
      if (nsp == 2) then
        call dpscop(w(oro),w(oh),nr(is),1+nr1*nlml,1,1d0)
        call dpdot(w(oh),w(orwgt),nr1,q1(2))
      endif
      call strxsu(nlml,nxi,lxi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      call defrr(os,          nlml*nlmbx)

c ------ start loop over neighbors to ib ---------------
      do 20 jb=1,nbas
      if(jb == ib .and. ixp(1) == 0) goto 20
      js=ips(jb)
      call dpdist(pos(1,ib),pos(1,jb),3,dd)
      if(alat*dd > rint(js)+r1) goto 20
      dr(1)=alat*(pos(1,jb)-pos(1,ib))
      dr(2)=alat*(pos(2,jb)-pos(2,ib))
      dr(3)=alat*(pos(3,jb)-pos(3,ib))
      j1=ioff(jb)+1
c ---------- loop over energies, add bessel tails ----
      call defrr(ohl,nlmp*nxi(js))
      call rstr0(nxi(js),lxi(1,js),exi(1,js),nlmp,1,dr(1),dr(2),dr(3),
     .  lmxl(is),1,w(ohl),w(ohl))
      do 11 je=1,nxi(js)
      e=exi(je,js)
      lb=lxi(je,js)
      nlmb=(lb+1)**2
      if(nlmb > nlmbx) call rx('rovlas: increase nlmbx')
      ojj=ojkl(lb)
      call hstrux(e,nlml,nlmb,nlmp,npow,je,je,w(oikl),w(ojj),w(oip),
     .  w(ocf),w(ohl),w(os))
      call ropbes(w(orofi),e,lmaxl,w(oy),w(oh),w(oxi),nr1,1)
      call rovlps(nlml,nlmb,w(os),nr1,rhoi(j1),w(oxi),w(oh),
     .  nri,nsp,w(orofi),w(oro),iov)
  11  j1=j1+nlmb
      call rlse(ohl)
  20  continue
      call dpdot(w(oro),w(orwgt),nr1,qtot)
      if (nsp == 2) then
        call dpscop(w(oro),w(oh),nr(is),1+nr1*nlml,1,1d0)
        call dpdot(w(oh),w(orwgt),nr1,qtot(2))
      endif
      if (ipr >= 35) then
        write(6,787) ib,srfpi*q1(1),srfpi*(qtot(1)-q1(1)),srfpi*qtot(1)
        if (nsp == 2) write(6,788) srfpi*q1(2),srfpi*(qtot(2)-q1(2)),
     .    srfpi*qtot(2),srfpi*(qtot(1)+qtot(2))
  787   format(7x,i4,3f12.6)
  788   format(7x,4x,3f12.6,'  sum=',f12.6)
      endif
      call rlse(orofi)
  10  continue
      call tcx('rovlas')
      end
      subroutine rovlps(nlml,nlmb,s,nr1,rhoi,xi,h,nri,nsp,rofi,rhol,iov)
C sign 1 for overlap, -1 for decompose
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhol(nr1,nlml,nsp),rofi(1),xi(nr1,0:1),s(nlml,nlmb),
     .   rhoi(nri,nsp),h(1)

C      do  iii = 1, nr1
C        write(*,100) rhol(iii,1,1)
C      enddo
C  100 format (6f12.6)
C      print *, '** In rovlps **'

      lmaxl=ll(nlml)
      do 1 i=1,nr1
    1 h(i)=rofi(i)
      ilml=0
      do 14 l=0,lmaxl
        do 2 i=1,nr1
    2   h(i)=h(i)*rofi(i)
        do 14 m=-l,l
        ilml=ilml+1
        do 16 isp=1,nsp
        sum=0d0
        do 15 ilmb=1,nlmb
C        if (rhoi(ilmb,isp) /= 0) then
C          call awrit3('ilmb=%i, isp=%i, rhoi=%d',' ',120,6,
C     .                ilmb,isp,rhoi(ilmb,isp))
C        endif
   15   sum=sum+s(ilml,ilmb)*rhoi(ilmb,isp)
        if (iov < 0) sum = -sum
        do 3 i=1,nr1
C            call awrit6('ilml=%i, l=%i, rhol=%d, sum=%d, h=%d, xi=%d',
C     .        ' ',120,6,ilml,l,rhol(i,ilml,isp),sum,h(i),xi(i,l))
    3   rhol(i,ilml,isp)=rhol(i,ilml,isp)+sum*h(i)*xi(i,l)
   16 continue
   14 continue
      end

