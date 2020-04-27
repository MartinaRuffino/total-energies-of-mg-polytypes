      subroutine vesmom(rmt,z,nr,a,lmxl,nbas,ips,orho,nvi,ioffv0,qmom)
c  Makes multipole moments of sphere densities.
c  28 Dec 94 spin pol (MvS)
      parameter( nrx=801, nlmx=49 )
      implicit real*8 (a-h,p-z), integer (o)
      integer intopt,nglob
      dimension rmt(1),z(1),nr(1),a(1),lmxl(1),ips(1),
     .   orho(nbas),qmom(nvi),ioffv0(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      fpi=16d0*datan(1d0)
      atepi=2d0*fpi
      y0=1d0/dsqrt(fpi)
      call defrr(orofi,    nrx)
      call defrr(orwgt,    nrx)
      call defrr(oh,       nrx)
      call dpzero(qmom,    nvi)
      if(ipr >= 50) write(6,221)
      intopt = 10*nglob('lrquad')

c ------ start loop over sites ----------
      do 10 ib=1,nbas
      is=ips(ib)
      j1=ioffv0(ib)+1
      lmaxl=lmxl(is)
      nlm=(lmaxl+1)**2
      oro=orho(ib)
      call radmsh(rmt(is),a(is),nr(is),w(orofi))
      call radwgt(intopt,rmt(is),a(is),nr(is),w(orwgt))
      call xxpotm(lmaxl,nr(is),w(orofi),w(orwgt),w(oro),w(oh),qmom(j1))
      qmom(j1)=qmom(j1)-z(is)*y0
      if(ipr >= 50) write(6,220) j1,ib,1,qmom(j1),qmom(j1)/y0,z(is)
      if(ipr >= 50) then
        do 11 ilm=2,nlm
        j=j1+ilm-1
  11    if(dabs(qmom(j)) > 1d-6) write(6,220) j,ib,ilm,qmom(j)
        endif
  220 format(i12,2i6,2f12.6,f8.2)
  221 format(/' vesmom:   j    ib   ilm      qmom',8x,'q',11x,'Z')
  10  continue
      call rlse(orofi)
      end
      subroutine xxpotm(lmaxl,nr1,rofi,rwgt,rho,h,qmom)
C- Numerically integrates density to make qmom
C  qmom(L) = r**l rho(r,L) r**2 dr: NB rho contains rho(r,L)*r**2
      implicit real*8 (a-h,p-z), integer (o)
      dimension rofi(1),rwgt(nr1),rho(nr1,1),h(nr1),qmom(1)
      nsp = lsp()+1
      do 1 i=1,nr1
    1 h(i)=rwgt(i)
      ilm=0
C     df=1d0
      do 10 l=0,lmaxl
C       df=df*(2*l+1)
        do 11 m=-l,l
          ilm=ilm+1
          sum=0d0
          do 2 i=1,nr1
    2     sum=sum+h(i)*rho(i,ilm)
          if (nsp == 2) then
            klm = ilm + (lmaxl+1)**2
            do 22 i=1,nr1
   22       sum=sum+h(i)*rho(i,klm)
          endif
   11     qmom(ilm)=sum
        do 3 i=1,nr1
    3   h(i)=h(i)*rofi(i)
   10 continue
      end
