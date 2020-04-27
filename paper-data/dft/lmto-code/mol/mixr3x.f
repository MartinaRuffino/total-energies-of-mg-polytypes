      subroutine mixr3x(beta,scale,nbas,ips,rmt,nr,a,lmxl,
     .   orhut,orhos,nri,nvi,rhoix,rhoi,rho0x,rho0,diff)
C- Mix density
      parameter( ncmx=20, nlmx=49, npwr=5 )
      implicit real*8 (a-h,p-z), integer (o)
      integer intopt,nglob
      dimension ips(1),rmt(1),nr(1),a(1),lmxl(1),orhut(1),orhos(1),
     .   rhoi(nri),rhoix(nri),rho0(1),rho0x(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      diff = 0
      do 10 ib=1,nbas
      is=ips(ib)
      nr1=nr(is)
      nlml=(lmxl(is)+1)**2
      ornew=orhut(ib)
      orold=orhos(ib)
C ... scale only interstitial now
C     call dpcopy(w(ornew),w(ornew),1,nr1*nlml,scale)
c ------ this part to compare old and new density ----
      call defrr(orofi,   nr1)
      call defrr(orwgt,   nr1)
      call radmsh(rmt(is),a(is),nr1,w(orofi))
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rmt(is),a(is),nr1,w(orwgt))
      print *, 'old density ...'
      call rhinfo(nr1,nlml,nsp,w(orold),w(orofi),w(orwgt),a(is))
      print *, 'new density ...'
      call rhinfo(nr1,nlml,nsp,w(ornew),w(orofi),w(orwgt),a(is))

      call defrr(oh,    nr1)
      call mixrhl(beta,w(ornew),w(orold),nr1,nlml,a(is),rmt(is),
     .   diffi,w(oh))
      diff = diff + diffi**2
  10  call rlse(oh)
      diff = dsqrt(diff/nbas)
c -------- mix interstitial density -------------
      do 20 i=1,nri
   20 rhoix(i)=scale*rhoix(i)
      do 21 i = 1, nvi
   21 rho0x(i)=scale*rho0x(i)
      top=0d0
      do 22 i=1,nri
      top=dmax1(top,dabs(rhoix(i)-rhoi(i)))
   22 rhoi(i)=beta*rhoix(i)+(1d0-beta)*rhoi(i)
      top0=0d0
      do 23 i=1,nvi
      top0=dmax1(top0,dabs(rho0x(i)-rho0(i)))
   23 rho0(i)=beta*rho0x(i)+(1d0-beta)*rho0(i)
      write(6,727) top, top0
  727 format(' top,top0=',2f12.6)

      end
