      subroutine rodecz(nxi,lxi,exi,n0,rmt,rint,a,nr,lmxl,
     .  nbas,alat,pos,ips,rhoi,ioff,cg,jcg,indxcg,cy,ixp,dxp,orho)
c  Decompose density into atom-centered densities
      implicit real*8 (a-h,p-z), integer (o)
      integer intopt,nglob
      dimension exi(n0,1),nxi(1),ips(1),lxi(n0,1),ioff(1),
     .  rhoi(1),pos(3,1),rint(1),cg(1),jcg(1),indxcg(1),cy(1),
     .  rmt(1),lmxl(1),a(1),nr(1),orho(1),dr(3),ojkl(0:20),ixp(1),dxp(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('rodecz')
      pi=4d0*datan(1d0)
      srfpi=dsqrt(4d0*pi)
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
      call dpdot(w(oro),w(orwgt),nr1,q1)
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
      if(nlmb > nlmbx) call rx('rodecz: increase nlmbx')
C|    call mstrux(e,dr,w(os),nlml,nlmb,nlml,cg,indxcg,jcg,cy)
      ojj=ojkl(lb)
C      call nstrux(e,dr,nlml,nlmb,nlmp,npow,w(oikl),w(ojj),w(oip),
C     .   w(ocf),cy,w(os))
C      ih = nlmp*(je-1)
C      call hstrux(e,nlml,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
C     .  w(ocf),w(ohl),w(os))
      call hstrux(e,nlml,nlmb,nlmp,npow,je,je,w(oikl),w(ojj),w(oip),
     .  w(ocf),w(ohl),w(os))
      call ropbes(w(orofi),e,lmaxl,w(oy),w(oh),w(oxi),nr1,1)
      call rodcm1(nlml,nlmb,w(os),nr1,rhoi(j1),w(oxi),w(oh),
     .   w(orofi),w(oro))
  11  j1=j1+nlmb
      call rlse(ohl)
  20  continue
      call rlse(orofi)
      call dpdot(w(oro),w(orwgt),nr1,qtot)
      if(ipr >= 35) write(6,787) srfpi*q1,srfpi*(qtot-q1),srfpi*qtot
  787 format(' rodecz:  q1,q2,qtot=',3f12.6)
  10  continue
      call tcx('rodecz')
      end

c ------ sub rodcm1 ------------------------
      subroutine rodcm1(nlml,nlmb,s,nr1,rhoi,xi,h,rofi,rhol)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhol(nr1,1),rofi(1),xi(nr1,0:1),s(nlml,nlmb),
     .   rhoi(1),h(1)
      lmaxl=ll(nlml)
      do 1 i=1,nr1
  1   h(i)=rofi(i)
      ilml=0
      do 14 l=0,lmaxl
        do 2 i=1,nr1
  2     h(i)=h(i)*rofi(i)
        do 14 m=-l,l
        ilml=ilml+1
        sum=0d0
        do 15 ilmb=1,nlmb
  15    sum=sum+s(ilml,ilmb)*rhoi(ilmb)
        do 3 i=1,nr1
  3     rhol(i,ilml)=rhol(i,ilml)-sum*h(i)*xi(i,l)
  14    continue
      end
