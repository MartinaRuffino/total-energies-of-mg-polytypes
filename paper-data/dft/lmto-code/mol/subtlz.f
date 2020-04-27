      subroutine subtlz(nxi,lxi,exi,n0,rint,nbas,alat,pos,ips,ioff,
     .   cg,jcg,indxcg,cy,pos1,r1,nr1,nr2,rofi,rwgt,nlml,nlmx,ixp,dxp,
     .   vl,vxc1,rl1,dvxc2,zeta,
     .   zetxc,ib,rhoi,poti,pot0,ioffv0,vval,qmom,for1,for2,for4)
c  Subtracts integrals of xi-tails times pot from zeta,zetxc.
c  Contributions to force:       f1= d integral( rho * eps )
c  f2=integral ( ves * drho)     f3=integral( rho * dves )
c  f4,f5: force from change in qval and vval
c  f6 is derivative of interstitial charge; gets scaled with vint later
c  f7=force from interstitial xc ring integrals  (15.4.92)
      implicit real*8 (a-h,p-z), integer (o)
      parameter( nlmxx=49 )
      dimension exi(n0,1),nxi(1),ips(1),lxi(n0,1),ioff(1),f7(3),f6(3),
     .   pos(3,1),rint(1),cg(1),jcg(1),indxcg(1),cy(1),dr(3),pos1(3),
     .   vl(nr1,1),vxc1(nr1,1),rofi(1),rwgt(1),zeta(1),zetxc(1),f1(3),
     .   rhoi(1),f2(3),poti(1),pot0(1),ioffv0(1),rl1(nr1,1),f3(3),f4(3),
     .   sum(nlmxx),vval(1),qmom(1),f5(3),xval(nlmxx),xmom(nlmxx),
     .   dvxc2(nr2,1),for1(3,nbas),for2(3,nbas),for4(3,nbas),ojkl(0:20),
     .   ixp(1),dxp(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      lmaxl=ll(nlml)
      nrwk=max0(nr1,nr2)
      is=ips(ib)
      i0=ioffv0(ib)+1
      i1=ioff(ib)+1
      call defrr(oxi,         nrwk*(lmaxl+1))
      call defrr(oh,          nrwk)
      call defrr(oy,          nrwk)
c ... setup for strux
      call strxsu(nlmx,nxi,lxi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      call defrr(os,          nlmx*nlmbx)
      call defrr(ogs,         nlmx*nlmbx*3)
      call defrr(osd,         nlmx*nlmbx)
      call defrr(ogd,         nlmx*nlmbx*3)
c ... one-center electrostatic terms
      call stles1(r1,nxi(is),lxi(1,is),exi(1,is),poti(i1),pot0(i0),
     .  zeta(i1))

c ------ start loop over other atoms ---------------
      do 20 jb=1,nbas
      js=ips(jb)
      call dpdist(pos1,pos(1,jb),3,dd)
      if(alat*dd < 1d-4 .and. ixp(1) == 0) goto 20
      dr(1)=alat*(pos(1,jb)-pos1(1))
      dr(2)=alat*(pos(2,jb)-pos1(2))
      dr(3)=alat*(pos(3,jb)-pos1(3))
      do 5 m=1,3
      f7(m)=0d0
      f1(m)=0d0
      f2(m)=0d0
      f3(m)=0d0
      f4(m)=0d0
      f5(m)=0d0
  5   f6(m)=0d0
      lmax=-1
      do 2 je=1,nxi(js)
  2   lmax=max0(lmax,lxi(je,js))
c ---------- loop over non-zero xi energies --------
      if(alat*dd <= rint(js)+r1) then
      j1=ioff(jb)+1
      kr=nr1+1
      lmxh = ll(nlmx)+lmax+1
      nlmh = (lmxh+1)**2
      call defrr(ohl,nlmh*nxi(js))
      call defrr(ohd,nlmh*nxi(js))
      call rstr0(nxi(js),lxi(1,js),exi(1,js),nlmh,1,dr(1),dr(2),dr(3),
     .  ll(nlmx)+1,0,w(ohl),w(ohd))
C     lmax=-1
      do 11 je=1,nxi(js)
      e=exi(je,js)
      lb=lxi(je,js)
C     lmax=max0(lmax,lb)
      nlmb=(lb+1)**2
      if(nlmb > nlmxx) call rx('subtlz: increase nlmbx')
      ojj=ojkl(lb)
C     call nstrpg(e,dr,nlmx,nlmb,nlmp,npow,w(oikl),w(ojj),w(oip),w(ocf),
C    .   cy,w(os),w(osd),w(ogs),w(ogd))
      ih = nlmh*(je-1)
      call hstrpg(e,nlmx,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .  w(ocf),w(ohl),w(ohd),w(os),w(osd),w(ogs),w(ogd))
c  2-center electrostatic terms
      call stles2(r1,nxi(is),lxi(1,is),exi(1,is),nlmx,nlmb,w(os),
     .   w(osd),e,poti(i1),poti(j1),pot0(i0),zeta(i1),zeta(j1),
     .   rhoi(i1),rhoi(j1),w(ogs),w(ogd),f2)
c  zetxc and force from rho shifting against eps
      call stlset(e,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,vxc1,sum)
      call stlxc1(nlmx,nlml,nlmb,sum,w(os),-1d0,zetxc(j1))
      call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhoi(j1),f1)
      call stlset(e,lmaxl,nr2,w(oxi),w(oh),w(oy),rofi(kr),rwgt(kr),
     .   dvxc2,sum)
      call stlxc1(nlmx,nlml,nlmb,sum,w(os),1d0,zetxc(j1))
      call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhoi(j1),f7)
c  zeta and force from rho shifting against vl
      call stlset(e,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,vl,sum)
      call stlxc1(nlmx,nlml,nlmb,sum,w(os),-1d0,zeta(j1))
      call stlfor(nlmx,nlml,nlmb,sum,w(ogs),rhoi(j1),f2)
c  force from ves shifting against rho
      call stlset(e,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,rl1,sum)
      call stlfor(nlmx,nlml,nlmb,sum,w(ogs),poti(j1),f3)
c  force from vval and qmom, also deriv of sphere charge
      call stlqi1(e,r1,lmaxl,vval,qmom,xval,xmom,x0)
      call stlfor(nlmx,nlml,nlmb,xval,w(ogs),rhoi(j1),f4)
      call stlfor(nlmx,nlml,nlmb,xmom,w(ogs),poti(j1),f5)
      call stlqi2(nlmx,nlml,nlmb,x0,w(ogs),rhoi(j1),f6)
  11  j1=j1+nlmb
      call rlse(ohl)
      endif
c ---------- functions with energy zero -------
      j0=ioffv0(jb)+1
      nlmb=(lmax+1)**2
      ojj=ojkl(lmax)
C      call nstrug(0d0,dr,nlmx,nlmb,nlmp,0   ,w(oikl),w(ojj),w(oip),
C     .   w(ocf),cy,w(os),w(ogs))
      call defrr(ohl,nlmh*nxi(js))
      call rstr0(1,lmax,0d0,nlmh,1,dr(1),dr(2),dr(3),ll(nlmx)+1,1,
     .           w(ohl),w(1))
      call hstrug(0d0,nlmx,nlmb,nlmp,0,0,w(oikl),w(ojj),w(oip),
     .  w(ocf),w(ohl),w(os),w(ogs))
      call rlse(ohl)

c  2-center electrostatic terms
      call stles0(r1,nxi(is),lxi(1,is),exi(1,is),nlmx,nlmb,w(os),
     .   pot0(j0),zeta(i1),rhoi(i1),w(ogs),f2)
c  force from ves shifting against rho
      call stlset(0d0,lmaxl,nr1,w(oxi),w(oh),w(oy),rofi,rwgt,rl1,sum)
      call stlfor(nlmx,nlml,nlmb,sum,w(ogs),pot0(j0),f3)
c  force from vval and qmom
      call stlqi1(0d0,r1,lmaxl,vval,qmom,xval,xmom,x0)
      call stlfor(nlmx,nlml,nlmb,xmom,w(ogs),pot0(j0),f5)

c ---------- add together force contributions ---------------
      if(ipr >= 60) write(6,809) ib,jb,f1,f2,f3,f4,f5,f6
  809 format(' ib,jb=',2i3,'   f1=',3f12.6/13x,'   f2=',3f12.6/13x,
     .   '   f3=',3f12.6/13x,'   f4=',3f12.6/13x,'   f5=',3f12.6
     .   /13x,'   f6=',3f12.6)
      if(ipr >= 60) write(6,899) f7
  899 format(13x,'   f7=',3f12.6)
      do 25 m=1,3
      for1(m,ib)=for1(m,ib)-f1(m)   +f7(m)
      for1(m,jb)=for1(m,jb)+f1(m)   -f7(m)
      fes=0.5d0*(f2(m)+f3(m)+f4(m)-f5(m))
      for2(m,ib)=for2(m,ib)-fes
      for2(m,jb)=for2(m,jb)+fes
      for4(m,ib)=for4(m,ib)-f6(m)
  25  for4(m,jb)=for4(m,jb)+f6(m)
  20  continue

      call rlse(oxi)
      end
c ------ stlset: set up integrals bessels*potential ----
      subroutine stlset(e,lmaxl,nr1,xi,h,y,rofi,rwgt,vxc1,sum)
      implicit real*8 (a-h,p-z), integer (o)
      dimension vxc1(nr1,1),rofi(1),xi(nr1,0:1),sum(1),
     .   y(1),h(1),rwgt(1)
      call ropbes(rofi,e,lmaxl,y,h,xi,nr1,1)
      do 1 i=1,nr1
  1   h(i)=rofi(i)*rwgt(i)
      ilm=0
      do 16 l=0,lmaxl
        do 2 i=1,nr1
  2     h(i)=h(i)*rofi(i)
        do 14 m=-l,l
          ilm=ilm+1
          sam=0d0
          do 3 i=1,nr1
  3       sam=sam+h(i)*xi(i,l)*vxc1(i,ilm)
  14    sum(ilm)=sam
  16    continue
      end
c ------ stlxc1: add to zetxc ---------------
      subroutine stlxc1(nlmx,nlml,nlmb,sum,s,fac,zetxc)
      implicit real*8 (a-h,p-z), integer (o)
      dimension sum(1),s(nlmx,nlmb),zetxc(1)
      do 15 jlm=1,nlmb
      do 15 ilm=1,nlml
  15  zetxc(jlm)=zetxc(jlm)+fac*s(ilm,jlm)*sum(ilm)
      end
c ------ stlfor: add to force -----------------
      subroutine stlfor(nlmx,nlml,nlmb,sum,gs,rhoi,f)
      implicit real*8 (a-h,p-z), integer (o)
      dimension sum(1),rhoi(1),f(3),gs(nlmx,nlmb,3)
      do 15 jlm=1,nlmb
      do 15 ilm=1,nlml
      f(1)=f(1)+gs(ilm,jlm,1)*rhoi(jlm)*sum(ilm)
      f(2)=f(2)+gs(ilm,jlm,2)*rhoi(jlm)*sum(ilm)
      f(3)=f(3)+gs(ilm,jlm,3)*rhoi(jlm)*sum(ilm)
  15  continue
      end
c ------ stlqi2: only for derivative of instl charge ----
      subroutine stlqi2(nlmx,nlml,nlmb,x0,gs,rhoi,f)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhoi(1),f(3),gs(nlmx,nlmb,3)
      do 15 jlm=1,nlmb
      f(1)=f(1)-gs(1,jlm,1)*rhoi(jlm)*x0
      f(2)=f(2)-gs(1,jlm,2)*rhoi(jlm)*x0
      f(3)=f(3)-gs(1,jlm,3)*rhoi(jlm)*x0
  15  continue
      end
c ------ stlqi1: setup for force terms from qmom and vval ------
      subroutine stlqi1(e,r1,lmaxl,vval,qmom,xval,xmom,x0)
      implicit real*8 (a-h,p-z), integer (o)
      dimension xval(1),xmom(1),phi(0:20),psi(0:20),vval(1),qmom(1)
      call bessl(e*r1*r1,lmaxl+1,phi,psi)
      ilm=0
      do 14 l=0,lmaxl
      bval=phi(l)*r1**l
      bmom=phi(l+1)*r1**(2*l+3)
      if(l == 0) x0=bmom*3.544907702d0
      do 14 m=-l,l
      ilm=ilm+1
      xval(ilm)=vval(ilm)*(bmom/r1**l)
      xmom(ilm)=qmom(ilm)*(bval/r1**l)
  14  continue
      end
c ------ stles2: add 2-center terms to zeta -----
      subroutine stles2(r1,nxi,lxi,exi,nlmx,nlmb,s,sd,e,poti,potj,
     .   pot0,zeti,zetj,rhoi,rhoj,gs,gd,f2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension poti(1),potj(1),s(nlmx,nlmb),zeti(1),zetj(1),
     .   exi(1),lxi(1),fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20),
     .   sd(nlmx,nlmb),pot0(1),rhoi(1),rhoj(1),gs(nlmx,nlmb,3),
     .   gd(nlmx,nlmb,3),f2(3)
      i=0
      lmax=-1
      do 10 ie=1,nxi
      lx=lxi(ie)
      lmax=max0(lmax,lx)
      ex=exi(ie)
      call wronkj(ex,e,r1,lx,fkk,fkj,fjk,fjj)
c ... for ex /= e
      if(dabs(e-ex) > 1d-8) then
        ilm=0
        do 11 l=0,lx
        fac=fkj(l)+1d0/(e-ex)
        do 11 m=-l,l
        i=i+1
        ilm=ilm+1
        do 14 jlm=1,nlmb
        xxx=-fac*(rhoj(jlm)*poti(i)+rhoi(i)*potj(jlm))
        f2(1)=f2(1)+xxx*gs(ilm,jlm,1)
        f2(2)=f2(2)+xxx*gs(ilm,jlm,2)
  14    f2(3)=f2(3)+xxx*gs(ilm,jlm,3)
        do 11 jlm=1,nlmb
        zetj(jlm)=zetj(jlm)+s(ilm,jlm)*fac*poti(i)
  11    zeti(i  )=zeti(i  )+s(ilm,jlm)*fac*potj(jlm)
c ... for e == ex
      else
        ilm=0
        do 12 l=0,lx
        fac=fkj(l)
        do 12 m=-l,l
        i=i+1
        ilm=ilm+1
        do 16 jlm=1,nlmb
        xxx=-(rhoj(jlm)*poti(i)+rhoi(i)*potj(jlm))
        f2(1)=f2(1)+xxx*(fac*gs(ilm,jlm,1)+0.5d0*gd(ilm,jlm,1))
        f2(2)=f2(2)+xxx*(fac*gs(ilm,jlm,2)+0.5d0*gd(ilm,jlm,2))
  16    f2(3)=f2(3)+xxx*(fac*gs(ilm,jlm,3)+0.5d0*gd(ilm,jlm,3))
C ...   This code causes problems when ib=jb.
C        do 12 jlm=1,nlmb
C        sx=s(ilm,jlm)*fac+sd(ilm,jlm)*0.5d0
C        zetj(jlm)=zetj(jlm)+sx*poti(i)
C  12    zeti(i  )=zeti(i  )+sx*potj(jlm)
C ...   Use this instead
        wk = 0d0
        do 13 jlm=1,nlmb
        sx=s(ilm,jlm)*fac+sd(ilm,jlm)*0.5d0
        zetj(jlm)=zetj(jlm)+sx*poti(i)
   13   wk       =wk       +sx*potj(jlm)
        zeti(i)  =zeti(i)+wk
   12   continue
      endif
  10  continue
c ... terms from v0
      call wronkj(0d0,e,r1,lmax,fkk,fkj,fjk,fjj)
      ilm=0
      do 21 l=0,lmax
      fac=fkj(l)+1d0/e
      do 21 m=-l,l
      ilm=ilm+1
      do 18 jlm=1,nlmb
      xxx=-fac*rhoj(jlm)*pot0(ilm)
      f2(1)=f2(1)+xxx*gs(ilm,jlm,1)
      f2(2)=f2(2)+xxx*gs(ilm,jlm,2)
  18  f2(3)=f2(3)+xxx*gs(ilm,jlm,3)
      do 21 jlm=1,nlmb
  21  zetj(jlm)=zetj(jlm)+s(ilm,jlm)*fac*pot0(ilm)
      end
c ------ stles1: add 1-center terms to zeta -----
      subroutine stles1(r1,nxi,lxi,exi,poti,pot0,zeti)
      implicit real*8 (a-h,p-z), integer (o)
      dimension poti(1),zeti(1),exi(1),lxi(1),pot0(1),
     .   fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20)
      i0=0
      do 10 ie=1,nxi
      li=lxi(ie)
      j0=0
      lmax=-1
      do 11 je=1,nxi
      lj=lxi(je)
      lmax=max0(lmax,lj)
      lx=min0(li,lj)
      call wronkj(exi(ie),exi(je),r1,lx,fkk,fkj,fjk,fjj)
      do 20 ilm=1,(lx+1)**2
      l=ll(ilm)
  20  zeti(ilm+i0)=zeti(ilm+i0)+poti(ilm+j0)*fkk(l)
  11  j0=j0+(lj+1)**2
      lx=min0(li,lmax)
      call wronkj(exi(ie),0d0,r1,lx,fkk,fkj,fjk,fjj)
      do 24 ilm=1,(lx+1)**2
      l=ll(ilm)
  24  zeti(ilm+i0)=zeti(ilm+i0)+pot0(ilm)*fkk(l)
  10  i0=i0+(li+1)**2
      end
c ------ stles0: add 2-center terms to zeta, e=0 terms ----
      subroutine stles0(r1,nxi,lxi,exi,nlmx,nlmb,s,pot0,zeti,
     .   rhoi,gs,f2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension pot0(1),s(nlmx,nlmb),zeti(1),gs(nlmx,nlmb,3),rhoi(1),
     .   exi(1),lxi(1),fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20),f2(3)
      i=0
      do 10 ie=1,nxi
      lx=lxi(ie)
      ex=exi(ie)
      call wronkj(ex,0d0,r1,lx,fkk,fkj,fjk,fjj)
        ilm=0
        do 11 l=0,lx
        fac=fkj(l)-1d0/ex
        do 11 m=-l,l
        i=i+1
        ilm=ilm+1
        do 18 jlm=1,nlmb
        xxx=-fac*rhoi(i)*pot0(jlm)
        f2(1)=f2(1)+xxx*gs(ilm,jlm,1)
        f2(2)=f2(2)+xxx*gs(ilm,jlm,2)
  18    f2(3)=f2(3)+xxx*gs(ilm,jlm,3)
        do 11 jlm=1,nlmb
  11    zeti(i)=zeti(i)+s(ilm,jlm)*fac*pot0(jlm)
  10  continue
      end
