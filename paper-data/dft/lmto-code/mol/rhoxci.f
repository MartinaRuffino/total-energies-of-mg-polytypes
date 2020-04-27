      subroutine rhoxci(nbas,nspec,pos,rsm,rcut,rint,radd,nxi,lxi,
     .   exi,n0,ips,dx,dy,dz,atr,nc0,nc1,nc2,c1,c2,ioff,cy,nri,rhoi,
     .   rhoep,rhomu,vxci,dvxci,f)
c  Makes xc integrals for smooth interstitial rho.
c  Density is made for one z-plane at a time.
c  The xc-potential for 5 succesive planes is saved in omu to
c  permit evaluation of the derivative with respect to z.
      implicit real*8 (a-h,p-z), integer (o)
      dimension rsm(1),rcut(1),rint(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   rhoi(nri),vxci(nri),pos(3,1),ioff(1),cy(1),ips(1),nc1(1),
     .   nc2(1),c1(1),c2(1),ipmu(-2:2),dvxci(nri,3),f(3,nbas)
      real w(1)
      common /w/ w
      call getpr(ipr)
c ------ get center of the global mesh and z limits -------
      call ropxc0(dx,dy,dz,pos,nbas,ips,rint,cenx,ceny,cenz,m1,m2)
      if(ipr >= 20) write(6,903) cenx,ceny,cenz,m1,m2,dx,dy,dz
  903 format(/' rhoxci:  center=  ',3f9.3,'     m1,m2=',2i5
     .       /'          dx,dy,dz=',3f9.3)
      if(m1 > m2) call rx('rhoxci: empty mesh')
      sum1=0d0
      sum2=0d0
      call dpzero(vxci,   nri)
      call dpzero(dvxci,  nri*3)
c ------ get maximal dimensions in x and y directions -----
      kb= 10000
      kt=-10000
      lb= 10000
      lt=-10000
      do 10 m=m1,m2
      z=m*dz+cenz
      call roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,k1,k2,l1,l2)
      kb=min0(kb,k1)
      kt=max0(kt,k2)
      lb=min0(lb,l1)
  10  lt=max0(lt,l2)
      if(ipr >= 20) write(6,736) kb,kt,lb,lt,m1,m2
  736 format(' max extent of mesh:  (x)',2i5,'   (y)',2i5,'   (z)',2i5)
      nx0=kt-kb+1
      ny0=lt-lb+1
      write(6,936) nx0,ny0
  936 format(' nx0=',i5,'   ny0=',i5)
      call defrr(omu,   nx0*ny0*5)

c ------ make density,vxc for the first 3 planes ------
      do 14 i=-2,2
  14  ipmu(i)=i+3
      call xmuzer(1,nx0,ny0,w(omu))
      call xmuzer(2,nx0,ny0,w(omu))
      do 20 m=m1,m1+2
      ii=m-m1+3
      call xmuzer(ii,nx0,ny0,w(omu))
      z=m*dz+cenz
      call roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,k1,k2,l1,l2)
      nx=k2-k1+1
      ny=l2-l1+1
      if(nx <= 0.or.ny <= 0) goto 20
      call defrr(orhot,     nx*ny)
      call dpzero(w(orhot), nx*ny)
      call ropmro(dx,dy,dz,pos,nbas,rsm,rcut,rint,radd,nxi,lxi,exi,n0,
     .   atr,nc0,nc1,c1,nc2,c2,ips,ioff,cy,rhoi,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,w(orhot))
      call defrr(oh,    nx*ny)
      call defrr(ou2,   nx*ny)

      call rx('update call to ropxcc')
C     call ropxcc(nx*ny,w(orhot),w(oh),w(ou2),sum1,sum2)
      call xmucop(ii,k1,k2,l1,l2,w(orhot),kb,kt,lb,lt,w(omu))
      call rlse(orhot)
  20  continue

c ------ start big loop over z-planes ---------
      do 80 m=m1,m2
      z=m*dz+cenz
      call roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,k1,k2,l1,l2)
      nx=k2-k1+1
      ny=l2-l1+1
c --------- make integrals over vxc -------------------
      if(nx > 0.and.ny > 0) then
      ii=ipmu(0)
C|    write(6,963) m,ipmu(0)
C|963 format(' m=',i5,'   ipmu(0)=',i5)
      call defrr(ovxc,     nx*ny)
      call defrr(odvx,     nx*ny)
      call defrr(odvy,     nx*ny)
      call defrr(odvz,     nx*ny)
      call xmuder(ipmu,kb,kt,lb,lt,w(omu),k1,k2,l1,l2,w(ovxc),
     .    dx,dy,dz,w(odvx),w(odvy),w(odvz))
      call ropvr7(dx,dy,dz,pos,nbas,rsm,rcut,rint,radd,nxi,lxi,exi,n0,
     .   atr,nc0,nc1,c1,nc2,c2,ips,ioff,cy,vxci,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,w(ovxc),w(odvx),w(odvy),w(odvz),
     .   dvxci(1,1),dvxci(1,2),dvxci(1,3))
      call rlse(ovxc)
      endif

c --------- shuffle pointers to the mu-planes -----------
      ipnew=ipmu(-2)
      do 15 ii=-2,2
      ipmu(ii)=ipmu(ii)+1
  15  if(ipmu(ii) == 6) ipmu(ii)=1
C|    write(6,150) ipmu,ipnew
C|150 format(' ipmu=',5i5,'   ipnew=',i5)
c --------- make density and vxc in z-plane 3 higher up -------
      call xmuzer(ipnew,nx0,ny0,w(omu))
      mm=m+3
      if(mm <= m2) then
      z=mm*dz+cenz
      call roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,k1,k2,l1,l2)
      nx=k2-k1+1
      ny=l2-l1+1
      call defrr(orhot,     nx*ny)
      call dpzero(w(orhot), nx*ny)
      call ropmro(dx,dy,dz,pos,nbas,rsm,rcut,rint,radd,nxi,lxi,exi,n0,
     .   atr,nc0,nc1,c1,nc2,c2,ips,ioff,cy,rhoi,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,w(orhot))
      call defrr(oh,    nx*ny)
      call defrr(ou2,   nx*ny)
      call rx('update call to ropxcc')
C     call ropxcc(nx*ny,w(orhot),w(oh),w(ou2),sum1,sum2)
      call xmucop(ipnew,k1,k2,l1,l2,w(orhot),kb,kt,lb,lt,w(omu))
      call rlse(orhot)
      endif

  80  continue
      call rlse(omu)
c --------- multiply integrals by dx*dy*dz -------------
      sum1=sum1*dx*dy*dz
      sum2=sum2*dx*dy*dz
      sum3=0d0
      do 40 i=1,nri
      vxci(i)=vxci(i)*(dx*dy*dz)
  40  sum3=sum3+rhoi(i)*vxci(i)
      q=sum1
      rhoep=sum2
      rhomu=sum3
      if(ipr >= 20) write(6,732) kb,kt,lb,lt,m1,m2,q,rhoep,rhomu
  732 format(' max extent of mesh:  (x)',2i5,'   (y)',2i5,'   (z)',2i5
     .   /' q=',f12.6,'   rhoep=',f12.6,'   rhomu=',f12.6)

c --------- forces -------------------
      if(ipr >= 30) write(6,289)
      xx=dx*dy*dz
      do 30 ib=1,nbas
      is=ips(ib)
      i1=ioff(ib)+1
      i2=ioff(ib+1)
      f(1,ib)=0d0
      f(2,ib)=0d0
      f(3,ib)=0d0
      do 31 i=i1,i2
C|      top=dmax1(dabs(rhoi(i)),dabs(dvxci(i,1)),dabs(dvxci(i,2)))
C|      top=dmax1(top,dabs(dvxci(i,3)))
C|      if(top > 1d-6) write(6,971) i,rhoi(i),(xx*dvxci(i,m),m=1,3)
C|971   format(i6,4f13.6)
      do 31 m=1,3
  31  f(m,ib)=f(m,ib)-rhoi(i)*dvxci(i,m)*xx
      if(ipr >= 30) write(6,288) ib,(f(m,ib),m=1,3)
  288 format(i6,3f13.6)
  289 format(/'    ib     xc-force from smooth density')
  30  continue



      end
c ------- subs to handle array mu ---------
      subroutine xmucop(ii,k1,k2,l1,l2,vxc1,kb,kt,lb,lt,vxc2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension vxc1(k1:k2,l1:l2),vxc2(kb:kt,lb:lt,5)
      do 10 l=l1,l2
      do 10 k=k1,k2
  10  vxc2(k,l,ii)=vxc1(k,l)
          end
      subroutine xmuzer(ii,nx0,ny0,vxc2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension vxc2(nx0,ny0,5)
      call dpzero(vxc2(1,1,ii),nx0*ny0)
          end
c ------- xmuder: make mu and gradient for one z-plane ----
      subroutine xmuder(ipmu,kb,kt,lb,lt,vxc2,k1,k2,l1,l2,vxc,
     .    dx,dy,dz,dvx,dvy,dvz)
      implicit real*8 (a-h,p-z), integer (o)
      dimension vxc2(kb:kt,lb:lt,5),vxc(k1:k2,l1:l2),ipmu(-2:2),
     .   dvx(k1:k2,l1:l2),dvy(k1:k2,l1:l2),dvz(k1:k2,l1:l2)
      nx=k2-k1+1
      ny=l2-l1+1
c  copy over vxc
      ii0=ipmu(0)
      do 1 l=l1,l2
      do 1 k=k1,k2
  1   vxc(k,l)=vxc2(k,l,ii0)
c  x-differentiation
      w1=8d0/(12d0*dx)
      w2=1d0/(12d0*dx)
C|       w1=1d0/(2d0*dx)
C|       w2=0d0
      call dpzero(dvx,nx*ny)
      do 10 l=l1,l2
      do 12 k=k1+2,k2
  12  dvx(k,l)=dvx(k,l)+vxc(k-2,l)*w2
      do 13 k=k1+1,k2
  13  dvx(k,l)=dvx(k,l)-vxc(k-1,l)*w1
      do 14 k=k1,k2-1
  14  dvx(k,l)=dvx(k,l)+vxc(k+1,l)*w1
      do 15 k=k1,k2-2
  15  dvx(k,l)=dvx(k,l)-vxc(k+2,l)*w2
  10  continue
c  y-differentiation
      w1=8d0/(12d0*dy)
      w2=1d0/(12d0*dy)
      call dpzero(dvy,nx*ny)
      do 20 k=k1,k2
      do 22 l=l1+2,l2
  22  dvy(k,l)=dvy(k,l)+vxc(k,l-2)*w2
      do 23 l=l1+1,l2
  23  dvy(k,l)=dvy(k,l)-vxc(k,l-1)*w1
      do 24 l=l1,l2-1
  24  dvy(k,l)=dvy(k,l)+vxc(k,l+1)*w1
      do 25 l=l1,l2-2
  25  dvy(k,l)=dvy(k,l)-vxc(k,l+2)*w2
  20  continue
c  z-differentiation
      w1=8d0/(12d0*dz)
      w2=1d0/(12d0*dz)
C|      w1=1d0/(2d0*dz)
C|      w2=0d0
      call dpzero(dvz,nx*ny)
      im2=ipmu(-2)
      im1=ipmu(-1)
      ip1=ipmu(1)
      ip2=ipmu(2)
      do 30 l=l1,l2
      do 30 k=k1,k2
      dvz(k,l)=w2*(vxc2(k,l,im2)-vxc2(k,l,ip2))
     .        +w1*(vxc2(k,l,ip1)-vxc2(k,l,im1))
  30  continue
      end
