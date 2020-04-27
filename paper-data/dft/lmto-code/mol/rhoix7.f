      subroutine rhoix7(nbas,nspec,pos,rsm,rint,nxi,lxi,exi,n0,
     .   ips,dx,dy,dz,ioff,cy,cg,jcg,indxcg,nri,rhoi,rhoep,rhomu,vxci,f)
c  Makes xc integrals for smooth interstitial rho.
c  Density is made for one z-plane at a time.
c  The xc-potential for 7 succesive planes is saved in omu to
c  permit evaluation of the derivative with respect to z.
Cu Updates
Cu   21 Jul 07 (S. Lozovoi) New gradient corrections
      implicit real*8 (a-h,p-z), integer (o)
      dimension rsm(1),rint(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   rhoi(nri,*),vxci(nri,*),pos(3,1),ioff(1),cy(1),ips(1),
     .   ipmu(-3:3),f(3,nbas),
     .   cg(1),jcg(1),indxcg(1)
      real w(1)
      common /w/ w
      call tcn('rhoix7')
      call getpr(ipr)
c ------ get center of the global mesh and z limits -------
      call ropxc0(dx,dy,dz,pos,nbas,ips,rint,cenx,ceny,cenz,m1,m2)
      if(ipr >= 20) write(6,903) cenx,ceny,cenz,m1,m2,dx,dy,dz
  903 format(/' rhoix7:  center=  ',3f9.3,'     m1,m2=',2i5
     .       /'          dx,dy,dz=',3f9.3)
      if(m1 > m2) call rx('rhoix7: empty mesh')
      sum1=0d0
      sum2=0d0
      ntot=0
      nsp=lsp()+1
      call dpzero(vxci,   nri*nsp)
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
      call defrr(omu,   nx0*ny0*7*nsp)

c ------ make density,vxc for the first 4 planes ------
      do 14 i=-3,3
  14  ipmu(i)=i+4
      call ymuzer(1,nx0,ny0,w(omu))
      call ymuzer(2,nx0,ny0,w(omu))
      call ymuzer(3,nx0,ny0,w(omu))
      do 20 m=m1,m1+3
      ii=m-m1+4
      call ymuzer(ii,nx0,ny0,w(omu))
      z=m*dz+cenz
      call roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,k1,k2,l1,l2)
      nx=k2-k1+1
      ny=l2-l1+1
      if(nx <= 0.or.ny <= 0) goto 20
      call defrr(orhot,    -nx*ny*nsp)
      if(lxcg() > 0) then
        call defrr(ogrhot,    -nx*ny*3*nsp)
        call defrr(og2rht,    -nx*ny*5*nsp)
        call defrr(oggrht,    -nx*ny*nsp)
      else
        ogrhot = 1
        og2rht = 1
        oggrht = 1
      endif
      call rxcmro(dx,dy,dz,pos,nbas,rsm,rint,nxi,lxi,exi,n0,ips,ioff,
     .   nri,rhoi,cg,jcg,indxcg,cenx,ceny,cenz,z,k1,k2,l1,l2,w(orhot),
     .   w(ogrhot),w(og2rht),w(oggrht),ntot)
      call defrr(oh,    nx*ny)
      call defrr(ou2,   nx*ny*nsp)
      call defrr(orhos, nx*ny)
      if(lxcg() > 0) then
        call defrr(oagrh, nx*ny*(2*nsp-1))
        call defrr(ogragr,nx*ny*(2*nsp-1))
        call defrr(ogrgr, nx*ny)
      else
        oagrh  = 1
        ogragr = 1
        ogrgr  = 1
      endif
      call ropxcc(nx*ny,w(orhos),w(orhot),
     .   w(ogrhot),w(og2rht),w(oggrht),w(oagrh),w(ogragr),w(ogrgr),
     .   w(oh),w(ou2),sum1,sum2,nx*ny)
      call ymucop(ii,k1,k2,l1,l2,w(orhot),kb,kt,lb,lt,w(omu))
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
C|    write(6,963) m,ipmu(0)
C|963 format(' m=',i5,'   ipmu(0)=',i5)
      call defrr(ovxc,     nx*ny*nsp)
      call defrr(odvx,     nx*ny*nsp)
      call defrr(odvy,     nx*ny*nsp)
      call defrr(odvz,     nx*ny*nsp)
      call ymuder(ipmu,kb,kt,lb,lt,w(omu),k1,k2,l1,l2,w(ovxc),
     .   dx,dy,dz,w(odvx),w(odvy),w(odvz))
      call rxcvrf(dx,dy,dz,pos,nbas,rsm,rint,nxi,lxi,exi,n0,
     .   ips,ioff,cy,vxci,cenx,ceny,cenz,z,k1,k2,l1,l2,
     .   w(ovxc),w(odvx),w(odvy),w(odvz),nri,rhoi,f)
      call rlse(ovxc)
      endif

c --------- shuffle pointers to the mu-planes -----------
      ipnew=ipmu(-3)
      do 15 ii=-3,3
      ipmu(ii)=ipmu(ii)+1
  15  if(ipmu(ii) == 8) ipmu(ii)=1
C|    write(6,150) ipmu,ipnew
C|150 format(' ipmu=',5i5,'   ipnew=',i5)
c --------- make density and vxc in z-plane 3 higher up -------
      call ymuzer(ipnew,nx0,ny0,w(omu))
      mm=m+4
      if(mm <= m2) then
      z=mm*dz+cenz
      call roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,k1,k2,l1,l2)
      nx=k2-k1+1
      ny=l2-l1+1

      call defrr(orhot,    -nx*ny*nsp)
      if(lxcg() > 0) then
        call defrr(ogrhot,    -nx*ny*3*nsp)
        call defrr(og2rht,    -nx*ny*5*nsp)
        call defrr(oggrht,    -nx*ny*nsp)
      else
        ogrhot = 1
        og2rht = 1
        oggrht = 1
      endif
      call rxcmro(dx,dy,dz,pos,nbas,rsm,rint,nxi,lxi,exi,n0,ips,ioff,
     .   nri,rhoi,cg,jcg,indxcg,cenx,ceny,cenz,z,k1,k2,l1,l2,w(orhot),
     .   w(ogrhot),w(og2rht),w(oggrht),ntot)

      call defrr(oh,    nx*ny)
      call defrr(ou2,   nx*ny*nsp)
      call defrr(orhos, nx*ny)
      if(lxcg() > 0) then
        call defrr(oagrh, nx*ny*(2*nsp-1))
        call defrr(ogragr,nx*ny*(2*nsp-1))
        call defrr(ogrgr, nx*ny)
      else
        oagrh  = 1
        ogragr = 1
        ogrgr  = 1
      endif

      call ropxcc(nx*ny,w(orhos),w(orhot),
     .   w(ogrhot),w(og2rht),w(oggrht),w(oagrh),w(ogragr),w(ogrgr),
     .   w(oh),w(ou2),sum1,sum2,nx*ny)
      call ymucop(ipnew,k1,k2,l1,l2,w(orhot),kb,kt,lb,lt,w(omu))
      call rlse(orhot)
      endif

  80  continue
      call rlse(omu)
c ------- multiply integrals by dx*dy*dz -------------
      sum1=sum1*dx*dy*dz
      sum2=sum2*dx*dy*dz
c
      rhomu=0.d0
      do isp=1,nsp
        sum3=0d0
        do i=1,nri
          vxci(i,isp)=vxci(i,isp)*(dx*dy*dz)
          sum3=sum3+rhoi(i,isp)*vxci(i,isp)
        enddo
        rhomu = rhomu + sum3
      enddo
      q=sum1
      rhoep=sum2
c     rhomu=sum3
      if(ipr >= 20) write(6,732) q,rhoep,rhomu
  732 format(' q=',f12.6,'   rhoep=',f12.6,'   rhomu=',f12.6)
      if(ipr >= 30) write(6,782) ntot
  782 format(' total points evaluated:',i10)
c ------- print forces -------------
      if(ipr >= 40) then
      write(6,289)
      do 30 ib=1,nbas
  30  write(6,288) ib,(f(m,ib),m=1,3)
  288 format(i6,3f13.6)
  289 format(/'    ib     xc-force from smooth density')
      endif

      call tcx('rhoix7')
      end
c ------- subs to handle array mu ---------
      subroutine ymucop(ii,k1,k2,l1,l2,vxc1,kb,kt,lb,lt,vxc2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension vxc1(k1:k2,l1:l2,*),vxc2(kb:kt,lb:lt,7,*)
c     call tcn('ymucop')
      nsp=lsp()+1
      do isp=1,nsp
        do 10 l=l1,l2
          do 10 k=k1,k2
  10        vxc2(k,l,ii,isp)=vxc1(k,l,isp)
      enddo
c     call tcx('ymucop')
          end
      subroutine ymuzer(ii,nx0,ny0,vxc2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension vxc2(nx0,ny0,7,*)
c     call tcn('ymuzer')
      nsp=lsp()+1
      do isp=1,nsp
        call dpzero(vxc2(1,1,ii,isp),nx0*ny0)
      enddo
c     call tcx('ymuzer')
          end
c ------- ymuder: make mu and gradient for one z-plane ----
      subroutine ymuder(ipmu,kb,kt,lb,lt,vxc2,k1,k2,l1,l2,vxc,
     .    dx,dy,dz,dvx,dvy,dvz)
      implicit real*8 (a-h,p-z), integer (o)
      dimension vxc2(kb:kt,lb:lt,7,*),vxc(k1:k2,l1:l2,*),ipmu(-3:3),
     .   dvx(k1:k2,l1:l2,*),dvy(k1:k2,l1:l2,*),dvz(k1:k2,l1:l2,*)
c     dsinh(x)=0.5d0*(dexp(x)-dexp(-x))
c     call tcn('ymuder')
      nsp=lsp()+1
      nx=k2-k1+1
      ny=l2-l1+1
      call dpzero(dvx,nx*ny*nsp)
      call dpzero(dvy,nx*ny*nsp)
      call dpzero(dvz,nx*ny*nsp)
      ii0=ipmu(0)
c------start big loop over spins----------------------------------
      do isp=1,nsp
c  copy over vxc for each spin
      do 1 l=l1,l2
      do 1 k=k1,k2
  1   vxc(k,l,isp)=vxc2(k,l,ii0,isp)
c  x-differentiation
      w1= 3d0/(4d0*dx)
      w2=-3d0/(20d0*dx)
      w3= 1d0/(60d0*dx)
      do 10 l=l1,l2
      do 17 k=k1+3,k2
  17  dvx(k,l,isp)=dvx(k,l,isp)-vxc(k-3,l,isp)*w3
      do 12 k=k1+2,k2
  12  dvx(k,l,isp)=dvx(k,l,isp)-vxc(k-2,l,isp)*w2
      do 13 k=k1+1,k2
  13  dvx(k,l,isp)=dvx(k,l,isp)-vxc(k-1,l,isp)*w1
      do 14 k=k1,k2-1
  14  dvx(k,l,isp)=dvx(k,l,isp)+vxc(k+1,l,isp)*w1
      do 15 k=k1,k2-2
  15  dvx(k,l,isp)=dvx(k,l,isp)+vxc(k+2,l,isp)*w2
      do 16 k=k1,k2-3
  16  dvx(k,l,isp)=dvx(k,l,isp)+vxc(k+3,l,isp)*w3
  10  continue
c  y-differentiation
      w1= 3d0/(4d0*dy)
      w2=-3d0/(20d0*dy)
      w3= 1d0/(60d0*dy)
      do 20 k=k1,k2
      do 28 l=l1+3,l2
  28  dvy(k,l,isp)=dvy(k,l,isp)-vxc(k,l-3,isp)*w3
      do 22 l=l1+2,l2
  22  dvy(k,l,isp)=dvy(k,l,isp)-vxc(k,l-2,isp)*w2
      do 23 l=l1+1,l2
  23  dvy(k,l,isp)=dvy(k,l,isp)-vxc(k,l-1,isp)*w1
      do 24 l=l1,l2-1
  24  dvy(k,l,isp)=dvy(k,l,isp)+vxc(k,l+1,isp)*w1
      do 25 l=l1,l2-2
  25  dvy(k,l,isp)=dvy(k,l,isp)+vxc(k,l+2,isp)*w2
      do 27 l=l1,l2-3
  27  dvy(k,l,isp)=dvy(k,l,isp)+vxc(k,l+3,isp)*w3
  20  continue
c  z-differentiation
C|    w1= 8d0/(12d0*dz)
C|    w2=-1d0/(12d0*dz)
      w1= 3d0/(4d0*dz)
      w2=-3d0/(20d0*dz)
      w3= 1d0/(60d0*dz)
      im3=ipmu(-3)
      im2=ipmu(-2)
      im1=ipmu(-1)
      ip1=ipmu(1)
      ip2=ipmu(2)
      ip3=ipmu(3)
      do 30 l=l1,l2
      do 30 k=k1,k2
      dvz(k,l,isp)=w3*(vxc2(k,l,ip3,isp)-vxc2(k,l,im3,isp))
     .        +w2*(vxc2(k,l,ip2,isp)-vxc2(k,l,im2,isp))
     .        +w1*(vxc2(k,l,ip1,isp)-vxc2(k,l,im1,isp))
  30  continue
c------end of the loop over spins----------------------------------
      enddo
c     call tcx('ymuder')
      end
