      subroutine ropvr7(dx,dy,dz,pos,nbas,rsm,rcut,rint,radd,nxi,lxi,
     .   exi,n0,atr,nc0,nc1,c1,nc2,c2,ips,ioff,cy,vxci,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,vxc,dvx,dvy,dvz,dvxcx,dvxcy,dvxcz)
c  loops over all atoms, makes integrals vxc*h in one xy plane
      implicit real*8 (a-h,p-z), integer (o)
      dimension rsm(1),rcut(1),rint(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   vxci(1),pos(3,1),ipd(1),ioff(1),cy(1),ips(1),vxc(1),
     .   c2(nc0,2,n0,1),nc2(1),c1(nc0,n0,n0,1),nc1(1),
     .   dvx(1),dvy(1),dvz(1),dvxcx(1),dvxcy(1),dvxcz(1)
      real w(1)
      common /w/ w
c --------- start loop over atoms ---------------
      do 10 ib=1,nbas
      is=ips(ib)
      posx=pos(1,ib)
      posy=pos(2,ib)
      zz=z-pos(3,ib)
      ri=rint(is)
      rc=rcut(is)
      npmx=4d0*(ri*ri-zz*zz)/(dx*dy)
      npmx=max0(npmx,40)
c --------- gather x,y coordinates and vxc ---------
      call defrr(ox,      npmx)
      call defrr(oy,      npmx)
      call defrr(oz,      npmx)
      call defrr(ovxc,    npmx)
      call defrr(odvx,    npmx)
      call defrr(odvy,    npmx)
      call defrr(odvz,    npmx)
      call defrr(ox1,     npmx)
      call defrr(oy1,     npmx)
      call defrr(ov1,     npmx)
      call ropgt1(posx,posy,cenx,ceny,dx,dy,ri,rc,zz,k1,k2,l1,l2,
     .   vxc,np,mp,w(ox1),w(oy1),w(ov1),w(ox),w(oy),w(oz),w(ovxc))
      call ropgt1(posx,posy,cenx,ceny,dx,dy,ri,rc,zz,k1,k2,l1,l2,
     .   dvx,np,mp,w(ox1),w(oy1),w(ov1),w(ox),w(oy),w(oz),w(odvx))
      call ropgt1(posx,posy,cenx,ceny,dx,dy,ri,rc,zz,k1,k2,l1,l2,
     .   dvy,np,mp,w(ox1),w(oy1),w(ov1),w(ox),w(oy),w(oz),w(odvy))
      call ropgt1(posx,posy,cenx,ceny,dx,dy,ri,rc,zz,k1,k2,l1,l2,
     .   dvz,np,mp,w(ox1),w(oy1),w(ov1),w(ox),w(oy),w(oz),w(odvz))
      if(np == 0) goto 10
      call rlse(ox1)
      if(np > npmx) call rx('ropvro: np gt npmx')
c --------- make and add density in xy-plane --------
      call defrr(or2,     np)
      call defrr(oq,      np)
      ir=ioff(ib)+1
      call ropvh7(nxi(is),lxi(1,is),exi(1,is),rsm(is),rcut(is),rint(is),
     .   radd,atr,vxci(ir),cy,n0,nc0,nc1(is),c1(1,1,1,is),nc2(is),
     .   c2(1,1,1,is),w(ox),w(oy),w(oz),w(ovxc),w(or2),w(oq),np,mp,
     .   w(odvx),w(odvy),w(odvz),dvxcx(ir),dvxcy(ir),dvxcz(ir))
      call rlse(ox)
  10  continue
      end
