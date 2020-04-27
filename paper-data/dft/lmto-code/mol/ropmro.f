      subroutine ropmro(dx,dy,dz,pos,nbas,rsm,rcut,rint,radd,nxi,lxi,
     .   exi,n0,atr,nc0,nc1,c1,nc2,c2,ips,ioff,cy,rhoi,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,rhot)
c  loops over all atoms, makes density in one xy-plane
      implicit real*8 (a-h,p-z), integer (o)
      dimension rsm(1),rcut(1),rint(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   rhoi(1),pos(3,1),ipd(1),ioff(1),cy(1),ips(1),rhot(1),
     .   c2(nc0,2,n0,1),nc2(1),c1(nc0,n0,n0,1),nc1(1)
      real w(1)
      common /w/ w
      msum=0d0
      nsum=0d0
      nat=0
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
c --------- gather x,y coordinates --------------
      call defrr(ox,      npmx)
      call defrr(oy,      npmx)
      call defrr(oz,      npmx)
      call defrr(ox1,     npmx)
      call defrr(oy1,     npmx)
      call ropgth(posx,posy,cenx,ceny,dx,dy,ri,rc,zz,
     .   np,mp,w(ox1),w(oy1),w(ox),w(oy),w(oz))
      if(np == 0) goto 10
      call rlse(ox1)
      nsum=nsum+np
      msum=msum+mp
      nat=nat+1
      if(np > npmx) call rx('ropmro: np gt npmx')
c --------- make and add density in xy-plane --------
      call defrr(orho,    npmx)
      call defrr(or2,     np)
      call defrr(oq,      np)
      ir=ioff(ib)+1
      call roprho(nxi(is),lxi(1,is),exi(1,is),rsm(is),rcut(is),rint(is),
     .   radd,atr,rhoi(ir),cy,n0,nc0,nc1(is),c1(1,1,1,is),nc2(is),
     .   c2(1,1,1,is),w(ox),w(oy),w(oz),w(orho),w(or2),w(oq),np,mp)
      call ropadd(posx,posy,cenx,ceny,dx,dy,ri,rc,zz,np,mp,w(orho),
     .    k1,k2,l1,l2,rhot)
      call rlse(ox)
  10  continue
      if(iprint() >= 50) write(6,817) z,k1,k2,l1,l2,nat,msum,nsum
  817 format('   z=',f9.3,'   x:',2i5,'    y:',2i5,'   nat',i3,2i8)
      end
