      subroutine ropexc(nbas,nspec,pos,rsm,rcut,rint,radd,nxi,lxi,
     .   exi,n0,ips,dx,dy,dz,atr,nc0,nc1,nc2,c1,c2,ioff,cy,nri,rhoi,
     .   rhoep,rhomu,vxci)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rsm(1),rcut(1),rint(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   rhoi(nri),vxci(nri),pos(3,1),ioff(1),cy(1),ips(1),
     .   nc1(1),nc2(1),c1(1),c2(1)
      real w(1)
      common /w/ w
c ------ get center of the global mesh and z limits -------
      call ropxc0(dx,dy,dz,pos,nbas,ips,rint,cenx,ceny,cenz,m1,m2)
      if(iprint() >= 20) write(6,903) cenx,ceny,cenz,m1,m2,dx,dy,dz
  903 format(/' ropexc:  center=  ',3f9.3,'     m1,m2=',2i5
     .       /'          dx,dy,dz=',3f9.3)
      if(m1 > m2) call rx('ropexc: empty mesh')
      sum1=0d0
      sum2=0d0
      call dpzero(vxci,   nri)
      kb= 1000
      kt=-1000
      lb= 1000
      lt=-1000

c ------ start big loop over z-planes ---------
      do 80 m=m1,m2
      z=m*dz+cenz
      call roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,k1,k2,l1,l2)
      nx=k2-k1+1
      ny=l2-l1+1
      if(nx <= 0.or.ny <= 0) goto 80
      kb=min0(kb,k1)
      kt=max0(kt,k2)
      lb=min0(lb,l1)
      lt=max0(lt,l2)

c --------- make total density in this plane ----------
      call defrr(orhot,     nx*ny)
      call dpzero(w(orhot), nx*ny)
      call ropmro(dx,dy,dz,pos,nbas,rsm,rcut,rint,radd,nxi,lxi,exi,n0,
     .   atr,nc0,nc1,c1,nc2,c2,ips,ioff,cy,rhoi,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,w(orhot))
c --------- this part to check ropmro -----------
C|    call ropcro(dx,dy,dz,pos,nbas,rsm,rcut,rint,nxi,lxi,exi,n0,atr,
C|   .   nc0,nc1,c1,nc2,c2,ips,ioff,cy,rhoi,cenx,ceny,cenz,
C|   .   z,k1,k2,l1,l2,w(orhot))

c --------- add to exc, overwrite rhot with vxc ------
C|    call ropvxc(w(orhot),nx*ny,sum1,sum2)

      call rx('update call to rhopxcc')
C      call defrr(oh,    nx*ny)
C      call defrr(ou2,   nx*ny)
C      call ropxcc(nx*ny,w(orhot),w(oh),w(ou2),sum1,sum2)
C      call rlse(oh)

c --------- make integrals over vxc -------------------
      call ropvro(dx,dy,dz,pos,nbas,rsm,rcut,rint,radd,nxi,lxi,exi,n0,
     .   atr,nc0,nc1,c1,nc2,c2,ips,ioff,cy,vxci,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,w(orhot))

      call rlse(orhot)
  80  continue
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
      if(iprint() >= 20) write(6,732) kb,kt,lb,lt,m1,m2,q,rhoep,rhomu
  732 format(' max extent of mesh:  (x)',2i5,'   (y)',2i5,'   (z)',2i5
     .   /' q=',f12.6,'   rhoep=',f12.6,'   rhomu=',f12.6)

      end
