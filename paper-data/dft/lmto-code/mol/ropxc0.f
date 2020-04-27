      subroutine ropxc0(dx,dy,dz,pos,nbas,ips,rint,cenx,ceny,cenz,m1,m2)
C-  some setup for ropexc: center and z-limits
Cu Updates
Cu   21 Jul 07 (S. Lozovoi) New implementation of GGA
      implicit real*8 (a-h,p-z), integer (o)
      dimension rint(1),pos(3,1),ips(1)
c  get center of the global mesh
      cenx=0d0
      ceny=0d0
      cenz=0d0
      do 1 ib=1,nbas
      cenx=cenx+pos(1,ib)/nbas
      ceny=ceny+pos(2,ib)/nbas
  1   cenz=cenz+pos(3,ib)/nbas
c  get limits along z-direction
      m1=100
      m2=-100
      do 2 ib=1,nbas
      is=ips(ib)
      call roplmz(pos(1,ib),pos(2,ib),pos(3,ib),cenx,ceny,cenz,
     .   dx,dy,dz,rint(is),mx1,mx2,ier)
      if(ier == 0) then
        m1=min0(m1,mx1)
        m2=max0(m2,mx2)
        endif
  2   continue
      end
c --------- roplxy: max limits in x,y directions for some z -----
      subroutine roplxy(z,dx,dy,dz,cenx,ceny,pos,nbas,ips,rint,
     .    k1,k2,l1,l2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rint(1),ips(1),pos(3,1)
      l1=100
      l2=-100
      k1=100
      k2=-100
      do 3 ib=1,nbas
      is=ips(ib)
      zz=z-pos(3,ib)
      call roplmy(pos(1,ib),pos(2,ib),cenx,ceny,dx,dy,rint(is),zz,
     .   lx1,lx2,ier)
      if(ier == 0) then
        l1=min0(l1,lx1)
        l2=max0(l2,lx2)
        endif
      call roplmy(pos(2,ib),pos(1,ib),ceny,cenx,dy,dx,rint(is),zz,
     .   kx1,kx2,ier)
      if(ier == 0) then
        k1=min0(k1,kx1)
        k2=max0(k2,kx2)
        endif
  3   continue
      end
c --------- ropgth: gather/sort points within circle -------------
      subroutine ropgth(posx,posy,cenx,ceny,dx,dy,r,rc,zz,
     .   np,mp,x1,y1,x,y,z)
      implicit real*8 (a-h,p-z)
      dimension x(1),y(1),x1(1),y1(1),z(1)
      np=0
      mp=0
      call roplmy(posx,posy,cenx,ceny,dx,dy,r,zz,l1,l2,ier)
      if(ier /= 0) return
      do 10 l=l1,l2
        yy=l*dy+ceny-posy
        call roplmx(posx,cenx,dx,r,rc,zz,yy,k1,k2,j1,j2)
        do 11 k=k1,j1-1
          np=np+1
          y1(np)=yy
  11      x1(np)=k*dx+(cenx-posx)
        do 15 k=j1,j2
          mp=mp+1
          y(mp)=yy
  15      x(mp)=k*dx+(cenx-posx)
        do 12 k=j2+1,k2
          np=np+1
          y1(np)=yy
  12      x1(np)=k*dx+(cenx-posx)
  10    continue
      do 20 i=1,np
      x(i+mp)=x1(i)
  20  y(i+mp)=y1(i)
      np=np+mp
      do 30 i=1,np
  30  z(i)=zz
      end
c --------- rxcgth: gather points within circle -------------
      subroutine rxcgth(posx,posy,cenx,ceny,dx,dy,r,zz,np,x,y,z)
      implicit real*8 (a-h,p-z)
      dimension x(1),y(1),x1(1),y1(1),z(1)
      np=0
      call roplmy(posx,posy,cenx,ceny,dx,dy,r,zz,l1,l2,ier)
      if(ier /= 0) return
      do 10 l=l1,l2
        yy=l*dy+ceny-posy
        call rxclmx(posx,cenx,dx,r,zz,yy,k1,k2)
        do 11 k=k1,k2
          np=np+1
          y(np)=yy
  11      x(np)=k*dx+(cenx-posx)
  10  continue
      do 30 i=1,np
  30  z(i)=zz
      end
c --------- ropadd: add density from one atom to rhot -----
      subroutine ropadd(posx,posy,cenx,ceny,dx,dy,r,rc,zz,np,mp,rho,
     .    k1,k2,l1,l2,rhot)
      implicit real*8 (a-h,p-z)
      dimension rho(np),rhot(k1:k2,l1:l2)
      ip=mp
      jp=0
      call roplmy(posx,posy,cenx,ceny,dx,dy,r,zz,lx1,lx2,ier)
      if(ier /= 0) return
      do 10 l=lx1,lx2
        yy=l*dy+ceny-posy
        call roplmx(posx,cenx,dx,r,rc,zz,yy,kx1,kx2,jx1,jx2)
        do 11 k=kx1,jx1-1
          ip=ip+1
  11      rhot(k,l)=rhot(k,l)+rho(ip)
        do 15 k=jx1,jx2
          jp=jp+1
  15      rhot(k,l)=rhot(k,l)+rho(jp)
        do 12 k=jx2+1,kx2
          ip=ip+1
  12      rhot(k,l)=rhot(k,l)+rho(ip)
  10    continue
C|    write(6,762) ip,jp,mp,np
C|762 format(' end of loop:  ip,jp=',2i5,'   mp,np=',2i5)
      if(jp /= mp.or.ip /= np) call rx('inconsistency in ropadd')
      end
c --------- rxcadd: add density from one atom to rhot -----
      subroutine rxcadd(posx,posy,cenx,ceny,dx,dy,r,zz,np,rho,
     .    grho,g2rho,ggrho,
     .    k1,k2,l1,l2,rhot,grhot,g2rht,ggrhot)
c      implicit real*8 (a-h,p-z)
      implicit none
      integer np,k1,k2,l1,l2
      integer nsp,lsp,lxcg
      integer ip,l,isp,ia,k,ier,lx1,lx2,kx1,kx2
      double precision posx,posy,cenx,ceny,dx,dy,r,zz
      double precision rho(np,*),rhot(k1:k2,l1:l2,*)
      double precision grho(np,3,*),grhot(k1:k2,l1:l2,3,*)
      double precision g2rho(np,5,*),g2rht(k1:k2,l1:l2,5,*)
      double precision ggrho(np,*),ggrhot(k1:k2,l1:l2,*)
      double precision yy
      call tcn('rxcadd')
      nsp = lsp()+1
      ip=0
      call roplmy(posx,posy,cenx,ceny,dx,dy,r,zz,lx1,lx2,ier)
      if(ier /= 0) return
      do 10 l=lx1,lx2
        yy=l*dy+ceny-posy
        call rxclmx(posx,cenx,dx,r,zz,yy,kx1,kx2)
        do 11 k=kx1,kx2
          ip=ip+1
          do 12 isp=1,nsp
            rhot(k,l,isp)=rhot(k,l,isp)+rho(ip,isp)
c ---- if GGA, make also density derivatives -----
            if(lxcg() > 0) then
              ggrhot(k,l,isp)=ggrhot(k,l,isp)+ggrho(ip,isp)
              do ia = 1,3
                grhot(k,l,ia,isp)=grhot(k,l,ia,isp)+grho(ip,ia,isp)
              enddo
              do ia = 1,5
                g2rht(k,l,ia,isp)=g2rht(k,l,ia,isp)+g2rho(ip,ia,isp)
              enddo
            endif
  12      continue
  11    continue
  10  continue
      if(ip /= np) call rx('inconsistency in rxcadd')
      call tcx('rxcadd')
      end
c --------- ropvxc: add to sum1,sum2, overwrite rhot with vxc ----
      subroutine ropvxc(rhot,np,sum1,sum2)
      implicit real*8 (a-h,p-z)
      dimension rhot(np)
      do 10 i=1,np
  10  sum1=sum1+rhot(i)
      do 11 i=1,np
      ro=dmax1(rhot(i),0d0)
      call evxc(ro,0.5d0*ro,exc,vxc)
C|       exc=1d0
C|       vxc=1d0
      sum2=sum2+ro*exc
  11  rhot(i)=vxc
C|11  rhot(i)=exc
      end
c --------- ropgt1: gather/sort points and vxc -----------
      subroutine ropgt1(posx,posy,cenx,ceny,dx,dy,r,rc,zz,
     .    kd1,kd2,ld1,ld2,rhot,np,mp,x1,y1,v1,x,y,z,v)
      implicit real*8 (a-h,p-z)
      dimension x(1),y(1),x1(1),y1(1),v(1),v1(1),z(1),
     .   rhot(kd1:kd2,ld1:ld2)
      np=0
      mp=0
      call roplmy(posx,posy,cenx,ceny,dx,dy,r,zz,l1,l2,ier)
      if(ier /= 0) return
      do 10 l=l1,l2
        yy=l*dy+ceny-posy
        call roplmx(posx,cenx,dx,r,rc,zz,yy,k1,k2,j1,j2)
        do 11 k=k1,j1-1
          np=np+1
          v1(np)=rhot(k,l)
          y1(np)=yy
  11      x1(np)=k*dx+(cenx-posx)
        do 15 k=j1,j2
          mp=mp+1
          v(mp)=rhot(k,l)
          y(mp)=yy
  15      x(mp)=k*dx+(cenx-posx)
        do 12 k=j2+1,k2
          np=np+1
          v1(np)=rhot(k,l)
          y1(np)=yy
  12      x1(np)=k*dx+(cenx-posx)
  10    continue
      do 20 i=1,np
      v(i+mp)=v1(i)
      x(i+mp)=x1(i)
  20  y(i+mp)=y1(i)
      np=np+mp
      do 30 i=1,np
  30  z(i)=zz
      end
c --------- rxcgt1: gather points and vxc -----------
      subroutine rxcgt1(posx,posy,cenx,ceny,dx,dy,r,zz,
     .    kd1,kd2,ld1,ld2,rhot,np,x,y,z,v)
      implicit real*8 (a-h,p-z)
      dimension x(1),y(1),v(1),z(1),rhot(kd1:kd2,ld1:ld2)
      np=0
      call roplmy(posx,posy,cenx,ceny,dx,dy,r,zz,l1,l2,ier)
      if(ier /= 0) return
      do 10 l=l1,l2
        yy=l*dy+ceny-posy
        call rxclmx(posx,cenx,dx,r,zz,yy,k1,k2)
        do 11 k=k1,k2
          np=np+1
          v(np)=rhot(k,l)
          y(np)=yy
  11      x(np)=k*dx+(cenx-posx)
  10    continue
      do 30 i=1,np
  30  z(i)=zz
      end
c ------------ roplmz: limits along z-direction --------
      subroutine roplmz(posx,posy,posz,cenx,ceny,cenz,dx,dy,dz,r,
     .   m1,m2,ier)
      implicit real*8 (a-h,p-z)
      k0=idnint((posx-cenx)/dx)
      l0=idnint((posy-ceny)/dy)
      xx0=k0*dx+cenx-posx
      yy0=l0*dy+ceny-posy
      w2=r*r-xx0*xx0-yy0*yy0
      ier=1
      if(w2 >= 0d0) then
        w=dsqrt(w2)
        m1=idnint((posz-cenz-w)/dz+0.5d0)
        m2=idnint((posz-cenz+w)/dz-0.5d0)
        ier=0
        endif
      end
c ------------ roplmy: limits along y-direction --------
      subroutine roplmy(posx,posy,cenx,ceny,dx,dy,r,zz,l1,l2,ier)
      implicit real*8 (a-h,p-z)
      k0=idnint((posx-cenx)/dx)
      xx0=k0*dx+cenx-posx
      w2=r*r-zz*zz-xx0*xx0
      ier=1
      if(w2 >= 0d0) then
        w=dsqrt(w2)
        l1=idnint((posy-ceny-w)/dy+0.5d0)
        l2=idnint((posy-ceny+w)/dy-0.5d0)
        ier=0
        endif
      end
c ------------ roplmx: limits along x-direction --------
      subroutine roplmx(posx,cenx,dx,r,rc,zz,yy,k1,k2,j1,j2)
      implicit real*8 (a-h,p-z)
      k1=1
      k2=0
      j1=k2+1
      j2=k2
      w2=r*r-zz*zz-yy*yy
      if(w2 <= 0d0) return
      w=dsqrt(w2)
      k1=idnint((posx-cenx-w)/dx+0.5d0)
      k2=idnint((posx-cenx+w)/dx-0.5d0)
      j1=k2+1
      j2=k2
      w2=rc*rc-zz*zz-yy*yy
      if(w2 > 0d0) then
        w=dsqrt(w2)
        i1=idnint((posx-cenx-w)/dx+0.5d0)
        i2=idnint((posx-cenx+w)/dx-0.5d0)
        if(i2 >= i1) then
          j1=i1
          j2=i2
          endif
        endif
      end
c ------------ rxclmx: limits along x-direction --------
      subroutine rxclmx(posx,cenx,dx,r,zz,yy,k1,k2)
      implicit real*8 (a-h,p-z)
      k1=1
      k2=0
      w2=r*r-zz*zz-yy*yy
      if(w2 <= 0d0) return
      w=dsqrt(w2)
      k1=idnint((posx-cenx-w)/dx+0.5d0)
      k2=idnint((posx-cenx+w)/dx-0.5d0)
      end
