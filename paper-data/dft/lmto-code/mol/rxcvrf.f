      subroutine rxcvrf(dx,dy,dz,pos,nbas,rsm,rint,nxi,lxi,exi,n0,
     .   ips,ioff,cy,vxci,cenx,ceny,cenz,
     .   z,k1,k2,l1,l2,vxc,dvx,dvy,dvz,nri,rhoi,f)
c  loops over all atoms, makes integrals vxc*h in one xy plane
      implicit real*8 (a-h,p-z), integer (o)
      dimension rsm(1),rint(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   vxci(nri,*),pos(3,1),ipd(1),ioff(1),cy(1),ips(1),
     .   rhoi(nri,*),f(3,1),vxc(k1:k2,l1:l2,*),
     .   dvx(k1:k2,l1:l2,*),dvy(k1:k2,l1:l2,*),dvz(k1:k2,l1:l2,*)

      real w(1)
      common /w/ w
      call tcn('rxcvrf')
      nsp=lsp()+1
      xxx=dx*dy*dz
c --------- start loop over atoms ---------------
      do 10 ib=1,nbas
      is=ips(ib)
      posx=pos(1,ib)
      posy=pos(2,ib)
      zz=z-pos(3,ib)
      ri=rint(is)
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
      do isp=1,nsp
C******** debug print **************
c            write(6,*) ' ************ rxcvrf ****'
c            write(6,*) ' nbas,nxi,n0,nk =',nbas,nxi,n0,nk
c            write(6,*) ' k1,k2,l1,l2 =',k1,k2,l1,l2
c            write(6,*) ' vxc ='
c          do  i = 1, nsp
c            write(6,*) ' isp =',i
c            do il=l1,l2,10
c               write(6,'(i3,100g12.4)') il,(vxc(ik,il,i),ik=k1,k2,10)
c            enddo
c          enddo
C******** debug print **************
        call rxcgt1(posx,posy,cenx,ceny,dx,dy,ri,zz,k1,k2,l1,l2,
     .   vxc(k1,l1,isp),np,w(ox),w(oy),w(oz),w(ovxc))
        call rxcgt1(posx,posy,cenx,ceny,dx,dy,ri,zz,k1,k2,l1,l2,
     .   dvx(k1,l1,isp),np,w(ox),w(oy),w(oz),w(odvx))
        call rxcgt1(posx,posy,cenx,ceny,dx,dy,ri,zz,k1,k2,l1,l2,
     .   dvy(k1,l1,isp),np,w(ox),w(oy),w(oz),w(odvy))
        call rxcgt1(posx,posy,cenx,ceny,dx,dy,ri,zz,k1,k2,l1,l2,
     .   dvz(k1,l1,isp),np,w(ox),w(oy),w(oz),w(odvz))
        if(np == 0) goto 10
        if(np > npmx) call rx('rxcvrf: np gt npmx')
c --------- make and add density in xy-plane --------
        call defrr(or2,     np)
        call defrr(oq,      np*2)
        ir=ioff(ib)+1
        call rxcvhf(nxi(is),lxi(1,is),exi(1,is),rsm(is),
     .     vxci(ir,isp),cy,n0,w(ox),w(oy),w(oz),w(ovxc),w(or2),w(oq),
     .     np,w(odvx),w(odvy),w(odvz),rhoi(ir,isp),xxx,f(1,ib))
      enddo
      call rlse(ox)
  10  continue
      call tcx('rxcvrf')
      end
