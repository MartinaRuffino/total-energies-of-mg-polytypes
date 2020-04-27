      subroutine rxcvhf(nxi,lxi,exi,rsm,vxci,cy,
     .   n0,x,y,z,v,rsq,q,n,dvx,dvy,dvz,rhoi,xxx,f)
c  for one atom and xy-plane, adds to vxc integrals
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),exi(1),cy(1),vxci(1),x(1),y(1),z(1),v(1),
     .   q(n,2),rsq(n),ioff(10),dvx(1),dvy(1),dvz(1),rhoi(1),f(3)
      real w(1)
      common /w/ w

      call tcn('rxcvhf')
      lmax=-1
      ioff(1)=0
      do 1 ie=1,nxi
      ioff(ie+1)=ioff(ie)+(lxi(ie)+1)**2
  1   lmax=max0(lmax,lxi(ie))
c ------- initialize rsq -----------------
      do 2 i=1,n
  2   rsq(i)=x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
c ... deb
c           print *,' ****** rxcvhf *******'
c           write(6,*) ' nxi,n,n0 =',nxi,n,n0
c            write(6,*) ' rhoi ='
c            do ip=1,nlm
c               write(6,'(100g12.4)') (rhoi(ir),ir=1,n,10)
c            enddo
c ... deb
c ------- make all radial parts xi -------
      call defrr(oxi,   n*(lmax+1)*nxi)
      call defrr(owk,   (4+lmax)*n)
      call defrr(oidx,  2*n)
      call hansr(rsm,0,lmax,nxi,lxi,exi,rsq,n,n,w(oidx),w(owk),
     .  000,w(oxi))
      call rlse(owk)
c ------- dimension arrays -----------------
      call defrr(of,   n)
      call defrr(ogcm, n)
      call defrr(ogsm, n)
      call defrr(ocm,  n)
      call defrr(osm,  n)
      call defrr(oh,   n)
c ------- start big loop over m ------------
      do 10 m=0,lmax
      call ropcsm(m,n,x,y,w(oh),w(ocm),w(osm))
      call dpzero(w(ogcm),   n)
      call dpzero(w(ogsm),   n)
      do 11 l=m,lmax
        call ropqlm(m,l,n,rsq,z,q,kk)
        call ropxvf(m,l,nxi,lxi,vxci,cy,lmax,ioff,w(oxi),w(ocm),
     .   v,q(1,kk),w(oh),w(of),w(ogcm),n,rhoi,001)
        call ropxvf(m,l,nxi,lxi,vxci,cy,lmax,ioff,w(oxi),w(osm),
     .   v,q(1,kk),w(oh),w(of),w(ogsm),n,rhoi,002)
  11    continue
c  add to forces ....
      call rxxvhf(n,w(ocm),w(osm),w(ogcm),w(ogsm),w(oh),
     .   dvx,dvy,dvz,xxx,f)
  10  continue
      call rlse(oxi)

      call tcx('rxcvhf')
      end
c ------- sub rxxvhf: add to forces --------------
      subroutine rxxvhf(n,cm,sm,gcm,gsm,h,dvx,dvy,dvz,xxx,f)
      implicit real*8 (a-h,p-z), integer (o)
      dimension cm(n),sm(n),gcm(n),gsm(n),h(n),dvx(n),dvy(n),
     .  dvz(n),f(3)
      do 10 i=1,n
  10  h(i)=gcm(i)*cm(i)+gsm(i)*sm(i)
      sx=0d0
      sy=0d0
      sz=0d0
      do 11 i=1,n
      sx=sx+h(i)*dvx(i)
      sy=sy+h(i)*dvy(i)
  11  sz=sz+h(i)*dvz(i)
      f(1)=f(1)-xxx*sx
      f(2)=f(2)-xxx*sy
      f(3)=f(3)-xxx*sz
      end
