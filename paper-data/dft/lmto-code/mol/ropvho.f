      subroutine ropvho(nxi,lxi,exi,rsm,rcut,rint,radd,atr,vxci,cy,
     .    n0,nc0,nc1,c1,nc2,c2,x,y,z,v,rsq,q,n,mp)
c  for one atom and xy-plane, adds to vxc integrals
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),exi(1),cy(1),vxci(1),x(1),y(1),z(1),v(1),
     .   c2(nc0,2,nxi),q(n,2),rsq(n),ioff(10),c1(nc0,n0,nxi)
      real w(1)
      common /w/ w
      lmax=-1
      ioff(1)=0
      do 1 ie=1,nxi
      ioff(ie+1)=ioff(ie)+(lxi(ie)+1)**2
  1   lmax=max0(lmax,lxi(ie))
c ------- initialize rsq -----------------
      do 2 i=1,n
  2   rsq(i)=x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
c ------- make all radial parts xi -------
      call defrr(oxi,   n*(lmax+1)*nxi)
      call defrr(oee,   n)
      call rophsm(nxi,lxi,exi,rsm,lmax,rcut,rint+radd,atr,
     .   n0,nc0,nc1,c1,nc2,c2,rsq,w(oee),w(oxi),n,mp)
C|    call rophsd(nxi,lxi,exi,rsm,lmax,rsq,w(oxi),n)
      call rlse(oee)
c ------- dimension arrays -----------------
      call defrr(of,   n)
      call defrr(og,   n)
      call defrr(ocm,  n)
      call defrr(osm,  n)
      call defrr(oh,   n)
c ------- start big loop over m ------------
      do 10 m=0,lmax
      call ropcsm(m,n,x,y,w(oh),w(ocm),w(osm))
      call dpzero(w(of),  n)
      call dpzero(w(og),  n)
      do 11 l=m,lmax
        call ropqlm(m,l,n,rsq,z,q,kk)
        call ropxfv(m,l,nxi,lxi,vxci,cy,lmax,ioff,w(oxi),
     .   w(ocm),w(osm),v,q(1,kk),w(oh),w(of),w(og),n)
  11    continue
  10  continue
      call rlse(oxi)
      end