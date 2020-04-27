      subroutine ropy0m(n,x,z,lmax,nd,yl,rsq)
C- Makes unnormalized spheric harmonic polynomials (vectorizes).
C  Position (l*(l+1))/2+m+1 holds element ylm(l,m) (loop
C  order l=0,lmax; m=0,l).
C  Also returns squares of points in rsq.
      implicit real*8 (a-h,p-z), integer (o)
      dimension x(1),z(1),yl(nd,1),rsq(1)
      real w(1)
      common /w/ w
      call defrr(ocm,  n)
      call defrr(osm,  n)
      call defrr(oq,   n*2)
      call defrr(oh,   n)
      call defrr(oy,   n)
      call dpzero(w(oy),n)
      do 2 i=1,n
    2 rsq(i)=x(i)*x(i)+z(i)*z(i)
c ------- start big loop over m ------------
      do 10 m=0,lmax
      call ropcsm(m,n,x,w(oy),w(oh),w(ocm),w(osm))
      do 11 l=m,lmax
        lm = (l*(l+1))/2+m+1
        call ropqlm(m,l,n,rsq,z,w(oq),kk)
        call ropy0x(lm,kk,n,w(oq),w(ocm),nd,yl)
   11 continue
   10 continue
      call rlse(ocm)
      end
c --------- sub ropyxx --------
      subroutine ropy0x(lm,kk,n,q,cm,nd,yl)
      implicit real*8 (a-h,p-z), integer (o)
      dimension q(n,2),cm(n),yl(nd,1)
      do 1 i=1,n
    1 yl(i,lm)=cm(i)*q(i,kk)
      end
