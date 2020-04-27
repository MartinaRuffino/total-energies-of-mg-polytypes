      subroutine ropylm(n,x,y,z,lmax,nd,yl,rsq)
C- Makes unnormalized spheric harmonic polynomials (vectorizes).
C- Also returns squares of points in rsq.
      implicit real*8 (a-h,p-z), integer (o)
      dimension x(1),y(1),z(1),yl(nd,1),rsq(1)
      real w(1)
      common /w/ w
      call defrr(ocm,  n)
      call defrr(osm,  n)
      call defrr(oq,   n*2)
      call defrr(oh,   n)
      do 2 i=1,n
  2   rsq(i)=x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
c ------- start big loop over m ------------
      do 10 m=0,lmax
      call ropcsm(m,n,x,y,w(oh),w(ocm),w(osm))
      do 11 l=m,lmax
        call ropqlm(m,l,n,rsq,z,w(oq),kk)
        call ropylx(m,l,kk,n,w(oq),w(ocm),w(osm),nd,yl)
  11    continue
  10  continue
      call rlse(ocm)
      end
c --------- sub ropyxx --------
      subroutine ropylx(m,l,kk,n,q,cm,sm,nd,yl)
      implicit real*8 (a-h,p-z), integer (o)
      dimension q(n,2),cm(n),sm(n),yl(nd,1)
      lav=l*(l+1)+1
      do 1 i=1,n
  1   yl(i,lav+m)=cm(i)*q(i,kk)
      if(m == 0) return
      do 2 i=1,n
  2   yl(i,lav-m)=sm(i)*q(i,kk)
      end
