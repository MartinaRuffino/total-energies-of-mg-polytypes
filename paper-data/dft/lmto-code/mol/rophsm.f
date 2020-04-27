      subroutine rophsm(nxi,lxi,exi,rsm,lmax,r1,r2,atr,
     .   n0,nc0,nc1,c1,nc2,c2,rsq,ee,xi,n,m)
c  Makes vectors of smoothed Hankel functions, l=0..lx.
c  n=number of points, m=number of points within rcut.
c  List must be sorted so that points 1..m come first.
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(nxi),exi(nxi),c2(nc0,2,nxi),
     .   rsq(n),xi(n,0:lmax,nxi),ee(n),c1(nc0,n0,nxi)
      real w(1)
      common /w/ w
      k=n-m
      r0=0d0
      a2=1d0/(rsm*rsm)
      do 1 i=m+1,n
  1   ee(i)=dexp(-a2*rsq(i))
      call defrr(oww,    3*n)
c ------ make transformed variables -------
      call defrr(ou1,    m)
      call defrr(ou2,    k)
      call roptra(r0,r1,atr,m,rsq,     w(ou1))
      call roptra(r1,r2,atr,k,rsq(m+1),w(ou2))
c ------ start loop over energies ---------
      do 50 ie=1,nxi
      lx=lxi(ie)
      if (lx > lmax) stop 'in rophsm'
      e=exi(ie)
      call rophs1(lx,c1(1,1,ie),nc0,nc1,w(oww),w(ou1),
     .   xi(1,0,ie),n,m)
      call rophs2(e,rsm,lx,c2(1,1,ie),nc0,nc2,rsq(m+1),
     .   ee(m+1),w(oww),w(ou2),xi(m+1,0,ie),n,k)
  50  continue
      call rlse(oww)
      end
