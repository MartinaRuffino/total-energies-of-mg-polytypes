      subroutine ylmrt0(lmax,nlm,xmat,p,cy)
c  Initializes xmat and p for use in ylmrot.
      parameter( lmx=6,  nmmx=2*lmx+1 )
      implicit real*8 (a-h,p-z), integer(o)
      dimension p(3,1),xmat(nlm,nlm),yl((lmx+1)**2),cy(1),
     .   ipiv(nmmx),a(nmmx,nmmx),p0(3,nmmx),b(nmmx,nmmx)
      data p0/ 0.020,-.025,-.118,0.419,-.538,0.513,0.245,-.717,-.600,
     .         -.056,0.224,-.309,-.034,-.180,0.207,-.351,-.614,0.950,
     .         -.782,-.134,-.308,0.568,0.716,-.457,-.528,-.927,-.562,
     .         -.856,-.443,0.267,-.111,0.794,0.598,-.985,-.144,-.617,
     .         0.678,0.400,-.617 /
      nmmax=2*lmax+1
      if(lmax > lmx) call rx('ylmrt0:  lmax gt lmx')
      do 30 ip=1,nmmax
      fc=1.d0/dsqrt(p0(1,ip)**2+p0(2,ip)**2+p0(3,ip)**2)
      p(1,ip)=p0(1,ip)*fc
      p(2,ip)=p0(2,ip)*fc
  30  p(3,ip)=p0(3,ip)*fc
      do 29 j=1,nlm
      do 29 i=1,nlm
  29  xmat(i,j)=0.d0
c -------------------------------
      do 10 l=0,lmax
      nm=2*l+1
      ilm0=l*l
      do 1 i=1,nm
      call sylm(p(1,i),yl,l,p2)
      do 1 m=1,nm
  1   a(m,i)=cy(m+ilm0)*yl(m+ilm0)
      call leqr2f(a,nm,nmmx,ipiv,yl)
      do 2 j=1,nm
      do 3 i=1,nm
  3   b(i,j)=0.d0
  2   b(j,j)=1.d0
      call leqr2s(a,nm,nmmx,b,nm,ipiv)
      do 4 j=1,nm
      do 4 i=1,nm
  4   xmat(i+ilm0,j+ilm0)=b(i,j)
  10  continue
      return
      end
c -------------------------------------------------
      subroutine ylmrt1(lmax,nlm,rotpar,rmat,xmat,p,cy)
c  makes matrix which transform the sph. harmonics under rotation.
c  def:  ylm(m1,rot(p))= sum rmat(m1,m2)*ylm(m2,p)   for each l
c  or:   ylm(l,rot-1(p))=sum rmat(m,l)*ylm(m,p)    (is equivalent).
c  calls a subroutine rotpnt(p1,p2,rotpar) which makes p2=rot(p1)
      parameter( lmx=6,  nmmx=2*lmx+1, lx=10)
      implicit real*8 (a-h,p-z), integer(o)
      dimension p(3,1),xmat(nlm,nlm),yl((lx+1)**2),cy(1),
     .  a(nmmx,nmmx),rmat(nlm,nlm),pp(3),rotpar(1)
      if (lmax > lx) call rx('ylmrt1: lmax gt lx')
      do 1 j=1,nlm
      do 1 i=1,nlm
  1   rmat(i,j)=0.d0
      nmmax=2*lmax+1
c ----- make matrix with transformed yl-values ----
      do 20 i=1,nmmax
      call rotpnt(p(1,i),pp,rotpar)
      call sylm(pp,yl,lmax,p2)
      do 20 l=0,lmax
      ilm0=l*l
      nm=2*l+1
      do 20 m=1,nm
  20  if(i <= nm) rmat(m+ilm0,i+ilm0)=cy(m+ilm0)*yl(m+ilm0)
c ----- multiply by xmat ---------
      do 21 l=0,lmax
      nm=2*l+1
      ilm0=l*l
      do 22 j=1,nm
      do 22 i=1,nm
  22  a(i,j)=rmat(i+ilm0,j+ilm0)
      do 23 j=1,nm
      do 23 i=1,nm
      sum=0.d0
      do 24 k=1,nm
  24  sum=sum+a(i,k)*xmat(k+ilm0,j+ilm0)
  23  rmat(i+ilm0,j+ilm0)=sum
  21  continue
c -----------------------------------------
c|    do 40 l=0,lmax
c|    ilm0=l*l
c|    nm=2*l+1
c|    do 40 i=1,nm
c|40  write(6,400) (rmat(i+ilm0,j+ilm0),j=1,nm)
c|400 format(1x,9f8.5)
      return
      end
