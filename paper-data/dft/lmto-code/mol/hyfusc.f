      subroutine hyfusc(mmax,ndim,nrhs,nx1,lx1,ex1,nx2,lx2,ex2,
     .  nalf,rmt1,rmt2,b,ndimx,nrhsx)
C- HYF analog of ftcusc
      implicit real*8 (a-h,p-z), integer (o)
      dimension nrhs(0:20),ndim(0:20),ex1(1),lx1(1),ex2(1),lx2(1),
     .   b(ndimx,nrhsx,0:mmax,nalf),fac(0:20)
      call getpr(ipr)
      if(ipr >= 70) write(6,200)
  200 format(/' ftcusc:  undo xi scaling')

C ------- start loop over m --------
      do 30 m=0,mmax
      ii=0

c ... first part
      do 31 ie=1,nx1
      call hrmsav(rmt1,lx1(ie),ex1(ie),fac)
      if(ipr >= 80) write(6,300) ex1(ie),rmt1,(fac(l),l=0,lx1(ie))
  300 format(' e',f9.3,'   r',f9.4,'   fac',7f9.4)
      do 31 l=m,lx1(ie)
      ii=ii+1
      xx=1d0/fac(l)
      do 36 ialf=1,nalf
      do 36 jrhs=1,nrhs(m)
  36  b(ii,jrhs,m,ialf)=b(ii,jrhs,m,ialf)*xx
  31  continue

c ... second part
      do 32 ie=1,nx2
      call hrmsav(rmt2,lx2(ie),ex2(ie),fac)
      if(ipr >= 80) write(6,300) ex2(ie),rmt2,(fac(l),l=0,lx2(ie))
      do 32 l=m,lx2(ie)
      ii=ii+1
      xx=1d0/fac(l)
      do 37 ialf=1,nalf
      do 37 jrhs=1,nrhs(m)
  37  b(ii,jrhs,m,ialf)=b(ii,jrhs,m,ialf)*xx
  32  continue

  30  continue

      end
