      subroutine ftcsrt(mmax,nx1,lx1,ex1,nx2,lx2,ex2,lp1,lp2,
     .   rmt1,rmt2,lscal,lsym,ndim,nrhs,ba,bb,ndimx,nrhsx,ndimbx,
     .   nalf,cofa,ncof)
c  sort into columns by ialf
      implicit real*8 (a-h,p-z), integer (o)
      dimension ndim(0:20),nrhs(0:20),cofa(ncof,nalf),lx1(1),
     .   lx2(1),ex1(1),ex2(1),
     .   ba(ndimx,nrhsx,0:mmax,nalf),bb(ndimbx,nrhsx,0:mmax)
      call getpr(ipr)
c ------- sort coeffs into ba by ialf ----------------
      do 10 m=0,mmax
      nd=ndim(m)
      nr=nrhs(m)
      do 11 ialf=1,nalf
      ia0=(ialf-1)*nd
      do 11 ir=1,nr
      do 11 id=1,nd
  11  ba(id,ir,m,ialf)=bb(id+ia0,ir,m)
  10  continue
c ------- undo xi-scaling ---------------------
      if(lscal == 1) then
      if(ipr >= 30) write(6,*) 'ftcsrt:  undo xi-scaling'
      call pshpr(1)
      do 30 ialf=1,nalf
      call ftcusc(mmax,ndim,nrhs,nx1,lx1,ex1,nx2,lx2,ex2,
     .   rmt1,rmt2,ba(1,1,0,ialf),ndimx,nrhsx)
  30  continue
      call poppr
      endif

c ------- call tcfsrt for each ialf -----------
      do 20 ialf=1,nalf
      call tcfsrt(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,ndim,nrhs,
     .    ba(1,1,0,ialf),ndimx,nrhsx,cofa(1,ialf))
  20  continue

      return
      end
