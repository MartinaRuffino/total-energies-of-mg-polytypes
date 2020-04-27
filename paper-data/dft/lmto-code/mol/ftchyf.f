      subroutine ftchyf(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,d0,ndist,adec,nalf,
     .  nb1,nb2,ri1,ri2,ifi)
C- Global fit to all distances simultaneously
C  lscal=1: scale to improve condition of overlap matrix
      implicit real*8 (a-h,p-z), integer (o)
      dimension lx1(8),lx2(8),ex1(8),ex2(8),dist(100),ndim(0:20),
     .   nrhs(0:20),ndimb(0:20),f(2000),wist(100),err1(100),err2(100)
      character*1 cl(0:8)
      data cl /'s','p','d','f','g','5','6','7','8'/
      real w(1)
      common /w/ w
      lscal=1
      call getpr(ipr)
      if(ipr >= 20) write(6,330) ndist,nb1,nb2,lscal
  330 format(/' ftchyf:  ndist=',i3,'   nb1,nb2=',2i4,'    lsc=',
     .   i1)
c ----- dimensions ----------
      call ftcdim(lx1,ex1,nx1,lx2,ex2,nx2,mmax,lmax,ep1,lp1,
     .   ep2,lp2,rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,lsym,ndim,nrhs,
     .   ndimx,nsdmx,nrhsx,ncof,nxi,000)
      nsdmbx=0
      ndimbx=0
      do 3 m=0,mmax
      ndimb(m)=ndim(m)*nalf
      nsdb=(ndimb(m)*(ndimb(m)+1))/2
      ndimbx=max0(ndimbx,ndimb(m))
  3   nsdmbx=max0(nsdmbx,nsdb)
      if(ipr >= 30) write(6,300) nsdmbx,(ndimb(m),m=0,mmax)
  300 format(' nsdmbx=',i6,'  ndimb=',12i5)
      call defrr(obb,    ndimbx*nrhsx*(mmax+1))
      call defrr(osb,    nsdmbx*(mmax+1))

c ----- set distances and f -----------
      if(ndist > 100)       call rx('ftchyf: ndist too big')
      if(ndist*nalf > 2000) call rx('ftchyf: increase dim of f')
      call ftcsch(ndist,nalf,d0,adec,dist,wist,f)
CL      write(71,710) ndist,adec,d0,dist(ndist),nalf,lscal
  710 format(' ndist',i3,'   a',f7.3,'   d0',f8.4,'   dn',f8.4,
     .   '   nalf',i3,'   lsc=',i1)

C ----- make matrices bb and sb for big fit problem ----
      call defrr(ob,   ndimx*nrhsx*(mmax+1))
      call defrr(os,   nsdmx*(mmax+1))
      call ftchad(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .   rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,nalf,ndist,dist,wist,f,
     .   nb1,nb2,ri1,ri2,mmax,lmax,lsym,ndim,nrhs,ndimx,nrhsx,nsdmx,
     .   ndimbx,nsdmbx,w(ob),w(os),w(obb),w(osb),lscal)
      call rlse(ob)

c ----- solve big lsqfit problem --------------------
C|    call ftcsev(mmax,ndimb,nsdmbx,w(osb),001)
      call ftcslv(mmax,w(osb),w(obb),ndimb,nrhs,ndimbx,nsdmbx,nrhsx)
      call rlse(osb)

c ----- sort coeffs into cof by ialf-columns ---------
      call defrr(ocof,   ncof*nalf)
      call defrr(oba,    ndimx*nrhsx*(mmax+1)*nalf)
      call ftcsrt(mmax,nx1,lx1,ex1,nx2,lx2,ex2,lp1,lp2,rmt1,rmt2,
     .   lscal,lsym,ndim,nrhs,w(oba),w(obb),ndimx,nrhsx,ndimbx,
     .   nalf,w(ocof),ncof)
      call dpcopy(w(ocof),w(ocof),1,ncof,2d0)
      call rlse(oba)

c ----- get error of final tables ----------
      jrep=3
      if(ndist == 1) jrep=1
      do 70 irep=1,jrep
      if(irep == 1) d=d0
      if(irep == 2) d=d0 * 1.2
      if(irep == 3) d=d0 * 1.5
      x=dexp(adec*(d0-d))
      y=2d0*x-1d0
      call defrr(ocuf,  ncof)
      call defrr(owk,   ncof*2)
      call ropecs(y,nalf,ncof,w(owk),w(ocof),w(ocuf))
      mb1=15
      mb2=21
      jpr=0
      if(irep == 1) jpr=1
      call pshpr(jpr)
      call ftcchk(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,rmt1,rmt2,
     .  rsmp1,rsmp2,rsm1,rsm2,d,mb1,mb2,ri1,ri2,w(ocuf),err1)
      call poppr
      call rlse(ocuf)
c ----- get error for direct fit for comparison ----
      call pshpr(1)
      npmx=2*nb1*nb2+600
      call defrr(oxp,    npmx)
      call defrr(owp,    npmx)
      call defrr(ozp,    npmx)
      call holint(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
     .  zc1,zc2)
      lchk=1
      call ftcgen(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,d,nb1,nb2,w(oxp),w(ozp),
     .  w(owp),np,np,zc1,zc2,ri1,ri2,000,001,err2,000)
      call poppr
      call rlse(oxp)
      nn=(lp1+1)*(lp2+1)
      if(ipr >= 20) write(6,872) ((cl(i),cl(j),i=0,lp1),j=0,lp2)
      if(ipr >= 20) write(6,870) d,(err2(i),i=1,nn)
      if(ipr >= 20) write(6,871)   (err1(i),i=1,nn)
  870 format(' d=',f6.3,'  dir ',9f6.2:/(14x,9f6.2))
  871 format(9x,        '  hyf ',9f6.2:/(14x,9f6.2))
  872 format(/' rms errs  (%):',9(3x,a1,'*',a1):/(14x,9(3x,a1,'*',a1)))
CL      write(71,872) ((cl(i),cl(j),i=0,lp1),j=0,lp2)
CL      write(71,870) d,(err2(i),i=1,nn)
CL      write(71,871)   (err1(i),i=1,nn)
  70  continue
c ----- output into file ifi -------------------------
      call hyfout(rmt1,rmt2,rsm1,rsm2,d0,adec,adec,ndist,nalf,mmax,
     .   ncof,nb1,nb2,ri1,ri2,nxi,lp1,ep1,lp2,ep2,lsym,lx1,ex1,nx1,
     .   lx2,ex2,nx2,w(ocof),err1,ifi)

      call rlse(obb)
      end
