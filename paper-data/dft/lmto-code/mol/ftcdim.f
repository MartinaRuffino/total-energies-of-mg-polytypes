      subroutine ftcdim(lx1,ex1,nx1,lx2,ex2,nx2,mmax,lmax,ep1,lp1,
     .   ep2,lp2,rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,lsym,ndim,nrhs,
     .   ndimx,nsdmx,nrhsx,ncof,nxi,lpr)
c  Get dimensions of the m-spaces, as well as number of rhs.
c  ncof=total number of coefficients.
      implicit real*8 (a-h,p-z), integer (o)
      dimension ndim(0:20),nrhs(0:20),lx1(1),lx2(1),ex1(1),ex2(1)
      fuzz=1.d-10
      mmax=lp1+lp2
c ------ determine whether case is symmetric --------
      lsym=1
      if(lp1 /= lp2.or.nx1 /= nx2)   lsym=0
      if(dabs(rsm1-rsm2) > fuzz)    lsym=0
      if(dabs(rsmp1-rsmp2) > fuzz)  lsym=0
      if(dabs(rmt1-rmt2) > fuzz)    lsym=0
      if(dabs(ep1-ep2) > fuzz)      lsym=0
      do 4 i=1,min0(nx1,nx2)
  4   if(lx1(i) /= lx2(i).or.dabs(ex1(i)-ex2(i)) > fuzz) lsym=0
c ------ print out specifications for functions -----
      if(lpr >= 1) then
      write(6,991) mmax,rmt1,rmt2
  991 format(' ftcdim:   mmax=',i3,'    rmt=',2f10.5)
      write(6,990) 1,(lx1(i),i=1,nx1),(200,i=nx1+1,6),
     .   (ex1(i),i=1,nx1)
      write(6,990) 2,(lx2(i),i=1,nx2),(200,i=nx2+1,6),
     .   (ex2(i),i=1,nx1)
      write(6,992) lp1,ep1,rsm1,lp2,ep2,rsm2,lsym
  990 format(' xi',i1,':   l= ',6(1x,i1),'   e=',6f7.2)
  992 format(' phi1:  lmx=',i3,'    e=',f8.3,'    rsm=',f7.3
     .  /' phi2:  lmx=',i3,'    e=',f8.3,'    rsm=',f7.3,/' lsym=',i3)
      endif
c ------ get number of rhs ----------------
      do 19 m=0,mmax
  19  nrhs(m)=0.d0
      nf1=0
      do 20 l1=0,lp1
      do 20 m1=0,l1
      nf1=nf1+1
      nf2=0
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 21 l2=0,ltop
      mtop=l2
      if(lsym == 1.and.l2 == l1) mtop=m1
      do 21 m2=0,mtop
      nf2=nf2+1
      mp=m1+m2
      mm=iabs(m1-m2)
      if(mp <= mmax) nrhs(mp)=nrhs(mp)+1
      if(mm /= mp.and.mm <= mmax) nrhs(mm)=nrhs(mm)+1
  21  continue
  20  continue
c ------ get nxi and lmax ------------------------
      nxi=0
      lmax=-1
      do 28 i=1,nx1
      lmax=max(lmax,lx1(i))
  28  nxi=nxi+(lx1(i)+1)**2
      do 29 i=1,nx2
      lmax=max(lmax,lx2(i))
  29  nxi=nxi+(lx2(i)+1)**2
c ------ get ndim, nsdm and printout -------------
      ndimx=0
      nsdmx=0
      nrhsx=0
      nrhst=0
      ncof=0
      do 10 m=0,mmax
      nd=0
      do 12 i=1,nx1
  12  if(m <= lx1(i)) nd=nd+(lx1(i)-m+1)
      do 13 i=1,nx2
  13  if(m <= lx2(i)) nd=nd+(lx2(i)-m+1)
      nsdm=(nd*(nd+1))/2
      ndim(m)=nd
      if(nd == 0) call rx('ftcdim: an m-space has dimension zero')
      ndimx=max0(ndimx,ndim(m))
      nsdmx=max0(nsdmx,nsdm)
      nrhsx=max0(nrhsx,nrhs(m))
      nrhst=nrhst+nrhs(m)
  10  ncof=ncof+ndim(m)*nrhs(m)
      if(lpr == 1) then
      write(6,941) (ndim(m),m=0,mmax)
      write(6,942) (nrhs(m),m=0,mmax)
  941 format(' ndim=',10i5)
  942 format(' nrhs=',10i5)
      write(6,950) ndimx,nsdmx,nrhsx,nrhst,ncof,nxi
  950 format(/' max ndim:',i7,'       max nsdm:',i7
     .       /' max nrhs:',i7,'       tot nrhs:',i7
     .       /' total number of non-zero coeffs:',i7
     .       /' number of 3-d fit functions, xi:',i7)
      endif
c ------ write to log ----------
      if(lpr == 0) return
CL      write(71,710) ep1,ep2,lp1,lp2
  710 format(' ------- ftcdim:  eph=',2f9.3,'   lph=',2i3,' -------')
CL      write(71,891) mmax,(ndim(m),m=0,mmax)
  891 format(' mmax',i3,'   ndim',10i4)
CL      write(71,890) 1,(lx1(i),i=1,nx1),(200,i=nx1+1,6),
CL     .   (ex1(i),i=1,nx1)
CL      write(71,890) 2,(lx2(i),i=1,nx2),(200,i=nx2+1,6),
CL     .   (ex2(i),i=1,nx1)
CL      write(71,892) rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,lsym
  890 format(' xi',i1,'   l= ',6(1x,i1),'   e=',6f7.2)
  892 format(' rmt',2f7.3,'   rsm',4f7.3,'   lsym=',i1)
      end
