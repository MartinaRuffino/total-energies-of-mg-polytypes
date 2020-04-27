      subroutine hyfevg(dr,r1,e1,nlm1,r2,e2,nlm2,x1,x2,g1,g2,rmat,
     .   nrot,nf1,nf2,nxi1,lxi1,exi1,nxi2,lxi2,exi2,tspec,tdata,itab)
c  Evaluates unrotated 2-center expansion. First searches through
c  tspec to find the correct table, then interpolates the tabulated
c  coefficients. Expansion coefficients for all nlm1*nlm2 products
c  are returned in x1 and x2. Also returns rotation matrix in rmat.
c  hyfevg does the same as hyfevl, but also returns the derivatives
c  of the fit coefficients respective to dr in g1 and g2.
      implicit real*8 (a-h,p-z), integer (o)
      dimension tspec(100,1),lxi1(1),lxi2(1),exi1(1),exi2(1),
     .   jxi1(10),jxi2(10),tdata(1),dr(3),lx1(10),lx2(10),
     .   x1(nf1,nlm1,nlm2),x2(nf2,nlm2,nlm1),dir(3),rmat(nrot,nrot),
     .   g1(nf1,nlm1,nlm2),g2(nf2,nlm2,nlm1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      l1=ll(nlm1)
      l2=ll(nlm2)
c ------- locate the correct table -----------------
      call hyfloc(r1,l1,e1,r2,l2,e2,nxi1,lxi1,exi1,
     .   nxi2,lxi2,exi2,tspec,itbl,jxi1,jxi2)
      if(itbl == 0) call rx('hyfevg:  no suitable table found''')
      it=iabs(itbl)
      ndt0=idnint(tspec(25,it))
      if(ipr >= 80) write(6,650) r1,r2,l1,l2,e1,e2,itbl,ndt0
  650 format(' hyfget:   r=',2f8.3,'   l=',2i3,'   e=',2f9.4,
     .   /' itbl=',i3,'    data offset=',i10)
      itab=itbl
c ------- get specifications of table ---------
      d0  =       tspec( 3,it)
      a1  =       tspec( 5,it)
      nalf=idnint(tspec( 7,it))
      mmax=idnint(tspec( 8,it))
      ncof=idnint(tspec( 9,it))
      nxi =idnint(tspec(12,it))
      lsym=idnint(tspec(13,it))
      lp1 =idnint(tspec(15,it))
      lp2 =idnint(tspec(17,it))
      nx1 =idnint(tspec(19,it))
      nx2 =idnint(tspec(20,it))
      do 14 i=1,nx1
  14  lx1(i)=idnint(tspec(40+i,it))
      do 15 i=1,nx2
  15  lx2(i)=idnint(tspec(60+i,it))
c ------- printout ---------------------
      if(ipr >= 80) then
      write(6,600) d0,a1,nalf,ncof,mmax,nxi,lsym
  600 format(' d0=',f8.4,'   a1=',f8.4,'   nalf=',i3,'   ncof=',i6
     .   /' mmax=',i3,'    nxi=',i5,'    lsym=',i2)
      write(6,601) lp1,nx1,(lx1(k),k=1,nx1)
      write(6,601) lp2,nx2,(lx2(k),k=1,nx2)
  601 format(' lp=',i2,'   nx=',i2,'   lx=',8i5)
      endif
c ------- some setup -------------
      mlm1=(lp1+1)**2
      mlm2=(lp2+1)**2
      lrot=max0(lp1,lp2)
      do 4 i=1,nx1
  4   lrot=max0(lrot,lx1(i))
      do 5 i=1,nx2
  5   lrot=max0(lrot,lx2(i))
      do 22 m=1,3
      dir(m)=dr(m)
  22  if(itbl < 0) dir(m)=-dr(m)
      call defrr(oa,   nxi*mlm1*mlm2)
      call defrr(od,   nxi*mlm1*mlm2)
c ------- evaluate interpolation at distance d ------------
      call defrr(ocof,    ncof)
      call defrr(odof,    ncof)
      d=dsqrt(dr(1)**2+dr(2)**2+dr(3)**2)
      call hyfipd(d,d0,a1,tdata(ndt0),w(ocof),w(odof),ncof,nalf)
C|    call hyfipl(d,d0,a1,tdata(ndt0),w(ocof),ncof,nalf)
c  primitive diff...
C|    hh=0.01d0
C|    call defrr(odf1,    ncof)
C|    call hyfipl(d+hh,d0,a1,tdata(ndt0),w(odof),ncof,nalf)
C|    call hyfipl(d-hh,d0,a1,tdata(ndt0),w(odf1),ncof,nalf)
C|    call dpadd(w(odof),w(odf1),1,ncof,-1d0)
C|    call dpcopy(w(odof),w(odof),1,ncof,1d0/(2*hh))
C|    call rlse(odf1)

      call tcfxpn(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,w(ocof),
     .   w(oa),nxi,mlm1,mlm2)
      call tcfxpn(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,w(odof),
     .   w(od),nxi,mlm1,mlm2)
      call rlse(ocof)
c ------- get matrix to rotate into direction of dir ------------
      if((lrot+1)**2 > nrot) call rx('hyfevg: nrot too small')
      call tcfrt1(dir,rmat,lrot,nrot)
c ------- copy over into x1,x2 -----------------
      if(itbl > 0) then
      call xxfevl(w(oa),nxi,mlm1,mlm2,x1,x2,nf1,nf2,nlm1,nlm2,
     .   nxi1,nxi2,lxi1,lxi2,nx1,nx2,lx1,lx2,jxi1,jxi2)
      call xxfevl(w(od),nxi,mlm1,mlm2,g1,g2,nf1,nf2,nlm1,nlm2,
     .   nxi1,nxi2,lxi1,lxi2,nx1,nx2,lx1,lx2,jxi1,jxi2)
      else
      call xxfevl(w(oa),nxi,mlm1,mlm2,x2,x1,nf2,nf1,nlm2,nlm1,
     .   nxi2,nxi1,lxi2,lxi1,nx1,nx2,lx1,lx2,jxi2,jxi1)
      call xxfevl(w(od),nxi,mlm1,mlm2,g2,g1,nf2,nf1,nlm2,nlm1,
     .   nxi2,nxi1,lxi2,lxi1,nx1,nx2,lx1,lx2,jxi2,jxi1)
      endif
      call rlse(oa)
      end
