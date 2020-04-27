      subroutine hyfevl(dr,r1,e1,nlm1,r2,e2,nlm2,x1,x2,rmat,nrot,
     .   nf1,nf2,nxi1,lxi1,exi1,nxi2,lxi2,exi2,tspec,tdata)
C  Evaluates unrotated 2-center expansion. First searches through
c  tspec to find the correct table, then interpolates the tabulated
c  coefficients. Expansion coefficients for all nlm1*nlm2 products
c  are returned in x1 and x2. Also returns rotation matrix in rmat.
      implicit real*8 (a-h,p-z), integer (o)
      dimension tspec(100,1),lxi1(1),lxi2(1),exi1(1),exi2(1),
     .   jxi1(10),jxi2(10),tdata(1),dr(3),lx1(10),lx2(10),
     .   x1(nf1,nlm1,nlm2),x2(nf2,nlm2,nlm1),dir(3),rmat(nrot,nrot)
      real w(1)
      common /w/ w
      call getpr(ipr)
      l1=ll(nlm1)
      l2=ll(nlm2)
c ------- locate the correct table -----------------
      call hyfloc(r1,l1,e1,r2,l2,e2,nxi1,lxi1,exi1,
     .   nxi2,lxi2,exi2,tspec,itbl,jxi1,jxi2)
      if(itbl == 0) call rx('hyfevl:  no suitable table found')
      it=iabs(itbl)
      ndt0=idnint(tspec(25,it))
      if(ipr >= 80) write(6,650) r1,r2,l1,l2,e1,e2,itbl,ndt0
  650 format(' hyfevl:   r=',2f8.3,'   l=',2i3,'   e=',2f9.4,
     .   /' itbl=',i3,'    data offset=',i10)
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
c ------- evaluate interpolation at distance d ------------
      call defrr(ocof,    ncof)
      d=dsqrt(dr(1)**2+dr(2)**2+dr(3)**2)
      call hyfipl(d,d0,a1,tdata(ndt0),w(ocof),ncof,nalf)
      call tcfxpn(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,w(ocof),
     .   w(oa),nxi,mlm1,mlm2)
      call rlse(ocof)
c ------- get matrix to rotate into direction of dir ------------
      if((lrot+1)**2 > nrot) call rx('hyfevl: nrot too small')
      call tcfrt1(dir,rmat,lrot,nrot)
c ------- copy over into x1,x2 -----------------
      if(itbl > 0) then
      call xxfevl(w(oa),nxi,mlm1,mlm2,x1,x2,nf1,nf2,nlm1,nlm2,
     .   nxi1,nxi2,lxi1,lxi2,nx1,nx2,lx1,lx2,jxi1,jxi2)
      else
      call xxfevl(w(oa),nxi,mlm1,mlm2,x2,x1,nf2,nf1,nlm2,nlm1,
     .   nxi2,nxi1,lxi2,lxi1,nx1,nx2,lx1,lx2,jxi2,jxi1)
      endif
      call rlse(oa)

      return
      end
c -------- sub xxfget --------------------------
      subroutine xxfevl(a,nxi,mlm1,mlm2,x1,x2,nf1,nf2,nlm1,nlm2,
     .   nxi1,nxi2,lxi1,lxi2,nx1,nx2,lx1,lx2,jxi1,jxi2)
c  copies the final tcf over to arrays x1,x2.
      implicit real*8 (a-h,p-z), integer (o)
      dimension x1(nf1,nlm1,nlm2),x2(nf2,nlm2,nlm1),jxi1(10),jxi2(10),
     .   a(nxi,mlm1,mlm2),lx1(1),lx2(1),lxi1(1),lxi2(1),j1(10),j2(10)
      call dpzero(x1,nf1*nlm1*nlm2)
      call dpzero(x2,nf2*nlm1*nlm2)
c -------- set up offsets for x1 and x2 -------
      j1(1)=0
      do 1 j=2,nxi1
  1   j1(j)=j1(j-1)+(lxi1(j-1)+1)**2
      j2(1)=0
      do 2 j=2,nxi2
  2   j2(j)=j2(j-1)+(lxi2(j-1)+1)**2
c -------- start looping over coeffs in a -----------
      do 10 ilm1=1,nlm1
      do 10 ilm2=1,nlm2
      ii0=0
      do 5 i=1,nx1
      j=jxi1(i)
      jj0=j1(j)
      nlmi=(lx1(i)+1)**2
      do 6 ilm=1,nlmi
  6   x1(jj0+ilm,ilm1,ilm2)=a(ii0+ilm,ilm1,ilm2)
  5   ii0=ii0+nlmi
      do 7 i=1,nx2
      j=jxi2(i)
      jj0=j2(j)
      nlmi=(lx2(i)+1)**2
      do 8 ilm=1,nlmi
  8   x2(jj0+ilm,ilm2,ilm1)=a(ii0+ilm,ilm1,ilm2)
  7   ii0=ii0+nlmi

  10  continue
      return
      end
