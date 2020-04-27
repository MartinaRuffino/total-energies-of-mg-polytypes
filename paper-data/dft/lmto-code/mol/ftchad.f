      subroutine ftchad(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .   rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,nalf,ndist,dist,wist,f,
     .   nb1,nb2,ri1,ri2,mmax,lmax,lsym,ndim,nrhs,ndimx,nrhsx,nsdmx,
     .   ndimbx,nsdmbx,b,s,bb,sb,lscal)
C- Adds together matrices sb,bb for hyperfit
C  lscal=1: scale to imporve condition; undo after solving lsq problem
      implicit real*8 (a-h,p-z), integer (o)
      dimension ex1(8),ex2(8),cx(300),dist(1),wist(1),
     .   ndim(0:20),nrhs(0:20),lx1(8),lx2(8),lx(20),
     .   s(nsdmx,0:mmax),b(ndimx,nrhsx,0:mmax),f(ndist,nalf),
     .   sb(nsdmbx,0:mmax),bb(ndimbx,nrhsx,0:mmax)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call sxlmnc(cx,10)
      if(iprint() >= 30) write(6,560) ndist,nalf
  560 format(/' ftchad:  ndist=',i3,'   nalf=',i3)
      call dpzero(bb,ndimbx*nrhsx*(mmax+1))
      call dpzero(sb,nsdmbx*(mmax+1))

c ------ start loop over distances ----------
      do 50 idist=1,ndist
      d=dist(idist)
      wgt=wist(idist)
      if(iprint() >= 31) write(6,550) idist,d,wgt
  550 format('          idist=',i3,'   d=',f11.6,'    wgt=',f11.6)
      jpr=ipr-20
      if(idist == 1) jpr=ipr
      call pshpr(jpr)

C ------ make mesh ------------
      npmx=2*nb1*nb2+600
      call defrr(oxp,    npmx)
      call defrr(owp,    npmx)
      call defrr(ozp,    npmx)
C|    call bisinl(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
C|   .  zc1,zc2)
      call holint(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
     .  zc1,zc2)

C ------ set up phi's and xi's -----
      nlm = ((lmax+1)*(lmax+2))/2
      nlm1 = ((lp1+1)*(lp1+2))/2
      nlm2 = ((lp2+1)*(lp2+2))/2
      call defrr(oxi,    np*nlm*(nx1+nx2))
      call defrr(oph1,   np*nlm1)
      call defrr(oph2,   np*nlm2)
      call ftcxip(lx1,ex1,nx1,lx2,ex2,nx2,rmt1,rmt2,rsmp1,rsmp2,
     .  rsm1,rsm2,lp1,ep1,lp2,ep2,cx,d,lscal,lmax,nx,lx,zc1,zc2,
     .  np,w(oxp),w(ozp),w(owp),w(oph1),w(oph2),w(oxi))

C ------ make the 2-center fit integrals -------
      call defrr(owk,   np)
      call dpzero(b,  ndimx*nrhsx*(mmax+1))
      call dpzero(s,  nsdmx*(mmax+1))
      call ftcnrm(mmax,lmax,lx,nx,w(owk),w(owp),np,w(oxi),nlm,cx,
     .  w(oph1),w(oph2),lp1,lp2,lsym,s,b,ndimx,nsdmx,nrhsx)

c ------ add to bb ---------
      do 10 m=0,mmax
      nd=ndim(m)
      nr=nrhs(m)
      do 11 ialf=1,nalf
      fac=f(idist,ialf)*wgt
      ia0=(ialf-1)*nd
      do 11 ir=1,nr
      do 11 id=1,nd
  11  bb(id+ia0,ir,m)=bb(id+ia0,ir,m)+fac*b(id,ir,m)
  10  continue

c ------ add to sb ----------------------
      do 20 m=0,mmax
      nd=ndim(m)
      do 20 ialf=1,nalf
      ia0=(ialf-1)*nd
      do 20 jalf=1,ialf
      ja0=(jalf-1)*nd
      fac=wgt*f(idist,ialf)*f(idist,jalf)
      is=0
      do 21 i=1,nd
      do 21 j=1,i
        is=is+1
        ia=i+ia0
        ja=j+ja0
        iii=(ia*(ia-1))/2+ja
        sb(iii,m)=sb(iii,m)+fac*s(is,m)
        if(jalf /= ialf.and.j /= i) then
          ia=j+ia0
          ja=i+ja0
          iii=(ia*(ia-1))/2+ja
          sb(iii,m)=sb(iii,m)+fac*s(is,m)
          endif
  21    continue
  20  continue

      call rlse(oxp)
      call poppr
  50  continue
      end
