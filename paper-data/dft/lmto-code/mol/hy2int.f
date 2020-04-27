      subroutine hy2int(mmax,lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,lmxl,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,
     .  nr1,nth1,rad1,wrad1,w1,phr1,nr2,nth2,rad2,wrad2,w2,phr2,
     .  npmx,lp1,ep1,lp2,ep2,lsym,s,b,ndim,nrhs,ndimx,
     .  nsdmx,nrhsx,sb,bb,ndimbx,nsdmbx,dist,wist,f,ndist,nalf)
C- Makes integrals for hyf (sb=overlap and bb=rhs)
      implicit none
      integer lp1,lp2,lsym,mmax,nb1,nb2,nalf,ndimbx,ndimx,ndist,nrhsx,
     .  nr1,nr2,nsdmbx,nsdmx,nth1,nth2,nx1,nx2,npmx,lmxl
      integer ndim(0:20),nrhs(0:20),lx1(6),lx2(6),lx(20)
      double precision ep1,ep2,rmt1,rmt2,zc1,zc2,
     .  rad1(1),rad2(1),wrad1(1),wrad2(1),rsm1,rsm2,rsmp1,rsmp2,w1,w2,
     .  ex1(1),ex2(1),dist(ndist),s(nsdmx,0:mmax),b(ndimx,nrhsx,0:mmax),
     .  sb(nsdmbx,0:mmax),bb(ndimbx,nrhsx,0:mmax),wist(1),f(ndist,nalf),
     .  cx(300),phr1(nr1,0:1),phr2(nr2,0:2)
      integer iprint,oph1,oph2,oidx,oxi,owk,oyl,owp,oxp,ozp,np,nbisi,
     .  i,ia,ia0,ialf,id,idist,iii,ir,is,ix,j,ja,ja0,jalf,lmax,
     .  m,nadd,nd,nlm,nlmp1,nlmp2,nr,nsd,nx,ioff1,ioff2
      double precision d,fac,wgt
      real w(1)
      common /w/ w

C --- Setup ---
      if (iprint() >= 30) print 555, ndist
  555 format(/' --------- hy2int:   ndist=',i3,' ---------')
      if (iprint() >= 60) then
        do 1 m=0,mmax
          nd=ndim(m)
          nsd=(nd*(nd+1))/2
          write(6,970) m,nd,nsd,nrhs(m)
  970     format(' m=',i2,'   ndim=',i5,'   nsd=',i5,'   nrhs=',i5)
    1   continue
      endif
      call sxlmnc(cx,12)
      call dpzero(sb,nsdmbx*(mmax+1))
      call dpzero(bb,ndimbx*nrhsx*(mmax+1))
      lmax = 0
      do  10  ix = 1, nx1
   10 lmax = max(lmax,lx1(ix))
      do  11  ix = 1, nx2
   11 lmax = max(lmax,lx2(ix))
      nlm = ((lmax+1)*(lmax+2))/2

C --- Start loop over distances --------------------
      do 50 idist=1,ndist
      if (idist > 1 .and. idist < ndist) call pshpr(iprint()-10)
      d=dist(idist)
      wgt=wist(idist)
      if(iprint() >= 40) write(6,550) idist,d,wgt
      if(iprint() > 40) write(6,551) (f(idist,i),i=1,nalf)
  550 format(' hy2int: idist=',i3,'   d=',f9.6,'   wgt=',f9.6)
  551 format(' f= ',6f12.6)

C --- Bispherical mesh, allocate arrays ---------------
      call defrr(oxp,    npmx)
      call defrr(owp,    npmx)
      call defrr(ozp,    npmx)
      call hy2msh(nr1,nth1,rad1,wrad1,w1,nr2,nth2,rad2,wrad2,w2,nadd,
     .  rmt1,rmt2,d,nb1,nb2,w(oxp),w(ozp),w(owp),npmx,nbisi,np,zc1,zc2)
      call defrr(oxi,    np*nlm*(nx1+nx2))
      nlmp1 = ((lp1+1)*(lp1+2))/2
      call defrr(oph1,   np*nlmp1)
      nlmp2 = ((lp2+1)*(lp2+2))/2
      call defrr(oph2,   np*nlmp2)
      call defrr(owk,    np)
      call defrr(oidx,   np*2)
      call defrr(oyl,    np*max(6,nlm,nlmp1,nlmp2))
      ioff1 = nbisi
      ioff2 = nbisi + nr1*nth1

C --- Tabulate xi's and phi's for this mesh ----------
      if (iprint() >= 60) call tm('start smoothed h')
      call defrr(owk,    np)
      call defrr(oidx,   np*2)
      call defrr(oyl,    np*max(6,nlm))
      stop 'update call to hy2xi'
C      call hy2xi(np,w(oxp),w(ozp),w(owp),1,lmax,nlm,
C     . zc1,nx1,lx1,ex1,nlmp1,rsmp1,rsm1,lp1,ep1,phr1,ioff1,nr1,nth1,
C     . zc2,nx2,lx2,ex2,nlmp2,rsmp2,rsm2,lp2,ep2,phr2,ioff2,nr2,nth2,
C     . cx,w(oidx),w(owk),w(oyl),w(oph1),w(oph2),nx,lx,w(oxi))
      call rlse(owk)

C --- Make the 2-center fit integrals at this distance -------
      if (iprint() >= 60) call tm('start tcfint')
      call tcfint(mmax,lmax,lx,nx,w(owp),np,nbisi,lmxl,w(oxi),nlm,cx,
     .  w(oph1),w(oph2),lp1,lp2,lsym,s,b,ndimx,nsdmx,nrhsx)
      call rlse(oxp)

C --- Start loop over m, first add to bb -------
      if (iprint() >= 60) call tm('start add to bb')
      do 20 m=0,mmax
      nd=ndim(m)
      nr=nrhs(m)
      do 21 ialf=1,nalf
      fac=f(idist,ialf)*wgt
      ia0=(ialf-1)*nd
      do 21 ir=1,nr
      do 21 id=1,nd
   21 bb(id+ia0,ir,m)=bb(id+ia0,ir,m)+fac*b(id,ir,m)
C --- Next add to sb ----------------------
      do 30 ialf=1,nalf
      ia0=(ialf-1)*nd
      do 30 jalf=1,ialf
      ja0=(jalf-1)*nd
      fac=wgt*f(idist,ialf)*f(idist,jalf)
      is=0
      do 31 i=1,nd
      do 31 j=1,i
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
   31   continue
   30 continue
   20 continue
      if (idist > 1 .and. idist < ndist) call poppr
   50 continue
      end
