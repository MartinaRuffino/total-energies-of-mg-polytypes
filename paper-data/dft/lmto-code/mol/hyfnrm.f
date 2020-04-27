      subroutine hyfnrm(mmax,lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,
     .  lscal,rmt1,rmt2,ri1,ri2,rsmp1,rsmp2,rsm1,rsm2,
     .  lp1,ep1,np1,lp2,ep2,np2,idxj,nph,lsymj,s,b,
     .  ndimx,nsdmx,nrhsx,dist,wist,ndist)
C- Makes normal matrix for hyfgen (s=overlap and b=rhs)
      implicit none
      integer lp1(1),lp2(1),lsymj(1),mmax,nb1,nb2,nph,idxj(2,1),np1,np2,
     .  nx1,nx2,ndimx,ndist,nrhsx,nsdmx,lscal
      integer lx1(6),lx2(6),lx(20)
      double precision ep1,ep2,ri1,ri2,zc1,zc2,
     .  rsm1,rsm2,rsmp1,rsmp2,rmt1,rmt2,
     .  ex1(1),ex2(1),dist(ndist),
     .  s(nsdmx,0:mmax,ndist),b(ndimx,nrhsx,0:mmax,ndist),
     .  wist(1),cx(300)
      integer iprint,oph1,oph2,oxi,owp,oxp,ozp,np,
     .  idist,lmax,
     .  nlm,nx,lmxp1,lmxp2,nlmp1,nlmp2,nbisi

      double precision d,wgt
      real w(1)
      common /w/ w

C --- Setup ---
      call tcn('hyfnrm')
      call sxlmnc(cx,12)

C --- Start loop over distances ---
      do 50 idist=1,ndist
      if (idist > 1 .and. idist < ndist) call pshpr(iprint()-20)
      d  = dist(idist)
      wgt= wist(idist)
      if(iprint() >= 30) write(6,550) idist,ndist,d,wgt
  550 format(' hyfnrm: idist=',i3,' of',i3,'   d=',f9.6,'   wgt=',f9.6)

      call hyfxim(lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,
     .  ri1,ri2,rsmp1,rsmp2,rsm1,rsm2,rmt1,rmt2,lscal,1,
     .  np,nbisi,lp1,ep1,np1,lp2,ep2,np2,cx,d,
     .  lmax,lmxp1,lmxp2,nx,lx,zc1,zc2,oph1,oph2,owp,oxi,oxp,ozp)
      nlm = ((lmax+1)*(lmax+2))/2
      nlmp1 = ((lmxp1+1)*(lmxp1+2))/2
      nlmp2 = ((lmxp2+1)*(lmxp2+2))/2

C --- Make the 2-center fit integrals at this distance -------
      call tcfnrm(mmax,lmax,lx,nx,w(owp),np,nbisi,w(oxi),nlm,
     .  w(oph1),w(oph2),lp1,lp2,nlmp1,nlmp2,idxj,nph,lsymj,
     .  s(1,0,idist),b(1,1,0,idist),ndimx,nsdmx,nrhsx)
      call rlse(oxp)
      if (idist > 1 .and. idist < ndist) call poppr
   50 continue
      call tcx('hyfnrm')
      end
      subroutine hyfxim(lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,
     .  ri1,ri2,rsmp1,rsmp2,rsm1,rsm2,rmt1,rmt2,lscal,lsw,
     .  np,nbisi,lp1,ep1,np1,lp2,ep2,np2,
     .  cx,d,
     .  lmax,lmxp1,lmxp2,nx,lx,zc1,zc2,oph1,oph2,owp,oxi,oxp,ozp)
C-  Generate mesh and tabulate phi's and xi's there
      implicit none
      double precision cx(1),d,ep1,ep2,ex1(1),ex2(1),ri1,ri2,
     .  rsm1,rsm2,rsmp1,rsmp2,zc1,zc2,rmt1,rmt2
      integer lmax,lp1(1),lp2(1),lx1(1),lx2(1),lscal,lsw,
     .  oidx,oph1,oph2,owk,owp,oxi,oxp,oyl,ozp,np1,np2,
     .  np,nbisi,nb1,nb2,nx,nx1,nx2,lx(20)
      integer iprint,nlmp1,nlmp2,nlm,npmx,ix,lmxp1,lmxp2
      integer ndCRAY
C#ifndef CRAY-DP
      parameter (ndCRAY=1)
C#elseC
C      parameter (ndCRAY=2)
C#endif
      real w(1)
      common /w/ w

C --- Setup ---
      lmax = 0
      do  10  ix = 1, nx1
   10 lmax = max(lmax,lx1(ix))
      do  11  ix = 1, nx2
   11 lmax = max(lmax,lx2(ix))
      nlm = ((lmax+1)*(lmax+2))/2
      lmxp1 = 0
      do  12  ix = 1, np1
   12 lmxp1 = max(lmxp1,lp1(ix))
      lmxp2 = 0
      do  14  ix = 1, np2
   14 lmxp2 = max(lmxp2,lp2(ix))
      nlmp1 = ((lmxp1+1)*(lmxp1+2))/2
      nlmp2 = ((lmxp2+1)*(lmxp2+2))/2
      if (nlmp1 > nlm .or. nlmp2 > nlm) call rx('hyfxim: nlmp gt nlm')

C --- Bispherical mesh, allocate arrays ---------------
      nbisi = 2*nb1*nb2
      npmx =  nbisi + 600
      call defrr(oxp,    npmx)
      call defrr(owp,    npmx)
      call defrr(ozp,    npmx)
      call holint(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
     .  zc1,zc2)
      call defrr(owk, np*3)
      call defrr(oidx,np)
      call hsrint(w(oxp),w(ozp),w(owp),w(owk),w(owk),w(oidx),nbisi)
      call rlse(owk)
      if (npmx < np) call rx('hyfxim: np gt npmx')

C --- Tabulate xi's and phi's for this mesh ----------
      call defrr(oxi,    np*nlm*(nx1+nx2))
      call defrr(oph1,   np*nlmp1*np1)
      call defrr(oph2,   np*nlmp2*np2)
      call defrr(owk,    np*ndCRAY)
      call defrr(oyl,    np*max(6,nlm)*ndCRAY)
      call defrr(oidx,   np*2)
      call hyfxi(np,w(oxp),w(ozp),w(owp),lscal,lsw,lmax,nlm,rmt1,rmt2,
     .  zc1,nx1,lx1,ex1,lmxp1,nlmp1,rsmp1,rsm1,lp1,ep1,np1,
     .  zc2,nx2,lx2,ex2,lmxp2,nlmp2,rsmp2,rsm2,lp2,ep2,np2,
     .  cx,w(oidx),w(owk),w(oyl),w(oph1),w(oph2),nx,lx,w(oxi))
      call rlse(owk)
      end
      subroutine hyfxi(np,xp,zp,wp,lscal,lws,lmax,nlm,rmt1,rmt2,
     .  zc1,nx1,lx1,ex1,lmxp1,nlmp1,rsmp1,rsm1,lp1,ep1,np1,
     .  zc2,nx2,lx2,ex2,lmxp2,nlmp2,rsmp2,rsm2,lp2,ep2,np2,
     .  cx,idx,wk,yl,ph1,ph2,nx,lx,xi)
C- Makes phi's and xi's.
C  lscal=1: scale to make rmsavg xi=1 at rmt.
C  lws=1:   fold in weights; wp are DESTROYED on output.
      implicit none
      integer nx1,lx1(1),lp1(1),np1,lmxp1,nlmp1,np,lscal,lws,nlm,
     .        nx2,lx2(1),lp2(1),np2,lmxp2,nlmp2,lx(1),nx
      double precision xp(1),zp(1),wp(1),
     .  cx(1),wk(1),idx(1),yl(1),tol,zc,xi(np,nlm,1),
     .  ex1(1),ep1,zc1,rsmp1,rsm1,ph1(np,nlmp1,1),rmt1,
     .  ex2(1),ep2,zc2,rsmp2,rsm2,ph2(np,nlmp2,1),rmt2
      integer ie,ip,l,m,lm,lmax

C --- Concatenate lxi for xi's centered on two sites ---
      call tcn('hyfxi')
      nx = nx1+nx2
      do  5  ie = 1, nx1
    5 lx(ie) = lx1(ie)
      do  6  ie = 1, nx2
    6 lx(ie+nx1) = lx2(ie)

      tol = 1d-10
      zc = -zc1
      do  10  ip = 1, np
   10 zp(ip) = zp(ip) + zc

C --- Smoothed Hankels xi and phi centered at origin ---
      call sol0sr(rsm1,lmax,nlm,nx1,lx1,ex1,xp,zp,np,np,
     .  tol,cx,idx,wk,yl,10,xi)
      call sol0sr(rsmp1,lmxp1,nlmp1,np1,lp1,ep1,xp,zp,np,np,
     .  tol,cx,idx,wk,yl,0,ph1)
      if (lscal == 1) call hyfxsc(nx1,lx1,ex1,lmax,rmt1,np,nlm,1,xi)

      zc = zc1 - zc2
      do  25  ip = 1, np
   25 zp(ip) = zp(ip) + zc

C --- xi and phi centered off origin ---
      call sol0sr(rsm2,lmax,nlm,nx2,lx2,ex2,xp,zp,np,np,
     .  tol,cx,idx,wk,yl,10,xi(1,1,nx1+1))
      call sol0sr(rsmp2,lmxp2,nlmp2,np2,lp2,ep2,xp,zp,np,np,
     .  tol,cx,idx,wk,yl,0,ph2)
      if (lscal == 1)
     .  call hyfxsc(nx2,lx2,ex2,lmax,rmt2,np,nlm,1,xi(1,1,nx1+1))

C --- Scale each xi by sqrt(wp) and each phi by sqrt(sqrt(wp)) ---
      zc = zc2
      do  34  ip = 1, np
   34 zp(ip) = zp(ip) + zc

      if (lws == 0) return
      do  40  ip = 1, np
   40 wp(ip) = dsqrt(dabs(wp(ip)))

      do  50  ie = 1, nx
      do  50  l = 0, lx(ie)
      do  50  m = 0, l
      lm = ((2*lmax-m+1)*m)/2 + l + 1
      do  50  ip = 1, np
   50 xi(ip,lm,ie) = wp(ip)*xi(ip,lm,ie)

      do  60  ip = 1, np
   60 wp(ip) = dsqrt(wp(ip))

      do  70 ie = 1, np1
      do  70  l = 0, lp1(ie)
      do  70  m = 0, l
      lm = (l*(l+1))/2+m+1
      do  70  ip = 1, np
   70 ph1(ip,lm,ie) = wp(ip)*ph1(ip,lm,ie)

      do  80 ie = 1, np2
      do  80  l = 0, lp2(ie)
      do  80  m = 0, l
      lm = (l*(l+1))/2+m+1
      do  80  ip = 1, np
   80 ph2(ip,lm,ie) = wp(ip)*ph2(ip,lm,ie)

      call tcx('hyfxi')
      end
