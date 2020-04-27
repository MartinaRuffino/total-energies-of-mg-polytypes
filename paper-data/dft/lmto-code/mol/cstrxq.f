      subroutine cstrxq(nxi,lxi,exi,ecut,q,tau,nrx,nlmx,lmax,wk,yl,
     .  a,alat,rlat,nkr,dlat,nkd,vol,cy,dl,dlp)
C- Reduced crystal strx, e < 0
C nlmx dimensions dl,dlp; nrx dimensions wk,yl.  Set lmax = max(lxi).
C wk must be dimensioned at least nrx*(2*lmax+9)
C Planned but not implemented:
C Sums for energies deeper than ecut are done directly in real space.
      implicit none
      integer nxi,lmax,nkr,nkd,nlmx,nrx,lxi(1)
      double precision alat,a,vol,ecut,exi(1),tau(3),q(3),
     .  wk(nrx,2*lmax+9),yl(nrx,(lmax+1)**2),cy(1),rlat(3,1),dlat(3,1),
     .  qdotr,y0,a2,pi,sp,gam,tpiba,xx,tpi,dl(2,nlmx,1),dlp(2,nlmx,1)
      integer ie,ir,ilm,l,m,ir1,nkd1,lc,ls
      complex*16 sm1

      real w(1)
      common /w/ w

      if (ecut < -1d-6) call rx('cstrxq: ecut not implemented')
      pi = 4*datan(1d0)
      y0 = 1/dsqrt(4*pi)
      a2 = a*a
      call dpzero(dl, 2*nlmx*nxi)
      call dpzero(dlp,2*nlmx*nxi)

C --- Energy-independent setup for Q-space part ---
      tpi   = 2*pi
      tpiba = 2*pi/alat
      gam = 0.25d0/a2
      do  10  ir = 1, nkr
        wk(ir,2) = tpiba*(q(1) + rlat(1,ir))
        wk(ir,3) = tpiba*(q(2) + rlat(2,ir))
        wk(ir,4) = tpiba*(q(3) + rlat(3,ir))
   10 continue
      call ropylm(nkr,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
      do  12  ir = 1, nkr
        sp = alat*(wk(ir,2)*tau(1) + wk(ir,3)*tau(2) + wk(ir,4)*tau(3))
        sm1 = cdexp(dcmplx(-gam*wk(ir,1), sp))
        wk(ir,2) = -dble(sm1)
        wk(ir,3) = -dimag(sm1)
   12 continue

C --- Q-space part of reduced strx for each energy ---
      call xstrxq(nxi,lxi,exi,a,vol,alat,nkr,wk,yl,nlmx,
     .  wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),dl,dlp)

C --- Energy-independent setup for real-space part ---
      if (nrx < max(nkd,nkr)) call rx('cstrxq: nrx < nkd,nkr')
      if (nlmx < (lmax+1)**2) call rx('cstrxq: nlmx < (lmax+1)**2')
      do  20  ir = 1, nkd
        wk(ir,2) = alat*(tau(1)-dlat(1,ir))
        wk(ir,3) = alat*(tau(2)-dlat(2,ir))
        wk(ir,4) = alat*(tau(3)-dlat(3,ir))
   20 continue
      call ropylm(nkd,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
      do  22  ir = 1, nkd
      qdotr = tpi*(q(1)*dlat(1,ir)+ q(2)*dlat(2,ir)+ q(3)*dlat(3,ir))
      wk(ir,4) = dcos(qdotr)
      wk(ir,5) = dsin(qdotr)
   22 wk(ir,2) = y0*dexp(-wk(ir,1)*a2)

C --- D-space part of reduced strx for each energy: xi(*,l)=wk(*,l+7)---
      ir1 = 1
      nkd1 = nkd
      if (wk(1,1) < 1d-6) then
        ir1 = 2
        nkd1 = nkd-1
      endif
      do  30  ie = 1, nxi
        if (lxi(ie) > lmax) call rx('cstrxq: lxi > lmax')
        call hansr4(wk,lxi(ie),nrx,nkd,exi(ie),1/a,wk(1,2),
     .    wk(1,3),wk(1,6))
C ...   Make cos(qR)*chi(aR), sin(qR)*chi(aR)
        lc = 6
        ls = lxi(ie)+8
        do  34  l = 0, lxi(ie)+1
        do  35  ir = ir1, nkd
   35   wk(ir,l+ls) = wk(ir,5)*wk(ir,l+lc)
        do  36  ir = ir1, nkd
   36   wk(ir,l+lc) = wk(ir,4)*wk(ir,l+lc)
   34   continue
        ilm = 0
        do  32  l = 0, lxi(ie)
        do  32  m = -l, l
          ilm = ilm+1
          call dpdot(yl(ir1,ilm),wk(ir1,l+lc+1),nkd1,xx)
          dl(1,ilm,ie)  = cy(ilm)*(dl(1,ilm,ie) + xx)
          call dpdot(yl(ir1,ilm),wk(ir1,l+ls+1),nkd1,xx)
          dl(2,ilm,ie)  = cy(ilm)*(dl(2,ilm,ie) + xx)
          call dpdot(yl(ir1,ilm),wk(ir1,l+lc),nkd1,xx)
          dlp(1,ilm,ie) = cy(ilm)*(dlp(1,ilm,ie) + xx/2)
          call dpdot(yl(ir1,ilm),wk(ir1,l+ls),nkd1,xx)
          dlp(2,ilm,ie) = cy(ilm)*(dlp(2,ilm,ie) + xx/2)
   32   continue
C ...   Site diagonal terms
        if (ir1 == 2) then
          dl(1,1,ie)  = dl(1,1,ie)  + cy(1)*wk(1,7)
          dlp(1,1,ie) = dlp(1,1,ie) + cy(1)*wk(1,6)
        endif
   30 continue

      end
      subroutine xstrxq(nx,lxi,exi,a,vol,alat,n,qsq,yl,nlmx,
     .  cs,sn,cs2,sn2,cs3,sn3,dl,dlp)
C- Reciprocal part of reduced strx
      implicit none
      integer lxi(1),nx,nlmx,n
      double precision alat,cof,e,gam,pi,vol
      double precision exi(1),a,qsq(n),yl(n,1),
     .  cs(n),sn(n),cs2(n),sn2(n),cs3(n),sn3(n),
     .  dl(2,nlmx,1),dlp(2,nlmx,1),c,s,xx,c2,c3,s2,s3
      integer i,ie,ilm,l,lx,m

      gam = 0.25d0/a**2
      pi = 4*datan(1d0)

C --- For each exi ne 0 and ilm do ---
      do  20  ie = 1, nx
        e = exi(ie)
        lx = lxi(ie)
        do  22  i = 1, n
          xx = 1/(e-qsq(i))
          c = cs(i)*xx
          cs2(i) = c
          cs3(i) = c*xx
          s = sn(i)*xx
          sn2(i) = s
          sn3(i) = s*xx
   22   continue
        cof = -4d0*pi*dexp(gam*e)/vol
        ilm = 0
        do  20  l = 0, lx, 2
          cof = -cof
          do  24  m = -l, l
            ilm = ilm+1
            call dpdot(yl(1,ilm),cs2,n,c2)
            dl(1,ilm,ie) = dl(1,ilm,ie) + cof*c2
            call dpdot(yl(1,ilm),cs3,n,c3)
            dlp(1,ilm,ie) = dlp(1,ilm,ie) + cof*(gam*c2 - c3)
            call dpdot(yl(1,ilm),sn2,n,s2)
            dl(2,ilm,ie) = dl(2,ilm,ie) + cof*s2
            call dpdot(yl(1,ilm),sn3,n,s3)
            dlp(2,ilm,ie) = dlp(2,ilm,ie) + cof*(gam*s2 - s3)
   24     continue
          if (l+1 <= lx) then
            do  26  m = -l-2, l
              ilm = ilm+1
              call dpdot(yl(1,ilm),sn2,n,s2)
              dl(1,ilm,ie) = dl(1,ilm,ie) + cof*s2
              call dpdot(yl(1,ilm),sn3,n,s3)
              dlp(1,ilm,ie) = dlp(1,ilm,ie) + cof*(gam*s2 - s3)
              call dpdot(yl(1,ilm),cs2,n,c2)
              dl(2,ilm,ie) = dl(2,ilm,ie) - cof*c2
              call dpdot(yl(1,ilm),cs3,n,c3)
              dlp(2,ilm,ie) = dlp(2,ilm,ie) - cof*(gam*c2 - c3)
   26       continue
          endif
   20 continue
      end
C      subroutine fmain
C      implicit none
C      integer nkdmx,nkrmx,lmxx,nlmx,nxx
C      parameter (nkdmx=1000, nkrmx=1000, lmxx=8,nlmx=(lmxx+1)**2, nxx=6)
C      double precision rlat(3,nkdmx), dlat(3,nkrmx),
C     .  wk(nkdmx*(2*lmxx+9)),yl(nkdmx*nlmx),cy(300)
C      complex*16 dl(nlmx,nxx),dlp(nlmx,nxx),dl0(nlmx),dlp0(nlmx)
C      double precision rb(3,3),qb(3,3),rb0(3,3),qb0(3,3),tau(3),q(3),
C     .  vol,a,as,tol,alat,ecut,exi(10),e,dmx,dpmx
C      integer lmax,nkd,nkr,lxi(10),ie,nrx,nxi,ilm
Cc     data rb0 /0d0,.5d0,.5d0,.5d0,0d0,.5d0,.5d0,.5d0,0d0/
C      data rb0 /0d0,.47d0,.52d0,.5d0,0d0,.5d0,.8d0,.5d0,0d0/
C      data tau /.1d0,.4d0,.2d0/  q /.2d0,-.1d0,.27d0/
C
C      common /static/ rlat,dlat,wk,yl,dl,dlp
C      real w(100000)
C      common /w/ w
C
C      call wkinit(100000)
C      as = 2
C      tol = 1d-6
C      lmax = 8
C      alat=7.6d0
C      call sylmnc(cy,12)
CC ... wk a temporary array here
C      call lattc(as,tol,alat,alat,rb0,1d0,1d0,1d0,1d0,rb,qb,
C     .   lmax,vol,a,dlat,nkd,rlat,nkr,nkdmx,nkrmx,wk)
C
C      nxi = 4
C      lxi(1) = 5
C      exi(1) = -3
C      lxi(2) = 7
C      exi(2) = -0.2
C      lxi(3) = 8
C      exi(3) = -8
C      lxi(4) = 1
C      exi(4) = -1.4
C
CC      nxi = 1
CC      exi(1) = -0.2
CC      exi(2) = -0.2
CC      lxi(1) = 2
CC      lxi(2) = 2
C      print *, 'tau=',tau
C      print *, '  q=',q
C      print *, 'set tau,q='
C      read(*,*) tau,q
C
C      lmax = -1
C      do  10  ie = 1, nxi
C   10 lmax = max(lmax,lxi(ie))
C      nrx = max(nkd,nkr)
C      call cstrxq(nxi,lxi,exi,ecut,q,tau,nrx,nlmx,lmax,wk,yl,
C     .  a,alat,rlat,nkr,dlat,nkd,vol,cy,dl,dlp)
C
C      print 334
C  334 format(/' Comparison of cstrxq with rcnsl:'/
C     .  11x,'dl',12x,'error',11x,'dlp',12x,'error')
C      do  20  ie = 1, nxi
C       e = exi(ie)
C       call rcnsl(e,q,tau,a,lmax,alat,rlat,nkr,dlat,nkd,vol,cy,dl0,dlp0)
C       print *, ' --- e=', sngl(e), '  lxi=',lxi(ie)
C       dmx  = 0
C       dpmx = 0
C       do  24  ilm = 1, (lxi(ie)+1)**2
C
C       if (cdabs(dl0(ilm))+cdabs(dl(ilm,ie))+
C     .     cdabs(dlp0(ilm))+cdabs(dlp(ilm,ie)) > 1d-18) then
C         print 333, ilm, dble(dl0(ilm)), dble(dl(ilm,ie)-dl0(ilm)),
C     .                   dble(dlp0(ilm)),dble(dlp(ilm,ie)-dlp0(ilm))
C         print 336,      dimag(dl0(ilm)), dimag(dl(ilm,ie)-dl0(ilm)),
C     .                   dimag(dlp0(ilm)),dimag(dlp(ilm,ie)-dlp0(ilm))
C       endif
C       dmx = max(dmx,cdabs(dl(ilm,ie)-dl0(ilm)))
C       dpmx= max(dpmx,cdabs(dlp(ilm,ie)-dlp0(ilm)))
C  333  format(i4,2(1pe18.10,1pe12.3),2x)
C  336  format(4x,2(1pe18.10,1pe12.3),2x)
C   24 continue
C      print 335, dmx,dpmx
C  335 format(' max error',12x,1pe12.3,18x,1pe12.3)
C   20 continue
C
C      end
