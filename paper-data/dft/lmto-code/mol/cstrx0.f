      subroutine cstrx0(nxi,lxi,exi,ecut,tau,nrx,nlmx,lmax,wk,yl,
     .  a,alat,rlat,nkr,dlat,nkd,vol,cy,dl,dlp)
C- Reduced crystal strx, Q=0
C nlmx dimensions dl,dlp; nrx dimensions wk,yl.  Set lmax = max(lxi).
C wk must be dimensioned at least wk(max(nkd,nkr),max(lmax+5,7))
C For now, if exi=0 it must occur last.  dlp not defined for exi=0
C Planned but not implemented:
C Sums for energies deeper than ecut are done directly in real space.
      implicit none
      integer nxi,lmax,nkr,nkd,nlmx,nrx,lxi(1)
      double precision alat,a,vol,ecut,exi(1),tau(3),
     .  wk(nrx,lmax+5),yl(nrx,(lmax+1)**2),
     .  dl(nlmx,1),dlp(nlmx,1),cy(1),rlat(3,1),dlat(3,1),
     .  y0,a2,pi,sp,gam,tpiba,xx
      integer ie,ir,ilm,l,m,ir1,nkd1
      complex*16 sm1

      real w(1)
      common /w/ w

      if (ecut < -1d-6) call rx('cstrx0: ecut not implemented')
      pi = 4*datan(1d0)
      y0 = 1/dsqrt(4*pi)
      a2 = a*a
      call dpzero(dl, nlmx*nxi)
      if (exi(1) < -1d-6) call dpzero(dlp,nlmx*nxi)

C --- Energy-independent setup for Q-space part ---
      tpiba = 2*pi/alat
      gam = 0.25d0/a2
      do  10  ir = 1, nkr
        wk(ir,2) = tpiba*rlat(1,ir)
        wk(ir,3) = tpiba*rlat(2,ir)
        wk(ir,4) = tpiba*rlat(3,ir)
   10 continue
      call ropylm(nkr,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
      do  12  ir = 1, nkr
        sp = alat*(wk(ir,2)*tau(1) + wk(ir,3)*tau(2) + wk(ir,4)*tau(3))
        sm1 = cdexp(dcmplx(-gam*wk(ir,1), sp))
        wk(ir,2) = -dble(sm1)
        wk(ir,3) = -dimag(sm1)
   12 continue

C --- Q-space part of reduced strx for each energy ---
      call xstrx0(nxi,lxi,exi,a,vol,alat,nkr,wk,yl,nlmx,
     .  wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),dl,dlp)

C --- Energy-independent setup for real-space part ---
      if (nrx < max(nkd,nkr)) call rx('cstrx0: nrx < nkd,nkr')
      if (nlmx < (lmax+1)**2) call rx('cstrx0: nlmx < (lmax+1)**2')
      do  20  ir = 1, nkd
        wk(ir,2) = alat*(tau(1)-dlat(1,ir))
        wk(ir,3) = alat*(tau(2)-dlat(2,ir))
        wk(ir,4) = alat*(tau(3)-dlat(3,ir))
   20 continue
      call ropylm(nkd,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
      do  22  ir = 1, nkd
   22 wk(ir,2) = y0*dexp(-wk(ir,1)*a2)

C --- D-space part of reduced strx for each energy: xi(*,l)=wk(*,l+5)---
      ir1 = 1
      nkd1 = nkd
      if (wk(1,1) < 1d-6) then
        ir1 = 2
        nkd1 = nkd-1
      endif
      do  30  ie = 1, nxi
        if (lxi(ie) > lmax) call rx('cstrx0: lxi > lmax')
        if (dabs(exi(ie)) > 1d-6) then
          call hansr4(wk,lxi(ie),nrx,nkd,exi(ie),1/a,wk(1,2),
     .      wk(1,3),wk(1,4))
        else
          call hansr5(wk,lxi(ie),nrx,nkd,1/a,wk(1,2),wk(1,3),wk(1,5))
        endif
        ilm = 0
        do  32  l = 0, lxi(ie)
        do  32  m = -l, l
          ilm = ilm+1
          call dpdot(yl(ir1,ilm),wk(ir1,l+5),nkd1,xx)
          dl(ilm,ie)  = cy(ilm)*(dl(ilm,ie) + xx)
          if (exi(ie) < -1d-6) then
            call dpdot(yl(ir1,ilm),wk(ir1,l+4),nkd1,xx)
            dlp(ilm,ie) = cy(ilm)*(dlp(ilm,ie) + xx/2)
          endif
   32   continue
C ...   Site diagonal terms
        if (ir1 == 2) then
          dl(1,ie)  = dl(1,ie)  + cy(1)*wk(1,5)
          if (exi(ie) < -1d-6)
     .    dlp(1,ie) = dlp(1,ie) + cy(1)*wk(1,4)
        endif
   30 continue

      end
      subroutine xstrx0(nx,lxi,exi,a,vol,alat,n,qsq,yl,nlmx,
     .  cs,sn,cs2,sn2,cs3,sn3,dl,dlp)
C- Reciprocal part of reduced strx for Q=0
      implicit none
      integer lxi(1),nx,nlmx,n
      double precision alat,cof,e,gam,pi,vol
      double precision exi(1),a,qsq(n),yl(n,1),
     .  cs(n),sn(n),cs2(n),sn2(n),cs3(n),sn3(n),dl(nlmx,1),dlp(nlmx,1),
     .  c,s,xx,c2,c3,s2,s3
      integer i,ie,ilm,l,lx,m,nx1

      gam = 0.25d0/a**2
      pi = 4*datan(1d0)

C --- For each exi ne 0 and ilm do ---
      nx1 = nx
      if (dabs(exi(nx)) < 1d-6) nx1 = nx-1
      do  20  ie = 1, nx1
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
            dl(ilm,ie) = dl(ilm,ie) + cof*c2
            call dpdot(yl(1,ilm),cs3,n,c3)
            dlp(ilm,ie) = dlp(ilm,ie) + cof*(gam*c2 - c3)
   24     continue
          if (l+1 <= lx) then
            do  26  m = -l-2, l
              ilm = ilm+1
              call dpdot(yl(1,ilm),sn2,n,s2)
              dl(ilm,ie) = dl(ilm,ie) + cof*s2
              call dpdot(yl(1,ilm),sn3,n,s3)
              dlp(ilm,ie) = dlp(ilm,ie) + cof*(gam*s2 - s3)
   26       continue
          endif
   20 continue
      if (nx1 == nx) return

C --- Case e = 0 ---
      lx = lxi(nx)
      do  122  i = 2, n
        xx = -1/qsq(i)
        cs2(i) = cs(i)*xx
        sn2(i) = sn(i)*xx
  122 continue
      cof = -4d0*pi/vol
      ilm = 0
      do  120  l = 0, lx, 2
        cof = -cof
        do  124  m = -l, l
          ilm = ilm+1
          call dpdot(yl(2,ilm),cs2(2),n-1,c2)
          dl(ilm,ie) = dl(ilm,ie) + cof*c2
  124   continue
        if (l+1 <= lx) then
          do  126  m = -l-2, l
            ilm = ilm+1
            call dpdot(yl(2,ilm),sn2(2),n-1,s2)
            dl(ilm,ie) = dl(ilm,ie) + cof*s2
  126     continue
        endif
  120 continue
      dl(1,nx) = dl(1,nx) - 4d0*pi/vol*gam
      end
C      subroutine fmain
C      implicit none
C      integer nkdmx,nkrmx,lmxx,nlmx,nxx
C      parameter (nkdmx=1000, nkrmx=1000, lmxx=8,nlmx=(lmxx+1)**2, nxx=6)
C      double precision rlat(3,nkdmx), dlat(3,nkrmx),
C     .  wk(nkdmx*(lmxx+5)),yl(nkdmx*nlmx),
C     .  dl(nlmx,nxx),dlp(nlmx,nxx),dl0(nlmx),dlp0(nlmx),cy(300)
C      double precision rb(3,3),qb(3,3),rb0(3,3),qb0(3,3),tau(3),
C     .  vol,a,as,tol,alat,ecut,exi(10),e,dmx,dpmx
C      integer lmax,nkd,nkr,lxi(10),ie,nrx,nxi,ilm
Cc     data rb0 /0d0,.5d0,.5d0,.5d0,0d0,.5d0,.5d0,.5d0,0d0/
C      data rb0 /0d0,.47d0,.52d0,.5d0,0d0,.5d0,.8d0,.5d0,0d0/
C      data tau /.1d0,.4d0,.2d0/
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
C      nxi = 5
C      lxi(1) = 5
C      exi(1) = -3
C      lxi(2) = 7
C      exi(2) = -0.2
C      lxi(3) = 8
C      exi(3) = -8
C      lxi(4) = 1
C      exi(4) = -1.4
C      lxi(5) = 5
C      exi(5) = 0
C
CC      nxi = 1
CC      lxi(1) = 2
C      print *, 'tau=',tau
C      print *, 'set tau='
C      read(*,*) tau
C
C      lmax = -1
C      do  10  ie = 1, nxi
C   10 lmax = max(lmax,lxi(ie))
C      nrx = max(nkd,nkr)
C      call cstrx0(nxi,lxi,exi,ecut,tau,nrx,nlmx,lmax,wk,yl,
C     .  a,alat,rlat,nkr,dlat,nkd,vol,cy,dl,dlp)
C
C      print 334
C  334 format(/' Comparison of cstrx0 with rcsk0:'/
C     .  11x,'dl',12x,'error',11x,'dlp',12x,'error')
C      do  20  ie = 1, nxi
C        e = exi(ie)
C        if (e == 0) then
C         call rcnsl0(tau,a,lmax,alat,rlat,nkr,dlat,nkd,vol,cy,dl0)
C         call dpzero(dlp(1,ie),nlmx)
C         call dpzero(dlp0,nlmx)
C        else
C         call rcsk0(e,tau,a,lmax,alat,rlat,nkr,dlat,nkd,vol,cy,dl0,dlp0)
C        endif
C        print *, ' --- e=', sngl(e), '  lxi=',lxi(ie)
C        dmx  = 0
C        dpmx = 0
C        do  24  ilm = 1, (lxi(ie)+1)**2
C        if (dabs(dl0(ilm))+dabs(dl(ilm,ie))+
C     .      dabs(dlp0(ilm))+dabs(dlp(ilm,ie)) > 1d-18)
C     .      print 333, ilm, dl0(ilm), dl(ilm,ie)-dl0(ilm),
C     .                      dlp0(ilm),dlp(ilm,ie)-dlp0(ilm)
C        dmx = max(dmx,dabs(dl(ilm,ie)-dl0(ilm)))
C        dpmx= max(dpmx,dabs(dlp(ilm,ie)-dlp0(ilm)))
C  333   format(i4,2(1pe18.10,1pe12.3),2x)
C   24 continue
C      print 335, dmx,dpmx
C  335 format(' max error',12x,1pe12.3,18x,1pe12.3)
C   20 continue
C
C      end
