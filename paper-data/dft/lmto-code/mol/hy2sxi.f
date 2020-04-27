      subroutine hy2sxi(np,xyz,nlm,
     .  tau1,nx1,lx1,ex1,nlmp1,rsmp1,rsm1,lph1,eph1,phr1,ioff1,nr1,nth1,
     .  tau2,nx2,lx2,ex2,nlmp2,rsmp2,rsm2,lph2,eph2,phr2,ioff2,nr2,nth2,
     .  cy,idx,wk,yl,ph1,ph2,nx,lx,xi)
C- Tabulates smoothed hankels xi's and phi's on a mesh
C  phi's are smoothed hankels with smoothing radii rsmp, except radial
C  functions for points noff+1...ns are taken from phr.
      implicit none
      integer nx1,lx1(1),lph1,nlm,nlmp1,nx2,lx2(1),lph2,nlmp2,lx(1),nx,
     .  np,nr1,nr2,ioff1,ioff2,nth1,nth2
      double precision xyz(np,3),phr1(nr1,0:1),phr2(nr2,0:1),
     .  cy(1),wk(1),idx(1),yl(np,1),tol,zc,xi(np,nlm,1),
     .  ex1(1),eph1,rsmp1,rsm1,ph1(np,1),tau1(3),
     .  ex2(1),eph2,rsmp2,rsm2,ph2(np,1),tau2(3)
      integer ie,ip,l,m,lm,ir,ith

C --- Concatenate lxi for xi's centered on two sites ---
      nx = nx1+nx2
      do  5  ie = 1, nx1
    5 lx(ie) = lx1(ie)
      do  6  ie = 1, nx2
    6 lx(ie+nx1) = lx2(ie)

      tol = 1d-10
      do  10  m = 1, 3
      zc = -tau1(m)
      do  10  ip = 1, np
   10 xyz(ip,m) = xyz(ip,m) + zc

C --- Smoothed Hankels xi and phi centered at origin ---
      call solhsr(rsm1,nlm,nx1,lx1,ex1,xyz(1,1),xyz(1,2),xyz(1,3),np,np,
     .  tol,cy,idx,wk,yl,xi)
      call solhsr(rsmp1,nlmp1,1,lph1,eph1,xyz(1,1),xyz(1,2),xyz(1,3),
     .  np,np,tol,cy,idx,wk,yl,ph1)

C --- Augmented functions ph1 for points ioff1+1 ... ioff1+nr1*nth1 ---
      lm = 0
      do  20  l = 0, lph1
      do  20  m = -l, l
        lm = lm+1
        ip = ioff1
        do  22  ir  = 1, nr1
        do  22  ith = 1, nth1
        ip = ip+1

C        if (ir >= 1) print 345,
C     .    ir,ith,lm,xyz(ip,1),xyz(ip,2),xyz(ip,3),
C     .    yl(ip,lm),phr1(ir,l), yl(ip,lm)*phr1(ir,l),
C     .    ph1(ip,lm)-yl(ip,lm)*phr1(ir,l)
C  345   format(3i3,7f10.5)
C       print *, ir,ith,lm,ph1(ip,lm),ph1(ip,lm)-yl(ip,lm)*phr1(ir,l)
   22   ph1(ip,lm) = yl(ip,lm)*phr1(ir,l)
   20 continue

      do  25  m = 1, 3
      zc = tau1(m) - tau2(m)
      do  25  ip = 1, np
   25 xyz(ip,m) = xyz(ip,m) + zc

C --- xi and phi centered off origin ---
      call solhsr(rsm2,nlm,nx2,lx2,ex2,xyz(1,1),xyz(1,2),xyz(1,3),np,np,
     .  tol,cy,idx,wk,yl,xi(1,1,nx1+1))
      call solhsr(rsmp2,nlmp2,1,lph2,eph2,xyz(1,1),xyz(1,2),xyz(1,3),
     .  np,np,tol,cy,idx,wk,yl,ph2)

C --- Augmented functions ph2 for points ioff2+1 ... ioff2+nr2*nth2 ---
      lm = 0
      do  30  l = 0, lph2
      do  30  m = -l, l
        lm = lm+1
        ip = ioff2
        do  32  ir  = 1, nr2
        do  32  ith = 1, nth2
        ip = ip+1
C       print *, ir,ith,lm,dsqrt(xyz(ip,1)**2+xyz(ip,2)**2+xyz(ip,3)**2)
C    .    ,sngl(ph2(ip,lm)),sngl(yl(ip,lm)*phr2(ir,l)),sngl(phr2(ir,l))
   32   ph2(ip,lm) = yl(ip,lm)*phr2(ir,l)
   30 continue

C --- Restore xyz ---
      do  34  m = 1, 3
      zc = tau2(m)
      do  34  ip = 1, np
   34 xyz(ip,m) = xyz(ip,m) + zc

      end
