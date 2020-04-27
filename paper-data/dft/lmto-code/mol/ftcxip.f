      subroutine ftcxip(lx1,ex1,nx1,lx2,ex2,nx2,rmt1,rmt2,
     .  rsmp1,rsmp2,rsm1,rsm2,lp1,ep1,lp2,ep2,cx,d,lscal,
     .  lmax,nx,lx,zc1,zc2,np,xp,zp,wp,ph1,ph2,xi)
C-  allocate and tabulate phi's and xi's on mesh
      implicit real*8 (a-h,p-z), integer (o)
      dimension cx(1),ex1(1),ex2(1),xp(1),zp(1),lx1(1),lx2(1),
     .   lx(20),ph1(1),ph2(1),xi(1)
      real w(1)
      common /w/ w
      np1=1
      np2=1
      nlm = ((lmax+1)*(lmax+2))/2
      nlm1 = ((lp1+1)*(lp1+2))/2
      nlm2 = ((lp2+1)*(lp2+2))/2
      if (iprint() >= 60) call tm('start smoothed h')
      call defrr(owk,    np)
      call defrr(oidx,   np*2)
      call defrr(oyl,    np*max(6,nlm,nlm1,nlm2))
      stop 'ftxcip update call to hyfxi'
C      call hyfxi(np,xp,zp,wp,lscal,lmax,nlm,rmt1,rmt2,
C     .  zc1,nx1,lx1,ex1,lp1,nlm1,rsmp1,rsm1,lp1,ep1,np1,
C     .  zc2,nx2,lx2,ex2,lp2,nlm2,rsmp2,rsm2,lp2,ep2,np2,
C     .  cx,w(oidx),w(owk),w(oyl),ph1,ph2,nx,lx,xi)

      call rlse(owk)
      end
