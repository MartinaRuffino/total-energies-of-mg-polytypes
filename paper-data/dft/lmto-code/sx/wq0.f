      subroutine wq0(nn,p0,p1,p2,v0,v1,v2,w0,w1,w2)
C- Screened potential ws := v0*(1- p0*v0)^-1 for q -> 0
C Inputs
C   nn
C   v0,v1 (from vbare.f)
C   v2    (from madmtq.f)
C Outputs
C   w0,w1,w2
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nn
      double precision p0(nn,nn),p1(nn,nn),p2(nn,nn),v0(nn,nn),v1(nn,nn)
      double precision w0(nn,nn),w1(nn,nn),w2(nn,nn),eps,v2
C ... Dynamically allocated local arrays
      real(8), allocatable :: px(:),vx(:),w1l(:),w2l(:)
C ... Local parameters
      double precision q1,q2
      integer iprint,lgunit

      q1 = .005d0
      q2 = 2d0*q1
      allocate(px(nn*nn*2))
      allocate(vx(nn*nn*2))
      allocate(w1l(nn*nn*2))
      allocate(w2l(nn*nn*2))
      call pofq(nn,px,p0,p1,p2,q1)
      call vofq(nn,vx,v0,v1,v2,q1)
      call wstat(nn,px,vx,w1l)
      call pofq(nn,px,p0,p1,p2,q2)
      call vofq(nn,vx,v0,v1,v2,q2)
      call wstat(nn,px,vx,w2l)
      call wofq(nn,w1l,w2l,w0,w1,w2,q1,q2)
C      if (iprint() >= 50) then
C        call o1r(w0,nn,'  w0q0  ')
C        call o1r(w1,nn,'  w1q0 ')
C        call o1r(w2,nn,'  w2q0  ')
C      endif
      eps = v2/w2(1,1)
      if (iprint() >= 20) call awrit1(
     .  '%N WQ0: q=0 dielectric constant = %,1;1d',' ',80,lgunit(1),eps)
      write(lgunit(2),101) eps
  101 format(/'  eps=',f8.1)
      deallocate(px,vx,w1l,w2l)
      end

      subroutine pofq(nn,px,p0,p1,p2,q)
C Outputs
C   px
C ----------------------------------------------------------------------
      implicit none
      integer nn,i,j
      double precision q,p0(nn,nn),p1(nn,nn),p2(nn,nn),px(nn,nn,2)
      do  10  i = 1, nn
      do  10  j = 1, nn
        px(j,i,1) = p0(j,i) + q**2*p2(j,i)
        px(j,i,2) = q*p1(j,i)
   10 continue
      end

      subroutine vofq(nn,vx,v0,v1,v2,q)
C Inputs
C   nn,q
C   v0,v1,v2
C Outputs
C   vx
C ----------------------------------------------------------------------
      implicit none
      integer nn,i,j
      double precision v2,q,v0(nn,nn),v1(nn,nn),vx(nn,nn,2)
      do  10  i = 1, nn
      do  10  j = 1, nn
        vx(j,i,1) = v0(j,i) + v2/q**2
        vx(j,i,2) = v1(j,i)/q
   10 continue
      end

      subroutine wofq(nn,zw1,zw2,w0,w1,w2,q1,q2)
C Outputs
C   w0,w1,w2
C ----------------------------------------------------------------------
      implicit none
      integer nn,i,j
      double precision q1,q2,w0(nn,nn),w1(nn,nn),w2(nn,nn)
      double precision zw1(nn,nn,2),zw2(nn,nn,2),x

      x = 1d0/(1d0/q2**2 - 1d0/q1**2)
      do  10  i = 1, nn
      do  10  j = 1, nn
        w1(j,i) = zw1(j,i,2)*q1
        w2(j,i) = (zw2(j,i,1) - zw1(j,i,1))*x
        w0(j,i) = zw2(j,i,1) - w2(j,i)/q2**2
   10 continue
      end

