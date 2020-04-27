      subroutine radgkg(mode,r,rsml,kmax,lmax,n0,g,dg,ddg)
C- Make radial parts g_kl of G_kL, and optionally radial derivative
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 generate g only
Ci         :>0 generate
Ci             dg = dg/dr and
Ci             ddg = 1/r d^2 (r*gkl) / dr^2  - l(l+1)/r^2 gkl
Ci                 = radial part of laplacian of 3-dimensional g_kl YL
Ci         :10s digit
Ci         :0  return radial parts of G_kL
Ci         :>0 return radial parts of G_kL divided by r**l
Ci   r     :radius
Ci   rsml  : vector of l-dependent smoothing radii of Gaussians
Ci         : EITHER must be specified for 0..lmax
Ci         : OR     rsml(0) < 0. Implies rsml(l) = -rsml(0) for all l
Ci   kmax  :make g_kl for k=0...kmax
Ci   lmax  :make g_kl for l=0...lmax
Ci   n0    :leading dimension of g
Co Outputs
Co   g     :radial part of generalized gaussians G_kL; see Remarks
Co         :radgkg produces g such that
Co         :  G_kL = g(k,l) r**l Y_L  (10s digit mode = 0) or:
Co         :  G_kL = g(k,l) Y_L  (10s digit mode = 1)
Co   dg    :dg/dr (10s digit mode 1) or d(g/r^l)/dr (10s digit mode 0)
Co   ddg   :  G_k+l,L
Cr Remarks
Cr   Definition:  G_kL = (lap)**k YL(-grad) g0 = g_kl*r^l*Y_L
Cr   Energy factor beta = dexp(e/(4*a*a)) is not included here.
Cr   See J. Math. Phys. {\bf 39},3393 (1998), Section V:
Cr
Cr     G_kL  = phi_kl (a^2/pi)^(3/2) exp(-a^2 r^2) r^l Y_L
Cr           = phi_kl r^l Y_L g_00
Cr     G_0L  = (2a^2)^l r^l Y_L g_00
Cr
Cr   where (Eq. 12.3)
Cr
Cr    phi_kl = (-1)^k (2a^2)^(k+l) 2^k k! L_k^(l+1/2)(a^2 r^2)
Cr    L_k    = generalized Laguerre polynomial
Cr    phi_0l = (2a^2)^l
Cr
Cr   also
Cr
Cr                          a^l (2l+1)!!
Cr      p_kl = phi_kl ------------------------ r^l
Cr                    (2a^2)^(k+l) (2k+2l+1)!!
Cr
Cr Gradient from general properties of functions yl(grad) f
Cr   g_kl'  = l/r g_kl - g_kl+1
Cr Alternative relation for gradient derived from solghg
Cr   g0' = - cz(l+1)(2*l+3)/(2*l+1)g(0,l+1) - cz(l)(2*l+1)/(2*l+3)g(1,l-1)
Cr where cz(l) = l/[(2*l-1)*(2*l+1)]
Cr Rewrite as
Cr   g0' = -1/(2l+1)[(l+1)g0_l+1 - lg0_l-1]
Cr
Cr Laplacian:
Cr   From the relation of the solid G_kl, lap G_0L = G_1L
Cr   we obtain radial part of Laplacian,
Cr     1/r d^2 (r*g0_l) / dr^2 = g1_l + l(l+1)g0_l/r^2
Cr   Laplacian of 3-dimensional G_0L is
Cr     lap(g0_l YL) = g1_l YL
Cr                  = lap(g0_l) YL + g0_l lap(YL)
Cr                  = [1/r d^2 (r*g0_l) / dr^2 - l(l+1)g0_l/r^2] YL
Cu Updates
Cu   28 Apr 11 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,kmax,lmax,n0
      double precision r,rsml(0:lmax)
      double precision g(0:n0,0:lmax),dg(0:n0,0:lmax),ddg(0:n0,0:lmax)
C ... Local parameters
      logical lrsm0
      integer l,k,mode0,mode1
      double precision gkll(0:n0+1,0:lmax+1),rl,facl,rsm,rsx,tol
      parameter (tol=1d-15)
C     double precision x,y

      if (kmax < 0 .or. lmax < 0) return
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
C     Negative 1st smoothing radius => constant rsm
      lrsm0 = rsml(0) <= 0

C      pi = 4d0*datan(1d0)
C      a = 1d0/rsm
C      ta2 = 2*a*a
C      g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0

      rsm = dabs(rsml(0))
      rsx = -9999
C ... G__kL/(r**l Y_L), Y_L = spherical harmonic
      do  l = lmax, 0, -1
        if (.not. lrsm0) rsm = abs(rsml(l))
        if (dabs(rsm-rsx) > tol) then
          if (mode0 == 0) then
            call radgkl(r,rsm,kmax,l,n0+1,gkll)
          else
            call radgkl(r,rsm,kmax+1,l+1,n0+1,gkll)
          endif
          rsx = rsm
        endif

C   --- Copy to g and make dg, possibly scaling by r^l ---
        if (mode1 == 0) then
          rl = r**l
          facl = l
        else
          rl = 1
          facl = 0
        endif
        do  k = 0, kmax
          g(k,l) = rl*gkll(k,l)
        enddo
        if (mode0 /= 0) then
          do  k = 0, kmax
            dg(k,l)  = g(k,l)*facl/r - r*rl*gkll(k,l+1)
            ddg(k,l) = rl*gkll(k+1,l)
          enddo
        endif
      enddo

C ...  This section below for alternate formula for dg, k=0 only
C      if (mode0 == 0) return
C
CC ... scale gkll by r^l
C      rl = 1
C      do  l = 0, lmax1
C        do  k = 0, kmax+1
C          gkll(k,l) = rl*gkll(k,l)
C        enddo
C        rl = rl*r
C      enddo
C
CC ... grad g0l = -1/(2l+1)[(l+1)g_l+1 - lg_l-1]
C      dg(0) = -gkll(0,1)
C      do  l = 1, lmax
CC       cz = (l+1)/dsqrt(dble((2*l+1)*(2*l+3)))
CC       dg(l) = -(gkll(0,l+1)*cz*dsqrt((2*l+3)/dble(2*l+1)))
C        dg(l) = -gkll(0,l+1)*(l+1)/dble(2*l+1)
C
CC       cz = l/dsqrt(dble((2*l-1)*(2*l+1)))
CC       dg(l) = dg(l) -cz*gkll(1,l-1) / dsqrt((2*l+1)/dble(2*l-1))
C        dg(l) = dg(l) - l*gkll(1,l-1)/dble(2*l+1)
C
C      enddo
C
C      if (mode1 /= 0) then
C        rl = 1
C        do  l = 1, lmax1
C          rl = rl*r
C          dg(l) = dg(l)/rl - l/r*g(0,l)
C        enddo
C      endif

      end
C      subroutine fmain
C      implicit none
C      integer n0,kmax,lmax,k,l
C      parameter (n0=4,kmax=3,lmax=4)
C      double precision gkl0(0:n0,0:lmax),gkl(0:n0,0:lmax),
C     .  dgkl(0:n0,0:lmax),ddgkl(0:n0,0:lmax)
C      double precision dgkl0(0:lmax),rsml(0:lmax)
C      double precision r,rsm,pi,a,g0,lapgl0(0:lmax),dg0,ddg0
C
C      r = 1.1d0
C      rsm = 1.198d0
C
CC     r = 1; rsm=1
C      r=2.3d0 ; rsm=2*r/3
CC     r = 2*r
C
C      do  l = 0, lmax
C        rsml(l) = rsm*(1-1*l/20d0)
C      enddo
C
C      print *, 'r=',r,'   rsm=',rsm
C
C      pi = 4d0*datan(1d0)
C      a = 1d0/rsm
C      g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0
C
C  333 format(1x,a6,5f16.10)
C
C      call radgkl(r,rsm,kmax,lmax,n0,gkl0)
C      call radgkg(10,r,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C
C      print *, 'gkl/r^l'
C      do  k = 0, kmax
C        print 333, 'mode 1',(gkl(k,l), l=0,lmax)
C      enddo
C
C      print *, 'compare gkl/r^l to radgkl'
C      do  k = 0, kmax
C        print 333,'diff',(gkl(k,l)-gkl0(k,l), l=0,lmax)
C      enddo
C
C      print *, 'l-dependent gkl'
C      call radgkg(10,r,rsml,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      do  k = 0, kmax
C        print 333, 'l-dep',(gkl(k,l), l=0,lmax)
C      enddo
C
C      call radgkg(00,r,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print *
C      print *, 'gkl'
C      do  k = 0, kmax
C       print 333, 'mode 0',(gkl(k,l), l=0,lmax)
C      enddo
C
C      print *
C      print *, 'compare gkl to radgkl'
C      do  k = 0, kmax
C       print 333, 'diff',(gkl(k,l)-gkl0(k,l)*r**l, l=0,lmax)
C      enddo
C
C      print *
C      print *, 'd(g0/r^l)/dr, compared to num diff'
C      call radgkg(10,r+1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      dgkl0 = gkl(0,:)
C      call radgkg(10,r-1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      dgkl0 = (dgkl0 - gkl(0,:))/2d-4
C      print 333, 'num', dgkl0
C      call radgkg(11,r,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'dg0', dgkl(0,:)
C      dg0 = -2*a*a*r*g0
C      print 333, 'l=0', dg0
C      print *, 'l-dependent d(g0/r^l)/dr'
C      call radgkg(11,r,rsml,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'l-dep', dgkl(0,:)
C
C      print *
C      print *, 'grad g0 = dg0/dr, compared to num diff'
C      call radgkg(0,r+1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      dgkl0 = gkl(0,:)
C      call radgkg(0,r-1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      dgkl0 = (dgkl0 - gkl(0,:))/2d-4
C      print 333, 'num', dgkl0
C      call radgkg(1,r,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'dg0', dgkl(0,:)
C      dg0 = -2*a*a*r*g0
C      print 333, 'l=0', dg0
C      call radgkg(1,r,rsml,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'l-dep', dgkl(0,:)
C
C      if (kmax > 0) then
C      dgkl = -9999
C      print *
C      print *, 'grad g1 = dg1/dr, compared to num diff'
C      call radgkg(0,r+1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      dgkl0 = gkl(1,:)
C      call radgkg(0,r-1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      dgkl0 = (dgkl0 - gkl(1,:))/2d-4
C      print 333, 'num', dgkl0
C      call radgkg(1,r,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'dg1', dgkl(1,:)
C      call radgkg(1,r,rsml,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'l-dep', dgkl(1,:)
C      endif
C
C      print *
C      print *, 'nabla^2 g0 = 1/r d^2(r*g0)/dr^2, compared to num diff'
C      call radgkg(1,r+1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
CC     print 333, 'xx', dgkl
C      lapgl0 = (r+1d-4)*gkl(0,:)
C      call radgkg(1,r-1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
CC     print 333, 'xx', dgkl(0,:)
C      lapgl0 = lapgl0 + (r-1d-4)*gkl(0,:)
C      call radgkg(1,r,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      lapgl0 = (lapgl0 - 2*(r)*gkl(0,:))/(1d-4)**2/r
CC      do  l = 1, lmax
CC        lapgl0(l) = lapgl0(l) - l*(l+1)*gkl(0,l)/r**2
CC      enddo
C      print 333, 'num', lapgl0
C      if (kmax > 0) then
C      print 333, 'ddg', (gkl(1,l)+l*(l+1)*gkl(0,l)/r**2, l=0,lmax)
C      endif
CC     1/r d^2(r.g0)/dr^2 = d^2(g0)/dr^2 + 2/r d(g0)/dr
C      dg0 = -2*a*a*r*g0
C      ddg0 = -2*a*a*g0 + (-2*a*a*r)**2*g0
C      ddg0 = ddg0 + 2/r*dg0
C      print 333, 'l=0', ddg0
C      print *, '... 3D laplacian'
C      print 333, 'num', (lapgl0(l)-l*(l+1)*gkl(0,l)/r**2, l=0,lmax)
C      if (kmax > 0) then
C      print 333, 'g1l', gkl(1,:)
C      endif
C      print 333, 'ddg', ddgkl(0,:)
C      print *, '... - 3D laplacian / g'
C      if (kmax > 0) then
C      print 333, 'k.e.', (-gkl(1,l)/gkl(0,l), l=0,lmax)
C      endif
C      call radgkg(1,r,rsml,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'l-dep', (-ddgkl(0,l)/gkl(0,l), l=0,lmax)
C
CC     print *, '!! scale rsm by .85 to check l-dep matches f orbital'
CC     rsm = rsm * .85d0
C
C      if (kmax > 0) then
C      print *
C      print *, 'nabla^2 g1 = 1/r d^2(r*g1)/dr^2, compared to num diff'
C      call radgkg(1,r+1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      lapgl0 = (r+1d-4)*gkl(1,:)
C      call radgkg(1,r-1d-4,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      lapgl0 = lapgl0 + (r-1d-4)*gkl(1,:)
C      call radgkg(1,r,-rsm,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      lapgl0 = (lapgl0 - 2*(r)*gkl(1,:))/(1d-4)**2/r
C      print 333, 'num', lapgl0
C      if (kmax > 1) then
C      print 333, 'ddg', (gkl(2,l)+l*(l+1)*gkl(1,l)/r**2, l=0,lmax)
C      endif
C      print *, '... 3D laplacian'
C      print 333, 'num', (lapgl0(l)-l*(l+1)*gkl(1,l)/r**2, l=0,lmax)
C      if (kmax > 1) then
C      print 333, 'g2l', (gkl(2,l), l=0,lmax)
C      endif
C      print 333, 'ddg', (ddgkl(1,l), l=0,lmax)
C      print *, '... - 3D laplacian / g'
C      if (kmax > 1) then
C      print 333, 'k.e.', (-gkl(2,l)/gkl(1,l), l=0,lmax)
C      endif
C      call radgkg(1,r,rsml,kmax,lmax,n0,gkl,dgkl,ddgkl)
C      print 333, 'l-dep', (-ddgkl(1,l)/gkl(1,l), l=0,lmax)
C      endif
C
C      end
