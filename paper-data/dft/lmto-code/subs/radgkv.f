      subroutine radgkv(nrx,kmx,r2,rsm,nr,kmax,lmax,gkl)
C- Make radial parts g_kl of G_kL divided by r**l for a vector of points
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrx   :max number of points r at which g_kl are evaluated;
Ci          leading dimension of gkl
Ci   kmx   :max power of Laplacian; second dimension of gkl
Ci   r2    :radius**2
Ci   rsm   :smoothing radius
Ci   nr    :actual number of points at which the g_kl are evaluated
Ci   kmax  :make g_kl for k=0...kmax
Ci   lmax  :make g_kl for l=0...lmax
Co Outputs
Co   gkl   :radial part of generalized gaussians G_kL; see Remarks
Co         :radgkv produces g such that G_kL(r) = gkl(r;k,l) r**l Y_L
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
Cr   radgkv is just a vectorized version of radgkl except that on input
Cr   it takes r**2 rather than r
Cb Bugs
Cu Updates
Cu   15 Jan 09 (S. Lozovoi) created from radgkl.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: nrx,kmx,nr,kmax,lmax
      double precision, intent(in) :: r2(nr),rsm
      double precision, intent(out) :: gkl(nrx,0:kmx,0:lmax)
C ... Local parameters
      integer il,ik
      double precision pi,a2,ta2,x,y,x0

      if(kmax. lt. 0) call rxi(' radgkv: negative kmax = ',kmax)
      if(lmax. lt. 0) call rxi(' radgkv: negative lmax = ',lmax)
      if(kmax. gt. kmx) call rx(' radgkv: kmax exceeds kmx')

      pi = 4d0*datan(1d0)
      a2 = 1d0/(rsm*rsm)
      ta2 = 2d0*a2
      x0 = (a2/pi)**1.5d0

C --- Do explicitly for k=0,1; See Eqs. 5.14 and 5.15 ---
      do  il = 0, lmax
        x = ta2**il*x0
        gkl(1:nr,0,il) = x * dexp(-a2*r2(1:nr))
      enddo
      if (kmax == 0) return

      do  il = 0, lmax
        y = dble(3+2*il)
        gkl(1:nr,1,il) = ta2*(ta2*r2(1:nr)-y) * gkl(1:nr,0,il)
      enddo
      if (kmax == 1) return

C --- Recursion for higher k; see Eq. 5.19 ---
      do  il = 0, lmax
        do  ik = 2, kmax
          x = dble(2*(ik-1)*(2*ik+2*il-1))*ta2
          y = dble(4*ik+2*il-1)
          gkl(1:nr,ik,il) = ta2*(ta2*r2(1:nr)*gkl(1:nr,ik-1,il) -
     .                  y*gkl(1:nr,ik-1,il) - x*gkl(1:nr,ik-2,il))
        enddo
      enddo

      end
