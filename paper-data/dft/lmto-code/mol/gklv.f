      subroutine gklv(rsm,rsq,nrx,kmx,nr,lmax,kmax,g)
C- Make radial parts g_kl of G_kL, divided by r**l (vectorises).
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsq     :radius**2 (!)
Ci   rsm     :smoothing radius
Ci   nr      :number of points
Ci   lmax    :make g_kl for l=0...lmax
Ci   kmax    :make g_kl for k=0...kmax
Ci   nrx,kmx :leading dimensions of g
Co Outputs
Co   g     :radial part of generalized gaussians G_kL; see Remarks
Co         :radgkl produces g such that G_kL = g(k,l) r**l Y_L
Cr Remarks
Cr   Definition:  G_kL = (lap)**k Y_L(-grad) g0 = g_kl*r^l*Y_L
Cr   Energy factor beta = dexp(e/(4*a*a)) is _not_ included here.
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
Cr                          a^l (2l+1)!!
Cr      p_kl = phi_kl ------------------------ r^l
Cr                    (2a^2)^(k+l) (2k+2l+1)!!
Cr
Cu Updates
Cu 8 Aug 2006 (S. Lozovi) First written based on radgkl(subs) and gausr(mol)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nrx,kmx,nr,lmax,kmax
      double precision rsq(nr),rsm,g(nrx,0:kmx,0:lmax)
C ... Local parameters
      integer l,k,ir
      double precision pi,a2,ta2,x,y,c,cl

      if (kmax < 0 .or. lmax < 0) return
      call tcn('gklv')

      pi = 4d0*datan(1d0)
      a2 = 1d0/rsm**2
      ta2 = 2d0*a2

C --- Do l = k = 0 explicitly ---
      c = (a2/pi)**1.5d0
      do ir = 1, nr
        g(ir,0,0) = c * dexp(-a2*rsq(ir))
      enddo

C --- Do explicitly for k=0,1; See Eqs. 5.14 and 5.15 ---
      do  l = 1, lmax
        do  ir = 1, nr
          g(ir,0,l) = ta2 * g(ir,0,l-1)
        enddo
      enddo

      if (kmax == 0) then
        call tcx('gklv')
        return
      endif

      do  l = 0, lmax
        c = ta2**(l+1)
        cl = 2*l + 3d0
        do  ir = 1, nr
          g(ir,1,l) = c*(ta2*rsq(ir)-cl) * g(ir,0,0)
        enddo
      enddo

      if (kmax == 1) then
        call tcx('gklv')
        return
      endif

C --- Recursion for higher k; see Eq. 5.19 ---
      do  l = 0, lmax
        do  k = 2, kmax
          x = 2*(k-1)*(2*k+2*l-1)
          y = 4*k+2*l-1
          do  ir = 1, nr
            g(ir,k,l) = ta2*((ta2*rsq(ir)-y)*g(ir,k-1,l) -
     .                  x*ta2*g(ir,k-2,l))
          enddo
        enddo
      enddo

      call tcx('gklv')
      end
