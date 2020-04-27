      subroutine vecpkl(r,rsm,nr,kmax,lmax,nrx,k0,wk,lrl,p,gp)
C- Vector of p_kl polynomials, or r^l p_kl
C ----------------------------------------------------------------------
Ci Inputs
Ci   r     :vector of points
Ci   rsm   :smoothing radius
Ci   nr    :number of points
Ci   kmax  :make p for k=0..kmax
Ci   lmax  :make p for l=0..lmax
Ci   nrx   :leading dimension of p
Ci   k0    :second dimension of p
Ci   wk    :work array of length nr
Ci   lrl   :if 1s digit = 0, returns p_kl; otherwise returns r^l p_kl
Ci         :if 10s digit nonzero, returns gp; otherwise gp is not touched.
Co Outputs
Co   p     :radial part of spherical polynomials P_kL; see Remarks
Co   gp    :radial derivative of p from l=0..lmax-1 (depending on lrl).
Cr Remarks
Cr   The p_kl = are defined in Eq. 12.7 J. Math. Phys. 39, 3393 (1988).
Cr   Note: Methfessel's gamma, a, and rs are related as gamma = 1/4/a^2 = rs^2/4
Cr   They are proportional to polynomials phi_kl  [Eq 12.7]
Cr      p_kl(r) := a**l (2l+1)!! / (2a**2)^(k+l) / (2k+2l+1)!! [r^l phi_kl(r)]
Cr   which are related to the generalized Laguerre polynomials  [Eq 12.3]
Cr      phi_kl(r) = (-)^l * (2a^2)^(k+l) * 2^k * k! * L_k^(l+1/2)(a^2 r^r)
Cr   The generalized Gaussians G_kL can be written in terms of the phi_kL as [Eq 12.2]
Cr      G_kL(a;r) = (a^2/pi)^(3/2) * exp(-a^2 r^2) * [phi_kl(r) r^l Y_L(rhat)]
Cr   Using the relation G_pL := (nabla)^p G_L [Eq 5.9] we have
Cr     [(nabla) G_kL]/G_kL = G_k+1,L/G_kL = phi_k+1,l / phi_kl = p_k+1,l / p_kl
Cr
Cr   P_kL = p_kl*Y_L are polyonomials orthogonal in the following sense [Eq 12.8]
Cr                                           (4a^2)^k a^l k! (2l+1)!!
Cr     int P_kL G_k'L' = delta_kk'*delta_ll'  ----------------------
Cr                                                     4pi
Cr   For any function F which can be integrated with the G_kL, viz
Cr     I_kL = int G_kL(a;r) F(r),
Cr   The P_kL can be uses as a basis to expand F.  Write
Cr     F(r) = f_L(r) r^l Y_L(rhat)       [Eq 12.10]
Cr   Then
Cr     f_L(r) r^l = sum_k a_kL p_kL(r)   [Eq 12.15]
Cr   with
Cr     a_kL = 4 pi [ (2a^2)^k * a^l 2^k * k! (2l+1)!!]^-1 I_kL    [Eq 12.16]
Cr   The a_kL that expand a smoothed Hankel function around a remote site:
Cr     H_L'(a;E,r-R') = sum_kL a_kL(a,E,R-R') P_kL(a,r-R)
Cr   are given in [Eq. 12.18].
Cr
Cr  *Practical calculation of the p_kl
Cr   Combining eqns 12.7 and 5.19 in that paper, we obtain
Cr    p_kl = a**l / (2a**2)^(k+l) (2l+1)!! / (2k+2l+1)!! phi_kl
Cr    p_0l = a**l
Cr    p_1l = a**l (2*(ar)**2/(2l+3) - 1)
Cr    p_kl = [(2*(ar)**2 - (4k+2l-1))p_k-1,l - 2(k-1)p_k-2,l] / (2k+2l+1)
Cu Updates
Cu   22 Aug 01 bug fix for gp when kmax=0
Cu   25 Jan 00 veckl generates gp as well as p.
C ----------------------------------------------------------------------
      implicit none
      integer nr,kmax,lmax,nrx,k0,lrl
      double precision r(nrx),wk(nr),rsm,p(nrx,0:k0,0:*),
     .  gp(nrx,0:k0,0:*)
      integer i,l,k
      double precision a,xx,xx2,xx3

      if (kmax < 0 .or. lmax < 0) return
      if (kmax > k0) call rx('vecpkl: kmax gt k0')
      if (rsm <= 0) call rx('vecpkl: rsm <= 0')
      a = 1d0/rsm

C --- Set wk = 2*a**2*r**2 ---
      xx = 2*a*a
      do  i = 1, nr
        wk(i) = xx*r(i)**2
      enddo

C --- Do explicitly for k=0,1 ---
      do  l = 0, lmax
      xx = a**l
        do  i = 1, nr
          p(i,0,l) = xx
        enddo
      enddo

      if (kmax > 0) then
        do  l = 0, lmax
          xx = a**l
          xx2 = 1/dble(2*l+3)
          do  i = 1, nr
            p(i,1,l) = xx*(wk(i)*xx2-1d0)
          enddo
        enddo
      endif

C --- Recursion for higher k ---
      do  k = 2, kmax
        xx3 = 2*(k-1)
        do  l = 0, lmax
          xx2 = (4*k+2*l-1)
          xx = 1/dble(2*k+2*l+1)
          do  i = 1, nr
            p(i,k,l) = xx*((wk(i)-xx2)*p(i,k-1,l) - xx3*p(i,k-2,l))
          enddo
        enddo
      enddo

C --- Radial derivative of p ---
      if (mod(lrl/10,10) /= 0) then

C  ... Set wk = 2*a**2*r**2
        xx = 2*a*a
        do  i = 1, nr
          wk(i) = xx*r(i)
        enddo

        do  k = 0, kmax
          do  l = 0, lmax-1
            xx2 = dble(2*k+2*l+3)/(a*(2*l+3))
            do  i = 1, nr
              gp(i,k,l) = wk(i)*(p(i,k,l) - xx2*p(i,k,l+1))
            enddo
          enddo
        enddo

      endif

C --- Scale by r^l if lrl nonzero ---
      if (mod(lrl,10) == 0) return

      do  i = 1, nr
        wk(i) = 1
      enddo
      do  l = 1, lmax

C   ... gP scales as  r*l gP +  l*r^l-1 P
        if (mod(lrl/10,10) /= 0 .and. l < lmax) then
          do  k = 0, kmax
            do  i = 1, nr
              gp(i,k,l) = wk(i)*r(i)*gp(i,k,l)+l*wk(i)*p(i,k,l)
            enddo
          enddo
        endif

        do  i = 1, nr
          wk(i) = wk(i)*r(i)
        enddo
        do  k = 0, kmax
          do  i = 1, nr
            p(i,k,l) = p(i,k,l)*wk(i)
          enddo
        enddo

      enddo

      end
