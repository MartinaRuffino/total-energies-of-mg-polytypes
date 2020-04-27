      subroutine hxpbl(ph,pg,q,rsmh,rsmg,eh,kmax,nlmh,nlmg,k0,ndim,
     .  cg,indxcg,jcg,cy,s_lat,c)
C- Coefficients to expand smooth bloch hankels centered at ph
C  into a sum of polynomials P_kL centered at pg.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ph    :Function is centered at position ph; see Remarks
Ci   pg    :Function is expanded about position pg; see Remarks
Ci   q     :wave number for Bloch sum
Ci   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
Ci         :rsmh must be specified for 1..ll(nlmh)+1
Ci   rsmg  :smoothing radius of gaussian
Ci   eh    :vector of l-dependent energies of smoothed Hankel
Ci         :eh must be specified for 1..ll(nlmh)+1
Ci   kmax  :polynomial cutoff
Ci   nlmh  :L-cutoff for smoothed Hankel functions being expanded
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   k0    :leading dimension of coefficient array c
Ci   ndim  :second dimension of coefficient array c
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   c     :cofficients c(k,M,L); see Remarks
Co         :here k=0..kmax, M=1..nlmg, L=1..nlmh
Cr Remarks
Cr   Expansion is:  H_L(r-ph) = sum(kM) c(k,M,L) * P_kM(r-pg)
Cr   See J. Math. Phys. {\bf 39},3393 (1998), Eq. 12.17
Cr
Cr   As rsmg -> 0, expansion turns into a Taylor series of H_L.
Cr   As rsmg increases the error is spread out over a progressively
Cr   larger radius, thus reducing the error for larger r while
Cr   reducing accuracy at small r.
Cb Bugs
Cb   Doesn't properly handle case rsmh<0
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 May 00 Made rsmh,eh l-dependent
Cu   24 Apr 00 Adapted from nfp hxp_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k0,kmax,ndim,nlmg,nlmh,jcg(*),indxcg(*)
      double precision eh(*),rsmg,rsmh(*),ph(3),pg(3),cg(*),cy(*),q(3)
      double complex c(0:k0,ndim,nlmh)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer ndim1,ndim2,ktop0,ilmg,ilmh,k,l,ll,lmaxg,m,nm
      double precision a,dfact,eg,fac,factk,fpi
      complex(8),allocatable:: s(:,:,:)

      if (nlmg == 0 .or. nlmh == 0) return
      fpi = 16d0*datan(1d0)
C     Memory allocation
      ndim1 = nlmg
      ndim2 = nlmh
      ktop0 = kmax
      allocate(s(ndim1,ndim2,0:ktop0))

      if (kmax > ktop0) call rxi('hxpbl: increase ktop0, need',kmax)
      if (nlmg > ndim1) call rxi('hxpbl: increase ndim1, need',nlmg)
      if (nlmh > ndim2) call rxi('hxpbl: increase ndim2, need',nlmh)
      if (nlmg > ndim1) call rxi('hxpbl: increase ndim1, need',nlmg)

C ... Integrals of gaussians and smoothed Hankels
      eg = 0d0
      call ghibl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,ndim1,ndim2,cg,indxcg,jcg,cy,s_lat,s)

C ... Scale to get coefficients of the P_kL
      a = 1d0/rsmg
      lmaxg = ll(nlmg)

      ilmg = 0
      dfact = 1d0
      do  l = 0, lmaxg
        nm = 2*l+1
        do  m = 1, nm
          ilmg = ilmg+1
          factk = 1d0
          do  k = 0, kmax
            fac = fpi / ((4*a*a)**k * a**l * factk * dfact)
            do  ilmh = 1, nlmh
              c(k,ilmg,ilmh) = s(ilmg,ilmh,k)*fac
            enddo
            factk = factk*(k+1)
          enddo
        enddo
        dfact = dfact*(2*l+3)
      enddo

      deallocate(s)
      end
