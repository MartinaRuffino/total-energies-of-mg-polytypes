      subroutine ggugbl(p1,p2,rsm1,rsm2,nlm1,nlm2,ndim1,ndim2,
     .  s_lat,s,ds)
C- Estatic energy integrals between Bloch gaussians, and gradients.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat awald tol vol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg indxcg jcg cy qlv dlv
Cio    Passed to:  gfigbl fklbl gklbl
Ci Inputs
Ci   p1    :first center
Ci   p2    :second center
Ci   rsm1  :smoothing radius of Gaussians at p1
Ci   rsm2  :smoothing radius of Gaussians at p2
Ci   e1    :energy  of Gaussians at p1
Ci   e2    :energy  of Gaussians at p2
Ci   nlm1  :L-max for  Gaussians at p1
Ci   nlm2  :L-max for  Gaussians at p2
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Co Outputs
Co   s     :integrals between Bloch Gaussians
Co   ds    :gradient of s; see Remarks
Cr Remarks
Cr   Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   22 Apr 00 Adapted from nfp ggug_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nlm1,nlm2,ndim1,ndim2
      double precision rsm1,rsm2,p1(3),p2(3)
      double complex s(ndim1,ndim2),ds(ndim1,ndim2,3)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer kmax,kdim,ilm2,ilm1

      kmax = 0
      kdim = 0
      call gfigbl(p1,p2,rsm1,rsm2,nlm1,nlm2,kmax,ndim1,ndim2,kdim,
     .   s_lat%cg,s_lat%indxcg,s_lat%jcg,s_lat%cy,s_lat,s,ds)

      do  ilm2 = 1, nlm2
        do  ilm1 = 1, nlm1
          s(ilm1,ilm2) = 2d0*s(ilm1,ilm2)
          ds(ilm1,ilm2,1) = 2d0*ds(ilm1,ilm2,1)
          ds(ilm1,ilm2,2) = 2d0*ds(ilm1,ilm2,2)
          ds(ilm1,ilm2,3) = 2d0*ds(ilm1,ilm2,3)
        enddo
      enddo

      end
