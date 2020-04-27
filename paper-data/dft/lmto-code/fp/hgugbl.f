      subroutine hgugbl(p1,p2,rsm1,rsm2,e1,nlm1,nlm2,ndim1,ndim2,s_lat,
     .  s,ds)
C- Estat energy integrals between Bloch Hankels and gaussians, and grad
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  vol alat plat qlat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg indxcg jcg cy qlv dlv
Cio    Passed to:  hhigbl phhigb hklbl gklbl fklbl hsmbl
Ci Inputs
Ci   p1    :first center
Ci   p2    :second center
Ci   rsm1  :smoothing radius of Hankels at p1
Ci   rsm2  :smoothing radius of gaussians at p2
Ci   e1    :energy  of Hankels at p1
Ci   nlm1  :L-max for  Hankels at p1
Ci   nlm2  :L-max for  gaussians at p2
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Co Outputs
Co   s     :integrals between Bloch Hankels and gaussians
Co   ds    :gradient of s; see Remarks
Cr Remarks
Cr   Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   22 Apr 00 Adapted from nfp hgug_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nlm1,nlm2,ndim1,ndim2
      double precision rsm1,rsm2,p1(3),p2(3),e1
      double complex s(ndim1,ndim2),ds(ndim1,ndim2,3)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer kmax,kdim,ilm2,ilm1
      double precision q(3),e2,vol
      data q /0d0,0d0,0d0/

      vol = s_lat%vol
      kmax = 0
      kdim = 0
      e2 = 0d0

      call hhigbl(0,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,ndim1,ndim2,
     .   kdim,s_lat%cg,s_lat%indxcg,s_lat%jcg,s_lat%cy,s_lat,s,ds)

      do  ilm2 = 1, nlm2
        do  ilm1 = 1, nlm1
          s(ilm1,ilm2) = 2d0*s(ilm1,ilm2)
          ds(ilm1,ilm2,1) = 2d0*ds(ilm1,ilm2,1)
          ds(ilm1,ilm2,2) = 2d0*ds(ilm1,ilm2,2)
          ds(ilm1,ilm2,3) = 2d0*ds(ilm1,ilm2,3)
        enddo
      enddo

      end
