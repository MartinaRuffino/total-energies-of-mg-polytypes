      subroutine ghigbl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,
     .  ndim1,ndim2,k0,cg,indxcg,jcg,cy,s_lat,s,ds)
C- Block of integrals between smooth hankels and gaussians with some
c  power of the laplace operator, and gradients.
C ----------------------------------------------------------------------
Ci Inputs
Ci   pg    :Function it expansed is at pg; see Remarks
Ci   ph    :Function is centered at ph; see Remarks
Ci   q     :wave number for Bloch sum
Ci   rsmg  :smoothing radius of gaussian
Ci   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
Ci         :rsmh must be specified for 1..ll(nlmh)+1
Ci   eg    :gkL scaled by exp(e*rsm**2/4)
Ci   eh    :vector of l-dependent energies of smoothed Hankel
Ci         :eh must be specified for 1..ll(nlmh)+1
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   nlmh  :L-cutoff for smoothed Hankel functions
Ci   kmax  :polynomial cutoff
Ci   ndim1 :leading dimension of s,ds
Ci   ndim2 :second dimension of s,ds
Ci   k0    :third dimenson of s,ds
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals of gaussian and Hankels; see Remarks
Co   ds    :gradients of s; see Remarks
Cr Remarks
Cr   s(L,M,k) contains integral of G_L^*(r-pg) (lap)^k H_M(r-ph)
Cr  ds((L,M,k) contains grad s wrt ph; take negative for grad wrt pg.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   25 May 00 Made rsmh,eh l-dependent
Cu   25 May 00 Adapted from nfp ghig_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k0,kmax,ndim1,ndim2,nlmg,nlmh,jcg(*),indxcg(*)
      double precision rsmg,rsmh(*),eg,eh(*)
      double precision ph(3),pg(3),cg(*),cy(*),q(3)
      double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer icg,icg1,icg2,ii,ilg,ilh,ilm,ilm1,ilm2,indx,ip,jlm,k,ktop,
     .  ktop0,l1,l2,lg,lh,ll,lm,lmaxg,lmaxh,lmaxx,m,nlm0,nlmx
      double precision ee,fac,gamg,gamh,rsmx,dr(3),e,rsm
      complex(8),allocatable:: hkl(:,:),dhkl(:,:,:)

      if (nlmh == 0 .or. nlmg == 0) return

C ... rsmh- and eh- independent setup
      do  m = 1, 3
        dr(m) = pg(m)-ph(m)
      enddo
      lmaxh = ll(nlmh)
      lmaxg = ll(nlmg)
      lmaxx = lmaxg+lmaxh
      nlmx = (lmaxx+1)**2
      ktop = max0(lmaxg,lmaxh)+kmax

      ktop0 = ktop+1
      nlm0  = (lmaxx+2)**2
      allocate( hkl(0:ktop0,nlm0),dhkl(0:ktop0,nlm0,3))

      do  k = 0, kmax
      do  jlm = 1, nlmh
      do  ilm = 1, nlmg
        s(ilm,jlm,k) = dcmplx(0d0,0d0)
        ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
        ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
        ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)
      enddo
      enddo
      enddo

C --- Loop over sequences of l with a common rsm,e ---
      l2 = -1
      do  l1 = 0, lmaxh
        if (l1 <= l2) cycle
        call gtbsl2(l1,lmaxh,eh,rsmh,l2)
        rsm  = rsmh(l1+1)
        e    = eh(l1+1)
        if (rsm <= 0 .or. e > 0) cycle
        ilm1 = l1**2+1
        ilm2 = (l2+1)**2
        lmaxx= lmaxg+l2
        nlmx = (lmaxx+1)**2
        gamh = 0.25d0*rsm*rsm
        gamg = 0.25d0*rsmg*rsmg
        rsmx = 2d0*dsqrt(gamg+gamh)
        ktop = max0(lmaxg,l2)+kmax
        call hklgbl(dr,rsmx,e,q,ktop,nlmx,ktop0,nlm0,cy,s_lat,
     .    hkl,dhkl)

C   ... Combine with Clebsch-Gordan coefficients
        ee = dexp(gamg*(eg-e))
        do  ilg = 1, nlmg
        lg = ll(ilg)
        do  ilh = ilm1, ilm2
          lh = ll(ilh)
          ii = max0(ilg,ilh)
          indx = (ii*(ii-1))/2 + min0(ilg,ilh)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  icg = icg1, icg2
            ilm = jcg(icg)
            lm = ll(ilm)
            k = (lg+lh-lm)/2
            fac = ee*(-1d0)**lg*cg(icg)
            do  ip = 0, kmax
              s(ilg,ilh,ip)    = s(ilg,ilh,ip)    + fac*hkl(k+ip,ilm)
              ds(ilg,ilh,ip,1) = ds(ilg,ilh,ip,1) + fac*dhkl(k+ip,ilm,1)
              ds(ilg,ilh,ip,2) = ds(ilg,ilh,ip,2) + fac*dhkl(k+ip,ilm,2)
              ds(ilg,ilh,ip,3) = ds(ilg,ilh,ip,3) + fac*dhkl(k+ip,ilm,3)
            enddo
          enddo
        enddo
        enddo
      enddo

      deallocate(hkl,dhkl)

      end
