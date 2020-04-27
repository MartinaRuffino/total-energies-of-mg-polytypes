      subroutine hhibl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,
     .  ndim1,ndim2,cg,indxcg,jcg,cy,s_lat,s)
C- Integrals between smooth Hankels with k-th power of Laplace operator.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1's digit (not implemented: always vectors)
Ci         :0 rsm1,rsm2,e1,e2 are scalars
Ci         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
Ci         :10's digit
Ci         :1: do not calculate strux for any (l1,l2) pairs for
Ci         :   if rsm(l1) or rsm(l2) is zero.
Ci   p1    :first center
Ci   p2    :second center
Ci   q     :wave number for Bloch sum
Ci   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
Ci   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
Ci   e1    :energies of smooth Hankels at p1 (l-dependent)
Ci   e2    :energies of smooth Hankels at p2 (l-dependent)
Ci   nlm1  :L-cutoff for functions at p1
Ci   nlm2  :L-cutoff for functions at p2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimension of s
Ci   ndim2 :second dimension of s
Ci   cg    :Clebsch Gordan coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals; see Remarks
Cr Remarks
Cr   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cr   Row L corresponds to p1 and col M corresponds to p2.
Cr   Strux s(L,M) is computed for dr=p1-p2
Cr   See JMP 39, 3393, Section 10
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   8 Jun 00  Added 10s digit mode.  New arg list.
Cu   19 May 00 Made rsm1,e1,rsm2,e2 l-dependent.  Elements for which
Cu             rsm1 =0 or rsm2 = 0 are not computed.
Cu   18 May 00 Adapted from nfp hhi_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,jcg(*),indxcg(*),nlm1,nlm2,kmax,ndim1,ndim2
      double precision p1(3),p2(3),q(3),cg(*),cy(*),
     .  rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*)
      double complex s(ndim1,ndim2,0:kmax)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer m,lmx1,lmx2,ll,l1,l2,k,jlm,ilm,lm11,lm21,lm12,lm22,l1t,l2t
      double precision dr(3)

      if (nlm1 == 0 .or. nlm2 == 0) return
      call tcn('hhibl')

      do  m = 1, 3
        dr(m) = p1(m)-p2(m)
      enddo
      lmx1 = ll(nlm1)
      lmx2 = ll(nlm2)

      do  k = 0, kmax
        do  jlm = 1, nlm2
          do  ilm = 1, nlm1
            s(ilm,jlm,k) = dcmplx(0d0,0d0)
          enddo
        enddo
      enddo

      l1t = -1
      do  l1 = 0, lmx1
        if (l1 <= l1t) cycle
        call gtbsl2(l1,lmx1,e1,rsm1,l1t)
C       l1t = l1

        l2t = -1
        do  l2 = 0, lmx2
          if (l2 <= l2t) cycle
          call gtbsl2(l2,lmx2,e2,rsm2,l2t)
C         l2t = l2

          lm11 = l1**2+1
          lm12 = (l1t+1)**2
          lm21 = l2**2+1
          lm22 = (l2t+1)**2
          if (mode/10 == 1 .and. rsm1(l1)*rsm2(l2) == 0) cycle
          call phhibl(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12,
     .      lm21,lm22,kmax,ndim1,ndim2,cg,indxcg,jcg,cy,s_lat,s)
        enddo
      enddo

      call tcx('hhibl')

      end
      subroutine phhibl(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,kmax,
     .  ndim1,ndim2,cg,indxcg,jcg,cy,s_lat,s)
C- Integrals between smooth Hankels with k-th power of Laplace operator.
C ----------------------------------------------------------------------
Ci Inputs
Ci   dr    :p1-p2
Ci   q     :wave number for Bloch sum
Ci   rsm1  :smoothing radius of Hankels at p1
Ci   rsm2  :smoothing radius of Hankels at p2
Ci   e1    :energy of smooth Hankels at p1
Ci   e2    :energy of smooth Hankels at p2
Ci   nlm1  :L-cutoff for functions at p1
Ci   nlm2  :L-cutoff for functions at p2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimension of s
Ci   ndim2 :second dimension of s
Ci   cg    :Clebsch Gordan coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals; see Remarks
Cr Remarks
Cr   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 May 00 Adapted from nfp hhi_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer jcg(*),indxcg(*),mlm1,nlm1,mlm2,nlm2,kmax,ndim1,ndim2
      double precision dr(3),q(3),cg(*),cy(*),rsm1,
     .  rsm2,e1,e2
      double complex s(ndim1,ndim2,0:kmax)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer nlm0,ktop0,icg,icg1,icg2,ii,ilm,ilm1,ilm2,indx,ip,k,
     .  ktop,l1,l2,ll,lm,lmax1,lmax2,lmaxx,nlmx
      parameter( nlm0=100, ktop0=10 )
      double precision fpi,e,fac,fac1,fac2,gam1,gam2,gamx,rsmx
      double complex hkl1(0:ktop0,nlm0),hkl2(0:ktop0,nlm0),
     .  hsm(nlm0),hsmp(nlm0)

      fpi = 16d0*datan(1.d0)
      gam1 = 0.25d0*rsm1*rsm1
      gam2 = 0.25d0*rsm2*rsm2
      gamx = gam1+gam2
      rsmx = 2d0*dsqrt(gamx)
      lmax1 = ll(nlm1)
      lmax2 = ll(nlm2)
      lmaxx = lmax1+lmax2
      nlmx = (lmaxx+1)**2
      ktop = max0(lmax1,lmax2)+kmax
      if (nlmx > nlm0) call rxi('increase nlm0 in hhibl need',nlmx)
      if (ktop > ktop0) call rx('hhibl: increase ktop0')

C ... Set up functions for connecting vector p1-p2
      if (dabs(e1-e2) > 1d-5) then
        fac1 = dexp(gam2*(e2-e1))/(e1-e2)
        fac2 = dexp(gam1*(e1-e2))/(e2-e1)
        call hklbl(dr,rsmx,e1,q,ktop,nlmx,ktop0,cy,s_lat,hkl1)
        call hklbl(dr,rsmx,e2,q,ktop,nlmx,ktop0,cy,s_lat,hkl2)
        do  ilm = 1, nlmx
          do  k = 0, ktop
            hkl1(k,ilm) = fac1*hkl1(k,ilm) + fac2*hkl2(k,ilm)
          enddo
        enddo
      else
        e = .5d0*(e1+e2)
        call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop0,cy,s_lat,hkl2)
        call hsmbl(dr,rsmx,e,q,lmaxx,cy,s_lat,hsm,hsmp)
        do  ilm = 1, nlmx
          hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
          do  k = 1, ktop
            hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
          enddo
        enddo

      endif

C ... Combine with Clebsch-Gordan coefficients
      do  ilm1 = mlm1, nlm1
      l1 = ll(ilm1)
      do  ilm2 = mlm2, nlm2
        l2 = ll(ilm2)
        ii = max0(ilm1,ilm2)
        indx = (ii*(ii-1))/2 + min0(ilm1,ilm2)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1)-1
        do  icg = icg1, icg2
          ilm = jcg(icg)
          lm = ll(ilm)
          k = (l1+l2-lm)/2
          fac = fpi*(-1d0)**l1*cg(icg)
          do  ip = 0, kmax
            s(ilm1,ilm2,ip) = s(ilm1,ilm2,ip)+fac*hkl1(k+ip,ilm)
          enddo
        enddo
      enddo
      enddo

      end
