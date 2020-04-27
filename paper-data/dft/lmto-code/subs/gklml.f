      subroutine gklml(p,rsml,e,kmax,nlm,k0,gkl)
C- k,L-dependent gaussians
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsml  : vector of l-dependent smoothing radii of Hankels
Ci         : EITHER must be specified for lmin..lmax
Ci         : OR     rsml(0) < 0. Implies rsml(l) = -rsml(0) for all l
Ci   e     :used to scale gkL by exp(e*rsml**2/4)
Ci   kmax  :polynomial cutoff
Ci   nlm   :L-cutoff for gkl
Ci   k0    :leading dimension of gkl
Co Outputs
Co   gkl   :generalized Gaussians
Cr Remarks
Cr   G_kL = laplace^k G_L; see Eq. 5.9, J Math Phys 39, 3393
Cu Updates
Cu   18 May 11 smoothing radius made l-dependent
Cu   24 Jan 07 Adapted from gklbl.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer k0,kmax,nlm
      double precision e,rsml(0:*),p(3)
      double precision gkl(0:k0,nlm)
C ... Local parameters
      logical lrsm0

      integer ilm,k,il,l,ll,m,lmax
      double precision r1,r2,cfac,rsm,rsx
      double precision wk(0:kmax,0:nlm-1),yl(nlm)
      double precision tol
      parameter (tol=1d-15)

      if (nlm <= 0) return

C      Not necessary
C      do  ilm = 1, nlm
C        do  k = 0, kmax
C          gkl(k,ilm) = 0d0
C        enddo
C      enddo

C     Negative 1st smoothing radius => constant rsml
      lrsm0 = rsml(0) <= 0

C     Spherical harmonic polynomials
      lmax = ll(nlm)
      call ropyln(1,p(1),p(2),p(3),lmax,1,yl,r2)
      r1 = dsqrt(r2)

C     Radial parts gkl / r^l
      rsm = dabs(rsml(0))
      rsx = -1d3
      do  l = lmax, 0, -1
        if (.not. lrsm0) rsm = dabs(rsml(l))
        if (dabs(rsm-rsx) > tol) then
          call radgkl(r1,rsm,kmax,l,kmax,wk)
C         Energy scaling
          if (e /= 0) then
            cfac = dexp(0.25d0*e*rsm*rsm)
            do  il = 0, l
              do  k = 0, kmax
                wk(k,il) = wk(k,il)*cfac
              enddo
            enddo
          endif
          rsx = rsm
        endif
      enddo

C     Contract radial with angular parts
      ilm = 0
      do  l = 0, lmax
        do  m = -l, l
          ilm = ilm+1
          do  k = 0, kmax
C           gkl(k,ilm) = gkl(k,ilm) + yl(ilm)*wk(k,l)
            gkl(k,ilm) = yl(ilm)*wk(k,l)
          enddo
        enddo
      enddo

      end

C     test rsm-dependent call
C      subroutine fmain
C      integer lmax,kmax
C      parameter (lmax=3,kmax=2)
C      double precision p(3),rsml(0:lmax),rsm1,eh
C      double precision gkl(0:kmax+1,(lmax+1)**2)
C      double precision gkl0(0:kmax+1,(lmax+1)**2)
C      integer l,ilm
C
C      data p /-2.4854171306613027d0,
C     .        -2.2983342613226059d0,
C     .        -2.1112513919839087d0/
C
C      rsm1 = 2.3d0 * 2.0d0/3
C      do  l = 0, lmax
C        rsml(l) = rsm1*(1-1d0*l/20d0)
C      enddo
C      eh = -0.5d0
C
C      do  l = lmax, 0, -1
C        call gklml(p,-rsml(l),eh,kmax,(l+1)**2,kmax+1,gkl0)
C      enddo
C      call gklml(p,rsml,eh,kmax,(lmax+1)**2,kmax+1,gkl)
C
C      do  ilm = 1, (lmax+1)**2
C        print 333, ilm,
C     .    gkl(0,ilm), gkl(0,ilm)-gkl0(0,ilm),
C     .    gkl(1,ilm), gkl(1,ilm)-gkl0(1,ilm),
C     .    gkl(2,ilm), gkl(2,ilm)-gkl0(2,ilm)
C  333   format(i4,6f15.10)
C      enddo
C
C      end
