      subroutine hdxpml(job,ph,pg,eh,rsmh,c01,rsmg,pmax,kmin,kmax,nlmh,
     .  nlmg,k0,ndim,p0,cg,indxcg,jcg,cy,sr,si,sdr)
C- Strux s expanding sm hankels H_pL(r-ph) and related functions
C  into sum(kM) s(k,M,p,L) * P_kM(r-pg) and possibly Bessel functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :1s digit fixes which structure constants s to return:
Ci         :0, expansion of H_pL
Ci         :1, real part of H_pL only
Ci         :2, expansion of [H_pL - H_pL(rsm=0)] (always real)
Ci         :3, not allowed
Ci         :4, expansion of hsm + c01(0)*g0l + c01(1)*g1l
Ci         :5, same as 4, real part only
Ci         :6, expansion of hsm + c01(0)*g0l + c01(1)*g1l - hsm(rsm=0))
Ci             If ALSO kmin<0,
Ci             s(-1,:,0,:) holds coff to Bessel expansion of hsm(rsm=0)
Ci         :if 10s digit > 0, same as above but in addition expands energy
Ci                            derivative of sm. Hankels
Ci   ph    :Function is centered at ph; see Remarks
Ci   pg    :origin of the expansion; see Remarks
Ci   eh    :vector of l-dependent energies of smoothed Hankel
Ci         :eh must be specified for 0..ll(nlmh)
Ci   c01   :Gaussian coefficients for generalized envelope functions
Ci   rsmh  :vector of l-dependent smoothing radii of smoothed Hankel
Ci         :EITHER must be specified for 0..ll(nlmh)
Ci         :OR rsmh(0) = const <0. Implies rsmh(l) = -const for all l
Ci   rsmg  :vector of l-dependent smoothing radii for polynomials
Ci         :same convention as for rsmh
Ci   pmax  :max power of the Laplace operator (lap)^p H_L := H_pL
Ci         :NB pmax should be 0 or 1 if 1st digit job is 4,5,6.
Ci   kmax  :P_kL polynomial k-cutoff
Ci   nlmh  :L-cutoff for smoothed Hankel functions being expanded
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   kmin,k0:upper and lower bounds dimensioning coff array s; see job=6
Ci   ndim,p0:dimensions coefficient arrays sr and si
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   sr,si :real and imaginary parts of cofficients s(k,M,p,L); see Remarks
Co         :here k=0..kmax, M=1..nlmg, p=0..pmax, L=1..nlmh
Co         :if job 1s digit = 0, si = 0 for negative energies
Co         :if job 1s digit > 0, si is not referenced
Co   sdr   :real part of s-dot(k,M,p,L), ie the coefficients of
Co         :H-dot expansion into P_kL
Cr Remarks
Cr   Expansion is:  H_pL(r-ph) = sum(kM) s(k,M,p,L) * P_kM(r-pg)
Cr   See J. Math. Phys. {\bf 39}, 3393 (1998), Section XII.
Cr   For p=0 case, see JMP, Eq. 12.17.
Cr   Recurrence relation for p>0 case: see JMP Eqs. 6.33 and 5.11
Cr
Cr   Same form for H-dot:
Cr          H-dot_pL(r-ph) = sum(kM) s-dot(k,M,p,L) * P_kM(r-pg)
Cr
Cr   As rsmg -> 0, the expansion turns into a Taylor series of H_L.
Cr   As rsmg increases the error is spread out over a progressively
Cr   larger radius, thus reducing the error for larger r while
Cr   reducing accuracy at small r.
Cb Bugs
Cb   Does not handle the case pg = ph
Cb   Does not handle zero smoothing case rsmh(l) or rsmg(l) = 0
Cb   H-dot option is implemented only for job = 11
Cu Updates
Cu   04 Apr 15 Fixed bug in call to ghgiml
Cu   19 Aug 11 New option to return strux expanding (hsm + c01*g0)
Cu             Bug fix for job0=2 or 6
Cu   15 May 08 (S. Lozovoi) H-dot added
Cu   17 Jul 07 (S. Lozovoi) Updated to accomodate changes in ghiml.f
Cu   23 Jan 07 (S. Lozovoi) Adapted from hxpbl.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,p0,k0,pmax,kmin,kmax,ndim,nlmg,nlmh,jcg(*),indxcg(*)
      double precision eh(0:*),rsmh(0:*),rsmg(0:*),c01(0:1,0:*)
      double precision ph(3),pg(3),cg(*),cy(*)
      double precision sr(kmin:k0,ndim,0:p0,nlmh)
      double precision si(kmin:k0,ndim,0:p0,nlmh)
      double precision sdr(kmin:k0,ndim,0:p0,nlmh)
C ... Dynamically allocated arrays
      real(8),allocatable:: slr(:,:,:,:),sli(:,:,:,:),sd(:,:,:,:)
C     real(8),allocatable:: sg(:,:,:,:)
C ... Local parameters
      logical ldot
      integer n0
      parameter (n0=10)
      integer job0,jobl,ilmg,ilmh,k,p,il,ll,lmaxg,lmaxh,m,nm,i
      double precision dfact,fac,factk,fpi,rs,rs2,facrs
      double precision rsmhl(0:n0),rsmgl(0:n0),dp(3)
C     double precision cglm(0:1,nlmh)

      if (nlmg == 0 .or. nlmh == 0) return
      job0 = mod(job,10)
      ldot = (mod(job/10,10) > 0)
      if (ldot .and. (job0 /= 1)) call rx('hdxpml: H-dot option'//
     .  ' is currently implemented only with mod(job,10) = 1')
      call sanrg(.true.,job0,0,6,' hdxpml','job')
      fpi = 16d0*datan(1d0)

C ... Local copies of rsmg and rsmh; const rsm flagged by rsm(0)<0
      lmaxh = ll(nlmh)
      if (rsmh(0) < 0d0) then
        call dvset(rsmhl(0),1,lmaxh+1,-rsmh(0))
      else
        call dcopy(lmaxh+1,rsmh(0),1,rsmhl(0),1)
      endif
      lmaxg = ll(nlmg)
      if (rsmg(0) < 0d0) then
        call dvset(rsmgl(0),1,lmaxg+1,-rsmg(0))
      else
        call dcopy(lmaxg+1,rsmg(0),1,rsmgl(0),1)
      endif

C ... Only the connecting vector matters
      do i = 1, 3
        dp(i) = pg(i) - ph(i)
      enddo

C --- Integrals of gaussians * sm. hankels ---
      jobl = job
      if (job0 >= 4) jobl = job-4
C     Until in ghgiml for job0 = 4,5,6 has been fixed,
C     comment next line and uncomment the jobl block below
      jobl = job
      allocate(slr(0:k0,ndim,0:p0,nlmh),sli(0:k0,ndim,0:p0,nlmh))
      if (ldot) allocate(sd(0:k0,ndim,0:p0,nlmh))
      if (.not. ldot) allocate(sd(1,1,1,1))
      call ghgiml(jobl,dp,eh,rsmhl,c01,rsmgl,nlmh,nlmg,pmax,kmax,
     .  k0,ndim,p0,cg,indxcg,jcg,cy,slr,sli,sd)
C      if (jobl < job) then
C        p2 = 2
C        allocate(sg(0:k0,ndim,0:p2,nlmh))
C        call ggiml(dp,rsmhl,rsmgl,nlmh,nlmg,2,kmax,k0,ndim,p2,cg,indxcg,jcg,cy,sg)
C        ilmh = 0
C        do  lh = 0, lmaxh
C          do  mh = -lh, lh
C            ilmh = ilmh+1
C            cglm(0,ilmh) = c01(0,lh)
C            cglm(1,ilmh) = c01(1,lh)
C            do  p = 0, pmax
C            do  k = 0, kmax
C            do  ilmg = 1, (lmaxg+1)**2
C              slr(k,ilmg,p,ilmh) = slr(k,ilmg,p,ilmh) +
C     .          cglm(0,ilmh)*sg(k,ilmg,p,ilmh) + cglm(1,ilmh)*sg(k,ilmg,p+1,ilmh)
C            enddo
C            enddo
C            enddo
C          enddo
C        enddo
C        deallocate(sg)
C      endif

C ... Scale integrals to get coefficients for the P_kL expansion
      ilmg = 0
      dfact = 1d0
      do  il = 0, lmaxg
        nm = 2*il+1
        dfact = dfact*nm
        rs = rsmgl(il)
        rs2 = rs*rs
        do  m = 1, nm
          ilmg = ilmg+1
          factk = 1d0
          facrs = rs**il
          do  k = 0, kmax
            fac = fpi*facrs / (factk*dfact)
            do  p = 0, pmax
            do  ilmh = 1, nlmh
              sr(k,ilmg,p,ilmh) = slr(k,ilmg,p,ilmh)*fac
c             if (ldot) sdr(k,ilmg,p,ilmh) = sd(k,ilmg,p,ilmh)*fac
            enddo
            enddo
c       ... make H-dot if required
            if (ldot) then
              do  p = 0, pmax
              do  ilmh = 1, nlmh
                sdr(k,ilmg,p,ilmh) = sd(k,ilmg,p,ilmh)*fac
              enddo
              enddo
            endif
C       ... Imaginary part is required only if job=0 and e positive
            if (job0 == 0 .or. job0 == 4) then
              do ilmh = 1, nlmh
                if (eh(ll(ilmh)) > 0) then
                  do p = 0, pmax
                    si(k,ilmg,p,ilmh) = sli(k,ilmg,p,ilmh)*fac
                  enddo
                endif
              enddo
            endif
            factk = 4*(k+1)*factk
            facrs = rs2*facrs
          enddo
        enddo
      enddo

      deallocate(slr,sli)
      if (ldot) deallocate(sd)

C ... Strux for expansion around connecting vector ph-pg
      if (job0 == 6 .and. kmin < 0) then
        allocate(slr(ndim,nlmh,1,1))
        call mstrx2(eh,ph-pg,nlmg,nlmh,ndim,cg,indxcg,jcg,cy,1,slr,slr)
        sr(-1,1:nlmg,0,1:nlmh) = slr(1:nlmg,1:nlmh,1,1)
        deallocate(slr)
      endif

      end
