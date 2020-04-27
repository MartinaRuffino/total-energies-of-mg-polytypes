      subroutine ghgiml(job,dp,elh,rsmh,c01,rsmg,nlmh,nlmg,pmax,kmax,
     .  ktop0,ndimg,ptop0,cg,indxcg,jcg,cy,sr,si,sdr)
C- Integrals between energy-independent Gaussians and smoothed Hankels
C  centered at Rg and Rh, respectively:  \int G_kL(r-Rg)*H_pM(e,r-Rh) d^3 r
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :1s digit fixes which integral to return:
Ci         :0, [Gaussians G_kL] * [Hankels H_pM]
Ci         :   The G_kL do NOT contain energy factor; see JMP Eq. 12.2
Ci         :1, Return only the real part of the preceding integral
Ci         :2, [G_kL] * [H_pM - H_pM(unsmoothed)]
Ci         :3, not allowed
Ci         :4, [G_kL] * [H_0M + c01(0,m)*G_0M + c01(1,m)*G_1M],  p=0
Ci         :   [G_kL] * [H_pM],                                  p>0
Ci         :5, Return only the real part of the preceding integral
Ci         :6, Same as 4, but subtract int [G_kL] *H_pM(unsmoothed)
Ci         :10s digit
Ci         :> 0, same also expands energy derivative
Ci         :  NB: can use only with 1s digit job = 1 for now.
Ci   dp    :Connecting vector between the sites: dp = Rg - Rh (see Remarks)
Ci   elh,rsmh:vectors of l-dependent energies and smoothing radii of
Ci         :smoothed Hankels; must be specified for 0..ll(nlmh)
Ci   rsmg  :vector of l-dependent smoothing radii of the Gaussians associated
Ci         :with polynomials P_kL centred at the origin of the expansion
Ci         :rsmg must be specified for 0..ll(nlmg)
Ci   nlmh  :L-cutoff for the Hankels
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   pmax  :max power of the laplace operator for Hankels to be expanded
Ci   kmax  :polynomial cutoff
Ci   ktop0,ndimg,ptop0:leading dimensions of sr and si
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   sr,si :real and imaginary part of s, s = sr + i*si (see Remarks)
Co         :for negative energies si set to 0
Co         :if job 1s digit > 1, si is not referenced
Co   sdr   :real part of s-dot (integrals of Gaussians with H-dot rather
Co         : than with H)
Cl Local variables
Cl   yl    :unnormalised sph. harmonic polynomials
Cr Remarks
Cr   1. Integrals are defined as:
Cr      s(k,L,p,M) = \int G_kL(e=0,rsg;r-Rg) * H_pM(eh,rsh;r-Rh) d^3r
Cr              = (-)^l \sum_N C_LMN exp(-rsg^2*eh/4) H_qN(eh,rs;Rg-Rh),  (*)
Cr   where q = k + p + (l+m-n)/2; rs = sqrt(rsg^2 + rsh^2),
Cr   Rg = pg (centre of expansion) and
Cr   Rh = ph (centre of the sm. Hankel to be expanded)
Cr   See MvS's notes on smooth Hankels.
Cr
Cr   If Hankel energies eh are not l-dependent, the exponent can be moved
Cr   outside the summation sign in (*)
Cr
Cr   2. For the difference between smoothed and bare Hankels (job 1s digit = 2),
Cr      H_qN(eh,rs;Rg-Rh) in (*) are replaced with:
Cr      [H_qN(eh,rs;Rg-Rh) - H_qN(eh,rsg;Rg-Rh)]
Cr
Cr   3. For the energy derivatives of sm. Hankels (job 10s digit > 0), s-dot
Cr      is found as:
Cr      s-dot(k,L,p,M) = \int G_kL(e=0,rsg;r-Rg) * H-dot_pM(eh,rsh;r-Rh) d^3 r
Cr      = (-)^l \sum_N C_LMN exp(-rsg^2*eh/4)
Cr        [ H-dot_qN(eh,rs;Rg-Rh) - (-rsg^2/4)*H_qN(eh,rs;Rg-Rh) ]  (**)
Cr
Cb Bugs
Cb   only job = 11 is currently implemented for energy derivatives of
Cb   sm. Hankels. Hence, in particular, sdi is not yet included into
Cb   the argument list.
Cb
Cu Updates
Cu   04 Apr 15 Bug fix for job0 = 4,5,6
Cu   19 Sep 11 New option to return strux expanding (hsm+c01(0)*g0+c01(1)*g1)
Cu   15 May 08 energy derivatives of sm. Hankels added
Cu   16 Jul 07 positive energies, laplacians, and (smooth - bare) Hankel
Cu             options included (S. Lozovoi)
Cu   23 Jan 07 Adapted from ghibl.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlmg,nlmh,pmax,kmax,ktop0,ndimg,ptop0
      integer job,jcg(*),indxcg(*)
      double precision rsmg(0:*),rsmh(0:*),elh(0:*),c01(0:1,0:*)
      double precision dp(3),cg(*),cy(*)
      double precision sr(0:ktop0,ndimg,0:ptop0,nlmh)
      double precision si(0:ktop0,ndimg,0:ptop0,nlmh)
      double precision sdr(0:ktop0,ndimg,0:ptop0,nlmh)
C ... Local parameters
      integer job0
      integer nlm0,ktopm,lmx0
      integer icg,icg1,icg2,ii,ilg,ilh,ilm,indx,mh,mg,
     . k,ip,ik,ktop,lg,lh,ll,lmin,lm,lmaxg,lmaxh,lmaxx,ione
      parameter (lmx0=12, nlm0=(lmx0+1)**2, ktopm=31)
      logical lpos,ldot
      double precision cc,chg
      double precision rmxx,rmx1,rmx2,dp1,dp2,eh,ehx,beta
      double precision rmg,rmg2,rmg25,rmgx
C     double precision rsmv(-1:lmx0)
      double precision gkl(0:ktopm,0:lmx0),gkl0(0:ktopm,0:lmx0),yl(nlm0)
      double precision xi(-1:lmx0,0:ktopm),xi0(-1:lmx0,0:ktopm)
      double precision fi(-1:lmx0,0:ktopm),xid(0:lmx0,0:ktopm)
      double precision fpi,tol
      data tol/1.d-16/

      if (nlmh == 0 .or. nlmg == 0) return
      job0 = mod(job,10)
      ldot = (mod(job/10,10) > 0)
      if (ldot .and. (job0 /= 1)) call rx('ghdiml: H-dot option'//
     .  ' is implemented only for mod(job,10) = 1')
      call sanrg(.true.,job0,0,6,' ghgiml','job')

C ... rsm-independent setup
      if (ldot) then
        lmin = -1
      else
        lmin = 0
      endif
      fpi = 16d0*datan(1d0)
      lmaxh = ll(nlmh)
      lmaxg = ll(nlmg)
      lmaxx = lmaxg + lmaxh
      if (job0 >= 4) then
        ktop = max0(lmaxg,lmaxh) + max(pmax,2) + kmax
      else
        ktop = max0(lmaxg,lmaxh) + pmax + kmax
      endif

      if (lmaxx > lmx0) call rxi('ghdiml: increase lmx0 to',lmaxx)
      if (ktop > ktopm) call rxi('ghdiml: increase ktopm to',ktop)

C ... Make sph. harmonics for connecting vector dp
      call sylm(dp,yl,lmaxx,dp2)
      dp1 = dsqrt(dp2)
      do  ilm = 1, (lmaxx+1)**2
        yl(ilm) = yl(ilm)*cy(ilm)
      enddo

C ... Initialise s
      call dpzero(sr(0,1,0,1),(ktop0+1)*ndimg*(ptop0+1)*nlmh)
      if (job0 == 0) call dpzero(si(0,1,0,1),(ktop0+1)*ndimg*(ptop0+1)*nlmh)
      if (ldot) call dpzero(sdr(0,1,0,1),(ktop0+1)*ndimg*(ptop0+1)*nlmh)
      call dpzero(xi,(lmx0+2)*(ktop+1))
      call dpzero(xi0,(lmx0+2)*(ktop+1))

C --- Outer loops over lh and lg  ---
      rmgx = -1.d2
      rmxx = -1.d2
      ione = -1
      do  lg = 0, lmaxg
        ione = -ione
        rmg = rsmg(lg)
        rmg2 = rmg*rmg
        rmg25 = rmg2*0.25d0
        ehx = 1.d20
        do  lh = 0, lmaxh
          eh = elh(lh)
          lpos = (eh > 0d0)
          rmx2 = rmg2 + rsmh(lh)**2
C     ... if total rsm or energy changed, update smoothed Hankels
          if (dabs(eh-ehx)+dabs(rmx2-rmxx) > tol) then
            ehx = eh
            rmxx = rmx2
            beta = fpi*dexp(0.25d0*eh*rmx2)
            chg = dexp(-0.25d0*eh*rmg2)
            rmx1 = dsqrt(rmx2)
C           rsmv(0) = -rmx1
C       ... radial Hankels / r^l
C           call hansmz(dp1,eh,rsmv,xi(0,0),fi(0,0),lmaxx)
            call hansrz(-rmx1,lmin,lmaxx,eh,dp2,1,1,10,
     .                  xi(lmin,0),fi(lmin,0))
C       ... make H-dot as: xi-dot_l = 0.5*xi_{l-1} [JMP](7.5)
            if (ldot) then
              do  ii = 0, lmaxx
                xid(ii,0) = 0.5d0*xi(ii-1,0)
              enddo
            endif

C     ...   make laplacians H_k+1,L = -e*H_kL - 4*pi*G_kL
            if (ktop > 0) then
              call radgkl(dp1,rmx1,ktop-1,lmaxx,ktopm,gkl)
              do  ik = 1, ktop
                do  ii = 0, lmaxx
                  xi(ii,ik) = -eh*xi(ii,ik-1) - beta*gkl(ik-1,ii)
C     ...         if ldot, also make H-dot_pL:
C                 H-dot_p+1,L = - e*H-dot_pL - H_pL - 4*pi*(rs^2/4)*G_pL
                  if (ldot) xid(ii,ik) = -eh*xid(ii,ik-1) - xi(ii,ik-1)
     .                                   - 0.25d0*rmx2*beta*gkl(ik-1,ii)
                enddo
              enddo
              if (job0 == 0 .and. lpos) then
                do  ik = 1, ktop
                  do  ii = 0, lmaxx
                    fi(ii,ik) = -eh*fi(ii,ik-1)
                  enddo
                enddo
              endif
            endif
C     ...   if job0 = 2 or 6, update bare Hankels if necessary and
C           replace xi --> (xi-xi0)
            if (job0 == 2 .or. job0 == 6) then
              if (dabs(rmg-rmgx) > tol) then
                rmgx = rmg
C               rsmv(0) = -rmg
C               call hansmz(dp1,eh,rsmv,xi0,fi,lmaxx)
                call hansrz(-rmg,lmin,lmaxx,eh,dp2,1,1,10,xi0(lmin,0),fi(lmin,0))
                if (ktop > 0) then
                  call radgkl(dp1,rmg,ktop-1,lmaxx,ktopm,gkl0)
                  do  ik = 1, ktop
                    do  ii = 0, lmaxx
                      xi0(ii,ik) = -eh*xi0(ii,ik-1) - beta*gkl0(ik-1,ii)
                    enddo
                  enddo
                endif
              endif
              call daxpy((lmx0+2)*(ktop+1),-1d0,xi0,1,xi,1)
            endif
          endif

C     --- Combine with Clebsch-Gordan coefficients
          ilh = lh*lh
          do  mh = -lh, lh
            ilh = ilh + 1
            ilg = lg*lg
            do  mg = -lg, lg
              ilg = ilg + 1
              ii = max0(ilg,ilh)
              indx = (ii*(ii-1))/2 + min0(ilg,ilh)
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1)-1
              do  icg = icg1, icg2
                ilm = jcg(icg)
                lm = ll(ilm)
                k = (lg+lh-lm)/2
                cc = ione*chg*cg(icg)*yl(ilm)
                do  ip = 0, pmax
                  do  ik = 0, kmax
                    sr(ik,ilg,ip,ilh) = sr(ik,ilg,ip,ilh) + xi(lm,k+ik+ip)*cc
                    if (ldot)
     .              sdr(ik,ilg,ip,ilh) = sdr(ik,ilg,ip,ilh) +
     .                     (xid(lm,k+ik+ip)-rmg25*xi(lm,k+ik+ip))*cc
                  enddo
                enddo
C               Add gaussian parts
                if (job0 >= 4) then
                  cc = ione*cg(icg)*yl(ilm)
                  do  ip = 0, pmax
                  do  ik = 0, kmax
                    sr(ik,ilg,ip,ilh) = sr(ik,ilg,ip,ilh) +
     .                                  cc*(c01(0,lh)*gkl(k+ik+ip,lm) +
     .                                      c01(1,lh)*gkl(k+ik+ip+1,lm))
                  enddo
                  enddo
                endif
                if (job0 == 0 .and. lpos) then
                  do  ip = 0, pmax
                    do  ik = 0, kmax
                      si(ik,ilg,ip,ilh) = si(ik,ilg,ip,ilh) +
     .                                    fi(lm,k+ik+ip)*cc
                    enddo
                  enddo
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

      end

