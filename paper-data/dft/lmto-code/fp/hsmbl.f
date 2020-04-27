      subroutine hsmbl(p,rsm,e,q,lmax,cy,s_lat,hsm,hsmp)
C- Bloch-sum of smooth Hankel functions and energy derivatives, by Ewald
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat awald tol vol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qlv dlv
Cio    Passed to:  *
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsm   :smoothing radius
Ci   e     :energy of smoothed Hankel
Ci   q     :wave number for Bloch sum
Ci   lmax  :l-cutoff for hsm
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   hsm   :Bloch-summed smoothed Hankels
Co   hsmp  :Energy derivatives of hsm
Cr Remarks
Cr  Bloch sum is defined as (J. Math. Phys. 39,3393 (1998), Eq. 3.16)
Cr    Hsm^b(q,r) = sum_T exp(iq.T) Hsm(r-T)
Cr  Hankel functions are evaluated by Ewald summation.
Cr  p in units of alat, qlv in units of 2 pi / alat
Cr  See also hsmq for a vectorized version.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   26 Jan 07 Works with positive energies e
Cu   1 May 00 Adapted from nfp hsmbl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lmax
      double precision rsm,e,q(3),p(3),cy(*)
      double complex hsm(*),hsmp(*)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer nkd,nkq,ilm,l,m
      double precision alat,p1(3),plat(3,3),qlat(3,3),sp,pi,rwald,awald,
     .  asm,tol,vol
      double complex cfac,phase

C ... Shorten p, multiply by corresponding phase later
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      call shorbz(p,p1,plat,qlat)
C      print *, '!!'
C      p1 = p
      pi = 4d0*datan(1d0)
      sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
      phase = dcmplx(dcos(sp),dsin(sp))

      awald = s_lat%awald
      tol = s_lat%tol
      vol = s_lat%vol
      nkd = s_lat%nkd
      nkq = s_lat%nkq
      rwald = 1d0/awald
      asm = 1d0/rsm

      if (rsm < rwald) then
        call hsmblq(p1,e,q,awald,lmax,alat,s_lat%qlv,nkq,vol,hsm,hsmp)
        call hsmbld(p1,rsm,e,q,awald,lmax,alat,s_lat%dlv,nkd,hsm,hsmp)
      else
        call hsmblq(p1,e,q,asm,lmax,alat,s_lat%qlv,nkq,vol,hsm,hsmp)
      endif

C ... Multiply by phase to undo shortening
      cfac = dcmplx(0d0,1d0)*phase
      ilm = 0
      do  l = 0, lmax
        cfac = cfac*dcmplx(0d0,-1d0)
        do  m = 1, 2*l+1
          ilm = ilm+1
          hsm(ilm)  = cfac*cy(ilm)*hsm(ilm)
          hsmp(ilm) = cfac*cy(ilm)*hsmp(ilm)
        enddo
      enddo

C     call zprm('hsm',2,hsm,ilm,ilm,1)
      end

      subroutine hsmblq(p,e,q,a,lmax,alat,qlv,nkq,vol,dl,dlp)
C- Bloch sum, smooth Hankel for smoothing radius 1/a, done in recip space
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p (units of alat)
Ci   e     :energy of smoothed Hankel
Ci   q     :wave number for Bloch sum, in units of 2*pi/a
Ci   a     :Ewald smoothing parameter
Ci   lmax  :l-cutoff for dl,dlp
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   qlv   :reciprocal lattice vectors, units of 2pi/alat
Ci   nkq   :number of r.l.v.
Ci   vol   :cell volume
Co Outputs
Co   dl    :k-summed smoothed Bloch hankels hsm(1/a,e,r)
Co         :dl is not yet scaled by cy(ilm)
Co   dlp   :energy derivative of dl
Co         :dlp is not yet scaled by cy(ilm)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkq,lmax
      double precision a,alat,e,vol
      double precision q(3),p(3),qlv(3,nkq)
      double complex dl(1),dlp(1)
C ... Local parameters
      integer lmxx,nlm,ilm,ir
      parameter (lmxx=11)
      double precision r(3),tpi,gamma,fpibv,tpiba,scalp,r2,den0,den1
      double precision yl((lmxx+1)**2)
      double complex eiphi

      if (lmax > lmxx) call rxi('hsmblq: increase lmxx to',lmax)
      tpi = 8d0*datan(1d0)
      gamma = .25d0/(a*a)
      fpibv = 2d0*tpi/vol
      tpiba = tpi/alat
      nlm = (lmax+1)**2
      do  ilm = 1, nlm
        dl(ilm) = dcmplx(0d0,0d0)
        dlp(ilm) = dcmplx(0d0,0d0)
      enddo
      do  ir = 1, nkq
        r(1) = tpiba*(q(1)+qlv(1,ir))
        r(2) = tpiba*(q(2)+qlv(2,ir))
        r(3) = tpiba*(q(3)+qlv(3,ir))
        scalp = alat*(r(1)*p(1)+r(2)*p(2)+r(3)*p(3))
        eiphi = dcmplx(dcos(scalp),dsin(scalp))
        call sylm(r,yl,lmax,r2)
        den0 = dexp(gamma*(e-r2))/(r2-e)
        den1 = den0/(r2-e)
        do  ilm = 1, nlm
          dl(ilm)  = dl(ilm) + yl(ilm)*eiphi*den0
          dlp(ilm) = dlp(ilm)+ yl(ilm)*eiphi*den1
        enddo
      enddo
      do  ilm = 1, nlm
        dl(ilm)  = fpibv*dl(ilm)
        dlp(ilm) = fpibv*dlp(ilm) + gamma*dl(ilm)
      enddo
      end

      subroutine hsmbld(p,rsm,e,q,a,lmax,alat,dlv,nkd,dl,dlp)
C- Adds real space part of reduced structure constants (ewald).
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :connecting vector, in units of a
Ci   rsm   :sm Hankel smoothing radius, in a.u.
Ci   e     :sm Hankel energy, in a.u.
Ci   q     :Bloch wave number, in units of 2*pi/a
Ci   a     :Ewald parameter
Ci   lmax  :l-cutoff for dl,dlp
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   dlv   :direct  lattice vectors, units of alat
Ci   nkd   :number of d.l.v
Co Outputs
Co   dl    :Difference in Bloch summed hsm(rsm,e,r) and hsm(1/a,e,r)
Co   dlp   :Energy derivative of dl
Cl Local variables
Cl         :
Cr Remarks
Cr  Let asm=1/rsm, r = p-T, r1=|r|, and Yl = Ylm(r).
Cr  Routine makes difference betw/ hsm(rsm,e,r) and hsm(1/a,e,r)
Cr    dl(ilm)  = sum_T yl(ilm)*(chi2(l)-chi1(l))*exp(iq.T)
Cr    dlp(ilm) = sum_T yl(ilm)*(chi2(l-1)-chi1(l-1))*exp(iq.T)/2
Cr    where chi1 = chi(a,kap,r), chi2 = chi(asm,kap,r)
Cr  chi is obtained recursively from chi0(0), chi(-1).
Cr  Formulas for chi(0) and chi(-1) explicitly: where a'= a or asm
Cr  Case e>0, kap=dsqrt(e)
Cr    expikr = exp(i*kap*r1)
Cr    zuplus = zerfc(i*kap/2/a'+r1*a')*expikr
Cr    chi(0) = expikr/r1 - dble(zuplus)/r1
Cr    chi(-1) = expikr/i*kap + dimag(zuplus)/kap
Cr  Case e<0, akap=dsqrt(-e)
Cr    uplus = derfc(.5d0*akap/a'+r1*a')/emkr
Cr    umins = derfc(.5d0*akap/a'-r1*a')*emkr
Cr    chi(0) = 0.5d0*(umins-uplus)/r1
Cr    chi(-1) = (umins+uplus)/(2*akap)
Cu Updates
Cu   10 May 07 New protections to handle error functions underflow
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmxx,ilm,ir,l,lmax,m,nkd,nm
      double precision q(3),p(3),dlv(3,nkd)
      double complex dl(1),dlp(1)
C ... Local parameters
      logical lpos,lzero
      parameter (lmxx=11)
      double precision a,a2,akap,alat,asm,asm2,cc,ccsm,derfc,e,emkr,gl,
     .  qdotr,r1,r2,rsm,srpi,ta,ta2,tasm,tasm2,umins,uplus,tpi,kap
      double precision yl((lmxx+1)**2),chi1(-1:10),chi2(-1:10),r(3)
C     double complex zchi1(-1:10),zchi2(-1:10)
      double complex cfac,zikap,expikr,zuplus,zerfc
      double precision srfmax,fmax
      parameter (srfmax=16d0, fmax=srfmax*srfmax)

      if (lmax > lmxx) call rxi('hsmbld: increase lmxx to',lmax)
      tpi = 8d0*datan(1d0)
      srpi = dsqrt(tpi/2d0)
      lpos = e > 0d0
      if (lpos) then
C     Methfessl's 'akap' is sqrt(-e) = i kap by standard kap=sqrt(e)
        kap = dsqrt(e)
        zikap = dcmplx(0d0,1d0)*kap
      else
        akap = dsqrt(-e)
      endif
      ta = 2d0*a
      a2 = a*a
      ta2 = 2d0*a2
      cc = 4d0*a2*a*dexp(e/(ta*ta))/srpi

      asm = 1d0/rsm
      tasm = 2d0*asm
      asm2 = asm*asm
      tasm2 = 2d0*asm2
      ccsm = 4d0*asm2*asm*dexp(e/(tasm*tasm))/srpi

      do  ir = 1, nkd
      r(1) = alat*(p(1)-dlv(1,ir))
      r(2) = alat*(p(2)-dlv(2,ir))
      r(3) = alat*(p(3)-dlv(3,ir))
      call sylm(r,yl,lmax,r2)
      r1 = dsqrt(r2)
      lzero = .false.

C --- Make the xi's from -1 to lmax ---
      if (r1 < 1d-6) then
        do  l = 1, lmax
          chi1(l) = 0d0
          chi2(l) = 0d0
        enddo
        if (lpos) then
C         do  l = 1, lmax
C           zchi1(l) = 0d0
C           zchi2(l) = 0d0
C         enddo
          zuplus = zerfc(zikap/ta)
          chi1(0) = -zikap*zuplus
     .      +2/dsqrt(4*datan(1d0))*a * exp(0.25_8*e/a2)
          zuplus = zerfc(zikap/tasm)
          chi2(0) = -zikap*zuplus
     .      +2/dsqrt(4*datan(1d0))*asm * exp(0.25_8*e/asm2)
          chi1(-1) = -zerfc(-zikap/ta)/zikap
          chi2(-1) = -zerfc(-zikap/tasm)/zikap
        else
          chi1(0) = ta*dexp(e/(2d0*ta2))/srpi - akap*derfc(akap/ta)
          chi1(-1) = derfc(akap/ta)/akap
          chi2(0) =tasm*dexp(e/(2d0*tasm2))/srpi - akap*derfc(akap/tasm)
          chi2(-1) = derfc(akap/tasm)/akap
        endif
      else
        if (lpos) then
          expikr = exp(zikap*r1)
          zuplus = zerfc(zikap/ta+r1*a)*expikr
          chi1(0) = expikr/r1 - dble(zuplus)/r1
          chi1(-1) = expikr/zikap + dimag(zuplus)/kap
C          zchi1(0) = expikr/r1 - dble(zuplus)/r1
C          zchi1(-1) = expikr/zikap + dimag(zuplus)/kap
        else
C         If large, these are -log(uplus),-log(umins); then chi->0
Cr        If both (akap*rsm/2 -/+ r/rsm) >> 1, we have
Cr        -log u(+/-) -> (akap*rsm/2 -/+ r/rsm)^2 -/+ akap*r
Cr                    =  (akap*rsm/2)^2 + (r/rsm)^2 >> 1
Cr         u(+/-)     -> exp[-(akap*rsm/2)^2 - (r/rsm)^2] -> 0
Cr        Also if akap*r >> 1,   chi < dexp(-akap*r1) -> 0
          emkr = dexp(-akap*r1)
          if (.5d0*akap/a+r1*a > srfmax .and.
     .        .5d0*akap/a-r1*a > srfmax .or. akap*r1 > fmax) then
            lzero = .true.
C            uplus = derfc(.5d0*akap/a+r1*a)/emkr
C            umins = derfc(.5d0*akap/a-r1*a)*emkr
C            print *, 'approx',(.5d0*akap/a+r1*a)**2 - akap*r1
C            print *, '-log up',-dlog(uplus)
C            print *, 'approx', (.5d0*akap/a-r1*a)**2 + akap*r1
C            print *, '-log um', -dlog(umins)
C            print *, r1*a
C            stop
          else
            uplus = derfc(.5d0*akap/a+r1*a)/emkr
            umins = derfc(.5d0*akap/a-r1*a)*emkr
            chi1(0) = 0.5d0*(umins-uplus)/r1
            chi1(-1) = (umins+uplus)/(2d0*akap)
          endif
        endif
        if (lzero) then
          do  l = -1, lmax
            chi1(l) = 0
          enddo
          lzero = .false.
        else
          gl = cc*dexp(-a2*r2)/ta2
          do  l = 1, lmax
            chi1(l) = ((2*l-1)*chi1(l-1)-e*chi1(l-2)-gl)/r2
            gl = ta2*gl
          enddo
        endif

C        chi2 is complex; so is chi1, but the imaginary part
C        is the same, so the difference is real
        if (lpos) then
          zuplus = zerfc(zikap/tasm+r1*asm)*expikr
          chi2(0) = expikr/r1 - dble(zuplus)/r1
          chi2(-1) = expikr/zikap + dimag(zuplus)/kap
C          zchi2(0) = expikr/r1 - dble(zuplus)/r1
C          zchi2(-1) = expikr/zikap + dimag(zuplus)/kap
        elseif (.5d0*akap/asm+r1*asm > srfmax .and.
     .        .5d0*akap/asm-r1*asm > srfmax .or. akap*r1 > fmax)then
          lzero = .true.
        else
          uplus = derfc(.5d0*akap/asm+r1*asm)/emkr
          umins = derfc(.5d0*akap/asm-r1*asm)*emkr
          chi2(0) = 0.5d0*(umins-uplus)/r1
          chi2(-1) = (umins+uplus)/(2d0*akap)
        endif
        if (lzero) then
          do  l = -1, lmax
            chi2(l) = 0
          enddo
          lzero = .false.
        else
          gl = ccsm*dexp(-asm2*r2)/tasm2
          do  l = 1, lmax
            chi2(l) = ((2*l-1)*chi2(l-1)-e*chi2(l-2)-gl)/r2
            gl = tasm2*gl
          enddo
        endif
      endif
      qdotr = tpi*(q(1)*dlv(1,ir)+q(2)*dlv(2,ir)+q(3)*dlv(3,ir))
      cfac = cdexp(dcmplx(0d0,qdotr))
      ilm = 0
      do  l = 0, lmax
        nm = 2*l+1
        do  m = 1, nm
          ilm = ilm+1
          dl(ilm) = dl(ilm) + yl(ilm)*(chi2(l)-chi1(l))*cfac
          dlp(ilm) = dlp(ilm)+yl(ilm)*0.5d0*(chi2(l-1)-chi1(l-1))*cfac
        enddo
        cfac = cfac*dcmplx(0d0,1d0)
      enddo
      enddo

      end
