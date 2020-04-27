      subroutine fpopm(s_ctrl,s_spec,s_lat,s_optic,nlmax,nrsmo,ndham,nphimx,
     .  isp,nsp,nspc,nbas,nev,ausp,cPkLq,opt,dsdk)
C- Contribution to electric dipole optical matrix elements from spheres
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham loptc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ips
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa kmxt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  fpopm_sfac
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  optme nlg ocrng unrng
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rgrad rgrade optme ocrng unrnge ocrng unrng
Cio    Passed to:  opt_nstate fpopm_sfac
Ci Inputs
Ci   nlmax :dimensions ausp,cPkLq
Cr   nrsmo :dimensions cPkLq ... max number of radial functions to expand envelopes
Ci   ndham :dimensions ausp,cPkLq
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nbas  :size of basis
Ci   nev   :number of eigenvectors
Ci   ausp  :coefficients to (phi,phidot) projection of w.f (makusq)
Ci   cPkLq :coefficients to (PkL) projection of w.f (makusq)
Ci   opt   :1s digit
Ci         :passed as 1s digit to gradme (>0 if cPkLq are in spherical harmonics)
Ci         :10s digit
Ci         :0 do not add nonlocal sigma
Ci         :1 add nonlocal sigma
Ci   dsdk  :k-gradient of nonlocal sigma
Co Outputs
Co   optme(i,j,m) (an element in s_optic).  Adds the local part of
Co         :<i|grad|j> connecting occ i with unocc j, polarization m
Co         :The total <i|grad|j> has three terms:
Co         : + matrix element of the true w.f. (calculated here)
Co         : - matrix element of the local rep of the envelope w.f. (calculated here)
Co         : + matrix element of the envelope w.f. (calculated in rsibl)
Co         : If ldsdk=T,  dsdk is also added.
Cr Remarks
Cr   Adapted from asaopm.f
Cu Updates
Cu   14 Jun 18 Reworked with updated gradme and real harmonics
Cu   21 Aug 14 Reworked to reflect the proper three-fold representation
Cu             of optical matrix elements
Cu   22 Jun 14 Extended to the noncollinear case
Cu   11 Feb 14 Some redesign for general case
Cu   15 Nov 13 (Ben Kaube) adapted from the ASA asaopm for full potential
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer opt,nsp,ndham,nphimx,nlmax,nrsmo,isp,nbas,nspc,nev
      double complex   ausp(nlmax,ndham*nspc,nphimx,nspc,nbas)
      double complex   cPkLq(nlmax,ndham*nspc,nrsmo,nspc,nbas)
      double complex   dsdk(ndham,ndham,3,nsp)
C ... For structures
!       include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_optic):: s_optic
C ... Dynamically allocated local arrays
      integer, allocatable :: lmxa(:),kmxa(:),ipc(:),nrfn(:),nlmi(:)
      complex(8), allocatable :: ccp(:,:,:,:,:)
C     real(8), allocatable :: gradm(:,:,:,:,:,:,:)
      complex(8), allocatable :: zgradm(:,:,:,:,:,:,:)
      complex(8), pointer :: optme(:,:,:)
C ... Local parameters
C     logical :: debug=.true.
C     integer,parameter:: nkap0=4
      logical ldsdk    !T if dsdk is to be added to optme
      integer job,i,ispc,iv(6),ksp,lham,loptic,morder,nrg,nlg,nsgrp,
     .  nfilo,nfiup,nemlo,nemup,nfilm,nempm,ivb,icb,ipol,iq,nkaph
      logical leps
      real(8) :: alat,frms,frmsi,maxfac,maxfaci,xx
      integer ncf,ncfi
      complex(8) :: cf,cfmean,cfmeani
      real(8), parameter :: pi = 4*datan(1d0)
      complex(8), parameter :: im = (0d0,1d0)
      procedure(integer) :: isw
      data cfmean /(0d0,0d0)/, frms /0d0/, maxfac /0d0/, ncf /0/
      save cfmean,frms,ncf,maxfac
C ... External calls
      external gradme,optdme,pvfpopm,rx,tcn,tcx
      procedure(integer) :: nglob

C     print *, '!! Interstitial part only'; return

      call opt_nstate(s_optic,3,nev,nfilo,nfiup,nemlo,nemup,nfilm,nempm)
      if (nfilm*nempm == 0) return

      call tcn('fpopm')

      alat = s_lat%alat
      optme => s_optic%optme
      if (size(optme,1)/=nfilm .or. size(optme,2)/=nempm .or.
     .    size(optme,3)/=3) call rx('fpopm: array mismatch')

      allocate(lmxa(nbas),kmxa(nbas),ipc(nbas),nrfn(nbas),nlmi(nbas))
      call spec2class(s_spec,nbas,s_ctrl%ips,'lmxa',1,lmxa,iv)
      call spec2class(s_spec,nbas,s_ctrl%ips,'kmxt',1,kmxa,iv)
      do  i  = 1, nbas
        ipc(i) = i
        nlmi(i) = (1+lmxa(i))**2
        nrfn(i) = nphimx ! for now: should be site-dependent
      enddo

      nlg = s_optic%nlg
      iv(1:5) = shape(s_optic%rgrad)
      i = nint(sqrt(dble(iv(1))))
      if (i /= nphimx) call rx('locpot: s_optic%rgrad has odd dimension')
C     if (nrg*nrg /= iv(1)) call rx('locpot: s_optic%rgrad has odd dimension')

      lham = s_ctrl%lham
      loptic = s_ctrl%loptc
      nsgrp = s_lat%nsgrp
      morder = mod(opt,10)  ! Passed as 1s digit to gradme
      leps  = loptic > 0 .and. mod(loptic,10) > 0
      if (.not. leps) goto 999

C --- Augmented wave functions ---
C ... Copy ausp to ccp; return ccp in spherical harmonics
      allocate(ccp(nlg**2,nbas,nspc,nphimx,nev)) !; call dpzero(ccp,2*size(ccp))
      call pvfpopm(nphimx,nbas,nlg,nlmax,ndham,lmxa,nspc,nev,ausp,ccp)

C ... Matrix element phi grad phi etc from radial matrix elements
      allocate(zgradm(nlg**2,nlg**2,3,nphimx,nphimx,nsp*nspc,nbas))
      call dpzero(zgradm,2*size(zgradm))
C     Ylm always in real harmonics
      do  ispc = 1, nspc
        ksp = max(ispc,isp)
        i = 30 + morder
C       i = i+ 300  ! Return gradient of x+iy, x-iy
C       print *, '!! debugging'; i = 331
        call gradme(i,nphimx,nphimx,ksp,nsp,nspc,1,nbas,nlg-1,lmxa,nlmax,
     .    s_optic%rgrad,zgradm,zgradm)
      enddo

C ... Electric dipole matrix elements from augmented w.f.
C     Uncomment for debugging
C     allocate(optme(nfilm,nempm,3)); call dpzero(optme,2*nfilm*nempm*3)
      call optdme(1,nphimx,nrfn,nlg,nbas,isp,nsp,nspc,nev,ipc,
     .  zgradm,nfilo,nfiup,nemlo,nemup,-1d0,ccp,optme)

      if (mod(opt/10,10) > 0) then
        ispc = min(isp,nspc)
        do  ipol = 1, 3
          do  ivb = nfilo, nfiup
            do  icb = nemlo, nemup
              optme(ivb,icb,ipol) = optme(ivb,icb,ipol)
     .          + (alat*im/(2d0*pi))*dsdk(ivb,icb,ipol,ispc)
            enddo
          enddo
        enddo
      endif

C     if (debug) call zprm('optme (1)',2,optme,nfilm,nfilm,nempm)
      deallocate(ccp,zgradm)

C --- Contribution from local projection of smooth envelope functions ---
      nkaph = nglob('nkaph')
      do  i  = 1, nbas
        nrfn(i) = kmxa(i)+1+nkaph ! number of envelope function types joined to (u,s,phiz)
      enddo

C ... Copy cPkLq to ccp; return ccp in spherical harmonics
      allocate(ccp(nlg**2,nbas,nspc,nrsmo,nev))
      call pvfpopm(nrsmo,nbas,nlg,nlmax,ndham,lmxa,nspc,nev,cPkLq,ccp)

C ... Matrix element cPkL from radial matrix elements
      allocate(zgradm(nlg**2,nlg**2,3,nrsmo,nrsmo,nsp*nspc,nbas))
      call dpzero(zgradm,2*size(zgradm))
      do  ispc = 1, nspc
        ksp = max(ispc,isp)
        i = 30 + morder
C       i = i+ 300  ! Return gradient of x+iy, x-iy
        call gradme(i,nrsmo,nrsmo,ksp,nsp,nspc,1,nbas,nlg-1,lmxa,nlmax,
     .    s_optic%rgrade,zgradm,zgradm)
      enddo

C ... Subtract electric dipole matrix elements from local part of envelope functions
C     allocate(optme(nfilm,nempm,3)); call dpzero(optme,2*3*nfilm*nempm)
      call optdme(1,nrsmo,nrfn,nlg,nbas,isp,nsp,nspc,nev,ipc,
     .  zgradm,nfilo,nfiup,nemlo,nemup,1d0,ccp,optme)
C     if (debug) call zprm('optme(1)-optme(2)',2,optme,nfilm,nfilm,nempm)

      deallocate(ccp,zgradm)
      deallocate(lmxa,kmxa,ipc,nrfn,nlmi)

  999 continue
      call tcx('fpopm')
      return

      entry fpopm_sfac(job,s_lat,s_optic,ldsdk,dsdk,iq,ndham,isp,nsp,nspc,nbas,nev)

      if (.not. ldsdk) return

C ... Accumulate
      if (mod(job,2) == 1) then

        call opt_nstate(s_optic,3,nev,nfilo,nfiup,nemlo,nemup,nfilm,nempm)

        cfmeani = (0d0,0d0); frmsi = 0d0; maxfaci = 0d0; ncfi = 0
        ispc = min(isp,nspc)
        do  ipol = 1, 3
          do  ivb = nfilo, nfiup
            do  icb = nemlo, nemup
              if (abs(s_optic%optme(ivb,icb,ipol)) >= 1d-5) then
                cf = 1 + (s_lat%alat*im/(2d0*pi))*dsdk(ivb,icb,ipol,ispc)/s_optic%optme(ivb,icb,ipol)
                cfmeani = cfmeani + cf
                maxfaci = max(maxfaci,abs(cf))
                frmsi = frmsi + dconjg(cf)*cf
                ncfi = ncfi+1
              endif
            enddo
          enddo
        enddo
        cfmean = cfmean + cfmeani
        maxfac = max(maxfac,maxfaci)
        frms = frms + frmsi
        ncf = ncf + ncfi
      endif

C ... Printout
      if (mod(job/2,2) == 1) then
        if (iq. gt. 0 .and. mod(job,2) == 1) then
          if (ncfi > 0) then
          xx = dsqrt(max(frmsi/ncfi-abs(cfmeani/ncfi)**2,0d0))
          call info5(20,0,0,
     .      ' fpopm,  iq=%i: mean scale f = %s,(%2,3;3d) +/- %;3d RMS for %i elements.  max = %;3d',
     .      iq,cfmeani/ncfi,xx,ncfi,maxfaci)
          endif
        elseif (ncf > 0 .and. iq. eq. 0) then
          xx = dsqrt(max(frms/ncf-abs(cfmean/ncf)**2,0d0))
          call info5(20,0,0,
     .      ' fpopm : mean scale f = %s,(%2,3;3d) +/- %;3d RMS for %i elements.  max = %;3d',
     .      cfmean/ncf,xx,ncf,maxfac,5)
        endif

      endif

      end

      subroutine pvfpopm(nfn,nbas,nlg,nlmax,ndham,lmxa,nspc,nev,aus,ccp)
C- Copies aus to ccp
C ----------------------------------------------------------------------
Ci Inputs
Ci   nfn   :number of radial functions, dimensions ccp and aus
Ci   nbas  :size of basis
Ci   nlg   :1 + l-cutoff for partial waves
Ci   nlmax :dimensions aus
Ci   ndham :dimensions aus
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nev   :actual number of eigenvectors generated
Ci   aus  :Coefficients to copy (and rotate to spherical harmonics)
Co Outputs
Ci   ccp   :Essentially aus, but coefficients ordered differently
Cu Updates
Cu   31 Dec 13 Adapted from asa pvopm
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nev,nbas,nlg,nlmax,ndham,nspc,nfn,lmxa(nbas)
      double complex ccp(0:nlg*nlg-1,nbas,nspc,nfn,nev),
     .  aus(0:nlmax-1,ndham*nspc,nfn,nspc,nbas)
C ... Local parameters
      integer ib,ibas,ll,jsp,ifn,nli

C --- For each each eigenvector, do ---
      do  ib = 1, nev

C --- ccp <- aus transformed to regular spherical harmonics Y_lm ---
      do  ibas = 1, nbas
      do  ifn = 1, nfn
      do  jsp= 1, nspc
        nli = lmxa(ibas)+1
        forall (ll = 0:nli**2-1) ccp(ll,ibas,jsp,ifn,ib) = aus(ll,ib,ifn,jsp,ibas)
        forall (ll = nli**2:nlg**2-1) ccp(ll,ibas,jsp,ifn,ib) = 0
      enddo                     ! noncollinear spins
      enddo                     ! kind of partial wave
      enddo                     ! sites
      enddo ! loop over eigenvectors
      end
