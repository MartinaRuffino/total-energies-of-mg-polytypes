      subroutine hklbl(p,rsm,e,q,kmax,nlm,k0,cy,s_lat,hkl)
C- Bloch-sums of k,L-dependent smooth Hankel functions.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat awald tol vol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qlv dlv
Cio    Passed to:  gklbl
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsm   :smoothing radius
Ci   e     :energy of smoothed Hankel
Ci   q     :wave number for Bloch sum
Ci   kmax  :polynomial cutoff
Ci   nlm   :L-cutoff for hkl
Ci   k0    :leading dimension of hkl
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   hkl   :Bloch-summed smoothed Hankels
Cr Remarks
Cr   H_kL = laplace^k H_L
Cr   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   24 Apr 00 Adapted from nfp hkl_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k0,kmax,nlm
      double precision e,rsm,q(3),p(3),cy(*)
      double complex hkl(0:k0,nlm)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer ilm,job,k,ll,lmax,nkd,nkq,nrx
C     parameter (nlm0=144)
      double precision alat,awald,fpi,pi,sp,tol,vol,plat(3,3),
     .  qlat(3,3),p1(3)
      double complex hsm(nlm),hsmp(nlm),phase,gklsav,gklnew
C     double complex hsmx(nlm),hsmpx(nlm)

      real(8),allocatable:: wk(:),yl(:)
      double precision faca
      parameter (faca=1d0)

      if (nlm == 0) return
C     if (nlm > nlm0) call rxi('increase nlm0 in hklbl need',nlm)

      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      lmax = ll(nlm)

C ... Shorten p, multiply by corresponding phase later
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      call shorbz(p,p1,plat,qlat)
      sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
      phase = dcmplx(dcos(sp),dsin(sp))

C ... Use standard routines
C     call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,s_lat,hkl)
C     call hsmbl(p1,rsm,e,q,lmax,cy,s_lat,hsm,hsmp)

C ... Alternatively, use vectorized equivalents (about 2x faster)
      awald = s_lat%awald
      tol = s_lat%tol
      vol = s_lat%vol
      nkd = s_lat%nkd
      nkq = s_lat%nkq
      nrx = max(nkd,nkq)
      allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
      call hsmq(1,0,ll(nlm),e,rsm,0000,q,p1,nrx,nlm,wk,yl,
     .  awald,alat,s_lat%qlv,nkq,s_lat%dlv,nkd,vol,hsm,hsmp)
      if (rsm > faca/awald) then
        call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,s_lat,hkl)
      else
        job = 2
        call gklq(lmax,rsm,q,p1,e,kmax-1,k0,alat,s_lat%dlv,nkd,nrx,
     .    yl,wk,job,hkl)
      endif
      deallocate(wk,yl)

C  ... for debugging
C       do  ilm = 1, nlm
C        if (cdabs(hsm(ilm)-hsmx(ilm)) > 1d-10) then
C          if (cdabs(hsm(ilm)/hsmx(ilm)-1) > 1d-8) then
C     .      .and. cdabs(hsm(ilm)) > .1d0
C     .      .and. cdabs(hsm(ilm)) < 1.d0) then
C          print *, 'error at ilm=',ilm
C          print *, hsm(ilm)
C          print *, hsmx(ilm)
C          print *, '---'
C          stop
C        endif
C        endif
C       enddo

C --- H_kL by recursion ---
      do  ilm = 1, nlm
      gklsav = hkl(0,ilm)
      hkl(0,ilm) = hsm(ilm)
      do  k = 1, kmax
        gklnew = hkl(k,ilm)
        hkl(k,ilm) = -e*hkl(k-1,ilm) - fpi*gklsav
        gklsav = gklnew
        enddo
      enddo

C ... Put in phase to undo shortening
      if (sp /= 0) then
        do  ilm = 1, nlm
          do  k = 0, kmax
            hkl(k,ilm) = phase*hkl(k,ilm)
          enddo
        enddo
      endif

      end
