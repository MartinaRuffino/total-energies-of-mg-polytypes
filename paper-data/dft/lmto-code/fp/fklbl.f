      subroutine fklbl(p,rsm,kmax,nlm,k0,cy,s_lat,fkl)
C- Bloch sum of smooth Hankels for e=0 and q=(0,0,0).
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
Ci   kmax  :polynomial cutoff
Ci   nlm   :L-cutoff for gkl
Ci   k0    :leading dimension of gkl
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   fkl   :Bloch-summed Hankels for q=0 and e=0
Cr Remarks
Cr   For (k=0,l=0) f equals the limit of hklbl minus the avg value.
Cr   For all other cases f is the limit of hklbl as e goes to zero.
Cu Updates
Cu   02 Jul 15 Remove fixed dimensioning nlm0
Cu   10 Nov 11 Begin migration to f90 structures
Cu   24 Apr 00 Adapted from nfp fkl_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer kmax,nlm,k0
      double precision p(3),cy(*),rsm
      double complex fkl(0:k0,nlm)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      real(8),allocatable:: wk(:),yl(:)
C ... Local parameters
      integer lmax,ll,nkd,nkq,nrx,job,ilm,k
C     parameter ( nlm0=196 )
      double precision q(3),alat,plat(3,3),qlat(3,3),p1(3)
      double precision faca,fpi,y0,e,vol,awald,tol
      double complex fsm(nlm),gklsav,gklnew
      parameter (faca=1d0)

      if (nlm == 0) return
      fpi = 16d0*datan(1d0)
      y0 = 1d0/dsqrt(fpi)
C     if (nlm > nlm0) call rx('fklbl: increase nlm0')
      lmax = ll(nlm)
      e = 0d0
      q(1) = 0d0
      q(2) = 0d0
      q(3) = 0d0

C ... Use standard routines
C     call gklbl(p,rsm,e,q,kmax-1,nlm,k0,cy,s_lat,fkl)
C     call fsmbl(p,rsm,lmax,cy,slat, fsm)

C ... Alternatively, use vectorized equivalents (about 2x faster)
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      awald = s_lat%awald
      tol = s_lat%tol
      vol = s_lat%vol
      nkd = s_lat%nkd
      nkq = s_lat%nkq
      call shorbz(p,p1,plat,qlat)
      nrx = max(nkd,nkq)
      allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
      call hsmqe0(lmax,rsm,0,q,p1,nrx,nlm,wk,yl,
     .  awald,alat,s_lat%qlv,nkq,s_lat%dlv,nkd,vol,fsm)
      if (rsm > faca/awald) then
        call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,s_lat,fkl)
      else
        job = 2
        call gklq(lmax,rsm,q,p1,e,kmax-1,k0,alat,s_lat%dlv,nkd,nrx,yl,wk,job,fkl)
      endif
      deallocate(wk,yl)

C ... Upward recursion in k: mainly sets fkl = -4*pi * g(k-1,l)
      do  ilm = 1, nlm
        gklsav = fkl(0,ilm)
        fkl(0,ilm) = fsm(ilm)
        do  k = 1, kmax
          gklnew = fkl(k,ilm)
          fkl(k,ilm) = -fpi*gklsav
          gklsav = gklnew
        enddo
      enddo

C ... Add extra term to F(k=1,l=0)
      fkl(1,1) = fkl(1,1) + fpi*y0/vol

      end
