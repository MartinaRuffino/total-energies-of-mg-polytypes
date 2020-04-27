      subroutine optsf(mode,nspx,nkp,efermi,evnl,evloc,
     .  nfilo,nfiup,nemlo,nemup,optmt,optmc)
C- Scale the optical matrix for nonlocal self energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 scale optmt
Ci         :2 scale optmc
Ci         :3 scale combination
Cl   nspx  :number of independent spin channels
Cl         :nspx=nsp unless noncollinear (nspc=2), in which case nspx=1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   evnl  :QSGW eigenvalues (For correction factor with NL potential)
Ci   evloc :LDA  eigenvalues (For correction factor with NL potential)
Ci   nfilo,nfiup :range of occupied states which have elements in optmt
Ci   nemlo,nemup :range of unocc states which have elements in optmt
Cl Local variables
Cl   fac   :NL matrix element scale factor calculated from evnl and evloc
Co Outputs
Co   optmt(*,ikp): <i|grad|j>^2 scaled by fac^2
Co   optmc(*,ikp): <i|grad|j> scaled by fac
Co   nsmallf     : number of k points where fac is < 1
Cr Remarks
Cr   Details of the method is provided in: Phys. Rev. B 48, 11789
Cu Updates
Cu   17 Jun 17 (B. Cunningham) optionally scales optmc
Cu   29 Mar 16 (P. Azarhoosh) restructure for reduced dimensions in evnl, evloc
Cu   26 Aug 14 Redesigned for the spin polarized and noncollinear cases
Cu   12 May 14 (Pooya Azarhoosh) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nspx,nkp,nfilo,nfiup,nemlo,nemup
      double precision optmt(3,nfilo:nfiup,nemlo:nemup,nspx,nkp),
     .                 evnl(nfilo:nemup,nspx,nkp),evloc(nfilo:nemup,nspx,nkp),
     .                 efermi
C ... Local parameters
      integer i,j,isp,nsmallf,nbigf,ikp,n
      double precision fac,tol1,tol2,dlength,maxfac,om1,om2,fmean,frms
      parameter (tol1=1d-8,tol2=1d-3,maxfac=5d0)
      complex(8) optmc(3,nfilo:nfiup,nemlo:nemup,nspx,nkp)
      nsmallf = 0; nbigf = 0; fmean = 0; frms = 0; n = 0

      if (mode == 0) return

C --- Loop over all k and spin  ---
      do  ikp = 1, nkp
      do  isp = 1, nspx

C   ... For each (occ, unocc) pair do
        do  i = nfilo, nfiup
        do  j = nemlo, nemup

          if (dlength(3,optmt(1,i,j,isp,ikp),1) < tol1) cycle ! 0 matrix elmement
          om1 = evloc(i,isp,ikp)-efermi
          om2 = evloc(j,isp,ikp)-efermi
          if(om1 <= 0.and.om2 <= 0 .or. om1 >= 0.and.om2 >= 0) cycle !Skip if same side of Fermi level
          if (om2-om1 < tol2) cycle ! Too similar: scale factor nonsensical
          fac = (evnl(j,isp,ikp)-evnl(i,isp,ikp))/
     .          (evloc(j,isp,ikp)-evloc(i,isp,ikp))
          fac = min(fac,maxfac)  ! Scale factor not reasonable
          if (fac < 1) nsmallf= nsmallf+1
          if (fac == maxfac) nbigf= nbigf+1
          if (fac < 1 .or. fac == maxfac) then
            call info5(50,0,0,' optsf: occ,unocc=%i,%i  iq=%i'//
     .        '  fac=%d',i,j,ikp,fac,0)
          endif
          fmean = fmean + fac**2; n = n+1
          frms  = frms  + fac**4
          if (mod(mode,2) == 1) then
            optmt(1:3,i,j,isp,ikp) = optmt(1:3,i,j,isp,ikp)*fac**2
          endif
          if (mod(mode/2,2) == 1) then
            optmc(1:3,i,j,isp,ikp) = optmc(1:3,i,j,isp,ikp)*fac
          endif
        enddo
        enddo
      enddo
      enddo
      if (n > 0) then
        frms = frms/n
        fmean = fmean/n
        frms = dsqrt(max(frms-fmean**2,0d0))
        call info5(20,0,0,
     .    ' optsf: mean scale f = %;3d +/- %;3d RMS.'//
     .  '  f<1 for %i pairs; f>%d for %i pairs',
     .    fmean,frms,nsmallf,maxfac,nbigf)
      endif
      end

      subroutine optpmme(mode,ndhamx,nspx,nkp,iprm,eval,nfilo,nfiup,nemlo,nemup,optmt,optmc)
C- Reorder the optical matrix elements according to a permuatation list
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 permute optmt
Ci         :2 permute optmc
Ci         :4 permute eval
Ci         :Any combination is allowed
Ci   ndhamx:leading dimension of eval
Cl   nspx  :number of independent spin channels
Cl         :nspx=nsp unless noncollinear (nspc=2), in which case nspx=1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   iprm  :permutation matrix
Ci   eval  :eigenvalues
Ci   nfilo,nfiup :range of occupied states which have elements in optmt
Ci   nemlo,nemup :range of unocc states which have elements in optmt
Co   optmt : matrix elements of the square of the velocity operator
Co   optmc : matrix elements of the velocity operator
Co Outputs
Co   optmt(:,:,:,ikp): permuted, dependending on mode
Co   optmc(:,:,:,:,ikp): permuted, dependending on mode
Co   eval(:,:,:) permuted, dependending on mode
Cr Remarks
Cu Updates
Cu   24 Mar 18 (MvS) Adapted from optsf.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ndhamx,nspx,nkp,nfilo,nfiup,nemlo,nemup,iprm(nkp)
      double precision optmt(3,nfilo:nfiup,nemlo:nemup,nspx,nkp),
     .                 eval(ndhamx,nspx,nkp)
      complex(8) optmc(3,nfilo:nfiup,nemlo:nemup,nspx,nkp)
C ... Dynamically allocated arrays
      real(8), allocatable :: evll(:,:,:),optmtl(:,:,:,:,:)
      complex(8), allocatable :: optmcl(:,:,:,:,:)
C ... Local parameters
      integer ikp

      if (mode == 0) return

C     n=0 for evals, n=1 for optmt, n=2 for optmc

      if (btest(mode,0)) then ! 1s bit
        allocate(evll(ndhamx,nspx,nkp))
        do  ikp = 1, nkp
            call dcopy(ndhamx*nspx,eval(1,1,iprm(ikp)),1,evll(1,1,ikp),1)
        enddo
        call dcopy(ndhamx*nspx*nkp,evll,1,eval,1)
        deallocate(evll)
      end if

      if (btest(mode,1)) then ! 2s bit
        allocate(optmtl(3,nfilo:nfiup,nemlo:nemup,nspx,nkp))
        do  ikp = 1, nkp
            call dcopy(size(optmt(:,:,:,:,1)),optmt(1,nfilo,nemlo,1,iprm(ikp)),1,optmtl(1,nfilo,nemlo,1,ikp),1)
        enddo
        call dcopy(size(optmtl),optmtl,1,optmt,1)
        deallocate(optmtl)
      end if

      if (btest(mode,2)) then ! 4s bit
        allocate(optmcl(3,nfilo:nfiup,nemlo:nemup,nspx,nkp))
        do  ikp = 1, nkp
            call zcopy(size(optmc(:,:,:,:,1)),optmc(1,nfilo,nemlo,1,iprm(ikp)),1,optmcl(1,nfilo,nemlo,1,ikp),1)
        enddo
        call zcopy(size(optmcl),optmcl,1,optmc,1)
        deallocate(optmcl)
      end if

!
!       do  n = 0, 2
!         k = 2**n
!
!         if (iand(mode,1) == k) allocate(evll(ndhamx,nspx,nkp))
!         if (iand(mode,2) == k) allocate(optmtl(3,nfilo:nfiup,nemlo:nemup,nspx,nkp))
!         if (iand(mode,4) == k) allocate(optmcl(3,nfilo:nfiup,nemlo:nemup,nspx,nkp))
!
! C   --- Loop over all k and spin  ---
!         do  ikp = 1, nkp
!           if (iand(mode,1) == k) then
!             call dcopy(size(eval(:,:,1)),eval(1,1,iprm(ikp)),1,evll(1,1,ikp),1)
!           endif
!           if (iand(mode,2) == k) then
!             call dcopy(size(optmt(:,:,:,:,1)),optmt(1,nfilo,nemlo,1,iprm(ikp)),1,optmtl(1,nfilo,nemlo,1,ikp),1)
!           endif
!           if (iand(mode,4) == k) then
!             call zcopy(size(optmc(:,:,:,:,1)),optmc(1,nfilo,nemlo,1,iprm(ikp)),1,optmcl(1,nfilo,nemlo,1,ikp),1)
!           endif
!         enddo
!
!         if (iand(mode,1) == k) then
!           call dcopy(size(eval),evll,1,eval,1)
!           deallocate(evll)
!         endif
!         if (iand(mode,2) == k) then
!           call dcopy(size(optmtl),optmtl,1,optmt,1)
!           deallocate(optmtl)
!         endif
!         if (iand(mode,4) == k) then
!           call zcopy(size(optmcl),optmcl,1,optmc,1)
!           deallocate(optmcl)
!         endif
!
!       enddo

      end
