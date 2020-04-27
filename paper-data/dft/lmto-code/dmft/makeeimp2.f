      subroutine makeeimp2(s_dmft,ndsig,ncix,ldcix,ndhamx,nlo,nhi,nsp,nspx,ef,mode,
     .        evlk,dmftu,sigbarinftyk,siginpinfty,sigdc_infty,eimpm)
C- Compute the local k-resolved impurity level according to formula 5.14 of Sponza's notes
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft
Ci     Elts read:  nzsig icix l iasig
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   ndsig  :maximum number of matrix elements in any one DMFT block
Ci   ncix   :dimensioning parameter, number of correlated l/site channels
Ci   ldcix  :dimensioning parameter, max number of correlated orbitals (lms character)
Ci   ndhamx :dimensioning parameter, largest hamiltonian dimension
Ci   nlo    :first eigenstate in hybridization window
Ci   nhi    :last  eigenstate in hybridization window
Ci   nsp    :dimensioning parameter,2 for spin-polarized case, otherwise 1
Ci   nspx   :2 for spin-polarized case and collinear moments; otherwise 1
Ci   ef     :chemical potential (updated for current iteration)
Ci   mode   :decide which definition of Eimp to be used.
Ci   evlk   :eigenvalues from the QSGW/LDA Hamiltonian
Ci   dmftu  :projector to local correlated Hilbert space
Ci   sigbarinftyk :high-frequency limit of sigbar in lattice representation {ijk}
Ci   siginpinfty  :high-frequency limit of sig,inp in local representation {LL'}
Ci   sigdc_infty  :sigDC_{LL'} at omega=infty.
Co Outputs
Co   eimpm :impurity levels in local representation
Cr Remarks
Cu Updates
Cu   18 May 18 (MvS) Enable equivalent cix blocks with spin flip
Cu   30 Apr 18 (MvS) First cut at equivalent cix blocks
Cu   17 Jan 18 (MvS) First cut at SO coupled case
Cu   27 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu    9 Nov 15 (Sponza) new version to get sigbar in input
Cu   22 Aug 15 (Pisanti) spin polarized
Cu   19 Oct 14 (Pisanti) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_dmft)::  s_dmft
C ... Passed parameters
      integer,    intent(in)  :: ndsig,ncix,ldcix,ndhamx,nlo,nhi,nsp,nspx,mode
      real(8),    intent(in)  :: ef
      real(8),    intent(in)  :: evlk(ndhamx,nsp)
      complex(8), intent(in)  :: dmftu(nlo:nhi,ldcix,nsp,ncix)
      complex(8), intent(in)  :: sigbarinftyk(nlo:nhi,nlo:nhi,nspx)
      complex(8), intent(in)  :: sigdc_infty(*)
      real(8),    intent(in)  :: siginpinfty(2*ndsig)
      complex(8), intent(out) :: eimpm(ldcix,ldcix,nsp,ncix)
C ... Dynamically allocated local arrays
      complex(8),allocatable :: tmp1(:,:),evldelta(:,:,:)
C ... Local parameters
      integer :: isp,ksp,i,j,nl1,icix,cix,cixi,ichan,m1,m2,isigind,nzsig,nev,nicix(ncix)
      complex(8) :: eimptmp(ldcix,ldcix,nsp)
C ... External calls
      external dpzero,zgemm

      nev = nhi-nlo+1
      allocate(evldelta(nlo:nhi,nlo:nhi,nspx))
      nzsig = s_dmft%nzsig

C --- Eigenvalues (plus sigbar) in lattice representation ---
      call dpzero(evldelta,2*size(evldelta))
      do  isp = 1, nspx
        do  i = nlo, nhi
          evldelta(i,i,isp) = evlk(i,isp)
          if (mode==1) then
            do  j = nlo, nhi
              evldelta(i,j,isp) = evldelta(i,j,isp) + sigbarinftyk(i,j,isp)
            enddo
          endif
        enddo
C       call zprm('sigbarinftyk',2,sigbarinftyk(nlo,nlo,isp),nhi-nlo+1,nhi-nlo+1,nhi-nlo+1)
      enddo

C --- Add local projection of evldelta to eimpm ---
      allocate(tmp1(ldcix,nev))
      do  isp = 1, nsp
        ksp = min(isp,nspx)
        do  cix = 1, ncix
          icix = iabs(s_dmft%icix(cix))
          nl1 = 2*s_dmft%l(icix)+1

          call zgemm('C','N',nl1,nev,nev,(1d0,0d0),
     .               dmftu(nlo,1,isp,cix),nev,
     .               evldelta(nlo,nlo,ksp),nev,(0d0,0d0),
     .               tmp1,ldcix)

          call zgemm('N','N',nl1,nl1,nev,(1d0,0d0),
     .               tmp1,ldcix,
     .               dmftu(nlo,1,isp,cix),nev,(1d0,0d0),
     .               eimpm(1,1,isp,cix),ldcix)

        enddo
      enddo
      deallocate(tmp1)

C --- Subtract impurity self-energy (in Rydberg) and chemical potential
C     For sake of simplicity this is done for all k-points and then everything is divided by nkfbz
      do  ichan = 1, nzsig
        m1  = s_dmft%iasig(1,ichan)
        m2  = s_dmft%iasig(2,ichan)
        cix  = s_dmft%iasig(5,ichan)
        isp = s_dmft%iasig(3,ichan)
C        if (s_dmft%icix(cix) < 0 .and. nspx == 2) isp = nspx+1-isp
        isigind = s_dmft%iasig(4,ichan)
        icix = iabs(s_dmft%icix(cix))

C       Exact definition for the 1-shot DMFT scheme (eq. 5.14)
        if (mode==1) eimpm(m1,m2,isp,cix) = eimpm(m1,m2,isp,cix)
     .    - dcmplx(siginpinfty(2*isigind-1),siginpinfty(2*isigind))

C       Approximate definition for the 1-shot DMFT scheme (eq.5.16)
        if (mode==2) eimpm(m1,m2,isp,cix) = eimpm(m1,m2,isp,cix) - sigdc_infty(isigind)

        if (m1 == m2 .and. (mode == 1 .or. mode == 2)) eimpm(m1,m2,isp,cix) = eimpm(m1,m2,isp,cix) - ef

      enddo

C     Average equivalent cix sites
      call ineqcix(s_dmft,ncix,nicix)
      do  cixi = 1, ncix            ! Loop over all inequivalent cixi
        if (nicix(cixi) >= 0) cycle ! nicix(:)<0 => first of equivalent cix
        icix = iabs(s_dmft%icix(cixi)) ! inequivalent cix index
        j = ldcix**2
        call dpzero(eimptmp,2*j*nsp)
        do  cix = 1, ncix           ! Loop over cix equivalent to cixi: gather average
          if (iabs(s_dmft%icix(cix)) /= icix) cycle ! Skip blocks not equivalent to this one
          do  isp = 1, nsp
            ksp = isp
            if (s_dmft%icix(cix) < 0 .and. nsp == 2) ksp = nsp+1-isp ! spin flip for AFM symmetry
            call daxpy(2*j,1/dble(iabs(nicix(cixi))),eimpm(1,1,ksp,cix),1,eimptmp(1,1,isp),1)
          enddo
        enddo
        do  cix = 1, ncix           ! Loop over cix equivalent to cixi: distribute average
          if (iabs(s_dmft%icix(cix)) /= icix) cycle ! Skip blocks not equivalent to this one
          do  isp = 1, nsp
            ksp = isp
            if (s_dmft%icix(cix) < 0 .and. nsp == 2) ksp = nsp+1-isp ! spin flip for AFM symmetry
            call dcopy(2*j,eimptmp(1,1,ksp),1,eimpm(1,1,isp,cix),1)
          enddo
        enddo
      enddo

C      do  cix = 1, ncix
C        call yprmi('eimpm, cix=%i',cix,0,3,eimpm(1,1,1,cix),0,ldcix,ldcix,ldcix*nsp)
C      enddo

      end subroutine makeeimp2
