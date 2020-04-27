      subroutine dmftevec(s_dmft,s_pade,ndham,nomg,nsp,nspc,nlohi,ldcix,iq,nkp,
     .  procid,kpproc,omega,dmftu,evlk,emu,sigdc,siginp,evecl,evecr,udeval)
C- Find eigenvectors of QP+DMFT hamiltonian for a range of frequencies, in QP basis
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft
Ci     Elts read:  nzsig ndsig iasig
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ndsig nicix iasig
Cio    Passed to:  embed_sigma compressproj
Cio  s_pade
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lpade npade zomega padcof
Cio    Passed to:  mksigpade
Ci Inputs
Ci   ndham  :dimensioning parameter, largest hamiltonian dimension
Ci   nomg   :number of frequency points
Ci   nsp    :2 for spin-polarized case, otherwise 1
Ci   nspc   :2 if spin-up and spin-down channels are coupled; else 1.
Ci   iq     :k-point index, used
Ci   nkp    :number of irreducible k-points (bzmesh.f)
Ci   nlohi  :range of bands included in projector
Ci   ldcix  :dimension of largest cix block, for dimensioning
Ci   procid :index to kpproc
Ci   kpproc :perform calculations for frequencies [kpproc(procid),kpproc(procid+1)]
Ci   omega  :frequency mesh
Ci   dmftu  :projector to local correlated Hilbert space
Ci   sigdc  :Double counting self-energy in local representation (can be freq.-dependent)
Ci   siginp :Impurity self-energy in local representation {LL'}
Ci   emu    :energy - chemical potential : used to interpolate sigma
Co Outputs
Co   evecl  :Left DMFT eigenfunction in the basis of single-particle states
Co   evecr  :Right DMFT eigenfunction in the basis of single-particle states
Co   udeval
Cr Remarks
Cr
Cu Updates
Cu   01 Jun 18 Not tested yet
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ndham,nomg,nsp,nspc,iq,nkp,ldcix,nlohi(2),procid,kpproc(0:procid+1)
      double precision emu,evlk(ndham,nsp),omega(nomg),siginp(nomg,*)
      complex(8), intent(in) :: dmftu(nlohi(1):nlohi(2),ldcix,*),sigdc(*)
      complex(8), intent(out) :: evecl(ndham,ndham,nomg,nsp),evecr(ndham,ndham,nomg,nsp),udeval(ndham*nspc,nomg,nkp,nsp)
C ... For structures
      type(str_dmft)::  s_dmft
      type (str_mpade):: s_pade
C ... Dynamically allocated local arrays
      complex(8), allocatable :: sigbark(:,:,:) ! DMFT Sigma in basis of 1-particle eigenfunctions, one frequency
      complex(8), allocatable :: udhamk(:,:,:)
      complex(8), allocatable :: udevalko(:) ! updated eigenvalues
      real(8), allocatable :: sigpade(:,:,:) ! Pade extrapolation of Self-energy
      complex(8), allocatable :: eveclko(:,:), evecrko(:,:) ! left and right eigenvectors of Hbar at given energy, spin, k-point
C ... Local parameters
      integer iomg,nspx,ndhamx,isp,nevn
C     double precision
      procedure(integer) :: mpipid

      nspx = nsp / nspc         ! 1 if nsp=1 or if noncollinear
      ndhamx = ndham * nspc
      nevn = nlohi(2)-nlohi(1)+1

      allocate(sigbark(nevn,nevn,nspx),udhamk(ndhamx,ndhamx,nspx),udevalko(ndhamx))
      allocate(sigpade(2,s_dmft%ndsig,s_dmft%nicix))
      allocate(eveclko(ndhamx,ndhamx),evecrko(ndhamx,ndhamx))

C ... Create and diagonalise the omega-dependent Hamiltonian
      do  iomg = 1, nomg
        if (iomg<kpproc(procid) .or. iomg>=kpproc(procid+1)) cycle

C       if (iomg > 60 .and. iomg < nomg) cycle

C       Construct Sigmabar in the lattice space
        call dpzero(sigbark,2*size(sigbark))
        call mksigpade(s_pade,1,nomg,s_dmft%ndsig,dcmplx(emu,omega(iomg)),siginp(iomg,1),sigpade)
        call embed_sigma(s_dmft,nlohi(1),nlohi(2),ldcix,nsp,nspx,dmftu,sigdc,sigpade,sigbark)

C       Construct udhamk =  evlk (all states) + sigbark (window) for spins 1:nspx
        call dpzero(udhamk,2*size(udhamk))
        call makehbar2(nevn,nlohi(1)-1,nspx,ndhamx,evlk,sigbark,udhamk)

C       Diagonalize one-particle+DMFT hamiltonian in basis of single-particle states
        do  isp = 1, nspx
          call dpzero(udevalko,2*size(udevalko))
          call dpzero(evecrko,2*size(evecrko))
          call dpzero(eveclko,2*size(eveclko))
          call diagoham_full(ndhamx,udhamk(1,1,isp),udevalko,eveclko,evecrko)
          evecl(:,:,iomg,isp) = eveclko
          evecr(:,:,iomg,isp) = evecrko
          udeval(:,iomg,iq,isp) = udevalko
        enddo
      enddo                     ! frequency loop
      deallocate(udhamk,udevalko,sigbark,sigpade,eveclko,evecrko)

C      iomg = kpproc(procid); isp = 1
C      call info5(1,0,0,'procid,iq,isp,iomg %4:1,4i evecl %3;11,6D evecr %3;11,6D udeval %3;11,6D',
C     .  [procid,iq,isp,iomg],sum(evecl(:,:,iomg,isp)),sum(evecr(:,:,iomg,isp)),sum(udeval(:,iomg,iq,isp)),5)

      end
