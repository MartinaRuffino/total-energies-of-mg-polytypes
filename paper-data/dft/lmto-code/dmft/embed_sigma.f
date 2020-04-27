      subroutine embed_sigma(s_dmft,nlo,nhi,ldcix,nsp,nspx,
     .                       dmftu,sigdc,siginp,sigbarijk)
C- Compute Sigmabar_{ijk}(w) as the embedding of (Sig.inp(w) - Sig.dc(w))_{LL'}
C -----------------------------------------------------------------------------
Cio Structures
Cio  s_dmft
Ci     Elts read:  nzsig ndsig iasig
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  compressproj
Ci Inputs
Ci   ndham  :dimensioning parameter, largest hamiltonian dimension
Ci   ldcix  :dimensioning parameter, max number of correlated orbitals (lms character)
Ci   nicix  :dimensioning parameter, number of correlated l/site channels
Ci   nzsig  :dimensioning parameter, total number of elements of self-energy matrix
Ci          :in any block
Ci   nsp    :dimensioning parameter,2 for spin-polarized case, otherwise
Ci   nspx   :2 for spin-polarized case and collinear moments; otherwise 1
Ci   dmftu  :projector to local correlated Hilbert space
Ci   sigdc  :Double counting self-energy in local representation (can be freq.-dependent)
Ci   siginp :Impurity self-energy in local representation {LL'}
Ci   sigbarijk : (siginp-sigdc) embedded in the lattice, basis of one-body crystal eigenstates.
Cl Local variables
Cb Bugs
Cb   This code does not consider rotation of sigma, I think
Cr Remarks
Cr  This can handle also energy-dependent sigdc
Cr  The use of embed_sigma requires each call to be stored in a very
Cr    large matrix sigbar(:,:,:,:,:). Consider writing on disk.
Cu Updates
Cu   17 Jan 18 (MvS) First cut at SO coupled case
Cu   27 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu      Fal 15 Sponza first created
C -----------------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!       include 'structures.h'
      type(str_dmft)         :: s_dmft
C ... Passed parameters
      integer,    intent(in) :: nlo,nhi,ldcix,nsp,nspx
      complex(8), intent(in) :: dmftu(nlo:nhi,ldcix,*)
      complex(8), intent(in) :: sigdc(*)
      real(8),    intent(in) :: siginp(*)
      complex(8), intent(out):: sigbarijk(nlo:nhi,nlo:nhi,nsp)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: strans(:,:,:) ! Projector
      complex(8), allocatable :: sigma(:) ! sigma with double counting removed
      complex(8), allocatable :: sigmabar(:) !
C ... Local parameters
      integer :: ichan,isigind,nzsig,ndsig,iv,jv,isp

      call tcn('embed_sigma')

      nzsig = s_dmft%nzsig
      ndsig = s_dmft%ndsig
      allocate(strans(nzsig,nlo:nhi,nlo:nhi),sigma(ndsig),sigmabar(nzsig))

C --- Compress the projector to 'embed' the self-energy ---
      call compressproj(s_dmft,dmftu,nlo,nhi,ldcix,nsp,nzsig,strans)

C      print "('checksum strans',4f15.10)",sum(strans(:,:,:))
C      call zprm('strns(3)',2,strans(3,:,:),nhi-nlo+1,nhi-nlo+1,nhi-nlo+1)

C --- Make sigmabar = siginp - sigdc ---
      call dpzero(sigma,2*size(sigma))
      do  ichan = 1, nzsig
        isigind = s_dmft%iasig(4,ichan)
        sigmabar(ichan) = dcmplx(siginp(2*isigind-1),siginp(2*isigind)) - sigdc(isigind)
C       print *, ichan, isigind, sigmabar(ichan)
      enddo

C --- Embed sigmabar into lattice space (iv,jv,spin)  per each k and omega ---
      call dpzero(sigbarijk,2*size(sigbarijk(:,:,1:nspx)))
      do  iv = nlo, nhi                ! over bands-1
        do  jv = nlo, nhi              ! over bands-2
          do  ichan = 1, nzsig
            isp = min(nspx,s_dmft%iasig(3,ichan))
C            cix = s_dmft%iasig(5,ichan)
C            if (s_dmft%icix(cix) < 0 .and. nspx == 2) isp = nspx+1-isp
            sigbarijk(iv,jv,isp) = sigbarijk(iv,jv,isp) + strans(ichan,iv,jv)*sigmabar(ichan)
          enddo
        enddo
      enddo

C      call zprm('sigbarijk',2,sigbarijk,nhi-nlo+1,nhi-nlo+1,nhi-nlo+1)
C      print "('checksum sigbarijk',4f15.10)",sum(sigbarijk(:,:,1)),sum(sigbarijk(:,:,nspx))
C      stop

      deallocate(strans,sigma,sigmabar)
      call tcx('embed_sigma')

      end subroutine embed_sigma
