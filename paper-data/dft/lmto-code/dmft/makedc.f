      subroutine makedc(s_dmft,dmftu,ndham,nlohi,orbmax,jsp,nsp,nindo,nicix,sigqp,sigdcmk)
C- Compute the local DC correction obtained as a projection of the QSGW self-energy (see L.Sponza's notes)
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft :struct for dmft interface; see structures.h
Ci Inputs
Ci   dmftu  :projector to local correlated Hilbert space
Ci   ndham  :dimensioning parameter, largest hamiltonian dimension
Ci   nlohi  :parameter controlling first and last band to include
Ci   orbmax  :maximum number of correlated channels in a subblock: dimensions dmftu,sidcmk
Ci   jsp    :spin index for eigenvalues
Ci   nsp    :dimensioning parameter,2 for spin-polarized case, otherwise 1
Ci   nindo  :number of orbitals per channel
Ci   nicix  :dimensioning parameter, number of correlated l/site channels
Ci   sigqp  :QSGW self-en obtained from subtraction of Vxc and Hartree en
Co Outputs
Co   sigdcmk :double-counting self-energy k-resolved on local space
Cr Remarks
Cu Updates
Cu   30 Mar 2015 P Pisanti first created
C ----------------------------------------------------------------------
!      use mod_paolo
      use structures
      implicit none
!      include 'structures.h'
      type(str_dmft)::  s_dmft
C ... Passed parameters
      integer, intent(in) :: nicix,orbmax,ndham,nlohi(2),jsp,nsp
!      real(8), intent(in) :: broad
      integer, intent(in) :: nindo(nicix)
      complex(8), intent(in) :: sigqp(ndham,ndham,nsp)
      complex(8), intent(in)  :: dmftu(ndham,orbmax,nicix)
      complex(8), intent(out) :: sigdcmk(orbmax,orbmax,nicix,nicix)
C ... Dynamically allocated local arrays
      complex(8),allocatable :: tmp1(:,:)
C ... Local parameters
      integer :: iorb1,iorb2,nind1,nind2,i,j,nev0,nev
      complex(8) :: sigqpij(ndham,ndham)
C ... External calls
      external dpzero,zgemm

      nev0 = nlohi(1)-1
      nev = nlohi(2)-nev0
      allocate(tmp1(orbmax,nev))
C --- DC sigma (eigenfunction basis) from QSGW self-en in energy subspace ---
      call dpzero(sigqpij,2*size(sigqpij))
      do  i = 1, nev
        do j = 1, nev
          sigqpij(i,j) = sigqp(i+nev0,j+nev0,jsp)
!          if (abs(dimag(sigqpj(i,j))) < broad) then
!           sigqpj(i,j) = dcmplx(dble(sigqpj(i,j)),broad)! Set minimum broadening for correlated orbitals
!          endif
        enddo
      enddo
C --- Local projections ---
      do iorb1 = 1, nicix
        nind1 = nindo(iorb1)
        do iorb2 = 1, nicix
          nind2 = nindo(iorb2)
C     ... sigdcloc ... zgemm call generating tmp1 : a very expensive way to make eval * dmftu ... rewrite as level 2 BLAS
          call zgemm('C','N',nind1,nev,nev,(1d0,0d0),dmftu(1,1,iorb1),ndham,sigqpij,ndham,(0d0,0d0),tmp1,orbmax)
          call zgemm('N','N',nind1,nind2,nev,(1d0,0d0),tmp1,orbmax,dmftu(1,1,iorb2),ndham,(0d0,0d0),sigdcmk(1,1,iorb1,iorb2),orbmax)
        enddo
      enddo
      deallocate(tmp1)
      end subroutine makedc

