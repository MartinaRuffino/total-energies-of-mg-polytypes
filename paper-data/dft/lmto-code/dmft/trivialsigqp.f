      subroutine trivialsigqp(ndham,nlohi,jsp,nsp,broad,evlk,sigqp)
C- Set the QSGW selfen (which is supposed to be non-diagonal) trivially to diagonal epsilon_qp = e_hartree + sigma_qsgw
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndham  :dimensioning parameter, largest hamiltonian dimension
Ci   nlohi  :
Ci   jsp    :
Ci   nsp    :2 for spin-polarized case, otherwise 1
Ci   broad  :Gaussian broadening, in Ry
Ci   evlk   :eigenvalue of QSGW hamiltonian
Co Outputs
Ci   sigqp   : QSGW self-en obtained from subtraction of Vxc and Hartree en
Cr Remarks
Cu Updates
Cu    31 Mar 2015
C ----------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: ndham,nlohi(2),jsp,nsp
      real(8), intent(in) :: evlk(ndham,nsp)
      real(8), intent(in) :: broad
C ... Local parameters
      integer :: i,nev0,nev
      complex(8) :: sigqp(ndham,ndham,nsp)
C ... External calls
      external dpzero
      nev0 = nlohi(1)-1
      nev = nlohi(2)-nev0
      call dpzero(sigqp,2*size(sigqp))
      do  i = 1, nev
          sigqp(i+nev0,i+nev0,jsp) = dcmplx(dble(evlk(i+nev0,jsp)),broad)! Set minimum broadening for correlated orbitals
      enddo


      end subroutine trivialsigqp
