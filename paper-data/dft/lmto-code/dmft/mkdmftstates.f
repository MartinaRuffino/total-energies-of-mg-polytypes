      subroutine mkdmftstates(ndhamx,wtkp,evald,rot,evec0,evecd,dmftwt)
C- Return DMFT eigtenvectors or density matrix in orbital basis, and possibly weights.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndhamx:rank of the hamiltonian
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   evald :eigenvalues of the DMFT density-matrix in basis of one-particle eigenstates
Ci   rot   :eigenvectors of DMFT hamiltonian in basis of eigenstates of one-particle hamiltonian
Ci         :or alternatively eigenvectors of the DMFT density matrix; see Remarks
Ci   evec0 :eigenvectors of one-particle hamiltonian, orbital basis
Co Outputs
Co   evecd :eigenvectors of DMFT hamiltonian, orbital basis
Co   dmftwt:DMFT density matrix weighted by qp weight
Cr Remarks
Cr   If rot are the eigenvectors of the DMFT density matrix,
Cr   W^{DMFT}_{ijk} of Haule's PRB 81, 195107 (2010) and
Cr   evald are the corresponding eigenvalues, this routine returns
Cr   eigenvectors that substitute for one-particle eigenvectors in the
Cr   generation of the charge density
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer,    intent(in)  :: ndhamx
      real(8),    intent(in)  :: wtkp,evald(ndhamx)
      complex(8), intent(in)  :: rot(ndhamx,ndhamx),evec0(ndhamx,ndhamx)
      complex(8), intent(out) :: evecd(ndhamx,ndhamx)
      real(8)   , intent(out) :: dmftwt(ndhamx)
C ... Dynamically allocated local arrays
      integer, allocatable :: iprm(:)
C ... Local parameters
      integer :: i
C ... External calls
      external dvheap,dvprm,zgemm

      call zgemm('n','n',ndhamx,ndhamx,ndhamx,(1d0,0d0),evec0,ndhamx,rot,ndhamx,(0d0,0d0),evecd,ndhamx)

C ... Make dmftwt
      do i = 1, ndhamx
        dmftwt(i) = dble(evald(i))*wtkp
      enddo

C ... Order the eigenvectors by increasing wtkp
      allocate(iprm(ndhamx))
      call dvheap(1,ndhamx,dmftwt,iprm,0d0,100)
      forall (i=1:ndhamx) iprm(i) = ndhamx+1-iprm(i)
      call dvprm(1,ndhamx,dmftwt,[wtkp],iprm,2)
      call dvprm(2*ndhamx,ndhamx,evecd,[wtkp],iprm,2)
      deallocate(iprm)

      end

C      subroutine mkdmftstatesz(ndhamx,wtkp,rot,evec0,evecd,evald,dmftwt)
CC- Same as mkdmftstates but eigenvalues are complex
C      implicit none
C      integer,    intent(in)  :: ndhamx
C      real(8),    intent(in)  :: wtkp
C      complex(8), intent(in)  :: rot(ndhamx,ndhamx),evec0(ndhamx,ndhamx),evald(ndhamx)
C      complex(8), intent(out) :: evecd(ndhamx,ndhamx)
C      real(8)   , intent(out) :: dmftwt(ndhamx)
C      integer :: i,j
C      complex(8) :: tmp(ndhamx,ndhamx)
C
C      evecd=cmplx(0d0,0d0)
C      do i=1,ndhamx
C        dmftwt(i)=dble(evald(i))*wtkp
C      enddo
C
C      call zgemm('n','n',ndhamx,ndhamx,ndhamx,(1d0,0d0),evec0,ndhamx,rot,ndhamx,(0d0,0d0),evecd,ndhamx)
C
C      end subroutine mkdmftstates
