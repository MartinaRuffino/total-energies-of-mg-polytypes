      subroutine makehbar2(nev,nev0,nspx,ndhamx,evlk,sigbark,updthamk)
C- Compute the updated hamiltonian k-resolved for each energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   nev     :number of states in the hybridization window
Ci   nev0    :offset to first state in hybridization window
Ci   nspx    :dimensioning parameter,2 for collinear, spin-polarized case, otherwise 1
Ci   ndhamx  :dimensioning parameter, largest hamiltonian dimension
Ci   evlk    :eigenvalues (from the QSGW/LDA Hamiltonian)
Ci           :evlk(nev0+1,isp) corresponds to sigbark(1,1,isp)
Ci   sigbark :embedded (sig.inp-sigDC) in lattice representation (ijk)
Ci           :sigbark(1:nev,1:nev,1:nspx) is restricted to the nev
Ci           :elements within the hybridization window.
Co Outputs
Co  updthamk : updated hamiltonian in lattice representation (ijk)
Cr Remarks
Cr  The hamiltonian is constructed for each k-point and each energy,
Cr   so it is possible to use frequency-dependent eigenvalues.
Cu Updates
Cu   17 Jan 18 (MvS) First cut at SO coupled case
Cu   22 Aug 15 (L. Sponza) first created
Cu   09 Nov 15 new version with sigbark passed as input
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer,    intent(in)   :: nev,nev0,nspx,ndhamx
      real(8),    intent(in)   :: evlk(ndhamx,nspx)
      complex(8), intent(in)   :: sigbark(nev,nev,nspx)
      complex(8), intent(inout):: updthamk(ndhamx,ndhamx,nspx)
C ... Local parameters
      integer :: isp,i,j

C     call tcn('makehbar2')

C     Construct the hamiltonian
C     It involves ALL eigenvalues, not only those in [lo,hi] window
      do  isp = 1, nspx
        forall (i = 1:ndhamx) updthamk(i,i,isp) = evlk(i,isp)
        do  i = 1, nev
          do  j = 1, nev
            updthamk(nev0+i,nev0+j,isp) = updthamk(nev0+i,nev0+j,isp) + sigbark(i,j,isp)
         enddo
        enddo
      enddo

C     call tcx('makehbar2')
      end subroutine makehbar2
