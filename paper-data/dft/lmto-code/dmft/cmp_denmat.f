         subroutine cmp_denmat(nbd,omg,kT,eik,cljk,crik,eik0,crik0,cljk0,
     .                         mu,Wij,mode)
C- Computes expression W^{DMFT}_{ijk} of Haule's PRB 81, 195107 (2010) and diagonalises it
C ----------------------------------------------------------------------
Ci   nbd   : number bands and spins
Ci   omg   : Matsubara frequency, imaginary part (omg is Im(omega)>0 plane
Ci   kT    : First Matsubara frequency
Ci   eik   : eigenvalues of DMFT hamiltonian at this omg
Ci   cljk  : left  eigenvector of DMFT hamiltonian at omg
Ci   crik  : right eigenvector of DMFT hamiltonian at omg
Ci   eik0  : eigenvalues of updated hamiltonian
Ci   crik0 : left  eigenvectors of DMFT hamiltonian at the last Matsubara frequency
Ci   cljk0 : right eigenvectors of DMFT hamiltonian at the last Matsubara frequency
Ci   mu    : chemical potential
Ci   mode  : 'sumfrq' contribution from this frequency
Ci         : 'fermif' contribution from endpoint
Co Outputs
Co   Wij   : DMFT density matrix element.  See p9 in PRB 81, 195107 (2010)
Cr Remarks
Cr  We are supposed to enter this subroutine ONLY if
Cr  frequencies are imaginary. That is if omf==cmplx(0.0,1.0)
Cr
Cu Updates
Cu   Written by Lorenzo Sponza
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character(len=6), intent(in) :: mode
      integer,   intent(in)   :: nbd ! number bands and spins
      real(8),   intent(in)   :: omg ! current freq. if sumfrq, last freq if fermif
      real(8),   intent(in)   :: kT ! frequencies (on Im(omega)>0.0 plane)
      complex(8),intent(in)   :: eik(nbd),eik0(nbd) ! eigenvalues of updated hamiltonian
      complex(8),intent(in)   :: cljk(nbd),cljk0(nbd) ! left  eigenvectors of updated hamiltonian (at omg)
      complex(8),intent(in)   :: crik(nbd),crik0(nbd) ! right eigenvectors of updated hamiltonian (at omg)
      real(8),   intent(in)   :: mu ! chemical potential
      complex(8),intent(inout):: Wij ! DMFT density matrix element
C ... Dynamically allocated local arrays
C ... Local parameters
      real(8),   parameter    :: pi=acos(-1d0)
      complex(8),parameter    :: ii=cmplx(0d0,1d0)
      integer  :: ib
      procedure(real(8)) :: fermifn

      select case(mode)
       case('sumfrq')

        do ib=1,nbd
          Wij = Wij + 2*kT*dble(crik(ib) *conjg(cljk(ib)) /(ii*omg+mu-eik(ib))
     .                  - dble(crik0(ib)*conjg(cljk0(ib)))/(ii*omg+mu-dble(eik0(ib))))
        enddo

       case('fermif')
        do ib=1,nbd
          Wij = Wij + dble(crik0(ib)*conjg(cljk0(ib)))*fermifn(dble(eik0(ib))-mu,kT)
        enddo

       case default
         call rx('Wrong mode in cmpdenmat')

      end select

      end subroutine cmp_denmat
