         subroutine cmpvalcharg_matsub4(nbd,nspx,omg,kT,eik,eikI,V,mu,N,dN,mode)
C- Electronic density from diagonal part of G on Matsubara frequencies
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbd   :Number of states to sum over
Ci   nspx  :2 for spin-polarized case, collinear moments; otherwise 1
Ci   omg   :omega for the given site (disordered local moment)
Ci   kT    :temperature
Ci   eik   :eigenvalues of hamiltonian at this omega
Ci   eikI  :eigenvalues of updated hamiltonian at omega = 0
Ci   v     :scalar potential
Ci   mu    :chemical potential
Ci   mode  :one of:
Ci         :'sumfrq'
Ci         :'fermif'
Co Outputs
Co   N     :Electron number at chemical potential
Co   dN    :DOS at chemical potential
Cs Command-line switches
Cr Remarks
Cr This subroutine computes expression (36) of Haule's PRB 81, 195107 (2010).
Cr So it computes the electronic density N of a diagonal Green's function
Cr  written on a grid of imaginary frequencies.
Cr
Cr Moreover it computes the derivative of this formula with respect to V
Cr  (dropping the 1/eik0 term and correspondingly the fermi function).
Cr This is the analytical derivative dN_an and the negative of the DOS
Cr
Cr Finally it takes the density of states near mu as Im[G(omega0)]
Cr  where omega0 is the closest energy to the real axis, that is o(1).
Cr This is the sampled derivative dN_s.
Cr
Cr N and derivative of N (or N and DOS near mu) are used by the rfalsi
Cr  algorithm
Cr
Cr NOTE: We are supposed to enter this subroutine ONLY if
Cr       frequencies are imaginary. That is if omf==cmplx(0.0,1.0)
Cu Updates
Cu   17 Jan 18 (MvS) First cut at SO coupled case
Cu        2016 Written by Lorenzo Sponza
C ----------------------------------------------------------------------
      implicit none
      real(8),   parameter    :: pi=acos(-1d0)
      complex(8),parameter    :: ii=cmplx(0d0,1d0),cZERO=cmplx(0d0,0d0),cONE=cmplx(1d0,0d0)
      integer,   intent(in)   :: nbd,nspx       ! number bands and spins
      real(8),   intent(in)   :: omg            ! current freq. if sumfrq, last freq if fermif
      real(8),   intent(in)   :: kT             ! frequencies (on Im(omega)>0.0 plane)
      complex(8),intent(in)   :: eik(nbd,nspx)  ! eigenvalues of updated hamiltonian (at given omg)
      complex(8),intent(in)   :: eikI(nbd,nspx) ! eigenvalues of updated hamiltonian (at omg==o(no))
      real(8),   intent(in)   :: V              ! scalar potential (might be 0.0)
      real(8),   intent(in)   :: mu             ! chemical potential
      character(len=6), intent(in) :: mode      ! either sum over frequencies or (arctans+fermi)
      real(8),   intent(inout):: N(nspx)        ! electron number at given V
      real(8),   intent(inout):: dN(nspx)       ! derivative of N wrt V
c     local variables
      integer  :: ib,is
      real(8)  :: fermifn,derfermi  ! arctan and fermi terms
      real(8) :: nloc, dnloc
      complex(8) :: t1, t2
C     double precision xx

      call tcn('cmpvalcharg_matsub4')

      select case(mode)
C   ... Called within the frequency loop.
        case('sumfrq')
          do is = 1, nspx
            nloc = 0; dnloc = 0
            do ib = 1, nbd
             t1 = 1/(ii*omg+mu-eik(ib,is)-v); t2 = 1/(ii*omg+mu-real(eiki(ib,is))-v)
             nloc = nloc + real(t1 - t2)
             dnloc = dnloc + real(t1*t1 - t2*t2)
            enddo
            n(is) = n(is) + 2*kt*nloc
            dn(is) = dn(is) + 2*kt*dnloc
          enddo


C   ... Called outside the frequency loop.
        case('fermif')
          do  is = 1, nspx
            do  ib = 1, nbd
              n(is)  = n(is)  + fermifn(dble(eiki(ib,is))-mu+v,kt)
              dn(is) = dn(is) + derfermi(dble(eiki(ib,is))-mu+v,kt)
            enddo
          enddo

        case default
          call rx('Wrong mode in cmpvalcharg_matsub4')

      end select

      call tcx('cmpvalcharg_matsub4')

      end subroutine cmpvalcharg_matsub4

      function fermifn(x,kt)
      implicit none
      real(8) :: fermifn, x,kt
      real(8) :: arg
      arg=x/kt
      if(arg<-100d0) then
        fermifn=1d0
      else if(arg>100d0) then
        fermifn=0d0
      else
        fermifn=1d0/(1d0 + exp(arg))
      endif
      end function fermifn

      function derfermi(x,kt)
      implicit none
      real(8) :: derfermi, x,kt
      real(8) :: arg
      arg=x/kt
      if(abs(arg)>100d0)  then
        derfermi=0
      else
        derfermi=-exp(arg)/(kt*(1d0+exp(arg))**2)
      endif
      end function derfermi
