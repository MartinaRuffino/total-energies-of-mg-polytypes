      subroutine setne(s_bz,s_ctrl,nzne,vne,lnoneq)
C- Set up parameters for non-equilibrium mode
C---------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  semsh
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lpgf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Co Outputs
Co nzne      :number of energy points along non-equilibrium contour
Co vne       :difference in fermi energies of right and left leads
Co           :vne=ef(R)-ef(L); may be <0
Co lnoneq    :=T for non-equilibrium mode (semsh(2)/100=1 and lpgf=1,5)
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   10 Feb (S.Faleev) First created
C---------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nzne
      double precision vne
      logical lnoneq
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
C ... Local parameters
      integer lpgf
      double precision semsh(10)
      real(8),parameter :: tol = 1d-8 ! Potential less than tol is treated as equilibrium
      logical, parameter :: T=.true., F=.false.

      semsh = s_bz%semsh
      lpgf = s_ctrl%lpgf(1)

      nzne = 0
      vne = 0d0
      lnoneq = F
      if ( mod(nint(semsh(2))/100,10) == 1
     .     .and.  (lpgf == 1 .or. lpgf == 5) ) then
         if (lpgf == 1) nzne = nint(semsh(7))
         if (nzne < 0) call rx('setne: illegal nzne < 0')
         vne  = semsh(8)
C      abs(vne) < tol considered to be an equilibrium case
         if (abs(vne) < tol .and. nzne /= 0)
     .      call rx('setne: abs(vne)<tol and nzne!=0')
         if (abs(vne) >= tol ) lnoneq = T
      endif

      end
