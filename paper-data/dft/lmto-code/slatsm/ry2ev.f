      real(8) function ry2eV(opt)
C- Convert Rydberg to eV
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt : 0 => return 1 Ry, in eV units
Ci       : 1 => return 1 Ha, in eV units
Ci       : 10 => return 1 eV, in Ry units
Ci       : 11 => return 1 eV, in Ha units
Co Outputs
Co   ry2eV
Cu Updates
Cu   06 Aug 18 First created: use NIST standard for ratio eV/Ry
Cu             https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev|search_for=rydberd
C ----------------------------------------------------------------------
      implicit none
      integer opt
      real(8), parameter :: eVbyRy = 13.605 693 009 d0

      select case (opt)
        case (1);     ry2eV = eVbyRy*2
        case (10);    ry2eV = 1/eVbyRy
        case (11);    ry2eV = 1/(2*eVbyRy)
        case default; ry2eV = eVbyRy
      end select

      end
C      subroutine fmain
C      procedure(real(8)) :: ry2eV
C
C      print *, ry2eV(0)
C      print *, ry2eV(1)
C      print *, ry2eV(10)
C      print *, ry2eV(11)
C      end
