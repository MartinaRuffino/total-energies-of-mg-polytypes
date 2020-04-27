      subroutine gfzkin(mode,z,efre,zkin)
C- Estimate an effective kinetic energy for free electrons
Ci mode  1 use zkin = z - efre
C  this will evolve ...
      implicit none
      integer mode
      double precision z(2),efre,zkin(2)

C ... for now
      call sanrg(.true.,mode,1,1,'gfzkin:','mode')

      zkin(1) = z(1) - efre
      zkin(2) = z(2)
      end
