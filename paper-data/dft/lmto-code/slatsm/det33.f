      real(8) function det33(am)
C- Calculates the determinant of a 3X3 matrix
      implicit none
      real(8),intent(in) :: am(3,3)
      det33= am(1,1)*am(2,2)*am(3,3)
     .       -am(1,1)*am(3,2)*am(2,3)
     .       -am(2,1)*am(1,2)*am(3,3)
     .       +am(2,1)*am(3,2)*am(1,3)
     .       +am(3,1)*am(1,2)*am(2,3)
     .       -am(3,1)*am(2,2)*am(1,3)
      end

