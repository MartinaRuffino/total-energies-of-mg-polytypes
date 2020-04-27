      subroutine grpop(v,v1,g,i)
      implicit none
      double precision g(3,3,*),v(3),v1(3)
      integer i
      v1(1) = g(1,1,i)*v(1) + g(1,2,i)*v(2) + g(1,3,i)*v(3)
      v1(2) = g(2,1,i)*v(1) + g(2,2,i)*v(2) + g(2,3,i)*v(3)
      v1(3) = g(3,1,i)*v(1) + g(3,2,i)*v(2) + g(3,3,i)*v(3)
      end

      subroutine zgrpop(v,v1,g,i)
      implicit none
      double precision g(3,3,*)
      double complex v(3),v1(3)
      integer i
      v1(1) = g(1,1,i)*v(1) + g(1,2,i)*v(2) + g(1,3,i)*v(3)
      v1(2) = g(2,1,i)*v(1) + g(2,2,i)*v(2) + g(2,3,i)*v(3)
      v1(3) = g(3,1,i)*v(1) + g(3,2,i)*v(2) + g(3,3,i)*v(3)
      end
