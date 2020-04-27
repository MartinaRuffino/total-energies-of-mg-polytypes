      subroutine cross(a,b,c)
C Example:
C dval a1=0 a2=1 a3=1 b1=1 b2=0 b3=1 c1=1 c2=1 c3=0 a2*b3-a3*b2 a3*b1-a1*b3 a1*b2-a2*b1
      implicit none
      double precision a(3),b(3),c(3)
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
      end
      double precision function tripl(a,b,c)
      implicit none
      double precision a(3),b(3),c(3)
      tripl = a(1)*b(2)*c(3) + a(2)*b(3)*c(1) + a(3)*b(1)*c(2)
     .       -a(3)*b(2)*c(1) - a(2)*b(1)*c(3) - a(1)*b(3)*c(2)
      end
      subroutine mkqlat(plat,qlat,vol0)
C- Reciprocal of a lattice vector
      implicit none
      double precision plat(3,3),qlat(3,3),tripl,vol0
      integer i,k
      call cross(plat(1,2),plat(1,3),qlat)
      call cross(plat(1,3),plat(1,1),qlat(1,2))
      call cross(plat(1,1),plat(1,2),qlat(1,3))
      vol0 = tripl(plat,plat(1,2),plat(1,3))
      do 3 i=1,3
      do 3 k=1,3
    3 qlat(i,k)=qlat(i,k)/vol0
      end
