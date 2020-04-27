      subroutine rotpnt(p,q,g)
      implicit none
      double precision p(3),q(3),g(3,3),h(3)
      integer i,j
      do  i = 1, 3
        h(i) = 0d0
        do  j = 1, 3
          h(i) = h(i) + g(i,j)*p(j)
        enddo
      enddo
      do  i = 1, 3
        q(i) = h(i)
      enddo
      end
