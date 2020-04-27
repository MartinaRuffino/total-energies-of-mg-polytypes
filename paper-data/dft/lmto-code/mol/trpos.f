      subroutine trpos(nbas,pos,i1,a)
C- Transform the pos vectors into plotting coordinates and print out
C ----------------------------------------------------------------------
Ci Inputs:
Ci
Co Outputs:
Co
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nbas,i1
      double precision pos(3,nbas),a(3,3)
C Local Variables
      integer i, j, ib
      double precision tpos(3),xpos(3)

      write (*,*) '  Transformed atom coordinates '

      do  ib = 1, nbas

C --- shift atom i1 to the origin ---
        do  i = 1, 3
          xpos(i) = pos(i,ib) - pos(i,i1)
        enddo

C --- rotate ---
        do  j = 1, 3
          tpos(j) = xpos(1)*a(1,j) + xpos(2)*a(2,j) + xpos(3)*a(3,j)
        enddo
        write (*,1) (tpos(i),i=1, 3)
    1   format (3f10.6)

      enddo

      end
