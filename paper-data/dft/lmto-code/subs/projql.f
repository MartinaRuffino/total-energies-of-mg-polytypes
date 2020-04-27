      subroutine projql(qlat)
C- Project 3D reciprocal lattice vecs on surface (2D recip lattice vecs)
C ----------------------------------------------------------------------
Ci Input qlat(3,3)  3D reciprocal lattice vectors
Co Output qlat(3,3) 2D reciprocal lattice vectors
Cr if bi i=1,3 are 3D vectors
Cr    Bi i=1,2 are 2D vectors
Cr    Bi(u)=bi(u)-(Bi.B3) B3(u)/(B3.B3) for u=x,y,z
C ----------------------------------------------------------------------
      implicit none
      double precision qlat(3,3),b3b3,bib3
      integer i,j,stdo,iprint
      procedure(integer) :: nglob

      b3b3 = qlat(1,3)**2 + qlat(2,3)**2 + qlat(3,3)**2
      do  i = 1, 2
        bib3 = qlat(1,i)*qlat(1,3) + qlat(2,i)*qlat(2,3)
     .       + qlat(3,i)*qlat(3,3)
        do  j = 1, 3
          qlat(j,i) = qlat(j,i) - qlat(j,3)*bib3/b3b3
        enddo
      enddo
C ... Ignore qlat(3,.) by making it zero
C      do  22  j = 1, 3
C   22 qlat(j,3)=0d0
      if (iprint() >= 20) then
      stdo = nglob('stdo')
      write(stdo,*) ' 2D-reciprocal vectors:'
        do  i = 1, 2
          write (stdo,1) (qlat(j,i),j=1,3)
    1     format(3F12.6)
        enddo
      endif
      end
