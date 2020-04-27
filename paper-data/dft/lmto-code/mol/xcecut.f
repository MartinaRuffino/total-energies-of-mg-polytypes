      subroutine xcecut(alat,plat,n1,n2,n3,ecut)
C- Makes an energy cutoff, equal to the largest primitive RLV
      implicit none
      integer i,j,n1,n2,n3,n123(3)
      double precision alat,plat(3,3),glat(3,3),ecut,sum(3),vol,dmax1

      n123(1) = n1
      n123(2) = n2
      n123(3) = n3
      call dinv33(plat,2,glat,vol)
      ecut = 0d0
      do  10  j = 1, 3
      do  20  i = 1, 3
   20 glat(i,j) = glat(i,j)*((n123(j)+1)/2)/alat
   10 call dpdot(glat(1,1),glat(1,1),3,sum(j))
      ecut = 1.3d0*dmax1(sum(1),sum(2),sum(3))

      end
