      logical function latvec(n,tol,qlat,vec)
C- Checks whether a set of vectors are lattice vectors
C ----------------------------------------------------------------------
Ci Inputs:
Ci   n     :number of vectors
Ci   tol   :tolerance
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   vec   :double-precision vector
Co Outputs:
Co   latvec:T if all vectors are lattice vectors within spec'd tol
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer n
      double precision qlat(3,3),vec(3,n),tol
C Local parameters:
      integer i,m
      double precision vdiff

      latvec = .false.
      do  i = 1, n
       do  m = 1, 3
         vdiff = vec(1,i)*qlat(1,m) +
     .           vec(2,i)*qlat(2,m) +
     .           vec(3,i)*qlat(3,m)
         if (dabs(vdiff-dnint(vdiff)) > tol) return
       enddo
      enddo
      latvec = .true.
      end
