      subroutine dsumdf(n,scal,a1,ofa1,l1,a2,ofa2,l2)
C- Returns scaled sum and difference of two vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   n    :number elements to scale and combine
Ci   scal :scale sum and difference by scal; see Outputs
Ci   a1   :first vector
Ci   ofa1 :offset to first entry in a1
Ci   l1   :skip length in a1
Ci   a2   :second vector
Ci   ofa2 :offset to first entry in a2
Ci   l2   :skip length in a2
Co Outputs
Co   a1   :a1 <- scal*(a1+a2)
Co   a2   :a2 <- scal*(a1-a2)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer n,l1,l2,ofa1,ofa2
      double precision scal, a1(1), a2(1)
C Local parameters
      real(8), allocatable :: a(:)

C --- a1-a2-> a  a1+a2 -> a1;  a -> a2 ---
      allocate(a(n))
      call dcopy(n,a1(1+ofa1),l1,a,1)  ! a <- a1
      call daxpy(n,-1d0,a2(1+ofa2),l2,a,1) ! a <- a1-a2
      call daxpy(n,1d0,a2(1+ofa2),l2,a1(1+ofa1),l1) ! a1 <- a1+a2
      call dcopy(n,a,1,a2(1+ofa2),l2) ! a2 <- a1-a2
      deallocate(a)

      if (scal == 1) return
      call dscal(n,scal,a1(1+ofa1),l1)
      call dscal(n,scal,a2(1+ofa2),l1)

      end
