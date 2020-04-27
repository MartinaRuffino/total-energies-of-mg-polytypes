
   module tbprl
      use mpi
      use mod_ctx
      implicit none

      private

      type tbc_t
!          integer :: wrld ! linear global communicator
         type(cart_ctx_t) :: c3d, c2d, c1d ! 3D, 2D, 1D communicators with cartesian topology
         integer, allocatable :: esamap(:) ! linear process -> atoms map for electrostatics L+1
         integer, allocatable :: escount(:) !, esdispls(:) ! sizes and displacements used in allgatherv
         integer, allocatable :: kmap(:), kcount(:) ! linear process block -> k-points map, (kcomm -> k-point block)
                                          ! the process blocks are expected to be rectangular and prefferably square
         integer, allocatable :: d2amap(:), d2count(:) !, d2displs(:) ! linear process -> atom map, counts and displacements
                                                                   ! for communication within each plane(kcomm) of c3d.
         integer :: desc(9)
         logical :: sl
         ! blacs' context, an object similar to mpi's communicator.
                                      ! row and column in blacs' grid, blocksize, descriptor, blocked ldim

         integer :: realsz(2) ! real sizes of the matrices (ldim)
         integer :: blcksz(2) ! block sizes in each dimension
         integer :: blocks(2) ! numbers of blocks
         integer :: lcblks(2) ! numbers of blocks per process in each dimension
         integer :: loclsz(2) ! local padded matrix sizes
         integer :: globsz(2) ! global padded matrix sizes


         integer, allocatable :: neighmap(:), neighidx(:)

         integer, allocatable :: nlclus(:), lclus(:,:)

      end type

      public :: tbc_t

   end module tbprl
