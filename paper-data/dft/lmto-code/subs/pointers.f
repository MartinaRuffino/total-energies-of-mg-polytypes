      module pointers
      use structures, only : str_symops, str_rhat
      implicit none
      private

!     Generate permanent pointers declared as public variables. This preserves the pointers.
!     Used by routines that need to allocate a permanent pointer
!     Example:
!     allocate(p_d1(n)); s_str%alph => p_d1
!---
!     This 'feature' likely amounts to a bug in the compiler, thus it only works by accident.
!     The sane way to perform the action is to declare the pointer (or allocatable object) argument properly and simply allocate
!     DMT

!     Generic 1d integer pointer
      integer, pointer :: p_i1(:)
!     Generic 2d integer pointer
      integer, pointer :: p_i2(:,:)
!     Generic 1d double precision pointer
      real(8), pointer :: p_d1(:)
!     Generic 2d double precision pointer
      real(8), pointer :: p_d2(:,:)
!     Generic 1d double complex pointer
      complex(8), pointer :: p_z1(:)
!     Generic 2d double complex pointer
      complex(8), pointer :: p_z2(:,:)
!     Generic symmetry structure pointer
      type(str_symops),pointer :: p_sym(:)
      type(str_rhat),pointer ::  p_rhat(:)

      public :: p_i1,p_i2,p_d1,p_d2,p_z1,p_z2,p_sym,p_rhat
      end module pointers
