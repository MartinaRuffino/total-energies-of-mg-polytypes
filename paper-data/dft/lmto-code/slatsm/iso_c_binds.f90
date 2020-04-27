

  module iso_c_binds

! Comment for compatibility with old open64 (older than ~2010?)
   use, intrinsic :: iso_c_binding

! On 32bit platform "long" is 4 bytes!
!    integer, parameter :: c_long = 8
!
!
!    integer, parameter :: c_char = 1
!    integer, parameter :: c_short = 2
!    integer, parameter :: c_int = 4
!    integer, parameter :: c_long_long = 8
!    integer, parameter :: c_float = 4
!    integer, parameter :: c_double = 8
!    integer, parameter :: c_long_double = 16
!
!    character, parameter :: c_new_line = achar(10)
!    integer, parameter :: c_null_ptr = 0
!    character, parameter :: c_null_char = achar(0)

  end module iso_c_binds
