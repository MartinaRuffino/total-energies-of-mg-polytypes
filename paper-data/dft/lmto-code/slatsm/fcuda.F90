   module fcuda

   implicit none

   private


   interface

   function cublasInit() bind(c, name='cublasInit') result(r)
      use iso_c_binding, only : c_int
      integer(c_int) :: r
   end function cublasInit

   function cublasShutdown() bind(c, name='cublasShutdown') result(r)
      use iso_c_binding, only : c_int
      integer(c_int) :: r
   end function cublasShutdown

!#ifdef! CUDAD
!   subroutine dsyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info) bind(c,name='magma_dsyevd')
!      use iso_c_binding, only : c_char, c_int, c_double
!      character(kind=c_char), intent(in), value :: jobz, uplo
!      integer(c_int), intent(in), value         :: n, lda, lwork, liwork
!! if f77 w joy
!      integer                                   :: a(*), work(*)
!      real(c_double)                            :: w(*)
!! else
!!             real(c_double), intent(inout)             :: a(*)
!!             real(c_double), intent(out)               :: w(*), work(*)
!! end joy
!      integer(c_int), intent(out)               :: iwork(*), info
!   end subroutine dsyevd
!
!   subroutine dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info) bind(c,name='magma_dsygvd')
!      use iso_c_binding, only : c_char, c_int, c_double
!      character(c_char), intent(in), value :: jobz, uplo
!      integer(c_int), intent(in), value :: itype, n, lda, ldb, lwork, liwork
!! if f77 w joy
!      integer                                   :: a(*), b(*), work(*)
!      real(c_double)                            :: w(*)
!! else
!!             real(c_double), intent(inout)             :: a(*)
!!             real(c_double), intent(out)               :: w(*), work(*)
!! end joy
!      integer(c_int), intent(out)               :: iwork(*), info
!   end subroutine dsygvd
!#endif

   subroutine fcuda_get_device_count(n) bind(c, name='cudaGetDeviceCount')
      use iso_c_binding, only : c_int
      integer(c_int), intent(out) :: n
   end subroutine fcuda_get_device_count

   subroutine fcuda_set_device(i) bind(c, name='cudaSetDevice')
      use iso_c_binding, only : c_int
      integer(c_int), intent(in), value :: i
   end subroutine fcuda_set_device

   subroutine fcuda_get_device(i) bind(c, name='cudaGetDevice')
      use iso_c_binding, only : c_int
      integer(c_int), intent(out) :: i
   end subroutine fcuda_get_device

   subroutine fcuda_get_local_dev_props(devid, name, major, minor,  &
   & totalGlobalMem, multiProcessorCount, computeMode) bind(c, name='fcuda_get_local_dev_props')
      use iso_c_binding, only : c_int, c_char
      implicit none
      integer, intent(out) :: devid, major, minor, totalGlobalMem, multiProcessorCount, computeMode
      character(c_char), intent(out) :: name
   end subroutine fcuda_get_local_dev_props

!    function device_is_usable(devid) bind(c, name='device_is_usable') result(r)
!       use iso_c_binding, only : c_int
!       integer(c_int), intent(in), value :: devid
!       integer(c_int) :: r
!    end function

   end interface

   public :: fcublasInit, fcublasShutdown
   public :: fcuda_get_device_count, fcuda_set_device, fcuda_get_device, &
            & fcuda_get_local_dev_props !, device_is_usable


!#ifdef! CUDAD
!   public :: dsyevd, dsygvd
!#endif

   contains

   function fcublasInit() result(ierr)
      integer :: ierr
      ierr = 0
!#ifdef! CUDA
!      ierr = cublasInit()
!#endif
   end function fcublasInit

   function fcublasShutdown() result(ierr)
      integer :: ierr
      ierr = 0
!#ifdef! CUDA
!      ierr = cublasShutdown()
!#endif
   end function fcublasShutdown



   end module fcuda