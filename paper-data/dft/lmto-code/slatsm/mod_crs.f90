

   module mod_crs
! Module for handling column resolved sparse matrices.
! Undocumented and so far only used to test the ls++ library

      implicit none

      private

      integer, parameter :: dp = 8

      type crs_t
         integer :: d(0:1)
         integer, allocatable :: cidx(:), rptr(:)
         real(dp), allocatable :: vals(:)
      end type crs_t

      public :: crs_t, tdense2crs_naive, print_crs, del_crs, tdense_print_as_crs, &
                                read_crs, read_crs_as_dense, tcrs2dense_naive

   contains


   subroutine tdense2crs_naive(a,lda,n,m, b)

      integer, intent(in) :: lda, n, m
      real(dp), intent(in) :: a(0:lda-1,0:m-1)
      type(crs_t), intent(out) :: b

      integer :: i,j,k
      logical :: first_found

      b % d = [m,n] ! transposed

      allocate(b % rptr(0:m))

      k = 0
      do j = 0, m - 1
         first_found = .false.
         do i = 0, n - 1
            if (a(i,j) == 0.0_dp) cycle
            if (.not. first_found) then
               b % rptr(j) = k
               first_found = .true.
            end if
            k = k + 1
         end do
      end do
      b % rptr(m) = k

      allocate(b % cidx(0:k-1))
      allocate(b % vals(0:k-1))

      k = 0
      do j = 0, m - 1
         do i = 0, n - 1
            if (k == b % rptr(j+1)) exit
            if (a(i,j) == 0.0_dp) cycle
            b % cidx(k) = i
            b % vals(k) = a(i,j)
            k = k + 1
         end do
      end do

   end subroutine tdense2crs_naive


   subroutine print_crs(fln, a)
      character(len=*), intent(in) :: fln
      type(crs_t), intent(in) :: a

      integer :: i, u

      open(unit=u, file=fln, action='write')
!     open(newunit=u, file=fln, action='write')

      write(u, '(" rows ", i0, /, " cols ", i0)') a % d

      write(u, '(" rowpointers ", i0)') a % d (0) + 1
      write(u, '(x, 16(x,i0))') a % rptr(0 : a % d (0))

      write(u, '(/, " indices ", i0)') a % rptr(a % d (0))
      write(u, '(x, 16(x,i0))') a % cidx(0 : a % rptr(a % d (0)) - 1)

      write(u, '(/, " values ", i0)') a % rptr(a % d (0))
      write(u, '(x, 16(x,es20.12))', advance='no') a % vals(0 : a % rptr(a % d (0)) - 1)

      write(u,'(/)')

      close(u)

   end subroutine print_crs



   subroutine tdense_print_as_crs(a,lda,n,m, fln)

      integer, intent(in) :: lda, n, m
      real(dp), intent(in) :: a(0:lda-1,0:m-1)
      character(len=*), intent(in) :: fln

      type(crs_t) :: b

      call tdense2crs_naive(a,lda,n,m, b)
      call print_crs(fln, b)
      call del_crs(b)
   end subroutine tdense_print_as_crs


   subroutine read_crs(fln, a)
      character(len=*), intent(in) :: fln
      type(crs_t), intent(out) :: a

      character(len=16) :: c
      integer :: tmpi
      integer :: i, u

      call del_crs(a)

      open(unit=u, file=fln, action='read')
!     open(newunit=u, file=fln, action='read')

      read(u, *) c(1: 4), a % d(0); if (c(1:4) /= 'rows') call rx('READ_CRS: trouble reading the rows line')
      read(u, *) c(1: 4), a % d(1); if (c(1:4) /= 'cols') call rx('READ_CRS: trouble reading the cols line')

      read(u, *) c(1:11), i
      if (c(1:11) /= 'rowpointers' .or. i /= a%d(0)+1) call rx('READ_CRS: trouble reading the rowpointers line')
      allocate(a % rptr(0 : a % d (0)))
      read(u, *) a % rptr(0 : a % d (0))

      read(u, *) c(1: 7), i
      if (c(1:7) /= 'indices' .or. i /= a % rptr(a % d (0))) call rx('READ_CRS: trouble reading the indices line')
      allocate(a % cidx(0 : a % rptr(a % d (0)) - 1))
      read(u, *) a % cidx(0 : a % rptr(a % d (0)) - 1)

      read(u, *) c(1: 6), i
      if (c(1: 6) /= 'values' .or. i /= a % rptr(a % d (0))) call rx('READ_CRS: trouble reading the values line')
      allocate(a % vals(0 : a % rptr(a % d (0)) - 1))
      read(u, *) a % vals(0 : a % rptr(a % d (0)) - 1)

      close(u)

   end subroutine read_crs


   subroutine tcrs2dense_naive(a, b)
       type(crs_t), intent(in) :: a
       real(8), allocatable, intent(out) :: b(:,:)

       integer :: i, j, k

       allocate(b(0: a % d(1)-1,0: a % d(0)-1))

       b = 0.0_8

       do j = 0, a % d(0)-1
            do k = a % rptr(j), a % rptr(j+1)-1
                i = a % cidx(k)
                b(i,j) = a % vals(k)
            end do
       end do

   end subroutine tcrs2dense_naive

   subroutine read_crs_as_dense(fln, a)
       character(len=*), intent(in) :: fln
       real(8), allocatable, intent(out) :: a(:,:)


       type(crs_t) :: b

       call read_crs(fln, b)
       call tcrs2dense_naive(b, a)
       call del_crs(b)

   end subroutine read_crs_as_dense


   subroutine del_crs(a)
      type(crs_t), intent(inout) :: a

      a % d = 0
      if (allocated(a%rptr)) deallocate(a%rptr)
      if (allocated(a%cidx)) deallocate(a%cidx)
      if (allocated(a%vals)) deallocate(a%vals)

   end subroutine del_crs

   end module mod_crs
