
    module meshes

    implicit none

    private

    integer, parameter :: p = 8

    public :: khmesh

    contains

    subroutine khmesh(n,d, w, nl, wmax)
! Generate Kristjan Haule's "LinLogMesh"
    integer, intent(in) :: n                  ! number of points in a single level (constant step)
    real(p), intent(in) :: d                  ! periferal (largest) step
    real(p), intent(out), allocatable :: w(:) ! generated frequency mesh
    integer, intent(inout) :: nl              ! number of levels, if -1 on input then it will be caculated
    real(p), intent(in), optional :: wmax     ! the highest frequency
!     If wmax not given nl has to be > 0, if wmax given nl will be reset and on output w(2*nl*n) == wmax == -w(0)

    integer :: l, i

    if (present(wmax)) then
        nl = int(nint(log(1+wmax/(d*n))/log(2.0_p)))
    end if

! This is very very unlikely but nevertheless..
    if (nl < 1) call rx('khmesh: nl < 1, increase wmax or reduce n or d')

    if (allocated(w)) deallocate(w)
    if (.not. allocated(w)) allocate(w(0:2*nl*n))

    w(nl*n) = 0
    do l = 0, nl-1
        w(1 + (l + nl)*n : n + (l + nl)*n) =  [(w((l+nl)*n) + i*d*2**l,i=1,n)]
    end do

    do i = 1, nl*n
        w(nl*n-i) = -w(nl*n+i)
    end do

!     do i = 0, nl*n*2
!         print *, i, w(i+1)
!     end do
    end subroutine khmesh

    end module meshes