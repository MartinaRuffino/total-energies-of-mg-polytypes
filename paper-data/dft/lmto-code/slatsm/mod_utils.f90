    module mod_utils
!     use mod_const
    implicit none

    private
! to be removed when transitioning to mod_const
    integer, parameter :: dp = 8

    public :: xor, log2i, mergesort

    contains

    logical function xor(p, q)
        logical, intent(in) :: p, q
        xor = (p .or. q) .and. .not. (p .and. q)
    end function xor

    integer function log2i(i) result(p)
        integer, intent(in) :: i
!         endian independent (using 2**p instead of ishift)
        p = 0
        do while (2**p < i)
            p = p+1
        end do
    end function log2i

    subroutine mergesort(a, k, n)
! stable key sort in n*log(n) time with memory usage of n. https://en.wikipedia.org/wiki/Merge_sort

! length of input
    integer, intent(in) :: n
! input to be key sorted
    real(dp), intent(in) :: a(0:n-1)
! the sorting key sought, indexing from 0
    integer, intent(out), target :: k(0:n-1)

! local vars
    integer :: w, i, p, pmax
    integer, target :: t(0:n-1)
    integer, pointer :: c(:), d(:)

    pmax = log2i(n)-1

    if (mod(pmax,2) == 0) then
        forall (i=0:n-1) t(i) = i
    else
        forall (i=0:n-1) k(i) = i
    end if

    do p = 0, pmax
        if (xor(mod(p,2) == 0, mod(pmax,2) == 0)) then
            c => k
            d => t
        else
            c => t
            d => k
        end if

        w = 2**p
        do i = 0, n-1, 2*w
            call mergepieces(a, i, min(i+w, n), min(i+2*w, n), c, d)
        end do
    end do

    contains
        subroutine mergepieces(v, ileft, iright, iend, a, b)

        implicit none
        integer, intent(in) :: ileft, iright, iend
        integer, intent(inout) :: a(0:iend-1), b(0:iend-1)
        real(dp), intent(in) :: v(0:iend-1)

        integer :: i, j, k

        i = ileft
        j = iright
        do k = ileft, iend-1
            if (i >= iright .or. (j < iend .and. v(a(i)) > v(a(j)))) then
                b(k) = a(j)
                j = j + 1
            else
                b(k) = a(i)
                i = i + 1
            end if
        end do
        end subroutine mergepieces
    end subroutine mergesort

    end module mod_utils

    subroutine reorder(k, n, a, lda, m)
    implicit none
! Order the n column vectors of length m, held in array a wich leading dimension lda, according to the key k.

! double precision size, to be removed when transitioning to mod_const
    integer, parameter :: dp = 8

! numner of vectors, length of k
    integer, intent(in) :: n
! k, key to guide the ordering
    integer, intent(in) :: k(0:n-1)
! leading dimension of a
    integer, intent(in) :: lda
! vector length (m <= lda)
    integer, intent(in) :: m
! array containing the vectors
    real(dp), intent(inout) :: a(lda,0:n-1)
!     type(*), intent(inout) :: a(0:*)

    integer :: i, i0
    logical :: ordered(0:n-1)
    real(dp) :: t(m)

    forall (i=0:n-1) ordered(i) = i == k(i)

    do i0 = 0, n-1
        if (ordered(i0)) cycle
        call dcopy(m,a(1,i0),1,t,1)
        i = i0
        do while (i0 /= k(i))
            call dcopy(m,a(1,k(i)),1,a(1,i),1)
            ordered(i) = .true.
            i = k(i)
        end do
        call dcopy(m,t,1,a(1,i),1)
        ordered(i) = .true.
    end do

    end subroutine reorder