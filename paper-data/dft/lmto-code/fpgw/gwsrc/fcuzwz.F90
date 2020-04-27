
    module fcugw

    implicit none

#ifdef CUD
    interface

    subroutine init_cublas() bind(c)
    end subroutine init_cublas

    subroutine stop_cublas() bind(c)
    end subroutine stop_cublas

    subroutine fcuzwz(m, n, k, w, ldw, z, ldz, s, lds, same_z) bind(c)
        use iso_c_binding, only : c_int, c_double
        integer(c_int), intent(in), value :: m, n, k, ldw, ldz, lds, same_z
        complex(c_double), intent(in), dimension(*) :: w, z
        complex(c_double), intent(out), dimension(*) :: s
    end subroutine fcuzwz

    end interface
#endif

    end module fcugw



