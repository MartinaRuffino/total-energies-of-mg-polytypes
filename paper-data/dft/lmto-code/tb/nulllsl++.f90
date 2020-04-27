    subroutine lsdm_cheb_c(p, beta, mu, accu, thresh, r, c, h, rho) bind(c)
! extern "C" void lsdm_cheb_c(int p, double beta, double mu, double accu, double thresh,
!                             int r, int c, double *h, double *rho /*, MPI_Comm *comm */)
        use iso_c_binding, only : c_int, c_double
        implicit none
        integer(c_int), intent(in), value :: p, r, c
        real(c_double), intent(in), value :: beta, mu, accu, thresh
        real(c_double), intent(in) :: h(*)
        real(c_double), intent(inout) :: rho(*)

        call rx('LS++ not enabled at compile time.')
    end subroutine lsdm_cheb_c
