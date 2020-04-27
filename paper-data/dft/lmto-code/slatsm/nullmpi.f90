    module mpi

! MPI replacement for when getting a real MPI library is too much hassle.
! This simulates the availability of MPI environment running with 1 process.
!
! Dimitar Pashov <d.pashov@gmail.com>

    use iso_c_binding, only : c_size_t

    implicit none

    logical, private :: fmpi_on = .false., fmpi_off = .false.

    integer :: mpi_real8      = 0, &
               mpi_complex16  = 0, &
               mpi_real       = 0, &
               mpi_logical    = 0, &
               mpi_integer    = 0, &
               mpi_character  = 0, &
               mpi_max        = 0, &
               mpi_min        = 0, &
               mpi_sum        = 0, &
               mpi_prod       = 0, &
               mpi_land       = 0, &
               mpi_band       = 0, &
               mpi_lor        = 0, &
               mpi_bor        = 0, &
               mpi_lxor       = 0, &
               mpi_bxor       = 0, &
               mpi_maxloc     = 0, &
               mpi_minloc     = 0, &
               mpi_datatype_null = 0, &
               mpi_order_fortran = 0, &
               mpi_request_null  = 0, &
               mpi_max_processor_name = 0, &
               mpi_distribute_cyclic  = 0


    integer, parameter :: mpi_address_kind = c_size_t
    integer, parameter :: mpi_status_size   = 1
! These are special opaque values defined by the standard, only tentatively represented by 'integer' type in the real implementations.
! More appropriate (stongly typed) representation is provided in the mpi_f08 module

    integer :: mpi_bottom          = 0, &
               mpi_status_ignore(mpi_status_size)   = 0, &
               mpi_errcodes_ignore = 0, &
               mpi_in_place        = 0, &
               mpi_argv_null       = 0, &
               mpi_argvs_null      = 0, &
               mpi_statuses_ignore(mpi_status_size) = 0, &
               mpi_unweighted      = 0, &
               mpi_any_source      = 0

    integer, parameter :: &
        mpi_comm_world =  0, &
        mpi_comm_null  = -1, &
        mpi_info_null  = -1


    integer(8), private :: nullmpi_clock_init, nullmpi_clock_rate, nullmpi_clock_max
    real(8), private :: nullmpi_inv_clock_rate

    contains

    subroutine mpi_init(ierror)
        integer :: ierror
        ierror = 0
        fmpi_on = .true.

        call system_clock(count_rate = nullmpi_clock_rate, count_max = nullmpi_clock_max)
        call system_clock(nullmpi_clock_init)
        nullmpi_inv_clock_rate = 1.0_8/nullmpi_clock_rate
    end subroutine mpi_init

    function mpi_wtime() result(r)
        integer(8) :: ticks
        real(8) :: r
        call system_clock(ticks)
        r = real(ticks - nullmpi_clock_init, 8)*nullmpi_inv_clock_rate
    end function mpi_wtime

    function mpi_wtick() result(r)
        real(8) :: r
        r = nullmpi_inv_clock_rate
    end function mpi_wtick

    subroutine mpi_initialized(flag, ierror)
        logical :: flag
        integer :: ierror
        flag = fmpi_on
        ierror = 0
    end subroutine mpi_initialized

    subroutine mpi_finalized(flag, ierror)
        logical :: flag
        integer :: ierror
        flag = fmpi_off
        ierror = 0
    end subroutine mpi_finalized


    subroutine mpi_finalize(ierror)
        integer :: ierror
        fmpi_off = .true.
        ierror = 0
    end subroutine mpi_finalize


    subroutine mpi_barrier(comm, ierror)
        integer :: comm, ierror
        ierror = 0
    end subroutine mpi_barrier


    subroutine mpi_comm_rank(comm, rank, ierror)
        integer :: comm, rank, ierror
        rank = 0
        ierror = 0
    end subroutine mpi_comm_rank


    subroutine mpi_comm_size(comm, size, ierror)
        integer :: comm, size, ierror
        size = 1
        ierror = 0
    end subroutine mpi_comm_size


    subroutine mpi_bcast(buffer, count, datatype, root, comm, ierror)
        type(*) :: buffer(..)
        integer ::  count, datatype, root, comm, ierror
        ierror = 0
    end subroutine mpi_bcast


    subroutine  mpi_comm_split(comm, color, key, newcomm, ierror)
        integer :: comm, color, key, newcomm, ierror
        ierror = 0
        newcomm = comm
    end subroutine mpi_comm_split


    subroutine mpi_get_processor_name(name, resultlen, ierror)
        character*(*) ::  name
        integer :: resultlen, ierror
        ierror = 0
        resultlen = 0
        name = ''
    end subroutine mpi_get_processor_name

    subroutine mpi_allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
        integer ::  count, datatype, op, comm, ierror
        type(*) :: sendbuf(..), recvbuf(..)
        ierror = 0
    end subroutine mpi_allreduce

    subroutine mpi_reduce(sendbuf, recvbuf, count, datatype, op, root, comm, ierror)
        integer ::  count, datatype, op, root, comm, ierror
        type(*) :: sendbuf(..), recvbuf(..)
        ierror = 0
    end subroutine mpi_reduce

    subroutine mpi_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
        type(*) ::  sendbuf (..), recvbuf (..)
        integer ::  sendcount, sendtype, recvcount, recvtype, comm, ierror
        ierror = 0
    end subroutine mpi_allgather

    subroutine mpi_allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcount, displs, recvtype, comm, ierror)
        type(*) ::  sendbuf (..), recvbuf (..)
        integer ::  sendcount, sendtype, recvcount(*), displs(*), recvtype, comm, ierror
        ierror = 0
    end subroutine mpi_allgatherv

    subroutine mpi_comm_free(comm, ierror)
        integer ::  comm,  ierror
        ierror = 0
    end subroutine mpi_comm_free

    subroutine mpi_type_contiguous(count, oldtype, newtype, ierror)
        integer :: count, oldtype, newtype, ierror
        ierror = 0
    end subroutine mpi_type_contiguous

    subroutine mpi_type_commit(datatype, ierror)
        integer :: datatype, ierror
        ierror = 0
    end subroutine mpi_type_commit


    subroutine mpi_type_free(datatype, ierror)
        integer :: datatype, ierror
        ierror = 0
    end subroutine mpi_type_free



    subroutine mpi_cartdim_get(comm, ndims, ierror)
    !    DO NOT set ndims to anything since it could be anything positive even with 1 process.
        integer ::  comm, ndims, ierror
        ierror = 0
    end subroutine mpi_cartdim_get


    subroutine mpi_cart_get(comm, maxdims, dims, periods, coords, ierror)
        integer :: comm, maxdims, dims(*), coords(*), ierror
        logical :: periods(*)
        ierror = 0
        dims(1:maxdims) = 1
        periods(1:maxdims) = .false.
        coords(1:maxdims) = 0
    end subroutine mpi_cart_get



    subroutine mpi_dims_create(nnodes, ndims, dims, ierror)
        integer :: nnodes, ndims, dims(*), ierror
        ierror = 0
    end subroutine mpi_dims_create

    subroutine mpi_cart_create(comm_old, ndims, dims, periods, reorder, comm_cart, ierror)
        integer :: comm_old, ndims, dims(*), comm_cart, ierror
        logical :: periods(*), reorder
        ierror = 0
        comm_cart = comm_old
    end subroutine mpi_cart_create


    subroutine mpi_cart_sub(comm, remain_dims, comm_new, ierror)
        integer :: comm, comm_new, ierror
        logical :: remain_dims(*)
        ierror = 0
        comm_new = comm
    end subroutine mpi_cart_sub


    subroutine mpi_cart_rank(comm, coords, rank, ierror)
        integer :: comm, coords(*), rank, ierror
        ierror = 0
        rank = 0
    end subroutine mpi_cart_rank


    subroutine mpi_cart_coords(comm, rank, maxdims, coords, ierror)
        integer :: comm, rank, maxdims, coords(*), ierror
        ierror = 0
        coords(1:maxdims) = 0
    end subroutine mpi_cart_coords

    subroutine mpi_gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
        type(*) ::    sendbuf(..), recvbuf(..)
        integer :: sendcount, sendtype, recvcount, recvtype, root
        integer :: comm, ierror
        ierror = 0
    end subroutine mpi_gather

    subroutine mpi_gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
        type(*) :: sendbuf(..), recvbuf(..)
        integer :: sendcount, sendtype, recvcounts(*), displs(*)
        integer :: recvtype, root, comm, ierror
        ierror = 0
    end subroutine mpi_gatherv

    subroutine mpi_scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
        type(*) :: sendbuf(..), recvbuf(..)
        integer :: sendcounts(*), displs(*), sendtype
        integer :: recvcount, recvtype, root, comm, ierror
        ierror = 0
    end subroutine mpi_scatterv

    subroutine mpi_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, ierror)
        type(*) :: sendbuf(..), recvbuf(..)
        integer :: sendcounts(*), sdispls(*), sendtypes(*)
        integer :: recvcounts(*), rdispls(*), recvtypes(*)
        integer :: comm, ierror
        ierror = 0
    end subroutine mpi_alltoallw



    subroutine mpi_type_create_darray(sz, rank, ndims, array_of_gsizes, array_of_distribs, &
                                            & array_of_dargs, array_of_psizes, order,oldtype, newtype, ierror)
        integer ::  sz, rank, ndims, array_of_gsizes(*), array_of_distribs(*), &
                & array_of_dargs(*), array_of_psizes(*), order, oldtype, newtype, ierror
        ierror = 0
    end subroutine mpi_type_create_darray


    subroutine mpi_sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status, ierror)
        type(*) :: buf(..)
        integer :: count, datatype, dest, sendtag
        integer :: source, recvtag, comm
        integer :: status(mpi_status_size), ierror
        ierror = 0
    end subroutine mpi_sendrecv_replace


    subroutine mpi_irecv(buf, count, datatype, source, tag, comm, request, ierror)
        type(*) ::   buf(..)
        integer   count, datatype, source, tag, comm, request, ierror
        ierror = 0
    end subroutine


    subroutine mpi_irsend(buf, count, datatype, dest, tag, comm, request, ierror)
        type(*) ::  buf(..)
        integer   count, datatype, dest, tag, comm, request, ierror
        ierror = 0
    end subroutine mpi_irsend


    subroutine mpi_waitall(count, array_of_requests, array_of_statuses, ierror)
        integer :: count, array_of_requests(*)
        integer :: array_of_statuses(mpi_status_size,*), ierror
        ierror = 0
    end subroutine mpi_waitall


    subroutine mpi_abort(comm, errorcode, ierror)
        integer :: comm, errorcode, ierror
        ierror = 0
    end subroutine mpi_abort


    subroutine mpi_recv(buf, count, datatype, source, tag, comm, status, ierror)
        type(*) :: buf(..)
        integer   count, datatype, source, tag, comm
        integer   status(mpi_status_size), ierror
        ierror = 0
    end subroutine mpi_recv

    subroutine mpi_send(buf, count, datatype, dest, tag, comm, ierror)
        type(*) :: buf(..)
        integer   count, datatype, dest, tag, comm, ierror
        ierror = 0
    end subroutine mpi_send

    subroutine mpi_isend(buf, count, datatype, dest, tag, comm, request, ierror)
        type(*) :: buf(..)
        integer   count, datatype, dest, tag, comm, request, ierror
        ierror = 0
    end subroutine mpi_isend


    subroutine mpi_info_create(info, ierror)
        integer   info, ierror
        ierror = 0
    end subroutine mpi_info_create

    subroutine mpi_info_free(info, ierror)
        integer   info, ierror
        ierror = 0
    end subroutine mpi_info_free


    end module mpi
