  subroutine cexit(pv, ps)
!     use iso_c_binding
    use mpi
    use fcuda, only: fcublasShutdown
    use posix

    implicit none

    integer :: err, ip, sz
    logical :: flag
    integer, intent(in) :: pv, ps

    integer, pointer :: nullptr => null()

    if (ps /= 0 ) then

#ifdef TRACE
!--- Print backtrace by writing/dereferencing a nullpointer
      if (pv /= 0) then
         call mpi_comm_rank(mpi_comm_world, ip, err)
         call mpi_comm_size(mpi_comm_world, sz, err)
         write(0, "('Induced backtrace from process ',i0,'/',i0,':')") ip,sz
         flush(0)
         nullptr = 343
         print *, nullptr
! FIXME: This is very crude and should be replaced by properly generated backtrace
      end if
#endif

      err = fcublasShutdown()
      call blacs_exit(1)
      if (pv == -1) then
        call mpi_comm_size(mpi_comm_world, sz, err)
        if (sz > 1) call mpi_abort(mpi_comm_world,pv,err)
      endif
      call mpi_finalized(flag,err)
      if (.not. flag) call mpi_finalize(err)
      call exit(pv)
    end if
  end subroutine cexit

  program cmain
    use mpi
    use fcuda, only : fcublasInit

    implicit none

    integer :: err

!     call loc_dev_init_insane();
    call mpi_init(err)
#ifdef TESTCOV
    call testcov_init(TESTCOV)
#endif
    call loc_dev_init(err)
    err = fcublasInit()
    call fmain()
    call cexit(0,1)
  end program cmain


