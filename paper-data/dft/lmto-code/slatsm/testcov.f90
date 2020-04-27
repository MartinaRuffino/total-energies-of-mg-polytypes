      subroutine testcov_init(strip)
! Send gcda files from different processes to different paths to avoid a race condition.
      use mpi
      implicit none

      integer, intent(in) :: strip
      integer :: procid, err
      character(len=24) :: var

      call mpi_comm_rank(mpi_comm_world, procid, err)

      write(var, '("GCOV_PREFIX=cov/",i8.8)') procid
      call ptenv(trim(var))

      if (strip > 0) then
        write(var, '("GCOV_PREFIX_STRIP=",i0)') strip
        call ptenv(trim(var))
      end if

      end subroutine testcov_init
