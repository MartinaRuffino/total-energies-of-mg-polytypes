C#define OPENMP
      subroutine fmain
      implicit none
      integer n,omppid,mpipid

C#ifdef OPENMP
      n = omppid(0)
      print *, 'omppid for number of processors outside brackets:', n
C$OMP PARALLEL PRIVATE(n)
      n = omppid(0)
      print *, 'omppid for number of processors:', n
C$OMP END PARALLEL

      n = omppid(1)
      print *, 'omppid(1), outside openmp PARALLEL:', n

C$OMP PARALLEL PRIVATE(n)
      n = omppid(1)
      print *, 'omppid(1) for processor id:', n
!$OMP do
      do  n = 1, 10
      print *, 'omppid(1) in do loop: i, procid =', n,omppid(1)
      enddo
!$OMP end do

C$OMP END PARALLEL
C#endif

      n = mpipid(0)

      if (n > 1) then
      print *, 'processor number, id:', n,mpipid(1)
      call MPI_FINALIZE(n)
      endif

      end
