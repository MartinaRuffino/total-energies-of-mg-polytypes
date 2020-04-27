      integer function omppid(mode)
C- Returns OPENMP thread id
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return number of threads
Ci         :1 return thread id
Ci         :<2 call OMP_SET_NUM_THREADS(-mode)
Ci         :Otherwise, return 0
Co Outputs
Co   omppid:thread id or number of threads (see mode)
Cr Remarks
C    OMP PARALLEL directives must bracket call to this routine, e.g.
C    C$OMP PARALLEL PRIVATE(n)
C    C      n = omppid(0)
C    C$OMP END PARALLEL
Cu Updates
Cu   24 Nov 05 Added mode 2
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
      implicit none
      integer mode
C#ifdefC OPENMP
C      integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
C#endif

C#ifdefC OPENMP
C      if (mode == 0) then
C        omppid = OMP_GET_NUM_THREADS()
C      else if (mode == 1) then
C        omppid = OMP_GET_THREAD_NUM()
C      else if (mode < 0) then
C        call OMP_SET_NUM_THREADS(-mode)
C        omppid = OMP_GET_THREAD_NUM()
C      else
C        omppid = 0
C      endif
C#else
      omppid = 0
C#endif

      end
C test
C      subroutine fmain
C      implicit none
C      integer n,omppid
C
CC     n = 0
C      n = omppid(-3)
C
CC#ifdef OPENMP
C      n = omppid(0)
C      print *, 'omppid for number of threads outside brackets:', n
CC$OMP PARALLEL PRIVATE(n)
C      n = omppid(0)
C      print *, 'omppid for number of threads:', n
CC$OMP END PARALLEL
C
C      n = omppid(1)
C      print *, 'omppid(1), outside openmp PARALLEL:', n
C
CC$OMP PARALLEL PRIVATE(n)
C      n = omppid(1)
C      print *, 'omppid(1) for processor id:', n
C!$OMP do
C      do  n = 1, 10
C      print *, 'omppid(1) in do loop: i, thread id =', n,omppid(1)
C      enddo
C!$OMP end do
C
CC$OMP END PARALLEL
CC#endif
C
C      end
