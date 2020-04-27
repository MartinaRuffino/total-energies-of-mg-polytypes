      subroutine gettime(datim)
      implicit none
      character datim*(*)

      datim = ' '
C#ifdefC SGI
C      call fdate(datim)
C#else
      call ftime(datim)
C#endif
      end
