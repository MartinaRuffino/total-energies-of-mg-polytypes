      subroutine spchar(i,ch)
C- Returns special characters
C ----------------------------------------------------------------------
Ci Inputs
Ci   i=1  return backward slash in ch
Ci   i=2  return tab in ch
Ci   i=3  return newline in ch
Ci Outputs
Co   ch
Cu Updates
Cu   04 Aug 13 D. Pashov replaced with architecture-independent implementation
C ----------------------------------------------------------------------
      implicit none
      integer i
      character ch
      character, parameter :: chs(3)=(/achar(92), achar(9), achar(10)/)

      ch = chs(i)

      end
C#ifdefC TEST
C      character ch*1
C      call spchar(1,ch)
C      print 333, 'backslash',ch
C      call spchar(2,ch)
C      print 333, 'tab',ch
C      call spchar(3,ch)
C      print 333, 'newline',ch
C  333 format(a,' :',a1,':')
C      end
C#endif
