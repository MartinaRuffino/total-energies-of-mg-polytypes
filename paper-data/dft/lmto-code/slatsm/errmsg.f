      subroutine errmsg (messg, iopt)
C- Write error message to standard error device
C ----------------------------------------------------------------
Ci Inputs
Ci   iopt: 0, return without message printed
Ci         1, return with message printed
Ci         2, stop with message printed
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      character*(*) messg

c      if (iopt /= 0) write(i1mach(4),*) messg
      if (iopt /= 0 .and. iprint() >= 40) write(i1mach(4),*) messg
      if (iopt < 2) return
      end
