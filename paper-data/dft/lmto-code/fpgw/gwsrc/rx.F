      subroutine rx(string)
c  routine for error exit
      character*(*) string
      write(6,892) string
      write(6,890)
  890 format(' ---- Error exit')
  892 format(/' ---- ',a)
      stop
      end
