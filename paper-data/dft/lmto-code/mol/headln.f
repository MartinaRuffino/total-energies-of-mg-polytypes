      subroutine headln(name)
C- Puts a heading line into protocol
      character*(*) name
      write(6,300) name
  300 format(' -------------------------  START ',a,
     .       ' -------------------------')
      end
