      subroutine fmain
C- Tests idalloc
      integer idalloc,i

C     i =  idalloc(' ',0,1,1)
      i =  idalloc('abc',10+2,19000,3900)

c     i = idasiz(1)
      i =  idalloc('c',0+2,433,300)  ! Doesn't do anything --- too small

      i =  idalloc('c',10+1,4333,3*300)
      i =  idalloc('c',10+4,4333,3*300)
      i =  idalloc('c',10+3,4333,3*300) ! No effect; already deallocated

      i =  idalloc('bc',0+0,1,1)
      i =  idalloc('abcd',10+0,1,1)
      i =  idalloc('a',10+0,1,1)
      call pshpr(0)
      i =  idalloc('abcde',10+1,12345,54321)
      i =  idalloc('abcde',10+4,12345,54321)
      call poppr


      i =  idalloc('bc',10+1,1000,2000)
      i =  idalloc('a',20+1,80000,2000)
      i =  idalloc('abc',10+4,19000,3900)
      i =  idalloc(' ',10,1,1)
      print "(' idalloc returned',i6)", i

      end
