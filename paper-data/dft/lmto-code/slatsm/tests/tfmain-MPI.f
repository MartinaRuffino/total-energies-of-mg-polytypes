C     Tests whether fmain is reached, and whether command line args are read, MPI
C     To test, invoke with e.g.
C     mpirun -n 4 ./a.out first second third
      subroutine fmain
      implicit none
      integer nargf,iarg,nargs,nproc,ichk(5)
      character*40 strn

      integer mpipid,procid

      nproc = mpipid(0)
      procid = mpipid(1)
      print 333, ' entered fmain, pid= ',procid,
     .  ' total number of proc=',nproc
  333 format(a,i4,a,i4)

      ichk(1) = mpipid(0)
      ichk(2) = mpipid(3)
      ichk(3) = mpipid(0)
      ichk(4) = mpipid(4)
      ichk(5) = mpipid(0)
      print "(' checking modes in mpipid',5i3)", ichk

      nargs = nargf()
      print 333, ' nargf found', nargs, ' arguments, pid=',procid
      do  iarg = 0, nargs-1
        call getarf(iarg,strn)
        print 334, procid, iarg, strn
  334   format(' procid=',i4,' iarg=',i4,2x,a)
      enddo

      end
