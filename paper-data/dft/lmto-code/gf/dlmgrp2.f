      subroutine dlmgrp2(nclass,nth,icdl,igrp)
C- Reset the grp2 array for DLM angle classes
      implicit none
C ... Passed parameters
      integer nclass,nth(nclass),icdl(nclass),igrp(*)
C ... Local parameters
      integer ic,ith,ifi,nangl,lgunit,isum,iprint
      logical flag
      integer mpipid,procid,master
      character*80 outs
      double precision xx

C ... Assign new group numbers to DLM classes
      procid = mpipid(1)
      master = 0

      flag = .false.
      do  ic = 1, nclass
        if (igrp(ic) == 0 .or. nth(ic) < 2) cycle
        flag = .true.
        do  ith = 1, nth(ic)
          igrp(icdl(ic) + ith - 1) = igrp(ic) * (ith + 200)
        enddo
      enddo
      nangl = isum(nclass,nth,1)

      if (flag .and. procid == master .and. iprint() >= 20) then
        print *
        outs = ' Class  grp2'
        do  ic = 1, 1   ! Don't write to log file
          ifi = lgunit(ic)
          if (ic == 1) write (ifi,901)
          call arrprt(outs,'%,5i%,6i','Ii',nclass+nangl,ifi,7,0,'  | ',
     .      xx,igrp,xx,xx,xx,xx,xx,xx)
  901     format (' DLMGRP2: New groups assigned - please check:')
          if (mpipid(0) > 1) exit ! this looks quite shady, hopefully we never come to it ...
        enddo
      endif

      end
