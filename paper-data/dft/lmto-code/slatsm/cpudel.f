      subroutine cpudel(unit,strn,delt)
C- incremental cup time, in seconds
C ----------------------------------------------------------------------
Ci Inputs:
Ci   unit>=0: printout of time to unit; else no printing
Ci   unit<-1: delt not calculated
Co Outputs:
Co   delt: incremental cpu time since last call, seconds
Cr Remarks
Cr   Uses cpusec
C ----------------------------------------------------------------------
C Passed parameters
      character*(*) strn
      integer unit
      double precision delt
C Local parameters
      character*1 timeu, outs*80
      double precision cpusec,told,tnew
      save told
      data told /0d0/

      if (unit < -1) return
      tnew = cpusec()
      delt = tnew - told
      told = tnew
      timeu = 's'
      if (tnew > 60) then
        timeu = 'm'
        tnew = tnew/60
        if (tnew > 60) then
        timeu = 'h'
        tnew = tnew/60
        endif
      endif
      if (unit >= 0 .and. tnew > 0d0) then
        outs = ' '
        write(outs,333) strn
        call awrit2('%a  %1;3,3g%53ptotal:  %1,3;3g'//timeu,
     .    outs,len(outs),-unit,delt,tnew)
      endif
C      if (unit >= 0 .and. tnew > 0d0)
C     .  write(unit,333) strn,delt,tnew,timeu
  333 format(' cpudel',a25,'  time(s):',g10.3,'  total:',f8.3,a1)
      end
C SUPERSEDED by D. Pashov's implementation
C      double precision function cpusec()
CC- process cputime, in seconds
CC ----------------------------------------------------------------------
CCi Inputs:
CCi   none
CCo Outputs:
CCo   returns cpu time, in seconds
CCr  Remarks
CCr    On the Apollo: time (in 4 microsecond units) is
CCr    (time(1)*65536 + time(2))*65536 + time(3)
CC ----------------------------------------------------------------------
CC#ifdefC APOLLO
CC      integer time(3)
CC      integer*2 t2(3)
CC      call proc1_$get_cput(t2)
CC      time(1) = t2(1)
CC      time(2) = t2(2)
CC      time(3) = t2(3)
CC      if (time(1) < 0) time(1) = time(1) + 65536
CC      if (time(2) < 0) time(2) = time(2) + 65536
CC      if (time(3) < 0) time(3) = time(3) + 65536
CC      cpusec = 4d-6 *
CC     .  ((dble(time(1))*65536 + dble(time(2)))*65536 + dble(time(3)))
CC#elseifC AIX
CC      integer mclock
CC      cpusec = dble(mclock())/100
CC#elseifC CRAY
CC      cpusec = second()
CC#elseifC VANILLA
CC      cpusec = 0
CC#elseifC DTIME
CC      real tarray(2),dtime,etime,tsave
CC      data tsave /0d0/
CC      tsave = tsave + dtime(tarray)
CC      cpusec = tsave
CC#else
C      real tarray(2),timet
CC#ifdefC POWER-PC
CC      real etime_
CC      timet = etime_(tarray)
CC#else
C      real etime
C      timet = etime(tarray)
CC#endif
C      timet = tarray(1)
C      cpusec = dble(timet)
CC#endif
C      end

      function cpusec()
        implicit none
        real(8) :: cpusec
        integer(8),save :: c0 = 0, ncalls = 0, cr = 0, cm = 0
        integer(8) :: c

        if (ncalls /= 0  ) then
          call system_clock(c)
        else
          call system_clock(c0, count_rate = cr, count_max = cm)
          c = c0
        endif

        cpusec = (real(c,8) - real(c0,8))/real(cr,8)

        ncalls = ncalls + 1

      end function cpusec
