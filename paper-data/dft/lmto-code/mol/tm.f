      subroutine tm(string)
C- Print out CPU time and time from last call to tm with message
      character*(*) string
      double precision x
      integer i1mach,iprint
      if (string == ' ') then
        call cpudel(-1,string,x)
      else
        call cpudel(i1mach(2),string,x)
      endif
c      save ncall,cpu0,i0
c      data ncall /0/, cpu0/0.0/
c      call vttime(i1,i2)
c      if(ncall == 0) i0=i1
cc|         write(6,*) 'tm  i0=',i0,'  i1=',i1
c      cpu=(i1-i0)/100.0
c      write(6,890) cpu,cpu-cpu0,string
c  890 format(' tm:',f12.4,f12.4,2x,a)
c      ncall=ncall+1
c      if(ncall >= 1) cpu0=cpu
c      return
      end
