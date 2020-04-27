C#define FORTRAN2003
      subroutine getarf(iarg,strn)
C- Returns a command-line argument
C ----------------------------------------------------------------------
Ci Inputs
Ci   iarg   :index to command-line argument (0 for executable)
Co Outputs
Co   strn   :string containing iargth command-line argument
Cu Updates
Cu        2011 (DP) use the standard f2003 routine
Cu   17 Jul 01 Written to accomodate problems with HP f90 compiler
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iarg
      character(len=*) strn

C     Get_command_argument is a Fortran 2003 standard.
C#ifdef FORTRAN2003
      call get_command_argument(iarg, strn)
C#elseifC HASGTARGC
CC ... This uses gtargc defined in fmain.c, and requires that the linker
CC     passes command-line arguments to the C main entry point.
C      call gtargc(iarg,strn)
C#elseifC HASGETARG
C      call getarg(iarg,strn)
C#endif
      end

      function nargf()
C- Returns number of command-line argument, including command name
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Co   nargf  :index to the last command line argument.
Cu Updates
Cu    3 Aug 04 First written
C ----------------------------------------------------------------------
      implicit none
      integer nargf

C     Get_command_argument_count is a Fortran 2003 standard.  Otherwise use C main
C#ifdef FORTRAN2003
      nargf = command_argument_count() + 1
C#elseC
CC ... This uses nargc defined in fmain.c, and requires that the linker
CC     passes command-line arguments to the C main entry point.
C#ifdefC HASNARGC
C      integer nargc
C      nargf = nargc()
CC
CC ... This works if a fortran-callable substitute for C iargc available
CC     (e.g. DEC fort)
C#elseifC HASIARGC
C      integer iargc
C      external iargc
C      nargf = iargc() + 1
C
CC ... Another possible fortran-callable substitute for C iargc
CC     (e.g. Intel ifort)
C#elseifC HASNARGS
C      integer nargs
C      external nargs
C      nargf = nargs()
C#endifC
C#endif

      end

      subroutine ftime(datim)
C- Returns date and time
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Co   datim  : string containing date and time.
Co          : On unix systems, this is a 26-character string.
Cu Updates
Cu        2011 (DP) use the standard f2003 routine
Cu    5 Aug 04 First written
C ----------------------------------------------------------------------
      implicit none
      character datim*(*)
      character(len=18) dt

      CALL DATE_AND_TIME(dt(1:8),dt(9:18))

      datim = dt(9:10)//':'//dt(11:12)//':'//dt(13:14)//' '
     & //dt(7:8)//'.'//dt(5:6)//'.'//dt(1:4)

      end

      double precision function delwc()
C- Returns change in wall clock time, in seconds, since last call
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Co   delwc: returns change in wall clock time since last call
Cr Remarks
Cr   Calls dwtime, which is machine dependent.
Cu Updates
Cu   16 May 13 First written
C ----------------------------------------------------------------------
      implicit none
      real(8),save :: tstart = -1d0
      real(8) :: dwtime,tnow

      tnow = dwtime()
      if (tstart == -1) tstart = tnow
      delwc = tnow - tstart
      tstart = tnow

      end

      double precision function dwtime()
C- Returns wall clock time, in seconds, as double precision number
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Co   returns wall clock time, as a floating point number, relative to
Co   start of the month
Cr Remarks
Cr   Uses DATE_AND_TIME, an F90 standard
Cu Updates
Cu   24 Dec 08 First written
C ----------------------------------------------------------------------
      implicit none
      integer itim(8)
      character datim*26

      CALL DATE_AND_TIME (datim(1:8),datim(9:18),datim(19:23),itim)
      dwtime = 24*3600*itim(3) +
     .           3600*itim(5)+60*itim(6)+itim(7)+dble(itim(8))/1000
      end

C#ifdefC TESTF
C      program test
C      integer i
C      character strn*50
C      i = nargf()
C      print *, 'no command line arguments + 1 = ', i
C      call getarf(1,strn)
C      print *, '1st command line argument = ', strn
C      end
C#endif
C#ifdefC TESTC
C      subroutine fmain
C      integer i
C      character strn*50
C      i = nargf()
C      print *, 'no command line arguments + 1 = ', i
C      call getarf(1,strn)
C      print *, '1st command line argument = ', strn
C      end
C#endif

C#ifdefC TDWTIME
C      subroutine fmain
C      implicit none
C      double precision dwtime,tstart,tnow
C
C      tstart = dwtime()
C      print *, 'start',tstart
C      call sleep(5)
C      tnow = dwtime()
C      print *, 'after sleep',tnow
C      call rx0('done')
C
C      end
C
C      subroutine sleep1(j)
C      integer i
C
C      print *, 'starting sleep'
C      do  i = 1, 3000000
C        call sleep2(i,j,l)
C      enddo
C      print *, 'end sleep',l
C      end
C
C      subroutine sleep2(i,j,l)
C      integer i,j,l
C      l = 0
C      do  k = 1, 100*j
C        l = l + j + (-1)**k
C      enddo
C      end
C#endif
