      subroutine getjobmode(mode)
C- Returns # to command-line switch --job=#, and NULLI if it does not exist
      implicit none
      integer mode
      procedure(logical) :: cmdopt
      character*160 strn
      integer,parameter :: NULLI=-99999

      mode = NULLI
      if (cmdopt('--job=',6,0,strn)) then
        if (len_trim(strn) <= 6) return
        read(strn(7:),*) mode
      endif

      end

      integer function cmdoptsw(cmdstr,substr,term)
C- Determines whether a whether a substr is embedded in a given command-line argument
C ----------------------------------------------------------------
Ci Inputs
Ci   cmdstr: Find command-line beginning with cmdstr
Ci   substr: string embedded within command-line argument.
Ci   term  : terminator string.  If not zero length, the character immediately
Ci         : following substr must be a character in term.
Co Outputs
Co   cmdoptsw returns:
Co          -1 if no command line argument begins with cmdstr
Co           i if command line argument beginning with cmdstr contains substr
Co             where i is the index to the element in the command-line argument
Co             substr may have zero length
Co           0 otherwise
Cr Remarks
Cr   Example:  cmdstr = '--opt' and substr = ':write'
Cr   See also wordsw
Cu Updates
Cu   24 Mar 17 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) cmdstr,substr,term
C ... Local parameters
      integer i
      character*256 strn
      procedure(logical) :: cmdopt

      cmdoptsw = -1
      if (.not. cmdopt(cmdstr,len(cmdstr),0,strn)) return
      i = index(strn,substr)
      if (len(substr) == 0) i = 1
      cmdoptsw = i
      if (i == 0 .or. len(term) == 0) return  ! If no match or no terminator
      i = i+len(substr)
      i = scan(strn(i:i),term)
      if (i == 0) cmdoptsw = 0  ! Terminator not found

      end
      logical function cmdopt(argstr,strln,nargs,outstr)
C- Determines whether a command-line argument supplied, and its argument
C ----------------------------------------------------------------
Ci Inputs
Ci   argstr,strln: command-line string to search; search to strln chars
Ci   nargs:        number of arguments associated with argstr
Co Outputs
Co   cmdopt: true if argument found, else false
Co   outstr (nargs>0 only) nth string after string argstr
Cr Remarks
Cu    3 Aug 04 Changed call to nargc with call to nargf
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*(*) argstr,outstr
      integer nargs,strln
C Local parameters
      integer iarg,nargf,idum,nxarg,ioff,jarg
      character*160 strn

      cmdopt = .false.
      nxarg = nargf()

C --- For each command line argument, check for match ---
      do  iarg  = 0, nxarg-1-nargs
        call getarf(iarg,outstr)
C   ... If not on the command-line, try 'extra' arguments
        if (outstr(1:strln) /= argstr(1:strln)) cycle
        cmdopt = .true.
        do  jarg = 1, nargs
          ioff = len_trim(outstr)
          call getarf(iarg+jarg,outstr(2+ioff:))
        enddo
        return
      enddo
C      if (nargf() > iarg+nargs) then
C        call getarf(iarg+nargs,outstr)
C      elseif (nargf() > iarg+nargs+nxarg) then
C        call rx('bug in CMDOPT')
C      else
C        call gcmdop(iarg-nargf()+1,outstr)
C      endif

      end

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

      call get_command_argument(iarg, strn)
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

C     Get_command_argument_count is a Fortran 2003 standard.
      nargf = command_argument_count() + 1

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
Cr   Fortran has no standard for this function, so it is
Cr   machine-dependent.  It may be that the compiler supplies this
Cr   function in a system call.  If this approach works, it is the
Cr   preferred method.  If no scheme is supplied, time is always 0
Cu Updates
Cu   24 Dec 08 First written
C ----------------------------------------------------------------------
      implicit none
!       integer itim(8)
!       character datim*26
      double precision cpusec

! C#ifdef DATE_AND_TIME
!       CALL DATE_AND_TIME (datim(1:8),datim(9:18),datim(19:23),itim)
!       dwtime = 24*3600*itim(3) +
!      .           3600*itim(5)+60*itim(6)+itim(7)+dble(itim(8))/1000
! C#elseC
! C      dwtime = cpusec()
! C#endif
        dwtime = cpusec()
      end

      function cpusec()
!       return wall clock time in seconds.
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


!       double precision function cpusec()
! C- process cputime, in seconds
! C ----------------------------------------------------------------------
! Ci Inputs:
! Ci   none
! Co Outputs:
! Co   returns cpu time, in seconds
!
!       call cpu_time(cpusec)
!
!       end function cpusec

      integer function omppid(mode)
C- Returns OPENMP thread id
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return number of threads
Ci         :1 return thread id
Ci         :<2 call OMP_SET_NUM_THREADS(-mode)
Ci         :Otherwise, return 0
Co Outputs
Co   omppid:thread id or number of threads (see mode)
Cr Remarks
C    OMP PARALLEL directives must bracket call to this routine, e.g.
C    C$OMP PARALLEL PRIVATE(n)
C    C      n = omppid(0)
C    C$OMP END PARALLEL
Cu Updates
Cu   24 Nov 05 Added mode 2
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
      implicit none
      integer mode
#ifdef OPENMP
      integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      if (mode == 0) then
        omppid = OMP_GET_NUM_THREADS()
      else if (mode == 1) then
        omppid = OMP_GET_THREAD_NUM()
      else if (mode < 0) then
        call OMP_SET_NUM_THREADS(-mode)
        omppid = OMP_GET_THREAD_NUM()
      else
        omppid = 0
      endif

#else
      omppid = 0
#endif

      end
