C#define unix
      subroutine finits(job,fcn,fcargs,iarg)
C- Machine and compiler-dependent inits for standard FORTRAN startup
C ----------------------------------------------------------------
Ci Inputs
Ci   job:  0, no command-line arguments;
Ci         1, no switches (extens only);
Ci         2, [-vnam=val ...] [-pr#] [switches] extens;
Ci         3, call fcn for switches first
Ci         10s digit:
Ci         1, allow ctrl.extens in place of extens
Ci   fcn,fcargs:  see Remarks
Co Outputs
Co   iarg: argument to extens, if found (0 if not)
Cr Remarks
Cr   finits parses command line arguments, ignoring args that
Cr   begin with "-", to find file extension.  When job=1,
Cr   switches -vnum are taken to be variable defs.
Cr   Calling prog can set own switches through external
Cr   fcn(iarg,fcargs).  fcn should return iarg as last arg parsed.
Cu Updates
Cu   15 Jan 10 Added 10s digit job
Cu   16 Jun 06 Added MPI Wall clock time printout
Cu    3 Aug 04 Changed call to nargc with call to nargf
C ----------------------------------------------------------------
C#ifndef LINUXF
      use posix
      use :: iso_c_binding, only : c_null_char
C#endif
      use mpi
      implicit none
C Passed parameters
      integer job,iarg
      double precision fcargs(1)
      external fcn
C#ifdef unix
      logical lsequ,lext
      integer i,fext,nargf,n,it(5),iv(5),a2vec,job0,job1,index,dot_idx
      character strn*256
      character(len=256) :: cstrn
      character*20 extns
C#endif
      double precision tstart,delwc,dwtime
C     common /finitp/ tstart

      integer master,procid
      double precision starttime, endtime
      parameter (master = 0)

      common /mpifinits/ starttime, endtime
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, i )
      if (procid == master) then
        starttime = MPI_WTIME()
      endif

      call awrit0('%u',' ',80,0)
      call initqu(.false.)
      tstart = delwc()              ! Initialize wall clock time
      tstart = dwtime()             ! For total wall clock time

C --- For Lahey F77L, open standard output as list-directed ---
C#ifdefC F77LAHEY
C      open(unit=*,carriage control='LIST')
C#endif

C --- Handle floating point exceptions in the IBM VM environment ---
C#ifdefC IBM_VM
C      call errset(208,999,-1)
C#endif

C --- Command line arguments and extension ---
      job0 = mod(job,10)
      job1 = mod(job/10,10)
      if (job0 == 0) return
      iarg = 0
      if (job0 == 3) call fcn(fcargs,iarg)
      lext = .false.
   10 iarg = iarg+1
      if (nargf() > iarg) then
        call getarf(iarg,strn)
        extns = strn
        if (job0 >= 2) then
C     ... -pr encountered ... set verbosity
          if (lsequ(strn,'-pr',2,' ',n)) then
            i = 3
            n = a2vec(strn,len(strn),i,2,', ',2,2,5,it,iv)
            do  12  i = 1, n
   12       call sprt(i-1,iv(i))
          endif
C     ... v encountered ... parse variables
          if (lsequ(strn,'-v',2,' ',n)) then
            i = 2
            call parsyv(strn,len(strn),999,0,i)
C           call shosyv(0,0,0,i1mach(2))
          endif
        endif
        if (lsequ(strn,'-',1,' ',n)) goto 10
        if (.not. lext) then
          dot_idx = index(strn,'.')
          if (job1 == 1 .and. strn(1:5) == 'ctrl.') then
            i = fext(strn(5:))
          elseif (job1 == 2 .and. dot_idx > 1) then
            i = fext(strn(dot_idx:))
C#ifndef LINUXF
            if (dot_idx > 5) then
                cstrn(1:dot_idx-4) = strn(1:dot_idx-5)//c_null_char
                i = chdir(cstrn)
                if (i /= 0)
     &         call rxi('Could not change path with chdir, err code:',i)
            end if
C#endif
          elseif (dot_idx == 1) then
            i = fext(extns)
          else
            i = fext('.'//extns)
          endif
          lext = .true.
          goto 10
        endif
      endif
      iarg = iarg-1
      end

      subroutine fexit(retval,iopt,strng,args)
C- Machine and compiler-dependent program termination
C ----------------------------------------------------------------------
Ci Inputs
Ci   retval:  return value passed to operating system
Ci   iopt decomposed into 3 one-digit numbers.
Ci   digit
Ci     1:  0: do not print string on exit;
Ci         9: print strng as Exit(retval): 'strng'
Ci      else: exit, using strn as a format statement and args a vector
Ci            of  c  double precision arguments
Ci    10:   0: do not print cpu time, else do
Ci   100:   0: do not print work array usage, else do
Co Outputs
Cr Remarks
Cu Updates
Cu   13 May 17 prints out number of warnings logged (see logwarn)
Cu   03 Jul 03 Open file tmp with fopnT, for local directories (MPI)
C ----------------------------------------------------------------------
      use mpi
      implicit none
C Passed parameters
      integer retval,iopt
      character*(*) strng
      double precision args(*),arg2(*),arg3(*)
C Local parameters
      integer i,i2,scrwid,fopnT,ientry
      parameter (scrwid=120)
      double precision tnew,wtime
      character*1 timeu
      character(len=256) :: strn
      character(len=26) :: datim
      character(len=20) :: hostnm
      integer master,procid
      parameter (master=0)
      procedure(logical) :: isopen
      procedure(integer) :: fopn,fhndl,iprint,getdig,i1mach,mpipid
      procedure(real(8)) :: cpusec

      double precision starttime, endtime
      common /mpifinits/ starttime, endtime

      ientry = 1
      goto 5

      entry fexit3(retval,iopt,strng,args,arg2,arg3)
      ientry = 3
      goto 5

      entry fexit2(retval,iopt,strng,args,arg2)
      ientry = 2

    5 continue
      procid = mpipid(1)
      if (procid == master) then
      i = getdig(iopt,0,10)
      if (i /= 0) then
        if (i == 9) then
          i = 0; call logwarn(i,'')
          if (i > 0) then
            call awrit2('%x Exit %i [%i%-1j severe warning%?;(n>1);s;;] '//strng,strn,len(strn),0,retval,i)
          else
            strn = ' Exit %i '//strng
          endif
          if (iprint() > 0) then
C            call strip(strn,i,i2)
C            i2 = min(i2,scrwid)
            call awrit1(trim(strn),' ',scrwid,i1mach(2),retval)
          endif
C          print 345, retval, strng(1:min(len(strng),scrwid-9))
C  345     format(' Exit',i3,' ',a)
        else
          if (iprint() > 0) then
             if (ientry == 3) call awrit3(strng,strn,-scrwid,i1mach(2),args,arg2,arg3)
             if (ientry == 2) call awrit2(strng,strn,-scrwid,i1mach(2),args,arg2)
             if (ientry == 1) call awrit1(strng,strn,-scrwid,i1mach(2),args)
          endif
C         strn = strng
C         print strn, retval, (args(j), j=1,i)
        endif
      endif

      i = getdig(iopt,1,10)
      if (i /= 0 .and. iprint() >= 10 .and. cpusec() /= 0) then
        timeu = 's'
        tnew = cpusec()
C       tstart = dwtime()-tstart
        endtime = mpi_wtime()
        wtime = endtime-starttime
        if (tnew > 3600) then
          timeu = 'm'
          tnew = tnew/60
          wtime = wtime/60
          if (tnew > 200) then
            timeu = 'h'
            tnew = tnew/60
            wtime = wtime/60
          endif
        endif

        datim = ' '
        call ftime(datim)
        hostnm = ' '
        call gtenv('HOST',hostnm)
        call word(hostnm,1,i,i2)
        i2 = max(i,i2)
C       write(i1mach(2),10) tnew,timeu,wtime,timeu,datim,hostnm(i:i2)
        strn = ' '
        call awrit2(' CPU time:  %,3;3d'//timeu//'   Wall clock %,3;3d'//timeu//
     .    '  at  '//trim(datim)//'  on  '//hostnm(i:i2),strn,len(strn),0,tnew,wtime)
        call awrit0('%a',strn,-len(strn),-i1mach(2))
        if (fhndl('LOG') >= 0)  then
          if (isopen(fhndl('LOG'),.false.))
     .      call awrit0('%a',strn,-len(strn),-fhndl('LOG'))
        endif
      endif
      endif

      i = fopnT('TMP' ,-1,0,2)
      if (i >= 0) call dfclos(i)
      if (fhndl('TMP') >= 0) call dfclos(fopn('TMP'))
      if (procid == master) then
      if (getdig(iopt,2,10) > 0 .and. iprint() > 0) call wkinfo
      endif

C   11 call tclev(hostnm,i)
C      if (i >= 0) then
C        call tcx(hostnm)
C        goto 11
C      endif
      call tcprt(i1mach(2))

C#ifndef IBM_VM
      call closea
C#endif
C      if ( procid == master ) then
C        endtime = MPI_WTIME()
C        call awrit2('%N Wall-clock time: %;3ds. Resolution %;9ds%N',
C     .              ' ',128,i1mach(2),endtime-starttime,MPI_WTICK())
C      endif
!      call MPI_FINALIZE(i)

C#ifdef unix
      call cexit(retval,1)
C#endif
      stop
      end
      subroutine rx(string)
C- Error exit
      implicit none
      character*(*) string

      call fexit(-1,119,string,0d0)
      end
      subroutine rx0(string)
C- Normal exit
      implicit none
      character*(*) string

      call fexit(0,119,string,0d0)
      end
      subroutine rx1(string,arg)
C- Error exit, with a single argument
      implicit none
      character*(*) string
      double precision arg
      character*120 outs
      outs = '%N Exit -1 '//string
      call fexit(-1,111,outs,arg)
      end
      subroutine rx2(string,arg1,arg2)
C- Error exit, with two arguments
      implicit none
      character*(*) string
      double precision arg1,arg2
      character*120 outs
      outs = '%N Exit -1 '//string
      call fexit2(-1,111,outs,arg1,arg2)
      end
      subroutine rxi(string,arg)
C- Error exit, with a single integer at end
      implicit none
      character*(*) string
      double precision arg
      character*120 outs
      outs = '%N Exit -1 '//string//' %i'
      call fexit(-1,111,outs,arg)
      end

      subroutine rxs(string,msg)
C- Error exit with extra string message
      implicit none
      character*(*) string,msg
      character*120 outs
      integer i
      outs = string // trim(msg)
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      subroutine rxs2(string,msg,msg2)
C- Error exit with extra string messages
      implicit none
      character*(*) string,msg,msg2
      character*120 outs
      integer i
      outs = string // msg // msg2
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      subroutine rxs4(string,msg,msg2,msg3,msg4)
C- Error exit with extra string messages
      implicit none
      character*(*) string,msg,msg2,msg3,msg4
      character*120 outs
      integer i
      outs = string // msg // msg2 // msg3 // msg4
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      subroutine rxx(test,string)
C- Test for error exit
      implicit none
      logical test
      character*(*) string

      if (test) call rx(string)
      end

C#ifndefC unix
C      subroutine ftimex(datim)
C      character*(*) datim
C      datim = ' '
C      end
C#endif
CC tests finits and fexit
C      subroutine fmain
CC      implicit none
C      integer iarg
C      double precision arg1(10)
C      external cmdarg
C
C      call finits(3,cmdarg,arg1,iarg)
CC     call fexit (-1,111,' Test fexit, one argument:%2:2;6d',arg1)
C      call fexit2(-1,111,' Test fexit, two args:%2:2;6d : %i',arg1,987)
C      end
C      subroutine cmdarg(cmargs,iarg)
CC- Command line arguments special to lmtoft
CC ----------------------------------------------------------------
CCi Inputs
CCi   cmargs(1,2): exi (if -e)
CCi         (3):   lmxf for getoro (i.e. n in -gn)
CCo Outputs
CCi   iarg
CCr Remarks
CC ----------------------------------------------------------------
CC      implicit none
CC Passed parameters
C      integer iarg
C      double precision cmargs(3)
CC Local variables
C      logical lsequ,a2bin
C      character strn*40
C      integer i,n,nargf
C
CC --- Command line arguments and extension ---
C      iarg = 0
C   10 iarg = iarg+1
C      if (nargf() > iarg) then
C        call getarf(iarg,strn)
C        i = 2
C        if (lsequ(strn,'-v',i,' ',n)) call parsyv(strn,40,999,0,i)
C        if (lsequ(strn,'-e',i,' ',n)) then
C          iarg = iarg+2
C          if (nargf() <= iarg) goto 99
C          call getarf(iarg-1,strn)
C          i = 0
C          if (.not. a2bin(strn,cmargs,4,0,' ',i,-1)) goto 99
C          call getarf(iarg,strn)
C          i = 0
C          if (.not. a2bin(strn,cmargs,4,1,' ',i,-1)) goto 99
C        endif
C        if (lsequ(strn,'-g',i,' ',n)) then
C          if (strn(3:3) /= ' ') then
C            i = 2
C            if (.not. a2bin(strn,cmargs,4,2,' ',i,-1)) goto 99
C          endif
C        endif
C        if (lsequ(strn,'-',1,' ',n)) goto 10
C      endif
C      iarg = iarg-1
C      return
C   99 call rx('error parsing switches')
C      end
