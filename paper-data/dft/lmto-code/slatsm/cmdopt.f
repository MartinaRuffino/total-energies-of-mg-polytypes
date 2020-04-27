      integer function cmdoptsw(cmdstr,sep,substr,term)
C- Determines whether a whether a substr is embedded in a given command-line argument
C ----------------------------------------------------------------
Ci Inputs
Ci   cmdstr: Find command-line beginning with cmdstr.
Ci         : The character following cmdstr acts as a separator between substring arguments.
Ci   sep   : Other separators preceding substr, to be appended to the separator
Ci         : following cmdstr.  Can be empty, in which case it is not used.
Ci   substr: string embedded within strn
Ci   term  : terminator string.  If not zero length, the character immediately
Ci         : following substr must be a character in term.
Co Outputs
Co   cmdoptsw returns:
Co          -1 if command-line argument not found
Co           i if argument contains sep//substr
Co             where i is the index to the element in strn
Co             AND term has no characters
Co             OR  character following substr is among characters in term
Co           0 otherwise
Cr Remarks
Cr   Example: command-line argument is 'readsek;ib=3:6'
Cr            cmdstr = 'readsek' sep = ';' and substr = 'ib' and term is '='
Cr   See also wordsw
Cu Updates
Cu    1 Apr 19 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) cmdstr,substr,term,sep
C ... Local parameters
      integer i,itrm
      character strn*256,dc*1
      procedure(logical) :: cmdopt
      procedure(integer) :: wordsw

      cmdoptsw = -1
      if (.not. cmdopt(cmdstr,len(cmdstr),0,strn)) return
      i = len(cmdstr)+1
      dc = strn(i:i)
      cmdoptsw = min(wordsw(strn,dc//sep,substr,term,itrm),1)

      end

      integer function cmdoptswx(cmdstr,substr,term)
C- Determines whether a whether a substr is embedded in a given command-line argument
C ----------------------------------------------------------------
Ci Inputs
Ci   cmdstr: Find command-line beginning with cmdstr
Ci   substr: string embedded within command-line argument.
Ci   term  : terminator string.  If not zero length, the character immediately
Ci         : following substr must be a character in term.
Co Outputs
Co   cmdoptswx returns:
Co          -1 if no command line argument begins with cmdstr
Co           i if command line argument beginning with cmdstr contains substr
Co             where i is the index to the element in the command-line argument
Co             substr may have zero length
Co           0 otherwise
Cr Remarks
Cr   Example:  cmdstr = '--opt' and substr = ':write'
Cr   See also wordsw
Cu Updates
Cu    1 Apr 19 Renamed to cmdoptswx (was cmdoptsw)
Cu   24 Mar 17 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) cmdstr,substr,term
C ... Local parameters
      integer i
      character*256 strn
      procedure(logical) :: cmdopt

      cmdoptswx = -1
      if (.not. cmdopt(cmdstr,len(cmdstr),0,strn)) return
      i = index(strn,substr)
      if (len(substr) == 0) i = 1
      cmdoptswx = i
      if (i == 0 .or. len(term) == 0) return  ! If no match or no terminator
      i = i+len(substr)
      i = scan(strn(i:i),term)
      if (i == 0) cmdoptswx = 0  ! Terminator not found

      end
      logical function cmdopt(argstr,strln,nargs,outstr)
C- Determines whether a command-line argument supplied, and its argument
C ----------------------------------------------------------------
Ci Inputs
Ci   argstr: command-line string to search
Ci   strln : search to strln chars.
Ci         : If strnl == NULLI use len(argstr) for strln
Ci   nargs : number of arguments associated with argstr
Co Outputs
Co   cmdopt: true if argument found, else false
Co   outstr: (nargs>0 only) nth string after string argstr
Cr Remarks
Cu   03 May 19 When strln=NULLI, use length(argstr) as substitute
Cu    3 Aug 04 Changed call to nargc with call to nargf
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*(*) argstr,outstr
      integer nargs,strln
C Local parameters
      integer iarg,idum,nxarg,strlnl
      character*160 strn
      integer, parameter :: NULLI=-99999
      procedure(logical) :: lsequ
      procedure(integer) :: nargf

      cmdopt = .false.
      strlnl = strln
      if (strln == NULLI) strlnl = len(argstr)

      iarg = 0
   10 iarg = iarg+1
C ... A usual command-line argument
      if (nargf() > iarg) then
        call getarf(iarg,strn)
C ... If not on the command-line, try 'extra' arguments
      else
        call ncmdop(nxarg)
*       print *, nxarg,iarg-nargf()
        if (nxarg <= iarg-nargf()) return
        call gcmdop(iarg-nargf()+1,strn)
      endif
      if (.not. lsequ(strn,argstr,strlnl,' ',idum)) goto 10
      cmdopt = .true.
      outstr = ' '
      if (nargf() > iarg+nargs) then
        call getarf(iarg+nargs,outstr)
      elseif (nargf() > iarg+nargs+nxarg) then
        call rx('bug in CMDOPT')
      else
        call gcmdop(iarg-nargf()+1,outstr)
      endif

      end

      subroutine acmdop(strn,lstr,opt)
C- Append strn to 'extra' command options
C ----------------------------------------------------------------
Ci Inputs
Ci   strn:  (acmdop,opt=0) is appended to internal cmdarg
Ci   lstr:  length of input string
Ci   opt:   0 append strn
Ci          1 print out table
Ci   iarg   (gcmdop) retrieve argument n
Co Outputs
Co   n       (ncmdop) number of arguments in list
Co   n       (gcmdop) number of arguments in list
Cr Remarks
C ----------------------------------------------------------------
      implicit none
      integer lstr,opt
      character*1 strn(1), sout*(*)
C Local variables
      integer mxarg,lcmd
C#ifdefC SUN-ULTRA
C      parameter (mxarg=2000,lcmd=2048)
C#else
      parameter (mxarg=2000,lcmd=20000)
C#endif
      integer marker(0:mxarg),nxarg,i1,i2,it,ia,i,lgunit,n,iarg
      character*(lcmd) cmdarg, ch*3
      save nxarg,marker
      data cmdarg /' '/ ch /' "'''/, nxarg /0/

*     print *, (strn(i), i=1,lstr)
      marker(0) = 1
      if (opt == 1) goto 100

      i2 = -1
   10 continue
      i1 = i2+1
C --- Continue until all command arguments exhausted ---
      call skipbl(strn,lstr,i1)
      if (i1 < lstr) then
        ia = marker(nxarg)
C   ... Find i2 : points to past last char of the argument
        i2 = i1
   12   i1 = i2
        call chrps2(strn,ch,3,lstr,i2,it)
        call strcop(cmdarg(ia:),strn(i1+1),i2-i1,' ',i)
        ia = ia+i
        if (ia >= lcmd) call rx('acmdop: increase lcmd')
*       print *, cmdarg(1:ia)
C   ... A quote encountered ... continue copying string
        if (it > 1) then
          i2 = i2+1
          i1 = i2
          call chrpos(strn,ch(it:it),lstr,i2)
          call strncp(cmdarg,strn,ia,i1+1,i2-i1)
          ia = ia+i2-i1
          if (ia >= lcmd) call rx('acmdop: increase lcmd')
*          print *, cmdarg(1:ia)
          i2 = i2+1
          if (i2 <= lstr) goto 12
        endif
C   ... End of this argument ... start on another
        nxarg = nxarg+1
        if (nxarg > mxarg) call rx('acmdop: increase mxarg')
        marker(nxarg) = ia
        goto 10
      endif
      if (opt == 0) return

C  --- Printout of command line arguments ---
  100 continue
      call awrit1(' acmdop:  %i file command switches:',
     .  ' ',80,lgunit(1),nxarg)
      do  20  i = 1, nxarg
   20 print 333, i, cmdarg(marker(i-1):marker(i)-1)
  333 format(i4,2x,'"',a,'"')
      return

      entry ncmdop(n)
      n = nxarg
      return

      entry gcmdop(iarg,sout)
      if (iarg > nxarg) then
        sout = ' '
        return
      endif
      sout = cmdarg(marker(iarg-1):marker(iarg)-1)
      return

      end
