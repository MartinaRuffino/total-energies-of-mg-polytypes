      subroutine readrssw(strn,is,irs)
C- Find a command-line argument --rs
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Cio  irs    : On input, if irs(1) < 0, to not scan irs
Cio         : Otherwise, read irs if command line argument is found
Co Outputs
Co   strn  : Command line argument
Co   is    : If argument found, returns index to last terminator in strn
Co         : If argument not found, returns -1
Cr Remarks
Cr   Written to find command-line argument --rs,
Cr   avoiding conflict with other arguments beginning with --rs
Cu Updates
Cu   15 Aug 17 First created
C ----------------------------------------------------------------------
      implicit none
      integer is,irs(5)
      character strn*(*)
      character dc*1
      integer i,ix(10)
      logical lfail
      procedure(logical) :: cmdstr
      procedure(integer) :: a2vec,wordsw

      is = -1
      i = 1
      do while (cmdstr(i,strn))
        i = i+1
        if (strn(1:4) /= '--rs') cycle
        dc = strn(5:5)
        lfail = dc >= 'A' .and. dc <= 'Z' .or.
     .          dc >= 'a' .and. dc <= 'z' .or.
     .          dc >= '0' .and. dc <= '9'
        if (lfail) cycle  ! --rs is the beginning of a different argument
        is = scan(strn,dc,BACK=.true.); i = is
        if (irs(1) < 0) exit
        if (a2vec(strn,len(strn),i,2,', ',2,2,5,ix,irs) < 1) call rxs('readrssw failed to parse ', strn)
        exit
      enddo
      end
