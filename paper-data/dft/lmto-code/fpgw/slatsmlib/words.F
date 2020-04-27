      subroutine words(str,nw)
C- Count blank-delimited words in str
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Co Outputs
Co   nw    :number of blank-delimited words in str
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) str
      integer nw
C ... Local parameters
      integer i1,i2,i0,i

      nw = 0
      i1 = 0
      i2 = 0
      i0 = 1
   99 do  10  i = i0, len(str)
        if(str(i:i) /= ' ') then
           i1 = i
           goto 90
        endif
   10 continue
      return
   90 nw = nw+1
      do  20  i = i1,len(str)
        if(str(i:i) == ' ') then
          i2 = i
          goto 91
        endif
   20   continue
      return
   91 i0 = i2
      goto 99
      end
      subroutine word(str,iw,j1,j2)
C- Returns j1,j2 so that str(j1:j2) is the iw-th word from beginning
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Ci   iw    :find iw-th word
Co Outputs
Co   j1    :str(j1:j2) is iw-th word
Co   j2    :-//-
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) str
      integer iw,j1,j2
C ... External calls
      external nword
      j1 = 1
      call nword(str,iw,j1,j2)
      end

      subroutine nword(str,iw,j1,j2)
C- Returns j1,j2 so that str(j1:j2) is the iw-th word from current pos
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Ci   iw    :find iw-th word
Ci   j1    :start search from str(j1:)
Co Outputs
Co   j1    :str(j1:j2) is iw-th word
Co   j2    :-//-
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iw,j1,j2
      character*(*) str
C ... Local parameters
      integer nw,i1,i2,i0,i
      nw = 0
      i1 = 0
      i2 = 0
      i0 = j1
      j2 = -1
  99  do  10  i = i0, len(str)
C   ... skip until nonblank char
        if(str(i:i) /= ' ') then
          i1 = i
          goto 90
        endif
  10  continue
      return
C   ... skip until a blank char
  90  nw = nw+1
      if (nw == iw) j1 = i1
      do  20  i = i1, len(str)
        if(str(i:i) == ' ') then
          i2 = i
          goto 91
        endif
   20 continue
C ... We have reached the end of the string
      if (nw == iw) j2 = len(str)
      return
C ... cleanup: exit if word sought, else try again
   91 i0 = i2
      if (nw == iw) then
        j2 = i2-1
        return
      endif
      goto 99
      end

      integer function wordsw(strn,sep,substr,term,itrm)
C- Determines whether a whether a substring is embedded in a given string
C ----------------------------------------------------------------
Ci Inputs
Ci   strn  : string to search
Ci   sep   : separator preceding strn.  Can be empty, in which case it is not used.
Ci   substr: string embedded within strn
Ci   term  : terminator string.  If not zero length, the character immediately
Ci         : following substr must be a character in term.
Co Outputs
Co   wordsw: returns:
Co           i if strn contains sep//substr
Co             where i is the index to the element in strn
Co             AND term has no characters
Co             OR  character following substr is among characters in term
Co           0 otherwise
Co  itrm   : index to termination of substr
Cr Remarks
Cr   Example:  strn = 'readsek;ib=3:6' and sep = ';' and substr = 'ib'
Cr   See also cmdoptsw
Cu Updates
Cu   15 Apr 17 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) strn,substr,term,sep
      integer itrm
C ... Local parameters
      integer i,n

      if (len(sep) > 0) then
        i = index(strn,sep(1:1)//substr)
        n = 1
      else
        i = index(strn,substr)
        n = 0
      endif
      wordsw = i
      itrm = i+len(substr)+n
      if (i == 0 .or. len(term) == 0 .or. itrm > len(substr)) return  ! If no match or no terminator
      if (scan(strn(itrm:itrm),term) == 0) wordsw = 0  ! Terminator not found

      end

C      subroutine fmain
C      implicit none
C      integer i,j
C      procedure(integer) :: wordsw
C
C      i = wordsw('readsek;ib=3:6',';','jb','',j)
C      print 333, "wordsw('readsek;ib=3:6',';','jb','',j) returned", i, " should return 0"
C  333 format(1x,a,i3,a)
C      i = wordsw('readsek;ib=3:6',';','ib','+-',j)
C      print 333, "wordsw('readsek;ib=3:6',';','ib','+-',j) returned", i, " should return 0"
C      i = wordsw('readsek;ib=3:6',';','ib','+=-',j)
C      print 333, "wordsw('readsek;ib=3:6',';','ib','+=-',j) returned", i, " should return 8"
C      i = wordsw('readsek;ib=3:6',';','ib','',j)
C      print 333, "wordsw('readsek;ib=3:6',';','ib','',j) returned", i, " should return 8"
C      end
