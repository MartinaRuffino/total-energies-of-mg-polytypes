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
Cu   26 Mar 19 Allow sep to consist of more than one character
Cu   15 Apr 17 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) strn,substr,term,sep
      integer itrm
C ... Local parameters
      integer i,n

      n=0
      if (len(sep) > 0) then
        do while (n < len(sep))
          n = n+1
          i = index(strn,sep(n:n)//substr)
          if (i > 0) exit
        enddo
        n = 1
      else
        i = index(strn,substr)
      endif
      wordsw = i
      itrm = i+len(substr)+n
      if (i == 0 .or. len(term) == 0 .or. itrm > len(strn)) return  ! If no match or no terminator
      if (scan(strn(itrm:itrm),term) == 0) wordsw = 0  ! Terminator not found

      end

      integer function parg2(strn,priorsep,substr,argsep,term,cast,nmin,n,itrm,ivec,dvec)
C- If a tag is embedded in a given string, parse for arguments
C ----------------------------------------------------------------
Ci Inputs
Ci   strn  : string to search
Ci priorsep: separator preceding strn.  Can be empty, in which case it is not used.
Ci   substr: string embedded within strn
Ci   argsep: characters that separate arguments without terminating
Ci   term  : terminator string.  If not zero length, the character immediately
Ci         : following substr must be a character in term.
Ci   cast  : 0=logical, 2=int, 3=real, 4=double
Ci   nmin  : minimum number of arguments to parse to return without error
Ci   n     : number of arguments to parse
Co Outputs
Co   itrm  : index to terminator that ended parsing
Co   ivec  : arguments returned in ivec if cast=2
Co   dvec  : arguments returned in dvec if cast=4
Co   parg  : Case n>0: number of arguments parsed.
Co         : Case n=0: position of substr if found, 0 otherwise
Cr Remarks
Cu Updates
Cu   15 Apr 17 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) strn,substr,priorsep,argsep,term
      integer itrm,cast,n,nmin,ivec(n)
      double precision dvec(n)
C ... Local parameters
      double precision xv(100)
      integer i,k,m,iv(100)
      procedure(integer) :: a2vec,wordsw

      if (n > 100) call rx('increase size of xvec in parg2')

      if (n == 0) then ! substr with no arguments
        i = wordsw(strn,priorsep,substr,term,itrm)
        parg2 = i+1
        return
      endif
      i = wordsw(strn,priorsep,substr,'',itrm) ! substr followed by arguments
      if (i == 0) then
        itrm = 0
        parg2 = 0
        return
      endif

      m = itrm-1
      k = a2vec(strn,len_trim(strn),m,cast,argsep//term,len(argsep)+len(term),len(argsep)+1,n,iv,xv)
      if (k < nmin .and. nmin /= 0) call rx('parg2: failed to parse '//trim(strn(i+1:)))
      itrm = m+1
      parg2 = k

      select case (cast)
        case (0:1,3,5:); call rx('parg not ready for this cast')
        case (2); call icopy(k,xv,1,ivec,1)
        case (4); call dcopy(k,xv,1,dvec,1)
      end select

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
