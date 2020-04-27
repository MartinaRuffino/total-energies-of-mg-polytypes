      module m_toksw
C- Module to create token switches, and read contents of tokens
C ----------------------------------------------------------------------
Cr First, a list of tokens must be created, one list each for
Cr any program which uses this module.
Cr Suppose a particular program 'foo' reads tokens
Cr   STRUC_ALAT STRUC_PLAT SITE_ATOM_POS HAM_TOL ...
Cr Further suppose that HAM_TOL is "optional", i.e. no error should occur
Cr should token HAM_TOL be missing from the input.
Cr Create the token list corresponding to 'foo' as follows:
Cr    call nwsadd()  ! starts list for a new program
Cr    call tkadd(" foo::" ) ! Identifies foo as progn associated w/ this list
Cr    call tkadd(" STRUC_PLAT") ! Append this token to the list
Cr    call tkadd(" HAM_TOL~")   ! ~ => token can be optionally read
Cr    call tkadd(" STRUC_PLAT SITE_ATOM_POS") ! lump several tokens in 1 call
Cr    call tkadd(" OPTIONS_SHORBZ!") ! => token not to be read (used by e.g.
Cr                                        lmfa, which starts from lmf structure
Cr                                        and makes changes to it.  See tksw.
Cr    etc
Cr
Cr For examples, see how routine toksw_init in m_rdctrl.f bundles the above
Cr calls for a particular set of programs, e.g. "LMF", "LM", etc.
Cr
Cr After the lists have been created, function tksw(foo,tag) returns
Cr a value 0,1,2 indicating :
Cr   0 tok is sought by foo, with input optional
Cr   1 tok is sought by foo, with input required
Cr   2 tok is is not sought by foo
Cr   3 tok is present but input not sought by foo
Cr
Cr The input file reader contains commands to read the contents of all
Cr tokens any program using this input system will seek.
Cr It calls tksw to determine whether a particular token is be read
Cr by the given calling program.
Cr
Cr This module also contains module m_gtv, the module that read and
Cr parses tokens and tokens' contents.  The entry points are:
Cr   gtv_r8,gtv_r8v,gtv_i4,gtv_i4v,gtv_lg,gtv_char,gtv_none
Cr and are usually invoked by the interface gtv, e.g. this call:
Cr     call gtv('HEADER',tksw('LMF','HEADER'),header,note=
Cr    .  'Contents displayed at beginning of program execution')
Cr reads string 'header'
C ----------------------------------------------------------------------
      implicit none
C ... Module variables
      logical,private:: debug=.false.
      integer,parameter,private :: lenmax=3000,nswmx=25
      character(lenmax),private:: swtok(nswmx+1)
      integer,private ::  nsw=0, iswlast = 0

      contains

      subroutine clear_swtok(debug_in)
C- Resets all tokens in swtok and sets internal debug switche
C ----------------------------------------------------------------------
Ci Inputs
Ci   debug_in : sets internal debug switch for future calls
Co Outputs
Cm Module variables
Cm  swtok : is initialized
Cm  nsw   : is initialized to zero
Cm  debug : is assigned to debug_in
Cr Remarks
Cr
Cu Updates
Cu   16 Jan 09
C ----------------------------------------------------------------------
C ... Passed parameters
      logical:: debug_in
C ... Module variables
      debug = debug_in
      swtok = ' '
      nsw = 0
      end subroutine

      subroutine nswadd()
C- Flag that all tokens have been appended for current main program
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Cm Module variables
Cm  nsw : nsw is incremented, starting a new
Cm      : entry in character array swtok
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none

      nsw = nsw+1
      if (nsw>nswmx) call rx('m_toksw: enlarge nswmx')
      end subroutine

      subroutine tkadd(tok)
C- Append tok to current set in swtok
C ----------------------------------------------------------------------
Ci Inputs
Ci   tok  : Append 'tok' to swtok(nsw)
Ci        : Note: nsw must be set previously by calling nswadd()
Co Outputs
Co  swtok : tok is appended to (module string) swtok(nsw)
Cl Local variables
Cr Remarks
Cr   This routine builds up character array swtok for a particular
Cr   main program.
Cr
Cr   swtok contains a set of lists of tokens, one list
Cr   for each main program, which are read by the main program,  e.g.
Cr     swtok(1) = "LMF:: HAM_TOL~ STRUC_PLAT ..."
Cr     swtok(2) = "LM:: SPEC_ATOM_DV~ HAM_NONCOL~ ..."
Cr   The list for each main program contains a sequence of strings,
Cr   separated by a space, each string corresponding to the full name
Cr   of a token (incl parents).  Optionally the name may be appended
Cr   by a "~" which signifies the token can be optionally read.
Cr
Cr   swtok is built up by the following sequence:
Cr     call nwsadd()  ! starts new list for another program
Cr     call tkadd(" prgn::" )  ! prgn = main program name, e.g. "LMF"
Cr     call tkadd(" tok")      ! tok = token name, e.g. SITE_ATOM_POS
Cr     call tkadd(" tok")      ! etcetera
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*):: tok
C ... Local parameters
      integer:: lensw

      if (nsw > nswmx) call rx('m_toksw: enlarge nswmx')
      lensw = len(trim(swtok(nsw)))
      if ( lensw+30>lenmax ) call rx('m_toksw: enlarge lenmax')
      swtok(nsw) = trim(swtok(nsw))//trim(tok)

C     lenc = max(lenc,lensw)  ! lenc (to check the length of buff)
C     print *, 'toksw for nsw=',nsw,' is now:', trim(swtok(nsw))
      end subroutine

      integer function tksw(prgn,tag)
C- Return 0, 1, 2 or 3 for given program and tag
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgn : program name, to which a family of tags is associated (see Remarks)
Ci   tag  : full name of token (incl parents) e.g. SITE_ATOM_POS
Co Outputs
Ci   tksw : 0 if 'tag'~ is present => input optionally sought by 'prgn'
Ci        : 1 if 'tag'  is present => input is required by 'prgn'
Ci        : 2 if 'tag'  is missing => input not sought by 'prgn'
Ci        : 3 if 'tag'! is present => input not sought by 'prgn'
Ci        : The next two quantities are designed to work with the
Ci        : 'express mode' input style.
Ci        :-1 if 'tag'% is present => input is required by 'prgn'
Ci        :   tksw=-1 similar to tksw=1.  tksw=-1 flags that once
Ci        :   this tag has been parsed, it should be excluded from
Ci        :   subsequent parsing of the same tag (change '%' -> '!')
Ci        :   This occurs when a tag is aliased to another tag, and the
Ci        :   two tags  point to the same quantity.  Converting the
Ci        :   tag to an excluded tag after it has been parse avoids duplication.
Ci        :-2 if 'tag'& is present => input optionally sought by 'prgn'
Ci        :   tksw=-2 is similar to tksw=0, The difference is that -2
Ci        :   provides the same flag function as tksw=-1.
Ci        :-2 if 'tag'# is present => input optionally sought by 'prgn'
Ci        :   tksw=-2 is similar to tksw=0, The difference is that -2
Ci        :   provides the same flag function as tksw=-1.
Cr Remarks
Cr   From a list of tags belonging to 'prgn' (see swtok below), determine
Cr   whether the supplied tag is in the list.
Cr  *If it is not, tksw returns 2, implying tag is not associated with prgn
Cr  *Otherwise:  tksw returns 0, 1, or 3 depending whether the tag stands alone, or
Cr   or is appended by a '~' or by a '!' (see below, and also tksw, above).
Cr
Cr   In more detail, this routine searches character array swtok (which must be
Cr   generated in advance) to determine how a program is to treat data associated with
Cr   the given tag.
Cr
Cr   swtok contains a set of lists of tags (full name of token), one list for each
Cr   kind of program, which identify input   read by the main program,  e.g.
Cr     swtok(1) = "LMF:: HAM_TOL~ STRUC_PLAT ..."
Cr     swtok(2) = "LM:: SPEC_ATOM_DV~ HAM_NONCOL~ ..."
Cr   The list for each main program contains a sequence of strings,
Cr   separated by a space, each string corresponding to the full name
Cr   of a token (including parents).  If the tag corresponding to 'prgn' is present:
Cr   it may be appended by a "~" which signifies the token can be optionally read;
Cr   it may be appended by a "!" which signifies the token should not be read.
Cr   it may be appended by a "%" or '#' which flags that if token is be read
Cr   more than once, it should be ignored after the first parsing.
Cu Updates
Cu   14 Apr 16 Can now return tksw=-1 or -2
Cu   10 Oct 14 Can now return tksw=3
C------------------------------------------------------
      implicit none
      character*(*):: tag,prgn
      integer:: isw,ix,sw,i
      integer, parameter:: npgsz=30

C ... Find which of swtok corresponds to prgn
      do  isw = 1, nsw
        ix = index(swtok(isw)(1:npgsz),trim(prgn)//'::')
        if (ix/=0) then
          if (swtok(isw)(ix-1:ix-1) == ' ') goto 100
        endif
      enddo
      call rx('tksw: no switches for program '//prgn)
  100 continue

      iswlast = isw

C     Quick check : if no token found return
      tksw = 2
      i = index(swtok(isw),' '//trim(tag))
C     Special case tag corresponds to program nam
      if (i == 1) i = index(swtok(isw)(2:),' '//trim(tag))
      if (i == 0) return

C     Case token among required tokens
      if (index(swtok(isw),' '//trim(tag)//' ') /= 0) then
        sw = 1
        if (debug) print *,'tksw: '//trim(tag)//' required'

C     Case token among optional tokens
      elseif (index(swtok(isw),' '//trim(tag)//'~') /= 0) then
        sw = 0
        if (debug) print *,'tksw: '//trim(tag)//' optional'

C     Case token among explicitly excluded tokens
      elseif (index(swtok(isw),' '//trim(tag)//'!') /= 0) then
        sw = 3
        if (debug) print *,'tksw: '//trim(tag)//' excluded'

C     Case token among required tokens, alias category
      elseif (index(swtok(isw),' '//trim(tag)//'%') /= 0) then
        sw = -1
        if (debug) print *,'tksw: '//trim(tag)//' required (*)'

C     Case token among optional tokens, alias category
      elseif (index(swtok(isw),' '//trim(tag)//'&') /= 0) then
        sw = -2
        if (debug) print *,'tksw: '//trim(tag)//' optional (*)'

C     Case token among required tokens, optional in alias category
      elseif (index(swtok(isw),' '//trim(tag)//'#') /= 0) then
        sw = -3
        if (debug) print *,'tksw: '//trim(tag)//' required (*)'

C     It is some other tag
      else
        sw = 2

      endif

      tksw = sw
      end function tksw

      integer function tksw2(prgnam,tag)
C- Return 0, 1, 2 or for tags associated with one of two program names
C ----------------------------------------------------------------------
Ci Inputs
Ci  prgnam: A pair of two program names.
Ci        : If the program names are different, two lists are to be searched.
Ci        : The first takes precedence.
Ci        : Algorithm:
Ci        : Evaluate i = tksw(prgnam(1),tag)
Ci        : If i<>2, tag is present in prgnam(1); return i
Ci        : Or if the two program names are the same, return i
Ci        : Otherwise, return j = tksw(prgnam(2),tag)
Ci  tag  : full name of token (incl parents) e.g. SITE_ATOM_POS
Co Outputs
Ci  tksw2 : 0 if 'tag'~ is present => input optionally sought by 'prgn'
Ci        : 1 if 'tag'  is present => input is required by 'prgn'
Ci        : 2 if 'tag'  is missing => input not sought by 'prgn'
Ci        : 3 if 'tag'! is present => input excluded by 'prgn'
Ci        : The next two quantities are designed to work with the
Ci        : express, or simplified input style.  See tksw() for further description
Ci        :-1 'tag'% is present => input is required by 'prgn'
Ci        :-2 'tag'& is present => input optionally sought by 'prgn'
Co        :-3 'tag'# is present => input is required by 'prgn' but optional for alias
Co        :    (not implemented)
Cr Remarks
Cr   This routine performs the same function as tksw, with the extension that tags
Cr   associated with prgnam(2) will be parses if no tag is associated with prgnam(1)
Cu Updates
Cu   10 Oct 14 Can now return tksw=3
C------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*):: tag,prgnam(2)
C ... Local parameters
      integer:: sw
      logical, parameter :: debug=.false.

      if (debug) print *,'tksw2: search for tag '//trim(tag)//' among progs '//trim(prgnam(1))//' '//trim(prgnam(2))

      sw = tksw(prgnam(1),tag)
      if (sw == 2 .and. prgnam(1) /= prgnam(2)) sw = tksw(prgnam(2),tag)
      tksw2 = sw

      end function tksw2

      integer function tksws(prgnam,tag)
C- Return info concerning whether tag or subtags accessible to express input style
C ----------------------------------------------------------------------
Ci Inputs
Ci  prgnam: A pair of two program names.
Ci        : If the program names are different, two lists are to be searched.
Ci        : The first takes precedence.
Ci  tag   : full name of token (incl parents) e.g. SITE_ATOM_POS
Co Outputs
Co  tksws :  3 tag is excluded ... presumes no subtags will be parsed
Co        : -1 'tag'% is present => input is required by 'prgn'
Co        : -2 'tag'& is present => input optionally sought by 'prgn'
Co        : -3 'tag'# is present => input is required by 'prgn' but optional for alias
Co        :     (not implemented)
Co        : ... If none of the above apply, then:
Co        :  2 No tag or subtag is the type 'tag'% or 'tag'#
Co        : -4 a subtag of 'tag' of the type 'tag_subtag'% or 'tag_subtag'#
Cb Bugs
Cb   It may be that prgnam(1) and prgnam(2) have the same tag
Cb   and the latter is read in 'express' mode while the former is not.
Cb   This routine searchs tags from both prgnam(1) and prgnam(2)
Cb   and icorrectly signals that the tag is available for express mode.
Cb   This bug is obscure enough so it probably doesn't matter.
Cr Remarks
Cr   Return
Cr   associated with prgnam(2) will be parses if no tag is associated with prgnam(1)
Cu Updates
Cu   10 Oct 14 Can now return tksw=3
C------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*):: prgnam(2),tag
C ... Local parameters
      integer:: i,isw,ix,ioff,i1,i2,taglen
      logical, parameter :: debug=.false.
      integer, parameter:: npgsz=30
      character(1) :: c

      if (debug) print *,'tksws: search for tag '//trim(tag)//' among progs '//trim(prgnam(1))//' '//trim(prgnam(2))

      taglen = len_trim(tag)

C     Return cases tag is present and ends in '%' or '&' or. '#' or '!'
      tksws = tksw2(prgnam,tag)
      if (tksws == -1 .or. tksws == -2 .or. tksws == -3 .or. tksws == 3) return

C ... Search for any subtag ending in '%' or '&' or '#'
      tksws = 2 ! This the return value if no or subtag ending in '%' or '&' or '#'
      do  i = 1, 2      ! Search in lists for each name

        if (i == 2 .and. prgnam(2) == prgnam(1)) cycle ! Only one name

        do  isw = 1, nsw  ! Find the list associated with this name
          ix = index(swtok(isw)(1:npgsz),trim(prgnam(i))//'::')
          if (ix/=0) then
            if (swtok(isw)(ix-1:ix-1) == ' ') exit
          endif
          if (isw == nsw) call rx('tksws: no switches for program '//prgnam(i))
        enddo

        ioff = 1
C       Loop until no more tags are found
   10   continue ! Re-entry
        ix = index(swtok(isw)(ioff:),' '//trim(tag))
        if (ix == 0) cycle ! No more tags found for prgnam
        c = swtok(isw)(ioff+ix+taglen:ioff+ix+taglen)
C       if (debug) print *,'tksws: found '//swtok(isw)(ioff+ix:ioff+ix+taglen)

        if (c /= '_') then
          ioff = ioff+ix+taglen
          goto 10
        endif
C       call snot
        i1 = ioff+ix
        call nwordg(swtok(isw),1,' ',1,i1,i2)
        c = swtok(isw)(i2:i2)
        if (c == '%' .or. c == '&' .or. c == '#') then
          tksws = -4
          return
        endif
        ioff = i2+1
        goto 10

      enddo

      end function tksws

      subroutine tkexcludem(tag)
C- Convert a tag in swtok to EXCLUDE
      implicit none
      character*(*):: tag
      integer ix,iend
      logical, parameter :: debug=.false.

      do
        ix = index(swtok(iswlast),' '//trim(tag)//'%')
        if (ix == 0) ix = index(swtok(iswlast),' '//trim(tag)//'&')
        if (ix == 0) ix = index(swtok(iswlast),' '//trim(tag)//'#')
        if (ix == 0) ix = index(swtok(iswlast),' '//trim(tag)//'~')
        if (ix == 0) then
          ix = index(swtok(iswlast),' '//trim(tag)//' ')
          if (ix == 0) exit
          iend = ix+len_trim(tag)+1
          if (swtok(iswlast)(iend+1:iend+1) /= ' ') then ! Insert a space
            swtok(nswmx+1)(iend:) = swtok(iswlast)(iend:)
            swtok(iswlast)(iend+1:) = swtok(nswmx+1)(iend:)
          endif
        endif
        iend = ix+len_trim(tag)+1
        if (debug) print *,'tkexclude: convert tag '//swtok(iswlast)(ix:iend)//' to EXCLUDE'
        swtok(iswlast)(iend:iend) = '!'  ! Mark tag with "exclude" symbol
        exit
      enddo

      end subroutine tkexcludem

      integer function tkswp(prgnam,express,tag)
Co Outputs
Co  tkswp : case express >= 0 :
Co        :  return tksw2(prgnam,tag)
Co        : case express <  0 :
Co        :  return tksws(prgnam,tag)
      implicit none
C ... Passed parameters
      character*(*):: tag,prgnam(2)
      integer:: express
C     procedure(integer) :: tksw2
C ... Local parameters
C     integer:: sw,i,isw,ix,ioff,i1,i2,taglen

      if (express >= 0) then
        tkswp = tksw2(prgnam,tag)
      else
        tkswp = tksws(prgnam,tag)
      endif

      end function tkswp

      end module

      subroutine tkexclude(tag)
C- Convert a tag in swtok to EXCLUDE
C     tkexclude placed outside module to circumvent need for caller to read it.
C     tkexclude merely calls tkexcludem, inside m_toksw
      use m_toksw

      character*(*):: tag

      call tkexcludem(tag)
      end


      subroutine find_region(lev,instr,token,toks,tokt,ntk,itrm,eor,
     .  is,ie)
C- Finds substring corresponding to ntk-th token.
C ----------------------------------------------------------------------
Ci Inputs
Ci   lev   :token level, reserved for acceleration (not used now).
Ci   instr :input string, where token and contents are embedded
Ci   token :token to match in instr; see Remarks
Ci   toks: :(only used if 1000s digit of itrm is set)
Ci         :list of characters that must precede token to satisfy a
Ci         :match.  Eliminates strings matching token embedded within
Ci         :other string.  Example:
Ci         :Example: given token ATOM, and
Ci         :toks='[ ' matches ' ATOM' and '[ATOM] but not 'CATOM'
Ci   tokt: :terminator(s) to token
Ci         :Example: given token ATOM,
Ci         :tokt='=:' matches matches ATOM= and ATOM:
Ci   ntk   :search for ntk-th occurence of token
Ci   itrm  :(itrm=0) start-of-region begins with second character following
Ci         :         token (1st char is token terminator)
Ci         :         end-of-region is delimited by first occurence of string
Ci         :         contained in eor following start-of-region
Ci         :
Ci         :(itrm=1) the first nonblank character after the token and its
Ci         :         terminator must be '['.  Note: '[' may also server as
Ci         :         start-of-region begins with the character following it.
Ci         :         end-of-region is delimited by the matching ']'
Ci         :         Note: '[...]' pairs may be nested; it is an error for
Ci         :         '[' not to have a matching ']'.
Ci         :
Ci         :(itrm=2) If the first nonblank character after the token and its
Ci         :         terminator is '[', follow the syntax of itrm=1.
Ci         :         Otherwise, follow the syntax of itrm=0.
Ci         :(itrm=11)Identical to itrm=1, with the addition that '[' may
Ci         :         serve both as terminator and start-of-region marker
Ci         :(itrm=12)Identical to itrm=2, with the addition that '[' may
Ci         :         serve both as terminator and start-of-region marker
Ci
Ci         :Adding 100 to itrm causes find_region to move
Ci         :start-of-region marker to 1st nonblank character
Ci         :Adding 1000 to itrm turns on the pre-token matching;
Ci         :see toks above
Ci
Ci   eor   :string demarcating end-of-region.  Its use depends on
Ci         :how start-of-region was determined (see itrm)
Ci         :If start-of-region is NOT determined by '[',
Ci         :the first string after start-of-region that matches the
Ci         :contents of eor demarcates end-of-region
Ci         :If start-of-region IS determined by '[', eor is not used
Co Outputs
Co   is    : start-of-region containing token contents, i.e.
Co         : token(is:is) is the first character after the terminator
Co         : If token is not matched, is=-99999
Co         : If token is matched but no '[' follows and itrm=1, is=-99998
Co         : If token is matched and itrm=1, but the matching ']' terminator
Co         :          cannot be found, is=-99997
Co   ie    : end-of-region containing ntk-th occurence of token contents
Co         : This index is EITHER:
Co         : (itrm=0) index to end-of-region.
Co         : (itrm>0) start of ntk+1 th occurence of token.
Co         : In either case, when marker is not found, ie=end of instr
Co         : -99999  token not found
Co         : -99998  missing start-of-region '['
Co         : -99997  missing end-of-region ']'
Cl Local variables
Cl         :
Cr Remarks
Cr   Find pointers is and ie of instr(is:ie) that demarcate token
Cr   contents.  How start-of region and end-of-region are determined depends
Cr   on itrm; see above.
Cr   Examples:
Cr      token        tokt   itrm  Matches
Cr     '@CAT@HAM'    ' '     0    '@CAT@HAM '
Cr     'SIG'         ' =:'   1    'SIG= [' or  'SIG: [' or 'SIG ['
Cr     'ATOM'        '='     2    'ATOM='  or 'ATOM= ['
Cu Updates
Cu   07 Aug 07 Adapted from region_finder (TK)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character(*),intent(in) :: instr,token,toks,tokt,eor
      integer is,ntk,lev,ie,itrm
C ... Local parameters
      integer j,i0,k,lentok,litrm,nnest,isave,i1,iprint,itrm2
      logical :: debug
      character(1) :: cc

      debug = iprint() >= 110
      ie = len(instr)
      lentok = len(token)

C --- Find ntk-th occurence of token ---
C     After this do loop:
C     1.  token+terminator was located
C     2.  If terminator is '[', nnest=1 and litrm=1
      is = 1
      nnest = 0
      itrm2 = mod(itrm,100)
      litrm = itrm2
      do  j = 1, ntk
   10   continue
C       i0 = 0 if token not found
        i0 = index(instr(is:ie),token)
C       No match ; exit
        if ( i0==0 ) then
          is = -99999
          if (debug) write(*,333) token, instr(1:min(20,len(instr)))
  333     format(' find_region: token `',a,
     .      ''' not found in string ',a,' ...')
          return
        endif
        is = is + i0-1 + lentok
C   ... One of toks must precede token; otherwise no match
        if (itrm >= 1000 .and. is-lentok > 1) then
          do  k = 1, len(toks)
            if (instr(is-lentok-1:is-lentok-1) == toks(k:k)) goto 15
          enddo
          goto 10
   15     continue
        endif
C   ... Terminator must follow; otherwise no match
        if (itrm2 > 10) then                ! Special case TOKEN[...
          cc = adjustl(instr(is:))
          if (cc == '[') then
            is = is-1+index(instr(is:ie),cc)
            nnest = 1
            litrm = 1
            goto 20
          endif
        endif
        do  k = 1, len(tokt)
          if (instr(is:is) == tokt(k:k)) goto 20
        enddo
        goto 10
   20   continue
      enddo
      is = is+1

C --- Find is = start-of-region ---
C     If itrm=0, token terminator marks start-of-region
C     In this case, this branch is not executed.
C     If nnest>0, terminator was '[' which marks start-of-region
C     In this case, this branch is not executed.
C     In remaining cases, if the next nonblank character is '[',
C     it marks start-of-region
C     litrm is either 0 or 1 after this branch
      if (itrm2 > 0 .and. nnest == 0) then
        cc = adjustl(instr(is:ie))
        if (itrm2 == 1 .or. itrm2 == 11) then
          if (cc /= '[') then
            if (debug) write(*,"(a)") ' find_region: missing ''['' '//
     .        'after '//instr(is-lentok:is+1)//'...'
            is = -99998
            return
          else
            i0 = index(instr(is:ie),'[')
C           if (i0 == 0) call rx('bug in find_region')
            is = is + i0
          endif
        elseif (itrm2 == 2 .or. itrm2 == 12) then
          if (cc /= '[') then
            litrm = 0
          else
            i0 = index(instr(is:ie),'[')
            is = is + i0
            litrm = 1
          endif
        else
          call rxi('illegal value for itrm:',itrm)
        endif
      endif

C --- Find ie = end-of-region.  Action depends on litrm ---
      if (litrm == 0) then
        i0 = index(instr(is:ie),eor)
C       i0=0 => no eor was found => ie remains (length of instr)
        if (i0 /= 0) ie = is-1 + i0-1
C     Require that end-of-region correspond to ']' matching '['
      else
        isave = is
        nnest = 1
        do  while (nnest > 0)
          i0 = index(instr(is:ie),']')
          if (i0 == 0) then
            if (debug) write(*,"(a)") ' find_region: missing '']'' '//
     .        'after '//instr(isave-lentok:min(isave+10,ie))//'...'
            is = -99997
            return
          endif
          i1 = index(instr(is:ie),'[')
          if (i1 > 0 .and. i1 < i0) then
            is = is+i1
            nnest = nnest+1
          else
            is = is+i0
            nnest = nnest-1
          endif
        enddo
        ie = is-2
        is = isave
      endif

C ... Move start-of-region to first nonblank character
      if (ie > is .and. mod(itrm,1000) >= 100) then
        cc = adjustl(instr(is:ie))
        if (cc /= ' ') then
          i1 = index(instr(is:ie),cc)
C         print *, instr(is:ie)
          is = is+i1-1
C         print *, instr(is:ie)
        endif
      endif

C ... Printout
      if (debug)
     .  write(*,"(' find_region: contents of ', a,' : |',a,'|')")
     .  token, instr(is:ie)
      end subroutine

C      subroutine snot
C      end
