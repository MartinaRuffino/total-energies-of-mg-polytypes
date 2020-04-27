      module iolib_cb
      implicit none
! Cv   recoff: offset to first char in char array (used by getcat)
! Cv   reclen: record length as determined by file i/o
! Cv   nrecs:  number of records in categ
! Cv   maxlen: maximum size of a token or category label
! Cv   catbeg: offset to first character in char array of current category
! Cv           (set by getcat)
! Cv   catsiz: size of this category
! Cv   subsiz: size of this subcategory
! Cv   iend:   in call to partok, position of last char parsed
! Cv   ichoos,nchoos used by partok for data to be input in one of
! Cv           nchoos possible ways
! Cv   optio:  0 make no attempt to parse token, but print out token
! Cv             that would be sought and associated information.
! Cv           1 attempt to parse token and data
! Cv           2 print out results of token
! this is a messy way to go about but better than the common block nevertheless, especially with the mixup around noerr
      integer, parameter :: catl0=7,recl0=128
      integer :: recoff = 0, reclen = recl0, nrecs, maxlen = catl0, catbeg,
     .         catsiz, subsiz, iend, ichoos = 0, nchoos = 0, optio
      logical noerr

!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio

      end module iolib_cb


      subroutine getcat(recrd,cat,msg,iopt)
C- find a category in a record
C ----------------------------------------------------------------
Ci Inputs
Ci   recrd,cat,term
Ci   lopt: if error, aborts if unmatched category
Ci  .. from the common block:
Ci   optio: 0: show category that would have been sought
Ci          1: attempt to find category
Ci          2: print out category sought
Ci   recoff,reclen,nrecs,maxlen
Co Outputs
Co   catbeg,catsiz,subsiz,iend, noerr (see remarks)
Cr Remarks
Cr   on return, noerr is always true unless a category was actually
Cr   sought (optio=1) and none was found, or if iopt=2.
Cu Updates
Cu   20 Oct 06 MPI: only master node prints messages
C ----------------------------------------------------------------
      use iolib_cb
      implicit none
      integer iopt
      character*1 recrd(0:*),cat*(*),msg*(*)
      integer procid,master,mpipid
      character strn*80, term*1, strn2*80
      integer irecs,i,i1mach
      logical lopt
      logical lsequ
!       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .        iend,ichoos,nchoos,optio
!       logical noerr

!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio


      master = 0
      procid = mpipid(1)

      if (iopt == 2) then
        noerr = .false.
        return
      endif
      lopt = iopt /= 0
      i = len(cat)
      term = cat(i:i)
      noerr = .true.
      goto (1,2,3), 1+mod(optio,10)

C --- print message showing category sought ---
    1 continue
      if (procid == master) then
      strn = ' category  ' // cat
      if (.not. lopt) call awrit0('%a (optional)',strn,len(strn),0)
      strn2 = msg
      if (msg /= ' ') call awrit0('%a:  '//strn2,strn,len(strn),0)
      call awrit0(strn,' ',-len(strn),i1mach(2))
      endif
      return

C --- print out category ---
    3 continue
      if (procid == master) then
      print 11, cat
   11 format(1x,60a)
      endif
      return

C --- find a category ---
    2 continue
      catbeg = recoff
      catsiz = 0
      subsiz = 0
      irecs = nrecs
   10 continue
        irecs = irecs-1
        if (.not. lsequ(recrd(catbeg),cat,maxlen,term,iend)) then
          catbeg = catbeg + reclen
          if (irecs == 0)  then
            noerr = .false.
            if (lopt) then
              call msgio(strn,cat,' ',' ',3,1,.true.,i)
              call cexit(-1,1)
            endif
            return
          endif
          goto 10
        endif
C Found a category.  Next determine the length of the input string
      catsiz = catbeg
   20 continue
        irecs = irecs-1
        catsiz = catsiz + reclen
      if (recrd(catsiz) == term .and. irecs >= 0) goto 20
C Clean up and quit
      catsiz = catsiz - catbeg
      subsiz = catsiz
      end
      logical function scat(ifi,categ,term,lrewnd)
C- Scan file for a category
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit
Ci   categ :seek category 'categ'
Ci   term  :string that terminates category
Ci   lrewnd:true, rewind file before searching; false, do not
Co Outputs
Co   scat  :true if category found, false if not
Co         :file unit is placed after record for which category was found
Cr Remarks
Cr   This is a file version of subroutine getcat
C ----------------------------------------------------------------
      integer ifi
      logical lrewnd
      character categ,term
      character*72 a
C Local variables
      integer catl0,recl0
      parameter (catl0=7,recl0=72)
      integer i
      logical rdstrn,lsequ
      external rdstrn,lsequ
      if (lrewnd) rewind ifi
      scat = .false.
   10 if (.not. rdstrn(ifi,a,recl0,.false.)) return
      if (.not. lsequ(categ,a,catl0,term,i)) goto 10
      scat = .true.
      end
      subroutine subcat(categ,token,term,i)
C- find a subcategory within a category
C ----------------------------------------------------------------
Ci Inputs
Ci   categ,token,term
Ci   i:     offset to first character of category where search begins
Co Outputs
Co   subsiz set to less than range of new token
Cr Remarks
Cr   noerr always returns unchanged.  Does nothing unless optio=1.
C ----------------------------------------------------------------
C passed variables
      use iolib_cb
      implicit none
      integer i
      character*1 categ(0:*),term
      character*(*) token
      logical parstr

C local variables
      integer i2,j

! C Common block
!       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .        iend,ichoos,nchoos,optio
!       logical noerr
!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio

      goto (1,2,1), 1+mod(optio,10)

C --- seek token within a category ---
    2 subsiz = catsiz
      i2 = i
      if (parstr(categ,token,subsiz,maxlen,term,i2,j)) subsiz = i2

    1 return
      end
      subroutine msgio(strng,categ,token,term,m1,m2,lopt,i)
C- Messages for getcat and partok
C ----------------------------------------------------------------
Ci Inputs
Ci   strng,categ,token,term,m1,m2,lopt,i
Co Outputs
Co   string is printed to std out
Co   i: length of printed string
Cr Remarks
Cr   String is composed as: msg(m1) msg(m2)
Cr   If the first character in <categ> is not blank, the string
Cr   'category <categ>' is appended to msg(m1)
Cr   If the first character in <token> is not blank, the string
Cr   '  token <token>' is prepended to msg(m2)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*1 strng(72)
      character*1 categ(*), token(*), term
      integer i,m1,m2
      logical lopt
C Local parameters
      integer imsg(11),k
      character*1 msgs(1)
      character*70 msg2
      equivalence (msg2,msgs)
      data imsg /1,2,17,28,41,52,1,1,1,1,1/
      data msg2
     ./'$read error in $unmatched $  is missing$(optional)$  of cast $'/
      call strcop(strng,msgs(imsg(m1)),72,'$',i)

      if (categ(1) /= ' ') then
        call strcat(strng,100,'$','category: $',100,'$',i)
        call strcat(strng,100,'$',categ,100,' ',i)
      endif
      if (token(1) /= ' ') then
        strng(i) = '$'
        call strcat(strng,100,'$','  token  $',100,'$',i)
        call strcat(strng,100,'$',token,100,term,i)
      endif
      call strcop(strng(i+1),msgs(imsg(m2)),72,'$',k)
      i = i+k-1
      if (lopt) print 10, (strng(k), k=1,i)
   10 format(' ',72a1)
      return
      end
      integer function partok(categ,token,term,result,chr,
     .                  count,cast,istart,iopt)
C- Reads a vector of numbers found after a token
C ----------------------------------------------------------------
Ci Inputs
Ci   categ: string to parse for token.
Ci          Parsing is between istart and internal variable subsiz.
Ci   token: token, eg 'nabc='
Ci   term:  term(1:1) = last char of token.
Ci          If len(term)=1, expressions are terminated by spaces.
Ci          If len(term)>1, partok will use any of the term(2:*)
Ci          the separator in parsing the list of elements; see Remarks
Ci   result:(for character cast only, result is integer) max str. size
Ci          26 May 1997: sign result <0: read into all characters
Ci          to abs(result), not exceeding catsiz.
Ci   chr:   (all casts except char case) description of token's function
Ci   count: maximum number of elements to read (if count passed as <0
Ci          count<0 => partok will read exactly -count elements
Ci                     or abort w/ error)
Ci          count=0 => see Outputs, result
Ci   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
Ci   istart:offset to first character of category where search begins
Ci   iopt   ones digit:
Ci          0  token is optional
Ci          1  token is required
Ci          2, do nothing but set noerr to F, return partok=-1
Ci          hundreds digit (for character tokens, nterm>1)
Ci          1: skip over any blanks following a terminator.
Ci          2: Token consists of all characters < min(result,subcat)
Ci   ... the following are in the common block
Ci   optio: ones digit
Ci          0: show token that would have been sought
Ci          1: attempt to find token and read contents
Ci          2: print out token and contents
Ci          tens digit
Ci          0: require ' ' preceding token
Ci          1: allow any of term to precede token
Ci          2: allow any of term to terminate token argument, in addition
Ci             to being separators of element list
Ci          3: combination of 1 and 2
Co Outputs
Co   result:array into which elements are read.
Co          if count=0 partok makes no attempt to parse for an element.
Co          if also cast=0 partok returns noerr in result
Co   chr:   result for case cast=1
Co   noerr: see remarks
Co    iend: offset to first character past match.
Co  partok: number of elements actually read (-1 if no attempt to read
Co          an element -- if optio is 0 or 2 or iopt is 2.)
Cr Remarks
Cr   on return, noerr is always false unless a token was actually
Cr   sought (optio=1) and found (note that this differs from getcat).
Cr
Cr   ichoos is automatically incremented on each call if nchoos /= 0
Cr   ichoos and nchoos are automatically reset to 0 if a match is found
Cr   or if ichoos=nchoos
Cr
Cr   Caution when using multiple-terminators ie len(term)>0:
Cr   partok attempts to parse the expression for each term(2:len(term)),
Cr   and will use any one that succeeds.  To avoid ambiguity, DO NOT use
Cr   terminators that may be part of an allowed arithmetic expression.
C ----------------------------------------------------------------
      use iolib_cb
      implicit none
      integer cast
      character*1 categ(0:*),term*(*)
      character*(*) chr,token
      integer count,result(*),istart,iopt
      logical parstr,a2bin,errflg,lopt
C Local variables
      integer procid,master,mpipid
      logical ldum,last
      integer j,j0,k,idum,i1mach,nterm,iterm,ltoken,rcln0,ilen,ival,is,
     .  i,m
      parameter (rcln0=128)
      character*1 strn(rcln0),fm*(rcln0),ss*(rcln0),ct*1,ss2*256
      equivalence (ss, strn)
      character(len=8), parameter :: nmcast(0:4)
     .  = ['logical ','char    ','integer ','real    ','double  ']
      character, parameter :: cc(5) = ['l',' ','i','r','g']
      character option*11,lterm*80

! C common block
!       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .        iend,ichoos,nchoos,optio
!       logical noerr
!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio

!       data nmcast /'logical','char','integer','real','double'/
!       data cc /'l',' ','i','r','g'/

      master = 0
      procid = mpipid(1)

      if (mod(iopt,10) == 2) then
        noerr = .false.
        partok = -1
        return
      endif
      call sanrg(.true.,mod(iopt,10),0,2,'partok','iopt')

      lopt = mod(iopt,10) /= 0
      errflg = count < 0
      idum = -1
      noerr = .false.
      if (nchoos /= 0) ichoos = ichoos+1
      goto (1,2,3), 1+mod(optio,10)

C --- Print message showing token sought, its cast and length ---
    1 continue
      if (procid == master) then
      option = ' (optional)'
      if (lopt) option = ' '
      if (nchoos /= 0 .and. ichoos /= nchoos) option = '    --- OR:'
      ss = '   token  ' // token
      if (count /= 0) then
        call awrit0('%a  of cast '//nmcast(cast),ss,len(ss),0)
        if (count == -99 .or. count == 99) then
           call awrit1('%a and length *',ss,len(ss),0,iabs(count))
         elseif (count > 1) then
           call awrit1('%a and length <= %i',ss,len(ss),0,iabs(count))
         elseif (count < -1) then
           call awrit1('%a and length %i',ss,len(ss),0,iabs(count))
         endif
      endif
      call awrit0('%a'//option,ss,len(ss),0)
      if (chr /= ' ' .and. cast /= 1) call awrit0('%a :',ss,len(ss),0)
      call awrit0('%a',ss,len(ss),-i1mach(2))
      if (chr /= ' ' .and. cast /= 1) then
        ss2 = '          '//chr//'%a'
        call awrit0(ss2,' ',-len(ss2),i1mach(2))
      endif
      endif
      goto 42

C --- Print out the contents of token ---
    3 if (nchoos == 0 .or. ichoos == 1 .and. count /= 0) then
      if (procid == master) then
        fm = ' '
        call strcop(fm(5:),token,len(token),term(1:1),j)
        if (cast == 1) then
          ss = chr
          call awrit0('%a '//ss,fm,-len(fm),-i1mach(2))
        else
          call awrit2('%a %n:1'//cc(cast+1),fm,len(fm),-i1mach(2),
     .      iabs(count),result)
        endif
      endif
      endif
      goto 42

C --- Seek token within a category ---
    2 iend = istart
      idum = 0

C ... Check for multiple terminators; flag with nterm>1
C ... Copy to lterm to get around SGI compiler bug
      nterm = len(term)
      lterm = term
      ltoken = min(len(token),maxlen)
      k = 1
      if (mod(optio/10,10) >= 2) k = nterm
C ... Find the token within categ(iend:subsiz); ldum=T if found
   21 continue
      do  28  iterm = 1, k
        ct = lterm(iterm:iterm)
        j0 = iend
        ldum = parstr(categ,token,subsiz,ltoken,ct,j0,j)
        if (ldum) goto 29
   28 continue
   29 iend = j0
C ... Case require ' ' preceding token
      if (mod(optio/10,2) == 0) then
        if (ldum .and. iend /= 0 .and. categ(iend-1) /= ' ') then
          iend = j
          goto 21
        endif
C ... Case require one of terminator to precede token
      elseif (mod(optio/10,2) == 1 .and. ldum .and. iend /= 0) then
        do  25  iterm = 2, nterm
          ct = lterm(iterm:iterm)
          if (categ(iend-1) == lterm(iterm:iterm)) goto 26
   25   continue
        iend = j
        goto 21
   26   continue
      endif

C ... Case count and cast both 0: put result of parse into result
      if (count == 0 .and. cast == 0) call lvset(result,1,1,ldum)
      if (ldum) then
        noerr = .true.
      else
        if (ichoos < nchoos .or. .not. lopt) goto 40
        call msgio(strn,categ,token,lterm(1:1),2,4,.true.,j)
        call cexit(-1,1)
      endif

C --- Parse for vector of elements ---
      idum = -1
   10 idum = idum+1
        if (idum == iabs(count)) goto 40
        call skipbl(categ,subsiz,j)
        j0 = j
        if (j < subsiz) then
C     ... Looking for a character string
          if (cast == 1) then
            chr = ' '
            ct = ' '
            is = 1
            ilen = ival(result,1)
C       ... Copy anything in quotation marks
            if (categ(j) == '"' .or. categ(j) == '''') then
              ct = categ(j)
              j = j+1
            endif
C       ... Copy up to min(result,end-of-category)
            if (mod(iopt/100,10) == 2 .and. ct == ' ') then
              call strncp(chr,categ,1,j+1,min(subsiz,ilen)-j)
              goto 10
            elseif (ct /= ' ' .or. nterm == 1) then
              call strcop(chr(is:),categ(j),min(ilen,subsiz-j),ct,k)
              if (chr(k:k) /= ct .and. k == subsiz-j) call rxs4(
     .          'partok token ',token,' missing quotation "',chr(1:k),
     .          '%a " ...')
              if (ct /= ' ') chr(k:k) = ' '
              is = k
              j = j + k
              goto 10
C       ... Copy up to first, or last terminator in token
            else
              i = j
   23         continue
              call chrps2(categ,term(2:),len(term)-1,ilen,j,k)
              last = .true.
              if (mod(iopt/100,10) == 1) last = k == len(term)-1
              j = j+1
C         ... Compress whitespace into one character
              m = i
              call skipbl(categ,ilen,i)
              if (i > m) i = i-1
              call strncp(chr,categ,is,i+1,j-i)
              is = is+j-i
              i = j
C             print *, chr(1:is)
              call skipbl(categ,ilen,j)
              if (j < ilen .and. .not. last) goto 23
            endif
            goto 10
          else
            if (nterm == 1) then
              if (a2bin(categ,result,cast,idum,' ',j,-1)) goto 10
            else
              k = j
              do  12  iterm = 2, nterm
                ct = lterm(iterm:iterm)
                j = k
                if (a2bin(categ,result,cast,idum,ct,j,-1)) goto 10
   12         continue
            endif
          endif
        endif

      if (idum < iabs(count) .and. errflg) then
        call msgio(strn,categ,token,lterm(1:1),2,1,.true.,k)
        print '('' parsing '',99a1)', (categ(k), k=j0,min(subsiz,j0+40))
        print 22, iabs(count), idum
   22   format(' sought',i3,' elements but found only',i3)
        call cexit(-1,1)
      endif

   40 iend = j
   42 if (noerr .or. optio == 2 .or. ichoos == nchoos) then
        nchoos = 0
        ichoos = 0
      endif

      partok = idum
      end
      subroutine partk0(recoff_a,reclen_a,nrecs_a,
     .  catlen_a,catbeg_a,catsiz_a,subsiz_a,optio_a,lget)
C- Initializes internal variables for partok
C ----------------------------------------------------------------
      use iolib_cb
      implicit none
      logical lget
      integer recoff_a,reclen_a,nrecs_a,catlen_a,catbeg_a,catsiz_a,subsiz_a,optio_a
      integer recl,catl,subs
C --- Defaults ---
      recl = reclen_a
      catl = catlen_a
      subs = subsiz_a
!       if (.not. lget) then
!         if (reclen < 0) call stlibv(2,recl,.true.)
!         if (catlen < 0) call stlibv(4,catl,.true.)
!         if (subsiz < 0) subs = catsiz
!         call stlibv(9,0,.false.)
!         call stlibv(10,0,.false.)
!       endif
!       call stlibv(1,recoff,lget)
!       call stlibv(2,recl,lget)
!       call stlibv(3,nrecs,lget)
!       call stlibv(4,catl,lget)
!       call stlibv(5,catbeg,lget)
!       call stlibv(6,catsiz,lget)
!       call stlibv(7,subs,lget)
!       call stlibv(12,optio,.false.)

      if (.not. lget) then
        if (reclen_a < 0) recl = reclen
        if (catlen_a < 0) catl = maxlen
        if (subsiz_a < 0) subs = catsiz_a
        ichoos = 0
        nchoos = 0
      endif

      if (lget) then
         recoff_a = recoff
         recl     = reclen
         nrecs_a  = nrecs
         catl     = maxlen
         catbeg_a = catbeg
         catsiz_a = catsiz
         subs     = subsiz
      else
         recoff = recoff_a
         reclen = recl
         nrecs  = nrecs_a
         maxlen = catl
         catbeg = catbeg_a
         catsiz = catsiz_a
         subsiz = subs
      end if

      optio = optio_a
      end
!       subroutine stlibv(i,ival,lget)
! C- Sets or retrieves an iolib common block variable
! C ----------------------------------------------------------------
!       use iolib_cb
!       implicit none
!       integer i,ival
!       logical lget
!
! ! C common block
! !       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
! !      .        iend,ichoos,nchoos,optio
! !       logical noerr
! !       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
! !      .               iend,ichoos,nchoos,noerr,optio
!       integer rectab(1)
!       equivalence (rectab,recoff)
!
!       if (lget) then
!         if (i == 11) then
!           ival = 0
!           if (noerr) ival = 1
!         else
!           ival = rectab(i)
!         endif
!       else
!         if (i == 11) then
!           noerr = ival /= 0
!         else
!           rectab(i) = ival
!         endif
!
!       endif
!       end
      integer function recln(i)
C- Returns (and optionally sets) the record length of ascii input files
C ----------------------------------------------------------------
      use iolib_cb
      implicit none
      integer i
! C For io routines ...
!       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .        iend,ichoos,nchoos,optio
!       logical noerr
!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio
      if (i > 0) reclen = i
      recln = reclen
      end
! C --- For io routines ---
!       block data drdfil
! C ----------------------------------------------------------------
!       implicit none
!       integer catl0,recl0
!       parameter (catl0=7,recl0=128)
! Cv Internal variables (iolib)
! Cv   recoff: offset to first char in char array (used by getcat)
! Cv   reclen: record length as determined by file i/o
! Cv   nrecs:  number of records in categ
! Cv   maxlen: maximum size of a token or category label
! Cv   catbeg: offset to first character in char array of current category
! Cv           (set by getcat)
! Cv   catsiz: size of this category
! Cv   subsiz: size of this subcategory
! Cv   iend:   in call to partok, position of last char parsed
! Cv   ichoos,nchoos used by partok for data to be input in one of
! Cv           nchoos possible ways
! Cv   optio:  0 make no attempt to parse token, but print out token
! Cv             that would be sought and associated information.
! Cv           1 attempt to parse token and data
! Cv           2 print out results of token
!
!       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .        iend,ichoos,nchoos,optio
!       logical noerr
!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio
!
!       data recoff /0/, reclen /recl0/, maxlen /catl0/,
!      .     ichoos /0/, nchoos /0/
!       end
