      module bin2a_cmn
        implicit none
        integer :: isw0 = 0
      end module bin2a_cmn

C --- Entry for setting up or determinining defaults ---
      subroutine bin2a0(is)
        use bin2a_cmn
        implicit none
        integer :: is
        if (is >= 0) isw0 = is
        if (is < 0) is = isw0
      end subroutine bin2a0

      subroutine bin2a(fmt,nblk,ndec,res,cast,count,mxlen,outstr,ip)
C- Converts number to ascii format, stripping leading blanks, trailing 0
C ----------------------------------------------------------------------
Ci Inputs:
Ci   fmt: cast=1:   holds the string to be appended to outstr
Ci        cast=0:   not used
Ci        cast=2:   not used, unless:
Ci                  ':n' outputs null string when res=nulli
Ci                  'i#' writes # characters, adding leading zeros if needed
Ci        cast=3,4: syntax X[#][:sw] where X is one of
Ci                  'd' to write in decimal representation
Ci                  'e' to write in exponential format
Ci                  'g' to use the smaller of 'd' and 'e'
Ci                  '(n.m)' fixed format, mimics fortran fmt (fn.m)
Ci                  'D' also mimics fortran fmt (Fn.m)
Ci                      D# => supplies n=#; arg ndec supplies m
Ci                  'F' fixed field, picking between d and e that
Ci                      F# => # is field width
Ci                      generates the most significant digits
Ci                  See Remarks for further description
Ci   nblk:  strip leading blanks, leaving a maximum of nblk
Ci   ndec:  (cast=3,4 only) retain a mininimum ndec digits after decimal
Ci          point, i.e. do not suppress trailing zeros to ndec.
Ci          ndec=0 does nothing.  ndec should not exceed precsn.
Ci          (cast=2 only): ndec specifies a field width
Ci   res:   binary value to be converted into ascii string
Ci   cast:  cast of res: 0=logical, 1=char, 2=int, 3=real, 4=double
Ci   count: res(count) is to be converted.  NB: count=0 for first entry
Ci   mxlen: maximum length of outstr
Cio Inputs/Outputs
Cio  ip:    on input, starting position in outstr for write
Cio         NB: ip=0 points to first character in string
Cio  ip:    on output, position of final character written to outstr
Co  Outputs
Co   outstr:binary res(count) written in ascii form to outstr(ip:..)
Cl Local variables
Cl  isw: (see pretty.f and Examples below)
Cl          1's digit: 0 => ndec and precsn are absolute
Cl                     (i.e. number of places to the right of '.').
Cl                     1 => ndec and precsn = no. significant digits
Cl         10's digit: 1 to suppress possible leading zero
Cl       isw uses the preset default isw0, unless specified
Cl       as described in Remarks below
Cr Remarks
Cr  *The string representation of floating point numbers is generated
Cr   by a "prettified" modification of the fortran write statement
Cr   (pretty.f), which includes suppression of trailing zeros and the
Cr   option to include or suppress the leading zero in decimal
Cr   fractions less than 1.  Floating-point formats include:
Cr     'd[n][:isw]' for decimal representation,
Cr     'e[n][:isw]' for exponential representation,
Cr     'g[n][:isw]' uses the minimum length of 'd' and 'e'
Cr     'D[n][:isw]' simulates the standard fortran format fn.m
Cr                  Here n follows D, ndec the role of m.  Or:
Cr     'Fn'         fixed field, picking between d and e that generates
Cr                  the most significant digits
Cr      (n.m)       also simulates the standard fortran format.
Cr
Cr  *Optional modifier 'n' is a number specifying how many decimals of
Cr   precision (n=6 if not specified). By default, n means:
Cr      for 'd' format, the absolute precision: i.e.
Cr        number of digits after the decimal point
Cr     for 'e' format, the relative precision , i.e.
Cr        number of digits printed
Cr     for 'D' format, it is the field width n in fortran format fn.m
Cr  *Optional modifier isw is a compound of the 1's and 10's digits.
Cr     1's digit isw can overwrite the default meaning of 'n' above.
Cr               isw=0 => n corresponds to absolute precision
Cr               isw=1 => n corresponds to relative precision
Cr       10's digit isw nonzero : suppress leading blanks.
Cr  *Entry bin2a0 allows the user to set the default of isw.
Cr  *Examples:
Cr     call bin2a('d2',1,3,1.234951d0,...)    => 1.23
Cr     call bin2a('d4',1,4,1.234951d0,...)    => 1.2350
Cr     call bin2a('d3:11',1,0,1.2349501d-6,4) => .00000123
Cr     call bin2a('e2',1,3,1.2499d7,...)      => 1.2e7
Cr     call bin2a('e5',1,5,1.2349510d7,...)   => 1.2350e7
Cr     call bin2a('e5:0',1,4,1.2349501d5,...) => 1.234950100e5
Cr     call bin2a('g',1,0,1.23d-5,...)        => 1.23e-5
Cr     call bin2a('g3:10',1,3,1.24996d-5,...) => .000
Cr     call bin2a('g4:10',1,4,1.24996d-5,...) => 1e-5
Cr     call bin2a('f4:10',1,4,1.24996d-5,...) => 1e-5
Cu Updates
Cu   31 Dec 17 Added i# option for integer formats
Cu   25 Sep 14 remove leading '-' sign from ostensibly zero numbers
Cu   02 Aug 07 Added :n outputs null string when res=nulli
C ----------------------------------------------------------------------
      use bin2a_cmn
      implicit none
C ... Passed parameters
      integer nblk,cast,count,ip,mxlen,ndec
CML   character*(*) outstr,fmt
      character*(mxlen) outstr
      character*(*) fmt
      double precision res(0:*)
C ... Local parameters
      logical ld,ls,lf,lnull,llnull,parstr
      integer i,j,k,lsmx,n1,n2,np,precsn,fw,ix(4),iv(4),p,isw,m,ndig,ndige
      parameter (lsmx=80)
      character*20 lfmt*20, fm*20, strn*(lsmx), strn2*(lsmx), ss*(lsmx)
      real rval
      double precision xx
      integer, parameter :: NULLI=-99999
      procedure(integer) :: iprint,a2vec,getdig

C     write(*,"('enter bin2a: cast,fmt=',i4,1x,a$)") cast,fmt

C --- Convert binary to ascii representation (log, int, char) ---
      lnull = .false.
      llnull = .false.
      goto (10,11,12,20,20), cast+1
      call rx('bin2a: bad cast')
   10 call bin2al('(L8)',res,count,strn2)
      goto 15
   11 strn2 = fmt
      goto 15
   12 continue
      call bin2ai('(I16)',res,count,strn2,lnull)
      if (fmt /= ' ' .and. fmt(1:min(len(fmt),2)) /= ':n') then
        if (fmt(1:2) /= '(i' .and. fmt(1:2) /= '(I' .and. fmt(1:1) /= 'i')
     .    call rx('bin2a cannot parse format '//fmt)
        j = 1; if (fmt(1:1) == '(') j = 2
        if (a2vec(fmt,len(fmt),j,2,': )',3,2,1,ix,iv) /= 1)
     .    call rx('bin2a cannot parse format '//fmt)
        ss = adjustl(strn2)
        k = len(trim(ss))
C       Handle overflow
        if (k > iv(1)) then
          strn2 = ' '
          forall (i=1:iv(1)) strn2(i:i) = '*'
          goto 15
        endif
        strn2 = ' '; forall (i=1:iv(1)) strn2(i:i) = '0'
        if (ss(1:1) == '-') then
          strn2(1:1) = '-'
          j = 2
        else
          j = 1
        endif
        m = iv(1)-k
        strn2(j+m:k+m) = ss(j:k)
      endif
      goto 15
C --- copy strn2 to strn with appropriate number of blanks ---
   15 continue
      i = 0
      call skipbl(strn2,lsmx,i)
      strn = ' '
C     If a field width specified, overwrite spillover with '*'
      if (ndec /= 0) then
        call skpblb(strn2,lsmx,j)
        j = j-ndec+1
        if (j > i) then
          strn2(j+1:j+ndec) = '****************'
        endif
        i  = j
      endif
      strn(1+nblk:lsmx) = strn2(i+1:lsmx)
      call skpblb(strn,lsmx,n1)
      n1 = n1+1
      if (lnull .and. fmt /= ' ') then
      i = 0
      if (parstr(fmt,':n',len(fmt)-1,2,'n',i,j)) then
        llnull = .true.
      endif
      endif
      goto 50

C --- Binary->ascii representation, floating-point ---
   20 continue
      if (cast == 3) xx = rval(res,count+1)
      if (cast == 4) xx = res(count)
      lnull = xx == dble(nulli)

C ... Determine appropriate format
      lfmt = fmt
      i = 0
      call skipbl(fmt,len(fmt),i)
      if (i >= len(fmt)) then
        lfmt = 'g'
      else
        lfmt = fmt(i+1:len(fmt))
      endif
      i = 0
      if (parstr(lfmt,':n',len(lfmt)-1,2,'n',i,j)) then
        lfmt(i+1:) = ' '
        llnull = .true.
      endif
C --- Do the conversion, floating point ---
      if (lfmt(1:1) == '(') then
        write(ss,lfmt) xx
        call pretty(ss,nblk,ndec,20,isw0,strn,n1)
      else
        strn  = ' '
        strn2 = ' '
        ld = .false.
        lf = .false.
        j = 0
C   ... i=1 => 'd'  i=2 =>  'e'  i=3 => 'g'
        call chrps2(lfmt,'degDF',5,len(lfmt),j,i)
        if (i <= 0) call rx('bin2a: bad format: '//lfmt)
        if (i == 5) then
          i = 3
          lf = .true.
        elseif (i == 4) then
          i = 1
          ld = .true.
        endif
C   ... Get precsn (or field width for D or F), in iv(1), sw in iv(2)
        j = j+1
        np = a2vec(lfmt,len(lfmt),j,2,': ',2,2,2,ix,iv)
        isw = 1 + isw0
        if (i == 1) isw = 0 + isw0
C   ... Simulated fortran format: precsn dictated by ndec
        fw = 0
        if (lf) then
          if (np <= 0) call rx('bin2a: bad format: '//lfmt)
          fw = iv(1)
        elseif (ld) then
          precsn = ndec
          fw = -1
          if (np >= 1) fw = iv(1)
C   ... if precsn explicit, use it
        elseif (np >= 1) then
          precsn = iv(1)
C   ... This is the default
        else
          precsn = 6
        endif
        if (np >= 2) isw = iv(2)
        if (isw >= 20) isw = mod(isw,10) + isw0
C  21   continue
C   ... p is the exponent
        p = 0
        if (xx /= 0) then
          p = int(dlog10(dabs(xx)))
          if (dabs(xx) < 1) p = p-1
        endif
C   ... fortran 'f' format
        if (i == 1 .or. i == 3) then
C     ... Estimate total width of format statement for fortran write
          if (lf) then
C       ... m is the space consumed by a '-' sign
            m = (1-int(dsign(1d0,xx)))/2
C       ... precsn = # rhs dec = field width - '-' - '.' - (p+1)
            precsn = fw - m - 1 - max(p+1,1)
C       ... Only works on some compilers
C            if (mod(isw,10) /= 0)
C     .      precsn = fw - m - 1 - max(p+1,0)
C       ... ndig = how many nonzero decimals printed
            ndig = max(precsn+p+1,0)
C       ... Exclude 'f' if it doesn't fit
            if (precsn < 0) then
              ndig = -1
C       ... Exclude 'e' if it does, and number isn't small
            else if (p > -2) then
              i = 1
            endif
C       ... Determine how many digits we get from 'e' format
            if (i /= 1) then
              write(ss,'(1pe20.0)') xx
C             print *, ss
C         ... We want at least 1 more digit than f format
              call pretty(ss,0,max(ndig+1,1),max(ndig+1,1),1,strn,j)
C         ... Tack on trailing 'e0' if pretty discarded it
              k = 0
              call chrpos(strn,'e',j,k)
              if (k >= j) j = j+2
C         ... How many decimals for 'e' format
              ndige = max(ndig+1,1) + fw - j
C         ... If pretty suppresses '.', add it back if ndige>1
              if (ndige > 1) then
                k = 0
                call chrpos(strn,'.',j,k)
                if (k >= j) ndige=ndige-1
              endif
C             print *, strn
            else
              ndige = ndig-1
            endif
C       ... Generate string for F format here.
            if (ndig < 0 .and. ndige < 0) then
              strn = ' '
              strn(nblk+1:nblk+fw) = '********************************'
              n1 = fw+nblk
              goto 50
            else if (ndig >= ndige) then
              i = 1
            else
              i = 2
              precsn = ndige
              goto 35
            endif
          elseif (.not. ld .or. (ld .and. fw == -1)) then
            fw = max(p+3,5) + precsn
            if (getdig(isw,0,10) == 1)
     .        fw = max(p+3,3) + max(precsn-p-1,0)
            fw = max(fw,10)
            if (fw > min(lsmx-2,99)) then
              strn = ' '
              strn(nblk+1:nblk+1) = '*'
              n1 = 1+nblk
              goto 35
            endif
          endif
          j = fw
C     ... Insert leading blanks
C         if (lf) then
          if (lf .or. ld) then
            j = j+nblk
          endif
          if (j >= 10) write(fm,'(''(f'',i2,''.'')') j
          if (j < 10) write(fm,'(''( f'',i1,''.'')') j
          k = j
C     ... Number of decimals for fortran write
          j = precsn
          if (.not. (ld .or. lf)) then
            if (getdig(isw,0,10) == 1) j = precsn-p-1
            j = max(j,0)
          endif
C         decimals can't exceed field width - 1
          j = max(min(k-1,j),0)
          if (j >= 10) write(fm(6:8),'(i2,'')'')') j
          if (j < 10) write(fm(6:7),'(i1,'')'')') j
          write(ss,fm) xx
          if (ld .or. lf) then
            if (nblk <= 0) then
            elseif (ss(1:nblk) /= ' ') then
              ss(1:k) =
     .          '*****************************************************'
              ss(1:nblk) = ' '
            endif
            strn = ss
            call skpblb(strn,lsmx,n1)
            n1 = n1+1

            k = 0
            call chrps2(strn,'-.0123456789',12,n1,k,j) ! First numerical character
            j = j-1
            ls = j == 0  ! ls=T if first character is a minus sign
C       ... Remove leading minus if string has no characters [1-9]
            if (ls) then
              m = 0
              call chrps2(strn,'123456789',9,n1,m,p) ! m,p destroyed here
              if (p == 0) then
                ls = .false.
                strn(k+1:k+1) = ' '
                call chrps2(strn,'.0123456789',11,n1,k,j) ! Remake k,j. Do not decrement j this time
              endif
            endif
            if (ls) call chrps2(strn,'.0123456789',11,n1,k,j)
C     ...   Case fraction should have a leading '0'
            if (j == 1 .and. getdig(isw,1,10) == 0) then
              if (ls .and. k > 1) strn(k-1:k) = '-0'
              if (.not. ls .and.k > 0) strn(k:k) = '0'
C     ...   Case fraction should have no leading '0'
            elseif (j == 2 .and. getdig(isw,1,10) /= 0) then
              if (ls) strn(k:k+1) = ' -'
              if (.not. ls)strn(k+1:k+1) = ' '
            endif
          else
!           print *, 'before pretty ...', fm, ss
            k = index(ss,'-')  ! Remove leading '-' if string is ostensibly zero
            if (k > 0) then
              m = 0
              call chrps2(ss,'123456789',9,len(ss),m,p) ! m,p destroyed here
              if (p == 0) ss(k:k) = ' '
            endif
            call pretty(ss,nblk,ndec,precsn,isw,strn,n1) ! Prettify number
          endif
   35     continue
        endif
C    .. fortran 'e' format
        if (i == 2 .or. i == 3) then
          j = p + precsn
          if (getdig(isw,0,10) == 1) j = precsn-1
          if (j > 22) then
            strn2 = ' '
            strn2(nblk+1:nblk+1) = '*'
            n2 = 1+nblk
            goto 45
          endif
          j = min(max(j,0),99)
          if (j >= 10) write(fm,'(''(1pe30.'',i2,'')'')') j
          if (j < 10) write(fm,'(''(1pe30.'',i1,'')'')') j
          write(ss,fm) xx
*         print *, 'before pretty ...', fm, ss
          j = ndec
          if (lf) j = precsn
          call pretty(ss,nblk,j,precsn,isw,strn2,n2)
C     ... Tack on trailing 'e0' if pretty discarded it
          j = 0
          call chrpos(strn2,'e',n2,j)
          if (j >= n2 .and. i == 2) then
            strn2(n2+1:n2+2) = 'e0'
            n2 = n2+2
          endif
C     ... Sometimes the '.' is suppressed; make fw right
          if (lf .and. n2 < fw+nblk) n2 = fw+nblk
   45     continue
        endif
*        if (i == 3)
*     .    print *, n1,n2,'compare |', strn(1:n1), '|', strn2(1:n2), '|'
        if (i == 2 .or. i == 3 .and.
     .    (n2 < n1 .or. strn(nblk+1:nblk+1) == '*')) then
          strn = strn2
          n1 = n2
        endif
      endif

C --- Copy to outstr ---
   50 continue
      n1 = max(n1,0)
      n2 = max(min(n1,mxlen-ip),0)
C     Handle null number: replace ascii string with 'NULL'
      if (lnull .and. llnull) then
        strn(1:n2) = ' '
        i = max(n2-3,1+nblk)
        strn(i:n2) = 'NULL'
      endif
      if (n2 > 0) outstr(ip+1:ip+n2) = strn(1:n2)
      ip = ip+n2
      if (ip == mxlen .and. n2 < n1) outstr(ip:ip) = '|'
      if (iprint() > 120) print '(1x,a,a)', 'bin2a:',outstr(1:ip)

      end
      subroutine bin2al(fmt,res,count,strn)
      character*(*) fmt, strn
      integer count
      logical res(0:*)
      write(strn,fmt) res(count)
      end
      subroutine bin2ai(fmt,res,count,strn,lnull)
C- Conversion of integer to ascii string
C ----------------------------------------------------------------------
Ci Inputs
Cl         :
Ci   fmt   : fortran format
Ci   res   : res(count) is converted
Ci   count : index to res: res(count) is converted
Co Outputs
Co   strn  : ascii representation of integer
Ci   lnull : true if res(count)=nulli, otherwise false
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) fmt, strn
      integer count
      integer res(0:*)
      logical lnull
C ... Local parameters
      integer nulli
      parameter (nulli=-99999)
      write(strn,fmt) res(count)
      lnull = res(count) == nulli
      end
