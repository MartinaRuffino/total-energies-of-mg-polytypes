C#define FORMATTED_READ
      logical function a2bin(instr,res,cast,count,term,j,jmaxi)
C- Convert ASCII to logical, integer, real, double
C ----------------------------------------------------------------
Ci Inputs
Ci   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
Ci   j:     offset to first character in string to read
Ci   term:  character that terminates the expression
Ci   jmaxi: j is not to exceed jmax
Ci          Also, jmaxi should not exceed (length-of-instr)-1
Ci   count'th element of array res is set to converted value
Co Outputs
Co   j is put past last character read
Cr Remarks
Cr   An ASCII string is converted to a double which is then cast to
Cr   the appropriate type and assigned to the count'th element of res.
Cr   A2BIN conforms to the 1977 ANSI FORTRAN standard.  Without using
Cr   recursion an ASCII string is interpreted with all normal operator
Cr   precidence and associativity through the use of a stack for
Cr   numbers and a stack for deferred operations.  The array OPRULE is
Cr   used to say when operations are deferred and when they are
Cr   executed.  The rules of precidence and associativity for all
Cr   operators are derived solely from the values in OPRULE.  Function
Cr   calls (sin, log, sqrt, ...) are parsed as unary operators.
Cr   OPRULE is initialized in a block data program and its contents
Cr   are this:
Cr        Current operator       Current operator is the most recently
Cr                             parsed operator from the ASCII string.
Cr T    ~ ^ * + < & | ? : ( )  Top operator is the operator on the top
Cr o  ~ F T T T T T T T T F T  of the operator stack.  If OPRULE is
Cr p  ^ F F T T T T T T T F T  .true. for a situation then the top
Cr    * F F T T T T T T T F T  operator is popped and executed.  This
Cr o  + F F F T T T T T T F T  continues until there are no more
Cr p  < F F F F T T T T T F T  operators or until OPRULE yields
Cr e  & F F F F F T T T T F T  .false., then the current operator is
Cr r  | F F F F F F T T T F T  pushed onto the operator stack.
Cr a  ? F F F F F F F F F F T  Parentheses are treated as operators.
Cr t  : F F F F F F F T T F T    A2BIN essentially does an on-the-fly
Cr o  ( F F F F F F F F F F F  conversion from infix to postfix and
Cr r  ) T T T T T T T T T T T  valuates the postfix expression as it is
Cr                             converted.
Cr
Cr   The operators are in order of precedence with associativity:
Cr 1> - (arithmetic negative), ~ (logical negative, .NOT.), and
Cr     functions abs(), exp(), log(), sin(), asin(), sinh(), cos(),
Cr      acos(), cosh(), tan(), atan(), tanh(), flor(), ceil(), erfc(),
Cr      and sqrt() (flor() rounds to the next lowest integer;
Cr      ceil() rounds up, erfc() is the the error function complement).
Cr 2> ^ (exponentiation)
Cr 3< * (times), / (divide), % (modulus)
Cr 4< + (add), - (subtract)
Cr 5< < ( < ); > ( > ); = ( == ); <> ( /= ); <= ( <= );  >= ( >= )
Cr 6< & (.and.)
Cr 7< | (.or.)
Cr 8&9  ?: conditional operators
Cr 10&11 parentheses
Cr
Cr   There is a granularity to the comparison operators (<, >, =, ...)
Cr   that is given by the variable MACHEP.  E.g., A=B if A and B are
Cr   within MACHEP of each other.  MACHEP is a regular variable; if
Cr   not set, it takes a default value of 1E-10.
Cr
Cr   The conditional operators work as follows, the expression
Cr   <logical>?<exp1>:<exp2> has the value <exp1> if <logical> is
Cr   .true.  or <exp2> if <logical> is .false.  If you put conditional
Cr   expressions inside any of the three expressions then A2BIN should
Cr   parse it in the common sense manner but parentheses should be
Cr   used just to be on the safe side.
Cr
Cr   The variables PTKTYP and PREVOP store the most recent token type
Cr   and operator. This is used by GETTOK to distinguish unary
Cr   arithmetic negation from binary subtraction, since they both use
Cr   the same operator.
Cr
Cr   PARENC is used to determine where the expression ends when the
Cr   terminator character is also an operator.
Cu Updates
Cu    3 Jun 10 New function besj0
Cu   09 Oct 09 a2bin can handle leading '+'
Cu   19 Dec 02 Added macro handling
C --------------------------------------------------------------------
      implicit none
      character*1 instr(0:*),term
      integer cast,count,j,jmaxi
      double precision res(0:*)
C ... Local parameters
      integer namlen,parenc
      integer opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok,maxsiz
      double precision numtok,numstk(0:32),machep
      parameter (maxsiz=72,namlen=40)
      double precision dum
      character*(namlen) vartok,macstr*256,strn*256
      logical   oprule(0:10,0:10)
      common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp,
     .                prevop,opstk,optok,parenc,oprule,vartok
      logical   gettok,fstpas,lerr
      integer   tj,ival,jmax,tm,k,t0,jmaxm
      save fstpas
      data fstpas /.true./

C --- Initialization: stacks empty, '-' must be unary negation ---
      a2bin= .false.
      tm = -1
      tj = j
      jmax = jmaxi
      if (jmax < 0) jmax = 999999
      numnum = 0
      opnum = 0
      parenc = 0
      machep = -1D0
      call numsyv(ival)
      if (ival == 0 .or. fstpas) then
        call getsyv('t',numtok,ival)
        if (ival == 0) call addsyv('t',1D0,ival)
        call getsyv('f',numtok,ival)
        if (ival == 0) call addsyv('f',0D0,ival)
        call getsyv('pi',numtok,ival)
        if (ival == 0) call addsyv('pi',4*datan(1d0),ival)
        fstpas = .false.
      endif
      toktyp = 3
      optok = 9
      call skipbl(instr,jmax,tj)
      if (instr(tj) == '+' .or. instr(tj) == ' ') tj = tj+1
C     No string: exit with error
      if (tj > jmax) return

C --- Get next token (variable, number, or operator) from ASCII string
   10 continue
      if (tm < 0) then
        t0 = tj
        if (.not. gettok(instr,tj,term,jmax)) goto 20
      else
        t0 = tm
        if (.not. gettok(macstr,tm,' ',jmaxm)) goto 20
C       end of macro expansion
        if (tm >= jmaxm) tm = -1
      endif
      goto (30,40,50,50,60),toktyp
      return
C ... handle variable token: get value and treat as number token
   30 call getsyv(vartok,numtok,ival)
      if (ival == 0) return
C ... handle number token: put it on the number stack
   40 numnum = numnum+1
      if (numnum > 32) call rx('A2BIN: expression too complex.')
      numstk(numnum) = numtok
!     if (tj == jmax) goto 20 ! Number ended with char parsed
      goto 10
C ... handle operator token: pop a few ops and push current op
   50 continue
C ... handle vector token by assigning ival to # stack
      if (toktyp == 4) then
        call getsvv(vartok,ival,0,1,1,dum)
        if (ival == 0) return
        numnum = numnum+1
        if (numnum > 32) call rx('A2BIN: expression too complex.')
        numstk(numnum) = ival
      endif
      if ((opnum == 0) .or. (.not. oprule(mod(opstk(opnum),16),
     .                                      mod(optok,16)))) goto 56
      opnum = opnum-1
C ... doop() expects the old op to already be popped, hence opnum+1
      call doop(opstk(opnum+1),lerr)
      if (lerr) return
      goto 50
   56 opnum = opnum+1
      if (opnum > 32) call rx('A2BIN: expression too complex.')
      if (optok == 9) parenc = parenc+1
      if (optok == 10) parenc = parenc-1
      opstk(opnum) = optok
      goto 10
C ... handle macro token
   60 continue
      do  k = t0, jmax
        strn(k-t0+1:k-t0+1) = instr(k)
C       print *, k, strn(1:40)
        if (instr(k) == ')') then
          tj = k+1
          strn(k-t0+2:k-t0+2) = ' '
          goto 61
        endif
      enddo
   61 continue
      k = 1
      call macevl(strn,macstr,k)
      if (k < 0) call rx('a2bin: could not parse macro')
      call word(macstr,1,k,jmaxm)
      tm = 0
      goto 10

C --- end expression, pop remaining operators ---
   20 if (opnum == 0) goto 25
      opnum = opnum-1
      call doop(opstk(opnum+1),lerr)
      if (lerr) return
      goto 20
C --- If expression passes the last few syntax checks, return it
   25 continue
C     if (numnum == 0) return
      if (numnum /= 1) return
C ... j,res are not modified unless there is a valid expression
      a2bin= .true.
      j = tj
      call setans(res,res,res,res,cast,count,numstk(1))
      end

      subroutine doop(op,lerr)
C- Perform the passed operation on the top 1, 2, or 3 numbers of the
C- number stack.
C -------------------------------------------------------------------
Ci Inputs
Ci   op:   operator code
Cr Remarks
Cr   lowest 4 bits of op are the operator class, next 28 bits are
Cr   the operator number.  The '?' and '(' operators do nothing
Cr   and should not be passed to DOOP in a syntactically correct
Cr   expression.  All the ')' operator does is pop the associated
Cr   '(' off the operator stack.
C -------------------------------------------------------------------
      implicit none
      integer   op
      integer   namlen,parenc
      integer   opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok
      double precision numtok,numstk(0:32),machep
      logical   oprule(0:10,0:10),lerr
      parameter (namlen=40)
      character*(namlen) vartok
      common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp,
     .                prevop,opstk,optok,parenc,oprule,vartok
      integer   t1,t2,ival,tvec
      double precision  v0,v1,log2d,derfc,dbesj0
      real ran1
      parameter (tvec=25)
C#ifdefC HANSMR
C      integer ixx,l
C      double precision e,rsm,xi(0:10),phi(0:10)
C#endif

C --- Get operator class and number ---
      lerr = .false.
      t1 = op/16
      t2 = mod(op,16)
C ... Make sure that there enough numbers on the number stack to do
C ... the operation
      if ((t2 == 0).and.(numnum == 0)) goto 10900
      if ((t2 == 8).and.(numnum < 3)) goto 10900
      if ((t2 > 0) .and. (t2 < 7) .and. (numnum < 2)) goto 10900
C ... Load the top 2 numbers into v0 and v1 to speed things up
      v0 = numstk(numnum)
      if (numnum /= 1) v1 = numstk(numnum-1)
C ... Except for unary ops and ')', pop v0 off the number stack
      if ((t2 /= 0) .and. (op /= 10)) numnum=numnum-1
C ... Array element is really a binop, so decrement numnum
      if ((t2 == 0) .and. (t1 == tvec)) numnum=numnum-1

C --- Go to appropriate operator class code
      goto (10000,10100,10200,10300,10400,10500,10600,10700,
     .          10800,10900,11000),t2+1

C --- Go to the appropriate unary op/function code
10000 goto (100,101,102,103,104,105,106,107,108,109,110,111,112,113,
     .          114,115,116,117,118,119,120,121,122,123,124,125),t1+1
C ... Arithmetic negation
  100 numstk(numnum) = -v0
      return
C ... Logical not
  101 numstk(numnum) = log2d(v0 == 0D0)
      return
C ... Functions
  102 numstk(numnum) = dabs(v0)
      return
  103 numstk(numnum) = dexp(v0)
      return
  104 numstk(numnum) = dlog(v0)
      return
  105 numstk(numnum) = dsin(v0)
      return
  106 numstk(numnum) = dasin(v0)
      return
  107 numstk(numnum) = dsinh(v0)
      return
  108 numstk(numnum) = dcos(v0)
      return
  109 numstk(numnum) = dacos(v0)
      return
  110 numstk(numnum) = dcosh(v0)
      return
  111 numstk(numnum) = dtan(v0)
      return
  112 numstk(numnum) = datan(v0)
      return
  113 numstk(numnum) = dtanh(v0)
      return
C ... flor()
  114 if (v0 /= int(v0)) then
        if (v0 < 0D0) then
          numstk(numnum) = int(v0)-1D0
        else
          numstk(numnum) = int(v0)
        endif
      endif
      return
C ... ceil()
  115 if (v0 /= int(v0)) then
        if (v0 > 0D0) then
          numstk(numnum) = int(v0)+1D0
        else
          numstk(numnum) = int(v0)
        endif
      endif
      return
C ... erfc()
  116 numstk(numnum) = derfc(v0)
      return
  117 continue
      if (v0 < 0D0) call rx('A2BIN: out of domain in SQRT')
      numstk(numnum) = sqrt(v0)
      return
C ... ran()  (to set seed, call ran(#), with # nonzero integer)
  118 continue
      if (v0 /= 0d0) call ran1in(nint(v0))
      numstk(numnum) = ran1()
      return

C ... functions not yet defined
  119 continue
      if (v0 < 0D0) call rx('A2BIN: out of domain in besj0')
      numstk(numnum) = dbesj0(v0)
      return

  120 continue
  121 continue
  122 continue
  123 continue
  124 continue
      call rx('doop: unknown function')

C ... Return vec[v1], element v0
  125 continue
C      ii1 = v1+dble(1)/2
C      ii2 = v0+dble(1)/2
C      call getsvv(' ',ii1,1,ii2,ii2,numstk(numnum))
      if (numnum < 1) then
        lerr = .true.
        return
      endif
      call getsvv(' ',nint(v1),1,nint(v0),nint(v0),numstk(numnum))
      return

C#ifdefC HANSMR
C  119 continue
C      if (v0 < 0D0) stop 'A2BIN: out of domain in HNSM.'
C      rsm = 0
C      call getsyv('l',rsm,ixx)
C      l = nint(rsm)
C      rsm = 0d0
C      call getsyv('rsm',rsm,ixx)
C      e = 0d0
C      call getsyv('e',e,ixx)
C      if (rsm == 0d0) then
C        call bessl(e*v0*v0,l,phi,xi)
C        xi(l) = xi(l)/v0**(l+l+1)
C      else
C        call hansmr(v0,e,1/rsm,xi,l)
C      endif
C      numstk(numnum) = xi(l)
C      return
C#endif

C --- Binary operators ---
C ... exponentiation
10100 numstk(numnum) = v1**v0
      return
C ... *, /, and %
10200 goto (200,201,202),t1+1
  200 numstk(numnum) = v1*v0
      return
  201 numstk(numnum) = v1/v0
      return
  202 numstk(numnum) = mod(v1,v0)
      return
C ... + and - (binary subtraction)
10300 continue
      if (t1 == 0) then
        numstk(numnum) = v1+v0
      else
        numstk(numnum) = v1-v0
      endif
      return
C ... Conditional operators
10400 if (machep < 0D0) then
        call getsyv('machep',machep,ival)
        if (ival == 0) machep = 1d-10
      endif
      goto (400,401,402,403,404,405),t1+1
  400 numstk(numnum) = log2d(v1 < v0)
      if (numstk(numnum) /= 0D0) goto 405
      return
  401 numstk(numnum) = log2d(v1 > v0)
      if (numstk(numnum) /= 0D0) goto 405
      return
  402 numstk(numnum) = log2d(abs(v1-v0) <= machep)
      return
  403 numstk(numnum) = log2d(v1 < v0)
      if (numstk(numnum) == 0D0) goto 402
      return
  404 numstk(numnum) = log2d(v1 > v0)
      if (numstk(numnum) == 0D0) goto 402
      return
  405 numstk(numnum) = log2d(abs(v1-v0) > machep)
      return
C ... Logical and
10500 numstk(numnum) = v1*v0
      return
C ... Logical or
10600 numstk(numnum) = abs(v0)+abs(v1)
      return
C ... '?' or '(' mean syntax error in expression
10700 continue
10900 continue
      call rx('A2BIN: expression syntax.')
C ... ')', pop the associated '(' and do nothing else
11000 continue
      lerr = ((opnum == 0) .or. (opstk(opnum) /= 9))
      opnum = opnum-1
      return
C ... Conditional expr: pop v1 and load either v0 or v1 onto the stack
10800 numnum = numnum-1
      if (numstk(numnum) == 0D0) then
        numstk(numnum) = v0
      else
        numstk(numnum) = v1
      endif
C ... Make sure that there is a '?' and pop it
      lerr = ((opnum == 0) .or. (opstk(opnum) /= 7))
      opnum = opnum-1

      end

      subroutine setans(resL,resI,resR,resD,cast,count,val)
C- cast val appropriately and copy to count'th element of res*
C -------------------------------------------------------------------
Ci Inputs
Ci   cast:   0=logical, 1=char, 2=int, 3=real, 4=double
Ci   count:  element to copy to
Ci   val:    value to copy in double form
Co Outputs
Co   resL:   logical array, count'th element is copied over
Co   resI:   integer array, count'th element is copied over
Co   resR:   real array, count'th element is copied over
Co   resD:   double array, count'th element is copied over
Cr Remarks
Cr   What's the purpose of casting to char? I don't know, so you
Cr   get an error.
C -------------------------------------------------------------------
      implicit none
      logical   resL(0:*)
      integer   resI(0:*),cast,count
      real      resR(0:*)
      double precision  resD(0:*),val

      goto (10,20,30,40,50),cast+1
   20 call rx('A2BIN: can''t cast result.')
   10 resL(count) = val /= 0D0
      return
   30 resI(count) = int(val)
      return
   40 resR(count) = val
      return
   50 resD(count) = val
      return
      end

      double precision function log2d(l)
C- Convert logical to double
C -------------------------------------------------------------------
Ci Inputs
Ci   l:   logical to be converted
Co Outputs
Co   Returns 1D0 for .true. and 0D0 for .false.
C -------------------------------------------------------------------
      implicit none
      logical   l
      log2d = 0D0
      if (l) log2d = 1D0
      return
      end
      logical function a2d(strn,lskp,j,jmax,term,res)
C- Convert ASCII string to number
C -------------------------------------------------------------------
Ci Inputs
Ci   strn:  string containing number
Ci   j:     location in string to start parsing number
Ci   jmax:  maximum value of j
Ci   lskp:  true, when a2d should skip blanks in number
Ci   term:  exit if a2d encounters character term
Co Outputs
Co   j:     location just after number in strn
Co   res:   double precision number
Co   a2d:   true if parse successful, false if not
Cr Remarks
Cr   a2d has a generic parser; however it is faster and more accurate
Cr   to use the Fortran internal read.  An unformatted read is
Cr   preferable, but does not conform to the ANSII 77 standard.  Do
Cr   Do ccomp -dUNFORMATTED_READ a2bin.f for this version
Cr   Some machines, such as the Cray, will correctly read an arbitrary
Cr   floating-point number using En.0 format; for this version do
Cr   comp -dFORMATTED_READ a2bin.f
Cr   Use the generic parser if your machine accepts neither of these.
Cr   ISFRAC and E1 are used to parse to the left of the decimal point.
Cr   OK is used to make sure that there is at least one decimal digit
Cr   in the mantissa and the exponent (if there is one).
C -------------------------------------------------------------------
      implicit none
      integer   i,j0,j,jmax,iprint
      logical   lskp
      character*1       t,strn(0:*),term
      double precision res,sgn
C#ifdef FORMATTED_READ
C#elseC
C      logical isfrac,neg,ok
C      integer e,e1
C#endif

C#ifdef FORMATTED_READ | UNFORMATTED_READ
      integer maxs,strs
      parameter (maxs=72)
      character*(maxs) strn2

C --- Find end of string; exit to 20 w/ strn(j-1) = last char of num ---
      a2d = .false.
      j0 = j
      if (lskp) call skipbl(strn,1000,j0)
      j = j0-1
   10 j = j+1
      if (jmax > 0 .and. j > jmax) goto 20
      if (j > j0+maxs) goto 99
      t = strn(j)
      if (t == term) goto 20
      if (t >= '0' .and. t <= '9') goto 10
      if (t == '.') goto 10
      if (t == 'd' .or. t == 'D' .or. t == 'e' .or. t == 'E') then
        if (lskp .and. strn(j) == ' ') call skipbl(strn,1000,j)
        if ((strn(j+1) == '+') .or. (strn(j+1) == '-')) j = j+1
        goto 10
      endif
      if (j == j0 .and. (t == '+' .or. t == '-')) goto 10
      if (lskp .and. (t == ' ' .or. t == achar(9))) goto 10
C --- Copy strn to strn2 for FORTRAN read ---
   20 continue
      if (lskp .and. j > j0 .and.
     .   (strn(j-1) == ' ' .or. strn(j-1) == achar(9))) then
        j = j-1
        goto 20
      endif
      strs = j-j0
      strn2 =  ' '
      call strcop(strn2(maxs+1-strs:maxs),strn(j0),strs,'z',i)
C#ifdef FORMATTED_READ
      read(strn2,'(g72.0)',err=99) res
C#elseC
C      read(strn2,*,err=99) res
C#endif
      a2d = .true.
      if (iprint() >= 130) print 333, strn2(maxs+1-strs:maxs), res
  333 format(' a2d: converted ', a,' to',g24.16)
      return
C --- Error handling ---
   99 continue
C#ifdefC APOLLO_BUG
C      if (iprint() >= 60) print 332, (strn(i), i=0, j-1)
C#else
      if (iprint() >= 60) print 332, (strn(i), i=j0, j-1)
C#endif
  332 format(' a2d: parse error for string ',72a1)
C#elseC  (not FORMATTED_READ | UNFORMATTED_READ)
CC --- Generic parser: initialization ---
C      sgn = 1
C      res = 0
C      e = 0
C      e1 = 0
C      isfrac = .false.
C      ok = .false.
C      j0 = j
C      if (lskp) call skipbl(strn,1000,j0)
C      j = j0
C      a2d = .false.
CC --- Parse mantissa ---
C   10 t = strn(j)
C      if (lskp .and. t == ' ') goto 20
C      if (t == '.') then
C        if (isfrac) return
C        isfrac = .true.
C        goto 20
C      elseif (j == j0 .and. t == '-') then
C        sgn = -1
C        goto 20
C      elseif (j == j0 .and. t == '+') then
C        goto 20
C      endif
C      if ((t < '0') .or. (t > '9')) goto 30
C      ok = .true.
C      res = res*10D0 + ichar(t)-ichar('0')
C      if (isfrac) e1 = e1-1
C   20 j = j+1
C      goto 10
CC --- Parse exponent ---
C   30 continue
C      if (.not. ok) return
C      ok = .false.
C      if ((t == 'd').or.(t == 'D').or.(t == 'e').or.(t == 'E')) then
C        j = j+1
C        if (lskp) call skipbl(strn,1000,j0)
C        if ((strn(j) == '+') .or. (strn(j) == '-')) then
C          neg = strn(j) == '-'
C          j = j+1
C          if (lskp) call skipbl(strn,1000,j0)
C        else
C          neg = .false.
C        endif
C   35   t = strn(j)
C        if ((t < '0') .or. (t > '9')) goto 38
C        ok = .true.
C        e = e*10 + ichar(t)-ichar('0')
C        j = j+1
C        if (lskp) call skipbl(strn,1000,j0)
C        goto 35
C   38   continue
C        if (.not. ok) return
C        if (neg) e= -e
C      endif
C      res = sgn*res*(10D0**(e+e1))
C      a2d = .true.
C
C      print *, 'a2d: converted ',(strn(i),i=j0,j-1), ' to ',res
C#endif FORMATTED_READ | UNFORMATTED_READ
      end

      logical function gettok(str,j,term,jmax)
C- Get next variable, number, or operator token
C -------------------------------------------------------------------
Ci Inputs
Ci   str:   string containing token
Ci   j,jmax:location in str to start parsing token, not to exceed jmax
Ci   term:  expression terminator character
Co Outputs
Co   j:     location just after token
Co   Returns whether a token was found
Co   toktyp: 1, a variable
Co           2, a number
Co           3, an operator or a function
Co           4, a vector
Co           5, a macro
Co           0, the 'assignment' operator (causes a2bin to return F)
Co   optok:  value  operator
Co           0+16   ~
Co           1      ^
Co           2      *
Co           2+16   /
Co           3      +
Co           3+16   -
Co           4      <
Co           4+16   >
Co           4+32   ==
Co           4+48   <=
Co           4+64   >=
Co           4+80   <>
Co           5      &
Co           6      |
Co           7      ?
Co           8      :
Co           9      (
Co           10     )
Co                  unop
Co           16*2   abs
Co           16*3   exp
Co           16*4   log
Co           16*5   sin
Co           16*6   asin
Co           16*7   sinh
Co           16*8   cos
Co           16*9   acos
Co           16*10  cosh
Co           16*11  tan
Co           16*12  atan
Co           16*13  tanh
Co           16*14  flor
Co           16*15  ceil
Co           16*16  erfc
Co           16*17  sqrt
Co           16*18  ran
Co           16*19  besj0
Co           16*25  vector element
Cr Remarks
Cr   If something that isn't some kind of token is encountered then
Cr   GETTOK returns .false.
C -------------------------------------------------------------------
      implicit none
      character*1 term,str(0:*)
      integer   j,jmax,namlen,maxsiz,parenc,nunop,k
      integer   opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok
      double precision numtok,numstk(0:32),machep,res(0:1)
      parameter (namlen=40,maxsiz=72,nunop=18)
      character*(namlen) vartok
      character*(maxsiz) strn
      logical   oprule(0:10,0:10)
      common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp,
     .                prevop,opstk,optok,parenc,oprule,vartok
      character*5 unop(0:nunop-1), unops(0:10), ctmp*1
      integer   i,s,t
      logical a2d,namchr
C ... Allowed characters in a name
      namchr(ctmp) =ctmp >= 'a' .and. ctmp <= 'z' .or. ctmp == '_'
     .         .or. ctmp >= 'A' .and. ctmp <= 'Z'
     .         .or. ctmp >= '0' .and. ctmp <= '9'
      data unops/'A~','^','*/%','+-','<>=','&','|','?',':','(',')'/
      data unop/'abs','exp','log','sin','asin','sinh','cos','acos',
     .          'cosh','tan','atan','tanh','flor','ceil','erfc','sqrt',
     .          'ran','besj0'/

C --- Check for terminator or end of expression
      gettok = .false.
      if (j > jmax) return
      if ((str(j) == term) .and. (parenc == 0)) then
        j = j+1
        return
      endif
C --- Save previous token type and operator code so that GETTOK can
C     distinguish unary negation from binary subtraction
      ptktyp = toktyp
      prevop = optok
      gettok = .true.
C --- Check for variable token first, variables are a letter followed
C     by a string of letters and digits
      if (((str(j) >= 'A') .and. (str(j) <= 'Z')) .or.
     .    ((str(j) >= 'a') .and. (str(j) <= 'z'))) then
        s = j
   10   j = j+1
        if (j > jmax) then
          if (s > jmax) j = j-1
          goto 11
        endif
        if (namchr(str(j))) goto 10
C        if (((str(j) >= 'A') .and. (str(j) <= 'Z')) .or.
C     .      ((str(j) >= 'a') .and. (str(j) <= 'z')) .or.
C     .      ((str(j) >= '0') .and. (str(j) <= '9'))) goto 10
   11   strn = ' '
        if (j-s > namlen) call rx('A2BIN: variable name too long.')
C --- Get variable (isn't FORTRAN string handling great?)
        call strcop(strn(1:j-s),str(s),j-s,' ',i)
        vartok = strn(1:j-s)
C   --- Check if variable is actually a macro
        if (str(j) == '(') then
          k = 0
          call macevl(strn,' ',k)
          if (k > 0) then
            toktyp = 5
            return
          endif
        endif
C --- Check if variable is actually a function or a vector
        call tokmat(vartok,unop,nunop,5,' ',i,s,.false.)
        if (i >= 0) then
          toktyp = 3
          optok = 16*(i+2)
        else
          call getsvv(vartok,i,0,1,1,res)
          if (i > 0) then
            toktyp = 4
            optok = 16*25
          else
            toktyp = 1
          endif
        endif
        return
C --- If it isn't a variable then check if it's a number
      else if ((str(j) == '.') .or.
     .         ((str(j) >= '0')  .and. (str(j) <= '9'))) then
        if (.not. a2d(str,.false.,j,jmax,term,numtok)) goto 99
        toktyp = 2
        return
C --- Otherwise it's an operator or a mistake
      else
        toktyp = 3
        do  30  i = 0, 10
          t = 0
C --- Check each operator class i for the operator
          call chrpos(unops(i),str(j),3,t)
          if (t /= 3) then
C --- Operator found, calculate the operator code
            optok = t*16 + i
            j = j+1
C           Valid special case unop is ')' and j=jmax+1
            if (j <= jmax+1 .and. i == 10) return
            if (j > jmax) goto 99
C --- See if it's unary negation
            if ((optok == 19) .and. (ptktyp == 3) .and. (prevop /= 10))
     .        optok = 0
C --- Distinguish '=' from '==': '=' only causes a2bin=F
            if (optok == 36) then
              if (str(j) /= '=') then
                toktyp = 0
                gettok = .true.
                return
              endif
              j = j+1
C --- Distinguish '<=' from '<'
            else if ((optok == 4) .and. (str(j) == '=')) then
              optok = 52
              j = j+1
              if (j > jmax) goto 99
C --- Distinguish '>=' from '>'
            else if ((optok == 20) .and. (str(j) == '=')) then
              optok = 68
              j = j+1
C --- Distinguish '<>' from '<'
            else if ((optok == 4) .and. (str(j) == '>')) then
              optok = 84
              j = j+1
            endif
            if (j > jmax) goto 99
            return
          endif
   30   continue
      endif
C --- If we fell through the op, check loop then an illegal character
C     is in the ASCII string.  Cause gettok to return .false.
   99 numnum = 0
      opnum = 0
      gettok = .false.
      end

      subroutine sstyle(style)
C- Sets style of expression parsing
C -------------------------------------------------------------------
Ci Inputs
Ci   style:   0=old style, 1=new style
Cr Remarks
Cr   Old style parsing means expressions are evaluated from left to
Cr   right with no operator precidence, but parenthesis work as
Cr   expected.  New style parsing means expressions are parsed with
Cr   algebraic precidence and associativity.  For old style parsing
Cr   the rule array OPRULE is changed to:
Cr
Cr        Current operator
Cr
Cr T    ~ ^ * + < & | ? : ( )
Cr o  ~ F T T T T T T T T F T
Cr p  ^ F T T T T T T T T F T
Cr    * F T T T T T T T T F T
Cr o  + F T T T T T T T T F T
Cr p  < F T T T T T T T T F T
Cr e  & F T T T T T T T T F T
Cr r  | F T T T T T T T T F T
Cr a  ? F F F F F F F F F F T
Cr t  : F T T T T T T T T F T
Cr o  ( F F F F F F F F F F F
Cr r  ) T T T T T T T T T T T
C -------------------------------------------------------------------
      implicit none
      integer style,i,j
      logical oprule(0:10,0:10),orule(0:10,0:10),nrule(0:10,0:10)
      integer opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok,namlen
      integer parenc
      double precision numtok,numstk(0:32),machep
      parameter (namlen=40)
      character*(namlen) vartok
      common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp,
     .                prevop,opstk,optok,parenc,oprule,vartok
      data nrule/10*.false.,2*.true.,9*.false.,4*.true.,7*.false.,
     .  5*.true.,6*.false.,6*.true.,5*.false.,7*.true.,
     .  4*.false.,8*.true.,3*.false.,8*.true.,.false.,.true.,
     .  .false.,8*.true.,.false.,.true.,.false.,.true.,
     .  10*.false.,10*.true.,.false.,.true./
      data orule/10*.false.,8*.true.,.false.,.true.,.false.,
     .  8*.true.,.false.,.true.,.false.,8*.true.,.false.,.true.,
     .  .false.,8*.true.,.false.,.true.,.false.,8*.true.,.false.,
     .  .true.,.false.,8*.true.,.false.,.true.,.false.,8*.true.,
     .  .false.,.true.,.false.,8*.true.,.false.,.true.,.false.,
     .  .true.,10*.false.,10*.true.,.false.,.true./

      if (style == 0) then
        do  10  i = 0, 10
          do  10  j = 0, 10
   10       oprule(i,j) = orule(i,j)
      else
        do  20  i = 0, 10
          do  20  j = 0, 10
   20       oprule(i,j) = nrule(i,j)
      endif
      return
      end

      block data da2bin
C- Block data for A2BIN
C -------------------------------------------------------------------
Cp Purpose
Cp   Initialize the array OPRULE once.
C -------------------------------------------------------------------
      integer opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok,namlen
      integer parenc
      double precision numtok,numstk(0:32),machep
      parameter (namlen=40)
      character*(namlen) vartok
      logical   oprule(0:10,0:10)
      common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp,
     .                prevop,opstk,optok,parenc,oprule,vartok
      data oprule/10*.false.,2*.true.,9*.false.,4*.true.,7*.false.,
     .            5*.true.,6*.false.,6*.true.,5*.false.,7*.true.,
     .            4*.false.,8*.true.,3*.false.,8*.true.,.false.,.true.,
     .            .false.,8*.true.,.false.,.true.,.false.,.true.,
     .            10*.false.,10*.true.,.false.,.true./
      end
