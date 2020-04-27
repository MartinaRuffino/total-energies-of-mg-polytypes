      subroutine query(instr,cast,var)
C- interactive flow control
C ----------------------------------------------------------------
Ci Inputs
Ci   strng: prompt string
Ci   cast:  <0, if nothing to change, otherwise
Ci          cast is  0=,logical, 2=int, 3=real 4=double
Co Outputs
Co   var:   query will change if requested
Cr Remarks
Cr   At the prompt, user enters either nothing, or one of
Cr     'S nnn', where nnn is number (or T or F for logical variable);
Cr     'V nnn', where nnn is the new verbosity;
Cr     'W' to toggle printing of work array;
Cr     'I' to turn off interactive prompt;
Cr     'A' to abort program execution
Cr     'T' to toggle
Cr   Documented here:
Cr   questaal.org/docs/input/inputfile/#interactive-mode
Cr   Example using file iact.ext:
Cr     s beta .3
Cu Updates
Cu   07 May 17 Updated noninteractive query from file read
Cu   11 Dec 14 T can be used as to 'T=#' to turn on timing of nested blocks
Cu   07 Jul 04 suppress interactive mode in MPI
Cu   19 May 04 Added getqu: returns value of interative mode
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) instr
      integer cast,var
      logical lset
C ... Local parameters
      integer lprmpt
      parameter (lprmpt=120)
      character*(lprmpt) rdstr, outstr
      logical lquery,lsequ,a2bin,lact,ltmp
      integer i,ivbs,iprint,j,i1mach,l,fxst,fopna,ifi,awrite
      double precision xx
      integer mpipid
      save lquery
      character*(1) rdstr2(2)
      equivalence (rdstr,rdstr2)
      external lsequ,a2bin,iprint

C     do nothing if MPI call
      if (mpipid(0) > 1) return

      ifi = 29
      lact = fxst('iact') /= 0
      goto 12
C --- Start of interactive loop (read from iact first pass only) ---
   10 continue
      if (lact) call dfclos(ifi)
      lact = .false.
   12 continue
      if (lact) then
        ifi = 29
        ifi = fopna('iact',ifi,1)
        rewind ifi
      endif
      if (.not. lquery .and. .not. lact) return

      outstr = ' QUERY:'
      outstr(9:lprmpt) = instr
      call skpblb(outstr,lprmpt,l)
      if (cast >= 0) then
        l = awrite('%a (def=',outstr,lprmpt,0,0,0,0,0,0,0,0,0)
        call bin2a(' ',0,0,var,cast,0,lprmpt,outstr,l)
        call awrit0('%a)',outstr,lprmpt,0)
      endif
      l = awrite('%a?',outstr,lprmpt,0,0,0,0,0,0,0,0,0)
      call cwrite(outstr,0,l,0)
  334 format(a60)
C --- Read string, exit if empty ---
      rdstr = ' '
      if (lact) then
        read(ifi,334,end=10,err=10) rdstr
        print *, trim(rdstr)
      else
        read(*,334) rdstr
      endif
      j = 0
      call skipbl(rdstr,lprmpt,j)
      if (j >= lprmpt) return
      j = j+1
      i = 0
C --- Handle query ---
   14 continue
      if (rdstr2(j) == '?') then
        outstr = ' (A)bort  (I)active  (V)erb  (C)pu  (T)iming  (W)ork'
        if (cast >= 0) call awrit0('%a  (S)et value',outstr,lprmpt,0)
        print '(a)',  outstr
        goto 10
      elseif (rdstr2(j) == 'S' .or. rdstr2(j) == 's') then
        if (cast < 0) then
          call info0(1,0,0,'no variable to set--continuing ...')
          goto 10
        endif
        if (lact) then
          i = index(rdstr(j:),instr)
          if (i == 0) then
            lact = .false.
            call fclose(ifi)
            goto 10
          endif
          j = j + i + len(instr) - 2
        endif
        if (.not. a2bin(rdstr,var,cast,0,' ',j,-1)) print *, 'conversion error'
        goto 10
      elseif (rdstr2(j) == 'V' .or. rdstr2(j) == 'v') then
        if (a2bin(rdstr2(j+1),ivbs,2,0,' ',i,-1))  then
          call setpr(ivbs)
        else
          print *, 'conversion error'
        endif
        goto 10
      elseif (rdstr2(j) == 'W' .or. rdstr2(j) == 'w') then
C        call wkprnt(2)
C        call wkchk('Called from query')
C        call wkinfo
        goto 10
      elseif (rdstr2(j) == 'I' .or. rdstr2(j) == 'i') then
        lquery = .not. lquery
      elseif (rdstr2(j) == 'A' .or. rdstr2(j) == 'a') then
        call rx0(outstr)
      elseif (rdstr2(j) == 'Q' .or. rdstr2(j) == 'q') then
        if (lact) call dfclos(ifi)
        call rx0(outstr)
      elseif (rdstr2(j) == 'T' .or. rdstr2(j) == 't') then
        call tc('tog')
        call tcget(ltmp)
        if (.not. ltmp) then
          print *, 'turn off trace'
          goto 10
        endif
        if (rdstr2(j+1) /= ' ') then
          if (a2bin(rdstr,i,2,0,' ',j,-1))  then
            print "(' turn on timing, nest=',i3)" , i
            call tcinit(i,i)
          else
            print *, 'could not parse integer following "t"'
          endif
        else
          print *, 'turn on trace'
        endif
        goto 10
      elseif (rdstr2(j) == 'C' .or. rdstr2(j) == 'c') then
        call cpudel(i1mach(2),'called from query:   ',xx)
        goto 10
C ... string not empty, but no keyword recognized.
      elseif (.not. lact) then
        if (cast < 0) return
        j = 0
        outstr = ' '
        if (.not. a2bin(rdstr,var,cast,0,' ',j,-1)) then
          call awrit0(' ? input "'//rdstr//'%a" not recognized',
     .      outstr,lprmpt,i1mach(2))
          rdstr2(1) = '?'
          j = 1
          goto 14
        endif
        goto 10
      endif

      if (lact) goto 10
      return

      entry initqu(lset)
      lquery = lset
      return

      entry getqu(lset)
      lset = lquery

      end

      subroutine querym(instr,cast,var)
C- Interactive flow control, MPI version (all processors must call it)
C ----------------------------------------------------------------
Ci Inputs
Ci   strng: prompt string
Ci   cast:  <0, if nothing to change, otherwise
Ci          cast is  0=,logical, 2=int, 3=real 4=double
Co Outputs
Co   var:   query will change if requested
Cr Remarks
Cr   At the prompt, user enters either nothing, or one of
Cr     'Snnn', where nnn is number (or T or F for logical variable);
Cr     'Vnnn', where nnn is the new verbosity;
Cr     'I' to turn off interactive prompt;
Cr     'A' to abort program execution
Cu Updates
Cu   07 Jul 04 suppress interactive mode in MPI
Cu   19 May 04 Added getqu: returns value of interative mode
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*(*) instr
      integer cast,var
      logical lset
C Local parameters
      integer lprmpt
      parameter (lprmpt=120)
      character*(lprmpt) rdstr, outstr
      logical lquery,lsequ,a2bin,lact
      integer i,ivbs,iprint,j,i1mach,l,fxst,fopna,ifi,awrite
      double precision xx
C     for MPI
      integer procid, master, mpipid
      parameter (master = 0)
      character*(1) rdstr2(10)
      equivalence (rdstr,rdstr2)
      external lsequ,a2bin,iprint

      if (mpipid(0) == 1) then
        call query(instr,cast,var)
        return
      endif

      ifi = 29
      procid = mpipid(1)
      if (procid == master) then
        lact = fxst('iact') /= 0
      endif
      call mpibc1(lact,1,1,.false.,' ',' ')
      goto 12

C --- Start of interactive loop (read from iact first pass only) ---
   10 continue
      if (lact .and. procid == master) call dfclos(ifi)
      lact = .false.
   12 continue
      if (lact .and. procid == master) then
        ifi = 29
        ifi = fopna('iact',ifi,1)
        rewind ifi
      endif
      lquery = .false.
      if (.not. lquery .and. .not. lact) return

      if (procid == master) then
      outstr = ' QUERY:'
      outstr(9:lprmpt) = instr
      call skpblb(outstr,lprmpt,l)
      if (cast >= 0) then
        l = awrite('%a (def=',outstr,lprmpt,0,0,0,0,0,0,0,0,0)
        call bin2a(' ',0,0,var,cast,0,lprmpt,outstr,l)
        call awrit0('%a)',outstr,lprmpt,0)
      endif
      l = awrite('%a?',outstr,lprmpt,0,0,0,0,0,0,0,0,0)
      call cwrite(outstr,0,l,0)
  334 format(a60)
      endif

C --- Read string, exit if empty ---
      rdstr = ' '
      if (lact) then
        if (procid == master) then
C         print *, ' '
          read(ifi,334,end=10,err=10) rdstr
          print *, trim(rdstr)
        endif
        call mpibcc(rdstr,len(rdstr),.false.,' ',' ')
      else
        read(*,334) rdstr
      endif
      j = 0
      call skipbl(rdstr,lprmpt,j)
      if (j >= lprmpt) return
      j = j+1
      i = 0

C --- Handle query ---
   14 continue
      if (rdstr2(j) == '?') then
        outstr = ' (A)bort  (I)active  (V)erb  (C)pu  (T)iming  (W)ork'
        if (cast >= 0) call awrit0('%a  (S)et value',outstr,lprmpt,0)
        print '(a)',  outstr
        goto 10
      elseif (rdstr2(j) == 'S' .or. rdstr2(j) == 's') then
        if (cast < 0) then
          call info0(1,0,0,'no variable to set--continuing ...')
          goto 10
        endif
        if (lact) then
          i = index(rdstr(j:),instr)
          if (i == 0) then
            lact = .false.
            if (procid == master) call fclose(ifi)
            goto 10
          endif
          j = j + i + len(instr) - 2
        endif
        if (.not. a2bin(rdstr,var,cast,0,' ',j,-1))
     .    print *, 'conversion error'
        goto 10
      elseif (rdstr2(j) == 'V' .or. rdstr2(j) == 'v') then
        if (a2bin(rdstr2(j+1),ivbs,2,0,' ',i,-1))  then
          call setpr(ivbs)
        else
          print *, 'conversion error'
        endif
        goto 10
      elseif (rdstr2(j) == 'W' .or. rdstr2(j) == 'w') then
C        call wkprnt(2)
C        call wkchk('Called from query')
C        call wkinfo
        goto 10
      elseif (rdstr2(j) == 'I' .or. rdstr2(j) == 'i') then
        lquery = .not. lquery
      elseif (rdstr2(j) == 'A' .or. rdstr2(j) == 'a') then
        call rx0(outstr)
      elseif (rdstr2(j) == 'Q' .or. rdstr2(j) == 'q') then
        if (lact .and. procid == master) call dfclos(ifi)
        call rx0(outstr)
      elseif (rdstr2(j) == 'T' .or. rdstr2(j) == 't') then
        call tc('tog')
        goto 10
      elseif (rdstr2(j) == 'C' .or. rdstr2(j) == 'c') then
        call cpudel(i1mach(2),'called from query:   ',xx)
        goto 10
C ... string not empty, but no keyword recognized.
      elseif (.not. lact) then
        if (cast < 0) return
        j = 0
        outstr = ' '
        if (.not. a2bin(rdstr,var,cast,0,' ',j,-1)) then
          call awrit0(' ? input "'//rdstr//'%a" not recognized',
     .      outstr,lprmpt,i1mach(2))
          rdstr2(1) = '?'
          j = 1
          goto 14
        endif
        goto 10
      endif

      if (lact) goto 10
      return

      end
C      subroutine fmain
C      double precision x
C      logical l
C      integer procid,mpipid
C
C      procid = mpipid(1)
C      call initqu(.true.)
CC      call initqu(.false.)
CC     call query('n?',-1,l)
C      l = .true.
C      call query('l',0,l)
C      print *, 'procid=',procid, ' l=',l
C      call query('i=?',2,i)
C      print *, 'procid=',procid, ' i=',i
C      x = -2
C      call query('d=?',4,x)
C      print *, 'procid=',procid, ' d=',x
C      end
