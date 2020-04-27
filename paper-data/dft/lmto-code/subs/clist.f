      subroutine clist(lstyle,slist,clabl,z,nclass,nlist,list)
C- Generates a list of classes from a string specification
C ----------------------------------------------------------------
Ci Inputs
Ci   slist:  string specifying list of classes
Ci   lstyle: style of slist specification; see Remarks
Ci   nclass  number of classes.
Co Outputs
Co   nlist,list a list of classes satisfying specifications
Cr Remarks
Cr *Syntax of slist: depends on one of three styles (lstyle)
Cr
Cr *lstyle=1 : a list of integers; see mkilst.f for complete syntax.
Cr             Example: '1,4:6,11' generates a list of five numbers,
Cr             1,4,5,6,11.
Cr
Cr *lstyle=2 : the list is specified according to an expression.
Cr             The expression can involve the class index ic and
Cr             atomic number z.  Any class satisfying expression is
Cr             included in the list.  Example:  'ic<6&z==14'
Cr
Cr *lstyle=3 : is specifically for unix systems.  slist is a filename
Cr             with the usual unix wildcards, eg a[1-6].  'clist'
Cr             makes a system 'ls' call for that string; any class
Cr             which 'ls' finds is included in the list.
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lstyle,nlist,nclass,list(*)
      character*(*) slist
      character*8 clabl(nclass)
      double precision z(nclass)
C ... Dynamically allocated arrays
      integer, allocatable :: llst(:)
C ... Local parameters
      integer iv0,ic,ival,is,i,j,ls,ls1
      logical a2bin,sw
      character strn*120,filnam*72,cnam*72

      ls = len(slist)
      nlist = 0
C -- lstyle=1 ---
      if (lstyle == 1) then

      call mkils0(slist,nlist,i)
C     call defi(ollst, nlist)
      allocate(llst(nlist))
      call mkilst(slist,nlist,llst)
      if (nlist == 0) return
      call ishell(nlist,llst)
      list(1) = llst(1)
      j = 1
      do  i = 2, nlist
        if (llst(i) > list(j) .and. llst(i) <= nclass) then
          j = j+1
          list(j) = llst(i)
        endif
      enddo
      nlist = j
      deallocate(llst)
      return

C --- lstyle=2 ---
      elseif (lstyle == 2) then

      call numsyv(iv0)
      nlist = 0
        do  ic = 1, nclass
        call lodsyv('ic',1,dble(ic),ival)
        call lodsyv('z',1,z(ic),ival)
        is = 0
          if ( .not. (a2bin(slist,sw,0,0,slist(ls:ls),is,ls))) then
            call rxs('clist: failed to parse',slist)
          elseif (sw) then
            nlist = nlist+1
            list(nlist) = ic
C   ... Abort if a2bin can't parse expression
        endif
        enddo
      call clrsyv(iv0)
      return

C --- lstyle=3 ---
      elseif (lstyle == 3) then
      nlist = 0
      call skpblb(slist,ls,ls1)
      call ffnam(slist(1:ls1+1),filnam)
        do  ic = 1, nclass
        call pshpr(1)
        call ffnam(clabl(ic),cnam)
        call poppr
        call awrit0('%xls ' // filnam //'%a|grep -s '
     .    // cnam // '%a>/dev/null',strn,len(strn),0)
        call locase(strn)
        call fsystm(strn,j)
        if (j == 0) then
          nlist = nlist+1
          list(nlist) = ic
        endif
        enddo
      else
        call rxi('clist: bad style',lstyle)
        return
      endif

      end
