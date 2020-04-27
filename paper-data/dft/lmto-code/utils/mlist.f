C#define MLIST
      integer function mkilsd(strn,nlist,list)
C- Resolve list (ascii string) into a vector of integers
C ----------------------------------------------------------------------
Ci Inputs
Ci   strn  :string holding list of integers
Ci   nlist :maximum number of integers in list
Ci         :nlist<0 => just return mkilsd = size of list;
Ci                     do not populate list
Co Outputs
Co   mkilsd:size of list
Co   list  :list of integers (not returned if input nlist<0)
Cl Local variables
Cl  lcount :0 return list and mkilsd = size of list
Cl         :1 return mkilsd only
Cl         :2 return -mkilsd only
Cr Remarks
Cr   This routine has the same syntax as mkilst (see mkilst below for syntax),
Cr   with extensions 'dup=#' and 'seq=#1,#2,...'.  Either may be used multiple times.
Cr
Cr   'dup' duplicates the entire list to this point, adding # to the duplicate elements.
Cr   The following example uses dup twice.  The first usage doubles the list from 2
Cr   to 4 elements, the second from 5 to 10 elements.
Cr      7:9:2,dup=2,0,dup=10  is expanded into  7 9 9 11 0 17 19 19 21 10
Cr
Cr   'seq' replicates the entire list to the current point, for each #1, #2, ...
Cr   #1-list(1) is added to the first replication, #2-list(1) to the second, etc,
Cr   so that replicated list n begins with n.  For example
Cr      7:9:2,seq=16,25   is expanded into  7 9 16 18 25 27
Cr   The list is thus expanded by factor n, n being the number of integers following seq
Cr
Cr   To continue the list after a seq tag, use a semicolon to delimit the
Cr   close of the tag from the remainder of the list.
Cr
Cr   See also documentation at www.questaal.org/docs/numerics/integerlists/
Cu Updates
Cu   10 Dec 17 seq tag can be terminated by ';', it no longer need be the last argument
Cu   11 Sep 17 Added 'seq' tag
Cu   30 Nov 16 Add 'dup' tag
Cu   02 Feb 01 strn is now a character string
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer list(*),nlist
      character*(*) strn
      integer,allocatable :: nseq(:)
C ... Local parameters
      integer ip,ip2,nlst1,nlst2,i,k,m,it(1),iv(1),a2vec,lcount,ms

      ip2 = 0
      lcount = 0
      if (nlist < 0) lcount = 1
      mkilsd = -1
      nlst1 = 0
      if (len(strn) <= 0) return
C      ip = 1
C      call skipbl(strn,len(strn),ip)
C      ip = ip+1
      k = index(strn,'dup=')
      m = index(strn,'seq=')-2
      if (m < 0) m = len(strn)
      if (m > 0) then
        ms = index(strn(m+2:),';')
        if (ms == 0 .and. k > m) return ! seq must come at the end
       endif
      ip = 1  ! Running index to current position in strn

C ... Handle lists with dup inside
C     k should point to dup, if it is present in the remaining list
C     m should point to last char in strn preceding ,seq or ;seq
   10 continue
      do  while (k >= ip)
        ip2 = k+3
        i = a2vec(strn,len(strn),ip2,2,', ',2,1,1,it,iv)
        if (i < 0) return

C       (ip:k) brackets new list; gather before duplicating
        k = k-1
        if (strn(ip:k) /= ' ') then
          call mkils0(strn(ip:k),nlst2,i)
          if (nlst2 == -1) return
          if (nlst1+nlst2 > nlist .and. lcount /= 1) lcount = 2
          if (lcount == 0) then
            call mkilst(strn(ip:k),nlst2,list(nlst1+1))
          endif
          nlst1 = nlst1+nlst2
        endif

C       Double list: second half = duplicate of original + shift
        if (nlst1*2 > nlist .and. lcount /= 1) lcount = 2
        if (lcount == 0) then
          do  k = 1, nlst1
            list(k+nlst1) = list(k) + iv(1)
          enddo
        endif
        nlst1 = 2*nlst1

        ip = ip2+1
        k = index(strn(ip:m),'dup=') + ip-1

      enddo

C ... Remaining list, possibly terminated before seq
      if (strn(ip:) /= ' ' .and. ip <= m) then
        if (strn(ip:min(ip+2,m)) == 'seq') then
          nlst2 = 0  ! no elements between ip and the next seq
        else
          call mkils0(strn(ip:m),nlst2,ip2)
          if (nlst1+nlst2 > nlist .and. lcount /= 1) lcount = 2
          if (lcount == 0) then
            call mkilst(strn(ip:m),nlst2,list(nlst1+1))
          endif
        endif
        nlst1 = nlst1+nlst2
      endif

      mkilsd = nlst1
      if (lcount == 2) mkilsd = -nlst1
C     m = index(strn,'seq=')

C ... Handle seq
      m = index(strn(ip:),'seq=') + ip-1
      if (mkilsd < 0 .or. m < ip) return
      ip = m+4
      ms = index(strn(m:),';')
      if (ms == 0) then
        ms = len_trim(strn)
      else
        ms = ms+m-2
      endif
      call mkils0(strn(ip:ms),nlst2,ip2) ! get nlst2 = # duplications
      m = (nlst2+1)*nlst1 ! number of new elements
      if (nlist > 0 .and. m > nlist) m = -m
      mkilsd = m
!     if (nlist < 0 .or. m < 0) return
      if (m < 0) return  ! Won't fit; return with error

C     Fill out list with duplicated elements
      if (nlist > 0) then
        allocate(nseq(nlst2))
        call mkilst(strn(ip:ms),nlst2,nseq)
        do  ip2 = 1, nlst2
          k = ip2*nlst1
          list(k+1:k+nlst1) = list(1:nlst1) + nseq(ip2) - list(1)
        enddo
        deallocate(nseq)
      endif

C ... Exit if string exhausted
      if (ms == len_trim(strn)) return

      ip = ms+2
      nlst1 = m ! The augmented list because the root list
      k = index(strn(ip:),'dup=') + ip-1
      m = index(strn(ip:),'seq=')-2 + ip-1
      if (m < ip) m = len(strn)
      goto 10

      end

      subroutine mkilst(strn,nlist,list)
C- Resolve list (ascii string) into a vector of integers
C ----------------------------------------------------------------------
Ci Inputs
Ci   strn  :string holding list of integers
Co Outputs
Co   nlist :number of integers in list
Co         :nlist<0 => mkilst failed to parse list
Co   list  :list of integers
Cr Remarks
Cr   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
Cr   low:high:step
Cr   low, high, and step are integer expressions specifying the sequence
Cr     low, low+step, low+2*step, ... high.
Cr   If :step is missing, the step size defaults to 1.  If also :high
Cr   is missing,  the sequence reduces to a single integer. Thus,
Cr     '5+1'       becomes a single number, 6.
Cr     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
Cr     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
Cr   Sequences may be strung together separated by commas, eg
Cr     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
Cu Updates
Cu   02 Feb 01 strn is now a character string
C ----------------------------------------------------------------------
      implicit none
      integer list(*),nlist
      character*(*) strn
      integer it(512),iv(512),a2vec,ip,i,j,k

      ip = 0
      nlist = -1
      call skipbl(strn,len(strn),ip)
      k = a2vec(strn,len(strn),ip,2,',: ',3,3,100,it,iv)
      if (k < 1) return
      if (k >= 99) call rx('mkilst: increase size of iv')
      it(k+1) = 0
      iv(k+1) = iv(k)
C ... loop over all iv
      nlist = 0
      i = 0
   14 i = i+1
C ... Case iv => a single number
      if (it(i) /= 2) then
        nlist = nlist+1
        list(nlist) = iv(i)
C ... Case iv => n1:n2:n3
      elseif (it(i+1) == 2) then
        do  12  j = iv(i), iv(i+1), iv(i+2)
        nlist = nlist+1
   12   list(nlist) = j
        i = i+2
C ... Case iv => n1:n2
      else
        do  17  j = iv(i), iv(i+1)
        nlist = nlist+1
   17   list(nlist) = j
        i = i+1
      endif
      if (i < k) goto 14
      end

      subroutine mkilss(iopt,slst,nlist,list)
C- Resolve a list into a sorted vector of integers, paring duplicates
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt  :1s digit
Ci         : 0, do not sort list or pare duplicates
Ci         : 1, sort list, paring duplicates
Ci         :10s digit
Ci         : 0 return nlist<0 if fail to parse list
Ci         : 1 abort with error message if fail to parse list
Ci   slst  :string holding list of integers
Co Outputs
Co   nlist :number of integers in list
Co         :nlist<0 => mkilss failed to parse list
Co   list  :list of integers
Cr Remarks
Cr   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
Cr   low:high:step
Cr   low, high, and step are integer expressions specifying the sequence
Cr     low, low+step, low+2*step, ... high.
Cr   If :step is missing, the step size defaults to 1.  If also :high
Cr   is missing,  the sequence reduces to a single integer. Thus,
Cr     '5+1'       becomes a single number, 6.
Cr     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
Cr     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
Cr   Sequences may be strung together separated by commas, eg
Cr     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
Cu Updates
Cu   30 Nov 16 calls mkilsd, to allow ~dup and ~seq tags
Cu   20 Oct 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iopt,list(*),nlist
      character*(*) slst
C ... Local parameters
      integer i,j
      procedure(integer) :: mkilsd

C     call mkilst(slst,nlist,list)
      nlist = mkilsd(slst,-1,[j])
      if (nlist < 0 .and. (iopt/10 > 0))
     .  call rxs('mkilss: failed to parse integer list ',slst)
      if (nlist < 0) return
      nlist = mkilsd(slst,nlist,list)
      if (mod(iopt,10) == 0) return

C     Sort list and pare duplicates
      call ishell(nlist,list)
      j = 1
      do  16  i = 2, nlist
        if (list(i) > list(j)) then
          list(j+1) = list(i)
          j = j+1
        endif
   16 continue
      nlist = j
      end

      subroutine mkilssr(iopt,slst,nlist,list,irnge)
C- Resolve a list into a sorted vector of integers, paring duplicates, checking bounds
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt  :1s digit
Ci         : 0, do not sort list or pare duplicates
Ci         : 1, sort list, paring duplicates
Ci         :10s digit
Ci         : 0 return nlist<0 if fail to parse list
Ci         : 1 abort with error message if fail to parse list
Ci   slst  :string holding list of integers
Co Outputs
Co   nlist :number of integers in list
Co         :nlist<0 => mkilss failed to parse list
Co   list  :list of integers
Cr Remarks
Cr   Identical to mkilss except bounds are checked.
Cr   mkilssr aborts if any element is not within [irnge(1),irnge(2)]
Cu Updates
Cu   30 Nov 16 calls mkilsd, to allow ~dup tag
Cu   20 Oct 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iopt,list(*),nlist,irnge(2)
      character*(*) slst
C ... Local parameters
      integer i

      call mkilss(iopt,slst,nlist,list)
      do  i = 1, nlist
        if (list(i) < irnge(1) .or. list(i) > irnge(2))
     .    call rx2('mkilssr: element %i outside range [%s,%2i]',list(i),irnge)
      enddo

      end
      subroutine mkils0(strn,nlist,ip)
C- Like mkilst, but returns size of list and last char parsed in strn
      implicit none
      integer nlist
      character*(*) strn
      integer it(100),iv(100),a2vec,ip,i,j,k

      ip = 0
      nlist = -1
      call skipbl(strn,len(strn),ip)
      if (ip >= len(strn)) return
      k = a2vec(strn,len(strn),ip,2,',: ',3,3,100,it,iv)
      if (k < 1) return
      if (k >= 99) call rx('mkilst: increase size of iv')
      it(k+1) = 0
      iv(k+1) = iv(k)
C ... loop over all iv
      nlist = 0
      i = 0
   14 i = i+1
C ... Case iv => a single number
      if (it(i) /= 2) then
        nlist = nlist+1
C ... Case iv => n1:n2:n3
      elseif (it(i+1) == 2) then
        do  12  j = iv(i), iv(i+1), iv(i+2)
   12   nlist = nlist+1
        i = i+2
C ... Case iv => n1:n2
      else
        do  17  j = iv(i), iv(i+1)
   17   nlist = nlist+1
        i = i+1
      endif
      if (i < k) goto 14

      end
      subroutine ilst2a(list,nlist,strn)
C- Make an ascii represention of a list from a vector of integers
C ----------------------------------------------------------------------
Ci Inputs
Co   list  :list of integers
Co   nlist :number of integers in list
Co Outputs
Ci   strn  :string holding list of integers
Cr Remarks
Cr   ilst2a performs the inverse operation to mkilst,
Cr   making a compact ascii representation of a list
Cu Updates
Cu   13 Sep 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlist,list(nlist)
      character*(*) strn
C ... Local variables
      integer i,i1,ls,awrite,ip,j
C     integer iilst(9)
C     data iilst /1,3,5,2,4,6,8,10,99/

      strn = ' '
      ip = 1
      if (nlist == 0) return
      ls = len(strn)
      i1 = 1
      i = 0
   10 continue
        i = i+1
        if (i > nlist) then
          return
        endif
        if (ip > 1) then
          strn(ip:ip) = ','
          ip = ip+1
        endif

C   ... Near the end of list --- no sequences to check
        if (i > nlist-2) then
          ip = ip+awrite('%i',strn(ip:),ls,0,list(i),0,0,0,0,0,0,0)
          goto 10
        endif

C   ... Case linear sequence of at least three
        if (list(i+1)-list(i) == list(i+2)-list(i+1)) then
          i1 = i
          j = list(i+1)-list(i)
   20     continue
          if (i <= nlist-2) then
            i = i+1
            if (j == list(i+2)-list(i+1)) goto 20
          endif
          if (j == 1) then
            ip = ip+awrite('%i:%i',strn(ip:),ls,0,
     .        list(i1),list(i+1),j,0,0,0,0,0)
          else
            ip = ip+awrite('%i:%i:%i',strn(ip:),ls,0,
     .        list(i1),list(i+1),j,0,0,0,0,0)
          endif
          i = i+1
          goto 10
        endif

C   ... Case no linear sequence ... just append single element
        ip = ip+awrite('%i',strn(ip:),ls,0,list(i),0,0,0,0,0,0,0)
        goto 10
      end

C#ifdefC TEST
C      subroutine fmain
C      implicit none
C      character*60 strn
C      integer nlist,list(99),mkilsd
C
C      strn = '                 2,1'
C      call mkilst(strn,nlist,list)
C      print *, 'mkilst',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C
C      strn = '             22:33:3,5+1:8+2:2'
C      call mkilst(strn,nlist,list)
C      print *, 'mkilst',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C
C      strn = '             22:33:3,dup=100'
C      nlist = mkilsd(strn,99,list)
C      print *, 'mkilsd',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C
C      strn = '             22:33:3,dup=100,5+1:8+2:2'
C      nlist = mkilsd(strn,99,list)
C      print *, 'mkilsd',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C
C      strn = '      22:33:8,dup=100,5+1:8+2:2,seq=2,12'
C      nlist = mkilsd(strn,99,list)
C      print *, 'mkilsd',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C
C      strn = '      22:33:8,dup=100,5+1:8+2:2,seq=2,12;1000,dup=10000'
C      nlist = mkilsd(strn,99,list)
C      print *, 'mkilsd',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C
C      strn = '      9,25,30,seq=39:159:30;seq=9+1000'
C      nlist = mkilsd(strn,99,list)
C      print *, 'mkilsd',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C
C      strn = '             22:33:3,dup=10,30:40:2'
C      call mkilssr(11,strn,nlist,list,[20,50])
C      print *, 'mkilssr',trim(strn)
C      call info2(0,0,0,' %i  %-1j%n:1i',nlist,list)
C      print *, 'should match:'
C      print *, '12   22 25 28 30 31 32 34 35 36 38 40 41'
C
C      call mkilssr(11,strn,nlist,list,[22,39])
C      print *, 'should not be here'
C
C      end
C#endif

C#ifdef MLIST
C This is the main entry point for the mlist utility.
C Resolve list represented as ascii string into a vector of doubles,
C to be printed to standard out.
C Alternatively you can specify an integer list.  This calls mkilss
C and follows the standard Questaal syntax for integer lists (see mkilst.f)
C The floating point case is a restricted form, not enabling the seq and dup tags.
C
C Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
C   low:high:step
C Example: mlist 1.1:4:.2 generates the following:
C    1.1 1.3 1.5 1.7 1.9 2.1 2.3 2.5 2.7 2.9 3.1 3.3
C low, high, and step are integer expressions specifying the sequence
C   low, low+step, low+2*step, ... high.
C If :step is missing, the step size defaults to 1.  If also :high
C is missing,  the sequence reduces to a single integer. Thus,
C   '5+1'       becomes a single number, 6.
C   '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
C   '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
C Sequences may be strung together separated by commas, eg
C   '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
C Example:
C    mlist -sort -ndec=1 -fd3:0 2.6:3.8:.2,3.058
C returns
C    2.6 2.8 3.0 3.058 3.2 3.4 3.6 3.8
C
C Integer lists call mkilss; see (slatsm/mkilst.f or
C added flexibility, enabling the ~seq and ~dup tags.
C Examples:
C   mlist -i 1:11:2,seq=2
C returns
C   1 3 5 7 9 11 2 4 6 8 10 12
C whereas
C   mlist -i -sort 1:11:2,seq=2
C returns
C   1 2 3 4 5 6 7 8 9 10 11 12
C whereas
C   mlist -i -sort -sep:, 1:11:2,seq=2
C returns
C   1,2,3,4,5,6,7,8,9,10,11,12
C
C Bugs: should pass max len of strn, have max size of list
C Updates
C   05 Aug 18 Added ability to call mkilss
C             Added -n : returns number of elements in list
C ----------------------------------------------------------------------
      subroutine fmain
      implicit none
      double precision d1mach,xxv(256),cval(500),tmp
      character outstr*2000,first*200,fmt*40
      character*1 f2(20),sep
      logical lcount,lsort,ilst,ipack
      integer iarg,n,i,ip,it(1024),ncval,k,nblk,ndec,stdo
      equivalence (f2,first)
      procedure(logical) :: cmdstr,lsequ,lgunit
      procedure(integer) :: mkilsd,a2vec

C     fmt = '(3f15.10)'
      ndec = 0
      fmt = 'g'
      sep = ' '
      iarg = 1
      stdo = lgunit(1)
      lcount = .false.; lsort = .false.; ilst = .false.; ipack = .false.
   15 continue
      if (.not. cmdstr(iarg,first)) goto 20
      if (first(1:1) .ne. '-') goto 30

      if (first(1:2) .eq. '-f') then
        fmt = first(3:)
C        fmt = '('
C        call strcat(fmt,1,' ',f2(3),99,' ',n)
C        call strcat(fmt,99,' ',')',1,' ',n)
      elseif (first(1:5) .eq. '-sep:') then
        sep = first(6:6)
      elseif (first(1:6) .eq. '-ndec=') then
        ip = 6
        k = a2vec(first,len(first),ip,2,' ',1,1,1,it,ndec)
        call sanrg(.true.,ndec,0,10,'mlist','ndec')
      elseif (first(1:5) .eq. '-n') then
        lcount = .true.
      elseif (first(1:5) .eq. '-sort') then
        lsort = .true.
      elseif (first(1:4) .eq. '-il') then
        ilst = .true.
        ipack = .true.
      elseif (first(1:4) .eq. '-i') then
        ilst = .true.
      elseif (lsequ(first,'- ',2,' ',n)) then
        iarg = iarg+1
        if (.not. cmdstr(iarg,first)) goto 20
        goto 30
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15

C --- Resolve syntax for range, generate array of numbers ---
   30 continue

      if (ilst) then
        k = 10; if (lsort) k = k+1
        call mkilss(k,trim(first),ncval,it)
        if (lcount) then
          call awrit1('%i',' ',80,stdo,ncval)
          call fexit(0,0,'',0)
        endif

        if (ipack) then
          call ilst2a(it,ncval,outstr)
          call awrit0(trim(outstr),' ',len(outstr),stdo)
        else
          outstr = '%n:1i'
          if (sep /= ' ') outstr = '%s'//sep//'%ni'
          call awrit2(trim(outstr),' ',512,stdo,ncval,it)
        endif
        call fexit(0,0,'',0)
      endif

      ip = 0
      k = a2vec(first,len(first),ip,4,',: ',3,3,256,it,xxv)
      if (k .lt. 1) goto 20
C ... loop over all xxv
      ncval = 0
      i = 0
   14 i = i+1
C ... Case xxv => a single number
      if (it(i) .ne. 2) then
        ncval = ncval+1
        cval(ncval) = xxv(i)
C ... Case xxv => n1:n2:n3
      elseif (it(i+1) .eq. 2) then
        do  12  tmp = xxv(i), xxv(i+1)*(1+2*d1mach(3)), xxv(i+2)
        ncval = ncval+1
   12   cval(ncval) = tmp
        i = i+2
C ... Case xxv => n1:n2
      else
        do  17  tmp = xxv(i), xxv(i+1)*(1+2*d1mach(3))
        ncval = ncval+1
   17   cval(ncval) = tmp
        i = i+1
      endif
      if (i .lt. k) goto 14

      if (lsort) then
        call dvheap(1,ncval,cval,it,0d0,0)
      endif

      outstr = ' '
      ip = 0
      nblk = 1
      if (sep .ne. ' ') nblk = 0

      if (lcount) then
        call awrit1('%x%i',' ',80,stdo,ncval)
        return
      endif

      call bin2av(fmt,0,nblk,ndec,cval,4,0,ncval-1,sep,len(outstr),.false.,outstr,ip)
      call cwrite(outstr,0,ip,1)
      call fexit(0,0,'',0)

   20 print *, 'usage: mlist [-n -i -sort -ndec=# -ffmt -sep:c -] list'
      print 111, '       -n     return size of list, not the list itself'
      print 111, '       -i     list is an integer list; use mkilsd to enable dup and seq'
      print 111, '       -sort  sort the list (and duplicates pared with -i)'
      print 111, '       -ndec  specify minimum number of decimals to print out. Not used with -i.'
      print 111, '       -ffmt  formatting (default=''g''); see bin2a conventions.  Not used with -i.'
      print 111, '       -sep   specify character separating elements to be printed'
  111 format(1x,a)
      call cexit(1,1)

      end
C#endif
