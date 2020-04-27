C#define AWRITE
      integer function idalloc(nam,add,n1,n2)
C- Maintains a log of dynamical arrays allocated, and their sizes
C ----------------------------------------------------------------------
Ci Inputs
Ci   nam   :If nonblank, name of array to be managed
Ci         :If blank, special printout mode
Ci   add   :controls what routine does
Ci         :... Case nam is nonblank:
Ci         :0 return size of current allocation (MB)
Ci         :1 This array is to be allocated with dimension n1
Ci         :2 Same as 1, but stop with error if array already allocated
Ci         :3 This array is to be deallocated
Ci         :4 Same as 3, but stop with error if array is not already allocated
Ci         :Add a multiple of 10 to print actions to stdout
Ci         :... Case nam is blank
Ci         :0 return total size of all allocations (MB)
Ci         :1 return maximum allocation (MB)
Ci         :Add a multiple of 10 to print actions to stdout
Ci   n1    :(add = 1,2 only) new array to be allocated as arr(n1,n2)
Ci         :Size is assumed to be double precision
Ci         :For add=3,4 use same n1 or n1=0 to avoid check on size
Ci   n2    :(add = 1,2 only) new array to be allocated as arr(n1,n2)
Ci         :Size is assumed to be double precision
Ci         :For add=3,4 use same n1 or n1=0 to avoid check on size
Co Outputs
Co idalloc :array size
Cs Command-line switches
Cl Local variables
Cl  minsiz :(call idasiz to modify this parameter)
Cl         :The table is not modified for arrays with size smaller than minsiz.
Cl         :Note: minsize is MB.
Cr Remarks
Cr *A table of entries for dynamic allocations is kept.
Cr  Each entry is associated with 'nam'.
Cr  If idalloc is called with a name not in the list, it is added to the list.
Cr  Call idalloc with blank name to interrogate the current status.
Cr  A log of the largest cumulative allocation is also kept.
Cr *If the given array size is smaller than minsiz the routine exits before
Cr  the table is modified (minsiz is in MB)
Cr  Call idasiz to alter minsiz's value .
Cr  Each entry has a size:
Cr     size>0 => array has been allocated
Cr     size=0 => array has never been allocated
Cr     size<0 => array was allocated and deallocated
Cu Updates
Cu   13 May 13  First created
C ----------------------------------------------------------------------
      implicit none
      character nam*(*)
      integer add,n1,n2
C Local
      integer mxnam,strnsiz,idasiz
      parameter (mxnam=30,strnsiz=8)
      integer,save :: nnam=0,minsiz=10
      real(8),save :: asize(mxnam)
      real(8),save :: xx,mxalloc
      real(8),parameter :: MB=8d0/1d6
      character*(strnsiz*mxnam),save :: namlst
      character*(strnsiz) :: anam(mxnam),naml
      equivalence (anam,namlst)
      integer k,kk,m

C ... One-time initialization
      if (nnam == 0) then
        namlst = ' '
        mxalloc = 0
      endif

C ... Special case blank name: printout mode
      if (nam == ' ') then

        kk = 0; xx = 0 ! Loop counts number and cumulative allocation
        do  k = 1, nnam
          if (asize(k) > 0) then
            kk = kk+1
            xx = xx + asize(k)
          endif
        enddo
        idalloc = xx/1d6 * 8

        if (mod(add,10) == 1) then
          idalloc = mxalloc/1d6 * 8
        endif

        if (add >= 10) then

C#ifdef AWRITE
          call info5(10,1,0,' idalloc:%,3i arrays allocated for'//
     .      '%;9,1D MB.  Max alloc:  %i.',kk,xx*MB,nint(mxalloc*MB),0,0)
C#elseC
C          write(*,101) kk,xx*MB,mxalloc*MB
C#endif
  101     format(/' idalloc:',i3,' arrays allocated for',f9.1,' MB.',
     .      '  Max alloc:',f9.0)
          do  k = 1, nnam

C#ifdef AWRITE
            call info2(10,0,0,'   '//anam(k)//
     .        '%;9,1D%?;n<0;* when last allocated;;',asize(k)*MB,
     .        int(sign(1d0,asize(k))))
C#elseC
C            if (asize(k) > 0) then
C              write(*,102) anam(k),asize(k)*MB
C  102         format(3x,a,f9.1:a,' when last allocated')
C            elseif (asize(k) < 0) then
C              write(*,102) anam(k),asize(k)*MB,'*'
C            else
C              write(*,"(3x,a,7x,'--')") anam(k)
C            endif
C#endif
          enddo
        endif

        return
      endif

C ... Do not look at table unless dimension exceeds minimum
C     or if array has null dimension
      if (n1*n2 /= 0 .and. (n1*MB)*n2 < minsiz) then
        return
      endif

C ... Fast search for nam already in list
C     Exits loop with k matching element; k=0 if no match found
C     Jumps to label 10 if match found
      kk = 0
      do
        k = index(namlst(1+kk*strnsiz:),nam)
C       No match
        if (k == 0) exit
        kk = kk + (k-1)/strnsiz + 1
C       Doesn't count if embedded in a different name
        if (mod(k-1,strnsiz) /= 0) cycle
        naml = nam
C       A match if the entire names are the same
        if (naml == anam(kk)) then
          k = kk
          exit
        endif
      enddo

C ... Name not present AND size < minsiz ... too small
      if (k == 0 .and. (n1*MB)*n2 < minsiz) then
        return
      endif

C ... Name not present: add new name
      if (k == 0) then
        if (nnam == mxnam) call rx('idalloc: increase mxnam')
        nnam = nnam+1
        anam(nnam) = nam
        asize(nnam) = 0
        k = nnam
      endif

C ... Memory management log function on the k'th entry
      select case (mod(add,10))
        case (0)

        case (1)
          asize(k) = dble(n1)*n2
          xx = asize(k)
          if (add >= 10) then
C#ifdef AWRITE
            call info2(10,0,0,' idalloc: '//
     .      'allocate array '//anam(k)//'%;9,1D MB',xx*MB,0)
C#elseC
C            write(*,103) 'allocate',anam(k),xx*MB
C  103     format(' idalloc: ',a,' array ',a,f9.1,' MB')
C#endif
          endif

        case (2)
          if (asize(k) > 0)
     .    call rx('idalloc: array  '//anam(k)//'%a  already allocated')
          asize(k) = dble(n1)*n2
          xx = asize(k)
          if (add >= 10)  then
C#ifdef AWRITE
            call info2(10,0,0,' idalloc: '//
     .      'allocate array '//anam(k)//'%;9,1D MB',xx*MB,0)
C#elseC
C            write(*,103) 'allocate',anam(k),xx*MB
C#endif
          endif

        case (3)
          xx = abs(asize(k))
          if (add >= 10 .and. asize(k) > 0) then
C#ifdef AWRITE
            call info2(10,0,0,' idalloc: '//
     .      ' release array '//anam(k)//'%;9,1D MB',xx*MB,0)
C#elseC
C            write(*,103) ' release',anam(k),xx*MB
C#endif
          endif
          asize(k) = -xx

        case (4)
          xx = abs(asize(k))
          if (add >= 10) then
C#ifdef AWRITE
            call info2(10,0,0,' idalloc: '//
     .      ' release array '//anam(k)//'%;9,1D MB',xx*MB,0)
C#elseC
C            write(*,103) ' release',anam(k),xx*MB
C#endif
          endif
          if (asize(k) <= 0)
     .    call rx('idalloc: array  '//anam(k)//' not allocated')
          asize(k) = -xx

        case default
          call rx('invalid parameter add')

      end select

      idalloc = xx/1d6 * 8

C
      xx = 0            ! Keep track of max allocation
      do  k = 1, nnam
        if (asize(k) > 0) xx = xx + asize(k)
      enddo
      mxalloc = max(mxalloc,xx)
      return

      entry idasiz(m)
C- Set minimum array size minsiz

      minsiz = m

      end

      integer function allocvb()
C- Return the allocation verbosity
      implicit none
      integer iprint
      integer,save :: ivb=-1
      integer set_allocvb,i
      logical cmdopt
      character strn*(20)

C ... Never set ... generate the verbosity
      if (ivb == -1) then
        ivb = 0
        if (cmdopt('--pra',5,0,strn)) then
          ivb = 10
        endif
      endif

      if (iprint() > 0) then
        allocvb = ivb
      else
        allocvb = 0
      endif

      return

      entry set_allocvb(i)
      if (i >= 0) ivb = min(10*i,10)
      set_allocvb = ivb

      end

      logical function memcheck(name,errmsg,maxmem,lreqd)
C- Check whether dynamic allocations logged by idalloc fall within maximum
C ----------------------------------------------------------------------
Ci Inputs
Ci   name  : name of calling routine
Ci   errmsg:error message if insufficient memory
Ci   maxmem:maximum memory, in MB
Ci   lreqd :If T, aborts if maximum memory exceeded
Ci         :If F, returns with memcheck=F
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Aug 14  First created
C ----------------------------------------------------------------------
C ... Passed parameters
      character name*(*),errmsg*(*)
      integer maxmem
      logical lreqd
C ... Local parameters
      integer i,idalloc
      character lstrn*(256)

      memcheck = .true.
      if (maxmem <= 0) return
      i = idalloc(' ',0,1,1)
      if (i <= maxmem) return

      memcheck = .false.
      lstrn = ' ' // name // ':  need %i MB for '// errmsg //
     .  ' but limited to %i'
      call info2(10,1,0,lstrn,i,maxmem)
      if (.not. lreqd) return
      call rx(name//' aborting because of insufficient memory')

      end

C      subroutine fmain
CC- Tests idalloc
C      integer idalloc,i
C
CC     i =  idalloc(' ',0,1,1)
C      i =  idalloc('abc',10+2,19000,3900)
C
Cc     i = idasiz(1)
C      i =  idalloc('c',0+2,433,300)  ! Doesn't do anything --- too small
C
C      i =  idalloc('c',10+1,4333,3*300)
C      i =  idalloc('c',10+4,4333,3*300)
C      i =  idalloc('c',10+3,4333,3*300) ! No effect; already deallocated
C
C      i =  idalloc('bc',0+0,1,1)
C      i =  idalloc('abcd',10+0,1,1)
C      i =  idalloc('a',10+0,1,1)
C      call pshpr(1)
C      i =  idalloc('abcde',10+1,12345,54321)
C      i =  idalloc('abcde',10+4,12345,54321)
C      call poppr
C
C
C      i =  idalloc('bc',10+1,1000,2000)
C      i =  idalloc('a',20+1,80000,2000)
C      i =  idalloc('abc',10+4,19000,3900)
C      i =  idalloc(' ',10,1,1)
C      print "(' idalloc returned',i6)", i
C
C      end
