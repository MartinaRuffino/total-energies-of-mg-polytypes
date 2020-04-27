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
Ci         :Add a multiple of 10 to print table to stdout
Ci   n1    :(add = 1,2 only) new array to be allocated as arr(n1,n2)
Ci   n2    :(add = 1,2 only) new array to be allocated as arr(n1,n2)
Co Outputs
Co idalloc :
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr  An internal is associated with 'nam'.
Cr  If nam is not already in the current list, it is added to the list
Cu Updates
Cu   13 May 13  First created
C ----------------------------------------------------------------------
      implicit none
      character nam*(*)
      integer add,n1,n2
C Local
      integer mxnam,strnsiz
      parameter (mxnam=30,strnsiz=8)
      integer,save :: nnam=0
      real(8),save :: asize(mxnam)
      real(8),save :: xx,mxalloc
      character*(strnsiz*mxnam),save :: namlst
      character*(strnsiz) :: anam(mxnam),naml
      equivalence (anam,namlst)
      integer k,kk

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

          write(*,101) kk,xx/1d6*8,mxalloc/1d6*8
  101     format(/' idalloc:',i3,' arrays allocated for',f9.1,' MB.',
     .      '  Max alloc:',f9.0)
          do  k = 1, nnam
            if (asize(k) > 0) then
              write(*,102) anam(k),asize(k)/1d6*8
  102         format(3x,a,f9.1:a,' when last allocated')
            elseif (asize(k) < 0) then
              write(*,102) anam(k),asize(k)/1d6*8,'*'
            else
              write(*,"(3x,a,7x,'--')") anam(k)
            endif
          enddo

        endif

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
C ... Name not present: add new name
      if (k == 0) then
        if (nnam == mxnam) call rx('idalloc: increase mxnam')
        nnam = nnam+1
        anam(nnam) = nam
        asize(nnam) = 0
        k = nnam
      endif

C ... Do out memory management log function on the k'th entry
      select case (mod(add,10))
        case (0)

        case (1)
          asize(k) = dble(n1)*n2
          xx = asize(k)
C          if (add >= 10) call info2(10,0,0,' idalloc: '//
C     .      'allocate array '//anam(k)//'%a : %;1,1d MB',xx/1d6*8,0)
          if (add >= 10) write(*,103) 'allocate',anam(k),xx/1d6*8
  103     format(' idalloc: ',a,' array ',a,f9.1,' MB')

        case (2)
          if (asize(k) > 0)
     .    call rx('idalloc: array  '//anam(k)//' already allocated')
          asize(k) = dble(n1)*n2
          xx = asize(k)
          if (add >= 10) write(*,103) 'allocate',anam(k),xx/1d6*8

        case (3)
          xx = abs(asize(k))
C          if (add >= 10) call info2(10,0,0,' idalloc: '//
C     .      'release array '//anam(k)//'%a : %;1,1d MB',xx/1d6*8,0)
          if (add >= 10) write(*,103) ' release',anam(k),xx/1d6*8
          asize(k) = -xx

        case (4)
          xx = abs(asize(k))
C          if (add >= 10) call info2(10,0,0,' idalloc: '//
C     .      'release array '//anam(k)//'%a : %;1,1d MB',xx/1d6*8,0)
          if (add >= 10) write(*,103) ' release',anam(k),xx/1d6*8
          if (asize(k) <= 0)
     .    call rx('idalloc: array  '//anam(k)//'%a  not allocated')
          asize(k) = -xx

        case default
          call rx('idalloc: invalid parameter add')

      end select

      idalloc = xx/1d6 * 8

C
      xx = 0            ! Keep track of max allocation
      do  k = 1, nnam
        if (asize(k) > 0) xx = xx + asize(k)
      enddo
      mxalloc = max(mxalloc,xx)

      end

C      integer function allocvb
CC- Return the allocation verbosity
C      implicit none
C      integer iprint
C      integer,save :: ivb=-1
C      logical cmdopt
C      character strn*(20)
C
CC ... Never set ... generate the verbosity
C      if (ivb == -1) then
C        ivb = 0
C        if ( cmdopt('--time=',7,0,strn) .or. iprint() >= 50) then
C          ivb = 10
C        endif
C      endif
C
C      allocvb = ivb
C
C      end

C      subroutine fmain
CC- Tests idalloc
C      integer idalloc,i
C
CC     i =  idalloc(' ',0,1,1)
C      i =  idalloc('abc',10+2,19000,3900)
C      i =  idalloc('c',0+2,433,300)
C      i =  idalloc('c',10+1,433,3*300)
C      i =  idalloc('c',10+4,433,3*300)
CC     i =  idalloc('c',10+3,433,3*300)
C      i =  idalloc('bc',0+0,1,1)
C      i =  idalloc('abcd',10+0,1,1)
C      i =  idalloc('a',10+0,1,1)
C      i =  idalloc('abcde',10+0,1,1)
C
C      i =  idalloc('bc',10+1,1000,2000)
C      i =  idalloc('a',20+1,80000,2000)
C      i =  idalloc('abc',10+4,19000,3900)
C      i =  idalloc(' ',10,1,1)
C      print "(' idalloc returned',i6)", i
C
C      end
