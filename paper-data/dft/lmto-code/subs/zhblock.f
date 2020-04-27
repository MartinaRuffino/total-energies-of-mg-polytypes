      subroutine zhblock(sopts,ndimh,isp,h,s)
C- Scales specified off-diagonal elements of Hamiltonian and Overlap
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :character string containing list of rows,columns
Ci         :Options in sopts separated by delimiter, which is 1st char
Ci   ndimh :dimension of hamiltonian and overlap h,s
Ci   isp   :current spin channel (1 or 2), or 0 if spin off-diagonal
Ci         :Add 10s digit to flag that s is missing: scale out h only
Ci         :Add 100s digit to flag off-diagonal block => e=0
Ci   h     :Hamiltonian matrix
Ci   s     :overlap matrix
Cl Local variables
Cr Remarks
Cr   Off-diagonal elements of h and s, specified by sopts, are scaled.
Cr   sopts consists of a sequence of blocks; each block specifies
Cr   one or more of the following:
Cr      list of rows
Cr      list of columns
Cr      energy of diagonal element e (optional; if not specified, left untouched)
Cr      scale factor scl (optional; if not specified; set to zero)
Cr
Cr   Delimiters separate each specification.  The delimiter is the first
Cr   character in sopts.  We assume in the description below that it is '~'.
Cr   You may specify multiple blocks.  Separate blocks by a double-delimiter '~~'.
Cr
Cr   For a single block, a list of rows and columns must be specified;
Cr   but there is some flexibility in how this is accomplished.
Cr   There are the following choices:
Cr      1. ~rows=(integer-list)
Cr          For each i and j in integer-list, each h(i,j) is scaled by scl for i ne. j
Cr          For i=j, h and s are left untouched unless e is given.
Cr          In that case s(i,i) is left untouched and h(i,i) = e*s(i,i)
Cr      2. ~cols=(integer-list)
Cr          Has the same effect as 1.
Cr      3. ~rows=(integer-list1)~cols=(integer-list2)
Cr          Similar to 1, but i must belong to integer-list1, j to integer-list2.
Cr      4. ~(integer-list)
Cr          Similar to 3, but integer-list2 is the entire hamiltonian
Cr      5. ~(integer-list-1)@(integer-list-2)
Cr          This is a variant of version 4. The first integer list is applied in the
Cr          spin-1 case; the second in the spin-2 case.
Cr   Example 1:
Cr     ~e=-1.5,1.4~10:16
Cr   This is an instance of option 4.
Cr   h(10:16,:) and h(:,10:16), and s(10:16,:) and s(:,10:16) are zeroed out
Cr   except for the diagonal elements.  The diagonal elements of s are left untouched
Cr   The diagonal elements of h(i,i) for i=10:16 are set to -1.5Ry*s(i,i)
Cr   for the majority spin and +1.4Ry*s(i,i) for the minority spin.  As a result
Cr   7 eigenvalues will be -1.5 Ry for spin 1, and +1.4 Ry for spin 2.
Cr   Example 2:
Cr     ~rows=5,6,8~s=.5~~rows=7,9~s=0
Cr   Scales the off-diagonal elements between t2g orbitals (5,6,8) by 0.5,
Cr   and sets the off-diagonal elements between eg orbitals to zero.
Cr   Example 3:
Cr     ~rows=5,6,8~cols=7,9~s=.5
Cr   Scales the matrix elements connecting t2g and eg orbitals by 1/2, leaving
Cr   the intra-t2g and intra-eg matrix elements intact.
Cr   Example 5:
Cr    ~e=-0.5~40,41,77@40,77,78~~e=0.5~42:46,78:83@41:46,79:83'
Cr    Zeros out off-diagonal elements for spin 1:
Cr      40,41,77 at energy = -0.5   and 42:46,78:83 at energy 0.5
Cr    For spin 2, zeros out off-diagonal elements
Cr      40,77,78 at energy = -0.5   and 41:46,79:83 at energy 0.5
Cr    Thus for each spin, 3 electrons are uncoupled and shifted down and 11 uncoupled and shifted up
Cr    but the orbitals entering into the shift depend on spin.
Cr    This example is useful for f states on a system with 2 Nd atoms.
Cr    The first acquires spin +1 from the f for the first atom and spin -1 for the second.
Cr   Further notes.
Cr   *For each (i,j) modified, (j,i) is modified in a the same manner
Cr    so that h and s remain hermitian.
Cr   *e can be one or two numbers, corresponding to two spins.
Cr    If the second number is not specified e(1) is used for both spins
Cr   *scl can be one or two numbers, corresponding to hamiltonian and overlap
Cr    If neither number is specified, s=0 for both h and s
Cr    If the second number is not specified, elements in the overlap are set to zero.
Cr   *Caution: since both h and its hermitian conjugate are scaled,
Cr    elements of h affected that belong to both integer-list1 and integer-list2
Cr    get scaled twice.
Cr   *s is not touched at all if isp<0.
Cr   *For the syntax of integer-lists, see slatsm/mkilst.f.
Cr   *If no list is specified, the zhblock returns without modifying anything.
Cr
Cu Updates
Cu   12 Oct 17 Added spin-dependent orbital list
Cu   08 Jan 17 Adapted from hblock.f
Cu    1 Jun 10 hblock first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character sopts*(*)
      integer ndimh,isp
      double complex h(ndimh,ndimh),s(ndimh,ndimh)
C ... Dynamically allocated arrays
      integer, pointer :: irow(:),icol(:)
C ... Local parameters
      character dc*1
      logical lsete,leob,HaveS,ZeroE
      integer i,j,k,nrow,ncol,ix(3),j1,j2,ispl,j1x,j2x
      double precision e(3),scl(2)
      double complex cxx
      character strn*256
      procedure(integer) :: a2vec,mkilsd,iprint,isw

C     debugging
C     forall (i=1:ndimh,j=1:ndimh) h(i,j) = i+10*j
C     call zprm('hblock starting h',12,h,ndimh,ndimh,ndimh)
C     call zprm('hblock starting s',12,s,ndimh,ndimh,ndimh)

      lsete = .false.
      ispl = mod(isp,10)
      HaveS = mod(isp/10,10) == 0
      ZeroE = mod(isp/100,10) /= 0
      dc = sopts(1:1)
      if (dc == ' ') return
      j2 = 0
C ... Next block
   10 continue
      nrow = -1 ; ncol = -1; scl = 0; leob = .false.
C ... Loop until rows and columns are both available
      do
        if (leob) exit ! End of block has been flagged
        j2 = j2+1
        leob = sopts(j2:min(len(sopts),j2+1)) == dc//dc
        if (sopts(j2:j2) == dc) cycle
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 < j1) exit
        if (sopts(j1:j1+1) == 'e=')  then
          j1 = j1+1
          k = a2vec(sopts,j2,j1,4,dc//', ',3,3,2,ix,e)
          if (k <= 0) goto 999
          if (k == 2 .and. ispl == 2) e(1) = e(2)
          if (ZeroE) e(1) = 0
          lsete = .true.
          cycle
        else if (sopts(j1:j1+1) == 's=')  then
          j1 = j1+1
          k = a2vec(sopts,j2,j1,4,dc//', ',3,3,2,ix,scl)
          if (k <= 0) goto 999
          if (k == 1) scl(2) = 0
        else if (sopts(j1:j1+4) == 'rows=')  then
          nrow = mkilsd(sopts(j1+5:j2),-1,ix)
          if (nrow <= 0) goto 999
          allocate(irow(nrow))
          nrow = mkilsd(sopts(j1+5:j2),nrow,irow)
        else if (sopts(j1:j1+4) == 'cols=')  then
          ncol = mkilsd(sopts(j1+5:j2),-1,ix)
          if (ncol <= 0) goto 999
          allocate(icol(ndimh))
          ncol = mkilsd(sopts(j1+5:j2),ncol,icol)
C       Only remaining option is an integer list
        else
C         Check whether list is spin-dependent
          j1x = j1; j2x = j2
          j = index(sopts(j1:j2),'@')
          if (j == 0) then
          elseif (ispl == 1) then
            j2x = j1+j-2
          else
            j1x = j1+j
          endif
          nrow = mkilsd(sopts(j1x:j2x),-1,ix)
          if (nrow < 0) goto 999
          allocate(irow(nrow),icol(ndimh))
          nrow = mkilsd(sopts(j1x:j2x),nrow,irow)
          forall (j=1:ndimh) icol(j) = j
          ncol = ndimh
        endif
      enddo

      if (nrow < 0 .and. ncol < 0) return ! Nothing to do
      if (nrow < 0) then
        nrow = ncol
        irow => icol
      endif
      if (ncol < 0) then
        ncol = nrow
        icol => irow
      endif

      call info8(30,0,0,' zhblock spin %i: scaling %i rows, %i columns by factor %d.%?;n;  e=%d;;',
     .  ispl,nrow,ncol,scl,isw(lsete),e,7,8)
      if (iprint() >= 40) then
        call ilst2a(irow,nrow,strn)
        call info0(20,0,0,' rows : '//trim(strn))
        call ilst2a(icol,ncol,strn)
        call info0(20,0,0,' cols : '//trim(strn))
      endif

C     Scale the blocks
      do  k = 1, nrow
        i = irow(k)
        if (i < 1 .or. i > ndimh) call rxi('zhblock: illegal element :',i)
        cxx = h(i,i)
        forall (j=1:ncol) h(icol(j),i) = h(icol(j),i) * scl(1)
        forall (j=1:ncol) h(i,icol(j)) = h(i,icol(j)) * scl(1)
        h(i,i) = cxx
        if (HaveS) then
          cxx = s(i,i)
          forall (j=1:ncol) s(icol(j),i) = s(icol(j),i) * scl(2)
          forall (j=1:ncol) s(i,icol(j)) = s(i,icol(j)) * scl(2)
          s(i,i) = cxx
        endif
        if (lsete) then
          h(i,i) = e(1)*s(i,i)
        endif
      enddo

C     call zprm('hblock ending h',12,h,ndimh,ndimh,ndimh)
C     call zprm('hblock ending s',12,s,ndimh,ndimh,ndimh)

      if (associated(irow,icol)) then
        deallocate(irow)
      else
        deallocate(irow,icol)
      endif
      goto 10  ! Start next block

      return
  999 call rxs('zhblock: failed to parse options : ',sopts(j1:j2))

      end
