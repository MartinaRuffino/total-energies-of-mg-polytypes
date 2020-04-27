      subroutine thdist2(mode,a,b,nthreads,ls1,le1,ls2,le2,
     .  first,last,idx)
C- Determine loop indices for a pair of perfectly nested multithreaded do loops
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 load is evenly balanced between for each loop element
Ci         :  (compound of inner and outer loop)
Ci         :1 load for loop element j = a*j+b, j = index to inner loop
Ci   a,b   :  load for loop element j = a*j+b
Ci   nthreads : number of threads
Ci   ls1   :  starting index for outer loop
Ci   le1   :  ending index for outer loop
Ci   ls2   :  starting index for inner loop
Ci   le2   :  ending index for inner loop
Co Outputs (See Remarks)
Co   first :  first(i) = starting (compound) loop index for thread i
Co   last  :  last(i)  = ending   (compound) loop index for thread i
Co   idx   :  idx(1,k) = outer loop index for compound loop index k
Co         :  idx(2,k) = inner loop index for compound loop index k
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   Perfectly nested i,j loop = > starting indices of both are fixed.
Cr
Cr   thdist2 generates enables the nested loop to be treated as
Cr           a single large loop, so that the division of labor
Cr           can be optimized over the loops.
Cr   Calling thdist2 generates:
Cr     idx : recover i, j for a compound loop index element
Cr     first, last :
Cr           starting and ending points for compound loop to
Cr           be partitioned into threads
Cr
Cr   Use outputs from thdist2 in a loop like this:
Cr     C$OMP parallel private(m,k,i,j)
Cr           m = omppid(1)
Cr           do  k = first(m), last(m) ! compound index
Cr             i = idx(1,k)            ! outer loop index
Cr             j = idx(2,k)            ! inner loop index
Cr             ...
Cr           enddo
Cr     C$OMP end parallel
Cr
Cu Updates
Cu   26 May 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nthreads,ls1,le1,ls2,le2
      integer first(nthreads),last(nthreads)
      integer idx(2,(le1-ls1+1)*(le2-ls2+1))
      double precision a(*),b
C ... Required for f90
      INTERFACE
      subroutine thdist1(mode,a,b,nthreads,lstart,lend,first,last)
      integer mode,nthreads,lstart,lend
      integer first(nthreads),last(nthreads)
      real(8),target :: a(lstart:lend)
      double precision b
      end subroutine thdist1
      END INTERFACE
C ... Local parameters
      logical getload
      integer i,j,nblk,ii
      real(8),pointer :: lload(:)

      nblk = (le1-ls1+1) * (le2-ls2+1)
      first = 1; last = 0

C ... Case where there are no threads or no elements in outer loop
      if (nthreads <= 1 .or. nblk <= 0) then
        first(1) = 1
        last(1)  = nblk
        if (nblk == 0) return
      endif

C ... Construct separate loop indices from compound index
C     and load table, if mode>0
      ii = 0
      if (mode == 0) then
        allocate(lload(1))
      elseif (mode == 1) then
        allocate(lload(nblk))
      else
        call rx('thdist2: bad mode')
      endif
      getload = mode == 1 .and. nthreads > 1
      do  i = ls1, le1
        do  j = ls2, le2
          ii = ii+1
          idx(1,ii) = i
          idx(2,ii) = j
          if (getload) lload(ii) = a(1)*iabs(j)+b
        enddo
      enddo
      if (ii /= nblk) call rx('bug in thdist2')
      if (nthreads <= 1) return

C ... Distribute compound loop among threads
      if (mode == 0) then
        call thdist1(0,a,b,nthreads,1,ii,first,last)
      else
C       print *, 'sum load',sum(lload)
        call thdist1(2,lload,b,nthreads,1,ii,first,last)
      endif

      end

      subroutine thdist1(mode,a,b,nthreads,lstart,lend,first,last)
C- Determine starting, ending indices for a single multithreaded do loop
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 load is evenly balanced between for each loop element j,
Ci         :  j ranges between lstart:lend
Ci         :1 load for loop element j = a*|j|+b
Ci         :2 load for loop element j = a(j)
Ci   a     :  (mode=0) not used
Ci         :  (mode=1) load for loop element j = a*|j|+b
Ci         :  (mode=2) load for loop element j = a(j)
Ci   b     :  (mode=0) not used
Ci         :  (mode=1) load for loop element j = a*|j|+b
Ci         :  (mode=2) not used
Ci   nthreads : number of threads
Ci   lstart:  starting point in loop
Ci   lend  :  ending point in loop
Co Outputs
Co   first :  first(i) = starting loop index for thread i
Co   last  :  last(i)  = ending   loop index for thread i
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   26 May 13 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nthreads,lstart,lend
      integer first(nthreads),last(nthreads)
      real(8),target :: a(lstart:lend)
      double precision b
C ... Local parameters
      integer i,j,nmisfit,nx,nblk,nblki,nstart
      double precision load,loadj,lastlj,loadi,loadii,spread,bload,iload
      real(8),pointer :: lload(:)

      first = 1; last = 0

C ... Case where there are no threads or no elements in outer loop
      if (nthreads <= 1 .or. lend <= lstart) then
        first(1) = lstart
        last(1)  = lend
        return
      endif

      nblk = lend-lstart+1
      nblki = nblk/nthreads ! 0 if nthreads>nblk
      nmisfit = mod(nblk,nthreads)
      nstart = lstart
      nx = 1; if (nmisfit <= 0) nx=0

C ... Case load is evenly balanced among threads
      if (mode == 0 .or. nthreads >= nblk) then
        do  i = 1, nthreads
          first(i) = nstart
          last(i)  = nstart-1 + nblki + nx
          if (last(i) == lend) exit
          nmisfit  = nmisfit-1
          if (nmisfit <= 0) nx=0
          nstart = last(i)+1
        enddo
      elseif (mode == 1 .or. mode == 2) then
        if (mode == 1) then
          load = 0              ! Will be average load per thread
          allocate(lload(lstart:lend))
          do  j = lstart, lend
            lload(j) = a(lstart)*iabs(j) + b ! Load for this loop element
            load = load + lload(j)
          enddo
        else
          lload => a
          do  j = lstart, lend
            load = load + lload(j)
          enddo
        endif
        if (load <= 0) call rx('thdist1: nonsensical a,b')
        load = load/nthreads
        loadii = 0              ! total load accumulated up to thread i-1
        loadi  = 0              ! total load accumulated up to thread i
        loadj  = 0              ! total load up loop index j
        do  i = 1, nthreads
          first(i) = nstart     ! loop starting point for thread i
          iload = i*load        ! Target summed load up to thread i

C         Find terminus for this thread
          do  j = nstart, lend
            lastlj = loadj      ! Load from prior j
            loadj = loadj + lload(j) ! Total load including loop element j
            last(i) = j         ! In case this is the terminus for this thread

C           Cases where this j should be loop terminus for thread i
            if (j == lend .and. i == nthreads) exit
            if (loadj >= iload .or. nthreads-i >= lend-j) then
              if (j == nstart) exit
C             Check whether the prior j balances load better
C             Must be closer to target load for this i
              if (abs(loadj-iload) > abs(lastlj-iload)) then
C               Load for prior thread should meet target, or
C               load for this thread should be at least average load
                if (loadi >= iload-load .or. lastlj-loadi >= load) then
                  last(i) = j-1
                  loadj = loadj - lload(j)
                endif
              endif
              exit
            endif
          enddo
C         Check whether increment of (prior last, current first) improves balance
          spread = (loadj-loadi)-(loadi-loadii) ! spread in load betw/ (i-1),i
          if (i > 1 .and. spread > 0) then
            bload = lload(first(i)) ! load at the boundary point
C            print "(3i4,2x,6f8.0)", i,first(i),last(i),
C     .        loadi,loadj,loadj-loadi,loadi-loadii,spread,bload
            if (abs(spread-2*bload) < abs(spread)) then
C             print *, 'switch',i,first(i)
              first(i)  = first(i)+1
              last(i-1) = last(i-1)+1
              loadi = loadi+bload
              loadj = loadj-bload
            endif
          endif

          loadii = loadi
          loadi = loadj
          nstart = last(i)+1
        enddo
        if (mode == 1) deallocate(lload)
      else
        call rx('thdist1: bad mode')
      endif

C     if (last(nthreads) /= lend) call rx('bug in thdist1')

      end
