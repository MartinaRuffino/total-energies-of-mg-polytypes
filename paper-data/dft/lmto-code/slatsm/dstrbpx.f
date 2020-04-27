      subroutine dstrbpx(nloop,nproc,wt,kpproc)
C- Distribute loop counter contiguously over threads,
C- where each element has a unique execution time
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nloop: number of times the loop is repeated
Ci   wt   : CPU weighting time for each element in the loop
C         : wt(1) = -1 -> all elements get unit weight
Ci   nproc: number of processors
Co Outputs:
Co  kpproc: upper and lower loop bounds for each thread
Co        : kpproc(procid) = starting kpproc for thread procid;
Co        : kpproc(procid+1)-1 = ending kpproc, i.e.
Co        : do  i = kpproc(procid), kpproc(procid+1)-1
Cr Remarks
Cr   For nested loops, create a separate index array, as in the example below
Cr   This example distributes a pair of nested loops:
Cr   the outer loop runs i=1:nloop and the inner j=i:nloop
Cr   But any nested loop does equally well where provided elements are independent
C      nnest = (nloop*(nloop+1))/2     ! Number of elements in the nested loop
C      allocate(index(nnest,2))
C      call dstrbpx(nnest,nproc,[-1d0],kpproc)
C      ! Make the index array. Depends on the loop structure but inest
C      inest = 0
C      do  i = 1, nloop
C      do  j = i, nloop
C          inest = inest+1
C          index(inest,1) = i; index(inest,2) = j
C        enddo
C      enddo
C      if (inest /= nnest) call rx('bug in setting up nested loop')
C
C      do  inest = kpproc(procid), kpproc(procid+1)-1
C        i = index(inest,1); j = index(inest,2)
C        print *, procid, inest, i, j
C        Do block for element (i,j) here
C      enddo
Cu Updates
Cu   17 Oct 17 Allows wt(1)=-1 => unit weights
Cu   27 Jun 14 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: nloop,nproc,kpproc(0:nproc)
      real(8), target :: wt(nloop)
C ... Local parameters
      integer procid,i,j
      double precision wtot,dsum,wtarg,wcount
      real(8), pointer :: wtl(:)
      real(8), allocatable :: wtproc(:)

C ... Setup, including case nproc == 1
      if (nloop < nproc .or. nproc < 1) return

      if (wt(1) == -1) then
        allocate(wtl(nloop))
        call dvset(wtl,1,nloop,1d0)
      else
        wtl => wt
      endif

      wtot = dsum(nloop,wtl,1) ! Total weight for all loop elements
      allocate(wtproc(0:nproc-1))
      wtproc(0) = wtot; kpproc(0) = 1; kpproc(nproc) = nloop+1
      if (nproc == 1) goto 99

C ... Distribute load by finding points closest to average
      procid = 0        ! To be incremented in the loop below
      wtarg = wtot*(procid+1)/nproc ! Running target for ith processor
      wcount = 0        ! Running total weight for current procid
      do  i = 1, nloop-1
        wcount = wcount + wtl(i)  ! Accumulated weight to this loop element
        if (abs(wcount-wtarg) >= abs(wcount+wtl(i+1)-wtarg)) cycle
C       Loop element i has weight closest to target
        procid = procid+1
        kpproc(procid) = i+1
        wtarg = wtot*(procid+1)/nproc
      enddo
      if (procid /= nproc-1) call rx('bug in dstrbpx')

C ... Get weights for each processor
      do  procid = 0, nproc-1
        i = kpproc(procid); j = kpproc(procid+1)-1
        wtproc(procid) = dsum(j-i+1,wtl(i),1)
C#ifdefC DEBUG
C        print *, procid, i, j, sngl(wtproc(procid))
C#endif
      enddo
      wtarg = maxval(wtproc)

C ... Reduce procid with maximum weight, if possible
   10 continue
      call idmxmn(nproc,wtproc,1,i,j)
      procid = j-1; i = kpproc(procid)

C#ifdefC DEBUG
C      print *, 'max weight at', procid,sngl(wtproc(procid))
C#endif

      if (procid > 0) then !  Increment left boundary?
      if (wtproc(procid-1)+wtl(i) < wtproc(procid) .and.
     .    wtl(i) > 0) then
        kpproc(procid) = kpproc(procid)+1
        wtproc(procid-1) = wtproc(procid-1)+wtl(i)
        wtproc(procid) = wtproc(procid)-wtl(i)
        goto 10
      endif
      endif

      if (procid < nproc-1) then !  Increment right boundary?
      i = kpproc(procid+1)-1
      if (wtproc(procid+1)+wtl(i) < wtproc(procid) .and.
     .    wtl(i) > 0) then
        kpproc(procid+1) = kpproc(procid+1)-1
        wtproc(procid+1) = wtproc(procid+1)+wtl(i)
        wtproc(procid) = wtproc(procid)-wtl(i)
        goto 10
      endif
      endif

C#ifdefC DEBUG
C      print *, '!!'
C      do  procid = 0, nproc-1
C        i = kpproc(procid); j = kpproc(procid+1)-1
C        wtproc(procid) = dsum(j-i+1,wtl(i),1)
C        print *, procid, i, j, sngl(wtproc(procid))
C      enddo
C      call idmxmn(nproc,wtproc,1,i,j)
C      procid = j-1; i = kpproc(procid)
C      print *, 'max weight recomputed', procid,sngl(wtproc(procid))
C#endif

C ... For each procid, balance load with left boundary shift
      do  procid = 1, nproc-1
        i = kpproc(procid)
        if (wtl(i) <= 0) cycle
        if (wtproc(procid-1)+wtl(i) < wtproc(procid)) then
          kpproc(procid) = kpproc(procid)+1
          wtproc(procid-1) = wtproc(procid-1)+wtl(i)
          wtproc(procid) = wtproc(procid)-wtl(i)
        endif
      enddo

C#ifdefC DEBUG
C      print *, '!!'
C      do  procid = 0, nproc-1
C        i = kpproc(procid); j = kpproc(procid+1)-1
C        wtproc(procid) = dsum(j-i+1,wtl(i),1)
C        print *, procid, i, j, sngl(wtproc(procid))
C      enddo
C#endif

C ... For each procid, balance load with right boundary shift
      do  procid = nproc-2, 0, -1
        i = kpproc(procid+1)-1
        if (wtl(i) <= 0) cycle
        if (wtproc(procid+1)+wtl(i) < wtproc(procid)) then
          kpproc(procid+1) = kpproc(procid+1)-1
          wtproc(procid+1) = wtproc(procid+1)+wtl(i)
          wtproc(procid) = wtproc(procid)-wtl(i)
        endif
      enddo

C#ifdefC DEBUG
C      print *, '!!'
C      do  procid = 0, nproc-1
C        i = kpproc(procid); j = kpproc(procid+1)-1
C        wtproc(procid) = dsum(j-i+1,wtl(i),1)
C        print *, procid, i, j, sngl(wtproc(procid))
C      enddo
C      call idmxmn(nproc,wtproc,1,i,j)
C      procid = j-1; i = kpproc(procid)
C      print *, 'max weight recomputed', procid,sngl(wtproc(procid))
C#endif

   99 continue
      i = 30 ; if (nproc == 1) i = 60
      call info5(i,0,0,' dstrbpx: %i elts, %i processors.  '//
     .  'wtarg=%;3g  wmax=%;3g',nloop,nproc,wtot/nproc,maxval(wtproc),0)
      deallocate(wtproc)
      if (wt(1) == -1) deallocate(wtl)

      end

C     Test
C      subroutine fmain
C      implicit none
C      integer, parameter :: nproc = 4, nloop = 11
CC      integer, parameter :: nproc = 11, nloop = 355
Cc      integer, parameter :: nproc = 41, nloop = 982
C
CC      integer, parameter :: nproc = 49, nloop = 1015 ! increment left
CC      Should generate this output:
CC dstrbpx: 1015 elts, 49 processors.  wtarg=10.2  wmax=10.7
CC      1    21    43    64    85   108
CC    129   151   181   202   223   239
CC    257   276   296   322   343   363
CC    379   398   419   443   461   483
CC    500   518   538   558   581   604
CC    627   647   667   687   711   734
CC    751   769   791   810   831   848
CC    870   888   912   932   950   973
CC    994  1016
CC Change to unit weights.  Should get
CC  dstrbpx: 1015 elts, 49 processors.  wtarg=20.7  wmax=21
CC      1    22    42    63    84   105
CC    125   146   167   187   208   229
CC    250   270   291   312   332   353
CC    374   395   415   436   457   477
CC    498   519   540   560   581   602
CC    622   643   664   685   705   726
CC    747   767   788   809   830   850
CC    871   892   912   933   954   975
CC    995  1016
C
CC     integer, parameter :: nproc = 53, nloop = 1988 ! increment left,right
CC      integer, parameter :: nproc = 1, nloop = 11
C
C      real ran1
C      double precision wt(nloop)
C      integer i,j,inest,nnest,procid,kpproc(0:nproc)
C      integer, allocatable :: index(:,:)
C
CC --- Interface Required for f90 compatibility ---
C      interface
C      subroutine dstrbpx(nloop,nproc,wt,kpproc)
C      implicit none
CC ... Passed parameters
C      integer :: nloop,nproc,kpproc(0:nproc)
C      real(8), target :: wt(nloop)
C      end
C      end interface
C
C      goto 999
C
CC     Basic check
C      call ran1in(11)
C      do  i = 1, nloop
C        wt(i) = 0.5d0
C        wt(i) = ran1()
CC        wt(i) = nint(ran1())
C      enddo
CC     Random weights
C      call dstrbpx(nloop,nproc,wt,kpproc)
C      print 333, kpproc
CC     Unit weights
C      call dstrbpx(nloop,nproc,[-1d0],kpproc)
C      print 333, kpproc
C  333 format(6i6)
C
C  999 continue
C
C      nnest = (nloop*(nloop+1))/2
C      allocate(index(nnest,2))
C      call dstrbpx(nnest,nproc,[-1d0],kpproc)
C      inest = 0
C      do  i = 1, nloop
C      do  j = i, nloop
C          inest = inest+1
C          index(inest,1) = i; index(inest,2) = j
C        enddo
C      enddo
C      if (inest /= nnest) call rx('bug in setting up nested loop')
C
C      do  procid = 0, nproc-1
C      do  inest = kpproc(procid), kpproc(procid+1)-1
C        i = index(inest,1); j = index(inest,2)
C        print 333, procid, inest, i, j
C      enddo
C      enddo
C
C      end
