      subroutine ivheap(m,n,vecs,iprm,opts)
C- Heapsort array of integer vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   vecs(m,n): n vectors of length m are to be sorted
Ci   iprm: an integer work array of dimension n, or if vecs returned
Ci        in sorted order, an array of length m*n
Ci   opts: ones digit
Ci           0 vecs returned sorted.
Ci           1 vecs is unchanged; only iprm is returned
Ci         tens digit
Ci           0 vecs sorted
Ci           1 vecs sorted by increasing length
Ci         hundreds digit
Ci           1 equal vectors preserve their original order
Ci
Co Outputs
Co   iprm a permutation table that sorts array 'vecs'
Co   vecs may be changed, depending on opts
Cu Updates
Cu   16 Jul 09 vecs no longer need to be doubly dimensioned for opts=0
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer m,n,iprm(n),opts
      integer vecs(m,n)
C ... Local parameters
      integer di,dl
      integer l,ir,irra,i,j,mm,i1,i2,nn
      integer, allocatable :: iwk(:,:)
      logical norm

C     Required for f90 compatibility
      interface
      subroutine ivprm(m,n,vecs,wk,iprm,opt)
      implicit none
C ... Passed parameters
      integer m,n,iprm(n),opt
      integer vecs(m,n)
      integer, target :: wk(m,n)
      end
      end interface

      do  2  ir = 1, n
    2 iprm(ir) = ir
      if (n <= 1) return
      norm = mod(opts/10,10) /= 0
      l = n/2+1
      ir = n

C --- For each l = n/2+1, 1, -1 do ---
   10 continue
C ... "Hiring phase"
      if (l > 1) then
        l = l-1
        irra = iprm(l)
C ... "Retirement-and-promotion phase"
      else
        irra = iprm(ir)
        iprm(ir) = iprm(1)
*       call awrit3('ir%i: %n:1i',' ',180,6,ir,n,iprm)
        ir = ir-1
        if (ir == 1) then
          iprm(1) = irra
*         call awrit2('exit %n:1i',' ',180,6,n,iprm)
          goto 100
        endif
      endif

C ... Setup to sift down element irra to proper level
      i = l
      j = l+l

C --- Do while j <= ir ---
   20 if (j <= ir) then

C   ... Increment j if vecs(iprm(j+1)) > vecs(iprm(j))
        if (j < ir) then
          if (norm) then
            di = 0
            dl = 0
            do  24  mm = 1, m
            dl = dl + vecs(mm,iprm(j))**2
   24       di = di + vecs(mm,iprm(j+1))**2
            if (di-dl > 0) j = j+1
          else
            do  26  mm = 1, m
            if (abs(vecs(mm,iprm(j+1))-vecs(mm,iprm(j))) <= 0) goto 26
            if (vecs(mm,iprm(j+1))-vecs(mm,iprm(j)) > 0) j = j+1
            goto 28
   26       continue
   28       continue
          endif
        endif

C   ... Demote rra to its level
        if (norm) then
          di = 0
          dl = 0
          do  34  mm = 1, m
          dl = dl + vecs(mm,irra)**2
   34     di = di + vecs(mm,iprm(j))**2
          if (di-dl > 0) then
            iprm(i) = iprm(j)
*           call awrit4('%i,%i: %n:1i',' ',180,6,i,j,n,iprm)
            i = j
            j = j+j
C     ... This is rra's level; set j to terminate the sift-down
          else
            j = ir+1
          endif
        else
          do  36  mm = 1, m
C     ...   Skip over equal elements
            if (abs(vecs(mm,iprm(j))-vecs(mm,irra)) <= 0) goto 36
            if (vecs(mm,iprm(j))-vecs(mm,irra) > 0) then
              iprm(i) = iprm(j)
*             call awrit4('%i,%i: %n:1i',' ',180,6,i,j,n,iprm)
              i = j
              j = j+j
C     ... This is rra's level; set j to terminate the sift-down
            else
              j = ir+1
            endif
            goto 38
   36     continue
C     ... Case rra = vec(iprm(j))
          j = ir+1
   38     continue
        endif
        go to 20
      endif
C ... Put rra into its slot
      iprm(i) = irra
*     call awrit3('%i: %n:1i',' ',180,6,i,n,iprm)
      go to 10

C --- For equal vectors, restore original ordering ---
  100 continue
      if (mod(opts/100,10) == 0) goto 200
      i2 = 0
C ... Find i1,i2 = next range of equal numbers
  110 i1 = i2+1
      if (i1 > n) goto 200
  120 i2 = i2+1
      if (i2 > n) goto 130
      if (norm) then
        di = 0
        dl = 0
        do  124  mm = 1, m
        dl = dl + vecs(mm,iprm(i1))**2
  124   di = di + vecs(mm,iprm(i2))**2
        if (di-dl > 0) goto 130
      else
        do  126  mm = 1, m
  126   if (abs(vecs(mm,iprm(i2))-vecs(mm,iprm(i1))) > 0) goto 130
      endif
C ... vec(i1) = vec(i2) ; imcrement i2 and try again
      goto 120

C --- Sort iprm(i1)..iprm(i2) ---
  130 continue
      i2 = i2-1
      i1 = i1-1
      nn = i2-i1
      if (nn <= 1) goto 110
      l = nn/2+1
      ir = nn

C ... For each l = (i2-i1+1)/2+1, 1, -1 do
  140 continue
C ... "Hiring phase"
      if (l > 1) then
        l = l-1
        irra = iprm(l+i1)
C ... "Retirement-and-promotion phase"
      else
        irra = iprm(ir+i1)
        iprm(ir+i1) = iprm(1+i1)
        ir = ir-1
        if (ir == 1) then
          iprm(1+i1) = irra
          goto 110
        endif
      endif

C ... Setup to sift down element irra to proper level ...
      i = l
      j = l+l

C ... Do while j <= ir ...
  150 if (j <= ir) then

C   ... Increment j if iprm(j+i11) > iprm(j+i1))
        if (j < ir) then
          if (iprm(j+i1) < iprm(j+1+i1)) j = j+1
        endif
C   ... Demote irra to its level
        if (irra < iprm(j+i1)) then
          iprm(i+i1) = iprm(j+i1)
          i = j
          j = j+j
C   ... This is irra's level; set j to terminate the sift-down
        else
          j = ir+1
        endif
        go to 150
      endif
C ... Put rra into its slot
      iprm(i+i1) = irra
      go to 140

C --- Sort vecs ---
  200 continue
      if (mod(opts,10) == 0) then
        allocate(iwk(m,n))
        call ivprm(m,n,vecs,iwk,iprm,1)
        deallocate(iwk)
      endif
      end
      subroutine ivprm(m,n,vecs,wk,iprm,opt)
C- Permute an array of integer vectors according to iprm
C ----------------------------------------------------------------
Ci Inputs
Ci   m     :length of each vector
Ci   n     :number of vectors to be sorted
Ci   vecs  :the array vectors, dimensioned (m,n)
Ci   iprm  :a permutation table by which array vecs is reordered
Ci   opt:  :ones digit
Ci          0 return sorted array in wk
Ci          1 same as 0 but copy wk back to vecs
Ci          2 do not use wk but allocate internal array.
Ci            permuted array returned in vecs.
Ci          tens digit
Ci          0 assume all iprm are nozero
Ci          1 cull list, excluding rows for which iprm is zero
Ci            return in n number of elements returned
Ci Inputs/Outputs
Co   n     : if list is to be culled, n is returned as the number
Co         : of elements for which iprm>0.  Otherwise, n is not touched
Co   wk    :returns vecs in permuted order
Co   vecs  :wk may be copied back into vecs, depending on opt.
Cu Updates
Cu   19 Aug 17 switched logical lopt for integer opt; added 10s digit
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer m,n,iprm(n),opt
      integer vecs(m,n)
      integer, target :: wk(m,n)
C ... Local parameters
      integer i,j,k,l
      integer, pointer :: wkl(:,:)

      if (mod(opt,10) /= 2) then
        wkl => wk
      else
        allocate(wkl(m,n))
        call icopy(n*m,vecs,1,wkl,1)
      endif

      if (mod(opt/10,10) /= 0) then
        l = 0
        do  i = 1, n
          k = iprm(i)
          if (k <= 0) cycle
          l = l+1
          forall (j=1:m) wkl(j,l) = vecs(j,k)
        enddo
        n = l
      else
        do  i = 1, n
          k = iprm(i)
          forall (j=1:m) wkl(j,i) = vecs(j,k)
        enddo
      endif
      if (mod(opt,10) /= 0) then
        do  i = 1, n
          forall (j=1:m) vecs(j,i) = wkl(j,i)
        enddo
      endif

      if (mod(opt,10) == 2) deallocate(wkl)
      end
      subroutine ivdup(m,n,vecs,iprm,opt)
C- Remove duplicates from a table of integer vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   m     :length of each vector
Ci   n     :number of vectors to be sorted
Ci   vecs  :the array vectors, dimensioned (m,n)
Ci   opt   :ones digit
Ci           0 input vecs is ordered vecs(m,n)
Ci           1 input vecs is ordered vecs(n,m)
Ci   opt   :tens digit
Ci           0 return iprm only
Ci           1 return vecs sorted and purged
Co Outputs
Co   vecs  :may be changed, depending on opt
Co   n     :number of vectors after culling (only set if 10s digit opt is nonzero)
Cu Updates
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer m,n,opt,iprm(n)
      integer, target :: vecs(n,m)  ! Indices should be reversed if 1s digit opt=0
C ... Dynamically allocated local arrays
      integer, allocatable :: vecl(:,:)
C ... Local parameters
      integer i,ig,jg,kg

      allocate(vecl(m,n))
      if (mod(opt,10) == 1) then
        forall (i=1:m, ig=1:n) vecl(i,ig) = vecs(ig,i)
      else
        call icopy(m*n,vecs,1,vecl,1)
      endif
C     call ivshel(m,n,vecl,iprm,.false.)
      call ivheap(m,n,vecl,iprm,1)
      do  ig = 2, n
        do  i = 1, m
          if (vecl(i,iprm(ig)) /= vecl(i,iprm(ig-1))) goto 100
        enddo
        iprm(ig) = -1
  100   continue
      enddo

      if (mod(opt/10,10) == 0) goto 900

      if (mod(opt,10) == 0) then
        call rx('ivdup not ready for this branch')
c       call ivprm(m,n,vecs,vecs,iprm,12)
      else
        ig = 0
        do  jg = 1, n
          kg = iprm(jg)
          if (kg <= 0) cycle
          ig = ig+1
          forall (i=1:m) vecs(ig,i) = vecl(i,kg)
        enddo
      endif
      n = ig

  900 continue
      deallocate(vecl)
      end

