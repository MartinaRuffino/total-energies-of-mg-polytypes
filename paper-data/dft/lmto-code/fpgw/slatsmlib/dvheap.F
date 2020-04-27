      subroutine dvheap(m,n,vecs,iprm,tol,opts)
C- Heapsort array of double-precision vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   m     :length of each vector
Ci   n     :number of vectors to be sorted
Ci   vecs  :the array vectors, dimensioned (m,n)
Ci   tol   :numbers differing by less than tol are treated as equal.
Ci   opts  :ones digit
Ci           0 vecs returned sorted.
Ci           1 vecs is unchanged; only iprm is returned.
Ci          tens digit
Ci           0 vecs sorted
Ci           1 vecs sorted by increasing length
Ci         hundreds digit
Ci           1 equal vectors preserve their original order
Co Outputs
Co   iprm  :a permutation table that sorts array 'vecs'
Co   vecs  : may be changed, depending on opts
Cu Updates
Cu   08 Sep 12 Last argument turned into integer, with new iopt=2
Cu   20 Sep 09 vecs no longer need to be doubly dimensioned for opts=0
Cu   02 Sep 02 Added 100s digit switch
C ----------------------------------------------------------------
      implicit none
      integer m,n,iprm(n),opts
      double precision vecs(m,n),tol
      double precision di,dl
      integer l,ir,irra,i,j,mm,i1,i2,nn
      real(8), allocatable :: wk(:,:)
      logical norm

C     Required for f90 compatibility
      interface
      subroutine dvprm(m,n,vecs,wk,iprm,opt)
      implicit none
C ... Passed parameters
      integer m,n,iprm(n),opt
      real(8) vecs(m,n)
      real(8), allocatable, target :: wk(:,:)
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
            di = 0d0
            dl = 0d0
            do  24  mm = 1, m
            dl = dl + vecs(mm,iprm(j))**2
   24       di = di + vecs(mm,iprm(j+1))**2
            dl = dsqrt(dl)
            di = dsqrt(di)
            if (abs(di-dl) > tol) then
              if (di-dl > tol) j = j+1
            endif
          else
            do  26  mm = 1, m
            if (abs(vecs(mm,iprm(j+1))-vecs(mm,iprm(j))) <= tol)
     .          goto 26
            if (vecs(mm,iprm(j+1))-vecs(mm,iprm(j)) > tol) j = j+1
            goto 28
   26       continue
   28       continue
          endif
        endif

C   ... Demote rra to its level
        if (norm) then
          di = 0d0
          dl = 0d0
          do  34  mm = 1, m
          dl = dl + vecs(mm,irra)**2
   34     di = di + vecs(mm,iprm(j))**2
          dl = dsqrt(dl)
          di = dsqrt(di)
          if (di-dl > tol) then
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
            if (abs(vecs(mm,iprm(j))-vecs(mm,irra)) <= tol) goto 36
            if (vecs(mm,iprm(j))-vecs(mm,irra) > tol) then
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
        di = 0d0
        dl = 0d0
        do  124  mm = 1, m
        dl = dl + vecs(mm,iprm(i1))**2
  124   di = di + vecs(mm,iprm(i2))**2
        dl = dsqrt(dl)
        di = dsqrt(di)
        if (abs(di-dl) > tol) goto 130
      else
        do  126  mm = 1, m
  126   if (abs(vecs(mm,iprm(i2))-vecs(mm,iprm(i1))) > tol) goto 130
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
        call dvprm(m,n,vecs,wk,iprm,2)
      endif
      end
      subroutine dvprm(m,n,vecs,wk,iprm,opt)
C- Permute an array of double precision vectors according to iprm
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
Ci          0 assume all iprm are nonzero
Ci          1 cull list, excluding rows for which iprm is zero
Ci            return in n number of elements returned
Ci Inputs/Outputs
Co   n     : if list is to be culled, n is returned as the number
Co         : of elements for which iprm>0.  Otherwise, n is not touched
Co   wk    :returns vecs in permuted order
Co   vecs  :wk may be copied back into vecs, depending on opt.
Cu Updates
Cu   19 Aug 17 new 10s digit opt
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer m,n,iprm(n),opt
      real(8) vecs(m,n)
      real(8), allocatable, target :: wk(:,:)
C ... Local parameters
      integer i,j,k,l
      real(8), pointer :: wkl(:,:)

      if (mod(opt,10) /= 2) then
        wkl => wk
      else
        allocate(wkl(m,n))
        call dcopy(n*m,vecs,1,wkl,1)
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
