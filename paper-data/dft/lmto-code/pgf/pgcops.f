      subroutine pgcops(l0,linl,linr,loutl,loutr,opt,sin,sout)
C- Copies iface strux into bulk strux
C ----------------------------------------------------------------
Ci Inputs
Ci   sin(l0,l0+linl+linr): copy from sin one or more of three blocks,
Ci        depending on opt.
Ci   opt  1s digit
Ci        0  leave sout_00 undisturbed
Ci        1  copy sin_00 to sout_00
Ci        10s digit
Ci        0  leave sout_0L undisturbed
Ci        1  copy sin_0L to sout_0L
Ci        2  copy sin_0R+ to sout_0L
Ci        100s digit
Ci        0  leave sout_0R undisturbed
Ci        1  copy sin_0R to sout_0R
Ci        2  copy sin_0L+ to sout_0R
Ci        1000s digit
Ci        0  s is real
Ci        1  s is complex
Co Outputs
Co   sout(l0+linr+loutr): copy into sout blocks from sin
C ----------------------------------------------------------------
      implicit none
      integer l0,linl,linr,loutl,loutr,opt
      double precision sin(l0,l0+linl+linr,2),sout(l0,l0+loutl+loutr,2)
      integer i,j,kin,kout,ic
      double precision xx

      ic = 1
      xx = 1
C ... Re-entry for imaginary part ...
    2 continue

C --- Copy into sout_00 ---
      if (mod(opt,10) == 1) then
        do  10  j = 1, l0
        do  10  i = 1, l0
   10   sout(i,j,ic) = sin(i,j,ic)
      endif

C --- Copy into sout_0L ---
      if (mod(opt/10,10) == 1) then
        if (linl /= loutl) goto 999
        do  110  j = 1, loutl
        do  110  i = 1, l0
  110   sout(i,j+l0,ic) = sin(i,j+l0,ic)
      endif

      if (mod(opt/10,10) == 2) then
        if (linr /= loutl .or. l0 /= loutl) goto 999
        kin  = l0+linl
        kout = l0
        do  120  j = 1, loutl
        do  120  i = 1, l0
  120   sout(i,j+kout,ic) = xx*sin(j,i+kin,ic)
      endif

C --- Copy into sout_0R ---
      if (mod(opt/100,10) == 1) then
        if (linr /= loutr) goto 999
        kin  = l0+linl
        kout = l0+loutl
        do  210  j = 1, loutr
        do  210  i = 1, l0
  210   sout(i,j+kout,ic) = sin(i,j+kin,ic)
      endif

      if (mod(opt/100,10) == 2) then
        if (linl /= loutr .or. l0 /= loutr) goto 999
        kin  = l0
        kout = l0+loutl
        do  220  j = 1, loutr
        do  220  i = 1, l0
  220   sout(i,j+kout,ic) = xx*sin(j,i+kin,ic)
      endif

C --- Another pass for the imaginary part ---
      if (ic == 1 .and. mod(opt/1000,10) == 1) then
        ic = 2
        xx = -1
        goto 2
      endif

      return

  999 call rx('pgcops: bad input')
      end
C      subroutine fmain
C      implicit none
C      double precision sin(1000),sout(1000)
C      integer l0,linl,linr,loutl,loutr,opt
C
C      l0 = 3
C      linl = 3
C      linr = 6
C      loutl = 2
C      loutr = 3
C      opt = 1200
C
C      call init(sin,l0,l0,(l0+linl+linr)*2)
C      call dpzero(sout,1000)
C      call yprm('sin',2,sin,l0*(l0+linl+linr),l0,l0,l0+linl+linr)
C      call pgcops(l0,linl,linr,loutl,loutr,opt,sin,sout)
C      call yprm('sout',2,sout,l0*(l0+loutl+loutr),l0,l0,l0+loutl+loutr)
C
C      end
C
C      SUBROUTINE init ( c, ldc, n, m )
CC- Initialize arrays
C      integer ldc, n, l
C      double precision c( ldc, m )
C      do  10  i = 1, n
C      do  10  j = 1, m
C   10 c(i,j) = dble( i )/j
C
C      end
