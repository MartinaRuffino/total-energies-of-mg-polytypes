      subroutine linemin(xn,fn,slope,xtol,dxmn,dxmx,alf,wk,ir)
C- Estimate a point along a given line minimize a function of many variables
C ----------------------------------------------------------------------
Ci Inputs
Ci   xn    :current distance along the line relative to the starting point.
Ci         :At the start of a new line minimization, xn should be zero.
Ci         :In subsequent calls, xn is updated by linemin.
Ci   fn    :function value at current xn (see Remarks).
Ci         :At the start of th eline minimization, xn is assumed to be zero
Ci   slope :gradient of fn at xn=0
Ci   xtol  :If zero, xtol is not used.
Ci         :If nonzero, xtol is used as a second criterion for ending the current
Ci         :line minimization; see Remarks
Ci   dxmn  :it should be no smaller than the smallest xn differing from zero for which
Ci         :the change in fn is meaningful.
Ci         :The mulivariate algorithm calling linemin can sometimes converge faster
Ci         :dxmn is somewhat larger.
Ci   dxmx  :If dxmx>0, then dxmx caps the maximum change linemin will make in xn
Ci   alf   :a parameter used in this criterion:
Ci         :"Sufficient" reduction when fn(xn) < fn(0) + alf*xn*slope
Ci         :This criterion is used to determine whether the current line minimizaton should end.
Cio Inputs/Outputs
Cio ir:  Program flow control (see Remarks)
Cio      Input:  To start a new line minimization, set ir to zero.
Cio      In each iteration linemin will set ir, and that value should
Cio      be preserved while iterating within a line minimization.
Cio      On exit, linemin sets ir to direct subsequent program flow.
Cio        ir   meaning
Cio         0:  f=0.  linemin does nothing.
Cio        -1:  continue on the existing line with a new xn linemin supplies.
Cio             Evaluate f at the new xn and call linemin again with all
Cio             other parameters unchanged.
Cio         1:  fn has decreased sufficiently end this line minimization (see Remarks)
Cio             The calling program should start a new line minimization, or
Cio             quit if the result is sufficiently converged.
Cio         2:  xn falls below dxmn. xn is returned unchanged.
Cio             This condition may be normal; it may be a spurious local minimum when
Cio             linemin is used to find the zeros of a vector function.  See gradzr.f
Cio             A calling program should decide on one of the following:
Cio               (1) convergence has been reached
Cio               (2) a spurious local minimum was found
Cio                  (occurs when finding zeros of a function and |slope|/fn is very small)
Cio               (3) a new line minimization should be started
Cio               (4) there is a some other numerical problem with the function.
Cio         The following values indicate abnormal conditions. linemin does nothing.
Cio         6: nonsensical tolerance
Cio         7: input xn equals x1 or x2 (some value of xn from a prior call).
Cio         8: If the initial grad fn . h is nonnegative.  Implies h does not point downhill.
Cio  wk:  : work array of length 12.  Should not be altered by calling program.
Cio         wk holds:
Cio         1..3  xn from prior calls
Cio         4..6  fn from prior calls
Cio         7     ir from the prior call => indicates whether >=3 points on current line
Cio         8     xn corresponding to fn in wk(9)
Cio         9     smallest fn encountered in line min
Cio         10    fn at starting point xn=0, i.e. when ir=0
Cio               Note slope=-2f when finding zeros of a vector function F and fn = |F|^2/2
Cio         11    internal copy of slope, in case calling program destroys it in the
Cio               course of a line minimization.
Co Outputs
Co   xn:  next estimate for minimum point along current line min.
Co        Initial xn (ir=0) is assumed to be zero.
Cr Remarks:
Cr  *linemin is a kernel for Newtonian methods that minimize some function fn,
Cr   given their gradient (or when finding the zeros of a vector function F, the
Cr   Jacobian). It differs from rfalsi, another line minimization scheme that works
Cr   with less information (no gradient).  See gradzr.f.
Cr
Cr  *linemin can be used in one of two contexts:
Cr     (1) to minimize a function fn given gradients df/dxi and the Jacobian d^2f/dxi dxj
Cr     (2) to find zeros of an n-component vector function F of n variables.
Cr         In this case define scalar fn = |F|^2/2.  Define the gradient and Jacobian to be
Cr   Let F be a multivalued function of n variables.  Minimize scalar fn = |F|^2/2.
Cr   The Jacobian is J(i,j) = dF(i) / dx(j).  Define gradient g(i) = sum F(j) J(j,i).
Cr
Cr  *linemin works by making excursions along line x(i) = x0(i) + xn*h(i), where x0 and h(i)
Cr   are starting point and displacement.  x0 and h are unspecified here, but :
Cr    - for minimization methods (1), h is normally the gradient direction at x=x0.
Cr    - for zero-finding methods (2), h is normally the "conjugate" direction, h = J^-1 F.
Cr   linemin returns some new estimate xn given  fn from a prior xn (at the start of a
Cr   line minimization, xn is assumed to be zero). The calling program should generate fn
Cr   new at points x0+h(xn) (linemin returns xn), until linemin returns a non-negative ir.
Cr   At that point, the calling program must decide whether :
Cr    (a) the system is globally converged
Cr    (b) start a new line minimization
Cr    (c) manage abnormal exits.
Cr  *linemin's first estimate for xn is xn=1.
Cr
Cr  *convergence : linemin returns ir=1 when the minmization is "sufficiently" converged.
Cr   The minimum criterion for "sufficient" is that fn(xn) < fn(0) + alf*xn*slope.
Cr   If construction of the Jacobian is not too expensive, there is little advantage
Cr   to continue minimizing along a fixed direction.  So any small reduction in fn
Cr   is sufficient (this is the strategy Numerical Recipes' lnsrch uses).
Cr   However, it may be that generating the gradient (or Jacobian) is relatively expensive.
Cr   If calculating the gradint (Jacobian) is relatively expensive, it may be advantageous
Cr   to impose further conditions on the line minimization.
Cr   If you sent xtol to a postive number, linemin imposes an extra criterion before
Cr   terminating a line minimization: he change in xn relative to its prior value must
Cr   be less than xtol.  Thus, if 0<xtol<1, linemin will always look for at least two
Cr   points along the line apart from the starting point (xn=0) since it always chooses
Cr   for its first iteration xn=1.
Cr
Cr  Based on surbroutine lnsrch in Numerical Recipes.
Cu Updates
Cu   01 Oct 14 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ir
      double precision xn,slope,xtol,fn,dxmn,dxmx,alf,wk(0:11)
C ... Local parameters
      integer i,ipr,ir0,stdo,i1mach
      double precision f0,f1,f2,fx0,x0,x1,x2,rhs1,rhs2,a,b,xtry,disc,fpp,de1,del1
      double precision xnl(0:2),fnl(0:2),fp(2)
      equivalence (x0,xnl(0)),(x1,xnl(1)),(x2,xnl(2))
      equivalence (f0,fnl(0)),(f1,fnl(1)),(f2,fnl(2))

      stdo = i1mach(2)

C ... Point-independent setup
      call getpr(ipr)
C     ipr = 55
      if (fn == 0) then  ! OK for zero finding.  But for minimization !?
        ir = 0; return
      endif

C --- Starting iteration ---
      if (ir == 0) then
        call dpzero(wk,12)
        wk(10) = fn    ! Starting function value
        wk(11) = slope ! Starting slope
        call info5(-40,0,0,' linemin: new line: slope=%1;3g  f=%1;3g  xtol=%1;3g  dxmn=%1;3g   alf=%1;3g',
     .    slope,fn,xtol,dxmn,alf)
        ir = 8; if (slope >= 0d0) goto 998
        ir = -1
        xn = 1
        if (dxmx > 0) xn = min(dxmx,1d0)
        fnl = fn; xnl = 0
        wk(8) = 0; wk(9) = fn
        goto 998
      endif

C ... Recover local variables from work array
      call dcopy(3,wk(1),1,xnl,1)
      call dcopy(3,wk(4),1,fnl,1)
      fx0 = wk(10); slope = wk(11) ! Function and slope at xn=0
      ir0 = wk(7) ! ir from the prior iteration
      wk(7) = ir  ! ir of the prior iteration next time linemin is called
      if (fn < wk(9)) then ! Update the current point of minimum f
        wk(8) = xn; wk(9) = fn
      endif

C ... Push prior points onto stack: (x1,f1) -> (x2,f2); (x0,f0) -> (x1,f1); (xn,fn) -> (x0,f0)
      do  i = 1, 0, -1
        fnl(i+1) = fnl(i)
        xnl(i+1) = xnl(i)
      enddo
      fnl(0) = fn
      xnl(0) = xn

C ... Sanity checks
      ir = 8; if (slope >= 0) goto 998
      ir = 7; if (x1 == xn .or. x2 == xn) goto 998

C ... Possible spurious root
      if (xn < dxmn) then ! function minimum apparently within dxmn of origin
        ir = 2
        call info5(-30,0,0,' linemin ir=%,2i: current xn=%1,8;8g < dxmn=%1;4g',ir,xn,dxmn,0,0)
        goto 999
      endif

      ir = -1 ! Default return

C ... Check for sufficient function decrease
      if (fn <= fx0 + ALF*xn*slope) then ! This block executes if fn drops enough
        ir = 1

        if (xtol == 0 .or. abs(xn-x1) <= abs(xtol)) goto 998  ! Also meets xtol criterion
C       Some improvement in fn, but last change in xn exceeded tolerance
        if (fn > wk(9)) then ! If some prior point is better than this one, terminate line min there
          xn = wk(8); fn = wk(9); goto 998
        endif
C       Make a new estimate for xn on this line and check whether it is significantly different
        fpp = 2*(fn-fx0-xn*slope)/xn**2  ! From slope at origin, fn at current and previous point
        if (abs(xn*fpp/slope) < 0.1d0) goto 998 ! Possibly ill conditioned case; skip improvements
        xtry = -slope/fpp ! Estimate for all points, including first
        if (ir0 /= 0) then ! Three points available on this line  !!? test should be ir0 /= 0
          de1 = (x1 - x2)/2; del1 = (x1 + x2)/2 - x0; fp(1) = f1; fp(2) = f2
          call dfphi(de1,del1,de1,del1,1,f0,fp,.false.)
          if (abs(fp(1)/fp(2)) < abs(-slope/fpp-x0)) then ! Choose the more conservative change
            xtry = x0 - fp(1)/fp(2)
          endif
        endif
        if (abs(xn-xtry) < dxmn) goto 998 ! New xn too similar to warrant continuing this line
        ir = -1 ! Continue line min because xtol condition not met
C        if (dxmx > 0) then  !? move outside if-then block ?
C          xtry = xn + dsign(min(dabs(xtry-xn),dxmx),xtry-xn) ! Do not allow change in xn to exceed dxmx
C        endif
      elseif (ir0 == 0) then   ! First point after start of line min
        xtry = -slope/(2d0*(fn-fx0-slope))
      else                    ! Subsequent points
        rhs1 = fn-fx0-xn*slope
        rhs2 = f1-fx0-x1*slope
        a = (rhs1/xn**2-rhs2/x1**2)/(xn-x1)
        b = (-x1*rhs1/xn**2+xn*rhs2/x1**2)/(xn-x1)
        if (a == 0d0) then
          xtry = -slope/(2d0*b)
        else
          disc = b*b - 3d0*a*slope
          if (disc < 0d0) then
            xtry = xn/2
          else if (b <= 0d0) then
            xtry = (-b+dsqrt(disc))/(3d0*a)
          else
            xtry = -slope/(b+dsqrt(disc))
          end if
        end if
        if (xtry > xn/2) xtry = xn/2 ! lambda < 0.5 lambda_1
      endif
      if (dxmx > 0) then
        if (xtry /= xn + dsign(min(dabs(xtry-xn),dxmx),xtry-xn)) then
          xtry = xn + dsign(min(dabs(xtry-xn),dxmx),xtry-xn) ! Do not allow change in xn to exceed dxmx
        endif
      endif
      xn = max(xtry,0.1d0*xn) ! lambda >= 0.1 lambda_1

C --- Cleanup: printout current status and new point ---
  998 continue
      if (ir >= 5) then   ! error
        call info2(30,0,0,' linemin: ir=%,2i (error exit) current xn=%1,8;8g',ir,xn)
      else
        call info8(30,0,0,
     .    ' linemin: ir=%,2i x1=%1;4g f1=%1;4g%?#n>0# x2=%1;4g f2=%1;4g#%2j#: %?#n<0#seek#final# xn=%1,8;8g',
     .    ir,x0,f0,ir,x1,f1,ir,xn)
      endif

C ... Exit: repack local variables into wk
  999 continue
      call dcopy(3,xnl,1,wk(1),1)
      call dcopy(3,fnl,1,wk(4),1)

      end subroutine linemin

