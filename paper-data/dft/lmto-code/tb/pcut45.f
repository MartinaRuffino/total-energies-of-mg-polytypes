      subroutine pcut45(n,x,r1,r2,val,slo,curv,e,de,dde)
C- Evaluate cutoff polynomial P_n and its first derivative P'_n at x
C ----------------------------------------------------------------------
Ci Inputs
Ci   n      :order of the cutoff polynomial P_n(x) (n=4 or n=5, see Remarks)
Ci   x      :P_n, P'_n and P"_n are evaluated at x
Ci   r1, r2 :cutoff interval within which function f(x) to cut off is replaced
Ci           (or multiplyed) by P_n(x)
Ci   val    : f(r1)
Ci   slo    : f'(r1)
Ci   curv   : f"(r1)
Co Outputs
Co   e      : P_n(x)
Co   de     : P'_n(x)
Co   dde    : P"_n(x)
Cr Remarks
Cr   1. The purpose of cutoff polynomial P_n(x) of order n is to smoothly cut
Cr      off given function f(x) (e.g., pair potential, hopping integrals)
Cr      to zero between x = r1 and x = r2
Cr
Cr   2. Polynomials of odd and even order.
Cr      Odd n = 2k+1 implies that P_{2k+1} matches continuously up to k-th
Cr      derivatives f(x) at r1 and 0 at x = r2
Cr      Even n = 2k implies that P_{2k} matches continuously up to k-th
Cr      derivatives V0(x) at r1 but only (k-1) derivatives = 0 at x = r2
Cr
Cr   3. Augmentative and multiplicative cutoffs
Cr      f(x) can either be augmented (replaced) by or multiplied by P_n(x)
Cr      within the interval r1 < x < r2.
Cr      augmentative cutoff:
Cr                         f(x)  <=>     Pn(x)   <=>  0
Cr                               x=r1            x=r2
Cr      multiplicative cutoff:
Cr                         f(x)  <=>  Pn(x)*f(x) <=>  0
Cr                               x=r1            x=r2
Cr      However, the multiplicative cutoff polynomial is equivalent to the
Cr      augmentative cutoff polynomial for the particular case f(x) = 1.
Cr      So, we do not consider it as a seperate case.
Cr
Cr   4. Currently only n=5 and n=4 are implemented, ie P_n(x) matches
Cr      value, slope, and curvature of f(x) at r1, whereas
Cr      P_n(r2) = P'_n(r2) = 0 (n = 4)
Cr      P_n(r2) = P'_n(r2) = P"_n(r2) = 0 (n = 5)
Cr
Cr   5. Explicit formulae.
Cr      n = 5:
Cr      P_5(x) = P_2(x) * (x-r2)^3
Cr      where P2(x) = a*(x-r1)^2 + b*(x-r1) + c
Cr         a = [0.5*f"(r1)*(r1-r2)^2 -3*f'(r1)*(r1-r2) + 6*f(r1)]/ (r1-r2)^5
Cr         b = [f'(r1)*(r1-r2) - 3*f(r1)]/ (r1-r2)^4
Cr         c = f(r1) / (r1-r2)^3
Cr      n = 4:
Cr      P_4(x) = P_2(x) * (x-r2)^2
Cr      where P2(x) = a*(x-r1)^2 + b*(x-r1) + c
Cr         a = [0.5*f"(r1)*(r1-r2)^2 -2*f'(r1)*(r1-r2) + 3*f(r1)]/ (r1-r2)^4
Cr         b = [f'(r1)*(r1-r2) - 2*f(r1)]/ (r1-r2)^3
Cr         c = f(r1) / (r1-r2)^2
Cu Updates
Cu   19 Apr 11 (SL) first created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer, intent(in) ::  n
      double precision, intent(in)  :: x,r1,r2,val,slo,curv
      double precision, intent(out) :: e,de,dde
C Local paramters
      double precision rr,xr1,xr2
      double precision a,b,c,pnorm,p2,dp2,ddp2

      rr = r1 - r2
      xr1 = x - r1
      xr2 = x - r2

      c = val*rr*rr
      if (n == 5) then
        pnorm = rr**(-5)
        a = (0.5d0*curv*rr - 3d0*slo)*rr + 6d0*val
        b = (slo*rr - 3d0*val)*rr
      elseif (n == 4) then
        pnorm = rr**(-4)
        a = (0.5d0*curv*rr - 2d0*slo)*rr + 3d0*val
        b = (slo*rr - 2d0*val)*rr
      else
        call rxi(' pcut45: order of polynomial should be either 4 or 5')
      endif

      p2 = pnorm*(c + xr1*(b + xr1*a))
      dp2 = pnorm*(b + xr1*2d0*a)
      ddp2 = pnorm*2d0*a

C ... P_n(x) = P_2(x-r1) * (x-r2)^(n-2), n = 4,5   ==>
C     P'_n(x) = [P'_2(x-r1) * (x-r2) + (n-2) * P_2(x-r1)] * (x-r2)^(n-3)
Cc    P"_n(x) = (n-3)P'_n(x)/(x-r2) + (x-r2)^(n-3)[(n-1)P'_2(x-r1) + (x-r2)P"_2(x-r1)]
C     P"_n(x) = [P"_2(x-r1) * (x-r2)^2 + 2*(n-2)*P'_2(x-r1) * (x-r2) + (n-2)*(n-3)*P_2(x-r1)]
c               * (x-r2)^(n-4)
      e = p2 * xr2**(n-2)
      de = (xr2*dp2 + float(n-2)*p2) * xr2**(n-3)
c     dde = float(n-3)*de/xr2 + xr2**(n-3)*(float(n-1)*dp2+xr2*ddp2)
c If n = 4 and  x = r2 then (x-r2)**(n-4) can be udefined. Deal with this case explicitly
      dde = (xr2*xr2*ddp2+float(2*(n-2))*xr2*dp2+float((n-2)*(n-3))*p2)
      if (n /= 4) dde = dde * xr2**(n-4)

      end
