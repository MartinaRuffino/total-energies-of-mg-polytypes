      double precision function cdf(x)
C- Cumulative normal distribution`
C     phi(x) = (1+erf(x/sqrt(2)))/2 = (2-erfc(x/dsqrt(2d0)))/2
C     dphi/dx = dexp(-x*x/2)/sqrt(2*pi)
      implicit none
C ... Passed parameters
      double precision x,derfc

      cdf = (2-derfc(x/dsqrt(2d0)))/2 ! = (1+erf(x/sqrt(2)))/2
      end

      double precision function cdfix(p)
C- Inverse of the cumulative distribution, with Newton-Raphson improvement
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Cumulative distribution phi
Co Outputs
Co   cdfix :inverse x of phi
Cr Remarks
Cr   The cumulative distribution phi is calculated by
Cr   polynomial approximation; see cdfi below.
Cr   cdfix refines that solution by a Newton-Raphson step.
Cr   The step improves cdfi's precision to:
Cr     1d-14 for -4<cdfix<4, (corresponding to p and 1-p > 1d-5)
Cr     1d-12 for -5<cdfix<5, (corresponding to p and 1-p > 1d-6)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      double precision p
C ... Local parameters
      double precision x,cdf,cdfi,dpdx
C     double precision dp
      real(8),parameter :: pi = 3.1415926535897931d0

      x = cdfi(p)
C     dp = cdf(x)-p ! error in p
      dpdx = dexp(-x*x/2)/sqrt(2*pi) ! d(phi)/dx
      cdfix = x - (cdf(x)-p)/dpdx
      end

      double precision function cdfi(p)
C- Rational approximation to inverse of the cumulative distribution
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Cumulative distribution phi
Co Outputs
Co   cdfi  :inverse x of cumulative distribution phi; see Remarks
Cl Local variables
Cl         :
Cr Remarks
Cr   The cumulative distribution phi is defined as:
Cr     phi(x) = 1/sqrt(2pi) int_-infty^x exp(-t^2/2) dt
Cr            = 1/2[1+erf(x/sqrt(2))], erf(x) = error function
Cr   It has these limiting cases:
Cr     phi(-infty) = 0
Cr     phi(0)      = 1/2
Cr     phi(infty)  = 1
Cr
Cr   The Q-function is defined as Q(x) = 1-phi(x) = phi(-x)
Cr
Cr   This algorithm is found to have a relative precision of 1d-9
Cr   for -6<cdfi<6, corresponding to range where p and 1-p > 1d-9
Cr   It is not reliable for p or 1-p < 1d-9.
Cr   Rational function coefficients were found on the web;
Cr   no further attempt was made to optimize them.
Cr   See cdfix for refinement of solution
Cu Updates
Cu   29 Jan 14  algorithm adapted from a web page
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      double precision p
C ... Local parameters
      double precision x,plow,phigh,r,q
C  Coefficients for rational approximations.
      real(8),parameter :: a(6) =
     .(/-3.969683028665376d+01,
     .   2.209460984245205d+02,
     .  -2.759285104469687d+02,
     .   1.383577518672690d+02,
     .  -3.066479806614716d+01,
     .   2.506628277459239d+00/)
      real(8),parameter :: b(5) =
     .(/-5.447609879822406d+01,
     .   1.615858368580409d+02,
     .  -1.556989798598866d+02,
     .   6.680131188771972d+01,
     .  -1.328068155288572d+01/)
      real(8),parameter :: c(6) =
     .(/-7.784894002430293d-03,
     .  -3.223964580411365d-01,
     .  -2.400758277161838d+00,
     .  -2.549732539343734d+00,
     .   4.374664141464968d+00,
     .   2.938163982698783d+00/)
      real(8),parameter :: d(4) =
     .(/ 7.784695709041462d-03,
     .   3.224671290700398d-01,
     .   2.445134137142996d+00,
     .   3.754408661907416d+00/)

      plow  = 0.02425d0 ! Region I
      phigh = 1 - plow  ! Region III

C  ... Region I
      if (p < plow) then
        q = sqrt(-2*dlog(p))
        x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) /
     .      ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1)

C ... Region II
      elseif (p <= phigh) then
        q = p - 0.5d0
        r = q*q
        x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q /
     .      (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1)

C ... Region III
      else
        q = sqrt(-2*dlog(1-p))
        x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) /
     .       ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1)
      endif
      cdfi = x

C ... Correct using Newton Raphson
C      px = cdf(x)
C      dpdx = dexp(-x*x/2)/sqrt(2*pi)
C      cdfi =

      end

C     Test precision of inverse phi
C      subroutine fmain
C      implicit none
C      double precision x,p,z,cdfi,cdfix,cdf
CC     double precision erf
C
C      print *, cdf(-1.2d0) + cdf(1.2d0) - 1; stop
C
C      do  x = -6, 6, .2d0
CC       p = (1+erf(x/dsqrt(2d0)))/2
C        p = cdf(x)
C        z = cdfi(p)
C        print *, sngl(p),sngl(1-p),sngl(x),sngl(z),sngl(x/z-1)
C        z = cdfix(p)
C        print *, sngl(p),sngl(1-p),sngl(x),sngl(z),sngl(x/z-1)
C
C      enddo
C      end
