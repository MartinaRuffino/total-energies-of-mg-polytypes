      subroutine idmxmn(n,a,incx,imin,imax)
C- Finds maximum and minimum absolute value of a double precision array
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of points to evaluate
Ci   a     :array of points
Ci   incx  :spacing between points to consider
Co Outputs
Co   imin  :points to lowest
Co   imax
Cl Local variables
Cl         :
Cr Remarks
Cr   Points a(1),a(1+incx),...,a(1+incx*n) are compared
Cr
Cu Updates
Cu   25 Aug 04 First created from iyamax
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,incx
      double precision a(n)
C ... Local parameters
      integer i,ix,imin,imax
      double precision smin,smax

      imin = 0
      imax = 0
      if (n < 1 .or. incx <= 0) return
      imin = 1
      imax = 1
      if (n == 1) return

      ix = 1
      smin = dabs(a(1))
      smax = dabs(a(1))
      ix = ix + incx
      do  10  i = 2, n
        if (dabs(a(ix)) <= smax) go to 5
        imax = i
        smax = dabs(a(ix))
    5   continue

        if (dabs(a(ix)) >= smin) go to 6
        imin = i
        smin = dabs(a(ix))
    6   continue
        ix = ix + incx
   10 continue
      end
