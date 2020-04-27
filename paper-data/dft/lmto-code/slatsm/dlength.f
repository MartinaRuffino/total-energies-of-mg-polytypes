      double precision function dlength(n,dx,incx)
C- Compute the length of a vector of dimension n
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :dimension of vector
Ci   dx    :dx vector
Ci   incx  :spacing between elements
Co Outputs
Co  dlength:sqrt(dx . dx)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   08 Feb 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,incx
      double precision dx(n)
C ... Local parameters
      double precision ddot

      dlength = dsqrt(ddot(n,dx,incx,dx,incx))

      end
