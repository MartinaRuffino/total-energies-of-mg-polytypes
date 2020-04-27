      subroutine gibbswts(nthet,cs,beff,wt)
C- Calculate and normalize DLM weights from effective field Legendre coefs
C ----------------------------------------------------------------------
Ci Inputs
Ci   nthet  :number of angles
Ci   cs     :cosines of theta angles
Ci   beff   :Ledendre coefficients for the effective field
Cio  wt     :weights (Gaussian on input, Gibbs on output)
Cr Remarks
Cb Bugs
Cu Updates
Cu   18 Apr 12 (Belashchenko) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nthet
      double precision cs(nthet),beff(10),wt(nthet)
C ... Local parameters
      integer ith,i
      double precision bf,wk(nthet),plegn,norm

c     print *,'GIBBSWTS:',(wt(i),i=1,nthet)
      norm = 0
      do  ith = 1, nthet
        bf = 0
        do  i = 1, 10
          bf = bf + beff(i)*plegn(i,cs(ith))
        enddo
        wk(ith) = dexp(-bf)
        norm = norm + wk(ith)*wt(ith)
      enddo
c     print *,'WK:',(wk(i),i=1,nthet)

      wt(:) = wk(:)*wt(:)/norm


      end
