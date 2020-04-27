      subroutine rotm(norb,theta,pmat)
C- Rotate complex matrix in spin space around the y axis
C ----------------------------------------------------------------------
Ci Inputs
Ci   norb  :dimension of the matrix
Ci   theta :angle by which to rotate
Ci   pmat  :angle by which to rotate
Co Outputs
Co   pmat  :matrix to rotate
Cu Updates
Cu   16 Apr 12 (Belashchenko) First created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer norb
      double precision theta
      double complex pmat(2*norb,2*norb)
C Local parameters:
      integer n1,n2
      double precision xc,xs
      double complex uu,ud,du,dd

      if (theta == 0d0) return

      xc = dcos(theta/2)
      xs = dsin(theta/2)
      do  n2 = 1, norb
        do  n1 = 1, norb
          uu = pmat(n1,n2)
          ud = pmat(n1,n2+norb)
          du = pmat(n1+norb,n2)
          dd = pmat(n1+norb,n2+norb)
          pmat(n1,n2)           = xc*xc*uu + xs*xs*dd - xc*xs*(ud+du)
          pmat(n1,n2+norb)      = xc*xs*(uu-dd) + xc*xc*ud - xs*xs*du
          pmat(n1+norb,n2)      = xc*xs*(uu-dd) + xc*xc*du - xs*xs*ud
          pmat(n1+norb,n2+norb) = xs*xs*uu + xc*xc*dd + xc*xs*(ud+du)
        enddo
      enddo
      end

