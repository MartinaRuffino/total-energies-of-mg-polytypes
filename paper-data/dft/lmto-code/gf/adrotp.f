      subroutine adrotp(mode,norb,p1,p2,theta,pmat)
C- Rotate matrix with pfun on the diagonal and record in pmat
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 add to pmat
Ci         :1 zero out pmat before adding P
Ci   norb  :dimension of the matrix
Ci   p1    :vector of spin-up diagonal elements of non-rotated matrix
Ci   p2    :vector of spin-dn diagonal elements of non-rotated matrix
Ci   theta :angle by which to rotate
Co Outputs
Co   pmat  :matrix to record output in
Cu Updates
Cu   09 Apr 12 (Belashchenko) First created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer mode,norb
      double precision theta
      double complex p1(norb),p2(norb),pmat(2*norb,2*norb)
C Local parameters:
      integer n,ns
      double precision xc,xs
      double complex ps,pa

      if (mode == 1) call dpzero(pmat,norb*norb*2*2*2)

      if (theta /= 0d0) then
        xc = dcos(theta)
        xs = dsin(theta)
        do  n = 1, norb
          ps = (p1(n) + p2(n))/2
          pa = (p1(n) - p2(n))/2
          ns = n + norb
          pmat(n,n)   = pmat(n,n)   + ps + pa*xc
          pmat(n,ns)  = pmat(n,ns)       + pa*xs
          pmat(ns,n)  = pmat(ns,n)       + pa*xs
          pmat(ns,ns) = pmat(ns,ns) + ps - pa*xc
        enddo
      else
        do  n = 1, norb
          ns = n + norb
          pmat(n,n)   = pmat(n,n)   + p1(n)
          pmat(ns,ns) = pmat(ns,ns) + p2(n)
        enddo
      endif
      end

