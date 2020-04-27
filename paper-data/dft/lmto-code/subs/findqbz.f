      subroutine findqbz(nqbz,qbz,plat,tol,qpi,iq)
C- Find iq in qbz that matches qpi
C ----------------------------------------------------------------------
Ci Inputs
Ci   nqbz  :number of k-points
Ci   qbz   :vectors of k-points
Ci   plat  :primitive lattice vectors, in units of alat
Ci         :(reciprocal transpose of qlat)
Ci   tol   :allowed tolerance in qp mismatch
Ci   qpi   :k-point to find
Co Outputs
Co   iq    :point iq that matches qbz.  No match -> iq = -1
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Jul 12 Updated to shorten vectors
Cu   08 Sep 05 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nqbz,iq
      double precision qbz(3,nqbz),qpi(3),plat(3,3),tol
C ... Local parameters
      integer i
      logical latvec
C     double precision ql(3)

      do  i = 1, nqbz
        if (latvec(1,tol,plat,qpi(:)-qbz(:,i))) then
C        call shorbz(qpi(:)-qbz(:,i),ql,qlat,plat)
C        if (abs(qpi(1)-qbz(1,i)) < tol) then
C        if (abs(qpi(2)-qbz(2,i)) < tol) then
C        if (abs(qpi(3)-qbz(3,i)) < tol) then
          iq = i
          return
C        endif
C        endif
        endif
      enddo
      iq = -1

      end
