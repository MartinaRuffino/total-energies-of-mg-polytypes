      subroutine shorigv(q,qs,plat,napw,igapw,igapws)
C- Shift APW G vectors to correspond to change in q
C ----------------------------------------------------------------------
Ci Inputs
Ci   q     :wave number before shift
Ci   qs    :wave number after shift
Ci   plat  :primitive lattice vectors, in units of alat
Ci   napw  :number of augmented PWs in basis
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   igapws:shifted igapw
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Feb 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      double precision q(3),qs(3),plat(3,3)
      integer :: napw,igapw(3,napw),igapws(3,napw)
C ... Local parameters
      integer i,ig,idq(3)
      double precision x0,tol
      parameter (tol=1d-8)

      do  i = 1, 3
C       x0 is projection of q-qs along qlat
        x0 = (q(1)-qs(1))*plat(1,i) +
     .       (q(2)-qs(2))*plat(2,i) +
     .       (q(3)-qs(3))*plat(3,i)
        idq(i) = idnint(x0) ! should be integer multiples of qlat
        if (abs(idq(i)-x0) > tol) call rx('shoigv: bug in shorbz')
      enddo

      do  ig = 1, napw
        igapws(:,ig) = idq + igapw(:,ig)
      enddo

      end
