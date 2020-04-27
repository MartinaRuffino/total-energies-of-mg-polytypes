      subroutine ocnock(s_optic,iq,evals,nevl,Ef,esmr)
C- Estimate range of (occ,unocc) states for a given set evals
C ----------------------------------------------------------------------
Cio Structures
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  window
Co     Stored:     unrng ocrng
Co     Allocated:  *
Cio    Elts passed:ocrng unrng
Cio    Passed to:  *
Ci Inputs
Ci   iq    :k-point index (printout only)
Ci   evals :Eigenvalues for this kp
Ci   nevl  :number of evals calculated
Ci   Ef    :Fermi energy
Ci   esmr  :smearing -- fixes spread in partial occupancy around Ef
Co Outputs
Co   s_optic%ocrng(3:4) and s_optic%unrng(3:4):
Co         :reduced range of occupied and unoccupied states
Co         :zero values are to be interpreted as "no further restriction"
Cr Remarks
Cr
Cu Updates
Cu   24 Aug 14 Redesigned
Cu   06 Mar 14 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_optic):: s_optic
C ... Passed parameters
      integer iq,nevl
      double precision evals(nevl),Ef,esmr
C ... Local parameters
      integer ie
      double precision omin,omax
      real(8), parameter :: tol = 1d-6
C ... External calls
      external info5

C ... Setup
      s_optic%unrng(3:4) = 0
      s_optic%ocrng(3:4) = 0
      if (s_optic%alltrans /= 0) return  ! No additional constraints

C     omin = s_optic%window(1)
      omin = 0  ! Ef may not be well defined in insulators
      omax = s_optic%window(2)

      do  ie = 1, nevl
C        print *, ie, evals(ie)-ef,omax,
C     .    (evals(ie)+omax < Ef-esmr), s_optic%ocrng(3)
C        print *, ie, evals(ie)-ef,omin,
C     .    (evals(ie)-omin < Ef-esmr),s_optic%unrng(3)

C       If E+omax<Ef exclude from occ lower bound
        if (evals(ie)+omax < Ef-esmr) s_optic%ocrng(3) = ie+1
C       If E-omin<Ef exclude from unocc lower bound
        if (evals(ie)-omin < Ef-esmr) s_optic%unrng(3) = ie+1

      enddo

      do  ie = nevl, 1, -1

C        print *, ie, evals(ie)-ef,omin,
C     .    evals(ie)+omin > Ef+esmr,s_optic%ocrng(4)
C        print *, ie, evals(ie)-ef,omax,
C     .    (evals(ie)-omax > Ef+esmr),s_optic%unrng(4)

C       If E+omin>Ef exclude from occ upper bound
        if (evals(ie)+omin > Ef+esmr) s_optic%ocrng(4) = ie-1
C       If E-omax>Ef exclude from unocc upper bound
        if (evals(ie)-omax > Ef+esmr) s_optic%unrng(4) = ie-1

      enddo

C      if (lmet == 0) then  ! Possible problem for noninteger qval?
C        s_optic%ocrng(4) = min(int(qval+.999999d0),s_optic%ocrng(4))
C        s_optic%unrng(3) = max(int(qval+1),s_optic%unrng(3))
C      endif

C     print *, s_optic%ocrng(3:4),'  ',s_optic%unrng(3:4)
      s_optic%ocrng(3) = max(s_optic%ocrng(3),s_optic%ocrng(1))
      s_optic%ocrng(4) = min(s_optic%ocrng(4),s_optic%ocrng(2))
      s_optic%unrng(3) = max(s_optic%unrng(3),s_optic%unrng(1))
      s_optic%unrng(4) = min(s_optic%unrng(4),s_optic%unrng(2))
C     print *, s_optic%ocrng(3:4),'  ',s_optic%unrng(3:4)

      call info5(35,0,0,' optics, iq=%i: restrict transitions to : '
     .  //'occ=(%i,%i) unocc=(%i,%i)',iq,s_optic%ocrng(3),s_optic%ocrng(4),
     .  s_optic%unrng(3),s_optic%unrng(4))

      end

      subroutine opt_nstate(s_optic,mode,nev,nfilo,nfiup,nemlo,nemup,nfilm,nempm)
C- Returns range of occ,unocc states for a given k needed for optics
C ----------------------------------------------------------------------
Cio Structures
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ocrng unrng (see routine ocnock)
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 return nfilo,nfiup,nemlo,nemup
Ci         :2 return nfilm,nempm
Ci         :1 and 2 may be taken in combination
Ci   nev   :maximum available number of eigensates
Co Outputs
Co   nfilo,nfiup :Matrix elements for occupied bands in range nfilo:nfiup
Co   nemlo,nemup :Matrix elements for unoccupied bands in range nemlo:nemup
Co   nfilm :number of filled states
Co   nempm :number of empty states
Cr Remarks
Cr
Cu Updates
Cu   21 Aug 14 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_optic):: s_optic
C ... Passed parameters
      integer nev,mode,nfilo,nfiup,nemlo,nemup,nfilm,nempm
C ... Local parameters
      integer nfilox,nfiupx,nemlox,nemupx  ! Work with local variables in case

      nfilox = s_optic%ocrng(1); nfiupx = min(s_optic%ocrng(2),nev)
      nemlox = s_optic%unrng(1); nemupx = min(s_optic%unrng(2),nev)
      if (s_optic%ocrng(3) > 0) nfilox = max(nfilox,s_optic%ocrng(3))
      if (s_optic%unrng(3) > 0) nemlox = max(nemlox,s_optic%unrng(3))
      if (s_optic%ocrng(4) > 0) nfiupx = min(nemupx,s_optic%ocrng(4))
      if (s_optic%unrng(4) > 0) nemupx = min(nemupx,s_optic%unrng(4))
      nfilox = min(nfilox,nfiupx+1)
      nemupx = max(nemupx,nemlox-1)
      if (mod(mode,2) /= 0) then
        nfilo = nfilox; nfiup = nfiupx; nemlo = nemlox; nemup = nemupx
      endif
      if (mod(mode/2,2) /= 0) then
        nfilm = nfiupx-nfilox+1
        nempm = nemupx-nemlox+1
      endif

      end
