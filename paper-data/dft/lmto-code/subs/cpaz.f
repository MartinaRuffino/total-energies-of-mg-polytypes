      subroutine cpaz(nclass,z,ncomp,wts)
C- Nuclear charge for CPA classes overwritten by average of chemical constituents
C ----------------------------------------------------------------------
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   ncomp :ncomp(ic) total number of components for class ic
Ci   wts   :component weights
Cio Inputs/Outputs
Cio  z     :nuclear charge, by class.  On ouput CPA classes are modified
Cio        :NB this can be any class-dependent scalar, e.g. core charge
Cl Local variables
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nclass,ncomp(nclass)
      double precision z(*),wts(*)
C ... Local parameters
      integer ic,icomp,ioff
      double precision z0

      ioff = nclass
      do  ic = 1, nclass
        if (ncomp(ic) < 2) cycle
        z0 = 0d0
        do  icomp = 1, ncomp(ic)
          z0 = z0 + z(ioff+icomp)*wts(ioff+icomp)
        enddo
        ioff = ioff + ncomp(ic)
        z(ic) = z0
      enddo

      end

