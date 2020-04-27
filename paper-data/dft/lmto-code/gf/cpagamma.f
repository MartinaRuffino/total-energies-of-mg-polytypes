      subroutine cpagamma(ncomp,ic0,cpawt,nl,nsp,pp,gam)
C- Returns ASA potential parameter gamma averaged over CPA components
C ----------------------------------------------------------------------
Ci Inputs
Ci   ncomp :number of components for this class
Ci   ic0   :starting CPA class index
Ci   cpawt :weighting of CPA class
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pp    :potential parameters (atomsr.f)
Co Outputs
Co   gam   :potential parameter.  If CPA, CPA averaged weighting
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 Jan 16
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ncomp,ic0,nsp,nl
      double precision cpawt(ic0+ncomp),pp(6,nl,nsp,*),gam(nl,nsp)
C ... Local parameters
      integer ic,icomp,isp,il
      double precision swt

      if (ncomp == 1) then ! No CPA components to average
        gam(:,:) = pp(5,:,:,ic0)
        return
      endif

C     No averaging over CPA components ... just take the first one
C      gam(:,:) = pp(5,:,:,ic0)
C      return

      gam = 0
      do  isp = 1, nsp
        do  il = 1, nl
          swt = 0
          do  icomp = 1, ncomp  ! Loop over CPA components
            ic = ic0 + icomp - 1 ! Class for this CPA component
            gam(il,isp) = gam(il,isp) + pp(5,il,isp,ic)*cpawt(ic)
            swt = swt + cpawt(ic)
          enddo
          if (abs(swt-1) > 1d-6) call rx('cpagamma: bad CPA weighting')
        enddo
      enddo

      end
