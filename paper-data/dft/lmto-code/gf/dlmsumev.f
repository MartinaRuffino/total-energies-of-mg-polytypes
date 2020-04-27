      subroutine dlmsumev(nclass,nangl,nl,qnu,lgcorr,qcorr,pp,efermi,
     .  etrms)
C- Sum of eigenvalues, DLM
C ----------------------------------------------------------------------
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   nangl :total number of angles for all DLM sites
Ci   nl    :(global maximum l) + 1
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   lgcorr:include correction from dO/dz
Ci   qcorr :charge corrections?
Ci   pp    :potential parameters (atomsr.f)
Ci   efermi:Fermi energy
Co Outputs
Ci   etrms :sumev -> etrms(16,*) for each 'class'
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   02 Jun 12
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lgcorr
      integer nclass,nangl,nl
      double precision qnu(3,nl,2,*),pp(6,nl,2,*),efermi,etrms(22,*)
      double precision qcorr(2,nl,2,*)
C ... Local parameters
      integer ic,isp,l  ,stdo
      double precision q,sumev

C     stdo = nglob('stdo')

      do  ic = nclass + 1, nclass + nangl
        q = 0
        sumev = 0
        do isp = 1, 2
          do  l = 1, nl
            if (lgcorr) then
              sumev = sumev + (qnu(2,l,isp,ic) + qcorr(2,l,isp,ic))
     .                      + (pp(1,l,isp,ic) - efermi)
     .                      * (qnu(1,l,isp,ic) + qcorr(1,l,isp,ic))
            else
              sumev = sumev + qnu(2,l,isp,ic)
     .              + (pp(1,l,isp,ic) - efermi) * qnu(1,l,isp,ic)
            endif
            q = q + qnu(1,l,isp,ic)
          enddo
        enddo
C        write (stdo,500) ic,sumev,q
C  500   format('Class ',i3,': sumev = ',F10.6,', q = ',F10.6)
        etrms(16,ic) = sumev
      enddo

      end
