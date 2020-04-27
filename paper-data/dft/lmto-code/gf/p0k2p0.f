      subroutine p0k2p0(aintra,nRlc,lpdim,nsp,iq,nkp,P0k,P01)
C- Transfers the response matrix for 1 k-point; converts to cplx0 format
C ----------------------------------------------------------------------
Ci Inputs
Ci   nRlc  :dimension of P0k: contracted over m
Ci   aintra:F if to contract P0k over l and m
Ci         :T if to contract P0k over m
Ci   lpdim :dimension of P0
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of k-points for which P0 is available
Ci   iq    :current k-point
Ci   P0    :response matrix for all k-points, aka dgdk
Co Outputs
Co   P01   :response function for 1 k-point, contracted over spin,
Co         :possibly contracted over l.
Co         :Also, Re, Im part split.
Cl Local variables
Cr Remarks
Cu   07 Nov 07   Rewritten
Cr   Oct 15 2004 (T. Sandu) Remade
Cr   Feb 20 2004 (T. Sandu) First written
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical aintra
      integer nRlc,nsp,nkp,iq,lpdim
      double precision P0k(2,nRlc,nRlc,nkp,2)
C ... Local parameters
      integer ii,ij
      double precision P01(nRlc,nRlc,2)

      if (lpdim /= nRlc) call rx('for now, lpdim must be nRlc')

      do  ii = 1, nRlc
      do  ij = 1, nRlc
        P01(ii,ij,1) = (P0k(1,ii,ij,iq,1)+P0k(1,ii,ij,iq,nsp))/(3-nsp)
        P01(ii,ij,2) = (P0k(2,ii,ij,iq,1)+P0k(2,ii,ij,iq,nsp))/(3-nsp)
      enddo
      enddo
      end

C      subroutine p02p0k0(nRLc,nkp,iq,WSk,WS1)
CC- Transfers static screened potential for 1 k-point
CC ----------------------------------------------------------------------
CCi Inputs
CCi   nRLc  :dimension of P0
CCi   nkp   :number of k-points for which P0 is available
CCi   iq    :current k-point
CCi   WS1   :static screened potential for 1 k-point, one spin
CCo Outputs
CCo     WSk static screened potential for all k-points
CCl Local variables
CCr Remarks
CCr Made Feb 26 2004 (T. Sandu)
CC ----------------------------------------------------------------------
CC     implicit none
CC ... Passed parameters
C      integer nRLc,nkp,iq
C      double precision WSk(nRLc,nRLc,2,nkp)
CC ... Local parameters
C      integer ii,ij
C      double precision WS1(nRLc,nRLc,2)
C
C
C      do  100  ii = 1, nRLc
C      do  100  ij = 1, nRLc
C      WSk(ii,ij,1,iq) = WS1(ii,ij,1)
C      WSk(ii,ij,2,iq) = WS1(ii,ij,2)
C  100 continue
Cc      call yprm('WS1',2,WS1,nRLc*nRLc,nRLc,nRLc,nRLc)
Cc      call yprm('WSk',2,WSk(1,1,1,iq),nRLc*nRLc,nRLc,nRLc,nRlc)
C      end

      subroutine gfgwk2gw(nlma,lidim,nkp,iq,P0k,P01)
C- Transfers sigma for 1 k-point
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma,lidim  :dimensions of sigma
Ci   nkp   :number of k-points for which sigma is available
Ci   iq    :current k-point
Ci   P0    :sigma matrix for all k-points
Co Outputs
Co    P01 sigma for 1 k-point
Cl Local variables
Cr Remarks
Cr Made May 21 2004
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlma,lidim,nkp,iq
      double precision P0k(2,nlma,lidim,nkp)
C     character*(20) fmt
C ... Local parameters
      integer ii,ij
      double precision P01(nlma,lidim,2)
C     fmt = '(9f15.7)'

c      print *, 'gfgwkgw merge with p0k2p0'

      do  100  ij = 1, lidim
      do  100  ii = 1, nlma
      P01(ii,ij,1) = P0k(1,ii,ij,iq)
      P01(ii,ij,2) = P0k(2,ii,ij,iq)
  100 continue
      end
