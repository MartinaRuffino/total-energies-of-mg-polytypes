      subroutine slhhg(e,lmax,ndim,hl,ghl)
C- Gradients of solid Hankel functions at a point, with Hankels as input
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :Hankel energy
Ci   lmax  :ghl is computed for l<=lmax
Ci   ndim  :Leading dimension of ghl
Ci   hl    :Solid Hankel functions, required for l<=lmax+1
Co Outputs
Co   ghl   :Gradients of the hl
Cr Remarks
Cr   hl must be supplied to one l higher than lmax
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax,ndim
      double precision e,hl(ndim),ghl(ndim,3)
C ... Local parameters
      integer nlm,m,ilm,nlm1,kz,kx1,kx2,ky1,ky2
      double precision cz,cx1,cx2,cy1,cy2

      nlm = (lmax+1)**2
      if((lmax+2)**2 > ndim) call rx('solhg: ndim too small')

      do  20  m = 1, 3
      do  20  ilm = 1, nlm
  20  ghl(ilm,m) = 0d0

      nlm1 = lmax*lmax
      do  ilm = 1, nlm
        call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
        ghl(ilm,1) = ghl(ilm,1) - cx1*hl(kx1) - cx2*hl(kx2)
        ghl(ilm,2) = ghl(ilm,2) - cy1*hl(ky1) - cy2*hl(ky2)
        ghl(ilm,3) = ghl(ilm,3) - cz*hl(kz)
        if (ilm <= nlm1) then
          ghl(kx1,1) = ghl(kx1,1) + e*cx1*hl(ilm)
          ghl(kx2,1) = ghl(kx2,1) + e*cx2*hl(ilm)
          ghl(ky1,2) = ghl(ky1,2) + e*cy1*hl(ilm)
          ghl(ky2,2) = ghl(ky2,2) + e*cy2*hl(ilm)
          ghl(kz,3)  = ghl(kz,3)  + e*cz*hl(ilm)
        endif
      enddo
      end
