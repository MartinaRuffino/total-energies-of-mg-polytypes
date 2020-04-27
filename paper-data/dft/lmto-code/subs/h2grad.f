      subroutine h2grad(mode,n,e,lmax,hl,hl1,ghl)
C- Make gradients of (smoothed) Hankel functions Hl, given Hl to one higher l
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 Use e*hl for Laplacian of Hankel (valid for unsmoothed Hankels)
Ci         :1 Use -hl1 for Laplacian of Hankel (for sm-hankels, hl1 should be H_1L)
Ci   n     :1 => hl, hl1, ghl are real (molecular case)
Ci         :2 => hl, hl1, ghl are complex (Bloch sums)
Ci   e     :energy (used if mode=0)
Ci   lmax  :maximum l
Ci   hl    :Hankels (smoothed or unsmoothed, real or complex), supplied to lmax+1
Ci   hl1   :Laplacian Hl (used if mode=1; not needed since reduces to e*hl if mode=0)
Co Outputs
Ci   ghl   :gradient of Hankel functions
Cr Remarks
Cr   A general-purpose gradient calculator of a smothed Hankel at one point.
Cr   You must supply either hl and, Hankels are smoothed, hl1.
Cu Updates
Cu   10 Apr 16  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: mode,lmax,n
      real(8) :: e,hl(n,*),hl1(n,*),ghl(n,3,*)
C ... Local parameters
      integer :: nlm,ilm,nlm1,kz,kx1,kx2,ky1,ky2,i
      double precision cz,cx1,cx2,cy1,cy2,xx

      nlm = (lmax+1)**2
      call dpzero(ghl,n*3*nlm)
      nlm1 = lmax*lmax
      do  ilm = 1, nlm
        call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
        do  i = 1, n
        ghl(i,1,ilm) = ghl(i,1,ilm) - cx1*hl(i,kx1) - cx2*hl(i,kx2)
        ghl(i,2,ilm) = ghl(i,2,ilm) - cy1*hl(i,ky1) - cy2*hl(i,ky2)
        ghl(i,3,ilm) = ghl(i,3,ilm) - cz*hl(i,kz)
        enddo
        if (ilm <= nlm1) then
          do  i = 1, n
          xx = e*hl(i,ilm)   ! More generally for smoothed H, e*H -> lap H
          if (mode == 1) xx = -hl1(i,ilm)
          ghl(i,1,kx1) = ghl(i,1,kx1) + cx1*xx
          ghl(i,1,kx2) = ghl(i,1,kx2) + cx2*xx
          ghl(i,2,ky1) = ghl(i,2,ky1) + cy1*xx
          ghl(i,2,ky2) = ghl(i,2,ky2) + cy2*xx
          ghl(i,3,kz)  = ghl(i,3,kz)  + cz*xx
          enddo
        endif
      enddo

      end
