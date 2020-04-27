      subroutine solgsg(px,rsm,lmxg,pmax,p0,gkl,ggrad)
C- Solid generalized Gaussians G_pL = (lap)^p G_L and gradient at point
C ----------------------------------------------------------------------
Ci Inputs
Ci   px    :coordinates of X relative to the centre of the Gaussian
Ci   rsm   :vector of l-dependent smoothing radii of the function
Ci         :rsm<0 => rsm, eh are independent of l.  Use |rsm(0)| for all l
Ci         :rsm>0 => rsm,eh must be specified for 0..lmax
Ci   lmxg  :max angular momentum
Ci   pmax  :max power of the Laplace operator
Ci   p0    :leading dimension of gkl
Co Outputs
Co   gkl   : G_kL(X), p = 0,...,pmax, L = 1,...,(lmax+1)^2
Co   ggrad : \nabla G_L(X), L = 1,...,(lmax+1)^2
Cr Remarks
Cu Updates
Cu   23 Aug 11 Remove cy from argument list
Cu   23 Jul 08 handles px=0 (S. Lozovoi)
Cu   22 Jun 07 gradents added (S. Lozovoi)
Cu   31 Aug 06 first written (S. Lozovoi)
C ----------------------------------------------------------------------
      implicit none
C Input parameters
      integer lmxg,pmax,p0
      double precision px(3),rsm(0:lmxg)
C Output parameters
      double precision gkl(0:p0,*),ggrad(3,*)
C Local variables
      logical fixdrs
      integer ll,ilm,ilr,ip,il,nlmp,lmax
      double precision tol
      parameter (tol=1d-14)
      double precision rsm0,rsx,r1
      double precision gklloc(0:p0,0:lmxg+1)
      double precision gloc(0:p0,(lmxg+2)**2),yl((lmxg+2)**2)
      integer kz,kx1,kx2,ky1,ky2
      double precision cz,cx1,cx2,cy1,cy2

C ... Spherical harmonics (unit radius)
      if (pmax > p0)
     .  call rxi(' solgsg: pmax exceeds dim p0, pmax=',pmax)
      if (pmax. lt. 1) call rxi
     .  ('solgsg: pmax should be at least 1 for gradients, but is',pmax)
      lmax = lmxg
C     call sylm(px,yl,lmax+1,r1)
      call ropyln(1,px(1),px(2),px(3),lmax+1,1,yl,r1)
      r1 = dsqrt(r1)

C ... if r1 = 0, only l=0 term survives in G_pL and
C                     l=1 terms survive in \grad G_L
      if (dabs(r1) < tol .and. lmxg > 1) then
        lmax = 1
        call dpzero(gkl, (p0+1)*(lmxg+1)**2)
        call dpzero(ggrad, 3*(lmxg+1)**2)
      endif
C     lmax+1 required for gradients
      nlmp = (lmax+2)**2

C ... Start big loop over smoothing radii
C     Negative 1st smoothing radius => constant rsm
      fixdrs = rsm(0) <= 0
      rsm0 = dabs(rsm(0))
      rsx = -9999
C     Decrement from highest l to simplify treatment of l-dependent rsm
      do  ilr = lmax, 0, -1
        if (.not. fixdrs) rsm0 = rsm(ilr)
        if (dabs(rsm0-rsx) > tol) then
          rsx = rsm0

          call radgkl(r1,rsm0,pmax,ilr+1,p0,gklloc)

C     ... Solid Gaussians and laplacians for given rsm0 up to curent lmax+1
          nlmp = (ilr+2)**2
          do  ilm = 1, nlmp
            il = ll(ilm)
            do ip = 0, pmax
              gloc(ip,ilm) = gklloc(ip,il)*yl(ilm)
            enddo
          enddo

C     ... Gradients in G_L up to ilr = current lmax
          call dpzero(ggrad, 3*(ilr+1)**2)

          do  ilm = 1, (ilr+1)**2
            call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
            ggrad(1,ilm) = ggrad(1,ilm) - cx1*gloc(0,kx1)
     .                                  - cx2*gloc(0,kx2)
            ggrad(2,ilm) = ggrad(2,ilm) - cy1*gloc(0,ky1)
     .                                  - cy2*gloc(0,ky2)
            ggrad(3,ilm) = ggrad(3,ilm) - cz*gloc(0,kz)

            if (ilm <= ilr*ilr) then
              ggrad(1,kx1) = ggrad(1,kx1) - cx1*gloc(1,ilm)
              ggrad(1,kx2) = ggrad(1,kx2) - cx2*gloc(1,ilm)
              ggrad(2,ky1) = ggrad(2,ky1) - cy1*gloc(1,ilm)
              ggrad(2,ky2) = ggrad(2,ky2) - cy2*gloc(1,ilm)
              ggrad(3,kz)  = ggrad(3,kz)  - cz *gloc(1,ilm)
            endif
          enddo

c     ... Copy result to G_pL, up to current lmax only
          call dcopy((p0+1)*(ilr+1)**2,gloc(0,1),1,gkl(0,1),1)

        endif
      enddo

      end
