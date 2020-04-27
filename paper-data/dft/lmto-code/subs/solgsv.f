      subroutine solgsv(rsm,nrx,nl,lmax,pmax,nr,rp2,yl,gex)
C- Solid, energy-independent generalized Gaussians for a vector of points
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsm   :vector of l-dependent smoothing radii of Gaussians
Ci         :EITHER must be specified for lmin..lmax
Ci         :OR     rsmh(0) = const < 0. Implies rsmh(l) = -const for all l
Ci   nrx   :max number of points r at which the Gaussians are evaluated;
Ci         :leading dimension of yl and gex
Ci   nl    :max angular momentum + 1; second dimension of gex
Ci   lmax  :G_pL are evaluated up to l = lmax
Ci   pmax  :G_pL are evaluated up to p = pmax
Ci   nr    :actual number of points r at which the Gaussians are evaluated
Ci   rp    :coordinates of points r relative to the centre of Gaussian
Ci   rp2   :squared distance from the centre of Gaussian; rp2 = rp.rp
Ci   yl    :vector of spherical harmonics corresponding to rp (ropyln.f)
Co Outputs
Co   gex   : G_pL(r), p = 0,...,pmax, L = 1,...,(lmax+1)^2
Cr Remarks
Cr     G_pL := (lap)^p G_L
Cr     Note that p and l indecies in solid and radial Gaussians arrays are
Cr     interchanged
Cb Bugs
Cu Updates
Cu   24 Mar 09 (S. Lozovoi) created from solgsg.f
C ----------------------------------------------------------------------
      implicit none
C Input parameters
      integer nrx,nl,lmax,pmax,nr
      double precision rp2(nr),rsm(0:lmax)
      double precision yl(nrx,nl**2)
C Output parameters
      double precision gex(nrx,nl**2,0:pmax)
C Local variables
      integer n0
      parameter (n0 = 10)
      integer ierr
      integer ilm,ip,il,im
      double precision rsml(0:n0),rsm0,rsx
      double precision, allocatable :: gkl(:,:,:)
      double precision tol
      data tol/1.d-28/

      if(lmax. gt. n0)
     . call rxi(' solgsv: Increase n0 up to at least ',lmax)
      if(pmax. lt. 0) call rxi(' solgsv: negative pmax = ',pmax)

C ... Allocate gkl
      if (.not. allocated(gkl)) then
        allocate(gkl(nr,0:pmax,0:lmax), stat=ierr)
        if(ierr /= 0) then
          call rx0(' solgsv: allocation of gkl failed')
        else
          call info5(0,1,0,' solgsv: gkl(%i,0:%i,0:%i) allocated',
     .               nr,pmax,lmax,0,0)
        endif
      endif

C ... Handle negative smoothing radii
      if (rsm(0) < 0d0) then
        call dvset(rsml(0),1,lmax+1,-rsm(0))
      else
        call dcopy(lmax+1,rsm(0),1,rsml(0),1)
      endif

C ... Make radial parts g_kl
      call tcn('radial gaussians')
      rsx = -1d2
      do il = lmax, 0, -1
        rsm0 = rsml(il)
        if (dabs(rsm0-rsx) > tol) then
          call radgkv(nr,pmax,rp2,rsm0,nr,pmax,il,gkl)
          rsx = rsm0
        endif
      enddo
      call tcx('radial gaussians')

c ... make solid Gaussians and their Laplacians up to p = pmax
      call tcn('G = g*Y')
      ilm = 0
      do  il = 0, lmax
        do im = -il, il
        ilm = ilm + 1
          do ip = 0, pmax
            gex(1:nr,ilm,ip) = gkl(1:nr,ip,il)*yl(1:nr,ilm)
          enddo
        enddo
      enddo
      call tcx('G = g*Y')

C --- The bit below is unnecessary ---
c ... if r = 0, only l=0 term should survive in G_pL
c (might already be contained within yl)
c     if (lmax >= 1) then
c       do ir = 1, nr
c         if(dabs(rp2(ir)) < tol) then
c            do ip = 0, pmax
c              do ilm = 2, (lmax+1)**2
c                gex(ir,ilm,ip) = 0d0
c              enddo
c            enddo
c         endif
c       enddo
c     endif

      deallocate(gkl)
      end
