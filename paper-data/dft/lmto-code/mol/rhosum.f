      subroutine rhosum(nrx,nxi,lxi,nri,rhoi,lmx,ioff,
     .   h,gh,g2h,ggh,n,rho,grho,g2rho,ggrho)
C- Interstitial density in a xy plane due to a single atom (vectorizes).
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrx     :leading dimension of sm. Hankels, their derivatives: h, gh, ggh, and g2h
Ci   nxi     :number of Hankel energies
Ci   lxi     :lmax for each energy
Ci   nri     :total number of functions in cd basis: dimensions of rhoi
Ci   rhoi    :coefficients to density in atom-centered density basis
Ci   lmx     :max. angular momentum, defines second dimension of h, gh, etc.
Ci   ioff    :table of offsets to rhoi for each energy
Ci   h       :sm. solid Hankels for each point, lm, and energy
Ci   gh      :gradient of sm. solid Hankels (third index runs over x, y, z)
Ci   g2h     :matrix of second derivatives of sm. solid Hankels (third index:
Ci           : 1 - xx, 2 - xy, 3 - xz, 4 - yz, 5 - zz, see hdercg.f)
Ci   ggh     :Laplacian of sm. solid Hankels
Ci   yl      :Ylm(i,ilm), the (real) spherical harmonic polynomials
Ci   xi      :smoothed Hankels xi(1..nr, 0..lxi(ie), 1..nxi)
Ci   n       :actual number of points in the xy plane, n must be <= nrx
Co Outputs
Co   rho    :interstitial density due to a given atom
Co   grho   :gradient of interstitial density due to a given atom
Co   ggrho  :Laplacian of interstitial density due to a given atom
Co   g2rho  :matrix of second derivatives of interstitial density
Cl Local variables
Cr Remarks
Cr gh, g2h, and ggh are not referenced for LDA. Same for grho, g2rho, and ggrho.
Cu Updates
Cu   21 Mar 07 (S. Lozovi) Rewritten to include density derivatives
Cu   21 Apr 02 (S. Lozovi) First created
C ----------------------------------------------------------------------
      implicit none
      integer nrx,lmx,nxi,nri,n,lxcg
      integer lxi(nxi),ioff(nxi)
      integer ie,joff,lx,l,m,ilm,isp,nsp,lsp,ia
      double precision   h(nrx,(1+lmx)**2,nxi)
      double precision  gh(nrx,(1+lmx)**2,3,nxi)
      double precision ggh(nrx,(1+lmx)**2,nxi)
      double precision g2h(nrx,(1+lmx)**2,5,nxi)
      double precision rhoi(nri,*),rho(nrx,*),
     .  grho(nrx,3,*),g2rho(nrx,5,*),ggrho(nrx,*)
      double precision cc
      call tcn('rhosum')
      nsp = lsp()+1

      call dpzero(rho,  nrx*nsp)
c ---- if GGA, initialise the arrays for density derivatives -----
      if(lxcg() > 0) then
        call dpzero(ggrho,  nrx*nsp)
        call dpzero(grho,  nrx*3*nsp)
        call dpzero(g2rho, nrx*5*nsp)
      endif

c ---- start big loop over energies -----
      do ie = 1,nxi
        joff = ioff(ie)
        lx = lxi(ie)
c       nlm = (lx+1)**2
        ilm = 0
        do l = 0,lx
c         lav = l*(l+1)+1
          do m = -l,l
            ilm = ilm+1
            do isp = 1,nsp
c             cc=rhoi(lav+m+joff,isp)
              cc=rhoi(ilm+joff,isp)
              call daxpy(n,cc,h(1,ilm,ie),1,rho(1,isp),1)
c ---- if GGA, make also density derivatives -----
              if(lxcg() > 0) then
                call daxpy(n,cc,ggh(1,ilm,ie),1,ggrho(1,isp),1)
                do ia = 1,3
                  call daxpy(n,cc,gh(1,ilm,ia,ie),1,grho(1,ia,isp),1)
                enddo
                do ia = 1,5
                  call daxpy(n,cc,g2h(1,ilm,ia,ie),1,g2rho(1,ia,isp),1)
                enddo
              endif
c             do i = 1,n
c               rho(i,isp) = rho(i,isp)+cc*h(i,ilm,ie)
c             enddo
            enddo
          enddo
        enddo
c ---- end big loop over energies -----
      enddo

      call tcx('rhosum')
      end
