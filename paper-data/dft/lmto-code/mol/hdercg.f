      subroutine hdercg(rsm,nrx,nxi,lmxg,lxig,lmx,lxi,exi,xi,yl,
     .   kmax,cg,jcg,indxcg,n,gkl,h,h0,gh,g2h,ggh)
C- Solid smoothed Hankels and their spacial derivatives for a vector of points.
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsm   :smoothing radius
Ci   nrx   :leading dimension of yl, xi, and h's
Ci   nxi   :number of smoothed Hankel energies
Ci   lmxg  :angular-momentum dimension of xi,
Ci         :is also used to define the yl and h arrays
Ci         :lmxg is different from lmx as for derivatives we need lmaxg = lmax+2
Ci   lxig  :same as lxi but contains lmax's increased by 2 if GGA
Ci   lmx   :angular-momentum dimension of h0, gh, and g2h
Ci   lxi   :lmax for each energy
Ci   exi   :smoothed Hankel energies
Ci   xi    :smoothed Hankels xi(1..nr, 0..lxi(ie), 1..nxi)
Ci   yl    :Ylm(i,ilm), the (real) spherical harmonic polynomials
Ci   kmax  :max order of the Laplacian required ( = 0 if LDA and 2 if GGA)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   n     :actual number of points in the xy plane, n must be <= nrx
Ci   gkl   :radial part of generalised Gaussians G_kL (without energy prefactor)
Co Outputs
Co   h0    :solid smoothed Hankels H_lm = xi_l * Y_lm
Co   h     :same as h0 but containing Hankels up to lmax+kmax, a temporary array
Co   gh    :first derivatives of solid smoothed Hankels
Co   g2h   :second derivatives of solid smoothed Hankels
Co   ggh   :Lapalcian of solid smoothed Hankels
Cl Local variables
Cr Remarks
Cr The ordering of the second derivatives in g2h is xx,xy,xz,yz,zz (the third index)
Cr                                            i2 =   1  2  3  4  5
Cu Updates
Cu   08 Mar 07 (S. Lozovi) First written
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nrx,lmxg,lmx,nxi,kmax,n
      integer lxi(nxi),lxig(nxi)
      double precision rsm,exi(nxi)
      double precision yl(nrx,(lmxg+1)**2)
      double precision xi(nrx,0:lmxg,nxi),gkl(nrx,0:(kmax-1),0:lmxg)
      double precision h(nrx,0:kmax,(lmxg+1)**2)
      double precision h0(nrx,(lmx+1)**2,nxi)
      double precision gh(nrx,(lmx+1)**2,3,nxi)
      double precision g2h(nrx,(lmx+1)**2,5,nxi)
      double precision ggh(nrx,(lmx+1)**2,nxi)
c cg arrays
      integer jcg(1),indxcg(1)
      double precision cg(1)
C ... Local parameters
      integer ie,lx,lxg,nlm,ilm,l,m,i,k,i2
      integer kz,kx1,kx2,ky1,ky2
      integer lk,jlm,klm,ii,ik,m2(5)
      integer indx,icg1,icg2,icg
      double precision e,fpi,beta,th
      double precision cz,cx1,cx2,cy1,cy2
      double precision cc,c2(5)
      integer ll
      data m2 /2, -2, 1, -1, 0/

      call tcn('hdercg')
      if (nrx < n)
     . call rxi('hdercg: increase nrx, needed at least ',n)

c --------- initialize some constants
      th = 1d0/3d0
      fpi = 16d0*datan(1d0)
      if(kmax > 0) then
        cc = dsqrt(fpi/15.d0)
        c2(1) =  2.d0*cc
        c2(2) =  cc
        c2(3) =  cc
        c2(4) =  cc
        c2(5) =  2.d0*dsqrt(fpi/5.d0)
      endif

c --------- Start big loop over energies ---------
      do ie = 1, nxi
        call tcn('solid hankels')
        e = exi(ie)
        beta = fpi*dexp(0.25d0*e*rsm*rsm)
        lx=lxi(ie)
        nlm= (lx+1)**2
        lxg = lxig(ie)
        ilm = 0
        do l = 0, lxg
          do m = -l, l
            ilm = ilm+1
            do i = 1, n
              h(i,0,ilm) = xi(i,l,ie)*yl(i,ilm)
            enddo
c --- save Hankels up to lmax only ---
            if(l <= lx) call dcopy(nrx,h(1,0,ilm),1,h0(1,ilm,ie),1)
c ---- if GGA, also calculate laplacians: H_k+1,L = -e*H_kL - 4*pi*G_kL ---
            if(kmax > 0) then
              call tcn('laplacians')
              do k = 1,kmax
c We need Laplacians only up to lmax - 2*(k-1) = lmaxg - 2*k
                if(l <= lxg-2*k) then
                  do i = 1,n
                    h(i,k,ilm) =
     .                 -e*h(i,k-1,ilm) - beta*gkl(i,k-1,l)*yl(i,ilm)
                  enddo
                endif
              enddo
              call tcx('laplacians')
            endif
          enddo
        enddo
        call tcx('solid hankels')
c --- if GGA, make gradients and second derivatives of H_0,L ---
        if(kmax > 0) then
          call tcn('first derivatives')
c ------------- first derivatives ---------------
          do  ilm = 1, nlm
            call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
            do i = 1,n
              gh(i,ilm,1,ie) = gh(i,ilm,1,ie)-cx1*h(i,0,kx1)
     .                         -cx2*h(i,0,kx2)
              gh(i,ilm,2,ie) = gh(i,ilm,2,ie)-cy1*h(i,0,ky1)
     .                         -cy2*h(i,0,ky2)
              gh(i,ilm,3,ie) = gh(i,ilm,3,ie)-cz*h(i,0,kz)
            enddo
            if (ilm <= lx*lx) then
              do i = 1,n
                gh(i,kx1,1,ie) = gh(i,kx1,1,ie) - cx1*h(i,1,ilm)
                gh(i,kx2,1,ie) = gh(i,kx2,1,ie) - cx2*h(i,1,ilm)
                gh(i,ky1,2,ie) = gh(i,ky1,2,ie) - cy1*h(i,1,ilm)
                gh(i,ky2,2,ie) = gh(i,ky2,2,ie) - cy2*h(i,1,ilm)
                gh(i,kz,3,ie) =  gh(i,kz,3,ie)  - cz *h(i,1,ilm)
              enddo
            endif
          enddo
          call tcx('first derivatives')
c --------------  second derivatives ---------------
          call tcn('second derivatives')
c       call dpzero(g2h(1,1,1,ie), nrx*(lmx+1)**2*5)
        lk = 2
        do i2 = 1, 5
          jlm = lk*(lk+1)+m2(i2)+1
          cc = c2(i2)
          do ilm = 1, nlm
          l = ll(ilm)
            ii = max0(ilm,jlm)
            indx = (ii*(ii-1))/2 + min0(ilm,jlm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            do icg = icg1, icg2
              klm  = jcg(icg)
              ik = (l+lk-ll(klm))/2
              call daxpy(n,cc*cg(icg),h(1,ik,klm),1,g2h(1,ilm,i2,ie),1)
            enddo
          enddo
        enddo
c Some further adjustments
        do ilm = 1, nlm
          do i = 1,n
            g2h(i,ilm,5,ie) = (g2h(i,ilm,5,ie) +
     .                   h(i,1,ilm))*th
            g2h(i,ilm,1,ie) = (g2h(i,ilm,1,ie) +
     .                   h(i,1,ilm) -
     .                   g2h(i,ilm,5,ie))*0.5d0
          enddo
c ---  save laplacians up to lmax  ---
          call dcopy(nrx,h(1,1,ilm),1,ggh(1,ilm,ie),1)
        enddo
        call tcx('second derivatives')


c end GGA if-loop
        endif


c end big loop over energies
      enddo

      call tcx('hdercg')
      end
