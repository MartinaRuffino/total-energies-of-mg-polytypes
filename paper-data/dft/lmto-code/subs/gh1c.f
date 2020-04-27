      subroutine gh1c(ivl,ph,rsmh,nlmh,eh,c01,pg,rsmg,nlmg,kmax,rhc,
     .  cy,cg,indxcg,jcg,ndim,gvl,dgvl,hvl,dhvl)
C- L-resolved value and Laplacian of Gaussians/Hankels centered
C  at ph over a sphere centered at pg using 1c decomposition
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :1s digit identifies the functions hvl
Ci         :10s digit identifies the functions gvl
Ci         :See Outputs for a description of these functions
Ci         :100s digit
Ci         :0 slopes dgvl,dhvl are left untouched
Ci         :1 slopes dgvl,dhvl are returned
Ci   ph    :coordinates of the head site
Ci   rsmh  :l-dependent smoothing radii for Gaussians and Hankels
Ci         :centered at the head site
Ci   nlmh  :L-cutoff for Gaussians / sm. Hankels at the head site
Ci   eh    :eh(l) = l-dependent sm. Hankel energies
Ci         :Not referenced if ivl = 0
Ci   c01   :Gaussian coefficients for generalized envelope functions
Ci   pg    :coordinates of the expansion site
Ci   rsmg  :P_kL smoothing radius at the expansion site
Ci   nlmg  :L-cutoff for P_kL at the expansion site
Ci   nds   :leading dimensions of sg
Ci   kmax  :Gaussians or sm. Hankels are decomposed into
Ci         :polynomials P_kL up to k = kmax
Ci   rhc   :l-dependent augmentation radii at the expansion site
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordan coefficients
Ci   ndim  :Leading dimension of gvl and hvl
Co Outputs
Co   hvl   :L-decomposed value and Laplacian of e0-function.
Co         :What is returned depends on 1s digit of ivl:
Co         : ivl   e0 function                         call hdxpml with:
Co         :  0    Gaussian G1                                   -
Co         :  1    sm. Hankel, Hsm                               1
Co         :  2    energy derivative of Hsm                     11
Co         :  3    Hsm + c01(0)*G0 + c01(1)*G1                   5
Co         :  4    Hsm - Hsm(rsm->0)                             2
Co         :  5    Hsm + c01(0)*G0 + c01(1)*G1 - Hsm(rsm->0)     6
Co  dhvl   :slopes corresponding to hvl
Co   gvl   :L-decomposed value and Laplacian of 2nd e0-function
Co         :What is returned depends on 10s digit of ivl:
Co         : 10s   2nd function
Co         :  0    gvl is untouched
Co         :  1    Gaussian G0
Co         :  2    Gaussian G1
Co  dgvl   :slopes corresponding to gvl
Cl Local variables
Cl  spg    :Expansion coefficients of G_kL in terms of P_kL at ph.
Cl  shr    :Expansion coefficients of H_kL in terms of P_kL at ph,
Cl         :real part
Cl  shi    :Expansion coefficients of H_kL in terms of P_kL at ph,
Cl         :imaginary part
Cl  sdr    :Expansion coefficients of dot-H_kL in terms of P_kL at ph
Cb Bugs
Cb   need to pass lmxcg to check if CG arrays are sufficiently large
Cr Remarks
Cr   gvl(L',L,p) :L'th component of G_pL centered at ph
Cr   hvl(L',L,p) :the meaning depends on ivl
Cr                ivl = 0 L'th component of G_p+1L centered at ph,
Cr                        ie hvl(L',L,p) = \lap gvl(L',L,p)
Cr                ivl = 1 L'th component of Hsm_pL centered at ph
Cr                ivl = 2 L'th component of Hsm-dot_pL centered at ph
Cr
Cr   Irrespective of a function type (Gaussians, sm. Hankels, or Hsm-dot),
Cr   gvl and hvl are obtained through the polynomial expansion of function
Cr   in question at a neighboring site. If ph and pg sites coincide, then
Cr   the actual function rather than its polynomial decomposition is used.
Cu Updates
Cu   09 Oct 11 gh1c optionally returns function derivatives
Cu   22 Sep 11 Redesigned for other kinds of envelope functions
Cu   18 Dec 08 (S. Lozovoi) ivl 10s digit
Cu   05 May 08 (S. Lozovoi) Hs and Hs-dot added
Cu   27 Feb 08 (S. Lozovoi) First written
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ivl,nlmh,nlmg,kmax,ndim
      double precision ph(3),pg(3)
      double precision rsmh(0:*),rsmg,eh(0:*),rhc(0:*),c01(0:1,*)
      double precision gvl(ndim,ndim,0:*),hvl(ndim,ndim,0:*)
      double precision dgvl(ndim,ndim,0:*),dhvl(ndim,ndim,0:*)
      integer indxcg(*),jcg(*)
      double precision cy(*),cg(*)
C ... Dynamically allocated arrays
      real(8), allocatable :: spg(:,:,:,:),shr(:,:,:,:),shdr(:,:,:,:)
C ... Local parameters
      integer n0,pmaxx
      parameter (n0=10, pmaxx=2)
      integer ndimx,kmaxx
c     parameter (ndimx=(n0+1)**2, kmaxx=n0)
      integer lmaxh,lmaxg,ilmg,pmax,ll,jvsl,jhd,im,il,ivl0,ivl1,ivl2,kmin,ikh
      double precision pkl(0:kmax,0:n0),gpkl(0:kmax,0:n0)
      double precision gkl(0:pmaxx,0:n0),hkl(0:1,0:n0),dhkl(0:1,0:n0),
     .  rsml(0:n0),rhcl(0:n0),dums(ndim,ndim)
      double precision dgkl(0:pmaxx,0:n0),ddgkl(0:1,0:n0)
C     double precision hs(0:n0),dhs(0:n0),ddhs(0:n0)
C     double precision hsp(0:n0),dhsp(0:n0),ddhsp(0:n0)
      double precision rsx,dhx,dhm0,rsm0,sss,ehvl0,ehx
      double precision g0,g1,h0,h1,dg0,dg1,dh0,dh1
      double precision tol
      parameter (tol=1d-6)

      call tcn('gh1c')
      ndimx = ndim
      kmaxx = kmax
      lmaxg = ll(nlmg)
      lmaxh = ll(nlmh)
      ivl0 = mod(ivl,10)
      ivl1 = mod(ivl/10,10)
      ivl2 = mod(ivl/100,10)

      dgkl = 0
      dhkl = 0

c...deb
c       print *,' gh1c: Input'
c       print *,' gh1c: ivl, nlmh, nlmg, kmax = ',ivl,nlmh,nlmg,kmax
c       print *,' gh1c: rsmh = ',(rsmh(il), il = 0, lmaxh)
c       print *,' gh1c: rsmg = ',rsmg
c       print *,' gh1c: eh = ',(eh(il), il = 0, lmaxh)
c       print *,' gh1c: rhc  = ',(rhc(il), il = 0, lmaxg)
c       print *,' gh1c: ndim = ',ndim
c       print *,' gh1c: ph = ',ph
c       print *,' gh1c: pg = ',pg
c       print '(a,6f12.4,5x,3f12.4)',
c    .    ' gh1c: ph, pg, ph - pg = ',ph,pg,ph-pg
c       print *,' '
c...deb
C --- Checks ---
      if (max(nlmh,nlmg) > ndimx)
     .  call rxi('gh1c: nlm exceeds ndimx. nlm = ',max(nlmh,nlmg))
      if (max(lmaxh,lmaxg) > n0)
     .  call rxi('gh1c: lmax exceeds n0. lmax = ',max(lmaxh,lmaxg))
C     if (ndim > ndimx) call rxi('gh1c: ndim is bigger than ndimx. ndim = ',ndim)
C     if (nlmh > ndimx) call rxi('gh1c: nlmh is bigger than ndimx. nlmh = ',nlmh)
C     if (nlmg > ndimx) call rxi('gh1c: nlmg is bigger than ndimx. nlmg = ',nlmh)
C     if (kmax > kmaxx) call rxi('gh1c: kmax is bigger than n0. kmax = ',kmax)

C --- Setup ---
C     sm Hankel energies must be negative for now
      if (ivl0 /= 0) then
      do  il = 0, lmaxh
        if (eh(il) > 0) call rxi('gh1c: positive energy for l = ',il)
      enddo
      endif
C     Highest order Laplacian required for 1c expansion of Gaussians
      pmax = -1   ! None
      if (ivl1 > 0) pmax = 1  ! Both G0 Laplacian G0 are needed
      if (ivl0 == 0 .or. ivl1 == 2) pmax = 2  ! G1 and Lap G1 needed

C     Argument to pass to hdxpml for first e0 (see hdxpml, hanszd, ..)
      jhd = -1                  ! illegal value
      if (ivl0 == 1) jhd =  1 ! Hsm
      if (ivl0 == 2) jhd = 11 ! Hsm-ddot
      if (ivl0 == 3) jhd =  5 ! Hsm + c01(0)*G0 + c01(1)*G1
      if (ivl0 == 4) jhd =  2 ! Hsm - Hsm(rsm->0)
      if (ivl0 == 5) jhd =  6 ! Hsm + c01(0)*G0 + c01(1)*G1 - Hsm(rsm->0)
C     if (ivl0 == 6) jhd =  7 ! Bare Hankel
      if (jhd == -1 .and. ivl0 /= 0) call rx('gh1c: illegal ivl0')

      kmin = 0; ikh = 0
      if (jhd == 5) then  ! Use  PkL + Bessel instead of pure PkL
        kmin = -1; ikh = -1; jhd = 6
      endif

      if (allocated(spg)) deallocate(spg)
      if (allocated(shr)) deallocate(shr)
      if (allocated(shdr)) deallocate(shdr)

      allocate(spg(0:kmax,ndim,0:pmaxx,ndim))
      allocate(shr(0:kmax,ndim,0:1,ndim),shdr(0:kmax,ndim,0:1,ndim))

C --- Case ph \= pg => proceed with 1c expansion ---
C     Determine whether expansion is off-site or on-site
      sss = abs(ph(1)-pg(1)) + abs(ph(2)-pg(2)) + abs(ph(3)-pg(3))
      if (sss > tol) then

C   ... Make 1c coefficients for 1st and 2nd function
        if (ivl0 >= 1) then    ! Needed unless 1st e0 has no sm-Hankels
          call hdxpml(jhd,ph,pg,eh,rsmh,c01,-rsmg,1,kmin,kmax,nlmh,
     .      nlmg,kmaxx,ndimx,1,cg,indxcg,jcg,cy,shr,shr,shdr)
        endif
        if (pmax >= 0) then    ! All cases that require Gaussians
          call gxpml(ph,pg,rsmh,-rsmg,pmax,kmax,nlmh,nlmg,kmaxx,ndimx,
     .      pmaxx,cg,indxcg,jcg,cy,spg)
        endif

C   ... Values and Laplacians, 1st e0
        jvsl = 11
C       hvl <- value and Laplacian of G1
        if (ivl0 == 0) then
          call vsl(jvsl,rhc,rsmg,0,kmax,0,0d0,lmaxg,(lmaxh+1)**2,kmaxx,
     .      ndimx,pmaxx,spg(0,1,1,1),hvl,dums,hvl(1,1,1),pkl,gpkl)
C       hvl <- value and Laplacian of sm-Hdot
        elseif (ivl0 == 2) then
          call vsl(jvsl,rhc,rsmg,0,kmax,0,0d0,lmaxg,(lmaxh+1)**2,kmaxx,
     .      ndimx,1,shdr,hvl,dums,hvl(1,1,1),pkl,gpkl)
C       hvl <- value and Laplacian of functions involving hsm
        else
          call vsl(jvsl,rhc,rsmg,kmin,kmax,ikh,eh,lmaxg,(lmaxh+1)**2,kmaxx,
     .      ndimx,1,shr,hvl,dums,hvl(1,1,1),pkl,gpkl)
        endif
        if (ivl2 == 1) then
          call dcopy(ndim*ndim,dums,1,dhvl,1)
        endif

C   ... Values and Laplacians, 2nd e0
        jvsl = 10
        if (ivl1 == 0) then
C       gvl <- value and Laplacian of G0
        elseif (ivl1 == 1) then
          call vsl(jvsl,rhc,rsmg,0,kmax,0,0d0,lmaxg,(lmaxh+1)**2,kmaxx,
     .      ndimx,pmaxx,spg,gvl,dums,gvl(1,1,1),pkl,gpkl)
          jvsl = 10
C       gvl <- value and Laplacian of G1
        elseif (ivl1 == 2) then
          call vsl(jvsl,rhc,rsmg,0,kmax,0,0d0,lmaxg,(lmaxh+1)**2,kmaxx,
     .      ndimx,pmaxx,spg(0,1,1,1),gvl,dums,gvl(1,1,1),pkl,gpkl)
        endif
        if (ivl2 == 1 .and. ivl1 /= 0) then
          call dcopy(ndim*ndim,dums,1,dgvl,1)
        endif

C --- Case ph = pg => do explicitly (only diagonal terms needed)
      else

C   ... Handle negative smoothing radii and rmt
        if (rsmh(0) < 0d0) then
          call dvset(rsml(0),1,lmaxh+1,-rsmh(0))
        else
          call dcopy(lmaxh+1,rsmh(0),1,rsml(0),1)
        endif
        if (rhc(0) < 0d0) then
          call dvset(rhcl(0),1,lmaxh+1,-rhc(0))
        else
          call dcopy(lmaxh+1,rhc(0),1,rhcl(0),1)
        endif

C   ... Make on-site Gaussians for each l
        rsx = -1d2
        dhx = -1d2
        if (pmax >= 0) then
        do  il = lmaxh, 0, -1
          dhm0 = rhcl(il)
          rsm0 = rsml(il)
          if (dabs(rsm0-rsx)*0 + dabs(dhm0-dhx) > tol) then
C           call radgkl(dhm0,rsm0,pmax,il,pmaxx,gkl)
C           scale Gaussians by r^l
C            if (il >= 1) then
C              fac = 1
C              do  ill = 1, il
C                fac = fac*dhm0
C                do  ip = 0, pmax
C                  gkl(ip,ill) = gkl(ip,ill)*fac
C                enddo
C              enddo
C            endif
            call radgkg(0,dhm0,rsml,pmax,il,pmaxx,gkl,dgkl,ddgkl)
          endif
          rsx = rsm0
          dhx = dhm0
        enddo
        endif
C       if ivl0 > 0, also make on-site Hankels or H-dots
        if (ivl0 > 0) then
          rsx = -1d2
          dhx = -1d2
          ehx = 1d2
          do  il = lmaxh, 0, -1
            dhm0 = rhcl(il)
            rsm0 = rsml(il)
            ehvl0 = eh(il)
C           Branch executes only if parameters changed from prior il
            if (dabs(rsm0-rsx) + dabs(dhm0-dhx) + dabs(ehvl0-ehx) > tol) then
C              OLD: works for ivl=1,2 only
C              call hanszd(jhd+1,dhm0,ehvl0,-rsm0,il,
C     .          hs,dhs,ddhs,hsp,dhsp,ddhsp)
C              if (ivl0 == 1) then
C                do  ill = 0, il
C                  hkl(0,ill) = hs(ill)
C                  hkl(1,ill) = ddhs(ill)
C                enddo
C              else
C                do  ill = 0, il
C                  hkl(0,ill) = hsp(ill)
C                  hkl(1,ill) = ddhsp(ill)
C                enddo
C              endif
              im = 100 + 10*ivl2 + ivl0-1
              call sole0g(dhm0,-rsm0,ehvl0,il,im,1,c01,hkl,dhkl)
            endif
            rsx = rsm0
            dhx = dhm0
            ehx = ehvl0
          enddo
        endif

C   ... Fill in gvl and hvl
        ilmg = 0
        do  il = 0, lmaxh
          g0 = gkl(0,il)
          g1 = gkl(1,il)
          dg0 = dgkl(0,il)
          dg1 = dgkl(1,il)
          if (ivl0 == 0) then
            h0 = g1
            h1 = gkl(2,il)
            dh0 = dg1
            dh1 = dgkl(2,il)
          else
            h0 = hkl(0,il)
            h1 = hkl(1,il)
            dh0 = dhkl(0,il)
            dh1 = dhkl(1,il)
          endif
          do  im = -il, il
            ilmg = ilmg + 1
            if (ivl1 >= 1) then
            call dpzero(gvl(1,ilmg,0),nlmh)
            call dpzero(gvl(1,ilmg,1),nlmh)
            gvl(ilmg,ilmg,0) = g0
            gvl(ilmg,ilmg,1) = g1
            endif
            call dpzero(hvl(1,ilmg,0),nlmh)
            call dpzero(hvl(1,ilmg,1),nlmh)
            hvl(ilmg,ilmg,0) = h0
            hvl(ilmg,ilmg,1) = h1
            if (ivl2 == 1) then
            if (ivl1 >= 1) then
            call dpzero(dgvl(1,ilmg,0),nlmh)
            call dpzero(dgvl(1,ilmg,1),nlmh)
            dgvl(ilmg,ilmg,0) = dg0
            dgvl(ilmg,ilmg,1) = dg1
            endif
            call dpzero(dhvl(1,ilmg,0),nlmh)
            call dpzero(dhvl(1,ilmg,1),nlmh)
            dhvl(ilmg,ilmg,0) = dh0
            dhvl(ilmg,ilmg,1) = dh1
            endif

          enddo
        enddo

      endif

      call tcx('gh1c')
      end
