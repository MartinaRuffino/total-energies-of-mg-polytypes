      subroutine dfratm(s_site,s_spec,lx,ib1,ib2,s_rhat)
C- Allocate arrays for local atomic densities.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec rho1 rho2 rhoc rho1x rho2x rhocx
Co     Stored:     *
Co     Allocated:  rho1 rho2 rhoc rho1x rho2x rhocx
Cio    Elts passed:rho1 rho2 rhoc rho1x rho2x rhocx
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  nr lmxl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     rho1 rho2 rhoc
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   lx      :if 1 allocate s_site(:)%{rho1,rho2,rhoc}
Ci           :   2 allocate s_site(:)%{rho1x,rho2x,rhocx}
Ci           :   Add 4 to de-allocate, rather than allocate
Ci           :   8 associate s_rhat(:) with s_site(:)%{rho1,rho2,rhoc}
Ci           :  16 associate s_rhat(:) with s_site(:)%{rho1x,rho2x,rhocx}
Ci           :Switches may be taken in combination, except 4,8,16 are exclusive
Ci   ib1,ib2 :allocate arrays for sites ib1..ib2
Co Outputs
Cr Remarks
Cr   rhoat(1,ib):  true local density
Cr   rhoat(2,ib):  smooth local density
Cr   rhoat(3,ib):  core density
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Jun 98 adapted from nfp df_rhoat.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lx,ib1,ib2
C ... For structures
!      include 'structures.h'
      type(str_site) :: s_site(*)
      type(str_spec) :: s_spec(*)
      type(str_rhat) :: s_rhat(*)
C ... Local parameters
      integer ib,is,lmxl,nlml,nr,nsp,nspc,nglob,nspch
      double precision xx

      nsp = nglob('nsp')
      nspc = nglob('nspc')
      nspch = nsp*nspc

      do  ib = ib1, ib2
        is = s_site(ib)%spec
        nr = s_spec(is)%nr
        lmxl = s_spec(is)%lmxl
        nlml = (lmxl+1)**2
        if (lmxl > -1) then
          if (mod(lx/4,2) == 1 .and. mod(lx,2) == 1) then
            if (associated(s_site(ib)%rho1)) then
            deallocate(s_site(ib)%rho1,s_site(ib)%rho2,s_site(ib)%rhoc)
            endif
            cycle
          endif
          if (mod(lx/4,2) == 1 .and. mod(lx/2,2) == 1) then
            if (associated(s_site(ib)%rho1x)) then
            deallocate(s_site(ib)%rho1x,s_site(ib)%rho2x,
     .                 s_site(ib)%rhocx)
            endif
            cycle
          endif
          if (mod(lx,2) == 1) then
            call ptr_site(s_site,1,'rho1',ib,nr,nlml*nspch,xx)
            call ptr_site(s_site,1,'rho2',ib,nr,nlml*nspch,xx)
            call ptr_site(s_site,1,'rhoc',ib,nr,nsp,xx)
          endif
          if (mod(lx/2,2) == 1) then
            call ptr_site(s_site,1,'rho1x',ib,nr,nlml*nspch,xx)
            call ptr_site(s_site,1,'rho2x',ib,nr,nlml*nspch,xx)
            call ptr_site(s_site,1,'rhocx',ib,nr,nsp,xx)
          endif
          if (mod(lx/8,2) == 1) then
            s_rhat(ib)%rho1 => s_site(ib)%rho1
            s_rhat(ib)%rho2 => s_site(ib)%rho2
            s_rhat(ib)%rhoc => s_site(ib)%rhoc
          elseif (mod(lx/16,2) == 1) then
            s_rhat(ib)%rho1 => s_site(ib)%rho1x
            s_rhat(ib)%rho2 => s_site(ib)%rho2x
            s_rhat(ib)%rhoc => s_site(ib)%rhocx
          endif
        endif
      enddo
      end
