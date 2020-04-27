      subroutine gfg2gr(mode,s_site,nbas)
C- Convert the local g_00 to G_00 by energy scaling (Dirac or spin-orbit)
C-----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb
Co     Stored:     gc
Co     Allocated:  *
Cio    Elts passed:gcu dpfr ddpfr bxc thet gc
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :1 Dirac
Ci         :2 spin-orbit coupling
Ci         :10s digit
Ci         :1 apply rotation specified by constraining fields
Ci   nbas  :number of sites
Ci   gcu   :On input s_site(ib)%gcu is unscaled g, eg (P-S)^-1 in lms representation
Co Outputs
Co   gc    :Scaled gc is written to s_site(ib)%gc
Co         :In the Dirac case (mode=1) the representation is rotated from lms to lambda-mu
Cr Remarks
Cu Updates
Cu  08 Jun 14 (Belashchenko) Extended to the relativistic case
Cu  24 Oct 13 (Belashchenko) First created
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer mod0,mod1,nspc,ib,ncomp,norb,n2,nl,icp,mxorb,nglob,idum
      complex(8), allocatable, dimension(:,:,:) :: g,gkmu,dpfr,ddpfr
      real(8) eula(3)
      real(8), allocatable :: bxc(:,:),th(:,:)

      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)
      call sanrg(.true.,mod0,1,2,'gfg2gr:','mod0')
      call sanrg(.true.,mod1,0,1,'gfg2gr:','mod1')
      mxorb = nglob('mxorb')
      nspc = 2

      do  ib = 1, nbas
        ncomp = s_site(ib)%ncomp ; norb = s_site(ib)%norb ; n2 = 2*norb
        allocate(g(n2,n2,ncomp),dpfr(n2,n2,ncomp),ddpfr(n2,n2,ncomp))

C   ... Copy local g from s_site(ib)%gcu
        call zcopy(n2*n2*ncomp,s_site(ib)%gcu,1,g,1)

C   ... If Dirac, convert local g from lms rep to gkmu in kappa-mu rep
        if (mod0 == 1) then
          allocate(gkmu(n2,n2,ncomp))
C         call zprm('gfg2gr: g(lms) before scaling',2,g,n2,n2,n2)
          call mstokm(0,0,ncomp,norb,g,gkmu,idum)
          call zcopy(n2*n2*ncomp,gkmu,1,g,1)
C         call yprmi('gfg2gr: g(kmu) for ib=%i',ib,0,3,g,0,n2,n2,n2)
          deallocate(gkmu)
        endif

C   ... Scale to make G_{lambda1,lambda2} and record in s_site(ib)%gc
        call zcopy(norb*norb*4*ncomp,s_site(ib)%dpfr,1,dpfr,1)
        call zcopy(norb*norb*4*ncomp,s_site(ib)%ddpfr,1,ddpfr,1)
C       call zprm('dpfr',2,dpfr,n2,n2,n2)
        if (mod1 == 1) then
          allocate(bxc(3,ncomp),th(ncomp,2))
          call dcopy(3*ncomp,s_site(ib)%bxc,1,bxc,1)
          call dcopy(2*ncomp,s_site(ib)%thet,1,th,1)
          nl = sqrt(norb + 0.1)
        endif
        do icp = 1, ncomp
          g(:,:,icp) = ddpfr(:,:,icp) + matmul(transpose(dpfr(:,:,icp)),matmul(g(:,:,icp),dpfr(:,:,icp)))
C    ...  Rotate if required (constraining fields)
          if (mod1 == 1) then
            eula(1) = th(icp,2)
            eula(2) = datan(bxc(1,icp)); eula(3) = - eula(1)
            call rot_LandS(10,eula,nl,norb,1,g(:,:,icp))
          endif
        enddo
        if (allocated(bxc)) deallocate(bxc,th)
        call zcopy(norb*norb*4*ncomp,g,1,s_site(ib)%gc,1)

C       Debugging
C       allocate(gkmu(n2,n2,ncomp))
C       call mstokm(1,0,ncomp,norb,gkmu,g,idum)
C       call zprm('gfg2gr: g(lms) after scaling',2,gkmu,n2,n2,n2)
C       deallocate(gkmu)

        deallocate(g,dpfr,ddpfr)
      enddo

      end
