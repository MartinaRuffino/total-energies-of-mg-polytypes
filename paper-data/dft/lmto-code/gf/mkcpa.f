      subroutine mkcpa(s_ctrl,s_site,s_pot,s_ham,s_lat,lspec,izp)
C- Make CPA coherent potential
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham nbas nl lrel ldomg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lbxc
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb domg omg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:amr amrt alr pfr dpfr ddpfr thet cpawt
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  pfr pf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pf cp dpf ddpf pfr
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lncol iprmb offH
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:istab symgr ag
Cio    Passed to:  *
Ci Inputs
Ci   lspec :?
Ci   izp   :complex energy number for access to Omega
Co Outputs
Co   cp    :coherent potential, in iprm order, stored in s_pot%cp
Cu Updates
Cu   18 Jun 18 Supersede call to pfr2block with call to new pokepf
Cu   09 Jul 14 (Belashchenko) pivotp only for lrel=2, otherwise use s_site
Cu   19 Aug 13 (Belashchenko) Added the relativistic case
Cu   18 Dec 12 Completed migration to F90 structures
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
Cu   18 Jan 12 (Belashchenko) Major cleanup
Cu   15 May 08 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer izp
      logical lspec
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_site)::  s_site(*)
C ... Dynamically allocated arrays
      complex(8), dimension(:,:,:,:),allocatable :: pfi,dpfi,ddpfi,pmomg
C     complex(8), pointer :: pf(:)
C ... Local parameters
      integer nbas,nl,nsp,norb,nspc,ib,ncomp,optrot,lrel,iopt,iproff
      integer nsgrp,morder
      integer ldham(16),lham,ldim,lidim,lihdim,pfdim,ioffd,mxorb,nglob
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (nspc,ldham(4))
      equivalence (pfdim,ldham(9))
      logical lso,bittst
      complex(8), pointer :: p_omg(:), p_pf(:,:)


      mxorb = nglob('mxorb')
      nsp = nglob('nsp')
      lham = s_ctrl%lham
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      lrel = mod(s_ctrl%lrel,10)
      ldham = s_ham%ldham
      nsgrp = s_lat%nsgrp
      lso = bittst(s_ham%lncol,4)
      call lmorder(0,morder,[0],[0])

      ioffd = lihdim
C ... Copy potpars to diagonal of cp for non-DLM sites (NB: iprmb order only
C                                                       for lrel !=2 !)
C ... (assume that s_pot%cp has been initialized to zero when allocated)

      if (lrel /= 2 .and. .not. lso) call dcpdiag(lidim,nspc,pfdim,s_pot%pf,s_pot%cp)

C --- For each site, do
      do  ib = 1, nbas

        ncomp = s_site(ib)%ncomp ; norb = s_site(ib)%norb
        iproff = 1 + mxorb*(ib-1)
        if (ncomp == 1) then
          if (lspec) then
            if (lrel == 2) call rx('specfun not ready for lrel=2')
            call prepsf0(norb,nspc,s_site(ib),lidim,
     .        s_ham%iprmb(iproff:iproff+norb-1),pfdim,s_pot%dpf,
     .        s_pot%ddpf,s_site(ib)%amr,s_site(ib)%amrt,s_site(ib)%alr)
          endif

C         This branch does not yet work for rel=2 case
C         because s_site(ib)%pfr is stored in kappa-mu repsn.  Needed: rotation to lms
          if (lso .and. lrel /= 2) then
            call pokeg0(0,norb,s_ham%iprmb(iproff),lidim,lidim,lidim,nspc,2,
     .        s_site(ib)%pfr,s_pot%cp)
          elseif (lrel == 2) then
C       ... Copy pfr for non-CPA sites directly to cp (check!)
C           s_pot%cp has the rank of the full hamiltonian offset for site ib
C           via s_ham%iprmb(iproff) or via offp=mxorb*(ib-1), but only one of them
C           pfr2block has been superseded by pokepf
C            call pfr2block(100*morder+11,nl,s_ham%iprmb,0,lidim,
C     .        mxorb*(ib-1),pfdim,s_pot%pfr,lidim,1,s_pot%cp)
            call pokepf(1+10*morder,1,pfdim,lidim,0,lidim,s_ham%iprmb,
     .        mxorb*(ib-1),ib,ib,nspc,nsp,s_pot%pfr,s_pot%cp)

          endif
C         Debugging check: copy pfib gen by mkfrpt for sites 1 and 2 to files 1 and 2.
C         mc -f16f12.6 -1:16 -s0 -a z 1 -split a 1:nr+1:16 1:nc+1:16 -pop 2 -split b 1:nr+1:16 1:nc+1:16 -pop \
C                      a11 z a12 z z b11 z b12 a21 z a22 z z b21 z b22 '-sarray[4,4]' out.fept -- -px
C         call yprmi('s_pot%%%%cp after ib=%i',ib,0,3,s_pot%cp,0,nspc*lihdim,nspc*lihdim,nspc*lihdim)
          cycle
        endif
        if (s_ctrl%ldomg == 1) then ; p_omg => s_site(ib)%domg(:,izp)
        else ; p_omg => s_site(ib)%omg(:,izp) ; endif
        if (.not. lso) then
          allocate(pfi(norb,2,nspc,ncomp))
        else
          allocate(pfi(1,1,1,ncomp))
        endif
        if (lrel == 2) then
          p_pf => s_pot%pfr
        elseif (.not. lso) then
          p_pf => s_pot%pf
        endif
        if (lrel == 2) then
          call pivotp(ioffd,nspc,norb,ncomp,pfdim,p_pf,pfi)
        elseif (.not. lso) then
          call zcopy(norb*2*ncomp,s_site(ib)%pfr,1,pfi,1)
        endif
C ...   Fix the following for lso
        if (lspec) then
          allocate(pmomg(norb,nspc,norb,2))
          if (lrel == 2) then
            allocate(dpfi(norb,2,nspc,ncomp))
            allocate(ddpfi(norb,2,nspc,ncomp))
            call pivotp(ioffd,nspc,norb,ncomp,pfdim,s_pot%dpf,dpfi)
            call pivotp(ioffd,nspc,norb,ncomp,pfdim,s_pot%ddpf,ddpfi)
          elseif (.not. lso) then
            allocate(dpfi(norb,2,1,ncomp))
            allocate(ddpfi(norb,2,1,ncomp))
            call zcopy(norb*2*ncomp,s_site(ib)%dpfr,1,dpfi,1)
            call zcopy(norb*2*ncomp,s_site(ib)%ddpfr,1,ddpfi,1)
          else
            allocate(dpfi(1,1,1,ncomp),ddpfi(1,1,1,ncomp))
          endif
          iopt = 3
        else
          if (.not.allocated(pmomg)) allocate(pmomg(1,1,1,1))
          iopt = 1
        endif
        if (bittst(s_ctrl%lbxc,4)) iopt = iopt + 10
        call cpadlm(iopt,nspc,lrel,lso,s_site(ib),norb,pfi,p_omg,
     .    s_pot%cp,pmomg,lidim,s_ham%iprmb(iproff))
        if (lspec) then
          iopt = 0
          if (bittst(s_ctrl%lbxc,4)) iopt = 1
          call prepsf(iopt,ncomp,nspc,s_site(ib),s_site(ib)%thet,
     .      s_site(ib)%cpawt,norb,pfi,dpfi,ddpfi,p_omg,pmomg,
     .      s_site(ib)%amr,s_site(ib)%amrt,s_site(ib)%alr)
        endif
        deallocate(pfi)
        if (lspec) deallocate(pmomg,dpfi,ddpfi)
        ioffd = ioffd + norb*ncomp
      enddo

C --- Symmetrize the coherent potential
C ... Rotate site-diagonal cp, other options see gfibz
      if (s_ctrl%nccomp /= 0 .and. lrel /= 2) then
        optrot = 1
        if (bittst(lham,256)) optrot = optrot + 20
C       There may be a bug in roth for lidim > ldim - check!
C       (Update: There seems to be no bug but be careful!)
        if (lidim > ldim) optrot = optrot + 4000
        call symcp(optrot,s_pot%cp,lidim,nspc,nbas,nl,s_ham%offH,
     .    s_ham%iprmb,nsgrp,s_lat%istab,s_lat%symgr,s_lat%ag,2,1)
      else
c       if (izp == 1)
c    .    print *,'WARNING(mkcpa): cohP not symmetrized for lrel=2'
      endif

      end

      subroutine dcpdiag(lidim,nspc,pfdim,pf,cp)
C- Copy pf to diagonal part of coherent potential array
      integer lidim,nspc,pfdim,i
      double complex pf(pfdim,2),cp(lidim,nspc,lidim,2)
      if (pfdim < lidim) call rx('dcpdiag: dimensioning error')
      do  i = 1, lidim
        cp(i,1,i,1) = pf(i,1) ; cp(i,nspc,i,2) = pf(i,2)
      enddo
      end

      subroutine pivotp(ioff,nspc,norb,ncomp,pfdim,pf,pfi)
C- Copy potential parameters for one site from pf to pfi
      implicit none
      integer ioff,nspc,norb,ncomp,pfdim
      double complex pf(pfdim,2,nspc),pfi(norb,2,nspc,ncomp)

      integer i,n,nc

      pfi = dcmplx(0d0,0d0)
      i = ioff
      do nc = 1, ncomp
        do n = 1, norb
          i = i + 1
          pfi(n,:,:,nc) = pf(i,:,:)
        enddo
      enddo
      end


