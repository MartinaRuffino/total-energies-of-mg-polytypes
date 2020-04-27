      subroutine mkgint(mode,s_ctrl,s_site,s_pot,s_ham,izp,wz,dz,ferm)
C- Conditionally-averaged GF for all DLM sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbasp lrel ldomg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lbxc
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb ncomp omgn omg
Co     Stored:     *
Co     Allocated:  sfvrtx
Cio    Elts passed:domg gc gcu gcorr sfvrtx
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  pfr dpfr ddpfr pf dpf ddpf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cp
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lncol
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do not return any of : unscaled g, scaled G, gcorr
Ci         :1 Return scaled G in s_site()%gc for each component
Ci         :2 Return unscaled g in s_site()%gcu for each component
Ci         :4 Return Lloyd correction in s_site()%gc
Ci         :  Any combination is allowed
Ci         :10s digit
Ci         :0 do not calculate vertex
Ci         :1 calculate vertex
Ci            For vertex see K. Belashchenko's notes dated 18 June 2012
Ci            and PRB 81, 064410
Ci         :100s digit
Ci         :0 use s_site%omgn for interactor
Ci         :1 use s_site%omg for interactor
Ci         :1000s digit
Ci         :0 Do not make one-site exchange parameters J0
Ci         :1 Make J0, contracting by L (one J0 for each component)
Ci         :2 Make J0, retaining L resolution
Ci   izp   :complex energy index for access to Omega
Ci   wz    :integration weight for this energy point
Ci   dz    :finite-difference delta z to calculate dOmega/dz
Ci         :(used only when calculating Lloyd correction to sev)
Ci   ferm  :true, if point is at fermi energy
Co Outputs
Co   gc,gcorr,j0,sfvrtx are generated and stored in s_site, depending on mode
Cu Updates
Cu   17 Jan 13 Redesigned switches and added spin flip vertex
Cu   18 Dec 12 Completed migration to F90 structures
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
Cu   24 Dec 11 (Belashchenko) Complete rewrite
Cu   08 Dec 08 (P. Larson) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed parameters:
      integer mode,izp
      double precision dz
      double complex wz
      logical ferm
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_site)::  s_site(*)
C Local parameters:
      integer nbasp,norb,nspc,ib,ncomp,i,mkgmod,cpamod,lcgf
      integer ldham(16),lrel,lidim,lihdim,pfdim,ioffd
      equivalence (lidim,ldham(2)),(lihdim,ldham(3)),(pfdim,ldham(9))
      logical lsev,ldomg,lvertex,lomg,lgu,lso,bittst
      complex(8), pointer :: omega(:),p_pf(:,:),p_dpf(:,:),p_ddpf(:,:)
      complex(8), allocatable :: pfi(:),dpfi(:),ddpfi(:),pmomg(:,:,:,:)
      double precision xx
      procedure(integer) nglob


C --- Setup ---
      nbasp = s_ctrl%nbasp
      ldham = s_ham%ldham
      nspc = nglob('nspc')
      lrel = mod(s_ctrl%lrel,10)
      ldomg = s_ctrl%ldomg == 1
      lsev = mod(mode,10) >= 2
      lvertex = mod(mode/10,10) == 1
      lgu = mod(mod(mode,10)/2,2) /= 0
      lomg  = mod(mode/100,10) == 1
      lso = bittst(s_ham%lncol,4)
      lcgf = s_ctrl%lcgf
      mkgmod = mode
      cpamod = 2
      if (bittst(s_ctrl%lbxc,4)) then
        mkgmod = mkgmod + 10000
        cpamod = cpamod + 10
      endif
C     print *, '!!'; lvertex = .true.
      if (ldomg .and. lsev)
     .  call rx('mkgint: ldomg and lsev can not be combined')

C --- For each CPA site, do ---
      ioffd = lihdim ! offset to DLM part of pf for current site
      do  ib = 1, nbasp
        norb = s_site(ib)%norb
        ncomp = s_site(ib)%ncomp
        if (ncomp < 2) then   ! Not a CPA site, copy gii to gcu
          call zcopy(norb*nspc*norb*2,s_site(ib)%gii,1,s_site(ib)%gcu,1)
          cycle
        endif
C       Choose which omega supplied by s_site ("new" or "old")
        omega => s_site(ib)%omgn(:,izp)
        if (lomg) omega => s_site(ib)%omg(:,izp)
C       Extract potential parameters
        i = ncomp*norb*2*nspc
        allocate(pfi(i),dpfi(i),ddpfi(i))
        if (lrel == 2) then
          p_pf => s_pot%pfr; p_dpf => s_pot%dpfr; p_ddpf => s_pot%ddpfr
        else
          p_pf => s_pot%pf; p_dpf => s_pot%dpf; p_ddpf => s_pot%ddpf
        endif
        if (lrel == 2) then
          call pivotp(ioffd,1+lrel/2,norb,ncomp,pfdim,p_pf,pfi)
          call pivotp(ioffd,1+lrel/2,norb,ncomp,pfdim,p_dpf,dpfi)
          call pivotp(ioffd,1+lrel/2,norb,ncomp,pfdim,p_ddpf,ddpfi)
        elseif (.not. lso) then
          call zcopy(norb*2*ncomp,s_site(ib)%pfr,1,pfi,1)
          call zcopy(norb*2*ncomp,s_site(ib)%dpfr,1,dpfi,1)
          call zcopy(norb*2*ncomp,s_site(ib)%ddpfr,1,ddpfi,1)
        endif

C   ... Setup for vertex
        if (lvertex) then
          allocate(pmomg(norb,nspc,norb,2))
C    ...  s_pot%cp is not changed with cpadlm mode 2
          call cpadlm(cpamod,nspc,lrel,lso,s_site(ib),norb,pfi,omega,
     .      s_pot%cp,pmomg,lidim,0)
          call ptr_site(s_site,1,'sfvrtx',ib,norb*norb,ncomp*2,xx)
          call ptr_site(s_site,1,'vtxrel',ib,4*4*4*norb*norb,ncomp,xx)
        else  ! The DEC alpha compiler needs pointer to be allocated
C         nullify(pmomg)
          allocate(pmomg(1,1,1,1))
        endif

C   ... Calculate g = (P - Omega)^-1; return g scaled to G.
C       and optionally the vertex
        if (nspc == 2) lrel = 2
        call mkgcpa(mkgmod,nspc,lrel,lso,lcgf,s_site(ib),norb,izp,wz,ncomp,
     .    omega,s_site(ib)%domg(:,izp),pmomg,pfi,dpfi,ddpfi,
     .    s_site(ib)%gc,s_site(ib)%gcu,s_site(ib)%gcorr,
     .    s_site(ib)%sfvrtx,s_site(ib)%vtxrel,dz,ferm)
        deallocate(pfi,dpfi,ddpfi,pmomg)
        ioffd = ioffd + norb*ncomp
C       print *, 'mkgint',izp,sngl(s_site(1)%j0(1:2,1:2))
      enddo
      end

      subroutine mkgcpa(mode,nspc,lrel,lso,lcgf,s_site,norb,izp,wz,ncp,omg,
     .  domg,pmomg,pfi,dpf,ddpf,gc,gcu,gcorr,sfvrtx,sfvrtxrel,dz,ferm)
C- Conditionally averaged g=(P-Omega)^-1 and G for each DLM component
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :1s digit
Ci         :0 do not return any of : unscaled g, scaled G, gcorr
Ci         :1 Return scaled G in s_site()%gc for each component
Ci         :2 Return unscaled g in s_site()%gcu for each component
Ci         :4 Return Lloyd correction in s_site()%gc
Ci         :  Any combination is allowed
Ci         :10s digit
Ci         :0 do not calculate vertex
Ci         :1 calculate vertex (sets lvertex=T)
Ci            For vertex see K. Belashchenko's notes dated 18 June 2012
Ci            and PRB 81, 064410
Ci         :100s digit
Ci         :  Not used here
Ci         :1000s digit
Ci         :0 Do not make one-site exchange parameters J0
Ci         :1 Make J0, contracting by L (one J0 for each component)
Ci         :2 Make J0, retaining L resolution
Ci         :10000s digit
Ci         :0 Do not use constraining fields
Ci         :1 Apply constraining fields from s_site(ib)%bxc
Ci  norb   :number of lower+intermediate blocks for this site
Ci  izp    :number of current z-point
Ci  wz     :weight for the current z-point
Ci  ncp    :number of CPA components
Ci  omg    :Omega(z) for the given site
Ci  domg   :Omega(z+dz) for the given site (if lsev=T)
Ci  pfi    :potential parameters for given site
Ci  dpf    :derivatives of pf for scaling g->F
Ci  ddpf   :ddpf for g->G
Ci  dz     :finite-difference delta z to calculate dOmega/dz
Ci         :(used only when lsev=T)
Ci  pmomg  :P(coherent)-Omega
Ci         :(used only when lvrtx=.true.)
Co Outputs
Co  gc     :conditionally-averaged onsite scaled G for each component
Co         :gc(:,s1,:,s2,icp) = scaled gf G for component icp, spin (s1,s2)
Co  gcorr  :Lloyd correction to be used for single-particle energy
Co  j0     :Exchange interaction
Co  sfvrtx :Vertex gamma (see notes, K. Belashchenko 18 June 2012)
Co         :sfvrtx(:,:,1,icp) = gamma(up,dn) for component icp
Co         :sfvrtx(:,:,2,icp) = gamma(dn,up) for component icp
Co         :not used unless lvrtx = .true.
Co  sfvrtxrel :Vertex gamma for nspc=2 (spin,p,spin,norb,norb); to be passed to pasajjnoncol for J
Cl Local variables
Cl  lG     :T => calculate G from g
Cl  lgu    :T => make g
Cl  lsev   :If true, calculate correction to sumev from dOmega/dz
Cl  lvrtx  :If true, calculate spin flip vertex
Cr Remarks
Cr   Let d Pi_m = (Pi+ - Pi-)_m
Cr     J0    = sum_L J0_L
Cr     J0_L  = (I0_L - J00_L)
Cr     I0_L  =-(4pi)^-1 Im int dz [dPi+_m * (gai+ - gai-)_LL]
Cr     J00_L = (4pi)^-1 Im int dz sum_L' [dPi_L' gai+_L'L dPi_L gai-_LL']
Cr   Need to clean up the passing of arguments already contained in s_site
Cu Updates
Cu   17 Apr 12 (Belashchenko) Added constraining fields
Cu   24 Dec 11 (Belashchenko) Complete rewrite
Cu   08 Dec 08 (P. Larson) First created
C ----------------------------------------------------------------------
      use gcpa
      use structures
      implicit none
C ... Passed parameters
      integer norb,nspc,lrel,ncp,izp,mode,lcgf
      logical lso,ferm
      double precision dz
      complex(8), dimension(norb,2,1+lrel/2,ncp) :: pfi,dpf,ddpf
      complex(8), dimension(norb,nspc,norb,2)  :: omg,domg
      complex(8), dimension(norb,2,norb,2,ncp)   :: gc,gcu
      complex(8) gcorr(norb,ncp,2),wz
      complex(8) pmomg(norb,nspc,norb,2),sfvrtx(norb,1,norb,2,ncp)
      complex(8) sfvrtxrel(4,4,4,norb,norb,ncp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site
C ... Dynamically allocated arrays
      integer, allocatable :: idxkm(:,:,:)
      complex(8), allocatable :: pii(:,:,:)
C ... Local parameters
      logical lsev,lvertex,lG,lgu,lJresL,lbxc
      integer icp,n1,n2,nl,id,jd,nll,n3,ll
      real(8) fpi,pth,pi,eula(3),th(ncp,nspc),bxc(3,ncp),thi(nspc),
     .  j0(ncp,1,norb),tau(ncp,norb)
      complex(8), dimension(norb,2,norb,2) :: wk,wk1,wkrel,dodz,dmat
      double complex wk11(4,norb,norb),wk22(4,4,norb,norb)
      double complex pfloc(norb,2,norb,2)
      complex(8) dp1u,dp1d,dp2u,dp2d,delPi(norb),cj,cj00,cxx,tau0
      integer master,procid,mpipid
C ... Pointers
C     real(8), pointer :: wrk(:)
      master = 0
      procid = mpipid(1)

      pi = 4*datan(1d0) ; fpi = 4*pi
      lG = mod(mod(mode,10),2) /= 0
      lgu = mod(mod(mode,10)/2,2) /= 0
      lsev = mod(mode,10) >= 4
      lvertex = mod(mode/10,10) == 1
      lJresL = mod(mode/1000,10) == 2
      lbxc = mod(mode/10000,10) == 1

      if (lrel == 2 .or. lso) then
!        if (lG) call rx('mkgcpa: lrel=2 or SOC and lG should not occur')
!        if (.not. lgu) call rx('mkgcpa: lrel=2 or SOC must request gcu')
        nl = sqrt(norb + 0.1)
      endif
      if (nspc == 2) then
!        if (lvertex) call rx('lvertex not ready for nspc=2')
!        if (lJresL) call rx('lvertex not ready for lJresL')
      endif

      th = s_site%thet
      if (lbxc) then
        call dcopy(3*ncp,s_site%bxc,1,bxc,1)
      else
        bxc = 0
      endif

C --- For each component, do ---
      do  icp = 1, ncp
        thi = th(icp,:)
        pth = datan(bxc(1,icp)) ; thi(1) = thi(1) + pth
        if (lso) then
          call gc00(0,norb,nspc,lrel,lso,thi,pfr=s_site%pfr(:,icp),
     .      omg=omg,g00=wk)
        else
          call gc00(0,norb,nspc,lrel,lso,thi,pfi=pfi(:,:,:,icp),
     .      omg=omg,g00=wk)
        endif

        if (lcgf == 27 .or. nspc == 2) then
          call zcopy(norb*norb*4,s_site%pfr(:,icp),1,pfloc,1) ! use pfr above?
          nll = ll(norb)+1
          call zcopy(size(pfloc),pfloc,1,wkrel,1)
          allocate(idxkm(2,0:nll-1,2*(nll+1)))
C         call zprm('P, kappa-mu',2,pfloc,norb*2,norb*2,norb*2)
          call mstokm(1,nll,1,norb,pfloc,wkrel,idxkm) ! check the mode (1 or 2)
C         call zprm('P, lms',2,pfloc,norb*2,norb*2,norb*2)
          deallocate(idxkm)
          allocate(pii(4,norb,norb))
          do id = 1, norb
          do jd = 1, norb
           pii(1,id,jd) =  (pfloc(id,1,jd,1) + pfloc(id,2,jd,2))/2d0                        ! 0
           pii(2,id,jd) =  (pfloc(id,1,jd,2) + pfloc(id,2,jd,1))/2d0 !(0d0,0d0)!        !x
           pii(3,id,jd) =  (-pfloc(id,1,jd,2) + pfloc(id,2,jd,1))*(0d0,1d0)/2d0                !y
           pii(4,id,jd) =  (pfloc(id,1,jd,1) - pfloc(id,2,jd,2))/2d0                        !z
          enddo
          enddo
        endif

        if (nspc /= 2) then
        delPi(:) = pfi(:,1,1,icp) - pfi(:,2,nspc,icp)
        else
          do id = 1, norb
           delPi(id) = pfloc(id,1,id,1) - pfloc(id,2,id,2)
          enddo
        endif

C   --- Calculate the exchange parameters (wk now holds g_ai)
        if (mod(mode/1000,10) /= 0) then
C         call zprm('delPi',2,delPi,norb,norb,1)
          call dcopy(ncp*1*norb,s_site%j0,1,j0,1)
!          call dcopy(ncp*norb,s_site%tau,1,tau,1)
          cj = 0; cj00 = 0
!          tau0 = 0
          do  n2 = 1, norb
            cj = cj + delPi(n2)*(wk(n2,1,n2,1)-wk(n2,2,n2,2))
            do  n1 = 1, norb
              cxx = delPi(n1)*wk(n1,1,n2,1)*delPi(n2)*wk(n2,2,n1,2)
              cj = cj + cxx ; cj00 = cj00 + cxx
!              if (ferm) then
!                print*,'j0 ferm'
!                tau0 = tau0 + dble(delPi(n1))*dimag(wk(n1,1,n2,1))*dble(delPi(n2))*dimag(wk(n2,2,n1,2))
!              endif
            enddo
            if (lJresL) then
              j0(icp,1,n2) = j0(icp,1,n2) - dimag(cj*wz)/fpi
C             j0(icp,2,n2) = j0(icp,2,n2) + dimag(cj00*wz)/fpi
              cj = 0; cj00 = 0
            endif
          enddo
!          if (ferm) tau(icp,1) = tau(icp,1) - dble(tau0*wz)/fpi
          if (.not. lJresL) then
            j0(icp,1,1) = j0(icp,1,1) - dimag(cj*wz)/fpi
C           j0(icp,2,1) = j0(icp,2,1) + dimag(cj00*wz)/fpi
          endif
          call dcopy(ncp*1*norb,j0,1,s_site%j0,1)
        endif

C   --- Vertex ---
C  For debugging: make gam(up,dn), gam(dn,up) using mc calculator
C  mc gii_up -i g1up -x delPi1 -v2dia -x gii_dn -i g1dn -x -x
C  mc gii_dn -i g1dn -x delPi1 -v2dia -x gii_up -i g1up -x -x
C  mc gii_up -i g2up -x delPi2 -v2dia -x gii_dn -i g2dn -x -x
C  mc gii_dn -i g2dn -x delPi2 -v2dia -x gii_up -i g2up -x -x
        if (lvertex) then
C         wk1(1:2,1) <- (Pi-Omega)*(Pa-Omega)^-1 for spins 1:2
C         For nspc = 2 the following needs to be rethought and redone

         if (lcgf /= 27) then
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      pmomg,norb,wk,2*norb,(0d0,0d0),wk1(1,1,1,1),2*norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      pmomg(1,1,1,2),norb,wk(1,2,1,2),2*norb,
     .      (0d0,0d0),wk1(1,2,1,1),2*norb)

C          call zprm('gii^-1 gai (up)',2,wk1(1,1,1,1),norb,norb,norb)
C          call zprm('gii^-1 gai (dn)',2,wk1(1,1,2,1),norb,norb,norb)

C         wk1(1:2,2) <- [(Pi-Omega) * (Pa-Omega)^-1] * delPi for spins 1,2
          do  n2 = 1, norb
            do  n1 = 1, norb
              wk1(n1,1,n2,2) = wk1(n1,1,n2,1) * delPi(n2)
              wk1(n1,2,n2,2) = wk1(n1,2,n2,1) * delPi(n2)
            enddo
          enddo

C          call zprm('... * delPi (up)',2,wk1(1,1,1,2),norb,norb,norb)
C          call zprm('... * delPi (dn)',2,wk1(1,1,2,2),norb,norb,norb)

C         wk1(1:2,1) <- (Pa-Omega)^-1*(Pi-Omega) for spins 1,2
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk,2*norb,pmomg,norb,(0d0,0d0),wk1(1,1,1,1),2*norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk(1,2,1,2),2*norb,pmomg(1,1,1,2),norb,
     .      (0d0,0d0),wk1(1,2,1,1),2*norb)

C         sfvrtx(1)<- wk1(1,2) * wk1(2,1); sfvrtx(2)<- wk1(2,2) * wk1(1,1)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk1(1,1,1,2),2*norb,wk1(1,2,1,1),2*norb,(0d0,0d0),
     .      sfvrtx(1,1,1,1,icp),norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk1(1,2,1,2),2*norb,wk1(1,1,1,1),2*norb,(0d0,0d0),
     .      sfvrtx(1,1,1,2,icp),norb)

C          print *, icp
C          call zprm('gam(up,dn)',2,sfvrtx(1,1,1,icp),norb,norb,norb)
C          call zprm('gam(dn,up)',2,sfvrtx(1,1,2,icp),norb,norb,norb)

         else ! lcgf = 27

C         wk11(1:4,1) <- (Pi-Omega)*(Pa-Omega)^-1 for spins upup,updown,downup,downdown

!          if (nspc /= 2) then
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      pmomg(:,1,:,1),norb,wk(:,1,:,1),norb,(0d0,0d0),wk11(1,:,:),norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      pmomg(:,1,:,2),norb,wk(:,1,:,2),norb,(0d0,0d0),wk11(2,:,:),norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      pmomg(:,2,:,1),norb,wk(:,2,:,1),norb,(0d0,0d0),wk11(3,:,:),norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      pmomg(:,2,:,2),norb,wk(:,2,:,2),norb,(0d0,0d0),wk11(4,:,:),norb)

C         wk1(1:4(p),4(spin)) <- [(Pi-Omega) * (Pa-Omega)^-1] * delPi(0,x,y,z) for all spins
          do  n2 = 1, norb
            do  n1 = 1, norb
            do  n3 = 1, 4
              wk22(n3,1,n1,n2) = wk11(n3,n1,n2) * pii(1,n1,n2)
              wk22(n3,2,n1,n2) = wk11(n3,n1,n2) * pii(2,n1,n2)
              wk22(n3,3,n1,n2) = wk11(n3,n1,n2) * pii(3,n1,n2)
              wk22(n3,4,n1,n2) = wk11(n3,n1,n2) * pii(4,n1,n2)
            enddo
            enddo
          enddo

C         wk11(1:4,1) <- (Pa-Omega)^-1*(Pi-Omega) for spins upup,updown,downup,downdown

          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk(:,1,:,1),norb,pmomg(:,1,:,1),norb,(0d0,0d0),wk11(1,:,:),norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk(:,1,:,2),norb,pmomg(:,1,:,2),norb,(0d0,0d0),wk11(2,:,:),norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk(:,2,:,1),norb,pmomg(:,2,:,1),norb,(0d0,0d0),wk11(3,:,:),norb)
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk(:,2,:,2),norb,pmomg(:,2,:,2),norb,(0d0,0d0),wk11(4,:,:),norb)

          do  n1 = 1, 4
          do  n2 = 1, 4
          do  n3 = 1, 4
          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
     .      wk22(n1,n2,:,:),norb,wk11(n3,:,:),norb,(0d0,0d0),
     .      sfvrtxrel(n1,n2,n3,:,:,icp),norb)
          enddo
          enddo
          enddo

C         sfvrtx(1)<- wk1(1,2) * wk1(2,1); sfvrtx(2)<- wk1(2,2) * wk1(1,1)
!          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
!     .      wk1(1,1,1,2),2*norb,wk1(1,2,1,1),2*norb,(0d0,0d0),
!     .      sfvrtx(1,1,1,1,icp),norb)
!          call zgemm('N','N',norb,norb,norb,(1d0,0d0),
!     .      wk1(1,2,1,2),2*norb,wk1(1,1,1,1),2*norb,(0d0,0d0),
!     .      sfvrtx(1,1,1,2,icp),norb)

         endif ! lcgf

        endif !vertex

        if (nspc == 2 .or. lcgf == 27) deallocate(pii)

!C ...        Lifetime using vertex
!        if (ferm .and. mod(mode/1000,10) /= 0 .and. nspc /= 2) then
!          call dcopy(ncp*norb,s_site%tau,1,tau,1)
!          tau0 = 0
!          do  n2 = 1, norb
!            do  n1 = 1, norb
!!              print*,'j0 ferm'
!              tau0 = tau0 + dble(sfvrtx(n1,1,n1,2,icp))*dimag(wk(n1,1,n2,1))*
!     .                dble(sfvrtx(n2,1,n2,1,icp))*dimag(wk(n2,2,n1,2))
!            enddo
!          enddo
!          tau(icp,1) = tau(icp,1) - dble(tau0*wz)/fpi
!          call dcopy(ncp*norb,tau,1,s_site%tau,1)
!          print*,'j0'
!        endif

C   --- Special branch for single-particle energy calculation ---
        if (lsev) then
C     ... Calculate U(dO/dz)U^-1
          dodz = 0
          if (nspc == 2) then
            call rx('mkgint: lsev does not work with nspc = 2')
            dodz(:,:,:,:) = (domg(:,:,:,:) - omg(:,:,:,:))/dz
            call rot_LandS(11,eula,nl,norb,1,dodz) ! this needs to be checked
          else
            dodz(:,1,:,1) = (domg(:,1,:,1) - omg(:,1,:,1))/dz
            dodz(:,2,:,2) = (domg(:,1,:,2) - omg(:,1,:,2))/dz
            call rotm(norb,-th(icp,1)-pth,dodz)
          endif
C     ... Multiply gc by it
          call zgemm('N','N',2*norb,2*norb,2*norb,(1d0,0d0),wk,2*norb,
     .      dodz,2*norb,(0d0,0d0),dmat,2*norb)
          do  n1 = 1, norb
            gcorr(n1,icp,1) = - dmat(n1,1,n1,1)
            gcorr(n1,icp,2) = - dmat(n1,2,n1,2)
          enddo
        endif

C   --- Make scaled G ---
        if (lG) then
        wk1 = wk !save copy

C   ... Convert g(icp) to G(icp) in the frame of Bxc
        do  n2 = 1, norb
          dp2u = dpf(n2,1,1,icp) ; dp2d = dpf(n2,2,1,icp)
          do  n1 = 1, norb
            dp1u = dpf(n1,1,1,icp) ; dp1d = dpf(n1,2,1,icp)
            wk(n1,1,n2,1) = wk(n1,1,n2,1) * sqrt(dp1u*dp2u)
            wk(n1,1,n2,2) = wk(n1,1,n2,2) * sqrt(dp1u*dp2d)
            wk(n1,2,n2,1) = wk(n1,2,n2,1) * sqrt(dp1d*dp2u)
            wk(n1,2,n2,2) = wk(n1,2,n2,2) * sqrt(dp1d*dp2d)
          enddo
          wk(n2,1,n2,1) = wk(n2,1,n2,1) + ddpf(n2,1,1,icp)
          wk(n2,2,n2,2) = wk(n2,2,n2,2) + ddpf(n2,2,1,icp)
        enddo

C   ... Rotate full G to the frame of M from that of Bxc
        if (pth /= 0d0) call rotm(norb,pth,wk)

C   ... Record in gc
        gc(:,:,:,:,icp) = wk
        wk = wk1 !restore
        endif

C   --- Make unscaled G ---
        if (lgu) then
C     ... Rotate g to the frame of M (with SOC it is done in gfg2gr)
          if (pth /= 0d0 .and. .not. lso .and. lrel /= 2)
     .      call rotm(norb,pth,wk)
          gcu(:,:,:,:,icp) = wk
        endif

      enddo
      end

