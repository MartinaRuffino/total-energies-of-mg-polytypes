      subroutine asajft(irep,s_ctrl,s_spec,s_site,s_ham,s_pot,s_lat,
     .  s_cpasitej,nsitej,nbcmpj,ib,icomp,icmp,nbcmpi,mxcomp,xycomp,vRLshf,
     .  offH,iprmb,nbas,nrc,pos,nkp,qp,iz,nzp,zp,wz,nk1,nk2,nk3,ipq,igstar,
     .  ifac,qb,JqRR,sumRL,sumR,qnu)
C- Magnetic exchange J(q) for a block of orbitals by the FT technique
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham lncol nl lgen3 lcgf ldlm nbas ldomg nbasp nspin
Ci                 nclass nspec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkcpa mkgint hamfbz mkpotf makpfz
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr z rmt lmxa hcr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkpotf makpfz
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  class ncomp norb domg omg omgn pnu spec clabel v0
Co     Stored:     *
Co     Allocated:  sfvrtx gcu
Cio    Elts passed: norb sfvrtx thet cpawt gcu domg bxc gc gcorr j0
Cio                ncomp dlmcl
Cio    Passed to:  mkcpa mkgint pasajq mkpotf makpfz pasajf
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham hord nlibu lgen3 lncol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb offH
Cio    Passed to:  mkcpa mkgint gfp0f1 pdglrs gfdpp gfg0g mkpotf makpfz
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  pf vmtz ves
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:vshft pp dpf ddpf pprel pf cp dddpf mxy qcorr sop
Cio    Passed to:  mkcpa mkgint gfp0f1 pdglrs gfdpp gfg0g pasajq mkpotf
Cio                makpfz
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:istab symgr ag
Cio    Passed to:  mkcpa
Cio  s_cpasitej
Ci     Elts read:  ib ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ncomp
Cio    Passed to:  pasajf
Ci Inputs
Ci   irep  :tells asajft about the nature of g written on disk.
Ci         :0 g is g=(P-S)^-1, AKA scattering path operator T
Ci         :    NB: beta-representation for g+ and g- should
Ci         :    be the same repsn, e.g. spin-averaged gamma rep
Ci         :1 g is proper Green's function matrix G
Ci         :  g will be subsequently scaled to g if s_ctrl%lcgf<20
Ci         :    and J is calculated in bare representation
Ci   nsitej:number of sites in source channel (column index)
Ci   nbcmpj:number of components associated with source channel
Ci         :Normal sites count as one component;
Ci         :CPA sites have multiple components.
Ci   ib    :site index
Ci   icomp0:CPA component within the site
Ci         :if 0, site is not a CPA site
Ci   icmp  :cumulative component for ib index, including all prior ones
Ci   nbcmpi:number of components associated with field channel
Ci         :Normal sites count as one component;
Ci         :CPA sites have multiple components.
Ci   mxcomp:(for dimensioning) max number of CPA components at any site
Ci   vRLshf:array of orbital- and spin-dependent potential shifts
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nbas  :size of basis
Ci   xycomp:3 for noncollinear for J_x,y,z(3*3), 1 otherwise
Ci   nrc   :nrc(i) = number of atoms in the ith class
Ci   pos   :basis vectors
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   qp    :irreducible k-points
Ci   iz    :energy index (used in CPA case to as index to Omega array)
Ci   nzp   :Last energy point, used to determine if at the Fermi level
Ci   zp    :complex energy
Ci   wz    :weights for complex energy integration
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   igstar:information needed for rotation qp to star of q
Ci         :bzmesh should be called with igstar(0)=-2
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   qb    :vectors of a microcell in the Brillouin zone
Co Outputs
Co   JqRR  : (lcgf = 10..14) exchange coupling between sites, r.s.
Co         : (lcgf = 20..26) susceptibility chi, q.s.
Co   sumRL :quantities needed for sum rules and on-site exchange;
Co         :see Remarks in pasajj below.
Co         :In the CPA case sumRL is resolved by:
Co         :CPA component icomp in the ib index
Co         :cumulative CPA component jccomp in the jb index.
Co   sumR  :L-contraction of sumRL, for current (cumulative) component
Co         :In the CPA case sumR is resolved by:
Co         :cumulative CPA component icmp in the ib index
Co         :cumulative CPA component jccomp in the jb index.
Co   qnu   :energy-weighted moments of the sphere charges
Cl Local variables
Cr Remarks
Cu Updates
Cu   18 Jun 18 Revise part of relativistic branch; synchronize with updated spherical harmonics
Cu   25 Jun 16 Bundle kp in gfg2gnc to speed energy scaling
Cu   25 Jun 16 (Vishina) extended to relativistic case, working version
Cu   05 Sep 15 Merge pvgfevcp with pvgfev
Cu   02 Sep 15 (Vishina) extended to relativistic case, first cut
Cu   05 Apr 13 Implement exact one-site rotation in CPA case,
Cu             Eq 37 in JMMM 67, 65 (1987)
Cu   17 Jan 13 First attempt at extending exchange calculations to CPA
Cu   18 Dec 12 Completed migration to f90 pointers
Cu   14 Nov 12 migrated f77 pointers to f90 pointers; setup for CPA
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Nov 10 use ALLOCATE to allocate large arrays
Cu   13 Jan 09 Longitudinal response function (lcgf=24)
Cu   16 Nov 07 Added vRLshf for LDA+U (diagonal-only)
Cu   21 Dec 01 Calls new mkpotf to make potential functions
Cu   05 Jun 01 New argument sumRL; sumR scaled to correspond to J0
Cu   15 Mar 00 extended to included i-waves in generation of gf
Cb Bugs    set up only for case nk1=k1, etc
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      double precision vRLshf(*)
      integer ib,irep,nbas,nkp,nk1,nk2,nk3,iz,icomp,icmp,nbcmpi,nsitej,nzp,xycomp
      integer mxcomp,nbcmpj
      integer ipq(1),ifac(3),igstar(0:*),iprmb(*),nrc(nbas)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,nbas+1)
C     Second spin in JqRR not used in transverse case
      double precision JqRR(nk1,nk2,nk3,xycomp**2,nbcmpi,nbcmpj,2)
      double precision pos(3,*),qb(3,3),zp(2),qnu(3,*),qp(3,nkp)
      double complex sumRL(11,xycomp,xycomp,0:mxcomp,0:mxcomp,*),
     .  sumR(11,xycomp,xycomp,0:mxcomp,0:mxcomp,*),wz
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_cpasite):: s_cpasitej(nsitej)
C ... Dynamically allocated arrays
      integer,allocatable:: ipc(:),ipdc(:)
      complex(8),allocatable:: gr(:,:,:,:,:,:,:),gc(:,:,:,:,:,:,:)
      complex(8),allocatable:: ghh(:,:),wk(:),gam2i(:,:),gam1i(:),mskm(:,:)
      complex(8), pointer :: pfbr(:),dpfbr(:),ddpfbr(:)
      complex(8), pointer :: pfloci(:,:,:,:)
C ... Local parameters
      logical lso,bittst
      integer hord,i,iopt,isp,k1,k2,k3,lcgf,ldh,ldham(16),ldim,ldimi,lrel,
     .  ldimx,ldpf,lgen3,lham,lidim,lihdim,lio,lncol,mxorb,nglob,nl,j,k,
     .  nlibu,nlma,nspc,offai,offc,offr,pfdim,jcomp,ibas,norbx,mode,offpi,offlmu
      integer, parameter :: nsp=2
      double precision plat(3,3),xx,vshft(nbas)
      double complex zsum
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (nspc,ldham(4)),(ldimx,ldham(5))
      equivalence (pfdim,ldham(9))
C ... For file I/O of gf
      integer clp(9,2),fopnx,ifi,iogfrs
      character*8 fnam
      integer MNC,MNB,MNT,MLB,MPL,MZP,MCD
      parameter (MNC=1,MNB=2,MNT=4,MLB=8,MPL=16,MCD=32,MZP=64)

C ... Setup
      call tcn('asajft')
      lham = s_ctrl%lham
      lrel = mod(s_ctrl%lrel,10)
      lncol = s_ctrl%lncol
      lso = bittst(lncol,4)
      nl = s_ctrl%nl
      lgen3 = s_ctrl%lgen3
      lcgf = s_ctrl%lcgf
      plat = s_lat%plat
      ldham = s_ham%ldham
      hord = s_ham%hord
      nlibu = s_ham%nlibu
      ldpf  = pfdim

C     Maximum norb among all sites: dimensions vertex
      norbx = 0
      do  i  = 1, nbas
        norbx = max(norbx,s_site(i)%norb)
      enddo

C      nclspd = s_ctrl%nclasp
C      if (ldlm /= 0) nclspd = s_ctrl%nclasp + s_ctrl%nccomp

!      call rxx(lncol/=0,'asajft not set up for noncollinear case')
      call fftz30(nk1,nk2,nk3,k1,k2,k3)
      if (allocated(ipc )) deallocate(ipc ); allocate(ipc(nbas))
      if (allocated(mskm)) deallocate(mskm); allocate(mskm(nl,2*nl))
      call sitepack(s_site,1,nbas,'class',1,ipc,xx)
C     vshft s_pot%vshft without the (-7:0) elements
      call dpscop(s_pot%vshft,vshft,nbas,9,1,1d0)
      if (lgen3 /= 0) then
        call rx('asajft not set up for 3rd generation')
      endif

C --- G for all qp in subblock connected to site ib ---
C ... Rewind file containing gf, read header
      call iinit(clp,9*2)
      clp(3,1) = lihdim*nspc
      clp(4,1) = lihdim*nspc
      clp(5,1) = lihdim*nspc
      mxorb = nglob('mxorb')
      fnam = 'gfqp'
C     do  isp = 1, 2
C      fnam = 'gfqp1'
C      if (isp == 2) fnam = 'gfqp2'
      ifi = fopnx(fnam,100,16+8+4+0,-1)
      rewind ifi
      call pshpr(1)
      lio = 2 + 10*(MNC*1+MNB*1+MNT*0+MLB*1+MPL*0+MCD*1+MZP*0)
      if (iogfrs(lio,0,0,' ',ifi,1,1,nbas,xx,zp,qp,plat,xx,clp,xx,xx,0,
     .  0,0,xx) /= 0) call rx('asajft: failed to read file header')
      call poppr
C     ifii(isp) = ifi
C     enddo

C --- Vertex for each CPA element in list ---
      if (s_ctrl%ldlm /= 0) then
        call mkcpa(s_ctrl,s_site,s_pot,s_ham,s_lat,.false.,iz)  ! Assemble s_pot%cp
C       Always use 100s digit mode=1 since s_site%omgn not available
        call mkgint(113,s_ctrl,s_site,s_pot,s_ham,iz,wz,0d0,iz==nzp)
      else                      ! Allocate trivial pointers to avoid compiler faults
        do  ibas = 1, nbas
          if (.not. associated(s_site(ibas)%sfvrtx))
     .      allocate(s_site(ibas)%sfvrtx(1,1))
        enddo
      endif

C ... Read g (irep=0) or G (irep=1) for subblock connected to site ib
      nlma = offH(4,1,ib+1) - offH(4,1,ib)
      offai = offH(4,1,ib)
      ldimi = offH(4,1,nbas+1)
      ldh = offH(3,1,nbas+1) - offH(3,1,ib)
      allocate(gr(nk1,nk2,nk3,nlma,nspc,ldimi,nsp))
      allocate(gc(nk1,nk2,nk3,ldimi,nspc,nlma,nsp))
      allocate(ghh(max(ldh,1),mxorb)); ghh = 0
      do  isp = 1, nsp
C       ifi = ifii(isp)
        offr = (isp-1)*nlma*ldimi ! spin-dependent initial offset to gr
        offc = (isp-1)*ldimi*nlma ! spin-dependent initial offset to gc
C        call gffbz(s_ctrl,nbas,offH,iprmb,ifi,ib,ib,440,pos,nkp,qp,zp,
C     .    nk1,nk2,nk3,nk1,nk2,nk3,ipq,plat,s_lat%istab,s_lat%symgr,
C     .    s_lat%ag,igstar,ifac,qb,offr,offc,ldh,gr,gc,ghh)
        if (lcgf == 20) then
          call rx('asajft not ready for lcgf = 20')
        endif
C        if (lcgf == 24 .and. isp == 2) cycle
C        if (lcgf == 25 .and. isp == 1) cycle
        call hamfbz(s_ctrl,s_ham,nbas,ifi,ib,ib,440,pos,nkp,qp,zp,
     .    nk1,nk2,nk3,nk1,nk2,nk3,ipq,plat,s_lat%istab,s_lat%symgr,
     .    s_lat%ag,igstar,ifac,qb,offr,offc,ldh,gr,gc,ghh)

C   ... Longitudinal susceptibility
        if (lcgf == 24) then
          call gfp0f1(s_ham,s_pot,0,31,.false.,offH,iprmb,ib,nbas,isp,2,
     .      nspc,wz,nk1,nk2,nk3,k1,k2,k3,ipq,igstar,nkp,nlma,offr,gr,gc,
     .      nbas,lidim,xx,xx,xx,xx,JqRR)
        endif

        if (nspc == 2) exit
      enddo

C ... Accumulate this energy's contribution to qnu for site ib
      if (lcgf < 20 .or. lcgf == 26 .or. lcgf == 27) then
        call pasajq(s_ctrl,s_pot,s_site,irep,nl,ib,lrel,nbas,iprmb,offH,
     .    nk1*nk2*nk3,zp,wz,s_pot%pp,vRLshf,s_pot%vshft,ipc,
     .    nrc,nspc,gr,nlma,ldimi,ghh,ldh,s_pot%dpf,s_pot%ddpf,pfdim,qnu)
      endif

C --- Exchange interactions ---
      if (lcgf < 20 .or. lcgf==26 .or. lcgf==27) then

C ... Bare potential functions P^0; prepare for scaling G->g
      if (irep == 1 .and. lcgf /= 26) then
C   ... Select order and kind of potential functions
        iopt = 0                  ! Second order P
        if (hord == 3) iopt = 1 ! Third order P
        if (hord == 1) iopt = 2 ! First order P
        if (hord == 4) iopt = 3 + 3400000 !  P computed from phi(z)
        iopt = iopt + 10000       ! P for bare representation (alpha=0)
C       iopt = iopt + 20000       ! P for gamma representation averaged over spins
CC       Make 0 repsn:       P     Pddot    sr(Pdot)
C        iopt = iopt + 10*(2**0  +  2**2  +  2**3) ! Make P, -1/2 P-dotdot/P-dot, sqrt(Pdot)

C       Make 0 repsn:       P     Pddot     Pdot
        iopt = iopt + 10*(2**0  +  2**2  +  2**1) ! Make P, -1/2 P-dotdot/P-dot, Pdot

        if (nlibu > 0) iopt = iopt + 10000000 ! Add LDA+U to shift C

        allocate(ipdc(nbas))
        call sitepack(s_site,1,nbas,'dlmcl',1,ipdc,xx)

C        call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,pfdim,ipc,ipdc,vRLshf,vshft,zp,
C     .    pfb,dpfb,ddpfb,xx,xx,xx,xx,xx)

C       Scalar relativistic P, Pdot, Pdotdot/Pdot
C        call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,pfdim,ipc,ipdc,
C     .    vRLshf,vshft,zp,pfb,dpfb,ddpfb,[xx],[xx],[xx],[xx],[xx])
        call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,pfdim,ipc,ipdc,
     .    vRLshf,vshft,zp,s_pot%pf,s_pot%dpf,s_pot%ddpf,[xx],[xx],[xx],[xx],[xx])
C       Relativistic P, sqrt(Pdot), Pdotdot/Pdot
        if (lrel == 2 .or. lso) then
          allocate(pfbr(ldpf*4))
          allocate(dpfbr(ldpf*4))
          allocate(ddpfbr(ldpf*4))
          call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,lrel,lso,pfdim,ipc,ipdc,
     .      vRLshf,vshft,zp,pfbr,dpfbr,ddpfbr,[xx],[xx],[xx],[xx],[xx])
          deallocate(pfbr,dpfbr,ddpfbr)
        endif

        deallocate(ipdc)
        if (s_ctrl%ldlm /= 0) then
!          call rx('asajft not ready for DLM  branch')
        endif

C ... g is given at the outset. Work in supplied representation
      else
C        pfb => s_pot%pf
C        allocate(dpfb(1,1),ddpfb(1,1)) ! Not used; anyway dpf should be sqrt(dpf)
      endif

      allocate(wk(2*k1*k2*k3))
      allocate(pfloci(norbx,2,norbx,2))
C     gamma in matrix form
      if (nspc == 2 .or. mxcomp /= 0) then
        allocate(gam2i(norbx,norbx))
        allocate(gam1i(1))
        mode = 12                                  ! Return gamma in matrix form from pfun in vector form
        if (lso .or. icomp /= 0) mode = mode-1   ! pfr is in matrix form
        if (icomp /= 0) mode = mode+200          ! Use 21 CPA vertex for gamma
        if (mod(s_ctrl%lrel,10) == 2) mode = mode+1000   ! convert to lms from (kappa,mu)
C     gamma in vector form
      else
        allocate(gam2i(1,1))
        allocate(gam1i(ldpf))
        mode = 2
C       call pvgfev(21,icomp,ib,offH,ldpf,s_pot%pf,s_site(ib)%norb,s_site(ib)%sfvrtx,gam1i)
      endif

C     Pup-Pdn, substituting vertex if icomp is a CPA site
      offpi = offH(1,1,ib)
      call pvgfevcp(s_site,mode,ib,icomp,s_pot%pf(1+offpi,1),ldpf,s_site(ib)%norb,norbx,
     .  s_site(ib)%sfvrtx,gam1i,gam2i,pfloci)

C      call getpfi(s_site,nl**2,2,ib,icomp,1,iprmb,ldim,lihdim,
C     .  ldpf,nl**2,s_pot%pf,pfi)

C ... Exchange coupling, by inverse of xi, or LW approx
      if (lcgf == 13) then
        if (icomp /= 0) call rx('pasajf2 not ready for CPA')
        call pasajf2(irep,nk1,nk2,nk3,k1,k2,k3,ib,nbas,offH,ldpf,s_pot%pf,
     .    s_pot%dpf,s_pot%ddpf,wz,wk,gr,gc,JqRR,sumRL)
      else
        call pasajf(s_ham,s_ctrl,s_site,s_cpasitej,s_pot,irep,nk1,nk2,nk3,k1,k2,k3,lcgf,nbas,ib,
     .    nspc,nlma,ldimi,norbx,icomp,icmp,nbcmpi,nsitej,mxcomp,nbcmpj,xycomp,
     .    gam2i,gam1i,pfloci,offH,iz==nzp,wz,wk,gr,gc,JqRR,sumRL)
      endif
      deallocate(wk,gam1i,gam2i)

      offlmu = (ib-1)*nl*nl*2
      do j = 1,xycomp
      do k = 1,xycomp
      do  jcomp = 0, mxcomp
        do  i = 1, 11
          sumR(i,j,k,icomp,jcomp,ib) = zsum(nlma,sumRL(i,j,k,icomp,jcomp,1+offai),11*(1+mxcomp)**2*xycomp*xycomp)
        enddo
      enddo
        enddo
      enddo
      endif

   99 continue
      deallocate(gr,gc,ghh,mskm)

C ... Print out sum rule quantities
C      do  jcomp = 0, mxcomp
C      call awrit4(' Sum rule: TpT = %,7;11D,%,7;11D; '//
C     .  '-(T+ - T-) = %,7;11D,%,7;11D',' ',90,6,
C     .    dble(sumR(1,icomp,jcomp,ib)), dimag(sumR(1,icomp,jcomp,ib)),
C     .   -dble(sumR(2,icomp,jcomp,ib)),-dimag(sumR(2,icomp,jcomp,ib)))
C      enddo

C      print *, 'exit asajft',
CC     .  sum(JqRR(1,1,1,:,:,1)),sum(JqRR(1,1,1,:,:,2))
C     .  (JqRR(3,1,1,1,2,1)),(JqRR(3,1,1,2,1,1))
C      stop

      call tcx('asajft')
      end

      subroutine pasajf(s_ham,s_ctrl,s_site,s_cpasitej,s_pot,irep,nk1,nk2,nk3,k1,k2,k3,
     .  lcgf,nbas,ib,nspc,nlma,ldimi,norbx,icomp,icmp,nbcmpi,nsitej,mxcomp,nbcmpj,xycomp,
     .  gam2i,gam1i,pfloci,offH,ferm,wz,wk,gr,gc,JRR,sumRL)
C- Kernel called by asajft to generate exchange interactions
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:norb sfvrtx
Cio    Passed to:  *
Cio  s_cpasitej
Ci     Elts read:  ib ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ncomp
Cio    Passed to:  *
Ci Inputs
Ci   ferm  :T, calculate lifetime (point is at fermi energy)
Ci   irep  :tells pasajf about the contents of gr,gc
Ci         :0 G is g=(P-S)^-1, AKA scattering path operator T
Ci         :  NB: epresentation for T+ and T- should be the same,
Ci         :  e.g. spin-averaged gamma rep
Ci         :  J is calculated in the given representation
Ci         :1 G is proper Green's function matrix
Ci         :  G will be scaled to T in the whatever repsn potential
Ci         :  functions pfg,dpfb,ddpfb are given.
Ci         :  J is calculated in that representation.
Ci         : Charge and pair exchange interactions should not be
Ci         : affected by representation, but components are.
Ci         : This alters, e.g. J_00 and pTpT.  Nevertheless the sum rule
Ci         : should continue to be satisfied.
Ci   nk1,nk2,nk3 k-divisions defining mesh for gr,gc
Ci   k1,k2,k3    dimensions for FFT, may be different from nk1..nk3
Ci   lcgf  :mode and opions for crystal Green's function code
Ci   nbas  :size of basis
Ci   ib    :calculated exchange coupling to sites ib
Ci   icomp :Current CPA component, or 0 if site not a CPA site
Ci   icmp  :Current combined site/CPA component corresponding to ib
Ci   nbcmpi:Number of combined site/CPA components (for dimensioning)
Ci   nsitej:Number of sites in site list
Ci   mxcomp:Maximum number of CPA components at any site (for dimensioning)
Ci   nbcmpj:>= 1+total number of CPA components in site list (for dimensioning)
Ci   gam2i :vertex for sitei, in matrix form
Ci   gam1i :Diagonal part of vertex, essentially (Pup-Pdn)
Ci         :CPA case: diagonal part of spin-21 vertex at site ib
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   ldpf  :dimensions potential functions P
Ci   pfb   :potential functions the repsn J is to be calculated in
Ci         :(see irep)
Ci   dpfb  :(used only when irep=1) sqrt-dot-pfb
Ci   ddpfb :(used only when irep=1) -1/2 pfb-dotdot/pfb-dot
Ci   pfbr   :potential functions (fully relativistic case)
Ci   wz    :weights for complex energy integration
Ci   wk    :work array
Ci   gr    :G_ij (see irep), i=field channel, j=source channel
Ci         :(?) which is field(?) compare dglrft, and gfree
Ci   gc    :G_ji
Co Outputs
Co   Jrr   :Contribution from one energy point to exchange J_RR';
Co          see pasajj below.  Jrr is resolved by site pairs.
Co         :In the CPA case :site => compound site/CPA component index
Co         :JRR(1,1,1,icmp,jcmp)
Co   sumRL :quantities needed for sum rules and on-site exchange;
Co         :see Remarks in pasajj below.
Co         :In the CPA case sumRL is resolved by:
Co         :1 CPA component icomp in the ib index
Co         :2 cumulative CPA component jccomp in the jb index.
Co         :  Component set to 0 for any site jb NOT a CPA site
Cl Local variables
Cl   gam2j :Crystal case: (Pup-Pdn)
Cl         :CPA case: diagonal part of spin-12 vertex at site jb
Cr Remarks
Cu Updates
Cu   02 Sep 15 (Vishina) extended to relativistic case
Cu   17 Jan 13 First attempt at extending exchange calculations to CPA
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical ferm
      integer irep,lcgf,nkap0,n0H,nspc
      parameter (nkap0=4,n0H=5)
      integer nk1,nk2,nk3,k1,k2,k3,nbas,nsitej,mxcomp,nbcmpj,norbx,xycomp,
     .  ib,nlma,ldimi,icomp,icmp,nbcmpi,offH(n0H,nkap0,nbas+1)
      double precision JRR(nk1,nk2,nk3,xycomp,xycomp,nbcmpi,*)
      integer, parameter :: nsp=2
      double complex wk(k1,k2,k3,2),wz
      double complex sumRL(11,xycomp,xycomp,0:mxcomp,0:mxcomp,*)
      double complex gr(nk1,nk2,nk3,nlma,nspc,ldimi,nsp)
      double complex gc(nk1,nk2,nk3,ldimi,nspc,nlma,nsp)
      double complex gam2i(norbx,norbx) ! Actually gam2i(ldpf,1) unless a vertex
      double complex pfloci(norbx,2,norbx,2)
      double complex gam1i(*)
C ... For structures
!      include 'structures.h'
      type(str_cpasite):: s_cpasitej(nsitej)
      type(str_ctrl)::  s_ctrl
      type(str_pot)::  s_pot
      type(str_site)::  s_site(*)
      type(str_ham)::   s_ham
C ... Dynamically allocated local arrays
      complex(8), pointer :: gam1j(:),gam2j(:,:),sumRLi(:,:,:,:)  !,pfi(:,:),mskm(:,:)
      complex(8), allocatable :: wkg(:,:,:,:)
      complex(8), pointer :: pflocj(:,:,:,:)
C ... Local parameters
      logical lso,bittst
      integer nkfbz,isp,nlmi,offpj,offpi,nlmj,ldim,offai,i,jb,jcomp,jcmp,jib,gamcomp,
     .   lhdim,nl,nglob,mode,mod1,ldpf,pfdim(2)
!      double complex pflocj(norbx,2,norbx,2)
      equivalence(ldpf,pfdim(1))
C     integer i1,i2,i3
      double precision fac

      call tcn('pasajf')
      nkfbz = nk1*nk2*nk3
      pfdim(1:2) = shape(s_pot%pf)  ! Recover ldpf from structure
C     nlma  = offH(4,1,ib+1)-offH(4,1,ib)
      ldim  = offH(1,1,nbas+1)
C     ldimi = offH(4,1,nbas+1)
      lhdim = ldimi + offH(3,1,nbas+1)
      nlmi  = offH(1,1,ib+1) - offH(1,1,ib) ! number of orbitals in lower block,  ib
      offpi = offH(1,1,ib) + 0*ldim
      nl = nglob('nl')
      lso = bittst(s_ctrl%lncol,4)

C ... Scale gr, gc to bare representation ... scale only lower block
      if (irep == 1 .and. lcgf /= 26) then
        allocate(wkg(nlma,nspc,ldimi,nspc))
C        call zprm('g345(1) before',2,gr(3,4,5,:,1,:,1:nsp),nlma,nlma,ldimi)

         do  isp = 1, 2
C        if (nspc == 2) then  ! This branch also works in the collinear case but it is slower
         if (nspc == 2 .or. .true.) then
          call sugfg2gnc(lso,mod(s_ctrl%lrel,10),1,2,s_ctrl%nccomp,.false.,mod1)
          call gfg2gnc(1000+mod1,nk1*nk2*nk3,s_pot,s_site,isp,nsp,nspc,offH,s_ham%iprmb,
     .      nlma,ldimi,ib,(/0,0/),gr(1,1,1,1,1,1,isp),nbas)
          call gfg2gnc(2000+mod1,nk1*nk2*nk3,s_pot,s_site,isp,nsp,nspc,offH,s_ham%iprmb,
     .      ldimi,nlma,ib,(/0,0/),gc(1,1,1,1,1,1,isp),nbas)
C          do  i3 = 1, nk3
C          do  i2 = 1, nk2
C          do  i1 = 1, nk1
C            call zcopy(nlma*nspc*ldimi*nspc,gr(i1,i2,i3,1,1,1,isp),nk1*nk2*nk3,wkg,1)
C           call gfg2gnc(1000+mod1,1,s_pot,s_site,isp,nsp,nspc,offH,s_ham%iprmb,nlma,ldimi,
C     .        ib,(/0,0/),wkg,nbas)
C            call zcopy(nlma*nspc*ldimi*nspc,wkg,1,gr(i1,i2,i3,1,1,1,isp),nk1*nk2*nk3)
C            call zcopy(nlma*nspc*ldimi*nspc,gc(i1,i2,i3,1,1,1,isp),nk1*nk2*nk3,wkg,1)
C           call gfg2gnc(2000+mod1,1,s_pot,s_site,isp,nsp,nspc,offH,s_ham%iprmb,ldimi,nlma,
C     .        ib,(/0,0/),wkg,nbas)
C            call zcopy(nlma*nspc*ldimi*nspc,wkg,1,gc(i1,i2,i3,1,1,1,isp),nk1*nk2*nk3)
C         enddo
C         enddo
C         enddo
        else
C         Scales only the lower block: rows 1:nlmi
          call pvgfe2(isp,nsp,nkfbz,offpi,0,nlmi,ldimi,ldpf,s_pot%dpf,s_pot%ddpf,nlma,ldimi,.true.,gr,gr)
          call pvgfe2(isp,nsp,nkfbz,0,offpi,ldimi,nlmi,ldpf,s_pot%dpf,s_pot%ddpf,ldimi,nlma,.true.,gc,gc)
        endif
        if (nspc == 2) exit  ! no second spin
        enddo
C        print *, 'cksum', sum(gr(1,1,1,:,:,:,:))
C        print *, 'cksum', sum(gr), sum(gc); stop

C       call zprm('g345(1) after',2,gr(3,4,5,:,1,:,1:nsp),nlma,nlma,ldimi)
        deallocate(wkg)
      endif

      offai = offH(4,1,ib)
      nlma  = offH(4,1,ib+1) - offH(4,1,ib)
      offpi = offH(1,1,ib)

      allocate(sumRLi(11,xycomp,xycomp,nlma))
      allocate(pflocj(norbx,2,norbx,2))

C     gamma in matrix form
      if (nspc == 2 .or. mxcomp /= 0) then
        allocate(gam2j(norbx,norbx))
        allocate(gam1j(1))
C     gamma in vector form
      else
        allocate(gam2j(1,1))
        allocate(gam1j(ldpf))
      endif

!     call getpfi(s_site,nl**2,2,ib,icomp,1,iprmb,ldim,lhdim,ldpf,nl**2,s_pot%pf,pfi) ! Not used

      jcmp = 0
      do  jib = 1, nsitej
      jb = s_cpasitej(jib)%ib
      do  jcomp = 0, s_cpasitej(jib)%ncomp
        if (jcomp == 0 .and. s_cpasitej(jib)%ncomp > 0) cycle
        if (jcomp == 0) then
          fac = 1
        else
          fac = s_site(jb)%cpawt(jcomp)
        endif
        jcmp = jcmp+1

C       gamma in matrix form
        if (nspc == 2 .or. mxcomp /= 0) then
          mode = 12                                  ! Return gamma in matrix form from pfun in vector form
          if (lso .or. jcomp /= 0) mode = mode-1   ! pfr is in matrix form
          if (jcomp /= 0) mode = mode+100          ! Use 21 CPA vertex for gamma
          if (mod(s_ctrl%lrel,10) == 2) mode = mode+1000   ! convert to lms from (kappa,mu)
C       gamma in vector form
        else
          mode = 2
        endif

        offpj = offH(1,1,jb)
        nlmj = offH(1,1,jb+1) - offH(1,1,jb)

C       Pup-Pdn, substituting vertex if jcomp is a CPA site
        call pvgfevcp(s_site,mode,jb,jcomp,s_pot%pf(1+offpj,1),ldpf,s_site(jb)%norb,norbx,
     .    s_site(jb)%sfvrtx,gam1j,gam2j,pflocj)

        call dpzero(sumRLi,2*11*nlma*xycomp*xycomp)
        if (lcgf == 10) then
          if (nspc == 2 .or. mxcomp /= 0) then
C           print*,'Entering pasajjcpa ib =',ib,'jb = ',jb
C            call pasajjcpa(nk1,nk2,nk3,nspc,mod(s_ctrl%lrel,10),s_site(jb)%norb,norbx,ldpf,k1,k2,k3,
C     .        nlmi,nlmj,0,ferm,wz,gam2i,gam2j,gam1i,gam1j,ib == jb,.true.,
C     .        gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),
C     .        s_pot%pprel(1,ib),s_pot%pprel(1,jb),
C     .        nlma,ldimi,nl,wk,JRR(1,1,1,icmp,jcmp),sumRLi)
            call pasajjcpa(nk1,nk2,nk3,nspc,0*mod(s_ctrl%lrel,10),norbx,ldpf,k1,k2,k3,
     .        nlmi,nlmj,0,ferm,wz,gam2i,gam2j,gam1i,gam1j,ib == jb,.true.,
     .        gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),
     .        nlma,ldimi,wk,JRR(1,1,1,1,1,icmp,jcmp),sumRLi)
          else
            call pasajj(nk1,nk2,nk3,k1,k2,k3,offpi,offpj,
     .        nlmi,nlmj,0,ferm,wz,gam1i,gam1j,ib == jb,.true.,
     .        gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),nlma,ldimi,wk,JRR(1,1,1,1,1,icmp,jcmp),sumRLi)
          endif
        elseif (lcgf == 27) then
          if (mod(s_ctrl%lrel,10) /= 2 .and. .not. lso) call rx('mode 27 is for rel=2 (relat. case)')
          if (nspc == 2 .or. mxcomp /= 0) then
C           print*,'Entering pasajjcpa ib =',ib,'jb = ',jb
C            call pasajjcpa(nk1,nk2,nk3,nspc,mod(s_ctrl%lrel,10),s_site(jb)%norb,norbx,ldpf,k1,k2,k3,
C     .        nlmi,nlmj,0,ferm,wz,gam2i,gam2j,gam1i,gam1j,ib == jb,.true.,
C     .        gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),
C     .        s_pot%pprel(1,ib),s_pot%pprel(1,jb),
C     .        nlma,ldimi,nl,wk,JRR(1,1,1,icmp,jcmp),sumRLi)
          if (s_ctrl%lrel == 2 .and. mxcomp == 0) then
            gamcomp = 1
            call pasajjnoncol(nk1,nk2,nk3,nspc,mod(s_ctrl%lrel,10),s_site(jb)%norb,norbx,ldpf,gamcomp,
     .        k1,k2,k3,nlmi,nlmj,0,ferm,wz,pfloci,pflocj,gam2i,gam2j,
     .        ib == jb,.true.,gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),s_site(ib)%vtxrel,
     .        s_site(jb)%vtxrel,JRR(1,1,1,1,1,icmp,jcmp),sumRLi)
          elseif (mxcomp /= 0) then
            gamcomp = norbx
            call pasajjnoncol(nk1,nk2,nk3,nspc,mod(s_ctrl%lrel,10),s_site(jb)%norb,norbx,ldpf,gamcomp,
     .        k1,k2,k3,nlmi,nlmj,0,ferm,wz,pfloci,pflocj,gam2i,gam2j,
     .        ib == jb,.true.,gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),s_site(ib)%vtxrel,
     .        s_site(jb)%vtxrel,JRR(1,1,1,1,1,icmp,jcmp),sumRLi)
          else
            call pasajjnoncol(nk1,nk2,nk3,nspc,mod(s_ctrl%lrel,10),s_site(jb)%norb,norbx,ldpf,gamcomp,
     .        k1,k2,k3,nlmi,nlmj,0,ferm,wz,s_site(ib)%pfr,s_site(jb)%pfr,gam2i,gam2j,
     .        ib == jb,.true.,gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),s_site(ib)%vtxrel,
     .         s_site(jb)%vtxrel,JRR(1,1,1,1,1,icmp,jcmp),sumRLi)
          endif
          else
            call pasajj(nk1,nk2,nk3,k1,k2,k3,offpi,offpj,
     .        nlmi,nlmj,0,ferm,wz,gam1i,gam1j,ib == jb,.true.,
     .        gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),nlma,ldimi,wk,JRR(1,1,1,1,1,icmp,jcmp),sumRLi)
          endif
        elseif (lcgf == 26) then
          call pasachi(nk1,nk2,nk3,nspc,k1,k2,k3,nlmi,nlmj,0,wz,ib == jb,
     .      gr(1,1,1,1,1,1+offpj,1),gc(1,1,1,1+offpj,1,1,1),nlma,ldimi,wk,
     .      JRR(1,1,1,1,1,icmp,jcmp),sumRLi)
        else
          if (icomp /= 0) call rx('asajft 511 not ready for CPA')
          call pasaj2(nk1,nk2,nk3,k1,k2,k3,offpi,offpj,nlmi,nlmj,
     .      mxcomp*(1+nbcmpj),wz,ldpf,s_pot%pf,ib == jb,.true.,
     .      gr(1,1,1,1+offpj*nlma,1,1,1),gc(1,1,1,1+offpj,1,1,1),nlma,ldimi,wk,
     .      JRR(1,1,1,1,1,icmp,jcmp),sumRL(1,1,1,icomp,jcomp,1+offai))
        endif

        do  i = 1, nlma
          sumRL(1,:,:,icomp,jcomp,i+offai) =
     .    sumRL(1,:,:,icomp,jcomp,i+offai) + fac*sumRLi(1,:,:,i)
C         Same for each jcomp ... scale by fac to make sum_jcomp = 1
          sumRL(2,:,:,icomp,0,i+offai) =
     .    sumRL(2,:,:,icomp,0,i+offai) +     fac*sumRLi(2,:,:,i)
          sumRL(3,:,:,icomp,jcomp,i+offai) =
     .    sumRL(3,:,:,icomp,jcomp,i+offai) + fac*sumRLi(3,:,:,i)
          sumRL(4,:,:,icomp,jcomp,i+offai) =
     .    sumRL(4,:,:,icomp,jcomp,i+offai) + fac*sumRLi(4,:,:,i)
          sumRL(5,:,:,icomp,jcomp,i+offai) =
     .    sumRL(5,:,:,icomp,jcomp,i+offai) + fac*sumRLi(5,:,:,i)
          sumRL(6,:,:,icomp,jcomp,i+offai) =
     .    sumRL(6,:,:,icomp,jcomp,i+offai) + fac*sumRLi(6,:,:,i)
          sumRL(7,:,:,icomp,jcomp,i+offai) =
     .    sumRL(7,:,:,icomp,jcomp,i+offai) + fac*sumRLi(7,:,:,i)
          sumRL(8,:,:,icomp,jcomp,i+offai) =
     .    sumRL(8,:,:,icomp,jcomp,i+offai) + fac*sumRLi(8,:,:,i)
          sumRL(9,:,:,icomp,jcomp,i+offai) =
     .    sumRL(9,:,:,icomp,jcomp,i+offai) + fac*sumRLi(9,:,:,i)
          sumRL(10,:,:,icomp,jcomp,i+offai) =
     .    sumRL(10,:,:,icomp,jcomp,i+offai) + fac*sumRLi(10,:,:,i)
          sumRL(11,:,:,icomp,jcomp,i+offai) =
     .    sumRL(11,:,:,icomp,jcomp,i+offai) + fac*sumRLi(11,:,:,i)
        enddo


      enddo
      enddo
      deallocate(gam1j,gam2j,sumRLi,pflocj)
      call tcx('pasajf')

      end

      subroutine pasajj(nk1,nk2,nk3,k1,k2,k3,offpi,offpj,idim,jdim,
     .  mxcomp,ferm,wz,gammai,gammaj,ldiag,lchk,tij,tji,ldi,ldj,wk,Jrr,sumRL)
C- ASA Magnetic exchange J(q) by FFT for a pair of sites or CPA components
C ----------------------------------------------------------------------
Ci Inputs
Ci   nk1,nk2,nk3:mesh of k-points
Ci   offpi,offpj:offsets to potential functions
Ci   idim,jdim  :size of G(i,j) subblock to calculate
Ci   ferm       :1 if fermi energy
Ci   wz         :energy integration weight
Ci   gammaj      :Crystal case: (Pup-Pdn)
Ci              :CPA case: spin12 vertex.
Ci   gammai      :Crystal case: (Pup-Pdn)
Ci              :CPA case: spin21 vertex.
Ci   tij,tji    :T-:matrix connecting orbitals ij and ji, aka gij,gji
Ci              :tij,tji are supplied on a uniform mesh of q points
Ci   ldiag      :ij is diagonal block, in which case tji is not
Ci               explicitly available, but must be computed from tij
Ci   lchk       :T if to generate sumRL(2,:,:,:) for sum rule
Ci   ldi,ldj    :dimensions of tij,tji:
Ci               tij = tij(nk1,nk2,nk3,ldi,ldj,nsp)
Ci               tji = tji(nk1,nk2,nk3,ldj,ldi,nsp)
Ci   wk         :d.c. work array dimensioned at least nk1*nk2*nk2*2
Ci
Co  Outputs
Co    Jrr       :Contribution from one energy point to integral for
Co               exchange interaction is accumulated connecting site R+T
Co               (pot fun offset offpi) to site R' (offset offpj)
Co               for all translation vectors T.
Co               Jrr is a contraction over L and L' of the orbital-resolved
Co               J_R+TL,RL' . Jrr is generated and returned in real space.
Co
Co    sumRL     :quantities needed for sum rules and on-site exchange.
Co              :interactions.  Real and imaginary parts retained.
Co              :In contrast to Jrr, sumRL is resolved by orbital L
Co              :sumRL(1,L=1..idim) = sum_TL' J_R+TL,RL', RHS of Eq(10) below
Co              :sumRL(2,id=1..idim) = 'one-site' rotation:
Co              :                       -1*LHS of Eq.(10) below
Co              :sumRL(3,id=1..idim) = onsite J_RL,RL, AKA J00
Cr  Remarks
Cr    These Remarks refer to the following papers:
Cr      (A) JMMM 67, 65 (1987), a paper on exchange interactions in CPA
Cr      (B) JMMM 200, 148 (1999) a review paper on exchange interactions
Cr          In that paper T = (P-S)^-1 is equivalent to g in this code.
Cr
Cr    This routine calculates the exchange interactions J_R+T,R' between
Cr    R+T and R'.  J is calculated in the RPA from the inverse Bloch sum
Cr    of scattering path operators g+_RL,R'L'(q) and g-_RL,R'L'(q) (AKA
Cr    T+_RL,R'L'(q) and T-_RL,R'L'(q)), as described below.  Thus the
Cr    pairwise exchange interactions are returned in real space, for
Cr    pairs (R+T,R') where (R,R') are in the unit cell and T is a
Cr    lattice translation vector.  J is calculated for a given R and R'
Cr    but for all T, using fast Fourier transform techniques.
Cr    J_R+T,R' is a contraction of the orbital-resolved exchange
Cr    J_R+TL,R'L' between orbitals L and L' at sites R+T and R'.
Cr
Cr    * In the CPA case, sites R and R' are generalized to components.
Cr      At either R or R' more than one component may occupy a single
Cr      site, so the pair exchanges are calculated for a compound
Cr       (R + component at R; R' + component at R')
Cr      g+_RL,RL' and g-_RL,RL' (AKA T+ and T-) are independent of component.
Cr      The component-dependence appears in vertices gammaj and gammai.
Cr      Note: gammaj and gammai are properly matrices of tensor rank 2
Cr      This routine approximates these arrays by the diagonal part.
Cr
Cr    * The orbital-specific character of the J is summed over
Cr      to obtain a single interaction connecting sites R+T to R'
Cr      Thus J_R+T,R' = sum_LL' J_R+TL,R'L' contracting over orbitals
Cr      L and L' at R and R', respectively.
Cr
Cr    * sumRL contracts J_R+TL,R'L' over T and also L' but retains
Cr      resolution by L at R.  sumRL then contains a net effective
Cr      field acting on channel RL from all sites R' as follows:
Cr      1. sumRL(1) sum_TL' J_R+TL,R'L'
Cr      sumRL(1) is the total effective field acting RL from all R'+T.
Cr      2. -sumRL(2) is the "total" field acting on RL.
Cr      This "one-site" rotation is called I_i in Table I of paper (B)
Cr      (caution: slightly different definition in the text! see Eq(38))
Cr      and is generated according Eq(10) below.
Cr      In the crystalline case sumRL(1) should match -sumRL(2).
Cr      In the CPA case, there is no exact sum rule.
Cr      3. sumRL(3) contains the onsite exchange sum_L' J_RL,RL', Eq.(7)
Cr      It is needed to remove the spurious self-interaction from two
Cr      orbitals rotating at the same site against each other.
Cr
Cr    The exchange interaction between orbital pairs L at R, L' at R'
Cr    is calculated as follows (R,R' confined to cell)
Cr      X_R+TL;R'L' = g+_R+TL;R'L' g-_R'L',R+TL                        (1)
Cr    X is a pairwise susceptibility-like object
Cr    Pairwise exchange:
Cr      pi J_R+TL;R'L' = int^E_f dz Im [ dP_RL X_R+TL;R'L' dP_R'L']    (2)
Cr    This equation is
Cr    In the crystalline case,
Cr      dP_RL = (P+_RL - P-_RL)/2                                      (3)
Cr    i.e. spin splitting of potential functions in channel RL.  In the
Cr    CPA case, dP_RL is approximated as the diagonal component of the
Cr    vertex (the correct CPA expression should use the matrix form)!
Cr
Cr    Eq(2) is equivalent to the convolution
Cr    pi J_RL,R'L'(q) = int^E_f dz int dk Im dP_L g+RL,R'L'(k) g-_R'L',RL(k+q)
Cr    Instead of direct convolution g+(k) * g-(k+q), g+ and g- are shifted
Cr    to R.S. by FFT.  Thus Eq(2) is used:  J is generated in real space.
Cr
Cr    The net pairwise interaction is contracted over individual orbitals:
Cr      J_R+T,R' = sum_LL' J_R+TL,R'L'                                 (4)
Cr    J_R+T,R' is what is returned in Jrr.
Cr    The cumulative "mean field" at R generated by R' is
Cr      I_R = sum_T J_R+T,R'                                           (5)
Cr    Note also that I_R can be written in k-space:
Cr      I_R = sum_T J_R+T,R' exp(iq.T), for q=0                        (6)
Cr    I_R includes a spurious on-site rotation of a spin against itself
Cr      J00 = J_R,R  (i.e. T=0 and R=R')                               (7)
Cr    Thus the "mean exchange field" at R generated by R' is
Cr      J_R = I_R - J00  (see pvgfed)                                  (8)
Cr    I_R can be computed another way by a sum rule, which is:
Cr      (g+_RL,RL - g-_RL,RL)/2 = - sum_TL' X_R+TL;R'L' dP_R'L'        (9)
Cr    which we scale by dP_L and integrate to extract J_R, aka J0
Cr      -1/pi Im int dz dP_RL (g+_RL - g-_RL) = sum_TL' J_R+TL,RL'    (10)
Cr    See JMMM 200, 148 (1999) Eq. 40 for a similar relation.
Cr    Note: in the CPA case, this sum rule is not satisfied.
Cr    In that case the LHS of Eq(10) is still exact; the RHS is not.
Cr    See JMMM 67, 65 (1987), Eq. 36.
Cr    For the LHS the true dP_RL is used, not gamma.
Cr
Cr    The following are accumulated into sumRL (parameters L-resolved)
Cr    sumRL(1) : the RHS of Eq(10), equivalent to -sumRL(2) for crystals
Cr    sumRL(2) : (ldiag=T only) (-1) * LHS of Eq(10) ('one-site' exchange)
Cr    sumRL(3) : (ldiag=T only) an L-resolved equivalent of J_00 Eq(7)
Cu Updates
Cu   17 Jan 13 Modifications for CPA: replace dP with gammai and gammaj
Cu   Aug 15 For CPA, rel and non-coll cases gammai, gammaj are now in the matrix form
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical ldiag,lchk,ferm
      integer offpi,offpj,idim,jdim,ldi,ldj,mxcomp,nk1,nk2,nk3,k1,k2,k3,nsp
      parameter (nsp=2)
      double precision Jrr(nk1,nk2,nk3)
      double complex wz,gammaj(*),gammai(ldj)
      double complex sumRL(11,0:mxcomp,0:mxcomp,idim)
      double complex tij(nk1,nk2,nk3,ldi,ldj,nsp),wk(k1,k2,k3,2)
      double complex tji(nk1,nk2,nk3,ldj,ldi,nsp)
C ... Local parameters
      integer id,jd,i1,i2,i3,nfbz
      complex(8), allocatable :: Imtij(:,:,:,:,:),Imtji(:,:,:,:,:)
      double complex wkso(k1,k2,k3,2)
      double complex DelPi,DelPj,wpi,wpi2,tpt,cx,tpti,tptiso,tptso
      double precision cImtij,cImtji

      call tcn('pasajj')
C ... Weight for integration 1/pi de ...
      wpi = wz/(4*datan(1d0))
      wpi2 = wpi
C ... Use the following if add symmetrized contribution to J
C     wpi2 = wpi/2
      nfbz = nk1*nk2*nk3

!      call zprm('gammai',2,gammai(1+offpi*0),idim,idim,1)
!      call zprm('gammaj',2,gammaj(1+offpj*0),jdim,jdim,1)
!      call zprm('tij',2,tij(2,3,4,1:idim,1:jdim,1),idim,idim,jdim)
!      call zprm('tji',2,tji(2,3,4,1:jdim,1:idim,2),jdim,jdim,idim)

!      return

C ... lifetime
      if (ferm .and. nsp == 2) then
       allocate(Imtij(nk1,nk2,nk3,idim,jdim))
       allocate(Imtji(nk1,nk2,nk3,jdim,idim))
       do  i3 = 1, nk3
       do  i2 = 1, nk2
       do  i1 = 1, nk1
        do id = 1,idim
        do jd = 1,jdim
         cImtij = dimag(tij(i1,i2,i3,id,jd,1))
         if (ldiag) then
           cImtji = dimag(tij(i1,i2,i3,jd,id,2))
         else
           cImtji = dimag(tji(i1,i2,i3,jd,id,2))
         endif
         Imtij(i1,i2,i3,id,jd) = cImtij*(1d0,0d0)
         Imtji(i1,i2,i3,jd,id) = cImtji*(1d0,0d0)
        enddo
        enddo
       enddo
       enddo
       enddo
      endif

C --- Loop over orbital pairs ---
      do  id = 1, idim ! Orbitals associated with this site
      tpt = 0
      if (ferm .and. nsp == 2) tptso = 0
      do  jd = 1, jdim

C   ... dP = (P^+ - P^-)/2
C       DelPi = (pfun(id+offpi,1) - pfun(id+offpi,2))/2
C       DelPj = (pfun(jd+offpj,1) - pfun(jd+offpj,2))/2
        DelPi = gammai(id+offpi*0)/2
        DelPj = gammaj(jd+offpj*0)/2

C   ... w1 = T^+, w2 = dPj T^-
        call dcopy(nfbz*2,tij(1,1,1,id,jd,1),1,wk,1)
        if (ldiag) then
          call zcopy(nfbz,tij(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
        else
          call zcopy(nfbz,tji(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
        endif
        call zscal(nfbz,delPj,wk(1,1,1,2),1)

C   ... FT T_ij(q) and p_j T_ji(q); overwrite w1 with product of result (=> R.S.)
        call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ... lifetime
        if (ferm .and. nsp == 2) then
         call dcopy(nfbz*2,Imtij(1,1,1,id,jd),1,wkso,1)
         if (ldiag) then
           call zcopy(nfbz,Imtji(1,1,1,jd,id),1,wkso(1,1,1,2),1)
         else
           call zcopy(nfbz,Imtji(1,1,1,jd,id),1,wkso(1,1,1,2),1)
         endif
C        call dscal(2*nfbz,dble(delPi),wkso(1,1,1,1),1)
         call dscal(2*nfbz,dble(delPj),wkso(1,1,1,2),1)
         call fftz3c(wkso,wkso(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
        endif

C       call zprm3(' wk %i',id+100+jd,wk,nk1,nk2,nk3)

C   ... Accumulate Jrr = J_R+T,R' = sum_LL' J_R+TL,R'L', Eq(4) in the notes
C       and sumRL(1) = dP_RL sum_TL' X_R+TL;R'L' dP_R'L', RHS of Eq(10)
C       and sumRL(3) = L-resolved version of J_00, Eq(7)
        tpti = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
        if (ferm .and. nsp == 2) tptiso = 0
        do  i3 = 1, nk3
        do  i2 = 1, nk2
        do  i1 = 1, nk1
          tpti = tpti + wk(i1,i2,i3,1)
          if (ferm .and. nsp == 2) tptiso = tptiso + wkso(i1,i2,i3,1)
          Jrr(i1,i2,i3) = Jrr(i1,i2,i3)+dimag(wpi2*delPi*wk(i1,i2,i3,1))
        enddo
        enddo
        enddo
        if (ldiag) sumRL(3,0,0,id) = sumRL(3,0,0,id) + wpi2*delPi*wk(1,1,1,1)
        tpt = tpt + delPi*tpti ! tpt = dP_RL sum_L' tpti_RL;R'L'
        if (ferm .and. nsp == 2) tptso = tptso + dble(delPi)*tptiso
C       tpt = tpt + 1*tpti

CC   ... This is the symmetrized, presumably equivalent contribution
CC   ... w1 = (Pi^+ - Pi^-)/2 T^-, w2 = (Pj^+ - Pj^-)/2 T^+
C        if (ldiag) then
C          do  40  i3 = 1, nk3
C          do  40  i2 = 1, nk2
C          do  40  i1 = 1, nk1
C            wk(i1,i2,i3,1) = tij(i1,i2,i3,id,jd,2)
C            wk(i1,i2,i3,2) = delPj*tij(i1,i2,i3,jd,id,1)
C   40     continue
C        else
C          do  42  i3 = 1, nk3
C          do  42  i2 = 1, nk2
C          do  42  i1 = 1, nk1
C            wk(i1,i2,i3,1) = tij(i1,i2,i3,id,jd,2)
C            wk(i1,i2,i3,2) = delPj*tji(i1,i2,i3,jd,id,1)
C   42     continue
C        endif
CC   ... Overwrite p_i T_ij with convolution of p_i T_ij and p_j T_ji
C        call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,nk1,nk2,nk3,30,-1)
C
CC   ... Accumulate 1/2 p_i T-_ij pj T+_ji into Jrr
C        do  44  i3 = 1, nk3
C        do  44  i2 = 1, nk2
C        do  44  i1 = 1, nk1
C        tpt = tpt + wk(i1,i2,i3,1)
C   44   Jrr(i1,i2,i3) = Jrr(i1,i2,i3) +dimag(wpi2*delPi*wk(i1,i2,i3,1))
C        tpt = tpt + delPi*tpti

      enddo
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
      sumRL(1,0,0,id) = sumRL(1,0,0,id) + wpi2*tpt
      if (ferm .and. nsp == 2) sumRL(7,0,0,id) = sumRL(7,0,0,id) +tptso*(4*datan(1d0))!/nfbz
      enddo
      if (ferm .and. nsp == 2) deallocate(Imtij,Imtji)

C     print *,'516 after id,jd loop', id,sumrl(1,0,0,1)

C --- Debugging ... See Eqns 4 and 5 above (ldiag only) ---
      if (.not. ldiag .or. .not. lchk) goto 999
C     tupdn = 0
      do  110  id = 1, idim
        jd = id

C   ... dP = (P^+ - P^-)/2
C       DelPi = (pfun(id+offpi,1) - pfun(id+offpi,2))/2
        DelPi = gammai(id+offpi*0)/2  ! For the closest match to sum rule
C       DelPi = (pfi(id,1) - pfi(id,2))/2  ! For exact J0 (JMMM 67, 65 (1987))
C        print *, id,cmplx((pfi(id,1) - pfi(id,2))/2),
C     .    cmplx(gammai(id+offpi*0)/2)

C   ... dP_i sum_q (T^+_RR(i,j;q) - T^-_RR(i,j;q))
        cx = 0
        do  130  i3 = 1, nk3
        do  130  i2 = 1, nk2
        do  130  i1 = 1, nk1
          cx = cx + tij(i1,i2,i3,id,jd,1) - tij(i1,i2,i3,id,jd,2)
  130   continue
        cx = cx/nfbz
C       print *, id, cmplx(cx), cmplx(delPi)

        sumRL(2,0,0,id) = sumRL(2,0,0,id) + wpi*delPi*cx/2
C       sumRL(4,0,0,id) = sumRL(4,0,0,id) + wpi*delPig*cx/2
  110 continue

  999 call tcx('pasajj')
      end

      subroutine pasajjcpa(nk1,nk2,nk3,nspc,lrel,norbx,ldpf,k1,k2,k3,
     .  idim,jdim,mxcomp,ferm,wz,gam2i,gam2j,gam1i,
     .  gam1j,ldiag,lchk,tij,tji,ldi,ldj,wk,Jrrcp,sumRLcp)
C- ASA Magnetic exchange J(q) by FFT for a pair of sites or CPA components
C ----------------------------------------------------------------------
Ci Inputs
Ci   nk1,nk2,nk3:mesh of k-points
Ci   nspc       : 2 for non-collinear case, otherwise 1
Ci   offpj:offsets to potential functions
Ci   idim,jdim  :size of G(i,j) subblock to calculate
Ci   ferm       :1 if fermi energy
Ci   wz         :energy integration weight
Ci   gam2j     :Crystal case: (Pup-Pdn)
Ci              :CPA case: diagonal part of (spin12) vertex.
Ci   gam2i     :Crystal case: (Pup-Pdn)
Ci              :CPA case: diagonal part of (spin21) vertex.
Ci   tij,tji    :T-:matrix connecting orbitals ij and ji, aka gij,gji
Ci              :tij,tji are supplied on a uniform mesh of q points
Ci   ldiag      :ij is diagonal block, in which case tji is not
Ci               explicitly available, but must be computed from tij
Ci   lchk       :T if to generate sumRL(2,:,:,:) for sum rule
Ci   ldi,ldj    :dimensions of tij,tji:
Ci               tij = tij(nk1,nk2,nk3,ldi,nspc,ldj,2)
Ci               tji = tji(nk1,nk2,nk3,ldj,nspc,ldi,2)
Ci   wk         :d.c. work array dimensioned at least nk1*nk2*nk2*2*nspc
Ci   qnur       :relativistic ASA energy moments
Co  Outputs
Co    Jrr       :Contribution from one energy point to integral for
Co               exchange interaction is accumulated connecting site R+T
Co               (pot fun offset offpi) to site R' (offset offpj)
Co               for all translation vectors T.
Co               Jrr is a contraction over L and L' of the orbital-resolved
Co               J_R+TL,RL' . Jrr is generated and returned in real space.
Co
Co    sumRL     :quantities needed for sum rules and on-site exchange.
Co              :interactions.  Real and imaginary parts retained.
Co              :In contrast to Jrr, sumRL is resolved by orbital L
Co              :sumRL(1,L=1..idim) = sum_TL' J_R+TL,RL', RHS of Eq(10) below
Co              :sumRL(2,id=1..idim) = 'one-site' rotation:
Co              :                       -1*LHS of Eq.(10) below
Co              :sumRL(3,id=1..idim) = onsite J_RL,RL, AKA J00
Co              :sumRL(4,id=1..idim) orb. exch.
Co              :sumRL(5,id=1..idim) orb-spin exch.
Co              :sumRL(6,id=1..idim) spin-orb exch.
Co              :sumRL(7,id=1..idim) reserved for the lifetime
Cr  Remarks
Cr    This routine is an extension to the ASA exchange maker pasajj where the
Cr    vertex plays the role of P+ - P- .  In the crystalline case without SO
Cr    coupling, P+ - P-  is diagonal, whereas P+ - P- (or the CPA vertex,
Cr    the analog of P+ - P-) is only site diagonal.
Cr    For lifetime (still in work):
Cr     Antropov Physica B 237-238 336 (1997):
Cr     lifetime Im J_ij = -pi*w*Tr_L [p*Im(T_up)*p*Im(T_down)]
Cr
Cr    This routine is closely patterned after pasajj, which see.
Cu Updates
Cu   02 Sep 15 (Vishina) extended to relativistic case
Cu   10 Aug 15 (Vishina) adapted from pasajj
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical ldiag,lchk,ferm
      integer idim,jdim,ldi,ldj,mxcomp,lrel,ldpf
      integer nk1,nk2,nk3,k1,k2,k3,nsp,nspc,norbx
      parameter (nsp=2)
      double precision Jrrcp(nk1,nk2,nk3)
      double complex wz,gam2j(norbx,norbx),gam2i(norbx,norbx),gam2jreal(norbx,norbx)
     .          ,gam2ireal(norbx,norbx),gam2iim(norbx,norbx),gam2jim(norbx,norbx)
!      double complex gam1j(1,ldpf),gam1i(ldpf,1) ! for gamma 64*1
      double complex gam1j(ldpf),gam1i(ldpf)
!      double complex gam2jcp(norb,norb),gam2icp(norb,norb)
      double complex sumRLcp(11,0:mxcomp,0:mxcomp,idim)
      double complex tij(nk1,nk2,nk3,ldi,nspc,ldj,2),wk(k1,k2,k3,2),wkso(k1,k2,k3,2),wksoi(k1,k2,k3,2)
      double complex tji(nk1,nk2,nk3,ldj,nspc,ldi,2)
C     double complex ppreli(5,nl,2*nl,2,2),pprelj(5,nl,2*nl,2,2)
C ... Dynamically allocated local arrays
C     integer, allocatable    :: idx(:,:,:)
      complex(8), allocatable :: tijPi(:,:,:,:,:),tjiPj(:,:,:,:,:)
      complex(8), allocatable :: tijPir(:,:,:,:,:,:),tjiPjr(:,:,:,:,:,:)
C     complex(8), allocatable :: tijkm(:,:,:,:,:),tjikm(:,:,:,:,:)
C     complex(8), allocatable :: wkkm(:,:)!,wkijr(:,:,:,:),wkjir(:,:,:,:)
      complex(8), allocatable :: Imtij(:,:,:,:,:),Imtji(:,:,:,:,:)
      complex(8), allocatable :: Imtijr(:,:,:,:,:),Imtjir(:,:,:,:,:)
      complex(8), allocatable :: Imtiji(:,:,:,:,:),Imtjii(:,:,:,:,:)
C ... Local parameters
      integer id,jd,i1,i2,i3,nfbz,l,ll,nli,nlj
      double precision JrrR(nk1,nk2,nk3),Jso(nk1,nk2,nk3),Jorb(nk1,nk2,nk3),
     .                  mi,mj,cImtij,cImtji
      double complex wkji(jdim,idim),wkji1(jdim,idim)
      double complex wkij(idim,jdim),wkij1(idim,jdim)
      double complex DelPi,DelPj,wpi,wpi2,cx,tptcp,tpticp,tptir,tptr,cxr,
     .           tptiso,tptso
      double complex sumRLr(11,0:mxcomp,0:mxcomp,idim)
      call tcn('pasajjcpa')

C     call zprm('gam2j',2,gam2j,norb,jdim,jdim)

C ... Weight for integration 1/pi de ...
      wpi = wz/(4*datan(1d0))
      wpi2 = wpi
C ... Use the following if add symmetrized contribution to J
C     wpi2 = wpi/2
      nfbz = nk1*nk2*nk3
!      JrrR = 0
!      print*,'Jrrcp(i1,i2,i3) entry',Jrrcp(1,1,1)
      if (lrel == 2) then

        JrrR(:,:,:) = Jrrcp(:,:,:)
        sumRLr(:,:,:,:) = sumRLcp(:,:,:,:)

C ...... Spin-spin exchange

        allocate(tijPir(nk1,nk2,nk3,idim,jdim,2))
        allocate(tjiPjr(nk1,nk2,nk3,jdim,idim,2))

C   ... Lifetime (in progress)................................................
C   ... Antropov Physica B 237-238 336 (1997)
C   ... lifetime Im J_ij = -pi*w*Tr_L [p*Im(T_up)*p*Im(T_down)]
        if (ferm) then
         allocate(Imtij(nk1,nk2,nk3,idim,jdim))
         allocate(Imtji(nk1,nk2,nk3,jdim,idim))
         allocate(Imtijr(nk1,nk2,nk3,idim,jdim))
         allocate(Imtjir(nk1,nk2,nk3,jdim,idim))
         allocate(Imtiji(nk1,nk2,nk3,idim,jdim))
         allocate(Imtjii(nk1,nk2,nk3,jdim,idim))
        endif

        do  i3 = 1, nk3
        do  i2 = 1, nk2
        do  i1 = 1, nk1
          tijPir(i1,i2,i3,1:idim,1:jdim,1) = tij(i1,i2,i3,1:idim,1,1:jdim,1)
          tijPir(i1,i2,i3,1:idim,1:jdim,2) = tij(i1,i2,i3,1:idim,2,1:jdim,2)
          tjiPjr(i1,i2,i3,1:jdim,1:idim,1) = tji(i1,i2,i3,1:jdim,1,1:idim,1)
          tjiPjr(i1,i2,i3,1:jdim,1:idim,2) = tji(i1,i2,i3,1:jdim,2,1:idim,2)
          if (ferm) then
           do id = 1,idim
           do jd = 1,jdim
            cImtij = dimag(tij(i1,i2,i3,id,1,jd,1))
            cImtji = dimag(tji(i1,i2,i3,jd,2,id,2))
            Imtij(i1,i2,i3,id,jd) = cImtij*(0d0,1d0)
            Imtji(i1,i2,i3,jd,id) = cImtji*(0d0,1d0)
           enddo
           enddo
!          dimag(Imtij(i1,i2,i3,1:jdim,1:idim)) = 0d0
!          dimag(Imtji(i1,i2,i3,1:jdim,1:idim)) = 0d0
          endif
        enddo
        enddo
        enddo

        do  id = 1, idim
        tptr = 0
        if (ferm) tptso = 0
        do  jd = 1, jdim

!          DelPi = gam1i(id+offpi)/2
!          DelPj = gam1j(jd+offpj)/2
          DelPi = gam1i(id)/2
          DelPj = gam1j(jd)/2
          l = ll(jd)

C   ... w1 = T^+, w2 = dPj T^-
          call dcopy(nfbz*2,tijPir(1,1,1,id,jd,1),1,wk,1)
          if (ldiag) then
            call zcopy(nfbz,tijPir(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
          else
            call zcopy(nfbz,tjiPjr(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
          endif
          call zscal(nfbz,delPj,wk(1,1,1,2),1)

C   ... FT T_ij(q) and p_j T_ji(q); overwrite w1 with product of result (=> R.S.)
          call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ... Lifetime (in progress)...
          if (ferm) then
            call dcopy(nfbz*2,Imtij(1,1,1,id,jd),1,wkso,1)
            if (ldiag) then
              call zcopy(nfbz,Imtij(1,1,1,jd,id),1,wkso(1,1,1,2),1)
            else
              call zcopy(nfbz,Imtji(1,1,1,jd,id),1,wkso(1,1,1,2),1)
            endif
            call zdscal(nfbz,dble(delPj),wkso(1,1,1,2),1)

            call fftz3c(wkso,wkso(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

            call dcopy(nfbz*2,Imtij(1,1,1,id,jd),1,wksoi,1)
            if (ldiag) then
              call zcopy(nfbz,Imtij(1,1,1,jd,id),1,wksoi(1,1,1,2),1)
            else
              call zcopy(nfbz,Imtji(1,1,1,jd,id),1,wksoi(1,1,1,2),1)
            endif
            call zdscal(nfbz,dimag(delPj),wksoi(1,1,1,2),1)

            call fftz3c(wksoi,wksoi(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          endif


C   ... Accumulate Jrr = J_R+T,R' = sum_LL' J_R+TL,R'L', Eq(4) in the notes
C       and sumRL(1) = dP_RL sum_TL' X_R+TL;R'L' dP_R'L', RHS of Eq(10)
C       and sumRL(3) = L-resolved version of J_00, Eq(7)
          tptir = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
          if (ferm) tptiso = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            tptir = tptir + wk(i1,i2,i3,1)
            if (ferm) tptiso = tptiso + dble(delPi)*wkso(i1,i2,i3,1)-dimag(delPi)*wksoi(i1,i2,i3,1)
            JrrR(i1,i2,i3) = JrrR(i1,i2,i3)+dimag(wpi2*delPi*wk(i1,i2,i3,1))
          enddo
          enddo
          enddo
          if (ldiag) sumRLr(3,0,0,id) = sumRLr(3,0,0,id) + wpi2*delPi*wk(1,1,1,1)
          tptr = tptr + delPi*tptir ! tpt = dP_RL sum_L' tpti_RL;R'L'
          if (ferm) tptso = tptso + tptiso ! p*Im t*p*Im t
        enddo
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
         sumRLr(1,0,0,id) = sumRLr(1,0,0,id) + wpi2*tptr
         if (ferm) sumRLr(7,0,0,id) = sumRLr(7,0,0,id) + tptso*(4*datan(1d0))/nfbz
        enddo

C   ... End Lifetime.............................................

C --- Debugging ... See Eqns 4 and 5 above (ldiag only) ---
        if (.not. ldiag .or. .not. lchk) goto 999
        do  id = 1, idim
          jd = id

C   ... dP = (P^+ - P^-)/2
!          DelPi = gam1i(id+offpi)/2  ! For the closest match to sum rule
          DelPi = gam1i(id)/2

C   ... dP_i sum_q (T^+_RR(i,j;q) - T^-_RR(i,j;q))
          cxr = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            cxr = cxr + tijPir(i1,i2,i3,id,jd,1) - tijPir(i1,i2,i3,id,jd,2)
          enddo
          enddo
          enddo
          cxr = cxr/nfbz
          sumRLr(2,0,0,id) = sumRLr(2,0,0,id) + wpi*delPi*cxr/2
        enddo

C ... ORBITAL exchange (in progress) and S-O .....

!        Jorb = 0
!        Jso = 0
!        nli = 0
!        do  id = 1, idim
!         tptr = 0
!         tptso = 0
!         nlj = 0
!         nli = nli + 1
!         mi = nli-ll(id)-1 ! m =-l,...,l
!!        print*,'id =',id,'nli =',nli,'mi = ',mi
!         do  jd = 1, jdim
!          nlj = nlj + 1
!          mj = nlj-ll(jd)-1
!!         print*,'jd =',jd,'nlj =',nlj,'mj = ',mj
!          DelPi = gam1i(id)/2
!          DelPj = gam1j(jd)/2
!C   ... w1 = T^+, w2 = M_j T^- (for orbit)
!          call dcopy(nfbz*2,tijPir(1,1,1,id,jd,1),1,wk,1)
!          if (ldiag) then
!            call zcopy(nfbz,tijPir(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
!          else
!            call zcopy(nfbz,tjiPjr(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
!          endif
!          call zscal(nfbz,mj,wk(1,1,1,2),1)
!          call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ... w1 = T^+, w2 = P T^- (for s-o)
!          call dcopy(nfbz*2,tijPir(1,1,1,id,jd,1),1,wkso,1)
!          if (ldiag) then
!            call zcopy(nfbz,tijPir(1,1,1,jd,id,2),1,wkso(1,1,1,2),1)
!          else
!            call zcopy(nfbz,tjiPjr(1,1,1,jd,id,2),1,wkso(1,1,1,2),1)
!          endif
!          call zscal(nfbz,DelPj,wkso(1,1,1,2),1)
!          call fftz3c(wkso,wkso(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)


!          tptir = 0
!          tptiso = 0
!          do  i3 = 1, nk3
!          do  i2 = 1, nk2
!          do  i1 = 1, nk1
!            tptir = tptir + wk(i1,i2,i3,1)
!            tptiso = tptiso + wkso(i1,i2,i3,1)
!            Jorb(i1,i2,i3) = Jorb(i1,i2,i3)+dimag(wpi2*mi*wk(i1,i2,i3,1))
!            Jso(i1,i2,i3) = Jso(i1,i2,i3)+dimag(wpi2*mi*wkso(i1,i2,i3,1))
!          enddo
!          enddo
!          enddo
!          if (ldiag) sumRLr(3,0,0,id) = sumRLr(3,0,0,id) + wpi2*mi*wk(1,1,1,1)
!          tptr = tptr + mi*tptir ! tpt = M sum_L' tMti_RL;R'L'
!          tptso = tptso + mi*tptiso ! tpt = M sum_L' tpti_RL;R'L'
!          if (nlj == 2*ll(jd)+1) nlj = 0
!         enddo ! jd
!C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
!         sumRLr(4,0,0,id) = sumRLr(4,0,0,id) + wpi2*tptr
!         sumRLr(5,0,0,id) = sumRLr(5,0,0,id) + wpi2*tptso
!         if (nli == 2*ll(id)+1) nli = 0
!        enddo !id
C ... End Orbital exchange (in progress) .....

        Jrrcp(:,:,:) = JrrR(:,:,:)
        sumRLcp(:,:,:,:) = sumRLr(:,:,:,:)
        deallocate(tijPir,tjiPjr)
        if (ferm) deallocate(Imtij,Imtji,Imtijr,Imtjir,Imtiji,Imtjii)
!        if (ferm .and. .false.) deallocate(Imtij,Imtji)

      else  ! Case lrel is not 2
!      endif ! lrel eq 2

C ... Non-relativistic part
      allocate(tijPi(nk1,nk2,nk3,idim,jdim))
      allocate(tjiPj(nk1,nk2,nk3,jdim,idim))
      if (ferm) then
       allocate(Imtij(nk1,nk2,nk3,idim,jdim))
       allocate(Imtji(nk1,nk2,nk3,jdim,idim))
       allocate(Imtijr(nk1,nk2,nk3,idim,jdim))
       allocate(Imtjir(nk1,nk2,nk3,jdim,idim))
       allocate(Imtiji(nk1,nk2,nk3,idim,jdim))
       allocate(Imtjii(nk1,nk2,nk3,jdim,idim))

       do  i3 = 1, nk3
       do  i2 = 1, nk2
       do  i1 = 1, nk1
         do id = 1,idim
         do jd = 1,jdim
          cImtij = dimag(tij(i1,i2,i3,id,1,jd,1))
          cImtji = dimag(tji(i1,i2,i3,jd,nspc,id,2))
          Imtij(i1,i2,i3,id,jd) = cImtij*(0d0,1d0)
          Imtji(i1,i2,i3,jd,id) = cImtji*(0d0,1d0)
         enddo
         enddo
       enddo
       enddo
       enddo
      endif ! ferm

      do  i3 = 1, nk3
      do  i2 = 1, nk2
      do  i1 = 1, nk1
C       tijPi(i1,i2,i3,1:idim,1:jdim) = matmul(gam2i(1:idim,1:idim),tij(i1,i2,i3,1:idim,1,1:jdim,1))/2
C       tjiPj(i1,i2,i3,1:jdim,1:idim) = matmul(gam2j(1:jdim,1:jdim),tji(i1,i2,i3,1:jdim,1,1:idim,2))/2

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)
        call zgemm('N','N',idim,jdim,idim,(.5d0,0d0),gam2i,norbx,wkij1,idim,(0d0,0d0),wkij,idim)
        tijPi(i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim) ! (.5d0,0d0) to account for (p-p)/2

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,nspc,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(.5d0,0d0),gam2j,norbx,wkji1,jdim,(0d0,0d0),wkji,jdim)
        tjiPj(i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        if (ferm) then

         gam2ireal(1:idim,1:jdim)=dble(gam2i(1:idim,1:jdim))*(1d0,0d0)
         gam2jreal(1:jdim,1:idim)=dble(gam2j(1:jdim,1:idim))*(1d0,0d0)
         gam2iim(1:idim,1:jdim)=dimag(gam2i(1:idim,1:jdim))*(1d0,0d0)
         gam2jim(1:jdim,1:idim)=dimag(gam2j(1:jdim,1:idim))*(1d0,0d0)

         wkij1(1:idim,1:jdim) = Imtij(i1,i2,i3,1:idim,1:jdim)
         call zgemm('N','N',idim,jdim,idim,(.5d0,0d0),gam2ireal,norbx,wkij1,idim,(0d0,0d0),wkij,idim)
         Imtijr(i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

         wkji1(1:jdim,1:idim) = Imtji(i1,i2,i3,1:jdim,1:idim)
         call zgemm('N','N',jdim,idim,jdim,(.5d0,0d0),gam2jreal,norbx,wkji1,jdim,(0d0,0d0),wkji,jdim)
         Imtjir(i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

         wkij1(1:idim,1:jdim) = Imtij(i1,i2,i3,1:idim,1:jdim)
         call zgemm('N','N',idim,jdim,idim,(.5d0,0d0),gam2iim,norbx,wkij1,idim,(0d0,0d0),wkij,idim)
         Imtiji(i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

         wkji1(1:jdim,1:idim) = Imtji(i1,i2,i3,1:jdim,1:idim)
         call zgemm('N','N',jdim,idim,jdim,(.5d0,0d0),gam2jim,norbx,wkji1,jdim,(0d0,0d0),wkji,jdim)
         Imtjii(i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

!!       cImtij = dimag(tijPi(i1,i2,i3,id,1,jd,1))
!!       cImtji = dimag(tji(i1,i2,i3,jd,nspc,id,2))
!        Imtij(i1,i2,i3,1:idim,1:jdim) = dimag(wkij(1:idim,1:jdim))*(0d0,1d0)
!        Imtji(i1,i2,i3,1:jdim,1:idim) = dimag(wkji(1:jdim,1:idim))*(0d0,1d0)
        endif

      enddo
      enddo
      enddo

!      call zprm('gam2i*tij',2,tijPi(2,3,4,1:idim,1:jdim),idim,idim,jdim)
!      call zprm('gam2j*tji',2,tjiPj(2,3,4,1:jdim,1:idim),jdim,jdim,idim)
!      return

      do  id = 1, idim ! Orbitals associated with this site
      tptcp = 0
      if (ferm) tptso = 0
      do  jd = 1, jdim
        call dcopy(nfbz*2,tijPi(1,1,1,id,jd),1,wk,1)
        call zcopy(nfbz,tjiPj(1,1,1,jd,id),1,wk(1,1,1,2),1)
        call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
        if (ferm) then
         call dcopy(nfbz*2,Imtijr(1,1,1,id,jd),1,wkso,1)
         call zcopy(nfbz,Imtjir(1,1,1,jd,id),1,wkso(1,1,1,2),1)
         call fftz3c(wkso,wkso(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
         call dcopy(nfbz*2,Imtiji(1,1,1,id,jd),1,wksoi,1)
         call zcopy(nfbz,Imtjii(1,1,1,jd,id),1,wksoi(1,1,1,2),1)
         call fftz3c(wksoi,wksoi(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
        endif
        tpticp = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
        if (ferm) tptiso = 0
        do  i3 = 1, nk3
        do  i2 = 1, nk2
        do  i1 = 1, nk1
          tpticp = tpticp + wk(i1,i2,i3,1)
          if (ferm) tptiso = tptiso + wkso(i1,i2,i3,1) - wksoi(i1,i2,i3,1)
          Jrrcp(i1,i2,i3) = Jrrcp(i1,i2,i3)+dimag(wpi2*wk(i1,i2,i3,1))
        enddo
        enddo
        enddo
        if (ldiag) sumRLcp(3,0,0,id) = sumRLcp(3,0,0,id)
     .         + wpi2*wk(1,1,1,1)
        if (ferm) tptso = tptso + tptiso
        tptcp = tptcp + tpticp
      enddo
      sumRLcp(1,0,0,id) = sumRLcp(1,0,0,id) + wpi2*tptcp
!      if (ferm) sumRLcp(7,0,0,id) = sumRLcp(7,0,0,id) + tptso/(4*datan(1d0))
      if (ferm) sumRLcp(7,0,0,id) = sumRLcp(7,0,0,id) + tptso*(4*datan(1d0))!/nfbz
      enddo
      if (ferm) deallocate(Imtij,Imtji,Imtijr,Imtjir,Imtiji,Imtjii)

C ... ORBITAL exchange (in progress) .....
      Jorb = 0
      Jso = 0
      nli = 0
      do  id = 1, idim
       tptcp = 0
       tptso = 0
       nlj = 0
       nli = nli + 1
       mi = nli-ll(id)-1
!      print*,'id =',id,'nli =',nli,'mi = ',mi
       do  jd = 1, jdim
          nlj = nlj + 1
          mj = nlj-ll(jd)-1
!         print*,'jd =',jd,'nlj =',nlj,'mj = ',mj
          call dcopy(nfbz*2,tij(1,1,1,id,1,jd,1),1,wk,1)
          call zcopy(nfbz,tji(1,1,1,jd,nspc,id,2),1,wk(1,1,1,2),1)
          call zdscal(nfbz,mj,wk(1,1,1,2),1)
          call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          call dcopy(nfbz*2,tij(1,1,1,id,1,jd,1),1,wkso,1)
          call zcopy(nfbz,tjiPj(1,1,1,jd,id),1,wkso(1,1,1,2),1)
          call fftz3c(wkso,wkso(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          tpticp = 0  ! tpti_RL;R'L' = sum TMT
          tptiso = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            tpticp = tpticp + wk(i1,i2,i3,1)
            tptiso = tptiso + wkso(i1,i2,i3,1)
            Jorb(i1,i2,i3) = Jorb(i1,i2,i3)+mi*dimag(wpi2*wk(i1,i2,i3,1))
            Jso(i1,i2,i3) = Jso(i1,i2,i3)+mi*dimag(wpi2*wkso(i1,i2,i3,1))
          enddo
          enddo
          enddo
!          if (ldiag) sumRLcp(3,0,0,id) = sumRLcp(3,0,0,id)
!     .         + wpi2*wk(1,1,1,1)
          tptcp = tptcp + mi*tpticp
          tptso = tptso + mi*tptiso
          if (nlj == 2*ll(jd)+1) nlj = 0
       enddo ! jd
       sumRLcp(4,0,0,id) = sumRLcp(4,0,0,id) + wpi2*tptcp
       sumRLcp(5,0,0,id) = sumRLcp(5,0,0,id) + wpi2*tptso
       if (nli == 2*ll(id)+1) nli = 0
      enddo ! id

C ... End ORBITAL exchange (in progress) .....

      do  i3 = 1, nk3
      do  i2 = 1, nk2
      do  i1 = 1, nk1

        tijPi(i1,i2,i3,1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)
     .                                - tij(i1,i2,i3,1:idim,nspc,1:jdim,2)
        wkij1(1:idim,1:jdim) = tijPi(i1,i2,i3,1:idim,1:jdim)
        call zgemm('N','N',idim,jdim,idim,(.5d0,0d0),gam2i,norbx,wkij1,idim,(0d0,0d0),wkij,idim)
        tijPi(i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

      enddo
      enddo
      enddo

C --- Debugging ... See Eqns 4 and 5 above (ldiag only) ---
      if (.not. ldiag .or. .not. lchk) goto 999
      do id = 1, idim
        jd=id
!      do  jd = 1, jdim

C   ... dP = (P^+ - P^-)/2
!        DelPi = gam2i(id,jd)/2

C   ...   dP_i sum_q (T^+_RR(i,j;q) - T^-_RR(i,j;q))
        cx = 0
        do i3 = 1, nk3
        do i2 = 1, nk2
        do i1 = 1, nk1
          cx = cx + tijPi(i1,i2,i3,id,jd)
        enddo
        enddo
        enddo
        cx = cx/nfbz
        sumRLcp(2,0,0,id) = sumRLcp(2,0,0,id) + wpi*cx/2
C       sumRL(4,0,0,id) = sumRL(4,0,0,id) + wpi*delPig*cx/2
!      enddo
      enddo
      deallocate(tijPi,tjiPj)
      endif ! lrel = 2

C   ... Lifetime (in work)...
C   ... Antropov Physica B 237-238 336 (1997)
C   ... lifetime Im J_ij = -pi*w*Tr_L [p*Im(T_up)*p*Im(T_down)]

!      if (ferm .and. .false.) then
!      allocate(tijPi(nk1,nk2,nk3,idim,jdim))
!      allocate(tjiPj(nk1,nk2,nk3,jdim,idim))
!      tptcp = 0; tpticp = 0
!        do  i3 = 1, nk3
!        do  i2 = 1, nk2
!        do  i1 = 1, nk1

!          wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)
!          call zgemm('N','N',idim,jdim,idim,(.5d0,0d0),gam2i,norb,wkij1,idim,(0d0,0d0),wkij,idim)
!          tijPi(i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

!          wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,2)
!          call zgemm('N','N',jdim,idim,jdim,(.5d0,0d0),gam2j,norb,wkji1,jdim,(0d0,0d0),wkji,jdim)
!          tjiPj(i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

!        enddo
!        enddo
!        enddo

!      call zprm('gam2i*tij',2,tijPi(2,3,4,1:idim,1:jdim),idim,idim,jdim)
!      call zprm('gam2j*tji',2,tjiPj(2,3,4,1:jdim,1:idim),jdim,jdim,idim)
!      return

!        do  id = 1, idim ! Orbitals associated with this site
!        tptcp = 0
!        do  jd = 1, jdim
!          call dcopy(nfbz*2,tijPi(1,1,1,id,jd),1,wk,1)
!          call zcopy(nfbz,tjiPj(1,1,1,jd,id),1,wk(1,1,1,2),1)
!          call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
!          tpticp = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
!          do  i3 = 1, nk3
!          do  i2 = 1, nk2
!          do  i1 = 1, nk1
!            tpticp = tpticp + wk(i1,i2,i3,1)
!            Jrrcp(i1,i2,i3) = Jrrcp(i1,i2,i3)+dimag(wpi2*wk(i1,i2,i3,1))
!          enddo
!          enddo
!          enddo
!          if (ldiag) sumRLcp(3,0,0,id) = sumRLcp(3,0,0,id)
!     .           + wpi2*wk(1,1,1,1)
!          tptcp = tptcp + tpticp
!        enddo
!        sumRLcp(7,0,0,id) = sumRLcp(7,0,0,id) + tptcp
!        print*,'PASAJJCPA: sum(4) = ',sumRLcp(7,0,0,id)
!        enddo

!       print*,'PASAJJCPA: sum(1)/wpi = ',sumRLcp(1,0,0,id)/wpi2
!        deallocate(tijPi,tjiPj)
!      endif

C   ... End lifetime ..................................


C       print*,'Jrrcp(i1,i2,i3)',Jrrcp(1,1,1),'JrrR(i1,i2,i3)',JrrR(1,1,1)



  999 call tcx('pasajjcpa')
      end


      subroutine pasaj2(nk1,nk2,nk3,k1,k2,k3,offpi,offpj,idim,jdim,
     .  nccomp,wz,ldpf,pfun,ldiag,lchk,tij,tji,ldi,ldj,wk,Jrr,sumRL)
C- ASA Longitudinal Magnetic exchange J(q) for a block of orbitals by FFT
C ----------------------------------------------------------------------
Ci Inputs
Ci   irep
Ci   nk1,nk2,nk3:mesh of k-points
Ci   offpi,offpj:offsets to potential functions
Ci   idim,jdim  :size of G(i,j) subblock to calculate
Ci   wz         :energy integration weight
Ci   pfun,ldpf  :potential functions, and leading dimension
Ci   tij,tji    :T-:matrix connecting orbitals ij and ji, aka bare gij,gji
Ci   ldiag      :ij is diagonal block, in which case tji is not
Ci               explicitly available, but must be computed from tij
Ci   lchk       :T if to check sum rule, and calculate sumRL
Ci   ldi,ldj    :dimensions of tij,tji:
Ci               tij = tij(nk1,nk2,nk3,ldi,ldj,nsp)
Ci               tji = tji(nk1,nk2,nk3,ldj,ldi,nsp)
Ci   wk         :d.c. work array dimensioned at least nk1*nk2*nk2*2
Ci
Co  Outputs
Co    Jrr       :Contribution from one energy point to integral for
Co               exchange interaction is accumulated connecting site R
Co               (pot fun offset offpi) to site R' (offset offpj)
Co               and those related to J' by a lattice vector.
Co
Co    sumRL     :quantities needed for sum rules and on-site exchange
Co              :interactions.  Real and imaginary parts retained.
Co              :sumRL(1,i=1..idim) =  sum_j J_ij computed by FT of (1)
Co              :sumRL(2,i=1..idim) = -sum_j J_ij computed by (5)
Co              :sumRL(3,i=1..idim) =  sum_j J_ij for j=1..idim (=J_00)
Co
Cr  Remarks
Cr
Cr    Adapted from pasajj ... Still experimental!
Cr    The exchange interaction is calculated from
Cr
Cr    pi J(q) = int^E_f de int dk Im Tr dP_i T+ij(k) dP_j T-_ji(k+q) (1)
Cr
Cr    pi J_i,j = int^E_f de Im Tr dP_i T+_i,j dP_j T-_j,i            (2)
Cr
Cr    where dP = (P+ - P-)/2, and i=i and j=R'L'                     (3)
Cr
Cr    There is a sum rule (where T_i = sum_k T_ii(k))
Cr
Cr      (T+_i - T-_i)/2 = sum_j T+_i,j dP_j T-_j,i                       (4)
Cr
Cr     which we integrate and scale by dP_i to extract on-site exchange
Cr
Cr      Im int de dP_RL (T+_RL - T-_RL) = sum_R'L J_RL,R'L'          (5)
Cr
Cr    When ldiag is set, the l.h.s. is evaluated in sumRL(2) for site R.
Cr    sumRL(1) accumulates the contribution to the  r.h.s for pairs
Cr    (R,R') where R' is some neighbor R' in the unit cell, and the
Cr    same neighbor translated by the nk1*nk2*nk3 lattice vectors.
Cr
Cr    The convolution T(k) * T(k+q) is done by FFT; result Jrr is
Cr    returned in real space.
C ----------------------------------------------------------------------
      implicit none
C Passed variables
      logical ldiag,lchk
      integer offpi,offpj,idim,jdim,ldi,ldj,nccomp
      integer ldpf,nk1,nk2,nk3,k1,k2,k3,nsp
      parameter (nsp=2)
      double precision Jrr(nk1,nk2,nk3)
      double complex wz,pfun(ldpf,nsp),sumRL(11*nccomp,idim)
      double complex tij(nk1,nk2,nk3,ldi,ldj,nsp),wk(k1,k2,k3,2)
      double complex tji(nk1,nk2,nk3,ldj,ldi,nsp)
C Local variables
      integer id,jd,i1,i2,i3,nfbz,ii
      double complex DelPi,DelPj,wpi,wpi2,tpt,cx,tpti

      call tcn('pasaj2')
C ... Weight for integration 1/pi de ...
      wpi = wz/(4*datan(1d0))
      wpi2 = wpi
C ... Use the following if add symmetrized contribution to J
      wpi2 = wpi/2
      nfbz = nk1*nk2*nk3

C --- Loop over orbital pairs ---
      do  id = 1, idim
      tpt = 0
      do  jd = 1, jdim

C   ... dP = (P^+ - P^-)/2
        DelPi = (pfun(id+offpi,1) - pfun(id+offpi,2))/2
        DelPj = (pfun(jd+offpj,1) - pfun(jd+offpj,2))/2

        do  ii = 1,2

C   ... w1 = T^+/-, w2 = dPj T^+/- : NB: change for case nk1 /= k1 etc
        call dcopy(nfbz*2,tij(1,1,1,id,jd,ii),1,wk,1)
        if (ldiag) then
          call zcopy(nfbz,tij(1,1,1,jd,id,ii),1,wk(1,1,1,2),1)
        else
          call zcopy(nfbz,tji(1,1,1,jd,id,ii),1,wk(1,1,1,2),1)
        endif
        call zscal(nfbz,delPj,wk(1,1,1,2),1)

C   ... Overwrite w1 with convolution of T_ij and p_j T_ji
        call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C      call zprm3(' wk %i',id+100+jd,wk,nk1,nk2,nk3)

C   ... Accumulate p_i T+/-_ij pj T+/-_ji into Jrr and sum T+/-_ij pj T+/-_ji
        tpti = 0
        do  i3 = 1, nk3
        do  i2 = 1, nk2
        do  i1 = 1, nk1
          tpti = tpti + wk(i1,i2,i3,1)
          Jrr(i1,i2,i3) = Jrr(i1,i2,i3)+dimag(wpi2*delPi*wk(i1,i2,i3,1))
        enddo
        enddo
        enddo
        if (ldiag) sumRL(3,id) = sumRL(3,id) + wpi2*delPi*wk(1,1,1,1)

        tpt = tpt + delPi*tpti
C       tpt = tpt + 1*tpti

        enddo


CC   ... This is the symmetrized, presumably equivalent contribution
CC   ... w1 = (Pi^+ - Pi^-)/2 T^-, w2 = (Pj^+ - Pj^-)/2 T^+
C        if (ldiag) then
C          do  40  i3 = 1, nk3
C          do  40  i2 = 1, nk2
C          do  40  i1 = 1, nk1
C            wk(i1,i2,i3,1) = tij(i1,i2,i3,id,jd,2)
C            wk(i1,i2,i3,2) = delPj*tij(i1,i2,i3,jd,id,1)
C   40     continue
C        else
C          do  42  i3 = 1, nk3
C          do  42  i2 = 1, nk2
C          do  42  i1 = 1, nk1
C            wk(i1,i2,i3,1) = tij(i1,i2,i3,id,jd,2)
C            wk(i1,i2,i3,2) = delPj*tji(i1,i2,i3,jd,id,1)
C   42     continue
C        endif
CC   ... Overwrite p_i T_ij with convolution of p_i T_ij and p_j T_ji
C        call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,nk1,nk2,nk3,30,-1)
C
CC   ... Accumulate 1/2 p_i T-_ij pj T+_ji into Jrr
C        do  44  i3 = 1, nk3
C        do  44  i2 = 1, nk2
C        do  44  i1 = 1, nk1
C        tpt = tpt + wk(i1,i2,i3,1)
C   44   Jrr(i1,i2,i3) = Jrr(i1,i2,i3) +dimag(wpi2*delPi*wk(i1,i2,i3,1))
C        tpt = tpt + delPi*tpti

      enddo
      sumRL(1,id) = sumRL(1,id) + wpi2*tpt
      enddo

C --- Debugging ... check sum rule (ldiag only) ---
      if (.not. ldiag .or. .not. lchk) goto 999
C     tupdn = 0
      do  110  id = 1, idim
        jd = id

C   ... dP = (P^+ - P^-)/2
        DelPi = (pfun(id+offpi,1) - pfun(id+offpi,2))/2

C   ... dP sum_q (T^+(q) - T^-(q))
        cx = 0
        do  130  i3 = 1, nk3
        do  130  i2 = 1, nk2
        do  130  i1 = 1, nk1
          cx = cx + tij(i1,i2,i3,id,jd,1) - tij(i1,i2,i3,id,jd,2)
  130   continue
        cx = cx/nfbz

        sumRL(2,id) = sumRL(2,id) + wpi*delPi*cx/2
  110 continue

  999 call tcx('pasaj2')
      end

      subroutine pasajf2(irep,nk1,nk2,nk3,k1,k2,k3,ib,nbas,offH,ldpf,
     .  pfb,dpfb,ddpfb,wz,wk,gr,gc,JRR,sumRL)
C- Kernel called by asajft
      implicit none
      logical lasttime
      integer irep,nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer nk1,nk2,nk3,k1,k2,k3,ldpf,nsp,ib,nbas,offH(n0H,nkap0,nbas)
      double precision JRR(nk1,nk2,nk3,nbas,nbas)
      parameter (nsp=2)
      double complex wk(k1,k2,k3,2),wz
      double complex pfb(ldpf,nsp),dpfb(ldpf,nsp),ddpfb(ldpf,nsp)
C     double complex gr(nk1,nk2,nk3,nlma,lidim,2)
C     double complex gc(nk1,nk2,nk3,lidim,nlma,2)
      double precision sumRL(22,*)
      double complex gr(nk1,nk2,nk3,*),gc(nk1,nk2,nk3,*)
C Local
      integer nkfbz,isp,nlmi,jb,offpj,iblk,offpi,nlmj,nlma,
     .  ldim,lidim,offg,offai
      integer, allocatable :: ipiv(:)
      real(8), allocatable :: work(:)
      complex(8), allocatable :: tuptdown(:),chi(:)

      call rx('pasajf2: handle irep')

      call tcn('pasajf2')
      nkfbz = nk1*nk2*nk3

      nlma  = offH(4,1,ib+1)-offH(4,1,ib)
      ldim  = offH(1,1,nbas+1)
      lidim = offH(4,1,nbas+1)

      allocate(tuptdown(nk1*nk2*nk3*ldim**2))
      allocate(chi(ldim**2))
      allocate(ipiv(ldim**2))
      allocate(work(ldim**2))
C     lhdim = lidim + offH(3,1,nbas+1)

C ... Scale gr, gc to bare representation ... scale only lower block
      offg = 0
      do  10  iblk = 1, 1
      nlmi  = offH(iblk,1,ib+1) - offH(iblk,1,ib)
      offpi = offH(iblk,1,ib) + (iblk-1)*ldim
      do  12  isp = 1, 2
        call pvgfe2(isp,nsp,nkfbz,offpi,0,nlmi,lidim,ldpf,dpfb,ddpfb,
     .    nlma,lidim,.true.,gr(1,1,1,1+offg),gr)
        call pvgfe2(isp,nsp,nkfbz,0,offpi,lidim,nlmi,ldpf,dpfb,ddpfb,
     .    lidim,nlma,.true.,gc(1,1,1,1+lidim*offg),gc)
   12 continue
      offg = offg + nlmi
   10 continue

      offai = offH(4,1,ib)
      offpi = offH(1,1,ib)
      do  20  jb = 1, nbas
        offpj = offH(1,1,jb)
        nlmj = offH(1,1,jb+1) - offH(1,1,jb)
        lasttime = .false.
        if (ib == nbas) then
          if (jb == nbas) then
            lasttime = .true.
          endif
        endif
        call pasajj2(nk1,nk2,nk3,k1,k2,k3,offpi,offpj,nlmi,nlmj,wz,ldpf,
     .    pfb,ib == jb,.true.,gr(1,1,1,1+offpj*nlma),gc(1,1,1,1+offpj),
     .    nlma,lidim,wk,JRR(1,1,1,1,1),sumRL(1,1+offai),tuptdown,
     .    chi,nbas,ldim,lasttime,offH(1,1,1))
   20 continue

      call tcx('pasajf2')
      end

      subroutine pasajj2(nk1,nk2,nk3,k1,k2,k3,offpi,offpj,idim,jdim,wz,
     .  ldpf,pfun,ldiag,lchk,tij,tji,ldi,ldj,wk,Jrr,sumRL,tuptdown,
     .  chi,nbas,ldim,lasttime,offH)
C- ASA Magnetic exchange J(q) for a block of orbitals by FFT
C ----------------------------------------------------------------------
Ci Inputs
Ci   nk1,nk2,nk3:mesh of k-points
Ci   offpi,offpj:offsets to potential functions
Ci   idim,jdim  :size of G(i,j) subblock to calculate
Ci   wz         :energy integration weight
Ci   pfun,ldpf  :potential functions, and leading dimension
Ci   tij,tji    :T-:matrix connecting orbitals ij and ji, aka bare gij,gji
Ci   ldiag      :ij is diagonal block, in which case tji is not
Ci               explicitly available, but must be computed from tij
Ci   lchk       :T if to check sum rule, and calculate sumRL
Ci   ldi,ldj    :dimensions of tij,tji:
Ci               tij = tij(nk1,nk2,nk3,ldi,ldj,nsp)
Ci               tji = tji(nk1,nk2,nk3,ldj,ldi,nsp)
Ci   wk         :d.c. work array dimensioned at least nk1*nk2*nk2*2
Ci
Co  Outputs
Co    Jrr       :Contribution from one energy point to integral for
Co               exchange interaction is accumulated connecting site R
Co               (pot fun offset offpi) to site R' (offset offpj)
Co               and those related to J' by a lattice vector.
Co
Co    sumRL     :quantities needed for sum rules and on-site exchange
Co              :interactions.  Real and imaginary parts retained.
Co              :sumRL(1,i=1..idim) =  sum_j J_ij computed by FT of (1)
Co              :sumRL(2,i=1..idim) = -sum_j J_ij computed by (5)
Co              :sumRL(3,i=1..idim) =  sum_j J_ij for j=1..idim (=J_00)
Co
Cr  Remarks
Cr
Cr    The exchange interaction is calculated from
Cr
Cr    pi J(q) = int^E_f de int dk Im Tr dP_i T+ij(k) dP_j T-_ji(k+q) (1)
Cr
Cr    pi J_i,j = int^E_f de Im Tr dP_i T+_i,j dP_j T-_j,i            (2)
Cr
Cr    where dP = (P+ - P-)/2, and i=i and j=R'L'                     (3)
Cr
Cr    There is a sum rule (where T_i = sum_k T_ii(k))
Cr
Cr      (T+_i - T-_i)/2 = sum_j T+_i,j dP_j T-_j,i                       (4)
Cr
Cr     which we integrate and scale by dP_i to extract on-site exchange
Cr
Cr      Im int de dP_RL (T+_RL - T-_RL) = sum_R'L J_RL,R'L'          (5)
Cr
Cr    When ldiag is set, the l.h.s. is evaluated in sumRL(2) for site R.
Cr    sumRL(1) accumulates the contribution to the  r.h.s for pairs
Cr    (R,R') where R' is some neighbor R' in the unit cell, and the
Cr    same neighbor translated by the nk1*nk2*nk3 lattice vectors.
Cr
Cr    The convolution T(k) * T(k+q) is done by FFT; result Jrr is
Cr    returned in real space.
C ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical ldiag,lchk,lasttime
      integer offpi,offpj,idim,jdim,ldim,ldi,ldj
      integer ldpf,nk1,nk2,nk3,k1,k2,k3,nsp,nbas,offH(nbas+1)
      parameter (nsp=2)
      double precision Jrr(nk1,nk2,nk3,nbas,nbas)
C     double precision work(ldim**2)
      double complex wz,pfun(ldpf,nsp),sumRL(11,idim)
      double complex tij(nk1,nk2,nk3,ldi,ldj,nsp),wk(k1,k2,k3,2)
      double complex tji(nk1,nk2,nk3,ldj,ldi,nsp),chi(ldim,ldim)
      double complex tuptdown(nk1,nk2,nk3,ldim,ldim),delta(ldim)
C ... Local parameters
      integer id,jd,nfbz,orbind1,orbind2,kpt1,kpt2,kpt3,ib,jb
      double complex DelPi,DelPj,wpi,wpi2,tpt
C     double precision cx,tpti

      call tcn('pasajj2')
C ... Weight for integration 1/pi de ...
      wpi = wz/(4*datan(1d0))
      wpi2 = wpi
C ... Use the following if add symmetrized contribution to J
C     wpi2 = wpi/2
      nfbz = nk1*nk2*nk3

C --- Loop over orbital pairs --- ---
      do  id = 1, idim
      tpt = 0
      do  jd = 1, jdim

C   ... dP = (P^+ - P^-)/2
        DelPi = (pfun(id+offpi,1) - pfun(id+offpi,2))/2
        DelPj = (pfun(jd+offpj,1) - pfun(jd+offpj,2))/2

C   ... w1 = T^+, w2 = dPj T^- : NB: change for case nk1 /= k1 etc
        call dcopy(nfbz*2,tij(1,1,1,id,jd,1),1,
     .    tuptdown(1,1,1,offpi+id,offpj+jd),1)
        if (ldiag) then
          call zcopy(nfbz,tij(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
        else
          call zcopy(nfbz,tji(1,1,1,jd,id,2),1,wk(1,1,1,2),1)
        endif
C        call zscal(nfbz,delPj,wk(1,1,1,2),1)

C   ... Overwrite w1 with convolution of T_ij and p_j T_ji
        call fftz3c(tuptdown(1,1,1,offpi+id,offpj+jd),wk(1,1,1,2),nk1,
     .  nk2,nk3,k1,k2,k3,20,-1)

C      call zprm3(' wk %i',id+100+jd,wk,nk1,nk2,nk3)

C   ... Accumulate pi T+_ij pj T-_ji into Jrr and sum T+_ij pj T-_ji
c        tpti = 0
c        do  i3 = 1, nk3
c        do  i2 = 1, nk2
c        do  i1 = 1, nk1
c          tpti = tpti + wk(i1,i2,i3,1)
c          Jrr(i1,i2,i3) = Jrr(i1,i2,i3)+dimag(wpi2*delPi*wk(i1,i2,i3,1))
c        enddo
c        enddo
c        enddo
c        if (ldiag) sumRL(3,id) = sumRL(3,id) + wpi2*delPi*wk(1,1,1,1)

c        tpt = tpt + delPi*tpti
C       tpt = tpt + 1*tpti

CC   ... This is the symmetrized, presumably equivalent contribution
CC   ... w1 = (Pi^+ - Pi^-)/2 T^-, w2 = (Pj^+ - Pj^-)/2 T^+
C        if (ldiag) then
C          do  40  i3 = 1, nk3
C          do  40  i2 = 1, nk2
C          do  40  i1 = 1, nk1
C            wk(i1,i2,i3,1) = tij(i1,i2,i3,id,jd,2)
C            wk(i1,i2,i3,2) = delPj*tij(i1,i2,i3,jd,id,1)
C   40     continue
C        else
C          do  42  i3 = 1, nk3
C          do  42  i2 = 1, nk2
C          do  42  i1 = 1, nk1
C            wk(i1,i2,i3,1) = tij(i1,i2,i3,id,jd,2)
C            wk(i1,i2,i3,2) = delPj*tji(i1,i2,i3,jd,id,1)
C   42     continue
C        endif
CC   ... Overwrite p_i T_ij with convolution of p_i T_ij and p_j T_ji
C        call fftz3c(wk,wk(1,1,1,2),nk1,nk2,nk3,nk1,nk2,nk3,30,-1)
C
CC   ... Accumulate 1/2 p_i T-_ij pj T+_ji into Jrr
C        do  44  i3 = 1, nk3
C        do  44  i2 = 1, nk2
C        do  44  i1 = 1, nk1
C        tpt = tpt + wk(i1,i2,i3,1)
C   44   Jrr(i1,i2,i3) = Jrr(i1,i2,i3) +dimag(wpi2*delPi*wk(i1,i2,i3,1))
C        tpt = tpt + delPi*tpti

      enddo
c      sumRL(1,id) = sumRL(1,id) + wpi2*tpt
      enddo
C ... If tuptdown is complete, invert, decorate with dP, contract over l
      if (lasttime) then
c  comment out matrix invert for now
c          do kpt1=1,nk1
c          do kpt2=1,nk2
c          do kpt3=1,nk3
c            chi=tuptdown(kpt1,kpt2,kpt3,:,:)
c            call zgetrf(ldim,ldim,chi,ldim,ipiv,check)
c            call zgetri(ldim,chi,ldim,ipiv,work,ldim**2,check)
c            tuptdown(kpt1,kpt2,kpt3,:,:)=chi
c          enddo
c          enddo
c          enddo
        do orbind1=1,ldim
          delta(orbind1)=(pfun(orbind1,1)-pfun(orbind1,2))/2
        enddo
        do  orbind1 = 1, ldim
          do  orbind2 = 1, ldim
            call zscal(nfbz,delta(orbind1)*delta(orbind2),
     .        tuptdown(1,1,1,orbind1,orbind2),1)
            call fftz3(tuptdown(1,1,1,orbind1,orbind2),nk1,nk2,nk3,
     .        k1,k2,k3,1,0,1)
          enddo
        enddo

        do  ib = 1, nbas
          do  jb = 1, nbas
C           write(*,*) "Jrr ib jb", Jrr(1,1,1,ib,jb), ib, jb
            do  orbind1 = (offH(ib)+1), (offH(ib+1)+1)
            do  orbind2 = (offH(jb)+1), (offH(jb+1)+1)
              do  kpt1 = 1, nk1
              do  kpt2 = 1, nk2
              do  kpt3 = 1, nk3
                Jrr(kpt1,kpt2,kpt3,ib,jb)=Jrr(kpt1,kpt2,kpt3,ib,jb)
     .            +dimag(tuptdown(kpt1,kpt2,kpt3,orbind1,orbind1)*wz)
              enddo
              enddo
              enddo
            enddo
            enddo
          enddo
        enddo
      endif
C --- Debugging ... check sum rule (ldiag only) --- ---
c      if (.not. ldiag .or. .not. lchk) goto 999
C     tupdn = 0
c      do  110  id = 1, idim
c        jd = id

C    ... dP = (P^+ - P^-)/2
c        DelPi = (pfun(id+offpi,1) - pfun(id+offpi,2))/2

C    ... dP sum_q (T^+(q) - T^-(q))
c        cx = 0
c        do  130  i3 = 1, nk3
c        do  130  i2 = 1, nk2
c        do  130  i1 = 1, nk1
c          cx = cx + tij(i1,i2,i3,id,jd,1) - tij(i1,i2,i3,id,jd,2)
c  130   continue
c        cx = cx/nfbz

C       tupdn = tupdn + delPi*cx/2
C       sumR(2) = sumR(2) + wpi*tupdn
c        sumRL(2,id) = sumRL(2,id) + wpi*delPi*cx/2
c  110 continue
C  999 continue
      call tcx('pasajj2')
      end

C      subroutine pvgfevcp(s_site,mode,icomp,ib,ldpf,pfun,norb,norbx,nl,nbas,lrel,nspc,indxsh,
C     .            sfvrtx,gamiicp,gamiirel,mskm)
CC- Vertex associated with a single site for exchange interactions
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :1s digit
CCi         :0 Not used here
CCi         :1 fill gamii with pfun(1)-pfun(2)
CCi         :10s digit
CCi         :0 make no contribution to gamii from a CPA vertex
CCi         :1 Copy diagonal part of CPA 12 vertex within component icomp, site ib
CCi         :  to appropriate place in gamii (nothing done if icomp<=0)
CCi         :2 Copy diagonal part of CPA 21 vertex within component icomp, site ib
CCi         :  to appropriate place in gamii (nothing done if icomp<=0)
CCi   icomp :component to CPA vertex. icomp=0 if not CPA
CCi   ib    :site containing CPA vertex
CCi   ldpf  :dimensions potential functions P
CCi   pfun  :vector of potential functions
CCi   norb  :dimensions sfvrtx
CCi   sfvrtx:spin flip vertex for all components in site ib
CCo Outputs
CCi   gamii :Pfun(up)-Pfun(dn) in the absence of vertices,
CCi         :CPA case: diagonal part of (spin) vertex is substituted
CCi         :at site ib for given component
CCl Local variables
CCl         :
CCr Remarks
CCr    For sites other than CPA sites, it is assumed that s_site%pfr has been filled
CCu Updates
CCu   02 Sep 15 (Vishina) extended to relativistic case
CCu   13 Jul 15 (Vishina) first created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,icomp,ib,ldpf,norb,norbx,nspc,lrel,nl,indxsh(*),nbas
C      double complex pfun(ldpf,2),sfvrtx(norb,norb,2,*)
C      double complex gamiicp(norbx,norbx)
C!      double complex gamiirel(norbx*2,norbx*2)
C      double complex gamiirel(ldpf)
CC ... For structures
C      include 'structures.h'
C      type(str_site)::  s_site(*)
CC ... Local parameters
C      logical ldsite
C      double complex wk(norb,nspc,norb,nspc),wk2(norb,2)
C      double complex wkkmu(norb*2,norb*2),wkms(ldpf,2,2)
C      double complex fun(nspc,nspc),xx
C      double complex wk3(norb,norb,nspc,nspc),wk4(norb,norb,nspc,nspc)
C      double complex mskm(0:nl-1,2*nl)
C      integer idx(2,0:nl-1,2*(nl+1))
C      integer i,j,mode1,ncomp,l,imu,i1,i2,pfdim,ibas,ib1,ib2
C      integer m11,m12,m21,m22,lmr11,lmr22,lmrd11,lmrd22,idxn,io
C
CC     mode0 = mod(mode,10)
C      ncomp = s_site(ib)%ncomp
C      mode1 = mod(mode/10,10)
C      pfdim = ldpf
C
C      if (icomp < 1 .and. nspc == 1) then    ! Non CPA branch, nonrelativistic.  P is a diagonal matrix
C        gamiicp = 0
C        call zcopy(norb*2,s_site(ib)%pfr,1,wk2,1)
C        do  i = 1, norb
C          gamiicp(i,i) = wk2(i,1) - wk2(i,2)
C        enddo
C
CC ... non-CPA branch
C      elseif (icomp < 1) then ! Non CPA branch
C
C        if (lrel == 2) then
C
C
C
C
C          gamiirel(:) = 0
C          call zcopy(norb*norb*4,s_site(ib)%pfr,1,wk,1)
C
C          print *, wk(1,1,1,1),wk(2,1,2,1),wk(3,1,3,1),wk(4,1,4,1),wk(5,1,5,1)
C          stop
C
C          call zcopy(ldpf*4,s_site(ib)%pfr,1,wkms,1)
C          do i = 1, ldpf
C            gamiirel(i) = wkms(i,1,1) - wkms(i,2,2)
C          enddo
C        else                    !lrel ne 2
C          if (mode1 == 1 .or. mode1 == 2) then
CC          print*,'pvgfevcp mode1 1 or 2'
C            call zcopy(norb*norb*4,s_site(ib)%pfr,1,wk,1)
C            do  i = 1, norb
C              do  j = 1, norb
C                gamiicp(i,j) = wk(i,1,j,1) - wk(i,2,j,2)
C              enddo
C            enddo
C          endif
C        endif                   ! lrel eq 2
CC ... CPA branch
C      else
C        if (mode1 == 1 .or. mode1 == 2) then
C          do  i = 1, norb
C            do  j = 1, norb
C              gamiicp(i,j) = sfvrtx(i,j,mode1,icomp)
C            enddo
C          enddo
C        endif
C      endif
C
C      end

      subroutine pvgfev(mode,icomp,ib,offH,ldpf,pfun,norb,sfvrtx,gamii)
C- Diagonal part of vertex for exchange interactions
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 make no contribution to gamii from pfun (crystalline atoms)
Ci         :1 fill gamii with pfun(1)-pfun(2)
Ci         :10s digit
Ci         :0 make no contribution to gamii from a CPA vertex
Ci         :1 Copy diagonal part of CPA 12 vertex within component icomp, site ib
Ci         :  to appropriate place in gamii (nothing done if icomp<=0)
Ci         :2 Copy diagonal part of CPA 21 vertex within component icomp, site ib
Ci         :  to appropriate place in gamii (nothing done if icomp<=0)
Ci   icomp :component to CPA vertex. icomp=0 if not CPA
Ci   ib    :site containing CPA vertex
Ci   ldpf  :dimensions potential functions P
Ci   pfun  :vector of potential functions
Ci   norb  :dimensions sfvrtx
Ci   sfvrtx:spin flip vertex for all components in site ib
Co Outputs
Ci   gamii :Pfun(up)-Pfun(dn) in the absence of vertices,
Ci         :CPA case: diagonal part of (spin) vertex is substituted
Ci         :at site ib for given component
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Jan 13 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,icomp,ib,ldpf,norb
      double complex sfvrtx(norb,norb,2,*),pfun(ldpf,2),gamii(ldpf)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,ib+1)
C ... Local parameters
      integer i,mode0,mode1,nlmi,offpi

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

      if (mode0 /= 0) then
        do  i = 1, ldpf
          gamii(i) = pfun(i,1) - pfun(i,2)
        enddo
      endif
      if (icomp < 1) return

      if (mode1 /= 0) then
        offpi = offH(1,1,ib)
        nlmi  = offH(1,1,ib+1) - offH(1,1,ib)
        if (mode1 == 1 .or. mode1 == 2) then
          do  i = 1, nlmi
            gamii(i+offpi) = sfvrtx(i,i,mode1,icomp)
          enddo
        endif
      endif

      end

      subroutine pasachi(nk1,nk2,nk3,nspc,k1,k2,k3,idim,jdim,mxcomp,wz,
     .  ldiag,Gij,Gji,ldi,ldj,wk,JrrR,sumRLr)
C- ASA Magnetic exchange J(q) by FFT for a pair of sites or CPA components
C ----------------------------------------------------------------------
Ci Inputs
Ci   nk1,nk2,nk3:mesh of k-points
Ci   nspc       : 2 for non-collinear case, otherwise 1
Ci   offpi,offpj:offsets to potential functions
Ci   idim,jdim  :size of G(i,j) subblock to calculate
Ci   ferm       :1 if fermi energy
Ci   wz         :energy integration weight
Ci   Gij,Gji    :G connecting orbitals ij and ji, aka Gij,Gji
Ci              :Gij,Gji are supplied on a uniform mesh of q points
Ci   ldiag      :ij is diagonal block, in which case Gji is not
Ci               explicitly available, but must be computed from Gij
Ci   ldi,ldj    :dimensions of tij,tji:
Ci               tij = tij(nk1,nk2,nk3,ldi,nspc,ldj,2)
Ci               tji = tji(nk1,nk2,nk3,ldj,nspc,ldi,2)
Ci   wk         :d.c. work array dimensioned at least nk1*nk2*nk2*2*nspc
Co  Outputs
Co    Jrr       :Contribution from one energy point to integral for
Co               exchange interaction is accumulated connecting site R+T
Co               (pot fun offset offpi) to site R' (offset offpj)
Co               for all translation vectors T.
Co               Jrr is a contraction over L and L' of the orbital-resolved
Co               J_R+TL,RL' . Jrr is generated and returned in real space.
Co
Co    sumRL     :quantities needed for sum rules and on-site exchange.
Co              :interactions.  Real and imaginary parts retained.
Co              :In contrast to Jrr, sumRL is resolved by orbital L
Co              :sumRL(1,L=1..idim) = sum_TL' G*G, RHS of Eq(10) below
Co              :sumRL(2,id=1..idim) = 'one-site' rotation:
Co              :                       -1*LHS of Eq.(10) below
Co              :sumRL(3,id=1..idim)  onsite chi, AKA chi00
Co              :sumRL(4,id=1..idim)  onsite chi, AKA chi00 ???
Co              :sumRL(5,id=1..idim)  not used
Co              :sumRL(6,id=1..idim)  not used
Co              :sumRL(7,id=1..idim)  not used
Cr  Remarks
Cr    This routine is closely patterned after pasajjcpa, which see.
Cu Updates
Cu   14 Jan 16 (Vishina) adapted from pasajjcpa
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical ldiag
      integer idim,jdim,ldi,ldj,mxcomp
      integer nk1,nk2,nk3,k1,k2,k3,nsp,nspc
      parameter (nsp=2)
      double precision JrrR(nk1,nk2,nk3)
      double complex wz
      double complex sumRLr(11,0:mxcomp,0:mxcomp,idim)
      double complex Gij(nk1,nk2,nk3,ldi,nspc,ldj,2),wk(k1,k2,k3,2),wkup(k1,k2,k3,2),wkdown(k1,k2,k3,2)
      double complex Gji(nk1,nk2,nk3,ldj,nspc,ldi,2)
C ... Dynamically allocated local arrays
C ... Local parameters
      integer id,jd,i1,i2,i3,nfbz,fac,fac1,fac2,md,nd,ld,kd,ldown,lup,ndown,nup
      real sq
      double precision Mi,Mj,li,lj,Lupi,Ldowni,Lupj,Ldownj,faci3,facj3,facij3
      double complex wpi,wpi2,ggir,ggr,faci,facj,facij,cxr

      call tcn('pasachi')

C ... Weight for integration 1/pi de ...
      wpi = wz/(4*datan(1d0))
      wpi2 = wpi
C ... Use the following if add symmetrized contribution to J
C     wpi2 = wpi/2
      nfbz = nk1*nk2*nk3

C ... Orbit-Orbit (LS*G*LS*G aka Lz*G*Lz*G, see PRB 90, 094406 (2014)) (B||z)
        do  md = 1, jdim
        do  nd = 1, idim
        do  ld = 1, jdim
        do  kd = 1, idim
          if ((md == 5 .and. nd == 9) .or. (md == 9 .and. nd == 5) .or. !1 <d|Lz|d>
     .       (md == 6 .and. nd == 8) .or. (md == 8 .and. nd == 6)) then

          if ((ld == 5 .and. kd == 9) .or. (ld == 9 .and. kd == 5) .or. !2
     .       (ld == 6 .and. kd == 8) .or. (ld == 8 .and. kd == 6)) then

           if (md == 5 .and. nd == 9) then
                fac1 = 2
           elseif (md == 9 .and. nd == 5) then
                fac1 = -2
           elseif (md == 6 .and. nd == 8) then
                fac1 = 1
           else
                fac1 = -1
           endif
           if (ld == 5 .and. kd == 9) then
                fac2 = -2 ! i*i=-1 gives the opposite sign, all of the elements have i
           elseif (ld == 9 .and. kd == 5) then
                fac2 = 2
           elseif (ld == 6 .and. kd == 8) then
                fac2 = -1
           else
                fac2 = 1
           endif
           fac = fac1*fac2
C   ...    w1 = G^++(ml), w2 = G^++(kn) (for <d_m|Lz|d_n> G_ml <d_l|Lz|d_k> G_kn)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1) !was nl
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1) ! was km
             else
               call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
           endif
C   ...    w1 = G^--(ml), w2 = G^--(kn)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,2,ld,2),1,wkdown,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             else
               call zcopy(nfbz,Gji(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,2),1,wkdown,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,2),1,wkdown(1,1,1,2),1)
           endif

C   ...    FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
           call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
           call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...    Accumulate chi = G^++(nl)*G^++(km) + G^--(nl)*G^--(km)
!          if (.not. ldiag) then
             sumRLr(3,0,0,md) = sumRLr(3,0,0,md) + fac*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))
!          elseif (ldiag .and. md == 5 .and. ld .and. 5 .and. kd == 9 .and. nd == 9) then
!             sumRLr(3,0,0,md) = sumRLr(3,0,0,md) + fac*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
!          elseif (ldiag .and. md == 8 .and. ld .and. 8 .and. kd == 6 .and. nd == 6) then
!             sumRLr(3,0,0,md) = sumRLr(3,0,0,md) + fac*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
!          endif
          endif !2
          endif !1

          if ((md == 2 .and. nd == 4) .or. (md == 4 .and. nd == 2)) then !11  <p|Lz|p>

          if ((ld == 2 .and. kd == 4) .or. (ld == 4 .and. kd == 2)) then !22

           if (md == 2 .and. nd == 4) then
                fac1 = 1
           else
                fac1 = -1
           endif
           if (ld == 2 .and. kd == 4) then
                fac2 = -1 ! i*i=-1 gives the opposite sign
           else
                fac2 = 1
           endif
           fac = fac1*fac2
C   ...    w1 = G^++(nl), w2 = G^++(km) (for <m|Lz|n> G_nl <l|Lz|k> G_km)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
             else
               call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
           endif
C   ...    w1 = G^--(nl), w2 = G^--(km)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,2,ld,2),1,wkdown,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             else
               call zcopy(nfbz,Gji(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,2),1,wkdown,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,2),1,wkdown(1,1,1,2),1)
           endif

C   ...    FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
           call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
           call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...    Accumulate chi = G^++(nl)*G^++(km) + G^--(nl)*G^--(km)
!          if (.not. ldiag) then
             sumRLr(3,0,0,md) = sumRLr(3,0,0,md) + fac*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))
!          elseif (ldiag .and. md == 2 .and. ld .and. 2 .and. kd == 4 .and. nd == 4) then
!             sumRLr(3,0,0,md) = sumRLr(3,0,0,md) + fac*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
!          endif
          endif !22
          endif !11

        enddo
        enddo
        enddo
        enddo

C   ...   End OO
!        if (.not. ldiag) then
C ... Spin-Orbit (L+S- & L-S+ parts)............................................................
        if (nsp /= 2) then
          print*,'pasachi: nsp must be eq 2 to calculate chi_so with L+,L-'
        else
        do  md = 1, jdim
        do  ld = 1, idim

          if (md >= 6 .and. md <= 9 .and. ld >= 5 .and. ld <= 8) then

           kd = ld + 1
           nd = md - 1
           if (md == 6) then
                faci = (-1.0d0,0.0d0)
           elseif (md == 7) then
                faci = sqrt(3.0d0)*(0.0d0,1.0d0)
           elseif (md == 8) then
                faci = sqrt(3.0d0)*(1.0d0,0.0d0)
           else
                faci = (1.0d0,0.0d0)
           endif
           if (ld == 5) then
                facj = (-1.0d0,0.0d0)
           elseif (ld == 6) then
                facj = sqrt(3.0d0)*(0.0d0,-1.0d0)
           elseif (ld == 7) then
                facj = sqrt(3.0d0)*(1.0d0,0.0d0)
           else
                facj = (1.0d0,0.0d0)
           endif

           facij = faci*facj
C   ...    w1 = G^-+(ml), w2 = G^-+(kn) (for <d_m|L-S+|d_n> G_ml <d_l|L+S-|d_k> G_kn)
           call dcopy(nfbz*2,Gij(1,1,1,md,2,ld,1),1,wkup,1) !was nl
           if (ldiag) then
             call zcopy(nfbz,Gij(1,1,1,kd,2,nd,1),1,wkup(1,1,1,2),1) ! was km
           else
             call zcopy(nfbz,Gji(1,1,1,kd,2,nd,1),1,wkup(1,1,1,2),1)
           endif
C   ...    w1 = G^+-(lm), w2 = G^+-(nk) (for <d_m|L+S-|d_n> G_lm <d_l|L-S+|d_k> G_nk)
           call dcopy(nfbz*2,Gij(1,1,1,ld,1,md,2),1,wkdown,1) !was nl
           if (ldiag) then
             call zcopy(nfbz,Gij(1,1,1,nd,1,kd,2),1,wkdown(1,1,1,2),1) ! was km
           else
             call zcopy(nfbz,Gji(1,1,1,nd,1,kd,2),1,wkdown(1,1,1,2),1)
           endif

C   ...    FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
           call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...    Accumulate chi = <d_m|L-S+|d_n>*G^+-(ml)*<d_l|L+S-|d_k>*G^-+(kn)
           if (.not. ldiag) then
           sumRLr(3,0,0,md) = sumRLr(3,0,0,md) + facij*wpi2*(wkup(1,1,1,1)+wkdown(1,1,1,1))
           endif
          endif !1
        enddo
        enddo
        endif

!       endif !ldiag


C ... Spin-Orbit (LS*G*LS*G aka Lz*G*Lz*G) (B||x)..................................................
        do  md = 1, jdim
        do  nd = 1, idim
        do  ld = 1, jdim
        do  kd = 1, idim
          if ((md == 6 .and. nd == 7) .or. (md == 7 .and. nd == 6) .or. !1111 <d|Lz|d>
     .       (md == 6 .and. nd == 9) .or. (md == 9 .and. nd == 6) .or.
     .       (md == 5 .and. nd == 8) .or. (md == 8 .and. nd == 5)) then

          if ((ld == 6 .and. kd == 7) .or. (ld == 7 .and. kd == 6) .or. !2111
     .       (ld == 6 .and. kd == 9) .or. (ld == 9 .and. kd == 6) .or.
     .       (ld == 5 .and. kd == 8) .or. (ld == 8 .and. kd == 5)) then

           if (md == 6 .and. nd == 7) then
                faci3 = -sqrt(3.0d0)
           elseif (md == 6 .and. nd == 9) then
                faci3 = -1
           elseif (md == 5 .and. nd == 8) then
                faci3 = -1
           elseif (md == 8 .and. nd == 5) then
                faci3 = 1
           elseif (md == 7 .and. nd == 6) then
                faci3 = sqrt(3.0d0)
           else
                faci3 = 1
           endif
           if (ld == 6 .and. kd == 7) then
                facj3 = sqrt(3.0d0) ! i*i=-1 gives the opposite sign, all of the elements have i
           elseif (ld == 6 .and. kd == 9) then
                facj3 = 1
           elseif (ld == 5 .and. kd == 8) then
                facj3 = 1
           elseif (md == 8 .and. nd == 5) then
                facj3 = -1
           elseif (md == 7 .and. nd == 6) then
                facj3 = -sqrt(3.0d0)
           else
                facj3 = -1
           endif
           facij3 = faci3*facj3
C   ...    w1 = G^++(ml), w2 = G^++(kn) (for <d_m|Lz|d_n> G_ml <d_l|Lz|d_k> G_kn)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1) !was nl
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1) ! was km
             else
               call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
           endif
C   ...    w1 = G^--(nl), w2 = G^--(km)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,2,ld,2),1,wkdown,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             else
               call zcopy(nfbz,Gji(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,2),1,wkdown,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,2),1,wkdown(1,1,1,2),1)
           endif

C   ...    FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
           call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
           call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...    Accumulate chi = G^++(nl)*G^++(km) + G^--(nl)*G^--(km)
!          if (.not. ldiag) then
             sumRLr(6,0,0,md) = sumRLr(6,0,0,md) + facij3*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))
!          elseif (ldiag .and. md == 7 .and. ld .and. 7 .and. kd == 6 .and. nd == 6) then
!             sumRLr(6,0,0,md) = sumRLr(6,0,0,md) + facij3*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
!          elseif (ldiag .and. md == 9 .and. ld .and. 9 .and. kd == 6 .and. nd == 6) then
!             sumRLr(6,0,0,md) = sumRLr(6,0,0,md) + facij3*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
!          elseif (ldiag .and. md == 5 .and. ld .and. 5 .and. kd == 8 .and. nd == 8) then
!             sumRLr(6,0,0,md) = sumRLr(6,0,0,md) + facij3*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
!          endif
          endif !2111
          endif !1111

          if ((md == 4 .and. nd == 3) .or. (md == 2 .and. nd == 2)) then !110  <p|Lx|p> B||x

          if ((ld == 4 .and. kd == 3) .or. (ld == 2 .and. kd == 2)) then !220

           if (md == 2 .and. nd == 2) then
                fac1 = 1
           else
                fac1 = -1
           endif
           if (ld == 2 .and. kd == 2) then
                fac2 = -1 ! i*i=-1 gives the opposite sign
           else
                fac2 = 1
           endif
           fac = fac1*fac2
C   ...    w1 = G^++(nl), w2 = G^++(km) (for <m|Lz|n> G_nl <l|Lz|k> G_km)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
             else
               call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,1),1,wkup,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,1),1,wkup(1,1,1,2),1)
           endif
C   ...    w1 = G^--(nl), w2 = G^--(km)
           if (nspc == 2) then
             call dcopy(nfbz*2,Gij(1,1,1,md,2,ld,2),1,wkdown,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             else
               call zcopy(nfbz,Gji(1,1,1,kd,2,nd,2),1,wkdown(1,1,1,2),1)
             endif
           else
             call dcopy(nfbz*2,Gij(1,1,1,md,1,ld,2),1,wkdown,1)
             call zcopy(nfbz,Gji(1,1,1,kd,1,nd,2),1,wkdown(1,1,1,2),1)
           endif

C   ...    FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
           call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
           call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...    Accumulate chi = G^++(nl)*G^++(km) + G^--(nl)*G^--(km)
!          if (.not. ldiag) then
             sumRLr(6,0,0,md) = sumRLr(6,0,0,md) + fac*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))
!          elseif (ldiag .and. md == 2 .and. ld .and. 2 .and. kd == 4 .and. nd == 4) then
!             sumRLr(3,0,0,md) = sumRLr(3,0,0,md) + fac*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
!          endif
          endif !220
          endif !110
        enddo
        enddo
        enddo
        enddo


C   ...  S*G*L*G ..............................................................
        if (nsp /= 2) then
          print*,'pasachi: nsp must be eq 2 to calculate chi_so'
        else
        if (.not. ldiag) then
        do  md = 1, idim

          if (md == 5 .or. md == 6 .or. md == 7 .or. md == 8) then !33
           ld = md + 1

           sq = dsqrt(3.0d0)
           if (md == 5) then
                fac1 = -2
           elseif (md == 6) then
                fac1 = -1
           elseif (md == 7) then
                fac1 = 0
           else
                fac1 = 1
           endif
C   ...    w1 = G^++(mn), w2 = G^-+(nl) (for <d_m|L+|d_l>G_mn(up-up)G_nl(down-up))
           do  nd = 1, idim
             call dcopy(nfbz*2,Gij(1,1,1,md,1,nd,1),1,wkup,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,nd,2,ld,1),1,wkup(1,1,1,2),1)
             else
             call zcopy(nfbz,Gji(1,1,1,nd,2,ld,1),1,wkup(1,1,1,2),1)
             endif
C   ...      w1 = G^++(nm), w2 = G^+-(ln) (for  G_nm(up-up)<d_m|L+_j|d_l>  G_ln(up-down))
             call dcopy(nfbz*2,Gij(1,1,1,nd,1,md,1),1,wkdown,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,ld,1,nd,2),1,wkdown(1,1,1,2),1)
             else
             call zcopy(nfbz,Gji(1,1,1,ld,1,nd,2),1,wkdown(1,1,1,2),1)
             endif

C   ...      FT L+G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
             call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
             call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...      Accumulate chi = L+G^++(mn)*G^-+(nl) + G^++(nm)L+*G^+-(ln)
             sumRLr(5,0,0,md) = sumRLr(5,0,0,md) + fac1*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
C   ...      *2 to account for L+G_up-down*Gup-up, etc.
           enddo !nd
          endif !33

          if (md == 6 .or. md == 7 .or. md == 8 .or. md == 9) then !44
           ld = md - 1

           if (md == 6) then
                fac1 = -1
           elseif (md == 7) then
                fac1 = 0
           elseif (md == 8) then
                fac1 = 1
           else
                fac1 = 2
           endif
C   ...    w1 = G^--(mn), w2 = G^+-(nl) (for <d_m|L-|d_l>G_mn(down-down)G_nl(up-down))
           do  nd = 1, idim
             call dcopy(nfbz*2,Gij(1,1,1,md,2,nd,2),1,wkup,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,nd,2,ld,1),1,wkup(1,1,1,2),1)
             else
             call zcopy(nfbz,Gji(1,1,1,nd,1,ld,2),1,wkup(1,1,1,2),1)
             endif
C   ...      w1 = G^--(nm), w2 = G^-+(ln) (for  G_nm(down-down)<d_m|L+_j|d_l>  G_ln(down-up))
             call dcopy(nfbz*2,Gij(1,1,1,nd,2,md,2),1,wkdown,1)
             if (ldiag) then
               call zcopy(nfbz,Gij(1,1,1,ld,1,nd,2),1,wkdown(1,1,1,2),1)
             else
             call zcopy(nfbz,Gji(1,1,1,ld,2,nd,1),1,wkdown(1,1,1,2),1)
             endif

C   ...      FT L+G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
             call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
             call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...      Accumulate chi = L+G^++(mn)*G^-+(nl) + G^++(nm)L+*G^+-(ln)
             sumRLr(5,0,0,md) = sumRLr(5,0,0,md) + fac1*wpi2*(wkup(1,1,1,1) + wkdown(1,1,1,1))*2
C   ...      *2 to account for L+G_up-down*Gup-up, etc.
           enddo !nd

          endif !44
        enddo !md
        endif !ldiag
        endif !nsp = 2

C   ... orbit-orbit (z-z part).........................................................
       Mi = -1
       li = 0
       do  id = 1, idim
       ggr = 0
       Mi = Mi + 1
       if (id == 2 .or. id == 5 .or. id == 10) then
         li = li + 1
         Mi = -li
       endif
!       if (id == 1 .or. id == 3 .or. id == 7 .or. id == 13) then
!       Mi = 0
!       elseif (id == 2 .or. id == 6 .or. id == 12) then
!       Mi = -1
!       elseif (id == 4 .or. id == 8 .or. id == 14) then
!       Mi = 1
!       elseif (id == 5 .or. id == 11) then
!       Mi = -2
!       else
!       Mi = 2
!       endif

       lj = 0
       Mj = -1
       do  jd = 1, jdim
         Mj = Mj + 1
         if (jd == 2 .or. jd == 5 .or. jd == 10) then
           lj = lj + 1
           Mj = -lj
         endif

C   ...  w1 = G^++, w2 = G^++ (ZZ part)
         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,jd,1,id,1),1,wkup(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,jd,1,id,1),1,wkup(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          call zcopy(nfbz,Gji(1,1,1,jd,1,id,1),1,wkup(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  w1 = G^--, w2 = G^-- (ZZ part)
         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,2,jd,2),1,wkdown,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,jd,2,id,2),1,wkdown(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,jd,2,id,2),1,wkdown(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,2),1,wkdown,1)
          call zcopy(nfbz,Gji(1,1,1,jd,1,id,2),1,wkdown(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  Accumulate chi = G^-- * G^-- + G^++ * G^++,
         ggir = 0
         do  i3 = 1, nk3
         do  i2 = 1, nk2
         do  i1 = 1, nk1
           ggir = ggir + wkup(i1,i2,i3,1) + wkdown(i1,i2,i3,1)
!           JrrR(i1,i2,i3) = JrrR(i1,i2,i3)+dimag(wpi2*wk(i1,i2,i3,1))
         enddo
         enddo
         enddo
!         if (ldiag) sumRLr(2,0,0,id) = sumRLr(2,0,0,id) + wpi2*wk(1,1,1,1)
         ggr = ggr + ggir*Mi*Mj ! tpt = dP_RL sum_L' tpti_RL;R'L'



C   ...  w1 = G^+-, w2 = G^-+
         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,2),1,wkup,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,jd,2,id,1),1,wkup(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,jd,2,id,1),1,wkup(1,1,1,2),1)
          endif

C   ...   FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...   w1 = G^-+, w2 = G^+-
          call dcopy(nfbz*2,Gij(1,1,1,id,2,jd,1),1,wkdown,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,jd,1,id,2),1,wkdown(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,jd,1,id,2),1,wkdown(1,1,1,2),1)
          endif

C   ...   FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
          call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...   Accumulate chi = G^-- * G^--,
          ggir = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            ggir = ggir + wkup(i1,i2,i3,1) + wkdown(i1,i2,i3,1)
          enddo
          enddo
          enddo
!         if (ldiag) sumRLr(2,0,0,id) = sumRLr(2,0,0,id) + wpi2*wk(1,1,1,1)
          ggr = ggr + ggir*Mi*Mj ! tpt = dP_RL sum_L' tpti_RL;R'L'
         endif ! nspc = 2

C   ...   NON-ZZ part (will be used for non-z components) ?..........................................

C   ...  w1 = G^++, w2 = G^++

         ldown = jd - 1
         lup = jd + 1
         ndown = id - 1
         nup = id + 1
         if (id <= (idim-1) .and. jd >= 2) then
         Lupi = dsqrt(li*(li+1)-Mi*(Mi+1))
         Ldownj = dsqrt(lj*(lj+1)-Mj*(Mj-1))
         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,ldown,1,nup,1),1,wkup(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,ldown,1,nup,1),1,wkup(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          call zcopy(nfbz,Gji(1,1,1,ldown,1,nup,1),1,wkup(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,2,jd,2),1,wkdown,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,ldown,2,nup,2),1,wkdown(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,ldown,2,nup,2),1,wkdown(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,2),1,wkdown,1)
          call zcopy(nfbz,Gji(1,1,1,ldown,1,nup,2),1,wkdown(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  Accumulate chi = G^-- * G^-- + G^++ * G^++,
         ggir = 0
         do  i3 = 1, nk3
         do  i2 = 1, nk2
         do  i1 = 1, nk1
           ggir = ggir + wkup(i1,i2,i3,1) + wkdown(i1,i2,i3,1)
         enddo
         enddo
         enddo
!         ggr = ggr + ggir*Lupi*Ldownj ! tpt = dP_RL sum_L' tpti_RL;R'L'
         endif

         if (id >= 2 .and. jd <= (jdim-1)) then
          Ldowni = dsqrt(li*(li+1)-Mi*(Mi-1))
          Lupj = dsqrt(lj*(lj+1)-Mj*(Mj+1))
         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,lup,1,ndown,1),1,wkup(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,lup,1,ndown,1),1,wkup(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          call zcopy(nfbz,Gji(1,1,1,lup,1,ndown,1),1,wkup(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,2,jd,2),1,wkdown,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,lup,2,ndown,2),1,wkdown(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,lup,2,ndown,2),1,wkdown(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,2),1,wkdown,1)
          call zcopy(nfbz,Gji(1,1,1,lup,1,ndown,2),1,wkdown(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  Accumulate chi = G^-- * G^-- + G^++ * G^++,
         ggir = 0
         do  i3 = 1, nk3
         do  i2 = 1, nk2
         do  i1 = 1, nk1
           ggir = ggir + wkup(i1,i2,i3,1) + wkdown(i1,i2,i3,1)
         enddo
         enddo
         enddo
!         ggr = ggr + ggir*Ldowni*Lupj ! tpt = dP_RL sum_L' tpti_RL;R'L'
         endif


       enddo
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
         sumRLr(4,0,0,id) = sumRLr(4,0,0,id) + wpi2*ggr
       enddo

C   ... spin-spin chi .................................................................
       do  id = 1, idim
       ggr = 0
       do  jd = 1, jdim

C   ...  w1 = G^++, w2 = G^--
         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,jd,2,id,2),1,wkup(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,jd,2,id,2),1,wkup(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,1),1,wkup,1)
          call zcopy(nfbz,Gji(1,1,1,jd,1,id,2),1,wkup(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  w1 = G^++, w2 = G^--
         if (nspc == 2) then
          call dcopy(nfbz*2,Gij(1,1,1,id,2,jd,2),1,wkdown,1)
          if (ldiag) then
            call zcopy(nfbz,Gij(1,1,1,jd,1,id,1),1,wkdown(1,1,1,2),1)
          else
            call zcopy(nfbz,Gji(1,1,1,jd,1,id,1),1,wkdown(1,1,1,2),1)
          endif
         else
          call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,2),1,wkdown,1)
          call zcopy(nfbz,Gji(1,1,1,jd,1,id,1),1,wkdown(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  Accumulate chi = G^++ * G^--,
         ggir = 0
         do  i3 = 1, nk3
         do  i2 = 1, nk2
         do  i1 = 1, nk1
           ggir = ggir + wkup(i1,i2,i3,1) + wkdown(i1,i2,i3,1)
           JrrR(i1,i2,i3) = JrrR(i1,i2,i3)+dimag(wpi2*wkup(i1,i2,i3,1))+dimag(wpi2*wkdown(i1,i2,i3,1))
         enddo
         enddo
         enddo
         if (ldiag) sumRLr(2,0,0,id) = sumRLr(2,0,0,id) + wpi2*wkup(1,1,1,1) + wpi2*wkdown(1,1,1,1)
         ggr = ggr + ggir ! tpt = dP_RL sum_L' tpti_RL;R'L'

       enddo
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
         sumRLr(1,0,0,id) = sumRLr(1,0,0,id) + wpi2*ggr
       enddo

       if (nspc == 2) then    ! Do I need that???
       do  id = 1, idim
       ggr = 0
       do  jd = 1, jdim

C   ...  w1 = G^+-, w2 = G^+-
         call dcopy(nfbz*2,Gij(1,1,1,id,1,jd,2),1,wkup,1)
         if (ldiag) then
          call zcopy(nfbz,Gij(1,1,1,jd,1,id,2),1,wkup(1,1,1,2),1)
         else
          call zcopy(nfbz,Gji(1,1,1,jd,1,id,2),1,wkup(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkup,wkup(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  w1 = G^-+, w2 = G^-+
         call dcopy(nfbz*2,Gij(1,1,1,id,2,jd,1),1,wkdown,1)
         if (ldiag) then
          call zcopy(nfbz,Gij(1,1,1,jd,2,id,1),1,wkdown(1,1,1,2),1)
         else
          call zcopy(nfbz,Gji(1,1,1,jd,2,id,1),1,wkdown(1,1,1,2),1)
         endif

C   ...  FT G_ij(q) and G_ji(q); overwrite w1 with product of result (=> R.S.)
         call fftz3c(wkdown,wkdown(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ...  Accumulate chi = G^++ * G^--,
         ggir = 0
         do  i3 = 1, nk3
         do  i2 = 1, nk2
         do  i1 = 1, nk1
           ggir = ggir + wkup(i1,i2,i3,1) + wkdown(i1,i2,i3,1)
           JrrR(i1,i2,i3) = JrrR(i1,i2,i3)+dimag(wpi2*wkup(i1,i2,i3,1))+dimag(wpi2*wkdown(i1,i2,i3,1))
         enddo
         enddo
         enddo
!         if (ldiag) sumRLr(2,0,0,id) = sumRLr(2,0,0,id) + wpi2*wkup(1,1,1,1) + wpi2*wkdown(1,1,1,1
         ggr = ggr + ggir ! tpt = dP_RL sum_L' tpti_RL;R'L'

       enddo
C      Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
!       sumRLr(1,0,0,id) = sumRLr(1,0,0,id) + wpi2*ggr
       enddo
       endif

C   ... For chi_0 G_up-G_down
       if (ldiag) then
       do  id = 1, idim
          jd = id
C   ... (G^+_RR(i,j;q) - G^-_RR(i,j;q))
          cxr = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            if (nspc == 2) then
              cxr = cxr + Gij(i1,i2,i3,id,1,jd,1) - Gij(i1,i2,i3,id,2,jd,2)
            else
              cxr = cxr + Gij(i1,i2,i3,id,1,jd,1) - Gij(i1,i2,i3,id,1,jd,2)
            endif
          enddo
          enddo
          enddo
          cxr = cxr/nfbz
          sumRLr(7,0,0,id) = sumRLr(7,0,0,id) + wpi2*cxr
       enddo
       endif


      call tcx('pasachi')
      end


      subroutine pasajjnoncol(nk1,nk2,nk3,nspc,lrel,norb,norbx,ldpf,gamcomp,k1,k2,k3,
     .  idim,jdim,mxcomp,ferm,wz,pfloci,pflocj,gam2i,gam2j,ldiag,lchk,tij,tji,sfvrtxreli,
     .  sfvrtxrelj,Jrrcp,sumRLcpz)
C- ASA Magnetic exchange J(q) by FFT for a pair of sites or CPA components
C ----------------------------------------------------------------------
Ci Inputs
Ci   nk1,nk2,nk3:mesh of k-points
Ci   nspc       : 2 for non-collinear case, otherwise 1
Ci   offpj:offsets to potential functions
Ci   idim,jdim  :size of G(i,j) subblock to calculate
Ci   ferm       :1 if fermi energy
Ci   wz         :energy integration weight
Ci   gam2j     :Crystal case: (Pup-Pdn)
Ci              :CPA case: diagonal part of (spin12) vertex.
Ci   gam2i     :Crystal case: (Pup-Pdn)
Ci              :CPA case: diagonal part of (spin21) vertex.
Ci   tij,tji    :T-:matrix connecting orbitals ij and ji, aka gij,gji
Ci              :tij,tji are supplied on a uniform mesh of q points
Ci   ldiag      :ij is diagonal block, in which case tji is not
Ci               explicitly available, but must be computed from tij
Ci   lchk       :T if to generate sumRL(2,:,:,:) for sum rule
Ci   wk         :d.c. work array dimensioned at least nk1*nk2*nk2*2*nspc
Ci   qnur       :relativistic ASA energy moments
Co  Outputs
Co    Jrr       :Contribution from one energy point to integral for
Co               exchange interaction is accumulated connecting site R+T
Co               (pot fun offset offpi) to site R' (offset offpj)
Co               for all translation vectors T.
Co               Jrr is a contraction over L and L' of the orbital-resolved
Co               J_R+TL,RL' . Jrr is generated and returned in real space.
Co
Co    sumRL     :quantities needed for sum rules and on-site exchange.
Co              :interactions.  Real and imaginary parts retained.
Co              :In contrast to Jrr, sumRL is resolved by orbital L
Co              :sumRL(1,L=1..idim) = sum_TL' J_R+TL,RL', RHS of Eq(10) below
Co              :sumRL(2,id=1..idim) = 'one-site' rotation:
Co              :                       -1*LHS of Eq.(10) below
Co              :sumRL(3,id=1..idim) = onsite J_RL,RL, AKA J00
Co              :sumRL(4,id=1..idim) orb. exch.
Co              :sumRL(5,id=1..idim) orb-spin exch.
Co              :sumRL(6,id=1..idim) spin-orb exch.
Co              :sumRL(7,id=1..idim) reserved for the lifetime
Cr  Remarks
Cr    This routine is an extension to the ASA exchange maker pasajj for the
Cr    noncollinear and relativistic cases.
Cr    For the formulas see
Cr     Antropov Physica B 237-238 336 (1997)
Cr    This routine is closely patterned after pasajj, which see.
Cu Updates
Cu   22 Dec 16 (Vishina) adapted from pasajj
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical ldiag,lchk,ferm
      integer idim,jdim,mxcomp,lrel,ldpf,gamcomp
      integer nk1,nk2,nk3,k1,k2,k3,nsp,nspc,norb,norbx
      parameter (nsp=2)
      double precision Jrrcp(nk1,nk2,nk3,3,3)
      double complex wz,pflocj(norbx,2,norbx,2),pfloci(norbx,2,norbx,2)
      double complex sumRLcpz(11,3,3,0:mxcomp,0:mxcomp,idim)
      double complex wk11(k1,k2,k3,2),wk12(k1,k2,k3,2),wk21(k1,k2,k3,2),wk22(k1,k2,k3,2)
      double complex wk00(k1,k2,k3,2),wkxx(k1,k2,k3,2),wkyy(k1,k2,k3,2),wkzz(k1,k2,k3,2)
      double complex tij(nk1,nk2,nk3,idim,2,jdim,2),tji(nk1,nk2,nk3,jdim,2,idim,2)
      double complex gam2j(gamcomp,gamcomp),gam2i(gamcomp,gamcomp)
      complex(8) sfvrtxreli(4,4,4,norb,norb),sfvrtxrelj(4,4,4,norb,norb)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: tijx(:,:,:,:,:),tjix(:,:,:,:,:),tijy(:,:,:,:,:),tjiy(:,:,:,:,:)
      complex(8), allocatable :: tijz(:,:,:,:,:),tjiz(:,:,:,:,:),tij0(:,:,:,:,:),tji0(:,:,:,:,:)
      complex(8), allocatable :: pii(:,:,:),pjj(:,:,:)
      complex(8), allocatable :: ptij(:,:,:,:,:,:),ptji(:,:,:,:,:,:)
      complex(8), allocatable :: Scalar(:,:,:,:,:,:)
C     complex(8), allocatable :: tijkm(:,:,:,:,:),tjikm(:,:,:,:,:)
C     complex(8), allocatable :: wkkm(:,:)!,wkijr(:,:,:,:),wkjir(:,:,:,:)
C     complex(8), allocatable :: Imtij(:,:,:,:,:),Imtji(:,:,:,:,:),pttotij(:,:,:),pttotji(:,:,:)
C ... Local parameters
      integer id,jd,i1,i2,i3,nfbz,l,ll,al,be
      double precision JrrR(nk1,nk2,nk3,3,3),Dz(nk1,nk2,nk3,3) ! cImtij,cImtji
      double complex wkji(jdim,idim),wkji1(jdim,idim),sigma(4,2,2),wk111
      double complex wkij(idim,jdim),wkij1(idim,jdim)
      double complex wpi,wpi2,tptir,tptr,cxr
      double complex tptirxx,tptrxx,tptiryy,tptryy,tptirzz,tptrzz
      double complex sumRLr(11,3,3,0:mxcomp,0:mxcomp,idim)
      call tcn('pasajjnoncol')

C     call zprm('gam2j',2,gam2j,norb,jdim,jdim)

C ... Weight for integration 1/pi de ...
      wpi = wz/(4*datan(1d0))
      wpi2 = wpi
C ... Use the following if add symmetrized contribution to J
C     wpi2 = wpi/2
      nfbz = nk1*nk2*nk3

      JrrR(:,:,:,:,:) = Jrrcp(:,:,:,:,:)
      sumRLr(:,:,:,:,:,:) = sumRLcpz(:,:,:,:,:,:)

      call dpzero(sigma,16)
C ... Pauli matrices; sigma(1)=I, sigma(2)=sigma_x, etc
!      sigma(1,1,1)=(1d0,0d0)
!      sigma(1,2,2)=(1d0,0d0)
!      sigma(2,1,2)=(1d0,0d0)
!      sigma(2,2,1)=(1d0,0d0)
!      sigma(3,1,2)=(0d0,-1d0)
!      sigma(3,2,1)=(0d0,1d0)
!      sigma(4,1,1)=(1d0,0d0)
!      sigma(4,2,2)=(-1d0,0d0)

!      sigma(1,1,1)=(.5d0,0d0)
!      sigma(1,2,2)=(.5d0,0d0)
!      sigma(2,1,2)=(.5d0,0d0)
!      sigma(2,2,1)=(.5d0,0d0)
!      sigma(3,1,2)=(0d0,-.5d0)
!      sigma(3,2,1)=(0d0,.5d0)
!      sigma(4,1,1)=(.5d0,0d0)
!      sigma(4,2,2)=(-.5d0,0d0)

C ...... Spin-spin exchange...................................................

      allocate(tij0(nk1,nk2,nk3,idim,jdim),tji0(nk1,nk2,nk3,jdim,idim))
      allocate(tijx(nk1,nk2,nk3,idim,jdim),tjix(nk1,nk2,nk3,jdim,idim))
      allocate(tijy(nk1,nk2,nk3,idim,jdim),tjiy(nk1,nk2,nk3,jdim,idim))
      allocate(tijz(nk1,nk2,nk3,idim,jdim),tjiz(nk1,nk2,nk3,jdim,idim))
      if (gamcomp == 1) then  ! non CPA
        allocate(ptij(7,nk1,nk2,nk3,idim,jdim),ptji(7,nk1,nk2,nk3,jdim,idim))
      else
        allocate(ptij(16,nk1,nk2,nk3,idim,jdim),ptji(16,nk1,nk2,nk3,jdim,idim))
      endif
      allocate(pii(4,idim,jdim),pjj(4,jdim,idim))
      allocate(Scalar(3,nk1,nk2,nk3,idim,jdim))
!      allocate(pttotij(nk1,nk2,nk3),pttotji(nk1,nk2,nk3))

C   ... Lifetime (in progress)................................................
C   ... Antropov Physica B 237-238 336 (1997)
C   ... lifetime Im J_ij = -pi*w*Tr_L [p*Im(T_up)*p*Im(T_down)]
!        if (ferm) then
!         allocate(Imtij(nk1,nk2,nk3,idim,jdim))
!         allocate(Imtji(nk1,nk2,nk3,jdim,idim))
!       endif
C ... Make p_x=p*sigma_x, etc
!        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)
!        call zgemm('N','N',idim,jdim,idim,(.5d0,0d0),gam2i,norbx,wkij1,idim,(0d0,0d0),wkij,idim)
!        tijPi(i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim) ! (.5d0,0d0) to account for (p-p)/2
      Scalar = 0

      if (gamcomp == 1) then  ! non CPA
      do id = 1, idim
      do jd = 1, jdim
        pii(1,id,jd) =  (pfloci(id,1,jd,1) + pfloci(id,2,jd,2))/2d0
        pjj(1,jd,id) =  (pflocj(jd,1,id,1) + pflocj(jd,2,id,2))/2d0
        pii(2,id,jd) =  (pfloci(id,1,jd,2) + pfloci(id,2,jd,1))/2d0 !(0d0,0d0)!
        pjj(2,jd,id) =  (pflocj(jd,1,id,2) + pflocj(jd,2,id,1))/2d0
        pii(3,id,jd) =  (-pfloci(id,1,jd,2) + pfloci(id,2,jd,1))*(0d0,1d0)/2d0
        pjj(3,jd,id) =  (-pflocj(jd,1,id,2) + pflocj(jd,2,id,1))*(0d0,1d0)/2d0
        pii(4,id,jd) =  (pfloci(id,1,jd,1) - pfloci(id,2,jd,2))/2d0
        pjj(4,jd,id) =  (pflocj(jd,1,id,1) - pflocj(jd,2,id,2))/2d0
      enddo
      enddo
C ... Make T_0,x,y,z

      do  i3 = 1, nk3
      do  i2 = 1, nk2
      do  i1 = 1, nk1
       do id = 1, idim
       do jd = 1, jdim
        tij0(i1,i2,i3,id,jd) =  (tij(i1,i2,i3,id,1,jd,1) + tij(i1,i2,i3,id,2,jd,2))/2d0
        tji0(i1,i2,i3,jd,id) =  (tji(i1,i2,i3,jd,1,id,1) + tji(i1,i2,i3,jd,2,id,2))/2d0
        tijx(i1,i2,i3,id,jd) =  (tij(i1,i2,i3,id,1,jd,2) + tij(i1,i2,i3,id,2,jd,1))/2d0
        tjix(i1,i2,i3,jd,id) =  (tji(i1,i2,i3,jd,1,id,2) + tji(i1,i2,i3,jd,2,id,1))/2d0
        tijy(i1,i2,i3,id,jd) =  (-tij(i1,i2,i3,id,1,jd,2) + tij(i1,i2,i3,id,2,jd,1))*(0d0,1d0)/2d0
        tjiy(i1,i2,i3,jd,id) =  (-tji(i1,i2,i3,jd,1,id,2) + tji(i1,i2,i3,jd,2,id,1))*(0d0,1d0)/2d0
        tijz(i1,i2,i3,id,jd) =  (tij(i1,i2,i3,id,1,jd,1) - tij(i1,i2,i3,id,2,jd,2))/2d0
        tjiz(i1,i2,i3,jd,id) =  (tji(i1,i2,i3,jd,1,id,1) - tji(i1,i2,i3,jd,2,id,2))/2d0
        if (ldiag) then
          Scalar(1,i1,i2,i3,id,jd) = Scalar(1,i1,i2,i3,id,jd) + pii(2,id,jd)*tijx(i1,i2,i3,id,jd)
          Scalar(2,i1,i2,i3,id,jd) = Scalar(2,i1,i2,i3,id,jd) + pii(3,id,jd)*tijy(i1,i2,i3,id,jd)
          Scalar(3,i1,i2,i3,id,jd) = Scalar(3,i1,i2,i3,id,jd) + pii(4,id,jd)*tijz(i1,i2,i3,id,jd)
        endif
       enddo
       enddo

!        ptij(1,i1,i2,i3,:,:) = matmul(pii(1,:,:),tij0(i1,i2,i3,:,:))
!       ptji(1,i1,i2,i3,:,:) = matmul(pjj(1,:,:),tji0(i1,i2,i3,:,:))
!       ptij(2,i1,i2,i3,:,:) = matmul(pii(2,:,:),tij0(i1,i2,i3,:,:))
!       ptji(2,i1,i2,i3,:,:) = matmul(pjj(2,:,:),tji0(i1,i2,i3,:,:))
!       ptij(3,i1,i2,i3,:,:) = matmul(pii(3,:,:),tij0(i1,i2,i3,:,:))
!       ptji(3,i1,i2,i3,:,:) = matmul(pjj(3,:,:),tji0(i1,i2,i3,:,:))
!       ptij(4,i1,i2,i3,:,:) = matmul(pii(4,:,:),tij0(i1,i2,i3,:,:))
!       ptji(4,i1,i2,i3,:,:) = matmul(pjj(4,:,:),tji0(i1,i2,i3,:,:))
!       ptji(5,i1,i2,i3,:,:) = matmul(pii(2,:,:),tijx(i1,i2,i3,:,:))
!       ptji(5,i1,i2,i3,:,:) = matmul(pjj(2,:,:),tjix(i1,i2,i3,:,:))
!       ptij(6,i1,i2,i3,:,:) = matmul(pii(3,:,:),tijy(i1,i2,i3,:,:))
!       ptji(6,i1,i2,i3,:,:) = matmul(pjj(3,:,:),tjiy(i1,i2,i3,:,:))
!       ptij(7,i1,i2,i3,:,:) = matmul(pii(4,:,:),tijz(i1,i2,i3,:,:))
!       ptji(7,i1,i2,i3,:,:) = matmul(pjj(4,:,:),tjiz(i1,i2,i3,:,:))

!       do id = 1, 2
!       do jd = 1, 2
!         ptij(1,i1,i2,i3,id,jd) = pii(1,id,jd)*tij0(i1,i2,i3,id,jd)
!         ptji(1,i1,i2,i3,jd,id) = pjj(1,jd,id)*tji0(i1,i2,i3,jd,id)
!         ptij(2,i1,i2,i3,id,jd) = pii(1,id,jd)*tijx(i1,i2,i3,id,jd)
!         ptji(2,i1,i2,i3,jd,id) = pjj(1,jd,id)*tjix(i1,i2,i3,jd,id)
!         ptij(3,i1,i2,i3,id,jd) = pii(1,id,jd)*tijy(i1,i2,i3,id,jd)
!         ptji(3,i1,i2,i3,jd,id) = pjj(1,jd,id)*tjiy(i1,i2,i3,jd,id)
!         ptij(4,i1,i2,i3,id,jd) = pii(1,id,jd)*tijz(i1,i2,i3,id,jd)
!         ptji(4,i1,i2,i3,jd,id) = pjj(1,jd,id)*tjiz(i1,i2,i3,jd,id)

        wkij1(1:idim,1:jdim) = tij0(i1,i2,i3,1:idim,1:jdim) !p0T0
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),pii(4,1:idim,1:idim),idim,wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(1,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji0(i1,i2,i3,1:jdim,1:idim)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),pjj(4,1:jdim,1:jdim),jdim,wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(1,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij0(i1,i2,i3,1:idim,1:jdim) !?? pxT0 ?? p(2, ) = p_x
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),pii(4,1:idim,1:idim),idim,wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(2,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji0(i1,i2,i3,1:jdim,1:idim)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),pjj(4,1:jdim,1:jdim),jdim,wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(2,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij0(i1,i2,i3,1:idim,1:jdim) !?? pyT0 ?? p(3, ) = p_y
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),pii(4,1:idim,1:idim),idim,wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(3,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji0(i1,i2,i3,1:jdim,1:idim)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),pjj(4,1:jdim,1:jdim),jdim,wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(3,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij0(i1,i2,i3,1:idim,1:jdim) !?? pzT0 p(4, ) = p_z
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),pii(4,1:idim,1:idim),idim,wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(4,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji0(i1,i2,i3,1:jdim,1:idim)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),pjj(4,1:jdim,1:jdim),jdim,wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(4,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

       wkij1(1:idim,1:jdim) = tijx(i1,i2,i3,1:idim,1:jdim) !pxTx p(2, ) = p_x
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),pii(4,1:idim,1:idim),idim,wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(5,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tjix(i1,i2,i3,1:jdim,1:idim)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),pjj(4,1:jdim,1:jdim),jdim,wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(5,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tijy(i1,i2,i3,1:idim,1:jdim) !pyTy p(3, ) = p_y
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),pii(4,1:idim,1:idim),idim,wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(6,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tjiy(i1,i2,i3,1:jdim,1:idim)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),pjj(4,1:jdim,1:jdim),jdim,wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(6,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tijz(i1,i2,i3,1:idim,1:jdim) !pzTz p(4, ) = p_z
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),pii(4,1:idim,1:idim),idim,wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(7,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tjiz(i1,i2,i3,1:jdim,1:idim)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),pjj(4,1:jdim,1:jdim),jdim,wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(7,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)
!       enddo
!       enddo

!          if (ferm) then
!          do id = 1,idim
!          do jd = 1,jdim
!           cImtij = dimag(tij(i1,i2,i3,id,1,jd,1))
!           cImtji = dimag(tji(i1,i2,i3,jd,2,id,2))
!           Imtij(i1,i2,i3,id,jd) = cImtij*(1d0,0d0)
!           Imtji(i1,i2,i3,jd,id) = cImtji*(1d0,0d0)
!          enddo
!          enddo
!         endif
      enddo !i3
      enddo !i2
      enddo !i1

      do al = 1,3 ! alpha=x,y,z
      do be = 1,3 ! beta=x,y,z

      do  id = 1, idim
        tptr = 0
!       if (ferm) tptso = 0

        do  jd = 1, jdim

          l = ll(jd)
          call dpzero(wk11,nk1*nk2*nk3*2*2)
          call dpzero(wk22,nk1*nk2*nk3*2*2)
C   ... (see Physica B 237-238 336 (1997))

C   ... wk1 = pT, wk2 = pT
          call dcopy(nfbz*2,ptij(al+4,:,:,:,id,jd),1,wk11,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(be+4,:,:,:,jd,id),1,wk11(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(be+4,1,1,1,jd,id),1,wk11(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)


C   ... wk1 = pT, wk2 = pT
          if (al == be) then
           call dcopy(nfbz*2,ptij(al+1,:,:,:,id,jd),1,wk22,1)
           if (ldiag) then
             call zcopy(nfbz,ptij(be+1,:,:,:,jd,id),1,wk22(1,1,1,2),1)
           else
             call zcopy(nfbz,ptji(be+1,1,1,1,jd,id),1,wk22(1,1,1,2),1)
           endif
C   ...    FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
           call fftz3c(wk22,wk22(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
          endif

          wk21 = - wk11 + wk22
!          call fftz3(wk21,nk1,nk2,nk3,k1,k2,k3,2,0,-1)


C   ... Lifetime (in progress)...
!          if (ferm) then
!            call dcopy(nfbz*2,Imtij(1,1,1,id,jd),1,wkso,1)
!            if (ldiag) then
!              call zcopy(nfbz,Imtij(1,1,1,jd,id),1,wkso(1,1,1,2),1)
!            else
!              call zcopy(nfbz,Imtji(1,1,1,jd,id),1,wkso(1,1,1,2),1)
!            endif
!            call zdscal(nfbz,dble(delPj),wkso(1,1,1,2),1)

!            call fftz3c(wkso,wkso(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
!         endif

C   ... Accumulate Jrr = J_R+T,R' = sum_LL' J_R+TL,R'L', Eq(4) in the notes
C       and sumRL(1) = dP_RL sum_TL' X_R+TL;R'L' dP_R'L', RHS of Eq(10)
C       and sumRL(3) = L-resolved version of J_00, Eq(7)
          tptir = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
!         if (ferm) tptiso = 0
!          call dpzero(wk12,nk1*nk2*nk3*2*2)
          wk111 = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
!           if (ldiag) then
!            tptir = tptir - ptij(al+4,i1,i2,i3,id,jd)*ptij(be+4,i1,i2,i3,jd,id)
!           wk12(i1,i2,i3,1) = wk12(i1,i2,i3,1) - ptij(al+4,i1,i2,i3,id,jd)*ptij(be+4,i1,i2,i3,jd,id)
!           if (i1 == 1 .and. i2 == 1 .and. i3 == 1) wk111 = wk111 - ptij(al+4,i1,i2,i3,id,jd)*ptij(be+4,i1,i2,i3,jd,id)
!          else
!            tptir = tptir - ptij(al+4,i1,i2,i3,id,jd)*ptji(be+4,i1,i2,i3,jd,id)
!          endif
!          if (al == be) then
!             if (ldiag) then
!               tptir = tptir + ptij(al+1,i1,i2,i3,id,jd)*ptij(be+1,i1,i2,i3,jd,id)
!              wk12(i1,i2,i3,1) = wk12(i1,i2,i3,1) + ptij(al+1,i1,i2,i3,id,jd)*ptij(be+1,i1,i2,i3,jd,id)
!              if (i1 == 1 .and. i2 == 1 .and. i3 == 1) wk111 = wk111  + ptij(al+1,i1,i2,i3,id,jd)*ptij(be+1,i1,i2,i3,jd,id)
!            else
!               tptir = tptir + ptij(al+1,i1,i2,i3,id,jd)*ptji(be+1,i1,i2,i3,jd,id)
!            endif
!          endif
            tptir = tptir + wk21(i1,i2,i3,1)
!           if (al == be) tptir = tptir + wk22(i1,i2,i3,1)
!!          if (ferm) tptiso = tptiso + wkso(i1,i2,i3,1)
           JrrR(i1,i2,i3,al,be) = JrrR(i1,i2,i3,al,be)+dimag(wpi2*wk21(i1,i2,i3,1))
!          if (al == be) JrrR(i1,i2,i3,al,be) = JrrR(i1,i2,i3,al,be) + dimag(wpi2*wk22(i1,i2,i3,1))
          enddo
          enddo
          enddo
          if (ldiag) then
!           call fftz3(wk12,nk1,nk2,nk3,k1,k2,k3,1,0,-1)
!           if (al == be) sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) + wpi2*ptij(al+1,1,1,1,id,jd)*ptij(be+1,1,1,1,jd,id)
!          sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) - wpi2*ptij(al+4,1,1,1,id,jd)*ptij(be+4,1,1,1,jd,id)
           sumRLr(3,al,be,0,0,id) = sumRLr(3,al,be,0,0,id) + wpi2*wk21(1,1,1,1)
!          sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) + wpi2*wk111
          endif

!          if (ldiag) sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) - wpi2*wk11(1,1,1,1)
!          if (ldiag .and. al == be) sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) + wpi2*wk22(1,1,1,1)
          tptr = tptr + tptir ! tpt = dP_RL sum_L' tpti_RL;R'L'
!         if (ferm) tptso = tptso + dble(delPi)*tptiso ! p*Im t*p*Im t
        enddo ! jd
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
         sumRLr(1,al,be,0,0,id) = sumRLr(1,al,be,0,0,id) + wpi2*tptr!/nfbz
!         if (ferm) sumRLr(7,0,0,id) = sumRLr(7,0,0,id) + wpi2*tptso
      enddo ! id

C   ... For Dzyaloshinski vector

C   ... Dzyaloshinski vector D = Int Im Tr[p_i T^0_ji p_j T_ji]
      do  id = 1, idim
        tptr = 0
!       if (ferm) tptso = 0

        do  jd = 1, jdim
          l = ll(jd)
C   ... wk1 = pT_0(11), wk2 = pT_x,y,z(11)
          call dcopy(nfbz*2,ptij(1,1,1,1,id,jd),1,wk11,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(be+1,1,1,1,jd,id),1,wk11(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(be+1,1,1,1,jd,id),1,wk11(1,1,1,2),1)
          endif

C   ... FT p_iT_ij(12,q) and p_jT_ji(21,q); overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          tptir = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
!         if (ferm) tptiso = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            tptir = tptir + wk11(i1,i2,i3,1)
!           if (ferm) tptiso = tptiso + wkso(i1,i2,i3,1)
            Dz(i1,i2,i3,be) = Dz(i1,i2,i3,be)+dimag(wpi2*wk11(i1,i2,i3,1))
          enddo
          enddo
          enddo
!          if (ldiag) sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) + wpi2*wk11(1,1,1,1) + wpi2*wk12(1,1,1,1) +
!     .          wpi2*wk21(1,1,1,1)+ wpi2*wk22(1,1,1,1)
          tptr = tptr + tptir ! tpt = dP_RL sum_L' tpti_RL;R'L'
!         if (ferm) tptso = tptso + dble(delPi)*tptiso ! p*Im t*p*Im t
        enddo
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
         sumRLr(8+be,be,be,0,0,id) = sumRLr(8+be,be,be,0,0,id) + wpi2*tptr
!         if (ferm) sumRLr(7,0,0,id) = sumRLr(7,0,0,id) + wpi2*tptso
      enddo
C   ... End Dzyaloshinski vector block


C   ...  For J_0 ---
      if (.not. ldiag .or. .not. lchk) goto 999
      if (al == be) then
      do  id = 1, idim
          jd = id

C   ... -p_iT_ii + pT_0pT_0
C   ... wk1 = pT0, wk2 = pT0
          call dcopy(nfbz*2,ptij(1,1,1,1,id,jd),1,wk11,1)
          call zcopy(nfbz,ptij(1,1,1,1,jd,id),1,wk11(1,1,1,2),1)

C   ... FT p_iT_ij(q) and p_jT_ji(q); overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ... dP_i sum_q (T^+_RR(i,j;q) - T^-_RR(i,j;q))
          cxr = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
!            cxr = cxr + wk11(i1,i2,i3,1) + wk12(i1,i2,i3,1) + wk21(i1,i2,i3,1) +
!     .          wk22(i1,i2,i3,1)! - Scalar(i1,i2,i3,id,1,jd,1) - Scalar(i1,i2,i3,id,2,jd,2)
!            cxr(1) = cxr(1) + ptij(4,i1,i2,i3,id,jd)
!            cxr(2) = cxr(2) + ptij(2,i1,i2,i3,id,jd)
!            cxr(3) = cxr(3) + ptij(3,i1,i2,i3,id,jd)
            cxr = cxr + Scalar(1,i1,i2,i3,id,jd) + Scalar(2,i1,i2,i3,id,jd) + Scalar(3,i1,i2,i3,id,jd)! '+' in front of Scalar (there's - in formula0 is since exasa adds - itself
!??            cxr = cxr + Scalar(al,i1,i2,i3,id,jd) ! '+' in front of Scalar (there's - in formula0 is since exasa adds - itself
          enddo
          enddo
          enddo
          cxr = cxr/nfbz
!          sumRLr(2,0,0,id,al,be) = sumRLr(2,0,0,id,al,be) + wpi*cxr
!          sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) - wpi2*wk11(1,1,1,1)
          sumRLr(2,al,be,0,0,id) = sumRLr(2,al,be,0,0,id) + wpi*cxr
!          sumRLr(2,0,0,id,al,be) = sumRLr(2,0,0,id,al,be) + wpi2*ptij(2,1,1,1,id,jd)
!     .      +  wpi2*ptij(3,1,1,1,id,jd) + wpi2*ptij(4,1,1,1,id,jd)
      enddo
      endif

      enddo ! beta
      enddo ! alpha


      else !CPA

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)   !gamma_uu,p0,uu*T_uu (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(1,1,1,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(1,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,1)   !gamma_uu,p0,uu*T_uu (j)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(1,1,1,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(1,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)   !gamma_dd,p0,uu*T_uu (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(4,1,1,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(2,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,1)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(4,1,1,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(2,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,2)   !gamma_uu,p0,dd*T_dd (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(1,1,4,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(3,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(1,1,4,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(3,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,2)   !gamma_dd,p0,dd*T_dd (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(4,1,4,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(4,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(4,1,4,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(4,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)



        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)   !gamma_uu,pz,uu*T_uu (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(1,4,1,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(5,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,1)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(1,4,1,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(5,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,1)   !gamma_dd,pz,uu*T_uu (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(4,4,1,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(6,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,1)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(4,4,1,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(6,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,2)   !gamma_uu,pz,dd*T_dd (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(1,4,4,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(7,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(1,4,4,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(7,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,2)   !gamma_dd,pz,dd*T_dd (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(4,4,4,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(8,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(4,4,4,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(8,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)



        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,2)   !gamma_ud,px,ud*T_ud (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(2,2,2,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(9,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(2,2,2,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(9,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,2)   !gamma_du,px,ud*T_ud (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(3,2,2,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(10,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(3,2,2,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(10,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,1)   !gamma_up,px,du*T_du (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(2,2,3,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(11,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,1)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(2,2,3,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(11,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,1)   !gamma_du,px,du*T_du (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(3,2,3,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(12,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,1)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(3,2,3,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(12,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)



        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,2)   !gamma_ud,py,ud*T_ud (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(2,3,2,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(13,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(2,3,2,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(13,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,1,1:jdim,2)   !gamma_du,py,ud*T_ud (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(3,3,2,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(14,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,1,1:idim,2)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(3,3,2,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(14,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,1)   !gamma_up,py,du*T_du (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(2,3,3,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(15,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,1)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(2,3,3,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(15,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)

        wkij1(1:idim,1:jdim) = tij(i1,i2,i3,1:idim,2,1:jdim,1)   !gamma_du,py,du*T_du (i)
        call zgemm('N','N',idim,jdim,idim,(1d0,0d0),sfvrtxreli(3,3,3,1:idim,1:idim),idim,
     . wkij1(1:idim,1:jdim),idim,(0d0,0d0),wkij,idim)
        ptij(16,i1,i2,i3,1:idim,1:jdim) = wkij(1:idim,1:jdim)

        wkji1(1:jdim,1:idim) = tji(i1,i2,i3,1:jdim,2,1:idim,1)
        call zgemm('N','N',jdim,idim,jdim,(1d0,0d0),sfvrtxrelj(3,3,3,1:jdim,1:jdim),jdim,
     . wkji1(1:jdim,1:idim),jdim,(0d0,0d0),wkji,jdim)
        ptji(16,i1,i2,i3,1:jdim,1:idim) = wkji(1:jdim,1:idim)


      do  id = 1, idim
        tptrxx = 0
        tptryy = 0
        tptrzz = 0
!       if (ferm) tptso = 0

        do  jd = 1, jdim
          l = ll(jd)
!          call dpzero(wk11,nk1*nk2*nk3*2*2)
!          call dpzero(wk12,nk1*nk2*nk3*2*2)
!          call dpzero(wk21,nk1*nk2*nk3*2*2)
!          call dpzero(wk22,nk1*nk2*nk3*2*2)
C   ... (see Physica B 237-238 336 (1997))

C   ... p0*T0*p0*T0,  pz*Tz*pz*Tz
C   ... wk1 = pT, wk2 = pT
          call dcopy(nfbz*2,ptij(5,:,:,:,id,jd),1,wk11,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(5,:,:,:,jd,id),1,wk11(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(5,1,1,1,jd,id),1,wk11(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          call dcopy(nfbz*2,ptij(6,:,:,:,id,jd),1,wk12,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(7,:,:,:,jd,id),1,wk12(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(7,1,1,1,jd,id),1,wk12(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk12,wk12(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

           call dcopy(nfbz*2,ptij(7,:,:,:,id,jd),1,wk21,1)
           if (ldiag) then
             call zcopy(nfbz,ptij(6,:,:,:,jd,id),1,wk21(1,1,1,2),1)
           else
             call zcopy(nfbz,ptji(6,1,1,1,jd,id),1,wk21(1,1,1,2),1)
           endif
C   ...    FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
           call fftz3c(wk21,wk21(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

           call dcopy(nfbz*2,ptij(8,:,:,:,id,jd),1,wk22,1)
           if (ldiag) then
             call zcopy(nfbz,ptij(8,:,:,:,jd,id),1,wk22(1,1,1,2),1)
           else
             call zcopy(nfbz,ptji(8,1,1,1,jd,id),1,wk22(1,1,1,2),1)
           endif
C   ...    FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
           call fftz3c(wk22,wk22(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          wk00 = wk11 + wk12 + wk21 + wk22
          wkzz = wk11 - wk12 - wk21 + wk22

C   ... px*Tx*px*Tx,
C   ... wk1 = pT, wk2 = pT
          call dcopy(nfbz*2,ptij(9,:,:,:,id,jd),1,wk11,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(9,:,:,:,jd,id),1,wk11(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(9,1,1,1,jd,id),1,wk11(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          call dcopy(nfbz*2,ptij(10,:,:,:,id,jd),1,wk12,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(11,:,:,:,jd,id),1,wk12(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(11,1,1,1,jd,id),1,wk12(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk12,wk12(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

           call dcopy(nfbz*2,ptij(11,:,:,:,id,jd),1,wk21,1)
           if (ldiag) then
             call zcopy(nfbz,ptij(10,:,:,:,jd,id),1,wk21(1,1,1,2),1)
           else
             call zcopy(nfbz,ptji(10,1,1,1,jd,id),1,wk21(1,1,1,2),1)
           endif
C   ...    FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
           call fftz3c(wk21,wk21(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

           call dcopy(nfbz*2,ptij(12,:,:,:,id,jd),1,wk22,1)
           if (ldiag) then
             call zcopy(nfbz,ptij(12,:,:,:,jd,id),1,wk22(1,1,1,2),1)
           else
             call zcopy(nfbz,ptji(12,1,1,1,jd,id),1,wk22(1,1,1,2),1)
           endif
C   ...    FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
           call fftz3c(wk22,wk22(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          wkxx = wk11 + wk12 + wk21 + wk22

C   ... py*Ty*py*Ty,
C   ... wk1 = pT, wk2 = pT
          call dcopy(nfbz*2,ptij(13,:,:,:,id,jd),1,wk11,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(13,:,:,:,jd,id),1,wk11(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(13,1,1,1,jd,id),1,wk11(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          call dcopy(nfbz*2,ptij(14,:,:,:,id,jd),1,wk12,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(15,:,:,:,jd,id),1,wk12(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(15,1,1,1,jd,id),1,wk12(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk12,wk12(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

           call dcopy(nfbz*2,ptij(15,:,:,:,id,jd),1,wk21,1)
           if (ldiag) then
             call zcopy(nfbz,ptij(14,:,:,:,jd,id),1,wk21(1,1,1,2),1)
           else
             call zcopy(nfbz,ptji(14,1,1,1,jd,id),1,wk21(1,1,1,2),1)
           endif
C   ...    FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
           call fftz3c(wk21,wk21(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

           call dcopy(nfbz*2,ptij(16,:,:,:,id,jd),1,wk22,1)
           if (ldiag) then
             call zcopy(nfbz,ptij(16,:,:,:,jd,id),1,wk22(1,1,1,2),1)
           else
             call zcopy(nfbz,ptji(16,1,1,1,jd,id),1,wk22(1,1,1,2),1)
           endif
C   ...    FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
           call fftz3c(wk22,wk22(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          wkyy = - wk11 + wk12 + wk21 - wk22

C   ... Lifetime (in progress)...
!          if (ferm) then
!            call dcopy(nfbz*2,Imtij(1,1,1,id,jd),1,wkso,1)
!            if (ldiag) then
!              call zcopy(nfbz,Imtij(1,1,1,jd,id),1,wkso(1,1,1,2),1)
!            else
!              call zcopy(nfbz,Imtji(1,1,1,jd,id),1,wkso(1,1,1,2),1)
!            endif
!            call zdscal(nfbz,dble(delPj),wkso(1,1,1,2),1)

!            call fftz3c(wkso,wkso(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)
!         endif

C   ... Accumulate Jrr = J_R+T,R' = sum_LL' J_R+TL,R'L', Eq(4) in the notes
C       and sumRL(1) = dP_RL sum_TL' X_R+TL;R'L' dP_R'L', RHS of Eq(10)
C       and sumRL(3) = L-resolved version of J_00, Eq(7)
          tptirxx = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
          tptiryy = 0
          tptirzz = 0
!         if (ferm) tptiso = 0
!          call dpzero(wk12,nk1*nk2*nk3*2*2)
          wk111 = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            tptirxx = tptirxx - wkxx(i1,i2,i3,1) + wk00(i1,i2,i3,1)
            tptiryy = tptiryy - wkyy(i1,i2,i3,1) + wk00(i1,i2,i3,1)
            tptirzz = tptirzz - wkzz(i1,i2,i3,1) + wk00(i1,i2,i3,1)
            JrrR(i1,i2,i3,1,1) = JrrR(i1,i2,i3,1,1)-dimag(wpi2*wkxx(i1,i2,i3,1))+dimag(wpi2*wk00(i1,i2,i3,1))
            JrrR(i1,i2,i3,2,2) = JrrR(i1,i2,i3,1,1)-dimag(wpi2*wkyy(i1,i2,i3,1))+dimag(wpi2*wk00(i1,i2,i3,1))
            JrrR(i1,i2,i3,3,3) = JrrR(i1,i2,i3,1,1)-dimag(wpi2*wkzz(i1,i2,i3,1))+dimag(wpi2*wk00(i1,i2,i3,1))
          enddo
          enddo
          enddo
          if (ldiag) then
           sumRLr(3,1,1,0,0,id) = sumRLr(3,1,1,0,0,id) - wpi2*wkxx(1,1,1,1) + wpi2*wk00(1,1,1,1)
           sumRLr(3,2,2,0,0,id) = sumRLr(3,2,2,0,0,id) - wpi2*wkyy(1,1,1,1) + wpi2*wk00(1,1,1,1)
           sumRLr(3,3,3,0,0,id) = sumRLr(3,3,3,0,0,id) - wpi2*wkzz(1,1,1,1) + wpi2*wk00(1,1,1,1)
          endif
          tptrxx = tptrxx + tptirxx ! tpt = dP_RL sum_L' tpti_RL;R'L'
          tptryy = tptryy + tptiryy
          tptrzz = tptrzz + tptirzz
!         if (ferm) tptso = tptso + dble(delPi)*tptiso ! p*Im t*p*Im t
        enddo ! jd
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
         sumRLr(1,1,1,0,0,id) = sumRLr(1,1,1,0,0,id) + wpi2*tptrxx !/nfbz
         sumRLr(1,2,2,0,0,id) = sumRLr(1,2,2,0,0,id) + wpi2*tptryy
         sumRLr(1,3,3,0,0,id) = sumRLr(1,3,3,0,0,id) + wpi2*tptrzz
!         if (ferm) sumRLr(7,0,0,id) = sumRLr(7,0,0,id) + wpi2*tptso
      enddo ! id

C   ... For Dzyaloshinski vector

C   ... Dzyaloshinski vector D = Int Im Tr[p_i T^0_ji p_j T_ji]
      do  id = 1, idim
!        tptrxx = 0
!        tptryy = 0
        tptrzz = 0
!       if (ferm) tptso = 0
        do  jd = 1, jdim
          l = ll(jd)

C   ... wk1 = pT, wk2 = pT
          call dcopy(nfbz*2,ptij(5,:,:,:,id,jd),1,wk11,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(5,:,:,:,jd,id),1,wk11(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(5,1,1,1,jd,id),1,wk11(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          call dcopy(nfbz*2,ptij(6,:,:,:,id,jd),1,wk12,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(7,:,:,:,jd,id),1,wk12(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(7,1,1,1,jd,id),1,wk12(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk12,wk12(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          call dcopy(nfbz*2,ptij(7,:,:,:,id,jd),1,wk21,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(6,:,:,:,jd,id),1,wk21(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(6,1,1,1,jd,id),1,wk21(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk21,wk21(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          call dcopy(nfbz*2,ptij(8,:,:,:,id,jd),1,wk22,1)
          if (ldiag) then
            call zcopy(nfbz,ptij(8,:,:,:,jd,id),1,wk22(1,1,1,2),1)
          else
            call zcopy(nfbz,ptji(8,1,1,1,jd,id),1,wk22(1,1,1,2),1)
          endif
C   ... FT p_iT_ij and p_jT_ji; overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk22,wk22(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

          wkzz = wk11 - wk12 + wk21 - wk22

          tptirzz = 0  ! tpti_RL;R'L' = sum_T X_R+TL;R'L' dP_R'L'
!         if (ferm) tptiso = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
            tptirzz = tptirzz + wkzz(i1,i2,i3,1)
!           if (ferm) tptiso = tptiso + wkso(i1,i2,i3,1)
            Dz(i1,i2,i3,3) = Dz(i1,i2,i3,3)+dimag(wpi2*wkzz(i1,i2,i3,1))
          enddo
          enddo
          enddo
          tptrzz = tptrzz + tptirzz ! tpt = dP_RL sum_L' tpti_RL;R'L'
!         if (ferm) tptso = tptso + dble(delPi)*tptiso ! p*Im t*p*Im t
        enddo
C     Accumulate sum_T,R',j J_R+T,R'(i,j)  (Still resolved by R,i)
         sumRLr(11,3,3,0,0,id) = sumRLr(11,3,3,0,0,id) + wpi2*tptrzz
!         if (ferm) sumRLr(7,0,0,id) = sumRLr(7,0,0,id) + wpi2*tptso
      enddo
C   ... End Dzyaloshinski vector block


C   ...  For J_0 ---
      if (.not. ldiag .or. .not. lchk) goto 999
      if (al == be) then
      do  id = 1, idim
          jd = id

C   ... -p_iT_ii + pT_0pT_0
C   ... wk1 = pT0, wk2 = pT0
          call dcopy(nfbz*2,ptij(1,1,1,1,id,jd),1,wk11,1)
          call zcopy(nfbz,ptij(1,1,1,1,jd,id),1,wk11(1,1,1,2),1)

C   ... FT p_iT_ij(q) and p_jT_ji(q); overwrite wk1 with product of result (=> R.S.)
          call fftz3c(wk11,wk11(1,1,1,2),nk1,nk2,nk3,k1,k2,k3,30,-1)

C   ... dP_i sum_q (T^+_RR(i,j;q) - T^-_RR(i,j;q))
          cxr = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1
!            cxr = cxr + wk11(i1,i2,i3,1) + wk12(i1,i2,i3,1) + wk21(i1,i2,i3,1) +
!     .          wk22(i1,i2,i3,1)! - Scalar(i1,i2,i3,id,1,jd,1) - Scalar(i1,i2,i3,id,2,jd,2)
!            cxr(1) = cxr(1) + ptij(4,i1,i2,i3,id,jd)
!            cxr(2) = cxr(2) + ptij(2,i1,i2,i3,id,jd)
!            cxr(3) = cxr(3) + ptij(3,i1,i2,i3,id,jd)
            cxr = cxr + Scalar(1,i1,i2,i3,id,jd) + Scalar(2,i1,i2,i3,id,jd) + Scalar(3,i1,i2,i3,id,jd)! '+' in front of Scalar (there's - in formula0 is since exasa adds - itself
!??            cxr = cxr + Scalar(al,i1,i2,i3,id,jd) ! '+' in front of Scalar (there's - in formula0 is since exasa adds - itself
          enddo
          enddo
          enddo
          cxr = cxr/nfbz
!          sumRLr(2,0,0,id,al,be) = sumRLr(2,0,0,id,al,be) + wpi*cxr
!          sumRLr(3,0,0,id,al,be) = sumRLr(3,0,0,id,al,be) - wpi2*wk11(1,1,1,1)
          sumRLr(2,al,be,0,0,id) = sumRLr(2,al,be,0,0,id) + wpi*cxr
!          sumRLr(2,0,0,id,al,be) = sumRLr(2,0,0,id,al,be) + wpi2*ptij(2,1,1,1,id,jd)
!     .      +  wpi2*ptij(3,1,1,1,id,jd) + wpi2*ptij(4,1,1,1,id,jd)
      enddo
      endif


      endif ! CPA


      Jrrcp(:,:,:,:,:) = JrrR(:,:,:,:,:)
      sumRLcpz(:,:,:,:,:,:) = sumRLr(:,:,:,:,:,:)
      deallocate(tij0,tji0,tijx,tjix,tijy,tjiy,tijz,tjiz,ptij,ptji,pii,pjj,Scalar)!,pttotij,pttotji)
!        if (ferm) deallocate(Imtij,Imtji)
!        if (ferm .and. .false.) deallocate(Imtij,Imtji)

  999 call tcx('pasajjnoncol')
      end
