      subroutine locpot(nbas,nsp,nspc,lso,lbf,lcplxp,v0beta,s_pot,s_ham,s_lat,
     .  s_optic,s_site,s_spec,s_rhat,s_atparms,qmom,vval,gpot0,rhobg,nlibu,
     .  lmaxu,vorb,lldau,nbf,bsite,ppnl,hab,vab,sab,qval,qsc,eterms,job)
C- Make the potential at the atomic sites and augmentation matrices.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  shfac bxcscali socscl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  eula lrsa
Co     Stored:     eula
Co     Allocated:  *
Cio    Elts passed:lrsa eula
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg jcg indxcg
Cio    Passed to:  *
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic nlg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:kmxax rgrad rgrade
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu v0 v1 pz
Co     Stored:     saxis
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z qc rg name a nr rmt rsma lmxa lmxl lmxb idu kmxt
Ci                 mxcst coreh stc lfoca rfoca ctail etail p pz coreq
Ci                 orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  corprm gtpcor uspecb
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rhoc rho2
Cio    Passed to:  *
Cio s_atparms
Ci    Elts read:  vvesat cpnves rhovv1 rhovv2 rhvxt1 rhvxt2 qloc aloc qlocc rhoexc rhoex rhoec rhovxc
Ci                focexc focex focec focvxc bfam rvepsv rvexv rvecv rvvxcv rveps rvvxc rhcv1 aloc rhcv1
Ci    Stored:     vvesat cpnves rhovv1 rhovv2 rhvxt1 rhvxt2 qloc aloc qlocc rhoexc rhoex rhoec rhovxc
Co                focexc focex focec focvxc bfam rvepsv rvexv rvecv rvvxcv rveps rvvxc rhcv1 aloc
Co    Allocated:
Cio   Elts passed:
Cio   Passed to:
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   lso   :0 no SO coupling
Ci         :1,3,4 Add L.S coupling to potential ppi
Ci         :2 Add LzSz only to potential ppi
Ci         :10s digit nonzero: print the SO weighted average of the electric field
Ci   lbf   :0 no external B field
Ci         :1 Add B-Sigma to potential ppi
Ci         :2 Add BzSigz only to potential ppi
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   v0beta:mixing beta for admixing sphere potential V0 defining phi,phidot
Ci         :Used only iif lfltwf=T
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   vval  :electrostatic potential at MT boundary; needed
Ci         :to compute matrix elements of local orbitals.
Ci   gpot0 :integrals of local gaussians * phi0~
Ci         :phi0~ is the estatic potential of the interstitial
Ci   job   :1s  digit
Ci         : 0 Indended for output density
Ci         :   Do not make core and augmentation matrices;
Ci         :   Do not change eterms(8,9,12)
Ci         : 1 make core and augmentation matrices
Ci         :   Also make eterms(8,9,12) needed for HF total energy
Ci         :10s digit
Ci         : 0 update potential used to define basis,
Ci         :   provided 1s digit job is also set
Ci         : 1 do not update potential
Ci         :   NB: Caller can also suppress update for specific species
Ci         :   through 4's bit of species->mxcst.
Ci         :100s digit
Ci         : 1 exclude exchange-correlation potential in making
Ci         :   matrix elements of the potential.
Ci         :   Also implies floating potential.
Ci         : 2 Use spin-averaged rho when making potential => Bxc=0
Ci         :1000s digit
Ci         :1 Generate valence-only rvvxcv,rvepsv,rvexv,rvecv, and rveps,rvvxc
Ci         :10000s digit
Ci         :1 write sphere density for site ib to file rhoMT.ib
Ci         :100000s digit
Ci         : 0 use density as is.
Ci         : 1 reset points of negative density to positive density
Ci             for purposes of calculating vxc
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci         : 3 like 1, but modifed rho is retained.
Ci  rhobg  :compensating background density
Ci  nlibu  : max number of lda+u blocks
Ci  lmaxu  : max l for U
Ci  vorb   : orbital dependent potential matrices
Ci  lldau  :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat beginning at dmats(*,lldau(ib))
Ci  nbf    :dimensions bsite
Ci  bsite  :external magnetic field for this site
Co Outputs
Co   eterms:quantities related to the core total energy
Co         :(8)  sumec = sum of foca=0 core eigenvalues
Co         :(9)  sumtc = sum of all core kinetic energies
Co         :(10) rhcvef1  = integral rhoc*(v1-2Z/r) for sites w/ lfoca=0
Co         :(12) sumt0 = sum of frozen core kinetic energies
Co   qval  :nominal total valence charge qv-z
Co   qsc   :semicore charge from local orbitals
Co   sig   :augmentation overlap integrals; see augmat
Co   tau   :augmentation kinetic energy integrals; see augmat
Co   ppi   :augmentation kinetic + potential integrals; see augmat
Co   ppnl  :NMTO-like potential parameters
Co   hab   :integrals of the ham. with true w.f.  See Remarks in augmat
Co   vab   :integrals of the pot. with true w.f.  See Remarks in augmat
Co   sab   :integrals of    unity with true w.f.  See Remarks in augmat
Co   ... Elements of s_atparms summed over all spheres:
Co    qv1    : true valence charge
Co    qv2    : sm   valence charge
Co    qcor1  : true core charge
Co    rhcv1  : int rhoc*(v1-2Z/r) for sites w/ lfoca=0
Co    qloc   : total valence charge rho1-rho2 in sphere
Co    qlocc  : total core charge in sphere
Co    atrue  : atrue(0) Magnetic moment of true density; atrue(1:3) direction
Co    aloc   : total valence magnetic moment in sphere+core moment for sites with core hole
Co    alocc  : core magnetic moment
Co
Co    vales1 : int n1_val ves1,  ves1~ = ves[n1] - 2*Z/r
Co    vales2 : int n2_val ves2,  ves2~ = ves[n2+gval+gcor+gnuc]
Co    vvesat : vales1 - vales2
Co           : Total electrostatic energy has a term vvesat/2
Co    cpnves : integral of core+nucleus times electrostatic potential
Co           : = int [ (rhoc-z) ves1~ - (rhocsm-rhonsm) ves2~ ]
Co           : where ves1~ is the estat potential of the true density
Co           : where ves2~ is the estat potential of the smooth density
Co           " =vales1 - vales2. Total electrostatic energy has a term cpnves/2
Co
Co    rep1   : int n1 * exc[n1], n1 = true density
Co    rep2   : int n2 * exc[n2], n2 = sm density (may be pert corr)
Co    rhoexc : rep1 - rep2
Co    rhoex  : exchange-only part of rhoexc
Co    rhoec  : correlation-only part of rhoexc
Co    rmu1   : int n1 * vxc[n1], n1 = true density
Co    rmu2   : int n2 * vxc[n2], n2 = sm density (may be pert corr)
Co    rhovxc : rmu1 - rmu2
Co    focvxc : (lfoc=2) int gcor and exc[n2], n2 without cores
Co    focexc : (lfoc=2) int gcor and (dvxc/dn2 * n2)
Co     focex : (lfoc=2) Exchange part of focexc
Co     focec : (lfoc=2) Correlation part of focexc
Co
Co    rhovv1 : int n1_val * v1~   where v1~ = ves[n1] + vxc[n1] + vext FIX
Co    rhovv2 : int n2_val * v2~ + sum_L qmom_L gpotb_L, v2~ = ves[n2] + vxc[n2] + vext FIX
Co
Co    rvepsv : (1000s digit = 1) like rhoeps, but for valence density only
Co    rvexv  : (1000s digit = 1) exchange    part of rvepsv
Co    rvecv  : (1000s digit = 1) correlation part of rvepsv
Co    rvvxcv : (1000s digit = 1) like rhovxc, but for valence density only
Co    rveps  : (1000s digit = 1) int rhov*exc(rhotot)
Co    rvvxc  : (1000s digit = 1) int rhov*vxc(rhotot)
Co
Co    rvext  : rvext(2,2) = int n(1,2)_(val,cor) * vext
Co    rhvxt1 : int n1_val * vext
Co    rhvxt2 : int n1_val * vext
Co    sgpote : sum_ilm qmom_ilm * int [gval * vext]_ilm
Co    bfam   : sum (bext * local moment)/2 --- for double counting.
Cl Local variables
Cl   lfltwf:T  update potential used to define basis
Cl   iblu  :index to current LDA+U block
Cl   idu   :idu(l+1)>0 => this l has a nonlocal U matrix
Cl   facso :scale L.S  matrix elements by facso. Should normally be 1.
Cl         :facso determined from shfac(2,ib), where ib = site index
Cl   bscal :scale Bxc field from true density by bscal(1), taken from shfac(2,ib)
Cl         :scale Bxc field from smooth density by bscal(2), taken from bxcscli
Cl         :Both should normally be 1.
Cl   nlg   :dimensions s_optic%rgrad
Cl   ns4 = 1 if local density is not spin polarized
Cl         2 if local spin density is collinear along z
Cl         4 if local spin density is to be rotated;
Cl   lrsa  :(from s_ham%lrsa)
Cl         1s digit governs how noncollinearity is treated.
Cl         0 Input density collinear along z.
Cl           Noncollinear output density generated to calculate energy only.
Cl           For purposes of this routine the density is collinear (nspc should be 1)
Cl         1 Rigid spin approximation.  Local input densities are
Cl           constrained to align with given Euler angles.
Cl           Output density treated as case 0, except spin is rigidly
Cl           rotated to axis given by site Euler angles.
Cl         2 Rigid spin approximation.  Like 1, except spin axis determined
Cl           by average magnetization.
Cl         3 Fully noncollinear density.
Cl           A local axis still must be defined to make the basis wave functions.
Cl           We use the same choice as 2.
Cl        10s digit determines how local z axis defining basis function is set.
Cl         0 Site spin quantization axis defining local density may be parallel to M or -M
Cl           whichever yields the largest projection.
Cl           Local moments can be positive or negative.
Cl         1 Site spin quantization axis defining local density is parallel to M
Cl           Local moments are positive by construction.
Cl           Note: there can be a difference between the two because
Cl           the core density is aligned with the quantization axis.
Cl           Small noncollinearity => use 0 to adiabatically move from collinear case.
Cl           Otherwise use 1 avoid discontinuous jumps with small changes in input conditions.
Cl           In this case starting magnetic moments should all be positive to avoid
Cl           ambiguities in the orientation (+M, -M) in the local rho.
Cl   gpotb :gpotb(ilm) = integral [gaussians(r,radius rg)_ilm * smooth ves
Cl                     = local analog of gpot0 generated by smves.
Cl   vnucl :true estat potential + vext at nucleus
Cl   buf   :Buffer to contain redirected standard output. :It is used if buffer_output is .true.
Cl         :In preparation for MPI parallelization of this routine.
Cl         :locpot and routines that call it write to stdout through calls to awrite,
Cl         :which has a special 'buffer' mode.
Cr Remarks
Cu Updates
Cu   29 May 17 Skips electric field gradient for sites with lmxa<2
Cu   17 Jan 16 option to render LDA Bxc = 0 (:100s digit job = 2)
Cu   12 Jul 15 New v0beta (modified argument list)
Cu   22 Jun 14 pass nspc to augmat (generates noncollinear optical matrix elmements)
Cu   13 Jan 14 (Ben Kaube) start optical matrix elements for PkL
Cu   15 Nov 13 (Ben Kaube) first cut at optical matrix elements
Cu             (MvS) first cut at the the noncollinear case
Cu   08 Sep 13 Returns parameters for core total energy in eterms
Cu   24 Jun 13 New facso, input via shfac, scaling L.S matrix elements
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   11 Mar 11 local densities can be reset to be positive definite
Cu   04 Jul 10 When writing density to file rhoMT, append core density
Cu   01 Apr 10 Add external B field to potential.  New argument list
Cu   10 Sep 08 Added electric field gradient (AxSv)
Cu   02 Jan 06 adds core magnetic moment to saloc
Cu   09 Nov 05 Convert dmat to complex form
Cu   06 Jul 05 Parameters for extended local orbitals are now input
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   12 Jun 05 Optionally write sphere densities to file
Cu   27 Apr 05 Added LDA+U  (Lambrecht)
Cu    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
Cu   15 Jul 04 First implementation of extended local orbitals
Cu   14 Mar 04 Makes rhov*exc and rhov*vxc
Cu   14 Mar 04 Makes rhov*exc and rhov*vxc
Cu   19 Sep 03 (ATP) Enabled partial core occupation
Cu   02 Oct 02 (WRL) Added background potential
Cu    9 May 02 Added species-specific freezing of potential
Cu    8 May 02 Added rhoex and rhoec (T. Miyake)
Cu    7 May 02 Added rvexv and rvecv (T. Miyake)
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   17 Sep 01 Returns qsc.  Altered argument list.
Cu   28 Aug 01 Extended to local orbitals.  Altered argument list.
Cu   15 Aug 01 Generates rvepsv and rvvxcv.  New argument list
Cu   18 Apr 01 Added ability to exclude exchange-correlation potential
Cu   20 Feb 01 Added ppnl to potential parameters generated
Cu   13 Jun 00 spin polarized
Cu    1 May 00 Adapted from nfp locpot.f
Cu   17 Jun 98 Adapted from DLN to run parallel on SGI
C ----------------------------------------------------------------------
      use structures
      use mpi
      use mod_adstrb
      implicit none
C ... Passed parameters
      integer nbas,nsp,nspc,job,n0,nppn,nab,lbf,lso,lcplxp
      integer nlibu,lmaxu,lldau(nbas),iblu,nbf
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      parameter (n0=10,nppn=12,nab=9)
      double precision qval,rhobg,v0beta,qmom(1),vval(1),
     .  hab(nab,n0,nsp,nbas),vab(nab,n0,nsp,nbas),sab(nab,n0,nsp,nbas),
     .  gpot0(*),ppnl(nppn,n0,nsp,nbas),bsite(nbas,nbf,3),eterms(22)

C ... For structures
!      include 'structures.h'
      type(str_pot)::   s_pot
      type(str_ham) ::  s_ham
      type(str_lat) ::  s_lat
      type(str_optic):: s_optic
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_rhat)::  s_rhat(*)
      type(str_atparms) ::  s_atparms,s_atparmsa
C ... Dynamically allocated local arrays
      real(8),allocatable:: rhol1(:),rhol2(:),rhoc(:,:),
     .  v1(:),v2(:),v1es(:),v2es(:),efg(:,:),zz(:),vext(:,:,:)
      real(8), pointer :: p_v0(:,:),p_v1(:,:)
      character, allocatable:: buf(:)*36000
C ... Local parameters
      character spid*8
      integer nrmx,nlmx,nkap0
      parameter (nrmx=5001, nlmx=81, nkap0=4)
      integer lh(nkap0),nkapi,nkape,nkaph,k,lx
      double precision eh(n0,nkap0),rsmh(n0,nkap0),ehl(n0),rsml(n0)
      logical lfltwf
      logical :: buffer_output = .false.
      integer i,ib,intopt,ipr,is,j1,jobl,kmax,lfoc,lmxa,lmxb,lmxl,
     .  loptic,ncore,nlg,nrg,nlml,nr,nrml,ns4,stdo,iv(10)
      double precision a,ceh,cofg,cofh,facso,pi,qc,qcorg,qcorh,qsc,qsca,qv,rfoc,rg,
     .  rmt,rsma,srfpi,xx,y0,z,alpha,beta,bscal(2)
      double precision rofi(nrmx),rwgt(nrmx),pnu(n0,2),bloc(nbf,3),pnz(n0,2),gpotb(81)
C     double complex magsp(2,2),umagsp(2,2)
C ... for sm. Hankel tails
      integer iltab(nbas)
      double precision rs3,vmtz
C ... for LDA+U
      integer idu(4)
C ... for core energy
      double precision smec,smtc,sumtc,sumec,stc0,sumt0
C ... for core hole
C     character chole*8
      character strn*9
      integer kcor,lcor
      double precision qcor(2),qc0,qsc0
      double precision rotl(3,3),saxis(3),eula(3)
      procedure(integer) :: iprint,isw,nglob,mpipid
      procedure(real(8)) :: dlength

      integer :: iproc, nproc, comm, err, nlms(0:nbas), iblus(0:nbas), l
      integer :: ab_mt, pn_mt, rg_mt, re_mt, tmp_mt
      integer, allocatable :: pmap(:), smap(:)
      real(8), allocatable :: tmp(:,:)
      type(pl_t) :: pl(0:nbas-1)

C --- Setup ---
      call tcn('locpot')
      stdo = nglob('stdo')
      ipr = iprint()
      allocate(buf(1)); buf(1) = ' '
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0 = 1d0/srfpi
      nkaph = nglob('nkaph')
      ns4 = nsp*nspc; if (mod(s_ham%lrsa,10) == 0) ns4 = nsp ! See "Local variables"
      k = nrmx*nlmx*2*ns4
      allocate(rhol1(k),rhol2(k),v1(k),v2(k),v1es(k),v2es(k))
      allocate(efg(5,nbas),zz(nbas))
      zz = 0
      if (nbf > 0) call dpzero(bloc,nbf*3)

      call info0(30,1,0,' locpot:')

C --- Start loop over sites ---
      nlms(0) = 0
      do ib = 1, nbas
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
!         if (lmxl == -1 .or. s_spec(is)%lmxa == -1) lmxl = 0
        nlml = (lmxl+1)**2
        nlms(ib) = nlms(ib-1) + nlml
      end do

      iblus(0) = 0
      do ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        idu = s_spec(is)%idu
        iblus(ib) = iblus(ib-1)
        do  l = 0, min(lmxa,3)
          if (mod(idu(l+1),10) /= 0) iblus(ib) = iblus(ib)+1
        end do
      end do

      comm = mpi_comm_world
      call mpi_comm_size(comm, nproc, err)
      call mpi_comm_rank(comm, iproc, err)
      allocate(pmap(0:nproc))
      call vbdist(nbas, nlms(1:nbas), nproc, pmap)

      buffer_output = nproc > 1
      if (buffer_output) then   ! Every node should print
        ipr = min(ipr,50)       ! Limit verbosity
        call mpibc1(ipr,1,2,.false.,'','')
      endif

      call iinit(iltab,nbas)
      qval = 0d0; qsc = 0d0; sumtc = 0d0; sumec = 0d0; sumt0 = 0d0

      do  ib = pmap(iproc)+1, pmap(iproc+1)
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
        p_v0 => s_site(ib)%v0
        p_v1 => s_site(ib)%v1
        pnz = s_site(ib)%pz
        z = s_spec(is)%z
        qc = s_spec(is)%qc
        rg = s_spec(is)%rg
        spid = s_spec(is)%name
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        rsma = s_spec(is)%rsma
        lmxa = s_spec(is)%lmxa
        lmxl = s_spec(is)%lmxl
        lmxb = s_spec(is)%lmxb
        zz(ib) = z
        if (lmxa == -1) cycle

        if (buffer_output) then
          call pshpr(ipr)
C         print *, 'setup buffer_output',iproc
          call awrit0('%&',' ',0,0)
        endif

        bscal(1) = 1+s_pot%shfac(1,ib)
        bscal(2) = s_pot%bxcscali
C       facso    = 1+s_pot%shfac(2,ib)  replace s_pot%shfac(2,ib) -> s_pot%socscl May 2015
        facso    = s_pot%socscl(ib)

        if (lbf > 0) call dcopy(nbf*3,bsite(ib,1,1),nbas,bloc,1)
        idu = s_spec(is)%idu
        kmax = s_spec(is)%kmxt
        i = s_spec(is)%mxcst
C       Float wave functions if:
C       100's digit job=1   OR:
C         4's bit mxcst=0   AND   10's digit job=0  AND  1's digit job=1
        lfltwf = mod(job/100,10) == 1 .or.
     .  mod(i/4,2) == 0 .and. mod(job/10,10) == 0 .and. mod(job,10) == 1

        call corprm(s_spec,is,qcorg,qcorh,qsca,cofg,cofh,ceh,lfoc,rfoc,z)
C       chole = s_spec(is)%coreh
        call gtpcor(s_spec,is,kcor,lcor,qcor)
        call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc0,qv,qsc0)
        if (qsc0 /= qsca .or. qc /= qc0-qsc0) then
          if (iprint() > 0)
     .    call awrit5(' is=%i qsc0=%d qsca=%d qc=%d qc0=%d',buf,120,
     .                stdo,is,qsc0,qsca,qc,qc0)
          call rxs('problem in locpot -- possibly low LMXA'//
     .      ', or orbital mismatch, species ',spid)
        endif
        qval = qval+qv+qsca
        qsc  = qsc+qsca

        nlml = (lmxl+1)**2
        nrml = nr*nlml
        j1 = nlms(ib-1)+1
        if (ipr >= 20) then
           strn = '%,3i   a='; if (nr>999) strn = '%,4i  a='
           call awrit8('%N site%,3i  z=%;5,1D  rmt=%;8,5D  nr='//trim(strn)//
     .       '%;5,3D  nlml=%,2i  rg=%;5,3D  Vfloat=%l',
     .       buf,100,stdo,ib,z,rmt,nr,a,nlml,rg,lfltwf)
           if (kcor /= 0) then
           if (qcor(1) /= 0 .or. qcor(2) /= 0) then
             if (ipr >= 30)
     .       call awrit5(' core hole:  kcor=%i  lcor=%i  qcor=%d  amom=%d',buf,100,stdo,kcor,lcor,qcor,qcor(2),0)
           endif
           endif
        endif
        call rxx(nr > nrmx,  'locpot: increase nrmx')
        call rxx(nlml > nlmx,'locpot: increase nlmx')
        call radmsh(rmt,a,nr,rofi)
        intopt = 10*nglob('lrquad')
        call radwgt(intopt,rmt,a,nr,rwgt)

C   ... Write true density to file rhoMT.ib
        if (mod(job/10000,10) == 1 .or. mod(job/10000,10) == 3) then
          call wrhomt(3,'rhoMT.','density',ib,s_rhat(ib)%rho1,
     .      s_rhat(ib)%rhoc,rofi,nr,nlml,nsp)
        endif

C   ... One-center expansion of s_pot%smvextc
        allocate(vext(nr,nlml,nsp))
        if (associated(s_pot%smvextc)) then
          call msh21c(ib,s_site,s_spec,s_lat,s_pot%smvextc,size(s_pot%smvextc),vext,vext,nr)
          jobl = 1000000
        else
          jobl = 0
          call dpzero(vext,size(vext))
        endif

C   --- Make potential and energy terms at this site ---
        call locpt2(job+jobl,z,rmt,rg,a,nr,nsp,ns4/nsp,cofg,cofh,ceh,rfoc,lfoc,
     .    nlml,qmom(j1),vval(j1),rofi,rwgt,s_rhat(ib)%rho1,s_rhat(ib)%rho2,
     .    s_rhat(ib)%rhoc,nbf,bscal,bloc,rhobg,vext,rhol1,rhol2,v1,v2,
     .    v1es,v2es,s_atparmsa,gpotb,efg(1,ib),buf)

        deallocate(vext)

C   ... Write true potential to file vtrue.ib
        if (mod(job/10000,10) == 2 .or. mod(job/10000,10) == 3) then
          call wrhomt(1,'vtrue.','potential',ib,v1,v1,rofi,nr,nlml,nsp)
        endif

C   ... Noncollinear local density
C       Define a spin quantization axis
        if (ns4 == 4) then  !Noncollinear case

C       Set local spin quantization axis according to the following
C       Let lten = 10s digit lrsa and l1 = 1s digit lrsa.
C       Let e = vector corresponding to given Euler angles E
C       Let M = Magnetization axis
C       Let S = Spin quantization axis
C         lten \      l1=1                            l1=2,3
C              -----------------------------------------------
C              |     Obtain S from e                  S = M/|M| or -M/|M|
C              |     M determined by                  lten determines which
C              |     alignment with e                 is chosen (see rotspv)
C              |
C          0   |     Use E as given                   S = M/|M|
C              |     S not determined by M
C              |
C          1   |     Use E as given                   S = +/- M/|M|
C              |     S not determined by M         => Mz>0 along spin axis
C              -------------------------------------------
          eula = s_ham%eula(ib,1:3) ! Initial value, subject to update
          select case (mod(s_ham%lrsa,10))
            case (1)
              call eua2rm(eula(1),eula(2),eula(3),rotl)
              saxis = rotl(3,:)

            case (2,3)
              saxis = s_atparmsa%atrue(1:3)/dlength(3,s_atparmsa%atrue(1),1)
              eula(1) = datan2(saxis(2),saxis(1))
              eula(2) = datan2(dsqrt(saxis(1)**2+saxis(2)**2),
     .                         saxis(3))
            case default
              call rx('locpot: illegal value of lrsa')

          end select

          alpha = datan2(s_atparmsa%atrue(2),s_atparmsa%atrue(1))
          beta = datan2(dsqrt(s_atparmsa%atrue(1)**2+s_atparmsa%atrue(2)**2),s_atparmsa%atrue(3))
          call info8(30,0,0,' Given alp,bet: %,6d %,6d'//
     .      '   from <M>: %,6d %,6d  Use: %,6d %,6d'//
     .      '%?#n# (*)##',
     .      s_ham%eula(ib,1),s_ham%eula(ib,2),alpha,beta,
     .      eula(1),eula(2),
     .      isw(mod(s_ham%lrsa,100) >= 12 .and. s_atparmsa%atrue(3) < 0),0)


C         Debugging
C         print *, 'EI',sngl(s_ham%eula(ib,1:3))
C         call eua2rm(s_ham%eula(ib,1),s_ham%eula(ib,2),0d0,rotl)
C         print *, 'VI',sngl(rotl(3,:))
C
C         print *, 'EM',sngl(alpha),sngl(beta)
C         call eua2rm(alpha,beta,0d0,rotl)
C         print *, 'VM',sngl(s_atparmsa%atrue(1:3)/dlength(3,s_atparmsa%atrue(1),1))
C
C         print *, "E'",sngl(eula)
C         call eua2rm(alpha,beta,0d0,rotl)
C         print *, "V'",sngl(saxis)

C         Update Euler angles in case they changed
          s_ham%eula(ib,1:3) = eula

C         Rotate potential to local quantization axis
          i = 1000*mod(s_ham%lrsa/10,10) + 544
C         Debugging
C          call rotspv(i,saxis,nr*nlml,nr*nlml,nr*nlml,
C     .      s_site(ib)%rho1,xx,s_site(ib)%rho1,xx)

C         call prmx('v1 before',v1,nr*nlml,nr,4)
          call rotspv(i,saxis,nr*nlml,nr*nlml,nr*nlml,v1,xx,v1,xx)
C         call prmx('v1 after',v1,nr*nlml,nr,4)

          s_site(ib)%saxis = saxis

        endif

C   ... Update the potential used to define basis set
        if (lfltwf) then
          call dscal(nr*nsp,1-v0beta,p_v0,1)
          do  i = 0, nsp-1
C           call dpscop(v1,p_v0,nr,1+nr*nlml*i,1+nr*i,y0)
            call dpsadd(p_v0,v1,nr,1+nr*i,1+nr*nlml*i,y0*v0beta)
          enddo
        endif
C       call prrmsh('v0',rofi,p_v0,nr,nr,nsp)

C   ... Store the potential used in mkrout to calculate the core
        do  i = 0, nsp-1
          call dpscop(v1,p_v1,nr,1+nr*nlml*i,1+nr*i,y0)
        enddo

C   ... Accumulate terms for LDA total energy
        s_atparms%vvesat = s_atparms%vvesat + s_atparmsa%vvesat
        s_atparms%cpnves = s_atparms%cpnves + s_atparmsa%cpnves
        s_atparms%rhovv1 = s_atparms%rhovv1 + s_atparmsa%rhovv1
        s_atparms%rhovv2 = s_atparms%rhovv2 + s_atparmsa%rhovv2
        s_atparms%rhvxt1 = s_atparms%rhvxt1 + s_atparmsa%rhvxt1
        s_atparms%rhvxt2 = s_atparms%rhvxt2 + s_atparmsa%rhvxt2
        s_atparms%qloc   = s_atparms%qloc  + s_atparmsa%qloc
        s_atparms%aloc   = s_atparms%aloc  + s_atparmsa%aloc
        s_atparms%qcor1  = s_atparms%qcor1 + s_atparmsa%qcor1
        s_atparms%qlocc  = s_atparms%qlocc + s_atparmsa%qlocc
        s_atparms%rhoexc = s_atparms%rhoexc + s_atparmsa%rhoexc
        s_atparms%rep1   = s_atparms%rep1  + s_atparmsa%rep1
        s_atparms%rep2   = s_atparms%rep2  + s_atparmsa%rep2
        s_atparms%rmu1   = s_atparms%rmu1  + s_atparmsa%rmu1
        s_atparms%rmu2   = s_atparms%rmu2  + s_atparmsa%rmu2
        s_atparms%vales1 = s_atparms%vales1  + s_atparmsa%vales1
        s_atparms%vales2 = s_atparms%vales2  + s_atparmsa%vales2
        s_atparms%qv1    = s_atparms%qv1  + s_atparmsa%qv1
        s_atparms%qv2    = s_atparms%qv2  + s_atparmsa%qv2
        s_atparms%rhoex  = s_atparms%rhoex  + s_atparmsa%rhoex
        s_atparms%rhoec  = s_atparms%rhoec  + s_atparmsa%rhoec
        s_atparms%rhovxc = s_atparms%rhovxc + s_atparmsa%rhovxc
        s_atparms%focexc = s_atparms%focexc + s_atparmsa%focexc
        s_atparms%focex  = s_atparms%focex  + s_atparmsa%focex
        s_atparms%focec  = s_atparms%focec  + s_atparmsa%focec
        s_atparms%focvxc = s_atparms%focvxc + s_atparmsa%focvxc
        s_atparms%rvext  = s_atparms%rvext + s_atparmsa%rvext
        s_atparms%atrue  = s_atparms%atrue + s_atparmsa%atrue
        s_atparms%bfam   = s_atparms%bfam + s_atparmsa%bfam
        s_atparms%sgpote = s_atparms%sgpote + s_atparmsa%sgpote
        if (mod(job/1000,10) == 1) then
          s_atparms%rvepsv = s_atparms%rvepsv + s_atparmsa%rvepsv
          s_atparms%rvexv  = s_atparms%rvexv  + s_atparmsa%rvexv
          s_atparms%rvecv  = s_atparms%rvecv  + s_atparmsa%rvecv
          s_atparms%rvvxcv = s_atparms%rvvxcv + s_atparmsa%rvvxcv
          s_atparms%rveps  = s_atparms%rveps  + s_atparmsa%rveps
          s_atparms%rvvxc  = s_atparms%rvvxc  + s_atparmsa%rvvxc
        endif
        if (lfoc == 0) then
          s_atparms%rhcvef1 = s_atparms%rhcvef1 + s_atparmsa%rhcvef1
          s_atparms%rhcvxt1 = s_atparms%rhcvxt1 + s_atparmsa%rvext(1,2)
        endif

C       Check for core moment mismatch ; add to total moment
        if (kcor /= 0) then
          if (dabs(qcor(2)-s_atparmsa%alocc) > 0.01d0) then
            if (iprint() >= 10) call awrit5(' (warning) core moment mismatch spec %i:'
     .        //'  input file=%;6d  core density=%;6;4d',
     .      buf,120,stdo,is,qcor(2),s_atparmsa%alocc,0,0)
          endif
          s_atparms%aloc(0) = s_atparms%aloc(0) + qcor(2)
        endif

C   --- Make augmentation matrices sig, tau, ppi ---
        if (mod(job,10) == 1) then

C     ... Smooth Hankel tails for local orbitals
          call dpzero(rsmh,n0*nkap0)
          call dpzero(eh,n0*nkap0)
          call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
          call uspecb(0,2,s_spec,is,is,lh,rsmh,eh,nkape)
          if (nkape > nkapi) then
            call dcopy(n0,rsmh(1,nkaph),1,rsml,1)
            call dcopy(n0,eh(1,nkaph),  1,ehl, 1)
          endif

C     ... Use effective potentials with modified xc
          if (mod(job/100,10) == 1) then
            call dcopy(nrml*nsp,v1es,1,v1,1)
            call dcopy(nrml*nsp,v2es,1,v2,1)
          endif

C          if (ib == 1) then
C          call prrmsh('v1 loc',rofi,v1,nr,nr,nlml*nsp)
C          call prrmsh('v2 loc',rofi,v2,nr,nr,nlml*nsp)
C          endif

          if (ipr >= 20) then
C            write(stdo,467) y0*(gpot0(j1)-gpotb(1))
C  467       format(' potential shift to crystal energy zero:',f12.6)
C            if (nsp == 2 .and. bscal(1) /= 1) then
C              call info2(20,0,0,' scale Bxc by %g',bscal,0)
C            endif
            call awrit1(' potential shift to crystal energy zero:%;12,6D',
     .        buf,120,stdo,y0*(gpot0(j1)-gpotb(1)))
            if (nsp == 2 .and. bscal(1) /= 1) then
              call awrit1(' scale Bxc by %g',buf,120,stdo,bscal)
            endif
          endif

          lx = 0
          if (mod(job/100,10) == 1) lx = 1
          loptic = s_optic%loptic
          if (loptic > 0) then
            i = ib; nlg = s_optic%nlg
            iv(1:5) = shape(s_optic%rgrad)
            nrg = nint(sqrt(dble(iv(1))))
            if (nrg*nrg /= iv(1)) call rx('locpot: s_optic%rgrad has odd dimension')
          else
            i = 1 ; nlg = 1; nrg = 1
          endif
          iblu = iblus(ib-1)
          call augmat(s_site(ib),z,rmt,rsma,lmxa,pnu,pnz,s_optic%kmxax,
     .      kmax,nlml,a,nr,nsp,nspc,lx,lso,lbf,loptic,nrg,nlg,facso,rofi,
     .      rwgt,s_lat%cg,s_lat%jcg,s_lat%indxcg,p_v0,v1,v2,gpotb,
     .      gpot0(j1),nkaph,nkapi,lmxb,lh,eh,rsmh,ehl,rsml,rs3,vmtz,
     .      lcplxp,lmaxu,vorb,lldau(ib),iblu,idu,nbf,bloc,
     .      ppnl(1,1,1,ib),hab(1,1,1,ib),vab(1,1,1,ib),sab(1,1,1,ib),
     .      s_optic%rgrad(1,1,1,1,i),s_optic%rgrade(1,1,1,1,i),buf)
        endif

C   ... Get sumec
        if (lfoc == 0) then
          allocate(rhoc(nr,nsp))
          call pshpr(ipr-11)
C         j1 = awrite('%-2&',' ',0,0,xx,xx,xx,xx,xx,xx,xx,xx) ! Returns value of buf
          if (buffer_output) call setpr(1) ! Suppress output from getcor in buffer mode
          call getcor(0,z,a,pnu,pnz,nr,lmxa,rofi,p_v1,kcor,lcor,
     .      qcor,smec,smtc,rhoc,ncore,0d0,0d0)
          call poppr
C     ... here if v0 used to make core... then effectively frozen
C         call p1kden(nr,rwgt,s_rout(ib)%rhoc,p_v1,p_v0,sum)
C         write(stdo,996) smec,sum,smec+sum
C 996     format('smec,sum,smec+sum=',3f14.6)
C         smec = smec+sum

          sumtc = sumtc + smtc
          sumec = sumec + smec
          deallocate(rhoc)
        else
          stc0 = s_spec(is)%stc
          sumtc = sumtc + stc0
          sumt0 = sumt0 + stc0
        endif

C       Debugging
        if (ns4 == 4) then
        i = 1000*mod(s_ham%lrsa/10,10) + 444
        call rotspv(i,saxis,nr*nlml,nr*nlml,nr*nlml,v1,xx,v1,xx)
C       call prmx('Restore v1',v1,nr*nlml,nr,4)
        endif

        if (buffer_output) then
C         print *, 'buffer_output',ib,nproc,iproc,len_trim(buf(1))
          call awrit0('%-2&',buf,1,stdo)
          call poppr
        endif

      enddo  ! loop over sites

      if (nproc > 1) then
! it will be better to use a custom mpi type for reducing all these in a single call.
         call mpi_allreduce(mpi_in_place, s_atparms%vvesat, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%cpnves, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhovv1, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhovv2, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhvxt1, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhvxt2, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%qloc  , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%aloc  , 4, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%qcor1 , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%qlocc , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhoexc, 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rep1  , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rep2  , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rmu1  , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rmu2  , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%vales1, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%vales2, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%qv1   , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%qv2   , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhoex , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhoec , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhovxc, 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%focexc, 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%focex , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%focec , 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%focvxc, 2, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rvext , 4, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%atrue , 4, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%bfam  , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%sgpote, 1, mpi_real8, mpi_sum, comm, err)

         if (mod(job/1000,10) == 1) then
           call mpi_allreduce(mpi_in_place, s_atparms%rvepsv , 1, mpi_real8, mpi_sum, comm, err)
           call mpi_allreduce(mpi_in_place, s_atparms%rvexv  , 1, mpi_real8, mpi_sum, comm, err)
           call mpi_allreduce(mpi_in_place, s_atparms%rvecv  , 1, mpi_real8, mpi_sum, comm, err)
           call mpi_allreduce(mpi_in_place, s_atparms%rvvxcv , 1, mpi_real8, mpi_sum, comm, err)
           call mpi_allreduce(mpi_in_place, s_atparms%rveps  , 1, mpi_real8, mpi_sum, comm, err)
           call mpi_allreduce(mpi_in_place, s_atparms%rvvxc  , 1, mpi_real8, mpi_sum, comm, err)
         endif

! lfoc can differ between atoms
         call mpi_allreduce(mpi_in_place, s_atparms%rhcvef1, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, s_atparms%rhcvxt1, 1, mpi_real8, mpi_sum, comm, err)


         call mpi_allreduce(mpi_in_place, qval , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, qsc  , 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, sumtc, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, sumec, 1, mpi_real8, mpi_sum, comm, err)
         call mpi_allreduce(mpi_in_place, sumt0, 1, mpi_real8, mpi_sum, comm, err)

         allocate(smap(0:nproc-1))
         smap = pmap(1:nproc) - pmap(0:nproc-1)

         if (mod(job,10) == 1) then
C           Setup MPI types to simplify call to mpi_allgatherv
            call mpi_type_contiguous( nab*n0*nsp, mpi_real8, ab_mt, err)  ! ab_mt = MPI type
            call mpi_type_contiguous(nppn*n0*nsp, mpi_real8, pn_mt, err)  ! pn_mt = vector MPI type

            call mpi_type_commit(ab_mt, err)
            call mpi_type_commit(pn_mt, err)

            call mpi_allgatherv(mpi_in_place, smap(iproc), ab_mt, hab, smap, pmap, ab_mt, comm, err)
            call mpi_allgatherv(mpi_in_place, smap(iproc), ab_mt, vab, smap, pmap, ab_mt, comm, err)
            call mpi_allgatherv(mpi_in_place, smap(iproc), ab_mt, sab, smap, pmap, ab_mt, comm, err)
            call mpi_allgatherv(mpi_in_place, smap(iproc), pn_mt, ppnl, smap, pmap, pn_mt, comm, err)

            call mpi_type_free(pn_mt, err)
            call mpi_type_free(ab_mt, err)

            if (s_optic%loptic > 0) then
              call mpi_type_contiguous(size(s_optic%rgrad)/nbas, mpi_real8, rg_mt, err)
              call mpi_type_contiguous(size(s_optic%rgrade)/nbas, mpi_real8, re_mt, err)

              call mpi_type_commit(rg_mt, err)
              call mpi_type_commit(re_mt, err)

              call mpi_allgatherv(mpi_in_place, smap(iproc), rg_mt, s_optic%rgrad , smap, pmap, rg_mt, comm, err)
              call mpi_allgatherv(mpi_in_place, smap(iproc), re_mt, s_optic%rgrade, smap, pmap, re_mt, comm, err)

              call mpi_type_free(rg_mt, err)
              call mpi_type_free(re_mt, err)
            end if

            if (mod(job/100,10) == 1) then
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sighhx;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % tauhhx;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % pihhx ;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sigkkx;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % taukkx;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % pikkx ;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sighkx;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % tauhkx;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % pihkx ;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
            else
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sighh;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % tauhh;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % pihh ;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sigkk;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % taukk;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % pikk ;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sighk;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % tauhk;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % pihk ;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
            end if

            if (mod(lso,10) /= 0) then
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sohh;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sokk;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
               do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % sohk;  end do
               call adstrb(nbas, pl, nproc, pmap, iproc, comm)
            end if
         end if

         if (mod(job/100,10) == 1 .or. mod(job/10,10) == 0 .and. mod(job,10) == 1) then
            do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % v0;  end do
            call adstrb(nbas, pl, nproc, pmap, iproc, comm)
         end if

         do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % v1;  end do
         call adstrb(nbas, pl, nproc, pmap, iproc, comm)

         if (ns4 == 4) then
! probably faster with noncontiguous mpi type avoiding the buffering in tmp.
            allocate(tmp(6,nbas))
            do ib = pmap(iproc)+1, pmap(iproc+1)
               tmp(1:3,ib) = s_site(ib)%saxis
               tmp(4:6,ib) = s_ham%eula(ib,1:3)
            end do

            call mpi_type_contiguous(6, mpi_real8, tmp_mt, err)
            call mpi_type_commit(tmp_mt, err)

            call mpi_allgatherv(mpi_in_place, smap(iproc), tmp_mt, tmp, smap, pmap, tmp_mt, comm, err)

            call mpi_type_free(tmp_mt, err)

            do ib = 1, nbas
               s_site(ib)%saxis   = tmp(1:3,ib)
               s_ham%eula(ib,1:3) = tmp(4:6,ib)
            end do

            deallocate(tmp)
         end if

         deallocate(smap)
      end if

      deallocate(pmap)

      eterms(10) = s_atparms%rhcvef1 + s_atparms%rhcvxt1
      if (mod(job,10) == 1) then ! Not if mkpot called for rhout
        eterms(8) = sumec
        eterms(9) = sumtc
        eterms(12)= sumt0
      endif

C     s_atparms%rvext(2,1) = s_atparms%rvext(2,1) + s_atparms%sgpote
C      print *, 's_atparms%rhvxt1,s_atparms%rhvxt2',s_atparms%rhvxt1,s_atparms%rhvxt2,s_atparms%sgpote
      s_atparms%rhvxt1 = s_atparms%rvext(1,1)
      s_atparms%rhvxt2 = s_atparms%rvext(2,1)
C      print *, 's_atparms%rhvxt1,s_atparms%rhvxt2',s_atparms%rhvxt1,s_atparms%rhvxt2,s_atparms%sgpote

      if (iprint() >= 40) then
C        write(stdo,251)
C        write(stdo,250)
C     .    s_atparms%rep1(1)+s_atparms%rep1(2),s_atparms%rep2(1)+s_atparms%rep2(2),
C     .    s_atparms%rhoexc(1)+s_atparms%rhoexc(2),
C     .    s_atparms%rmu1(1),s_atparms%rmu2(1),s_atparms%rhovxc(1)
C        if (nsp == 2) write(stdo,253)
C     .    s_atparms%rmu1(2),s_atparms%rmu2(2),s_atparms%rhovxc(2),
C     .    s_atparms%rmu1(1)+s_atparms%rmu1(2),s_atparms%rmu2(1)+s_atparms%rmu2(2),
C     .    s_atparms%rhovxc(1)+s_atparms%rhovxc(2)
C         write(stdo,252)
C     .    s_atparms%vales1,s_atparms%vales2,s_atparms%vales1-s_atparms%vales2,
C     .    s_atparms%rhovv1,s_atparms%rhovv2,s_atparms%rhovv1-s_atparms%rhovv2,
C     .    s_atparms%qv1,s_atparms%qv2,s_atparms%qloc
C         call info5(40,0,0,' core chg:%;15,6D%;15,6D%;15,6D',
C     .     s_atparms%qcor1,s_atparms%qcor1-s_atparms%qlocc,s_atparms%qlocc,0,0)
C         if (jobl == 1000000) then
C           write(stdo,258) s_atparms%rvext(1,1),s_atparms%rvext(2,1),
C     .       s_atparms%rvext(1,1) - s_atparms%rvext(2,1)
C         endif
C        if (nsp == 2) write(stdo,254) s_atparms%atrue(0)
C  251   format(/' Local terms accumulated over augmentation sites:'/18x,
C     .    'true',11x,'sm,loc',6x,'difference')
C  250   format(' rhoeps:  ',3f15.6/' rhomu:   ',3f15.6)
C  253   format(' spin2:   ',3f15.6/' total:   ',3f15.6)
C  252   format(' val*ves  ',3f15.6/' val*vef  ',3f15.6/' val chg: ',3f15.6)
C  258   format(' val*vext ',3f15.6)
C  254   format(' val mom: ', f15.6)
CC  255   format(a,           3f15.6,a,f11.6)
CC  256   format(a,           3f15.6,a,f10.6)
C  257   format(' core chg:',3f15.6)

        call awrit0('%N Local terms accumulated over augmentation sites:%N%18ftrue%11fsm,loc%6fdifference',
     .    buf,160,stdo)
        call awrit6(' rhoeps:  %;15,6D%;15,6D%;15,6D%N rhomu:   %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%rep1(1)+s_atparms%rep1(2),s_atparms%rep2(1)+s_atparms%rep2(2),
     .    s_atparms%rhoexc(1)+s_atparms%rhoexc(2),
     .    s_atparms%rmu1(1),s_atparms%rmu2(1),s_atparms%rhovxc(1))
        if (nsp == 2)
     .    call awrit6(' spin2:   %;15,6D%;15,6D%;15,6D%N total:   %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%rmu1(2),s_atparms%rmu2(2),s_atparms%rhovxc(2),
     .    s_atparms%rmu1(1)+s_atparms%rmu1(2),s_atparms%rmu2(1)+s_atparms%rmu2(2),
     .    s_atparms%rhovxc(1)+s_atparms%rhovxc(2))
        call awrit6(' val*ves  %;15,6D%;15,6D%;15,6D%N val*vef  %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%vales1,s_atparms%vales2,s_atparms%vales1-s_atparms%vales2,
     .    s_atparms%rhovv1,s_atparms%rhovv2,s_atparms%rhovv1-s_atparms%rhovv2)
        call awrit3(' val chg: %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%qv1,s_atparms%qv2,s_atparms%qloc)
        call awrit3(' core chg:%;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .     s_atparms%qcor1,s_atparms%qcor1-s_atparms%qlocc,s_atparms%qlocc)
        if (nsp == 2) call awrit1(' val mom: %;15,6D',buf,160,stdo,s_atparms%atrue(0))

      endif

C ... Electric field gradient
      if (ipr > 50) then
        call elfigr(nbas,stdo,zz,efg)
      endif
      deallocate(efg,zz)
      deallocate(rhol1,rhol2,v1,v2,v1es,v2es)
      call tcx('locpot')

      end

      subroutine elfigr(nc,stdo,z,efg1)
C- Computation of electric field gradient
C ----------------------------------------------------------------------
Ci Inputs
Ci   nc    :number of classes or sites
Ci   stdo  :standard output
Ci   z     :nuclear charge
Cio Inputs/Outputs
Cio  efg1  :input:  l=2 part of electric field at nucleus
Cio        :output: Electric field gradient
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Sep 08 Adapted from old FP (A Svane)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nc,stdo
      double precision z(nc),efg1(5,*)
C ... Local parameters
      integer i,ic,ifesn,j,n,ier,imax,ixx,iyy
      character*(1024) buf
      double precision v(3,3),vi(3,3),d(3),e(3),e2(3),tau(2,3)
      double precision tr(3,3),ti(3,3),da(3)
      double precision conv1,emmsfe,emmssn,Qfe,Qsn,pi,f0,s3,conv2,conv3,
     .  ax,bx,cx,dx,ex,dmax,eta,split
      procedure(real(8)) dlength
      data conv1 /162.1d0/
      data emmsfe,emmssn /4.8d-8,7.97d-8/
      data Qfe,Qsn /0.21d-28,-0.08d-28/

      pi = 4d0*datan(1d0)
C     remember sign:
      f0 = -dsqrt(15d0/16d0/pi)
      s3 = dsqrt(3d0)
C     Axel's conv2 is 486.3e19 Ry/e/au^2 -> V/m^2.
C     A more precise number is 13.6057/.52917725e-10/.52917725e-10 = 4.859e21 V/m^2
C     To convert E field in Ry/e/au to V/m, the conversion is 13.6057/.52917725e-10 = 2.5711e11
      conv2 = conv1*3d0

C      if (iprint() < 30) return

C      write(stdo,'(/," Electric Field Gradient:")')
C      write(stdo,'(" Site ",5x,"Tensor axes",5x,"esu/cm^2",2x,
C     .          "  V/m^2  ",4x,"eta",4x,"line splitting")')
C      write(stdo,'("      ",5x,"           ",5x," x10^13 ",2x,
C     .          "  x10^19 ",4x,"   ",4x,"    (mm/s)")')
C
      call awrit0('%N Electric Field Gradient:',buf,160,stdo)
      call awrit0(' Site %5fTensor axes%5fesu/cm^2%2f  V/m^2  %4feta%4fline splitting',
     .  buf,160,stdo)
      call awrit0('      %5f           %5f x10^13 %2f  x10^19 %4f   %4f    (mm/s)',buf,160,stdo)

C --- For each class or site, make and print electric field gradient ---
      do  ic = 1, nc
      if (z(ic) < 0.01d0) cycle
      if (dlength(5,efg1(1,ic),1) == 0) cycle
      ax = efg1(1,ic)

      ifesn = 1
      if (dabs(z(ic)-26d0) < 0.01d0) then
        conv3 = conv2*1.0d19*0.5d0*Qfe/emmsfe
      else if (dabs(z(ic)-50d0) < 0.01d0) then
        conv3 = conv2*1.0d19*0.5d0*Qsn/emmssn
      else
        ifesn = 0
        conv3 = 0d0
      endif

C      ax=efg1(1,ic)+efg2(1,ic)
C      bx=efg1(2,ic)+efg2(2,ic)
C      cx=efg1(3,ic)+efg2(3,ic)
C      dx=efg1(4,ic)+efg2(4,ic)
C      ex=efg1(5,ic)+efg2(5,ic)
      ax = efg1(1,ic)
      bx = efg1(2,ic)
      cx = efg1(3,ic)
      dx = efg1(4,ic)
      ex = efg1(5,ic)
      v(1,1) = 2*ex-2/s3*cx
      v(2,2) = -2*ex-2/s3*cx
      v(3,3) = 4/s3*cx
      v(1,2) = 2*ax
      v(1,3) = 2*dx
      v(2,3) = 2*bx
      v(2,1) = v(1,2)
      v(3,1) = v(1,3)
      v(3,2) = v(2,3)
      do  i = 1, 3
      do  j = 1, 3
        v(i,j) = f0*v(i,j)
        vi(i,j) = 0d0
        tr(i,j) = 0d0
        ti(i,j) = 0d0
      enddo
        tr(i,i) = 1d0
      enddo

      n = 3
      call htridi(n,n,v,vi,d,e,e2,tau)
      call imtql2(n,n,d,e,tr,ier)
      if (ier > 0) call rx(' ELFIGR : IER ne 0')
      call htribk(n,n,v,vi,tau,n,tr,ti)
      do  i = 1, 3
        da(i) = dabs(d(i))
      enddo
      dmax = 0d0
      imax = 0
      do  i = 1, 3
         if (da(i) > dmax) then
             dmax = da(i)
             imax = i
         endif
      enddo
c  d(imax) is field gradient (Vzz)
      ixx = mod(imax,3)+1
      iyy = mod(imax+1,3)+1
      eta = 0d0
      if(dabs(d(imax)) > 1.d-2) eta = dabs((d(ixx)-d(iyy))/d(imax))
c     do  i = 1, 3
c        write(stdo,98) conv1*d(i),conv2*d(i),(tr(j,i),j=1,3)
c     enddo
c     if (ifesn == 0) write(stdo,97) eta
      split = conv3*da(imax)*dsqrt(1d0+eta**2/3d0)
c     if (ifesn == 1) write(stdo,96) eta,split
c 98  format(3X,F12.4,' 10**13 esu/cm**2',3X,F12.4,' 10**19 V/m*2',
c    +         /,12X,3F12.6)
c 97  format(/' eta = ',F12.4/)
c 96  format(/' eta = ',F12.4,' line splitting = ',F12.4,' mm/s'/)

      do  i = 1, 3
        if (i == 1) then
C          write(stdo,'(i4,3x,3f6.3,2x,f8.2,2x,f8.2,5x,f6.3,5x,f8.5)')
C     .      ic,(tr(j,i),j=1,3),conv1*d(i),conv2*d(i),eta,split
          call awrit6('%,4i%3f%3;6,3D%2f%;8,2D%2f%;8,2D%5f%;6,3D%5f%;8,5D',buf,160,stdo,
     .      ic,tr(1,i),conv1*d(i),conv2*d(i),eta,split)
        else
C          write(stdo,'(4x,3x,3f6.3,2x,f8.2,2x,f8.2,5x)')
C     .      (tr(j,i),j=1,3),conv1*d(i),conv2*d(i)
          call awrit4('%7f%3;6,3D%2f%;8,2D%2f%;8,2D%?#n==3#%N##',buf,160,stdo,
     .      tr(1,i),conv1*d(i),conv2*d(i),i)
        endif
      enddo
      enddo

      end

      subroutine magvec(nr,nlml,nsp,nspc,rwgt,rho,amom)
C- Returns l=0 magnetic moment or noncollinear magnetization vector
      implicit none
C ... Passed parameters
      integer nr,nlml,nsp,nspc
      double precision rwgt(nr),rho(nr,nlml,nsp*nspc), amom(0:3)
C ... Local parameters
      double precision y0,ddot,dlength,az


      y0 = 1/dsqrt(16*datan(1d0))
      az = ddot(nr,rwgt,1,rho,1)-ddot(nr,rwgt,1,rho(1,1,nsp),1)

      if (nspc == 1) then
        amom(0) = az/y0
      else
        amom(1) =  2*ddot(nr,rwgt,1,rho(1,1,3),1)/y0
        amom(2) = -2*ddot(nr,rwgt,1,rho(1,1,4),1)/y0
        amom(3) =  az/y0
        amom(0) = dlength(3,amom(1),1)
      endif
      end


