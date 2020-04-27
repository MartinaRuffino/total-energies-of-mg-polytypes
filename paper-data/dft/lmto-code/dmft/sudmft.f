C#define H5
      subroutine sudmft(s_bz,s_ctrl,s_ham,s_lat,s_pot,s_spec,s_site,s_optic,s_gw,s_dmft,
     .  dmxp,iter,s_strn)
C- Driver linking lmf to DMFT, or to analyze properties of DMFT-generated Green's function
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  w nkp ef egap star zinsul wtkp lmet nkabc ntet n
Ci                 nevmx efmax fsmom ndos dosw def numq lshft nef qp
Ci                 sopertp range dos idtet ipq pdos wtkb swtk
Co     Stored:     nef egap numq wtkp wtkb ndos dosw ef def sopertp w n
Co                 dos idtet ipq pdos qp star swtk
Co     Allocated:  qp wtkb swtk
Cio    Elts passed: nkabc lshft ipq zinsul wtkp wtkb swtk lio sopertp
Cio                qp nef numq idtet egap n w def
Cio    Passed to:  bndfp rdsigm subzi addrbl sosite optinq optin2
Cio                optint iorsf bcast_strx iinit mpibc1
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lfp elind nl quit lekkl ldos lmet zbak plbnd lwsig
Ci                 ips lham loptc lfrce cllst clp clssl group initc ics
Ci                 idcc ipc ipcp mxcst nrc ncomp pgfsl pgfvl pgord
Ci                 pgplp spid dclabl clabl pos rmax
Co     Stored:     lwsig plbnd lfp cllst clp clssl group initc ics idcc
Co                 ipc ipcp ips mxcst nrc ncomp pgfsl pgfvl pgord pgplp
Co                 spid clabl dclabl pos rmax
Co     Allocated:  *
Cio    Elts passed: lcd lfp lbas lrs nspec lncol lbxc maxmem lwsig ldos
Cio                spid ips
Cio    Passed to:  bndfp grp2av suham2 rdsigm siged fpopm optinq vcdmel
Cio                rxes dfrce iorsf bcast_strx iinit mpibc1
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham nqsig ndham oveps ovncut nlmto lsig pmin pmax
Ci                 pnudef pwmode pwemin pwemax nbf lncol elind evals
Ci                 eseavr sigp eterms eula lrsa rsrnge rsstol nprs qsig
Ci                 offH ndhrs
Co     Stored:     lsig eterms lncol elind sigp ehf ehk eula ndhrs
Co                 eseavr nqsig iaxs
Co     Allocated:  nprs iaxs hrs qsig
Cio    Elts passed: lncol evals nlmto offH iprmb elind lrsa magf qsig
Cio                eula nprs iaxs hrs rsrnge
Cio    Passed to:  bndfp mkpot locpot suham2 rdsigm hft2rs siged hambls
Cio                hambl augmbl sopert3 mkhso blsig makusq addrbl
Cio                sosite mkrout mkehkf mkekin
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat nabc gmax npgrp nsgrp vol afmt pos ng
Ci                 kv gv awald tol nkd nkq tolft kv2 igv igv2 napw ag
Ci                 bgv cg cy dlv gvq indxcg ips0 istab jcg qlv symgr s_
Ci                 sym
Co     Stored:     napw igv igv2 ng gmax nabc alat plat ag bgv cg cy
Co                 dlv gv gvq indxcg ips0 istab jcg kv kv2 pos qlv
Co                 symgr s_sym
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed: plat pos istab symgr ag igv2 nsgrp vol ng nabc kv
Cio                gv ips0 bgv cy cg indxcg jcg qlv dlv alat qlat igv
Cio                napw
Cio    Passed to:  bndfp rhopos mkpot ioden2 symsmr smves vesft vesgcm
Cio                mshvmt symvvl ugcomp ggugbl gfigbl fklbl gklbl
Cio                hhugbl hhigbl phhigb hklbl hsmbl hgugbl smvxc2
Cio                vxcnlm smvxcm smcorm locpot msh21c suham2 rdsigm
Cio                siged suqlst sugvec hambls hambl augmbl bstrux hxpbl
Cio                ghibl hxpgbl ghigbl hklgbl smhsbl hhibl phhibl hsibq
Cio                mkhso makusq pusq1 fpopm addrbl fsmbl rsibl rlocbl
Cio                sosite sugvec0 vcdmel rxes mkrout symrat symrho
Cio                dfrce pvdf4 pvdf2 pvdf1 mkekin totfrc mixrho rhgcmp
Cio                iorsf bcast_strx iinit mpibc1
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  ppn qval smrout nlml nlma hab sab lcplxp vesrmt
Ci                 shfac bxcscali socscl rnew rhat ves aamom bxc cp
Ci                 ddpf dddpf ddpfr dlmwt dmatk dpf dpfr gibbs gma gmar
Ci                 grrme mad mxy palp papg pf pfnc pfr pmpol pnu pp
Ci                 pprel pti qc qcorr qnu qnur qpp qt rhos rhrmx sop
Ci                 thetcl vdif vintr vrmax vshft smpot smvextc smrho GFr
Co     Stored:     lcplxp qval ves aamom bxc cp ddpf dddpf ddpfr dlmwt
Co                 dmatk dpf dpfr gibbs gma gmar grrme hab sab mad mxy
Co                 palp papg pf pfnc pfr pmpol pnu pp ppn pprel pti qc
Co                 qcorr qnu qnur qpp qt rhat rnew rhos rhrmx socscl
Co                 sop shfac thetcl vdif vintr vrmax vshft smpot
Co                 smvextc smrho smrout GFr
Co     Allocated:  smrout hab sab
Cio    Elts passed: ppn qval smpot smvcnst lcplxp rnew hab sab smrout
Cio                smrho vesrmt rhat smvextc bxcscali v0beta
Cio    Passed to:  bndfp mkpot locpot suham2 mixrho iorsf bcast_strx
Cio                iinit mpibc1
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa z a nr rmt pz lmxl grp2 lmxb kmxt qc rg lfoca
Ci                 rfoca ctail etail stc p name rs3 eh3 vmtz orbp rsma
Ci                 idu mxcst coreh coreq ngcut idxdn ncomp idmod nxi
Ci                 exi chfa rsmfa rsmv kmxv rhoc
Co     Stored:     orbp ngcut idxdn name lmxa a nr rmt z lmxl kmxt p pz
Co                 qc lfoca rfoca coreh pb1 pb2 ctail etail stc nxi exi
Co                 chfa rsmfa rhoc
Co     Allocated:  rhoc
Cio    Elts passed:z rhoc
Cio    Passed to:  wlattc bndfp dfratm grp2av pgrp2a prrhat rhopos
Cio                siterhopos dfaugm lkdim mkpot rhomom corprm smves
Cio                vesgcm mshvmt symvvl ugcomp smvxcm smcorm smvxc4
Cio                elocp uspecb locpot gtpcor msh21c suham2 sugcut
Cio                rdsigm siged makidx nscpa dfqkkl suclst surho sumlst
Cio                hambls hambl augmbl bstrux smhsbl hsibq tbhsi
Cio                hsubblock hsibq2 hsibq4 mkhso makusq pusq1 mkorbm
Cio                mullmf mkpdos fpopm mkdmtu addrbl fsmbl fsmbpw rsibl
Cio                rsibl1 rlocbl sosite psymrq1 iorbtm mshn3p mchan
Cio                vcdmel rxes mkrout symrat symrho pnunew dfrce pvdf4
Cio                pvdf2 pvdf1 mkekin mixrho rhgcmp rhogkl ftlxp pvmix5
Cio                pvmix3 pvmix7 iorsf pvsms2 bcast_strx
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos pnu pz class clabel rho1 rho2 rhoc rho1x
Ci                 rho2x rhocx v0 v1 sighh sighk sigkk pihh pihk pikk
Ci                 sighhx sighkx sigkkx pihhx pihkx pikkx sohh sohk
Ci                 sokk qkkl qhkl qhhl eqkkl eqhkl eqhhl tauhh tauhk
Ci                 taukk force vel bxc cpawt omg omgn domg dmat gc gcu
Ci                 gcorr gii sfvrtx j0 pdos tauhhx tauhkx taukkx thet
Co     Stored:     saxis pnu pz force pos pos0 vel spec clabel bxc
Co                 cpawt omg omgn domg dmat gc gcu gcorr gii sfvrtx j0
Co                 pdos rho1 rho2 rhoc rho1x rho2x rhocx qhhl qhkl qkkl
Co                 eqhhl eqhkl eqkkl sighh sighk sigkk tauhh tauhk
Co                 taukk pihh pihk pikk sohh sohk sokk sighhx sighkx
Co                 sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx thet
Co                 v0 v1
Co     Allocated:  rho1 rho2 rhoc rho1x rho2x rhocx sigkk taukk sigkkx
Co                 taukkx pikk sokk pikkx sighk tauhk sighkx tauhkx
Co                 pihk sohk pihkx sighh tauhh sighhx tauhhx pihh sohh
Co                 pihhx qkkl eqkkl qhkl eqhkl qhhl eqhhl v0 v1
Cio    Elts passed: rho1 rho2 rhoc rho1x rho2x rhocx v0 qkkl qhkl qhhl
Cio                eqkkl eqhkl eqhhl pihhx tauhh pikkx taukk pihkx tauhk
Cio    Passed to:  readindmfl wlattc bndfp dfratm grp2av pgrp2a prrhat
Cio                rhopos siterhopos dfaugm lkdim mkpot rhomom smves
Cio                vesgcm mshvmt symvvl ugcomp smvxcm smcorm smvxc4
Cio                elocp locpot msh21c suham2 dfqkkl suclst surho
Cio                sumlst hambls hambl augmbl bstrux smhsbl hsibq
Cio                hsubblock hsibq2 hsibq4 mkhso makusq pusq1 mkorbm
Cio                mullmf mkpdos mkdmtu addrbl fsmbl fsmbpw rsibl
Cio                rsibl1 rlocbl sosite psymrq1 mshn3p mchan vcdmel
Cio                rxes mkrout symrat symrho pnunew dfrce pvdf4 pvdf2
Cio                pvdf1 mkekin totfrc mixrho rhgcmp rhogkl ftlxp
Cio                pvmix5 pvmix3 pvmix7 iorsf pvsms2 bcast_strx
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  nkabc eseavr lshft
Co     Stored:     eseavr
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bndfp rdsigm siged
Cio  s_dmft
Ci     Elts read:  lsigim ndsigi knorm ldadc dcmode ndsig nicix ncix
Ci                 nomg beta omg sig nlohi icix l ncatom ndim sigind
Ci                 iasig ib nzsig nzsigi
Co     Stored:     dcmode omg iproj lsigim gammac ncatom ncix nicix cf
Co                 sigind nzsigi iasig ndsig nzsig
Co     Allocated:  *
Cio    Elts passed: omg sig nicix ndim ndsig lsigim nomg beta iproj
Cio                gammac ib icix l qsplit ncatom sigind cf nzsigi iasig
Cio    Passed to:  readindmfl makeproj iodmftu embed_sigma compressproj
Cio                makegloc2 makeeimp2 makedelta3 sigp2sigij
Ci Inputs
Ci   dmxp  :vector of mixing parameters; see mixrho.f for dmxp(1..25)
Ci         :Additionally:
Ci         :dmxp(26)  is the iteration where the minimum RMS error was found
Ci         :dmxp(27)  is the minimum RMS error.
Ci         :dmxp(33)  is the Lindhard parameter
Ci   iter  :current iteration number
Co Outputs
Cl Local variables
Cl  iq,iqs :index to irreducible qp, and corresponding point in full BZ
Cl  lsig0  :if T, siginp is zero => no double counting
Cl  lmxax  :largest augmentation l cutoff, for dimensioning
Cl  ldcix  :dimension of largest cix block, for dimensioning
Cl  ndham  :dimensioning parameter, at least as large as largest
Cl         :hamiltonian dimension
Cl  nspx   :number of independent spin channels containing eigenvalues
Cl         :and total DOS; nspx=nsp unless nspc=2, in which case nspx=1
Cl  nspc   :1 in collinear case, 2 in noncollinear
Cl  ndimh  :Hamiltonian dimension for one spin and current k.
Cl         :In the APW case with pwmode>10 it may be less than ndham
Cl  ndimhx :ndimh*nspc
Cl  gpmode :grpt mode; see command line switches
Ci  emu    :energy - chemical potential : used to interpolate sigma
Cl  LW5    :set s_ctrl%lwsig=LW5 to cause bndfp to write evec file
Cl  LW8    :if given s_ctrl%lwsig==LW8, lmfdmft is to read use existing evecs file
Cs Command-line switches
Cs   --job=#     :Set DMFT job to #
Cs               :job 0 creates files lattice, class, dmftin
Cs               :job 1 execute mode
Cs   --ldadc=#   :Set double counting to #
Cs   --dos       :Makes DOS on the real axis (formula is not correct if Sigma depends on omega)
Cs   --makesigqp :--makesigqp[:mode=#] makes quasiparticle DMFT sigma,
Cs               :in format of QSGW sigma file
Cs   --novxc     :not used now
Cs   --udrs      :Generate and save output density
Cs   --gprt      :--gprt[~rdsigr=fn][~rdsigr=fn][~broad=#][~mode=#][~nom=#]
Cs               :Generate and write gloc to disk acccording to mode
Cs               : mode  write
Cs               :  0    (default) make gloc, delta, eimp
Cs               :  1    k-integrated gloc (this branch doesn't work; use mode=3)
Cs               :  2    k-resolved gkloc
Cs               :  3    both 1 and 2
Cs               :  4    k-resolved hkloc (has not been maintained)
Cs               :  5    both 1 and 2 in indmfl compressed format
Cs               :  8    k-resolved g_ij(k,omega), in one-particle basis
Cs               :  9    k-resolved sigma(k,omega), in one-particle basis
Cs               : 18    k-resolved g_ij(k,omega), diagonal part only, in one-particle basis
Cs               : 19    k-resolved sigma(k,omega), diagonal part only, in one-particle basis
Cs               : 20    Like 19, but effective diagonal sigma defined from diagonal of gij^-1
Cs               : rdsigr=fn reads a new DMFT (frequency, sigma) set from file fn
Cs               : Frequencies are assumed to be on the real axis.
Cs               : nom restricts the calculation to # frequency points
Cs               : rdsig=fn reads a new DMFT (frequency, sigma) set from file fn
Cs               : Frequencies are assumed to Matsubara frequencies.
Cs               : nom=# restricts the calculation to # frequency points
Cs               : wse writes the diagonal part of the self-energy to file se.
Cs               : Only valid with mode 9
Cs   --pade      : Analytically continue DMFT sigma to real axis by Pade approximant
Cr               : Example: --pade~nw=121~window=-4/13.6,4/13.6~icut=30,75
Cr Remarks
Cr   This is an alpha-test version of an interface to K. Haule's DMFT code.
Cr   Additional files required are:
Cr     indmfl  : information about correlated blocks
Cr     sig.inp : File with self-energies for correlated subblocks
Cr
Cr   One or more subsystems are projected and/or embedded, using DMFT to generate
Cr   the self-energy for one or more subsystems.  In this version, it is assumed that:
Cr   a. A subsystem is defined by the augmentation partial waves of a particular l, following Haule.
Cr   b. The reference (host) Green's function is noninteracting.
Cr
Cr   sudmft provides informtion needed to the DMFT self-energy solver. It requires:
Cr   1. a. Which frequencies the DMFT self-energy is defined on.  At present sudmft uses Matsubara frequencies.
Cr      The frequency mesh and DMFT self-energies are read/created by routine iodsig.
Cr      b. A choice of double-counting.  This information is read by routine iodmftdc.
Cr      c. A definition of the effective local interaction.
Cr   2. partitioning into one or more subsystems, which requires
Cr      a: information about the site and l- quantum number of the subsystem
Cr      The subsystem is defined, and which particular matrix elements are computed,
Cr      the number of subsystems (cix blocks) (both the full number and the number of inequivalent ones)
Cr      are set up by routine readindmfl, which requires indmfl.ext.  In future, this may change.
Cr      b: a definition of the projector
Cr      c: states to be included in the projector
Cr      d: The one-body hamiltonian from which the local G is projected, for each k
Cr      The projector is constructed in routine lmlproj.
Cr   3. The self-energy linking the bath to the local system, called the hybridization function delta.
Cr      The one-body hamiltonian of the local system is given by the local part of the
Cr      noninteracting reference hamiltonian and delta.
Cr      In the usual setup execute (job=1, no special branches) sudmft calls makedelta3 to make delta.
Cr      The interface to Haule's code also requires Eimp, which sudmft makes as well.
Cr
Cr   Once a self-energy has been generated, each block must be embedded into the reference hamiltonian
Cr   to obtain a crystalline Green's function G(omega,k).
Cr   Embedding is required to obtain the Fermi level and other quantities such as the charge density
Cr   (specified through --udrs) or other special-purpose objects (set through --gprt).
Cr   Routine embed_sigma embeds the local self-energy to make sigma in the basis of one-body eigenstates.
Cr
Cb Bugs
Cb   *Needs rewriting with multiple cix
Cu Updates
Cu   26 Jul 18 (MvS) revised management of LDA double counting
Cu   09 Feb 18 (MvS) New --gprt mode 20
Cu   17 Jan 18 (MvS) First cut at SO coupled case
Cu   28 Nov 17 (MvS) --gprt modes 2-4 work with MPI
Cu   30 Apr 17 (MvS) --gprt modes 18 and 19 can record for frequencies on Matsubara axis
Cu   15 Apr 17 (MvS) new --pade
Cu   31 Mar 17 (MvS) new --gprt:rdsig and --gprt:rdsigr
Cu   05 Oct 16 (Sponza) sudmft can print out gloc
Cu   09 Aug 16 (Sponza) sudmft can generate output charge density
Cu   27 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu   14 Feb 16 (MvS) New MPI version, parallelizing over omega
Cu   13 Feb 16 (MvS) New Pade interpolation of sigma for DOS
Cu   13 Jan 16 (MvS) new option --makesigqp to write embedded static sigma
Cu   07 Dec 15 (MvS) Improved printout in spin polarized case
Cu   06 Dec 15 (Sponza)  DOS maker generates spin polarized DOS
Cu   01 Dec 15 (Sponza)  Revised eimp maker
Cu   30 Oct 15 (Sponza)  Revised Fermi level finder
Cu   29 Oct 15 (Sponza)  Revisions for spin polarized case
Cu   29 Oct 15 (MvS) Added --dos switch, fermi level finder for insulators
Cu   23 Oct 15 (MvS) sudmft makes use of crystal symmetry
Cu   22 Aug 15 (Sponza)  Find Charge neutrality point
Cu   22 Aug 15 (Pisanti) spin polarized
Cu    1 May 15 (MvS) routines makes QSGW sigma for double counting
Cu   21 Apr 15 (Pisanti) intermediate stage, adding DMFT sigma and double counting
Cu   25 Nov 14 First created by P Pisanti.
C ----------------------------------------------------------------------
      use mpi
      use structures
C#ifdef H5
      use h5
C#endif
      use meshes, only : khmesh
      implicit none
C ... Passed parameters:
      integer iter
      double precision dmxp(33)
C ... For structures
!      include 'structures.h'
      type(str_bz)   :: s_bz
      type(str_ctrl) :: s_ctrl
      type(str_ham)  :: s_ham
      type(str_lat)  :: s_lat
      type(str_pot)  :: s_pot
      type(str_spec) :: s_spec(*)
      type(str_site) :: s_site(*)
      type(str_optic):: s_optic
      type(str_gw)   :: s_gw
      type(str_strn) :: s_strn(*)
      type(str_dmft) :: s_dmft
C     Local structures
      type (str_mpade):: s_pade
      type(str_gwse)::    s_gwse

C ... Dynamically allocated local arrays
      integer, allocatable :: ips(:),ipc(:),ipcx(:),lmxa(:),ipb(:),konft(:,:,:),iblst(:),nicixi(:)
      integer, allocatable :: ifblst(:) ! List of orbitals entering into color weights
      real(8), allocatable :: bas(:,:)  ! Basis vectors
      real(8), allocatable :: evlk(:,:) ! Eigenvalues for 1 k-point
      real(8), allocatable :: eigiq(:,:,:) ! Eigenvalues for all k-points, window of bands
      real(8), pointer     :: siginp(:,:) ! Matrix elements of sigma
      real(8), pointer     :: ppnl(:)   ! NMTO-like potential parameters
      real(8), allocatable :: qbyl(:)   ! l-decomposed site charge
      real(8), allocatable :: hbyl(:)   ! l-decomposed site eigenvalue sum
      real(8), pointer     :: omega(:)  ! Frequencies for DMFT.  Whether real or imaginary depends on s_dmft%lsigim
      real(8), allocatable :: omegap(:) ! Frequencies for interpolation to the real axis
      real(8), allocatable :: N(:,:)    ! total and k-contrib. to N. electrons band-by-band
      real(8), allocatable :: dNan(:,:) ! dN/dV where V is potential shift
      real(8), allocatable :: dos(:,:)  ! total DOS on the real axis
      real(8), allocatable :: sigpade(:,:,:) ! Pade extrapolation of Self-energy
      real(8), allocatable :: sig_ldadc(:) ! LDA double counting --- one number for each cix block
C     real(8), allocatable :: sigdc(:)  ! double counting --- one-to-one correspondence with delta
      real(8), allocatable :: wk(:)
      complex(8), allocatable :: zsp(:,:,:) ! Eigenvectors with both spins
      complex(8), allocatable :: ausp(:) ! val,slo of w.f. at MT sphere surface
      complex(8), allocatable :: dmftu(:,:,:,:) ! Factored form of projector
      complex(8), allocatable :: sigdc(:)
      complex(8), allocatable :: gloc(:,:,:,:,:)
      complex(8), allocatable :: gkloc(:,:,:,:,:,:)
      complex(8), allocatable :: sigbark(:,:,:) ! DMFT Sigma in basis of 1-particle eigenfunctions, one frequency
      complex(8), allocatable :: sigwq(:,:,:,:) ! diagonal DMFT Sigma in basis of 1-particle eigenfunctions
      complex(8), allocatable :: udhamk(:,:,:)
      complex(8), allocatable :: udeval(:,:,:,:), udevalko(:) ! updated eigenvalues
      complex(8), allocatable :: eimpm(:,:,:,:)
      complex(8), allocatable :: dmftdelta(:,:),dmftdeltah5(:,:),g0h5(:,:)
      complex(8), allocatable :: sigqp(:,:,:) ! QSGW selfen for DC correction to impurity selfen

C ... Local variables
      integer, parameter :: NULLI=-99999, PRTU=50, PRTG=70
      integer, parameter :: lSOf=4,lSzLz=32,lSzLzp=64,lSzLzp0=128
      integer, parameter :: n0=10, nppn=12, nkap0=4, LW5=5, LW8=8
      real(8), parameter :: eVbyRy = 13.60569193d0 !? 1Ry = 13.605 693 009(84) eV according to https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev|search_for=rydberd
      real(8), parameter :: pi = acos(-1d0)
      real(8),parameter  :: NULLR = -99999
      complex(8), parameter :: srm1=(0d0,1d0),NULLZ=(-99999D0,0d0)
      integer :: rankprocid                ! procid that captures sigbark(nomg), used to make eimp1
      integer, parameter :: rfalsi_limit=100
      integer, parameter :: lchk_sigbar=3    ! print check files on  P+E[sigbar(infty)] HARD CODED
      integer, parameter :: lchk_prjnrm=7    ! print check files on Normalized Projectors HARD CODED
      logical :: minimization_required       ! if .TRUE. charge neutrality not achieved
      integer :: rfalsi_niter                ! number of iteration of minimization procedure
      real(8) :: rfwk(12)                    ! working vector for rfalsi
      real(8) :: efshft                  ! Fermi level of interacting system = ef0-efshft
      real(8) :: efshfti(2)                  ! Instances of efshft in insulator mode

      logical lsig0,ldiaghbar
      integer n1,n2,n3,k1,k2,k3,nkabc(3),ngabc(3) ! For BZ
      equivalence (nk1,nkabc(1)), (nk2,nkabc(2)), (nk3,nkabc(3))
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer cix,i,iseh,ifi,ifis,ifiz,ig,ipass,ipr,iq,iqfbz,iqs,is,isp,j,j1,
     .  jfi,job,k,knorm,l,lchk,ldim,lfrce,lmxax,lrout,lrsig,lsig,lso,lpnu,dcmode,
     .  lwsig,maxcor(2),mnevl,mxint,napw,nat,nbas,ndham,ndhamx,ndimh,ndimhx,
     .  nevn,ncix,nicix,nl,nlmax,nlmto,nmcore,npgrp,nphimx,nqbnd,nfbn(4),
     .  nqsig,nsgrp,nsp,nspc,nspx,nspc3,ldcix,onesp,ovncut,stdo
      double precision alat,ebot,ef0,epsovl,esmear,etop,gcutb,gmax,qval,tim,ef,
     .  plat(3,3),pnu(n0,2),pnz(n0,2),qlat(3,3),qp(3),xv(20),ninsul,chempot,emu
      integer ifac(3)            ! used for iqstar
      double precision qb(9),qpr(3) ! used for iqstar
      complex(8) :: omfac        ! 1 or sqrt(-1) depending on s_dmft%lsigim
      integer lshft(3),nk1,nk2,nk3,nkp,nqi,nkfbz,nlohi(2),iomg,nkpan(100)
      integer nomg               ! Number of frequencies for DMFT mesh
      integer nomgx              ! Usually same as nomg, but if gpmode>0 nomgx = nomax_gpr
      double precision eseavr(2) ! Constant term for high-lying QSGW sigma
      character strn*160,fn*80
      character*(1024) fmt

C     For dos
      logical doslsw(9),lidos,ldig,lrange,ldpade
      integer dosisw(5),idfmt,npts
C     dosw(1:2) = energy window for DOS. dosw(3) preserves c.n. point.  dosw(4) = constant shift for writing DOS
      double precision dosw(4)
      real(8) :: Nnew(2), dNan_new(2) ! total density and derivative,resolved by spin
      equivalence (lidos,doslsw(1)),(ldig,doslsw(4)),(lrange,doslsw(3))
      equivalence (npts,dosisw(1)),(idfmt,dosisw(2))

C     For writing gkloc
      integer :: nomax_gpr      ! max. number of frequencies written to disk
      integer :: gpmode         ! Special modes for writing G or related properties; see description of --gprt above
      integer :: wgpmode        ! subset of gpmode where file I/O is done on the fly.  MPI is handled specially
      integer :: ifgk,jfgk,lascii ! File logical units for writing gpmode-related output, lascii => whether ascii or binary
      character*72 :: dc*1
      integer :: ii,jj,iv(12)   ! to parse arguments of gprt
      integer :: klist_u, qlist_u
      logical :: fullg
      logical :: quick
C     For density update
      logical lbin,lstar
      integer :: ib1,ib2,lekkl,numq,irs(5),lfrzw
      character fileid*68
      complex(8), allocatable :: evecl(:,:,:,:), evecr(:,:,:,:) ! left and right eigenvectors of Hbar
!     complex(8), allocatable :: eveclko(:,:), evecrko(:,:) ! left and right eigenvectors of Hbar at given energy, spin, k-point
      complex(8), allocatable :: denmat(:,:,:)    ! density matrix
      complex(8), allocatable :: dmftvec(:,:,:)   ! rotated QSGW wavefunctions
      real(8), allocatable    :: dmftwktb(:,:,:)  ! DMFT-weights of k-points
      complex(8), allocatable :: dmevec(:,:)      ! eigenvecs of density matrix
      complex(8), allocatable :: dmevecinv(:,:)   ! formerly eigenvecs of density matrix.  Now work array
      real(8), allocatable    :: dmeval(:,:,:)    ! eigenvals of density matrix
!     complex(8), allocatable :: dmeval(:,:,:)    ! eigenvals of density matrix
!     complex(8), allocatable :: dmevecinv(:,:,:,:) ! eigenvecs of density matrix
!     complex(8), allocatable :: dmevec(:,:,:,:)  ! eigenvecs of density matrix

      complex(8), allocatable :: evecq(:,:,:,:)   ! LDA/QSGW evecs, loaded for all q
      real(8), allocatable    :: evalq(:,:,:)     ! LDA/QSGW evals, read for all q

      complex(8), pointer     :: srout(:,:)
      real(8)                 :: sumqv(3,2),sumev(2,3)
      real(8), allocatable    :: dmatk(:),f(:),tso(:)


C     For embedding static sigma (obtained with mk_siginp-freq0.py)
      integer :: is0,ifs0 ! file sig.inp.f0 (input and output)
      real(8), allocatable :: sigf0(:,:,:)  ! dimensioned as sigpade

      integer :: khm_n, khm_l, khm_u
      real(8) :: khm_d, khm_m
      character(len=:), allocatable :: khm_usage
C     For chi0
      integer,allocatable ::i123(:,:)
      integer :: nq, nomgv, nomgw
      complex(8), allocatable :: chi0(:,:,:,:,:,:,:), chi0loc(:,:,:,:,:,:)


C     For MPI
      integer, allocatable :: kpproc(:) ! q loop blocking indices for multithreaded case
      integer, allocatable :: displs(:), recvcounts(:)
      procedure(integer) :: mpipid
      integer procid,master,nproc,err

C     Procedures
      procedure(logical) :: a2bin,cmdopt,isopen
      procedure(real(8)) :: ddot,dlength,cpusec,mpiquery
      procedure(integer) :: iorsf,iodsig,bitand,fopna,fopng,fopnx,fopnn,fopnig,fxst,iprint,
     .                      str_pack,a2vec,cmdoptswx,nglob,parg,wordsw,rdm,isw

      integer :: u ! unit to use in open(newunit=u, ...)

C --- Setup ---
      call tcn('sudmft')
      stdo = nglob('stdo')
      call getpr(ipr)
      nproc = mpipid(0)
      procid = mpipid(1)
      master = 0
      gcutb = -1                  ! Not needed for DMFT
      lchk  = 1                   ! Sanity checks
      ldim  = s_ham%ldham(1)
      knorm = s_dmft%knorm
      nullify(s_dmft%omg,s_dmft%sig)
      nullify(s_gwse%qp,s_gwse%sgq,s_gwse%eig,s_gwse%sigwq)
      lstar = .true.              ! Loop over star of k when looping 1 ... nkfbz

      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      ngabc = s_lat%nabc
      gmax = s_lat%gmax
      etop = NULLI
      ebot = -etop
      nmcore = 1000*(bitand(s_ctrl%lcd,16)/16)
      npgrp = s_lat%npgrp
      nsgrp = s_lat%nsgrp
      call fftz30(n1,n2,n3,k1,k2,k3)
      nbas = nglob('nbas')
      nsp  = nglob('nsp')
      nspc = nglob('nspc')    ! 2 for noncollinear case
      nspc3 = nglob('nspc3')  ! 2 if spins coupled perturbatively (lso=3)
      nspx = nsp / nspc         ! 1 if nsp=1 or if noncollinear
      nlmax = nglob('nlmax')
      lso =   isw(IAND(s_ham%lncol,4) /= 0)
     .    + 2*isw(IAND(s_ham%lncol,lSzLz) /= 0)
     .    + 3*isw(IAND(s_ham%lncol,lSzLzp) /= 0)
     .    + 4*isw(IAND(s_ham%lncol,lSzLzp0) /= 0)
!     lwvxc = .not. cmdopt('--novxc',7,0,strn)
      lsig = 1
C     lsig = 0
      nqsig = s_ham%nqsig
      allocate(ipb(nbas))
      allocate(ifblst(1))
      ndham = s_ham%ndham
      ndhamx = s_ham%ndham * nspc
      napw = 0  ! No PMT method yet
      esmear = s_bz%w
      nphimx = 3                ! Max number of partial waves/l = 3 : (phi, phidot, LO)
      ninsul = 0                ! A metal by default
      nnew = 0; dNan_new= 0     ! Initialize total charge and DOS for both spins
      lpnu = 0                  ! 1 updates linearization energies when making output rho

C ... No Pade initially
      emu = 0                   ! energy - mu : used to interpolate sigma
      ldpade = .false.; s_pade%lpade = .false.; nullify(s_pade%zomega,s_pade%padcof)

C ... For reduced hilbert space
      epsovl = s_ham%oveps
      ovncut = s_ham%ovncut
      mnevl = 9999 !for printout, smallest ham dimension
      nlmto = s_ham%nlmto

C ... Read job from command line
      job=-999; i = 0
      if (cmdopt('--job=',6,0,strn)) then
        i = 6
      elseif (cmdopt('-job=',5,0,strn)) then; i = 5
      endif
      if (i /= 0) then
        if (.not. a2bin(strn,job,2,0,' ',i,-1)) call
     .    rxs2('SUDMFT: failed to parse "',strn(1:30),'%a"')
      endif
      if (job == -999) then
        write(stdo,*) ' lmfdmft: input one of the following jobs:'
        write(stdo,*)
     .     '   0 : init mode; creates files lattice, class, sym.out'
        write(stdo,*) '   1 : DMFT execute mode'
        write(stdo,'(a)',advance='no') ' job? '
        read (5,*) job
      endif

      select case (job)
        case (0)
          call info0(20,1,0,
     .      ' ... sudmft job=0: create files lattice, class, dmftin')
        case (1)
          call info0(20,1,0,' ... sudmft job=1: make projectors')
        case DEFAULT; call rxi('SUDMFT: illegal job',job)
      end select

      fullg = cmdopt('--fullg',7,0,strn)
      quick= cmdopt('--quick',7,0,strn)
C ... Count number of atoms excluding floating orbital sites
      nat = 0
      do  i = 1, nbas
        is = s_site(i)%spec
        if (s_spec(is)%lmxa > -1) nat = nat + 1
        ipb(i) = nat
      enddo

C ... Assemble s_dmft <= information about DMFT-correlated orbitals
      if (procid == master) i = fxst('indmfl'); call mpibc1(i,1,2,.false.,'','')
      if (s_dmft%ncix == 0) then
         call rxx(i/=1,'sudmft: No indmfl file ... aborting')
      endif
      call readindmfl(s_site,s_dmft,s_lat,nsp)

      call info5(10,0,0,' indmfl file expects %i DMFT block(s) '//
     .  '(max dimension %i),  %i matrix elements',
     .  s_dmft%nicix,maxval(s_dmft%ndim),s_dmft%ndsig,0,0)
      call rxx(s_dmft%ndsig==0,'no matrix elements specified in indmfl file')
      nicix = s_dmft%nicix; ncix = s_dmft%ncix; ldcix = maxval(s_dmft%ndim)
      omfac = 1 ; if (s_dmft%lsigim) omfac = srm1

C ... Init mode
      if (job == 0) then

C   ... Make lmxa,bas,lmxax
        allocate(ips(nbas),lmxa(nbas),bas(3,nbas))
        call sitepack(s_site,1,nbas,'spec',1,ips,xv)
        do  i = 1, nbas
          lmxa(i) = s_spec(ips(i))%lmxa
          bas(1:3,i) = s_site(i)%pos
C          if (lmxa(i) > -1) then
C          call orbl(i,0,ldim,s_ham%iprmb,norb,ltab,ktab,xv,offl,i1)
C          if (gwcode == 0 .or. gwcode == 2) then
C            write(ifi,"(3i10)") i1
C          endif
C          endif
        enddo
        deallocate(ips)
        lmxax = mxint(nbas,lmxa)

C   ... Parameters needed for for class file
        allocate(ipc(nbas),ipcx(nbas))
        call sitepack(s_site,1,nbas,'class',1,ipc,xv)
        call pvsug1(nbas,lmxa,ipc,ipcx)
C       maxcor(1) = largest number of core radial functions in system
C       maxcor(2) = largest l for which a core state exists
        maxcor(1) = 0; maxcor(2) = -1
        allocate(konft(0:lmxax,nbas,nsp))
        do  i = 1, nbas
          is = s_site(i)%spec
          pnu = s_site(i)%pnu
          pnz = s_site(i)%pz
          do  isp = 1, nsp
            do  l  = 0, lmxa(i)
              konft(l,i,isp) = pnu(l+1,isp)
              if (mod(pnz(l+1,isp),10d0) < pnu(l+1,isp) .and.
     .            pnz(l+1,isp) > 0)
     .          konft(l,i,isp) = mod(pnz(l+1,isp),10d0)
              maxcor(1) = max(maxcor(1),konft(l,i,isp)-l-1)
              if (konft(l,i,isp)-l-1 > 0) then
                maxcor(2) = max(maxcor(2),l)
              endif
            enddo
          enddo
        enddo

        ifi = fopna('lattice',-1,2); rewind ifi
        jfi = fopna('class',-1,2); rewind jfi
        call wlattc(3,ifi,jfi,s_site,s_spec,alat,plat,nbas,
     .    nat,maxcor,gcutb,ipb,ipcx,lmxax,bas,lmxa,nsp,konft)
        call fclose(ifi); call fclose(jfi)

C  ... Setup structure to write dmftin file
C        ifi = fopna('dmftin',-1,2)
C        call iodmftp(2,s_dmft,-ifi)

C       call info0(10,1,0,' Files lattice, class, dmftin, written to disk')
        call rx0('finished job 0')

      endif

      if (nicix == 0) call rx0('no correlated orbitals')

C --- Create and/or read DMFT sigma and omega mesh ---
C     This could follow the band pass but doing it here
C     enables code to find inconsistencies early on.
      nomg = s_dmft%nomg
      if (associated(s_dmft%omg)) deallocate(s_dmft%omg)
      if (associated(s_dmft%sig)) deallocate(s_dmft%sig)
      allocate(s_dmft%omg(nomg),s_dmft%sig(nomg,2*s_dmft%ndsig))
      if (procid == master) then
        i   = fopnx('sig.inp',72,-1,-1)
        ifi = fopnx('sig.inp',2,0,-1); rewind ifi
        is0 = fopnx('sig.inp.f0',72,-1,-1)
      else
        ifi = 1
      endif
      call mpibc1(i,1,2,.false.,'','')
      call mpibc1(is0,1,2,.false.,'','')
C     If sig.inp and sig.inp.f0 do not exist, create sig.inp it and exit
      if (i /= 1 .and. is0 /= 1 ) then
        call info0(10,0,0,' Missing sigma file : create it and exit ...')
        forall (i = 1:nomg) s_dmft%omg(i) = pi*(2*i-1)/s_dmft%beta
        call dpzero(s_dmft%sig,size(s_dmft%sig))
        if (procid == master) i =
     .    iodsig(001,s_dmft%nomg,s_dmft%beta,s_dmft%ndsig,s_dmft%omg,s_dmft%sig,-ifi)
        call rx0('done writing template sigma')
C     To read sig.inp.f0 needs also sig.inp
C     That's stupid but I can't make it work otherwise... ???
      elseif (i /= 1 .and. is0 == 1) then
        call rx0('At the moment, embedding sig.inp.f0 requires also sig.inp.')
      elseif (i == 1 .and. is0 == 1) then
        call info0(20,1,0,' Found sig.inp.f0 to be embedded...')
      endif

      k = 131; if (nsp == 2 .and. cmdopt('--nmsig',7,0,strn)) k = k+400
      i = iodsig(k,s_dmft%nomg,s_dmft%beta,s_dmft%ndsig,s_dmft%omg,s_dmft%sig,ifi)
      omega => s_dmft%omg; siginp => s_dmft%sig; nomg = s_dmft%nomg
      lsig0 = dlength(size(siginp),siginp,1) < 1d-7
      if (lsig0) then
        call info0(30,0,0,' DMFT sigma is zero ... no double counting')
      elseif (cmdopt('--symsig',8,0,strn)) then
        call info2(30,0,0,' Symmetrizing sigma with %i group ops ...',s_lat%nsgrp,2)
        call symdmftsig(s_lat,s_dmft,nomg,s_dmft%sig)
      endif

      nlohi = s_dmft%nlohi
      if (.not. (nlohi(1) > 0 .and. nlohi(2) > nlohi(1)))
     .  call rx('sudmft : you must specify range of states to include in projector')

C --- Printout setup info ---
      call info8(20,1,0,
     .  ' SUDMFT:%?#(n>1)#%-1j %i blocks,## projector #%i with'//
     .  ' k-%?#n#resolved#integrated# norm.  %i bands in window (%i,%i)',
     .  nicix,s_dmft%iproj,knorm,nlohi(2)-nlohi(1)+1,nlohi(1),nlohi(2),0,0)
      if (s_dmft%lsigim) then
        call info5(20,0,0,'%9f%i Matsubara frequencies, interval (%g,%g) eV',
     .    nomg,omega(1)*eVbyRy,omega(nomg)*eVbyRy,0,0)
      else
        call info5(20,0,0,'%9f%i frequencies on real axis, interval (%d,%d) eV',
     .    nomg,omega(1)*eVbyRy,omega(nomg)*eVbyRy,0,0)
      endif

C ... Read --dc from command line
C      Lorenzo proposes dc term to be:
C      positive = constant double counting; negative = dynamical DC
C      The information is stored for now in sig_ldadc itself,
C      The variable sig_ldadc is assed assigned a value;
C      In future should be merged with sigdc (see makedelta3) and s_dmft%ldadc will be redundant.
       allocate(sig_ldadc(s_dmft%ndsig))   ! Should be replaced by sigdc
       call iodmftdc(s_dmft,sig_ldadc,nsp,ifi)

C  ... Set up double counting.  For now, just copy sig_ldadc to sigdc
       allocate(sigdc(s_dmft%ndsig)); call dpzero(sigdc,2*size(sigdc))
       if (lsig0 .and. (sig_ldadc(1) < 0d0)) call rx('sudmft: no sigdc') ! ?? SigDC computed otherwise, depending on the negative value of sig_ldadc
       if ((.not. lsig0) .and. (sig_ldadc(1) > 0d0)) then
         call dcopy(s_dmft%ndsig,sig_ldadc,1,sigdc,2)
         sig_ldadc = 0d0
      endif

C --- Analytic continuation of sigma to real axis by Pade ---
C     Abuse nomax_gpr = number of frequencies, xv = window, nlohi = icut
      if (cmdopt('--pade',6,0,strn)) then
        fmt = 'sig2'
        dc = strn(7:7)

        if (strn(8:10)=='fn=') then
        ii = wordsw(strn,dc,'fn=','',j)
          call nwordg(strn,0,dc,1,j,j1)
          fn = strn(j:j1)
          ifi = fopnx(fn,70,-1,-1)  ! read fn, name corresponds to true file fn.ext
          call info0(20,1,0,' pade interpolation: reading real-axis frequencies from file '//trim(fn))
          call rxx(ifi /= 1," missing file '"//trim(fn)//"'")
          ifi = fopna(trim(fn),-1,0)
          If (procid == master) then
            rewind ifi
            nomax_gpr = 0; j = 0
            i = rdm(ifi,0,0,' ',xv,nomax_gpr,j)
            if (i/=1 .or. nomax_gpr<=0 .or. j<0) call rx('no frequencies read')
            allocate(wk(0:nomax_gpr*j),omegap(0:nomax_gpr-1))
            rewind ifi
            if (rdm(ifi,0,nomax_gpr*j,' ',wk,nomax_gpr,j) /= 1) call rx('cannot read file!')
            forall (i = 0:nomax_gpr-1) omegap(i) = wk(i)/eVbyRy
            deallocate(wk)
          endif
        else if (strn(8:10)=='nw=') then
          ii = cmdoptswx('--pade','nw=','') + 2
          if (a2vec(strn,len_trim(strn),ii,2,dc//' ',2,3,1,iv,nomax_gpr) < 0) goto 911
          ii = cmdoptswx('--pade','window=','') + 6
          j = a2vec(strn,len_trim(strn),ii,4,', '//dc,3,2,2,iv,xv)
          allocate(omegap(0:nomax_gpr-1))
              omegap(0) = xv(1)
          do  iomg = 1, nomax_gpr-1
              omegap(iomg) = xv(1) + iomg*(xv(2)-xv(1))/(nomax_gpr-1)
          enddo
          if (j < 1) goto 911; if (nomax_gpr > 1 .and. j < 2) goto 911
        else if (strn(8:10)=='khm') then
          khm_usage = 'usage: --pade:khm:d=last_step:n=points_in_a_level:[l=levels|m=last_freq]'

          ii = cmdoptswx('--pade','d=','')
          if (ii<=0) call rx('"d" not passed, '//trim(khm_usage))
!           read(strn(ii+2:),*) khm_d
          j = a2vec(strn,len_trim(strn),ii+1,4,' '//dc,2,1,1,iv,khm_d)

          ii = cmdoptswx('--pade','n=','')
          if (ii<=0) call rx('"n" not passed, '//trim(khm_usage))
!           read(strn(ii+2:),*) khm_n
          j = a2vec(strn,len_trim(strn),ii+1,2,' '//dc,2,1,1,iv,khm_n)

          ii = cmdoptswx('--pade','l=','')
          khm_l = -1
!           if (ii>0) read(strn(ii+2:),*) khm_l
          if (ii>0) j = a2vec(strn,len_trim(strn),ii+1,2,' '//dc,2,1,1,iv,khm_l)


          ii = cmdoptswx('--pade','m=','')
!           if (ii>0) read(strn(ii+2:),*) khm_m
          if (ii>0) j = a2vec(strn,len_trim(strn),ii+1,4,' '//dc,2,1,1,iv,khm_m)

          if (ii>0 .and. khm_l /= -1) call rx('l and m are mutually exclusive, '//trim(khm_usage))
          if (ii<=0 .and. khm_l == -1) call rx('neither l nor m passed, '//trim(khm_usage))
          if (ii>0)  call khmesh(khm_n,khm_d, omegap, khm_l, khm_m)
          if (ii<=0) then
            call khmesh(khm_n,khm_d, omegap, khm_l)
            khm_m = -omegap(0)
          end if

          nomax_gpr = size(omegap)

C         Uses constructs ifort v11 does not recognize
C         open(newunit=khm_u,file='rmesh.dat',action='write')
C         write(khm_u,'(2(g0,x),i0,a)') khm_d, khm_m, khm_n,' # delta, ommax, Nd'
          khm_u = fopng('rmesh.dat',-1,8)
          call awrit3('%;16F %;16F %i  # delta, ommax, Nd',' ',80,khm_u,khm_d,khm_m,khm_n)
          do  iomg = 0, nomax_gpr-1
            call awrit2('%i %;16F',' ',80,khm_u,iomg,omegap(iomg))
            ! write(khm_u,'(i0,x,g0)') iomg, omegap(iomg) ! ifort v11 cannot compile
          end do
          close(khm_u)
          call info0(20,0,0,' Output frequency mesh written to rmesh.dat')

! this workaround is needed because the following routines work with Ry but then print out the results in eV
          omegap = omegap/eVbyRy
        endif
        ii = cmdoptswx('--pade','icut=','') + 4
        nlohi = [0,0]
        if (ii > 4) then
          if (a2vec(strn,len_trim(strn),ii,2,', '//dc,3,2,2,iv,nlohi) < 2) goto 911
        endif
        call supade(s_pade,0,nomg,nlohi,s_dmft%ndsig,omega,siginp)
C       Debugging check : efshfti used as a dummy
        call pade(1,dcmplx(0d0,omega(1)),s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof,efshfti)
        call info(20,0,0,' debugging check : Pade fit to first point:  err = (%,3g,%,3g)',
     .    siginp(1,1)-efshfti(1),siginp(1,2)-efshfti(2))
        allocate(sigf0(nomax_gpr,2,s_dmft%ndsig))
        do  i = 1, s_dmft%ndsig
          do  iomg = 0, nomax_gpr-1
            call pade(1,dcmplx(omegap(iomg),1d-6),s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),efshfti)
            sigf0(iomg+1,1:2,i) = efshfti(1:2)
          enddo
        enddo
        if (procid == master) then
          ifi = fopna(trim(fmt),-1,0)
          rewind ifi
          i = iodsig(201,nomax_gpr,0d0,s_dmft%ndsig,omegap,sigf0,-ifi)
        endif
        call rx0('wrote Pade-continued sigma to file '//trim(fmt))
  911   call rx('usage : --pade:nw=#[:icut#,#]:window=#,#')
      endif !pade



C --- Main loop to generate projectors ---
C     No LDA+U for now ... wait until redo with s_ldau
      s_ctrl%plbnd = 2          ! bndfp will not make output density
      if (s_ctrl%lwsig == LW8) then ! bndfp will not make evec file
!        s_ctrl%plbnd = 10       ! flags bndfp to quit after making potential ... needs checking
      else
        s_ctrl%lwsig = LW5      ! bndfp will write evals, evecs to file
      endif
      s_optic%loptic = 0        ! bndfp will not generate optics
      s_ctrl%lfp = s_ctrl%lfp - bitand(s_ctrl%lfp,2) ! Do not shorten q vectors
      nkp = s_bz%nkp
      allocate(s_ham%evals(ndham,nsp,nkp))
      allocate(s_pot%ppn(1))
      lrout = 0; lfrce = 0
      s_ctrl%lfp = s_ctrl%lfp - IAND(s_ctrl%lfp,4) + 4 ! Retain hab and sab
      s_ctrl%lfp = s_ctrl%lfp - IAND(s_ctrl%lfp,8) + 8 ! Update s_ctrl%elind
      if(.not. quick) then
      call bndfp(s_strn,s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot,s_bz,
     .  s_optic,s_gw,nbas,nsp,0,0,0,0,lrout,lfrce,0,xv,1,1,xv,xv,xv)
      ppnl => s_pot%ppn
      endif
      lrsig = s_ham%lsig
      nkfbz = s_bz%nkabc(1)*s_bz%nkabc(2)*s_bz%nkabc(3)

      if (lrsig /= 0 .and. s_dmft%dcmode /= 1) then
C        In future, if we need to hang onto sigma-vxc ...
C        s_hamx%nlmto = s_ham%nlmto
C        s_hamx%pwmode = s_ham%pwmode
C        s_hamx%ndhrs = s_ham%ndhrs
C        s_hamx%sigp = s_ham%sigp
C        s_hamx%eseavr = s_ham%eseavr
C        s_hamx%oveps = s_ham%oveps
C        s_hamx%ovncut = s_ham%ovncut
C        s_hamx%lrsa = s_ham%lrsa
C        s_hamx%iprmb => s_ham%iprmb
C        s_hamx%nprs => s_ham%nprs
C        s_hamx%iaxs => s_ham%iaxs
C        s_hamx%hrs => s_ham%hrs
        if (procid == master) ifi = fopna('sig0',-1,4)
        call rdsigm(lrsig,nbas,nsp,s_ham%nlmto,s_ctrl,s_spec,s_lat,s_ham,s_bz,s_gw,ifi,lwsig)
        call rxx(lwsig/=0,'sudmft: does not allow lwsig')
      endif

C --- Setup to make hybridization function or other quantities ---
      allocate(evlk(ndham,nsp),sigqp(ndham,ndham,nsp))
      call dpzero(sigqp,2*size(sigqp))
      allocate(dmftdelta(nomg,2*s_dmft%ndsig)) !hybridization written just for non-zero self-en elements
      allocate(dmftdeltah5(nomg,2*s_dmft%ndsig)) !hybridization written just for non-zero self-en elements
      allocate(g0h5(nomg,2*s_dmft%ndsig)) !hybridization written just for non-zero self-en elements
      nevn = nlohi(2)-nlohi(1)+1
      ef0 = s_bz%ef; if (s_bz%egap /= NULLI) ef0 = s_bz%ef + s_bz%egap/2
      nl = s_ctrl%nl

C --- Make static sigma from embedded static sig.inp ---
      if (cmdopt('--makesigqp',11,0,strn)) then
       dcmode = 0
       if (cmdopt('--makesigqp:mode=',17,0,strn)) then
        j1 = 0
        if (parg('mode=',2,strn(13:),j1,len(strn(13:)),' ',1,1,i,dcmode) <= 0)
     .    call rx('--makesigqp:mode= improperly formed')
       endif
C      Spin loop has to be done twice because the QSGW sigma file has it outside the k-loop
       if (procid == master) then
C       Sanity checks on embedded sig.inp
        ifi = fopng('sig.inp.f0.emb',-1,1); rewind ifi
        call readsiginp(1000+100*dcmode,ifi,i,j,k,0,0,qb,qb)
        if (i/=nkp .or. j/=nsp .or. k/=nevn)
     .    call rx('sudmft: file mismatch, sig.inp.f0.emb')
C       Write header to sigm file
        ifis = fopna('sigm1',-1,4); rewind ifis
        eseavr(1:2) = 0
        call iosigh(0,0,nsp,nspc,ndhamx,nlmto,s_bz%nkabc(1),s_bz%nkabc(2),s_bz%nkabc(3),nkp,0,0,0,0,-ifis,eseavr)
        allocate(zsp(ndhamx,ndhamx,2))
        do  isp = 1, nsp

C         Read header from evec file
          ifiz = fopna('evec',-1,4); rewind ifiz
          call iosigh(2,LW5,nsp,nspc,ndhamx,nlmto,nk1,nk2,nk3,nkp,nqi,lshft(1),lshft(2),lshft(3),ifiz,eseavr)
          if (ndhamx /= s_ham%ndham * nspc)
     .      call rx('sudmft: mismatch in evec file')
          rewind ifi  ! sig.inp
C     ... For each k and spin, read sig.inp make <z^-1 | sig.inp | z^-1>
          do  iq = 1, nkp
            do  i = 1, nsp
              read(ifiz) qp, ndimh
              ndimhx = ndimh * nspc
              call dpdump(evlk(1,i),ndhamx,ifiz)
              call dpdump(zsp(1,1,i),ndimhx**2*2,ifiz)
            enddo
            call readsiginp(1001+100*dcmode,ifi,nkp,nsp,nevn,nlohi(1),ndhamx,zsp,sigqp) ! sigqp for both spins
            write(ifis) qp
            call dpdump(sigqp(1,1,isp),ndimh**2*2,-ifis)
          enddo
        enddo
       endif
       call rx0('wrote embedded sigma (orbital basis) to file sigm1')
       deallocate(zsp)
      endif

C --- Make, normalize and save projectors in irreducible BZ ---
      allocate(s_dmft%Olapp(ldcix,ldcix,nsp,ncix))
      if(cmdopt('--stop=proj',11,0,strn)) then
         call lmlproj(111,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,ef0)
         call rx0('Hamiltonian save in h5')
      endif

      call info0(10,1,0,' Make and renormalize projectors ...')
      if(.not. quick) then
         ifi = fopna('proj',-1,4); rewind ifi ! Ensure proj file is rewound
         call lmlproj(11,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,ef0)
      endif
      call mpibc1(sigdc,2*size(sigdc),4,.false.,'','')




CC --- Make GW Sigloc ---
C      allocate(dmftu(nlohi(1):nlohi(2),ldcix,nsp,ncix))
C      ifi = fopna('proj',-1,4); rewind ifi
C      do  iqfbz = 1, nkfbz
C
CC   ... Read projector from disk
C        call iodmftu(s_dmft,.true.,ndham,nlohi(1),nlohi(2),ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifi)
C        if (ig == 1) iqs = 0  ! First point of new star
C        call iqstar(2,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,iv)
CC        call awrit6('%% rows %i cols %i  iq %,5i  irr %,5i  i123 %3,4i  qp %3:2,6;6d',' ',
CC     .    128,6,nomax_gpr,1+2*nsp*(2*s_dmft%l(1)+1),iqfbz,iq,iv,qpr)
C        print '(3i4,2x,3i3,3f12.6)', iq,ig,iqs,iv(1:3),qpr
C
CC        call dpzero(sigbark,2*size(sigbark))
CC        call embed_sigma(s_dmft,nlohi(1),nlohi(2),ldcix,nsp,nspx,dmftu,sigdc,siginp(iomg,:),sigbark)
CC
C
C      enddo
C      stop

C --- Find chemical potential or make DOS ---
      allocate(kpproc(0:nproc)) ! Distribution of processors for MPI
      kpproc(0) = 1; kpproc(1) = nomg+1 ! Serial case
      if (procid == master .and. nproc > 1) then
        allocate(N(nomg,1)); call dvset(N,1,nomg,1d0)
        call dstrbpx(nomg,nproc,N,kpproc)
        deallocate(N)
      endif
      call mpibc1(nkfbz,1,2,.false.,'','')
      call mpibc1(kpproc,nproc+1,2,.false.,'','')

      if (.not. s_dmft%lsigim) goto 100  ! Chemical potential finder works only for Matsubara frequencies now
      if (cmdopt('--mu=',5,0,strn)) then
        i = 0; i = parg('--mu=',4,strn,i,len(strn),', ',2,1,k,chempot)
        if (i < 0) call rxs('sudmft: failed to parse ',strn)
        call info2(10,1,0,' Reading chemical potential from command-line: %,6;6d',chempot,2)
        efshft = 0d0            ! 1010 block iteratively finds efshft where mu = ef0 - efshft
        goto 100
      endif
      call tcn('getmu')

      allocate(N(nkp,nsp),dNan(nkp,nsp))
      allocate(udeval(ndhamx,nomg,nkp,nspx))
      allocate(dmftu(nlohi(1):nlohi(2),ldcix,nsp,ncix))

      ldiaghbar = .true.        ! Require eigenvalues
      lrange = .false.          ! lrange will be set if DOS is sought after mu is found
      efshft = 0d0              ! 1010 block iteratively finds efshft where mu = ef0 - efshft
      rfwk(:) = 0d0
      minimization_required = .true.
C     ipass = controls flow for Fermi level finding in insulator mode
C     ipass =  0 => normal mode
C     ipass =  1 => 1st pass in insulator mode : find mu for c.n. point + s_bz%zinsul
C     ipass =  2 => 1st pass in insulator mode : find mu for c.n. point - s_bz%zinsul
      ipass =  0; if (s_bz%zinsul /= 0) ipass = 1
      call info2(10,1,0,' Find chemical potential mu%?#n# (insulator mode, dZ=%d)## ...',
     .  ipass,s_bz%zinsul)
      rfalsi_niter = 0
 1010 continue ! Reentry point for iterative charge neutrality search or DOS
      select case(ipass)
        case(0); qval = s_pot%qval
        case(1); qval = s_pot%qval + s_bz%zinsul
        case(2); qval = s_pot%qval - s_bz%zinsul
        case default; call rx('Wrong value for ipass')
      end select

      if (procid == master) then
        ifi = fopna('proj',-1,4); rewind ifi
        if (is0 == 1) then
         ifs0 = fopnx('sig.inp.f0.emb',2,0,-1); rewind ifs0
        endif
      endif
      call mpibc1(ifi,1,2,.false.,'','')

C ... Loop over full BZ
C     Note proj contains points within star for given iq, then over irr iq
      do  iqfbz = 1, nkfbz

C   ... Read projector from disk.  iodmftu returs iq, ig within i
        call iodmftu(s_dmft,.true.,ndham,nlohi(1),nlohi(2),ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifi)

C   ... Total DOS requires only irreducible kp
        if (ig == 1) then
          if (is0==1) then    ! just embedding sig.inp.f0
           allocate(sigbark(nevn,nevn,nspx))
           allocate(sigf0(2,s_dmft%ndsig,s_dmft%nicix))
           call read_sigf0(s_dmft%ndsig,s_dmft%nicix,sigf0)
           call embed_sigma(s_dmft,nlohi(1),nlohi(2),ldcix,nsp,nspx,dmftu,sigdc,sigf0,sigbark)
           call print_sigf0(nevn,nsp,sigbark,iqfbz,ifs0)
           deallocate(sigf0,sigbark)
          else ! working with sig.inp
           if (ldiaghbar .or. ldpade) then
            allocate(sigbark(nevn,nevn,nspx),udhamk(ndhamx,ndhamx,nspx),udevalko(ndhamx))
            allocate(sigpade(2,s_dmft%ndsig,s_dmft%nicix))
            do  isp = 1, nspx
              call dpzero(udeval(1,1,iq,isp),2*ndhamx*nomg)
            enddo
            do  iomg = 1, nomg
              if (iomg<kpproc(procid) .or. iomg>=kpproc(procid+1)) cycle

c         ... Construct Sigmabar(1:nev,1:nev,1:nsp) in the lattice space
              call dpzero(sigbark,2*size(sigbark))
              call mksigpade(s_pade,1,nomg,s_dmft%ndsig,dcmplx(emu,omega(iomg)),siginp(iomg,1),sigpade)
              call embed_sigma(s_dmft,nlohi(1),nlohi(2),ldcix,nsp,nspx,dmftu,sigdc,sigpade,sigbark)

!             if (iqfbz==1) write(*,'(i4,f12.6,12f14.6)') iomg, omega(iomg), (dble(sigbark(i,i,1)), i=31,40)

C         ... Construct updated hamiltonian (at each k-point and frequency)
              call dpzero(udhamk,2*size(udhamk))
              call makehbar2(nevn,nlohi(1)-1,nspx,ndhamx,evlk,sigbark,udhamk)
C         ... Diagonalize updated hamiltonian at each k-point and frequency
              call dpzero(udevalko,2*size(udevalko))
              do  isp = 1, nspx
                call diagoham(ndhamx,udhamk(1,1,isp),udevalko)
                udeval(:,iomg,iq,isp) = udevalko(:)
              enddo
            enddo               ! loop over iomg
C           print "('checksum udeval',1p4e18.8)",sum(udeval(:,:,iq,1)),sum(udeval(:,:,iq,nspx))

C           Combine udeval contributed by various processors
            do  isp = 1, nspx
              call mpibc2(udeval(1,1,iq,isp),2*ndhamx*nomg,4,3,.false.,'','')
            enddo
            deallocate(udhamk,udevalko,sigbark,sigpade)
           endif   ! ldiaghbar

C      ... Compute sum_over frequency-dependent contribution to the electron number for a given efshft
           N(iq,:)=0 ; dNan(iq,:)=0
           do  iomg = 1, nomg
             call cmpvalcharg_matsub4(ndhamx,nspx,omega(iomg),omega(1)/pi,
     .       udeval(:,iomg,iq,:),udeval(:,nomg,iq,:),efshft,ef0,N(iq,:),dNan(iq,:),'sumfrq')
           enddo
C      ... Analytic contribution
           call cmpvalcharg_matsub4(ndham,nsp,xv,omega(1)/pi,
     .     xv,udeval(:,nomg,iq,:),efshft,ef0,N(iq,:),dNan(iq,:),'fermif')
C          Weight by number of irreducible qp in the full BZ
           do  isp = 1, nspx
             N(iq,isp)    = N(iq,isp)   *s_bz%wtkp(iq)/nsp
             dNan(iq,isp) = dNan(iq,isp)*s_bz%wtkp(iq)/nsp
           enddo
          endif                 ! sig.inp or sig.inp.f0
        endif                   ! ig == 1
      enddo                     ! k-points in the full BZ
      if (is0==1) call rx0('File sig.inp.f0 embedded, and recorded in sig.inp.f0.emb')
      ldiaghbar = .false.       ! Eigenvalues now calculated

C ... Integrated charge
      do  isp = 1, nspx
        Nnew(isp)     = sum(N(:,isp))
        dNan_new(isp) = sum(dNan(:,isp))
      enddo

C ... Branch to tabulate total DOS on a mesh
      if (lrange) then
        if (nspc == 2) call rx('total DOS not implemented for noncollinear case')
        j1 = j1+1               ! Counter for current energy shift
        do  isp = 1, nsp
          dos(j1,isp) = -dNan_new(isp); if (lidos) dos(j1,isp) = Nnew(isp)
        enddo

        efshft = ef0 - dosw(1) - dble(j1)/dble(npts-1)*(dosw(2)-dosw(1))
        emu = dosw(1) + dble(j1)/dble(npts-1)*(dosw(2)-dosw(1)) - chempot

        if (mod(j1,10) == 1) then
          call info2(30,0,0,' ... finished energy point %i (%;0d sec)',j1,cpusec()-tim)
          tim = cpusec()
        endif
        if (j1 == npts) then
          i = 3; if (ldig) i = i+100
          if (procid == master) then
          call iodos(i,-fopnn('DOS'),dos,npts,nsp,npts,1,dosw(1)-dosw(4),dosw(2)-dosw(4),nsp,ef0-dosw(3)-dosw(4),idfmt)
          endif
          call rx0('done writing total DOS')
        endif
        goto 1010
      endif

      if (rfalsi_niter == 0) call info5(20,0,0,' Seek mu for %d electrons ... '//
     .  '%d electrons at Ef0=%d.  D(Ef0)=%;3d',qval,sum(Nnew),ef0,-sum(dNan_new),0)
C ... rfalsi returns an estimate for efshft.  Iterate until convergence
      if (minimization_required.and.rfalsi_niter<=rfalsi_limit) then
        call getVnew(efshft,qval,sum(Nnew),sum(dNaN_new),minimization_required,rfwk,rfalsi_niter)
        if (rfalsi_niter==rfalsi_limit) write(*,111)
  111   format(' + PREVENT INFNITE LOOP . Reliable Vnew not found')
        goto 1010
      endif
      call info5(20,0,min(ipass,1),
     .  ' mu%?#n==1#(+)##%-1j%?#n==2#(-)## = %d = Ef0-%d.  '//
     .  'Deviation from neutrality = %;3e  D(mu)=%;3,3d',
     .  ipass,ef0-efshft,efshft,sum(Nnew)-s_pot%qval,-sum(dNan_new))
      if (nspx == 2) call info5(20,0,0,
     .  ' Electron charge:  %;6,6d   moment:  %;6,6d spin resolved DOS:  %2:2;3,3d',
     .    Nnew(1)+Nnew(2),Nnew(1)-Nnew(2),-dNan_new,0,0)
      chempot = ef0-efshft

C     Flow control for insulator mode
      select case(ipass)
      case(1)
        efshfti(1) = efshft
        if (procid == master) then
          ifi = fopna('proj',-1,4); rewind ifi
        endif
        ipass = 2; rfalsi_niter = 0; efshft = 0; minimization_required = .true.
        goto 1010
      case(2)
        efshfti(2) = efshft; efshft = (efshfti(1)+efshfti(2))/2
        call info5(20,0,0,
     .    ' mu(+) = %d   mu(-) = %d   mubar = %d   gap = %d Ry = %d eV',
     .    ef0-efshfti(1),ef0-efshfti(2),ef0-efshft,efshfti(1)-efshfti(2),(efshfti(1)-efshfti(2))*eVbyRy)
        ipass = 0
        chempot = ef0-efshft
      end select

C ... Option to generate DOS
      if (cmdopt('--dos',5,0,strn)) then
        call iinit(dosisw,size(dosisw)); idfmt = 1
        call sudossw(strn,1+2+4+16,doslsw,dosisw,dosw,' ')
        dosw(4) = 0; if (doslsw(8)) dosw(4) = chempot ! Shift Fermi level to 0 when writing DOS
        dosw(1) = dosw(1) + dosw(4); dosw(2) = dosw(2) + dosw(4)
        if (.not. lrange .or. npts == 0)
     .    call rxs('sudmft: missing or improper dos options: ',strn)
      endif

C ... Setup to generate DOS
      if (lrange) then
        if (nspc == 2) call rx('total DOS not implemented for noncollinear case')
        tim = cpusec()
        ldpade = doslsw(9)
        dosw(3) = efshft
C        if (doslsw(8)) then
C          dosw(1) = dosw(1) + ef0-dosw(3); dosw(2) = dosw(2) + ef0-dosw(3)
C        endif
        call info8(20,1,0,' Make DOS for %i points from %i bands in'
     .  //' window (%;3d,%;3d) = mu + (%;3d,%;3d)  %?#n#with Pade interp##',
     .    npts,nevn,dosw(1)-dosw(4),dosw(2)-dosw(4),dosw(1)-chempot,dosw(2)-chempot,isw(ldpade),8)
        efshft = ef0 - dosw(1)  ! energy relative to noninteracting Fermi level
        emu = dosw(1)-chempot   ! -omega on real axis (negative of energy relative to mu)
        if (.not. ldpade) call info0(10,1,0,' *** Warning : You are calculating DOS without Pade interpolation!')

C       Pade coefficients DOS to be calculated with varying Re[omega]
        if (ldpade) then
          call supade(s_pade,0,nomg,[1,0],s_dmft%ndsig,omega,siginp)
C          Check
C           npade = s_pade%npade
C           call pade(1,dcmplx(0d0,omega(100)),npade,s_pade%zomega,npade,s_pade%padcof,xv)
C           print *, siginp(100,1,1) - xv(1), siginp(100,2,1) - xv(2)
C           stop
        endif
        j1 = 0  ! Counter for omega loop
        allocate(dos(npts,nsp))
        goto 1010
      endif
      deallocate(udeval,dmftu,N,dNan)
      call tcx('getmu')
  100 continue                  ! End of branch to get chemical potential

C --- Compute the charge density ---
C     Possibly some bugs.  Ni works fine (no sigma)
C     1. lscoq test case has strange weights, e.g. dmft weights 1 up to band 81 !?
C        Possibly range nlohi not properly taken into account
C     2. lscoq test case : MPI doesn't exactly tally with serial
      if (cmdopt('--udrs',6,0,strn)) then ! update the rs file.
        if (nspc == 2) call rx('--udrs not ready yet, noncollinear case')
        call tcn('getrho')
        call info0(10,1,0,' Compute updated density ...')
!       open relevant files
        ndhamx = s_ham%ndham * nspc

        if (procid == master) then
          ifi  = fopna('proj',-1,4); rewind ifi  ! projectors
          ifiz = fopna('evec',-1,4); rewind ifiz ! QSGW wavefunctions
          call iosigh(2,LW5,i,j,ndhamx,nlmto,nk1,nk2,nk3,nkp,nqi,lshft(1),lshft(2),lshft(3),ifiz,eseavr)
        endif
        if (ndhamx /= s_ham%ndham * nspc) call rx('sudmft : mismatch in evec file (update rho)')
        call mpibc1(i,1,2,.false.,'','')
        call mpibc1(j,1,2,.false.,'','')
        call mpibc1(nkabc,3,2,.false.,'','')
        call mpibc1(nkp,1,2,.false.,'','')
        call mpibc1(lshft,3,2,.false.,'','')

        lrout = 1; numq = 1
        lekkl = 0 ! Methfessel style band CG since 1-particle evals are not meaningful
        call dfqkkl(nbas,lekkl,s_site,s_spec,numq)
        call ptr_pot(s_pot,8+1,'smrout',k1*k2*k3,numq*nsp*nspc,0)
        srout => s_pot%smrout
        call surho(nbas,s_site,s_spec,0,0,lrout,lekkl,numq,k1,k2,k3,srout,0,xv,sumev,sumqv)

        allocate(udeval(ndham,nomg,nkp,nsp))
        call dpzero(udeval,2*size(udeval))
        allocate(dmftu(nlohi(1):nlohi(2),ldcix,nsp,ncix))
        allocate(evecr(ndham,ndham,nomg,nsp))
        allocate(evecl(ndham,ndham,nomg,nsp))
        allocate(dmevec(ndham,ndham))
        call dpzero(dmevec,2*size(dmevec))
        allocate(dmevecinv(ndham,ndham))
        allocate(dmeval(ndham,nsp,nkp))
        call dpzero(dmeval,size(dmeval))
        allocate(denmat(ndham,ndham,nsp))
        allocate(dmatk(ndham*nspc))
        allocate(f(1))
        allocate(tso(1))
        allocate(dmftvec(ndham,ndham,nsp))
        allocate(dmftwktb(ndham,nsp,nkp))

        allocate(evalq(ndham,nsp,nkp))
        allocate(evecq(ndham*nspc,ndham*nspc,nsp,nkp))
        call readevec(2,ndham,nsp,nspc,nkabc,nkp,lshft,nproc,evalq,evecq)

C       Loop over irreducible qp
        do  iqfbz = 1, nkfbz
         call iodmftu(s_dmft,.true.,ndham,nlohi(1),nlohi(2),ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifi)

         if (ig /= 1) cycle
C        if (iq == 3) call rx0('done')

         call dpzero(evecl,2*size(evecl)) ! Probably not needed
         call dpzero(evecr,2*size(evecr)) ! Probably not needed
         call dmftevec(s_dmft,s_pade,ndham,nomg,nsp,nspc,nlohi,ldcix,iq,nkp,
     .     procid,kpproc,omega,dmftu,evlk,emu,sigdc,siginp,evecl,evecr,udeval)

C         iomg = kpproc(procid); isp = 1
C         call info5(1,0,0,'Procid,iq,isp,iomg %4:1,4i evecl %3;11,6D evecr %3;11,6D udeval %3;11,6D',
C     .     [procid,iq,isp,iomg],sum(evecl(:,:,iomg,isp)),sum(evecr(:,:,iomg,isp)),sum(udeval(:,iomg,iq,isp)),5)


C        Last frequency must be distributed to all processors
         do  isp = 1, nspx
           call mpibc2(udeval(1,nomg,iq,isp),2*ndhamx,4,3,.false.,'','')
           call mpibc2(evecl(1,1,nomg,isp),2*ndhamx*ndhamx,4,3,.false.,'','')
           call mpibc2(evecr(1,1,nomg,isp),2*ndhamx*ndhamx,4,3,.false.,'','')
         enddo
C        End of the diagonalization. Right and left eigenvectors are used
C        to construct the density (actually the DMFT rotation matrices)

C   ...  Construct the (static) density matrix denmat (sum over Matsubara's)
         call dpzero(denmat,2*size(denmat))
         do  iomg = 1, nomg
           if (iomg<kpproc(procid) .or. iomg>=kpproc(procid+1)) cycle
           do  ib1 = 1, ndhamx
           do  ib2 = 1, ndhamx
             do  isp = 1, nspx

               call cmp_denmat(ndhamx,omega(iomg),omega(1)/pi,
     .                         udeval(1,iomg,iq,isp),evecl(1,ib2,iomg,isp),evecr(1,ib1,iomg,isp),
     .                         udeval(1,nomg,iq,isp),evecl(1,ib2,nomg,isp),evecr(1,ib1,nomg,isp),
     .                         chempot,denmat(ib1,ib2,isp),'sumfrq')
               if (iomg == nomg) then
               call cmp_denmat(ndhamx,NULLR,omega(1)/pi,
     .                         NULLZ,NULLZ,NULLZ,
     .                         udeval(1,nomg,iq,isp),evecr(1,ib2,nomg,isp),evecr(1,ib1,nomg,isp),
     .                         chempot,denmat(ib1,ib2,isp),'fermif')
               endif

             enddo ! spin
           enddo   ! ib2
           enddo   ! ib1
          enddo ! freqs

          call mpibc2(denmat,2*size(denmat),4,3,.false.,'','')

C          if (procid == 0) then
C            print *, 'writing denmat...'
C            ifi = fopna('denmat',-1,0)
C            call ywrm(0,' ',3,ifi,'(9f18.11)',denmat,1,ndhamx,ndhamx,ndhamx)
C          endif
C          call rx0('done')

C     ... Diagonalise denmat and get DMFT rotations (dmevec) and weights (dmeval)
          if (procid == master) then
          do  isp = 1, nspx

C            Use this if density matrix is non hermitian.  Use complex cast for dmeval
C            call diagdenmat(ndhamx,denmat(1,1,isp),dmevec,
C     .                      dmevecinv,dmevalz(1,isp,iq))
C             ! Useful subroutine for check.
C            !call testdenmat(ndhamx,dmevec,denmat(:,:,isp),dmeval(:,isp,iq),isp,iqfbz)
C
CC           Read QSGW wavefunctions and rotate them with mkdmftstates
C            read(ifiz) qp, ndimh   ! file is opened before the loop on q points
C            ndimhx = ndimh * nsp
C            call dpdump(evlk(1,isp),ndhamx,ifiz)
C            call dpdump(zsp(1,1,isp),ndhamx**2,ifiz)
C            call mkdmftstatesz(ndhamx,s_bz%wtkp(iq)/nsp,dmevec,
C     .                        zsp(1,1,isp),dmftvec(1,1,isp),
C     .                        dmevalz(1,isp,iq),dmftwktb(1,isp,iq),nspc)
C            call info2(1,1,0,'eigenvalues of density matrix spin %i:',isp,iq)
C            print '(i4,2f12.6,2x,i4,2f12.6)', ((i, dmeval(i,isp,iq), i+1, dmeval(i+1,isp,iq)), i=1,ndhamx,2)

            call zhevx(ndhamx,ndhamx,denmat(1,1,isp),xv,0,.true.,ndhamx,-NULLR,
     .        i,dmevecinv,.false.,dmeval(1,isp,iq),ndhamx,dmevec)

C           call zprm('evec from zhev',2,dmevec,ndhamx,ndhamx,ndhamx)

C           Debugging
C            call info2(1,1,0,'eigenvalues of density matrix spin %i:',isp,iq)
C            print '(4(i4,f12.6))', ((i+k, dmeval(i,isp,iq), k=0,3), i=1,ndhamx,4)

C           Rotate QSGW wavefunctions with mkdmftstates
            ndimh  = ndhamx        ! Assume all LMTOs for now
            call mkdmftstates(ndhamx,s_bz%wtkp(iq)/nsp,dmeval(1,isp,iq),dmevec,evecq(1,1,isp,iq),
     .        dmftvec(1,1,isp),dmftwktb(1,isp,iq))

C           Print out density weights, single-particle and interacting cases
C           Note : there is not a one-to-one correspondence QSGW/LDA and DMFT states
            call print_wtkb(iq,iqfbz,isp,ndhamx,nsp,nkp,ndhamx,s_bz%wtkp,chempot,evlk,s_bz%wtkb,dmftwktb)

C           Accumulate the output density
            s_bz%nef = 1
            call addrbl(s_site,s_spec,s_lat,s_ham,s_bz,s_optic,isp,nsp,
     .                  nspc,qp,ndham,ndimh,NULLI,lrout,4,0,0,dmftwktb,-1,s_bz%swtk,
     .                  iq,lfrce,0,0,k1,k2,k3,s_pot%smpot,s_pot%smvcnst,s_pot%lcplxp,numq,
     .                  qval,dmftvec(1,1,isp),evlk(1,isp),ndham,
     .                  chempot,0,esmear,NULLR,NULLR,NULLI,xv,
     .                  srout,sumqv,sumev,f,dmatk,tso)

          enddo     ! Loop over spins
          endif     ! master process

        enddo       ! Full BZ

        deallocate(dmftu,udeval,evecr,evecl,denmat)
        if (procid == master) then
          call fclose(ifi)
          call fclose(ifiz)
        endif

C   ... Restore q=0 lattice vectors
        call pshpr(10)
        call sugvec0(s_lat)
        call poppr

C   ... Make local parts of output density
        call dfratm(s_site,s_spec,16+2,1,nbas,s_pot%rnew)
        allocate(qbyl(n0*nsp*nbas))
        allocate(hbyl(n0*nsp*nbas))
        call mkrout(' ',s_site,s_spec,s_lat,s_ham,nbas,nsp,ldim,0,
     .    s_pot%rnew,s_pot%hab,s_pot%sab,qbyl,hbyl,lrout)

C   ... Symmetrize output density
        call symrho(s_site,s_spec,s_lat,7+10*lfrce,s_pot%smrout,s_pot%rnew,qbyl,hbyl,f)
C       call zprm3('smrho after symrho',0,srout,k1,k2,k3*nsp)
        if (iand(s_ctrl%lcd,256) /= 0) then
          i = 1
          call rhopos(i,s_site,s_spec,s_lat,s_pot%rnew,s_pot%smrout)
        endif

C   ... Mix input and output densities
        if (lpnu > 0) then
          xv(1:10) = s_ham%pmin; xv(11:20) = s_ham%pmax
          i = s_ham%pnudef
          if (i > 40) i = i-40  ! Strip off flag that retains lvs11 bug
          lfrzw = isw(IAND(s_ctrl%lbas,16)/=0)
          call pnunew(nbas,nsp,s_site,s_spec,xv,xv(11),10*i+lfrzw,
     .      s_pot%hab,s_pot%sab,qbyl,hbyl)
        endif

C   ... Mix input and output densities
        i = str_pack('mix',-2,s_strn,strn)
        call mixrho(s_site,s_spec,s_lat,s_pot,srout,s_pot%smrho,
     .    nsp,iter,strn,qval,s_ham%elind,k1,k2,k3,dmxp)
        dmxp(26) = iter; dmxp(27) = dmxp(11)

        deallocate(dmevec,dmeval,dmevecinv,srout,dmatk,tso,f,dmftvec,dmftwktb)

        if (procid == master .and. s_ctrl%quit /= 32) then
          call info2(10,1,0,' Writing restart file ... RMS DQ is now %1,3;3e',dmxp(27),0)

          irs(1) =     IAND(s_ctrl%lrs,7)
     .           +     8*isw(IAND(s_ctrl%lrs,256)/=0)
          irs(2) =     IAND(s_ctrl%lrs,8+16)/8 + 10*isw(IAND(s_ctrl%lrs,512)/=0)
          irs(3) = isw(IAND(s_ctrl%lrs,32)/=0)
          irs(4) = isw(IAND(s_ctrl%lrs,64)/=0)
          irs(5) = isw(IAND(s_ctrl%lrs,128)/=0)

          lbin = mod(irs(2),10) /= 2
          if (lbin) fileid = 'rst'
          if (.not. lbin) fileid = 'rsta'
          if (mod(irs(2),10) == 3) then
            call word(fileid,1,i,j)
            j = j+1
            fileid(j:j) = '.'
            call bin2a(' ',0,0,iter,2,0,len(fileid),fileid,j)
            if (lbin) ifi = fopng(fileid,-1,8+4)
            if (.not. lbin) ifi = fopng(fileid,-1,8)
            call info0(10,1,-1,' sudmft:  writing to restart file '//fileid)
          else
            if (lbin) ifi = fopna(fileid,-1,4)
            if (.not. lbin) ifi = fopna(fileid,-1,0)
          endif
          i = str_pack('jobid',-1,s_strn,strn)
          fileid = 'sudmft: ' // strn
          nat = nglob('nat')
          k = iorsf(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .      fileid,nbas,nat,s_ctrl%nspec,xv,iter,lbin,-ifi)
          call fclose(ifi)
        endif                     ! Master process writing restart file

        call rx0('lmdmft: rst file written')
      endif                     ! branch to generate output density

C --- Check for switches make Gloc (total and/or k-resolved) or G or Sigma ---
      i = 7
      nqbnd = 0                 ! Use the existing k-mesh for k list
      if (cmdopt('--gprt',i-1,0,strn)) then
        gpmode = 1              ! Default in the absence of mode=
        nomax_gpr = 0           ! Flags that it has not been set
        dc = strn(7:7)
        if (dc /= ' ') then
          do
           i = i+1  ! i points to start of next switch.  Preserve i in this loop!
           if (strn(i:i) == dc) cycle
           j = min(len(strn),i)
           call nwordg(strn,0,dc//' ',1,j,i) ! i now terminates switch and must be preserved
           if (i < j) exit
           if (.false.) then
C          Select what to calculate.  For documentation, see command-line switches
           elseif (strn(j:j+4) == 'mode=') then
            ii = j+4
            if (ii >= i) call rx('gprt: bad argument to mode')
            jj = a2vec(strn,i,ii,2,dc//' ',2,3,1,iv,gpmode)
            if (jj/=1) call rxs('gprt, cannot parse, ',strn(j:i))
            if (gpmode == 4) gpmode = -2
C          Override calculated chemical potential
           elseif (strn(j:j+2) == 'mu=') then
            ii = j+2
            if (ii >= i) call rx('gprt: bad argument to mode')
            jj = a2vec(strn,i,ii,4,dc//' ',2,3,1,iv,chempot)
            if (jj/=1) call rxs('gprt, cannot parse, ',strn(j:i))
C          Read in new k list
           elseif (strn(j:j+3) == 'band') then   ! 1=local, 2=k-resolved, 3=both, 4=k-resolved h, 8=full gk
             fn = strn(j+4:i)
C            Return nfbn,ifblst
             deallocate(ifblst); allocate(ifblst(ndhamx*3))
             call rdqlst(0,s_lat,trim(fn),onesp,nqbnd,xv,xv) ! count number of k-points
             if (nqbnd <= 0) call rx('improper qp list from '//strn(j:i))
C          Read in new self-energy
           elseif (strn(j:j+6) == 'rdsigr=' .or. strn(j:j+5) == 'rdsig=') then   ! Read omega, sigma from file
             k = 7; if (strn(j:j+5) == 'rdsig=') k = 6
             if (j+k > i) call rx('gprt: bad argument to rdsig')
             call info2(20,0,0,
     .         ' Read new frequency mesh (%?#(n==7)#real#imaginary# axis) '//
     .         'and sigma from file '//strn(j+k:i),k,2)
             if (k == 7) then  ! Real sigma
               s_dmft%beta = 0
               omfac = 1
               s_dmft%lsigim = .false.
               s_dmft%gammac = 0 ! fix makegloc2 before using broadening on real axis
C              Set broadening for r.s. G
               ii = wordsw(strn,dc,'broad=','',jj)
               if (ii /= 0) then
                 jj = jj-1
                 if (a2vec(strn,len_trim(strn),jj,4,', '//dc,3,3,1,iv,s_dmft%gammac) /= 1)
     .             call rx('gprt: bad argument to broad')
               endif
             endif

             if (procid == master) l = fxst(strn(j+k:i)); call mpibc1(l,1,2,.false.,'','')
             if (l /= 1) call rxs('--gprt failed to find missing file: ',strn(j+k:i))
C            Get new dimensions
             if (procid == master) ifi = fopna(strn(j+k:i),-1,1)
             jj = iodsig(151,nomg,s_dmft%beta,l,s_dmft%omg,s_dmft%sig,ifi)
             if (jj /= 0) call rx('failed to read sigm file')
             if (l /= s_dmft%ndsig) call rxs('mismatch in DMFT channels, file ',strn(j+k:i))
             call mpibc1(nomg,1,2,.false.,'','')
             s_dmft%nomg = nomg
             if (associated(s_dmft%omg)) deallocate(s_dmft%omg)
             if (associated(s_dmft%sig)) deallocate(s_dmft%sig)
             allocate(s_dmft%omg(nomg),s_dmft%sig(nomg,2*s_dmft%ndsig))
             jj = iodsig(131,s_dmft%nomg,s_dmft%beta,s_dmft%ndsig,s_dmft%omg,s_dmft%sig,ifi)
             if (procid == master) call fclose(ifi)
             omega => s_dmft%omg; siginp => s_dmft%sig

C            Update distribution of omegas
             kpproc(0) = 1; kpproc(1) = nomg+1 ! Serial case
             if (procid == master .and. nproc > 1) then
               allocate(N(nomg,1)); call dvset(N,1,nomg,1d0)
               kpproc(2:) = 0
               call dstrbpx(nomg,nproc,N,kpproc)
               deallocate(N)
             endif
             call mpibc1(kpproc,nproc+1,2,.false.,'','')
C            if (nproc > nomg) call rx2('sudmft : nproc (%i) must not be larger than nomg (%i)',nproc,nomg)


C          Set number of frequency points
           elseif (strn(j:j+3) == 'nom=') then ! print k-resolved
             if (j+4 > i) call rx('gprt: bad argument to nom')
             ii = j+3
             jj = a2vec(strn,i,ii,2,dc//' ',2,3,1,iv,nomax_gpr)
           else
            call rxs('gprt: failed to parse argument, ',strn(j:i))
           endif
          enddo
          if (nomax_gpr==0) nomax_gpr = nomg
          if (nomax_gpr>nomg) nomax_gpr = nomg
        endif

        if (iprint() >= 20 .and. procid == master) then
          if (mod(gpmode,2)==1 .and. gpmode<4) then
            write(*,'(x,a)') 'k-integrated gloc'
          else if (gpmode==2 .or. gpmode==3) then
            write(*,'(x,a)') 'k-resolved gloc'
          else if (gpmode==-2) then
            write(*,'(x,a)') 'k-resolved hloc'
          else if (gpmode==8) then
            write(*,'(x,a)') 'G(k) in 1-particle basis'
          else if (gpmode==9) then
            write(*,'(x,a)') 'Sigma(k) in 1-particle basis'
          else if (gpmode==18) then
            write(*,'(x,a)') 'diagonal G(k) in 1-particle basis'
          else if (gpmode==19) then
            write(*,'(x,a)') 'diagonal Sigma(k) in 1-particle basis'
          else if (gpmode==20) then
            write(*,'(x,a)') 'effective diagonal Sigma(k) in 1-particle basis'
          else if (gpmode==1) then
            write(*,'(x,a,f6.4)') 'Im#, Real# axis,  mu=',chempot
          else
            write(*,'(x,a,i0)') 'band mode, ',nqbnd,' qp'
          end if
        end if
!         call info5(20,1,0,' Create and write to disk: '//
!      .    '%?#(n%2==1&n<4)# k-integrated gloc##%-1j'//
!      .    '%?#(n==2|n==3)# k-resolved gloc##%-1j'//
!      .    '%?#(n==-2)# k-resolved hloc##%-1j'//
!      .    '%?#(n==8)# G(k) in 1-particle basis##%-1j'//
!      .    '%?#(n==9)# Sigma(k) in 1-particle basis##%-1j'//
!      .    '%?#(n==18)# diagonal G(k) in 1-particle basis##%-1j'//
!      .    '%?#(n==19)# diagonal Sigma(k) in 1-particle basis##%-1j%j'//
!      .    '%?#(n==20)# effective diagonal Sigma(k) in 1-particle basis##%-1j%j'//
!      .    '%?#(n==1)#, Im#, Real# axis,  mu=%,6;6d'//
!      .    '%?#(n)# (band mode, %-1j%i qp)',
!      .    gpmode,isw(s_dmft%lsigim),chempot,nqbnd,5)

        if (gpmode < -2 .or. gpmode > 20 .or.
     .      gpmode > 5 .and. gpmode < 8 .or.
     .      gpmode > 9 .and. gpmode < 18)
     .    call rxi('unexpected value for --gprt mode',gpmode)
        if (nqbnd > 0) then
          deallocate(s_bz%qp)
          s_bz%nkp = nqbnd
          nkp = nqbnd
          nkfbz = nqbnd
          lstar = .false.
          allocate(s_bz%qp(3,nqbnd))
          call rdqlst(3,s_lat,trim(fn),onesp,nqbnd,s_bz%qp,nkpan)
C         Get color weight, if any
          nfbn = 0
          call suqlst(s_lat,fn,7,ndhamx,0d0,nsp,xv,nfbn,ifblst,0,qp,onesp)
          deallocate(s_ham%evals)
          allocate(s_ham%evals(ndham,nsp,nkp))
        endif
      else
        gpmode=0 ; nomax_gpr=nomg
      endif ! gprt switch

C --- Create local quantities Delta, Eimp, possibly G or Sigma, depending on gpmode ---
      call tcn('delta')

C      optional print of sigbar
C      if(iprint() >= PRTG) then
C        allocate(chk_projsbark(ldcix,ldcix,nsp,ncix))
C        allocate(chk_projsbar(ldcix,ldcix,nsp,ncix))
C        call dpzero(chk_projsbar,2*size(chk_projsbar))
C      endif

C --- Change k mesh ---
      if (nqbnd > 0) then
        if (mod(gpmode,2)==1 .and. gpmode<4)
     .    call rx('cannot use band mode for k-integrated properties')

C   ... Eigenvectors of 1-particle hamiltonian on for list of k points
        ifiz = fopna('evec',-1,4); call fclose(ifiz)
        s_ctrl%plbnd = 3        ! bndfp will not integrate over BZ
        s_bz%nkabc = 0
C       bndfp will need to replace nkp, s_bz%qp
        call bndfp(s_strn,s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot,s_bz,
     .  s_optic,s_gw,nbas,nsp,0,0,0,0,lrout,lfrce,0,xv,1,1,xv,xv,xv)
        ppnl => s_pot%ppn
        lrsig = s_ham%lsig

C   ... Make the projectors
        call info0(10,1,0,' Make and renormalize projectors for new qp ...')
        ifi = fopna('proj',-1,4); rewind ifi ! Ensure proj file is rewound
        i = 0 ; if (s_dmft%knorm == 0) i = 20
C       print 777, 'sudmft call lmlproj',procid
        call lmlproj(i,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,ef0)
C       print 777, 'sudmft exit lmlproj',procid
      endif

C --- Setup including gpmode-specific parts ---
      allocate(dmftu(nlohi(1):nlohi(2),ldcix,nsp,ncix))


      call bzmsh00(s_lat%plat,s_bz%lshft,0,s_bz%nkabc,ifac,qb) ! Needed for iqstar
      nkabc = s_bz%nkabc
C ... default mode
      if (gpmode == 0) then
        call info0(10,1,0,' Make gloc, delta, eimp ...')
        allocate(eimpm(ldcix,ldcix,nsp,ncix))  ! Assumed to be collinear for now
        call dpzero(eimpm,2*size(eimpm))
      endif

      ifgk = 0; wgpmode = 0
      if (gpmode >= 18 .and. gpmode <= 20) then
        if (procid == master) then
          ifgk = fopna('se',-1,2); rewind ifgk
          allocate(iblst(nlohi(1):nlohi(2)))
          forall (i = nlohi(1):nlohi(2)) iblst(i) = i
          iseh = 100*0 + 20 + 4 ! in Ry, ASCII, write ommin, ommax, chemical potential
!         j = 1
          if (s_dmft%lsigim) then
            j = -nomax_gpr
            xv(1) = s_dmft%beta*eVbyRy
            xv(2) = 0
          else
            j = nomax_gpr
            xv(1) = s_dmft%omg(1)
            xv(2) = s_dmft%omg(nomax_gpr)
          endif
          call ioseh(iseh,0,nsp,nspx,nlohi(2)-nlohi(1)+1,nkp,j,xv(1),xv(2),chempot,
     .      nlohi(2)-nlohi(1)+1,iblst,0,-ifgk)
          deallocate(iblst)
        endif                     ! procid = master
      elseif (gpmode == 9) then
        if (nproc > 1) call rx1('--gprt mode=%i does not yet work with MPI, sorry',gpmode)
        wgpmode = gpmode
        fn = 'sigk'; lascii = 4      ! Out file name, binary mode
      elseif (gpmode == 8) then
        if (nproc > 1) call rx1('--gprt mode=%i does not yet work with MPI, sorry',gpmode)
        wgpmode = gpmode
        fn = 'gk'; lascii = 4        ! Out file name, binary mode
      elseif (mod(gpmode,8) > 1) then
        wgpmode = gpmode
        fn = 'gkloc'; lascii = 0        ! Out file name, ascii mode
      elseif (gpmode == -2) then
        wgpmode = gpmode
        fn = 'hkloc'; lascii = 0     ! Out file name, ascii mode
      endif
C     In the following modes for nproc>1, sudmft will write to processor-dependent file and
C     clean up afterwards.  Master node should append _0
      jfgk = 0                  ! jfgk > 0 => special management of gpmode output files (look for jfgk below)
      if (wgpmode > 0 .and. (gpmode>5)) then     ! Special file management for MPI case
        ifgk = fopna(fn,-1,lascii); rewind ifgk  ! File logical unit for sigk, gi, gkloc, or hkloc
        jfgk = ifgk             ! jfgk and ifgk are identical except for the master node in the MPI case

        if (gpmode == 3 .and. fullg .and. procid == master) then

           write(jfgk,'("#",5(x,i0),a,/,a)') nkfbz, 1, nomax_gpr, s_dmft%nzsig, 1
     &          ,"  # nkpt, nsymop, nom, cixdms norbitals"
     &          ,"#  0.00000  0.00000  0.00000    # actual position of correlated atoms in the unit cell" ! not really actual...
        else if (gpmode == 5 .and. procid == master) then
           !nzsig replaced by ndsig
           write(jfgk,'("#",5(x,i0),a,/,a)') nkfbz, 1, nomax_gpr, s_dmft%ndsig, 1
     &          ,"  # nkpt, nsymop, nom, cixdms norbitals"
     &          ,"#  0.00000  0.00000  0.00000    # actual position of correlated atoms in the unit cell" ! not really actual...
        end if

        if (nproc > 1) then
          if (procid == 0 .and. nproc > 1) then ! Rename master node file; preserve ifgk for cumulative results
            k = fopnx(fn,192,0,ifgk)          ! Full name of true output file.  fn should be preserved for merging
            ifgk = fopng(trim(fn)//'_0',-1,lascii) ! Append _0
            rewind ifgk
          endif
        endif
      endif

C ... Modes 18-20: Read qp levels from disk in full BZ; cull irr qp.  eigiq is relative to chemical potential
      if (gpmode >= 18 .and. gpmode <= 20) then
        ifi = fopna('proj',-1,4); rewind ifi

        allocate(eigiq(nlohi(1):nlohi(2),nkp,nspx)) ! qp levels, written to file se
        allocate(sigwq(nomax_gpr,nevn,nkp,nspx))    ! Diagonal sigma or Gk, written to file se
        call dpzero(sigwq,2*size(sigwq))

C   ... Read qp levels from disk in full BZ; cull irr qp. eigiq is relative to chemical potential
        do  iqfbz = 1, nkfbz
          call iodmftu(s_dmft,.true.,ndham,nlohi(1),nlohi(2),ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifi)
C         print 777, '1563 dmftu(nlohi(1)-1+3,1:3)',procid,iq,iqfbz,dmftu(nlohi(1)-1+3,1:3,1,1)
C          call zprm('dmftu',2,dmftu,nlohi(2)-nlohi(1),30,ncix)
          if (ig /= 1 .and. lstar) cycle
          eigiq(nlohi(1):nlohi(2),iq,1:nspx) = evlk(nlohi(1):nlohi(2),1:nspx) - chempot
        enddo
      else
C       Dummy array so compiler doesn't complain
        allocate(sigwq(1,1,1,1))
      endif

C     For the remainder of this routine, use nomgx as proxy for nomg.
C     The following cases apply nomgx may differ (elta is not made)
C     * Any case for which gpmode>1   (only one kind for now)

      nomgx = nomg
      if (iabs(gpmode) > 1) then
        nomgx = nomax_gpr
      endif
C     Update distribution of omegas
      kpproc(0) = 1; kpproc(1) = nomgx+1 ! Serial case
      if (procid == master .and. nproc > 1) then
        allocate(N(nomgx,1)); call dvset(N,1,nomgx,1d0)
        kpproc(2:) = 0
        call dstrbpx(nomgx,nproc,N,kpproc)
        deallocate(N)
      endif
      call mpibc1(kpproc,nproc+1,2,.false.,'','')
      if(cmdopt('--chiloc',8,0,strn)) then
         if (procid == master) then
            ef = ef0-efshft
            ifi = fopna('params_chi',-1,0)
            nomgv = s_dmft%nomf
            nomgw = s_dmft%nomb
            write(*,*)'chi0 calculated for ',nomgv,'fermionic frequencies'
            write(*,*)'chi0 calculated for ',nomgw,'bosonic frequencies'

            if(cmdopt('--qlist',7,0,strn)) then
               write(*,*) 'qlist will be read'
               if ( (nk1 /= nk2) .or. (nk1 /=nk3) ) then
                  call rx('qlist : not ready for non cubic k-mesh')
               endif

               call readqlist(0,nq,nk1,i123)
               allocate(i123(3,nq))
               call readqlist(1,nq,nk1,i123)

               call makechi0loc(1,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,
     .              nspx,nkfbz,ldcix,nsp,ncix,sigdc,ef,nomgv,nomgw,i123,nq,chi0,chi0loc)
            else
               call makechi0loc(0,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,
     .              nspx,nkfbz,ldcix,nsp,ncix,sigdc,ef,nomgv,nomgw,i123,0,chi0,chi0loc)
            endif
         endif
         call rx0('chi0 finished ')
      endif

C --- Make sigma, G, delta for all qp ---
C     gpmode=18-20 => restricted to the irreducible set
C     Band mode (nqband>0) => all points are included: no star available
      ifi     = fopna('proj',-1,4); rewind ifi
      allocate(gloc(ldcix,ldcix,nsp,ncix,nomgx))
      allocate(gkloc(ldcix,ldcix,nsp,ncix,nkfbz,nomgx))
      call dpzero(gloc,2*size(gloc))
      call dpzero(gkloc,2*size(gloc))

      if ( ((gpmode == 3  .and. fullg) .or.  gpmode == 5) .and. procid == 0) then
        if (any(s_bz % nkabc /= s_bz % nkabc(1))) call rx("sudmft: Only cubic BZ mesh accepted by suscept.py")
        klist_u = fopng('suscept.klist',-1,8)
        qlist_u = fopng('suscept.qlist',-1,8)
      end if
      allocate(i123(nkfbz,3))
      do  iqfbz = 1, nkfbz

         if (iqfbz >= 1 .and. nproc > 1) then
            i = fopnig(fileid,ifgk,lascii) ! Reopen ifgk
         endif

C     ... Read projector from disk
         call iodmftu(s_dmft,.true.,ndham,nlohi(1),nlohi(2),ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifi)
         qpr = qp

C     ... Irreducible qp only for these modes
         if (ig /= 1 .and. gpmode >= 18 .and. gpmode <= 20) cycle

C     ... Write k-point info to gkloc
         if (iabs(gpmode) > 1 .and. gpmode < 18) then
            if (ig == 1) iqs = 0 ! First point of new star
            if (lstar) call iqstar(2,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,iv)
C     Write full G in binary format.
C     First record : dimensioning parameters ndham,nomax_gpr,nkfbz,s_bz%nkabc,s_dmft%lsigim
            do i =1,3
               i123(iqfbz,i)=iv(i)
            enddo
C     Second record : omega(1:nomax_gpr)
C     Then follow pairs of records for each particular omega, k
C     First record : iomg,nlohi(1),nlohi(2),iqfbz,iq,iv,qpr (k-point information)
C     Second record : G(nlohi(1):nlohi(2),nlohi(1):nlohi(2))
            if (gpmode == 8) then
               if (iqfbz == 1) then
                  write(ifgk) ndham,nomax_gpr,nkfbz,s_bz%nkabc,isw(s_dmft%lsigim)
                  write(ifgk) omega(1:nomax_gpr)
               endif
            else
               if (procid == master) then
                  if ( ((gpmode == 3  .and. fullg) .or.  gpmode == 5) .and. procid == 0) then
                     write(klist_u, '(5(x,i0),x,f3.1)') iqfbz, iv(1:3)-1, s_bz%nkabc(1), 1.0
                     write(qlist_u, '(5(x,i0),x,f3.1)') iqfbz, iv(1:3)-1, s_bz%nkabc(1), 1.0
                  endif
               endif
            endif
         endif

C     --- Local Green's function and Eimp with updated chemical potential ---
C     call dpzero(gmk,2*size(gmk))
         allocate(sigbark(nevn,nevn,nspx))
!     do  iomg = 1, nomgx
!     if (iomg<kpproc(procid) .or. iomg>=kpproc(procid+1)) cycle
         do iomg = kpproc(procid), kpproc(procid+1)-1
            if (gpmode >= 8 .and. iomg > nomax_gpr) cycle

            if (gpmode == 8 .or.  gpmode== 9) then ! Write header for this omega, qp
               write(ifgk) iomg,nlohi(1),nlohi(2),iqfbz,iq,iv,qpr
            endif

c     ... Embed sigmabar in lattice representation
            call dpzero(sigbark,2*size(sigbark))
            call embed_sigma(s_dmft,nlohi(1),nlohi(2),ldcix,nsp,nspx,dmftu,sigdc,siginp(iomg,:),sigbark)


            if (gpmode == 9) then
               call rx('sudmft not ready for gpmode=9')
            elseif (gpmode == 19) then
               forall (i = 1:nevn) sigwq(iomg,i,iq,1:nspx) = sigbark(i,i,1:nspx)
               cycle
            endif

C     ... Add to local Green's function
            call makegloc2(s_dmft,ncix,ldcix,ndhamx,nlohi(1),nlohi(2),iq,nkp,nsp,nspx,iomg,
     .           ef0-efshft,s_dmft%gammac,evlk,omfac*omega(iomg),dmftu,sigbark,
     .           gpmode,nomax_gpr,gloc(1,1,1,1,iomg),sigwq,ifgk,gkloc(1,1,1,1,iqfbz,iomg))

         enddo                  ! frequency loop

C     --- Optional test. Hard-coded lchk_sigbar=3
         if (iprint() >= PRTG) then
            call rx('rewrite chk_sigbar')
C     call chk_sigbar(s_dmft,s_dmft%ndsig,nicix,orbmax,ndham,nevn,nsp,ndsigi(1),nkfbz,
C     .      dmftu,sigbark,siginp(nomgx,:,:),sigdc,chk_projsbark)
C     do  i = 1, nicix
C     chk_projsbar(1:ndsigi(i),1:ndsigi(i),:,:) =
C     .               chk_projsbar(1:ndsigi(i),1:ndsigi(i),:,:) +
C     .               chk_projsbark(1:ndsigi(i),1:ndsigi(i),:,:)/nkfbz
C     enddo
         endif

C     --- Accumulate Eimp as defined in eq 5.14 of Sponza's notes ---
         if (iabs(gpmode) == 0) then
C     Communicate to sigbark(nomgx) to head node
            call hunti(kpproc,nproc,nomgx,[-1],rankprocid)
            rankprocid = rankprocid-1
            call mpibc5(sigbark,size(sigbark),6,0,rankprocid,0)
C     print 777, '1668 before eeimp2',procid,rankprocid,0,sigbark(1,1,1),sum(sigbark)
            if (procid == master) then
               call makeeimp2(s_dmft,s_dmft%ndsig,ncix,ldcix,ndhamx,nlohi(1),nlohi(2),nsp,nspx,
     .              ef0-efshft,1,evlk,dmftu,sigbark,siginp(nomgx,:),sigdc,eimpm)
            endif
         endif
         deallocate(sigbark)
!     #Mark last term is a part of eq(3.10) but without DC subtraction (so it's just sum_k of projection of eps-mu)
!     #Lorenzo Correction: Eimp is now redefined as in eq5.14 of my notes: P(eval + E(Simp-DC)) - Simp-mu

      enddo                     ! k-points in the full BZ


      allocate(displs(0:nproc))
      allocate(recvcounts(0:nproc-1))

      displs = (kpproc-1)*ldcix*ldcix*nsp*ncix*nkfbz*2
      recvcounts = displs(1:nproc) - displs(0:nproc-1)


      if (procid == master) then
        call mpi_gatherv(mpi_in_place, 0, mpi_real8, gkloc, recvcounts, displs, mpi_real8, master, mpi_comm_world, err)
      else
        call mpi_gatherv(gkloc(1,1,1,1,1,kpproc(procid)), recvcounts(procid), mpi_real8, gkloc, recvcounts,
     &       displs, mpi_real8, master, mpi_comm_world, err)
      end if
c$$$      if (procid == master) then
c$$$         block
c$$$         integer :: nomgv,nomgw
c$$$         nomgv=50
c$$$         nomgw=1
c$$$         call  makechikh(s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,
c$$$     .        gkloc,ldcix,nsp,ncix,nkfbz,nomg,i123,nomgw,nomgv)
c$$$         end block
c$$$      endif
      deallocate(recvcounts)
      deallocate(displs)


C ... Combine sigma contributed by various processors
      if (nproc > 1 .and. gpmode >= 18 .and. gpmode <= 20) then
        call info0(20,0,-1,' ... Sharing sigma between processes...')
        xv(1) = mpiquery(0)  ! Time
        call mpibc2(sigwq,2*size(sigwq),4,3,.false.,'','')
        xv(2) = mpiquery(0)  ! Time
        call info2(20,0,0,' MPI broadcast took %;1d sec',xv(2)-xv(1),0)
      endif
C ... Combine gloc contributed by various processors
      if (gpmode < 18) then
        call mpibc2(gloc,2*ldcix*ldcix*nsp*ncix*nomgx,4,3,.false.,'','')
        call dscal(2*size(gloc),1d0/nkfbz,gloc,1)

      endif
C....... write gkloc  and gloc
      if ((procid == master) .and. (gpmode < 6)) then
         call info0(20,0,-1,' ... Writing gkloc and gloc...')
         if(fullg) gpmode=5
         if(gpmode==3) then
            call  printgprt(s_dmft,s_bz,s_ham,s_lat,1,ldcix,nkp,nspx,nsp,ncix,nkfbz,nomgx,gkloc,gloc)
         else if (gpmode==5) then
            if(cmdopt('--reduc',7,0,strn))then
               call  printgprt(s_dmft,s_bz,s_ham,s_lat,3,ldcix,nkp,nspx,nsp,ncix,nkfbz,nomgx,gkloc,gloc)
            else
               call  printgprt(s_dmft,s_bz,s_ham,s_lat,2,ldcix,nkp,nspx,nsp,ncix,nkfbz,nomgx,gkloc,gloc)
            endif
         endif
         call info0(20,0,-1,' ... gkloc finished...')
      endif                     !master
      if (procid == master .and. (fullg .or. (gpmode==5))) then
        write(klist_u, '("END")')
        write(qlist_u, '("END")')
        write(klist_u, '(3(3(x,f11.6),/))') 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0
        close(qlist_u)
        close(klist_u)
      end if



C ... Write qp, eigenvalues
      if (iabs(gpmode) > 1) then
         if (procid == master) then
            if (gpmode >= 18 .and. gpmode <= 20) then
C     call zprm('sigwq',2,sigwq,nomgx,nomgx,nevn*nkp)
               call iose2(iseh,0,0,nspx,nlohi(2)-nlohi(1)+1,nkp,nomax_gpr,s_bz%qp,eigiq,xv,sigwq,xv,xv,-ifgk)
C#ifdef H5
               open(file='dos.h5',unit=123)
               close(123,status='delete')

               call zdscal(size(sigwq),1/eVbyRy,sigwq,1)
               !recall : sigwq is g ....
               call h5_write('dos.h5:/g', sigwq, find_h5dtype(sigwq))
               call h5_write('dos.h5:/omega', eVbyRy*omega, find_h5dtype(omega))
C#endif
            endif
            call fclose(ifgk)
         endif

        if (gpmode < 8) call rx0(': '//trim(strn)//' generated file gkloc')
        if (gpmode == 8) call rx0(': '//trim(strn)//' generated file gk')
        if (gpmode == 18) call rx0(': '//trim(strn)//' wrote diagonal g to file se')
        if (gpmode == 19) call rx0(': '//trim(strn)//' wrote diagonal sigma to file se')
        if (gpmode == 20) call rx0(': '//trim(strn)//' wrote effective diagonal sigma to file se')
        call rx0(strn)
      endif

      if (procid == master) then
        call dscal(2*size(eimpm),1d0/nkfbz,eimpm,1)
C       call zprm('eimpm(1)',2,eimpm,ldcix,ldcix,ldcix)
C       call zprm('eimpm(nsp)',2,eimpm(1,1,nsp,1),ldcix,ldcix,ldcix)
      endif

c --- Optional test of sigbar at infinity
      if (iprint() >= PRTG) then
        call rx('rewrite chk_sigbar_prt')
C        call chk_sigbar_prt(s_dmft,orbmax,nicix,ndsigi(1),s_dmft%ndsig,
C     .    chk_projsbar,siginp(nomgx,:,:),sigdc,lchk_sigbar)
C        deallocate(chk_projsbark,chk_projsbar)
      endif

C --- Construct hybridization function ---
      if (procid == master) then
         call dpzero(dmftdelta,2*size(dmftdelta))
         call dpzero(dmftdeltah5,2*size(dmftdeltah5))
         call dpzero(g0h5,2*size(g0h5))
        call makedelta3(s_dmft,ldcix,nsp,ncix,nomgx,gloc,eimpm,omega,omfac,dmftdelta,dmftdeltah5,g0h5,sig_ldadc)
!       #Mark in last routine the DC is subtracked from Eimp to obtain eq.(3.10) and
!       sigma_imp subtracted from Delta, to obtain final hybridization
!       #Lorenzo this has been changed in makedelta2: Nothing is subtracted to Eimp.
      endif


C ... Write delta and eimp1 files, in eV
      if (procid == master) then
        call info0(30,0,0,' Writing files delta and eimp1 ...')
        ifi = fopna('delta',-1,2); rewind ifi
        i = iodsig(201,s_dmft%nomg,s_dmft%beta,s_dmft%ndsig,s_dmft%omg,dmftdelta,-ifi)
        call fclose(ifi)

        ifi = fopna('eimp1',-1,2); rewind ifi
        if (lsig0)  then
           i = 0; call bin2av('D11',0,0,6,eVbyRy*dble(sig_ldadc),4,0,size(sig_ldadc)-1,
     .    ',',len(fmt),.false.,fmt,i)
        write(ifi,543) 'Edc=[',fmt(1:i),']  # Double counting'


c$$$          call daxpy(s_dmft%ndsig,-1d0,sig_ldadc,1,ausp,2)
c$$$          sig_ldadc = 0d0
        else
          i = 0; call bin2av('D11',0,0,6,eVbyRy*dble(sigdc),4,0,size(sigdc)-1,
     .      ',',len(fmt),.false.,fmt,i)
          write(ifi,543) 'Edc=[',fmt(1:i),']  # Double counting'
        endif

  543   format(a,a,a)
        allocate(ausp(s_dmft%ndsig))  ! ausp is just a work array here
        allocate(nicixi(ncix))        ! Loop over inequivalent cix only
        call ineqcix(s_dmft,ncix,nicixi)
        do  cix = 1, ncix
          if (nicixi(cix) >= 0) cycle ! skip equivalent cix
          call sigp2sigij(s_dmft,102,1,nsp,1,cix,ldcix,ausp,ausp,eimpm(:,:,:,cix))
        enddo
        i = 0; call bin2av('D11',0,0,6,eVbyRy*dble(ausp+sigdc),4,0,size(ausp)-1,
     .    ',',len(fmt),.false.,fmt,i)
        write(ifi,543)
     .       'Eimp=[',fmt(1:i),'] # Eimp: P.(e+E.sbar)-sinp-mu+Edc'

        i = 0; call bin2av('D11',0,0,6,eVbyRy*dble(ausp-ausp(1)),4,0,size(ausp)-1,
     .       ',',len(fmt),.false.,fmt,i)
        write(ifi,543) 'Eimp=[',fmt(1:i),'] # Eimp: shifted by Eimp[0]'

        open(newunit=u, file='solver_input.h5')
        close(u,status='delete')

        call zdscal(size(dmftdeltah5), evbyry, dmftdeltah5, 1)
        call zdscal(size(g0h5), 1./evbyry, g0h5, 1)

      open(file='solver_input.h5',unit=123)
      close(123,status='delete')
      call h5_write('solver_input.h5:/eimp_wo_dc', eVbyRy*dble(ausp+sigdc),h5t_native_real8)
      call h5_write('solver_input.h5:/delta', dmftdeltah5, h5t_native_cmplx8)
      call h5_write('solver_input.h5:/g0', g0h5, h5t_native_cmplx8)
       call h5_write('solver_input.h5:/nsp', nsp, find_h5dtype(nsp))
       call h5_write('solver_input.h5:/ncix', ncix, find_h5dtype(ncix))
       call h5_write('solver_input.h5:/nicix', nicix, find_h5dtype(ncix))
       call h5_write('solver_input.h5:/nomg', s_dmft%nomg, find_h5dtype(nomg))
       call h5_write('solver_input.h5:/beta', s_dmft%beta, find_h5dtype(s_dmft%beta))
       call h5_write('solver_input.h5:/icix', s_dmft%icix, find_h5dtype(s_dmft%icix(1)))
       call h5_write('solver_input.h5:/l', s_dmft%l, find_h5dtype(s_dmft%l(1)))
       call h5_write('solver_input.h5:/sigind', s_dmft%sigind, find_h5dtype(s_dmft%sigind(1,1,1,1)))

C       if Sig=0, we need to shift Ed with double counting because it has not been done above.
        if (lsig0) then
          call daxpy(s_dmft%ndsig,-1d0,sig_ldadc,1,ausp,2)
          sig_ldadc = 0d0
        endif

        i = 0; call bin2av('D11',0,0,6,eVbyRy*dble(ausp),4,0,size(ausp)-1,',',len(fmt),.false.,fmt,i)
        write(ifi,543) 'Ed [',fmt(1:i),'] # Ed: Eimp-Edc for PARAMS'
        call h5_write('solver_input.h5:/eimp', eVbyRy*dble(ausp), find_h5dtype(dble(ausp(1))))

        i = 0; call bin2av('D11',0,0,6,-eVbyRy*dble(ausp(1)),4,0,0,' ',len(fmt),.false.,fmt,i)
        write(ifi,543) 'mu ',fmt(1:i),' # mu = -Eimp[0] for PARAMS'
        call fclose(ifi)

      endif

      call tcx('delta')
C ... Check self-consistent condition Gloc=Gimp (and optional print with  '--gprt')
!      call chksccondition(ldcix,nsp,ncix,nomgx,gloc,omega,cmdopt('--gprt',6,0,strn))
      call chk_selfcons(ldcix,nsp,ncix,nomgx,gloc)

      deallocate(omega)
      call tcx('sudmft')
      call rx0('done making DMFT hybridization function')

  777 format(a,3i5,10f14.7)
      end subroutine sudmft

      subroutine wlattc(mode,ifi,jfi,s_site,s_spec,alat,plat,nbas,nat,
     .  maxcor,GcutH,ipb,ipcx,lmxax,bas,lmxa,nsp,konf)
C- Write LATTC (or lattice) file
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 format for fpgw v033a5
Ci         :1 format for spex v02.05
Ci         :2 format for DMFT v1
Ci   ifi   :file handle for lattice data
Ci   jfi   :file handle for lmto data
Ci  s_site :struct for site-specific data; see structures.h
Ci    Elts read: spec pos pnu pz ov0 clabel
Ci    Stored:    *
Ci    Passed to: *
Ci  s_spec :struct for species-specific data; see structures.h
Ci    Elts read: lmxa pz p idxdn a nr z rmt
Ci    Stored:    *
Ci    Passed to: *
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci         :Written to disk (modes 0,1)
Ci   plat  :primitive lattice vectors, in units of alat
Ci         :Written to disk (modes 0,1)
Ci   nbas  :size of basis, including floating orbitals sites
Ci   nat   :size of true basis (exclude floating orbitals sites)
Ci         :Written to disk (modes 0,1)
Ci   maxcor:(1) = largest number of core radial functions of any site in system
Ci         :(2) = largest l for which a core state exists
Ci   GcutH :G cutoff for eigenfunctions
Ci         :Written to disk (mode 1)
Ci   ipb   :index to reduced basis (excluding floating orbitals)
Ci         :given site index including those orbitals
Ci         :Written to disk (modes 0,1)
Ci   z     :Nuclear charge
Ci         :Written to disk (mode 1)
Ci   ipcx  :class index: site ib belongs to class ipcx(ib)
Ci         :Written to disk (mode 1)
Ci   lmxax :global maximum of lmxa
Ci         :Written to disk (modes 0,1)
Ci   bas   :site positions, Cartesian coordinates, units of alat
Ci         :Written to disk (mode 1)
Ci   lmxa  :augmentation l-cutoff for each site
Ci         :Written to disk (modes 0,1)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci         :Written to disk (mode 1)
Ci   konf  :principal quantum numbers
Ci         :Written to disk (modes 0,1)
Co Outputs
Co   data is written to file logical unit ifi
Cl Local variables
Cl  maxradv:Maximumum number of radial functions describing partial waves
Cl         :for a given l (2 for phi,phidot; 3 if any local orbitals)
Cl         :Written to disk (mode 1)
Cr Remarks
Cr   fpgw mode writes to a single file (mixeed data)
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   07 Jul 09 New mode; altered argument list
Cu   05 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,ifi,jfi,nbas,nat,nsp,lmxax,maxcor(2)
      integer ipb(nbas),ipcx(nbas),lmxa(nbas),konf(0:lmxax,nbas,nsp)
      double precision plat(3,3),alat,GcutH,bas(3,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer n0,nrmx,maxrad
      parameter (n0=10,nrmx=5001)
      double precision qlat(3,3),xx,posp(3),tol,roundp,z,a,rmt,pz(n0,2)
      integer i,l,ib,is,ic,isp,nclasx,mxint,nlod,iclbas,nr,nrl,nglob,llo(n0),iplo(n0)
      character slabl*2,spid*8,fmt*80,fmt0*80,strn*100
C     integer,allocatable:: lmxa2(:),ipc2(:)

      tol = 1d-6
      nclasx = mxint(nbas,ipcx)

C      allocate(lmxa2(nat),ipc2(nat))
C      ia = 0
C      do  ib = 1, nbas
C        if (lmxa(ib) > -1) then
C          ia = ia+1
C          lmxa2(ia) = lmxa(ib)
C          ipc2(ia) = ipcx(ib)
C        endif
C      enddo

C --- spex format ---
      if (mode == 1 .or. mode == 3) then
        if (nglob('lpz') > 0) then
           maxrad = 3
         else
           maxrad = 2
         endif
        nlod = nglob('lpz')/1000
        call awrit8(' nat=%i nclass=%i nsp=%i alat=%;8d maxloc=%i maxradv=%i maxradc=%i maxcorl=%i',
     .     ' ',80,ifi,nat,nclasx,nsp,alat,nlod,maxrad,maxcor(1),maxcor(2))
        call awrit8('%x nat=%i nsp=%i nclass=%i lmxax=%i maxloc=%i '//
     .    'maxradv=%i maxradc=%i maxcorl=%i',strn,len(strn),0,
     .    nat,nsp,nclasx,lmxax,nlod,maxrad,maxcor(1),maxcor(2))
        call awrit1('%a gcuth=%;3d',strn,len(strn),-jfi,GcutH)

C   ... Lattice vectors are extended to 12 digits
        write(ifi,"(' plat:')")
        do  ic = 1, 3
          do  ib = 1, 3
            posp(ib) = roundp(plat(ib,ic),tol)
C           posp(ib) = plat(ib,ic)
          enddo
          write(ifi,"(3f18.12,'  (P',i1,')')") posp,ic
        enddo

C   ... Class data
C        write(ifi,"(' Class data:'/'  ic    Z    Symbol')")
        write(ifi,345)
        write(jfi,345)
  345   format(' Class data:')
        call awrit2('%x  ic    Z        rmt        a       nr   lmxa       konf'//
     .    '         %?#n#llo  %nfilo%-1j%nf#',strn,len(strn),jfi,nlod,3*nlod-3)
        call awrit2('%x  ic    Z    Symbol   rmt        a       nr   lmxa   %nfkonf'//
     .    '   %?#n#%-2j%nfllo  ilo##',strn,len(strn),ifi,3*lmxax/2,nlod)

C       fortran format for writing class line
        fmt0 = '(i4,f8.2,f12.7,f10.6,2i5,3x,a1,'
C       call awrit1('%a%?#n#%-1j%ii3,2x,%-1j%ii3,3x,##99i3)',fmt0,len(fmt0),0,nlod)
        call awrit2("%a%ii3,a1,1x,%?#n#%-1j%ii3,2x,%-1j%ii3,3x,##)",fmt0,len(fmt0),0,lmxax+1,nlod)
        fmt = '(i4,f8.2,3x,a,f12.7,f10.6,2i5,3x,a1,'
        call awrit2("%a%ii3,1x,a1,1x,%?#n#%-1j%ii3,2x,%-1j%ii3,3x##)",fmt,len(fmt0),0,lmxax+1,nlod)

        do  ic = 1, nclasx
          ib = iclbas(ic,ipcx,size(ipcx))
          if (ib == 0) cycle
          is = s_site(ib)%spec
          z = s_spec(is)%z
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          nrl = nr
          spid = s_site(ib)%clabel
          if (mode == 1) call mshnrl(1,rmt,a,nr,nrmx,nrl)
          pz = s_spec(is)%pz
          i = 0
          llo = -1
          iplo = 0
          do  l = 0, lmxa(ib)
            if (pz(l+1,1) /= 0) then
              i = i+1
              llo(i) = l
              iplo(i) = int(pz(l+1,1))
            endif
          enddo
          konf(lmxa(ib)+1:lmxax,ib,1) = 0
          if (i > nlod) call rx('wlattc: number of local orbitals > nlod')
C         print *, is, sngl(pz(1:3)),llo(1:4)
C         z = s_spec(is)%z
          call zslabl(-1,slabl,iabs(nint(z)))
C         write(jfi,"(i4,f8.2,f12.7,f10.6,i5,i5,3x,99i4)")
          write(jfi,fmt0) ic,z,rmt,a,nrl,lmxa(ib),'[',konf(0:lmxax,ib,1),']',
     .      (llo(i),i=1,nlod),(iplo(i),i=1,nlod)
C         write(ifi,"(i4,f8.2,3x,a)") ic,z,slabl
          write(ifi,fmt) ic,z,slabl,rmt,a,nrl,lmxa(ib),'[',konf(0:lmxax,ib,1),']',
     .      (llo(i),i=1,nlod),(iplo(i),i=1,nlod)
        enddo

C   ... Site data.  Positions in unit of plat ... expanded to 12 digits
C       qlat = (plat+)^-1
        call dinv33(plat,1,qlat,xx)
        write(ifi,"(' Site data:'/
     .    '  ib  ic',15x,'pos (units of plat)',22x,'spid')")
        do  ib = 1, nbas
C         call spacks(0,'site clabel',ssite,spid,ib,ib)
          spid = s_site(ib)%clabel
          if (lmxa(ib) > -1) then
C           posp+ = (plat)^-1 pos+
            call dgemm('T','N',3,1,3,1d0,qlat,3,bas(1,ib),3,0d0,posp,3)
            do  ic = 1, 3
              posp(ic) = roundp(posp(ic),tol)
            enddo
            write(ifi,"(2i4,3f18.12,3x,a)") ipb(ib),ipcx(ib),posp,trim(spid)
          endif
        enddo

        return

C --- fpgw format --
      else

        write(ifi,"(e24.16)") alat
        write(ifi,"(3e24.16)") plat(1:3,1)
        write(ifi,"(3e24.16)") plat(1:3,2)
        write(ifi,"(3e24.16)") plat(1:3,3)
        write(ifi,*) ' -1d10 ! dummy entry. True QpGcut_psi is in GWinput'
        write(ifi,*) ' ------------------------------------------- '
        write(ifi,"(2i4,' ! nbas lmxax (max l for argumentaion)')") nat,lmxax
        write(ifi,*) ' ------------------------------------------- '
        do  isp = 1, nsp
          write(ifi,"(' -- ibas lmxa konf(s) konf(p) konf(d)... ',' isp=',2i2)") isp
          do  ib = 1, nbas
            if (lmxa(ib) > -1) then
              write(ifi,"('   ',99i4)") ipb(ib),lmxa(ib),konf(0:lmxa(ib),ib,isp)
            endif
          enddo
        enddo

      endif

      end

      subroutine pvsug1(nbas,lmxa,ipc,ipcx)
C- Remakes class table, expunging classes with floating orbitals.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   lmxa  :augmentation l-cutoff
Ci   ipc   :class table: site ib belongs to class ipc(ib)
Co Inputs/Outputs
Co   ipcx  :expunged class table: classes with lmxa=-1 are expunged
Co         :and the remaining classes are sequentially renumbered
Co         :preserving the order of the remaining classes.
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Mar 07 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,lmxa(nbas),ipc(nbas),ipcx(nbas)
C ... Local parameters
      integer i,ibp,ic,ipskip,nelim
      integer prm(nbas)

      call ivshel(1,nbas,ipc,prm,.true.)
C      do  i = 1, nbas
C        prm(i) = prm(i)+1
C      enddo

C     nelim = number of classes that have been eliminated so far
      nelim = 0
C     ipskip is the number of the last class that was skipped.
C     Multiple occurrences of a skipped class must still only
C     reduce the net number of classes by one.
C     We must avoid incrementing nelim when multiple sites
C     correspond to a skipped class.
      ipskip = 0
C     Loop over sites in order of increasing class index, ibp
      do  i = 1, nbas
        ibp = prm(i)+1
        ic = ipc(ibp)
C       Test whether this class should be purged
        if (lmxa(ibp) < 0) then
          if (ipskip /= ic) nelim = nelim+1
          ipskip = ic
          ipcx(ibp) = -1
        else
          ipcx(ibp) = ic - nelim
        endif
      enddo
      end
