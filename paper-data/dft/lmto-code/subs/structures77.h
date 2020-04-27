! Structures used globally.
! Detailed descriptions are given with individual declarations.
!  Structure Contents
!  str_bz    Parameters for Brillouin zone integration
!  str_ctrl  for control parameters
!  str_gw    for GW and dielectric response
!  str_ham   for hamiltonian parameters
!  str_lat   for lattice structure
!  str_mix   Parameters for mixing charge density and other quantities
!  str_move  Parameters for molecular and spin dynamics
!  str_optic Parameters for optics
!  str_ordn  Parameters for order-N
!  str_pot   Parameters for density, potential, matrix elements
!  str_spec  Parameters for species data
!  str_site  Parameters for site data
!  str_str   Parameters for structure constants, lmto, tb and TCF
!  str_tb    Parameters for Levenberg-Marquardt fitting of bands, DOS
!  str_strn  Holds global strings
!  str_symops Information about symmetry operations
!  str_gwse  structure containing dynamical GW self-energy

!u Updates
!u   11 Nov 13 (Belashchenko) New parameters for fully relativstic branch
!u             (MvS) additional parameter for noncollinear magnetism
!u   02 May 12 New str_symops
!u   25 Apr 12 (Belashchenko) New parameters for chemical CPA

! --- Notes on symmetry operations ---
!r    Space group operation is defined by (g,a)
!r    g is the 3x3 point group part; a is the translation part.
!r    If g is in cartesian coordinates and b is some vector,
!r    then (g,a) rotates b into b':
!r      b' = g b + a
!r    If B and B' are corresponding vectors in units of P
!r    b = P B; b' = P B'
!r    Then
!r      P B' = g P B + a
!r      B' = P^-1 g P B + P^-1 a = G B + A where
!r      G = P^-1 g P  and A = P^-1 a
!r    In a crystal, P are primitive lattice vectors => G is integer
!r    For further description, see subs/mks_sym.f
      type str_symops
      sequence

!       Point group and translation defining space group operation
        real(8) :: rg(3,3),ag(3)
!       Points to which group operation is inverse of this one
        integer :: inv
!       prg = P^-1 g P, called G in Remarks
!       pag = P^-1 a, called A in Remarks
        integer :: prg(3,3)
        real(8) :: pag(3)

        integer, pointer  :: mib(:)
        real(8), pointer  :: tib(:,:),rotylm(:,:)

      end type str_symops

! --- GW self-energy structure ---
      type str_gwse
      integer   ::  mode     ! Contents of sigma
                             ! 1s digit 0 => real axis  1 => Matsubara axis
      integer   ::  nband    ! Number of bands in SE file
      integer   ::  nsp      ! Number of spins in SE
      integer   ::  nomg     ! Number of frequency points
      integer   ::  nkp      ! Number of qp
      integer   ::  nblst    ! Number of bands for which sigma is retained.  See iblst.
      real(8)   ::  ommin,ommax ! range of frequencies (real axis only)
      real(8)   ::  beta        ! Temperature (Matsubara axis only)

      integer, pointer :: iblst(:,:) ! iblst(1,1:nblst) List of bands for which sigma is read
                                     ! iblst(2,:)       List of bands for which sigma is on disk
      real(8), pointer :: qp(:,:)    ! qp qp(3,nqp) at which sigma is given
      real(8), pointer :: eig(:,:,:) ! QP levels eig(nband,nqp,nsp) (QSGW sense)
      real(8), pointer :: sxq(:,:,:) ! sxq(nband,nqp,nsp) corresponding Fock exchange potential
      real(8), pointer :: xcq(:,:,:) ! xcq(nband,nqp,nsp) corresponding LDA xc potential
      real(8), pointer :: sgq(:,:,:) ! sgq = sigwq at eigenvalue
      complex(8), pointer :: sigwq(:,:,:,:) ! sigwq(nomg,nband,nqp,nsp)

      end type str_gwse

! ... General purpose vector of integer work arrays
      type ip_wk_vec
      integer, pointer :: p(:)
      end type ip_wk_vec

! ... General purpose vector of 2D integer work arrays
      type ip_wk2_vec
      integer, pointer :: p(:,:)
      end type ip_wk2_vec

! ... General purpose vector of double precision work arrays
      type dp_wk_vec
      real(8), pointer :: p(:)
      end type dp_wk_vec

! ... General purpose vector of double complex work arrays
      type dc_wk_vec
      complex(8), pointer :: p(:)
      end type dc_wk_vec

! ... For vectorized structure constants; see strxsu.f
      type str_strxsu
      integer, pointer :: ip(:),ikl(:)
      real(8), pointer :: cf(:)
      type(ip_wk_vec) jkl(0:20)
      end type str_strxsu

! ... For Pade interpolation on the Matsubara frequencies
      type str_mpade
      logical :: lpade  ! .true. if contents have been set
      logical :: lnull  ! Not used
      integer :: mode   ! not used now
      integer :: npade  ! number of points entering into Pade
      complex(8), pointer :: zomega(:)        ! Energies used to make fit
      complex(8), pointer :: padcof(:,:,:)    ! Pade coefficients
      end type str_mpade

! --- str_bz ---
!- Parameters governing Brillouin zone integration
! ----------------------------------------------------------------
!r  element  purpose
!r  def      Uncertainty in Fermi level (rdctrl)
!r  dosw     Energy window over which DOS accumulated (rdctrl,bndfp)
!r  ef       Fermi level (rdctrl,iors,bndasa,bndfp)
!r  efmax    eigenvalues above efmax are not needed (rdctrl)
!r  egap     energy bandgap (bzwtsf)
!r  fsmom    Impose global external magnetic field (rdctrl)
!r           1st arg: field fixed by specifing moment M (fixed spin moment method)
!r           2nd arg: constant field addded
!r  lcond    Used by lmdos to calculate conductivity in place of DOS (rdctrl)
!r           lcond(1):
!r           0: Calculate DOS D(E) instead of conductivity (default)
!r           1: ballistic conductivity <v.DIR> (must also supply DIR)
!r           2: <v_i . v_j> (diffusive conductivity)
!r  sopertp  parameters used to estimate Fermi surface contribution to change in
!r           band energy from perturbed eigenvalues. Estimate by sampling.
!r           sopertp(1) = sampling Fermi level
!r           sopertp(2) = sampling sumev
!r           sopertp(3) = sum [delta w] sum of change in weights
!r           sopertp(4:6) hold terms involving 1st order perturbation t1
!r           sopertp(4) = sum [delta (w e)] change in weight-eval product, sampling
!r           sopertp(5) = sum [w (delta e)] weighted change in evals by sampling
!r           sopertp(6) = sum [w (delta e)] weighted change in evals by tetrahedron
!r           sopertp(7:9) hold terms involving 2nd order perturbation t1-t2
!r           sopertp(7) = Same as sopertp(4), but for t1-t2
!r           sopertp(8) = Same as sopertp(5), but for t1-t2
!r           sopertp(9) = Same as sopertp(6), but for t1-t2
!r  lio      controls which of various quantities are I/O to disk (rdctrl)
!r           1  read qp
!r           2  write qp
!r           4  write density matrix in k space
!r           8  write (density matrix - 1/2) in k space
!r           16 save evcs
!r           32 save dos weights
!r           NB: 4,32 out of date ; see ldos
!r  lmet     controls how band codes q-point integration weights (rdctrl)
!r           0 assume insulator: weight given by number of times irr
!r             q-point is represented in full BZ.
!r           Other modes deal with metals and weights are modified for
!r           states near the Fermi level
!r           1 save condensed evecs on disk (ASA only)
!r           2 read qp weights from disk from prior iteration; otherwise
!r             make and initial band pass only to generate the weights
!r           3 Always make 2 band passes
!r           4 Fermi level by tetrahedron, density by sampling, 3-point interpolation
!r           5 keep evecs in memory; reuse on second pass (lmf only)
!r  lmull    Turn on Mulliken analysis; used by tbe (rdctrl)
!r  lopt     a special mode (rdctrl)
!r           0 default mode (nothing special)
!r           2 BZMAP mode
!r  lshft    shift mesh to straddle q=0.  One component for each axis (rdctrl)
!r  lwtkb    describes status of wktb array (subzi,bzwtsf)
!r           0 weights are neither required nor available a priori
!r           1 weights are required a priori, and are read from disk
!r          -1 weights are required a priori, but were not read
!r           2 weights generated with constrained global moment
!r           4 Special interacting hamiltonian mode with the following generalization
!r             1. eigenvectors may not be eigenstates of a 1-body hamiltonian
!r             2. eigenvalues are not given; properties depending on eigenvalues
!r                are not generated.
!r             3. wtkb is supplied: wtkb and eigenvectors are synchronized,
!r                but they are not necessarily ordered by energy
!r  lswtk    Flags whether to make 'spin weights' array swtk (subzi)
!r          -2 do not make spin weights
!r           1 spin weights array allocated; make them
!r  n        Polynomial order for Methfessel-Paxton sampling (rdctrl)
!r  ndos     No. DOS points for sampling integration, and lmdos  (rdctrl)
!r  nevmx    eigenvectors above nevmx are not needed (rdctrl)
!r  nkabc    number of divisions in qp mesh (rdctrl,mqp,evalqp,sugws)
!r  nkp      Number of qp (mkqp,rdsigm)
!r  ntet     number of tetrahedra (mkqp)
!r  numq     number of trial Fermi levels by sampling scheme for interpolation (subzi)
!r  nef      number of Fermi levels, used e.g. for L.S. perturbation theory
!r  dos      pointer array for total DOS
!r  idtet    pointer array marking tetrahedron corners (mkqp)
!r  ipq      pointer array:maps point in full BZ to element in the irreducible set (mkqp,evalqp)
!r  pdos     pointer array for partial DOS
!r  qp       pointer array for qp (mkqp,evalqp)
!r  star     pointer array containing info on star of k (mkqp)
!r  wtkp     pointer array containing symop weighting of qp (mkqp,evalqp)
!r  wtkb     pointer array containing qp weights, wtkp*occupation  (subzi)
!r  range    number of FWHM for sampling integration (rdctrl)
!r  semsh    nz, modec, emin,emax, ecc,eps delta (rdctrl)
!r           See gf/emesh.f for detailed explanation of its contents
!r  swtk     'spin weights': diagonal part of  (z)^-1 sigmz z
!r           where z are eigenvectors, sigma is the Pauli spin matrix.
!r           Supplies information about spin moment in noncoll case.
!r  w        Line broadening for sampling integration (rdctrl,iors)
!r  zinsul   tolerance in deviation from c.n., used by Green's function technique for insulators
!r           Fermi level is average of levels for for (c.n. point +/- zinsul)
!r  zval     No. electrons to accumulate in BZ integration (rdctrl,getzv)
! ----------------------------------------------------------------
      type str_bz
!     sequence

      integer    ::   lio
      integer    ::   lmet
      integer    ::   lmull
      integer    ::   lopt
      integer    ::   lswtk
      integer    ::   lwtkb
      integer    ::   lshft(3)
      integer    ::   n
      integer    ::   ndos
      integer    ::   nevmx
      integer    ::   nkabc(3)
      integer    ::   nkp
      integer    ::   ntet
      integer    ::   numq
      integer    ::   nef
      integer    ::   inull ! add for 8-byte word boundary

      real(8)    ::   def
      real(8)    ::   dosw(2)
      real(8)    ::   ef
      real(8)    ::   efmax
      real(8)    ::   egap
      real(8)    ::   fsmom(2)
      real(8)    ::   lcond(4)
      real(8)    ::   sopertp(12)
      real(8)    ::   range
      real(8)    ::   semsh(10)
      real(8)    ::   w
      real(8)    ::   zinsul
      real(8)    ::   zval

! ... Integer pointer arrays
      integer, pointer::    idtet(:,:)
      integer, pointer::    ipq(:)
      integer, pointer::    star(:)

! ... Real pointer arrays
      real(8), pointer::    dos(:)
      real(8), pointer::    pdos(:)
      real(8), pointer::    qp(:,:)
      real(8), pointer::    wtkp(:)
      real(8), pointer::    wtkb(:)
      real(8), pointer::    swtk(:)

      end type str_bz

! --- str_ctrl ---
!- Parameters governing program flow and what codes calculate
! ----------------------------------------------------------------
!r  element purpose
!r  cllst   (O-N) list of sites belonging to each cluster (susite)
!r  clp     (O-N) cluster-related dimensioning parameters (suordn)
!r  clssl   (O-N) site ib belongs to cluster clssl(ib)
!r  dclabl  class label array (r8tos8 format) (mksym.f)
!r  defm    Lines for deformation when relaxing plat (rdctrl)
!r  elin    ccor linearization energy, 2C hamiltonian (rdctrl)
!r  group   (ASA) force equivalence of classes, even if not sym-related (clsprm)
!r  initc   Switches containing what info available, by class
!r          1 P,Q  2 pp  4 sop  8 vintra  16 pmpol  32 gradrm  64 pot
!r  ics     Index to which species each class belongs (mksym)
!r  idcc    (CPA) points to links between CPA and parent classes (mksym,dlmidx)
!r          idcc(ic=1:nclass) points to child classes:
!r          If class ic is NOT a CPA class, idcc(ic) points to itself
!r          If ic IS a CPA class, idcc(ic) points to the first of
!r          of the CPA or DLM classes associated with ic.
!r          idcc(i=nclass+1:nclass+nccomp) (CPA classes):
!r          index to the parent class for this CPA class
!r  ipc     index to which class each site belongs (mksym)
!r  ipcp    similar to ipc, but differs for pgf padded basis atoms(mksym)
!r  ips     index to which species each site belongs (susite)
!r  lasa    parameters governing execution of ASA code (rdctrl)
!r          1 Make V from P,Q  2 Make pp from given V  4 ccor  8 free atm
!r          16 map  32 nonspherical mpol moms 64 MT corr
!r          128 interpretation of sphere Q2; see newrho.f
!r          256 how ASA Q1,Q2 are accumulated; see makwts.f
!r          512 (spin pol) beta repsn = (gamma(1) + gamma(nsp))/2
!r          1024 (spin pol) rotate to beta repsn by scaling S, not P
!r          2048 not set if s_site(:)%pfr, s_site(:)%dpfr, s_site(:)%ddpfr stored as vector
!r               set     if s_site(:)%pfr, s_site(:)%dpfr, s_site(:)%ddpfr stored as matrix
!r          4096 Make phidot, phidotdot by outward integration only
!r  lbas    1 Hamiltonian has no screening transformation (rdctrl)
!r          2 Hamiltonian is nfp style
!r          16 freeze phi,phidot, and log derivative parameter P for all species
!r             In the ASA, freeze P only
!r  lbxc    switch to impose constraining fields through scaling
!r          of LDA XC potential.  Scaling is site-dependent and
!r          is read/written through file bxc (rdctrl).
!r          0 no modification of LDA XC potential
!r          1 Scale field as 1+bxc(ib), ib=site index
!r          2 Scale field as 1+bxc^2(ib), ib=site index
!r            bit lbxc can be 1 or 2 but not 3.
!r          4 impose constraining fields for vector DLM  (rdctrl)
!r            This automatically sets 1+2's bits lbxc to 2.
!r          8 (lmf) scale Bxc for smooth density as well as local density
!r         16 (lmf) spin-average potential so LDA part of Bxc=0
!r  lcd     1 freeze core  (rdctrl)
!r          2 non-self-consistent Harris
!r          4 represent full potential density on a uniform mesh
!r          8 represent full potential density via TCF
!          16 When calculating core, calculate w/ spin-averaged potential
!r         32 average l=0 input site rho in sites sharing a common GRP2
!r         64 average l=0 output site rho in sites sharing a common GRP2
!r         64 (molecules) XC potential by FFT
!r        128 force input site rho positive at all points
!r        256 force output site rho positive at all points
!r  lcgf    Green's function (rdctrl)
!r         >0:  Green's function calculation
!r          1:  Calculate G for irreducible k-points
!r          10: Transverse exchange interactions J(q) from MST, LW approx
!r          11: Read J(q) from disk and generate derivative properties
!r          12: Make G(R,R') (not implemented)
!r          13: Transverse exchange interactions J(q) from MST (not implemented)
!r          14: Longitudinal exchange interactions J(q) from MST, LW approx
!r          20: Transverse chi+- from ASA GF (not implemented)
!r          21: Read chi(q) from disk and print derivative properties
!r          24: chi++ from ASA GF
!r          25: chi-- from ASA GF
!r  ldlm    Switch for CPA (rdctrl)
!r          0  no CPA
!r          12 full CPA/DLM calculation: charge and omega both iterated
!r          32 restricted calculation: only omega iterated
!r          100s digit 1 => calculate Lloyd correction
!r  ldomg   (CPA) not documented (gfasa)
!r  ldos    1 make dos (rdctrl)
!r          2 generate weights for partial dos
!r          4 generate weights for m-decompos'n of pdos
!r          8 BZ map
!r  lfp     switches for full-potential (rdctrl)
!r          1 program uses full potential
!r          2 Have bndfp shorten q vectors before using them
!r          4 cause bndfp to retain hab and sab in s_pot%hab and s_pot%sab
!r          8 cause bndfp to update elind in s_ham%elind
!r  lfrce   How forces are calculated (rdctrl)
!r          0 do not calculate forces
!r          1 shift in FA density
!r          2 shift in core+nucleus
!r         10 added for shift to be screened
!r  lgen3   switches for third generation LMTO (rdctrl)
!r          1 3rd generation LMTO
!r          2 make structure constants by Ewald summation
!r  lham    1 (ASA) 2-center (rdctrl)
!r          1 (molecules) two-panel
!r          2 (ASA) 2-c + pert. corr
!r          4 (ASA) auto-down-fold
!r          8 (ASA) change rep interactively
!r         16 (ASA) suppress d hybridization
!r         32 (ASA) preserve ortho. evecs
!r         64 (ASA) save evecs to disk
!r        128 (ASA) gamma-rep
!r        256       use true spherical harmonics
!r        512       Read shfac file, if it exists
!r  lmet    1 metal
!r          2 tetrahedron (set in rdctrl)
!r          4 (GF) V-shift1 is zero
!r          8 (GF) V-shift2 is zero
!r         16 When using tetrahedron integration,
!r            replace weights for output density with sampling weights
!r  lncol   1 noncollinear magnetism (rdctrl)
!r          2 spin spirals
!r          4 spin-orbit coupling
!r          8 External magnetic field (B . sigma)
!r         16 mag. forces
!r         32 spin-orbit coupling, LzSz only
!r         64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
!r        128 spin-orbit coupling, entire L.S-LzSz pert
!r        256 External magnetic field, Bz Sigmaz only
!r        512 Resolve SO contribution to energy by site
!r       1024 Print out SO-weighted average of electric field
!r         NB: ham->lncol and ctrl->lncol should be duplicates
!r  lrsa      For the FP noncollinear case
!r         1s digit governs how noncollinearity is treated.
!r          0 Input density collinear along z.
!r            Noncollinear output density generated to calculate energy only.
!r            Only the z component is kept for self-consistency
!r          1 Rigid spin approximation.  Local input densities are
!r            constrained to align with supplied Euler angles.
!r            Output density treated as case 0, except spin is rigidly
!r            rotated to axis given by site Euler angles.
!r          2 Rigid spin approximation.  Like 1, except spin axis determined
!r            by average magnetization.
!r          3 Fully noncollinear density.
!r        10s digit determines how local z axis defining basis function is set.
!r          0 Site spin quantization axis defining local density may be parallel to M or -M
!r            whichever yields the largest projection.
!r            Local moments can be positive or negative.
!r          1 Site spin quantization axis defining local density is parallel to M
!r            This switch doesn't have effect unless 1s digit lrsa is at least 2.
!r            Local moments are positive by construction.
!r            Note: there can be a difference between the two because
!r            the core density is aligned with the quantization axis.
!r            Small noncollinearity => use 0 to adiabatically move from collinear case.
!r            Otherwise use 1 avoid discontinuous jumps with small changes in input conditions.
!r            In this case starting magnetic moments should all be positive to avoid
!r            ambiguities in the orientation (+M, -M) in the local rho.
!r            Alternatively you can set the 100s digit:
!r        ... 100s digit
!r          1 Take spin-average of core density
!r  loptc   used by optics package  (rdctrl); also copied to s_optic%loptic
!r           1 generate Im(eps) (spin 1 only, if LTET=3)
!r           2 generate Im(eps) spin 2 (LTET=3)
!r           8: Im eps for  nonequilibrium absorption : independent Ef for electrons and holes
!r           9: Im eps for  nonequilibrium emission : independent Ef for electrons and holes
!r          Add 20 for Second Harmonic Generation
!r          -1 generate joint DOS
!r          -2: generate joint DOS, spin 2 (use with LTET=3)
!r          -3: generate up-down joint DOS (LTET=3)
!r          -4: generate down-up joint DOS (LTET=3)
!r          -5: generate spin-up single DOS (LTET=3)
!r          -6: generate spin-dn single DOS (LTET=3)
!r          -8: Joint DOS for  nonequilibrium absorption : independent Ef for electrons and holes
!r          -9: Joint DOS for  nonequilibrium emission : independent Ef for electrons and holes
!r  lordn   1 Embedded GF  4 Vanderbilt  8 Embedded cluster (rdctrl)
!r  lpgf    control parameter in layer GF (rdctrl)
!r             lpgf(1)
!r          1: Generate diagonal GF and output density
!r          2: Find Fermi left- and right-bulk Ef
!r          3: Trace out band structure for spec'd energy mesh
!r             In this case make Im(z)=0.
!r          4: Ditto for right bulk.
!r          5: Calculate current through structure
!r          6: Special case of 5 ... not documented
!r          7: Calculate reflection at L and R interfaces
!r             lpgf(2), 1s digit
!r          0  Make GF using standard approach
!r          1  Use LU decomposition
!r             lpgf(2), 10s digit
!r          0  Make surface GF as though collinear, even if in nc mode (e.g. with SO coupling)
!r          1  Make surface GF as noncollinear
!r  lqp     1 do not add inversion 2 inverse iteration (rdctrl)
!r  lrel    1 scalar relativistic (rdctrl)
!r          2 fully relativistic
!r            Add 10 to specify that core is calculated fully relativistically
!r  lrs     controls input from restart file (rdctrl)
!r          0 overlap atomic densities
!r          1 Read from restart file
!r          2 Read from restart file, ascii mode
!r          3 Always overlap atomic densities, even after relaxation step
!r          4 Read from restart file, invoke smshft
!r          8 Write new density to restart file
!r         16 Write new density to restart file, ascii format
!r         32 read site positions from input file
!r         64 read starting fermi level from input file
!r        128 read starting pnu level from input file
!r        256 rotate local density after reading
!r  lscr    controls screening of (nout-nin) for charge mixing  (rdctrl)
!r          0 do nothing
!r          1 Make P0(0)
!r          2 Screen output q and ves
!r          3 Screen output ves only
!r          4 Use model response to screen output q
!r            Add 1 to combine mode 1 with another mode
!r            Add 10*k to compute intra-site contribution to
!r            vbare each kth iteration
!r            Add 100*k to compute response function only
!r            each kth iteration
!r  lrquad  Type of quadrature for radial mesh (rdctrl)
!r          0 3-point (Simpson''s rule)
!r          1 5-point (Boole''s rule)
!r          2 7-point
!r  lstonr  second digit for graphical output (rdctrl)
!r  lstr    Note: no longer used; see str->lshow,str->lequiv
!r          1 print strux   2 find equiv strux
!r  lsx     1 Calculate screened exchange sigma  (rdctrl)
!r         10's digit nonzero make vintra to include in sigma
!r  ltb     Parameters for tight-binding program tbe  (rdctrl)
!r          1 overlap        2 crystal-field     4 ovlp+CF
!r          8 add ebarLL    16 forces           32 fij
!r         64 not used     128 pressure        256 evdisc
!r        512 pair pot    1024 TrH & local E  2048 local rho
!r       2^12 Hubbard U   2^13 No Madelung    2^14 wgt avg U
!r       2^15 L>0 estat   2^16 disc read incr 2^17 gamma-pt
!r  lves    controls how electrostatic potential is calculated (rdctrl)
!r          1 take ves as input
!r          2 Add site vshft to ves
!r  lxcf    parameter defining XC functional (rdctrl)
!r         *Case lxcf>1000 : use functionals from the libxc library.
!r          Then lxcf is a compound of 2 digits, so that two functionals
!r          may be invoked (exchange and correlation separately)
!r             lxcf/2**16 = index to the first functional in libxc
!r             mod(lxcf,2**16) = index to the 2nd functional in libxc
!r         *Case lxcf<1000 : lxcf corresponds to the following:
!r          1s digit:
!r          1 for Ceperly-Alder
!r          2 for Barth-Hedin (ASW fit)
!r          3 for PW91
!r          4 for PBE (same as PW91?)
!r          100s digit for GGAs
!r          0 for LSDA
!r          1 for LMH
!r          2 for PW91
!r          3 for PBE
!r          4 for PBE with Becke exchange
!r  lwsig   (bndfp) special modes for reading/writing
!r          evecs or self-energies to disk
!r          1  Rotates sigm to LDA basis; saves in file
!r          2  Similar to lwsig=1, except
!r             low- and high- energy blocks replaced by diagonal parts
!r          3  Rotates sigm to LDA basis and replaces low- and high-
!r             energy blocks as in lwsig=2.  But sigm is rotated back
!r             to orbital basis before being written to file
!r          4  Writes evals,evecs of LDA hamiltonian to file
!r          5  Writes evals,evecs of full hamiltonian to file
!r          6  Similar to 5, but occupation numbers are also saved.
!r          7  evecs of hamiltonian only (for reading in ASCII)
!r          8  special mode causing lmf to qp and evecs from evec file, rather than create them
!r          9  matrix elements of LDA vxc in LDA basis
!r  maxit   max. no.  iterations in self-consistency cycle (rdctrl)
!r  maxmem  If > 0, upper limit to dynamic memory allocation,
!r  mdprm   switches for molecular statics/dynamics (rdctrl)
!r  mdprm   arg 1: 1 new dynamics (rdctrl)
!r                 2 restart dynamics
!r                 4 relax with conjugate gradients
!r                 5 relax with variable metric
!r                 6 relax with Broyden
!r          arg 2: statics: switch
!r                 1 read hessian matrix
!r                 dynamics:
!r                   time step
!r          arg 3: (stat) relaxation x-tolerance
!r                 (dyn)  temperature
!r          arg 4: (stat) relaxation g-tolerance
!r                 (dyn)  thermostat relaxation time
!r          arg 5: (stat) step length
!r                 (dyn)  total MD time
!r          arg 6: (stat) Remove hessian after this many steps
!r                 (dyn)  barostat relaxation time
!r          arg 7: (stat) -
!r                 (dyn)  external pressure
!r  modep.. which dimensions are periodic (rdctrl)
!r  mxcst   mixing constraint (lmasa)
!r  nccomp  total number of extra CPA components for all classes (mksym)
!r          nccomp = sum (s_spec%ncomp) over all classes
!r  ndlmcl  number of classes subject to CPA/DLM treatment (mksym)
!r  nbas    size of basis (rdctrl)
!r  nbasp   size of padded basis, used by layer GF lmpg (susite)
!r  nclass  number of classes (mksym)
!r  nesabc  mesh for determining empty spheres (rdctrl)
!r  nitmv   max number of mol-dynamics iterations (rdctrl)
!r  nl      1 + maximum augmentation l (rdctrl)
!r  nmap    number of maps of ASA class parameters (rdctrl)
!r          (now permanently zero)
!r  ncl     (O-N) number of clusters
!r  ncomp   (CPA) total number of CPA components for class ic
!r  nclasp  (pgf) number of classes including PL -1,npl
!r  nofgl   (pgf) number of left  PL for off-diagonal GF
!r  nofgr   (pgf) number of right PL for off-diagonal GF
!r  npadl   (pgf) number of basis atoms in PL -2,-1,0
!r  npadr   (pgf) number of basis atoms in PL npl-1,npl,npl+1
!r  npl     (pgf) number of principal layers (susite)
!r  nrc     number of sites in each class (mksym)
!r          classes may be extended to include elements not in basis,
!r          e.g. for layer GF, and for CPA
!r  nsite   number of sites, related to number of atoms  (rdctrl,susite)
!r  nspec   number of species (rdctrl)
!r  nspin   number of spins (rdctrl)
!r  nvario  number of variables to output (rdctrl)
!r  omax1   sphere overlap constraints, type 1 (rdctrl)
!r  omax2   sphere overlap constraints, type 2 (rdctrl)
!r  lekkl   affects how the band center-of-gravity is calculated (mkrout).
!r          0 ekkl from density matrix * hab (Methfessel style)
!r          1 ekkl from energy-weighted density matrix (Kotani style)
!r  pgfsl   (PGF) pgfsl(ib) = PL to which site ib corresponds (susite)
!r  pgfvl   (PGF) pgfvl(ib) = PL potential group (susite)
!r  pgord   (PGF) table of indices that sort basis (pgfset.f)
!r  pgplp   (PGF) PL-dependent indices (pgfset.f)
!r  plbnd   Controls what FP code lmf makes in band pass (rdctrl,evalqp,sugws)
!r          0 band pass, integrates over BZ for density, total energy
!r          1 band pass makes bands at qp tabulated in file, written to disk
!r          2 like 0, but no output density.  Bands stored in s_ctrl%evals
!r          3 like 2, but no attempt to integrate over BZ
!r         10 make potential and exit
!r  pos     lattice positions
!r  quit    Switch to terminate execution at selected break points (rdctrl)
!r          stop execution after:
!r          1 show      2 atom    4 bands     8 ham     16 dos   32 rho    64 pot
!r  rmax    sphere radius, by class (clsprm)
!r  rmaxes  upper limit to ES radius when finding new empty spheres (rdctrl)
!r  rmines  lower limit to ES radius when finding new empty spheres (rdctrl)
!r  refinegstol
!r          Used by lmpg, to iterate gs to convergence within refinegstol(1)
!r          refinegstol(1) converge gs to this tolerance
!r          refinegstol(2) holds maximum final RMS difference in gs
!r          refinegstol(3) holds  RMS difference in gs, 1st iteration
!r  sclwsr  scale wsr until reaching this fractional vol (rdctrl)
!r          10s digit used for assymetric treatment of ES:
!r           0 ES and other sites are treated symmetrically
!r           1 all sites with z>0 are resized first; then
!r             all sites are resized.
!r           2 all sites with z>0 are resized first; then
!r             the ES sites only are resized
!r  sdmod   spin dynamics mode (rdctrl)
!r  sdprm   spin dynamics parameters (rdctrl)
!r  sdxsi   Bulgac and Kusnezov global deamons (rdctrl,susite)
!r  size    size of this structure (not used)
!r  smalit  parameters for small iterations (rdctrl)
!r  stpnrho flag indicating whether Dirac solver has relativistic moments Q
!r          If false, solver does not have qnur
!r          If true,  solver does have qnur
!r  tdlm    (DLM) spin temperature; a parameter if Weiss field? (rdctrl)
!r  tol     tolerances for self-consistency (rdctrl)
!r          tol(1) = charge tolerance
!r          tol(2) = charge tolerance (not used now)
!r          tol(3) = energy tolerance
!r          tol(4) = tolerance in max HF correction to force
!r          tol(5) = current value of HF correction to force
!r  wsrmax  constraint on size of largest WS sphere (rdctrl)
!r  zbak    background charge and MT correction parameters (rdctrl,aiocls)
!r          zbak(1) = uniform background charge included
!r                    in electrostatics but not in xc pot.
!r          zbak(2) = charge used to make MT correction
!r                    (ASA only)
!r spid     (character) table of species labels
!r clabl    (character) table of class labels
! ----------------------------------------------------------------
      type str_ctrl
      sequence

      integer    ::   lasa
      integer    ::   lbas
      integer    ::   lbxc
      integer    ::   lcd
      integer    ::   lcgf
      integer    ::   ldlm
      integer    ::   ldomg
      integer    ::   ldos
      integer    ::   lfp
      integer    ::   lfrce
      integer    ::   lgen3
      integer    ::   lham
      integer    ::   lmet
      integer    ::   lncol
      integer    ::   lrsa
      integer    ::   loptc
      integer    ::   lordn
      integer    ::   lpgf(2)
      integer    ::   lqp
      integer    ::   lrel
      integer    ::   lrs
      integer    ::   lscr
      integer    ::   lrquad
      integer    ::   lstonr(3)
      integer    ::   lsx
      integer    ::   ltb
      integer    ::   lves
      integer    ::   lxcf
      integer    ::   lwsig
      integer    ::   maxit
      integer    ::   maxmem
      integer    ::   modep(3)
      integer    ::   nbas
      integer    ::   nbasp
      integer    ::   nccomp
      integer    ::   ncl
      integer    ::   nclasp
      integer    ::   nclass
      integer    ::   ndlmcl
      integer    ::   nesabc(3)
      integer    ::   nitmv
      integer    ::   nl
      integer    ::   nmap
      integer    ::   nofgl
      integer    ::   nofgr
      integer    ::   npadl
      integer    ::   npadr
      integer    ::   npl
      integer    ::   nsite
      integer    ::   nspec
      integer    ::   nspin
      integer    ::   nvario
      integer    ::   lekkl
      integer    ::   plbnd
      integer    ::   quit
      integer    ::   sdmod
      integer    ::   smalit(2)

      logical    ::   stpnrho
!     integer    ::   inull ! ensure 8-byte word boundary

      real(8)    ::   defm(6)
      real(8)    ::   elin
      real(8)    ::   mdprm(7)
      real(8)    ::   omax1(3)
      real(8)    ::   omax2(3)
      real(8)    ::   rmaxes
      real(8)    ::   rmines
      real(8)    ::   refinegstol(3)
      real(8)    ::   sclwsr
      real(8)    ::   sdprm(5)
      real(8)    ::   sdxsi(4)
      real(8)    ::   tdlm
      real(8)    ::   tol(5)
      real(8)    ::   wsrmax
      real(8)    ::   zbak(2)

! ... Integer pointer arrays
      integer, pointer :: cllst  (:)
      integer, pointer :: clp    (:)
      integer, pointer :: clssl  (:)
      integer, pointer :: group  (:)
      integer, pointer :: initc  (:)
      integer, pointer :: ics    (:)
      integer, pointer :: idcc   (:)
      integer, pointer :: ipc    (:)
      integer, pointer :: ipcp   (:)
      integer, pointer :: ips    (:)
      integer, pointer :: mxcst  (:)
      integer, pointer :: nrc    (:)
      integer, pointer :: ncomp  (:)
      integer, pointer :: pgfsl  (:)
      integer, pointer :: pgfvl  (:)
      integer, pointer :: pgord  (:)
      integer, pointer :: pgplp  (:)

! ... Character pointer arrays
      character*8, pointer :: spid(:)
      character*8, pointer :: clabl(:)

! ... Real pointer arrays
      real(8), pointer :: dclabl (:)
      real(8), pointer :: pos    (:,:)
      real(8), pointer :: rmax   (:)

      end type str_ctrl

! --- str_dmft ---
!- Parameters for DMFT interface
! ----------------------------------------------------------------
!r Correlated orbitals are grouped into subblocks; each is calculated separately from the others.
!r Self-energy sigma contains nicix inequivalent correlated subblocks.
!r Subblock cix is an (ndim,ndim,nsp) matrix, where ndim is cix-dependent.
!r Typically there is a single l per subblock, in which case ndim for that subblock is 2*l+1.
!r Two or more sites may point to the same subblock (equivalent either with a spin exchange or not).
!r There are ncix total correlated subblocks among all sites, nicix inequivalent ones.
!r
!r Only a fraction of matrix elements of the correlated block are calculated
!r (typically but not necessarily the diagonal matrix elements).
!r All matrix elements of DMFT sigma and delta are stored contiguously in a single array in
!r a compressed form; see iodsig for I/O. readindmfl makes information to compress/decompress sig.
!r See symdmftsig for an instance of how sig is compressed/decompressed.
!r
!r  iproj    Index labelling projector type
!r  knorm    0 => Haule style k-independent normalization
!r           1 => French style k-dependent normalization
!r  lsigim   True if sigma is to be calculated on the imaginary axis
!r  lmaxd    Largest l of a correlated block
!r  dcmode   Double counting mode
!r  ncatm    Number of correlated atoms
!r  ncix     Number of correlated blocks of l channels at all sites in the lattice
!r  nicix    Number of inequivalent correlated blocks
!r  ndsig    Total number of indpendent matrix elements in all correlated subblocks, including spin
!r           ndsig dimensions DMFT delta and sigma
!r  nzsig    Total number of nonzero matrix elements in all correlated subblocks, including spin
!r           nzsig dimensions iasig
!r  nkabc    the k-mesh used to make and embed DMFT projectors
!r  beta     inverse temperature, Ry
!r  gammac   Broadening for correlated orbitals
!r  nlohi    first, last of states to be included in projector
!r           NB: window can alternatively be specified by wlohi
!r  wlohi    energy window of states to be included in projector, relative to Ef
!r           NB: window can alternatively be specified by nlohi
!r  nomg     Number of frequencies.  For Matsubara integration, this and beta specify mesh
!r  omg      Frequency mesh, stored as real numbers.
!r           If lsigim=T, true omg is sqrt(-1)*omg
!r  sig      DMFT sigma array, (called siginp in the Haule WEIN2K code)
!r           It contains entries only for the ndsig matrix elements to be calculated by DMFT.
!r           sig is complex but declared as a real matrix, following Haule convention:
!r           sig = sig(nomg,2,ndsig)
!r           Information needed to map sig to a particular (m1,m2) element is contained in iasig.
!r           The l quantum number is given by dmft%l and the site index by dmft%ib
!r  iasig    contains orbital indices, relative to start of augmentation l block,
!r           for each nonzero matrix matrix element in sig
!r           For now:
!r           iasig(1,1:nzsig) = row index m1, which ranges betw/ 1 and 2l+1
!r           iasig(2,1:nzsig) = column index m2, which ranges betw/ 1 and 2l+1
!r           iasig(3,1:nzsig) = spin index
!r           iasig(4,1:nzsig) = index to element in compressed form of DMFT sigma or delta
!r           iasig(5,1:nzsig) = cix index
!r  ... In each of the arrays below, one element for each correlated subblock
!r  ib       site associated with correlated subblock
!r  icix     icix(1:ncix)| points to inequivalent cix block: 1<|icix|<nicix
!r           icix(i) > 0 => i is equivalent to  icix(i)
!r           icix(i) < 0 => i is equivalent to -icix(i) with spin up,down exchanged
!r  ig       index of group symetry (symgr in str_lat) from cix to icix(cix)
!r  rmat      matrix rotation of cubic harmonic corresponding to ig
!r  ... In each of the arrays below, one element for each inequivalent correlated subblock
!r  l        l quantum numbers for this subblock
!r  ncl      number of l quantum numbers for this subblock
!r  qsplit   defines meaning of orbitals in subblock (taken from qsplit in Haule's CTQMC code)
!r            0  average GF, non-correlated
!r            1  |j,mj> basis, no symmetry, except time reversal (-jz=jz)
!r           -1  |j,mj> basis, no symmetry, not even time reversal (-jz=jz)
!r            2  real harmonics basis, no symmetry, except spin (up=dn)
!r           -2  real harmonics basis, no symmetry, not even spin (up=dn)
!r            3  t2g orbitals  only
!r           -3  eg orbitals   only
!r            4  |j,mj>, only l-1/2 and l+1/2
!r            5  axial symmetry in real harmonics
!r            6  hexagonal symmetry in real harmonics
!r            7  cubic symmetry in real harmonics
!r            8  axial symmetry in real harmonics, up different than down
!r            9  hexagonal symmetry in real harmonics, up different than down
!r           10  cubic symmetry in real harmonics, up different then down
!r           11  |j,mj> basis, non-zero off diagonal elements
!r           12  real harmonics, non-zero off diagonal elements
!r           13  J_eff=1/2 basis for 5d ions, non-magnetic with symmetry
!r           14  J_eff=1/2 basis for 5d ions, no symmetry
!r  sigind   maps individual elements of sigma into 2D array
!r  cf       matrix rotating orbitals in driver's representation to spherical harmonics
!r  csize    Number of independent matrix elements calculated for the subblock -> nisigi
!r  nzsigi   Number of nonzero matrix elements in subblock
!r           Indices for subblock cix span nzsigi(cix-1)+1:nzsigi(cix) for cix=1:ncix
!r           Note: nzsigi(cix-1)+1:nzsigi(cix) gives the total number of nonzero elements in the subblock
!r  ndim     (1:nicix) Size of sigind matrix for each independent cix block
!r  umode    Defines meaning U and of elements in array u
!r           1s digit
!r           0 => u(:,1) = Hubbard U, u(:,2) = Hubbard J
!r           1 => u(:,1) = F0, u(:,2) = F2, u(:,3) = F4
!r           2 => u(:,1) = screening length, Yukawa potential
!r           10s digit
!r           0 => density-density only
!r           1 => Full matrix form
!r  nomf     number of fermionic frequencies for the susceptibility
!r  nomb     number of bosonic   frequencies for the susceptibility						      						      
!r  The following are NO LONGER USED
!r  ic       class index
!r  ip       Points to 1st equivalent orbital in list
!r  idd      Double counting prescription, DMFT
!r           0 => Not treated at DMFT level
!r           1 => Around mean-field limit
!r           2 => Fully localized limit
!r           3 => D.C. determined by fixed nominal charge
!r  idm      for each l, a sequence of (2l+1) digits defining
!r           m-dependent handling of correlation
!r  u        Hubbard parameters.  Meaning depends on umode, above
!r  ndc      (dcmode=2) nominal charge that fixes double counting
! ----------------------------------------------------------------
      type str_dmft
      sequence

      integer :: iproj    ! index defining projector type.
      integer :: knorm    ! whether projectors are normalized at each k or not
      logical :: lsigim   ! True if frequencies are imaginary

      integer :: dcmode   ! Double counting style
      integer :: lmaxd    ! Largest l for DMFT-correlated orbitals
      integer :: ncatom   ! Number of correlated atoms
      integer :: ncix     ! Number of correlated blocks at all sites in the lattice
      integer :: nicix    ! Number of inequivalent correlated blocks
      integer :: ndsig    ! # of inequivalent matrix elements calculated by DMFT
      integer :: nzsig    ! total # of matrix elements calculated by DMFT
      integer :: nomg     ! Number of points for frequency mesh
      integer :: nlohi(2) ! First, last band used to calculate projector
      integer :: nkabc(3) ! the DMFT k-mesh for the bath
      integer :: nomb
      integer :: nomf						      

!     integer :: inull    ! pad for 8-byte word boundary

      real(8) :: beta     ! Inverse temperature, Ry
      real(8) :: wlohi(2) ! omegamin, omegamax for projector
      real(8) :: gammac   ! Broadening for correlated orbitals
!     real(8) :: omegamin ! Minimum frequency
!     real(8) :: omegamax ! Maximum frequency
      real(8), pointer  :: omg(:) ! Frequency mesh

! ... In each of the arrays below, one element for each correlated subblock
      integer, pointer   :: ib(:)      ! site associated with correlated subblock
      integer, pointer   :: icix(:)    ! which inequivalent block: this site belongs to
      integer, pointer   :: ig(:)
      real(8), pointer   :: rmat(:,:,:)
! ... In each of the arrays below, one element for each inequivalent correlated subblock
      integer, pointer   :: l(:)       ! l quantum number.  should have 2nd dimension in more l's than 1
      integer, pointer   :: ndim(:)    ! (1:nicix) Size of sigind matrix for each independent cix block
      integer, pointer   :: ncl(:)     ! Number of l's in this subblock
      integer, pointer   :: nzsigi(:)  ! index to last matrix element in this subblock
      integer, pointer   :: sigind(:,:,:,:) ! maps elements in sig into 2D array
!     integer, pointer   :: sigind(:,:,:) ! maps elements in sig into 2D array
      integer, pointer   :: iasig(:,:)  ! contains orbital and cix info for each ME in sig
      integer, pointer   :: qsplit(:)  ! specifies meaning of orbital basis in subblock
      real(8), pointer   :: sig(:,:)   ! Matrix elements of DMFT sigma
      integer, pointer   :: umode(:)   ! Defines meaning of elements in U,J
      complex(8),pointer :: cf(:,:,:)  ! Rotation of basis used by driver to spherical harmonics
      complex(8), pointer:: Olapp(:,:,:,:) ! Projector overlap matrix

!     Not used now
!     integer, pointer   :: ic(:)      ! class index
!     integer, pointer   :: ip(:)      ! Points to 1st equivalent orbital in table
!     integer, pointer   :: idd(:)     ! Double counting prescription, DMFT
!     integer, pointer   :: idm(:)     ! markers for m-dependent treatment of correlation
!     real(8), pointer   :: u(:,:) ! Hubbard parameters U,J,.. (remove?)
!     real(8), pointer   :: ndc(:) ! For DC determined by fixed nominal charge (remove?)

      end type str_dmft

! --- str_gw ---
!- GW and dielectric response; also used to hold DMFT parameters
! ----------------------------------------------------------------
!r  element purpose
!r  code     specify for which GW code to create input:
!r           0 fpgw (v033a5)
!r           1 spex (v02.04)
!r           2 fpgw (Sep 2012)
!r  delre    spacing on real axis for energy contour
!r           delre(1) = dw(0)
!r           delre(2) = omegac (coff to parabolic term)
!r  deltaw   width in finite diff for sigma energy derivative
!r  deltax   Im. E (broadening delta) in computation of x0
!r  ecuts    energy cutoff for calculating sigma_nn'
!r  ecutpb   (sugwin) include core states in product basis with energy > ecutpb(1)
!r           (sugwin) include core states in product basis with charge > ecutpb(2) ... not used
!r  eseavr   average diagonal sigma above energy cutoff
!r  gcutb    G-vector cutoff for basis
!r  gcutx    G-vector cutoff for response function
!r  gsmear   smearing in pole of Green's function
!r  lsigint  0 Standard mode
!r           1 Use eigenvectors from LMTO part of LDA hamiltonian
!              for sigm interpolation
!r  lgw      specifies mode for GW computation
!r           1s digit
!r           0 standard mode
!r           1 Faleev's fast energy integration
!r  lmxpb    default l-cutoff for product basis
!r  lshft    shift mesh to straddle q=0.  Component for each axis
!r  mksig    Mode for generating self-energy sigma matrix
!r           0 do not create sigma matrix
!r           1 Make sigma_nn'(E_f)
!r           2 Make sigma_nn'(E_n/2+E_n'/2)
!r           3 Make sigma_nn'(E_n)/2 + sigma_nn'(E_n')/2
!r           lsig>0 implies 1s digit lgw>0
!r  nband    number of bands in perturbation calculation
!r  nime     # energies along Im axis for energy contour
!r  nkabc    the number of k-points for GW calculation
!r  pb1      default product basis 1
!r  pb2      default product basis 2
!r  pbtol    product-basis overlap cutoff
!r  qoffp    qp offset parameter
!r  rdsig    parameters for how to read, interpolate SE matrix
!r  usebsew  whether to use Bethe-Salpeter W
! ----------------------------------------------------------------
      type str_gw
      sequence

      integer    ::   code
      integer    ::   lgw
!     integer    ::   lsigint
      integer    ::   lmxpb
      integer    ::   lshft(3)
      integer    ::   mksig
      integer    ::   nband
      integer    ::   nime
      integer    ::   nkabc(3)
      integer    ::   rdsig
      integer    ::   usebsew ! alignment ok now without inull below
      !integer    ::   inull ! ensure 8-byte word boundary

      character*8::   pb1
      character*8::   pb2

      real(8)    ::   delre(2)
      real(8)    ::   deltaw
      real(8)    ::   deltax
      real(8)    ::   ecuts
      real(8)    ::   ecutpb
      real(8)    ::   eseavr(2)
      real(8)    ::   gcutb
      real(8)    ::   gcutx
      real(8)    ::   gsmear
      real(8)    ::   pbtol(8)
      real(8)    ::   qoffp


      end type str_gw

!r --- str_ham ---
!- Hamiltonian parameters
! ----------------------------------------------------------------
!r  element  purpose
!r  alfsi    stability factor alfsi*(S+)*S (molecules)
!r  amgm     magnetization
!r  bandw    Maximum size of off-diagonal, packed storage
!r  basopt   Parameters used generation of LMTO basis parms
!r            1 = autob (used in autoset basis)
!r            1s   digit 1 or 3 (lmfa) Autogenerate RSMH,EH
!r                       2 or 4 (lmfa) Autogenerate RSMH,EH, RSMH2,EH2
!r                       1 or 3 (lmf)  Read RSMH,EH,RSMH from basp file
!r                       2 or 4 (lmf)  READ RSMH,EH, RSMH2,EH2 from basp file
!r            10s  digit 1 (lmfa) Find and estimate PZ from free atom wf
!r                         (lmf)  Read P from basp file
!r            100s digit 1 (lmfa) Estimate P from free atom wf
!r                         (lmf)  Read P from basp file
!r            2 = global cutoff to lmxb, 1st kappa
!r                (used in autogenerating parameters)
!r            3 = global cutoff to lmxb, 2nd kappa
!r                (used in autogenerating parameters)
!r            4 = rsmmx: maximum rsm, in units of rmt
!r                (used in autogenerating parameters)
!r            5 = ehmx : maximum eh, Ry
!r                (used in autogenerating parameters)
!r            6 = esprd: default spread in EH,EH2
!r                (used in autogenerating EH,EH2)
!r            7 = modeV: specifies a mode to modify potential
!r                (used in autogenerating EH,EH2)
!r            8 = vbar : parameter in V modification
!r                (used in autogenerating EH,EH2)
!r            9 = pqrmx : Set local orbital PZ set when q(r>rmt)>pqrmx
!r           10 = 1s  digit LMTO basis parameters
!r                 0 standard basis, typically spd
!r                 1 hyperminimal basis
!r                 2 hyperminimal basis + 1 higher l
!r                 3 same as 0
!r                 4 increment standard basis by adding 1 l
!r               10s digit
!r                 1 Set defaults for GW calculation
!r              100s digit
!r                 0 Traditional mode
!r                 1 Jerome's experimental mode
!r            ... when basopt(10)'s 100s digit set (function fabaspj)
!r           11 = deep energy for V(r)/ctp basis scheme
!r           12 = shallow energy for V(r)/ctp basis scheme
!r           13 = Hankel energy for V(r)/ctp basis scheme
!r  dabc     Spacing for real-space mesh (molecules)
!r  ehf      Harris-Foulkes energy
!r  ehk      Hohnberg-Kohn energy
!r  elind    Lindhard screening parameter
!r  eseavr   energy-weighted average of self-energy, high states
!r  eterms   terms making up the total energy.   For FP:
!r            1  =  ehf
!r            2  =  eks
!r            3  =  utot
!r            4  =  rhoves  int (rhoval x Ves) (true-sm)
!r            5  =  cpnves  int (core+nucleus x Ves) (true-sm)
!r            6  =  rhoeps  int (rho * exc)
!r            7  =  rhomu   int (rho * vxc)
!r            8  =  sumec   sum of foca=0 core eigenvalues
!r            9  =  sumtc   sum of all core kinetic energies
!r            10 =  rhcvef1 int (rhoc*vef1) for sites w/ lfoca=0
!r            11 =  valvef  int (rhov1*vef1 - rhov2*vef2)
!r            12 =  sumt0   sum of foca core energies
!r            13 =  sumev   band structure energy
!r            14 =          not used
!r            15 =  amom    system magnetic moment
!r            16 =  sumev   sum of single-particle evals
!r            17 =  rinvxt  input density * external pot
!r            18 =  rouvxt  output density * external pot
!r            19 =  rhosig  trace of self-energy over occ states
!r            20 =  external bfield * local moment
!r            21 =  constraining bfield * local moment
!r  etrms    Analog of eterms, by resolved by class and used by CPA (dlminit)
!r  evals    Holds eigenvalues of hamiltonian
!r  hord     order of polynomial approximation to ham.
!r           In 2nd gen LMTO, relevant only in GF context.
!r  jpo      For jigsaw puzzle orbitals
!r  kmto     envelope kinetic energies making up chord LMTO
!r  lasa     not used; see ctrl lasa
!r  ldham    vector describing hamiltonian dimensions:
!r           1: ldim   = dimension of lmto basis
!r           2: lidim  = ldim + size of downfolding block
!r           3: lidhim = lidim + size of higher block
!r           4: nspc   = number of coupled spins
!r                     = 1 unless noncollinear magnetism
!r           5: ldimc  = ldim * nspc
!r           6: lidimc = lidim * nspc
!r           7: lihdimc = lihdim * nspc
!r           8: nspx   = number of separate spin channels
!r                     = nsp if nspc is not 2; else 1
!r           9: pfdim  = leading dimension of potential function arrays.
!r                       Special dimensioning for CPA case; see makidx.f
!r  lgen3     shouldn't be used; just a copy of ctrl->lgen3
!r  lham      1 (ASA) 2-center
!r            1 (molecules) two-panel
!r            2 (ASA) 2-c + pert. corr
!r            4 (ASA) auto-down-fold
!r            8 (ASA) change rep interactively
!r           16 (ASA) suppress d hybridization
!r           32 (ASA) preserve ortho. evecs
!r           64 (ASA) save evecs to disk
!r          128 (ASA) gamma-rep
!r          256       use true spherical harmonics
!r          512       not used
!r         1024       not used
!r              NB: ham->lham and ctrl->lham should be duplicates
!r  lmaxu     dimensioning parameter for LDA+U potential
!r  lmxax     largest augmentation lmax in basis
!r  lncol     1 noncollinear magnetism (rdctrl)
!r            2 spin spirals
!r            4 spin-orbit coupling
!r            8 External magnetic field (B . sigma)
!r           16 mag. forces
!r           32 spin-orbit coupling, LzSz only
!r           64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
!r          128 spin-orbit coupling, entire L.S-LzSz pert
!r          256 External magnetic field, Bz Sigmaz only
!r          512 Resolve SO contribution to energy by site
!r         1024 Print out SO-weighted average of electric field
!r           NB: ham->lncol and ctrl->lncol should be duplicates
!r  lrsa    Applies to FP noncollinear case
!r          ... bits 0 and 1:
!r            0 Input density collinear along z.
!r              Noncollinear output density generated to calculate energy only.
!r              It is rotated to z before mixing with input density.
!r            1 Rigid spin approximation.  Local input densities
!r              constrained to align with supplied Euler angles.
!r              Output density treated as case 0, except rotation
!r              is to user-supplied axis given through Euler angles.
!r            2 Rigid spin approximation.  Like 1, execpt axis determined
!r              by average magnetization.
!r            3 Fully noncollinear density
!r          ... bit 2 determines how local z axis defining basis function is set.
!r            4 If set, local z axis determined by direction of |M|.
!r              Local moments are positive by construction.
!r              Otherwise local z axis determined by direction of M or -M,
!r              whichever yields the largest projection.
!r              Local moments can be positive or negative.
!r              Note: there can be a difference between the two because
!r              the core density is aligned with the quantization axis.
!r              If noncollinearity is small, set bit to 0 to
!r              adiabatically move from collinear case.
!r              Otherwise set bit to 4 to avoid discontinuous jumps in the axis
!r              with small changes in input conditions.
!r              User is advised also to set bit 8
!r            8 Take spin-average of core density
!r           NB: ham->lrsa and ctrl->lrsa should be duplicates
!r  lsig     parameters concerning self-energy
!r            1s digit
!r            0 do not read Sigma'//
!r            1 read sigma, add to potential
!r            Add 2 to symmetrize Sigma(T)
!r            Add 4 to include Re(Sigma(T)) only
!r            10s digit
!r            0 simple interpolation
!r            1 approx high and low sigma by diagonal (see rdsigm)
!r            3 interpolate sigma by LDA evecs ... no longer supported
!r  ltb      not used; see ctrl->ltb
!r  lxcf     parameter defining local XC functional.
!r           Not used: see ctrl->lxcf
!r  nbf      number of channels for magnetic field
!r  ndham    Largest dimension of lmto+PW basis
!r  ndhrs    dimension of site block for r.s. h or sigma
!r  ndofH    leading dimension to offH
!r  neula    number of channels for euler angles
!r           1 -> No euler angles defined
!r  nkaph    number of repetitions of one l-quant.n. in basis
!r  nlibu    Number of LDA+U blocks
!r  nlmto    Dimension of LMTO part of basis
!r  nmto     Polynomial order to NMTO
!r  npwmin   Estimate for minimum number of PWs to be calculated
!r           (PW basis is q-dependent; max size not known a priori)
!r  npwpad   Padding to be added to estimated max basis dimension
!r           and subtracted from min basis dimension
!r           (PW basis is q-dependent; max size not known a priori)
!r  nqsig    Number of k-points for which self-energy known
!r           (Used for interpolating between qpoints)
!r  bdots    mag. field in local coordinates, dotted with Pauli matrices (suham)
!r  eula     Euler angles (susite)
!r  hrs      r.s. h or sigma (rdsigm,bndasa)
!r  iaxs     neighbor table for self-energy (rdsigm)
!r  iprmb    permutation table arranging order of basis functions in hamiltonian (suham,makidx)
!r  lmxa     vector of augmentation l-cutoffs
!r  magf     external magnetic field (susite,iomagf)
!r  nprs     number-of-pairs table for self-energy (rdsigm)
!r  offH     table of offsets to hamiltonian matrix (suham, makidx, pgfasa)
!r  qsig     list of qp at which sigma can be computed
!r  oveps    When diagonalizing hamiltonian, discard part of hibert space
!r           corresponding to evals of overlap < oveps
!r  ovncut   When diagonalizing hamiltonian, discard part of hibert space
!r           corresponding to first ovncut evals of overlap
!r  pmax     global minimum allowed values for pnu
!r  pmin     global minimum allowed values for pnu
!r  pnudef   mode that controls default value of pnu
!r           1s digit
!r           0 Use version 6 defaults
!r           1 Use defaults tailored for LDA
!r           2 Use defaults tailored for GW
!r           10s digit for spin averaging pnz
!r           0 pnz (low) is allowed to float and is spin-dependent
!r           1 pnz (low) is allowed to float and is spin-averaged
!r           2 pnz is frozen
!r           4 add to retain pre-7.14 compatibility :
!r             pnz from spin 1 is used to semicore envelope parms (see loctsh)
!r             Default beginning 18 Oct 2018: pnz calculated from spin average of pnz
!r             Supersedes prior meaning of this switch, which additionally used spin-1 potential
!r  pwemax   High Energy cutoff for PW part of basis
!r  pwemin   Low Energy cutoff for PW part of basis
!r  pwmode   Controls PW part of basis
!r           0 => no PW part of basis
!r           1 => include PWs in basis
!r           2 => include only PWs in basis
!r  qpoff    qp offset when generating transformation of self-energy
!r  qss      spin spiral
!r  rdvext   Switch for reading parameters for external potential
!r  rsrnge   cutoff in connecting vector length for r.s. sigma
!r  rsstol   tolerance in Bloch sum error for r.s. sigma
!r  seref    Sum of reference energies
!r  sigp     parameters for approximating self-energy sigma
!r           arg 1: mode : specifies how to set its diagonal part
!r                  for states above the high-energy cutoff
!r                  0 constrain sigii to be > asig+bsig*e
!r                  1 constrain sigii to be = asig+bsig*e
!r                  2 constrain sigii to be > asig and < bsig
!r                  3 constraint same as mode 1.
!r                    Mode 3 differs in that the least-squares fit to
!r                    sigii (for informational purposes only, to help
!r                    estimate asig and bsig) is done for states between
!r                    efit and nmax or emax
!r                  4 constrain sigii to be eseavr, read from sigm file
!r           arg 2: nmin : sigma for states 1..nmin are approximated by sigii
!r           arg 3: emin : (used only if nmin<0)
!r                       : sigma for levels e<emin are approximated by sigii
!r           arg 4: nmax : sigma for levels i>nmax are approximated by
!r                         sigii AND constrained according to mode
!r           arg 5: emax : (used only if nmax<=0)
!r                       : sigma for levels e<emax are approximated by
!r                         sigii AND constrained according to mode
!r           arg 6: asig : constraint used to approximate
!r                         sigii = asig + E * bsig  or
!r                         asig < sigii < bsig
!r           arg 7: bsig : constraint used to approximate
!r                         sigii = asig + E * bsig  or
!r                         asig < sigii < bsig
!r           arg 8: efit : (mode 3) energy minimium (not used?)
!r                         for fitting asig and bsig
!r           arg 9: ? lwrite:T, write sigii to file (not used?)
!r           10 not used now
!r  size     size of this structure
!r  thrpv    3 PV for cell
!r  udiag    diagonal-only LDA+U
! ----------------------------------------------------------------
      type str_ham
      sequence

      integer    ::  bandw
      integer    ::  hord
      integer    ::  jpo
      integer    ::  lasa
      integer    ::  ldham  (16)
      integer    ::  lgen3
      integer    ::  lham
      integer    ::  lmaxu
      integer    ::  lmxax
      integer    ::  lncol
      integer    ::  lrsa
      integer    ::  lsig
      integer    ::  ltb
      integer    ::  lxcf
      integer    ::  nbf
      integer    ::  ndham
      integer    ::  ndhrs
      integer    ::  ndofH
      integer    ::  neula
      integer    ::  nkaph
      integer    ::  nlibu
      integer    ::  nlmto
      integer    ::  nmto
      integer    ::  npwmin
      integer    ::  npwpad
      integer    ::  nqsig
      integer    ::  ovncut
      integer    ::  pnudef
      integer    ::  pwmode
      integer    ::  rdvext
      integer    ::  udiag
!     integer    ::  inull ! ensure 8-byte word boundary

      real(8)    ::  alfsi
      real(8)    ::  amgm
      real(8)    ::  basopt (15)
      real(8)    ::  dabc   (3)
      real(8)    ::  ehf
      real(8)    ::  ehk
      real(8)    ::  elind
      real(8)    ::  entropy
      real(8)    ::  eseavr (2)
      real(8)    ::  eterms (22)
      real(8)    ::  kmto   (6)
      real(8)    ::  oveps
      real(8)    ::  pmax   (10)
      real(8)    ::  pmin   (10)

      real(8)    ::  pwemax
      real(8)    ::  pwemin
      real(8)    ::  qpoff  (3)
      real(8)    ::  qss    (4)
      real(8)    ::  rsrnge
      real(8)    ::  rsstol
      real(8)    ::  seref
      real(8)    ::  sigp   (10)
      real(8)    ::  thrpv

!     Pointer arrays
      integer, pointer::    iaxs(:)
      integer, pointer::    iprmb(:)
      integer, pointer::    lmxa(:)
      integer, pointer::    nprs(:)
      integer, pointer::    offH(:,:)

      real(8), pointer::    etrms(:)
      real(8), pointer::    eula(:,:)
      real(8), pointer::    hrs(:)
      real(8), pointer::    magf(:)
      real(8), pointer::    qsig(:,:)
      real(8), pointer::    evals(:,:,:)

      complex(8), pointer:: bdots(:)

      end type str_ham

!r --- str_lat ---
!- Lattice and parameters related to structure
! ----------------------------------------------------------------
!r  Element  Purpose
!r   ag      pointer to symmetry group translations
!r   alat    lattice parameter, in a.u.
!r   as      dimensionless Ewald smoothing parameter
!r   avw     average MT radius, in a.u.
!r   awald   Ewald smoothing parameter (lattic)
!r   bgv     phase factor sum for symmetrization of mesh rho (supot)
!r   cg      pointer to Clebsch Gordan coeffs (setcg)
!r   cy      pointer to Ylm normalization constants (setcg)
!r   dist    deformation parameters (lattic)
!r   dlv     pointer to direct lattice vector (lattic)
!r   gam     lattice shear parms: gam, gx,gy,gz
!r   gmax    cutoff gmax for Fourier transforms
!r   gv      pointer to list of F.T. G vectors (supot)
!r   gvq     pointer to q-dependent G-vectors (sugvec)
!r   indxcg  pointer to Clebsch Gordan indxcg (setcg)
!r   igv     pointer to q-dependent G-vectors, in units of qlat (sugvec)
!r   igv2    analog of igv, but transpose, used for APWs (sugvec)
!r   ips0    pointer to first vec in star (for symm mesh rho) (supot)
!r   istab   pointer to site permutations table for group ops
!r   jcg     pointer to Clebsch Gordan jcg (setcg)
!r   kv      pointer to indices in list of F.T. G vectors (supot->gvlst)
!r   kv2     Analog of kv, for APWs (sugvec)
!r   ldist   switch specifying what kind of dist (lattdf)
!r   lsym    0 no special symmetry considerations
!r           bits 0,1:
!r           1 Make symops including SOC with spin quantized along z
!r           2 Allow operations that take z to -z
!r           bit 2
!r           4 When symmetrizing over k in the noncollinear case,
!r             do not include spinor rotation
!r   lmxst   L=cutoff when setting tolerance in Ewald sums.  If 0 => use internal default
!r   nabc    no. divisions for F.T. mesh
!r   ng      no. G vectors
!r   napw    no. G vectors for APW
!r   ngq     -
!r   nkd     no. direct latt. vecs. for Ewald sum (lattic)
!r   nkdmx   dimensioning for arrays holding latt. vecs
!r   nkq     no. reciprocal latt. vecs. for Ewald sum (lattic)
!r   nkqmx   dimensioning for arrays holding latt. vecs
!r   npgrp   Number of point symmetry group operations
!r   nsgrp   Number of space symmetry group operations
!r   nsafm   Column in (symgr, ag, istab) for special AFM symmetry, if any
!r            0 => no special AFM symmetry
!r           >0 => special AFM symmetry, but both spins calculated
!r           <0 => special AFM symmetry, spin 2 derived from spin 1
!r   plat    lattice vectors, units of alat (lattic)
!r   plat0   lattice vectors before distortion (lattic)
!r   plat2   secondary lattice vecs used in various contexts
!r   plate   order-N
!r   platl   pgf (lattic)
!r   platr   pgf (lattic)
!r   pos     pointer to site positions (susite)
!r   qlat    reciprocal lattice vectors, units 2pi/a (lattic)
!r   qlv     pointer to Ewald reciprocal lattice vectors (lattic)
!r   rpad    truncate Ewald to rpad*rmax when lattice vector
!r           list has to be padded in order to include at
!r           least one lattice vector
!r   slat    superlattice vectors
!r   symgr   pointer to symmetry group rotation matrices
!r   tol     Ewald tolerance
!r   tolft   FT mesh tolerance
!r   vol     cell volume
! ----------------------------------------------------------------
      type str_lat
      sequence

      integer    ::  ldist
      integer    ::  lmxst = 0
      integer    ::  lsym
      integer    ::  nabc  (3)
      integer    ::  ng
      integer    ::  napw
      integer    ::  ngq
      integer    ::  ngvcc ! number of lattice vectors igvcc
      integer    ::  nkd
      integer    ::  nkdmx
      integer    ::  nkq
      integer    ::  nkqmx
      integer    ::  npgrp
      integer    ::  nsgrp
      integer    ::  nsafm
      integer    ::  inull ! ensure 8-byte word boundary

!     real(8)    ::  afmt  (3)
      real(8)    ::  alat
      real(8)    ::  as
      real(8)    ::  avw
      real(8)    ::  awald
      real(8)    ::  dist  (3,3)
      real(8)    ::  gam   (4)
      real(8)    ::  gmax

      real(8)    ::  plat    (3,3)
      real(8)    ::  plat0   (3,3)
      real(8)    ::  plat2   (3,3)
      real(8)    ::  plate   (3,3)
      real(8)    ::  platl   (3,3)
      real(8)    ::  platr   (3,3)
      real(8)    ::  qlat    (3,3)
      real(8)    ::  rpad
      real(8)    ::  slat    (3,3)
      real(8)    ::  tol
      real(8)    ::  tolft
      real(8)    ::  vol

! ... Integer pointer arrays
      integer, pointer ::  indxcg(:)
      integer, pointer ::  ips0(:)
      integer, pointer ::  istab(:,:)
      integer, pointer ::  jcg(:)
      integer, pointer ::  kv(:,:)
      integer, pointer ::  igv(:,:)
      integer, pointer ::  igv2(:,:)
      integer, pointer ::  igvcc(:,:)
      integer, pointer ::  kv2(:,:)

! ... Real pointer arrays
      real(8), pointer::   ag(:,:)
      complex(8),pointer:: bgv(:)
      real(8), pointer ::  cg(:)
      real(8), pointer ::  cy(:)
      real(8), pointer ::  dlv(:,:)
      real(8), pointer ::  gv(:,:)
      real(8), pointer ::  pos(:,:)
      real(8), pointer ::  gvq(:,:)
      real(8), pointer ::  qlv(:,:)
      real(8), pointer ::  symgr(:,:)

      type(str_symops),pointer ::  s_sym(:)

      end type str_lat

! --- str_mix ---
!- Parameters for mixing charge density and other quantities
! ----------------------------------------------------------------
!r  element purpose
!r  b       mixing beta
!r  bl      previous mixing beta
!r  bv      extra potential mixing
!r  elind   Lindhard energy for model screening
!r  fn      mixing file name
!r  kill    kill the mixing file after k iterations
!r  lxpot   decouple potential and the charge
!r          1: mix vin and v(qmix)
!r          2: mix vin and v(qout)
!r  mmix    maximum number to mix
!r  mode    1 Anderson 2 Broyden
!r  model   previous mixing mode
!r  n       Not used?
!r  nitu    max number of LDA+U itreations
!r  nmix    actual number mixed
!r  nsave   number of iterations to save on disk
!r  r..     expression for rmscst
!r  rms1    1st rms error
!r  rms2    2nd rms error
!r  tj..    Anderson t's
!r  tolu    tolerance for LDA+U
!r  umix    mixing parameter for LDA+U
!r  w..     Linear mixing weights
!r  wc      Broyden weight
! ----------------------------------------------------------------
      type str_mix
      sequence

      integer    ::   kill
      integer    ::   lxpot
      integer    ::   mmix
      integer    ::   mode
      integer    ::   model
      integer    ::   nitu
      integer    ::   nmix
      integer    ::   nsave
!     integer    ::   inull ! ensure 8-byte word boundary
!     integer    ::   n not used

      character*8::   fn
!     character*8::   r(3) !string not used

      real(8)    ::   b
      real(8)    ::   bl
      real(8)    ::   bv
      real(8)    ::   elind
      real(8)    ::   rms1
      real(8)    ::   rms2
      real(8)    ::   tj(10)
      real(8)    ::   tolu
      real(8)    ::   umix
      real(8)    ::   w(3)
      real(8)    ::   wc
      end type str_mix

! --- str_move ---
!- Parameters for molecular and spin dynamics
! ----------------------------------------------------------------
!r  element purpose
!r  ct      coefficients for global deamons thermostat
!r  gyro    gyromagnetic ratio (magnetic dynamics)
!r  kt      temperature, atomic units
!r  modt    mode for thermostat
!r  nmodt   number of thermostat modes
!r  prmint  Parameters for numerical integration
!r          For Bulirsch-Stoer integration:
!r          1:   mode: (1 for BS-integration)
!r          2:   rst:  1 for start new integration, 0 to continue
!r          3:   ts0:  minimum time step size
!r          4:   tol:  tolerance in integration errors
!r          5:   mx:   order of rational function extrapolation
!r          6:   mi:   number of midpoint rules
!r          7-17 nseq: sequence of no. midpoint divisions
!r          18:        offset to bs workspace
!r  size
!r  tnow    duration of simulation to this point, in units of ts
!r  ts      suggested time step size
!r  tsequ   initial equilibration time before stats accumulated
!r  tstot   duration of total simulation, in units of ts
!r          In SQS finder context, sign of tstot is used as a flag
!r          tstot<0 => flag to end simulation on minimum energy
! ----------------------------------------------------------------
      type str_move
      sequence

      integer    ::   modt(3)
      integer    ::   nmodt
!     integer    ::   inull ! ensure 8-byte word boundary

      real(8)    ::   ct(3)
      real(8)    ::   gyro
      real(8)    ::   kt
      real(8)    ::   prmint(20)
      real(8)    ::   tnow
      real(8)    ::   ts
      real(8)    ::   tsequ
      real(8)    ::   tstot
      end type str_move

! --- str_optic ---
!- Parameters for optics
! ----------------------------------------------------------------
!r  element  purpose
!r  axes     (abc) axes in xi^(abc) for nonlinear optics
!r  cll      core level l quantum number
!r  cln      core level n quantum number
!r  cls      1 or 2 for EELS or XANES
!r  clsite   site id for core level
!r  dw       energy mesh spacing
!r  esciss   Energy shift (scissors operator)
!r  gradrm   matrix elements of the radial part of grad, phi,phidot
!r  gradrsl  matrix elements of the radial part of grad, PkL
!r  optme    Optical matrix elements |<i | grad | j>|
!r           connecting occ states i to unocc states j for 1 k
!r  optmt    Optical matrix elements |<i | grad | j>|^2
!r           connecting occ states i to unocc states j
!r           optmt is symmetrized over group rotations
!r           Note: not used now
!r  ffmt     file format for writing optics
!r           1s digit for opt file
!r           0   6-digit floating format, e.g.
!r               0.011200    17.606809    17.606809    17.606809
!r           1   6-digit exponential form
!r               1.120000E-02  1.760681E+01  1.760681E+01  1.760681E+01
!r  gradrm   radial matrix elements of gradient operator and phi,phidot
!r  gradrsl  radial matrix elements of gradient operator and polynomials PkL
!r  gradsm   matrix elements the smoothed part of the eigenfunctions (not used)
!r  iq       q vector for chi(q), in multiples of qlat/nkabc
!r  imref    quasi-fermi levels: imref(1) for holes, imref(2) for electrons
!r  kmxax    Largest polynomial order: dimensions gradrsl
!r  kt       temperature for finite temperature case
!r  loptic   What optical property to calculate
!r           Should be same as s_ctrl%loptc
!r           1 generate Im(eps) (spin 1 only, if LTET=3)
!r           2 generate Im(eps) spin 2 (LTET=3)
!r           8: Im eps for  nonequilibrium absorption : independent Ef for electrons and holes
!r              Transition probability scaled by f(E-imrefp) * [1-f(E-imrefn)]
!r           9: Im eps for  nonequilibrium emission : independent Ef for electrons and holes
!r              Transition probability scaled by [1-f(E-imrefp)] * f(E-imrefn)
!r          Add 20 for Second Harmonic Generation
!r          -1 generate joint DOS
!r          -2: generate joint DOS, spin 2 (use with LTET=3)
!r          -3: generate up-down joint DOS (LTET=3)
!r          -4: generate down-up joint DOS (LTET=3)
!r          -5: generate spin-up single DOS (LTET=3)
!r          -6: generate spin-dn single DOS (LTET=3)
!r          -8: Joint DOS for  nonequilibrium absorption : independent Ef for electrons and holes
!r              Transition probability scaled by f(E-imrefp) * [1-f(E-imrefn)]
!r          -9: Joint DOS for  nonequilibrium emission : independent Ef for electrons and holes
!r              Transition probability scaled by [1-f(E-imrefp)] * f(E-imrefn)
!r  lpart    Resolve JDOS or Im eps into parts
!r           0  No decomposition
!r           1  Resolve eps or dos into (occ,unocc) parts
!r           2  Resolve eps or dos by k
!r           3  Both 1 and 2
!r          Add 10 to save eps as a binary file
!r  ltet     Integration mode for optics
!r           0 sampling
!r          -1 on-the-fly sampling
!r           1 original tetrahedron integration (no longer supported)
!r           2 tetrahedron, calling bzjnos
!r             In spin polarized cases, optics generated for both spins
!r           3 tetrahedron, calling tetwtq
!r             In spin polarized cases, optics generated for one spin only
!r  ltcnst   Controls how constraints are employed for restricting (occ,unocc) transitions
!r           Default occ and unocc ranges are stored in ocrng(1:2) and unrng(1:2) below
!r           Reduced occ and unocc ranges are stored in ocrng(3:4) and unrng(3:4) below
!r           0 routine ocnock adds further restrictions, storing them in ocrng(3:4) and unrng(3:4)
!r           1 No further restrictions on ocrng(3:4) and unrng(3:4) below
!r  alltrans Use Fermi function to further restrict range of (occ,unocc) transitions (set to F for loptc=8 or 9)
!r  mefac    If 1 and nonlocal sigm added to LDA hamiltonian,
!r           Add correction to velocity operator from dSig/dk
!r           If 2 and nonlocal sigm added to LDA hamiltonian,
!r           estimate correction change in velocity operator from LDA eigenvalues.
!r           See Phys. Rev. B59, 7486 Eq. 12.
!r  nchi2    number of (abc) for which to calculate xi^(abc)
!r  nlg      dimensioning parameter for radial gradient matrix elements
!r           (maximum number of l's)
!r  nfilm    number of filled states
!r  nempm    number of empty states
!r  ocut     Fermi level fuzz used for k-dependent occ,unocc cutoffs (see ocnock)
!r  ocrng    ocrng(1:2) default range of occupied bands for (occ,unocc) transitions
!r           ocrng(3:4) k-dependent possible narrowing of  ocrng(1:2) (see ocnock)
!r                      If zero, no narrowing in beyond ocrng(1:2)
!r  unrng    unrng(1:2) default range of unoccupied bands
!r           unrng(3:4) k-dependent possible narrowing of unrng(1:2) (see ocnock)
!r                      If zero, no narrowing in beyond unrng(1:2)
!r  Nmp      Methfessel-Paxton polynomial order (sampling integration)
!r  w        gaussian line broadening width or temperature (optint)
!r  window   energy window for epsilon
!r           (1) emin
!r           (2) emax
!r           (3) emin, possibly modified from initial value
!r           (4) emax, possibly modified from initial value
! ----------------------------------------------------------------
      type str_optic
      sequence

      integer    ::  axes   (3,6)
      integer    ::  cll
      integer    ::  cln
      integer    ::  cls    = 0
      integer    ::  clsite
      integer    ::  ffmt
      integer    ::  iq     (3)
      integer    ::  kmxax
      integer    ::  loptic
      integer    ::  lpart
      integer    ::  ltet
      integer    ::  alltrans
      integer    ::  mefac
      integer    ::  nchi2
      integer    ::  Nmp
      integer    ::  nlg
      integer    ::  nfilm
      integer    ::  nempm
      integer    ::  ocrng  (4)
      integer    ::  unrng  (4)
      integer    ::  inull ! ensure 8-byte word boundary

      real(8)    ::  dw     (3)
      real(8)    ::  esciss
      real(8)    ::  window (4)
      real(8)    ::  esmr
      real(8)    ::  kt
      real(8)    ::  imref  (2)
      real(8)    ::  w

      real(8), pointer :: rgrad(:,:,:,:,:) ! radial matrix elements of phi,phidot,phiz for each site
                                             ! dimensioned for now (nkaph**2,2,nlg,nsp,nbas)
      real(8), pointer :: rgrade(:,:,:,:,:)! radial matrix elements of P_kL
                                             ! dimensioned ((nkaph+1+kmxax)**2,2,nlg,nsp,nbas)
!     real(8), pointer :: optmt(:,:,:,:,:)   ! Square of optical matrix elements for all spin, k
      complex(8), pointer :: optme(:,:,:)    ! Optical matrix elementso for one spin, k

      end type str_optic

! --- str_ordn ---
!- Parameters for order-N
! ----------------------------------------------------------------
!r  element purpose
!r  efre    effective energy at which free electrons
!r          start propagating
!r  mode    order-N mode
!r  ncl     Number of clusters
!r  ndofh   leading dimension to offcH
!r  clp
!r  clssl
!r  iaxg    neighbor table for Green's function
!r  mapgv   maps GF pairs into iax table
!r  offcH   offset to table of Hamiltonian offsets
!r  rmaxg   ranges of Green's function:
!r          (1) = range for direct inversion
!r          (2) = range for iterative inversion
!r          (3) = range for free-electron GF
! ----------------------------------------------------------------
      type str_ordn
      sequence

      integer    ::   mode
      integer    ::   ncl
      integer    ::   ndofh
      integer    ::   inull ! ensure 8-byte word boundary

      real(8)    ::   efre
      real(8)    ::   rmaxg(3)

      integer, pointer :: iaxg(:)
      integer, pointer :: mapgv(:)
      integer, pointer :: clp(:)
      integer, pointer :: clssl(:)
      integer, pointer :: offcH(:)

      end type str_ordn

! --- str_rhat ---
!- Parameters for density data
! ----------------------------------------------------------------
!r  Element  Purpose
!r  ... Pointer  arrays
!r  rho1     true site density, rho1(nr,nlml,nsp)
!r  rho2     local repsn of smoothed density, rho2(nr,nlml,nsp)
!r  rhoc     core density, rhoc(nr,nsp)
!r  eqhhl    site-local energy-weighted density matrix, head-head contribution
!r  eqhkl    site-local energy-weighted density matrix, head-tail contribution
!r  eqkkl    site-local energy-weighted density matrix, tail-tail contribution
!r  qhhl     site-local density matrix, head-head contribution
!r  qhkl     site-local density matrix, head-tail contribution
!r  qkkl     site-local density matrix, tail-tail contribution
! ----------------------------------------------------------------
      type str_rhat
      sequence
      real(8), pointer    :: rho1(:,:)
      real(8), pointer    :: rho2(:,:)
      real(8), pointer    :: rhoc(:,:)
      real(8), pointer    :: eqhhl(:,:)
      real(8), pointer    :: eqhkl(:,:)
      real(8), pointer    :: eqkkl(:,:)
      real(8), pointer    :: qhhl(:,:)
      real(8), pointer    :: qhkl(:,:)
      real(8), pointer    :: qkkl(:,:)

      end type str_rhat


      type str_cpasite

      integer    ::   ib
      integer    ::   idcc
      integer    ::   ncomp

      real(8), pointer    :: angles(:)

      end type str_cpasite

! --- str_pot ---
!- Parameters for density, potential, matrix elements
! ----------------------------------------------------------------
!r  element purpose
!r  bfield  global magnetic field direction and amplitude (4)
!r  lcplxp  0 if augmentation ppi is real; 1 if ppi is complex
!r  nlma    number of augmentation channels in system
!r  nlml    (FP) total number of site density channels
!r  nrhos   (FP) total number of site density channels
!r  nrmxv   (ASA) [not used] Largest radial mesh dimension among all classes; dimensions v0
!r  pfdim   (GF) Dimensions global potential function arrays (also stored in s_ham%ldham(9))
!r          Note: not used, since dimension can be extracted from fortran shape function
!r  vconst  Constant estat potential shifts, used where
!r          the Fermi level is specified and the potential
!r          adjusts to it.
!r          vconst(1) = potential shift
!r          vconst(2) = potential shift of L end region (PGF)
!r          vconst(3) = potential shift of R end region (PGF)
!r  vmtz    ASA muffin-tin zero
!r  vmtz0   Input ASA muffin-tin zero
!r  ... dynamical arrays
!r  aamom   local ASA magnetic moments (susite;amagnc)
!r  bxc     constant magnetic field inside augm. sphere (clsprm,iomagf)
!r bxcscali scaling factor for smoothed part of XC b-field
!r  cp      Coherent potential P matrix for CPA, iprmb order (suhamz;mkcpa)
!r  del     (tbe) potential shifts.  No longer part of structure
!r  dlmwt   (clsprm; made in dlmwgts)
!r  dmatk   k-space density-matrix
!r  pf      ASA potential functions (suhamz;mkptfp)
!r  pfnc    noncollinear ASA potential functions (pgflu)
!r  dpf     energy derivative of potential functions (suhamz;mkptfp)
!r          Compare to s_site(:)%dpfr (stored either in matrix or vector form, orbital order)
!r  dddpf   3rd derivative of potential functions (suhamz;mkptfp)
!r  ddpf    2nd derivative of potential functions (suhamz;mkptfp)
!r  pfr     CPA or relativistic case potential functions : lrel=2 or SO coupling for all atoms (suhamz)
!r          s_pot%pfr holds P in lms representation, vector(s) of dimension pfdim
!r          Compare to s_site%pfr, which holds P for each site (and separate CPA species)
!r          s_pot%pfr is no longer made for CPA and SO cases; use instead s_site%pfr
!r  dpfr    fully relativistic analog of dpf (suhamz;mkptfp)
!r          Stored in vector form for the whole lattice, in downfolding order
!r          Compare to s_site(:)%dpfr (stored either in matrix or vector form, orbital order)
!r  ddpfr   fully relativistic analog of ddpf (suhamz;mkptfp)
!r  GFr     overlap integrals for the calculation of averages [Turek (6.141)] (4,4,2,2,l,imu,4)
!r  gibbs   gibbs weights for disordered local moments (clsprm)
!r  gma     gamma-alpha in potential functions (suhamz)
!r  gmar    (gamma-alpha)*P^alpha/P^gamma fully relativistic (suhamz)
!r  grrme   radial matrix elements of gradient (clsprm)
!r  hab     hamiltonian matrix elements of partial waves at each (potpus.f)
!r  sab     overlap matrix elements of partial waves at each (potpus.f)
!r  mad     Madelung matrix (allocated in supot)
!r  mxy     off-diagonal local moments for one site (clsprm)
!r  palp    ASA potential functions, alpha repsn (suhamz;mkptfp)
!r  papg    p^alpha/p^gamma fully relativistic (suhamz)
!r  pmpol   multipole integrals of phi,phidot, for nonspherical output rho, ASA (clsprm)
!r  pnu     P_nu: Methfessel's log derivative function (clsprm)
!r  pp      potential parameters (clsprm;potpar)
!r  ppn     NMTO generation potential parameters (bndasa)
!r  pprel   Dirac potential parameters (clsprm)
!r  pti     inverse potential functions (suham)
!r  qc      ASA sphere core charge (clsprm)
!r  qcorr   Lloyd corrections for single-electron energy
!r  qnu     Energy moments of the charge density (clsprm;bndasa,gfasa,pgfasa)
!r  qbyl    l-resolved charges (lmf)
!r  qnur    Energy moments of the charge density, for fully relativistic (clsprm;gfasa)
!r  qpp     ASA Multipole moments of the nonspherical charge (clsprm;asaddq,ioqpp)
!r  qt      ASA sphere total charge (clsprm)
!r  qval    total valence charge, including semicore states (bndfp)
!r  qsc     total charge from semicore states (local orbitals)
!r  rhos    spin-density matrix (clsprm)
!r  rhrmx   ASA sphere electron density at rmax (clsprm)
!r  shfac   shfac(:,i) are artificial scale factors in 1-electron hamiltonian H,
!r          that modify site i (or class i in the ASA)  (rdctrl)
!r          Normally all factors are 0, but these factors make it possible
!r          to probe how the hamiltonian varies with scalings.
!r          Data are read from file 'shfac', which can be edited by the user.
!r          if it does not exist the program will create it automatically from
!r          species data in s_spec%shfac
!r          Meaning of vector elements:
!r          shfac(1,i) = scaling factor for the XC field for site or class i
!r          shfac(2,:) = scale SO parameter lambda by 1+shfac(1,i)
!r          shfac(3,:) = scale channel bandwidth by 1+shfac(1,i)
!r  smpot   smoothed potential (supot;mkpot)
!r  smpoth  smoothed potential, Hartree part (bndfp;mkpot)
!r  smvextc smoothed external potential, stored in packed form (mkpot,ioden2)
!r  smvcnst additional constant added to smpot (mkpot
!r  smrho   smoothed density (made in mkrout; allocated in supot)
!r  socscl  scaling factors for spin-orbit parameters, by class
!r  sop     spin-orbit parameters (clsprm)
!r  thet    Angles for disordered local moments (clsprm)
!r  v0      (ASA only) spherical part of potential for all classes
!r  v0beta  mixing beta for admixing sphere potential V0 defining phi,phidot.  Usually 1
!r  vdif    difference ves(q,rmax) - vactual(rmax) (lmasa)
!r  ves     electrostatic potential at rmax, by class (used by ASA: clsprm, iors)
!r  vesrmt  electrostatic potential at rmax, by site (used by FP; made in mkpot)
!r  vintr   intra-atomic W for screening, monopole approx (clsprm,asasx)
!r  vrmax   total potential at rmax (clsprm)
!r  vshft   ASA constant potential shifts, mainly for GF codes (susite,iovshf)
! -----------------------------------------------------------------
      type str_pot
      sequence

      integer    ::   lcplxp
      integer    ::   nlma
      integer    ::   nlml
      integer    ::   nrhos
!     integer    ::   nrmxv = 0 ! Largest radial mesh dimension among all classes
!     integer    ::   pfdim ! Dimensions pf,dpf,pfr,dpfr, etc.  Not needed
!     integer    ::   inull ! ensure 8-byte word boundary

      real(8)    ::   bfield  (4)
      real(8)    ::   bxcscali
      real(8)    ::   smvcnst
      real(8)    ::   vconst  (3)
      real(8)    ::   qval
!     real(8)    ::   qsc
      real(8)    ::   v0beta
      real(8)    ::   vmtz
      real(8)    ::   vmtz0

! ... Pointer arrays
      real(8), pointer   :: aamom(:)
      real(8), pointer   :: bxc(:)
      real(8), pointer   :: dlmwt(:)
      real(8), pointer   :: GFr(:,:)
      real(8), pointer   :: gibbs(:)
      real(8), pointer   :: gma(:,:)
      real(8), pointer   :: grrme(:)
      real(8), pointer   :: hab(:,:)
      real(8), pointer   :: sab(:,:)
      real(8), pointer   :: mad(:)
      real(8), pointer   :: mxy(:)
      real(8), pointer   :: pmpol(:)
      real(8), pointer   :: pnu(:,:)
      real(8), pointer   :: pp(:,:)
      real(8), pointer   :: ppn(:)
      real(8), pointer   :: pprel(:,:)
      real(8), pointer   :: pti(:)
      real(8), pointer   :: qc(:)
      real(8), pointer   :: qcorr(:)
      real(8), pointer   :: qnu(:,:)
      real(8), pointer   :: qbyl(:,:)
      real(8), pointer   :: qnur(:,:)
      real(8), pointer   :: qpp(:)
      real(8), pointer   :: qt(:)
      real(8), pointer   :: rhos(:)
      real(8), pointer   :: rhrmx(:)
      real(8), pointer   :: shfac(:,:)
      real(8), pointer   :: socscl(:)
      real(8), pointer   :: sop(:,:)
      real(8), pointer   :: thetcl(:,:)
      real(8), pointer   :: v0(:,:)
      real(8), pointer   :: vdif(:)
      real(8), pointer   :: ves(:)
      real(8), pointer   :: vesrmt(:)
      real(8), pointer   :: vintr(:)
      real(8), pointer   :: vrmax(:,:)
      real(8), pointer   :: vshft(:)

      complex(8), pointer:: cp(:)
      complex(8), pointer:: dddpf(:,:)
      complex(8), pointer:: ddpf(:,:)
      complex(8), pointer:: ddpfr(:,:)
      complex(8), pointer:: dpf(:,:)
      complex(8), pointer:: dpfr(:,:)
      complex(8), pointer:: gmar(:,:)
      complex(8), pointer:: palp(:,:)
      complex(8), pointer:: papg(:,:)
      complex(8), pointer:: pf(:,:)
      complex(8), pointer:: pfm(:,:)   ! P in matrix form
      complex(8), pointer:: pfma(:,:)  ! P^alp in matrix form
      complex(8), pointer:: pfnc(:,:)  ! Seems to be used only by lmpg
      complex(8), pointer:: pfr(:,:)

      complex(8), pointer:: dmatk(:,:)

      complex(8), pointer:: smpot(:,:)
      complex(8), pointer:: smrho(:,:)
      complex(8), pointer:: smrout(:,:)
      complex(8), pointer:: smpoth(:,:)
      complex(8), pointer:: smvextc(:)
      complex(8), pointer:: smcv0(:,:)

      type(str_rhat),pointer ::  rhat(:),rnew(:)

      end type str_pot

! --- str_spec ---
!- Parameters for species data
! ----------------------------------------------------------------
!r  element  purpose
!r  a        a for mesh
!r  alpha    screening alphas
!r  beta     mixing beta
!r  chfa     coefficients to fit of free-atom density tails
!r  colxbs   color of atom (for graphics)
!r  coreh    core hole channel (char: eg '1s')
!r  coreq    coreq(1) = charge in core hole channel
!r           coreq(2) = moment in core hole channel (nsp=2 only)
!r  ctail    coefficients to fit of free-atom core tail by unsm. Hankel
!r  dv       constant shift in potential
!r  eh3      sm Hankel energy for high local orbitals
!r  ehvl     l-dependent sm Hankel energies for val-lap basis
!r  enu      linearization energies
!r  eref     reference energy
!r  etail    energy to fit of free-atom core tail
!r  etf      Thomas-Fermi screening length
!r  exi      Hankel energies for fit to c.d.
!r           For free atoms, fit to free-atom density tails.
!r  group    group
!r  grp2     (ASA) groups together classes that are not equivalent by symmetry.
!r           lmgf-CPA makes special use of array grp2 to make DLM classes equivalent
!r  hcr      kinky wave hard sphere radii, atomic units
!r  hsfitp   When generating c01 defining generalized e0 functions : e0(rsm,eh;r) = hsml + c01(0)*g0l + c01(1)*g1l
!r           c01 are determined by minimizing an objective function f, namely the average difference
!r            in K.E. between e0 and the bare Hankel h(eh;r) over an interval.  Specifically
!r              f = integral_rmin^infinity [nabla(e0(r)-hunsm(r))]^2*wt(r)]
!r           hsfitp(1) determines rmin through rmin = hsfitp(1)*rmt.
!r           In one minimization scheme (see Mode 4, cgmin.f), the objective function is instead
!r              f + hsfitp(2)*t^2 where t = (nabla (h+c0*g0+c1*g1) + ttarg*(h+c0*g0+c1*g1))|r=rMT
!r           Thus hsfitp(2) determines ratio of local to asymptotic error in K.E. to minimize.
!r           See e0parm and cgmin for further description of hsfitp (called fitp in cgmin, which see)
!r  idmod    idmod(l) controls how linearization energy is
!r           determined for an orbital of angular momentum l
!r           0 float to center of gravity of occupied part of band
!r           1 freeze to spec'd logarithmic derivative
!r           2 freeze to spec'd linearization energy
!r           3 float to natural center of band (C parameter)
!r           4 similar to idmod=1, but Pup, Pdn are adjusted to eliminate spin splitting in enu
!r  idu      (LDA+U, or DMFT) identifies l-channels with Hubbard U
!r           1s digit (LDA+U)
!r           0 no U
!r           1 D.C. around MF limit; see Petukhov, PRB 67, 153106 (2003)
!r           2 D.C. around FLL limit; see PRB 52, R5467 (1995)
!r           10s digit
!r           0 put U on phi only
!r           1 put U on phi,phidot,LO
!r           100s digit (LDA+DMFT)
!r           0 no U
!r           1 D.C. around MF limit; see Petukhov, PRB 67, 153106 (2003)
!r           2 D.C. around FLL limit; see PRB 52, R5467 (1995)
!r           3 D.C. with Haule's nominal charge constraint
!r  idm      (DMFT -- no longer used)
!r           One number for each nonzero l in 100s digit idu, defining m-specific treatment of correlation
!r           One number consists of (2l+1) digits corresponding to markers for each m
!r           digit=0 or 1: normal treatment of correlation
!r           digit=2 : do not treat as correlated
!r           digit=3 : treat at Hartree Fock level
!r  idxdn    idxdn
!r  iq1      (asa --fit) switches to freeze ASA C parameters
!r           (tbe --fit) switches to freeze TB site energies
!r  iq2      (asa --fit) switches to freeze ASA Delta parameters
!r           (tbe --fit) switches to freeze TB spin-orbit params
!r  jh       LDA+U J parameters for each l-channel
!r  kmxh     k cutoffs for head augmentation expansion
!r  kmxt     k cutoffs for tail augmentation expansion
!r  kmxv     k-cutoff for 1-center projection of free-atom rho
!r  lfoca    switch specifying treatment of core density
!r  lmxa     l cutoff for augmentation expansion
!r  lmxb     highest l in basis
!r  lmxf     fitting l
!r  lmxl     cutoff for local density
!r  lmxpb    l-cutoff for product basis
!r  lxi      l-cutoffs for interstital product basis
!r  mass     atomic mass
!r  mxcst    species-specific constraints of various types:
!r           bit function
!r           1   suppresses c.d. mixing for this species (ASA)
!r           2   Exclude this species when auto-resizing sphere radii
!r           4   Freeze augmentation w.f. for this species (FP)
!r  name     species name
!r  naug     number of k,L for augmentation
!r  ngcut    orbital-dependent G cutoffs (for nfp basis)
!r           ngcut(l+1,ik) = G cutoff for given l, kappa
!r  norb     number of orbitals for this species
!r  norp     number of parameters needed to define orbital
!r  nr       nr for mesh
!r  ntorb    number of orbital types in basis
!r  nxi      Number of energies in fit of free-atom density tails
!r  ncpa     (CPA) Number of chemical components
!r  ncomp    (CPA) (Number of chemical components)*(Number of DLM angles)
!r  nthet    (CPA) number of DLM angles for this species
!r  nang     (CPA) number of solid DLM angles for this species (relativistic)
!r  iscpa    (chemical CPA) iscpa(1:ncpa) : indices to parent species
!r  orbp     (lmf) parameters defining orbital shape
!r  rhoc     pointer to core density
!r  p        log derivative parameter
!r  pb1      product basis for linear response
!r  pb2      second product basis for linear response
!r  permt    Potential energy at hard core radius, resolved by l and spin
!r  pz       log derivative parameter, semicore states
!r  q        starting q's (charges)
!r  qc       core charge
!r  qpol     (tbe) 10 polarisability parameters
!r  radxbs   radius of atom (for graphics)
!r  rcfa     renormalization radius of free atom density, and width
!r  rcut     range over which true, smoothed Hankels differ---may
!r           be used in the calculation of XC potential
!r  rfoca    smoothing radius for frozen core overlap approx
!r  rg       rsm for gaussians to fix multipole moments
!r  rham     range of hamiltonian
!r  rint     range of interstitial product basis
!r  rmt      augmentation radius
!r  rs3      Lower bound to rsm for local orbital
!r  rsminl   Envelope smoothing radii constrained not to fall below rsminl
!r  rsma     rsm for augmentation expansion
!r  rsmfa    rsm to fit free atom density
!r  rsmv     rsm for one-center projection of free-atom density
!r           In TCF, smoothing radius to calculate exchange
!r           (artificially smoothed density)
!r  shfac    Artificial factors modifying in 1-electron hamiltonian H.
!r           Normally all factors are 0, but these factors make it possible
!r           to probe how the hamiltonian varies with scalings.
!r           s_spec%shfac are read through the input file.
!r           Data actually used by the pgram are read from file 'shfac' and kept
!r           s_pot%shfac. If 'shfac' does not exist and at least one of the species
!r           elements is nonzero 'shfac' is created and filled with species data.
!r           Otherwise, species data is not used.
!r           See str_pot for meaning of vector elements.
!r  stc      core kinetic energy
!r  stni     (tbe) Stoner I
!r  stnm     Stoner mmax
!r  uh       Hubbard U
!r  vmtz     Asymptotic potential for fitting functions at rmt
!r  vso      (tbe) spin-orbit parameters
!r  z        atomic number
! ----------------------------------------------------------------
      integer, private, parameter :: nkap0=4

      type str_spec

      sequence

      integer    ::   ntorb
      integer    ::   naug
      integer    ::   norb
      integer    ::   lmxa
      integer    ::   kmxh
      integer    ::   lmxl
      integer    ::   kmxt
      integer    ::   norp
      integer    ::   nr
      integer    ::   lfoca
      integer    ::   lmxb
      integer    ::   lmxf
      integer    ::   mxcst
      integer    ::   group
      integer    ::   grp2
      integer    ::   nxi
      integer    ::   ncpa
      integer    ::   nthet
      integer    ::   nang(2)
      integer    ::   ncomp
      integer    ::   iscpa(10)
      integer    ::   nbeff
      integer    ::   lmxpb
      integer    ::   lxi
      integer    ::   kmxv
      integer    ::   idmod(10)
      integer    ::   idxdn(10,nkap0)
      integer    ::   ngcut(10,nkap0)
      integer    ::   idu(4)
!     integer    ::   idm(4)
      integer    ::   iq1(10)
      integer    ::   iq2(10)
      integer    ::   ivso(10)
      integer    ::   inull ! ensure 8-byte word boundary

      character*8::   coreh
      character*8::   name
      character*8::   pb1
      character*8::   pb2

      real(8)    ::   z
      real(8)    ::   mass
      real(8)    ::   rmt
      real(8)    ::   rsmfa
      real(8)    ::   rsma
      real(8)    ::   rg
      real(8)    ::   rsmv
      real(8)    ::   coreq(2)
      real(8)    ::   a
      real(8)    ::   eref
      real(8)    ::   etf
      real(8)    ::   beta
      real(8)    ::   ctail
      real(8)    ::   etail
      real(8)    ::   stc
      real(8)    ::   stni
      real(8)    ::   stnm
      real(8)    ::   rham
      real(8)    ::   rfoca
      real(8)    ::   dv
      real(8)    ::   xcpa(10)
      real(8)    ::   beff(10)
      real(8)    ::   qc
      real(8)    ::   colxbs(3)
      real(8)    ::   radxbs
      real(8)    ::   rcut
      real(8)    ::   rint
      real(8)    ::   eh3
      real(8)    ::   rs3
      real(8)    ::   vmtz
      real(8)    ::   rcfa(2)
      real(8)    ::   p(10,2)
      real(8)    ::   q(10,2)
      real(8)    ::   alpha(10)
      real(8)    ::   hcr(10)
      real(8)    ::   exi(10)
      real(8)    ::   hsfitp(2)
      real(8)    ::   qpol(10)
      real(8)    ::   chfa(10,2)
      real(8)    ::   orbp(10,2,nkap0)
      real(8)    ::   enu(10)
      real(8)    ::   rsminl(10)
      real(8)    ::   permt(10,2)
      real(8)    ::   pz(10,2)
      real(8)    ::   uh(4)
      real(8)    ::   jh(4)
      real(8)    ::   ehvl(10)
      real(8)    ::   shfac(3)
      real(8)    ::   vso(4)

      real(8), pointer :: rhoc(:)

      end type str_spec

! --- str_site ---
!- Parameters for site data
! ----------------------------------------------------------------
!r  Element  Purpose
!r  amag     vector magnetization in sphere,
!r           defining quantization axis in noncollinear case
!r           amag(0) = size of magnetization
!r           amag(1:3) = magnetization as a vector
!r  bfield   site magnetic field (not used)
!r  cg0      coefficients to gaussian part of sm h basis
!r  clabel   string labelling class
!r  class    class index
!r  cli      (RGF) cluster index
!r  delta    (tbe) vector of on-site energy shifts
!r  dlmcl    class index
!r           For a normal site, site ib belongs to class dlmcl = ipc(ib) (mksym.f)
!r           For a CPA site, dlmcl is the first CPA class associated with ib
!r           There are ncomp classes associated with this site
!r  dpole    point charge dipole (molecules)
!r  eula     Euler angles (noncollinear magnetism)
!r  force    force (dynamics)
!r  mpole    point charge monopole (molecules)
!r  ncpa     (CPA) number of chemical CPA components
!r  ndelta   size of delta
!r  norb     Number of lower+intermediate+higher orbitals belonging to this site (setnorb.f) ...
!r  ncomp    (CPA) total number of components.  ncomp=1 if no CPA
!r  offh     hamiltonian offset for this site
!r  pl       (PGF) principal layer index
!r  plv      (PGF) principal layer potential index
!r  pnu      log derivative parameter
!r  pos      true coordinates of atom
!r  pos0     coordinates of atom before distortion
!r  pz       log derivative parameter, semicore states
!r  relax    (dynamics) flags which coordinates to relax
!r  rmaxs    Site-dependent radial cutoff for strux, in a.u.
!r  rsm      smoothing radius in sm h basis
!r  saxis    (FP) spin quantization axis
!r  sid      (RGF) Site id, used to merge coincident clusters
!r  spec     species index
!r  vel      velocity
!r  vshft    constant potential shift for this site
!r  ... Pointer  arrays
!r  bxc      exchange-correlation field
!r  cpawt    (CPA) CPA weights for each component
!r  domg     (CPA) Omega(z+dz) to differentiate Omega (see omg)
!r  gcorr    (CPA) Lloyd correction for single-particle energy
!r  gii      Unscaled on-site g for the given site.
!r                 gii is in lms representation, even in Dirac case
!r  sfvrtx   (CPA) site vertex for spin flips
!r  vtxrel   (CPA) site vertex for spin flips, non-collinear case
!r  j0       (CPA) sum of exchange interaction at a site
!r  tau      (CPA) spin lifetime at a site
!r  omg      (CPA) coherent interactor Omega (dlminit; gfomg)
!r  omgn     (CPA) output coherent interactor Omega
!r  alr      (CPA) Matrices alr,amr,amrt are used in CPA spectral function calculations
!r  amr      (CPA) See Eq. (4.26) in Turek's book
!r  amrt     (CPA) See specfun.f
!r  dmat     (CPA) Site Density-of-states (gfidos2)
!r  pfr      Potential functions for one site. Stored in orbital order as:
!r             nonrelativistic case : vector(s) of size norb (mkptfp)
!r             relativistic case    : norb*2 x norb*2 matrix in (mkfrfp or mkptso)
!r             fully relativistic pfr stored in mu,lambda,lambda repsn (mkfrfp)
!r           Compare to s_pot%pfr.  It holds P in lms representation, vector(s) or size norb
!r  dpfr     Potential functions pdot or sqrt(pdot), pdot = energy derivative of pfr
!r             dpfr is stored in the same manner as pfr , additionally:
!r             in vector form, dpfr is stored as pdot; in matrix form stored as sqrt(pdot)
!r  ddpfr    Potential functions -1/2*(P_dot_dot)*(P_dot)^(-1) with SO coupling (mkptso)
!r             ddpfr is stored in the same manner as pfr
!r  pfra     Potential functions for one site, alpha repsn
!r             Same structure as pfr
!r  gc       (CPA) conditionally-averaged onsite scaled G for each component (mkgint)
!r                 In the Dirac case gc is stored in the lambda-mu repsn
!r  gcu      (CPA) Unscaled version of gc (used when lrel=2 or SOC) (mkgint)
!r                 For sites NOT CPA sites, gcu and gii are the same.
!r                 gcu is in lms representation, even in Dirac case.
!r  gcorr    (CPA) Lloyd correction to be used for single-particle energy (mkgint)
!r  pdos     partial DOS
!r  rho1     One-center expansion of output density inside augmentation sphere
!r  rho2     One-center expansion of smoothed output density
!r  rhoc     Core density
!r  rho1x    auxiliary pointer for rho1 used e.g. for output density
!r  rho2x    auxiliary pointer for rho2 used e.g. for output density
!r  rhocx    auxiliary pointer for rhoc used e.g. for output density
!r  qhhl     Site-density matrix, head-head part (mkrout)
!r  qhkl     Site-density matrix, head-tail part (mkrout)
!r  qkkl     Site-density matrix, tail-tail part (mkrout)
!r  eqhhl    Energy-weighted site-density matrix, head-head part (mkrout)
!r  eqhkl    Energy-weighted site-density matrix, head-tail part (mkrout)
!r  eqkkl    Energy-weighted site-density matrix, tail-tail part (mkrout)
!r  sighh    Overlaps of (radial wave function)*(radial wave function)
!r  sighk    Overlaps of (radial wave function)*(PkL)
!r  sigkk    Overlaps of (PkL)*(PkL)
!r  tauhh    Kinetic energy of (radial wave function)*(radial wave function)
!r  tauhk    Kinetic energy of (radial wave function)*(PkL)
!r  taukk    Kinetic energy of (PkL)*(PkL)
!r  pihh     Potential matrix element of (radial wave function)*(radial wave function)
!r  pihk     Potential matrix element of (radial wave function)*(PkL)
!r  pikk     Potential matrix element of (PkL)*(PkL)
!r  sohh     L.S perturbation (radial wave function)*(radial wave function)
!r  sohk     L.S perturbation (radial wave function)*(PkL)
!r  sokk     L.S perturbation (PkL)*(PkL)
!r  sighhx   auxiliary pointer for sighh used e.g. for hamiltonian without exchange
!r  sighkx   auxiliary pointer for sighk used e.g. for hamiltonian without exchange
!r  sigkkx   auxiliary pointer for sigkk used e.g. for hamiltonian without exchange
!r  tauhhx   auxiliary pointer for tauhh used e.g. for hamiltonian without exchange
!r  tauhkx   auxiliary pointer for tauhk used e.g. for hamiltonian without exchange
!r  taukkx   auxiliary pointer for taukk used e.g. for hamiltonian without exchange
!r  pihhx    auxiliary pointer for pihh  used e.g. for hamiltonian without exchange
!r  pihkx    auxiliary pointer for pihk  used e.g. for hamiltonian without exchange
!r  pikkx    auxiliary pointer for pikk  used e.g. for hamiltonian without exchange
!r  thet     DLM angles
!r  v0       Spherical potential that defines wave functions
!r  v1       Sperical part of MT potential
! ----------------------------------------------------------------
      type str_site
      sequence

      integer    ::  class
      integer    ::  cli
      integer    ::  dlmcl
      integer    ::  ncpa
      integer    ::  ndelta
      integer    ::  norb
      integer    ::  ncomp
      integer    ::  offH
      integer    ::  pl
      integer    ::  plv
      integer    ::  relax  (3)
      integer    ::  sid
      integer    ::  spec
      integer    ::  inull ! ensure 8-byte word boundary

      character(len=8) ::  clabel

      real(8)    ::  amag   (0:3)
      real(8)    ::  bfield (3)
      real(8)    ::  cg0    (10)
      real(8)    ::  delta  (6)
      real(8)    ::  dpole  (3)
      real(8)    ::  eula   (3)
      real(8)    ::  force  (3)
      real(8)    ::  mpole
      real(8)    ::  pnu    (10,2)
      real(8)    ::  qnu    (3,10,2)
      real(8)    ::  pos    (3)
      real(8)    ::  pos0   (3)
      real(8)    ::  pz     (10,2)
      real(8)    ::  rmaxs
      real(8)    ::  rsm    (10)
      real(8)    ::  saxis  (3)
      real(8)    ::  vel    (3)
      real(8)    ::  vshft

! ... Pointer arrays
!     for CPA (ASA)
      real(8), pointer    :: bxc(:)
      real(8), pointer    :: j0(:,:)
      real(8), pointer    :: tau(:,:)
      real(8), pointer    :: cpawt(:)
      real(8), pointer    :: thet(:,:)
      complex(8), pointer :: domg(:,:)
      complex(8), pointer :: dmat(:,:)
      complex(8), pointer :: pfr(:,:)
      complex(8), pointer :: dpfr(:,:)
      complex(8), pointer :: ddpfr(:,:)
      complex(8), pointer :: pfra(:,:)
      complex(8), pointer :: gcu(:,:)
      complex(8), pointer :: gc(:,:)
      complex(8), pointer :: gcorr(:,:)
      complex(8), pointer :: gii(:,:)
      complex(8), pointer :: omg(:,:)
      complex(8), pointer :: omgn(:,:)
      complex(8), pointer :: sfvrtx(:,:)
      complex(8), pointer :: vtxrel(:,:)
      complex(8), pointer :: alr(:,:,:,:)
      complex(8), pointer :: amr(:,:,:,:)
      complex(8), pointer :: amrt(:,:,:,:)

!     Spherical part of potential
      real(8), pointer    :: v0(:,:)
      real(8), pointer    :: v1(:,:)

!     density-like objects
      real(8), pointer    :: pdos(:,:)
      real(8), pointer    :: rho1(:,:)
      real(8), pointer    :: rho2(:,:)
      real(8), pointer    :: rhoc(:,:)
      real(8), pointer    :: rho1x(:,:)
      real(8), pointer    :: rho2x(:,:)
      real(8), pointer    :: rhocx(:,:)

      real(8), pointer    :: qhhl(:,:)
      real(8), pointer    :: qhkl(:,:)
      real(8), pointer    :: qkkl(:,:)

      real(8), pointer    :: eqhhl(:,:)
      real(8), pointer    :: eqhkl(:,:)
      real(8), pointer    :: eqkkl(:,:)

!     Matrix elements
      real(8), pointer    :: sighh(:,:)
      real(8), pointer    :: sighk(:,:)
      real(8), pointer    :: sigkk(:,:)
      real(8), pointer    :: tauhh(:,:)
      real(8), pointer    :: tauhk(:,:)
      real(8), pointer    :: taukk(:,:)
      real(8), pointer    :: pihh(:,:)
      real(8), pointer    :: pihk(:,:)
      real(8), pointer    :: pikk(:,:)
      real(8), pointer    :: sohh(:,:)
      real(8), pointer    :: sohk(:,:)
      real(8), pointer    :: sokk(:,:)


      real(8), pointer    :: sighhx(:,:)
      real(8), pointer    :: sighkx(:,:)
      real(8), pointer    :: sigkkx(:,:)
      real(8), pointer    :: tauhhx(:,:)
      real(8), pointer    :: tauhkx(:,:)
      real(8), pointer    :: taukkx(:,:)
      real(8), pointer    :: pihhx(:,:)
      real(8), pointer    :: pihkx(:,:)
      real(8), pointer    :: pikkx(:,:)

      end type str_site

! --- str_str0 ---
!- Parameters for structure constants (subset of str containing allocatable arrays)
! ----------------------------------------------------------------
!r  element purpose
!r  a      :screening parameters, e.g. alpha in 2nd gen LMTO
!r  i      :neighbor table iax; see pairc.f
!r  n      :number of pairs in iax table preceding ib; see ntab in pairc.f
!r  s      :structure constants
! ----------------------------------------------------------------
      type str_str0

!     Neighbor tables
      integer, pointer :: i(:)
      integer, pointer :: n(:)
      integer, pointer :: ng(:)

!     For structure constants
      real(8), pointer :: a(:)
      real(8), pointer :: s(:)

      end type str_str0

! --- str_str ---
!- Parameters for structure constants, lmto, tb and TCF
! ----------------------------------------------------------------
!r  element purpose
!r  adec    (mol) log spacing in distance for TCF fit (input TCF_ADEC)
!r  amode   mode for calculating screening alpha's (input STR_ENV_MODE)
!r          1s digit
!r           0: alpha=alp0, (which is input)
!r           1: alpha=gamma, calculated from pot pars.
!r              NO LONGER SUPPORTED
!r           2: alpha calculated for hard spheres.
!r              NO LONGER SUPPORTED; see mktral
!r          10s digit 0: adot=Andersen's localized
!r           0: adot=Andersen's localized
!r           1: adot so Kanpur notes Eq 3.87 is zero
!r           2: adot calculated for hard spheres.
!r              NO LONGER SUPPORTED; see mktral
!r  delrx    Range of screened function beyond last site in cluster (input STR_DELRX)
!r  drwats   padding beyond cluster for Watson sphere (input STR_DRWATS)
!r  iinv     (mol) Parameters for iterative inversion of tb-strx (input STR_IINV)
!r           arg 1: number of sites for cutoff
!r           arg 2: number of iterations
!r           arg 3: allowed tolerance
!r           arg 4-5: not used now.
!r  ivl      (lmstr) Type of val-lap functions (input STR_ENV_VLFUN)
!r  kaps     Kinetic energies at which strux evaluated (NMTO) (input STR_ENV_EL)
!r  lequiv   look for equivalent connecting vectors
!r  lmaxw    lmax for Watson sphere
!r  lmem     memory handling
!r  loka     conventions for Hankel and Bessel functions
!r           0 Methfessel's conventions
!r           1 Andersen's 2nd gen LMTO conventions
!r  lshow    show strux
!r  mnnbr    minimum allowed number of neighbors in cluster
!r           (alternative to rmax)
!r  mxnbr    maximum number of neighbors in cluster (pairc)
!r  nalf     (mol) polynomial order for d-dependence of TCF fit coffs (input TCF_NALF)
!r  nbisi    (mol) mesh parameters for TCF (input TCF_NBISI)
!r  ncupl    (mol) reduced polynomial order in TCF (input TCF_NCUPL)
!r  ndust    (mol) no. iterations to improve on TCF (input TCF_NDUST)
!r  nkaps    number of energies for envelope functions (input STR_ENV_NEL)
!r  nttab    number of pairs in this table (pairc)
!r  nitab    number of pairs stored in strux (asastr)
!r           Typically same as nttab, but equivalent pairs may be mapped.
!r  nds      leading dimension of structure constant matrix s, sdot (asastr,rdstrx)
!r  alph     screening alpha's (asastr,rdstrx)
!r  adot     screening alpha's, e. derivative (asastr,rdstrx)
!r  iax      neighbor table (pairc,tbzint,lmce,rdstrx)
!r  npr      number-of-pairs table (pairc,rdstrx)
!r  s        strux (asastr,rdstrx)
!r  sdot     strux energy derivative (asastr,rdstrx)
!r  rfit     ?? (not used) default hard core radius
!r  rmax     range of strux (clsprm)
!r           tbe: range of hamiltonian (tbzint)
!r  rmaxg    range of val-lap strux (input STR_RVL/R)
!r  skmsh    parameters for 3rd gen LMTO energy interpolation
!r           not used?
!r  tolg     Tolerance in l=0 gaussians, which determines their range
!r  wztcf    (mol) weighting of z-points in TCF (input TCF_WZTCF)
!r  envtype  Indicates how new envelope functions e0 are defined (e0parm)
!r           The envelope function consists of a "base" which are constructed
!r           to be as "close as possible" to bare Hankels in the interstitial
!r           Each "head" is further modified to make the K.E. match a specified value (t0parm)
!r           1s digit
!r           0,1  Not used now
!r           2    e0(rsm,eh;r) = hsml + c01(0)*g0l + c01(1)*g1l
!r           10s digit
!r           2    adjust rsmh so K.E. deviation matches tolke (see e0parm.f)
!r           100s digit
!r           0    eh for each kappa is independent
!r           1    eh for kappa=1,2 are the same
!r  rsmh     Smoothing radii for e0 functions (e0parm)
!r  c01g     Coefficients to gaussian contribution to the "base" part of e0 functions (e0parm.f)
!r           There is a family of coefficients for each species.
!r  c01site  Coffs to the modification of the the "base" part e0 to fit the head (t0parm.f)
!r           There is a family of coefficients for each site.
! ----------------------------------------------------------------
      type str_str
      sequence

      integer    ::  amode
      integer    ::  ivl
      integer    ::  lequiv
      integer    ::  lmaxw
      integer    ::  lmem
      integer    ::  loka
      integer    ::  lshow
      integer    ::  mnnbr
      integer    ::  mxnbr
      integer    ::  nalf
      integer    ::  nbisi  (3)
      integer    ::  ncupl
      integer    ::  ndust
      integer    ::  nitab
      integer    ::  nkaps
      integer    ::  nttab
      integer    ::  nds
      integer    ::  envtype      ! Flags affecting how generalized envelope functions are defined
C     integer    ::  inull ! ensure 8-byte word boundary

      real(8)    ::  adec
      real(8)    ::  delrx
      real(8)    ::  drwats
      real(8)    ::  iinv   (3)
      real(8)    ::  kaps   (6)
      real(8)    ::  rfit
      real(8)    ::  rmax
      real(8)    ::  rmaxg
      real(8)    ::  skmsh  (6)
      real(8)    ::  tolg
      real(8)    ::  wztcf

!     Neighbor tables
!#ifndef LINUXF
      integer, pointer :: iax(:) => null()
      integer, pointer :: npr(:) => null()
!#else!
!      integer, pointer :: iax(:)
!      integer, pointer :: npr(:)
!#endif

!     For structure constants
!#ifndef LINUXF
      real(8), pointer :: adot(:) => null()
      real(8), pointer :: alph(:) => null()
      real(8), pointer :: s(:)    => null()
      real(8), pointer :: sdot(:) => null()
!#else!
!      real(8), pointer :: adot(:)
!      real(8), pointer :: alph(:)
!      real(8), pointer :: s(:)
!      real(8), pointer :: sdot(:)
!#endif

!     For screened structure constants with generalized envelope functions
      real(8), pointer :: rsmh(:,:)       ! Smoothing radii for e0 functions; see e0parm.f
      real(8), pointer :: c01g(:,:,:)     ! Coffs to gaussian part of e0 functions; see e0parm.f
      real(8), pointer :: c01site(:,:,:)  ! Modification of c01 to fit head K.E.; see t0parm.f

      end type str_str

! --- str_tb ---
!- Parameters for Levenberg-Marquardt fitting of bands, DOS
! ----------------------------------------------------------------
!r  element purpose
!r  alam    lambda variable, Levenburg-Marquardt (Numerical Recipes)
!r  alsc    scale factor for Levenburg-Marquardt lambda
!r  ebfit   energy window for fitting eigenvalues or DOS
!r          lower and upper limit of eigenvalues (DOS) to be fit
!r  fmode   mode for fitting parameters
!r          tbe:
!r          0: fit parameters used in tbham
!r          1: not implemented
!r          lm (see lmfitmr1)
!r          0: No fitting
!r          1: Fit ASA C and Delta
!r  nbfit   indices to first and last bands to be included
!r          in fit for all k-points
!r  nbfitf  index to first reference band to be fit
!r  ndos    number of energy points in ebfit (fitting DOS)
!r  rmfit   real space range to fit TB parameters (fmode=1)
!r  shft    shift bands by a constant
!r          0 no shift
!r          1 Find Ef and align with Ef in given bands
!r          2 Use Ef as specified in ctrl file, and
!r            align with Ef in given bands
!r          Future: 2nd argument specifies particular k
!r  wg      wg(1) gaussian broadening of fit DOS
!r          wg(2) gaussian broadening of reference DOS
!r  wt      weights for fitting bands or DOS to reference
!r          (1) defines functional form of weights
!r           0: unit weights
!r           1: Fermi function,
!r           2: Fermi-like gaussian cutoff
!r           3: Exponential weights
!r           4: Gaussian weights
!r          (2) energy shift (e.g. center for Fermi function)
!r          (3) energy scale (e.g. width of Fermi function)
! ----------------------------------------------------------------
      type str_tb
      sequence

      integer    ::  fmode
      integer    ::  nbfit  (2)
      integer    ::  nbfitf
      integer    ::  ndos
      integer    ::  shft   (2)
      integer    ::  inull ! ensure 8-byte word boundary

      real(8)    ::  alam
      real(8)    ::  alsc
      real(8)    ::  ebfit  (2)
      real(8)    ::  rmfit  (2)
      real(8)    ::  wg     (2)
      real(8)    ::  wt     (3)

      end type str_tb

! --- str_strn ---
!- Holds global strings
! ----------------------------------------------------------------
      type str_strn
!     character(len=1),allocatable :: amix(:)
!     character(len=1),allocatable :: blockh(:)
!     character(len=1),allocatable :: gemb(:)
!     character(len=1),allocatable :: gfopt(:)
!     character(len=1),allocatable :: jobid(:)
!     character(len=1),allocatable :: map(:)
!     character(len=1),allocatable :: mix(:)
!     character(len=1),allocatable :: mmham(:)
!     character(len=1),allocatable :: sxopt(:)
!     character(len=1),allocatable :: symg(:)
!     sequence

!#ifndef LINUXF
      character(len=1),pointer :: strn(:) => null()
!#else!
!      character(len=1),pointer :: strn(:)
!#endif

      end type str_strn

! --- Structure str_atparms contains integrals inside augmentation spheres ---
!     These terms are needed for total energy.  Used by mkpot, generated in locpot.
!     True local density n1 = n1_val + n1_core (valence + core parts)
!     Smooth local density n2, no core, not compensated by gaussians
!     Compensating gaussians gval making n2~-n2.  Multipole moments for n2~ and n1 are equal
!     gcor  : gaussians making up smooth core density
!     gnuc  : gaussians making up smooth nuclear density
!     ves1~ : True   estat potential = ves[n1_val + n1_core] - 2*Z/r
!     ves2~ : Smooth estat potential = ves[n2+gval+gcor+gnuc]
!     gpotb : integrals of gaussians (radius rg) times ves2~

      type str_atparms
      sequence

!     Charges and magnetic moments
      real(8)    :: qv1       = 0d0  ! true valence charge
      real(8)    :: qv2       = 0d0  ! sm   valence charge
      real(8)    :: qcor1     = 0d0  ! true core charge
!     real(8)    :: qcor2     = 0d0  ! sm   core charge
      real(8)    :: qloc      = 0d0  ! total valence charge rho1-rho2 in sphere
      real(8)    :: qlocc     = 0d0  ! total core charge in sphere
      real(8)    :: atrue(0:3)= 0d0  ! Magnetic moment of true density, and avg direction
      real(8)    :: aloc(0:3) = 0d0  ! total valence magnetic moment in sphere+core moment for sites with core hole
      real(8)    :: alocc     = 0d0  ! core magnetic moment

!     Integrals between density and electrostatic potential
      real(8)    :: vales1    = 0d0  ! int n1_val ves1,  ves1~ = ves[n1] - 2*Z/r
      real(8)    :: vales2    = 0d0  ! int n2_val ves2,  ves2~ = ves[n2+gval+gcor+gnuc]
      real(8)    :: vvesat    = 0d0  ! vales1 - vales2.  Total electrostatic energy has a term vvesat/2
      real(8)    :: sgpotb    = 0d0  ! integrals of compensating gaussians times ves2~
!     real(8)    ::  vcpn1    = 0d0  ! int n1_core * ves1~ - Z ves1~(r=0)
!     real(8)    ::  vcpn2    = 0d0  ! int (n2_core + sm nuc) * ves2~
      real(8)    :: cpnves    = 0d0  ! vcpn1 - vcpn2.  Total electrostatic energy has a term cpnves/2

!     Integrals between density and exchange-correlation potential
      real(8)    :: rep1(2)   = 0d0  ! int n1 * exc[n1], n1 = true density
      real(8)    :: rep2(2)   = 0d0  ! int n2 * exc[n2], n2 = sm density (may be pert corr)
      real(8)    :: rhoexc(2) = 0d0  ! rep1 - rep2
      real(8)    :: rhoex(2)  = 0d0  ! exchange-only part of rhoexc
      real(8)    :: rhoec(2)  = 0d0  ! correlation-only part of rhoexc
      real(8)    :: rmu1(2)   = 0d0  ! int n1 * vxc[n1], n1 = true density
      real(8)    :: rmu2(2)   = 0d0  ! int n2 * vxc[n2], n2 = sm density (may be pert corr)
      real(8)    :: rhovxc(2) = 0d0  ! rmu1 - rmu2
      real(8)    :: focexc(2) = 0d0  ! (lfoc=2) int gcor and xc energy density of n2
      real(8)    :: focvxc(2) = 0d0  ! (lfoc=2) int gcor and (dvxc/dn2 * n2)
      real(8)    :: focex(2)  = 0d0  ! (lfoc=2) Exchange part of focexc
      real(8)    :: focec(2)  = 0d0  ! (lfoc=2) Correlation part of focexc

!     Integrals between density and total potential
      real(8)    :: rhovv1    = 0d0  ! int n1_val * v1~   where v1~ = ves[n1] + vxc[n1]
      real(8)    :: rhovv2    = 0d0  ! int n2_val * v2~ + sum_L qmom_L gpotb_L, v2~ = ves[n2] + vxc[n2]
                                     ! int [(rho2*v2) + (n0~-n0)*gpotb] + focvxc
      real(8)    :: rhcvef1   = 0d0  ! int rhoc*(v1-2Z/r) for sites w/ lfoca=0

!     Integrals between valence-only density and exchange-correlation potential
!     Quantities are generated only when the appropriate switch is set
      real(8)    :: rvepsv    = 0d0  ! like rhoeps, but for valence density only
      real(8)    :: rvexv     = 0d0  ! exchange    part of rvepsv
      real(8)    :: rvecv     = 0d0  ! correlation part of rvepsv
      real(8)    :: rvvxcv    = 0d0  ! like rhovxc, but for valence density only
      real(8)    :: rveps     = 0d0  ! int rhov*exc(rhotot)
      real(8)    :: rvvxc     = 0d0  ! int rhov*vxc(rhotot)

!     Interaction of density with external fields
      real(8)    :: rvext(2,2)= 0d0  ! int n(1,2)_(val,cor) * vext
      real(8)    :: rhvxt1    = 0d0  ! int n1_val * vext
      real(8)    :: rhvxt2    = 0d0  ! int n2_val * vext
      real(8)    :: rhcvxt1   = 0d0  ! int n1_cor * vext for sites w/ lfoca=0
      real(8)    :: sgpote    = 0d0  ! integrals of compensating gaussians times vext
      real(8)    :: bfam      = 0d0  ! sum (bext * local moment)/2 --- for double counting.
                                     ! NB: only correct so far for collinear, l-independent bfield

      end type str_atparms

      type s_lgf
      complex(8), pointer :: gll(:,:)
      end type s_lgf
