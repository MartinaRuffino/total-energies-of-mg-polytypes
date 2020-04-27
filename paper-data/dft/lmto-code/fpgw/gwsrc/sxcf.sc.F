      subroutine sxcf_sc(qip,nq,itq,ntq,nbandmx,qlat,qibz,qbz,wbz,nstbz,irkip,nrkip,
     .  nqibz,nqbz,nctot,nl,nnc,nclass,natom,nband,nlmto,nsp,isp,
     .  iclass,mdim,nlnmv,nlnmc,nblochpmx,icore,ncore,imdim,ecore,
     .  ppb,nlnmx,mdimx,nbloch,wgt0,nq0i,q0i,ngrp,symgr,shtvg,tiat,miat,
     .  invg,il,in,im,nlnm,nn_,lx,nx_,nxx_,cgr,ppbrd,
     .  exchange,ifexsp,screen,cohtest,
     .  wklm,lxklm,iSigMode,nbmx,ebmx,niw,freqx,wx,dw,
     .  ef,esmr,diagonly,eftrue,
     o  zsec,coh)
C- Calculate full sigma_ij(e_i) = <i|Re[S](e_i)|j>
C ----------------------------------------------------------------------
Ci Inputs
Ci   qip    : list of irreducible k-points for which to calculate <i|sigma|j>
Ci   nq     : number of qp-points for which self-energy is made
Ci   itq    : list of states i which to calculate <i|sig>j>.
Ci          : Since this code makes full sigma, itq should be 1,2,3,...ntq
Ci          : REPLACE with [0] and remove itq
Ci   ntq    : dimensions zsec; maximum number of states for <i|sigm|j> calculated
Ci          : Eventually rename to: ndsec
Ci  nbandmx : nbandmx(iq) number of bands for matrix elements <i|sigma|j> at qp iq
Ci   qlat   : primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   qibz   : k-points for the irreducible BZ
Ci   qbz    : k-points for the full BZ
Ci   wbz    : sampling weight for qbz
Ci   nstbz  : (?) see genqbz
Ci   irkip  : If ip = index to some k where sigma is made and kx = index to some q in the IBZ,
Ci          : kr = irkip(k,ig,ip) = index to q in the FBZ corresponding to kx for group op ig
Ci   nrkip  : number of kp in FBZ for rotation ig and k in IBZ (?)
Ci   nqibz  : number of k-points in the irreducible BZ
Ci   nqbz   : number of k-points in the full BZ
Ci   nctot  : number of allowed core states
Ci   nl     : (global maximum l) + 1
Ci   nnc    : max number of core states
Ci   nclass : number of inequivalent classes
Ci   natom  : number of atoms
Ci   nband  : dimensioning parameter: maximum number of bands at any k.
Ci          : Eventually rename to: ndham
Ci   nlmto  : number of augmentation channels ... called ndima in sugw.f
Ci          : note: can be bundled in a structure, as it does not dimension anything.
Ci          : Eventually rename to: naugm
Ci   nsp    : 2 for spin-polarized case, otherwise 1
Ci   isp    : current spin channel (1 or 2)
Ci   iclass : the jth atom belongs to class iclass(j)
Ci   mdim   : number of product basis function for class i
Ci          : Eventually rename to: ntpba
Ci   nlnmv  : total number of valence partial waves for each site (for dimensioning)
Ci   nlnmc  : total number of core partial waves for each site (for dimensioning)
Ci nblochpmx: largest dimension of mixed product basis functions
Ci          : Eventually rename to: npbmbx
Ci   icore  : index for core states
Ci   ncore  : no. core states in each class
Ci   imdim  : imdim(i) = index to first product basis function for site i
Ci          : Used when all functions are strung together as a long vector
Ci          : Eventually rename to: ipbfn
Ci   ecore  : core eigenvalues
Ci   nlnmx  : maximum number of l,n,m
Ci   ppb    : matrix elements of product basis functions <phi(RLn) phi(RL'n') B(R,i)>
Ci          : ppb is not passed to sxcf, but it is made in ppbafp_v2 as CG * ppbrd
Ci   mdimx  : maximum value of mdim; third dimension of ppb
Ci          : Eventually rename to: mxnpba
Ci   nbloch : total number of local product basis functions within MTs
Ci          : Eventually rename to: npbmb
Ci   wgt0   : weight for each k-point near 0 (offset gamma method).  See q0iwgt3
Ci   nq0i   : number of offset points
Ci   q0i    : vector of offset k-points
Ci   ngrp   : number of space group operations
Ci   symgr  : rotation part of space group
Ci   shtvg  : translation part of space group, called delta in mptauof and ag in mksym
Ci   tiat   : group op ig maps ib -> ib''. ib'' equivalent to basis atom ib'=miat(ibas,ig)
Ci          : but shifted by translation vector tiat.  See mptauof.
Ci          : tiat must be an integer multiple of plat.  If shtvg changes by
Ci          : some lattice vector, tiat must change by the same amount.
Ci   miat   : miat(ibas,ig); is atom that ibs is rotated into (and shifted by tiat)
Ci   invg   : maps point group operation to its inverse
Ci   il     : l quantum number for Phi, for each function in product basis.  See m_genallcf_v3
Ci   in     : n for Phi, for each function in product basis.  See m_genallcf_v3
Ci   im     : m quantum number for Phi, for each function in product basis.  See m_genallcf_v3
Ci   nlnm   : number of (valence  + core; nlnmv and nlnmc contain separate valence, core info)
Ci   nn_    : dimensioning parameter max number of radial product functions for a given l
Ci          : Eventually rename to: ndphi
Ci   lx     : maximum l in product basis for each site.  In practice, lx=2*(nl-1)
Ci          : This array can be eliminated
Ci   nx_    : number of radial product basis functions for each l and site
Ci          : Eventually rename to: nrbl
Ci   nxx_   : Dimensions ppbrd; must be at least the maximum value of nx_
Ci   cgr    : Clebsch-Gordan coefficients rotated by symmetry group operations
Ci   ppbrd  : radial parts of product basis functions; read from files PPBRD_V2_*
Ci          : Eventually rename to: prbasr
Ci  exchange: T, Calculate the exchange part of the self-energy
Ci          : F, Calculate correlated part of the self-energy
Ci   ifexsp : file logical unit for writing the spectrum function for exchange.
Ci          : 0 => spectrum function not written
Ci   screen : for screened exchange (not implemented now)
Ci   cohtest: T, evaluate self-energy by screened exchange - coulomb-hole partitioning
Ci          : not implemented now
Ci   wklm   : coefficient w_L in Eqns. 28,33 of JPSP 83, 094711 (2014)
Ci   lxklm  : dimensions wklm
Ci  iSigMode: Mode determining form of QSGW approximation to Sigma
Ci   nbmx   : nbmx(1) : max number of bands <i|sigma|j> to calculate
Ci          : nbmx(2) : max number of bands to include in G, calculating G * W
Ci   ebmx   : ebmx(1) : exclude states above ebmx(1) in calculating <i|sigma|j>
Ci          : ebmx(2) : exclude states above ebmx(2) in G, calculating G * W
Ci   niw    : number of frequencies along the imaginary axis
Ci   freqx  : Imaginary frequencies on gaussian quadrature, points x between (0,1)
Ci          : Actual frequencies are Im omega = (1-x)/x; thus x=1/(1+omega)
Ci   wx     : Corresponding weights
Ci   dw     : spacing for frequency mesh along real axis
Ci   ef     : Fermi level
Ci   esmr   : broadening parameter for pole in G (?)
Ci  diagonly:
Ci   eftrue : T => fermi energy for valence exchange and correlation mode (added jan2013)
Ci          : ef is not the true fermi energy for core-exchange mode.
Co Outputs
Co   zsec   : self-energy S_ij= <i|Re[S](e_i)|j>
Co          : i and j are band indices.
Co          : zsec from im-axis integral part was originally symmetrized as
Co          :   wtt*.5d0*(   sum(zwzi(:,itp,itpp)) + !S_{ij}(e_i)
Co          :   dconjg( sum(zwzi(:,itpp,itp)) ))     !S_{ji}^*(e_j)= S_{ij}(e_j)
Co          : However, just the 1st term is calculated now:
Co          :   wtt* sum(zwzi(:,itp,itpp))   !S_{ij}(e_i)
Co          : This is OK because hqpe.sc.F does the symmetrization.
Co          : See also mkseci.F.
Co   coh    : For coulomb hole.  Not working now.
Cl Local variables
Cl  ngc     : number of G vectors for Coulomb interaction and
Cl          : number of interstitial product functions (q dependent)
Cl          : Eventually rename to: npwmb
Cl  ngb     : nbloch + ngc = total number of mixed product basis functions I,J in which <I|v|J> is calculated
Cl          : Eventually rename to: nmb
Cl  nw      : number of frequencies along the real axis
Cl          : Eventually rename to: nrw
Cl  nwxi    : the first frequency required for pole contribution to sigma
Cl          : and a dimensioning parameter for zw3, below
Cl          : Eventually rename to: nrwmin
Cl  nwx     : the last frequency required for pole contribution to sigma
Cl          : and a dimensioning parameter for zw3, below
Cl          : Eventually rename to: nrwmax
Cl  freqw   : Actual frequencies on Im-axis, freqw=(1-freqx)/freqx; see wintzg.
Cl  ua_auto : T, autodetermine Gaussian parameter a
Cl  ua      : Parameter in exp(-ua^2 w'^2) used in numerical integration of Sigma (wintzg)
Cl  ntqxx   : number of bands for <i|sigma|j>.
Cl  wtt     : Integration weight; see zsec above
Cl  cphikq  : Coefficients to partial waves (phi,phidot) for intermediate states at q-rk
Cl          : Eventually rename to: cphiqk
Cl  iatomp  : iatomp(i) : points to site i' which site i rotates into under group op
Cl  freqr   : frequencies along real axis
Cl          : OLD: freq_r
Cl  nt0p    : number of states below ef+ddw*esmr (ddw=10)
Cl          : Eventually rename to: noccp
Cl  nt0m    : number of states below ef-ddw*esmr (ddw=10)
Cl          : Eventually rename to: noccm
Cl  nt_max  : number of states contributing to pole in omega integration
Cl          : essentially similar to noccp ... not sure of the difference
Cl          : Eventually rename to: noccs
Cl  nbmax   : number of intermediate band states for SE
Cl          : = noccp-nctot in Exchange mode
Cl          : More complicated formula with nbmx,ebmx cutoff
Cl          : Eventually rename to: niband
Cl  nctot   : total number of core states to include in the calculation of G*W
Cl  nt      : number of intermediate states (n', Eq. 34) = nctot+niband
Cl          : For exchange mode, includes only occupied states
Cl          : Eventually rename to: nistate
Cl  nstate  : same as nt
Cl          : Eventually rename to: nistate
Cl  ntqxx   : number of states which <i|sigma|j> calculated for current q
Cl          : ntsec is limited to states which satisfy  omega-freqw<ebmx(1)  and  ntsec<=ndsec
Cl          : Eventually rename to: ntsec
Cl  qk      : q - qbz_kr
Cl  ekq     : evals at q-k
Cl  zmel   : Matrix elements of mixed basis, <psi(q,t') | psi(q-k,t) M(k,I)>, Eq. 34, PRB76, 165106
Cl          : Stored as zmel(igb(qrk), iocc(q), iunocc(q-qrk))
Cl          : Set mkmelt.F
Cl          : RENAME as zmel0
Cl  zmeltt  : zmeltt(it,itp,ivc) = S[I] zmel0(I,it,itp) * ppovlz(I,ivc)
Cl          : RENAME as zmelc
Cl  zmel    : zmel(ivc,it,itp) = S[I] zmel0(I,it,itp) * ppovlz(I,ivc)
Cl          : RENAME as zmelc (?) Note transpose of zmeltt
Cl  zw      : screened coulomb interaction for one frequency,
Cl          : read by direct access from file WVR.xxxx or WVI.xxxx
Cl  zw3     : screened coulomb interaction on the real frequency mesh.
Cl          : It is a symmetrized form of zw, and kept for a cluster of
Cl          : frequencies in the vicinity of Ef that are needed to calculated
Cl          : the contribution to sigma from the poles of G on the real axis.
Cl  ekq     : eigenvalues at q-k
Cr Remarks
Cr   The FP version of this code was written by T. Kotani (2000), adapted from
Cr   an ASA based GW code written by F. Aryasetiawan.
Cr   S. Faleev (2006) included the full sigma, used for QSGW
Cr   M. van Schilfgaarde (2013) redesigned for efficient implementation.
Cr   Martin Leuders (2014) developed a MPI version
Cr   It assumes that this routine is called on all nodes and parallelization is done internally
Cr
Cr  Main steps of algorithm for correlation part of self-energy
Cr    A.  Make mixed basis overlap matrix ppovlz (Made first, because fixed through rotations of k)
Cr        The full ppovlz is ppovlz(nmb,nmb)
Cr        MPI0 : Keep full matrix
Cr        MPI1 : Partition as ppovlz(I1:I2,J1:J2)
Cr
Cr    B.  Make matrix elements zmel0 (mkmelt)
Cr        The full zmel0 is zmel0(nmb,nistate,ntsec)
Cr        MPI0: Partition as zmel0(nmb,nilo:nihi,ntsec)
Cr        MPI1: Partition as zmel0(I1:I2,nilo:nihi,ntsec)
Cr
Cr    C.  Scale zmel0 by ppovlz : zmelc(it,itp,J) = S[I] zmel0(i,it,itp) * ppovlz(I,J)
Cr        The full zmelc is zmelc(nistate,ntsec,nmb)
Cr        MPI0: Partition as zmelc(nilo:nihi,ntsec,nmb)
Cr        MPI1: Partition as zmelc(nilo:nihi,ntsec,J1:J2)
Cr              zmelc and ppovlz are stored in slices => multiply them with BLACS call
Cr
Cr        Below this point MPI0 and MPI1 are the same
Cr
Cr    D.  Make zwz0 = zmelc*(W(0)-v)*zmelc (matzwz3)
Cr        Full matrix is zwz0(nistate,ntsec,ntsec))
Cr        MPI1: Partition as zwz0(nilo:nihi,ntsec,ntsec))
Cr
Cr    ... Im axis contribution to self-energy, steps E,F :
Cr        [note that there is a parallelizable loop over frequencies for steps E,F]
Cr    E.  Make zwz = zmel*(W(iw)-v)*zmel  (matzwz3; this is the most expensive step)
Cr        MPI1: Partition as zwz(nilo:nihi,ntsec,ntsec))
Cr
Cr    F.  Add contribution to zsec for this frequency (mkseci)
Cr        Uses zwz and zwz0 ... contributions from each intermediate state are independent
Cr
Cr    ... Re axis contribution to self-energy from poles of G in window between omega and Ef
Cr    G.  Add contribution to zsec for poles on real axis (mksecp)
Cr
Cr  The following comes from the original document for sce.f (correlation case), from F.Aryasetiawan.
Cr  The current GW code is very different now, but the documentation is kept as a reference.
Cr  92.02.24
Cr  93.10.18 from sec.f modified to take into account equivalent atoms
Cr  calculates the correlated part of the self-energy SE
Cr  SEc(q,t,t') = <psi(q,t) |SEc| psi(q,t'>
Cr  SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
Cr
Cr  The zeroth order Green function
Cr  G(r,r';w)   = S[occ]   psi(kn,r) psi(kn,r')^* /(w-e(kn)-i*delta)
Cr              + S[unocc] psi(kn,r) psi(kn,r')^* /(w-e(kn)+i*delta)
Cr
Cr  The screened coulomb potential
Cr  Wc(r,r';w)  = W(r,r';w) - v(|r-r'|)
Cr              = < [r1,r2] v(|r-r1|) X(r1,r2;w) v(|r2-r'|) >
Cr  W(r,r';w)   = < [r''] eps^-1(r,r'';w) v(|r''-r'| >
Cr  eps^-1      = 1 + v*X = inverse dielectric matrix
Cr  eps         = 1 - v*X0 in RPA
Cr
Cr  Expand Wc(r,r';w) in mixed basis M = {PW, B} where B are atom-centered functions
Cr  Wc(r,r';w)  = S[k=FBZ] S[i,j=1,nbloch]
Cr                M(k,I,r) Wc(k,w)(I,J) M(k,J,r')^*
Cr  Wc(k,w)(I,J) are  the matrix elements of Wc in M
Cr
Cr  SEc(q,t,t') = S[k=FBZ] S[n=occ] S[i,j=1,nbloch]
Cr         <psi(q,t) |psi(q-k,n) M(k,I)> <M(k,J) psi(q-k,n) |psi(q,t')>
Cr         (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)-i*delta)>
Cr
Cr              + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
Cr         <psi(q,t) |psi(q-k,n) M(k,I)> <M(k,J) psi(q-k,n) |psi(q,t')>
Cr         (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)+i*delta)>
Cr
Cr --------------------------------------------
Cr  contribution from the poles of G: SEc(pole)
Cr --------------------------------------------
Cr
Cr  for w <= ef
Cr  SEc(pole) = - S[k=FBZ] S[n=occ] S[i,j=1,nbloch]
Cr         <psi(q,t) |psi(q-k,n) M(k,I)> <M(k,J) psi(q-k,n) |psi(q,t')>
Cr              Wc(k,e(q-k,n)-w)(i,j) theta(e(q-k,n)-w)
Cr
Cr  for w > ef
Cr  SEc(pole) = + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
Cr         <psi(q,t) |psi(q-k,n) M(k,I)> <M(k,J) psi(q-k,n) |psi(q,t')>
Cr              Wc(k,w-e(q-k,n))(i,j) theta(w-e(q-k,n))
Cr
Cr  theta(x)  = 1   if x > 0
Cr            = 1/2 if x = 0
Cr            = 0   if x < 0
Cr
Cr  the analytic structure of GWc for w <= ef
Cr                                |
Cr                                |   o = pole of G
Cr                                ^   x = pole of Wc
Cr                                |
Cr                                |   ef-w
Cr                                |----<-----
Cr                                |          |
Cr                  o  o  o  o  o |o  o  o   ^
Cr                x  x  x  x  x  x|          |
Cr   -----------------------------|---->------------------------------
Cr                                |x  x  x  x  x  x  x  x
Cr                                |              o  o  o  o  o
Cr                                |       <----->
Cr                                ^        gap in insulator
Cr                                |
Cr                                |
Cr
Cr  the analytic structure of GWc for w > ef
Cr                                |
Cr                                |   o = pole of G
Cr                                |   x = pole of Wc
Cr                                |
Cr          gap in insulator      ^
Cr                 <----->        |
Cr       o  o  o  o               |
Cr          x  x  x  x  x  x  x  x|
Cr   ------------------------>----|-----------------------------------
Cr                    |           |x  x  x  x  x  x  x  x
Cr                    ^   o  o  o  o  o  o  o
Cr                    |           |
Cr                     ------<----|
Cr                        w-ef    |
Cr                                ^
Cr                                |
Cr
Cr  Integration along the real axis from -inf to inf is equivalent to
Cr  the integration along the path shown.
Cr
Cr  For integration along the imaginary axis (see wintz0, wintzg)
Cr    (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
Cr  The i*delta becomes irrelevant.
Cr  Omit k and basis index for simplicity and denote e(q-k,n) = e
Cr    wint = (i/2pi) < [w'=-inf,inf] Wc(w')/(w+w'-e) >
Cr
Cr  w' ==> iw', w' is now real and
Cr    wint = - (1/pi) < [w'=0,inf] Wc(iw') (w-e)/{(w-e)^2 + w'^2} >
Cr
Cr  Integral is done in practice by transforming to: x = 1/(1+w')
Cr  this leads to a denser mesh in w' around 0 for equal mesh x.
Cr  Desirable since Wc and the lorentzian are peaked around w'=0.
Cr  See routines wintz0 or wintzg for details.
Cr
Cr  This together with the contribution from the pole of G (s.u.)
Cr  gives the so called static screened exchange -Wc(0)
Cr
Cu Updates
Cu   05 May 13
C ----------------------------------------------------------------------
      use m_readqg
      use m_readeigen,only: readeval,readcphi,readgeig
      use keyvalue
C     use m_nbandgw,only: getnbandgw
      use mpi_mod, q_comm => comm_1, q_inter => inter_1,
     &             k_comm => comm_2, k_inter => inter_2,
     &             rank => mpi_rank
      use fcugw
      implicit none
C ... Passed parameters
      logical :: exchange,screen,cohtest
      integer :: ntq,natom,nqbz,nqibz,ngrp,nq,nw,niw,nband,nlmto,nq0i,nctot,isp,nsp,nlnmx,
     .  nbloch,nclass,nl,nnc,nblochpmx,nxx_,mdimx,ifexsp,lxklm
      integer mdim(natom),il(*),in(*),im(*),nlnm(*),nn_,lx(*),nx_(*),nbmx(2),itq(ntq),
     .  iatomp(natom),miat(natom,ngrp),irkip(nqibz,ngrp,nq),nrkip(nqibz,ngrp,nq),invg(ngrp)
      integer :: nbandmx(nq)    !sep2012 number of bands for calculated self-energy.
      integer :: nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)
      real(8) :: ef,eftrue,esmr,dw,qlat(3,3),tiat(3,natom,ngrp),tr(3,natom),
     &  wbz(nqbz),qibz(3,nqibz),qbz(3,nqbz),freqx(niw),wx(niw),
     &  ekq(nband),ekc(nctot+nband),coh(ntq,nq),wklm((lxklm+1)**2),cgr(*),ppbrd(*),qip(3,*)
      real(8) :: wgt0(nq0i,ngrp),symgr(3,3,ngrp),q0i(1:3,1:nq0i),shtvg(3,ngrp),ecore(nctot),ppb(*)
      complex(8) :: zsec(ntq,ntq,nq)
C ... Local parameters
      logical :: GaussSmear,eibz4sig,diagonly,ltmp,smbasis,ua_auto,test_symmetric_W,onceww
      logical :: debug=.false.,initp=.true.
      character(5):: charnum5
      integer GWversion,i,iSigMode,ia,iclose,idalloc,ifrcw,ifrcwi,ifvcoud,
     .  ifwd,igb1,igb2,invr,invrot,iopen,ip,iqend,iqini,irot,it,itend,itini,itp,
     .  itpp,ivc,iverbose,iw,iwp,ix,ixs,ixsmin,ixsmx,kr,kx,mrecl,nbcut,nbcutc,
     .  nbmax,ngb,ngb0,ngc,nlmtobnd,nocc,npm,nprecx,nrec,nrot,nstate,nstatetot,
     .  nstbz(nqbz),nt,nt0,nt0m,nt0p,nt_max,ntqxx,nwr0,nwx,nwxi,verbose
      integer ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi
      integer :: num_q_groups, num_k_groups, ipstart, ipend, ipstep, ikstart, ikend, ikstep

      real(8) :: ddw,det,efm,efp,esmrx,fpi,freqw,omg,pi,tpi,ua_,vc,we,weavx2,
     .  wexx,wfac,wfacx,wfacx2,wgtq0p,wtt
      real(8) :: wfaccut=1d-10
      real(8) :: symope(3,3),qxx(3),shftv(3),q(3),qk(3),qlatinv(3,3),expa_(niw),
     .  omega(ntq),ebmx(2),quu(3),qibzk(3),qrk(3),qvcoud(3),eq(nband)

      real(8),allocatable :: freqr(:),vcoud(:),vcoudsq(:)
      complex(8), allocatable :: cphiq(:,:),cphikq(:,:)
      complex(8), allocatable :: zmel(:,:,:),tmp(:,:,:) ! ,zmeltt(:,:,:)
      complex(8), allocatable :: vcoul(:,:),w3p(:,:,:),ppovl(:,:),ppovlz(:,:),zcousq(:,:),zwz0(:,:,:)
      complex(8), allocatable :: z1r(:,:),z2r(:,:),zw3(:,:,:),zw3x(:,:),w3pi(:,:),zw(:,:),zwz(:,:,:)

      complex(8) :: zfac

      integer(8) :: cr,cm,c1,c2

      integer :: nblk

      external GaussSmear

      call system_clock(count_rate = cr, count_max = cm)


C --- Setup ---
      pi  = 4d0*datan(1d0)
      tpi = 8d0*datan(1d0)
      fpi = 4d0*pi
      if (verbose()>=90) debug=.true. !; debug=.true.
C      if (bzcase()==1) then
        if (abs(sum(qibz(:,1)**2))/=0d0) call rx('sxcf assumes 1st qibz/=0 ')
        if (abs(sum( qbz(:,1)**2))/=0d0) call rx('sxcf assumes 1st qbz /=0 ')
C      endif
!     Next line was setup prior to newaniso
!     iqini=2; if (bzcase()==2) iqini=1; iqend= nqibz+nq0i
      iqini = 1
      iqend = nqibz
      if (mdimx /= maxval(mdim)) call rx('sxcf: wrong mdimx')
      if (sum(mdim(iclass(1:natom))) /= nbloch) call rx('mkmelt: product basis mismatch')

      if (rank == 0) call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
#ifdef USE_MPI
      call MPI_Bcast(nbcut,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
      nlmtobnd  = nlmto*nband
      nstatetot = nctot + nband
      nbcutc = nctot+nbcut
      if (nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
C     call dinv33(qlat,0,qlatinv,det)
      call minv33(qlat,qlatinv)

      allocate(cphiq(nlmto,nband), cphikq(nlmto,nband))

C --- Correlation case: read general real frequency mesh from file 'freq_r' ---
C     Note freqr(0)=0 (from feb2006)
C     We assume freqr(iw) == -freqr(-iw), so freqr(iw<0) is a little redundant
      if (.not.exchange) then
        if (rank == 0) then
          open(UNIT=3111,file='freq_r')
          read(3111,*) i,nwr0
        endif
#ifdef USE_MPI
        call MPI_Bcast(i, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(nwr0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
        nw = i-1
        allocate(freqr(nwr0:nw))
        if (rank == 0) then
          do iw = nwr0, nw
            read(3111,*) freqr(iw)
          enddo
          close(3111)
        endif
#ifdef USE_MPI
        call MPI_Bcast(freqr, nw-nwr0+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
#endif
        if (nwr0/=0) then
          if (nw/= -nwr0)       call rx("sxcf_fal3_scz: nw/=-nwr0")
          if (freqr(0)/=0d0)    call rx("sxcf_fal3_scz: freqr(0)/=0")
          if (sum(abs(freqr(1:nw)+freqr(-1:-nw:-1)))/=0) call rx("sxcf_fal3_scz: freqr /= -freqr")
        endif

        if (rank == 0) then
          ifwd = iopen('WV.d',1,-1,0)
          read (ifwd,*) nprecx,mrecl
          ifwd = iclose('WV.d')
        endif
#ifdef USE_MPI
        call MPI_Bcast(nprecx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(mrecl, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
      endif ! .not.exchange

! --- ip loop to make zsec(:,:,ip). Uses symmetry to reduce number of q ---
      num_q_groups = mpi_N1
      call mpi_make_groups(comm_world,nq,num_q_groups,q_comm,q_inter)

#ifdef USE_MPI
      write(*,'("MPI: rank = ",I4,", q_comm%ID = ",I4,", q_comm%group_ID = ",I4)')
     &          rank, q_comm%ID, q_comm%group_ID

      write(*,'("MPI: rank = ",I4,", q_inter%ID = ",I4,", q_inter%group_ID = ",I4)')
     &          rank, q_inter%ID, q_inter%group_ID
#endif

      num_k_groups = mpi_N2
      call mpi_make_groups(q_comm, iqend-iqini, num_k_groups, k_comm, k_inter)
      ipstart = 1 + q_comm%group_ID
      ipend   = nq
      ipstep  = num_q_groups
      do  1001 ip = ipstart, ipend, ipstep
        if (sum(irkip(:,:,ip))==0) cycle
        q(1:3) = qip(1:3,ip)
        if (q_comm%ID == 0) then
          write(*,"(' do 1001: iq=',i3,' of',i3,'  q=',3f12.6)") ip,nq,q
          call cputid(0)
        endif

        call readeval(q,isp,eq)
        do  i = 1, ntq
          omega(i) = eq(itq(i)) !sf 21May02 + 2d0*(dble(iw)-shtw)*deltaw
        enddo

        call readcphi(q,nlmto,isp,quu,cphikq)
        do  it = 1, ntq ! augmentation part of eigenfunction for itp at q
          itp = itq(it)
          cphiq(1:nlmto,it) = cphikq(1:nlmto,itp)
        enddo

!       gauss_img : interpolation gaussion for W(i \omega).
        if (q_comm%ID == 0) call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
#ifdef USE_MPI
        call MPI_Bcast(ua_,1,MPI_DOUBLE_PRECISION,0,q_comm%communicator,ierror)
#endif

        if (ua_ <= 0d0) then
          ua_auto = .true.
          write(6,"(' ua_auto=T')")
        else
          ua_auto = .false.
          do  ix = 1, niw
            freqw     = (1d0 - freqx(ix))/freqx(ix)
            expa_(ix) = exp(-(ua_*freqw)**2)
          enddo
        endif

        if (debug .and. q_comm%ID == 0) then
          do ix = 1, niw
            write(6,"(' expa=',10d13.6)") freqx(ix),expa_(ix)
          enddo
        endif

C   --- Big loop over irreducible k-points kx, index to q ---
        ikstart = iqini + k_comm%group_ID
        ikend   = iqend
        ikstep  = num_k_groups
        do 1100 kx = ikstart, ikend, ikstep
          if (sum(irkip(kx,:,ip))==0) cycle
          iverbose = (verbose()/30)*10
          if (kx /= iqini .or. ip > 1) iverbose = 0
          if (k_comm%ID == 0) then
            write (*,'(" do 1100 kx",i3,"  out of ",i3,"  k-points "$)') kx,iqend
            call cputid(0)
          endif
          if (kx <= nqibz) then
            qibzk = qibz(:,kx)
          else
            qibzk = 0d0
          endif

C     ... Get the number of mixed basis functions for kx
          call getngvec('QGcou',qibzk,qlatinv, quu,ngc)
          ngb = nbloch + ngc
          if (debug .and. k_comm%ID == 0) print *, 'sxcf: ngb=',ngb,nbloch

!     --- Read diagonalized Coulomb interaction from file Vcoud.ID  where ID=kx ---
!!        Note sep2012 takao
!!        Vcoud is a sequential file containing the coulomb matrix for qibzk.
!!        It is is written in hvccfp0.m.F.
!!        For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
          if (kx<=nqibz) qxx = qibzk
          if (kx>nqibz ) qxx = q0i(:,kx-nqibz)
          if (k_comm%ID == 0) ifvcoud = iopen('Vcoud.'//charnum5(kx),0,0,0)
          do
            if (k_comm%ID == 0) then ! only k-head reads
              read(ifvcoud) ngb0
              read(ifvcoud) qvcoud
            endif
#ifdef USE_MPI
            call MPI_Bcast(ngb0, 1, MPI_INTEGER, 0, k_comm%communicator, ierror)
            call MPI_Bcast(qvcoud,  3, MPI_DOUBLE_PRECISION, 0, k_comm%communicator, ierror)
#endif
            if (allocated(vcoud)) deallocate(vcoud)
            allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
            if (k_comm%ID == 0) then
              read(ifvcoud) vcoud
              read(ifvcoud) zcousq
            endif
#ifdef USE_MPI
            call MPI_Bcast(vcoud,  ngb0,    MPI_DOUBLE_PRECISION, 0, k_comm%communicator, ierror)
            call MPI_Bcast(zcousq, ngb0**2, MPI_DOUBLE_COMPLEX, 0, k_comm%communicator, ierror)
#endif

            if (sum(abs(qvcoud-qxx))<1d-6) goto 1133
          enddo ! done reading vcoud, zcousq

C     ... Sanity check
          if (sum(abs(qvcoud-qxx))>1d-6) then
            print *,'qvcoud =',qvcoud
            print *,'qxx=',qxx,kx
            call rx( 'sxcf: qvcoud/=qibz(:,kx) hvcc is not consistent')
          endif
 1133     continue
          if (ngb0/=ngb ) then  !sanity check
            write(0,'("[",I4,"] qxx = ",I5,", ngb0 = ",I6", ngb = ",I6)') rank,qxx,ngb0,ngb
#ifdef USE_MPI
            call MPI_abort(MPI_COMM_WORLD,-1,ierror)
#endif
            call rx( 'hsfp0.m.f:ngb0/=ngb')
          endif

C     --- Make overlap matrix (Step A, Remarks) ---
!!        <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!!        zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
!!        ppovlz= O Zcousq
!!        (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
          if (allocated(ppovlz)) deallocate(ppovlz)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
          if (k_comm%ID == 0) call readppovl0(qibzk,ngc,ppovl)
#ifdef USE_MPI
          call MPI_Bcast(ppovl,ngc*ngc,MPI_Double_Complex,0,k_comm%communicator,ierror)
#endif
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
C         ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
          print *, 'ngc, ngb', ngc, ngb
          call system_clock(c1)
          call zgemm('N','N',ngc,ngb,ngc,(1d0,0d0),ppovl,ngc,zcousq(nbloch+1,1),ngb0,(0d0,0d0),ppovlz(nbloch+1,1),ngb)
          call system_clock(c2)
          write(*,'("sxcf ppzc:",x,g20.12)') real(c2 - c1,8)/real(cr,8)
C         call zprm('ppovlz',ppovlz,ngb,ngb,ngb)
          deallocate(zcousq,ppovl)
          if (verbose()>=50) call cputid2('complete <I|v|J>',0)

!     ... Open WVR,WVI for correlation mode
          if (.not.exchange) then
            if (k_comm%ID == 0) then  ! MPI:  files are only opened on head node of the k_groups.
              ifrcw  = iopen('WVR.'//charnum5(kx),0,-1,mrecl)
              ifrcwi = iopen('WVI.'//charnum5(kx),0,-1,mrecl)
            endif
          endif

!     ... Find nrot = number of k-points kx is rotated into (for printout only)
          nrot = 0
          do irot = 1,ngrp
            if (kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if (kr == 0) cycle   ! next irot
              qrk = qbz(:,kr)
            else
              kr = -99999 !for sanity check
              qrk = 0d0
              if (wgt0(kx-nqibz,irot) == 0d0) cycle ! next irot
            endif
            nrot = nrot+1
          enddo

!     --- Loop over rotations of k to full BZ ---
!         rk (also called qrk) is rotated k; kr is index in qbz to rk
          do 1000 irot = 1, ngrp
            if (kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if (kr==0) cycle   ! next irot
              qrk = qbz(:,kr)
            else
              kr = -99999 !for sanity check
              qrk = 0d0
              if (wgt0(kx-nqibz,irot) == 0d0) cycle ! next irot
            endif
            if (debug) write(*,"('ip,kx irot=',3i5,' out of',2i4)") ip,kx,irot,iqend,nrot

!       ... Rotate atomic positions invrot*R = R' + T
            invr   = invrot(irot,invg,ngrp)  ! s_sym(ig)%inv
            symope = symgr(:,:,irot)         ! symgr = s_sym(ig)%rg
            shftv  = matmul(symope,shtvg(:,invr)) ! shtvg = s_sym(ig)%rg but may differ by some plat
            tr(1:3,:)= tiat(1:3,:,invr)        ! s_sym(invr)%tib(1:3,:) but may differ by shtvg - s_sym(ig)%rg
            iatomp = miat(:,invr)            ! s_sym(ig)%mib(:)
C           call system_clock(c1)
!       ... Get ppb = <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
            call ppbafp_v2(irot,ngrp,isp,nsp,
     i           il,in,im,nlnm,
     d           nl,nn_,nclass,nlnmx,
     i           mdimx,lx,nx_,nxx_,
     i           cgr,nl-1,
     i           ppbrd,
     o           ppb)
C           call system_clock(c2)
C           write(*,'("sxcf ppba:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

!       ... qk = q-rk. rk is inside 1st BZ, but not restricted to the irreducible BZ
            qk = q - qrk        ! qbz(:,kr)
!!          Read eigenvalue and Cphi at q-rk
            call readeval(qk,isp,ekq)
            call readcphi(qk,nlmto,isp,quu,cphikq)

!!          No. occupied (core+valence) and unoccupied states at q-rk
            ekc(nctot+1:nctot+nband) = ekq(1:nband)
            nt0 = nocc(ekc,ef,.true.,nstatetot)
            ddw = .5d0
            if (GaussSmear()) ddw = 10d0
            efp = ef+ddw*esmr
            efm = ef-ddw*esmr
            nt0p = nocc(ekc,efp,.true.,nstatetot)
            nt0m = nocc(ekc,efm,.true.,nstatetot)

c     tetra
c     if (tetra) then
c     irkp = ibzx(kp)  ! IBZ corresponding to qk
c     do i=1,nband
c     if ( wtetef(i,irkp) == 0d0 ) exit
c     enddo
c     if (sum(abs(wtetef(i:nband,irkp)))/=0d0)
c     &   call rx(' sxcf: wtetef sum err1')
c     nt0org = nocc (ekc,ef,.true.,nstatetot)
c     nt0 = i-1 + nctot
c     endif

!! nbmx(1),ebmx(1): restricts number of bands of <i|sigma|j> to calculate.
!! nbmx(2),ebmx(2): restrict number of intermediate bands in calculating G W
            if (exchange) then
              nbmax = nt0p-nctot
            else
              nbmax = nband  ! intermediate states => all states initially
              i = nocc(ekc,ebmx(2),.true.,nstatetot)-nctot; nbmax  = min(nband,nbmx(2),i)
              if (initp) then
                if (k_comm%ID == 0) print *,' nbmax=',nbmax
                initp=.false.
              endif
            endif
!! ntqxx is number of bands for <i|sigma|j> for this qp
            freqw = eftrue; if (IAND(mod(GWversion()/10,10),1) == 0) freqw = 0  ! For backwards compatibility
            ntqxx = nocc(omega-freqw,ebmx(1),.true.,ntq)
                                !bug -ef is added jan2013
                                !previous version do not give wrong results, but inefficient.
            ntqxx = min(ntqxx,nbmx(1))
            if (ntqxx<nband) then
              do i=ntqxx,1,-1   !reduce ntqxx when band tops are degenerate sep2012
                if (omega(i+1)-omega(i)<1d-2) then
                  ntqxx=i-1
                else
                  exit
                endif
              enddo
            endif
            nbandmx(ip)=ntqxx   !number of bands to be calculated Sep2012.

            if (debug) print *,' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
            nstate = nctot + nbmax ! = nstate for the case of correlation

!       --- Matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)> (Step B, Remarks) ---
!           including the phase factor exp(ik.T)
!           B(rot*k,r) = B(k,invrot*r)
!           Stored in zmel. See description of zmel in the header
            nt = nctot + nbmax  ! = include cores in intermediate states for the correlation case ... note nt identical to nstate above
            if (nt*ntqxx > 0) ia = idalloc('zmel',iverbose+2,ngb*nt,ntqxx*2)
            allocate(zmel(ngb,ntqxx,nt))
            ltmp = verbose()>=40 .and. kx == iqini .or. verbose()>=41
            ipblo = 1
            ipbhi = ngb
            it1lo = 1
            it1hi = nt
            it2lo = 1
            it2hi = ntqxx

!           Debugging check
!           ipblo=2; ipbhi=4 ! Slice of pb with 3 states from atom 1
!           ipblo=4; ipbhi=nbloch+3 ! some core states from atom 2; valence and some PW states
!           it1lo = 5; it1hi = 9 ! End of cores, some valence states
!           it2lo = 2; it2hi = 5 ! some unocc states
!            print *, '!!', ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi
!            deallocate(zmel); allocate(zmel(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi))
            call system_clock(c1)
            call mkmelt(1,isp,isp,qlat,q,qrk,qibzk,cphiq,cphikq,ppb,icore,ncore,symope,shftv,tr,
     .        nl,natom,nclass,iclass,iatomp,nlnmx,nlnmv,nlnmc,imdim,mdim,mdimx,nlmto,ngc,nbloch,
     .        1,ntqxx,itq,1,1,nbmax,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,nnc,ltmp,zmel)
            call system_clock(c2)
            write(*,'("sxcf melt:",x,g20.12)') real(c2 - c1,8)/real(cr,8)
!             stop
!           Debugging check
!           ipblo=4; ipbhi=nbloch+3 ! Slice of pb with 3 states from atom 1
!           it1lo = 5; it1hi = 9 ! End of cores, some valence states
!           it2lo = 2; it2hi = 5 ! some unocc states
!            call zprm('zmel from mkmelt',zmel(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi),
!     .        ipbhi-ipblo+1,ipbhi-ipblo+1,(it1hi-it1lo+1)*(it2hi-it2lo+1))

            if (debug) write(6,"('mkmelt sum zmel, abszmel=',4d23.12)") sum(zmel),sum(abs(zmel))

            if (nbcut/=0.and.(.not.exchange)) then
              do it = max(it1lo,nctot+1),min(it1hi,nctot+nbcut,nctot+nbmax)
                zmel(:,:,it) = 0
              enddo
            endif

            if (smbasis()) then
               call rx('smbasis not supported now')
C              ntp0= ntqxx
C              nn= nnr(kx)
C              no= nor(kx)
C              allocate( pomat(nn,no) )
C              pomat= pomatr(1:nn,1:no,kx)
C              if ( sum(abs(qibzk-qrr(:,kx)))>1d-10 .and.kx <= nqibz ) then
C                call rx('qibz/= qrr')
C              endif
C              if (no /= ngb.and.kx <= nqibz) then
C!!     A bit sloppy check only for kx<nqibz because qibze is not supplied...
C                write(6,"(' q  ngb  ',3d13.5,3i5)")  qibzk,ngb
C                write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
C                call rx('x0kf_v2h: POmat err no/=ngb')
C              endif
C
C              ngb = nn          ! Renew ngb !!!
C              allocate ( zmeltt  (nn, nctot+nbmax, ntp0) )
C                  call matm( pomat, dcmplx(rmelt,cmelt), zmeltt,
C     &                 nn, no, (nctot+nbmax)*ntp0 )
C                 Try something like (note nt = (nctot + nbmax))
C                  call ZGEMM ( "N", "N", nn, nt*ntp0, no, (1d0,0d0),
C     .              pomat, nn,zmel, no,(0d0,0d0), zmeltt, nn )
C
C                  deallocate(rmelt, cmelt)
C                  allocate( rmelt(ngb, nctot+nbmax, ntp0), !ngb is reduced.
C     &                 cmelt(ngb, nctot+nbmax, ntp0) )
C                  rmelt = dreal(zmeltt)
C                  cmelt = dimag(zmeltt)
C                  deallocate(zmeltt,pomat)
C     1113 continue
C            else
C             nn=ngb
C             no=ngb
            endif
            if (debug) write(6,"('sum of zmel=',4d23.16)") sum(zmel)
!! === End of matrix element generation ===

C       ... Integration weights
C            if (bzcase()==2) then
C              if (kx<=nqibz) then
C                wtt = wbz(kr)
C                if (nstbz(kr)/=0) wtt = wbz(kr)*(1d0-wgtq0p()/nstbz(kr))
C              elseif (kx>nqibz) then
C                wtt = wgt0(kx-nqibz,irot)
C              endif
C            else
              if (kx<= nqibz) then
                wtt = wbz(kr)
              else
                wtt = wbz(1)*wgt0(kx-nqibz,irot)
                if (abs(wbz(1)-1d0/dble(nqbz))>1d-10) call rx('sxcf:wbz(1) inconsistent')
              endif
C            endif
            if (eibz4sig()) then
              wtt = wtt*nrkip(kx,irot,ip) !25sep2012
            endif

C       --- Exchange self-energy ---
            if (exchange) then
Cr            (need to check following comments)
Cr             S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
Cr             v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
Cr
Cr             > z1p(j,n,t) = S[i=1,nbloch] <psi(q,t) | psi(q-rk,n) B(rk,i)> v(k)(i,j)
Cr
Cr             write(6,*)' vcoulsum=',sum(vcoul)
Cr             if (debug) write(6,*)'  sumz=',dcmplx(rmelt,cmelt),sum(vcoul)
Cr
Cr        --- screened exchange case ---
Cr            if (screen) then
Cr              allocate( zw (nblochpmx,nblochpmx))
Cr              ix = 1
Cr              ! write(*,*)(kx-2)*(nw_w+1)+ix
Cr              read(ifrcw,rec=((kx-2)*nw+ix)) zw ! Readin W(0) - v         !sf 22May02
Cr              !nw is number of frequency points in general mesh: freqr(nw), freqr(1)=0
Cr              vcoul = vcoul + zw(1:ngb,1:ngb) !c  screen test
Cr              deallocate(zw)
Cr            endif
              if (debug .and. (k_comm%ID == 0)) then
                do it = 1, nctot+nbmax
                  write(6,"('it sum(zmel(it)',i5,2f10.4)") it, sum(abs(zmel(:,1,it)))
                enddo
                write(6,*)'sumcheck ppovlz=',sum(abs(ppovlz(:,:)))
                write(6,*)'sumcheck vcoud=',sum(abs(vcoud(:)))
              endif

              if (diagonly) call rx('sxcf_sc not ready for diagonly')
! ! Apparently nstate and nt hold the same value
!               allocate(zmeltt(ngb,ntqxx,nt))
!
! C             zmeltt(it,itp,ivc) = S[J] zmel(J,it,itp) * ppovlz(J,ivc)
!               call zgemm('T','N',ngb,ntqxx*nt,ngb,(1d0,0d0),ppovlz,ngb,zmel,ngb,(0d0,0d0),zmeltt,ngb)
! C             call zprm('zmeltt, vectorized',zmeltt,ntqxx*(nctot+nbmax),ntqxx*(nctot+nbmax),ngb)
!               if (verbose()>=30 .and. kx == iqini .or. verbose()>=40) call cputid2(' complete zmel * overlap',0)
!               if (nstate*ntqxx > 0) then
!                 ia = idalloc('zmel',iverbose+4,ngb*nstate,ntqxx*2)
!               endif
!               deallocate(zmel)

! This should really be exported to a routine of its own, maybe zgemm_inplace taking an extra argument for the size of the desired buffer, or the buffer itself (probably less useful since it may be very large and only used a couple of times).
              if (nstate > 0) then
              nblk =  min(nstate,max(1,((4*2**23)/(ngb*ntqxx))*8)) ! made for up to ~4GiB and multiple of 8
              write(*,'("sxcf nblk, tmpsz, zmel:",x,i0,",",x,i0," MiB",x,f6.1," GiB")') nblk,
     &           (ngb*ntqxx*nblk)/2**16, real(int(ngb,8)*ntqxx*nstate,8)/2**26

              call system_clock(c1)
              allocate(tmp(ngb,ntqxx,nblk))
              do it = 1, nstate, nblk ! this complication is done to reduce excessive memory usage while not damaging performance
                  if (nstate - it + 1 < nblk) nblk = nstate - it + 1 ! last iteration
                  call zgemm('t','n',ngb,ntqxx*nblk,ngb,(1d0,0d0),ppovlz,ngb,zmel(1,1,it),ngb,(0d0,0d0),tmp,ngb)
                  call zcopy(ngb*ntqxx*nblk,tmp,1,zmel(1,1,it),1)
              end do
              deallocate(tmp)
              call system_clock(c2)
              write(*,'("sxcf orth:",x,g20.12)') real(c2 - c1,8)/real(cr,8)
              endif


C             Vectorized version of w3p_t,t',t'' = sum_c z_t,t',c v_c z*_t,t'',c
C             Write as standard matrix product of t-dependent arrays
C             w3pi(t)_t',t'' = w3p_t,t',t''
C             z1(t)_t',c     = z_t,t',c v_c
C             z2(t)_t",c     = z_t,t'',c
C             w3pi(t)_t',t'' = sum_c z1(t)_t',c z2*(t)_t'',c
C             w3pi(t)        = z1(t) z2+(t)

              print *, 'sxcf: ngb,ntqxx,nt', ngb,ntqxx,nt
              ltmp = kx == iqini .and. wklm(1) < 0
C             ltmp = kx == iqini

              if (.not. ltmp) then  ! Dimitar's fast branch
              print *, 'sxcf: living in the fast lane'

              allocate(vcoudsq(ngb))
              do ivc = 1, ngb
                if (vcoud(ivc) < 0.0d0) call rx("sxcf: negative Vcou eval")
                vcoudsq(ivc) = sqrt(vcoud(ivc))
              end do
              if (kx == iqini) then
                vcoudsq(1) = wklm(1) * fpi*sqrt(fpi)/wbz(kx) ! Replace with vbar
                if (vcoudsq(1) < 0.0d0) call rx("sxcf: negative wklm(1)")
                vcoudsq(1) = sqrt(vcoudsq(1))
              end if

              do it = 1, nt
!                 call zdscal(ntqxx, vcoudsq(ivc), zmeltt(ivc,1,it), ngb)
                forall (ivc=1:ngb,itp=1:ntqxx) zmel(ivc,itp,it) = conjg(zmel(ivc,itp,it))*vcoudsq(ivc)
!                 forall (ivc=1:ngb,itp=1:ntqxx) zmeltt(ivc,itp,it) = zmeltt(ivc,itp,it)*vcoudsq(ivc)
                wfac = -wtt
                if (it > nctot) wfac = wfac*wfacx(-1d99, ef, ekc(it), esmr)

                call zherk('u', 'c', ntqxx, ngb, wfac, zmel(1,1,it), ngb, 1d0, zsec(1,1,ip), ntq)
!                 call zgemm('c','n', ntqxx, ntqxx, ngb, cmplx(wfac,0d0,8), zmeltt(1,1,it), ngb, zmeltt(1,1,it), ngb,
!      &                    cmplx(1d0,0d0,8), zsec(1,1,ip), ntq)
              end do

              if (nstate*ntqxx > 0) ia = idalloc('zmel',iverbose+4,ngb*nstate,ntqxx*2)
              deallocate(zmel,vcoudsq)

              else   ! Slow branch

!               sum_it,ig z_it,itp1,ig * v_ig * z_it,itp2,ig * w_it
!               sum_it w_it * sum_ig z_it,itp1,ig * v_ig * z_it,itp2,ig : nt * (ntp*ng*ntp + npt*ng)
!               sum_ig v_ig * sum_it z_it,itp1,ig * w_it * z_it,itp2,ig : ng * (ntp*nt*ntp + npt*nt)

              allocate(w3p(nctot+nbmax,ntqxx,ntqxx))
              allocate(z1r(ntqxx,ngb),z2r(ntqxx,ngb),w3pi(ntqxx,ntqxx))
              vc = vcoud(1)
              if (kx == iqini) vcoud(1) = wklm(1) * fpi*sqrt(fpi)/wbz(kx)  ! Replace with vbar

              do  it = 1, nctot+nbmax
                do  ivc = 1, ngb
                  do  itp = 1, ntqxx
                    z1r(itp,ivc) = zmel(ivc,itp,it) * vcoud(ivc)
                    z2r(itp,ivc) = zmel(ivc,itp,it)
                  enddo
                enddo
C               Straight inline
!                call zgemm('N','C',ntqxx,ntqxx,ngb,(1d0,0d0),z1r,ntqxx,z2r,ntqxx,(0d0,0d0),w3pi,ntqxx)
C               call zprm('w3pi',w3pi,ntqxx,ntqxx,ntqxx)
C               job 11 a little faster, but harder to parallelize
                call zqsmpy(12,'N','C',ntqxx,ngb,z1r,ntqxx,z2r,ntqxx,(0d0,0d0),w3pi,ntqxx)
C               call zprm('w3pi',w3pi,ntqxx,ntqxx,ntqxx)
                do  itp = 1, ntqxx
                  do itpp = 1, ntqxx
                    w3p(it,itp,itpp) = w3pi(itp,itpp)
                  enddo
                enddo

              enddo



              vcoud(1) = vc
              deallocate(z1r,z2r,w3pi,zmel)
              if (verbose()>=40 .and. kx == iqini) call cputid2(' complete unweighted exchange spectral function. ',0)
C             call zprm('w3p',w3p,(nctot+nbmax),(nctot+nbmax),ntqxx**2)

!         ... Write the Spectrum function for exchange May. 2001
              if (ifexsp/=0) then
                do  it  = 1, nctot+nbmax
                  do itp = 1, ntqxx
                    if (k_comm%ID == 0) write(ifexsp,
     .                "(3i4, 3f12.4,1p,d23.15,d23.15)") ip,itp,it,qrk, ekc(it), -wtt*dble(w3p(it,itp,itp))
                  enddo
                enddo
              endif
              if (debug .and. (k_comm%ID == 0)) then
                do  it = 1,nctot+nbmax
                  do  itp = 1,ntqxx
                    write(6,"(' diagonal w3p =',2i4,2d14.6)") it,itp,w3p(it,itp,itp)
                  enddo
                enddo
              endif

!         ... Correct wfac = weights for valence by esmr
              do  it = nctot+1, nctot+nbmax
                wfac = wfacx(-1d99, ef, ekc(it), esmr) ! gaussian
                do itp = 1, ntqxx
                  w3p(it,itp,1:ntqxx) = wfac * w3p(it,itp,1:ntqxx)
                enddo
              enddo

!ML           This is the partial summation for zsec over kx and ip.
!ML           MPI:   At the end of the loops we need to reduce the results
              do itpp = 1, ntqxx
                do itp = 1, ntqxx !S[j=1,nbloch]  z1p(j,t,n) <B(rk,j) psi(q-rk,n) |psi(q,t')>
                  if (diagonly.and.(itpp/=itp)) cycle
C                 if (itpp>ntqxxd .and.itp/=itpp) cycle ! A test
                  zsec(itp,itpp,ip) = zsec(itp,itpp,ip) - wtt * sum(w3p(:,itp,itpp) )
                enddo
              enddo
              deallocate(w3p)
              endif ! Fast and slow branches

              cycle ! next irot
              if (verbose()>=40) call cputid2(' complete weighted exchange spectral function. ',0)

            endif !exchange

C       --- Correlation part of self-energy ---
!           The correlated part of the self-energy:
!           S[n=all] S[i,j=1,nbloch]
!           <psi(q,t) |psi(q-rk,n) B(rk,i)>
!            < [w'=0,inf] (1/pi) (w-e)/{(w-e)^2 + w'^2} Wc(k,iw')(i,j) >
!                                          <B(rk,j) psi(q-rk,n) |psi(q,t)>
!           e = e(q-rk,n), w' is real, Wc = W-v

!       ... Make matrix elements * overlap (Step C, Remarks)
C           <M Psi2| Psi1> -> <Psi1| Psi2 M >
            call system_clock(c1)
            call ztoccz(zmel,ngb,ngb,nstate*ntqxx)
            call system_clock(c2)
            write(*,'("sxcf ztoccz:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

            nblk =  min(nstate,max(1,((4*2**23)/(ngb*ntqxx))*8)) ! made for ~4GiB and multiple of 8
            write(*,'("sxcf nblk, tmpsz, zmel:",x,i0,",",x,i0," MiB",x,f6.1," GiB")') nblk,
     &         (ngb*ntqxx*nblk)/2**16, real(int(ngb,8)*ntqxx*nstate,8)/2**26

            call system_clock(c1)
            allocate(tmp(ngb,ntqxx,nblk))
            do it = 1, nstate, nblk ! this complication is done to reduce excessive memory usage while not damaging performance
                if (nstate - it + 1 < nblk) nblk = nstate - it + 1 ! last iteration
#ifdef Z3M
                call zgemm3m('c','n',ngb,ntqxx*nblk,ngb,(1d0,0d0),ppovlz,ngb,zmel(1,1,it),ngb,(0d0,0d0),tmp,ngb)
#else
                call zgemm('c','n',ngb,ntqxx*nblk,ngb,(1d0,0d0),ppovlz,ngb,zmel(1,1,it),ngb,(0d0,0d0),tmp,ngb)
#endif
                call zcopy(ngb*ntqxx*nblk,tmp,1,zmel(1,1,it),1)
            end do
            deallocate(tmp)
            call system_clock(c2)
            write(*,'("sxcf orth:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

            if (verbose()>=40) call cputid2(' complete zmel = zzmel * overlap+ .',0)
            if (debug) write(6,"('sum zmel, abszmel=',1p,4d23.16)") sum(zmel),sum(abs(zmel))
C           call zprm('zmel*ppovlz+, vectorized',zmel,ngb,ngb,nstate*ntqxx)

!       --- Contribution to SEc(qt,w) from integration along the imaginary axis ---
!           Loop over w' = (1-x)/x, frequencies in Wc(k,w'),
!           where {x} are gaussian points between (0,1)
            ia = idalloc('zwz0+zwz',iverbose+2,2*nstate*ntqxx,ntqxx*2)
            allocate(zw(nblochpmx,nblochpmx),zwz0(ntqxx,ntqxx,nstate))
            ix = 1 - nwr0 !at omega=0
            nrec=ix

            if (k_comm%ID == 0) read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v
#ifdef USE_MPI
            if (k_comm%size > 1) then
              call MPI_Bcast(zw,nblochpmx*nblochpmx, MPI_DOUBLE_COMPLEX, 0, k_comm%communicator, ierror)
            endif
#endif

            write(*,"('sxcf: ntq; ngb, ntqxx, nstate:',4(x,i8))") ntq, ngb, ntqxx, nstate

C       ... Make zwz0 = zmel*(W(0)-v)*zmel (Step D, Remarks)
            call system_clock(c1)
#ifdef CUD
            call fcuzwz(ngb, ntqxx, nstate, zw, nblochpmx, zmel, ngb, zwz0, ntqxx, 0)
#else
            call matzwz3(2,1,zw,zmel,nblochpmx,ntqxx,nstate,ngb,zwz0) ! zwz0 = zmel*(W(0)-v)*zmel
#endif
            if (debug) write(6,"('sum zwz0, abszwz0=',4d23.14)") sum(zwz0),sum(abs(zwz0))
            call system_clock(c2)
            write(*,'("sxcf zwz0:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

!             if (verbose()>=30 .and. kx == iqini .or. verbose()>=40) call cputid2(' complete z * W(0) * z. ',0)

!!          COH term test ----- The sum of the all states for zwz0 gives the delta function.
C           if (cohtest) then
C             do itp = 1,ntq
C               coh(itp,ip)  = coh(itp,ip) + .5d0*wtt*sum(dreal(zwz0(1:nstate,itp,itp)))
C             enddo
C             deallocate(zw,zwz0,zmel)
C             cycle
C           endif

            if (niw<1) call rx("sxcf:niw <1")
            npm = 1
            if (nwr0/=0) npm = 2 ! TimeReversal off case.

C       ... Frequency integration along the imaginary axis for SEc (Steps E,F Remarks)
            allocate(zwz(ntqxx,ntqxx,nstate))
            do  ix = 1, niw
              nrec = ix
              if (k_comm%ID == 0) then
!               write(0,'("[",I4,"] reading WVI.",I5.5," file at.",I5)') rank,kx,nrec
                read(ifrcwi,rec=nrec) zw ! direct access read Wc(ix) = W(ix) - v
              endif
#ifdef USE_MPI
              if (k_comm%size > 1) then
                call MPI_Bcast(zw,nblochpmx*nblochpmx, MPI_DOUBLE_COMPLEX, 0, k_comm%communicator, ierror)
              endif
#endif

C         ... Make zwz = zmel*(W(ix)-v)*zmel
              call system_clock(c1)
#ifdef CUD
              call fcuzwz(ngb, ntqxx, nstate, zw, nblochpmx, zmel, ngb, zwz, ntqxx, 1)
#else
              call matzwz3(2,npm,zw,zmel,nblochpmx,ntqxx,nstate,ngb,zwz)
#endif
              if (debug) write(6,"('sum zwz(ix), abszwz(ix)=',4d23.14)") sum(zwz(:,:,:)),sum(abs(zwz(:,:,:)))
              call system_clock(c2)
              write(*,'("sxcf zwz",i0,":",x,g20.12)') ix, real(c2 - c1,8)/real(cr,8)



              call system_clock(c1)
              call mkseci(ix,iSigMode,npm,nstate,ntq,ntqxx,nctot,ef,omega,ekc,
     .          GaussSmear(),esmr,niw,freqx,wx,ua_,expa_,wtt,zwz,zwz0,zsec(1,1,ip))
              call system_clock(c2)
              write(*,'("sxcf seci:",x,g20.12)') real(c2 - c1,8)/real(cr,8)
C             call zprm('zsec',zsec(1,1,ip),ntq,ntq,ntq)

            enddo ! ix

            if (verbose()>=20 .and. kx == iqini .or. verbose()>=40) call cputid2(' complete int z * W(omega) * z. ',0)
C       ia = idalloc(' ',iverbose,1,1)

C       Optionally call mkseci after all zwz have been calculated.
C       Then zwz must be allocated zwz(:,:,:,niw)
C       call mkseci(0,iSigMode,npm,nstate,ntq,ntqxx,nctot,ef,omega,ekc,
C    .    GaussSmear(),esmr,niw,freqx,wx,ua_,expa_,wtt,zwz,zwz0,zsec(1,1,ip))

            ia = idalloc('zwz0+zwz',iverbose+4,2*nstate*ntqxx,ntqxx*2)
            deallocate(zwz,zwz0)
            if (debug) write(6,"('sum zsec, abszsec',4d23.16)")sum(zsec(:,:,ip)),sum(abs(zsec(:,:,ip)))

            if (iSigMode==0 .or. iSigMode==4) then !no pole contribution for SE(e_f)
              ia = idalloc('zmel',iverbose+4,ngb*nstate,ntqxx*2)
              deallocate(zmel,zw)
              cycle ! to next rotation
            endif

C       ... Contribution to SEc(qt,w) from the poles of G (real axis) (Step G, Remarks)
C       ... Get ixsmin, ixsmx
            ixsmx = 0; ixsmin = 0
            do  3001 itp = 1, ntqxx
              omg  = omega(itp) !sf 23 May02
              if (omg < ef) then
                itini = 1
                itend = nt0p
              else
                itini = nt0m+1
                itend = nstate
              endif
              do  3011 it = itini, itend
                esmrx = esmr
                if (it<=nctot) esmrx = 0d0
                wfac = wfacx2(omg,ef, ekc(it),esmrx)
                if (GaussSmear()) then
                  if (wfac<wfaccut) cycle ! next it
                  we = .5d0*(weavx2(omg,ef,ekc(it),esmr)-omg)
                else
                  if (wfac==0d0) cycle ! next it
                  if (omg>=ef) we = max(.5d0*(omg-ekc(it)),0d0) ! positive
                  if (omg< ef) we = min(.5d0*(omg-ekc(it)),0d0) ! negative
                endif
                if (it<=nctot) then
                  if (wfac>wfaccut) call rx("sxcf: it<=nctot.and.wfac/=0")
                endif
                do  iwp = 1, nw
                  ixs = iwp
                  if (freqr(iwp)>abs(we)) exit
                enddo
!               This change is because G(omega-omg') W(omg') !may2006
!               if (ixs>ixsmx  .and. omg<=ef ) ixsmx  = ixs
!               if (ixs>ixsmin .and. omg> ef ) ixsmin = ixs
                if (ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
                if (ixs>ixsmin .and. omg< ef ) ixsmin = ixs
                wexx  = we
                if (ixs+1 > nw) then
                  write (*,*) ' nwr0 ixsmin',nwr0, ixsmin
                  write (*,*) ' wexx, dw ',wexx,dw
                  write (*,*) ' omg ekc(it) ef ', omg,ekc(it),ef
                  call rx(' sxcf 222: |w-e| out of range')
                endif
 3011         continue
 3001       continue !end of SEc w and qt -loop

            if (nwr0==0) then
              nwxi = 0
              nwx  = max(ixsmx+1,ixsmin+1)
            else
              nwxi = -ixsmin-1
              nwx  =  ixsmx+1
            endif
            if (nwx > nw) then
              call rx(' sxcf_fal3_sc nwx check : |w-e| > max(w)')
            endif
            if (nwxi < nwr0) then
              call rx(' sxcf_fal3_sc nwxi check: |w-e| > max(w)')
            endif
            if (verbose()>=40 .and. kx == iqini) call cputid2(' complete G-pole part of sigma',0)

C       ... Find nt_max = number of states contributing to pole in omega integration
            nt_max = nt0p !initial nt_max
            do  4001  itp = 1, ntqxx
              omg = omega(itp)
              if (omg > ef) then
                do  it = nt0m+1, nstate ! nt0m corresponds to efm
                  wfac = wfacx2(ef,omg, ekc(it),esmr)
                  if ((GaussSmear().and.wfac>wfaccut) .or. (.not.GaussSmear().and.wfac/=0d0)) then
                    if (it > nt_max) nt_max=it ! nt_max is unocc. state
                  endif         ! that ekc(it>nt_max)-omega > 0
                enddo           ! so it > nt_max does not contribute to omega pole integral
              endif
 4001       continue            !end of  w and qt -loop

C       ... Case Sigma_en =1,3,5
            if (iSigMode == 1 .or. iSigMode == 3 .or. iSigMode==5) then
            else
              call rx('sxcf_fal2_sc: calc contr to sigma from poles of G but iSigMode /= 1 3 5')
            endif

            call system_clock(c1)
            allocate(zw3(ngb,ngb,nwxi:nwx))
            write(*,'("sxcf w3 sz: ",i0," MiB")') size(zw3)/2**16

            do  ix = nwxi, nwx

              nrec= ix-nwr0+1
              if (k_comm%ID == 0) then
!                     write(0,'("[",I4,"] reading WVR.",I5.5," file at ",I5)') rank,kx,nrec
                read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v
              endif
#ifdef USE_MPI
              if (k_comm%size > 1) then
                call MPI_Bcast(zw,nblochpmx*nblochpmx,MPI_DOUBLE_COMPLEX, 0, k_comm%communicator, ierror)
              endif
#endif

C            Make W hermitian, because we need only hermitian part of Sigma_nn'

!$omp parallel do collapse(2)
              do  igb2 = 1, ngb
                do  igb1 = 1, ngb
                  zw3(igb1,igb2,ix) = .5d0*(zw(igb1,igb2)+dconjg(zw(igb2,igb1)) )
                enddo
              enddo
!$omp end parallel do

            enddo               ! ix
            deallocate(zw)
            call system_clock(c2)
            write(*,'("sxcf w read & c+c*:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

            if (test_symmetric_W().and.npm==2) then
              if (onceww(4)) print *,' test_symmetric_W()=',test_symmetric_W(),nwxi,nwx
              call system_clock(c1)
              allocate(zw3x(ngb,ngb))
              do ix = 1,min(abs(nwxi),nwx)
                zw3x = 0.5d0*(zw3(:,:,ix) + zw3(:,:,-ix))
                zw3(:,:, ix) = zw3x
                zw3(:,:,-ix) = zw3x
              enddo
              deallocate(zw3x)
              call system_clock(c2)
              write(*,'("sxcf iw avg:",x,g20.12)') real(c2 - c1,8)/real(cr,8)
            endif

            call system_clock(c1)
C           Add real-axis contribution to zsec(:,:,ip)
            call mksecp(iSigMode,npm,nstate,ntq,ntqxx,nt0p,nt0m,nt_max,
     .        nctot,ngb,ekc,ef,omega,GaussSmear,esmr,nw,wtt,nwxi,nwx,nwr0,
     .        freqr,zw3,zmel,zsec(1,1,ip))
            call system_clock(c2)
            write(*,'("sxcf secp:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

            ia = idalloc('zmel',iverbose+4,ngb*nstate,ntqxx*2)
            deallocate(zw3)
            deallocate(zmel)
            if (verbose()>=30 .and. kx == iqini .or. verbose()>=40) call cputid2(' complete zsec',0)

            if (debug) then
              print *,' end of do  2001 or 2002 '
              do  itp = 1, ntqxx
                write(6,'(" zsec=",i3,2d15.7)') itp,zsec(itp,itp,ip)
              enddo
            endif
            if (verbose()>=40) call cputid2(' complete irot',0)
 1000     continue              ! end do irot.  Was not parallelized, no reduction has to be performed for it.

          if (k_comm%ID == 0) ifvcoud = iclose('Vcoud.'//charnum5(kx))
          if (.not.exchange) then
            if (k_comm%ID == 0) then
              ifrcw  = iclose('WVR.'//charnum5(kx))
              ifrcwi = iclose('WVI.'//charnum5(kx))
            endif
          endif
          if (k_comm%ID == 0) then
            if (verbose()>=30) call cputid2(' end of kx loop',0)
          endif

 1100   continue                ! end of kx-loop

!ML   MPI reduction of zsec:
!ML
!ML      complex(8) :: zsec(ntq,ntq,nq)
!ML                                 /\
!ML                    ip  _________||
!ML
!ML      for each ip we need to sum all contributions of kx.
!ML      each k_comm group has already summed over it set of kx
!ML      ==> we need to reduce from all k-groups: use k_inter!
!ML
!ML      each q-group only knows its own ip sections of zsec
!ML      as zsec was set to 0 before calling sxcf_fal3, we can combine then
!ML      by an ALL_REDUCE over the q_inter communicators.
!ML
!ML   Possible optimization: don't stride through kx, but use blocks, so that data
!ML   can be more easily distributed.
!ML
#ifdef USE_MPI
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, zsec(:,:,ip), ntq*ntq,
     &        MPI_DOUBLE_COMPLEX, MPI_SUM,
     &        k_inter%communicator, ierror)
#endif

        if (verbose()>=10) call cputid2(' exit 1100 loop',0)
        ifvcoud =iclose('Vcoud')
        if (irot==1) write(6,"('  sum(abs(zsec))=',d23.15)") sum(abs(zsec))
ccccccccccccccccccccccccccccc
c     print *," end kount ip irot=", irot, ip, sum(kount(:,ip))
c     print *," sumkount 1=",irot,ip,sum(kount(:,:))
cccccccccccccccccccccccccccc
        if (allocated(vcoul))deallocate(vcoul)
 1001 continue                  ! end do ip

!      write(0,'("[",I4,"] ntq = ",I4)') rank,ntq
!
!      do ip=1,nq
!         do i=1,ntq
!            write(0,"('[',I4,'] BEFORE ALLREDUCE: ip isp i=' i4,i2,i3,'  Re(Sc) =',f12.4,'  Img(Sc) =',f12.4 )") !sf 21May02
!     &           rank, ip,isp,i,
!     &           dreal(zsec(i,i,ip)), !sf 21May02
!     &           dimag(zsec(i,i,ip)) !,zfac,zfac1,zfac2  !sf 21May02
!         enddo
!      enddo


#ifdef USE_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, zsec, ntq*ntq*nq,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM,
     &     q_inter%communicator, ierror)


      CALL MPI_ALLREDUCE(MPI_IN_PLACE, nbandmx, nq,
     &     MPI_INTEGER, MPI_SUM,
     &     q_inter%communicator, ierror)

#endif

!ML debug output
!      if (rank == 0) then
!
!         write(0,'("[",I4,"] After Reduce: ntq = ",I4)') rank,ntq
!
!         do ip=1,nq
!            do i=1,ntq
!               write(0,"('[',I4,'] AFTER ALLREDUCE: ip isp i=' i4,i2,i3,'  Re(Sc) =',f12.4,'  Img(Sc) =',f12.4 )") !sf 21May02
!     &              rank, ip,isp, i,
!     &              dreal(zsec(i,i,ip)), !sf 21May02
!     &              dimag(zsec(i,i,ip)) !,zfac,zfac1,zfac2  !sf 21May02
!            enddo
!         enddo
!
!      endif

      if (verbose()>=40) then
        itp = idalloc(' ',11,1,1)
      else
        itp = idalloc(' ',1,1,1)
        write(*,"(/' Max dynamic allocation sxcf_fal3_scz:',i7,' MB'/)") itp
      endif

      end
