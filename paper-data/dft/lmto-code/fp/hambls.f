      subroutine hambls(mode,s_site,s_spec,s_lat,s_ham,nbas,isp,ispc,q,
     .  k1,k2,k3,qsig,nqsig,smpot,vconst,lcplxp,lso,alfa,ndimh,h,s,dsdk,ismidb,nevl)
C- Hamiltonian and overlap of LDA+(sigma-vxc) for one k-point.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos sighh sighk sigkk pihh pihk pikk sighhx
Ci                 sighkx sigkkx pihhx pihkx pikkx
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hambl augmbl bstrux smhsbl hsibq hsubblock hsibq2
Cio                hsibq4
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt rsma pz name orbp ngcut
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hambl augmbl bstrux uspecb smhsbl hsibq tbhsi
Cio                hsubblock hsibq2 hsibq4
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat alat qlat vol awald tol nkd nkq gmax nabc ng kv
Ci                 kv2
Co     Stored:     *
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:cg indxcg jcg cy qlv dlv gv igv kv igv2
Cio    Passed to:  hambl augmbl bstrux hxpbl ghibl hklbl gklbl hxpgbl
Cio                ghigbl hklgbl smhsbl hhibl phhibl hsmbl hsibq sugvec
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  oveps ovncut ndhrs nprs sigp eseavr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb iaxs hrs
Cio    Passed to:  hambl sopert3
Ci Inputs
Ci   mode  :On input:
Ci         :1s digit specifies LDA matrix elements
Ci         :  0 compute LDA hamiltonian and overlap
Ci         :    If the self-energy is not to be added, return with h,s
Ci         :  1 Assume h,s already computed.
Ci         :    Self-energy added to h, if 10s digit set
Ci         :  2 Like 0, but also diagonalize h,s and return with mode=-1
Ci         :    In this h,s contain evals, evecs, respectively; see Remarks.
Ci         :    Note: if 10s digit mode is nonzero, this switch is implied.
Ci         :    Switches 1+2 may be combined.
Ci         :  4 Add 4 to 1s digit to return <i| ds/dk |j>
Ci         :    Can add 4 to 1 or to 2 (NEED CHECK, also for noncollinear case)!
Ci         :10s digit specifies self-energy matrix elements
Ci         :  0 Self-energy is not added
Ci         :  1 Add Bloch sum of r.s. sigma
Ci         :  2 Same as 1
Ci         :  4 If this bit is present, Assume r.s. sigma is real
Ci         :100s digit specifies interpolation mode
Ci         :  0 Simple bloch sum of sigma
Ci         :  1 approx high- and low-energy sigma with diagonal in LDA basis;
Ci         :    see Remarks
Ci         :  2 perturbation approach, diagonalizing exactly only
Ci         :    the first kappa block
Ci         :  3 Linear interpolation of sigma from neighboring points
Ci         :    sig(q) = sum_i wt_i A+_i sig(qpi) A_i
Ci         :    where A_i = z^LDA_qi (z_q^LDA)^-1
Ci         :1000s digit nonzero flags quantities other than h and s are
Ci         :      to be returned in arrays h and s
Cl         :    0 no special treatment; h and s are hamiltonian, overlap
Cl         :      For the following, mode is returned as -1
Ci         :    1 evecs of the LDA hamiltonian in array s;
Ci         :      evals of the LDA hamiltonian in array h, first column
Ci         :      sigma not calculated.
Ci         :    2 return raw sigma, orbital basis, in s and
Ci         :      1st col of array h, LDA evals (real vector)
Ci         :    3 return in array s, sigma in LDA basis and
Ci         :      1st col of array h, LDA evals (real vector)
Ci         :    4 return in array s, modified sigma in LDA basis and
Ci         :      1st col of array h, LDA evals (real vector)
Ci         :      (Same as 1000s digit mode 3, but
Ci         :      low- and high- energy blocks replaced by diagonal parts)
Ci         :    5 return in array s, modified sigma in orbital basis
Ci         :      (Same as 1000s digit mode 4, but sigma is in orbital basis)
Ci         :    6 Returns evecs of full hamiltonian in array s,
Ci         :      evals of the full hamiltonian in array h (NOT IMPLEMENTED)
Ci         :(Digit has a special meaning unique to interpolation mode 3:)
Ci         :    number of interpolation points
Ci         :    0 -> use default (4 points)
Co         :On output:
Co         :mode = -1 => hamiltonian was diagonalized internally
Ci   nbas  :size of basis
Ci   isp   :spin index
Ci   ispc  :index where to store h,s. (for noncollinear mode)
Ci   q     :Bloch vector (k-point)
Ci   k1,k2,k3 dimensions of smpot
Ci   qsig  :q-points at which sigma can be computed
Ci   nqsig :number of qsig
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   vconst:additional constant potential
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   osig,otau,oppi  augmentation matrices
Ci   alfa  :add alfa * overlap to hamiltonian
Ci         :This is for stability in evals.  Preferably alfa=0
Ci   ndimh :dimension of hamiltonian and overlap h,s
Co Outputs
Co   mode  :if mode returns <0, hambls finds (approximate)
Co         :evals and evecs internally; see also description of h,s
Co   h     :Hamiltonian matrix, unless 1s digit or 1000s digit mode is nonzero.
Co         :In that case mode is returned as -1, and
Co         :what is returned in h depends on 1000s digit mode, e.g.
Co         :if it is 1, h is diagonalized internally, and on exit
Co         :h is returned as a diagonal matrix, whose elements are eigenvalues;
Co         :s contains eigenvectors.  The method for doing this, and the
Co         :motivation for it, are described in Remarks.
Co   s     :overlap matrix, unless 1s digit or 1000s digit mode is nonzero.
Co         :See description for h.
Co  ismidb :upper limit to middle block of sigma
Co  nevl   :(If mode returned -1): dimension of reduced h
Cl Local variables
Cl   lwsig :is 1000s digit of mode.  See Remarks for further description
Cl         :0 no special treatment
Cl         :1 hambls does not generate hamiltonian, but returns
Cl         :  LDA eigenvectors in place of s
Cl         :  LDA eigenvalues in the first column of h
Cl         :  In this case, sigm is not generated or used
Cl         :2 hambls does not generate hamiltonian, but returns
Cl         :  sigm in orbital basis in place of s
Cl         :  LDA eigenvalues in the first column of h
Cl         :3 hambls does not generate hamiltonian, but returns
Cl         :  sigm rotated to lda basis in place of h
Cl         :  LDA eigenvalues in the first column of h
Cl         :4 Same as lwsig=3, but
Cl         :  low- and high- energy blocks replaced by diagonal parts
Cl         :5 sigm returned in orbital basis, but with
Cl         :  low- and high- energy blocks replaced by diagonal parts
Cl         :  Same as lwsig=4, but sigm is returned in orbital basis
Cl   mxitrp:maximum no qp to include in interpolation for self-energy
Cl   sigp  :parameters for approximating self-energy sigma.  sigma
Cl         :is approximated by its diagonal part sigii for energies
Cl         :below a low-energy cutoff (specified nmin or emin) and
Cl         :above a low-energy cutoff (specified nmax or emax).
Cl         : arg 1: (mode) specifies how to set diagonal part sigii
Cl         :        for states above the high-energy cutoff nmax or emax
Cl         :        0 constrain sigii to be > asig+bsig*e
Cl         :        1 constrain sigii to be = asig+bsig*e
Cl         :        2 constrain sigii to be > asig and < bsig
Cl         :        3 constraint same as mode 1.
Cl         :          Mode 3 differs in that the least-squares fit to
Cl         :          sigii (for informational purposes only, to help
Cl         :          estimate asig and bsig) is done for states between
Cl         :          efit and nmax or emax
Cl         :        4 constrain sigii to be asig
Cl         : arg 2: nmin : nmin=0 => not used
Cl         :               nmin<0 => use emin in place of nmin
Cl         :               Usage depends on whether mode is 0 or 1
Cl         :               mode = 0: for states 1..nmin, off-diagonal
Cl         :               parts of sig(1:nmin,1:nmin) are zeroed out.
Cl         :               mode = 1: sigd(1..nmin) is filled with emin
Cl         : arg 3: emin : (used only if nmin<0)
Cl         :               Usage depends on whether mode is 0 or 1
Cl         :               mode = 0: for states e_i<emin, off-diagonal
Cl         :               parts of sig(1:i,1:i) are zeroed out.
Cl         :               mode = 1: sigd(1..nmin) is filled with emin
Cl         : arg 4: nmax : sigma for levels i>nmax are approximated by
Cl         :               sigii AND constrained according to mode, arg 1
Cl         : arg 5: emax : (used only if nmax<=0)
Cl         :             : sigma for levels e<emax are approximated by
Cl         :               sigii AND constrained according to mode, arg 1
Cl         : arg 6: asig : constraint used to approximate
Cl         :               sigii = asig + E * bsig  or
Cl         :               asig < sigii < bsig
Cl         : arg 7: bsig : constraint used to approximate
Cl         :               sigii = asig + E * bsig  or
Cl         :               asig < sigii < bsig
Cl         : arg 8: efit : fit sigii between efit and emax (mode 3)
Cl   napw  :number of PWs in APW part of basis
Cl   igapwl:PWs in units of reciprocal lattice vectors,
Cl         :possibly modified if q is shortened.
Cr Remarks
Cr   This routine can add Vxc(QSGW)-Vxc(LDA) to the LDA hamiltonian.  Vxc(QSGW)
Cr   is calculated on a fixed k-mesh.  To interpolate it to an arbitrary k is
Cr   challenging, particularly, because of problems with interpolation of the
Cr   high-energy parts of sigma.
Cr
Cr   Interpolation of the low-energy parts can be accomplished with a set of
Cr   short-ranged basis functions, as follows:
Cr
Cr   (1) Project matrix elements of sigm in the low-energy subspace onto the
Cr       short-ranged functions, approximate the high-energy part with
Cr       a constant, diagonal matrix as described below.
Cr       This is done at the end of the GW cycle; see hqpe.sc.m.F
Cr   (2) Make an inverse Bloch transform of this subspace into real space;
Cr       This is done in rdsigm.f.
Cr   (3) Evaluate sigm at arbitrary k in the subspace by a forward Bloch sum;
Cr   (4) Project back to the full hilbert space.
Cr       Steps (3) and (4) are accomplished here.
Cr
Cr   The sm-Hankels provide a short-range basis of functions, which also
Cr   well describe eigenfunctions in the low-energy subspace.  We use them
Cr   since they are the basis (or some portion in the PMT case) already.
Cr
Cr   1. Project out the low-energy subspace.  At the end of the GW step,
Cr   sigm is rotated from the orbital basis to the eigenfunction basis
Cr   (where H is diagonal).  sigm is ordered by eigenvalue, and partitioned:
Cr                 ( low-low    low-high )
Cr                 (                     )
Cr                 ( high-low   high-high)
Cr   The low-high and high-low parts are set to zero, and the high-high block is
Cr   replaced by a constant matrix (some weighted average of the diagonal
Cr   elements; see eseavr in hqpe.sc.m.F).  This is a good approximation because
Cr     (a) off-diagonal parts of sigm are generally small
Cr     (b) the high-energy parts only weakly affect the low-energy
Cr         parts (as is evident from perturbation theory)
Cr     (c) the high-energy parts are generally in the regime (>2Ry) where
Cr         states are strongly delocalized, and typically free-electron like.
Cr         Sigma is reasonably like a constant, diagonal matrix in that region.
Cr
Cr   3.  Obtain sigm(k) in the sm-Hankel subspace of basis funcitions
Cr       by doing a Bloch sum.  This is readily done; see step B below.
Cr
Cr   4.  Project back to the full hilbert space.  Step 1 is repeated, but this
Cr       time rotated to the basis of LDA eigenfunctions (which are known)
Cr       rather than the QSGW eigenfunctions.  The high-high part is replaced
Cr       with a constant, diagonal array (global constant eseavr); the resulting
Cr       sigm has the rank of the full hilbert space and rotated back to the
Cr       orbital basis.  This step accomplishes the desired results, while
Cr       surmounting two complications:
Cr
Cr           * eigenfunctions of H with sigm included are not known at
Cr             `arbitrary k, precluding the rotation into the eigenfunction
Cr             basis of H.  But LDA and QSGW subspace of eigenfunctions
Cr             between the bottom and 2 or 3 Ry above Ef are very similar,
Cr             so this approximation essentially no better or worse than
Cr             the original projection onto the sm-Hankel subspace.
Cr
Cr           * This interpolation procedure work well only if the sm-Hankel
Cr             functions are very short-ranged.  sm-Hankels for typical EH
Cr             (EH=-0.4) appear to be capable of interpolating sigma to about
Cr             2Ry or so, but no higher.  Thus the low-energy part of hilbert
Cr             space must be truncated to an energy somewhat lower than the rank
Cr             of the sm-Hankel basis would otherwise permit.
Cr
Cr   In detail, the procedure outlined aboved is accomplished as follows.
Cr
Cr   A. generate evals and evecs of the LDA hamiltonian
Cr
Cr   B. Make a Bloch sum sigm(k) of sigm(RS) in the LMTO subspace
Cr
Cr   C. Rotate (sigm-vxc)(k) to LDA basis:
Cr      sigm(LDA) = (Zlda)+ sigma-vxc (Zlda)
Cr      The rank of sigm(LDA) is the rank of the LMTO basis.
Cr
Cr   D. Overwrite sigm(k) with a diagonal array for high-lying
Cr      (and optionally low-lying) energies.  The untouched lower part
Cr      is typically 2-3 Ry or less, typically with significantly fewer
Cr      states than than the rank of the LMTO basis.  At all events it cannot
Cr      be larger than the LMTO basis.    The rank of sigm, included the
Cr      higher parts, is the rank of the full basis. This is the crucial
Cr      step that enables interpolation to proceed smoothly.
Cr
Cr   E. Add Hlda = Elda to sigm generate total hamiltonian in LDA basis
Cr
Cr   F. Diagonalize Hlda+sigm, LDA basis
Cr
Cr   G. Rotate evecs from LDA basis to orbital basis
Cr
Cr   Other quantities may be sought, depending on lswig.
Cr   lwsig:
Cr     1  Return LDA evecs in s, evals in h (stop at step A above).
Cr        sigm is not used in this case.
Cr
Cr     2  Return with unmodified sigma (orbital basis) in s,
Cr        (stop at step B above) and LDA evals in h
Cr
Cr     3  Return with unmodified sigma (lda basis) in s and LDA evals in h
Cr        (stop at step C above)
Cr
Cr     4  Return with modified sigma (lda basis) in s and LDA evals in h
Cr        (stop at step D above)
Cr
Cr     5  Return with modified sigma (orbital basis) in s and LDA evals in h
Cr        (after step D above rotate sig to orbial basis)
Cu Updates
Cu   11 Aug 14 Adapted for PMT method, interpolating using LMTO basis
Cu   08 Feb 13 Internally shortens q vector
Cu   30 Aug 12 Modifications to handle SO=3
Cu   10 Nov 11 Begin migration to f90 structures
Cu   17 Mar 11 Incorporate eseavr
Cu   30 Jan 10 Handles epsovl; new arg list
Cu   04 Jul 08 (T. Kotani) New PW addition to basis
Cu   27 Jun 08 Redesigned lwsig modes
Cu   26 May 07 Some preparation for rotation betw/ LDA, GW basis
Cu   16 Jan 07 Redefined functions of lwsig
Cu   20 Feb 05 Add switch (1000s digit mode) to return sigma in h
Cu    7 Feb 07 SO coupling made to with with self-energy
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   25 Jun 03 Update to incorporate sigp mode 3
Cu   24 May 03 New interpolation modes
Cu   10 Jan 03 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,isp,ispc,k1,k2,k3,lcplxp,ndimh,nqsig,lso,ismidb
      double precision qsig(3,nqsig),alfa
      double precision q(3),vconst
      double complex smpot(k1,k2,k3)
      double complex h(ndimh,ndimh,*),s(ndimh,ndimh,*),dsdk(ndimh,ndimh,3,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Dynamically allocated arrays
      integer,pointer :: igapwl(:,:)
      integer,allocatable:: ipm(:)
      real(8),pointer :: el(:)
      real(8),allocatable,target :: e(:),e0(:)
      complex(8),allocatable:: sigq(:,:)
      complex(8),allocatable, target:: zqi(:,:),z0(:,:)
      complex(8),pointer :: zl(:,:)
C     complex(8),allocatable:: hqi(:,:,:),sqi(:,:,:)
C ... Local parameters
      character outs*256
      logical ldsdk
      double precision plat(3,3),qlat(3,3),wk(ndimh*11),qs(3),
     .  swt,emin,emax,asig,bsig,xx,sigp(10),eseavr(2),epsovl
      equivalence (emin,sigp(3)),(emax,sigp(5))
      equivalence (asig,sigp(6)),(bsig,sigp(7))
      integer i,j,lrsig,ipr,ndim1,hreal,nl,mxorb,nttabs,nsp,ndhrs,
     .  mode0,moditp,n123(4),off1,nspc,ndimhx,lwsig,ovncut,nevl,napw,nlmto
      logical, parameter :: T=.true., F=.false.
      integer, parameter :: NULLI=-99999
      real(8), parameter :: tol=1d-7
      complex(8), parameter :: zer=(0d0,0d0),one=(1d0,0d0)
C     integer iqfbz(8)
C     integer iprm(nqsig),ipiv(ndimh),nev,nsim,mxitrp,evl(ndimh)
C     double precision qloc(3,nqsig),wt(4),ddot,norm(4,4),qpi(3,4),dq,dinv
      procedure(logical) :: cmdopt
      procedure(real(8)) :: dlength
      procedure(integer) :: nglob

C     q(1) = -0.17500d0; q(2) =  0.17500d0; q(3) =  0.17500d0

C     print *, '!!y shift q'; q(1) = q(1)-.001d0; print *, sngl(q)

C --- Setup ---
      call tcn('hambls')
      mode0 = mod(mode,10); lrsig = mod(mode/10,10); moditp = mod(mode/100,10); lwsig = mod(mode/1000,10)
      ldsdk = mode0 >= 4; mode0 = mod(mode0,4)
      napw = s_lat%napw         ! Number of APW's for this k
      nlmto = s_ham%nlmto       ! Number of LMTO's ... independent of k
      nspc = nglob('nspc')
      ndimhx = ndimh * nspc

C     print *, '!!y debugging turn off ldsdk'; ldsdk = .false.

C     See above for mode0=2
      if (lrsig /= 0 .and. mode0 < 2) mode0 = mode0+2
      if (lwsig /= 0 .and. mode0 >= 2) mode0 = mode0-2
      ismidb = 0
      call sanrg(.true.,lwsig,0,5,'hambls:','lwsig')

C ... Shorten q; shift APW G vectors to correspond
      igapwl => s_lat%igv2
      call shorps(1,s_lat%qlat,(/72,2,2/),q,qs)
      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        allocate(igapwl(3,napw))
        call shorigv(q,qs,s_lat%plat,napw,s_lat%igv2,igapwl)
      endif

C --- LDA hamiltonian ---
      if (mod(mode0,2) == 0) then
        j = 0
        if (lso == 1 .or. lso == 3) j = 10
        call hambl(j,s_site,s_spec,s_lat,s_ham,isp,qs,k1,k2,k3,smpot,vconst,lcplxp,alfa,
     .    ndimh,napw,igapwl,h(1,1,ispc),s(1,1,ispc),h(1,1,3))
C        call zprm('hamiltonian, LDA',12,h,ndimh,ndimh,ndimh)
C        call zprm('overlap',12,s,ndimh,ndimh,ndimh)
C       nevl = ndimh
      endif
C     call zprm('hso',2,h(1,1,3),ndimh,ndimh,ndimh)

C ... Store LDA eigenvectors, eigenvalues in s and h and exit
      if (lwsig == 1) then

        allocate(zqi(ndimh,ndimh))
        allocate(z0(ndimh,ndimh))
        allocate(e(ndimh))

        call zhevx(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),1,T,ndimh,9d9,i,z0,F,e,ndimh,zqi)
C       call zprm('evecs, LDA',2,zqi,ndimh,ndimh,ndimh)
C       i = fopna('out',-1,0)
C       call ywrm(0,'lda evecs',2,i,'(9f20.10)',zqi,1,ndimh,ndimh,ndimh)
C       call ywrm(0,'lda eval',1,i,'(9f20.10)',e,1,ndimh,ndimh,1)
C       call rx0('done')

        call dcopy(ndimh**2*2,zqi,1,s(1,1,ispc),1)
        call dcopy(ndimh,e,1,h,1)

        deallocate(z0,zqi,e)
        goto 999
      endif

C ... Exit if no further elaborations: add self-energy and/or diagonalize
      if (mode0 < 2 .and. lwsig == 0) goto 999

C --- Setup for Bloch transform ---
C     Case skip parts having to do with sigma
      if (lrsig == 0) then
        call dcopy(ndimh,e,1,e0,1)
        call zcopy(ndimh**2,zqi,1,s(1,1,ispc),1)
        goto 98
      endif
C     if (napw > 0) call rx('sigm not implemented for PWs')
      ndhrs = s_ham%ndhrs
C     hreal=1 -> sig(T) is real
      hreal = 0; if (lrsig >= 4) hreal = 1
      nl = nglob('nl'); mxorb = nglob('mxorb'); nsp  = nglob('nsp')
      call sanrg(.true.,isp,1,nsp,'hambl:','isp')
      nttabs = s_ham%nprs(nbas+1)
      plat = s_lat%plat
      call dinv33(plat,1,qlat,swt)
      sigp = s_ham%sigp; eseavr = s_ham%eseavr
      if (sigp(1) == 4) then
        asig = eseavr(isp); bsig = 0
        call rxx(asig == NULLI,
     .    'hambls: average SE not defined ... cannot use SIGP_MODE=4')
      endif

C --- Standard interpolation of sigma ---
      if (moditp == 0 .and. lwsig == 0) then
C           iprmb   transp   add    complex s    c*16
        i = 100000 + 4000  + 100 + 40*(1-hreal) + 10
        call bloch(i,qs,nl,plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,
     .    s_ham%hrs,ndhrs,isp,nsp,nlmto,nlmto,0,ndimh,0,ndimh,0,
     .    h(1,1,ispc),xx,xx)
C       call zprm('h+sigm, orb bas',2,h(1,1,ispc),ndimh,nlmto,nlmto)
        goto 999
      endif

C --- Diagonalize LDA hamiltonian ---
      allocate(zqi(ndimh,ndimh),z0(ndimh,ndimh))
      allocate(e(ndimh),e0(ndimh))
      zl => zqi; el => e

C ... PMT method ... require e,z from LMTO part of basis only
      if (s_ham%nlmto /= ndimh .and. lrsig /= 0) then
        allocate(zl(ndimh,ndimh),sigq(ndimh,ndimh),el(ndimh))
        call zcopy(ndimh**2,h(1,1,ispc),1,zqi,1)
        call zcopy(ndimh**2,s(1,1,ispc),1,sigq,1)
        call zhevx(nlmto,ndimh,zqi,sigq,1,T,nlmto,9d9,nevl,z0,F,el,ndimh,zl)
        deallocate(sigq)
C       call zprm('evecl',2,zl,ndimh,nlmto,nlmto)
      endif

C ... zqi <- evecs, evals of LDA hamiltonian = Z.lda
C     This is step A in Remarks
C     call zprm('LDA hamiltonian, LDA',12,h,ndimh,ndimh,ndimh)
C     call zprm('overlap',12,s,ndimh,ndimh,ndimh)
      epsovl = s_ham%oveps
      ovncut = s_ham%ovncut
      if (epsovl == 0 .and. ovncut == 0) then
        call zhevx(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),1,T,ndimh,9d9,nevl,z0,F,e,ndimh,zqi)
      else
        nevl = -1
        e0(1) = 99999d0
        call zhevo(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),ndimh,9d9,epsovl,ovncut,nevl,i,e,e0,ndimh,zqi)
      endif
!     call zprm('evecs, LDA',2,zqi,ndimh,ndimh,ndimh)
!     call prmx('evals, LDA',e,ndimh,ndimh,1)
C     i = fopna('out',-1,0)
C     call ywrm(0,'lda evecs',2,i,'(9f20.10)',zqi,1,ndimh,ndimh,ndimh)
C     call ywrm(0,'lda eval',1,i,'(9f20.10)',e,1,ndimh,ndimh,1)
C     call rx0('done')

C --- Interpolate sigma keeping higher energy LDA eigenvectors ---
      if (moditp == 1) then

C       In this branch, sigq is just a work array
        allocate(sigq(ndimh,ndimh))

C   ... Bloch transform sigm(RS)-sigm(k).  This is step B in Remarks
C       Call to bloch : 104110=perm orb, transpose, no add, c*16
C           iprmb   transp  add    complex s    c*16
        i = 100000 + 4000 + 000 + 40*(1-hreal) + 10
C       i = i + 3000000 +[0d0,0d0,-.001d0]
        call bloch(i,qs,nl,plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,s_ham%hrs,
     .    ndhrs,isp,nsp,nlmto,nlmto,0,ndimh,0,ndimh,0,s(1,1,ispc),xx,xx)
C        if (cmdopt('--zhblock',9,0,outs)) then
C          call zhblock(outs(10:),ndimh,110,s(1,1,ispc),s)
C        endif
!       call zprm('raw sigm(orb)',2,s(1,1,ispc),ndimh,nlmto,nlmto)

C       print *, '!!y debugging ldsdk zero out sigma'; s(:,:,ispc)=0

C        i = fopna('out',-1,0)
C        print *, 'hambls',ndhrs,qs
C        call ywrm(0,'sig(k)',3,i,'(9f20.10)',s(1,1,ispc),1,ndimh,nlmto,
C     .    nlmto)
C        call rx0('done')
        if (ldsdk) then
          do  j  = 1, 3
            call bloch(i+j*1000000,qs,nl,plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,s_ham%hrs,
     .        ndhrs,isp,nsp,nlmto,nlmto,0,ndimh,0,ndimh,0,dsdk(1,1,j,ispc),xx,xx)
C           print *, j; call zprm('raw dsdk, orb bas',2,dsdk(1,1,j,ispc),ndimh,nlmto,nlmto)
          enddo
        endif

C       LDA evals in h
        if (lwsig >= 2 .and. lwsig <= 4) then
          call dcopy(ndimh,e,1,h,1)
        endif

C       Return with raw sigm (orbital basis) in s, LDA evals in h
        if (lwsig == 2) goto 99

C   ... Rotate sigm to LDA basis; step C in Remarks.
C       s <- (Z.lda)+ sigma-vxc Z.lda = sigma-vxc in LDA basis.
C       call zprm('LDA evecs',2,z0,ndimh,nlmto,nlmto)
        call phmbls(3,ndimh,nlmto,e0,xx,z0,s(1,1,ispc),zl,zl,s(1,1,ispc))
C       call zprm('sigma, LDA basis',2,s,ndimh,nlmto,nlmto)
C       Rotate corresponding ds/dk
        if (ldsdk) then
          do  j  = 1, 3
            call phmbls(1,ndimh,nlmto,e0,xx,z0,dsdk(1,1,j,ispc),zl,zl,dsdk(1,1,j,ispc))
C           print *, j; call zprm('raw dsdk, lda basis',2,dsdk(1,1,j,ispc),ndimh,nlmto,nlmto)
          enddo
        endif

C       Return with raw sigm (LDA basis) in s, LDA evals in h
        if (lwsig == 3) then
          if (.not. associated(zl,zqi)) deallocate(el,zl)
          goto 99
        endif

C   ... Overwrite z+ sigma z with approx for low,high-lying states
C       This is step D in Remarks
        do  i = nlmto+1, ndimh
          el(i) = el(nlmto)
        enddo
C       Zero out dsig/dk for these states
        if (ldsdk) then
          wk(1:7) = sigp(1:7); sigp(1) = 4; asig = 0
          do  j  = 1, 3
            call psigd(0,ndimh,dsdk(1,1,j,ispc),el,sigp,n123,[xx])
C           print *, j; call zprm('approx dsdk, LDA basis',2,dsdk(1,1,j,ispc),ndimh,nlmto,nlmto)
          enddo
          sigp(1:7) = wk(1:7)
        endif
        call psigd(0,ndimh,s(1,1,ispc),el,sigp,n123,wk)
        ismidb = n123(3)
C       call zprm('approx sigma, LDA basis',2,s,ndimh,ndimh,ndimh)

        if (.not. associated(zl,zqi)) deallocate(el,zl)

C       Return with modified sigm (LDA basis) in s, LDA evals in h
        if (lwsig == 4) goto 99

C   ... Return with modified sigm (orb basis) in s, evals in h
        if (lwsig == 5) then
          if (ldsdk) call rx('hambls not ready for lwsig=5 and dsdk')
!         call zprm('modified sigm(lda)',2,s(1,1,ispc),ndimh,ndimh,ndimh)
C         Rotate from LDA to orbital basis
          call phmbls(33,ndimh,ndimh,e0,xx,z0,s(1,1,ispc),zqi,zqi,s(1,1,ispc))
C          print *, '!!1'
!          call zprm('modified sigm(orb)',2,s(1,1,ispc),ndimh,ndimh,ndimh)
          call dcopy(ndimh,e,1,h,1)
          goto 99
        endif

C   ... s <- Z.lda+ (h(lda) + sigma) Z.lda = e(lda) + Z.lda+ sigma Z.lda
C       This is h+sigma(approx) in LDA basis, step E in Remarks
C       Initialize e0 with diagonal term of h+sigma_approx in LDA basis
C       For low and high-lying states, e0 is already the final answer
C       The middle eigenvalues e0 will be updated later.
C       dsdk is not modified in this step.
        do  i = 1, ndimh
          s(i,i,ispc) = s(i,i,ispc) + e(i)
          e0(i)  = s(i,i,ispc)
        enddo

C   ... s <- evecs of h + sig in LDA basis, for states in middle block
C       This is step F described in Remarks
C       Middle block has range (off1+1:off1+ndim1,off1+1:off1+ndim1)
C       dsdk is not modified in this step.
        off1  = n123(2)
        ndim1 = n123(3)-n123(2)
        call zhevx(ndim1,ndimh,s(1+off1,1+off1,ispc),xx,0,T,ndim1,9d9,i,
     .    sigq,F,e0(1+off1),ndimh,z0)
C       call zprm('evecs of h+sig, LDA basis',2,z0,ndimh,ndim1,ndim1)
C       call prmx('evals',e0,ndimh,ndimh,1)

C   ... Rotate dsdk (LDA basis) to dsdk (eigenfunction basis)
        if (ldsdk) then
          call dpzero(sigq,2*size(sigq))  ! should not be necessary
          do  j  = 1, 3
            call phmbls(1,ndimh,ndim1,e0,xx,sigq,dsdk(1+off1,1+off1,j,ispc),z0,z0,dsdk(1+off1,1+off1,j,ispc))
!           print *, j; call zprm('approx dsdk, e basis',2,dsdk(1,1,j,ispc),ndimh,nlmto,nlmto)
          enddo
        endif

C   ... Rotate LDA evecs(orbital basis) to evecs of h+sig, orbital basis
C       This is step G described in Remarks
C       Given zqi = LDA evecs (orbital basis)
C       and   z0  = h+sig evecs, LDA basis
C       Then elda = zqi+  Hlda zqi  =>  e = z0+ zqi+ (Hlda+sig) zqi z0
C       Initialize with entire Z.lda; overwrite evecs off1+1..off1+ndim1
        call zcopy(ndimh**2,zqi,1,s(1,1,ispc),1)
        call zgemm('N','N',ndimh,ndim1,ndim1,one,zqi(1,off1+1),
     .    ndimh,z0,ndimh,zer,s(1,off1+1,ispc),ndimh)
C       call zprm('evecs of h+sig+LzSz, orbital basis',2,s,ndimh,ndimh,ndimh)

C   ... Reorder evecs,evals in ascending eval order
        allocate(ipm(ndimh))
        if (ldsdk) then
          do  j  = 1, 3
            call phmbls(16,ndimh,ndimh,e0,ipm,sigq,xx,z0,z0,dsdk(1,1,j,ispc))
          enddo
        endif
        call phmbls(8+4,ndimh,ndimh,e0,ipm,sigq,xx,s(1,1,ispc),s(1,1,ispc),xx)
        deallocate(ipm)

C --- Interpolate sigma via evecs h+sig, 1st LMTO kappa block ---
C     Define     z^1 = evecs of 1st LMTO kappa block
C     subscript _L -> LDA evec ; subscript S -> GW evec
C     This branch estimates z = z_L z_L^1 S^1 Z^1_S
C     where S^1 is overlap matrix of 1st LMTO kappa block.
      elseif (moditp == 2) then

        call rx('interp mode not fully checked')

CC#ifdef F90
C        allocate(hqi(ndimh,ndimh,1))
CC       In this branch, sigq is just a work array
C        allocate(sigq(ndimh,ndimh))
C        allocate(sqi(ndimh,ndimh,1))
C        allocate(zqi(ndimh,ndimh))
C        allocate(z0(ndimh,ndimh))
C        allocate(z0g(ndimh,ndimh))
C        allocate(e(ndimh))
C        allocate(e0(ndimh))
C        allocate(ipm(ndimh))
CC#endif
C
CC   ... z0 <- Evecs of small LDA hamiltonian = Z1.lda
CC       call zprm('h before permutation',2,h,ndimh,ndimh,ndimh)
CC       call zprm('s before permutation',2,s,ndimh,ndimh,ndimh)
C        call permhk(0,nbas,ndimh,s_ham%offH,h(1,1,ispc),hqi,ndim1)
CC       call zprm('h after permutation',2,hqi,ndimh,ndimh,ndimh)
C        call permhk(0,nbas,ndimh,s_ham%offH,s(1,1,ispc),sqi,ndim1)
CC       z0 used here as a work array
C        call zhevx(ndim1,ndimh,hqi,sqi,1,T,ndim1,9d9,i,sigq,F,e,ndimh,
C     .    z0)
CC       call zprm('LDA evecs of small h',2,zqi,ndimh,ndim1,ndim1)
C
CC   ... zqi <- Evecs of large LDA hamiltonian = Z.lda
CC       Also e = LDA eigenvalues
CC       call zprm('h before permutation',2,h(1,1,ispc),ndimh,ndimh,ndimh)
C        call permhk(0,nbas,ndimh,s_ham%offH,h(1,1,ispc),hqi,ndim1)
CC       call zprm('h after permutation',2,hqi,ndimh,ndimh,ndimh)
C        call permhk(0,nbas,ndimh,s_ham%offH,s(1,1,ispc),sqi,ndim1)
CC       z0 used here as a work array
C        call zhevx(ndimh,ndimh,hqi,sqi,1,T,ndimh,9d9,i,sigq,F,e,ndimh,
C     .    zqi)
CC        call zprm('LDA evecs of entire h',2,zqi,ndimh,ndimh,ndimh)
CC        call permhk(4,nbas,ndimh,s_ham%offH,zqi,hqi,ndim1)
CC        call zprm('permuted z',2,hqi,ndimh,ndimh,ndimh)
CC        zqi = 0
CC        call permhk(5,nbas,ndimh,s_ham%offH,zqi,hqi,ndim1)
CC        call zprm('restore permuted z',2,zqi,ndimh,ndimh,ndimh)
CC
CC   ... sigq <- Evecs of small GW hamiltonian = Z1.gw
CC       Call to bloch : 104110=perm orb, transpose, add, c*16
C        i = 100000 + 4000 + 40*(1-hreal) + 100 + 10
C        call bloch(i,qs,nl,plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,
C     .    s_ham%hrs,ndhrs,isp,nsp,ndimh,ndimh,0,ndimh,0,ndimh,0,
C     .    h(1,1,ispc),w,w)
CC       call zprm('h+sigm-vxc',2,h(1,1,ispc),ndimh,ndimh,ndimh)
C        call permhk(0,nbas,ndimh,s_ham%offH,h(1,1,ispc),hqi,ndim1)
C        call permhk(0,nbas,ndimh,s_ham%offH,s(1,1,ispc),sqi,ndim1)
C        call zhevx(ndim1,ndimh,hqi,sqi,1,T,ndim1,9d9,i,sigq,F,e0,ndimh,
C     .    z0g)
CC       call zprm('GW evecs of small h',2,sigq,ndimh,ndim1,ndim1)
C
CC   ... sigq <- (Z1.lda)+ S11 Z1.sig
C        call permhk(0,nbas,ndimh,s_ham%offH,s(1,1,ispc),sqi,ndim1)
C        call zgemm('C','N',ndim1,ndim1,ndim1,one,z0,ndimh,sqi,ndimh,
C     .    zer,hqi,ndimh)
C        call zgemm('N','N',ndim1,ndim1,ndim1,one,hqi,ndimh,z0g,ndimh,
C     .    zer,sigq,ndimh)
CC        call zprm('z1(lda)+ S z1(sig)',2,sigq,ndimh,ndim1,ndim1)
C
CC   ... z0 <- Z.lda (Z1.lda)+ S11 Z1.sig
C        call zgemm('N','N',ndimh,ndim1,ndim1,one,zqi,ndimh,sigq,ndimh,
C     .    zer,z0,ndimh)
C        call zmscop(0,ndimh,ndimh-ndim1,ndimh,ndimh,0,ndim1,0,ndim1,
C     .    zqi,z0)
C
CC   ... zqi <- Z.lda (Z1.lda)+ S11 Z1.sig in original orbital order
C        call permhk(5,nbas,ndimh,s_ham%offH,zqi,z0,ndim1)
CC       call zprm('z in original orbital order',2,zqi,ndimh,ndimh,ndimh)
C
CC   ... subtract sigma from h to restore LDA h; also s <- z+ sig z
CC       Call to bloch : 104110=perm orb, transpose, no add, c*16
C        i = 100000 + 4000 + 40*(1-hreal) + 000 + 10
C        call bloch(i,qs,nl,plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,
C     .    s_ham%hrs,ndhrs,isp,nsp,ndimh,ndimh,0,ndimh,0,ndimh,0,sigq,w,w)
C        call daxpy(ndimh**2*2,-1d0,sigq,1,h(1,1,ispc),1)
C        call phmbls(3,ndimh,ndimh,e0,ipm,z0,sigq,zqi,zqi,s(1,1,ispc))
C        if (ddot(3,qs,4,qs,4) < tol**2 .or. .true.) then
C          call phmbl3(2,ndimh,0,e,e0,sigp,fasig,fbsig)
CC          call info5(30,0,0,' hambls: fit sigma for e>%d :'//
CC     .      ' sig = %,4;4d + %,4;4d * elda (q=%3:1,5;5d)',
CC     .      emax,fasig,fbsig,qs,0)
C        endif
C
CC   ... Overwrite z+ sigma z for high-lying energies
C        call psigd(0,ndimh,s(1,1,ispc),e,sigp,n123,wk)
C
CC   ... h <- z+ h(lda) z + z+ sigma z
C        call phmbls(1,ndimh,ndimh,e0,ipm,z0,h(1,1,ispc),zqi,zqi,h(1,1,ispc))
C        call daxpy(ndimh**2,1d0,s(1,1,ispc),1,h(1,1,ispc),1)
C
CC   ... h <- z+ h z ; s <- z ; reorder evecs,h in ascending eval order
C        call phmbls(31-1,ndimh,ndimh,e,ipm,sigq,h(1,1,ispc),zqi,zqi,
C     .    h(1,1,ispc))
C        call zcopy(ndimh**2,zqi,1,s(1,1,ispc),1)
C        mode = -1
C
CC#ifdef F90
C        deallocate(hqi,sigq,sqi,zqi,z0,z0g,e,e0,ipm)
CC#endif

C --- Interpolate sigma using eigenvectors ---
      elseif (moditp == 3) then

        if (moditp == 3) then
C          if (qplin(qs,plat,nqfbz,qfbz,0,iqfbz,wt) < 0) then
            call rx1('hambls: wt gen failed for qp=%3:1,5;5d',qs)
C          endif
        endif

CC ... Get qpi = four k-points in BZ of sigma matrix nearest to qp
C      do  i = 1, nqsig
C        do  j = 1, 3
C          qpi(j,1) = qsig(j,i) - qs(j)
C        enddo
C        call shorbz(qpi,qloc(1,i),qlat,plat)
C      enddo
C      call dvheap(3,nqsig,qloc,iprm,tol,11)
C      do  i = 1, 4
C        call dcopy(3,qloc(1,iprm(i)),1,qpi(1,i),1)
C      enddo
C
CC --- If qpi(1) = q, Bloch transform sigma at q and exit ---
CC     if (ddot(3,qpi,4,qpi,4) < tol**2) then
C      if (ddot(3,qpi,1,qpi,1) < tol**2 .or. moditp == 0) then
C
CC       Standard interpolation.  Too bad it doesn't work well!
CC       call zprm('h(lda)',2,h,ndimh,ndimh,ndimh)
CC       call zprm('s(lda)',2,s,ndimh,ndimh,ndimh)
CC       Call to bloch : 104110=perm orb, transpose, add, c*16
C        i = 100000 + 4000 + 40*(1-hreal) + 100 + 10
C        call bloch(i,qs,nl,plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,
C     .    s_ham%hrs,ndhrs,isp,nsp,nlmto,nlmto,0,ndimh,0,ndimh,0,
C     .    h(1,1,ispc),xx,xx)
CC       call zprm('h+sigm-vxc',2,h(1,1,ispc),ndimh,ndimh,ndimh)
C        goto 999
C      endif
C
CC ... Find weights for each of four points:
CC     nsim = # independent pts; wt = weights for points
CC  99 continue
C      do  i = mxitrp, 1, -1
C        nsim = i
CC       Transpose of qpi into norm
C        call dmcpy(qpi,1,3,norm(1,2),4,1,4,3)
C        call dvset(norm,1,4,1d0)
CC       call dcopy(16,qpi,1,norm,1)
CC       call dvset(norm(1,2),1,12,1d0)
CC       call prmx('norm',norm,4,4,4)
C        j = dinv(' ',i,4,norm)
CC       If inversion failed, reduce the number of points
C        if (j /= 0) goto 10
CC       If any weight exceeds 1 or swt<>1, reduce the number of points
C        swt = 0
C        do  j = 1, i
C          if (norm(1,j) > 5) then
C            goto 10
C          endif
C          swt = swt + norm(1,j)
C        enddo
C        if (abs(swt-1) > tol) goto 10
CC       All tests passed ... use this set of points
C        goto 20
C   10   continue
C      enddo
C   20 continue
CC     call prmx('inverse',norm,4,nsim,nsim)
CC     norm(1,1..nsim) are weights
C
CC     for debugging
CC     qpi(1,1) = 0.75d0 - qs(1)
CC     qpi(2,1) = 0.00d0 - qs(2)
CC     qpi(3,1) = 0.00d0 - qs(3)
C
CC     Assign weights to each qp; compute sum-of-weights as check
C      swt = 0
C      j = nsim
C      do  i = nsim, 1, -1
C        wt(i) = norm(1,i)
C        if (abs(wt(i)) < tol .and. j == i) j=j-1
C        if (j >= i) swt = swt + wt(i)
C      enddo
C      nsim = j
C      if (abs(swt-1) > tol) then
C        call info5(0,0,0,' hambls failed for q=%3:1,5;5d, '//
C     .  'wt=%n:1,5;5d',qs,nsim,wt,0,0)
CC       goto 99
C        call rx1('hambls:  interpol. weights'//
C     .  ' should sum to 1; actual sum=%d',swt)
C      endif
C
C      allocate(hqi(ndimh,ndimh,nsim))
C      allocate(sqi(ndimh,ndimh,nsim))
C      allocate(zqi(ndimh,ndimh))
C      allocate(z0(ndimh,ndimh))
C
CC ... Inverse z0^-1, where z0 = eigenvectors at desired q
C      call zcopy(ndimh**2,h(1,1,ispc),1,hqi,1)
C      call zcopy(ndimh**2,s(1,1,ispc),1,sqi,1)
CC      call zprm('h(lda)',2,hqi,ndimh,ndimh,ndimh)
CC      call zprm('s(lda)',2,sqi,ndimh,ndimh,ndimh)
C      call zhev(ndimh,hqi,sqi,T,T,ndimh,9d9,nev,wk,F,-1,evl,z0)
CC     call zprm('z(lda)',2,z0,ndimh,ndimh,ndimh)
C      call zgetrf(ndimh,ndimh,z0,ndimh,ipiv,i)
C      if (i == 0) call zgetri(ndimh,z0,ndimh,ipiv,hqi,ndimh**2,i)
C      call rxx(i /= 0,'hambls: failed to invert eigenvectors')
CC     call zprm('z(lda)^-1',2,z0,ndimh,ndimh,ndimh)
C
CC ... For each qpi(1..nsim), add wt(qpi) * sigma(qpi)
C      call info5(36,0,0,' Interp SE for q=%3:1,5;5d, '//
C     .  'wt=%n:1,5;5d',qs,nsim,wt,0,0)
C
C      do  j = 1, nsim
C        dq = dsqrt(ddot(3,qpi(1,j),1,qpi(1,j),1))
C        call daxpy(3,1d0,qs,1,qpi(1,j),1)
C        call info5(50,0,0,' qi(%i) = %3:-1,6;6d  |qi-q|=%,6;6d',
C     .    j,qpi(1,j),dq,0,0)
C        call hambl(0,s_site,s_spec,s_lat,s_ham,isp,qpi(1,j),
C     .    k1,k2,k3,smpot,vconst,lcplxp,alfa,ndimh,napw,
C     .    igapwl,hqi,sqi,xx)
CC       call zprm('h(lda)',2,hqi,ndimh,ndimh,ndimh)
CC       call zprm('s(lda)',2,sqi,ndimh,ndimh,ndimh)
C        call zhev(ndimh,hqi,sqi,T,T,ndimh,9d9,nev,wk,F,-1,evl,zqi)
CC       call zprm('z(lda)',2,zqi,ndimh,ndimh,ndimh)
C
CC       hqi <- z(qpi) z(q)^-1
C        call zgemm('N','N',ndimh,ndimh,ndimh,(1d0,0d0),zqi,
C     .    ndimh,z0,ndimh,(0d0,0d0),hqi,ndimh)
C
CC       sqi <- sigma(qpi)
CC       Call to bloch : 104010=perm orb, transpose, c*16
C        i = 100000 + 4000 + 40*(1-hreal) + 10
C        call bloch(i,qpi(1,j),nl,plat,mxorb,s_ham%iprmb,1,nttabs,
C     .    s_ham%iaxs,s_ham%hrs,ndhrs,isp,nsp,nlmto,nlmto,0,ndimh,0,
C     .    ndimh,0,sqi,xx,xx)
CC       call zprm('sig at qpi',2,sqi,ndimh,ndimh,ndimh)
C
CC       zqi <- (z(qpi) z(q)^-1)+ sigma(qpi)
C        call zgemm('C','N',ndimh,ndimh,ndimh,(1d0,0d0),
C     .    hqi,ndimh,sqi,ndimh,(0d0,0d0),zqi,ndimh)
C
CC       sqi <- (z(qpi) z(q)^-1)+ sigma(qpi) (z(qpi) z(q)^-1)
C        call zgemm('N','N',ndimh,ndimh,ndimh,(1d0,0d0),
C     .    zqi,ndimh,hqi,ndimh,(0d0,0d0),sqi,ndimh)
C
CC       call zprm('sig(from qpi)',2,sqi,ndimh,ndimh,ndimh)
C
CC       Add wt(qpi) sigma(translated from qpi) to hamiltonian
C        call daxpy(ndimh**2*2,wt(j),sqi,1,h(1,1,ispc),1)
C
C      enddo
C
C      deallocate(hqi,sqi,zqi,z0)
      else
        call rxi('unknown interpolation mode %i',moditp)
      endif

C ... At this point, e0 has evals, s(:,:,ispc) has evecs of hamiltonian
   98 continue
      if (allocated(sigq)) deallocate(sigq)

C     Overwrite h with eigenvalues
      call dpzero(h(1,1,ispc),ndimh**2*2)
      do  i = 1, ndimh
        h(i,i,ispc) = e0(i)
      enddo
C      call zprm('h+sig in evec repsn',2,h(1,1,ispc),ndimh,ndimh,ndimh)
!      call prmx('evals of h+sig',e0,ndimh,ndimh,1)
!      call zprm('evecs of h+sig',2,s(1,1,ispc),ndimh,ndimh,ndimh)
C      call phmbls(32+1,ndimh,nlmto,e0,xx,z0,h,s,s,dsdk)  ! Overwrites s, dsdk
C      call zprm('h+sig in orbital representation',2,dsdk,ndimh,ndimh,ndimh)
C      j = 3
C      call phmbls(32+1,ndimh,nlmto,e0,xx,z0,dsdk(1,1,j,1),s,s,dsdk(1,1,j,1))   ! Comment out h+sig, orb
C      call zprm('dsdk in orbital representation',2,dsdk(1,1,j,1),ndimh,ndimh,ndimh); stop

C ... Add off-diagonal spin-orbit hamiltonian to h; re-diagonalize
C     h(:,:,1..2) = ++ and -- blocks; h(:,:,3) = L+S- block
      if (lso == 4) call rx('hambls not ready for lso=4')
      if ((lso == 1 .or. lso == 3) .and. ispc == 2) then
        if (ldsdk) then
          stop 're-diagonalized dsdk'
        endif
        call sopert3(10*lso+1,s_ham,ndimh,h,s,xx)
      endif

C ... Cleanup
   99 mode = -1
      deallocate(zqi,z0,e,e0)

      call getpr(ipr)
      if (ipr >= 90 .and. mode >= 0) then
        call info(0,0,0,'h and s for q=%3:1;6,6d',qs,0)
        call zprm('h+sig',2,h(1,1,ispc),ndimh,ndimh,ndimh)
      endif
C     i = fopna('out',-1,0)
C     print *, 'hambls, h'
C     call ywrm(0,'h(k)',2,i,'(9f20.10)',h(1,1,ispc),1,ndimh,ndimh,
C    .  ndimh)
C     call rx0('done')

  999 continue
      call tcx('hambls')

      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        deallocate(igapwl)
      endif

C     For some crazy reason, the return is needed to compile using
C     pgf90 6.0-2 on x86-64 Linux  !!
      return
      end
      subroutine blsig(mode,nbas,s_ham,isp,nsp,nspc,plat,q,
     .  lwtkb,zval,iq,wtkp,wtkb,ndimh,zq,sig)
C- sigma-vxc or sum_occ <sigma-vxc> for one k-point.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit specifies what to do:
Ci         :  0 Make sigma and store in sig
Ci         :  1 Make sum_states wtkb(iq) (zq+)*sig*zq
Ci         :10s digit specifies self-energy matrix elements
Ci         :  0 No self-energy to compute
Ci         :  1 r.s. sigma is real
Ci         :  2 r.s. sigma is complex
Ci   nbas  :size of basis
Ci     Elts read: ndhrs s s s oindxo
Ci     Stored:
Ci     Passed to:
Ci   isp   :generate self-energy for this spin
Ci         :Not used in the noncollinear case
Ci   plat  :primitive lattice vectors
Ci   q     :Bloch vector (k-point)
Ci   wtkp  :q-point weights from symmetry operations
Ci   lwtkb :(used for 1s digit mode=1)
Ci         :0 assign weight 1 for first zval states, 0 otherwise
Ci         :1 weights are required a priori, and are read from disk
Ci   wtkb  :integration weights, needed if lwtkb is 1
Ci   ndimh :dimension of hamiltonian and self-energy sig
Ci   zq    :evecs (used for 1s digit mode=1)
Cio Inputs/Outputs
Cio  sham  :struct for parameters defining hamiltonian; see routine uham
Cio    Elts read: ndhrs s s s oindxo eterms
Cio    Stored:    eterms  (1s digit mode=1)
Co Outputs
Co   sig   :(1s digit mode = 0) self-energy matrix
Co         :(1s digit mode = 1) <zq sig zq>
Co         :Note: in noncollinear case, passed sig MUST be
Co         :dimensioned sig(ndimh,nspc,ndimh,nspc) even though
Co         :dimension appears as sig(ndimh,ndimh) in this routine.
Cl Local variables
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   11 Jan 05 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,iq,isp,ndimh,nsp,nspc,lwtkb
      double precision plat(3,3),q(3),wtkb(ndimh,nsp,*),wtkp(iq),zval
      double complex zq(ndimh,nspc,ndimh,nspc),sig(ndimh,ndimh,2)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
C ... Local parameters
      integer i,lrsig,ipr,jsp,ksp
      integer hreal,nl,mxorb,nttabs,ndhrs,nevec,ndimhx,nlmto
      logical, parameter :: T=.true., F=.false.
      double precision :: eterms(22),rhosig,swgt=0,xx
      equivalence (eterms(19),rhosig)
      complex(8),allocatable:: wk(:,:,:,:)
      double precision zsz(ndimh*nspc),ewgt(ndimh*nspc)
      procedure(integer) :: nglob

      lrsig = mod(mode/10,10)
      if (lrsig == 0) return
      call tcn('blsig')

C ... Setup for Bloch transform
      ndhrs = s_ham%ndhrs
      nlmto = s_ham%nlmto
C     hreal=1 -> sig(T) is real
      hreal = 0
      if (lrsig >= 4) hreal = 1
      nl = nglob('nl')
      mxorb = nglob('mxorb')
      nttabs = s_ham%nprs(nbas+1)
      call sanrg(.true.,isp,1,nsp,'blsig:','isp')
      ndimhx = ndimh*nspc

      allocate(wk(ndimh,nspc,ndimh,nspc))

C --- Bloch sum sigma-vxc, orbital basis ---
C     Call to bloch : 104110=perm orb, transpose, no add, c*16
      i = 100000 + 4000 + 40*(1-hreal) + 000 + 10
      do  jsp = 1, nspc
      ksp = isp
      if (nspc == 2) ksp = jsp
C     collinear case:     ksp = isp, jsp=1
C     noncollinear case:  ksp=1..2,  jsp=1..2
      call bloch(i,q,nl,plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,
     .  s_ham%hrs,ndhrs,ksp,nsp,nlmto,nlmto,0,ndimh,0,ndimh,0,
     .  sig(1,1,jsp),xx,xx)
      enddo
C     Convert sig to full matrix with zero off-diagonal spin blocks
      if (nspc == 2) then
        call sopert(1,ndimh,nspc,wk,sig,sig)
      endif

      call getpr(ipr)
      if (ipr >= 100 .and. mode >= 0) then
        call info(0,0,0,'sig for q=%3:1;6,6d',q,0)
        call zprm('sig',2,sig,ndimhx,ndimhx,ndimhx)
      endif

C ... ksp = current spin index (isp or 1, in noncollinear case)
      ksp = isp
      if (nspc == 2) ksp = 1

C --- Compute double-counting contribution to total energy ---
      if (mod(mode,10) == 1) then
C  ...  zsz <- diagonal part of <sigma-vxc>
        call phmbls(3,ndimhx,ndimhx,zsz,xx,wk,sig,zq,zq,sig)
C  ...  Make weights for d.c. sum
        call dpzero(ewgt,ndimhx)
        if (lwtkb == 0) then
          nevec = zval/2d0 + 0.500001d0
          if (nspc == 2) nevec = 2*nevec
          do  i = 1, nevec
            ewgt(i) = abs(wtkp(iq))/nsp
          enddo
        else
          do  i = ndimhx, 1, -1
            nevec = i
            if (abs(wtkb(i,ksp,iq)) > 1d-8) goto 12
          enddo
   12     continue
          do  i = 1, nevec
            ewgt(i) = abs(wtkb(i,ksp,iq))
          enddo
        endif
C  ...  Contribution to d.c. sum for this qp
        eterms = s_ham%eterms
        if (iq == 1 .and. ksp == 1) then
          rhosig = 0
          swgt = 0
        endif

        do  i = 1, nevec
          swgt = swgt + ewgt(i)
          rhosig = rhosig + ewgt(i)*zsz(i)
        enddo
        s_ham%eterms = eterms
C  ...  Sanity check
C       print *, 'iq=',iq,'swgt=',swgt,'rhosig=',rhosig
      endif

      deallocate(wk)

      call tcx('blsig')
      end
      subroutine permhk(mode,nbas,ldh,offH,h,hp,lblock)
C- Permutation of H by kappa blocks
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 copy permuted h into hp.  Both rows, cols permuted.
Ci         :1 copying in reverse permutation order hp to h.
Ci            mode=1 is inverse operation of mode 0.
Ci            Thus this sequence of calls:
Ci              call permhk(0,nbas,ldh,offH,h,hp,lblock)
Ci              h = 0
Ci              call permhk(1,nbas,ldh,offH,h,hp,lblock)
Ci            copies h to hp, permuting hp; then initializes h and
Ci            restores the original h from hp.
Ci         :2 copy permuted h into hp.  Only cols are permuted.
Ci         :3 copying in reverse permutation order hp to h.
Ci            mode=3 is inverse operation of mode 2.
Ci            Thus this sequence of calls:
Ci              call permhk(2,nbas,ldh,offH,h,hp,lblock)
Ci              h = 0
Ci              call permhk(3,nbas,ldh,offH,h,hp,lblock)
Ci            copies h to hp, permuting hp; then initializes h and
Ci            restores the original h from hp.
Ci         :4 copy permuted h into hp.  Only rows are permuted.
Ci         :5 copying in reverse permutation order hp to h.
Ci            mode=5 is inverse operation of mode 4
Ci            Thus this sequence of calls:
Ci              call permhk(4,nbas,ldh,offH,h,hp,lblock)
Ci              h = 0
Ci              call permhk(5,nbas,ldh,offH,h,hp,lblock)
Ci            copies h to hp, permuting hp; then initializes h and
Ci            restores the original h from hp.
Ci   nbas  :size of basis
Ci   ldh :dimension of h
Ci   offH  :table of offsets to unpermuted hamiltonian (makidx.f)
Ci         :offH(m,ik,ib) = offset to subblock (ik,ib), where
Ci         :m=1 for lower orbitals, m=2 for intermediate, etc.
Ci         :offH(m,ik+1,ib)-offH(m,ik,ib) = size of block
Cio Inputs/Outputs
Cio  h     :mode=0: (input) matrix to be permuted
Cio        :mode=2: (input) matrix to be permuted
Cio        :mode=1: (output) matrix to be restored from hp.
Cio        :        Both rows and columns are (reverse) permuted
Cio        :mode=3: (output) matrix to be restored from hp
Cio        :        Only columns are (reverse) permuted.
Cio  hp    :mode=0: (output) matrix copied from h where both
Cio                 rows and columns are permuted.
Cio        :mode=2: (output) matrix copied from h where only
Cio                 columns are permuted.
Cio        :mode=1: (input) matrix used to restore h
Cio        :mode=3: (input) matrix used to restore h
Co Outputs
Co   lblock:dimension of lower block
Cu Updates
Cu   17 May 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,ldh,lblock
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offh(n0H,nkap0,nbas+1)
      double complex h(ldh,ldh),hp(ldh,ldh)
C ... Local parameters
      integer nblk,m,i,ib,ikap

      integer iprm(nbas,nkap0),offs(nkap0,nbas+1)
      double precision xx

      nblk = nkap0*nbas
      m = 1
      call icopy(nblk+1,offH(m,1,1),n0H,offs,1)
      i = 0
      lblock = 0
      do  ib = 1, nbas
        do  ikap = 1, nkap0
          i = i+1
          iprm(ib,ikap) = i
        enddo
        lblock = lblock + offs(2,ib) - offs(1,ib)
      enddo

      if (mode == 0) then
        call zpmblk(nblk,iprm,offs,xx,ldh,h,8120,xx,ldh,ldh,hp,i)
      elseif (mode == 1) then
        call zpmblk(nblk,iprm,xx,offs,ldh,hp,8202,xx,ldh,ldh,h,i)
      elseif (mode == 2) then
        call zpmblk(nblk,iprm,offs,xx,ldh,h,120,xx,ldh,ldh,hp,i)
      elseif (mode == 3) then
        call zpmblk(nblk,iprm,xx,offs,ldh,hp,202,xx,ldh,ldh,h,i)
      elseif (mode == 4) then
        call zpmblk(nblk,iprm,offs,xx,ldh,h,4120,xx,ldh,ldh,hp,i)
      elseif (mode == 5) then
        call zpmblk(nblk,iprm,xx,offs,ldh,hp,4202,xx,ldh,ldh,h,i)
      else
        call rxi('permhk: bad mode',mode)
      endif

      end

      subroutine phmbl3(mode,ndimh,nsmidb,e,de,sigp,a,b)
C- Fit coffs a,b in linear fit a+b*e(lda) = de[e(lda)]
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 initialize normal matrix
Ci         :2 add points to normal matrix
Ci         :4 solve normal matrix for coefficients a and b
Ci         :8 print out nsmidb
Ci         :NB: these switches may be taken in combination
Ci   ndimh :number of eigenvalues
Ci   nsmidb:smallest value of nmax encountered
Ci         :(only used in printout, mode=8)
Ci   e     :list of (LDA) eigenvalues
Ci   de    :list of self-energies
Ci   sigp  :parameters for approximating self-energy sigma.  sigma
Ci         :is approximated by its diagonal part sigii for energies
Ci         :below a low-energy cutoff (specified nmin or emin) and
Ci         :above a low-energy cutoff (specified nmax or emax).
Ci         :NB: this routine uses parts of sigp only to determine
Ci         :    which states are to be included in fit
Ci         : arg 1: mode : specifies how to set its diagonal part
Ci         :        for states above the high-energy cutoff nmax or emax
Ci         :        0 constrain sigii to be > asig+bsig*e
Ci         :        1 constrain sigii to be = asig+bsig*e
Ci         :        2 constrain sigii to be > asig and < bsig
Ci         :        3 constraint same as mode 1.
Ci         :          Mode 3 differs in that the least-squares fit to
Ci         :          sigii (for informational purposes only, to help
Ci         :          estimate asig and bsig) is done for states between
Ci         :          efit and nmax or emax
Ci         : arg 2: nmin : sigma for states 1..nmin are approximated
Ci         :               by sigii
Ci         : arg 3: emin : (used only if nmin<0)
Ci         :             : sigma for levels e<emin are approximated
Ci         :               by sigii
Ci         : arg 4: nmax : sigma for levels i>nmax are approximated by
Ci         :               sigii AND constrained according to mode
Ci         : arg 5: emax : (used only if nmax<=0)
Ci         :             : sigma for levels e<emax are approximated by
Ci         :               sigii AND constrained according to mode
Ci         : arg 6: asig : constraint used to approximate
Ci         :               sigii = asig + E * bsig
Ci         : arg 7: bsig : constraint used to approximate
Ci         :               sigii = asig + E * bsig
Ci         : arg 8: efit : fit sigii between efit and emax
Ci         : arg 9: lwrite:write sigii to 'sigii' file
Co Outputs
Co   a     :(4s bit mode set) constant coefficient
Co   b     :(4s bit mode set) linear coefficient
Cr Remarks
Cr
Cu Updates
Cu    4 Nov 03 Write SE to sigii file also for states not in fit
Cu   25 Jun 03 Update to incorporate sigp mode 3
Cu   21 May 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ndimh,nsmidb
      double precision e(ndimh),de(ndimh),a,b,sigp(10)
C ... Local parameters
      logical lhigh,lfit
      integer i,nmin,nmax,modsgp,fopna,ifi,lwrite
      double precision norm(2,2),rhs(2),tnorm(2,2),emin,emax,elo,ehi,
     .  efit
      save norm,rhs,elo,ehi
      data norm /4*0d0/,rhs/2*0d0/,elo/9d0/,ehi/-9d0/

      if (mod(mode,2) == 1) then
        call dpzero(norm,4)
        call dpzero(rhs,2)
      endif
      if (mode == 1) return

      modsgp = nint(sigp(1))
      nmin   = nint(sigp(2))
      emin   = sigp(3)
      nmax   = nint(sigp(4))
      emax   = sigp(5)
C     asig   = sigp(6)
C     bsig   = sigp(7)
      efit   = sigp(8)
      lwrite = nint(sigp(9))

      if (mod(mode/2,2) == 1) then
        if (lwrite == 1) then
          ifi = fopna('sigii',-1,0)
          write(ifi,'(''#  state    E(lda)      dE(sig)'')')
        endif
        do  i = 1, ndimh
          lhigh = (nmax <= 0 .and. e(i) > emax) .or.
     .            (nmax > 0 .and. i > nmax)
          if (modsgp == 3) then
            lfit = e(i) > efit .and. .not. lhigh
          else
            lfit = e(i) > efit
          endif
          if (lfit) then
            norm(1,1) = norm(1,1) + 1
            norm(1,2) = norm(1,2) + e(i)
            norm(2,1) = norm(2,1) + e(i)
            norm(2,2) = norm(2,2) + e(i)**2
            rhs(1) = rhs(1) + de(i)
            rhs(2) = rhs(2) + e(i)*de(i)
            elo = min(elo,de(i))
            ehi = max(ehi,de(i))
            if (lwrite == 1) then
              ifi = fopna('sigii',-1,0)
              write(ifi,'(i6,2f12.6)') i,e(i),de(i)
            endif
          else
            if (lwrite == 1) then
              ifi = fopna('sigii',-1,0)
              write(ifi,'(''#'',i5,2f12.6)') i,e(i),de(i)
            endif
          endif
        enddo
      endif

      if (mod(mode/8,2) == 1 .or. mod(mode/4,2) == 1)
     .  call info0(30,0,0,' ')


      if (mod(mode/8,2) == 1) then
        call info2(30,0,0,' hambls: smallest nmax encountered for '//
     .    'sigm middle block = %i',nsmidb,0)
      endif

      if (mod(mode/4,2) == 1) then
        if (norm(1,1) /= 0) then
          call dinv22(norm,tnorm)
          a = tnorm(1,1)*rhs(1) + tnorm(1,2)*rhs(2)
          b = tnorm(2,1)*rhs(1) + tnorm(2,2)*rhs(2)
          call info5(30,0,0,' hambls: sig(low,high) = %,4;4d,%,4;4d'//
     .      '  fit : %,4;4d + %,4;4d * E(lda) (%i points)',
     .      elo,ehi,a,b,nint(norm(1,1)))
          if (lwrite == 1) then
            ifi = fopna('sigii',-1,0)
            write(ifi,345) modsgp,nmin,emin,nmax,emax,a,b,efit
  345       format('# mode nmin    emin   nmax    emax',8x,
     .        'asig      bsig',8x,'efit'/2i5,f10.4,i5,f10.4,2x,2f10.4,
     .        2x,f10.4)
          endif
        endif
      endif

      end
C ... checked :
C  against num. differentiation: dsigorbraw[xyz]
C  dsigorbrawx -> dsigldarawx checked against mc *(1 below) Also y,z
C  dsigldarawx -> dsigldax checked against mc *(2 below).  Also y,z
C  dsigorbx -> dsigldax checked against mc *(3 below).  Also y,z
C  Internally consistent: dsigorbrawx dsigldarawx dsigorbx dsigldax dsigqsgwx
C
C Show hlda is diagonal in lda eigenfunction basis
C   mc -f9f15.10 zlda -cc -t hlda -x zlda -x  -sub 1,6,1,6
C   mc zlda -cc -t hlda -x zlda -x -p -v2dia -v2dia -- -px:5
C Show hlda+sigorb gives QSGW evals
C   mc -f9f15.10 hlda sigorb -+ -p -p slda -tog -gevl  eqsgw -- -px:6
C Make zqsgw by hand (tricky to compare to hambls)
C   mc -f9f16.10 hlda sigorb -+ -p -p slda -tog -gevc > zqsgw
C Show zqsgw+ (hlda+sigorb) zqsgw = eqsgw
C   mch -f9f15.10 zqsgw -cc -t hlda sigorb -+  zqsgw -x -x eqsgw -v2dia -- -abs -max:g
C Make zlda2qsgw in full hilbert space
C   mc -f9f16.11 zlda2qsgw -1:29 -a one one -s0 -ccat one -s0 one -ccat -rcat > zlda2qsgwfull
C Show (zlda zlda2qsgwfull)+ (hlda+sigorb) (zlda zlda2qsgwfull) accomplishes the same thing
C   mch -f9f15.10 zlda zlda2qsgwfull -x -a z z -cc -t hlda sigorb -+  z -x -x eqsgw -v2dia -- -abs -max:g
C Show sigqsgw = zlda2qsgwfull+ siglda zlda2qsgwfull
C   mch -f9f15.10 zlda2qsgwfull -a z z -cc -t siglda z -x -x  sigqsgw -- -abs -max:g -px
C *(1) dsigx in lda eigenfunction basis, compare to what hambls generates
C   mc -f9f15.10 zlda -cc -t dsigorbrawx -x zlda -x dsigldarawx -- -abs -max:g -px
C *(2) dsigx in lda eigenfunction basis, stripped, compare to what hambls generates
C   mc -f9f15.10 dsigldarawx dsigldax -- -sub 1,29,1,29 -px
C   mc -f9f15.10 dsigldax -split s 1,30,59 1,30,59 s12 -px
C   mc -f9f15.10 dsigldax -split s 1,30,59 1,30,59 s22 -px
C (3) rotate dsigx,y,z in orbital basis to LDA basis, show agreement
C   mc -f9f15.10 zlda -cc -t dsigorbx -x zlda -x dsigldax -- -abs -max:g -px
C   mc -f9f15.10 zlda -cc -t dsigorby -x zlda -x dsiglday -- -abs -max:g -px
C   mc -f9f15.10 zlda -cc -t dsigorbz -x zlda -x dsigldaz -- -abs -max:g -px
C Inverse: rotate dsigx in LDA basis to orbital basis, show agreement
C   mc -f9f15.10 zlda -i -a z z -cc -t dsigldax -x z -x dsigorbx  -- -abs -max:g -px
C Compare dsigqsgwx from zlda2qsgwfull+ dsigldax zlda2qsgwfull ... disagreement in phase
C   mc -f9f15.10 zlda2qsgwfull -a z z -cc -t dsigldax z -x -x  dsigqsgwx -- -abs -max:g
C This is because some evecs changed phase in the course of the development: signs are opposite
C   mc -f9f15.10 zlda2qsgwfull -a z z -cc -t dsigldax z -x -x dsigqsgwx -rowl 25 -coll 17
C But evals compare ...
C   mc -f9f15.10 zlda2qsgwfull -a z z -cc -t dsigldax z -x -x -herm -evl dsigqsgwx -herm -evl -- -px
C   mc -f9f15.10 zlda2qsgwfull -a z z -cc -t dsiglday z -x -x -herm -evl dsigqsgwy -herm -evl -- -px
C   mc -f9f15.10 zlda2qsgwfull -a z z -cc -t dsigldaz z -x -x -herm -evl dsigqsgwz -herm -evl -- -px
C rotate dsigx in orbital basis to QSGW basis.  Can compare evales only
C   mc -f9f15.10 zqsgw -cc -t dsigorbx -x zqsgw -x -herm -evl dsigqsgwx -herm -evl -- -px
C   mc -f9f15.10 zqsgw -cc -t dsigorby -x zqsgw -x -herm -evl dsigqsgwy -herm -evl -- -px
C   mc -f9f15.10 zqsgw -cc -t dsigorbz -x zqsgw -x -herm -evl dsigqsgwz -herm -evl -- -px
C... old
C hlda in QSGW eigenfunction basis
C mc z -cc -t hlda -x z -x  -sub 1,5,1,5
C Show hlda + sig in QSGW eigenfunction basis is diagonal
C mc z -cc -t hlda -x z -x  sigzqsgw  -+ -p -v2dia -v2dia -- -px:5
C mc -f9f15.10 out out.si -- -s1/.002 > snot
C mc -f9f15.10 dhqsgwx dhldax --
C mch zqsgw -cc -t hqsgw -x zqsgw -x -p -v2dia -v2dia -- -px:5
