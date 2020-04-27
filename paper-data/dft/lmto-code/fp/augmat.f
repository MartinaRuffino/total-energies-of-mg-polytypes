      subroutine augmat(s_site,z,rmt,rsma,lmxa,pnu,pnz,kmxax,kmax,nlml,
     .  a,nr,nsp,nspc,lx,lso,lbf,loptic,nrg,nlg,facso,rofi,rwgt,cg,jcg,
     .  indxcg,v0,v1,v2,gpotb,gpot0,nkaph,nkapi,lmxh,lh,eh,rsmh,ehl,
     .  rsml,rs3,vmtz,lcplxp,lmaxu,vorb,lldau,iblu,idu,nbf,bloc,ppnl,
     .  hab,vab,sab,rgrad,rgrade,buf)
C- Make augmentation matrices sig,tau,pi for one site
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  sighh tauhh pihh sighhx tauhhx pihhx sohh sigkk
Ci                 taukk pikk sigkkx taukkx pikkx sokk sighk tauhk pihk
Ci                 sighkx tauhkx pihkx sohk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   z     :nuclear charge
Ci   rmt   :augmentation sphere radius
Ci   rsma  :augmentation smoothing radius
Ci   lmxa  :augmentation L-cutoff
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnz   :boundary conditions for second p.q.n (local orbital).
Ci          10s digit controls how local orbital included in hamiltonian
Ci   kmxax :dimensions rgrade (for optics)
Ci   kmax  :polynomial cutoff
Ci   nlml  :L-cutoff for density, potential on the radial mesh
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci         :  Used here to make off-diagonal part of optical matrix elements
Ci   lx    :0 poke results into s_site%{sig,tau,ppi}
Ci         :1 poke results into s_site%{sigx,taux,ppix}
Ci   lso   :0 no SO coupling
Ci         :1,3,4 Add L.S coupling to potential ppi
Ci         :2 Add LzSz only to potential ppi
Ci   lbf   :0 no external B field
Ci         :1 Add B-Sigma to potential ppi
Ci         :2 Add BzSigz only to potential ppi
Ci   loptic: calculate rgrad if loptic = 1
Ci   nrg   :dimensions rgrad, radial quantum number
Ci   nlg   :dimensions rgrad, l quantum number
Ci   facso :scale L.S  matrix elements by facso. Should normally be 1.
Ci   nbas  :size of basis
Ci   rwgt  :radial mesh weights
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   v0    :spherical potential that defines phi and phidot
Ci   v1    :true nonspherical potential, seen by augmented functions, excluding
Ci         :nuclear contribution 2*Z/r
Ci   v2    :smooth nonspherical potential, seen by unaugmented functions
Ci         :estat part from n2~ = n2 + compensating gaussians + pseudocore charge
Ci   gpotb :integrals of local gaussians times local smooth ves;
Ci         :see Remarks
Ci   gpot0 :integrals of local gaussians * phi0~
Ci         :phi0~ is the estatic potential of the interstitial;
Ci         :see Remarks
Ci   nkaph :number of orbital types for a given L quantum no. in basis
Ci         :dimensions sig,tau,ppi
Ci   nkapi :number of valence envelope function types per l q.n.
Ci   lmxh  :largest l in basis; must be >= max(lh)
Ci   lh    :list of nkaph l-cutoffs in basis
Ci   eh    :energy of smoothed Hankel
Ci   rsmh  :smoothing radii of smoothed hankel for each l, energy.
Ci         :rsmh(l+1,ik) for must be zero for all ik=1..nkaph unless
Ci         :it channel is a valence or extended local orbital.
Ci     ... The following four parameters are used to extrapolate
Ci         quantities outside the MT radius.  Extrapolation is used to
Ci         correct matrix elements for local orbitals.
Ci   ehl   :energy of smoothed Hankel tail for local orbital
Ci   rsml  :smoothing radius for smoothed Hankel tail of local orbital
Ci   rs3   :smoothing radius for extrapolation of MT potential
Ci   vmtz  :muffin-tin zero: subtracted from V in the fitting procedure.
Ci         :The asymptotic form of V-vmtz is taken to be zero.
Ci  lcplxp :0 if ppi is real; 1 if ppi is complex
Ci   ...   The following are LDA+U-related
Ci   lmaxu :dimensioning parameter for U matrix
Ci   vorb  :orbital dependent potential matrices
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci          U on site ib with dmat beginning at dmats(*,lldau(ib))
Ci   idu   :idu(l+1)=1 => this l has a nonlocal U matrix
Ci   nbf   :dimensions bloc
Ci   bloc  :external magnetic field for this site
Ci   buf   :Buffer to contain redirected standard output.
Ci         :This routine uses awrite to write output to stdout.
Ci         :It gets redirected to buf if buffer mode is turned on
Cio Inputs/Outputs
Cio  iblu  :index to current LDA+U block
Cio        :On input, index to last LDA+U block that was accessed
Cio        :iblu will be incremented to from blocks at this site
Cio  rofi  :radial mesh points.  On input, rofi(1..nr) are made.
Cio        :if V is to be extrapolated outside its MT sphere, to
Cio        :V(1..nrbig), rofi(nr+1,nrbig) are also generated
Cio        :Thus MUST be dimensioned at least rofi(1..nrbig)
Cio        :nrbig is internally generated, but will not
Cio        :exceed parameter nrx defined vxtrap.
Co Outputs
Co   ppnl  :NMTO potential parameters
Co   hab   :matrix elements of the ham. with true w.f.  See Remarks.
Co   vab   :matrix elements of the pot. with true w.f.  See Remarks.
Co   sab   :matrix elements of    unity with true w.f.  See Remarks.
Co   augmentation matrices (see s_site)
Co  rgrad  :radial matrix elements of gradient operator (phi/phidot)
Co  rgrade :radial matrix elements of gradient operator (envelopes)
Cr Remarks
Cr   This subroutine implements the computation of matrices
Cr   sigma, tau, ppi that comprise the local (augmented) part of
Cr   the hamiltonian in the full-potential method described in
Cr      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
Cr      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
Cr      ed. (Springer-Verlag, Berlin) 2000.
Cr   See discussion in Eqns 20-29.
Cr
Cr  In this method, the integral of a product of two augmented functions
Cr  F~i and F~j, with corresponding (smooth) envelopes Fi and Fj, and
Cr  their one-center expansions F^i and F^j, is calculated as
Cr
Cr     int_all F~i F~j  =  int_all Fi Fj  + int_MT (F~i F~j - F^i F^j)
Cr
Cr  Augmented F~i matches continuously and differentiably to Fi.
Cr  The one-center expansion F^i of Fi is identical except the former
Cr  is truncated to a finite L-cutoff in its on-center expansion.
Cr
Cr     Fi  = sum_L=0..infinity Fi_L
Cr     F^i = sum_L=0..lmxa     Fi_L
Cr
Cr  In the following description, we assume one atom per cell and
Cr  one orbital per L channel to simplify notation: i->L and j->L'.
Cr  If FL is further approximated by a polynomial expansion inside
Cr  one augmentation sphere, viz
Cr
Cr     F~L = sum_kL CkL P~kL                                 (cf Eq 16)
Cr     F^L = sum_kL CkL PkL
Cr
Cr  the integral becomes
Cr
Cr     int_all F~L F~L' = int_all FL FL' + sum CkL sig_kk'l Ck'L' (20)
Cr
Cr  where
Cr
Cr     sig_kk'l = int_MT (P~kL P~k'L' - PkL Pk'L')              (21)
Cr
Cr  is independent of the shape of any particular augmentation
Cr  function --- thus according to (20) one can compute sig once and
Cr  for all for a set of (kL) pairs, and evaluate the augmentation `on
Cr  the fly' when the coefficients C are known.  There is a
Cr  corresponding expression for kinetic energy and a (slightly more
Cr  complicated) form for pot. matrix elements pi:
Cr
Cr     int F~L V F~L'  = int FL V0~ FL' + ppi_LL'
Cr
Cr  where
Cr
Cr     ppi_LL' = int (F~L V1 F~L' - F^L V2~ F^L')
Cr             + sum_M Q_LL'M * int [(V0~ - V2~) G_M]
Cr
Cr     Q_LL'M  = int ( F~L F~L' - F^L F^L') r^m Y_M
Cr
Cr  As before, we can further approximate Fi~ and Fi^ as in Eq. 16; then
Cr  each of these terms can be expressed as a linear combination of
Cr  matrix elements ppi_kk'LL' independent of the shape of any
Cr  particular augmentation function:
Cr
Cr     ppi_LL' = sum_kk' C*_kL C_k'L'  ppi_kk'LL'
Cr
Cr  where
Cr
Cr     ppi_kk'LL' = int (P~kL V1 P~k'L' - PkL V2~ Pk'L')
Cr                + sum_M Q_kk'LL'M * int [(V0~ - V2~) G_M]
Cr
Cr     Q_kk'LL'M  = int ( P~kL P~k'L' -  PkL Pk'L') r^m Y_M
Cr
Cr  In the Springer book the augmentation for all spheres is presumed
Cr  to proceed along the lines of (20).  For numerical reasons it has
Cr  been found necessary to separate out heads from tails and compute
Cr  integrals of involving heads directly.  This is because for the
Cr  heads the P_kL expansion must be expanded to very high k, so
Cr  high that the expansion sometimes has numerical problems with
Cr  convergence.  Only tail-tail matrix elements are computed according
Cr  to (20).  For 1-center (head-head) augmentation integrals, the
Cr  integral (Fi~ Fj~ - Fi Fj) is explicitly computed; head-tail
Cr  integrals are mixed explicit repsn in one index and polynomial
Cr  expansions in the other; the (original) poly-poly integrals are
Cr  kept and used for the 3-center expansions.  Consequently, three
Cr  kinds of integrals are kept.
Cr           sigkk             sighk             sighh
Cr         P~P~ - PP        H~P~ - HP        H~H~ - HH
Cr  with corresponding integrals tau and ppi.
Cr
Cr  In this version, augmented functions are linear combinations
Cr  of radial wave functions u and s that match continuously and
Cr  differentiably onto the envelope functions at rmt.
Cr  u and s are linear combinations of and phi,phidot defined as:
Cr  u has val=1, slo=1 at rmax, s has val=0, slo=1 at rmax.
Cr
Cr  Local orbitals : local orbitals occur in one of these types:
Cr  1. val,slo=0 at rmt                         (10s digit pnz=0)
Cr  2. val,slo matches sm hankel                (10s digit pnz=1)
Cr  3. val,slo matches sm hankel, perturbative  (10s digit pnz=2)
Cr  In all cases tails are not expanded about other sites.  Thus local
Cr  orbitals involve only one- and two-center terms.
Cr
Cr  Consider the overlap matrix where Fi is a local orbital:
Cr    int_all F~i F~j  =  int_all Fi Fj  + int_MT (F~i F~j - F^i F^j)
Cr
Cr  Case 1: Fi = F~i = F^i.  Then no error occurs in the L truncation
Cr  of F~j because
Cr     int_MT Fi Fj - F^i F^j = int_MT Fi (Fj-F^j) = 0
Cr  This is because the second factor is zero for l<=lmxa, and the
Cr  first factor is zero for l>lmxb.  This is similarly true for the
Cr  kinetic energy matrix elements.
Cr
Cr  Cases 2 and 3: The true envelope function and its one-center
Cr  expansion are identical. Also we have Fi=F~i for r<rmt, and Fi=F^i
Cr  for r>rmt; and the value and slope of F~i matches F^i at rmt.
Cr  There is a also exact cancellation in the overlap
Cr     int Fi Fj - F^i F^j = int_r<rmt F~i (Fj-F^j) +
Cr                           int_r>rmt F^i (Fj-F^j)
Cr  In each term, the second factor is zero for l<=lmxa, and the first
Cr  factor is zero for l>lmxb.
Cr
Cr  For the potential matrix element there is an additional
Cr  complication and an approximation.  Let F~i be a local orbital
Cr  with angular momentum L and consider a one-center expansion of the
Cr  matrix element a one-center expansion of the matrix element
Cr
Cr     int F~i V F~j  = F~i sum_L'' V_L'' sum_L' F~L'
Cr
Cr  If F~j has a finite Y_L expansion at this site, then the integral
Cr  is again exact provided that L'' is expanded at least up to the
Cr  difference in L's between F~i F~j (the l-cutoff for V is lmxl),
Cr  because F~i as only a single L in its Ylm expansion; therefore
Cr  only the L projection of the product (sum_L'' V_L'' sum_L' F~L')
Cr  makes any contribution to the integral.
Cr
Cr  But when F~j has higher L components (as it does when it is
Cr  centered at another site) but is truncated in its L expansion,
Cr  there is an error because there is a missing projection onto L of
Cr  the following:
Cr
Cr    sum_L'' sum_L'>La  V_L'' F~L' ~ sum_L'>La V_(L-L') F~L'
Cr
Cr  Since both factors are small (V_(L-L') is small for L-L'>0, and
Cr  F~L' is supposed to be small for L'>La), and moreover the product
Cr  V_(L-L') F~L is largest near the MT radius, at which point the
Cr  local orbital's value and slope are both zero, the integral of a
Cr  local orbital with this product clearly is very rapidly convergent
Cr  in its L-cutoff.
Cr
Cr  The complication arises that, if no smoothed analog is to be
Cr  computed, it is no longer true that the local representation of
Cr  the true potential can be shifted by an arbitrary harmonic
Cr  function r^m Y_M, as in the usual 'three-fold' representation of
Cr  the potential.  Therefore, in considering these matrix elements,
Cr  the local expansion of V must be written as
Cr
Cr    V(r) = sum_L V_L(r) Y_L,   V_L(r) = V1_l(r) + VVAL_L/R^l r^l
Cr
Cr  where the VVAL_L is the L projection of the e.s. potential at the
Cr  MT radius R.
Cr
Cr  *Documentation of hab,vab,sab and sig,tau,ppi :
Cr   hab,vab,sab are generated in potpus and are matrix elements of the
Cr   triplet of wave functions (u,s,gz) in the spherical part of V(r)
Cr   u and s are linear combinations of and phi,phidot defined as:
Cr     u has val=1, slo=1 at rmax;   s has val=0, slo=1
Cr   There may additionally be local orbitals gz, specified by nonzero
Cr   pnz, of one of two types, as described in potpus.f
Cr   Orbitals of the first type are confined to the augmentation sphere;
Cr   they have no corresponding smooth part.  Orbitals of the second type
Cr   extend into the interstitial by attaching a Hankel tail.
Cr
Cr  *Structure of sig,tau,ppi, and hab,vab,sab.
Cr   All these arrays contain matrix elements between combinations
Cr   of valence states and local orbitals and are linear combinations of
Cr      int (ul,sl,gz,r*P_k)_l {h,v,1} (ul,sl,gz,r*P_k)_l'
Cr
Cr   There is some duplication in the two sets of matrix elements,
Cr   which can be confusing.  The following points help to explain
Cr   their relationship.
Cr
Cr  *{h,v,s}ab are matrix elements of products of true wave functions
Cr   in the spherical V only, in the (ul,sl,gz) form discussed above.
Cr
Cr  *sig,tau,ppi are matrix elements of products of augmented wave
Cr   functions, with the smooth part subtracted off.  ppi contains
Cr   matrix elements of the full potential, not merely the spherical
Cr   part. {h,v,s}ab are used to help assemble sig,tau,ppi.
Cr
Cr  *With the exception of ppi, all matrix elements are diagonal in l,
Cr   and independent of m; require only one l index.
Cr   In the case of ppi, the full LL' matrix is required.
Cr   Contributions from the nonspherical part of V are calculated with
Cr   the large component of the radial wave functions only.
Cr
Cr  *Matrix elements with a local orbital of the first type do not have
Cr   a smooth contribution subtracted.
Cr
Cr  *sig,tau,ppi come in three flavors as discussed above, e.g.
Cr   H~H~-HH; H~P~-HP; P~P~-PP, there is a separate array for each:
Cr     sigkk(Pk,Pk) is dimensioned (1+kmax,1+kmax)
Cr     sighk(H, Pk) is dimensioned (nkaph,1+kmax)
Cr     sighh(H, H)  is dimensioned (nkaph,nkaph)
Cr   Here nkapi is 1 or 2 depending on whether there are one
Cr   or two kinds of envelope functions in the basis.
Cr
Cr   Thus the sig matrix elements are dimensioned
Cr     sig{1,2,3}(nf1,nf2,0..lmax,1..nsp)
Cr   with nf{1,2} = nkaph or 1+kmax, depending on which sig
Cr   and lmax = lmxa for those involving P's and lmxh for those not.
Cr
Cr * Documentation of LDA+U
Cr   the orbital dependent potential is of the form Vmm' for each LDA+U
Cr   block (fixed l and site), which makes matrix elements
Cr   <phi_m|Vnonlocal|phi_m'>
Cr   To add to ppi, these need to be first rotated to (u,s) basis which is
Cr   done in vlm2us whose output is vumm array for this site (possibly
Cr   for different l's). vumm is passed to gaugm.
Cr
Cr * Noncollinear extension
Cr   As a first cut at a noncollinear extension, the present
Cr   implementation makes the rigid spin approximation: that is the
Cr   quantization axis is constrained to a single vector inside the
Cr   unit sphere.  On this local axis the potential and moments are
Cr   diagonal in spin.  Implementation is particularly simple because
Cr   the matrix elements are diagonal in spin, as in the collinear
Cr   case.  At the end they they need to be rotated to the global
Cr   quantization axis.
Cr
Cr   A full noncollinear extension would be rather complicated.  A
Cr   simple generalization would be to extend the potential vum to a
Cr   noncollinear form; the difference in the (locally) collinear and
Cr   noncollinear can be treated as a perturbation, analogous to the
Cr   (spatially) nonspherical parts of the potential.  The
Cr   perturbation is given by the difference v1-v0/Y0, for both spin
Cr   and spatial components.
Cr
Cr   Since the basis set is constructed for a (locally) collinear
Cr   density, moments qum remain spin diagonal (I think).  They are
Cr   used only in a term contributing to the Hartree potential, Eq. 27
Cr   in the book chapter.  Without the rigid spin approximation it
Cr   needs to be established that term qm * (gpot0-gpotb) be correct
Cr   for each spin.  For l=0 it should be true since the magnetic axis
Cr   is chosen so the charge is diagonal in spin.  Not sure about
Cr   l>0.
Cu Updates
Cu   24 Aug 15 (Scott McKecknie) SO weighted average of electric field
Cu   22 Jun 14 augmat generates noncollinear optical matrix elements
Cu   11 Feb 14 Some redesign of optical matrix elements for general case
Cu   13 Jan 14 (Ben Kaube) start optical matrix elements for PkL
Cu   15 Nov 13 (Ben Kaube) returns radial part of optical matrix elements
Cu   05 Jul 13 tso was held in ppiz(...,2), now a separate array
Cu   24 Jun 13 New facso
Cu   01 Apr 10 Add external B field to potential
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   27 Apr 05  LDA+U (Lambrecht)
Cu   24 Dec 04 (A. Chantis) ppi matrix elements for full L.S
Cu    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   15 Jul 04 (Chantis) radial integrals for spin-orbit coupling
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   27 Aug 01 Extended to local orbitals.  Altered argument list.
Cu   20 Feb 01 Added ppnl to potential parameters generated
Cu   13 Jun 00 spin polarized
Cu   17 May 00 Adapted from nfp augmats.f
C ----------------------------------------------------------------------
      use structures
      implicit none
      integer, parameter :: n0=10,nkap0=4,nppn=12,nab=9
C ... Passed parameters
      integer, intent(in) :: lmxa, kmax, kmxax, nlml, nr, nsp, nspc
     . , nkaph, nkapi, lmxh, lx, lso, lbf, lcplxp, nbf, loptic
     . , lmaxu, lldau, iblu, idu(4), nrg, nlg
      complex(8), intent(in) :: vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      real(8), intent(in) :: z,rmt,rsma,a,facso,rs3,vmtz
      integer, intent(in) :: jcg(*),indxcg(*),lh(nkap0)
      real(8), intent(in) :: v0(nr,nsp), v1(nr,nlml,nsp), v2(nr,nlml,nsp)
     . , cg(*), pnu(n0,nsp), pnz(n0,nsp), eh(n0,nkaph), rsmh(n0,nkaph)
     . , ehl(n0), rsml(n0)
      real(8), intent(inout) :: rofi(nr), rwgt(nr), bloc(nbf,3)
     . , ppnl(nppn,n0,2), gpot0(nlml), gpotb(nlml), hab(nab,n0,nsp)
     . , vab(nab,n0,nsp), sab(nab,n0,nsp)
      real(8), intent(inout) :: rgrad(nrg,nrg,2,nlg,nsp*nspc)
      real(8), intent(inout) :: rgrade(-nkaph:kmxax,-nkaph:kmxax,2,nlg,nsp*nspc)
      real(8) :: rmee(1)
      character(len=*), intent(inout) :: buf(*)
      type(str_site), intent(in) ::  s_site
C ... Dynamically allocated local arrays
      real(8), pointer :: sig(:,:),tau(:,:),ppi(:,:),tso(:,:)
      real(8), allocatable :: fhp(:,:,:), grv1(:,:)
C ... Local parameters
      integer k,ll,lmxl,nlma,nlmh,i
      double precision pi,y0,pp(n0,2,5)
      integer lmxx
      parameter (lmxx=25)
      integer lxa(0:kmax),lsol
      double precision vdif(nr*nsp),sodb(nab,n0,nsp,2,2),
     .  qum((lmxa+1)**2*(lmxx+1)*6*nsp),vum((lmxa+1)**2*nlml*6*nsp),
     .  fh(nr,0:lmxh,nkap0),xh(nr*(lmxh+1)*nkap0),
     .  vh((lmxh+1)*nkap0),fp(nr,0:lmxa,0:kmax),
     .  dh((lmxh+1)*nkap0),xp(nr,0:lmxa,0:kmax),
     .  vp((lmxa+1)*(kmax+1)),dp((lmxa+1)*(kmax+1))
!     .  grv1(nr,3)
      double complex vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)

      call tcn('augmat')
      pi   = 4d0*datan(1d0)
      y0   = 1d0/dsqrt(4d0*pi)
      nlma = (lmxa+1)**2
      lmxl = ll(nlml)
      if (lmxl > lmxx) call rxi('augmat: increase lmxx to',lmxl)
      lsol = mod(lso,10)
      if (lso == 0) allocate(tso(1,1))


C --- Make hab,vab,sab and potential parameters pp ---
      do  i = 1, nsp
        call dpcopy(v1(1,1,i),vdif(1+nr*(i-1)),1,nr,y0)
        call dpadd(vdif(1+nr*(i-1)),v0(1,i),1,nr,-1d0)
      enddo

C --- SCOTT: get gradient of v1 for SO-weighted electric field ---
      if (lso >= 10) then
        allocate(grv1(nr,3))
        call gradv1(nlml, z, nr, nsp, nspc, rofi, v1, grv1)
      else
        allocate(grv1(1,1))
      endif

C     NB:rofi MUST be dimensioned rofi(1..nrx) if V is to be extrapol.
      call potpus(z,rmt,lmxa,v0,vdif,grv1,a,nr,nsp,nspc,lso+lbf,loptic,
     .  nrg,nlg,facso,rofi,pnu,pnz,ehl,rsml,rs3,vmtz,nab,n0,pp,ppnl,hab,
     .  vab,sab,sodb,rgrad,buf)

C --- Moments and potential integrals of ul*ul, ul*sl, sl*sl ---
      call momusl(z,rmt,lmxa,pnu,pnz,rsml,ehl,lmxl,nlml,a,nr,nsp,rofi,
     .  rwgt,v0,v1,qum,vum)

C --- Set up all radial head and tail functions, and their BC's ---
      do  k = 0, kmax
        lxa(k) = lmxa
      enddo
C     Debugging check; see description of rsmh above
      do  k = 1, lmxh+1
        if (pnz(k,1) /= 0 .and. mod(pnz(k,1),100d0) < 10) then
          if (rsmh(k,nkapi+1) /= 0)
     .    call rx1('augmat: illegal value for rsmh',rsmh(k,nkapi+1))
        endif
      enddo

      nlmh = (lmxh+1)**2
      call dpzero(fh,nr*(lmxh+1)*nkap0)
      call fradhd(nkaph,eh,rsmh,lh,lmxh,nr,rofi,fh,xh,vh,dh)
      call fradpk(kmax,rsma,lmxa,nr,rofi,fp,xp,vp,dp)

C ... LDA+U: rotate vorb from (phi,phidot) to (u,s) for all l with U
C     at this site and store in vumm
      if (lldau /= 0) then
        call vlm2us(lmaxu,rmt,idu,lmxa,iblu,vorb,ppnl,sab,vumm)
      endif
C ... Optics: matrix elements of envelope fns + PkL
      if (loptic > 0) then
        allocate(fhp(nr,0:lmxa,-nkaph:kmxax));call dpzero(fhp,size(fhp))
        do  k = 1, nkaph  ! Head functions
          fhp(:,0:lmxh,-k) = fh(:,0:lmxh,k)
        enddo
        do  k = 0, kmax  ! Polynomials
          fhp(:,:,k) = fp(:,:,k)
        enddo
        call rgrme(0,kmxax+1+nkaph,kmax+1+nkaph,kmax+1+nkaph,1,1,nlg-1,
     .    lmxa,nr,rofi,rwgt,fhp,fhp,rgrade)
        deallocate(fhp)
        if (nsp == 2) then
          k = size(rgrade(:,:,:,:,1))
          do  i = 2, nsp*nspc
            call dpscop(rgrade,rgrade,k,1,1+k*(i-1),1d0)
          enddo
        endif
      endif

C --- Augmentation matrices for cases P*P, H*H, H*P ---

C ... Hsm*Hsm
      sig => s_site%sighh
      tau => s_site%tauhh
      ppi => s_site%pihh
      if (lx /= 0) then
        sig => s_site%sighhx
        tau => s_site%tauhhx
        ppi => s_site%pihhx
      endif
      if (lsol /= 0) tso => s_site%sohh
      call gaugm(nr,nsp,lsol,lbf,lcplxp,rofi,rwgt,lmxa,lmxl,nlml,v2,
     .  gpotb,gpot0,hab,vab,sab,sodb,qum,vum,cg,jcg,indxcg,
     .  nkaph,nkapi,lmxh,lh,fh,xh,vh,dh,
     .  nkaph,nkapi,lmxh,lh,fh,xh,vh,dh,
     .  lmxh,sig,tau,nlmh,nlmh,ppi,tso,
     .  lmaxu,vumm,lldau,idu,nbf,bloc)

C ... Pkl*Pkl
      sig => s_site%sigkk
      tau => s_site%taukk
      ppi => s_site%pikk
      if (lx /= 0) then
        sig => s_site%sigkkx
        tau => s_site%taukkx
        ppi => s_site%pikkx
      endif
      if (lsol /= 0) tso => s_site%sokk
      call gaugm(nr,nsp,lsol,lbf,lcplxp,rofi,rwgt,lmxa,lmxl,nlml,v2,
     .  gpotb,gpot0,hab,vab,sab,sodb,qum,vum,cg,jcg,indxcg,
     .  kmax+1,kmax+1,lmxa,lxa,fp,xp,vp,dp,
     .  kmax+1,kmax+1,lmxa,lxa,fp,xp,vp,dp,
     .  lmxa,sig,tau,nlma,nlma,ppi,tso,
     .  lmaxu,vumm,lldau,idu,nbf,bloc)

C ... Hsm*Pkl
      sig => s_site%sighk
      tau => s_site%tauhk
      ppi => s_site%pihk
      if (lx /= 0) then
        sig => s_site%sighkx
        tau => s_site%tauhkx
        ppi => s_site%pihkx
      endif
      if (lsol /= 0) tso => s_site%sohk
      call gaugm(nr,nsp,lsol,lbf,lcplxp,rofi,rwgt,lmxa,lmxl,nlml,v2,
     .  gpotb,gpot0,hab,vab,sab,sodb,qum,vum,cg,jcg,indxcg,
     .  nkaph,nkapi,lmxh,lh,fh,xh,vh,dh,
     .  kmax+1,kmax+1,lmxa,lxa,fp,xp,vp,dp,
     .  lmxh,sig,tau,nlmh,nlma,ppi,tso,
     .  lmaxu,vumm,lldau,idu,nbf,bloc)

      if (lsol == 0) deallocate(tso)
      deallocate(grv1)

      call tcx('augmat')
      end subroutine augmat
