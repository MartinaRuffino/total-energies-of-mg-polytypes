      subroutine strg1c(ivl,nkap,nds,nlmy,nl,iat,ixi,ips,rsm,c01,ehl,
     .  rl,kmx,rsma,alat,plat,pos,iax,ntab,ntabg,cy,cg,indxcg,jcg,s,
     .  slj,dslj,lslj,dlslj)
C- 1c decomposition of screened basis functions and their Laplacians
C  via the polynomial expansion (for testing)
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :Identifies the functions used in the screened basis
Ci         :ivl=0,1,2 use double-kappa functions F_0 and F_1
Ci         :Other ivl correspond to 1-kappa functions F_1
Co         : ivl   F1 function
Co         :  0    Gaussian G1
Co         :  1    sm. Hankel, Hsm
Co         :  2    energy derivative of Hsm
Co         :  3    Hsm + c01(0)*G0 + c01(1)*G1
Co         :  4    Hsm - Hsm(rsm->0)
Co         :  5    Hsm + c01(0)*G0 + c01(1)*G1 - Hsm(rsm->0)
Co         : ivl   F0 function
Co         :  0    Gaussian G0
Co         :  1    Gaussian G0
Co         :  2    Gaussian G10
Co         : else  not used
Ci   nkap  :number of functions per channel entering into strux
Ci   nds   :leading dimensions of s
Ci   nlmy  :leading dimensions of slj and lslj
Ci   nl    :Dimensions rsm,rl,ehl,c01; must be at least max lmxh+1
Ci   iat   :site R at which the screened basis functions are defined
Ci   ixi   :Indicates site R' for one-center expansion.  ixi is index
Ci         :to pair table relative to start of cluster at iat
Ci         :Thus ixi=1 points to head; ixi=2 points to 1st neighbor, etc
Ci   ips   :index to which species each site belongs
Ci   rsm   :Gaussian/Hankel smoothing radii for each species and
Ci          angular momentum
Co   c01   :Gaussian coefficients for generalized envelope functions (ivl>2)
Ci   ehl   :Energies for sm. Hankels where they enter into envelope functions
Ci         :Not referenced if ivl = 0
Ci   rl    :l-dependent radius at which to evaluate the expansion
Ci   kmx   :kmx(is) is kmax for decomposition of Gaussians or sm. Hankels
Ci         :into polynomials P_kL at a site with species index is
Ci   rsma  :P_kL smoothing radius (one for each species)
Ci   alat  :lattice constant
Ci   plat  :primitive lattice translation vectors in units of alat
Ci   pos   :basis vectors
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairc.f)
Ci   ntabg :ntabg(ib) no. of neighbors in cluster centered at ib (pairg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordan coefficients
Ci   s     :Structure constant matrix defining the screened basis; see Remarks
Co Outputs
Co   slj   :slj(L',L,i) is L'-component of 1c expansion of U_iRL around site ixi
Co          i = 0 and i = 1 correspond to value and Laplacian parts
Co          of the value-Laplacian unit basis derived for site iat
Co   dslj  :radial derivative of slj
Co   lslj  :lslj(L',L,i) = \lap slj(L',L,i)
Co  dlslj  :radial derivative of lslj
Cl Local variables
Cl  gvl(L',L,p) :L'th component of F_0 centered at ph; F_0 defined by ivl, above
Cl              :p=0 is function; p=1 is Laplacian
Cl  hvl(L',L,p) :L'th component of F_1 centered at ph; F_1 defined by ivl, above
Cl              :p=0 is function; p=1 is Laplacian
Cl  nlmb   :Number of L channels at iat = R
Cl  nlmix  :Number of L channels at R"
Cl  nlma   :Number of L channels at R'
Cb Bugs
Cb   need to pass lmxcg to check if CG arrays are large enough
Cr Remarks
Cr   Functions ivl = 0,1,2, with nkap=2
Cr   Screened basis U_iRL is the unit value-Laplacian set of functions having
Cr   the property that their values/laplacians at sphere R' are zero except
Cr   when indices match:
Cr     \lap^i' U_iRL(R'L') = \delta_ii' * \delta_LL' * \delta_RR'
Cr
Cr   U_iRL are constructed from a double set of functions F_0, F_1 as
Cr      U_iRL = \sum_i"R"L" s_i"R"L",iRL * F_i"R"L"                  (*)
Cr   where the sum over i" runs over 0,1, R" and L" are sites and L channels
Cr   in a given cluster, and R is always the center of the cluster.
Cr   In the code, indices to s are stored as:
Cr      U_iRL = \sum_i"R"L" s(L",L,i",i,ix) * F_i"R"L"               (*)
Cr
Cr   F_0 are always the Gaussians G0, whereas F_1 are either G1, Hs, or Hs-dot
Cr   depending on switch ivl (see above).
Cr
Cr   Program evaluates U_iRL(R'L') as array slj, and \lap U_iRL(R'L') as
Cr   array lslj with R = iat and R' = ixi using (*) and 1c decomposition of F_i.
Cu Updates
Cu   06 Apr 15 Some bug fixes
Cu   09 Oct 11 Returns function gradients (WARNING: never checked!)
Cu             Change hcr(:,1:nspec) -> rl(:) -- radius just for aug. sphere.
Cu   25 Sep 11 Workaround GNU compiler bug
Cu   23 Sep 11 Extended to new kinds of envelope functions (see sole0g)
Cu   18 Dec 08 (S. Lozovoi) switch to whether to use dgemm
Cu   05 May 08 (S. Lozovoi) Hs and Hs-dot added
Cu   27 Feb 08 (S. Lozovoi) First written
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ivl,iat,nds,nlmy,nl,ixi,nkap,niax
      parameter (niax=10)
      integer iax(niax,*),ntab(*),ntabg(*)
      integer ips(*),kmx(*)
      double precision alat,plat(3,3),pos(3,*),c01(0:1,0:nl-1,*)
      double precision rl(0:nl-1),rsm(nl,*),ehl(nl,*),rsma(*)
      double precision s(nds,nds,nkap,nkap,*)
      double precision slj(nlmy,nlmy,nkap),lslj(nlmy,nlmy,nkap)
      double precision dslj(nlmy,nlmy,nkap),dlslj(nlmy,nlmy,nkap)
      integer indxcg(*),jcg(*)
      double precision cy(*),cg(*)
C ... Dynamically allocated local arrays
      real(8), pointer :: e0(:,:,:)
C ... Local parameters
      integer ib,ikap,il,isg,isi,isj,it,itt,jb,jkap,jvl,kmax,lix,ll,
     .  nclus,nlma,nlmb,nlmix,offR,isum
      integer n0,nclusm,ndim0,pmx
      parameter (n0=10, nclusm=200, ndim0=(n0+1)**2, pmx=1)
      double precision pclus(3,nclusm)
      double precision rsmix(0:n0),eix(0:n0)
C     25 Sep 11 dimension gvl,hvl 1:pmx+1 to work around GNU compiler bug
      real(8),target:: gvl(ndim0,ndim0,1:pmx+1),hvl(ndim0,ndim0,1:pmx+1)
      real(8),target:: dgvl(ndim0,ndim0,1:pmx+1),
     .                 dhvl(ndim0,ndim0,1:pmx+1)

      nclus = ntabg(iat)
C ... offR and isg are the offsets to iax table and s matrix, respectively
      offR = ntab(iat)
      isg = 0
      if (iat > 1) isg = isum(iat-1,ntabg,1)

C --- Checks ---
      if (nds > n0**2) call rxi('strg1c: nds is bigger than n0**2. nds = ',nds)
      if (nclus <= 0) call rxi('strg1c: empty cluster encountered for iat = ',iat)
      if (nclus > ntab(iat+1)-ntab(iat))
     .  call rxi('strg1c: size of the cluster is too big. nclus = ',nclus)
      if (nclus > nclusm)
     .  call rxi('strg1c: size of the cluster exceeds the maximum. nclus = ',nclus)

C --- ixi-independent setup ---
C ... Head of the cluster
C     ic = iax(1,offR+1)
      nlmb = iax(9,offR+1)

c ... find coordinates of all atoms in the cluster
      do  it = 1, nclus
        call acoord(iat,it,alat,plat,pos,iax,ntab,pclus(1,it))
      enddo

c ... Expansion site
      ib = iax(2,offR+ixi)
      isi = ips(ib)
      nlma = iax(9,offR+ixi)
C     lh = ll(nlma)
      kmax = kmx(isi)

      call dpzero(slj,nkap*nlmy*nlmy)
      call dpzero(lslj,nkap*nlmy*nlmy)
      call dpzero(dslj,nkap*nlmy*nlmy)
      call dpzero(dlslj,nkap*nlmy*nlmy)

C --- Loop over all sites in the cluster ---
      do  it = 1, nclus
C       jb = site index to cluster atom, R" in Remarks
        itt = it+offR
        jb = iax(2,itt)
        isj = ips(jb)
C       Number of L channels from R"
        nlmix = iax(9,itt)
        lix = ll(nlmix)
        do  il = 0, lix
          rsmix(il) = rsm(il+1,isj)
        enddo
        if (ivl /= 0) then
          do  il = 0, lix
            eix(il) = ehl(il+1,isj)
          enddo
        endif

C   ... 1c-decomposition of F_0 and F_1 (nkap=2), or F_1 (nkap=1) around ixi
        jvl = ivl
        if (nkap == 2) jvl = jvl+10
        call gh1c(100+jvl,pclus(1,it),rsmix,nlmix,eix,c01(0,0,isj),pclus(1,ixi),
     .    rsma(isi),nlma,kmax,rl,cy,cg,indxcg,jcg,ndim0,gvl,dgvl,hvl,dhvl)

C   ... Combine with s to make slj and lslj
C       Index-site correspondence: iat --> L, ixi --> Lp, it --> ig
        do  ikap = 1, nkap
        do  jkap = 1, nkap
          if (jkap == 1) e0 => gvl
          if (jkap == 2 .or. nkap == 1) e0 => hvl
          call dgemm('N','N',nlma,nlmb,nlmix,1d0,e0(1,1,0+1),ndim0,
     .      s(1,1,jkap,ikap,isg+it),nds,1d0,slj(1,1,ikap),nlmy)
          call dgemm('N','N',nlma,nlmb,nlmix,1d0,e0(1,1,1+1),ndim0,
     .      s(1,1,jkap,ikap,isg+it),nds,1d0,lslj(1,1,ikap),nlmy)
          if (jkap == 1) e0 => dgvl
          if (jkap == 2 .or. nkap == 1) e0 => dhvl
          call dgemm('N','N',nlma,nlmb,nlmix,1d0,e0(1,1,0+1),ndim0,
     .      s(1,1,jkap,ikap,isg+it),nds,1d0,dslj(1,1,ikap),nlmy)
          call dgemm('N','N',nlma,nlmb,nlmix,1d0,e0(1,1,1+1),ndim0,
     .      s(1,1,jkap,ikap,isg+it),nds,1d0,dlslj(1,1,ikap),nlmy)
        enddo
        enddo

C ... End loop over cluster atoms
      enddo

      end
