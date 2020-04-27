      subroutine prodbasa(mode,s_lat,s_site,s_spec,ig,nl,nsp,nbas,nat,ndrphi,
     .  nrphiv,nrphic,nrpb,ntpba,ndlmto,rprbme,prbasme)
C- Matrix elements of partial waves and product basis for all augmented sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:symgr
Cio    Passed to:  makecgr
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt nr
Cio    Passed to:  *
Ci Inputs
Ci  mode   :1s digit
Ci         :0 exclude all states core states
Ci         :1 include core states
Ci         :2 include valence states
Ci         :3 combination of 1+2
Ci         :10s digit : effective dimension for prbasme
Ci         :0 Kotani convention: prbasme =  prbasme(ndlmto,ndlmto,nsp,ntpba(0),iat)
Ci         :  ntpba(0) must not be < max PB dimension for one site
Ci         :1 functions from different sites are strung together in a long chain
Ci         :  ntpba(iat) is the PB dimension for site iat
Ci  ig     :index to group operation
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  nbas   :size of basis
Ci  nat    :number of sites in basis with augmentation spheres
Ci  ndrphi :(dimensioning) global max number of core + valence radial functions of a particular l
Ci         :Formerly named nn in old GW code
Ci  nrphiv :number of (valence partial waves) for a particular l and site
Ci         :Formerly named nindxv in old GW code
Ci  nrphic :number of (core partial waves) for a particular l and site
Ci         :Formerly named nindxc in old GW code
Ci  ndrpb  :(dimensioning) global max number of radial product basis functions of a particular l
Ci  nrpb   :cumulative number of radial product basis functions for a particular l and site.
Ci         :nrpb is the analog of nrphiv, but for the product basis and it is cumulative.
Ci         :Functions nrpb(l,iat):nrpb(l+1,iat) are the family of radial B functions
Ci         :stored in rprodb for quantum number l and site iat.
Ci         :In former GW code, nrpb(l+1,iat)-nrpb(l,iat) was stored in nxx
Ci  ntpba  :ntpba(iat) = number of full product basis functions at site iat.
Ci         :ntpba formerly named nblocha or mdim in old GW code
Ci  ndlmto :(dimensioning) global max no. partial waves (l,n,m) in LMTO basis at any site
Ci  rprbme :radial matrix elements < gtoto gtoto B> for all sites
Co Outputs
Co  prbasme:full matrix elements CG * rprbme for all sites
Cl Local variables
Cr Remarks
Cr   States are ordered by m,n,l,core followed by m,n,l,val
Cr   This adopts the sequence of the old GW code; see mtoindex
Cu Updates
Cu   23 Aug 18 Adapted from old GW ppbafp_v2
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Passed parameters
      integer mode,nl,nsp,nbas,nat,ndrphi,ig,ndlmto
      integer :: nrphiv(0:nl-1,nat),nrphic(0:nl-1,nat),ntpba(0:nat)
      integer :: nrpb(0:2*(nl-1),nat+1) !dimensioned as a single vector to avoid compiler complaints
      real(8), intent(in) :: rprbme(0:nl-1,ndrphi,0:nl-1,ndrphi,0:2*(nl-1),*)
      real(8), intent(out) :: prbasme(ndlmto,ndlmto,*)
C ... Dynamically allocated arrays
      real(8), allocatable :: cgr(:,:,:)
C ... Local parameters
      logical addcore,addval
      integer,parameter:: NULLI = -99999
      integer indxlm(3,ndrphi*nl*nl),mtoindx(ndrphi,nl*nl),nblr(0:2*(nl-1))
      integer ib,iat,is,isp,mxrbli,ntpbi,nmtoi,offpb,offpr,npbi

      call info5(20,1,0,
     .  ' prodbas:  make prbasme = (phi phi | full PB) for %i sites   group op %i   dimen (%i,%-1j%i,%ix%i)',
     .  nat,ig,ndlmto,ntpba,nat)

      if (nsp == 2) call rx('prodbas not spin polarized yet')

      allocate(cgr(nl**2,nl**2,(2*nl-1)**2))
      call makecgr(s_lat,nl-1,s_lat%symgr(1,ig),cgr)

      addcore = mod(mod(mode,10),2) == 1
      addval  = mod(mod(mode,10)/2,2) == 1
      iat = 0; offpr = 0; offpb = 0; isp = 1
      do  ib = 1, nbas
        is = s_site(ib)%spec
        if (s_spec(is)%lmxa < 0) cycle
        iat = iat+1
        nmtoi = 0; npbi = 0
        call nrpb2nbl(03,nl,iat,nrpb,nblr,mxrbli,ntpbi) ! nblr,mxrbli,ntpbi needed for prodbasi

        if (addcore) call mtoindex(nl,nrphic(0,iat),-1,ndrphi,indxlm,mtoindx,nmtoi)
        if (addval)  call mtoindex(nl,nrphiv(0,iat),nrphic(0,iat),ndrphi,indxlm,mtoindx,nmtoi)
        if (nmtoi > ndrphi*nl*nl) call rx('dimensioning error in prodbasa')

        call prodbasi(nl,ndrphi,nmtoi,ndlmto,ntpbi,nblr,indxlm,cgr,rprbme(0,1,0,1,0,1+offpr),
     .    prbasme(1,1,1+offpb))
C       print *, 'sumcheck prodbas',iat,sngl(sum(prbasme(:,:,offpb+1:offpb+ntpbi)))
        call info2(30,0,0,' sumcheck prodbas %,4i %;10F',iat,sum(prbasme(:,:,offpb+1:offpb+ntpbi)))
        offpr = offpr + mxrbli
        offpb = offpb + ntpba(0)
        if (mod(mode/10,2) /= 0) offpb = offpb + ntpba(iat) - ntpba(0)
      enddo

      deallocate(cgr)

      end
      subroutine prodbasi(nl,ndrphi,nmtoi,ndlmto,ntpbi,nblr,indxlm,cgr,rprbme,prbasme)
C- Assemble full product basis for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci  nl     :(maximum l) + 1
Ci  ndrphi  :max number of core + valence radial functions of a particular l (for dimensioning)
Ci         :Formerly named nn in old GW code
Ci  nmtoi  :number of partial waves with unique l,n,m for this site
Ci  ndlmto :dimensions prbasme: must not be < nmtoi
Ci  ntpbi  :Number of product basis functions (used for sanity checks)
Ci  nblr   :number of radial product basis functions for each l
Ci         :Formerly named nxx in old GW code
Ci  indxlm :for each partial wave, index to n,l,m quantum number (mtoindex)
Ci         :indxlm(1) = index to which radial function
Ci         :indxlm(2) = index to l quantum number
Ci         :indxlm(3) = index to m quantum number
Ci         :            index to L = l*l + l + m + 1
Ci  cgr    :Clebsch-Gordan coefficients <lm3 | lm1 lm2>, possibly rotated (makecgr)
Ci  rprbme :radial part of matrix elements <gtoto gtoto | B>
Co Outputs
Co  prbasme:full matrix elements CG * rprbme
Co         :prbasme(j,k,i) = <phi(j) phi(k) B(i)>, with j,k,i each a compound l,n,m index
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   23 Aug 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,ndrphi,nmtoi,ndlmto,ntpbi,nblr(0:2*(nl-1)),indxlm(3,nmtoi)
      real(8), intent(in) :: cgr(nl**2,nl**2,(2*nl-1)**2)
      double precision rprbme(0:nl-1,ndrphi,0:nl-1,ndrphi,0:2*(nl-1),maxval(nblr))
      real(8), intent(out) :: prbasme(ndlmto,ndlmto,ntpbi)
C ... Local parameters
      integer i,lb,nb,mb,lmb,i2,l2,n2,m2,lm2,i1,l1,n1,m1,lm1

C     For each PB matrix element i, do
      i = 0
      do  lb = 0, 2*(nl-1)
        do  nb = 1, nblr(lb)
          do  mb  = -lb, lb
            i    = i+1 ! Current product basis function
            lmb  = lb*lb + lb + mb + 1

C           Loop over unique partial wave pairs (l2,n2,mp; l1,n1,m1) and spin
            do  i2 = 1, nmtoi
              l2  = indxlm(1,i2)
              n2  = indxlm(2,i2)
              m2  = indxlm(3,i2)
              lm2 = l2*l2 + l2 + m2 + 1

              do  i1 = 1, nmtoi
                l1  = indxlm(1,i1)
                n1  = indxlm(2,i1)
                m1  = indxlm(3,i1)
                lm1 = l1*l1 + l1 + m1 + 1

                prbasme(i1,i2,i) = cgr(lm1,lm2,lmb)*rprbme(l1,n1,l2,n2,lb,nb)
              enddo
C             print *, i2,i,sngl(sum(prbasme(:,i2,1,i)))
            enddo
C           print *, 'done',i,sngl(sum(prbasme(:,:,i,1)))
          enddo
        enddo
      enddo

      if (i /= ntpbi) call rx('prodbasi: dimensioning mismatch')

      end

