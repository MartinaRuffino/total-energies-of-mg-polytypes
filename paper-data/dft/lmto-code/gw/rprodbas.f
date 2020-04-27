      subroutine rgprodbasa(mode,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,nrphi,lcuta,
     .  nocc,nunocc,tolbas,nrmx,gtota,
     .  nradpb,ndrpb,nrpb,ntpba,rprodb,nrgpb,ntgpba,rgprodb)
C- Generate radial part of product basis for all augmented sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa a rmt nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt nr
Cio    Passed to:  *
Ci Inputs
Ci  mode   :1s digit
Ci         : should be 0 for now.   Special modes not yet implemented
Ci         :10s digit.
Ci         : 0 return only dimensioning parameters nradpb, nrpb, ntpba
Ci         : 1 Turn ntpba into cumulative number of functions for site iat
Ci         : 2 Also return rprodb and/or grprodb
Ci         : These these ops may apply to p.b. or to grad p.b., or to both
Ci         : For gradients, nrpb -> nrgpb, ntpba -> ntgpba, rprodb -> rgprodb
Ci         :100s digit tells which kind of functions operations in 10s digit apply to
Ci         : 0 do nothing
Ci         : 1 regular product basis
Ci         : 2 grad+ = dg/dr - (l+1)/r*g  and grad+ = dg/dr + l/r*g
Ci         : 4 same as 2, but return (grad phi2) * phi1 instead of phi2 * (grad phi1)
Ci         : 6 grad dg(l)/r for g(l)   for (only one kind is made)
Ci         : 8 same as 6, but  (grad phi2) * phi1 instead of phi2 * (grad phi1)
Ci         : 1 may be taken in combination with 2, 4, 6, or 8
Ci         :1000s digit
Ci         : 1 check completeness of raw product basis (commented out for now)
Ci         : 2 check completeness of truncated product basis
Ci         : 4 check completeness of truncated product basis, verbose output
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  nbas   :size of basis
Ci  nat    :number of sites in basis with augmentation spheres
Ci  ndrphi :(dimensioning) max number of core + valence radial functions of a particular l
Ci         :Formerly named nn in old GW code
Ci  nrphi  :number of (core states + valence partial waves) for a particular l
Ci         :nrphi formerly named nindx in old GW code
Ci  lcuta  :exclude partial waves with l>lcuta from participating in prod. basis.
Ci  nocc   :noccv(l,ic) = 1 if to include as 1st fn in product basis, 0 otherwise
Ci  nunocc :nunoccv(l,ic) = 1 if to include as 2nd fn in product basis, 0 otherwise
Ci  tolbas :product basis tolerance, by l.  See Remarks in rprodbasi
Ci  nrmx   :leading dimension of gtota, rprodb
Ci         :nrmx formerly named nrofi in old GW code
Ci  gtota  :(r * valence partial waves), or combination of them and core states
Co Outputs
Co  nradpb :total number of radial product basis functions generated.  See Remarks
Co         :nradpb(1) for regular product basis functions
Co         :nradpb(2) for grad+ functions or grad functions
Co         :nradpb(3) for grad- functions
Co         :Which are returned depends on 100s digit mode, which see for meaning of grad+ and grad-
Co  ndrpb  :maximum number of radial functions for a given l at any site = maxval(nblr)
Co         :ndrpb(1) for regular product basis functions
Co         :ndrpb(2) for grad+ functions or grad functions
Co         :ndrpb(3) for grad- functions
Co         :Which are returned depends on 100s digit mode, which see for meaning of grad+ and grad-
Co  ntpba  :ntpba(iat) = number of full product basis functions at site iat.
Co         :ntpba formerly named nblocha or mdim in old GW code
Co  ntgpba :analog to ntpba, for gradient funnctions
Co  nrpb   :cumulative number of radial product basis functions for a particular l and site.
Co         :Functions nrpb(l,iat):nrpb(l+1,iat) are the family of radial B functions
Co         :stored in rprodb for quantum number l and site iat.
Co         :In former GW code, nrpb(l+1,iat)-nrpb(l,iat) was stored in nxx
Co  nrgpb  :analog to nrpb, for gradient functions
Co  rprodb :radial parts of orthonormalized product basis functions B for all l and sites
Co         :rprodb is constructed from r * (gtota1/r) * (gtota2/r)
Co         :Note that rprodb = r * true B(r).
Co  rgprodb:analog to rprodb, for gradient funnctions
Cl Local variables
Cl  off    :Offset to index in rprodb for current site
Cl         :rprodb is stored in rprodb(:,off+1).  If off<0, rprodb is not made
Cl  nblraw :number of radial functions for each l before reduction by orthogonalization
Cl  nblr   :number of radial functions for each l
Cl         :Formerly named nxx in old GW code
Cr Remarks
Cr   Calls rprodbasa for either product basis, gradient product basis, or combination
Cr   For the product basis, nradpb(1),ndrpb(1),nrpb,ntpba,rprodb are returned
Cr   The gradient product basis may contain two types of gradients
Cr   In that case nradpb(2;3),ndrpb(2;3),nrgpb,ntgpba,rgprodb are returned
Cr   For the simple gradient only one set is returned.
Cr
Cr   See rprodbasa for a detailed description of what is generated.
Cu Updates
Cu   20 Sep 18 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Passed parameters
      integer :: mode,nl,nsp,nbas,nat,ndrphi,nrmx,ndrpb(3),nradpb(3),lcuta(nat)
      integer :: nrphi(0:nl-1,nat),nocc(0:nl-1,ndrphi,nat),nunocc(0:nl-1,ndrphi,nat)
      integer, target :: nrpb(0:2*(nl-1),nat+1),ntpba(0:nat)
      integer, target :: nrgpb(0:2*(nl-1),nat+1,2),ntgpba(0:nat,2)
      real(8) :: rprodb(nrmx,*),rgprodb(nrmx,*)
      real(8) :: tolbas(0:2*(nl-1)),gtota(nrmx,0:nl-1,ndrphi,nsp,nat)
C ... Local parameters
      integer k,mode01,mode2,mode3,opt,offg

      mode01 = mod(mode,100)      ! Passed to rprodbasa as 1s+10s digit opt
      mode2 = mod(mode/100,10)    ! selects whether p.b. or gradient, or both
      mode3 = mod(mode/1000,10)   ! Passed to rprodbasa as 1000s digit opt
      if (mode2 == 0) return      ! No PB specified

      do  k = 1, 3
        if (k == 1 .and. mod(mode2,2) == 1) then ! Regular product basis
          opt = mode01 + 1000*mode3
          call rprodbasa(opt,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,nrphi,lcuta,
     .      nocc,nunocc,tolbas,nrmx,gtota,
     .      nradpb(1),ndrpb(1),nrpb,ntpba,rprodb)
        elseif (k > 1 .and. mode2/2 > 0) then ! gradient product basis
          if (mode2 >= 6 .and. k == 3) cycle ! simple gradient has only one type
          select case (mode2/2)
          case (1); opt = k     ! 2 if k==2;  3 if k==3
          case (2); opt = k + 4 ! Same as case 1, but make (grad phi) phi
          case (3); opt = 1     ! Simple gradient
          case (4); opt = 1 + 4 ! Same as case 3, but make (grad phi) phi
          end select
          opt = mode01 + 1000*mode3 + 100*opt
          offg = 0; if (k==3 .and. mod(mode/10,10) > 0) offg = nradpb(2)
          call rprodbasa(opt,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,nrphi,lcuta,
     .      nocc,nunocc,tolbas,nrmx,gtota,
     .      nradpb(k),ndrpb(k),nrgpb(0,1,k-1),ntgpba(0,k-1),rgprodb(1,1+offg))
        endif
      enddo

      end

      subroutine rprodbasa(opt,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,nrphi,lcuta,
     .  nocc,nunocc,tolbas,nrmx,gtota,nradpb,ndrpb,nrpb,ntpba,rprodb)
C- Generate radial part of product basis for all augmented sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa a rmt nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt nr
Cio    Passed to:  *
Ci Inputs
Ci  opt    :1s digit
Ci         : should be 0 for now.   Special modes not yet implemented
Ci         :10s digit
Ci         : 0 return only dimensioning parameters nradpb, nrpb, ntpba
Ci         : 1 Turn ntpba into cumulative number of functions for site iat
Ci         : 2 Also return rprodb
Ci         :100s digit
Ci         : 0 regular product basis; no gradient
Ci         : 1 Use dg(l)/r for g(l)   for (only one kind is made)
Ci         : 2 Use dg/dr - (l+1)/r*g  for gradient of the first kind (grad+)
Ci         : 3 Use dg/dr + l/r*g      for gradient of the second kind (grad-)
Ci         : Add 4 to return (grad phi2) * phi1 instead of phi2 * (grad phi1)
Ci         :1000s digit
Ci         : 1 check completeness of raw product basis (commented out for now)
Ci         : 2 check completeness of truncated product basis
Ci         : 4 check completeness of truncated product basis, verbose output
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  nbas   :size of basis
Ci  nat    :number of sites in basis with augmentation spheres
Ci  ndrphi  :max number of core + valence radial functions of a particular l (for dimensioning)
Ci         :Formerly named nn in old GW code
Ci  nrphi  :number of (core states + valence partial waves) for a particular l
Ci         :nrphi formerly named nindx in old GW code
Ci  lcuta  :exclude partial waves with l>lcuta from participating in prod. basis.
Ci  nocc   :noccv(l,ic) = 1 if to include as 1st fn in product basis, 0 otherwise
Ci  nunocc :nunoccv(l,ic) = 1 if to include as 2nd fn in product basis, 0 otherwise
Ci  tolbas :product basis tolerance, by l.  See Remarks in rprodbasi
Ci  nrmx   :leading dimension of gtota, rprodb
Ci         :nrmx formerly named nrofi in old GW code
Ci  gtota  :(r * valence partial waves), or combination of them and core states
Co Outputs
Co  nradpb :total number of radial product basis functions generated.  See Remarks
Co         :See 1000s digit opt for meaning of grad+ and grad-
Co  ndrpb :maximum number of radial functions for a given l at any site = maxval(nblr)
Co  mxnpba :maximum number of full product functions at any site = maxval(ntpba)
Co  nrpb   :cumulative number of radial product basis functions for a particular l and site
Co         :Functions nrpb(l,iat):nrpb(l+1,iat) are the family of radial B functions
Co         :stored in rprodb for quantum number l and site iat.
Co         :In former GW code, nrpb(l+1,iat)-nrpb(l,iat) was stored in nxx
Co  rprodb :radial parts of orthonormalized product basis functions B for all l and sites
Co         :rprodb is constructed from r * (gtota1/r) * (gtota2/r)
Co         :Note that rprodb = r * true B(r).
Cl Local variables
Cl  off    :Offset to index in rprodb for current site
Cl         :rprodb is stored in rprodb(:,off+1).  If off<0, rprodb is not made
Cl  ntpba  :ntpba(iat) = number of full product basis functions at site iat.
Cl         :One element was formerly named nblocha in old GW code
Cl  nblraw :number of radial functions for each l before reduction by orthogonalization
Cl  nblr   :number of radial functions for each l
Cl         :Formerly named nxx in old GW code
Cr Remarks
Cr   This routine calls rprodbasi for every site with an augmentation sphere.
Cr   For each site, all pairs of participating partial waves (see nocc, nunocc)
Cr   are initially created.  Their overlap matrix is computed and diagonalized.
Cr
Cr   The product basis is orthonormalized and functions with eigenvalues
Cr   of the overlap < tolbas are removed from the basis.
Cr
Cr   Remarks concerning gradients
Cr   A product basis constructed from r g2(l')/r grad g1(l)/r generates two kinds of radial functions:
Cr   because grad g1(l)/r Y_lm generates f_l+1 Y_l+1,m  with f proportional to dg/dr - (l+1)g/r
Cr                             and       f_l-1 Y_l-1,m  with f proportional to dg/dr + (l)g/r
Cr   See https://www.questaal.org/docs/numerics/spherical_harmonics/#gradients-of-spherical-harmonics
Cr
Cr   The raw number of product basis functions and their corresponding gradients differ because
Cr   gi*gj is symmetry in ij, while gi grad gj is not.
Cr   Moreover, the orthogonalization step need not reduce the number in an equivalent manner.
Cr
Cr   Call tree:
Cr     do  ib = 1, nbas
Cr       call rprodbasi <- product basis for a particular site
Cr       In rprodbasi
Cr         call prodphi <- products of partial waves and/or gradients
Cr         construct product basis, possibly grad PB
Cr         orthonormalize PB, possibly grad PB
Cr     enddo
Cu Updates
Cu   19 Aug 18 First cut at gradients
Cu   12 Dec 17 Adapted from old GW basnfp_v2
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Passed parameters
      integer opt,nl,nsp,nbas,nat,ndrphi,nrmx,ndrpb,nradpb
      integer :: nrphi(0:nl-1,nat),nocc(0:nl-1,ndrphi,nat),nunocc(0:nl-1,ndrphi,nat)
      integer :: lcuta(nat),nrpb(0:2*(nl-1),nat+1),ntpba(0:nat)
      real(8) :: tolbas(0:2*(nl-1)),gtota(nrmx,0:nl-1,ndrphi,nsp,nat)
      real(8) :: rprodb(nrmx,*)
C ... Local parameters
CC     interface serves no purpose, but gfortran requires it because of 'target' in rprodbasi
C      interface
C      integer function rprodbasi(mode,nl,nsp,nrmx,nr,ri,rwt,ndrphi,nrphi,
C     .  nocc,nunocc,phipb,tolbas,lcut,off,nblraw,nblr,rprodb,fiterr)
CC ... Passed parameters
C      integer nl,ndrphi,nrmx,nr,nsp,lcut,mode,off
C      integer nocc(0:nl-1,ndrphi),nunocc(0:nl-1,ndrphi),nrphi(0:nl-1),nblraw(0:2*(nl-1)),nblr(0:2*(nl-1))
C      real(8) :: ri(nrmx),rwt(nrmx),tolbas(0:2*(nl-1)),fiterr(2)
C!     real(8) :: phime(nrmx,0:nl-1,ndrphi,nsp)
C      real(8),target :: phipb(nrmx,0:nl-1,ndrphi,nsp)
C      real(8) :: rprodb(nrmx,*)
C      end function rprodbasi
C      end interface
      real(8) :: rofi(nrmx),rwgt(nrmx),xx(1),err(2)
      integer,parameter:: NULLI = -99999
      integer i,ib,iat,is,lb,nrawl,nrawlm,nredl,nredlm,nnzero,off,mode,mode1,mode2,
     .  nrpbi,ntpb,nblraw(0:2*(nl-1)),nblr(0:2*(nl-1)),ntpbwk(nat)
      procedure(integer) :: GWinputcat,fopng,a2vec,nglob,rprodbasi

C ... Setup
      do  lb = 2*(nl-1), 1, -1
        if (tolbas(lb) /= tolbas(lb-1)) exit
        nnzero = lb
      enddo
      mode1 = mod(opt/100,10)   ! 100s digit opt -> 10s digit mode
      mode2 = mod(opt/1000,10)  ! 1000s digit opt -> 100s digit mode
      if (mod(mode1,4) == 0) mode1 = 0

      call info8(20,1,0,' prbas: for'//
     .  '%?#n# grad##%-1j%?#n==2#+##%-1j%?#n==3#-## phi'//
     .  '%?#n# grad##%-1j%?#n==2#+##%-1j%?#n==3#-## phi'//
     .  ' for %i sites  nrphi=%i  lmax=%i  tol=%s,%n;3e',
     .  mod(mode1/4,2)*mod(mode1,4),(1-mod(mode1/4,2))*mod(mode1,4),
     .  nat,ndrphi,nl-1,nnzero,tolbas,8)

      iat = 0; nradpb = 0; nrpbi = 0; ndrpb = 0; ntpb = 0; ntpba(0) = 0
      call info0(30,0,0,' site raw nrad ntot  reduced nrad ntot  nrad by l%55p mean, max err')

C --- For each site, make product basis ---
      mode = 10*mode1 + 100*mode2
      off = -1 ; if (mod(opt/10,10) > 1) off = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        if (s_spec(is)%lmxa < 0) cycle
        iat = iat+1
        call rmeshprm(2+4,s_spec,nrmx,is,is,20,xx,xx,xx,xx,rofi,rwgt)
        i = rprodbasi(mode,nl,nsp,nrmx,s_spec(is)%nr,rofi,rwgt,
     .    ndrphi,nrphi(0,iat),nocc(0,1,iat),nunocc(0,1,iat),gtota(1,0,1,1,iat),
     .    tolbas,lcuta(iat),off,nblraw,nblr,rprodb,err)
        nradpb = nradpb + i

C   ... Bookkeeping and printout
        nrawl = 0; nrawlm = 0; nredl = 0; nredlm = 0; nnzero = 0
        do  lb = 0, 2*(nl-1)
          if (nblr(lb) > 0) nnzero = lb+1
          nrawl  = nrawl  + nblraw(lb)
          nrawlm = nrawlm + nblraw(lb)*(2*lb+1)
          nredl  = nredl  + nblr(lb)
          nredlm = nredlm + nblr(lb)*(2*lb+1)
          nrpb(lb,iat) = nrpbi
          nrpbi = nrpbi + nblr(lb)
          ndrpb = max(ndrpb,nblr(lb))
        enddo
        nrpb(0,iat+1) = nrpbi
        ntpbwk(iat) = nredlm
        ntpb = ntpb + nredlm

        call info8(30,0,0,' %,3i   %2,6i%,13i%,6i  %s,%ni%?#n>=100#%52p%2:1;7,7F##',
     .    ib,[nrawl,nrawlm],nredl,nredlm,nnzero,nblr,mode,err)
      enddo
      call info2(20,1,0,' System has %i (radial) %i (total) B functions',nrpb(0,nat+1),ntpb)
      if (iat /= nat) call rx('rprodbas : atom mismatch')

C ... Fill out ntpba
      do  iat = 1, nat
        ntpba(iat) = ntpbwk(iat)
        if (mod(opt/10,2) == 1) ntpba(iat) = ntpbwk(iat) + ntpba(iat-1)
      enddo

      end

      subroutine mergephivc(mode,nl,nsp,nat,ndrphiv,ndrphic,ndrphi,nrmx,
     .  nrphiv,nrphic,noccv,nunoccv,noccc,nunoccc,nocc,nunocc,nrphi,gvala,gcora,gtota)
C- Merges core and valence occ number tables into one large table.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 nocc,nunocc are not touched
Ci         :1 merge nrphi(v,c) -> nrphi; nocc(v,c) -> nocc; nunocc(v,c) -> nunocc
Ci         :10s digit
Ci         :1 merge (gval, gcore) -> gtota
Ci   nl    :(global maximum l) + 1 --- here merely a dimensioning parameter
Ci   nat   :number of sites in basis with augmentation spheres
Ci  ndrphiv:(dimensioning) : max number of radial valence partial waves of a particular l
Ci         :ndrphiv formerly named nnv in old GW code
Ci  ndrphic:(dimensioning) : max number of core levels of a particular l
Ci         :ndrphic formerly named nnc in old GW code
Ci   ndrphi:(dimensioning) : max number of radial (core states + valence partial waves)
Ci   nrmx  :number of radial mesh points (needed if partial waves are copied)
Ci   nrphiv:number of partial waves for this l and site
Ci         :nrphiv formerly named nphitv in old GW code
Ci   nrphic:number of core levels for this l and site
Ci         :nrphic formerly named nphitc in old GW code
Ci     ... If indices are copied, the following are needed
Ci   noccv :nocc for valence states (see nocc below)
Ci  nunoccv:nunocc for valence states (see nunocc below)
Ci   noccc :nocc for core states (see nocc below)
Ci  nunoccc:nunocc for core states (see nunocc below)
Ci     ... If waves are copied, the following are needed
Ci   ... If partial waves are copied (10s digit mode set)
Ci   gvala :valence partial waves
Ci   gcora :core partial waves
Co Outputs
Co   ... If indices are copied:
Co   nrphi :number of core levels + partial waves for a particular l and site
Co   nocc  :nocc(l,i) = 0 or 1 for i = 1..ndrphi
Co         :1 => include as left function when assembling all
Co         :     possible combinations for product basis
Co         :0 => do not include
Co         :nocc is made only if 1s digit mode is 1
Co   nunocc:nocc(l,i) = 0 or 1 for i = 1..ndrphi
Co         :1 => include as right function when assembling all
Co         :     possible combinations for product basis
Co         :0 => do not include
Co         :nunocc is made only if 1s digit mode is 1
Co   gtota :combined core and valence partial waves
Co         :gtota is made only if 10s digit mode set
Cu Updates
Cu   14 Dec 17  Taken from GW code reindx
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,nsp,ndrphi,ndrphiv,ndrphic,nat,nrmx
      integer nrphiv(0:nl-1,nat),nrphic(0:nl-1,nat),
     .  noccv(0:nl-1,ndrphiv,nat),noccc(0:nl-1,ndrphic,nat),
     .  nunoccv(0:nl-1,ndrphiv,nat),nunoccc(0:nl-1,ndrphic,nat),
     .  nocc(0:nl-1,ndrphi,nat),nunocc(0:nl-1,ndrphi,nat),nrphi(0:nl-1,nat)
      real(8) :: gvala(nrmx,0:nl-1,ndrphiv,nsp,nat)
      real(8) :: gcora(nrmx,0:nl-1,ndrphic,nsp,nat)
      real(8) :: gtota(nrmx,0:nl-1,ndrphi,nsp,nat)

C ... Local parameters
      integer ib,l,n,ncore,nval,ir,isp

      do  ib = 1, nat
        do  l = 0, nl-1
          ncore = nrphic(l,ib)
          nval  = nrphiv(l,ib)

C     ... Copy indices
          if (mod(mode,10) == 1) then

            nrphi(l,ib)= ncore + nval
            if (ncore+nval > ndrphi) call rx('mergephivc: ncore+nval > nrphi')

            do  n = 1, ncore
              nocc(l,n,ib)   = noccc(l,n,ib)
              nunocc(l,n,ib) = nunoccc(l,n,ib)
            enddo

            do  n = 1,nval
              nocc(l,ncore+n,ib)   = noccv(l,n,ib)
              nunocc(l,ncore+n,ib) = nunoccv(l,n,ib)
            enddo
          endif

C     ... Copy wave functions
          if (mod(mode/10,10) == 1) then
            do  isp = 1, nsp
            do  n = 1, ncore
              forall (ir=1:nrmx) gtota(ir,l,n,isp,ib) =  gcora(ir,l,n,isp,ib)
            enddo
            do  n = 1, nval
              forall (ir=1:nrmx) gtota(ir,l,n+ncore,isp,ib) =  gvala(ir,l,n,isp,ib)
            enddo
            enddo
          endif

        enddo                   ! loop over l
      enddo                     ! loop over sites

      end

      integer function maxnphi(nrphiv,nrphic,nl,nat)
C- Finds the maximum number of (core states + valence partial waves), for dimensioning
C  Taken from old GW, maxnn
      implicit none
C ... Passed parameters
      integer nl,nat
      integer nrphiv(nl*nat),nrphic(nl*nat)
C ... Local parameters
      integer i

      maxnphi = -1
      do i = 1, nl*nat
        maxnphi = max(maxnphi,nrphiv(i)+nrphic(i))
      enddo
      end

      integer function nodnum(f,n)
C- Count the number of nodes in a partial wave
C Taken from old GW, which was taken from basn.f
      implicit none
      integer n
      real(8):: f(n)
      integer i
      nodnum = 0
      do i = 2, n-1
        if (f(i)*f(i+1) < 0) nodnum=nodnum+1
      enddo
      end
