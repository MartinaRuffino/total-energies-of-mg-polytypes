      subroutine mtodim(mode,s_site,s_spec,nl,nbas,nat,nrphiv,nrphic,nlnmaug)
C- Return dimensioning of partial waves in LMTO basis
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
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
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci  mode   :1s digit
Ci         :0 exclude all states core states
Ci         :1 include core states
Ci         :2 include valence states
Ci         :3 combination of 1+2
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  nbas   :size of basis
Ci  nat    :number of sites in basis with augmentation spheres
Ci  nrphiv :number of (valence partial waves) for a particular l
Ci         :Set nrphiv(0,1) = -1 if to exclude valence states
Ci         :Formerly named nindxv in old GW code
Ci  nrphic :number of (core partial waves) for a particular l
Ci         :Set nrphic(0,1) = -1 if to exclude core states
Ci         :Formerly named nindxc in old GW code
Co Outputs
Co  nlnmaug:total number of partial waves in LMTO basis by site, including all l,m,radial waves
Co         :These include all orbitals phi for which <phi phi | prbas> are to be made (see prodbasa)
Co         :naugl(iat,1) = number of core + valence
Co         :naugl(iat,2) = number belonging to core only
Co         :Formerly nlnm (aka mnl) contained naugl(:,1);  nlnmv,nlnmc divided nlnm into valence, core
Cu Updates
Cu   22 Aug 18 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Passed parameters
      integer mode,nl,nbas,nat
      integer :: nrphiv(0:nl-1,nat),nrphic(0:nl-1,nat),nlnmaug(nat,2)
C ... Local parameters
      integer,parameter:: NULLI = -99999
      integer ib,iat,is,indexi,iwk(1)
      logical addcore,addval

      addcore = mod(mod(mode,10),2) == 1
      addval  = mod(mod(mode,10)/2,2) == 1

      iat = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        if (s_spec(is)%lmxa < 0) cycle
        iat = iat+1
        indexi = 0
        if (addcore) call mtoindex(-nl,nrphic(0,iat),-1,1,iwk,iwk,indexi)
        nlnmaug(iat,2) = indexi
        if (addval) call mtoindex(-nl,nrphiv(0,iat),nrphic(0,iat),1,iwk,iwk,indexi)
        nlnmaug(iat,1) = indexi
      enddo
      end

      subroutine mtoindex(nl,nrphi,nphi0,ndrphi,indxlm,mtoindx,offlm)
C- Indices for LMTO basis functions or product basis functions for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci  nl     :abs(nl) = 1 + maximum l for LMTO basis.
Ci         :Set nl<0 if not to accumulate indxlm, mtoindx (only offlm is modified)
Ci         :If mtoindex is called in the product basis context, use 2*(nl-1) for this arg
Ci  nrphi  :number of partial waves for this l and site
Ci         :This routine does nothing if nrphi(0)<0
Ci         :Formerly called nindx
Ci         :If mtoindex is called in the product basis context, use nrpb for this arg
Ci  nphi0  :Offsets to be added to nrphi when assembling indxlm(2,:)
Ci         :Not used if nphi0(0)<0.  Use when appending valence states to cores.
Ci  ndrphi :dimensioning parameter, must be no smaller than maxval(nrphi)
Ci         :Not used if nl<0
Cio Inputs/Outputs
Cio  offlm :offset to next partial wave.
Cio        :offlm is incremented by the number of waves with unique l,n,m
Cio        :where l and m are angular moment, and n is index to radial part
Cio        :Set to 0 to start at beginning; see Remarks
Co Outputs
Co         :Partial waves are ordered in a list by : m(inner), n(middle), l(outer)
Co         :for a particular (l,m), L = l*l + l + m + 1
Co         :mtoindx contains index to element in list given (n,l,m)
Co         :indxlm yields (n,l,m) for a particular element in list.
Co  indxlm :(returns only if nl>0) three values for each partial wave
Co         :indxlm(1,mtoindx) = index to radial function n
Co         :indxlm(2,mtoindx) = l quantum number
Co         :indxlm(3,mtoindx) = m quantum number
Co         :Former names in old GW code are the arrays  il, im, in
Co  mtoindx:(returns only if nl>0).  For given l,n,m, points to all
Co         :mtoindx(n,lm) = position in list of partial wave with q.n. (m,n,l)
Co         :mtoindx formerly named ilnm in old GW code
Cr Remarks
Cr   This routine constructs indices to a list of augmented partial waves
Cr   in the order that they will be used in assembling the full product basis.
Cr
Cr   To start the list from the beginning, set offlm=0 before calling mtoindex.
Cr   For core states only, set nrphi to nrphic
Cr   For valence states only, set nrphi to nrphiv
Cr   To combine in one group, pass nrphi as nrphic+nrphiv
Cr   To join valence following core (the convention of the old GW code),
Cr   make two calls to mtoindex; see mtodim above for example.
Cr
Cr   For now, follows Kotani convention: ordering is m(inner), n(middle), l(outer)
Cu Updates
Cu   22 Aug 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,ndrphi,offlm
      integer nrphi(0:abs(nl)-1),nphi0(0:abs(nl)-1)
      integer indxlm(3,ndrphi*nl*nl),mtoindx(ndrphi,nl*nl)
C ... Local parameters
      integer l,l2,n,m,lm
      logical lnphi0

      if (nrphi(0) < 0) return

      lnphi0 = nphi0(0) >= 0
      do  l = 0, abs(nl)-1
        l2 = l*l
        do  n = 1, nrphi(l)
          do  m = 1, 2*l+1
            offlm = offlm + 1
            if (nl<0) cycle

            lm = l2 + m
            indxlm(1,offlm) = l
            indxlm(2,offlm) = n
            if (lnphi0) indxlm(2,offlm) = indxlm(2,offlm) + nphi0(l)
            indxlm(3,offlm) = m - l - 1
            mtoindx(n,lm)= offlm
            if (n>ndrphi) call rx('mtoindx: increase dimensioning')
          enddo
        enddo
      enddo

      end

C      subroutine prbindex(nl,nrpb,ndrpb,indxpb,pbindx,offpb)
CC- Indices for product basis functions for one site
CC ----------------------------------------------------------------------
CCi Inputs
CCi  nl     :abs(nl) = 1 + maximum l for product basis.
CCi         :Set nl<0 if not to accumulate indxpb, pbindx (only offpb is modified)
CCi         :If mtoindex is called in the product basis context, use 2*(nl-1) for this arg
CCi  nrpb   :number of partial waves for particular l
CCi  ndrpb  :dimensioning parameter, must be no smaller than maxval(nrpb)
CCi         :Not used if nl<0
CCio Inputs/Outputs
CCio  offpb :offset to next partial wave.
CCio        :offpb is incremented by the number of waves with unique l,n,m
CCio        :where l and m are angular moment, and n is index to radial part
CCio        :Set to 0 to start at beginning; see Remarks
CCo Outputs
CCo         :Partial waves are ordered in a list by : m(inner), n(middle), l(outer)
CCo         :for a particular (l,m), L = l*l + l + m + 1
CCo         :pbindx contains index to element in list given (n,l,m)
CCo         :indxpb yields (n,l,m) for a particular element in list.
CCo  indxpb :(returns only if nl>0) three values for each partial wave
CCo         :indxpb(1,pbindx) = index to radial function n
CCo         :indxpb(2,pbindx) = l quantum number
CCo         :indxpb(3,pbindx) = m quantum number
CCo         :Former names in old GW code are the arrays  il, im, in
CCo  pbindx:(returns only if nl>0).  For given l,n,m, points to all
CCo         :pbindx(n,lm) = position in list of partial wave with q.n. (m,n,l)
CCo         :pbindx formerly named ilnm in old GW code
CCr Remarks
CCr   This routine constructs a list of augmented partial waves
CCr   in the order that they will be used in assembling the full product basis.
CCr
CCr   To start the list from the beginning, set offpb=0 before calling mtoindex.
CCr
CCr   For now, follows Kotani convention: ordering is m(inner), n(middle), l(outer)
CCu Updates
CCu   22 Aug 18 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nl,ndrpb,offpb
C      integer nrpb(0:abs(nl)-1),indxpb(3,ndrpb*nl*nl),pbindx(ndrpb,nl*nl)
CC ... Local parameters
C
C      print *, 'replacement for pbindex.  Not needed yet'
C
CC     call mtoindex(nl,nrpb,-1,ndrpb,indxpb,pbindx,offpb)
C      end
