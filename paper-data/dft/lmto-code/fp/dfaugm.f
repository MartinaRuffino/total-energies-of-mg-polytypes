      subroutine dfaugm(nbas,lcplxp,lso,lx,s_site,s_spec)
C- Allocate augmentation matrices sigma,tau,pi for all atoms
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  sigkk taukk sigkkx taukkx pikk pikkx sighk tauhk
Co                 sighkx tauhkx pihk pihkx sighh tauhh sighhx tauhhx
Co                 pihh pihhx
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   lcplxp:if nonzero allocate space for complex ppi
Ci   lso   :flag for SO coupling
Ci   lx    :if 0, allocate sig,tau,pi
Ci         :else, allocate sigx,taux,pix
Co Outputs
Cr Remarks
Cr   Pointers are specified as osig(itype,ibas) where
Cr     {sig,tau,ppi}kk: case < Pkl | {sig,tau,ppi} | Pkl >
Cr     {sig,tau,ppi}hk: case < Pkl | {sig,tau,ppi} | Hsm >
Cr     {sig,tau,ppi}hh: case < Hsm | {sig,tau,ppi} | Hsm >
Cr   sig and tau are diagonal in m, ppi is full matrix
Cr   Thus integral (P~_kL P~_k'L' - P_kL P_k'L') is diagonal in LL',
Cr       sig(nf1,nf2,0..lmax) with lmax the l-cutoff
Cr   For sig(Pkl,Pkl), nf1=nf2==1+kmax; lmax=lmxa
Cr   For sig(Hsm,Pkl), nf1=nkaph and nf2=1+kmax; lmax=lmxh
Cr   For sig(Hsm,Hsm), nf1=nf2=nkaph; lmax = lmxh
Cl Local variables
Cl   nkapi :number of envelope function types per l q.n. for spec is2
Cl   nkaph :number of orbital types for a given L quantum no. in basis
Cu Updates
Cu   05 Jul 13 Split off SO part of hamiltonian --- separate from ppi
Cu   30 Aug 12 Modifications for SO=3
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle lmxa=-1 -> no allocation
Cu   29 Jun 05 Adapted to store SO in ppi separately from ppi
Cu    1 Sep 04 Adapted to handle complex ppi.  so folded into ppi
Cu   29 Jun 04 (A. Chantis) memory allocation for LzSz matrix elements.
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   25 Aug 01 Extended to local orbitals
Cu   11 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp df_augm.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,lcplxp,lso,lx !,osig(3,nbas),otau(3,nbas),oppi(3,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,is,kmax,lmxa,lmxh,nelt1,nelt2,nglob,nlma,nlmh,
     .  nsp,nkaph,nspc3
      double precision xx

      nsp = nglob('nsp')
C     nspc = nglob('nspc')     ! number of coupled spins
      nspc3 = nglob('nspc3')   ! number of coupled spins, incl perturbatively
      nkaph = nglob('nkaph')
      call rxx((lso == 1.or.lso == 3.or.lso == 4).and.nspc3 /= 2,
     .  'dfaugm: inconsistent parameters')

C --- Loop over sites, allocating sig,tau,pi for each site ---
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        lmxh = s_spec(is)%lmxb
        kmax = s_spec(is)%kmxt
        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
        if (lmxa == -1) cycle

C   ... Case Pkl*Pkl
        nelt1 = (kmax+1)*(kmax+1)*(lmxa+1)*nsp
        nelt2 = (kmax+1)*(kmax+1)*nlma*nlma*nsp*nspc3
        if (lx == 0) then
          call ptr_site(s_site,1,'sigkk',ib,nelt1,1,xx)
          call ptr_site(s_site,1,'taukk',ib,nelt1,1,xx)
        else
          call ptr_site(s_site,1,'sigkkx',ib,nelt1,1,xx)
          call ptr_site(s_site,1,'taukkx',ib,nelt1,1,xx)
        endif
        if (lx == 0 .and. lcplxp == 0) then
          call ptr_site(s_site,1,'pikk',ib,1,nelt2,xx)
        elseif (lx == 0 .and. lcplxp /= 0) then
          call ptr_site(s_site,1,'pikk',ib,2,nelt2,xx)
          if (lso /= 0) then
            call ptr_site(s_site,1,'sokk',ib,2,nelt2,xx)
          else ! Dummy array
            call ptr_site(s_site,1,'sokk',ib,2,1,xx)
          endif
        elseif (lx == 1 .and. lcplxp == 0) then
          call ptr_site(s_site,1,'pikkx',ib,1,nelt2,xx)
        else
          call ptr_site(s_site,1,'pikkx',ib,2,nelt2,xx)
        endif

C   ... Case Hsm*Pkl
        if (lmxh > lmxa) call rx('dfaugm: lmxh > lmxa unexpected')
        nelt1 = nkaph*(kmax+1)*(lmxh+1)*nsp
        nelt2 = nkaph*(kmax+1)*nlmh*nlma*nsp*nspc3
        if (lx == 0) then
          call ptr_site(s_site,1,'sighk',ib,nelt1,1,xx)
          call ptr_site(s_site,1,'tauhk',ib,nelt1,1,xx)
        else
          call ptr_site(s_site,1,'sighkx',ib,nelt1,1,xx)
          call ptr_site(s_site,1,'tauhkx',ib,nelt1,1,xx)
        endif
        if (lx == 0 .and. lcplxp == 0) then
          call ptr_site(s_site,1,'pihk',ib,1,nelt2,xx)
        elseif (lx == 0 .and. lcplxp /= 0) then
          call ptr_site(s_site,1,'pihk',ib,2,nelt2,xx)
          if (lso /= 0) then
            call ptr_site(s_site,1,'sohk',ib,2,nelt2,xx)
          else ! Dummy array
            call ptr_site(s_site,1,'sohk',ib,2,1,xx)
          endif
        elseif (lx == 1 .and. lcplxp == 0) then
          call ptr_site(s_site,1,'pihkx',ib,1,nelt2,xx)
        else
          call ptr_site(s_site,1,'pihkx',ib,2,nelt2,xx)
        endif

C   ... Case Hsm*Hsm
        nelt1 = nkaph*nkaph*(lmxh+1)*nsp
        nelt2 = nkaph*nkaph*nlmh*nlmh*nsp*nspc3
        if (lx == 0) then
          call ptr_site(s_site,1,'sighh',ib,nelt1,1,xx)
          call ptr_site(s_site,1,'tauhh',ib,nelt1,1,xx)
        else
          call ptr_site(s_site,1,'sighhx',ib,nelt1,1,xx)
          call ptr_site(s_site,1,'tauhhx',ib,nelt1,1,xx)
        endif
        if (lx == 0 .and. lcplxp == 0) then
          call ptr_site(s_site,1,'pihh',ib,1,nelt2,xx)
        elseif (lx == 0 .and. lcplxp /= 0) then
          call ptr_site(s_site,1,'pihh',ib,2,nelt2,xx)
          if (lso /= 0) then
          call ptr_site(s_site,1,'sohh',ib,2,nelt2,xx)
          else
          call ptr_site(s_site,1,'sohh',ib,2,1,xx)
          endif
        elseif (lx == 1 .and. lcplxp == 0) then
          call ptr_site(s_site,1,'pihhx',ib,1,nelt2,xx)
        else
          call ptr_site(s_site,1,'pihhx',ib,2,nelt2,xx)
        endif

      enddo

      end
