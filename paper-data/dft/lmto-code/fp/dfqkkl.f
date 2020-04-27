      subroutine dfqkkl(nbas,lekkl,s_site,s_spec,numq)
C- Allocates arrays to accumulate local output density
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  qkkl eqkkl qhkl eqhkl qhhl eqhhl
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
ci   lekkl :Make eq{hh,hk,kk}l as well as q{hh,hk,kk}l
Ci   numq  :number of Fermi levels for which to accumulate c.d.
Co Outputs
Cl Local variables
Cl   nkapi :number of envelope function types per l q.n. for spec is2
Cl   nkaph :number of orbital types for a given L quantum no. in basis
Cr Remarks
Cu Updates
Cu   22 Nov 12 Replace qkkl with structures
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle lmxa=-1 -> no allocation
Cu   15 Jun 05 Allocation for noncollinear case
Cu   25 Aug 01 Extended to local orbitals
Cu   15 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp df_qkkl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,lekkl,numq
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,is,kmax,lmxa,lmxh,nelt1,nglob,nlma,nlmh,nsp,
     .  nspc,nelt3,nelt2,nkaph,n2,nhk
      double precision xx

      nsp = nglob('nsp')
      nspc = nglob('nspc')
      nkaph = nglob('nkaph')

C --- Loop over sites, allocating qkkl for each site ---
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        lmxh = s_spec(is)%lmxb
        kmax = s_spec(is)%kmxt
        if (lmxa == -1) goto 10

        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
        n2 = numq*nsp*nspc
        nhk = numq*nsp*(2*nspc-1)

C   ... Case Pkl*Pkl
        nelt1 = (kmax+1)*(kmax+1)*nlma*nlma
        call ptr_site(s_site,1,'qkkl',ib,nelt1,n2,xx)
        if (lekkl == 1) call ptr_site(s_site,1,'eqkkl',ib,nelt1,n2,xx)

C   ... Case Pkl*Hsm
        nelt2 = (kmax+1)*nkaph*nlma*nlmh
        call ptr_site(s_site,1,'qhkl',ib,nelt2,nhk,xx)
        if (lekkl == 1) call ptr_site(s_site,1,'eqhkl',ib,nelt2,nhk,xx)

C   ... Case Hsm*Hsm
        nelt3 = nkaph*nkaph*nlmh*nlmh
        call ptr_site(s_site,1,'qhhl',ib,nelt3,n2,xx)
        if (lekkl == 1) call ptr_site(s_site,1,'eqhhl',ib,nelt3,n2,xx)

c|        write(6,836) nelt1,nelt3,nelt2
c|  836   format('   nelt=',3i6)
   10   continue
      enddo

      end
