      subroutine surho(nbas,s_site,s_spec,lmet,ldos,lrout,
     .  lekkl,numq,k1,k2,k3,smrout,ndos,dos,sumev,sumqv)
C- Initialize output density coeffs, dos, and eigenvalue sum
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qkkl qhkl qhhl eqkkl eqhkl eqhhl
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   lmet  :nonzero => assume metal and init dos
Ci   ldos  :nonzero => init dos
Ci   lrout :nonzero => init qkkl,smrout
C    lekkl :nonzero => init eqkkl
Ci   numq  :number of trial Fermi levels
Ci   k1..3 :dimensions smrout
Ci   smrout:output density
Ci   ndos  :number of energy mesh points
Ci   dos   :density of states
Ci   sumev :band sum of eigenvalues
Ci   sumqv :band sum of charges
Co Outputs
Co   See Remarks
Cr Remarks
Cr   The following steps are taken:
Cr     1.  Local arrays qkkl are zeroed
Cr     2.  smrout is zeroed
Cr     3.  DOS and eigenvalue sums are zeroed.
Cu Updates
Cu   22 Nov 12 Replace qkkl with structures
Cu   10 Nov 11 Begin migration to f90 structures
Cu   27 Jul 08 (T. Kotani) added eqkkl
Cu   02 Jan 06 sumqv doubled for two spin channels
Cu   01 Jul 05 handle lmxa=-1 -> no local arrays
Cu   15 Jun 05 qkkl zeroed for noncollinear case
Cu   23 Jan 01 Added lrout switch
Cu   19 Jun 00 spin polarized
Cu    5 May 00 Adapted from nfp init_density
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k1,k2,k3,ldos,lmet,lrout,lekkl,nbas,ndos,numq
C     integer oqkkl(3,nbas),oeqkkl(3,nbas)
      double precision sumev(2,3),sumqv(3,2),
     .  dos(ndos,2)
      double complex smrout(k1,k2,k3,numq)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,i,is,kmax,lmxa,lmxh
      integer nglob,nsp,nspc,nkaph,n2

      nsp = nglob('nsp')
      nspc = nglob('nspc')
      nkaph = nglob('nkaph')

C ... Reset coeffs of local density
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        lmxh = s_spec(is)%lmxb
        kmax = s_spec(is)%kmxt
        if (lmxa == -1) goto 10

C       nlma = (lmxa+1)**2
C       nlmh = (lmxh+1)**2
C       nelt1 = (kmax+1)*(kmax+1)*nlma*nlma
C       nelt2 = (kmax+1)*nkaph*nlma*nlmh
C       nelt3 = nkaph*nkaph*nlmh*nlmh
C       nhk = numq*nsp*(2*nspc-1)
        n2 = numq*nsp*nspc

        if (lrout /= 0) then
          call dpzero(s_site(ib)%qkkl, size(s_site(ib)%qkkl))
          call dpzero(s_site(ib)%qhkl, size(s_site(ib)%qhkl))
          call dpzero(s_site(ib)%qhhl, size(s_site(ib)%qhhl))
        endif
        if (lrout /= 0 .and. lekkl /= 0) then
          call dpzero(s_site(ib)%eqkkl, size(s_site(ib)%eqkkl))
          call dpzero(s_site(ib)%eqhkl, size(s_site(ib)%eqhkl))
          call dpzero(s_site(ib)%eqhhl, size(s_site(ib)%eqhhl))
        endif

   10   continue
      enddo

C ... Reset array for smooth density
      if (lrout /= 0) call dpzero(smrout,2*k1*k2*k3*n2)

C ... DOS and eigenvalue sum
      do  i = 1, numq
        sumev(1,i) = 0d0
        sumev(2,i) = 0d0
        sumqv(i,1) = 0d0
        sumqv(i,2) = 0d0
      enddo
      if (ldos /= 0 .or. lmet > 0) call dpzero(dos,2*ndos*nsp)

      end
