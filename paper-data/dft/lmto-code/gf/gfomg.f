      subroutine gfomg(s_ctrl,s_site,s_pot,s_ham,izp)
C- Calculate CPA coherent interactor Omega for all DLM sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbasp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:omgn
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cp
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Ci Inputs
Ci   izp   :complex energy number for access to Omega
Co Outputs
Co   Omega updated for all DLM sites, stored in s_site->omgn
Cu Updates
Cu   19 Aug 13 (Belashchenko) Added relativistic case (2 x 2)
Cu   19 Dec 12 conversion to F90 pointers complete
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
Cu   10 Jan 12 (Belashchenko) All redone
Cu   08 Dec 08 (P. Larson) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed parameters:
      integer izp
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_site)::  s_site(*)
C Local parameters:
      integer nbas,nspc,ib
      integer ldham(16),lidim,mxorb,nglob
      equivalence (lidim,ldham(2))

      mxorb = nglob('mxorb')
      nspc = nglob('nspc')

      nbas = s_ctrl%nbas
      ldham = s_ham%ldham

C ... Omega = (P - g^-1)^-1
      do  ib = 1, nbas
        if (s_site(ib)%ncomp < 2) cycle
        call omedlm(s_site(ib)%gii,nspc,s_site(ib)%norb,lidim,s_pot%cp,
     .    s_ham%iprmb(1+mxorb*(ib-1):s_site(ib)%norb+mxorb*(ib-1)),s_site(ib)%omgn(:,izp))
      enddo

      end

      subroutine omedlm(gii,nspc,norb,lidim,cp,iprm,omega)
C- Make coherent potential interactor for one site, for ASA GF
C ----------------------------------------------------------------------
Ci Inputs
Ci   gii   :Green function
Ci   norb  :number of orbitals belonging to current site
Ci   lidim :dimension of gii
Ci   cp    :Coherent potential
Ci   iprm  :permutations ordering orbitals for given site (makidx.f)
Co Outputs
Co   omega :coherent interactor Omega for current site
Cr Remarks
Cr   Omega = cp - gloc^-1, where gloc is on-site part of gii:
Cr           i.e. sum_q gii(q) and R=R'.
Cu Updates
Cu   09 Jul 14 (Belashchenko) further simplified, now calls pokeg0
Cu   18 Jun 14 (Belashchenko) simplified, omedlmr eliminated
Cu   19 Dec 12 conversion to F90 pointers complete
Cu   10 Jan 12 (Belashchenko) Rewritten
Cu   08 Dec 08 (P. Larson) First created
C ----------------------------------------------------------------------
      implicit none
      integer norb,nspc,lidim,iprm(norb)
      complex(8) cp(lidim,nspc,lidim,2)
      complex(8),dimension(norb,nspc,norb,2) :: gii,omega
      complex(8),dimension(norb,nspc,norb,2) :: wk
      real(8), allocatable :: wrk(:)
      integer isp,ndim,ierr

C ... wk = on-site g in complex format
      wk = gii

      ndim = norb*nspc
      allocate(wrk(66*ndim))
      do  isp = 1, 3-nspc
        call zqinv('N',wk(1,1,1,isp),ndim,-66*ndim,ndim,wrk,ndim,ierr)
        if (ierr /= 0) call rx1('omedlm: inversion,ierr=%i',ierr)
      enddo
      deallocate(wrk)

C ... Omega = cp - g^-1
      call pokeg0(1,norb,iprm,lidim,lidim,lidim,nspc,2,omega,cp)
      omega = omega - wk

      end


