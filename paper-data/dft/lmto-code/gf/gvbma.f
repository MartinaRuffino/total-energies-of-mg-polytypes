      subroutine gvbma(s_site,s_pot,mode,nl,nsp,nbas,lhdim,indxsh,pp,bma)
C- Create a vector beta-alpha
C ----------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  dlmcl ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dlmwt
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1 make bma for bare representation (beta=0)
Ci          2 make bma for gamma representation (beta=gamma)
Ci          3 make bma for gamma representation averaged over spins
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   lhdim:size of lower + downfolding block
Ci   idpc  :class indexo
Ci         :For a normal site, site ib belongs to class ipc(ib) (mksym.f)
Ci         :For a CPA site, idpc is the first CPA class associated with ib
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci          (makidx.f)
Ci   pp    :potential parameters (atomsr.f)
Co Outputs
Co   bma   :gamma-alpha for transformation to new representation
Cr Remarks
Cr   Mapping of composite RL index to this vector is same as that of
Cr   eigenvectors permutated according to downfolding rules in indxsh
Cu Updates
Cu   04 Jan 13 Standard class index ipc replaced with ipdc for CPA
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nl,nsp,nbas,lhdim
      integer indxsh(*)
      double precision pp(6,nl,nsp,*),bma(lhdim,nsp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(nbas)
      type(str_pot)::   s_pot
C ... Local parameters
      integer ibas,ic,l,m,lmr,isp,ncomp
      double precision alpha,beta,cpagam(nl,nsp)

      do  isp = 1, nsp
      lmr = 0
      do  ibas = 1, nbas
C       ic = ipdc(ibas)
        ic = s_site(ibas)%dlmcl         ! Parent class

C       Returns gamma ppar, or if a CPA site, CPA-averaged gamma
        ncomp = s_site(ibas)%ncomp  ! Number of CPA components
        call cpagamma(ncomp,ic,s_pot%dlmwt,nl,nsp,pp,cpagam)

        do  l = 0, nl-1
          if (indxsh(lmr+1) > lhdim) then
            lmr = lmr + 2*l+1
            cycle
          endif
          do  m = -l, l
            lmr = lmr+1
            alpha = pp(6,l+1,isp,ic)
            if (mode == 0) beta = alpha
            if (mode == 1) beta = 0
            if (mode == 2) beta = cpagam(l+1,isp)
            if (mode == 3) beta = (cpagam(l+1,1)+cpagam(l+1,nsp))/2
            bma(indxsh(lmr),isp) = beta - alpha
          enddo
        enddo
      enddo
      enddo

C     call prmx('beta-alpha',bma,lhdim,lhdim,nsp)
      end
