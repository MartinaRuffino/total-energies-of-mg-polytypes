      subroutine gfgii(s_ctrl,s_site,s_ham,gii)
C- Extract the site-diagonal blocks of gii and record in s_site
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
Cio    Elts passed:gii
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Ci Inputs
Ci   gii   :q-summed g (g_RR' for R,R' inequivalent sites in the basis)
Ci         :Generated in gfibz
Co Outputs
Co   gii stored for all sites in s_site->gii
Cu Updates
Cu   18 Jun 14 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed parameters:
      complex(8) gii(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_site)::  s_site(*)
C Local parameters:
      integer nspc,ib
      integer ldham(16),lidim,mxorb,nglob
      equivalence (lidim,ldham(2))

      mxorb = nglob('mxorb')
      nspc = nglob('nspc')

      ldham = s_ham%ldham

C ... Extract the site-diagonal block of gii and record in s_site
      do  ib = 1, s_ctrl%nbasp
        call pokeg0(1,s_site(ib)%norb,s_ham%iprmb(1+mxorb*(ib-1)),lidim,lidim,lidim,
     .    nspc,2,s_site(ib)%gii,gii)
      enddo

      end

