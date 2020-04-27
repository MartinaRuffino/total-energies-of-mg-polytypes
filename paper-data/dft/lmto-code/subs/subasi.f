      subroutine subasi(s_ctrl,s_spec,s_ham)
C- Read some parameters defining hamiltonian indpendent of potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   s_ctrl (not used)
Ci   s_spec (not used)
Ci   s_ham (not used)
Co Outputs
Co   Parameters defining basis are set in sspec
Co   Global parameters nkaph and mxorb are set
Cr Remarks
Cr   This routine generates energy-independent hamiltonian setup.
Cr  *It generates and packs a table of hamiltonian offsets offH,
Cr   orbital permutation indices oindxo.
Cr
Cr  *For the ASA 2nd generation LMTO:
Cr   Extract order of potential function from gfopts
Cr   Transform pp's to alpha representation
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
C ... Local parameters
C      double precision dglob,xx
C      integer fopna,ifi,nkaph,nkmax,nlmax,
C     .  nspec,lfp

C     lfp = mod(s_ctrl%lfp,2)
C     nbasp = s_ctrl%nbasp
C     nspec = s_ctrl%nspec

C --- Get the maximum L-cutoff ---
C      Moved to rdccat.f
C      lmxax = -1
C      do  is = 1, nspec
C        lmxax = max(lmxax,igetss('spec lmxa',is,sspec))
C      enddo
C      nlmax = (lmxax+1)**2
C      xx = dglob('nkaph',1d0,1)
C      xx = dglob('mxorb',dble(nlmax),1)

      end
