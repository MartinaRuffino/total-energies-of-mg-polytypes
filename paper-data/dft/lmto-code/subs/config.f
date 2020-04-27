      subroutine config(pnu,lmax,z,konfig,lmaxc)
C- Returns principal quantum numbers from pnu, estimating those unknown
C ----------------------------------------------------------------
Ci Inputs
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   lmax  :maximum l for which pnu is supplied
Ci   z     :nuclear charge
Co Outputs
Co   konfig:estimated principal quant. no. of valence state for l>lmax
Co         :Core orbitals are specified by:
Co         :  1, 2, ..., konf(0)-1 for s
Co         :  2, 3, ..., konf(1)-1 for p
Co         :  3, 4, ..., konf(2)-1 for d, and so on.
Co   lmaxc :largest l for which there is a core state
Cr Remarks
Cu Updates
Cu   08 Feb 01 generate config for l=0..8
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer lmax,konfig(0:8),lmaxc
      double precision pnu(0:*),z
C Local parameters
      integer l
      integer LI,NA,K,RB,CS,FR,CU,AG,AU,HF
      parameter (LI=3,
     .           NA=11, K=19, RB=37, CS=55, FR=87,
     .           CU=29, AG=47, AU=79, HF=72)

C --- Calculate lmaxc ---
      lmaxc = lmax
      if (z >= LI) lmaxc = max(lmax,0)
      if (z >= NA) lmaxc = max(lmax,1)
      if (z >= CU) lmaxc = max(lmax,2)
      if (z >= HF) lmaxc = max(lmax,3)

C --- Estimate konfig ---
      do  l = 0, 8
        konfig(l) = l+1
      enddo
      if (z >= LI) konfig(0) = 2
      if (z >= NA) konfig(0) = 3
      if (z >= K) konfig(0) = 4
      if (z >= RB) konfig(0) = 5
      if (z >= CS) konfig(0) = 6
      if (z >= FR) konfig(0) = 7
      konfig(1) = max(konfig(0),2)

      if (z >= CU) konfig(2) = 4
      if (z >= AG) konfig(2) = 5
      if (z >= AU) konfig(2) = 6

      if (z >= HF) konfig(3) = 5

C --- Override konfig with given pnu ---
      do  l = 0, lmax
        konfig(l) = pnu(l)
      enddo
      end
