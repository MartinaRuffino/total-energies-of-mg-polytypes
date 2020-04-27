      subroutine suham2(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_strn)
C- Further hamiltonian setup after potential generated
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lfp
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  ng alat tolft
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb pz name orbp ngcut
Co     Stored:     orbp ngcut
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb sugcut
Ci   s_site (not used)
Ci   s_ham (not used)
Ci   s_pot (not used)
Ci   s_strn (not used)
Ci Inputs
Co Outputs
Cl Local variables
Cr Remarks
Cr   This routine completes energy-independent hamiltonian setup,
Cr   which requires information not when suham is called.
Cr   The basis, e.g. local orbitals, depend on the potential
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   16 Aug 04 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_strn) :: s_strn(*)
C ... Local parameters
      integer nspec,lfp,nglob
      integer ng,i
      double precision alat,tolgv,dum

C --- Setup ---
      lfp = mod(s_ctrl%lfp,2)
      nspec = nglob('nspec')

C --- FP setup ---
      if (lfp /= 0) then
        call uspecb(0,0,s_spec,1,nspec,dum,dum,dum,i)
C       If any local orbitals ...
        if (i >= 100) then
          ng = s_lat%ng
          alat = s_lat%alat
          tolgv = s_lat%tolft
          call sugcut(2,nspec,s_spec,alat,ng,s_lat%gv,tolgv)
        endif
      endif

      end
