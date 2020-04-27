      subroutine defspc(s_spec)
C- Sets defaults for species data
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rmt rg rsma rfoca rsmfa
Co     Stored:     rg rsma rfoca rsmfa rsmv
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Cu Updates
Cu   06 Sep 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)

      double precision rg,rsma,rfoca,rsmfa,rmt
      integer nspec,is,nglob

      nspec = nglob('nspec')
      do  is = 1, nspec

        rmt = s_spec(is)%rmt
        rmt = s_spec(is)%rmt
        rg = s_spec(is)%rg
        rsma = s_spec(is)%rsma
        rfoca = s_spec(is)%rfoca
        rsmfa = s_spec(is)%rsmfa

        rmt = s_spec(is)%rmt
        rg = s_spec(is)%rg
        rsma = s_spec(is)%rsma
        rfoca = s_spec(is)%rfoca
        rsmfa = s_spec(is)%rsmfa

        if (rg == 0) rg    = -1
        if (rg < 0) rg    = -rg*0.25d0*rmt
        if (rsma == 0) rsma  = -1
        if (rsma < 0) rsma  = -rsma*0.4d0*rmt
        if (rfoca == 0) rfoca = -1
        if (rfoca < 0) rfoca = -rfoca*0.4d0*rmt
        if (rsmfa == 0) rsmfa = -1
        if (rsmfa < 0) rsmfa = -rsmfa*0.5d0*rmt

        s_spec(is)%rg = rg
        s_spec(is)%rsma = rsma
        s_spec(is)%rfoca = rfoca
        s_spec(is)%rsmfa = rsmfa
        s_spec(is)%rsmv = rmt*.5d0

C        print 333, is,rg,rsma,rfoca,rsmfa
C  333   format(i3,4f12.6)
      enddo
C     stop
      end
