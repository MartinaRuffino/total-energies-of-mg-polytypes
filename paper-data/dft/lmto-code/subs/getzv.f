      subroutine getzv(nclass,nrc,z,qc,s_bz,s_ctrl,zval)
C- Calculate and save number of valence electrons
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  zval
Co     Stored:     zval
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  zbak
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lasa
Cio    Passed to:  *
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   nrc   :nrc(ic) = number of sites belonging to class ic
Co Outputs
Co   zval  :valence charge
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nclass
      integer nrc(nclass)
      double precision qc(nclass),z(nclass)
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
C ... Local parameters
      double precision zval,zbak(2)
      integer ic
      integer,parameter:: NULLI = -99999
      procedure(integer) :: iprint,lgunit,isw

      zval = s_bz%zval
      if (zval == 0 .or. zval == NULLI) then
        zval = 0
        do  ic = 1, nclass
          zval = zval+nrc(ic)*(z(ic)-qc(ic))
        enddo
      endif
      s_bz%zval = zval

C --- Printout ---
      zbak = s_ctrl%zbak
      if (iprint() >= 30) then
        call awrit4('%N GETZV:  %d valence electrons'//
     .    '%?#n#  zbak=%d  qbak=%d##',' ',80,lgunit(1),zval,
     .    isw(IAND(s_ctrl%lasa,64) /= 0),zbak,zbak(2))
      endif
      if (iprint() >= 20) then
        call awrit4(' zval %d'//'%?#n#  zbak %;2d  qbak %;2d##',
     .    ' ',80,lgunit(2),zval,isw(IAND(s_ctrl%lasa,64) /= 0),
     .    zbak,zbak(2))
      endif

      end
