      subroutine gtpcor(s_spec,is,kcore,lcore,qcore)
C- Unpacks parameters related to partial core occpation
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  coreh coreq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   is    :species for which to unpack kcore,lcore,qcore
Co Outputs
Co   kcore  :p.q.n for occupation
Co   lcore  :l quantum for occupation
Co   qcore  :core charge and magnetic moment
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer is,kcore,lcore
      double precision qcore(2)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      character*8 ch

      kcore = 0
      lcore = -1
      qcore(1) = 0
      qcore(2) = 0
      ch = s_spec(is)%coreh
      if (ch == ' ') return
      qcore = s_spec(is)%coreq
      read (ch,'(i1)') kcore
      if (ch(2:2) == 's' .or. ch(2:2) == 'S') lcore = 0
      if (ch(2:2) == 'p' .or. ch(2:2) == 'P') lcore = 1
      if (ch(2:2) == 'd' .or. ch(2:2) == 'D') lcore = 2
      if (ch(2:2) == 'f' .or. ch(2:2) == 'F') lcore = 3
      end
