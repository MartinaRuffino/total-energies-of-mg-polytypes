      subroutine intdrv(s_ctrl,s_optic,s_bz,eband,nbmax,nsp,nspc,efermi,
     .  idtet,vol,nfilm,nempm)
C- Driver for SHG Integrals (nonlinear optics)
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  loptc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  dw window ocrng unrng esciss nchi2 axes
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp ntet
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs:
Ci   eband :energy bands for irreducible part of BZ
Ci   nbmax :leading dimension of eband
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   efermi:Fermi level
Ci   idtet :(0,i) no. of tetrahedra of the i'th kind
Ci         :(1-4,i) identifies the i'th tetrahedron in terms of
Ci         :four irreducible k-points
Ci   vol   :volume
Ci   nfilm,nempm: dimensions optmt
Co Outputs:
Co   Im(eps) or joint density of states written to file 'optdf'
Cr Remarks
Cr   Adapted from bzints to make joint density of states or Im(eps)
Cr   All energy differences between states below ef and states
Cr   above ef+emin are summed, and integrated over the BZ
Cr   Treatment near the critical points (ef and ef+emin) handled crudely
Cr   Optics package adapted from Sergey Rashkeev with Walter Lambrecht,
Cr   which was adapted from an earlier version by V. Antropov.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   02 Mar 01 Added scissors operator
Cu   20 Dec 00 (wrl) extended to noncollinear case
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbmax,nsp,nspc,idtet(0:4,*),nfilm,nempm
      double precision efermi,vol
      double precision eband(nbmax,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_optic):: s_optic
      type(str_bz)::    s_bz
C ... Local parameters
      logical bittst
      integer nfilo,nfiup,nemlo,nemup,nkp,ntet,npts
      integer ifio,fopn,loptic
      integer ocrng(2),unrng(2)
      integer nkabc(3),n1,n2,n3
      double precision optrng(4),emin,emax
      equivalence (emin,optrng(1)),(emax,optrng(2))
      double precision esciss,dw(3)
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      integer nchi2, iabc(3,6)

C ... Setup
      loptic = s_ctrl%loptc
      if (bittst(loptic/10,4)) then
        ifio = fopn('SHGD')
        rewind ifio
      endif

      dw = s_optic%dw
      optrng = s_optic%window
      ocrng = s_optic%ocrng(1:2)
      unrng = s_optic%unrng(1:2)
      npts = nint((emax-emin)/dw(1)) + 1
      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      ntet = s_bz%ntet
      esciss = s_optic%esciss
      nchi2 = s_optic%nchi2
      iabc = s_optic%axes

C ... Call Integration Subroutines
      if (bittst(loptic/10,2)) then
        call shsint(n1,n2,n3,nkp,nsp,ntet,idtet,vol,esciss,nfilo,
     .  nfiup,nemlo,nemup,nbmax,eband,nchi2,iabc)
      elseif (bittst(loptic/10,4)) then
        call shdint(ifio,n1,n2,n3,nkp,nsp,emin,emax,efermi,npts,
     .  ntet,idtet,vol,esciss,nfilo,nfiup,nemlo,nemup,nbmax,eband)
      endif

C ... Close Files
      if (bittst(loptic/10,4)) then
        call fclose(ifio)
      endif

      end

