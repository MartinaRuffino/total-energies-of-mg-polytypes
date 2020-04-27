      subroutine asaopm(s_ctrl,s_pot,s_lat,s_optic,nl,isp,nsp,nspc,
     .  nclass,nbas,lmxa,ldimx,nev,qp,ikp,nkp,eband,nbmax,accob,nfilm,
     .  nempm,optmt,velmt)
C- Makes electric dipole and second harmonic optical matrix elements
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas lham loptc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ipc
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:grrme
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:symgr
Cio    Passed to:  *
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ocrng unrng nchi2 axes esciss
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nclass:number of inequivalent classes
Ci   nbas  :size of basis
Ci   ldimx :hamiltonian lower block dimension, (2x dim, noncoll case)
Ci   nev   :actual number of eigenvectors generated
Ci   qp    :k-point
Ci   ikp   :k-point index, used to specify where to poke optmt,velmt
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   eband :energy bands; alias eb (sec*.f)
Ci   nbmax :maximum number of bands
Cio Inputs/Outputs
Co   accob :Contains decomposition of norm (aka dnpp; see makwts.f)
Co         :accob(:,ib,:,1,1:2,:) projection of phi onto site ib (real,imag)
Co         :accob(:,ib,:,2,1:2,:) projection of phid onto site ib (real,imag)
Co         :accob is transformed to spherical harmonics in cmplx*16 format
Co Outputs
Co   optmt(*,ikp) ASA <i|grad|j> connecting occ i with unocc j
Co   velmt(*,ikp) ASA <i|grad|i>
Cr Remarks
Cr   Optics package adapted from Sergey Rashkeev with Walter Lambrecht,
Cr   which was adapted from an earlier version by V. Antropov.
Cu Updates
Cu   30 Dec 13 Redesign to use with FP package
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   22 Feb 04 Redefined flags in loptc switch
Cu   20 Nov 02 (jek) extended to metals
Cu   20 Dec 00 (wrl) extended to noncollinear case
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nsp,nclass,ldimx,nbmax,ikp,nkp,nfilm,nempm,isp,
     .  nbas,nspc,nev,lmxa(nclass)
      double precision qp(3,*),eband(nbmax)
      double precision accob(0:nl*nl-1,nbas,nspc,2,2,nev)
      double precision optmt(3,nfilm,nempm,nsp/nspc,nkp),
     .                 velmt(3,nfilm,nsp/nspc,nkp)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_optic):: s_optic
C ... Dynamically allocated local arrays
      real(8), allocatable :: ccp(:)
      real(8), allocatable :: gradm(:,:,:,:,:,:,:)
C ... Local parameters
      integer bitand,i,iabc(3,6),icheck,leps,lham,loptic,lshg,
     .  lshgd,nchi2,nemlo,nemup,nfilo,nfiup,nsgrp,ocrng(2),unrng(2)
      logical lsph,bittst
      double precision esciss
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))

      call rx('OPTINQ: OPTICS package not installed')

      end
