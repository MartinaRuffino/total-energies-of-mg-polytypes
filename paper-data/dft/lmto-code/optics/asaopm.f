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
Ci   lmxa  :augmentation l-cutoff, by class
Ci   ldimx :hamiltonian lower block dimension, (2x dim, noncoll case)
Ci   nev   :actual number of eigenvectors generated
Ci   qp    :k-point
Ci   ikp   :k-point index, used to specify where to poke optmt,velmt
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   eband :energy bands; alias eb (sec*.f)
Ci   nbmax :maximum number of bands
Ci   nfilm :(no occ states) dimensions optmt and velmt
Ci   nempm :(no unocc states) dimensions optmt and velmt
Cio Inputs/Outputs
Co   accob :Contains decomposition of norm (aka dnpp; see makwts.f)
Co         :accob(:,ib,:,1,1:2,:) projection of phi onto site ib (real,imag)
Co         :accob(:,ib,:,2,1:2,:) projection of phid onto site ib (real,imag)
Co         :accob is transformed to spherical harmonics in cmplx*16 format
Co Outputs
Co   optmt(m,i,j,ikp) |<i|grad|j>| connecting occ i with unocc j, polarization m
Co   velmt(m,i,ikp) i=j part of optmt
Cr Remarks
Cr   Optics package adapted from Sergey Rashkeev with Walter Lambrecht,
Cr   which was adapted from an earlier version by V. Antropov.
Cu Updates
Cu   18 Dec 14 Adjust for new modes loptic=8,9,-8,-9
Cu   21 Aug 14 Revisited with simplified treatment of range in occ, unocc states
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
      integer, allocatable :: nrfn(:),nlmi(:)
      real(8), allocatable :: ccp(:)
C     real(8), allocatable :: gradm(:,:,:,:,:,:,:)
      complex(8), allocatable :: zgradm(:,:,:,:,:,:,:)
      complex(8), allocatable :: optme(:,:,:)
C ... Local parameters
      logical lsph,leps
      integer i,iabc(3,6),icheck,lham,loptic,lshg,lshgd,morder,nchi2,
     .  nemlo,nemup,nfilo,nfiup,nsgrp,ocrng(4),unrng(4)
      double precision esciss
      procedure(logical) :: bittst
      procedure(integer) :: isw,bitand
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))

      call tcn('asaopm')
      ocrng = s_optic%ocrng
      unrng = s_optic%unrng
      nchi2 = s_optic%nchi2
      iabc = s_optic%axes
      nbas = s_ctrl%nbas
      lham = s_ctrl%lham
      loptic = s_ctrl%loptc
      esciss = s_optic%esciss
      nsgrp = s_lat%nsgrp

      lsph  = bittst(lham,256)
      leps  = loptic > 0 .and. mod(loptic,10) > 0
      lshg  = bitand(loptic/10,2)
      lshgd = bitand(loptic/10,4)
      if (.not. leps .and. lshg == 0 .and. lshgd == 0) return

C     ASA always has 2 radial functions for each site
C     And use fixed lmxa = nl-1 for each site
      allocate(nrfn(nbas),nlmi(nbas))
      do  i  = 1, nbas
        nrfn(i) = 2
        nlmi(i) = nl**2
      enddo

C --- Transform accob to spherical harmonics ---
      allocate(ccp(nspc*nl*nl*nbas*4))
C     call pvasaopm(lsph,nbas,nl,nspc,nev,ccp,accob)
C     if (ikp == 1) print *, '!! call to pvasaopm suppresses SH'
      call pvasaopm(.true.,nbas,nl,nspc,nev,ccp,accob)
      deallocate(ccp)

      if (nspc==2) call rx('asaopm: grrme not available for nc case')

C --- Matrix element phi grad phi etc from radial matrix elements ---
      i = 3*nl**4*nsp*nclass
      allocate(zgradm(nl**2,nl**2,nsp*nspc,nclass,3,2,2))

C     Debugging: compatibility with old code
C     opt=301 => ordering same as ogradme but grady -> -grady
C     opt=302 => m order reversed from ogradme but grady -> -grady
C     In both cases, if1,if2 order is reversed relative to ogradme
C      call gradme(301,2,2,isp,nsp,nspc,1,nclass,nl-1,lmxa,nl**2,
C     .  s_pot%grrme,zgradm,zgradm)

C --- Electric dipole matrix elements ---
      if (leps) then

C       For spherical gradients:
        call lmorder(0,morder,iabc,iabc)  ! iabc just dummy
        morder = 2-morder
        if (.not. lsph) morder = 0
        i =  30 + morder
C       i = i+ 300  ! Return gradient of x+iy, x-iy
        call gradme(i,2,2,isp,nsp,nspc,1,nclass,nl-1,lmxa,nl**2,
     .    s_pot%grrme,zgradm,zgradm)

        allocate(optme(nfilm,nempm,3)); call dpzero(optme,2*nfilm*nempm*3)
C       Use opt=1 for compatibility with old code from Vladimir
        call optdme(1,2,nrfn,nl,nbas,isp,nsp,nspc,nev,s_ctrl%ipc,
     .    zgradm,nfilo,nfiup,nemlo,nemup,1d0,accob,optme)

        call symdme(1,s_lat%symgr,nsgrp,
     .    nfilo,nfiup,nemlo,nemup,nfilo,nfiup,nemlo,nemup,
     .    1d0,optme,optmt(1,1,1,isp,ikp),optmt)
        deallocate(optme)

      endif


C --- Second harmonic generation ---
      if (lshg /= 0) then
        if (.not. lsph) call rx('SHARM needed for SHG ... please set OPTICS_SHARM=t')
        allocate(optme(nfilo:nemup,nfilo:nemup,3)); call dpzero(optme,2*size(optme))

        call lmorder(0,morder,iabc,iabc)  ! iabc just dummy
        morder = 2-morder
        i =  30 + morder
        call gradme(i,2,2,isp,nsp,nspc,1,nclass,nl-1,lmxa,nl**2,
     .    s_pot%grrme,zgradm,zgradm)

C       Use opt=1 for compatibility with old code from Vladimir
        call optdme(1,2,nrfn,nl,nbas,isp,nsp,nspc,nev,s_ctrl%ipc,
     .    zgradm,nfilo,nemup,nfilo,nemup,-1d0,accob,optme)
        call optshs(ikp,nkp,nchi2,iabc,
     .    optme,s_lat%symgr,nsgrp,esciss,nfilo,nfiup,nemlo,nemup,nbmax,eband)
        deallocate(optme)
       endif

       if (lshgd /= 0) then
         call optshd(2,2,nl,nbas,isp,nsp,ikp,nkp,ldimx,s_ctrl%ipc,nchi2,
     .     iabc,nclass,zgradm,icheck,
     .     s_lat%symgr,nsgrp,qp,esciss,nfilo,nfiup,nemlo,nemup,nbmax,
     .     eband,accob)
      endif

      deallocate(zgradm,nrfn,nlmi)

      call tcx('asaopm')

      end

      subroutine pvasaopm(lsph,nbas,nl,nspc,nev,ccp,accob)
C- Transforms accob to complex*16 format, optionally to spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci   lsph  :F transform to spherical harmonics
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nev   :actual number of eigenvectors generated
Ci   ccp   :work array dimensioned at least nl**2*nbas*nspc*2
Cio Inputs/Outputs
Cio  accob :Contains decomposition of norm (aka dnpp; see makwts.f)
Cio        :Indices are permuted and possible rotation to spherical harmonics
Cio        :On input:
Cio        :  accob(:,ibas,:,1,1:2,ib) projection of phi from state ib onto site ibas (real,imag)
Cio        :  accob(:,ibas,:,2,1:2,ib) projection of phid, state ib,   onto site ibas (real,imag)
Cio        :accob is overwritten in cmplx*16 format, true Y_lm:
Cio        :accob(lm,ibas,isp,:,1:2,:) -> accob(1:2,lm,ibas,isp,:,:)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lsph
      integer nl,nev,nbas,nspc
      double precision accob(0:nl*nl-1,nbas,nspc,2,2,nev)
      double complex ccp(0:nl*nl-1,nbas,nspc,2)
C ... Local parameters
      integer i,ib,ibas,l1,ll,lmm,lpp,m1,jsp
      double complex cccos,ccsin
C ... External calls
      external dcopy

C --- For each each eigenvector, do ---
      do  ib = 1, nev

C --- ccp <- accob transformed to regular spherical harmonics Y_lm ---
      do  ibas = 1, nbas
      do  i = 1, 2 ! phi and phidot
      do  jsp= 1, nspc
        if (.not. lsph) then
          do  l1 = 0, nl-1
          do  m1 = 0, l1
            if (m1 == 0) then
              ll = l1*l1 + l1 + m1
              ccp(ll,ibas,jsp,i) =
     .          dcmplx(accob(ll,ibas,jsp,i,1,ib),accob(ll,ibas,jsp,i,2,ib))
            else
              lpp = l1*l1 + l1 + m1
              lmm = l1*l1 + l1 - m1
              ccsin = dcmplx(accob(lmm,ibas,jsp,i,1,ib),accob(lmm,ibas,jsp,i,2,ib))
              cccos = dcmplx(accob(lpp,ibas,jsp,i,1,ib),accob(lpp,ibas,jsp,i,2,ib))
C             This is standard definition of spherical harmonics, m ordered l:-l
              ccp(lmm,ibas,jsp,i) = (cccos - dcmplx(0d0,1d0)*ccsin)/dsqrt(2d0)
              ccp(lpp,ibas,jsp,i) = (cccos + dcmplx(0d0,1d0)*ccsin)/dsqrt(2d0)*(-1)**m1
C             This is definition of spherical harmonics used until 2018
              ccp(lmm,ibas,jsp,i) = (cccos + dcmplx(0d0,1d0)*ccsin)/dsqrt(2d0)
              ccp(lpp,ibas,jsp,i) = (cccos - dcmplx(0d0,1d0)*ccsin)/dsqrt(2d0)*(-1)**m1
            endif
          enddo
          enddo
        else
          do  ll = 0, nl**2-1
            ccp(ll,ibas,jsp,i) = dcmplx(accob(ll,ibas,jsp,i,1,ib),accob(ll,ibas,jsp,i,2,ib))
          enddo
        endif
      enddo                     ! noncollinear spins
      enddo                     ! kind of partial wave
      enddo                     ! sites

C --- Overwrite accob with result in spherical harmonics ---
      call dcopy(nspc*nl*nl*nbas*4,ccp,1,accob(0,1,1,1,1,ib),1)

      enddo
      end
