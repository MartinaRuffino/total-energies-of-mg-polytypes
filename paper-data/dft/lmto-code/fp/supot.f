C#define FP
      subroutine supot(mode,s_ctrl,s_lat,s_pot)
C- Initialization for potential
C ----------------------------------------------------------------------
Cio Structures
C    ... For ASA
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: nbas nbasp nspin
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat vol awald nkd nkq
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:pos dlv qlv
Cio    Passed to: *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read: *
Co     Stored:    *
Co     Allocated: mad
Cio    Elts Passed:mad
Cio    Passed to: *
C    ... for FP
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nbasp nspin
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lcd
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat vol awald nkd nkq plat nabc gmax ng qlat kv kv2
Co     Stored:     nabc gmax ng
Co     Allocated:  ips0 bgv gv kv igv kv2 igv2
Cio    Elts passed:pos dlv qlv nsgrp symgr ag gv ips0 bgv igv igv2
Cio    Passed to:  sugvec0 sugvec
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  mad smrho smpot
Cio    Elts passed:mad rhat rnew
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1 make Madelung matrix for monopoles (ASA)
Co Outputs
Cr Remarks
Cr   Permanent arrays are created to handle the potential, depending
Cr   on the mode
Cu Updates
Cu   08 Aug 17 lmf phases for special AFM symmetrization
Cu   05 May 12 Repackaged G vector generation into sugvec
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
!      type(str_ham)::   s_ham
C ... Local parameters
C#ifdef FP
      integer, allocatable :: ipca(:),istab(:,:)
      real(8), allocatable :: g(:,:),ag(:,:)
      integer ng,ngabc(3),n1,n2,n3,k1,k2,k3,ngmx
      double precision gmax,plat(3,3)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C#endif
      integer nbas,nbasp,nbaspp,nsp,nspc,i
C     integer npadl,npadr
      double precision xx
      procedure (integer) :: iprint,idalloc,allocvb,nglob

      call info0(50,1,0,' supot : potential setup ...')

      nbas = s_ctrl%nbas
C     npadl = s_ctrl%npadl
C     npadr = s_ctrl%npadr
      nbasp = s_ctrl%nbasp
      nsp = s_ctrl%nspin
      nspc  = nglob('nspc')

      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) == 1) then
        call ptr_pot(s_pot,1,'mad',nbaspp**2,0,xx)
        call madmat(nbaspp,s_lat%pos,s_lat%awald,s_lat%alat,abs(s_lat%vol),
     .    s_lat%dlv,s_lat%nkd,s_lat%qlv,s_lat%nkq,s_pot%mad)
C        call smadmat(nbaspp,s_lat%pos,s_lat%awald,s_lat%alat,abs(s_lat%vol),
C     .    s_lat%dlv,s_lat%nkd,s_lat%qlv,s_lat%nkq,s_lat%cg,s_lat%indxcg,s_lat%jcg,s_pot%mad)
C       call madmat2D(nbas,npadl,npadr,s_lat,s_pot%mad)
      endif

C#ifdef FP
C --- Setup for FT charge density, potential representation ---
      if (IAND(s_ctrl%lcd,4) /= 0) then
        plat = s_lat%plat
        ngabc = s_lat%nabc
        gmax = s_lat%gmax

C   ... Generate energy cutoff gmax or n1..n3
        call pshpr(iprint()-0)
        call mshsiz(s_lat%alat,plat,0,gmax,ngabc,ngmx)
        call poppr
        call fftz30(n1,n2,n3,k1,k2,k3)
        s_lat%nabc = ngabc
        s_lat%gmax = gmax
        if (ngmx <= 0) call rx('supot: missing G vectors.  Check HAM_GMAX.')

C   ... Make list of lattice vectors within cutoff
        call sugvec0(s_lat)
        ng = s_lat%ng
        if (ng /= ngmx) call rx('supot: bug in gvlist')

        i = idalloc('irho+v',allocvb()+2,k1*k2*k3,nsp)
        call ptr_pot(s_pot,8+1,'smrho',k1*k2*k3,nsp*nspc,xx)
        call ptr_pot(s_pot,8+1,'smpot',k1*k2*k3,nsp*nspc,xx)
C       allocate(s_pot%rhat(nbaspp),s_pot%rnew(nbaspp))
        call ptr_pot(s_pot,1,'rhat',nbaspp,1,xx)
        call ptr_pot(s_pot,1,'rnew',nbaspp,1,xx)

C   ... Setup for symmetrization and possible AFM symmetrization
        i = 1; if (iabs(s_lat%nsafm) /= 0) i = 2
        call ptr_lat(s_lat,1,'ips0',ng*i,0,0,0)
        call ptr_lat(s_lat,1,'bgv',ng*i,0,0,0)
        call sgvsym(s_lat%nsgrp,s_lat%symgr,s_lat%ag,ng,s_lat%gv,s_lat%ips0,s_lat%bgv)
C   ... Setup for AFM symmetrization
        if (iabs(s_lat%nsafm) /= 0) then
          allocate(ipca(nbas),istab(nbas,2),g(9,2),ag(3,2))
          call info0(20,0,0,'%9f... setup for AFM symmetrization')
          call suafmsym(s_lat,nbas,ipca,istab,g,ag)
          call sgvsym(-2,g,ag,ng,s_lat%gv,s_lat%ips0(1+ng),s_lat%bgv(1+ng))
        endif
      endif
C#endif
      end
