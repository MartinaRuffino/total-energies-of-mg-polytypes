      subroutine smcorm(nbas,s_site,s_spec,s_lat,ng,gv,cgh1,cgh2,
     .  lfoc1,lfoc2)
C- For foca, add together density of smoothed part of core
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lfoca rfoca qc z ctail etail stc lmxb p pz rmt rg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  corprm
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cy
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Co Outputs
Co   cgh1  :Portion of smoothed core that is treated directly
Co   cgh2  :Portion of smoothed core that is treated perturbatively
Co   lfoc1 :returned nonzero if any site lfoca is direct (1)
Co   lfoc2 :returned nonzero if any site lfoca is perturbative
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   02 Jul 05  skip sites for which cofh=0
Cu    1 May 00  Adapted from nfp smc_mesh.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ng,nbas,lfoc1,lfoc2
      double precision gv(ng,3)
      double complex cgh1(ng),cgh2(ng)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Local parameters
      integer k0,nlmx,kmax,ib,is,lfoc,i
      double precision tau(3),v(3),alat,vol,qcorg,qcorh,qsc,cofg,cofh,
     .  ceh,rfoc,z
C     integer stdo,iprint,lgunit
      parameter (k0=3, nlmx=25)
      double complex gkl(0:k0,nlmx)

C     stdo = lgunit(1)
      alat = s_lat%alat
      vol = s_lat%vol
      kmax = 0

C --- Accumulate FT of smooth-Hankel foca heads ---
      call dpzero(cgh1,  2*ng)
      call dpzero(cgh2,  2*ng)
      lfoc1 = 0
      lfoc2 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
C       qc = qcorg+qcorh
C        if (iprint() >= 50) write(stdo,351) qc,lfoc,qcorg,qcorh
C  351   format(' qc=',f12.6,'   lfoc',i2,'   qcorg,qcorh',2f12.6)

        if (cofh /= 0) then
        if (lfoc == 1) then
          lfoc1 = 1
          do  i = 1, ng
            v(1) = gv(i,1)
            v(2) = gv(i,2)
            v(3) = gv(i,3)
            call hklft(v,rfoc,ceh,tau,alat,kmax,1,k0,s_lat%cy,gkl)
            cgh1(i) = cgh1(i) + cofh*gkl(0,1)/vol
          enddo
        else if (lfoc == 2) then
          lfoc2 = 1
          do  i = 1, ng
            v(1) = gv(i,1)
            v(2) = gv(i,2)
            v(3) = gv(i,3)
            call hklft(v,rfoc,ceh,tau,alat,kmax,1,k0,s_lat%cy,gkl)
            cgh2(i) = cgh2(i) + cofh*gkl(0,1)/vol
          enddo
        endif
        endif
      enddo

C      write(stdo,846) cgh1(1)*vol,cgh2(1)*vol
C  846 format(' cgh(1) times vol',4f12.6)

      end




