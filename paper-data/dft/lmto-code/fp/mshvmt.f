      subroutine mshvmt(nbas,s_site,s_spec,s_lat,ng,gv,cv,vval)
C- Makes potential at MT surfaces given potential on a uniform mesh
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rmt lmxl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   cv    :smooth potential in reciprocal space, in glist form
Ci         :generate it from:
Ci         :call gvgetf(ng,1,kv,k1,k2,k3,smpot,cv)
Co Outputs
Co   vval  :coffs to YL expansion of es potential at MT boundary
Co         :for each site, computed from the mesh density.
Cr Remarks
Cr   A PW exp(i.q.r) has a one-center expansion at radius r
Cr      sum_L C_L Y_L(r) where C_L = 4 pi i^l j_l(|rq|) Y_L(q)
Cr   Routine symvvl symmetrizes the vval generated here.
Cb Bugs
Cb   Possible to make ves for sites with lmxl=-1, which tells
Cb   value of ves at point.  However, vval doesn't have the
Cb   space allocated.  So skip for now
Cu Updates
Cu   16 Jun 16 Modifications for new ability to add sm external potential
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxl=-1
Cu   22 Aug 01 Newly created.
C ----------------------------------------------------------------------
      use structures
      use mpi

      implicit none
C ... Passed parameters
      integer nbas,ng
      double precision gv(ng,3),vval(*)
      double complex cv(ng)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      double precision,allocatable:: phil(:,:),yl(:,:)
      double precision,allocatable:: gv2(:,:),agv(:),cgp(:),sgp(:)
C ... Local parameters
      integer i,ib,is,lmxx,nlmx,iv0,lmxl,nlm,ngabc(3),
     .  n1,n2,n3,m,ilm,l,ipr
      double precision alat,pi,tpiba,tau(3),rmt,fac,plat(3,3)
C     double precision q(3),qlat(3,3)
      double complex vvali,fprli
C     procedure(logical) :: cmdopt
      procedure(integer) :: cmdoptswx
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      parameter (lmxx=8, nlmx=(lmxx+1)**2)

      integer :: iproc, nproc, comm, err, nlms(0:nbas)
      integer, allocatable :: pmap(:), smap(:), nmap(:)

      call tcn('mshvmt')

      allocate(phil(ng,0:lmxx),yl(ng,nlmx))
      allocate(gv2(ng,3),agv(ng),cgp(ng),sgp(ng))
      call getpr(ipr)
      pi = 4d0*datan(1d0)
      alat = s_lat%alat
      plat = s_lat%plat
      ngabc = s_lat%nabc
      tpiba = 2*pi/alat

C --- YL(G)*G**l, agv=|G| for each g ---
      call dpcopy(gv,gv2,1,3*ng,tpiba)
      call ropyln(ng,gv2(1,1),gv2(1,2),gv2(1,3),lmxx,ng,yl,agv)
      do  i = 1, ng
        agv(i) = sqrt(agv(i))
      enddo

      nlms(0) = 0
      do ib = 1, nbas
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
!         if (lmxl == -1) lmxl = 0
        nlm = (lmxl+1)**2
        nlms(ib) = nlms(ib-1) + nlm
      end do

      comm = mpi_comm_world
      call mpi_comm_size(comm, nproc, err)
      call mpi_comm_rank(comm, iproc, err)

      allocate(pmap(0:nproc))
      call vbdist(nbas, nlms(1:nbas), nproc, pmap)

C --- For each ib in nbas, do ---
      do ib = pmap(iproc)+1, pmap(iproc+1)
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        rmt = s_spec(is)%rmt
        lmxl = s_spec(is)%lmxl
        if (lmxl == -1) cycle
        nlm = (lmxl+1)**2
        if (nlm > nlmx) call rxi('mshvmt: increase nlmx to',nlm)
C       Add a negligibly small amount to rmt to handle case rmt=0
        rmt = rmt+1d-32

C   --- j_l(|rmt*q|)/rmt**l for each G and l=0..lmax ---
        i = 1
C       Without i=11, will not evolve correctly in the large agv*r limit
        if (cmdoptswx('--newvesmt','','') > 0) i = 11
        call ropbes(agv,rmt**2,lmxl,cgp,sgp,phil,ng,i)

C   ... Phases exp(-i.G.tau), fast version
C       call suphs0(plat,ng,gv,gv2)
C       call dinv33(plat,1,qlat,fac)
C       call dpzero(q,3)
C       call suphas(q,tau,ng,gv2,n1,n2,n3,qlat,cgp,sgp)
C   ... Phases calculated straightforwardly.  Fast enough not to matter.
        call dscal(3,alat,tau,1)
        do  i = 1, ng
          fac = -(tau(1)*gv2(i,1)+tau(2)*gv2(i,2)+tau(3)*gv2(i,3))
          cgp(i) = dcos(fac)
          sgp(i) = dsin(fac)
        enddo

        iv0 = nlms(ib-1)
C   --- Sum_G 4*pi*(i*rmt)**l j_l(|rmt*G|)/(rmt*G)**l YL(G) G**l ---
C       call dpzero(vval(iv0+1),nlm)
        ilm = 0
        fprli = 4*pi
        do  l  = 0, lmxl
          do  m = -l, l
            ilm = ilm+1
            vvali = 0
            do  i = 2, ng
              vvali = vvali +(phil(i,l)*yl(i,ilm))*(cv(i)*dcmplx(cgp(i),-sgp(i)))
            enddo
            vval(ilm+iv0) = fprli*vvali
          enddo
          fprli = fprli*(0d0,1d0)*rmt
        enddo

C   ... Printout
C        if (ipr > 1) then
C          do  ilm = 1, nlm
C            if (ilm == 1) then
C              write(stdo,650) ib,ilm,vval(ilm+iv0)
C            elseif (dabs(vval(ilm+iv0)) > 1d-6) then
C              write(stdo,651)    ilm,vval(ilm+iv0)
C            endif
C  650              format(i4,i6,2f12.6)
C  651                     format(4x,i6,f12.6)
C          enddo
C        endif

      enddo

      deallocate(phil,yl)
      deallocate(gv2,agv,cgp,sgp)

      if (nproc > 1) then
        allocate(nmap(0:nproc))
        allocate(smap(0:nproc-1))

        nmap = nlms(pmap)
        smap = nmap(1:nproc) - nmap(0:nproc-1)

        call mpi_allgatherv(mpi_in_place, smap(iproc), mpi_real8, vval, smap, nmap, mpi_real8, comm, err)

        deallocate(smap)
        deallocate(nmap)
      end if

      deallocate(pmap)

      call tcx('mshvmt')
      end

      subroutine symvvl(nbas,s_site,s_spec,s_lat,vval,vrmt)
C- Symmetrizes the potential at the MT boundary.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   s_site :struct for site-specific data; see structures.h
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: *
Ci   s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxl
Ci     Stored:    *
Ci     Passed to: *
Ci   s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: plat nsgrp symgr ag
Ci     Stored:    *
Ci     Passed to: *
Co Outputs
Cio  vval  :On input,  unsymmetrized potential
Cio        :On output, elements of potential for sites in the same
Cio        :class are symmetrized.
Co Outputs
Co   vrmt  :spherical average of potential (i.e. Y0*vval(l=0)) returned
Co         :for each site.
Cr Remarks
Cr   This routine symmetrizes any vector of the same structure as vval.
Cb Bugs
Cb   Possible to make ves for sites with lmxl=-1, which tells
Cb   value of ves at point.  However, vval doesn't have the
Cb   space allocated.  So skip for now
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxl=-1
Cu   23 Aug 01 Newly created.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas
      double precision vval(*),vrmt(nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      real(8), allocatable :: qwk(:)
      real(8), allocatable :: sym(:)
C ... Local parameters
      integer ic,ib,ilm,mxint,nclass,ipa(nbas),nrclas,stdo,iv0
      integer ipc(nbas),ips(nbas),lmxl(nbas)
      double precision pos(3,nbas),posc(3,nbas),plat(3,3),pi,y0,xx
      integer nlml,lgunit,ipr,jpr,ngrp,nn,iclbas

      stdo = lgunit(1)
      call getpr(ipr)
      plat = s_lat%plat
      ngrp = s_lat%nsgrp
      call sitepack(s_site,1,nbas,'class',1,ipc,xx)
      call sitepack(s_site,1,nbas,'spec',1,ips,xx)
      call sitepack(s_site,1,nbas,'pos',3,xx,pos)
      nclass = mxint(nbas,ipc)
      do  ib = 1, nbas
        lmxl(ib) = s_spec(ips(ib))%lmxl
      enddo

      do  ic = 1, nclass

C   ... Make nrclas,ipa,posc
        call psymr0(lmxl,ic,nbas,ipc,pos,posc,ipa,nrclas)
        if (nrclas > 0) then
        ib = iclbas(ic,ipc,nbas)
C       if (ib == 0) cycle  ! shouldn't be needed
        if (lmxl(ib) > -1) then
        nlml = (lmxl(ib)+1)**2
        if (ipr >= 50) write(stdo,800) ic,nrclas,nlml
  800   format(' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)

        allocate(qwk(nlml),sym(nlml*nlml*nrclas))
        call symqmp(nrclas,nlml,nlml,plat,posc,ngrp,s_lat%symgr,
     .    s_lat%ag,qwk,ipa,sym,vval,nn)
        deallocate(qwk,sym)
       endif
       endif
      enddo

C ... Extract vrmt = l=0 term for each site, and printout
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      if (ipr >= 45) write(stdo,221)
  221 format(/' site class  ilm      vval',6x,'ves(rmax)')
      iv0 = 0
      do  ib = 1, nbas
        if (lmxl(ib) == -1) goto 10

        nlml = (lmxl(ib)+1)**2
        vrmt(ib) = vval(1+iv0)*y0

C   ... Printout
        ic = ipc(ib)
        jpr = 0
        if (ipr > 60) jpr = 2
        if (ib == iclbas(ic,ipc,nbas)) then
          if (ipr >= 45) jpr = 1
          if (ipr >= 50) jpr = 2
        endif
        if (jpr > 0) then
          do  ilm = 1, nlml
            if (ilm == 1) then
              write(stdo,650) ib,ic,ilm,vval(ilm+iv0),vrmt(ib)
            elseif (dabs(vval(ilm+iv0)) > 1d-6  .and. jpr > 1) then
              write(stdo,651)    ilm,vval(ilm+iv0)
            endif
  650       format(i4,2i6,2f12.6)
  651       format(10x,i6,f12.6)
          enddo
        endif

        iv0 = iv0 + nlml
   10  continue
      enddo

      end
