      subroutine ugcomp(nbas,s_site,s_spec,s_lat,qmom,gpot0,hpot0,ugg,f)
C- Part of the smooth estatic energy from compensating G's alone.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl rg lfoca rfoca qc z ctail etail stc lmxb p pz
Ci                 rmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  corprm
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat awald tol vol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg indxcg jcg cy qlv dlv
Cio    Passed to:  ggugbl gfigbl fklbl gklbl hhugbl hhigbl phhigb hklbl
Cio                hsmbl hgugbl
Ci Inputs
Ci   nbas  :size of basis
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Cio Inputs/Outputs
Cio  Let n0  = smooth potential without compensating gaussians
Cio      n0~ = smooth potential with compensating gaussians
Cio    phi0  = ves[n0]
Cio    phi0~ = ves[n0~]
Cio    g_RL  = gaussian in RL channel
Cio    h_R   = l=0 sm hankel in RL channel, for core density
Cio  Then:
Cio  gpot0 :On input, integrals g_RL * phi0
Cio        :On output, integrals g_RL * phi0~
Cio  hpot0 :On input, integrals h_R * phi0
Cio        :On output, integrals h_R * phi0~
Co Outputs
Co   ugg   :electrostatic energy integral [n0~-n0]*[phi0~-phi0]
Ci   f     :contribution to forces is added
Cr Remarks
Cu Updates
Cu   02 Jul 15 Remove fixed dimensioning ndim
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   22 Apr 00 Adapted from nfp ugcomp
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters

C#ifdefC MPE
C      include "mpef.h"
C#endif
      integer procid, master, numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
C     character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1),lgunit
C     double precision starttime, endtime
      logical mlog,cmdopt
      character*120 strn
C#ifdefC MPE
CC Event numbers:
C      include "events.ins"
C#endif
      integer nbas
      double precision qmom(*),gpot0(*),f(3,nbas),hpot0(nbas),ugg
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      real(8),allocatable :: xgpot0(:,:)
      complex(8),allocatable :: s(:,:),ds(:,:,:)
C ... Local parameters
      integer ndim0,i,ib,ilm1,ilm2,is,iv0,jb,js,jv0,nvl,l1,l2,
     .  lfoc1,lfoc2,ll,lmax1,lmax2,m,nlm1,nlm2
      parameter (ndim0=2)
      double precision ceh1,ceh2,cof1,cof2,cofg1,cofg2,cofh1,cofh2,fpi,
     .  pi,qcorg1,qcorg2,qcorh1,qcorh2,qsc1,qsc2,qm1,qm2,rg1,rg2,rh1,
     .  rh2,srfpi,y0,z1,z2
      double precision df(0:20),ff(3),tau1(3),tau2(3)
      double complex s0(ndim0,ndim0),
     .   ds0(ndim0,ndim0,3),wk(ndim0,ndim0),dwk(ndim0,ndim0,3)
C ... For parallel threads
      integer nlmx,npmx,ip,mp,nbmx
C#ifndef SGI-PARALLEL
      parameter (npmx=1, nbmx=256)
      double precision xf(3,nbas,npmx),xhpot0(nbas,npmx),xugg(npmx)
C#elseC
C      parameter (npmx=32)
C      double precision xf(3,nbas,npmx),xhpot0(nbas,npmx),xugg(npmx)
C#endif
      integer , dimension(:), allocatable :: bproc
      double precision , dimension(:), allocatable :: buffer
      integer nvl0,iiv0(nbas)

      call tcn('ugcomp')
      call mpi_comm_rank( mpi_comm_world, procid, ierr )
      call mpi_comm_size( mpi_comm_world, numprocs, ierr )
      master = 0
      if (numprocs > 1) then
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      mlog = cmdopt('--mlog',6,0,strn)
      endif

      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      srfpi = dsqrt(fpi)
      y0 = 1d0/srfpi
      nlmx = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        nlmx = max(nlmx,(s_spec(is)%lmxl+1)**2)
      enddo
      allocate(s(nlmx,nlmx),ds(nlmx,nlmx,3))

C ... Setup array iiv0 = (vector of iv0 for parallel); allocate work arrays
      mp = 1

      if (numprocs > 1) then
      call setofl(0,s_site,s_spec,nbas,nvl0,iiv0)
      if (nlmx*nbas < nvl0) call rx('ugcomp: increase nlmx')
      endif

      if (npmx < mp) call rxi('ugcomp: increase npmx, needed',mp)

C --- Loop over sites where charge lump making pot is centered ---
      ugg = 0d0
      iv0 = 0
      ip = 1
      call dpzero(xugg, mp)
      allocate(xgpot0(nlmx*nbas,npmx))
      call dpzero(xgpot0, nlmx*nbas*mp)
      call dpzero(xf, 3*nbas*mp)
      call dpzero(xhpot0, nbas*mp)
      allocate (bproc(0:numprocs))
      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_UGCOMP,procid,"ugcomp")
C#endif
      call dstrbp(nbas,numprocs,1,bproc(0))
      else
      bproc(0:1) = [1,nbas+1]
      endif
      do  ib = bproc(procid), bproc(procid+1)-1
        if (numprocs > 1) then
        if (mlog .and. ib == bproc(procid)) then
          call gettime(datim)
          call awrit4(' ugcomp '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting atoms %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        iv0 = iiv0(ib)
        end if
        is = s_site(ib)%spec
        tau1 = s_site(ib)%pos
        lmax1 = s_spec(is)%lmxl
        rg1 = s_spec(is)%rg
        call corprm(s_spec,is,qcorg1,qcorh1,qsc1,cofg1,cofh1,ceh1,lfoc1,rh1,z1)
        nlm1 = (lmax1+1)**2

C   ... Loop over sites where charge lump sees the potential
        if (lmax1 > -1) then
        jv0 = 0
        do  jb = 1, nbas
          js = s_site(jb)%spec
          lmax2 = s_spec(js)%lmxl
          if (lmax2 < 0) cycle
          tau2 = s_site(jb)%pos
          rg2 = s_spec(js)%rg
          call corprm(s_spec,js,qcorg2,qcorh2,qsc2,cofg2,cofh2,ceh2,lfoc2,rh2,z2)
          nlm2 = (lmax2+1)**2
          if (numprocs > 1) jv0 = iiv0(jb)

          if (nlm1 > nlmx) call rxi('ugcomp: nlmx < nlm1=',nlm1)
          if (nlm2 > nlmx) call rxi('ugcomp: nlmx < nlm2=',nlm2)
          call ggugbl(tau1,tau2,rg1,rg2,nlm1,nlm2,nlmx,nlmx,s_lat,s,ds)

          ff(1) = 0d0
          ff(2) = 0d0
          ff(3) = 0d0
          do  ilm1 = 1, nlm1
            l1 = ll(ilm1)
            qm1 = qmom(iv0+ilm1)
            if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
            cof1 = qm1*fpi/df(2*l1+1)
            do  ilm2 = 1, nlm2
              l2 = ll(ilm2)
              qm2 = qmom(jv0+ilm2)
              if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
              cof2 = qm2*fpi/df(2*l2+1)
              xugg(ip) = xugg(ip) + cof1*cof2*s(ilm1,ilm2)
              xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip)
     .                            + s(ilm1,ilm2)*cof1*fpi/df(2*l2+1)
C         ... Forces
              ff(1) = ff(1) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,1)
              ff(2) = ff(2) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,2)
              ff(3) = ff(3) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,3)
            enddo
          enddo

C     --- Additional h*h, h*g, g*h terms for foca ---
          if (lfoc1 > 0 .or. lfoc2 > 0) then
            call hhugbl(0,tau1,tau2,rh1,rh2,ceh1,ceh2,1,1,ndim0,ndim0,
     .        s_lat,wk,dwk,s0,ds0)
            xugg(ip) = xugg(ip) + cofh1*s0(1,1)*cofh2
            xhpot0(jb,ip) = xhpot0(jb,ip) + cofh1*s0(1,1)
            ff(1) = ff(1) + 0.5d0*cofh1*cofh2*ds0(1,1,1)
            ff(2) = ff(2) + 0.5d0*cofh1*cofh2*ds0(1,1,2)
            ff(3) = ff(3) + 0.5d0*cofh1*cofh2*ds0(1,1,3)

            call hgugbl(tau1,tau2,rh1,rg2,ceh1,1,nlm2,nlmx,nlmx,s_lat,s,ds)
            do  ilm2 = 1, nlm2
              l2 = ll(ilm2)
              qm2 = qmom(jv0+ilm2)
              if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
              cof2 = qm2*fpi/df(2*l2+1)
              xugg(ip) = xugg(ip) + cofh1*s(1,ilm2)*cof2
              ff(1) = ff(1) + 0.5d0*cofh1*cof2*ds(1,ilm2,1)
              ff(2) = ff(2) + 0.5d0*cofh1*cof2*ds(1,ilm2,2)
              ff(3) = ff(3) + 0.5d0*cofh1*cof2*ds(1,ilm2,3)
              xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip)
     .                            + s(1,ilm2)*cofh1*fpi/df(2*l2+1)
            enddo

            call hgugbl(tau2,tau1,rh2,rg1,ceh2,1,nlm1,nlmx,nlmx,s_lat,s,ds)
            do  ilm1 = 1, nlm1
              l1 = ll(ilm1)
              qm1 = qmom(iv0+ilm1)
              if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
              cof1 = qm1*fpi/df(2*l1+1)
              xugg(ip) = xugg(ip) + cof1*s(1,ilm1)*cofh2
              ff(1) = ff(1) - 0.5d0*cof1*cofh2*ds(1,ilm1,1)
              ff(2) = ff(2) - 0.5d0*cof1*cofh2*ds(1,ilm1,2)
              ff(3) = ff(3) - 0.5d0*cof1*cofh2*ds(1,ilm1,3)
              xhpot0(jb,ip) = xhpot0(jb,ip) + cof1*s(1,ilm1)
            enddo
          endif

          if (jb /= ib) then
            do  m = 1, 3
              xf(m,ib,ip) = xf(m,ib,ip) - ff(m)
              xf(m,jb,ip) = xf(m,jb,ip) + ff(m)
            enddo
          endif

          jv0 = jv0+nlm2
        enddo
        iv0 = iv0+nlm1
        endif
      enddo
      if (numprocs > 1) then
      nvl = nvl0
      else
      nvl = iv0
      endif
      deallocate(s,ds)

C ... Assemble data from separate threads
      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_UGCOMP,procid,"ugcomp")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
C#endif
      allocate(buffer(1:nvl), stat=ierr)
      call MPI_ALLREDUCE(xgpot0,buffer,nvl,
     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' ugcomp '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce gpot0 nvl=%i',' ',256,lgunit(3),
     .        procid,numprocs,nvl)
      endif
      call daxpy(nvl,1d0,buffer,1,gpot0,1)
      deallocate(buffer, stat=ierr)

      allocate(buffer(1:nbas), stat=ierr)
      call MPI_ALLREDUCE(xhpot0,buffer,nbas,
     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' ugcomp '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce hpot0 nbas=%i',' ',256,lgunit(3),
     .        procid,numprocs,nbas)
      endif
      call daxpy(nbas,1d0,buffer,1,hpot0,1)
      deallocate(buffer, stat=ierr)

      allocate(buffer(1:3*nbas), stat=ierr)
      call MPI_ALLREDUCE(xf,buffer,3*nbas,
     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' ugcomp '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce f 3nbas=%i',' ',256,lgunit(3),
     .        procid,numprocs,3*nbas)
      endif
      call daxpy(3*nbas,1d0,buffer,1,f,1)
      deallocate(buffer, stat=ierr)

      call MPI_ALLREDUCE(xugg,ugg,1,
     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit2(' ugcomp '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce ugg',' ',256,lgunit(3),
     .        procid,numprocs)
      endif
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
C#endif
      deallocate(bproc, stat=ierr)
      else
      do  ip = 1, mp
        do  ib = 1, nbas
          f(1,ib) = f(1,ib) + xf(1,ib,ip)
          f(2,ib) = f(2,ib) + xf(2,ib,ip)
          f(3,ib) = f(3,ib) + xf(3,ib,ip)
          hpot0(ib) = hpot0(ib) + xhpot0(ib,ip)
        enddo
        do  i = 1, nvl
          gpot0(i) = gpot0(i) + xgpot0(i,ip)
        enddo
        ugg = ugg + xugg(ip)
      enddo
      endif

      deallocate(xgpot0)
      call tcx('ugcomp')
      end
