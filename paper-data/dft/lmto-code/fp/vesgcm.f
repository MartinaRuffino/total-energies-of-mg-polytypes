      subroutine vesgcm(job,s_site,s_spec,s_lat,nbas,qmom,ng,
     .  cv,cgsum,f,gpot0,hpot0,qsmc,zsum,vrmt)
C- Adds contribution from compensating gaussians to smooth estat potential
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl rg rmt lfoca rfoca qc z ctail etail stc lmxb p
Ci                 pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  corprm
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   job   :0 Make cgsum, gpot0, hpot0, qsmc, zsum; add to f
Ci         :1 Also add potential from n0~-n0 into cv, and make vrmt
Ci   nbas  :size of basis
Ci   qmom  :multipole moments of on-site the local density:
Ci         :integral r^l (rho1-rho2) + l=0 contr. from core spillout
Ci   ng    :number of G-vectors
Co Inputs/Outputs
Cio  cv    :On input, smpot phi[n0], stored in packed form (gvlist)
Cio        :On output, estat potential from compensating (gaussians + hankels) is added
Co Outputs
Co   cgsum :G-rep of compensating gaussian density, in packed form
Co   f     :contribution from compensating gaussians is added to force
Co   gpot0 :vector of integrals [ (compensating gaussians g_RL) * phi0]
Co         :Here phi0 = phi[n0], where n0 = sm rho w/out gaussians.
Co         :Note : The integrals for total energy are [ g_RL * phi0~] where
Co         :phi0~ - phi0  = potential from the compensating gaussians
Co         :g_RL * [phi0~ -phi0] is computed elsewhere (in routine ugcomp) when
Co         :computes g_RL * [phi0~-phi0] from structure constants to improve accuracy.
Co   hpot0 :integrals of smooth Hankels * phi0 (part of core contribution)
Co   qsmc  :smoothed core charge
Co   zsum  :total nuclear charge
Co   vrmt  :(job=1) potential at MT boundary
Cl Local variables
Cl   qcorg :gaussian part of the pseudocore charge; see corprm.f
Cl   qcorh :hankel part of the pseudocore charge; see corprm.f
Cr Remarks
Cr   Local charge consists of a sum of gaussians that compensate for
Cr   the difference in multipole moments of true and smooth local charge
Cr   and a contribution from the smooth core charge.
Cr     g(qmpol) + g(qcore-z) + h(ncore)
Cu Updates
Cu   02 Jul 16 cgsum(1) now returned as (compensating gaussian charge)/vol
Cu   16 Jun 16 New job; shorten argument list
Cu   02 Jul 15 Removed fixed dimension nlmx
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
Cu   20 Apr 01 Generates vrmt
Cu   23 Apr 00 Adapted from nfp ves_gcomp.f
C ----------------------------------------------------------------------
      use structures
      use mpi

      implicit none
C ... Passed parameters
      integer job,nbas,ng
      double precision qsmc,zsum
      double precision qmom(*),f(3,nbas),gpot0(*),hpot0(nbas),vrmt(nbas)
      double complex cv(ng),cgsum(ng)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      real(8),allocatable :: cof(:)
      complex(8),allocatable :: gkl(:,:),cg1(:)
C ... Local parameters
      integer stdo,k0,i,ib,ilm,is,iv0,kb,l,lgunit,ll,lmxl,m,nlm,lfoc
      parameter (k0=3)
      double precision alat,ceh,cofg,cofh,g2,pi,qc,qcorg,qcorh,qsc,rfoc,
     .  rg,sum1,sum2,sum3,tpiba,vol,xx,y0,z
      double precision df(0:20),tau(3),v(3),ddot,gvr,rmt,fac,gvb
      double complex img

      integer :: iproc, nproc, comm, err, fg_mt, nlms(0:nbas), j
      integer, allocatable :: pmap(:), smap(:), nmap(:)
      real(8) :: fg(4,nbas), fgav(4)


      call tcn('vesgcm')

      call stdfac(20,df)
      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      alat = s_lat%alat
      vol = s_lat%vol
      tpiba = 2*pi/alat

C     call zprm('cv',2,cv,ng,ng,1)

C --- Accumulate FT of Gaussian + Hankel density for listed vectors ---
C     and make integrals g*phi0, h*phi0
      allocate(cg1(ng))
      call dpzero(cgsum,2*ng)

      qsmc = 0d0
      zsum = 0d0
      fgav = 0

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

      do  ib = pmap(iproc)+1, pmap(iproc+1)
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        lmxl = s_spec(is)%lmxl
        rg = s_spec(is)%rg
        fg(:,ib) = 0
        if (lmxl == -1) cycle
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        qc = qcorg+qcorh
        qsmc = qsmc+qc
        zsum = zsum+z
        nlm = (lmxl+1)**2
        allocate(gkl(0:k0,nlm),cof(nlm))
        ilm = 0
        iv0 = nlms(ib-1)
        do  l = 0, lmxl
          do m = -l,l
            ilm = ilm+1
            cof(ilm) = qmom(ilm+iv0)*4*pi/df(2*l+1)
            gpot0(ilm+iv0) = 0d0
          enddo
        enddo
        hpot0(ib) = 0d0
        cof(1) = cof(1) + 4*pi*y0*(qcorg-z)

        call dpzero(cg1,2*ng)
        do  i = 1, ng
          v = s_lat%gv(i,1:3)
          call gklft(v,rg,0d0,tau,alat,0,nlm,k0,s_lat%cy,gkl) ! FT of compensating gaussians g_L at tau
          do  ilm = 1, nlm
            cg1(i) = cg1(i) + cof(ilm)*gkl(0,ilm)/vol
            gpot0(ilm+iv0) = gpot0(ilm+iv0) + dconjg(cv(i))*gkl(0,ilm) ! g_L * V[n0]
          enddo
          call hklft(v,rfoc,ceh,tau,alat,0,1,k0,s_lat%cy,gkl) ! FT of hankels h_L at tau for pseudo-nuc
          cg1(i) = cg1(i) + cofh*gkl(0,1)/vol
          hpot0(ib) = hpot0(ib) + dconjg(cv(i))*gkl(0,1) ! h_L * V[n0]
        enddo
        call dpadd(cgsum,cg1,1,2*ng,1d0)  ! cgsum += qmom_L*V[n0]*g_L + cofh_0*V[n0]*h_0

C   ... Multiply factors into gpot0
        do  ilm = 1, nlm
          l = ll(ilm)
          gpot0(ilm+iv0) = gpot0(ilm+iv0)*4*pi/df(2*l+1)
        enddo

C   ... Force of smooth density on the compensating gaussians
        do  i = 1, ng
          xx = -dimag(dconjg(cv(i))*cg1(i))
          do j = 1, 3
            fg(j,ib) = fg(j,ib) + xx*s_lat%gv(i,j)
          end do
        enddo
        fg(:,ib) = fg(:,ib)*vol*tpiba

        fgav = fgav + fg(:,ib)

        deallocate(gkl,cof)
      enddo

      deallocate(cg1)

      if (nproc > 1) then
        allocate(smap(0:nproc-1))
        smap = pmap(1:nproc) - pmap(0:nproc-1)

        call mpi_type_contiguous(4, mpi_real8, fg_mt, err)
        call mpi_type_commit(fg_mt, err)

!         call mpi_allreduce(mpi_in_place, fgav, 1, fg_mt, mpi_sum, comm, err) ! mpi_sum not defined for custom types
        call mpi_allreduce(mpi_in_place, fgav, 4, mpi_real8, mpi_sum, comm, err)
      end if

      fgav = fgav/nbas

      do ib = pmap(iproc)+1, pmap(iproc+1)
        fg(:,ib) = fg(:,ib) - fgav
      end do

      if (nproc > 1) then
        call mpi_allgatherv(mpi_in_place, smap(iproc), fg_mt, fg, smap, pmap, fg_mt, comm, err)
        call mpi_allgatherv(mpi_in_place, smap(iproc), mpi_real8, hpot0, smap, pmap, mpi_real8, comm, err)

        call mpi_type_free(fg_mt, err)

        allocate(nmap(0:nproc))

        nmap = nlms(pmap)
        smap = nmap(1:nproc) - nmap(0:nproc-1)

        call mpi_allgatherv(mpi_in_place, smap(iproc), mpi_real8, gpot0, smap, nmap, mpi_real8, comm, err)

        deallocate(nmap)
        deallocate(smap)

        call mpi_allreduce(mpi_in_place, cgsum, ng, mpi_complex16, mpi_sum, comm, err)
        call mpi_allreduce(mpi_in_place, zsum, 1, mpi_real8, mpi_sum, comm, err)
        call mpi_allreduce(mpi_in_place, qsmc, 1, mpi_real8, mpi_sum, comm, err)
      end if

      do ib = 1, nbas
        f(1:3,ib) = f(1:3,ib) + fg(1:3,ib)
      end do

      if (job == 0) then
        deallocate(pmap)
        return
      end if

C --- Add 8pi/G**2 * (FT gaussian+Hankel density) into smpot ---
      cgsum(1) = cgsum(1) - (qsmc-zsum)/vol
C      call info5(40,1,0,' vesgcm: compensating density'//
C     .  ' Qval = %;6,6d  core-nuc = %;6;6d  tot = %;6;6d',
C     .  dble(cgsum(1)*vol),qsmc-zsum,dble(cgsum(1)*vol)+qsmc-zsum,0,0)
      do  i = 2, ng
        g2 = tpiba*tpiba*(s_lat%gv(i,1)**2+s_lat%gv(i,2)**2+s_lat%gv(i,3)**2)
        cv(i) = cv(i) + (8*pi)*cgsum(i)/g2
      enddo

C --- Electrostatic potential at rmt ---
      call dpzero(vrmt,nbas)
      img = (0d0,1d0)
      do ib = pmap(iproc)+1, pmap(iproc+1)
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        call dscal(3,alat,tau,1)
        rmt = s_spec(is)%rmt
C       Add a negligibly small amount to rmt to handle case rmt=0
        rmt = rmt+1d-32
        do  i = 2, ng
          v(1) = s_lat%gv(i,1)*tpiba
          v(2) = s_lat%gv(i,2)*tpiba
          v(3) = s_lat%gv(i,3)*tpiba
          g2 = dsqrt(ddot(3,v,1,v,1))
          gvr = g2*rmt
          fac = sin(gvr)/gvr
          gvb = v(1)*tau(1) + v(2)*tau(2) + v(3)*tau(3)
          vrmt(ib) = vrmt(ib) + dble(cv(i)*fac*exp(img*gvb))
        enddo
      enddo

      if (nproc > 1) then
! doto: this is likely unneccessary since each process probably only accesses its own vrmt entries. lone mpi_reduce in mkpot or smves might work fine..
        allocate(smap(0:nproc-1))
        smap = pmap(1:nproc) - pmap(0:nproc-1)

        call mpi_allgatherv(mpi_in_place, smap(iproc), mpi_real8, vrmt, smap, pmap, mpi_real8, comm, err)

        deallocate(smap)
      end if

      deallocate(pmap)

      call tcx('vesgcm')
      end
