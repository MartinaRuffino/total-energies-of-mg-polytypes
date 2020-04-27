      subroutine ovlpfa(s_site,s_lat,s_ham,ib1,nbas,nxi,nxi0,exi,hfc,
     .  rsmfa,ng,ngmx,gv,cv)
C- Set up Fourier coeffs to overlap the smooth part of FA densities.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   nxi   :number of Hankel functions per l for c.d. basis, by species
Ci   nxi0  :leading dimension of hfc
Ci   exi   :smoothed Hankel energies for each function in c.d. basis, by species
Ci   hfc   :coefficients to smoothed Hankels
Ci   rsmfa :Hankel smoothing radius
Ci   ng    :number of G-vectors
Ci   ngmx  :leading dimension of gv
Ci   gv    :list of reciprocal lattice vectors G (glist.f)
Co Outputs
Co   cv    :Fourier coefficients to smooth density
Cl Local variables
Cl   ns4 = 1 if local density is not spin polarized
Cl         2 if local spin density is collinear along z
Cl         4 if local spin density is to be rotated;
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 May 07 parallelized (MPI)
Cu   01 Jul 05 Zero-radius empty spheres skip as having no local part
Cu   13 Jun 00 spin polarized
Cu   24 Apr 00 Adapted from nfp ovlpfa.f
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
      integer ib1,nbas,nxi(*),nxi0,ng,ngmx
      double precision gv(ngmx,3),rsmfa(*),exi(nxi0,*),
     .  hfc(nxi0,2,*)
      complex(8),target :: cv(ng,4)  ! Maybe be as many as 4 if nsp=nspc=2
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Local parameters
      integer stdo,ipr,iprint,lgunit,ib,is,nx,i,ixi,ig,nsp,nspc,nglob,
     .  isp,ns4
      double precision v(3),pi,y0,alat,vol,tpiba,sum(2),px,py,pz,pos(3),
     .  sam(2),e,cof,rsm,gam,v2,aa,scalp
      complex(8),pointer :: cvl(:,:)
      double complex phase,cv0(2)
      equivalence (px,pos(1)),(py,pos(2)),(pz,pos(3))
      integer, dimension(:),allocatable :: kpproc
      integer ierr
C     integer ifi,fopna
      integer procid, master, numprocs
      logical mlog,cmdopt
      character strn*120

      call tcn('ovlpfa')
      ipr  = iprint()
      stdo = lgunit(1)
      nsp  = nglob('nsp')
      nspc  = nglob('nspc')
      ns4 = nsp*nspc; if (mod(s_ham%lrsa,10) == 0) ns4 = nsp ! See "Local variables"
      pi   = 4d0*datan(1d0)
      y0   = 1d0/dsqrt(4*pi)
      alat = s_lat%alat
      vol = s_lat%vol
      tpiba = 2*pi/alat
      call dpzero(cv,2*ng*ns4)
C     Accumulate into cvl for one site; add to cv later
C     Enables site-specific rotation of cvl to noncollinear axis
      if (ns4 == 4) then
        allocate(cvl(ng,ns4))
      else
        cvl => cv
      endif

      call mpi_comm_rank(mpi_comm_world, procid, ierr)
      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      master = 0
      if (numprocs > 1) mlog = cmdopt('--mlog',6,0,strn)

C --- Loop over sites ---
      sum(1) = 0d0
      sum(2) = 0d0
      call info0(31,1,0,' ovlpfa: overlap smooth part of FA densities')
      allocate (kpproc(0:numprocs))
      if (numprocs > 1) then
        call pshpr(ipr-10)
        call dstrbp(nbas,numprocs,1,kpproc(0))
        call poppr
        ipr = 0
      else
        kpproc(0:1) = [ib1,nbas+1]
      endif
      cv0 = 0
      do ib = kpproc(procid), kpproc(procid+1)-1
        is = s_site(ib)%spec
        pos = s_site(ib)%pos
        nx = nxi(is)

        if (ns4 == 4) then
          call dpzero(cvl,2*ng*ns4)
        endif

C   ... Loop over Hankels at this site
        sam(1) = 0
        sam(2) = 0
        do  isp = 1, nsp
          do  ixi = 1, nx
            e = exi(ixi,is)
            cof = hfc(ixi,isp,is)
            rsm = rsmfa(is)
            gam = 0.25d0*rsm**2
            sam(isp) = sam(isp) - cof*y0*4d0*pi*dexp(gam*e)/e
C       ... Loop over reciprocal lattice vectors
            do  ig = 1, ng
              v(1) = gv(ig,1)*tpiba
              v(2) = gv(ig,2)*tpiba
              v(3) = gv(ig,3)*tpiba
              v2 = v(1)**2+v(2)**2+v(3)**2
              aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
              scalp = -alat*(px*v(1)+py*v(2)+pz*v(3))
              phase = dcmplx(dcos(scalp),dsin(scalp))
              cvl(ig,isp) = cvl(ig,isp) + cof*aa*phase*y0/vol
            enddo
          enddo
          sum(isp) = sum(isp) + sam(isp)
        enddo
        if (ipr > 30 .and. nx > 0) then
          call awrit6(' site%,4i  spec%,3i  pos%3;8,4D'//
     .      '  Qsmooth %,6;6d%?#n==2#  mom %,5;5d##',
     .      ' ',90,stdo,ib,is,pos,sam(1)+sam(2),nsp,sam(1)-sam(2))
          if (ipr >= 40) then
          write(stdo,700) 'energy:',(exi(i,is),i=1,nx)
          write(stdo,700) 'coeff:',(hfc(i,1,is),i=1,nx)
          if (nsp == 2) write(stdo,700) 'spin2:',(hfc(i,2,is),i=1,nx)
          write(stdo,'(1x)')
          endif
  700   format(2x,a7,6f11.3)
        endif

        cv0(1:nsp) = cvl(1,1:nsp)
        if (ns4 == 4) then
          if (s_ham%neula/=1) call rx('not ready for l-dependent eula')
          if (ipr >= 40)
     .      write(stdo,700) 'Euler:',(s_ham%eula(ib,ixi),ixi=1,3)
          call rotspv(613,s_ham%eula(ib,1:3),ng,ng,ng,cvl,cvl,cvl,cvl)
C         call rotspu(0,1,1,1,1,s_ham%eula(ib,1),1,u)
C         call rotspv(13,u,ng,ng,ng,cvl,cvl,cvl,cvl)
C         call zprm('cv-nc',2,cvl,ng,ng,4)
          call daxpy(2*ng*ns4,1d0,cvl,1,cv,1)
        endif

      enddo

C ... Combine cv from separate threads
      if (numprocs > 1) then
      deallocate(kpproc, stat=ierr)
      call mpibc2(cv,ng*ns4*2,4,3,mlog,'ovlpfa','cv')
      call mpibc2(sum,2,4,3,mlog,'ovlpfa','sum')
      call mpibc2(cv0,nsp,6,3,mlog,'ovlpfa','cv0')

C     Debugging
C     ifi = fopna('out',-1,0)
C     call ywrm(0,'cv',3,ifi,'(9f20.10)',cv,1,ng,ng,1)
C     call rx0('done')
      endif

      ipr  = iprint()
      call info5(31,0,0,' total smooth Q = %,6;6d'//
     .  '%?#n==2#  moment = %,5;5d#%j#%?#n>40#  FT (0,0,0) = %,6;6d',
     .  sum(1)+sum(2),nsp,sum(1)-sum(2),ipr,
     .  (cv0(1)+cv0(nsp))*vol/(3-nsp))

C     if (.not. associated(cvl,cv)) deallocate(cvl)
      if (ns4 == 4) deallocate(cvl)

      call tcx('ovlpfa')

      end
