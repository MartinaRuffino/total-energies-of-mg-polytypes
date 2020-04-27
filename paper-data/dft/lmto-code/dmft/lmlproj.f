      subroutine lmlproj(mode,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,ef0)
C- Make projectors from LMTO eigenfunctions
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  star
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lshft nkabc ipq
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ndham ldham lsig pwmode
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:offH iprmb
Cio    Passed to:  makusq
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  napw alat qlat vol plat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed: plat pos istab symgr ag napw igv2 qlat cg indxcg
Cio                jcg cy qlv dlv
Cio    Passed to:  makusq pusq1 bstrux hxpbl ghibl hklbl gklbl hxpgbl
Cio                ghigbl hklgbl
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  ppn
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rsma lmxa lmxl kmxt a nr rmt lmxb pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  makusq uspecb pusq1 bstrux
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  makusq pusq1 bstrux
Cio  s_dmft
Ci     Elts read:  ndsigi nicix ncix nlohi knorm icix l dcmode ib
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ndim iproj
Cio    Passed to:  makeproj
Ci Inputs
Ci    mode: 1s digit
Ci        : 0 => qp is generated for the given nkp qp only.
Ci        : 1 => qp is generated for the full BZ (sets lfbz; see Local variables)
Ci        : 10s digit
Ci        : 0 => allocate overlap matrix of projectors locally.
Ci        :      Note: when s_dmft%knorm=0, this is the only mode
Ci        :      since the overlap matrix is made for each k.
Ci        : 1 => use s_dmft%Olapp as memory address for overlap matrix
Ci        : 2 => s_dmft%Olapp is projector overlap matrix, and it is already
Ci        :      already generated
Ci        :  100s digit
Ci        : 0 => write proj in binary format
Ci        : 1 => write proj in binary+hdf5 format
Co Outputs
Ci   sigdc:
Cs Command-line switches
Cl Local variables
Cl   lstar : F => qp is generated for the given nkp points only.
Cl         : T => qp is generated for the full BZ.
Cl         :      Requires that qp are on a regular mesh, set by
Cl         :      s_bz%nkabc.  points in the full BZ are generated
Cl         :      by rotation to the star of k
Cl   lolap : 10s digit of mode
Cr Remarks
Cr   Requires that evecs reside on disk in file 'evec'
Cu Updates
Cu   17 Jan 18 (MvS) First cut at SO coupled case
Cu   28 Sep 17 (MvS) MPI parallel
Cu   16 Apr 17 Adapted from sudmft
C ----------------------------------------------------------------------
      use structures
      use h5
      use mpi
      implicit none
C ... Passed parameters:
      integer mode
      real(8) :: ef0
C ... For structures
!     include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_dmft)::  s_dmft
      type(h5space) :: pfs_h5,efs_h5, ms_h5,s_h5
      type(h5dataset) :: pdset, edset
      type(h5file) :: pf

C ... Dynamically allocated local arrays
      integer, allocatable :: kpproc(:) ! For MPI distribution of processes
      real(8), allocatable :: evlk(:,:) ! Eigenvalues for 1 k-point
      real(8), allocatable :: evlkh5(:,:) ! Eigenvalues for 1 k-point
      real(8), pointer     :: ppnl(:)   ! NMTO-like potential parameters
      real(8), allocatable :: rmat(:)   ! Rotation matrix for rotating eigenfunctions
      complex(8), allocatable, target :: h(:,:) ! Work space for call to hambls
      complex(8), pointer :: z(:,:)     ! Eigenvectors
      complex(8), allocatable :: ausp(:,:,:,:,:)  ! val,slo of w.f. at MT sphere surface for projector
!     integer, pointer :: ndsigi(:) ! Points to s_dmft%ndsigi
      complex(8), pointer     :: dmftu(:,:,:,:),dmftuq(:,:,:,:,:,:) ! Factored form of projector
      complex(8), allocatable :: dmftUi(:,:,:)   ! Locally dimensioned form of dmftu needed for file I/O
      complex(8), pointer     :: Olapm(:,:,:,:)  ! Local version of overlap matrix
      complex(8), pointer     :: Olapmk(:,:,:,:) ! Local version of k-resolved overlap matrix (knorm=1)
      complex(8), pointer     :: projemb(:,:)    ! projected embedded identity.

!     complex(8), allocatable :: sig_padcof(:,:,:,:)  ! Pade coefficients for interpolating DOS
      complex(8), pointer :: zr(:,:) !Rotated eigenvectors for symmetry reduced BZ mesh lmfdmft

C ... Local variables
      integer, parameter :: NULLI=-99999, PRTU=50, PRTG=70
      integer, parameter :: lSOf=4,lSzLz=32,lSzLzp=64,lSzLzp0=128
      integer, parameter :: n0=10, nppn=12, nkap0=4, LW5=5
      real(8), parameter :: ry2ev = 13.60569193d0
      real(8), parameter :: pi = acos(-1d0)
      real(8),parameter  :: NULLR = -99999
      complex(8), parameter :: srm1=(0d0,1d0),NULLZ=(-99999D0,0d0)
      integer, parameter :: rfalsi_limit=100
      integer, parameter :: lchk_sigbar=3    ! print check files on  P+E[sigbar(infty)] HARD CODED
      integer, parameter :: lchk_prjnrm=7    ! print check files on Normalized Projectors HARD CODED
      integer procid,master,nproc            ! For MPI

      character fn*128,extmpi*6
      logical lstar
      integer cix,icix,i,ifi,ifiz,ig,ipass,iq,iqs,isp,istar,j,jsp,ksp,k,ldim,lrsig,nbas,ndham,
     .  ndhamx,ndimh,ndimhx,nev,nevn,ncix,nl,nl1,nlmax,nphimx,nsp,nspc,lolap
      integer rank                  ! MPI rank, used when merging projector files
      integer nstarq                ! number of stars at this irreducible qp
      integer ldcix                 ! dimensioning parameter for projector matrices
      integer jfi                   ! MPI case, file handle for proj file
      integer oldrank               ! MPI rewrite of proj file; to track of when file changes
      double precision qp(3),xv(20)
      integer :: ifac(3), i123(3)               ! used for iqstar
      double precision qb(9),qpr(3) ! used for iqstar
      complex(8) :: facalpha,facbeta
C     complex(8) :: chk(1000,2)      ! for debugging sanity check

      integer lshft(3),nk1,nk2,nk3,nkabc(3),nkp,nqi,nkfbz,nlohi(2)
      equivalence (nk1,nkabc(1)), (nk2,nkabc(2)), (nk3,nkabc(3))
      double precision eseavr(2) ! Constant term for high-lying QSGW sigma
      integer :: lcix(s_dmft % ncix)
      procedure(logical) :: a2bin,cmdopt
      procedure(real(8)) :: ddot,dlength,cpusec
      procedure(integer) :: nglob,fopna,fopng,fopnx,iprint,mpipid

      !     for hdf5
      logical :: save_h5
      integer :: comm, iqfbz
      real(8) :: lstk(3,s_bz%nkabc(1)*s_bz%nkabc(2)*s_bz%nkabc(3))
      integer,allocatable :: wtk(:)
      integer :: u ! unit to use in open(newunit=u, ...)
!     integer :: debug = 0

C ... Setup
      call tcn('lmlproj')
      procid = mpipid(1); master = 0
      nl = nglob('nl')
      nbas = nglob('nbas')
      nsp  = nglob('nsp')
      nspc = nglob('nspc')      ! 2 for noncollinear case
      ndham = s_ham%ndham
      ndhamx = s_ham%ndham * nspc
      nphimx = nglob('nphimx')  ! number of envelope function types joined to (u,s,phiz)
      ldim  = s_ham%ldham(1)
      ldcix = maxval(s_dmft%ndim)
      nlmax = nglob('nlmax')
      if (s_lat%napw /= 0) call rx('lmlproj : no PMT method yet')
      ncix = s_dmft%ncix
      nlohi = s_dmft%nlohi
      nevn  = nlohi(2)-nlohi(1)+1
      lrsig = s_ham%lsig
      lstar = mod(mode,10) == 1
      lolap = mod(mode/10,10)
      if (s_dmft%knorm == 1) lolap = 0
      ppnl => s_pot%ppn
      allocate(evlk(ndham,nsp))
      allocate(evlkh5(nevn,nsp))
      allocate(z(ndhamx,ndhamx),h(ndhamx,ndhamx))
      save_h5=.false.
      if (mod(mode/100,10) == 1) then
         open(newunit=u, file='solver_input.h5')
         close(u,status='delete')
         save_h5=.true.
      endif

      if (lolap == 0) then
        allocate(Olapm(ldcix,ldcix,nsp,ncix))
      else
        Olapm => s_dmft%Olapp
      endif
      if (lolap < 2) then
        call dpzero(Olapm,2*size(Olapm))
      endif
      allocate(Olapmk(ldcix,ldcix,nsp,ncix)) ! debugging check to confirm that U U+ is orthonormal
      if (size(Olapm) /= size(Olapmk)) call rx('s_dmft%Olapp improperly dimensioned')
      call dpzero(Olapmk,2*size(Olapmk))
      if (lstar) call bzmsh00(s_lat%plat,s_bz%lshft,0,s_bz%nkabc,ifac,qb) ! Needed for iqstar

C     Open evecs file, read in qp data. Header is read from master only.
      ifiz = fopna('evec',-1,4); rewind ifiz
      if (procid == master) then
        call iosigh(2,LW5,i,j,k,k,nkabc(1),nkabc(2),nkabc(3),nkp,nqi,lshft(1),lshft(2),lshft(3),ifiz,eseavr)
      endif
      call mpibc1(i,1,2,.false.,'','')
      call mpibc1(j,1,2,.false.,'','')
      call mpibc1(k,1,2,.false.,'','')
      call mpibc1(nkabc,3,2,.false.,'','')
      call mpibc1(nkp,1,2,.false.,'','')
      call mpibc1(lshft,3,2,.false.,'','')

      if (save_h5) then
         if (procid==master) then
            open(file='proj.h5',unit=123)
            close(123,status='delete')
            call h5_write('h_lm.h5:/ldcix', ldcix, find_h5dtype(ldcix))
            call h5_write('h_lm.h5:/ncix',  ncix, find_h5dtype(ncix))
            call h5_write('h_lm.h5:/nsp',    nsp, find_h5dtype(nsp))
            call h5_write('h_lm.h5:/nkfbz', [s_bz%nkabc(1)*s_bz%nkabc(2)*s_bz%nkabc(3)], find_h5dtype(s_bz%nkabc(1)))

            call h5_write('h_lm.h5:/nkabc', s_bz%nkabc, find_h5dtype(s_bz%nkabc))
            call h5_write('h_lm.h5:/plat', s_lat%plat, find_h5dtype(s_lat%plat))
            call h5_write('h_lm.h5:/qlat', s_lat%qlat, find_h5dtype(s_lat%plat))
            call h5_write('h_lm.h5:/alat', s_lat%alat, find_h5dtype(s_lat%alat))
            do cix = 1, ncix
               lcix(cix) = s_dmft%l(iabs(s_dmft%icix(cix)))
            enddo
            call h5_write('h_lm.h5:/lcix', lcix, find_h5dtype(lcix))


            allocate(wtk(s_bz%nkabc(1)*s_bz%nkabc(2)*s_bz%nkabc(3)))
            call bzmshp('makeprojcrpa',1,s_bz%nkabc,s_bz%lshft,s_lat%plat,s_lat%symgr,
     .           s_lat%nsgrp,s_lat%nsgrp,.true.,0,0,qb,nkfbz,lstk,wtk,lstk,0)
            call h5_write('h_lm.h5:/lstk', lstk, find_h5dtype(lstk(1,1)))


         endif
         comm = mpi_comm_world
         call ms_h5 % init( [nevn,ldcix])
         call pfs_h5 % init( [nevn,ldcix,nsp,ncix,nkabc(1)*nkabc(2)*nkabc(3)])
         call pf % init('h_lm.h5',comm=comm)
         call pdset % init(pf, '/proj', dtype = find_h5dtype(srm1),  space = pfs_h5)
         call efs_h5 % init( [nevn,nsp,nkabc(1)*nkabc(2)*nkabc(3)])
         call edset % init(pf, '/ek', dtype = h5t_native_real8,  space = efs_h5)

      endif


C     Distribute k over processors
      nproc = mpipid(0)
      allocate (kpproc(0:nproc))
      if (nproc > 1) then
        call info0(30,1,0, ' ... Start MPI k-loop for projectors')
        if (nproc > nkp) call rxi('MPIK job cannot allocate more processors than nkp =',nkp)
        call dstrbp(nkp,nproc,1,kpproc(0))
      else
        kpproc(0:1) = [1,nkp+1]
      end if

C     Debugging
C      allocate(evlq(ndham,nsp,nkp))
C      allocate(zq(ndham*nspc,ndham*nspc,nsp,nkp))
C      call readevec(2,ndham,nsp,nspc,nkabc,nkp,lshft,nproc,evlq,zq)
C      call rx0('done')

!      print *, '!!'; nproc = 2


C --- Make projectors for each point in the full Brillouin zone ---
C     Pass 1 : assemble projector overlap.  Skipped if lolap is 2
C     Pass 2 : Normalize projectors and write to disk; make double counting
      do  ipass  = 1, 2

      if (ipass == 1 .and. lolap == 2) cycle

C ... Read file header to evec.  Read also checks match on nsp,nspc,ndham
      ifiz = fopna('evec',-1,4); rewind ifiz
      if (procid == master) then  ! No distinction between ndham and nlmto for now
        call iosigh(3,LW5,nsp,nspc,ndham,ndham,nk1,nk2,nk3,nkp,nqi,lshft(1),lshft(2),lshft(3),ifiz,eseavr)
      endif
      nkfbz = nk1*nk2*nk3       ! nkfbz=0 flags that points are not defined on a regular mesh
      if (nkfbz == 0 .and. lstar) call rx1('sudmft: expected regular k-mesh but found nkabc=%s,(%3i)',nkabc)

C ... Initial setup assuming no stars of k
      allocate(dmftuq(nevn,ldcix,nsp,ncix,1,1))
      call dpzero(dmftuq,2*size(dmftuq))
      ig = 1

      do  iq = kpproc(procid), kpproc(procid+1)-1

C       if (iq == 16) call snot

C       print 777, 'starting q loop',procid,iq
        do  isp = 1, nsp
          if (nkp /= nkfbz .and. lstar) call info0(PRTU,0,0,' ')
          read(ifiz) qp, ndimhx, nev
          ndimh = ndimhx / nspc
          call dpdump(evlk(1,isp),nev,ifiz)
          call dpdump(z,ndimhx**2*2,ifiz)

C         print 777, 'read evec q,evlk,z',procid,iq,0,qp,evlk(nlohi(1),1),z(1,1)

C     ... For each point in the star of k, make unnormalized projectors and add to the overlap
          iqs=0; istar=0; qpr=qp; zr => z
          do  while (.true.)    ! loop over star of qp
            if (lstar) then     ! Generate dmftu for entire star
              if (iqs==0 .and. isp==1) then ! First point: count number of points in star
                call iqstar(-1,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,nstarq,qpr,[0])
                deallocate(dmftuq)
                allocate(dmftuq(nevn,ldcix,nsp,ncix,nstarq,1))
                call dpzero(dmftuq,2*size(dmftuq))
              endif
              call iqstar(2,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,i123)
              if (iqs == 0) exit  ! No more qpr in star
              ig = s_bz%star(iqs+1)
              if (ig > 1) then
                zr => h
                allocate(rmat(nl**4*2))
                call rotwf(140,nl,nbas,nspc,s_lat%pos,s_ham%offH,s_ham%iprmb,s_lat%istab(1,ig),
     .            s_lat%symgr(1,ig),s_lat%ag(1,ig),qpr,rmat,0,ldim,ldim,nev,z,zr)
                deallocate(rmat)
              endif
            endif

C      ...  Unnormalized DMFT projectors
            istar = istar+1
            dmftu => dmftuq(:,:,:,:,istar,1)
            allocate(ausp(nlmax,ndhamx,nphimx,nsp,nbas)); call dpzero(ausp,2*size(ausp))
            call makusq(1,s_site,s_spec,s_lat,s_ham,nbas,nbas,0,nlmax,1,ndham,nphimx,ndimh,
     .        s_lat%napw,s_lat%igv2,nlohi(2),nsp,nspc,nsp,isp,1,qpr,zr,ppnl,ausp,xv)
C           In noncollinear case, isp=1 always => need internal jsp=1..2
C           ksp is the current spin index in both cases:
C           ksp = isp  in the collinear case
C               = jsp in the noncollinear case
C           jsp=1 for independent spins, and spin index when nspc=2
            do  jsp = 1, nspc
              ksp = min(max(jsp,isp),nsp)
              call makeproj(s_dmft,nbas,nsp,ksp,ndhamx,ldcix,nlohi(1),nlohi(2),nlmax,nphimx,ausp,
     .          ppnl,nppn,n0,iq,qpr,dmftu)
C             print 777, 'unnorm dmftu',procid,iq,istar,dmftu(1,1,ksp,1)
            enddo
            deallocate(ausp)

C       ... Collect dmftu for both spins before remaining block proceeds
            if (isp < nsp .and. .not. lstar .and. nspc == 1) exit
            if (isp < nsp .and. nspc == 1) cycle

C       ... Accumulate projector overlap
            if (s_dmft%knorm == 0 .and. lolap < 2) then
              facalpha = 1d0/nkfbz; facbeta = 1
            else
              facalpha = 1; facbeta = 0
            endif

            if (ipass == 1) then
               do  cix  = 1, ncix
                  icix = iabs(s_dmft%icix(cix))
                  nl1 = 2*s_dmft%l(icix)+1 ! 2*l+1
                  if (nl1 > ldcix) call rx('bug in lmlproj')
                  do  jsp = 1, nsp
C     call zprm('U',2,dmftu(1,1,jsp,cix),nevn,nevn,nl1)
                     call zgemm('C','N',nl1,nl1,nevn,facalpha,dmftu(1,1,jsp,cix),nevn,
     .                    dmftu(1,1,jsp,cix),nevn,facbeta,Olapm(1,1,jsp,cix),ldcix)
C     if (iq == 1) then
C     print *, iq, ig, jsp, sngl(qpr)
C     call zprm('O',2,Olapm(1,1,jsp,cix),ldcix,nl1,nl1)
C     endif
                  enddo         ! spin
               enddo            ! cix
            endif
            if (s_dmft%knorm == 1) then
               do  cix = 1, ncix
                  icix = iabs(s_dmft%icix(cix))
                  nl1 = 2*s_dmft%l(icix)+1 ! 2*l+1
                  do  jsp = 1, nsp
                     call zpotrf('U',nl1,Olapm(1,1,jsp,cix),ldcix,j) ! Cholesky decomposition
                     if (j /= 0) call rx('zpotrf info /= 0')
                     call ztrtri('U','N',nl1,Olapm(1,1,jsp,cix),ldcix,j) ! inverse of upper triangular
                     if (j /= 0) call rx('ztrtri info /= 0')
                     do  j = 1, nl1-1 ! Zero out lower triangle
                        forall (k = j+1:nl1) Olapm(k,j,jsp,cix) = 0
                     end do
C     Check : mc -f10f15.10 out.lsco -p -t -x -i ovl --
C     call zprm('L',2,Olapm(1,1,jsp,cix),ldcix,nl1,nl1)
                  enddo         ! spin
               enddo            ! cix
            endif

C     ...   Normalize projectors, make double counting and write to proj file
C           Writes to file 'proj' 3 records: header, evlk, dmftUi for each each cix, star in iq, irr iq
C           Let Olap = LL+ and Projector = U U+.  Olap stores L^-1
C           Write renormalized projector U~ = L^-1 U  So that U~ U~+ = L^-1 U (L^-1 U)+
            if (ipass == 1) allocate(projemb(ldcix,ldcix)) ! debugging check to confirm that U U+ is orthonormal
            if (ipass == 2 .or. s_dmft%knorm == 1) then
              ifi = fopna('proj',-1,4)
              if (procid == 0 .and. nproc > 1) then  ! Write to proj.ext_0 for head node
                 jfi = ifi                 ! file id for proj.ext, written after all processes write to their own files
                cix = fopnx(fn,192,0,ifi) ! Full file name project.ext
                ifi = fopng(trim(fn)//'_0',-1,4) ! Write to proj.ext_0
              endif
              do  cix = 1, ncix
                icix = iabs(s_dmft%icix(cix))
                nl1 = 2*s_dmft%l(icix)+1 ! 2*l+1
                allocate(dmftUi(nevn,nl1,nsp))  ! Allocate here to ensure contiguous spin 1,2 blocks
                do  jsp = 1, nsp
                  call zgemm('N','N',nevn,nl1,nl1,(1d0,0d0),dmftu(1,1,jsp,cix),nevn,
     .              Olapm(1,1,jsp,cix),ldcix,(0d0,0d0),dmftUi(1,1,jsp),nevn)

C                 Test: assemble overlap of renormalized U (should be unit matrix)
                  if (iprint() >= PRTG) then
                    if (ipass == 2) projemb => olapmk(:,:,jsp,cix)
                    call zgemm('C','N',nl1,nl1,nevn,facalpha,dmftUi(1,1,jsp),nevn,
     .                dmftUi(1,1,jsp),nevn,facbeta,projemb,ldcix)
                    if (ipass == 1) then
                      do  i = 1, nl1
                        projemb(i,i) = projemb(i,i) - 1
                      enddo
                      call info5(10,0,0,' deviation from unity, renormalized overlap iq=%i cix=%i : %;3g',
     .                  iq,cix,dlength(2*size(projemb),projemb,1),0,0)
C                    else
C                      if (jsp == 1) print *, iq,olapmk(1,1,1,1)
                    endif
                  endif
                enddo           ! jsp

C               print 777, 'writing to proj',procid,iq,ig,dmftui(3,1:3,1)
                write(ifi) cix,nevn,nl1*nsp,qpr,iq,ig,qp ! Write to proj file
                call dpdump(evlk,ndham*nsp,-ifi)
                call dpdump(dmftUi,2*nevn*nl1*nsp,-ifi)
                if ((save_h5) .and. (ipass==2) .and. (lstar)) then
                   comm = mpi_comm_world
!     iqfbz=(i123(3)-1)*nkabc(1)*nkabc(2)+nkabc(1)*(i123(2)-1)+i123(1)
                   iqfbz=(i123(1)-1)*nkabc(3)*nkabc(2)+nkabc(3)*(i123(2)-1)+i123(3) ! FOR TRIQS OR GW code
!     Inverse nsp and ncix order to make the compression ncix*ldcix easier. Most of quantity are diagonal in spin. in cix, not always.
                   call ms_h5 % resize( [nevn,ldcix,nsp])
                   call pfs_h5 % resize ([nevn,ldcix,nsp,ncix,nkabc(1)*nkabc(2)*nkabc(3)])
                   call pfs_h5 % select_hyperslab(start=[0,0,0,cix-1,iqfbz-1], count=[nevn,ldcix,nsp,1,1])
                   call pdset % write(dmftUi(1,1,1), mem_space = ms_h5, file_space =  pfs_h5)
                   if(cix == 1) then
                      call ms_h5 % resize( [nevn,nsp])
                      call efs_h5 % resize ([nevn,nsp,nkabc(1)*nkabc(2)*nkabc(3)])
                      call efs_h5 % select_hyperslab(start=[0,0,iqfbz-1], count=[nevn,nsp,1])
                      evlkh5(:,:)=  evlk(s_dmft%nlohi(1):s_dmft%nlohi(2),:) - ef0
                      call dscal(nsp*nevn,ry2ev,evlkh5,1)
                      call edset % write(evlkh5(1,1), mem_space = ms_h5, file_space =  efs_h5)
                   endif

                endif
                deallocate(dmftUi)
             enddo              ! cix

C         ... Double counting from QSGW static self-energy
              if (lrsig /= 0 .and. s_dmft%dcmode /= 1) then
                call rx('update QSGW double counting')
C                if (nsp == 2) call rx('update for spin polarized case')
C                allocate(sigdcmk(orbmax,orbmax,nicix,nicix))
C                call hambls(5000+10*lrsig,s_site,s_spec,s_lat,s_ham,nbas,isp,1,qp,k1,k2,k3,s_ham%qsig,nqsig,
C     .            s_pot%smpot,s_pot%smvcnst,s_pot%lcplxp,lso,0d0,ndimh,h,sigqp(1,1,isp),xv,i,nev)
C                call makedc(s_dmft,dmftu,nevn,nlohi,orbmax,isp,nsp,ndsigi(1),nicix,sigqp,sigdcmk)
C!               #Mark the result of last routine is k-res DC
C                if (.not. s_dmft%ldadc) then
C                  call rx('patch this dc sigma')
CC                 sigdc = sigdc + sigdcmk/nkfbz
C                endif
C                deallocate(sigdcmk)
              endif
!             #Mark last term is eq(3.13), provided sigqp is sigma_qsgw or e_hartree, it's up to you

            endif
            if (ipass == 1) deallocate(projemb)  ! Used for printout
            if (.not. lstar) exit
          enddo                 ! Loop over star of q
          if (nspc == 2) exit   ! No second spin in noncollinear case
       enddo                    ! Loop over spins
      enddo                     ! Loop over irreducible qp
      deallocate(dmftuq)

C     Debugging : print overlap matrix
C      do  cix = 1, ncix
C        icix = iabs(s_dmft%icix(cix))
C        nl1 = 2*s_dmft%l(icix)+1
C        do  jsp = 1, nsp
C          if (ipass == 1) call zprm('Olap',2,Olapm(1,1,jsp,cix),ldcix,nl1,nl1)
C          if (ipass == 2) call zprm('Olap~ ',2,Olapmk(1,1,jsp,cix),ldcix,nl1,nl1) ! Should be identity
C        enddo
C      enddo

      if (s_dmft%knorm == 1) exit    ! No 2nd pass if normalization is k-dependent

      if (iprint() >= PRTG .and. ipass == 2) then
        do  cix = 1, ncix
          icix = iabs(s_dmft%icix(cix))
          nl1 = 2*s_dmft%l(icix)+1
          do  jsp = 1, nsp
            projemb => olapmk(:,:,jsp,cix)
            do  i = 1, nl1
              projemb(i,i) = projemb(i,i) - 1
            enddo
            call info2(10,0,0,' deviation from unity, renormalized overlap cix=%i : %;3g',
     .        cix,dlength(2*size(projemb),projemb,1))
          enddo
        enddo
      endif
C ... Cholesky decompose overlap: Rewrite Olapm with (L+)^-1  where L L+ = overlap
      if (ipass == 1) then
        call mpibc2(olapm,2*size(olapm),4,3,.false.,'','')
        do  cix = 1, ncix
          icix = iabs(s_dmft%icix(cix))
          nl1 = 2*s_dmft%l(icix)+1
          do  jsp = 1, nsp
            call zpotrf('U',nl1,Olapm(1,1,jsp,cix),ldcix,j) ! Cholesky decomposition
            if (j /= 0) call rx('zpotrf info /= 0')
            call ztrtri('U','N',nl1,Olapm(1,1,jsp,cix),ldcix,j)
            if (j /= 0) call rx('ztrtri info /= 0')
            do  j = 1, nl1-1    ! Zero out lower triangle
              do  k = j+1, nl1
                Olapm(k,j,jsp,cix) = 0
              end do
            end do
C           Check : mc -f10f15.10 out.lsco -p -t -x -i ovl --
C           call zprm('Lx',2,Olapm(1,1,jsp,cix),ldcix,nl1,nl1)
          enddo                 ! spin
        enddo                   ! cix
      endif

      enddo                     ! Passes 1 and 2
      if (save_h5) then
         call pdset%close()
         call edset%close()
         call pf%close()
      endif

      call fclr(' ',ifi)
      if (mpipid(2) /= 0) call rx('mpi barrier returned with error')
      deallocate(evlk,z,h,Olapmk)
      if (lolap == 0) deallocate(Olapm)

C --- Multithreaded case: merge projectors into a single file ---
      if (nproc > 1 .and. procid == 0) then
        allocate(dmftUi(nevn,64,2)) ! Allocate array larger than anything actually needed
        allocate(evlk(ndham,nsp))

C   ... Loop over k-points serially, reading iq-dependent proj file
        oldrank = -1
        do  iq = 1, nkp

          call hunti(kpproc,nproc,iq+1,[-1],rank)  ! returns rank = procid+1 for this iq
          rank = rank-1
          if (oldrank /= rank) then  ! New file
            call awrit1('%x_%i',extmpi,len(extmpi),0,rank)
            if (oldrank /= -1) call fclr(' ',ifi)
            ifi = fopng(trim(fn)//trim(extmpi),-1,4)
C           print 777, 'switching to new file: '//trim(fn)//trim(extmpi),rank,iq,ifi
            oldrank = rank
          endif

          iqs=0; istar=0
          do  while (.true.)    ! loop over star of qp
            if (lstar) then  ! Generate dmftu for entire star
              if (iqs==0) then ! First point: count number of points in star
                call iqstar(-1,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,nstarq,qpr,[0])
C               istarq(iq) = nstarq
              endif
              call iqstar(1,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,[0])
              if (iqs == 0) exit  ! No more qpr in star
            else
              istar = 0
            endif

            istar = istar+1

            do  cix = 1, ncix
              icix = iabs(s_dmft%icix(cix))
              nl1 = 2*s_dmft%l(icix)+1

              read(ifi) i,nevn,j,qpr,k,ig,qp ! read from file proj.ext_n
              if (i/=cix .or. j/=nl1*nsp .or. k/=iq) then
                print 777, 'cix,dim,iq ', i,j,k
                print 777, 'file values',cix,nl1*nsp,iq,qp
  777           format(a,3i5,10f14.7)
                call rx('bug in lmlproj')
              endif
              call dpdump(evlk,ndham*nsp,ifi)
              call dpdump(dmftUi,2*nevn*j,ifi)
C              print 777, 'reading from proj',rank,k,ig,dmftui(3,1:3,1)

              write(jfi) i,nevn,j,qpr,k,ig,qp ! write to file proj.ext
              call dpdump(evlk,ndham*nsp,-jfi)
              call dpdump(dmftUi,2*nevn*j,-jfi)

C            chk(iq,2) = chk(iq,2) + sum(dmftUi)

            enddo               ! cix
            if (.not. lstar) exit
          enddo                 ! Loop over star of q

C        print 777, 'sanity check',procid,iq,ig,chk(iq,1),chk(iq,2)-chk(iq,1)

        enddo                   ! Loop over iq
        call fclr(' ',ifi)      ! close file generated by last procid
        deallocate(dmftUi)

      endif ! MPI to serial conversion of proj file

      call fclr(' ',ifiz)        ! close evec file

      end
C      subroutine snot
C      return
C      end
