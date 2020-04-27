      subroutine bndtb(s_ctrl,s_bz,s_lat,s_spec,tbc,nbas,nl,nspc,nsp,&
         &  nsp1,lmx,idxdn,nclass,ips,ipc,nrclas,pv,force,zval,mull,   &
         &  wtbzm,nsite,iax, &
         &  npr,xyzfrz,hrs,h0,dh,dhcf,vso,hso,srs,ds,dscf,pot0,    &
         &  efermi,sumev,entrpy,f,thrpv,esite,rho,rhoc,rhon,zos, &
         &  drhos,ldim, idxsh)
!C- k-integration of tight-binding bands
!C ----------------------------------------------------------------------
!Ci Inputs
!Ci   nbas, no. of atoms in the basis; nl, number of l's;
!Ci   nspc (2 for coupled spins) nsp (2 for spin polarised TB)
!Ci   nps1 is for dimensioning arrays (equivalent to nspx in LMASA)
!Ci   lmx, max. l for classes;
!Ci   if idxdn(l,i) > 1 then this L for ith class is excluded from basis;
!Ci   nclass, no. of classes - atoms in the same class are
!Ci   symmetry-related; the jth atom belongs to class ipc(j);
!Ci   nrclas(i), no. of atoms in the ith class;  force: TRUE for
!Ci   forces; nbmax, max no. of bands;
!Ci   zval, total number of valence electrons
!Ci   mull, Mulliken DOS switch, see mkwttb;
!Ci   npln,nwmx,nqmx,nw,ew,nqpl,vxpl,vypl,xpl,ypl,zpl,wtbzm:
!Ci   parameters for BZ maps (integrated DOS vs k), see getbzp and bzmio
!Ci   nsite:   total number of neighbors in all clusters;
!Ci   iax, neighbor lists; npr, see tbham;
!Ci   hrs,h0,dh: real-space hamiltonian, and dh(x,y,z,r) its derivatives
!Ci   vso:    table of spin-orbit parameters and hso is
!Ci   hso:    the spin-orbit hamiltonian
!Ci   srs,ds: real-space overlap matrix, and ds(x,y,z,r) its derivatives
!Ci   dhcf,dscf: the crystal field ham. and overlap derivatives
!Ci   pot0: monopole potential at site R (eq. 7.81, Finnis)
!Co Outputs
!Co   iwk, number of orbitals for each atom; eband, energy bands; efermi;
!Co   sumev, sum occupied e'vals; f(3,nbas), force on each atom from the
!Co   entrpy, "entropy term" actually TS
!Co   band structure; thrpv, 3PV from the band structure; esite, band
!Co   energies for each atom;
!Co   rho, s,p, and d Mulliken charges;
!Co   rholm, {lm} decomposed Mulliken charges;
!Co   If STONER=T zos returns the local number of states and index points
!Co   to the d-electron components.
!Co   rho, rhoc and rhon are, respectively, local charges by l,
!Co   eigenvector products for TB-L and TB+U (see tbfrce and tbesel)
!Co   drhos: (c*_RL dS_RLR'L'/dR c_R'L' + cc) s-c TB+ovlp
!Co          built in tbfrce
!Co   ldim is passed back up for dimensioning drhos in tbesel
!Cf  Files:
!Cf    BAND file has structure (compatibility with ASA)
!Cf    1.   nl 1 1 nkp ldim*nspc nfstg
!Cf    For each qpt, the following records:
!Cf    2.   nchan  nev (if nfstg nonzero)
!Cf         eband      (if 1s   digit of nfstg nonzero)
!Cf         doswt      (if 10s  digit of nfstg 2)
!Cf    3.   efermi, 0d0
!Cm MPI
!Cm   Parallelise over the k-loop. Each process computes contributions
!Cm   to efermi, sumev, f, thrpv, esite, rho. eband is assembled at
!Cm   different k by different processors in slices of a 3D communicator.
!Cm   Each slice is a 2D communicator and preferably deals with single k-point.
!Cm
!Cu Updates
!Cu        2013 (DP) Major rewrite to handle distributed matrices and process topologies
!Cu   10 Nov 11 Begin migration to f90 structures
!Cu   04 Jun 08 (ATP) Modificiations for gamma-only ; faster execution
!Cu   10 Apr 02 Bug fix --- nfstg now set correctly
!Cu   15 Feb 02 (ATP) Added MPI parallelization
!Cu   03 Feb 01 use ldos instead of lbzio
!C ----------------------------------------------------------------------

      use mpi
      use tbprl
      use mod_ctx, only : cart_ctx_t
      use mod_tbdos
      use structures

      implicit none
      integer, parameter :: dp = 8

!C Passed parameters
      type(str_ctrl) ::  s_ctrl
      type(str_bz) ::    s_bz
      type(str_lat) ::   s_lat
      type(str_spec) ::  s_spec(*)
      type(tbc_t), intent(in), target :: tbc

      procedure(real(dp)) :: ddot

      integer nbas,nl,nsp,nspc,nsp1,nclass,nbmax,mull,npln,nqmx,nsite,ldim, idxsh(*)
      integer lmx(*),idxdn(0:nl-1,*),ips(*),ipc(nbas),nrclas(*), iax(*),npr(*),iwk(nbas)
      real(dp) :: zval,wtbzm,efermi,sumev,entrpy,thrpv,egap
      real(dp) :: hrs(*),h0(*),dh(*),vso(*),hso(*),      &
         &  srs(*),ds(*),f(3,nbas),esite(nbas,nsp1),    &
         &  rho(nl,nsp1,nbas),rholm(nl*nl,nsp1,nbas),rhoc(nl*nl,nl*nl,nbas),   &
         &  rhon(nl*nl,nl*nl,nbas,nsp1),zos(*),dhcf(*),dscf(*),drhos(3,nbas,nbas), pot0(nbas)
      logical pv,force,xyzfrz(3)
!C Local variables
      integer parg,ip,it,ldos,lncol,lstnr(3),ltb,getef,nste(2)
      integer i,j,k,l,ikp,nfilet,nfileb,isite,        &
         &  lihdim,ipr,iprx,nev,fopn,fopnT,iprint,icheck,ibas,indx,nchan,  &
         &  ic,ifi,nevhi,nevlo,nat0,nlm,mull0,mull1,fopno,i1mach,lm,ll,    &
         &  imfil,ndos0,ndos,nfstg,isp,ispx,ldimx,nbmaxx,ik0,ikstep, nevmx0,iomoms
      integer lmet,nkabc(3),n1,n2,n3,nkp,ntet,mpsord,norder,nevmx
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      integer nsgrp,oistab,oag,osymgr, lgunit, mode
      real(dp), parameter :: enull=99999d0
      real(dp), allocatable :: eband(:,:,:), wtkb(:,:,:)
      real(dp), allocatable :: zll(:,:,:,:,:)
      integer :: inulla(1)

!C --- for reconstructing RS hamiltonian
      integer :: metal, tetra
      real(dp) :: dum,del,evlo,qp(3),dval,dosw(2),plat(3,3),alat, srnge,swidth,efmax,dosef(2),qval(2),cv, &
         & gd(5,5,s_lat % nsgrp), tmpr
      logical :: leig,doswt,fsym,charge,rl,trh,lso,donly,lmol, &
         & lgamma,lgamsp,lov,bzmp,stoner,cryf,ocryf,rstrt,swnos,TBU,rn,bittst,cmdopt,excite,sc
      character*80 outs
      character(len=45), parameter :: multyp(0:3,0:1) = reshape([ &
            &  ' Partial DOS, resolved by class and l        ',   &
            &  ' Bond charge, resolved by class pairs        ',   &
            &  ' Bond charge, resolved by {class, l} pairs   ',   &
            &  ' Bond charge, resolved by {class, l, m} pairs',   &
            &  ' Partial DOS, resolved by atom, l, and m     ',   &
            &  ' Bond charge, resolved by atom pairs         ',   &
            &  ' Bond charge, resolved by {atom, l} pairs    ',   &
            &  ' Bond charge, resolved by {atom, l, m} pairs '], [4,2]  )

      integer err, kbid, nl2, nl4, typesz, kpstep
      logical :: mpiav
!       logical, parameter :: aprox = .true.
      logical, parameter :: aprox = .false.

      type(cart_ctx_t), pointer :: c3d, c2d, c1d
      real(8) :: tube_r8(8)
!       integer, allocatable :: wcounts(:), wdispls(:)
      integer :: eval_mt, rho_mt, rhoc_mt, rholm_mt, e_mt, f_mt, drhos_mt, wts_mt
      integer :: d2cnt, ia
      integer(8) :: memstr

      call tcn('bndtb')

      inulla = 0

      c3d => tbc % c3d
      c2d => tbc % c2d
      c1d => tbc % c1d

      kbid = c1d % id

      mpiav = c3d % sz > 1


!       call mpiprn(shortname,length)
!       mlog = cmdopt('--mlog',6,0,outs)



!       call mpi_barrier(mpi_comm_world, ierr)
!       print "(/,'p bndtb: ', i0/)", pid



!C --- Setup ---
      ltb = s_ctrl % ltb
      ldos = s_ctrl % ldos
      lncol = s_ctrl % lncol
      lstnr = s_ctrl % lstonr
!C     call upack('bz nkabc nkp ntet oidtet lmet',sbz,nkabc,nkp,ntet,
!C    .  oidtet,lmet)
      nkabc = s_bz % nkabc
      nkp = s_bz % nkp
      ntet = s_bz % ntet
!       oidtet = s_bz % oidtet
      lmet = s_bz % lmet
!C     call upack('bz ndos dosw oqp owtkp',sbz,ndos0,dosw,oqp,owtkp,0)
      ndos0 = s_bz % ndos
      dosw = s_bz % dosw
!       oqp = s_bz % oqp
!C     call upack('bz n range w nevmx efmax',sbz,mpsord,srnge,
!C    .  swidth,nevmx,efmax)
      mpsord = s_bz % n
      srnge = s_bz % range
      swidth = s_bz % w
      nevmx = s_bz % nevmx
      efmax = s_bz % efmax
      mpsord = isign(1,mpsord) * mod(iabs(mpsord),100)
!C     call upack('lat plat alat',slat,plat,alat,0,0,0)
      plat = s_lat % plat
      alat = s_lat % alat
!C     call upack('lat nsgrp oistab oag osymgr',slat,nsgrp,oistab,oag,
!C    .  osymgr,0)
      nsgrp = s_lat % nsgrp
!       oistab = s_lat % oistab
!       oag = s_lat % oag
!       osymgr = s_lat % osymgr

      tetra = 0; if (ntet > 0) tetra = 1
      metal = 0; if (lmet /= 0) metal = 1
      stoner = lstnr(1) /= 0
!C      doswt  = bittst(ldos,2) .or. bittst(ldos,4)
      doswt  = (mull >= 0)
      trh    = bittst(ltb,2**10)
      charge = bittst(ltb,2**11) .or. bittst(ltb,2**13) .or. bittst(ltb,2**15)
      rl     = bittst(ltb,2**15)
      lso    = bittst(lncol,4)
      donly  = bittst(ltb,2**6)
      cryf   = bittst(ltb,2)
      lov    = bittst(ltb,1)
      ocryf  = bittst(ltb,4)
      bzmp   = bittst(ldos,8)
!C     T if just gamma point
      lgamma = bittst(ltb,2**17)
      lmol   = bittst(ltb,2**18)
      if (lmol) lgamma = .true.

!C spin polarised gamma only, or molecule:
      lgamsp = lgamma .and. nsp == 2

      if (lgamsp) then
        call rxx(tetra /= 0,'BNDTB: for spin pol at gamma use TETRA=F')
      endif

      TBU = bittst(ltb,2**13)
      rn = TBU
      sc = rl .or. rn
      trh = trh .or. rl .or. rn
      if (stoner) doswt = .true.
      rstrt = cmdopt('-cont',5,0,outs)
      fsym = nsgrp > 1 .and. (.not. lmol)
      leig = force .or. charge .or. trh .or. pv
      swnos = ndos0 < 0
      ndos = iabs(ndos0)
      nfstg = 1
      if (doswt) nfstg = 11
      call rxx(lso .and. nspc /= 2,'BNDTB: must set NSPIN=2 for SO=T')
!C --- Set up permutation matrix for Bloch ---
      lihdim = nbas * nl**2
      ldimx = ldim*nspc
!       nbmaxx = nbmax*nspc
!       lgamma = .false.
!       print *, 'reminder to fix lgamma in bndtb and secmtb'

      nl2 = nl*nl
      nl4 = nl2*nl2

      typesz = 1
      if (.not. lgamma) typesz = 2

      if (rl)  rhoc = 0.0_dp
      if (lov .and. force .and. sc) drhos = 0.0_dp
      if (rn)  rhon = 0.0_dp


!       print *, 'c1d:', c1d%id, c1d%crd, c1d%lrt
!       print *, 'c2d:', c2d%id, c2d%crd, c2d%lrt
!       print *, 'c3d:', c3d%id, c3d%crd, c3d%lrt
!
!       call mpi_barrier(mpi_comm_world, err)
!       stop

!C --- Get number of orbitals on each atom ---
      if (leig) then
        iwk = 0
        icheck = 0
        do ibas = 1, nbas
          do l = 0, nl-1
            indx = idxdn(l,ipc(ibas))
            if (indx > 1) cycle
            icheck = icheck + 2*l + 1
            iwk(ibas) = iwk(ibas) + 2*l + 1
          end do
        end do
        call rxx(icheck /= ldim,'TBZINT: icheck /= ldim')
      endif

!C --- Determine verbosity for secmtb ---
      ipr=0
      if (iprint() >= 30) ipr=1
      if (iprint() >= 35) ipr=2
      if (iprint() > 40) ipr=3



!       nevmx0 = (int(zval) + 1)/2
!       if (nevmx == 0) then
!         if (metal /= 0) then
!           nevmx = min(max(nevmx0+nevmx0/2,nl*nl),ldim)
!         else
!           nevmx = min(nevmx0+1,ldim) ! to have one unoccup. state for gap estimation
!         endif
!         nevmx = nspc*min(nevmx,nbmax)
!       endif
!       nevmx0 = nspc*min(nevmx0,nbmax)
!       if (iprint() > 10 .and. nevmx < nevmx0 .and. leig) &
!          & call awrit2(' ***WARNING*** nevmx=%i < %i, expect errors in integrals',' ',80,i1mach(2),nevmx,nevmx0)

      nevmx = ldim
      iprx = ipr
      if (mpiav) iprx = 0

      thrpv = 0d0

      if (force) f = 0.0_dp
      if (trh) esite = 0.0_dp
      if (trh .or. charge) rho = 0.0_dp
      if (trh .or. charge) rholm = 0.0_dp




      if (.not. aprox) then

      memstr = typesz * tbc%loclsz(1) * tbc%loclsz(2) * nsp1 * (tbc % kmap(kbid+1) - tbc % kmap(kbid)) * 8
      memstr = memstr/2**20
      if (memstr > 24) call info2(10,1,1,' BNDTB: allocating  %i MiB memory for evecs at root',memstr,0)
      allocate(zll(typesz, tbc%loclsz(1), tbc%loclsz(2), nsp1, tbc % kmap(kbid)+1 : tbc % kmap(kbid+1)))
!          allocate(zll(typesz, tbc % globsz(1), tbc % globsz(2), nsp1, tbc % kmap(kbid)+1 : tbc % kmap(kbid+1)))
!       print *, 'g,l sizes',   tbc % globsz, tbc % loclsz

      zll = 0.0_dp
      qp = 0.0_dp

!C...   initialize eigenvalues to enull (see, eg, subs/efrng2.f)
      if (c3d%lrt) then
         allocate(eband(nevmx,nsp1, nkp))
      else
         allocate(eband(nevmx,nsp1, tbc % kmap(kbid)+1 : tbc % kmap(kbid+1)))
      endif

      eband = enull


      do ikp = tbc % kmap(kbid)+1, tbc % kmap(kbid+1)
         if (.not. lgamma) call dpscop(s_bz%qp,qp,3,3*ikp-2,1,1d0)
         do isp = 1, nsp
            call secmtb(s_ctrl,tbc,plat,nbas,nl,nspc,nsp,isp,lmx,  &
               & ipc,idxsh,ldimx,nevmx,efmax,ikp,nkp,qp,nsite,     &
               & iax,npr,hrs,vso,hso,srs,pot0,rl,iprx,leig,nev,    &
               & zll(1,1,1,isp,ikp),eband(1,isp,ikp))
         end do
      end do


      if (c2d%lrt .and. tbc % c1d % sz > 1) then
         call mpi_type_contiguous(nevmx*nsp1, mpi_real8, eval_mt, err)
         call mpi_type_commit(eval_mt, err)
         if (c1d%lrt) then
            call mpi_gatherv(mpi_in_place, 0, eval_mt, &
                           & eband, tbc % kcount, tbc % kmap, eval_mt, c1d % rt, c1d % comm, err)
         else
            call mpi_gatherv(eband(1,1,tbc%kmap(kbid)+1), tbc % kcount(kbid), eval_mt, &
                           & eband, inulla, inulla, eval_mt, c1d % rt, c1d % comm, err)
         end if
         call mpi_type_free(eval_mt, err)
      end if

      if (c3d%lrt) then
         allocate(wtkb(nevmx,nsp1,nkp))
      else
         allocate(wtkb(nevmx,nsp1, tbc % kmap(kbid)+1 : tbc % kmap(kbid+1)))
      endif

      wtkb = 0.0_8

      if (c3d%lrt) then
!C --- BZ integration for fermi level, band sum and qp weights ---
!C      if (.not. donly .and. (.not. lgamma .or. lgamsp)) then
         norder = mpsord
!C --- enforce metal sampling even if filled bands encountered ---
!C        norder = sign(1,mpsord)*(abs(mpsord) + 100)
!C        call rxx(norder < 0,
!C     .    'F-D not implemented for spin pol. Restart with N>=0')

         call bzwts(nevmx,nevmx,nsp1,nspc,n1,n2,n3,nkp,ntet,s_bz%idtet,  &
            & zval,metal,tetra,norder,ndos,swidth,srnge,s_bz%wtkp,eband,   &
            & efmax,efermi,sumev,wtkb,dosef,qval,entrpy,egap)


!       if (cmdopt('--wvecs',7,0,outs)) then
!          call pshprt(1)
!          call splwts(nkp,ldim,nevmx,nsp1,w(owtkp),eband,norder,      &
!                   &  swidth,efermi,metal,sumev,wtkb,qval,entrpy,dosef,cv)
!          call popprt
!          call iowts(nsp1,ldim,nevmx,eband,wtkb,efermi,qval)
!       endif
         tube_r8 = [efermi,sumev,dosef,qval,entrpy,egap]
      end if

      call mpi_bcast(tube_r8, 6, mpi_real8, c3d % rt, c3d % comm, err )

      if (.not. c3d%lrt) then
         efermi = tube_r8(1)
         sumev  = tube_r8(2)
         dosef  = tube_r8(3)
         qval   = tube_r8(4)
         entrpy = tube_r8(5)
         egap   = tube_r8(6)
      end if

      ! this needs to be in the 'exact' (not approx) branch and before zll is overwritten in tbfrce!
      if (cmdopt('--dos',5,0,outs)) call plotdos(tbc, nevmx, nsp1, s_bz%wtkp, eband, zll, ldim, efermi, &
                    lov, lgamma, nl, nbas, nsite, lmx, ipc, idxsh, iax, npr, srs, plat, s_bz%qp, vso, hso)

      if (metal /= 0 .and. c3d%sz > 1) then

         call mpi_type_contiguous(nevmx*nsp1, mpi_real8, wts_mt, err)
         call mpi_type_commit(wts_mt, err)
         if (c2d%lrt  .and. tbc % c1d % sz > 1) then
            if (c1d%lrt) then
               call mpi_scatterv(wtkb, tbc%kcount, tbc % kmap, wts_mt, &
                     & mpi_in_place, 0, wts_mt, c1d%rt, c1d%comm, err)
            else
               call mpi_scatterv(wtkb, inulla, inulla, wts_mt, wtkb, tbc%kcount(c1d%id), wts_mt, c1d%rt, c1d%comm, err)
            end if
         end if
         call mpi_bcast(wtkb(1,1,tbc%kmap(kbid)+1), tbc%kcount(c1d%id), wts_mt, c2d%rt, c2d%comm, err)
         call mpi_type_free(wts_mt, err)


!          call mpi_type_contiguous(nevmx*nsp1, mpi_real8, wts_mt, err)
!          call mpi_type_commit(wts_mt, err)
!          if (c3d%lrt) then
!             allocate(wcounts(0:c3d%sz-1), wdispls(0:c3d%sz-1))
!          else
!             allocate(wcounts(c3d%id:c3d%id), wdispls(c3d%id:c3d%id))
!          endif
!          wcounts(c3d%id) = tbc % kmap(c1d%id+1) - tbc % kmap(c1d%id)
!          wdispls(c3d%id) = tbc % kmap(c1d%id)
!
!          if (c3d%lrt) then
!             call mpi_gather(mpi_in_place, 1, mpi_integer, wcounts, 1, mpi_integer, c3d%rt, c3d%comm, err)
!             call mpi_gather(mpi_in_place, 1, mpi_integer, wdispls, 1, mpi_integer, c3d%rt, c3d%comm, err)
!             call mpi_scatterv(wtkb, wcounts, wdispls, wts_mt, mpi_in_place, 0, wts_mt, c3d%rt, c3d%comm, err)
!          else
!             call mpi_gather(wcounts, 1, mpi_integer, 0, 0, mpi_integer, c3d%rt, c3d%comm, err)
!             call mpi_gather(wdispls, 1, mpi_integer, 0, 0, mpi_integer, c3d%rt, c3d%comm, err)
!             call mpi_scatterv(wtkb, 0, 0, wts_mt, wtkb, wcounts(c3d%id), wts_mt, c3d%rt, c3d%comm, err)
!          endif
!          deallocate(wcounts, wdispls)
!          call mpi_type_free(wts_mt, err)
      end if
      end if


      if (iprint() >= 30) then
         if (lgamma) then
            print *,' SECMTB: hamiltonian real '
         endif
         j = min(9,ldim)
         if (iprint() >= 35) j = nevmx
         ikp = 1
         kpstep = 100
         if (nkp <= 100) kpstep = 10
         if (iprint() >= 40) kpstep = 1
         do while (ikp <= nkp)
!          do ikp = 1, nkp, 99
            call dpscop(s_bz%qp,qp,3,3*ikp-2,1,1d0)
            do isp = 1, nsp
               write(lgunit(1),'()')
               call awrit3(' SECMTB:  kpt %i of %i, k=%3:2,5;5d',' ',80,lgunit(1),ikp,nkp,qp)
               write(lgunit(1),'(255(9f8.4:/))') (eband(i,isp,ikp), i=1,j)
               if (ipr >= 2) call awrit5(' nev, nevmx, nevec, ldim=  %i  %i  %i  %i  efmax= %1;5d', &
                                                      & ' ',80,i1mach(2),nevmx,nevmx,nevmx,ldim,efmax)
            end do
            if (ikp == 1 .and. kpstep /= 1) ikp = 0
            ikp = ikp + kpstep
         end do
      endif

      mode = 0
      if (metal /= 0) mode = 1
      do ikp = tbc % kmap(kbid) + 1, tbc % kmap(kbid+1)
         if (.not. lgamma) call dpscop(s_bz%qp,qp,3,3*ikp-2,1,1d0)
         do isp = 1, nsp
            call tbfrce(tbc, mode,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp,   &
                     & lmx,ipc,idxdn,idxsh,iwk,ldimx,nevmx,zval,qp,   &
                     & abs(s_bz%wtkp(ikp))/nsp,wtkb(1,isp,ikp),mpsord,           &
                     & swidth,metal /= 0,efermi,nsite,iax,npr,xyzfrz,              &
                     & eband(1,isp,ikp),zll(1,1,1,isp,ikp),tbc%loclsz(1),typesz, &
                     & hrs,h0,dh,dhcf, vso,hso,srs,ds,dscf,charge,rl,rn,trh,pv,force,lso, &
                     & lov,cryf,ocryf, sumev,entrpy,f,thrpv,esite(1,isp), &
                     & rho,rholm,rhoc,rhon(1,1,1,isp),drhos)
         end do
      end do

      deallocate(zll, eband, wtkb)

!     sync: f drhos esite rho rholm rhoc thrpv


      if (c3d%sz > 1) then
         call tcn('netw')

!   First allreduce the local pieces across the kblock using the 1D communicators,
!   then allgather the pieces for different atoms within each kblock plane using the 2D communicators.
!   This shall be complete before the symmetrisation.
!       Beware mpi_sum is not defined for non-intrinsic(derived) types.

        if (c1d % sz > 1) then

         d2cnt = tbc%d2count(c2d%id)
         ia = tbc%d2amap(c2d%id)+1
         if (charge)call mpi_allreduce(mpi_in_place, rho(1,1,ia)   , nl*nsp1 *d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (charge)call mpi_allreduce(mpi_in_place, rholm(1,1,ia) , nl2*nsp1*d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (rl)    call mpi_allreduce(mpi_in_place, rhoc(1,1,ia)  , nl4     *d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (rn)    call mpi_allreduce(mpi_in_place, rhon(1,1,ia,1), nl4     *d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (rn .and. nsp==2) &
                    call mpi_allreduce(mpi_in_place, rhon(1,1,ia,2), nl4     *d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (trh)   call mpi_allreduce(mpi_in_place, esite(ia,1)   ,          d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (trh .and. nsp==2) &
                    call mpi_allreduce(mpi_in_place, esite(ia,2)   ,          d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (force) call mpi_allreduce(mpi_in_place, f(1,ia)       , 3       *d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
         if (force .and. lov .and. sc) &
                  & call mpi_allreduce(mpi_in_place, drhos(1,1,ia), 3*nbas   *d2cnt, mpi_real8, mpi_sum, c1d % comm, err)
        end if
   ! Utilising the mpi derived types we can use the same recvcnts and displs for different kinds of data in the allgather calls.
        if (c2d % sz > 1) then
         if (charge) call mpi_type_contiguous(nsp1*nl , mpi_real8, rho_mt  , err)
         if (charge) call mpi_type_contiguous(nsp1*nl2, mpi_real8, rholm_mt, err)
         if (rl .or. rn) call mpi_type_contiguous(nl4, mpi_real8, rhoc_mt , err)
!          if (trh)   call mpi_type_contiguous(nsp1   , mpi_real8, e_mt    , err)
         if (force) call mpi_type_contiguous(3      , mpi_real8, f_mt    , err)
         if (force .and. lov) call mpi_type_contiguous(3 * nbas, mpi_real8, drhos_mt, err)

         if (charge)call mpi_type_commit(rho_mt  , err)
         if (charge)call mpi_type_commit(rholm_mt, err)
         if (rl .or. rn) call mpi_type_commit(rhoc_mt, err)
!          if (trh)   call mpi_type_commit(e_mt, err)
         if (force) call mpi_type_commit(f_mt, err)
         if (force .and. lov) call mpi_type_commit(drhos_mt, err)

         if (charge) call mpi_allgatherv(mpi_in_place, 0,   rho_mt,rho  , tbc%d2count, tbc%d2amap,rho_mt , c2d%comm, err)
         if (charge) call mpi_allgatherv(mpi_in_place, 0, rholm_mt,rholm, tbc%d2count, tbc%d2amap,rholm_mt,c2d%comm, err)
         if (rl) call mpi_allgatherv(mpi_in_place, 0, rhoc_mt, rhoc, tbc%d2count, tbc % d2amap, rhoc_mt , c2d % comm, err)
         if (rn) call mpi_allgatherv(mpi_in_place, 0, rhoc_mt, rhon(1,1,1,1), tbc%d2count, tbc % d2amap, rhoc_mt ,c2d%comm,err)
         if (rn .and. nsp==2) &
            call mpi_allgatherv(mpi_in_place, 0, rhoc_mt, rhon(1,1,1,2), tbc%d2count, tbc % d2amap, rhoc_mt ,c2d%comm,err)
         if (trh) call mpi_allgatherv(mpi_in_place, 0, mpi_real8, esite(1,1), tbc%d2count, tbc % d2amap, mpi_real8, c2d%comm, err)
         if (trh .and. nsp==2) &
            call mpi_allgatherv(mpi_in_place, 0, mpi_real8, esite(1,2), tbc%d2count, tbc % d2amap, mpi_real8, c2d % comm, err)
         if (force) call mpi_allgatherv(mpi_in_place, 0, f_mt, f    , tbc%d2count, tbc % d2amap, f_mt, c2d % comm, err)
         if (force .and. lov .and. sc) &
            call mpi_allgatherv(mpi_in_place, 0, drhos_mt, drhos, tbc%d2count, tbc % d2amap, drhos_mt, c2d % comm, err)

         if (charge) call mpi_type_free(rho_mt  , err)
         if (charge) call mpi_type_free(rholm_mt, err)
         if (rl .or. rn) call mpi_type_free(rhoc_mt, err)
!          if (trh)   call mpi_type_free(e_mt, err)
         if (force) call mpi_type_free(f_mt, err)
         if (force .and. lov .and. sc) call mpi_type_free(drhos_mt, err)
        end if

! The above is unsightly but seems to be faster than the following
!       if (charge) call twolvl1d2d_rdst(tbc, mpi_real8, nsp1*nl , rho  )
!       if (charge) call twolvl1d2d_rdst(tbc, mpi_real8, nsp1*nl2, rholm)
!       if (rl)     call twolvl1d2d_rdst(tbc, mpi_real8, nl2*nl2 , rhoc )
!       if (trh)    call twolvl1d2d_rdst(tbc, mpi_real8, nsp1    , esite)
!       if (force)  call twolvl1d2d_rdst(tbc, mpi_real8, 3       , f    )
!       if (force .and. lov) &
!                 & call twolvl1d2d_rdst(tbc, mpi_real8, 3*nbas  , drhos)
!       if (rn)     call twolvl1d2d_rdst(tbc, mpi_real8, nl2*nl2 , rhon(1,1,1,1))
!       if (rn .and. nsp == 2) &
!                 & call twolvl1d2d_rdst(tbc, mpi_real8, nl2*nl2 , rhon(1,1,1,2))

         if (pv) call mpi_allreduce(mpi_in_place, thrpv, 1, mpi_real8, mpi_sum, c3d%comm, err)
!       call mpi_allreduce(mpi_in_place, entrpy, 1, mpi_real8, mpi_sum, c3d%comm, err)
         call tcx('netw')
      end if

!C --- Symmetrize forces ---
      if (fsym .and. force) call symfor(nbas,1,s_lat % symgr,nsgrp,s_lat%istab,f)

!C --- Symmetrise rho, rholm and esite ---
      if (fsym .and. (charge .or. trh)) call symre(lov,trh,nl,nsp1,nbas,nsgrp,s_lat%istab,rho,rholm,esite)


!C --- Symmetrise rhoc, rhon ---
      if (fsym .and. (rl .or. rn)) then
        call symtbd(nsgrp,s_lat % symgr,gd)
        if (rn) then
          do  isp = 1, nsp
            call symrtb(nl,1,nbas,nsgrp,s_lat % symgr,gd,s_lat%istab,rhon(1,1,1,isp))
          enddo
        endif
        if (rl) call symrtb(nl,1,nbas,nsgrp,s_lat % symgr,gd,s_lat%istab,rhoc)
      endif

      call tcx('bndtb')


!       write(420,'(2(9(9(x,f16.8),/),/))') rhoc
      end subroutine bndtb


!
!       do isp = 0, nsp-1
!           write(653, '("isp:", i0)') isp+1
!           do isite = 0, nsite-1
!               tmpr = ddot(nl4, hrs(isite*nl4 +isp*nsite*nl4+1), 1, hrs(isite*nl4 +isp*nsite*nl4+1), 1)
!               if (tmpr < 1d-6) cycle
!
!               write(653, '("  isite:", i0)') isite+1
!               do j = 0, nl2-1
!                   write(653,'("    ")', advance='no')
!                   do i = 0, nl2-1
!                       tmpr = hrs(i+j*nl2+isite*nl4 +isp*nsite*nl4+1)
!                       if (abs(tmpr) < 1e-8 ) tmpr = abs(tmpr)
!                       write(653, '(1x,f12.8)', advance='no') tmpr
!                   end do
!                   write(653,'("")')
!               end do
!           end do
!       end do
!
!
!
!       write(654, '("isp:", i0)') 1
!       do isite = 0, nsite-1
!           tmpr = ddot(nl4, srs(isite*nl4+1), 1, srs(isite*nl4+1), 1)
!           if (tmpr < 1d-6) cycle
!
!           write(654, '("  isite:", i0)') isite+1
!           do j = 0, nl2-1
!               write(654,'("    ")', advance='no')
!               do i = 0, nl2-1
!                   tmpr = srs(i+j*nl2+isite*nl4+1)
!                   if (abs(tmpr) < 1e-8 ) tmpr = abs(tmpr)
!                   write(654, '(1x,f12.8)', advance='no') tmpr
!               end do
!               write(654,'("")')
!           end do
!       end do
!
!
!
! ! for convergence conparisons (for Jorge)
!       integer, save :: counter = 0, iou(5), countmax = 10
!       character(len=100), save :: fmt_f, fmt_pgm
!       real(dp), allocatable, save ::  olde(:), oldz(:,:)
!       real(dp), allocatable :: diff(:)
!       real(dp), parameter :: pi = acos(-1.0_dp) !3.141592653589793_8


!             write(320,'(i5,x,i5,12(x,f16.8))') ikp, isp, rho
!             write(330,'(i5,x,i5,18(x,f16.8))') ikp, isp, eband(:,isp,ikp)
!             write(340,'(i5,x,i5,18(x,f16.8))') ikp, isp, wtkb(:,isp,ikp)

!             write(341,'(18(x,f16.8))') wtkb(:,isp,ikp)

!             write(321,'(i5,x,i5,12(x,f16.8))') ikp, isp, rho
!             write(342,'(18(x,f16.8))') wtkb(:,isp,ikp)
!             write(449,'(2(9(9(x,f16.8),/),/))') rhoc
!             if (ikp == nkp .and. isp == 2) write(439,'(2(9(9(x,f16.8),/),/))') rhoc


!       write(245, *) ldim, eband
!       stop

! ! convergence comparisons (za Jorge)
!       print *, 'nevmx, ldimx:', nevmx, ldimx
!       if (nevmx < ldimx) stop 'ASK FOR ALL EIGPAIRS FOR THE COMPARISON, STOPPING..'
!       if (counter > countmax) stop 'more iterations that envisioned, the pgm images will not be correct. STOPPING..'
!
!       if (counter == 0) then
!          allocate(olde(ldimx), oldz(ldimx, ldimx))
!          write(fmt_f  , "('(i3, ',i0,'(x, f20.16))')") ldimx
!          write(fmt_pgm, "('( ',i0,'(x,i5))')") ldimx+2
!          open(newunit=iou(1), file='diff_eval', action='write')
!          open(newunit=iou(2), file='diff_supf', action='write')
!          open(newunit=iou(3), file='evals', action='write')
!          open(newunit=iou(4), file='diff_c_fi.pgm', action='write')
!          open(newunit=iou(5), file='diff_eval.pgm', action='write')
!
!          write(iou(4), '(a,/,a,/,i4,x,i4,/,i6)') 'P2', "# abs(angle(C, C_{prev})) ", ldimx+2, countmax+2, 65535
!          write(iou(4), fmt_pgm) (32768, i = 1, ldimx+2)
!
!          write(iou(5), '(a,/,a,/,i4,x,i4,/,i6)') 'P2', "# abs eval diff", ldimx+2, countmax+2, 65535
!          write(iou(5), fmt_pgm) (32768, i = 1, ldimx+2)
!
!
!       else
!          allocate(diff(ldimx))
!
!          diff = eband(1:ldimx,1,1) - olde
!          write(iou(1), trim(fmt_f)) counter, diff
!
!          write(iou(5), fmt_pgm) 32768, int(65535*abs(diff)/0.25_dp), 32768 ! the max value is under just under 0.25 so i pick this for normalisation
!
!          diff = acos(sum(zll(:,:,1) * oldz, dim=1))/pi
!          where (diff > 0.5_dp) diff = diff - 1.0_dp
!
!          write(iou(4), fmt_pgm) 32768, int(65535*abs(diff)/0.5_dp), 32768 ! could have gone 1/0.5 -> 2
!
!          do i = 1, ldimx
!             diff(i) = diff(i) + real(i, dp)
!          end do
!          write(iou(2), trim(fmt_f)) counter, diff
!
!          deallocate(diff)
!
!       end if
!
!
!       write(iou(3), trim(fmt_f)) counter, eband(1:ldimx,1,1)
!
!
!       if (counter == countmax) then
!          write(iou(4), fmt_pgm) (32768, i = 1, ldimx+2)
!          write(iou(5), fmt_pgm) (32768, i = 1, ldimx+2)
!
!          do i = 1, size(iou)
!             close(iou(i))
!          end do
!          deallocate(olde, oldz)
!       else
!          olde = eband(1:ldimx,1,1)
!          oldz = zll(:,:,1)
!       end if
!
!       counter = counter + 1
!
! ! end of conv comparison


!       print *, 'p diagd:', pid
!       call mpi_barrier(mpi_comm_world, ierr)
!       print *, 'p diagd cont:', pid


!       call wkchk('CRASH')
!             do ikp = tbc % kmap(kbid) + 1, tbc % kmap(kbid+1)
!                do isp = 1, nsp1
!                   do i = 1, nevmx
!                      write(1300+c3d%id, '(x,es10.2)', advance='no') wtkb(i,isp,ikp)
!                   end do
!                end do
!                write(1300+c3d%id,'(" ")')
!             end do
!          flush(1300+c3d%id); flush(1400+c3d%id)
!          call mpi_barrier(c3d%comm, err)
!       call wkchk('CRASH')
