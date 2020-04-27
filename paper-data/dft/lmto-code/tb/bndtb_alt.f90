subroutine bndtb_alt(s_ctrl,s_bz,s_lat,s_spec,tbc,nbas,nl,nspc,nsp, &
    nsp1,lmx,idxdn,nclass,ips,ipc,nrclas,pv,force,zval,mull,  &
    wtbzm,nsite,iax,                                                &
    npr,xyzfrz,hrs,h0,dh,dhcf,vso,hso,srs,ds,dscf,pot0,             &
    efermi,sumev,entrpy,f,thrpv,esite,rho,rhoc,rhon,zos,            &
    drhos,ldim, idxsh)
    ! The arguments shall be the same as those bndtb, bar the k-point, overlap and spin dependent things for now.
    use mpi
    use tbprl
    use mod_ctx, only : cart_ctx_t
    use structures

    implicit none
    integer, parameter :: dp = 8

    type(str_ctrl) ::  s_ctrl
    type(str_bz) ::    s_bz
    type(str_lat) ::   s_lat
    type(str_spec) ::  s_spec(*)
    type(tbc_t), intent(in), target :: tbc

    integer nbas,nl,nsp,nspc,nsp1,nclass,mull,npln,nqmx,nsite,ldim, idxsh(*)
    integer lmx(*),idxdn(0:nl-1,*),ips(*),ipc(nbas),nrclas(*), iax(*),npr(*),iwk(nbas)
    real(dp) :: zval,wtbzm,efermi,sumev,entrpy,thrpv,egap
    real(dp) :: hrs(*),h0(*),dh(*),vso(*),hso(*),                                         &
        srs(*),ds(*),f(3,nbas),esite(nbas,nsp1),                                          &
        rho(nl,nsp1,nbas),rholm(nl*nl,nsp1,nbas),rhoc(nl*nl,nl*nl,nbas),                  &
        rhon(nl*nl,nl*nl,nbas,nsp1),zos(*),dhcf(*),dscf(*),drhos(3,nbas,nbas), pot0(nbas)
    logical pv,force,xyzfrz(3)

!C Local variables
    integer parg,ip,it,ldos,lncol,lstnr(3),ltb,getef,nste(2)
    integer i,j,k,l,ikp,nfilet,nfileb,                                &
        lihdim,ipr,iprx,nev,fopn,fopnT,iprint,icheck,ibas,indx,nchan, &
        ic,ifi,nevhi,nevlo,nat0,nlm,mull0,mull1,fopno,i1mach,lm,ll,   &
        imfil,ndos0,ndos,nfstg,isp,ispx,ldimx,ik0,ikstep, nevmx0,iomoms
    integer lmet,nkabc(3),nkp,ntet,mpsord,norder,nevmx

    integer nsgrp,oistab,oag,osymgr, lgunit, mode

!C --- for reconstructing RS hamiltonian
    integer ohrs,oors,ontab,opos,oclabl,oalpha
    real(dp) :: dum,del,evlo,qp(3),dval,dosw(2),plat(3,3),alat,       &
        srnge,swidth,efmax,dosef(2),qval(2),cv, gd(5,5,s_lat % nsgrp)
    logical :: metal,tetra,leig,doswt,fsym,charge,rl,trh,lso,donly,lmol, &
        lgamma,lgamsp,lov,bzmp,stoner,cryf,ocryf,rstrt,swnos,TBU,rn,bittst,cmdopt,excite
    character(len=80) :: outs
    character(len=45), parameter :: multyp(0:3,0:1) = reshape([ &
        ' Partial DOS, resolved by class and l        ',        &
        ' Bond charge, resolved by class pairs        ',        &
        ' Bond charge, resolved by {class, l} pairs   ',        &
        ' Bond charge, resolved by {class, l, m} pairs',        &
        ' Partial DOS, resolved by atom, l, and m     ',        &
        ' Bond charge, resolved by atom pairs         ',        &
        ' Bond charge, resolved by {atom, l} pairs    ',        &
        ' Bond charge, resolved by {atom, l, m} pairs '], [4,2]  )

    integer err, kbid, nl2, nl4, typesz
    logical :: mpiav
!       logical, parameter :: aprox = .true.
    logical, parameter :: aprox = .false.

    type(cart_ctx_t), pointer :: c3d, c2d, c1d
    real(8) :: tube_r8(8)
!       integer, allocatable :: wcounts(:), wdispls(:)
    integer :: eval_mt, rho_mt, rhoc_mt, rholm_mt, e_mt, f_mt, drhos_mt, wts_mt

    real(8), allocatable :: hk(:,:), frho(:,:), eband(:), zll(:,:)
    real(8) :: ef, gap, accu, thresh, lsxx_b
    integer :: nst, lsxx_p
    character(len=4) :: lsprog

    call tcn('bndtb_alt')


c3d => tbc % c3d
c2d => tbc % c2d
c1d => tbc % c1d

    kbid = c1d % id

    mpiav = c3d % sz > 1

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

    tetra  = ntet > 0
!     metal  = lmet /= 0
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

    if (.not. lgamma) call rx('BNDTB_ALT: ALT branch only works for Gamma point or molecule')

!C spin polarised gamma only, or molecule:
    if (nsp == 2) call rx('BNDTB_ALT: nsp=2 not implemented in the ALT branch')

    TBU = bittst(ltb,2**13)
    rn = TBU
    if (rl .or. rn) trh = .true.
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
!       lgamma = .false.
!       print *, 'reminder to fix lgamma in bndtb and secmtb'

    typesz = 1
    if (.not. lgamma) typesz = 2

    if (rl)  rhoc = 0.0_dp
    if (lov .and. force) drhos = 0.0_dp
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



    isp  = 1
    nsp  = 1
    nsp1 = 1
    nspc = 1

    if (c2d % lrt) print '(" ldim: ", i0)', ldim
    allocate(hk(ldim,ldim), frho(ldim,ldim), eband(ldim), zll(ldim,ldim))

    call tcn('blt0')
    call ztbloch(lgamma,qp,nl,1,1,1,nbas,plat,lmx,ipc,idxsh,nsite,iax,npr,hrs,vso,hso,.false., &
                                                                   ldim,hk,ldim,typesz)
    call tcx('blt0')


    if (c2d % lrt) then
        if (.not. cmdopt('--bgc',5,1,outs)) then
            call tcn('eig')
!         This is just to get proper ef and gap for each sc step. shall not be used in principle
            call h2rho_eig(ldim, hk, zval, efermi, gap, frho, sumev)
            call tcx('eig')
            print '(" band gap centre & width: ", f12.6, x, f12.6)', efermi, gap
        else
            read(outs, *) efermi
            print "(' band gap centre passed: ', f12.6)", efermi
            sumev = 0.0d0
            gap = 1.0d0
        end if
    end if

    if (c2d % sz > 1) then
        call mpi_bcast(efermi, 1, mpi_real8, c2d % rt, c2d % comm, err)
!         call mpi_bcast(gap   , 1, mpi_real8, c2d % rt, c2d % comm, err)
!         call mpi_bcast(sumev , 1, mpi_real8, c2d % rt, c2d % comm, err)
!         call mpi_bcast(frho  , ldim*ldim, mpi_real8, c2d % rt, c2d % comm, err)
    end if

    lsprog = 'msi'
    accu = 1d-12
    thresh = 0d0
    lsxx_p =  100000
    lsxx_b =  50



    if (c2d % lrt) then
        if (cmdopt('--linscl',8,1,outs)) read(outs, *) lsprog
        if (cmdopt('--lstol',7,1,outs)) read(outs, *) accu
        if (cmdopt('--lstol',7,2,outs)) read(outs, *) thresh
        if (cmdopt('--lsxxp',7,1,outs)) read(outs, *) lsxx_p
        if (cmdopt('--lsxxb',7,1,outs)) read(outs, *) lsxx_b ! beta may be estimated from the gap but will leave it for another time
        print "(' accuracy, threshold, ls++ beta, ls++ p: ', 3(f12.6,x),i0)", accu, thresh, lsxx_b, lsxx_p
    end if

    if (c2d % sz > 1) then
!     wasteful and latentful
        call mpi_bcast(lsprog, len(lsprog), mpi_character, c2d % rt, c2d % comm, err)
        call mpi_bcast(accu  , 1, mpi_real8, c2d % rt, c2d % comm, err)
        call mpi_bcast(thresh, 1, mpi_real8, c2d % rt, c2d % comm, err)
        call mpi_bcast(lsxx_b, 1, mpi_real8, c2d % rt, c2d % comm, err)
        call mpi_bcast(lsxx_p, 1, mpi_integer, c2d % rt, c2d % comm, err)
    end if

    call tcn('h2rho')
    if (trim(lsprog) == 'msi') then
        if (c2d % lrt) print '(" matrix sign iteration branch")'
        if (c2d % lrt) call h2rho_msi(ldim, hk, zval, efermi, gap, accu, thresh, frho)
        if (c2d % sz > 1) call mpi_bcast(frho, ldim*ldim, mpi_real8, c2d % rt, c2d % comm, err)
    else if (trim(lsprog) == 'ls++') then
        if (c2d % lrt) print '(" LS++ branch")'
        call h2rho_lsx(ldim, hk, zval, efermi, gap, accu, thresh, lsxx_p, lsxx_b, frho)
    end if
    call tcx('h2rho')

    deallocate(hk)

    thrpv = 0d0
    if (force) f = 0.0_dp
    if (trh) esite = 0.0_dp
    if (trh .or. charge) rho = 0.0_dp
    if (trh .or. charge) rholm = 0.0_dp



    call rho2frce(tbc,lmol,lgamma,plat,nbas,nl,nsp,isp,lmx,ipc,idxsh, &
        iwk,ldim,qp,nsite,iax,npr,typesz,h0,dh,charge,      &
        rl,trh,pv,force,f,thrpv,esite(1,isp),rho,rholm,rhoc,frho)

    nl2 = nl*nl
    nl4 = nl2*nl2

    if (c3d%sz > 1) then
        call tcn('netw')

! By utilising the mpi derived types the same recvcnts and displs can be used for different kinds of data in the allgather calls.
        if (c2d % sz > 1) then
        if (charge) call mpi_type_contiguous(nsp1*nl , mpi_real8, rho_mt  , err)
        if (charge) call mpi_type_contiguous(nsp1*nl2, mpi_real8, rholm_mt, err)
        if (rl .or. rn) call mpi_type_contiguous(nl4, mpi_real8, rhoc_mt , err)
        if (trh)   call mpi_type_contiguous(nsp1   , mpi_real8, e_mt    , err)
        if (force) call mpi_type_contiguous(3      , mpi_real8, f_mt    , err)
        if (force .and. lov) call mpi_type_contiguous(3 * nbas, mpi_real8, drhos_mt, err)

        if (charge)call mpi_type_commit(rho_mt  , err)
        if (charge)call mpi_type_commit(rholm_mt, err)
        if (rl .or. rn) call mpi_type_commit(rhoc_mt, err)
        if (trh)   call mpi_type_commit(e_mt, err)
        if (force) call mpi_type_commit(f_mt, err)
        if (force .and. lov) call mpi_type_commit(drhos_mt, err)

        if (charge) call mpi_allgatherv(mpi_in_place, 0,   rho_mt,rho  , tbc%d2count, tbc%d2amap,rho_mt , c2d%comm, err)
        if (charge) call mpi_allgatherv(mpi_in_place, 0, rholm_mt,rholm, tbc%d2count, tbc%d2amap,rholm_mt,c2d%comm, err)
        if (rl) call mpi_allgatherv(mpi_in_place, 0, rhoc_mt, rhoc, tbc%d2count, tbc % d2amap, rhoc_mt , c2d % comm, err)
        if (rn) call mpi_allgatherv(mpi_in_place, 0, rhoc_mt, rhon(1,1,1,1), tbc%d2count, tbc % d2amap, rhoc_mt ,c2d%comm,err)
        if (rn .and. nsp==2) &
            call mpi_allgatherv(mpi_in_place, 0, rhoc_mt, rhon(1,1,1,2), tbc%d2count, tbc % d2amap, rhoc_mt ,c2d%comm,err)
        if (trh)   call mpi_allgatherv(mpi_in_place, 0, e_mt, esite, tbc%d2count, tbc % d2amap, e_mt, c2d % comm, err)
        if (force) call mpi_allgatherv(mpi_in_place, 0, f_mt, f    , tbc%d2count, tbc % d2amap, f_mt, c2d % comm, err)
        if (force .and. lov) &
            call mpi_allgatherv(mpi_in_place, 0, drhos_mt, drhos, tbc%d2count, tbc % d2amap, drhos_mt, c2d % comm, err)

        if (charge) call mpi_type_free(rho_mt  , err)
        if (charge) call mpi_type_free(rholm_mt, err)
        if (rl .or. rn) call mpi_type_free(rhoc_mt, err)
        if (trh)   call mpi_type_free(e_mt, err)
        if (force) call mpi_type_free(f_mt, err)
        if (force .and. lov) call mpi_type_free(drhos_mt, err)
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

    call tcx('bndtb_alt')


!       write(420,'(2(9(9(x,f16.8),/),/))') rhoc

    end subroutine bndtb_alt



    subroutine h2rho_lsx(n,h,q, ef,gap,accu,thresh, p, beta, rho)
!     Use the sign iteration. The polynomial is equivalent to the vanderbuilt one.
        implicit none
        integer, intent(in) :: n, p
        real(8), intent(in) :: h(n,n), q, ef, gap, accu, thresh, beta
        real(8), intent(out) :: rho(n,n)

        interface
            subroutine lsdm_cheb_c(p, beta, mu, accu, thresh, r, c, h, rho) bind(c)
! extern "C" void lsdm_cheb_c(int p, double beta, double mu, double accu, double thresh,
!                             int r, int c, double *h, double *rho /*, MPI_Comm *comm */)
                use iso_c_binding, only : c_int, c_double
                implicit none
                integer(c_int), intent(in), value :: p, r, c
                real(c_double), intent(in), value :: beta, mu, accu, thresh
                real(c_double), intent(in) :: h(*)
                real(c_double), intent(out) :: rho(*)
            end subroutine lsdm_cheb_c
        end interface

        call lsdm_cheb_c(p, beta, ef, accu, thresh, n, n, h, rho)

        call dscal(n*n,2.0d0,rho,1)

    end subroutine h2rho_lsx



    subroutine h2rho_eig(n,h,q, ef,gap,rho,sumev)

        implicit none

        integer, intent(in) :: n
        real(8), intent(in) :: h(n,n), q
        real(8), intent(out) :: rho(n,n), ef, gap, sumev

        integer :: nst, m, nz
        real(8) :: ev(n), z(n,n)

        call tcn('rho_eig')
        call unievp1(.false., .false., 1, 'v', 'a', 'u', n, h, n, rho, n, 0.0d0, 0.0d0, 0, n, m, nz, ev, z, n) ! rho here is just a placeholder
        nst = (int(q) + 1)/2
        ef = 0.5*(ev(nst+1) + ev(nst))
        gap =     ev(nst+1) - ev(nst)
        sumev = 2*sum(ev(:nst))
        call dgemm('n', 't', n, n, nst, 2.0d0, z, n, z, n, 0.0d0, rho, n)
        call tcx('rho_eig')
    end subroutine h2rho_eig



    subroutine h2rho_msi(n,h,q, ef,gap,accu,thresh, rho)
!     Use the sign iteration. The polynomial is equivalent to the vanderbuilt one.
!     This implementation is not actually linear scalling but cubic as it uses dense
!     matrices for simplicity. Extra tricks are applied to reproduce real sparcity with thresholding.
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: h(n,n), q, ef, gap, accu, thresh
        real(8), intent(out) :: rho(n,n)

        integer :: i

        rho = h
        do i = 1, n
            rho(i,i) = rho(i,i) - ef
        end do
        call matsign_itr(n, rho, 200, accu, thresh) !1d-12 dava syotvetstvashti si sili
        do i = 1, n
            rho(i,i) = rho(i,i) - 1.0_8
        end do
        rho = -rho
    end subroutine h2rho_msi


!     subroutine h2rho_msi_cu(n,h,q, ef,gap,rho)
! !     Same as h2rho_msi but faster on the GPU (bar the copy line below)
! !     Commented out to not bother with linking the cuda stuff (or a fake r2rho_host)
!         implicit none
!         integer, intent(in) :: n
!         real(8), intent(in) :: h(n,n), q, ef, gap
!         real(8), intent(out) :: rho(n,n)
!         interface
!             subroutine h2rho_cu(d, hx, ef, maxit, efilt) bind(c, name='h2rho_host')
!                 use iso_c_binding, only: c_int, c_double
!                 implicit none
!                 integer(c_int), intent(in), value :: d, maxit
!                 real(c_double), intent(in), value :: ef, efilt
!                 real(c_double), intent(inout) :: hx(*)
!             end subroutine h2rho_cu
!         end interface
!
!         rho = h
!         call h2rho_cu(n, rho, ef, 200, 1d-12)
!
!     end subroutine h2rho_msi_cu


    subroutine matsign_itr(d, x, maxit, e, thresh)
!       X must be rescaled by the conditioning number    X = A/||A|| (largest abs eigenvalue)
!       approximation to the sign function: sign(X) = X(X^2)^(-1/2)
!       Xn+1 = 1/2 Xn(3I-Xn^2)
!       Xinf -> sign(A)
!       use tbprl
      implicit none
      integer, parameter :: dp = 8

      real(dp), intent(in) :: e, thresh
      integer, intent(in) :: d, maxit
      real(dp), intent(inout) :: x(d,d)

      real(dp) :: c, ix2fnrm, ex2fnrm, esq
      real(dp), allocatable :: t(:,:,:), x2(:,:)
      integer :: i, j, k, it, d2, sparsity(2)
      logical :: last

      procedure(real(dp)) :: ddot

      d2 = d*d
      esq = sqrt(e)

      allocate(t(d,d,0:1),x2(d,d))



      it = 0
      j = 0
      k = 1
      do
         if (it /= 0) then
            j = mod(it  , 2)
            k = mod(it+1, 2)
            call dgemm('t','n',d,d,d,-1.0_dp,t(1,1,j),d,t(1,1,j),d,0.0_dp,x2,d)
         else
            call dgemm('t','n',d,d,d,-1.0_dp,x,d,x,d,0.0_dp,x2,d)
         end if

         where ( -thresh < x2 .and. x2 < thresh ) x2 = 0.0_dp
         sparsity(1) = count(-thresh < x2 .and. x2 < thresh)

         ex2fnrm = esq*ddot(d2, x2, 1, x2, 1)

         do i = 1,d
            x2(i,i) = x2(i,i) + 1.0_dp
         end do

         ix2fnrm = ddot(d2, x2, 1, x2, 1)

         last = ix2fnrm < ex2fnrm
         write(*,'(a,i3,a,2(x,f18.12),a,l)',advance='no') &
            & '  sign itr:',it,' ||I-X^2||F, efilt^0.5*||X^2||F:',ix2fnrm,ex2fnrm,' last:',last

         do i = 1,d
            x2(i,i) = x2(i,i) + 2.0_dp
         end do


         if (last) then
            call dgemm('t','n',d,d,d,0.5_dp,t(1,1,j),d,x2,d,0.0_dp,x,d)
            where ( -thresh < x .and. x < thresh ) x = 0.0_dp
            sparsity(2) = count( -thresh < x .and. x < thresh)
            write (*,'(" sparsity:", 2(x,f6.3), " %", 2(x,i0))') 100*real(sparsity,dp)/real(d2,dp), sparsity
            exit
         else
            if (it /= 0) then
               call dgemm('t','n',d,d,d,0.5_dp,t(1,1,j),d,x2,d,0.0_dp,t(1,1,k),d)
            else
               call dgemm('t','n',d,d,d,0.5_dp,x,d,x2,d,0.0_dp,t(1,1,k),d)
            end if
            where ( -thresh < t(:,:,k) .and. t(:,:,k) < thresh ) t(:,:,k) = 0.0_dp
            sparsity(2) = count( -thresh < t(:,:,k) .and. t(:,:,k) < thresh)
            write (*,'(" sparsity:", 2(x,f6.3), " %", 2(x,i0))') 100*real(sparsity,dp)/real(d2,dp), sparsity
            it  = it + 1
         end if

      end do


      deallocate(x2,t)

    end subroutine matsign_itr





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
!       print *, 'nbmax, ldimx:', nbmax, ldimx
!       if (nbmax < ldimx) stop 'ASK FOR ALL EIGPAIRS FOR THE COMPARISON, STOPPING..'
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
!                   do i = 1, nbmax
!                      write(1300+c3d%id, '(x,es10.2)', advance='no') wtkb(i,isp,ikp)
!                   end do
!                end do
!                write(1300+c3d%id,'(" ")')
!             end do
!          flush(1300+c3d%id); flush(1400+c3d%id)
!          call mpi_barrier(c3d%comm, err)
!       call wkchk('CRASH')
