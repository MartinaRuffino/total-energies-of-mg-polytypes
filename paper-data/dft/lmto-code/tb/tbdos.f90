    module mod_tbdos

    use mpi
    use tbprl
    use mod_ctx, only : cart_ctx_t

    implicit none

    private

    procedure(logical) :: cmdopt

    integer, parameter :: p = 8
    real(p), parameter :: pi = 3.14159265358979323846_p
    complex(p), parameter :: znul = 0
    complex(p), parameter :: zone = 1


    public :: plotdos

    contains

    subroutine plotdos(tbc, nev, nsp, wtkp, evls, z, ldim, ef, &
        lov, lgamma, nl, nbas, npairs, lmx, ipc, indxsh, iax, npr, sr, plat, qp, vso, hso)
! Generate total and (optionally) projected DOS
!
!       g_e  = sum_n g_ne
!       g_ie = sum_n z_in g_ne z'_ni = sum_n zz'_in im(g_ne)
!
! Known issues, to be fixed in the (hopefully near) future:
!
!   Scalapack is not supported yet error will be thrown if distributed z is passed, it is however written in a manner to allow easy upgrade.
!   No symmetrisation is not performed.
!   Only the basic --dos args are parsed and none of the --pdos ones are.
!   The Im(1/(e-i*eta))/pi delta function is hardcoded, will be better to use the integration method specifid in the ctrl file instead as bzwts does.
!
! Desirable features:
!   Write to easier to files easier to load with numpy.loadtxt and possibly hdf5 files.
!

! The eigenvalues are already gathered at the 3d root process but the respective parts at other processes are still there.

! c' s c = i
! s c = c1'
! c' s = c1
! s = c1' c1
!
! h c = s c e
! h = s c e c1
! h = c1' e c1
! s1 h = c e c1
! s1 h s1 = c e c'
! h s1 = s c e c'


! tb parallel configuration structure
        type(tbc_t), intent(in), target :: tbc
! number of eigenstates
        integer, intent(in) :: nev
! number of spins
        integer, intent(in) :: nsp
! basis size (used for the projected dos)
        integer, intent(in) :: ldim
! Fermi level, passed only to be printed in the dos files
        real(p), intent(in) :: ef
! BZ integration weights
        real(p), intent(in) :: wtkp(:)
! eigenvalues(nev,nsp,kpts)
! it is important to keep the implicit dimensioning because the k point slices differ between processes (and root holds copies of all of them)
        real(p), intent(in) :: evls(:,:,:)
! eigenvectors(typesz,ldim,nev,nsp,kpts), similarly to evls, k point slices differ between processes and implicit dims shall be kept. typesz is 1 for real/float and 2 for complex packed in 2 reals.
        real(p), intent(in), target :: z(:,:,:,:,:)
! says whether overlap is to be applied
        logical, intent(in) :: lov
! overlap matrix, used only if lov == .true., real space neigh pair structure

        logical, intent(in) :: lgamma
        integer, intent(in) :: nl, nbas, npairs
        integer, intent(in) :: lmx(0:*), ipc(nbas), indxsh(nl**2*nbas), iax(*), npr(*)
        real(p), intent(in) :: sr(*), plat(3,3), qp(3,*), vso(*), hso(*)

        character(len=128) :: outs
        character(len=16) :: ext
        logical :: txt2d
        integer :: n, k, s, v, i, u, kbid, err, typesz
        real(p) :: emin, emax, eta, eta2, ed, tmp
        real(p), allocatable :: e(:), tdos(:,:), pdos(:,:,:), zz(:,:), g(:,:), w(:)

        integer :: lsz, ldh
        real(p), allocatable :: sk(:,:,:)
        real(p), pointer :: zi(:,:,:)

        type(cart_ctx_t), pointer :: c3d, c2d, c1d

        procedure(integer) :: fextg

        if (tbc % sl) call rx('tb dos plotting does not work with scalapack yet')
! from here on assuming that only k-points are spread across processes.


!         c3d => tbc % c3d
!         c2d => tbc % c2d
        c1d => tbc % c1d

        if (c1d % lrt) then
!             write(*,'("DOS mesh size and smearing? ")', advance='no')
            n = 1001
            eta = 0.005_p
            emin = 0
            emax = 0
            txt2d = .false.
            call parse_dos_opts(n, emin, emax, eta, txt2d)
            err = fextg(ext)
        end if

        call mpi_barrier(c1d % comm, err) !not strictly necessary, just in case there is diff timeout in the bcast
!         call mpi_bcast(txt2d, 1, mpi_logical, c1d % rt, c1d % comm, err)
        call mpi_bcast(n, 1, mpi_integer, c1d % rt, c1d % comm, err)
        call mpi_bcast(eta, 1, mpi_real8, c1d % rt, c1d % comm, err)
        call mpi_bcast(emin, 1, mpi_real8, c1d % rt, c1d % comm, err)
        call mpi_bcast(emax, 1, mpi_real8, c1d % rt, c1d % comm, err)

        eta2 = eta * eta

        kbid = c1d % id

! window not given so figure it out
        if (emin == 0 .and. emax == 0) then
            emin = minval(evls(:, :, tbc % kmap(kbid) + 1 : tbc % kmap(kbid+1)))
            emax = maxval(evls(:, :, tbc % kmap(kbid) + 1 : tbc % kmap(kbid+1)))

            call mpi_allreduce(mpi_in_place, emin, 1, mpi_real8, mpi_min, c1d % comm, err)
            call mpi_allreduce(mpi_in_place, emax, 1, mpi_real8, mpi_min, c1d % comm, err)
        end if

        allocate(e(0:n-1))
        e = [(emin + real(i,8)*(emax-emin)/real(n,8), i = 0, n-1)]

        allocate(w(tbc % kmap(kbid) + 1 : tbc % kmap(kbid+1)))

        tmp = nsp*pi*sum(abs(wtkp(1:tbc%kmap(c1d%sz))))
        forall (k = tbc % kmap(kbid) + 1 : tbc % kmap(kbid+1)) w(k) = abs(wtkp(k)) / tmp

        allocate(tdos(0:n-1, nsp))
        tdos = 0

        do k = tbc % kmap(kbid) + 1, tbc % kmap(kbid+1)
            do s = 1, nsp
                do i = 0, n-1
                    tmp = 0
                    do v = 1, nev
                        ed =  e(i) - evls(v,s,k)
                        tmp = tmp + eta/(ed*ed + eta2)
                    end do
                    tdos(i,s) = tdos(i,s) + tmp * w(k)
                end do
            end do
        end do

        call mpi_reduce(mpi_in_place, tdos, n*nsp, mpi_real8, mpi_sum, c1d % rt, c1d % comm, err)
        if (c1d % lrt) then
            if (txt2d) then
                call write_2ddos_fmt('dos'//trim(ext), n, 1, nsp, e, tdos)
            else
                call write_pldos_fmt('dos'//trim(ext), emin, emax, n, 1, nsp, ef, 0.0_p, tdos)
            end if
        end if

        if (cmdopt('--pdos',6,0,outs) .or. cmdopt('--mull',5,0,outs)) then
            allocate(pdos(0:n-1, ldim, nsp))
            allocate(g(nev, 0:n-1))

            allocate(zz(ldim, nev))

            typesz = size(z, 1)

            if (lov) then
                ldh = tbc%loclsz(1)
                allocate(sk(typesz,ldh,tbc%loclsz(2)))
                allocate(zi(typesz,ldh,tbc%loclsz(2)))
            end if

            do k = tbc % kmap(kbid) + 1, tbc % kmap(kbid+1)
! s is spin independent so only using isp=1
                if (lov) call ztblochd(tbc, lgamma, qp(1,k), nl, 1, 1, 1, nbas, plat, lmx, ipc, indxsh, &
                                        npairs, iax, npr, sr, vso, hso, .false., ldim, sk, ldh, typesz)

                do s = 1, nsp
                    if (lov) call  zgemm('c','n',ldim,nev,ldim,zone,sk,ldh,z,ldh,znul,zi,ldh)
                    if (.not. lov) zi => z(:,:,:,s,k)

                    do v = 1, nev
                        do i = 1, ldim
                            zz(i,v) = sum(zi(1:typesz,i,v)*z(1:typesz,i,v,s,k))
                        end do
                    end do

                    do i = 0, n-1
                        do v = 1, nev
                            ed =  e(i) - evls(v,s,k)
                            g(v,i) = eta/(ed*ed + eta2)
                        end do
                    end do

                    tmp = 1; if (k == tbc % kmap(kbid) + 1) tmp = 0
                    call dgemm('t','t',n,ldim,nev,w(k),g,nev,zz,ldim,tmp,pdos(0,1,s),n)
                end do
            end do

            if (lov) deallocate(zi)
            if (lov) deallocate(sk)
            deallocate(zz)
            deallocate(g)

            call mpi_reduce(mpi_in_place, pdos, n*ldim*nsp, mpi_real8, mpi_sum, c1d % rt, c1d % comm, err)

            if (c1d % lrt) then
                if (txt2d) then
                    call write_2ddos_fmt('dosp'//trim(ext), n, ldim, nsp, e, pdos)
                else
                    call write_pldos_fmt('dosp'//trim(ext), emin, emax, n, ldim, nsp, ef, 0.0_p, pdos)
                end if
            end if

            deallocate(pdos)
        end if

! deallocating tdos so late to allow a sanity checks on pdos
        deallocate(tdos)
        deallocate(e)
        deallocate(w)

    end subroutine plotdos

    subroutine write_2ddos_fmt(fln, ne, nc, ns, e, dos)
! filename, the file will be truncated
        character(len=*), intent(in) :: fln
! number of energies
        integer, intent(in) :: ne
! number of spins
        integer, intent(in) :: ns
! number of chanels
        integer, intent(in) :: nc
! energy mesh
        real(p), intent(in) :: e(ne)
! DOS to write out
        real(p), intent(in) :: dos(ne,nc,ns)

        integer :: i, j, s, u

        open(newunit=u, file=trim(fln), action='write')
        write(u, "('#% rows ',i0,' cols ',i0)") ne, nc*ns+1
        do j = 1, ne
            write(u, '(es18.8)', advance='no') e(j)
            do s = 1, ns
                do i = 1, nc
                    write(u, '(x,es18.8)', advance='no') dos(j,i,s)
                end do
            end do
            write(u, '()')
        end do
        close(u)
    end subroutine write_2ddos_fmt

    subroutine write_pldos_fmt(fln, emin, emax, ne, nc, ns, ef, dlt, dos)
! Write dos in the pldos format v1
! https://www.questaal.org/docs/input/data_format/#the-dos-file

! number of energies
        integer, intent(in) :: ne
! number of spins
        integer, intent(in) :: ns
! number of chanels
        integer, intent(in) :: nc
! first energy
        real(p), intent(in) :: emin
! ending energy
        real(p), intent(in) :: emax
! Fermi energy
        real(p), intent(in) :: ef
! sampling broadening
        real(p), intent(in) :: dlt
! DOS to write out
        real(p), intent(in) :: dos(ne,nc,ns)
! filename, the file will be truncated
        character(len=*), intent(in) :: fln

        integer :: i, s, u

        open(newunit=u, file=trim(fln), action='write')
        write(u, '(2(x,f12.6),3(x,i0),2(x,f12.6),x,i0)') emin, emax, ne, nc, ns, ef, dlt, 1

        do i = 1, nc
            do s = 1, ns
                write(u, '(5(x,es18.8))') dos(:,i,s)
            end do
        end do

        close(u)
    end subroutine write_pldos_fmt

    subroutine parse_dos_opts(npts, emin, emax, eta, txt2d)
! try to parse DOS command line options --dos:rdm:npts=1001:window=-1,1
        integer, intent(inout) :: npts
        real(p), intent(inout) :: emin, emax, eta
        logical, intent(inout) :: txt2d

        logical :: ltmp
        integer :: i, j
        character(len=128) :: outs, t
        character :: s

        t = '--dos'
        i = len_trim(t)
        ltmp = cmdopt(t,i,0,outs)
        i = i + 1
        s = outs(i:i)
        do j = i, len_trim(outs)
            if (outs(j:j) == s) outs(j:j) = ' '
            if (outs(j:j) == ',') outs(j:j) = ' '
        end do

        i = token_val_idx(' npts=', outs)
        if (i > 0) read(outs(i:), *) npts

        i = token_val_idx(' window=', outs)
        if (i > 0) read(outs(i:), *) emin, emax

        i = token_val_idx(' eta=', outs)
        if (i > 0) read(outs(i:), *) eta

        txt2d = token_val_idx(' rdm', outs) > 0

    end subroutine parse_dos_opts

    integer function token_val_idx(tok, line)
        character(len=*), intent(in) :: tok, line
        token_val_idx = index(line, trim(tok))
        if (token_val_idx > 0) token_val_idx = token_val_idx + len_trim(tok)
    end function token_val_idx


    end module mod_tbdos
