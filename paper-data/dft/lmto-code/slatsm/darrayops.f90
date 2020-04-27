
   subroutine darray_scatter(mtype, a, c, desc, c2d)
!- Scatter a global array to local 2d block-cyclic distributed pieces.
!-------------------------------------------------------------------------------
!i Inputs:
!i   mtype : the basic MPI type of the data (only mpi_real8 tested)
!i   a     : global matrix to be scattered
!i   desca : descriptor for a, see unievp.f90 for details
!i   c2d   : 2D cartesian communicator
!i
!o Outputs:
!o   c     : local 2d block-cyclic distributed pieces, The actual size may be
!o            smaller than the container c for alignment requirements by the mpi_darray type and the mpi_alltolallw routine
!r Remarks:
!r  the local containers 'c' shall be of the same dimension on all processes

   use mpi
   use mod_ctx
   implicit none
   integer, parameter :: dp = 8
!         rather limited for now only 0,0 can send

   integer, intent(in) :: mtype

   real(dp), intent(in), target :: a(*)
   real(dp), intent(inout) :: c(*)
   integer, intent(inout) :: desc(9)
   type(cart_ctx_t), intent(in) :: c2d

   real(dp), pointer :: sp(:)
   real(dp), target :: nula(1)

   integer :: nb, lr, lc, nprow, npcol
   integer, parameter :: ndims = 2

   integer :: n, nproc, comm
   integer :: err, stat(mpi_status_size)

   integer :: pi, pj, ciproc, i, typesz
   integer(mpi_address_kind) :: lb, extent
   integer, dimension(ndims) :: globsz, distribs, blcksz, psizes, blocks, bl_prc, loclsz, realsz
   integer, dimension(:), allocatable :: dars, scnts, sdispls, rcnts, rdispls, rtypes, stypes, szs

   logical :: srv


   typesz = 1
   if (mtype == mpi_complex16)  typesz = 2

   srv = c2d % lrt
   comm   = c2d % comm
   nproc  = c2d % sz
   psizes = c2d % szs(0:1)

   realsz = desc(3:4)
   blcksz = desc(5:6)
   blocks = (realsz - 1) / blcksz + 1
   bl_prc = (blocks - 1) / psizes + 1
   loclsz = bl_prc * blcksz
   globsz = loclsz * psizes


   allocate(dars(0:nproc-1), scnts(0:nproc-1), sdispls(0:nproc-1), &
      & rcnts(0:nproc-1), rdispls(0:nproc-1), rtypes(0:nproc-1), stypes(0:nproc-1), szs(0:nproc-1))

   distribs = [mpi_distribute_cyclic, mpi_distribute_cyclic]

   if (srv) then
      do pi = 0, psizes(1)-1
         do pj = 0, psizes(2)-1
            call mpi_cart_rank(comm, [pi, pj], i, err)
            call mpi_type_create_darray(nproc, i, 2, globsz, distribs, blcksz, psizes, mpi_order_fortran, mtype, dars(i), err)
            call mpi_type_commit(dars(i), err)
!             call mpi_type_size(dars(i), szs(i), err)
!             szs(i) = szs(i)/dp
!             ln = (ns/psizes + min(mod(ns, psizes)/([pi,pj]+1), 1))*blcksz
!             if (product(ln) /= szs(i)) print *, 'size mismatch at pcoords', pi, pj, 'type size:', szs(i), ' est size:', product(ln)
         end do
      end do
   end if


!    ln = (ns/psizes + min(mod(ns, psizes)/(pc+1), 1))*blcksz
!          if (any( ln /= [lr*nb, lc*nb])) print *, 'size mismatch: ln, [lr*nb, lc*nb]:', ln, lr*nb, lc*nb

!          deallocate(c)
!          allocate(c(0:ln(1)-1, 0:ln(2)-1))
!    desc(9) = ln(1)
!
!    if ((.not. allocated(c))) then
!       allocate(c(ln(1), ln(2)))
!    else
!       if (any(size(c) /= ln)) then
!          deallocate(c)
!          allocate(c(ln(1), ln(2)))
!       end if
!    end if


   sdispls = 0
   rdispls = 0
   rtypes = mtype
   rcnts = 0
   rcnts(0) = product(loclsz)


   if (srv) then
      scnts = 1
      stypes = dars
      sp => a(1:typesz*product(globsz))
!             print *, 'a:', size(a,1), size(a,2)
!             call print_c(a, 99, 'a', 'f4.1', achar(10))
   else
      scnts = 0
      stypes = mtype
      sp => nula
   end if

!          print *, 'c:', size(c,1), size(c,2), ciproc


   call mpi_alltoallw(sp, scnts, sdispls, stypes, c, rcnts, rdispls, rtypes, comm, err)

   if (srv) then
      do i = 0, nproc-1
         call mpi_type_free(dars(i), err)
      end do
   end if

!          call print_c(c, 100+ciproc, 'c', 'f4.1', achar(10))


!          call mpi_finalize(err)
!          stop

   end subroutine darray_scatter




   subroutine darray_gather(mtype, a, c, desc, c2d)
!- Recreate the global array from the local 2d block-cyclic distributed local pieces.
!-------------------------------------------------------------------------------
!i Inputs:
!i   mtype : the basic MPI type of the data (only mpi_real8 tested)
!i   c     : local 2d block-cyclic distributed pieces, The actual size may be
!i            smaller than the container c for alignment requirements by the mpi_darray type and the mpi_alltolallw routine
!i   desca : descriptor for a, see unievp.f90 for details
!i   c2d   : 2D cartesian communicator
!i
!o Outputs:
!o   a     : global matrix to be asembled
!
!r Remarks:
!r  the local containers 'c' shall be of the same dimension on all processes


   use mpi
   use mod_ctx
   implicit none
   integer, parameter :: dp = 8
   type(cart_ctx_t), intent(in) :: c2d
   integer, intent(in) :: mtype
   real(dp), intent(out), target :: a(*)
   real(dp), intent(in) :: c(*)
   integer, intent(in) :: desc(9)

   real(dp), pointer :: sp(:)
   real(dp), target :: nula(1)

   integer :: nb, lr, lc, nprow, npcol
   integer, parameter :: ndims = 2

   integer :: n, nproc, comm
   integer :: err, stat(mpi_status_size)

   integer :: pi, pj, ciproc, iproc, i, typesz
   integer(mpi_address_kind) :: lb, extent
   integer, dimension(ndims) :: globsz, distribs, blcksz, psizes, blocks, bl_prc, loclsz, realsz

   integer, dimension(:), allocatable :: dars, scnts, sdispls, rcnts, rdispls, rtypes, stypes, szs

   logical :: periods(2), srv



   typesz = 1
   if (mtype == mpi_complex16) typesz = 2

   srv    = c2d % lrt
   comm   = c2d % comm
   nproc  = c2d % sz
   psizes = c2d % szs(0:1)

   realsz = desc(3:4)
   blcksz = desc(5:6)
   blocks = (realsz - 1) / blcksz + 1
   bl_prc = (blocks - 1) / psizes + 1
   loclsz = bl_prc * blcksz
   globsz = loclsz * psizes


   allocate(dars(0:nproc-1), scnts(0:nproc-1), sdispls(0:nproc-1), &
      & rcnts(0:nproc-1), rdispls(0:nproc-1), rtypes(0:nproc-1), stypes(0:nproc-1)) ! szs(0:nproc-1)
   distribs = [mpi_distribute_cyclic, mpi_distribute_cyclic]

   do pi = 0, psizes(1)-1
      do pj = 0, psizes(2)-1
         call mpi_cart_rank(comm, [pi, pj], i, err)
         call mpi_type_create_darray(nproc, i, ndims, globsz, distribs, blcksz, psizes, mpi_order_fortran, mtype, dars(i), err)
         call mpi_type_commit(dars(i), err)
!          call mpi_type_size(dars(i), szs(i), err)
!          szs(i) = szs(i)/dp
!          ln = (ns/psizes + min(mod(ns, psizes)/([pi,pj]+1), 1))*blcksz
!          if (product(ln) /= szs(i)) print *, 'wh size mismatch at pcoords', pi, pj, 'type size:', szs(i), ' est size:', product(ln)
      end do
   end do

!    ln = (ns/psizes + min(mod(ns, psizes)/(pc+1), 1))*blcksz
!          if (any( ln /= [lr*nb, lc*nb])) print *, 'wh size mismatch: ln, [lr*nb, lc*nb]:', ln, lr*nb, lc*nb

!
!          if ((.not. allocated(c))) then
!             allocate(c(ln(1), ln(2)))
!          else
!             if (any(size(c) /= ln)) then
!                deallocate(c)
!                allocate(c(ln(1), ln(2)))
!             end if
!          end if
!

   sdispls = 0
   rdispls = 0
   stypes = mtype
   scnts = 0
   scnts(0) = product(loclsz) !szs(iproc)
   rtypes = dars

   if (srv) then
      rcnts = 1
      sp => a(1:typesz*product(globsz))
!             print *, 'a:', size(a,1), size(a,2)
!             call print_c(a, 99, 'a', 'f4.1', achar(10))
   else
      rcnts = 0
      sp => nula ! looks like mpi does access this even if it does not do anything with it since count=0, value of null() seems to cause weird issues.
!       sp => null()
   end if

!          print *, 'c:', size(c,1), size(c,2), ciproc


   call mpi_alltoallw(c, scnts, sdispls, stypes, sp, rcnts, rdispls, rtypes, comm, err)

   do i = 0, nproc-1
      call mpi_type_free(dars(i), err)
   end do

!    print *, 'end gthr'
   end subroutine darray_gather










   subroutine darray_allgather(mtype, a, c, desc, c2d)
!- Recreate the global array from the local 2d block-cyclic distributed local pieces.
!-------------------------------------------------------------------------------
!i Inputs:
!i   mtype : the basic MPI type of the data (only mpi_real8 tested)
!i   c     : local 2d block-cyclic distributed pieces, The actual size may be
!i            smaller than the container c for alignment requirements by the mpi_darray type and the mpi_alltolallw routine
!i   desca : descriptor for a, see unievp.f90 for details
!i   c2d   : 2D cartesian communicator
!i
!o Outputs:
!o   a     : global matrix to be asembled
!
!r Remarks:
!r  the local containers 'c' shall be of the same dimension on all processes


   use mpi
   use mod_ctx
   implicit none
   integer, parameter :: dp = 8
   type(cart_ctx_t), intent(in) :: c2d
   integer, intent(in) :: mtype
   real(dp), intent(out), target :: a(*)
   real(dp), intent(in) :: c(*)
   integer, intent(in) :: desc(9)

   real(dp), pointer :: sp(:)
   real(dp), target :: nula(1)

   integer :: nb, lr, lc, nprow, npcol
   integer, parameter :: ndims = 2

   integer :: n, nproc, comm
   integer :: err, stat(mpi_status_size)

   integer :: pi, pj, ciproc, iproc, i, typesz
   integer(mpi_address_kind) :: lb, extent
   integer, dimension(ndims) :: globsz, distribs, blcksz, psizes, blocks, bl_prc, loclsz, realsz

   integer, dimension(:), allocatable :: dars, scnts, sdispls, rcnts, rdispls, rtypes, stypes, szs

   logical :: periods(2), srv

   typesz = 1
   if (mtype == mpi_complex16) typesz = 2

   srv = c2d % lrt
   comm   = c2d % comm
   nproc  = c2d % sz
   psizes = c2d % szs(0:1)

   realsz = desc(3:4)
   blcksz = desc(5:6)
   blocks = (realsz - 1) / blcksz + 1
   bl_prc = (blocks - 1) / psizes + 1
   loclsz = bl_prc * blcksz
   globsz = loclsz * psizes


   allocate(dars(0:nproc-1), scnts(0:nproc-1), sdispls(0:nproc-1), &
      & rcnts(0:nproc-1), rdispls(0:nproc-1), rtypes(0:nproc-1), stypes(0:nproc-1), szs(0:nproc-1))

   distribs = [mpi_distribute_cyclic, mpi_distribute_cyclic]

   do pi = 0, psizes(1)-1
      do pj = 0, psizes(2)-1
         call mpi_cart_rank(comm, [pi, pj], i, err)
         call mpi_type_create_darray(nproc, i, ndims, globsz, distribs, blcksz, psizes, mpi_order_fortran, mtype, dars(i), err)
         call mpi_type_commit(dars(i), err)
!          call mpi_type_size(dars(i), szs(i), err)
!          szs(i) = szs(i)/dp
!          ln = (ns/psizes + min(mod(ns, psizes)/([pi,pj]+1), 1))*blcksz
!          if (product(ln) /= szs(i)) print *, 'wh size mismatch at pcoords', pi, pj, 'type size:', szs(i), ' est size:', product(ln)
      end do
   end do

!    ln = (ns/psizes + min(mod(ns, psizes)/(pc+1), 1))*blcksz
!          if (any( ln /= [lr*nb, lc*nb])) print *, 'wh size mismatch: ln, [lr*nb, lc*nb]:', ln, lr*nb, lc*nb

!
!          if ((.not. allocated(c))) then
!             allocate(c(ln(1), ln(2)))
!          else
!             if (any(size(c) /= ln)) then
!                deallocate(c)
!                allocate(c(ln(1), ln(2)))
!             end if
!          end if
!

   sdispls = 0
   rdispls = 0
   stypes = mtype
!    scnts = 0
!    scnts(0) = product(loclsz) !szs(iproc)
   scnts = product(loclsz)
   rtypes = dars
   rcnts = 1

!    if (srv) then
!       rcnts = 1
!       sp => a(1:typesz*product(globsz))
! !             print *, 'a:', size(a,1), size(a,2)
! !             call print_c(a, 99, 'a', 'f4.1', achar(10))
!    else
!       rcnts = 0
!       sp => nula
!    end if

!          print *, 'c:', size(c,1), size(c,2), ciproc


   call mpi_alltoallw(c, scnts, sdispls, stypes, a, rcnts, rdispls, rtypes, comm, err)

   do i = 0, nproc-1
      call mpi_type_free(dars(i), err)
   end do


   end subroutine darray_allgather











!
!
!    subroutine twolvl1d2d_rdst(tbc, base, sz, data)
! !- 2 level restribution
! !-------------------------------------------------------------------------------
! !i Inputs:
! !i   tbc   : TB communicators and offsets for networking
! !i   base  : the base mpi type
! !i   sz    : size of the mpi_contiguous derived type in units of base
! !io  data : data to be reduced and assembled according to tbc % d2count and tbc % d2amap
! !
! !r Remarks:
! !r   First allreduce the local pieces across the kblock using the 1D communicators,
! !r   then allgather the pieces for different atoms within each kblock plane using the 2D communicators.
! !r   This shall be complete before the symmetrisation.
! !r       Beware mpi_sum is not defined for non-intrinsic types.
!
!       use tbprl
!       use mpi
!       use iso_c_binding, only : c_ptr
!
!       implicit none
!
!       type(tbc_t), intent(in) :: tbc
!       integer, intent(in) :: base, sz
!       type(c_ptr), intent(inout) :: data(*)
! !       real(8), intent(inout) :: rdata(*)
!
!
!       integer :: err, dtype
!
!       call mpi_allreduce(mpi_in_place, data, sz * tbc % d2count(tbc%c2d%id), base, mpi_sum, tbc % c1d % comm, err)
!       call mpi_type_contiguous(sz, base, dtype, err)
!       call mpi_type_commit(dtype, err)
!       call mpi_allgatherv(mpi_in_place, 0, dtype, data, tbc % d2count, tbc % d2amap, dtype, tbc % c2d % comm, err)
!       call mpi_type_free(dtype, err)
!
!    end subroutine twolvl1d2d_rdst
!
!    subroutine dbcast(a, c, desc)
!       implicit none
!       integer, parameter :: dp = 8
! !         rather limited for now only 0,0 can send
!       real(dp), intent(in) :: a(0:,0:)
!       real(dp), intent(out) :: c(0:,0:)
!       integer, intent(in) :: desc(9)
!
!       integer :: ctx, nprow, npcol, nb, ip, jp, il, jl, lr,lc, iprow, ipcol
!       integer :: i,j,n
!
!       ctx = desc(2)
!       n = desc(3)
!       nb = desc(5)
!       call fblacs_gridinfo(ctx, nprow, npcol, iprow, ipcol)
!       lr = ((n-1)/nb)/nprow+1
!       lc = ((n-1)/nb)/npcol+1
!
!       if (iprow == 0 .and. ipcol == 0) then
! !              generate the locals and   send out
!          do jp=0,npcol-1
!             do ip=0,nprow-1
!                if (ip /= iprow .or. jp /= ipcol) then
!                   c(0:lr*nb-1,0:lc*nb-1) = 0.0_dp
!                   do jl = 0, lc-1
!                      j = jl*npcol + jp
!                      if (j > (n-1)/nb) cycle
!                      do il = 0, lr-1
!                         i = il*nprow+ip
!                         if (i > (n-1)/nb) cycle
! !                                     print *, ip,jp, il,jl, i, j
!                         c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1) &
!                         & = a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1)
!                      end do
!                   end do
!                   call dgesd2d( ctx, lr*nb, lc*nb, c, lr*nb, ip, jp )
! !                             call print_c(c,200,'c'//i2s([ip,jp])//new_line(''),'f5.0')
!                end if
!             end do
!          end do
! !
! !        assemble the pieces that have to stay here
! !         this can be incorporated in the above by revesing the ip, jp counters e.g.: do ip = nprow-1,0,-1
!          c(0:lr*nb-1,0:lc*nb-1) = 0.0_dp
!          do jl = 0, lc-1
!             j = jl*npcol
!             if (j > (n-1)/nb) cycle
!             do il = 0, lr-1
!                i = il*nprow
!                if (i > (n-1)/nb) cycle
!                c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1) &
!                            & = a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1)
! !                 c(il*nb:(il+1)*nb-1,jl*nb:(jl+1)*nb-1) = a(i*nb:(i+1)*nb-1,j*nb:(j+1)*nb-1)
!             end do
!          end do
! !                 call print_c(c,200,'c'//i2s([0,0])//new_line(''),'f5.0')
!
!       else
!          call dgerv2d( ctx, lr*nb, lc*nb, c, lr*nb, 0, 0 )
!       end if
!
!    end subroutine dbcast
!
!
!
!
!
!
!
!
!    subroutine zbcast(a, c, desc)
!       implicit none
!       integer, parameter :: dp = 8
! !         rather limited for now only 0,0 can send
!       complex(dp), intent(in) :: a(0:,0:)
!       complex(dp), intent(out) :: c(0:,0:)
!       integer, intent(in) :: desc(9)
!
!       integer :: ctx, nprow, npcol, nb, ip, jp, il, jl, lr,lc, iprow, ipcol
!       integer :: i,j,n
!
!       ctx = desc(2)
!       n = desc(3)
!       nb = desc(5)
!       call fblacs_gridinfo(ctx, nprow, npcol, iprow, ipcol)
!       lr = ((n-1)/nb)/nprow+1
!       lc = ((n-1)/nb)/npcol+1
!
!       if (iprow == 0 .and. ipcol == 0) then
! !              generate the locals and   send out
!          do jp=0,npcol-1
!             do ip=0,nprow-1
!                if (ip /= iprow .or. jp /= ipcol) then
!                   c(0:lr*nb-1,0:lc*nb-1) = 0.0_dp
!                   do jl = 0, lc-1
!                      j = jl*npcol + jp
!                      if (j > (n-1)/nb) cycle
!                      do il = 0, lr-1
!                         i = il*nprow+ip
!                         if (i > (n-1)/nb) cycle
! !                                     print *, ip,jp, il,jl, i, j
!                         c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1) &
!                         & = a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1)
!                      end do
!                   end do
!                   call zgesd2d( ctx, lr*nb, lc*nb, c, lr*nb, ip, jp )
! !                             call print_c(c,200,'c'//i2s([ip,jp])//new_line(''),'f5.0')
!                end if
!             end do
!          end do
! !
! !        assemble the pieces that have to stay here
! !         this can be incorporated in the above by revesing the ip, jp counters e.g.: do ip = nprow-1,0,-1
!          c(0:lr*nb-1,0:lc*nb-1) = 0.0_dp
!          do jl = 0, lc-1
!             j = jl*npcol
!             if (j > (n-1)/nb) cycle
!             do il = 0, lr-1
!                i = il*nprow
!                if (i > (n-1)/nb) cycle
!                c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1) &
!                            & = a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1)
! !                 c(il*nb:(il+1)*nb-1,jl*nb:(jl+1)*nb-1) = a(i*nb:(i+1)*nb-1,j*nb:(j+1)*nb-1)
!             end do
!          end do
! !                 call print_c(c,200,'c'//i2s([0,0])//new_line(''),'f5.0')
!
!       else
!          call zgerv2d( ctx, lr*nb, lc*nb, c, lr*nb, 0, 0 )
!       end if
!
!    end subroutine zbcast
!
!
!
!    subroutine dwrite_home(a, c, desc)
!       implicit none
!       integer, parameter :: dp = 8
! !         rather limited for now only 0,0 can send
!       real(dp), intent(out) :: a(0:,0:)
!       real(dp), intent(in) :: c(0:,0:)
!       integer, intent(in) :: desc(9)
!
!       integer :: ctx, nprow, npcol, nb, ip, jp, il, jl, lr,lc,  iprow, ipcol
!       integer :: i,j,n
!
!       ctx = desc(2)
!       n = desc(3)
!       nb = desc(5)
!       call fblacs_gridinfo(ctx, nprow, npcol, iprow, ipcol)
!       lr = ((n-1)/nb)/nprow+1
!       lc = ((n-1)/nb)/npcol+1
!
!       if (iprow == 0 .and. ipcol == 0) then
!
!          a(0:lr*nb-1,0:lc*nb-1) = 0.0_dp
!
! !         this can be incorporated below by revesing the ip, jp counters e.g.: do ip = nprow-1,0,-1
!          do jl = 0, lc-1
!             j = jl*npcol
!             if (j > (n-1)/nb) cycle
!             do il = 0, lr-1
!                i = il*nprow
!                if (i > (n-1)/nb) cycle
!                a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1) &
!                      & = c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1)
! !                 c(il*nb:(il+1)*nb-1,jl*nb:(jl+1)*nb-1) = a(i*nb:(i+1)*nb-1,j*nb:(j+1)*nb-1)
!             end do
!          end do
! !                 call print_c(c,200,'c'//i2s([0,0])//new_line(''),'f5.0')
!
!
!          do jp=0,npcol-1
!             do ip=0,nprow-1
!                if (ip /= iprow .or. jp /= ipcol) then
!                   call dgerv2d( ctx, lr*nb, lc*nb, c, lr*nb, ip, jp )
!                   do jl = 0, lc-1
!                      j = jl*npcol + jp
!                      if (j > (n-1)/nb) cycle
!                      do il = 0, lr-1
!                         i = il*nprow+ip
!                         if (i > (n-1)/nb) cycle
!                         a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1) &
!                               & = c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1)
!                      end do
!                   end do
! !                             call print_c(c,200,'c'//i2s([ip,jp])//new_line(''),'f5.0')
!                end if
!             end do
!          end do
!
!
! !                 call print_c(a,300,'a'//new_line(''),'f5.0')
!
!       else
!          call dgesd2d( ctx, lr*nb, lc*nb, c, lr*nb, 0, 0 )
!       end if
!
!    end subroutine dwrite_home
!
!
!
!
!
!    subroutine zwrite_home(a, c, desc)
!       implicit none
!       integer, parameter :: dp = 8
! !         rather limited for now only 0,0 can send
!       complex(dp), intent(out) :: a(0:,0:)
!       complex(dp), intent(in) :: c(0:,0:)
!       integer, intent(in) :: desc(9)
!
!       integer :: ctx, nprow, npcol, nb, ip, jp, il, jl, lr,lc, iprow, ipcol
!       integer :: i,j,n
!
!       ctx = desc(2)
!       n = desc(3)
!       nb = desc(5)
!       call fblacs_gridinfo(ctx, nprow, npcol, iprow, ipcol)
!       lr = ((n-1)/nb)/nprow+1
!       lc = ((n-1)/nb)/npcol+1
!
! !       print *, 'lr, lc, nprow, npcol', lr, lc, nprow, npcol, iprow, ipcol
!
!
!       if (iprow == 0 .and. ipcol == 0) then
!
!          a(0:lr*nb-1,0:lc*nb-1) = 0.0_dp
!
! !         this can be incorporated below by revesing the ip, jp counters e.g.: do ip = nprow-1,0,-1
!          do jl = 0, lc-1
!             j = jl*npcol
!             if (j > (n-1)/nb) cycle
!             do il = 0, lr-1
!                i = il*nprow
!                if (i > (n-1)/nb) cycle
!                a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1) &
!                      & = c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1)
! !                 c(il*nb:(il+1)*nb-1,jl*nb:(jl+1)*nb-1) = a(i*nb:(i+1)*nb-1,j*nb:(j+1)*nb-1)
!             end do
!          end do
! !                 call print_c(c,200,'c'//i2s([0,0])//new_line(''),'f5.0')
!
!
!          do jp=0,npcol-1
!             do ip=0,nprow-1
!                if (ip /= iprow .or. jp /= ipcol) then
!                   call zgerv2d( ctx, lr*nb, lc*nb, c, lr*nb, ip, jp )
!                   do jl = 0, lc-1
!                      j = jl*npcol + jp
!                      if (j > (n-1)/nb) cycle
!                      do il = 0, lr-1
!                         i = il*nprow+ip
!                         if (i > (n-1)/nb) cycle
!                         a(i*nb:min((i+1)*nb,n)-1,j*nb:min((j+1)*nb,n)-1) &
!                               & = c(il*nb:il*nb+min((i+1)*nb,n)-i*nb-1, jl*nb:jl*nb+min((j+1)*nb,n)-j*nb-1)
!                      end do
!                   end do
! !                             call print_c(c,200,'c'//i2s([ip,jp])//new_line(''),'f5.0')
!                end if
!             end do
!          end do
!
!
! !                 call print_c(a,300,'a'//new_line(''),'f5.0')
!
!       else
!          call zgesd2d( ctx, lr*nb, lc*nb, c, lr*nb, 0, 0 )
!       end if
!
!    end subroutine zwrite_home
!
!
!
!    subroutine mkapidx(n, ipc, lmxl, lxl, idx)
! !    this exists just to get around 'the w'
!
!       integer, intent(in) :: n, ipc(n), lmxl(*)
!       logical, intent(in) :: lxl
!       integer, intent(out) :: idx(n)
! !       xl == false : B_{R,L;R',L'}
! !       xl == true : B_{R,L+1;R',L'}
!
!       integer :: i, xl
!
!       xl = 1
!       if (lxl) xl = xl+1
!
!       idx(1) = (lmxl(ipc(1))+xl)**2
!       do i = 2, n
!          idx(i) = idx(i-1) + (lmxl(ipc(i))+xl)**2
!       end do
!
!    end subroutine mkapidx
!
