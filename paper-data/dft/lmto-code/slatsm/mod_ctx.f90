!> Handle communicators without calling mpi inquiry functions all the time,
!> or passing loads of extra arguments through all functions.
!! DMT
   module mod_ctx

   use mpi

   implicit none

   private

   integer, parameter :: maxdims = 3


   type :: cart_ctx_t
      integer :: id, sz, comm, rt, nd
      integer, dimension(0:maxdims-1) :: crd, szs
      logical :: lrt, circ(0:maxdims-1)
   end type

   interface init_ctx
      module procedure init_cart_ctx_from_comm
   end interface init_ctx

   public :: cart_ctx_t, init_ctx, sub_cart_ctx, del_cart_ctx, cart_ctx_fill, inplace_selection_sort, init_cart_ctx_0

   contains

!> Auxiliary routine.
!> Sort "a" inplace by successively swapping the elements with the smallest left.
!> This is a badly scaling algorithm but shall not matter because the
!> length is not expected to be more than what one can count on one hand.
   subroutine inplace_selection_sort(n, a, d)
      integer, intent(in) :: n !< number of elements to sort
      integer, intent(inout) :: a(n) !< array to be sorted in place
      integer, intent(in), optional :: d !< order/direction: (default) > 0 ascending order, < 0 descending order

      integer :: i, t, p(1)


      do i = 1, n-1
         p = minloc(a(i:n)) + i - 1
         if (p(1) /= i) then
            t = a(i)
            a(i) = a(p(1))
            a(p(1)) = t
         end if
      end do


!       procedure pointer to the intrinsics minloc, maxloc does not work
      if (present(d)) then
         if (d < 0) then
            do i = 1, n/2
               t = a(i)
               a(i) = a(n-i+1)
               a(n-i+1) = t
            end do
         end if
      end if


   end subroutine inplace_selection_sort


!> Copy the properties of the cartesian communicator comm to the context ctx.
   subroutine cart_ctx_fill(ctx, comm, root)
      type(cart_ctx_t), intent(inout) :: ctx !< Context, the fields of which are to be filled
      integer, intent(in) :: comm !< mpi communicator with cartesian topology.
      integer, intent(in), optional :: root !< root id, not part of the mpi specification

      integer :: err

      ctx % comm = comm

      call mpi_comm_rank(ctx % comm, ctx % id, err)
      call mpi_comm_size(ctx % comm, ctx % sz, err)

      ctx % rt = 0
      if (present(root)) ctx % rt = root

      ctx % lrt = ctx % rt == ctx % id

      call mpi_cartdim_get(ctx % comm, ctx % nd, err)
      if (ctx % nd > maxdims) call rx("Error: cart_ctx % nd > maxdims")

      call mpi_cart_get(ctx % comm, ctx % nd, ctx % szs, ctx % circ, ctx % crd, err)


   end subroutine cart_ctx_fill



   subroutine init_cart_ctx_from_comm(c, nd, comm, root, sizes, circular, reorder, dims_order)

      type(cart_ctx_t), intent(inout) :: c
      integer, intent(in) :: nd
      integer, intent(in), optional :: comm, root, sizes(0:nd-1), dims_order(0:nd-1)
      logical, intent(in), optional :: circular(0:nd-1), reorder

      integer :: err
      integer :: comm_loc, root_loc, ccomm, rank_old, root_new

      logical :: reorder_loc

      if (nd > maxdims) call rx("init_cart_ctx_from_comm: nd > maxdims")

      comm_loc = mpi_comm_world
      if (present(comm)) comm_loc = comm

      call mpi_comm_size(comm_loc, c % sz, err)

      c % szs = 0
      if (present(sizes)) c % szs(0:nd-1) = sizes(0:nd-1)

      c % nd = nd
!       print *, c % sz, c % nd, c % szs
      call mpi_dims_create(c % sz, c % nd, c % szs, err)
      if (err /= 0) call rx('Inappropriate dimensions passed to mpi_dims_create')
      if (present(dims_order)) then
         call inplace_selection_sort(nd, c % szs)
         c % szs = c % szs ( dims_order )
      end if

      reorder_loc = .true.
      if (present(reorder)) reorder_loc = reorder

      c % circ = .false.
      if (present(circular)) c % circ = circular

      call mpi_cart_create(comm_loc, nd, c % szs, c % circ, reorder_loc, ccomm, err)

      root_loc = 0
      if (present(root)) root_loc = root

      root_new = root_loc
      if (reorder_loc) then
         call mpi_comm_rank(comm_loc, rank_old, err)
         if (rank_old == root_loc) call mpi_comm_rank(ccomm, root_new, err )
         call mpi_bcast(root_new, 1, mpi_integer, root_loc, comm_loc, err)
      end if

      call cart_ctx_fill(c, ccomm, root_new)

      if (c % lrt) print *, 'Process grid dims: ', c % szs(0:nd-1)
! !#ifndef LINUXF
!       flush(6); call mpi_barrier(comm, err)
!       print  "(' Process id: ', i5.5, '/', i0, ', coords:', 3(x,i5.5))", c % id, c % sz, c % crd(0:nd-1)
!       flush(6); call mpi_barrier(comm, err)
! !#endif
   end subroutine init_cart_ctx_from_comm


!> Partition the cartesian context c preserving the specified dimension
!> remain_dims resulting in multitude of contexts of lower dimensionality.
!> Each process obtains its respective new context.
!> The root processes of the resulting contexts are set appropriately using the parent context's root.
!> See mpi_cart_sub.
   function sub_cart_ctx(c, remain_dims) result (r)
      type(cart_ctx_t), intent(in) :: c !< Source context, serving as parent.
      logical, intent(in) :: remain_dims(0:) !< Dimensions set to .true. are preserved.
      type(cart_ctx_t) :: r !< result

      integer :: lcomm, err, root_crds(0:maxdims-1), new_root, i, idx(0:count(remain_dims)-1)

      call mpi_cart_sub(c % comm, remain_dims, lcomm, err)

      idx = pack([(i, i=0, c % nd-1)], remain_dims)
      call mpi_cart_coords(c % comm, c % rt, c % nd, root_crds, err)
      call mpi_cart_rank(lcomm, root_crds(idx), new_root, err)
      r % nd = size(idx)
      call cart_ctx_fill(r, lcomm, new_root)

   end function sub_cart_ctx

!> Finalizer (destructor), not yet used syntactically as such.
   subroutine del_cart_ctx(c)
      type(cart_ctx_t), intent(inout) :: c !< context to be cleared
      integer :: err

      call mpi_comm_free(c % comm, err)

      c % id   = -1
      c % sz   = -1
      c % rt   = 0
      c % nd   = 0
      c % crd  = -1
      c % szs  = 0
      c % lrt  = .false.
      c % circ = .false.

   end subroutine del_cart_ctx


   subroutine init_cart_ctx_0(c, nd)
      integer, intent(in) :: nd
      type(cart_ctx_t), intent(out) :: c

      if (nd > maxdims) call rx("init_cart_ctx_0: nd > maxdims")

      c % comm = -1

      c % id   = 0
      c % sz   = 1
      c % rt   = 0
      c % nd   = nd
      c % crd  = 0
      c % szs  = 1
      c % lrt  = .true.
      c % circ = .false.

   end subroutine init_cart_ctx_0


   end module mod_ctx
