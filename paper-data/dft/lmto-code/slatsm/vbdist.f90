   subroutine vbdist(n, idx, np, amap)
!- Distribute n variable size elements in np boxes as evenly as you can.
!-------------------------------------------------------------------------------
!i Inputs:
!i   n    : number of elements
!i   idx  : cumulative array , ie. idx(1) == size(1)
!i                                 idx(i) == idx(i-1) + size(i)
!i                                 idx(n) == the sum of the sizes of all elements
!i   np    : number of boxes (processes in the currently intended usage)
!
!o Output:
!o   amap  : cumulative array of evenly mapped elements
!o           box i is assigned the elements from (amap(i)+1) to (amap(i+1)) inclusive
!
!ir Remarks
!ir This may benefit from some improvements and sophistication.


      integer, intent(in) :: n, idx(n), np
      integer, intent(out) :: amap(0:np)

!       extra +1 for forces..

      real :: av, tg, dfp, df !, dfs(0:np), tgs(np)
      integer :: a, p

      av = real(idx(n))/real(np)

!       dfs(0) = 0
      amap(0) = 0
      a = 1
!       print *, 'p,a,idx,tg:', amap(0),0,0,0
      do p = 1, np
         tg = p*av

         dfp = huge(0.0)
         df = abs(tg - idx(a))
!          print *, 'p,tg,a,idx,dfp:',p,tg,a,idx(a),dfp
         do while (df <= dfp .and. a<n)
            dfp = df
            a = a + 1
            df = abs(tg - idx(a))
!             print *, '     a,idx,df:',a,idx(a),dfp,df
         end do
         if (a < n .or. df > dfp) a = a-1
!          if (a == amap(p-1)) then
!             if (a == n)
!          end if
         amap(p) = a
!          dfs(p) = abs(tg - idx(a))
!          print *, 'p,a,idx,tg:', p,a,idx(a),tg
      end do

!       tgs = (/(p*av, p=1,np)/)

!       print *, 'idx:             ', idx
!       print *, 'amap:', amap
!       print '(" dfs:    ",5(8x,f4.2))', dfs
!       print '(" tgs:             ",4(6x,f6.2))', tgs
!       print *, 'rdx:             ', idx(amap(1:))
! !       print *, 'asum1:', sum(abs(idx(amap(1:)) - tgs)), sum(dfs)
! !       print *, 'asum2:', sum(abs(idx([1,2]) - tgs))

   end subroutine vbdist


   module mod_adstrb
      use mpi
      implicit none
      private

!       pointer list type
      type pl_t
         real(8), pointer :: r(:,:)
      end type

      public :: pl_t, adstrb

   contains

      subroutine adstrb(na, pl, np, aip, ip, comm)
! Distribute atomic quantities held in pointer lists.
! For example the augmentation matrices in s_site(:). To use this, initialise an
! array (size na) of instances of pl_t and associate the member pointer "r" of each
! element to the respective target then call this routine. See locpot.f for examples
!   :
! use mod_adstrb
! implicit none
!   :
! type(pl_t) :: pl(0:nbas-1)
!   :
! do ib = 0, nbas - 1; pl(ib) % r => s_site(ib+1) % pikk ;  end do
! call adstrb(nbas, pl, nproc, pmap, iproc, comm)
! ...
!
!     na: number of atoms
!     pl: pointer list
!     np: number of processes
!     aip: atoms on process displacements
!     ip: process index
!     comm: communicator
!
!     Naming scheme for the 3 letter index variables:
!
!     elems: e
!     atoms: a
!     procs: p
!
!     sizes/counts: s
!     integral sizes/counts: i
!
!     thus 'asp' for example means atoms on process (count)
!          'aip' would then be atoms on process integrated (also can be seen as displacements)
!
! This routine assembles disjoint pieces poined to by the pl pointer list, then
! uses mpi_allgatherv and scrambles the pack back to the pointed data.
!
! It is tempting to simply use mpi_gatherv with buffer 0x0 and displacements
! table made of the pointers to the disjoint pieces but the standard limits the
! diplacements type to 32bit int.
!

      integer, intent(in) :: na, np, ip
      integer, intent(in) :: comm
      integer, intent(in) :: aip(0:np)

      type(pl_t), intent(inout) :: pl(0:na-1)

      integer :: ia, err, rsz

      integer :: eia(0:na), esa(0:na-1), &
                 eip(0:np), esp(0:np-1)

      real(8), allocatable :: buf(:)

      eia(0) = 0
      do ia = 1, na
        rsz = 0
        if (associated(pl(ia-1) % r)) rsz = size(pl(ia-1) % r)
        eia(ia) = eia(ia-1) + rsz
      end do

      esa = eia(1:na) - eia(0:na-1)

      allocate(buf(0:eia(na)-1))

      do ia = aip(ip), aip(ip+1)-1
         if (associated(pl(ia) % r)) call dcopy(esa(ia), pl(ia) % r, 1, buf(eia(ia)), 1)
      end do

      eip = eia(aip)
      esp = eip(1:np) - eip(0:np-1)

      call mpi_allgatherv(mpi_in_place, esp(ip), mpi_real8, buf, esp, eip, mpi_real8, comm, err)

      do ia = 0, na - 1
        if (associated(pl(ia) % r)) call dcopy(esa(ia), buf(eia(ia)), 1, pl(ia) % r, 1)
      end do

      deallocate(buf)

      end subroutine adstrb

   end module mod_adstrb

