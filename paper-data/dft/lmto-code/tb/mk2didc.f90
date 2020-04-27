   subroutine mk2didc(tbc, nbas, iax, niax, nsites, npr, ipc, idxdn, nl)
! 1. Distribute the atoms over the processes of the 2D communicators in contiguous blocks ballancing the approximate effort.
! 2. Find process local clusters or atoms to help tbfrce.f90. This is fairly experimental at this point.

      use tbprl
      use mpi

      implicit none

      type(tbc_t), intent(inout), target :: tbc
      integer, intent(in) :: nbas, nl, niax, nsites, iax(0:niax-1,nsites), npr(0:1,nbas), ipc(nbas), idxdn(0:nl-1,*)



      integer :: i, j, ij, l, k, pid, err
      integer, allocatable :: pnlm(:), pidx(:)
      logical, allocatable :: walked(:)
      integer, pointer :: amap(:), neighmap(:), neighidx(:)

      procedure(logical) :: cmdopt
      character(len=128) :: outs

      call tcn('mk 2D idc')

      amap => tbc % d2amap


      allocate(walked(nbas))

! find the packed offsets

      if (.not. allocated(tbc % neighmap)) allocate(tbc % neighmap(nbas+1))
      neighmap => tbc % neighmap

      neighmap(1) = 0
      do i = 1, nbas
         neighmap(i+1) = neighmap(i)
         walked = .false.
         do ij = 1, npr(0, i)
            j = iax(1,npr(1, i)+ij)
            if (walked(j)) cycle
            neighmap(i+1) = neighmap(i+1) + 1
            walked(j) = .true.
         end do
      end do


! set the indexee array
      if (allocated(tbc % neighidx)) deallocate(tbc % neighidx)
      allocate(tbc % neighidx(neighmap(nbas+1)))
      neighidx => tbc % neighidx
      k = 1
      do i = 1, nbas
         walked = .false.
         do ij = 1, npr(0, i)
            j = iax(1,npr(1, i)+ij)
            if (walked(j)) cycle
            neighidx(k) = j
            k = k + 1
            walked(j) = .true.
         end do
      end do
      deallocate(walked)
! ! test
!       do i = 1, nbas
!          do ij = neighmap(i)+1, neighmap(i+1)
!             j = neighidx(ij)
!             write(645, *) i,j
!          end do
!       end do
!       call rx('stop')

      if (.not. cmdopt('--leanrho',9,0,outs)) then
! find number of orbitals per atom
          allocate(pnlm(nbas))
          do i = 1, nbas
             pnlm(i) = 0
             do l = 0, nl-1
                 if (idxdn(l,ipc(i))<=1) pnlm(i) = pnlm(i) + 2*l + 1
             end do
          end do

! approximate the effort for each atom
          allocate(pidx(0:nbas))
          pidx(0) = 0
          do i = 1, nbas
             pidx(i) = pidx(i-1)
             do ij = neighmap(i)+1, neighmap(i+1)
                j = neighidx(ij)
                pidx(i) = pidx(i) + pnlm(j)
             end do
          end do
          deallocate(pnlm)

! distribute the atoms over the processes in contiguous blocks ballancing the approximate effort
          call vbdist(nbas, pidx(1:nbas), tbc % c2d % sz, amap)
          deallocate(pidx)
      else
          amap(0) = 0
          amap(1 : tbc % c2d % sz) = nbas
      endif

      tbc % d2count = amap(1:tbc % c2d % sz) - amap(0:tbc % c2d % sz-1)


! !       pid = tbc % c2d % id
!
! ! find process local clusters with fast uniqsort for indexes
! !  1) set the available indexes and count the unique elements
! !  2) allocate the obtained number and pack them in
!       allocate(walked(nbas))
!       if (.not. allocated(tbc%nlclus)) allocate(tbc%nlclus(0 : tbc % c2d % sz - 1))
!       if (.not. allocated(tbc%lclus )) allocate(tbc%lclus(nbas, 0 : tbc % c2d % sz - 1))
!
!       do pid = 0, tbc % c2d % sz - 1
!          walked = .false.
!          tbc % nlclus(pid) = 0
!          do ij = neighmap(amap(pid)+1)+1, neighmap(amap(pid+1)+1)
!             j = neighidx(ij)
!             if (.not. walked(j)) then
!                tbc % nlclus(pid) = tbc % nlclus(pid) + 1
!                walked(j) = .true.
!             end if
!          end do
!
! !          allocate( tbc % lclus( tbc % nlclus ) )
!          tbc % lclus(1:nbas,pid) = -9999
!          k = 1
!          do i = 1, nbas
!             if (walked(i)) then
!                tbc % lclus(k,pid) = i
!                k = k + 1
!             end if
!          end do
!       end do
!
!       deallocate(walked)
!
!
!
!
! !
! ! !  Find the necessary distributed blocks to receive
! ! !  [rmprcr, rmblcr, lcprcr, lcblcr]
! !       allocate(pidx(nbas+1))
! !       pidx(1) = 0
! !       do i = 1, nbas
! !          pidx(i+1) = pidx(i) + pnlm(i)
! !       end do
! !
! !
! !       do ij = 1, nlclus
! !          i = lclus(ij)
! !          bi = pidx(i)
! !          ei = pidx(i+1)-1
! !
! !
! !       end do
! !
! ! !
! ! ! ! walk over blocks of the 2D blockscyclic distributed array and record the ones required locally
! ! ! ! 'required' here means they contain atom of the local cluster.
! ! ! ! This is specifically made for the transposed eigenvectors so the walk is over the columns only.
! ! !
! ! !    do bj = 0, tbc % blocks(2)
! ! !
! ! !    end do
! !
! !
!
!
! !
! !       write(700+pid, *) 'nlclus:', nlclus
! !       write(700+pid, *) 'amap:  ', amap
! !       write(700+pid, *) 'neimap:', neighmap
! !       write(700+pid, *) 'neiidx:', neighidx(neighmap(amap(pid)+1)+1 : neighmap(amap(pid+1)+1))
! !       write(700+pid, *) 'lclus: ', lclus
! !
! !       flush(700+pid)
! !
! !       call mpi_barrier(tbc%c2d%comm,err)
! !
! !       call rx('stop')

      call tcx('mk 2D idc')

   end subroutine mk2didc


