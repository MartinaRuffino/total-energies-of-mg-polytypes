

subroutine comm_phys_split(comm, new_comm)
!    Split the input communicator to smaller nonoverlapping communicators
! each containing the processes of one whole node.
   use mpi
   implicit none
   integer, intent(in) :: comm
   integer, intent(out) :: new_comm
   integer :: root, rank, comm_size, maxnamelen, color, ierr, i, j
   character(len=mpi_max_processor_name), allocatable :: names(:)
   logical :: lroot, match

! //    this could probably be done better with hwloc ...
   root = 0;
   call mpi_comm_rank(comm, rank, ierr);

   lroot = rank == root;

   call mpi_comm_size(comm, comm_size, ierr);

   allocate(names(0:comm_size-1))

   do i = 0, comm_size-1
      do j = 1, mpi_max_processor_name
         names(i)(j:j) = ' '
      end do
   end do

   call mpi_get_processor_name(names(rank), maxnamelen, ierr);
   call mpi_allgather(mpi_in_place, 0, mpi_datatype_null, names, mpi_max_processor_name, mpi_character, comm, ierr);

   color = 0;
   match = .false.
   do while (color < comm_size .and. .not. match)
      match = names(rank)(1:mpi_max_processor_name) == names(color)(1:mpi_max_processor_name)
      color = color + 1
   end do
   deallocate(names)

   call mpi_comm_split(comm, color, rank, new_comm, ierr)
end subroutine comm_phys_split


subroutine loc_dev_init(ierr)
! Initialise the CUDA devices and associate each process with a local device without overloading.
! Here we rely on all devices being of acceptable capability. This is ensured by the shim library: devremapshim.c

   use mpi
   use fcuda
   use iso_c_binding, only : c_null_char
   implicit none

   integer, intent(out) :: ierr
   integer :: comm, rank, csize, total_ng, usable_ng, namelen, devid, grank, gsize, &
            & major, minor, totalGlobalMem, multiProcessorCount, computeMode, devnamelen
   logical :: lroot
   character(len=mpi_max_processor_name) :: name
   character(len=128) :: devname
   integer :: i

   call comm_phys_split(mpi_comm_world, comm)

   call mpi_comm_rank(comm, rank, ierr)

   call mpi_comm_size(comm, csize, ierr)

   call mpi_comm_free(comm, ierr)

   lroot = rank == 0

!    usable_ng = 0
   call fcuda_get_device_count (total_ng)



! This is now handled by the cudaDevRemap shim so we can call fcuda_set_device(rank) directly.
!    do i = 0, total_ng-1
!       if (1 == device_is_usable(i))  then
!          if (rank == usable_ng) call fcuda_set_device(i);
!          usable_ng = usable_ng + 1
!       end if
!    end do

   if (total_ng > 0) then
      call fcuda_set_device(mod(rank,total_ng))
      call mpi_get_processor_name(name, namelen, ierr);

      call fcuda_get_local_dev_props(devid, devname, major, minor, totalGlobalMem, multiProcessorCount, computeMode)
      devnamelen = index(devname, c_null_char)-1

      call mpi_comm_rank(mpi_comm_world, grank, ierr);
      call mpi_comm_size(mpi_comm_world, gsize, ierr);

      write(*, "('proc global: ',i0,'/',i0,', local: ',i0,'/',i0,' @ ',a, &
               & ': dev (',i0,') ',a,', v',i0,'.',i0,', ',i0,' MiB, ',i0,' SMs, mode: ',i0)") &
               & grank, gsize, rank, csize, name(1:namelen), devid, devname(1:devnamelen), major, &
               & minor, totalGlobalMem, multiProcessorCount, computeMode
   end if

end subroutine loc_dev_init



subroutine loc_dev_init_insane()
! Apparantly this has to be called before mpi_init
! There is a shoddy way to get the rank of a process but communication is still necessary afterwards
! This defies sanity so I am guessing there must be a good reason for it, and at the same time hoping it gets 'fixed'.

   use fcuda
   use iso_c_binding, only : c_null_char
   implicit none

   character(len=10) :: val
   integer :: rank, csize, usable_ng, namelen, devid, grank, gsize, &
            & major, minor, totalGlobalMem, multiProcessorCount, computeMode, devnamelen
   character(len=100) :: devname

   call get_environment_variable('MV2_COMM_WORLD_LOCAL_RANK',val)

!    if (len_trim(val) > 0)

   read(val, *) rank


   if (rank == 0) devid = 0
   if (rank == 1) devid = 2
   if (rank > 1) stop 'loc_dev_init_insane: comm_size hardcoded to 2 due to a bulshit MPI-CUDA library design'
   grank = rank
   csize = 2
   gsize = csize

   call fcuda_set_device(devid);

   call fcuda_get_local_dev_props(devid, devname, major, minor, totalGlobalMem, multiProcessorCount, computeMode)
   devnamelen = index(devname, c_null_char)-1
   write(*, "('proc global: ',i0,'/',i0,', local: ',i0,'/',i0,' @ ',a, &
            & ': dev (',i0,') ',a,', v',i0,'.',i0,', ',i0,' MiB, ',i0,' SMs, mode: ',i0)") &
            & grank, gsize, rank, csize, ' ', devid, devname(1:devnamelen), major, &
            & minor, totalGlobalMem, multiProcessorCount, computeMode

end subroutine loc_dev_init_insane

