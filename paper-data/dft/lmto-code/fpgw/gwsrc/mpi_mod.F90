      module mpi_mod

! Purpose of the module:
!
! Encapsulate all MPI book keeping, and provide corresponding dummies if
! compiled without MPI support
!
! Besides the top level (rank, size, MPI_COMM_WORLD, etc.) the module
! provides means to split communicators for nested parallelizations.


#ifdef USE_MPI
      include 'mpif.h'
#else
        integer MPI_comm_world
#endif


      type comm_struct

      integer :: size
      integer :: ID
      integer :: communicator
      integer :: group_ID
      integer :: num_of_siblings
      integer :: ID_of_head

      end type comm_struct

!.......................................................................


! Structures to keep account of different groups and communicators:
!
! Global data:
!
!
! Each processor needs to know it's id with respect to each (relevant) communicator
! (relevant means that the processor is part of the communicator)


      integer :: ierror

      integer :: mpi_size = 1,      mpi_rank = 0

      type(comm_struct) comm_world

! provide comm structures for 4 levels:

      type(comm_struct) comm_1,  comm_2,  comm_3,  comm_4
      type(comm_struct) inter_1, inter_2, inter_3, inter_4

! the lower level communicators don't know about it's siblings
! we also need inter-group communicators which connect the heads of
! each group at a given level

!---------------------------------------------------------------------------
! level 0:  (COMM_WORLD)  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | ...
!---------------------------------------------------------------------------
! level 1:  comm_1        |    group 0    |    group 1    |  group 2  | ..
!  ID                     | 0 | 1 | 2 | 3 | 0 | 1 | 2 | 3 | 0 | 1 | ...
!---------------------------------------------------------------------------
!           inter_1       | 0 |           | 1 |           | 2 |
!                             | 0 |           | 1 |           | 2 |
!                                 | 0 |           | 1 |           | 2 |
!                                     | 0 |           | 1 |           | 2 |
!---------------------------------------------------------------------------
! level 2:  comm_2 (1)    | gr 0  |  gr 1 |
!  ID                     | 0 | 1 | 0 | 1 |
!---------------------------------------------------------------------------
!           inter_2(1)    | 0 |   | 1 |
!                             | 0 |   | 1 |
!---------------------------------------------------------------------------
! level 2:  comm_2 (2)                    | gr 0  |  gr 1 |
!  ID                                     | 0 | 1 | 0 | 1 |
!---------------------------------------------------------------------------
!           inter_2(2)                    | 0 |   | 1 |
!                                             | 0 |   | 1 |
!---------------------------------------------------------------------------
! level 2:  comm_2 (2)                                    | gr 0  |  gr 1 |
!  ID                                                     | 0 | 1 | 0 | 1 |
!---------------------------------------------------------------------------
!           inter_2(3)                                    | 0 |   | 1 |
!                                                             | 0 |   | 1 |
!---------------------------------------------------------------------------
! ...       ...
!
!

! global storage of number of processors/groups for each level


      integer mpi_N1, mpi_N2, mpi_N3

!  integer, allocatable:: proclist(:,:)
!  real(kind(0.d0)),allocatable:: loadlist(:)

!.......................................................................

! some parametgers for BLACS / Scalapack

      integer DLEN_
      PARAMETER          ( DLEN_ = 9 )


      contains

      subroutine mpi_start

      implicit none

      mpi_rank=0
      mpi_size=1

      ierror = 0


#ifdef USE_MPI
      call MPI_Init(ierror)
      call MPI_Comm_Rank(MPI_Comm_World,mpi_rank,ierror)
      call MPI_Comm_Size(MPI_Comm_World,mpi_size,ierror)
#else
        MPI_Comm_World = 0
#endif

      comm_world%size            = mpi_size
      comm_world%ID              = mpi_rank
      comm_world%communicator    = MPI_Comm_World
      comm_world%group_ID        = 0
      comm_world%num_of_siblings = 1
      comm_world%ID_of_head      = 0
      comm_1%size            = 1
      comm_1%ID              = 0
      comm_1%communicator    = 0
      comm_1%group_ID        = 0
      comm_1%num_of_siblings = 1
      comm_1%ID_of_head    = 0

      comm_2%size            = 1
      comm_2%ID              = 0
      comm_2%communicator    = 0
      comm_2%group_ID        = 0
      comm_2%num_of_siblings = 1
      comm_2%ID_of_head      = 0

      comm_3%size            = 1
      comm_3%ID              = 0
      comm_3%communicator    = 0
      comm_3%group_ID        = 0
      comm_3%num_of_siblings = 1
      comm_3%ID_of_head      = 0

      comm_4%size            = 1
      comm_4%ID              = 0
      comm_4%communicator    = 0
      comm_4%group_ID        = 0
      comm_4%num_of_siblings = 1
      comm_4%ID_of_head      = 0

      end subroutine mpi_start

!----------------------------------------------------------------------------

      subroutine mpi_end

      implicit none

#ifdef USE_MPI
        call MPI_BARRIER  (MPI_COMM_WORLD, ierror)
        call MPI_FINALIZE (ierror)
#endif

      end subroutine mpi_end

!----------------------------------------------------------------------------

      function mpi_split_comm(base_comm, processor_list, new_comm, inter_comm)

      implicit none

      type(comm_struct), intent(IN)  :: base_comm
        integer ,          intent(IN)  :: processor_list(0:base_comm%size)
!  type(comm_struct), intent(OUT) :: new_comm, inter_comm
        type(comm_struct) :: new_comm, inter_comm

      integer mpi_split_comm

        integer color, ID, size

      color = processor_list(base_comm%ID) + 1


        mpi_split_comm = 0

#ifdef USE_MPI


        call MPI_COMM_Split (base_comm%communicator,color,0,new_comm%communicator,ierror)
!  call MPI_COMM_Split (base_comm%communicator,color,mpi_rank,new_comm%communicator,ierror)

        call MPI_COMM_RANK  (new_comm%communicator, ID,   ierror)
        call MPI_COMM_SIZE  (new_comm%communicator, size, ierror)

      new_comm%group_ID = color - 1
        new_comm%ID       = ID
        new_comm%size     = size

      mpi_split_comm = ierror

        write(*,*) "MPI-split: ", mpi_rank, ID, color-1

!  Create the inter-group communicator:

!  if (new_comm%ID == 0) then
!     color = 1
!  else
!     color = MPI_UNDEFINED
!  endif

      color = new_comm%ID  + 1

        call MPI_COMM_Split(base_comm%communicator,color,0,inter_comm%communicator,ierror)
!  call MPI_COMM_Split(base_comm%communicator,color,mpi_rank,inter_comm%communicator,ierror)

        if (inter_comm%communicator /= MPI_COMM_NULL) then
           call MPI_COMM_RANK(inter_comm%communicator,inter_comm%ID,ierror)
           call MPI_COMM_SIZE(inter_comm%communicator,inter_comm%size,ierror)
           inter_comm%group_ID = color - 1
        else
           inter_comm%size = 1
           inter_comm%ID   = 0
           inter_comm%ID   = MPI_UNDEFINED
        end if

      mpi_split_comm = ierror

#else
        ierror = 0
#endif

      return

      end function



!----------------------------------------------------------------------------

      subroutine mpi_mk_groups (base_size, loop_count, num_groups, processor_list)
!
! Generate groups for outer loops:
!
! processor_list(base_size-1)+1 is the size of the smallest group
!

      integer, intent(IN)  :: base_size, loop_count
        integer :: num, i
      integer, intent(OUT) :: processor_list(0:base_size-1)



      if(num_groups == 0) then

!        find largest number of outer loop groups:

           if(base_size > 1) then
        num = CEILING( loop_count / (1.0*base_size) )
           else
              num = 1
           endif

        num_groups = loop_count / num

        do while ( (mod(base_size, num_groups) /= 0) .or. (mod(loop_count, num_groups) /= 0))
            num_groups = num_groups - 1
        end do

      elseif(num_groups > base_size) then

        num_groups = base_size

      endif


      do i=0, base_size-1
        processor_list(i) = mod(i,num_groups)
      end do

      end subroutine mpi_mk_groups

!----------------------------------------------------------------------------

      subroutine mpi_make_groups (base_comm, loop_count, num_groups, new_comm, inter_comm)

      type(comm_struct), intent(IN)  :: base_comm
      integer,           intent(IN)  :: loop_count
      integer,           intent(INOUT) :: num_groups
      type(comm_struct), intent(OUT) :: new_comm, inter_comm

      integer                        :: processor_list(0:base_comm%size-1)
      integer                        :: ierr

      call mpi_mk_groups(base_comm%size, loop_count, num_groups, processor_list)

      ierr = mpi_split_comm(base_comm, processor_list, new_comm, inter_comm)

      new_comm%num_of_siblings = num_groups

      inter_comm%num_of_siblings = 1


#ifdef USE_MPI
      write(*,'( "MPI_make_groups: rank = ",I4,", ID = ",I4,", group_ID = ",I4 )') mpi_rank, new_comm%ID, new_comm%group_ID
#endif
      if (ierr /= 0) call rx('Problem creating MPI groups')

      end subroutine mpi_make_groups

!----------------------------------------------------------------------------

      subroutine find_divisors (size, nx, ny)

      integer, intent(in)  :: size
      integer, intent(out) :: nx, ny


      nx = INT(SQRT(REAL(size)))

      do while ( mod(size,nx) /= 0 )
        nx = nx - 1
      enddo

      ny = size/nx

      end subroutine find_divisors



      end module mpi_mod

