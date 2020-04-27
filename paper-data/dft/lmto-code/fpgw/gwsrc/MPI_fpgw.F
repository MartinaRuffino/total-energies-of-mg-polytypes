C      module m_mpi
C      implicit none
C      include "mpif.h"
C
C      integer :: mpi__info
C      integer :: mpi__size
C      integer :: mpi__rank
C      logical :: mpi__root
C      integer :: mpi__iini, mpi__iend
C      logical, allocatable ::  mpi__task(:)
C      integer,allocatable :: mpi__ranktab(:)
C      integer:: ista(MPI_STATUS_SIZE )
C
C      interface MPI__Send
C      module procedure
C     &  MPI__Send_i,  MPI__Send_iv,
C     &  MPI__Send_d,  MPI__Send_dv
C      end interface
C
C      interface MPI__Recv
C      module procedure
C     &  MPI__Recv_i,  MPI__Recv_iv,
C     &  MPI__Recv_d,  MPI__Recv_dv
C      end interface
C
C      contains
C
C
C      subroutine MPI__Initialize()
C      implicit none
C      character(1024*4) :: cwd, stdout
C      call getcwd(cwd)          ! get current working directory
C
C      call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
C      call MPI_Comm_rank( MPI_COMM_WORLD, mpi__rank, mpi__info )
C      call MPI_Comm_size( MPI_COMM_WORLD, mpi__size, mpi__info )
C
C      if( mpi__rank == 0 ) then
C         mpi__root = .true.
C      else
C         mpi__root = .false.
C      end if
C
C      if( mpi__root ) then
C         call chdir(cwd)        ! recover current working directory
C      endif
Cc! console-output from different nodes to different files
Cc      if( mpi__size > 1 ) then
Cc        if(mpi__root )write(*,*)'MPI console outputs to following files.'
Cc        write(*,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
Cc        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
Cc        open(unit=6,file=trim(stdout))
Cc        write(*,"(a,i3)")" ### console output for rank=",mpi__rank
Cc      endif
Cc      if(mpi__root ) then
Cc        close(unit=6)
Cc      endif
C      return
C      end subroutine MPI__Initialize
C
C      subroutine MPI__consoleout(idn)
C      implicit none
C      character(1024*4) :: cwd, stdout
C      character*(*):: idn
C! console-output from different nodes to different files
C      if( mpi__size > 1 ) then
C        if(mpi__root )write(*,*)'MPI console outputs to following files.'
C        write(*,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
C        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
C        open(unit=6,file=trim(stdout))
C        write(*,"(a,i3)")" ### console output for rank=",mpi__rank
C      endif
C      return
C      end subroutine MPI__consoleout
C
C
C      subroutine MPI__Barrier
C      implicit none
C
C      call MPI_Barrier( MPI_COMM_WORLD, mpi__info )
C
C      end subroutine MPI__Barrier
C
C      subroutine MPI__Finalize
C      implicit none
C
C      call MPI_Finalize ( mpi__info )
C
C      end subroutine MPI__Finalize
C
C
C
C      subroutine MPI__getRange( mpi__indexi, mpi__indexe, indexi, indexe )
C      implicit none
C
C      integer, intent(out) :: mpi__indexi, mpi__indexe
C      integer, intent(in)  :: indexi, indexe
C
C      integer, allocatable :: mpi__total(:)
C      integer              :: total
C      integer :: p
C
C      allocate( mpi__total(0:mpi__size-1) )
C
C      total = indexe-indexi+1
C      mpi__total(:) = total/mpi__size
C
C      do p=1, mod(total,mpi__size)
C         mpi__total(p-1) = mpi__total(p-1) + 1
C      end do
C
C      mpi__indexe=indexi-1
C      do p=0, mpi__rank
C         mpi__indexi = mpi__indexe+1
C         mpi__indexe = mpi__indexi+mpi__total(p)-1
C      end do
C      deallocate(mpi__total)
C
C      return
C      end subroutine MPI__getRange
C
C
C      subroutine MPI__Broadcast( data )
C      implicit none
C      integer, intent(inout) :: data
C
C      call MPI_Bcast( data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__info )
C
C      return
C      end subroutine MPI__Broadcast
C
C
C      subroutine MPI__send_d(data,dest)
C      implicit none
C      real(8):: data
C      integer :: n,dest,ierr
C      n=1
C      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
C      end subroutine MPI__send_d
C
C      subroutine MPI__recv_d(data,src)
C      implicit none
C      real(8):: data
C      integer :: n,src,ierr
C      n=1
C      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
C      end subroutine MPI__recv_d
C
C      subroutine MPI__send_dv(data,dest)
C      implicit none
C      real(8):: data(:)
C      integer :: n,dest,ierr
C      n=size(data)
C      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
C      end subroutine MPI__send_dv
C
C      subroutine MPI__recv_dv(data,src)
C      implicit none
C      real(8):: data(:)
C      integer :: n,src,ierr
C      n=size(data)
C      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
C      end subroutine MPI__recv_dv
C
C      subroutine MPI__send_i(data,dest)
C      implicit none
C      integer:: data
C      integer :: n,dest,ierr
C      n=1
C      call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
C      end subroutine MPI__send_i
C
C      subroutine MPI__recv_i(data,src)
C      implicit none
C      integer:: data
C      integer :: n,src,ierr
C      n=1
C      call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
C      end subroutine MPI__recv_i
C
C      subroutine MPI__send_iv(data,dest)
C      implicit none
C      integer:: data(:)
C      integer :: n,dest,ierr
C      n=size(data)
C      call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
C      end subroutine MPI__send_iv
C
C      subroutine MPI__recv_iv(data,src)
C      implicit none
C      integer:: data(:)
C      integer :: n,src,ierr
C      n=size(data)
C      call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
C      end subroutine MPI__recv_iv
C
C
C      subroutine MPI__REAL8send(data,n,dest)
C      implicit none
C      real(8):: data(n)
C      integer :: n,dest,ierr
C      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
C      end subroutine MPI__REAL8send
C
C      subroutine MPI__REAL8recv(data,n,src)
C      implicit none
C      real(8):: data(n)
C      integer :: n,src,ierr
C      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
C      end subroutine MPI__REAL8recv
C
C      subroutine MPI__DbleCOMPLEXsend(data,n,dest)
C      implicit none
C      complex(8):: data(n)
C      integer :: n,dest,ierr
C      call MPI_Send(data,n,MPI_COMPLEX16,dest,mpi__rank, MPI_COMM_WORLD,ierr)
C      end subroutine MPI__DbleCOMPLEXsend
C      subroutine MPI__DbleCOMPLEXrecv(data,n,src)
C      implicit none
C      complex(8):: data(n)
C      integer :: n,src,ierr
C      call MPI_Recv(data,n,MPI_COMPLEX16,src,src, MPI_COMM_WORLD,ista,ierr)
C      end subroutine MPI__DbleCOMPLEXrecv
C
C      subroutine MPI__AllreduceSum( data, sizex )
C      implicit none
C      integer, intent(in) :: sizex
C      complex(8), intent(inout) :: data(sizex)
C      complex(8), allocatable   :: mpi__data(:)
C
C      if( mpi__size == 1 ) return
C
C      allocate(mpi__data(sizex))
C      mpi__data = data
C
C      call MPI_Allreduce( mpi__data, data, sizex,
C     &     MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi__info )
C
C      deallocate( mpi__data )
C
C      return
C      end subroutine MPI__AllreduceSum
C
C
C      subroutine MPI__AllreduceMax( data, sizex )
C      implicit none
C      integer, intent(in) :: sizex
C      integer, intent(inout) :: data(sizex)
C      integer, allocatable   :: mpi__data(:)
C
C      if( mpi__size == 1 ) return
C
C      allocate(mpi__data(sizex))
C      mpi__data = data
C
C      call MPI_Allreduce( mpi__data, data, sizex,
C     &     MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpi__info )
C
C      deallocate( mpi__data )
C
C      return
C      end subroutine MPI__AllreduceMax
C
C
C      subroutine MPI__sxcf_rankdivider(irkip,irkip_all,nspinmx,nqibz,ngrp,nq)
C      implicit none
C
C      integer, intent(out) :: irkip    (nspinmx,nqibz,ngrp,nq)
C      integer, intent(in)  :: irkip_all(nspinmx,nqibz,ngrp,nq)
C      integer, intent(in)  :: nspinmx,nqibz,ngrp,nq
C
C      integer :: ispinmx,iqibz,igrp,iq
C      integer :: total
C      integer, allocatable :: vtotal(:)
C      integer :: indexi, indexe
C      integer :: p
C
C      if( mpi__size == 1 ) then
C         irkip = irkip_all
C         return
C      end if
C
Cc$$$      write(*,*) "irkip_all", mpi__rank
C      total = 0
C      do ispinmx=1, nspinmx
C         do iq=1, nq
C            do iqibz=1, nqibz
C               do igrp=1, ngrp
C                  if( irkip_all(ispinmx,iqibz,igrp,iq) /= 0 ) then
C                     total = total + 1
Cc$$$                     write(*,*) ispinmx, iq, iqibz, igrp, irkip_all(ispinmx,iqibz,igrp,iq)
C                  end if
C               end do
C            end do
C         end do
C      end do
C      write(6,"('MPI__sxcf_rankdivider: nspinmx,nqibz,ngrp,nq,total=',5i6)") nspinmx,nqibz,ngrp,nq,total
C
C      allocate( vtotal(0:mpi__size-1) )
C      vtotal(:) = total/mpi__size
C
C      do p=1, mod(total,mpi__size)
C         vtotal(p-1) = vtotal(p-1) + 1
C      end do
C
C      indexe=0
C      do p=0, mpi__rank
C         indexi = indexe+1
C         indexe = indexi+vtotal(p)-1
C      end do
C      deallocate(vtotal)
C
C      total = 0
C      irkip(:,:,:,:) = 0
C
Cc$$$      write(*,*) "irkip", mpi__rank
C      do iq=1, nq
C         do iqibz=1, nqibz
C            do igrp=1, ngrp
C               do ispinmx=1, nspinmx
C               !spin loop must be internal! Nov2013 (it was in outer loop--> allreduce problem)
C                  if( irkip_all(ispinmx,iqibz,igrp,iq) /= 0 ) then
C                     total = total + 1
C                     if( indexi<=total .and. total<=indexe ) then
C                        irkip(ispinmx,iqibz,igrp,iq) = irkip_all(ispinmx,iqibz,igrp,iq)
Cc           write(6,"(a,5i5)") "isp iq iqibz igrp irkip",
Cc    &                         ispinmx, iq, iqibz, igrp, irkip(ispinmx,iqibz,igrp,iq)
C                     end if
C                  end if
C               end do
C            end do
C         end do
C      end do
C
C      return
C      end subroutine MPI__sxcf_rankdivider
C
C!!
C      subroutine MPI__hx0fp0_rankdivider(iqxini,iqxend,nqibz)
C      implicit none
C      integer, intent(in) :: iqxini, iqxend, nqibz
C
C      integer :: iq
C      allocate( mpi__task(1:iqxend),mpi__ranktab(1:iqxend) )
C
C      if( mpi__size == 1 ) then
C         mpi__task(:) = .true.
C         mpi__ranktab(:) = mpi__rank
C         return
C      end if
C!!
C      mpi__task(:) = .false.
C      do iq=iqxini, iqxend
C        if(iq==1.or. iq>nqibz) then
C           mpi__ranktab(iq) = 0
C        else
C           mpi__ranktab(iq) = mod(iq,mpi__size-1)+1
C        endif
C!!
C        if( mpi__rank == 0 ) then
C            if( iq == 1 .or. iq>nqibz ) then
C               mpi__task(iq) = .true.
C            else
C               mpi__task(iq) = .false.
C            end if
C        else
C            if( iq == 1 .or. iq>nqibz ) then
C               mpi__task(iq) = .false.
C            else
C               if( mpi__rank == mod(iq,mpi__size-1)+1 ) then
C                  mpi__task(iq) = .true.
C               else
C                  mpi__task(iq) = .false.
C               end if
C            end if
C        end if
C      end do
C
C      return
C      end subroutine MPI__hx0fp0_rankdivider
C
C
C      subroutine MPI__hx0fp0_rankdivider2(iqxini,iqxend)
C      implicit none
C      integer, intent(in) :: iqxini, iqxend
C      integer :: iq,i
C      allocate( mpi__task(1:iqxend),mpi__ranktab(1:iqxend) )
C      mpi__task(:) = .false.
C      mpi__ranktab(1:iqxend)=999999
C      if( mpi__size == 1 ) then
C         mpi__task(:) = .true.
C         mpi__ranktab(:) = mpi__rank
C         return
C      end if
C      if(mpi__rank==0) write(6,*) "MPI_hx0fp0_rankdivider2:"
C      do iq=iqxini, iqxend
C         mpi__ranktab(iq) = mod(iq-1,mpi__size)  !rank_table for given iq. iq=1 must give rank=0
C         if( mpi__ranktab(iq) == mpi__rank) then
C            mpi__task(iq) = .true.               !mpi__task is nodeID-dependent.
C         endif
C         if(mpi__rank==0) then
C           write(6,"('  iq irank=',2i5)")iq,mpi__ranktab(iq)
C         endif
C      enddo
C      return
C      end subroutine MPI__hx0fp0_rankdivider2
C
C      end module m_mpi

      subroutine MPI_sxcf_rankdivider(nrank,nspinmx,nqibz,ngrp,nq,irkip_all, irank,irkip)
      implicit none
      integer:: nspinmx,nqibz,ngrp,nq, nrank,irank
      integer:: irkip(nspinmx,nqibz,ngrp,nq)
      integer:: irkip_all(nspinmx,nqibz,ngrp,nq)
      if(nrank==1.and.irank==0) then
         irkip=irkip_all
      else
         stop 'not implemented yet'
      endif
      end
