module driver

   use f90sockets, ONLY : open_socket, writebuffer, readbuffer
   implicit none

        ! SOCKET VARIABLES
        INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
        INTEGER socket, inet, port        ! socket ID & address of the server
        CHARACTER(LEN=1024) :: host

        ! SOCKET COMMUNICATION BUFFERS
        CHARACTER(LEN=12) :: header
        LOGICAL :: isinit=.false., hasdata=.false., verbose = .true.
        INTEGER cbuf, rid
        CHARACTER(LEN=2048) :: initbuffer      ! it's unlikely a string this large will ever be passed...
        DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

        ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
        INTEGER nat
        DOUBLE PRECISION mtxbuf(9), virial(9)

   contains

   subroutine driver_init(lroot)
         character(LEN=256) outs
         logical cmdopt, a2bin, lroot
         integer i

         port=12345
         inet=1
         if (cmdopt('--ipi-host=',11,0,outs)) then
           host = trim(outs(12:))
        endif
        if (cmdopt('--port=',7,0,outs)) then
           i=7
           if (.not.a2bin(outs,port,2,0,' ',i,len(outs)) ) WRITE(*,*) "Could not parse --port command line argument"
        endif
        if (cmdopt('--unix',6,0,outs)) then
           inet=0
        endif

        if (lroot) write(*,*) "Connecting to host: ", trim(host), " port: ", port
        host=trim(host)//achar(0)
        if (lroot) CALL open_socket(socket, inet, port, host)
        nat=1
    end subroutine

    subroutine driver_step(natoms, atoms, cell, pot, forces, stress, lroot)
          use mpi
          IMPLICIT NONE
          LOGICAL lroot
          INTEGER natoms
          REAL(8) atoms(3, natoms), cell(3,3), pot, forces(3,natoms), stress(3,3)
          REAL(8) ipipot
          INTEGER ierror, i
      do
          if (lroot) call readbuffer(socket, header, MSGLEN)
          call MPI_Bcast(header,MSGLEN,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)! broadcast MSGLEN chars to all nodes
          if (lroot) write(*,*) "Read header ", header
          if (trim(header) == "STATUS") then
            if (lroot) then  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
                if (hasdata) then
                   call writebuffer(socket,"HAVEDATA    ",MSGLEN)
                else if (isinit) then
                   call writebuffer(socket,"READY       ",MSGLEN)
                else
                   call writebuffer(socket,"NEEDINIT    ",MSGLEN)
                endif
            endif
         elseif (trim(header) == "INIT") then
            if (lroot) call readbuffer(socket, rid) ! actually this reads the replica id
            if (lroot) call readbuffer(socket, cbuf) ! length of parameter string -- ignored at present!
            if (lroot) call readbuffer(socket, initbuffer, cbuf)
            call MPI_Bcast(rid,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
            isinit=.true.
         elseif (trim(header)=="GETFORCE") then
         if (hasdata) then

            ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper

            DO i = 1, nat
                msgbuffer(3*(i-1)+1:3*i) = forces(:,i) *0.5d0
            ENDDO
            ipipot = pot * 0.5d0
            virial = reshape(stress,(/9/)) * 0.5d0
   !     x     " @ DRIVER MODE: Returning v,forces,stress "
            if (lroot) then
               call writebuffer(socket,"FORCEREADY  ",MSGLEN)
               call writebuffer(socket,ipipot)
               call writebuffer(socket,nat)
               call writebuffer(socket,msgbuffer,3*nat)
               call writebuffer(socket,virial,9)
               call writebuffer(socket,0)
            endif
            hasdata=.false.
            !imginit = .false.
         endif
      elseif  (trim(header)=="POSDATA") then
         if (lroot) then
            call readbuffer(socket, mtxbuf, 9)
            cell=transpose(RESHAPE(mtxbuf, (/3,3/)))
            call readbuffer(socket, mtxbuf, 9)
            call readbuffer(socket, nat)
         endif

         call MPI_Bcast(cell,9,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_Bcast(nat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
         if (nat/=natoms) then
            write(*,*) "Atom number mismatch between driver ", natoms, " and server ", nat
            call exit(-1)
         endif
         if (.not.allocated(msgbuffer)) then
            allocate(msgbuffer(3*nat))
         endif
         if (lroot) call readbuffer(socket, msgbuffer, nat*3)
         call MPI_Bcast(msgbuffer,3*nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         DO i = 1, nat
            atoms(:,i) =  msgbuffer(3*(i-1)+1:3*i)
         ENDDO
!         if (ionode) call readbuffer(socket, combuf, natwrap*3*8)
!         call MPI_Bcast(combuf,3*natwrap,MPI_DOUBLE,0,MPI_COMM_WORLD,ierror)
!         cell=reshape(transpose(cellh), (/9/))*0.52917721 !cell to Angstrom
!         call dcell(cell,celldata)
!         volm = celldata(10)

!         com3n = RESHAPE(combuf, (/ 3 , natwrap /) )*0.52917721 ! positions to Angstrom
         exit
      else
         write(*,*) "Unexpected message from wrapper", trim(header)
      endif
      enddo
    end subroutine
end module driver
