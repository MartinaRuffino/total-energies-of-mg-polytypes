      subroutine iostrt(nbas,nf,iter,bas,vel,eps,zeta,zacc,zsqacc,veps,
     .                  time,ifi,ierr)
C- I/O positions and velocities
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas; nf, number of degrees of freedom (set to 3*nbas-3 in initv);
Ci   iter, current MD iteration number; bas, vel, current positions and
Ci   velocities; eps, current volume strain; zeta, Hoover viscosity;
Ci   zacc, ditto accumulated over the whole MD run; zsqacc, zeta^2
Ci   accumulated over the whole MD run; veps, barostat viscosity;
Ci   time, current MD time; ifi, file handle
Co Outputs:
Co   ierr, output error code: 0 on successful read; 1 if read failed
Cr Remarks
C ----------------------------------------------------------------------
Cu Updates
Cu   06 Dec 12 (DMT) Revamp with bugfixes
Cu   07 Jul 09 (ATP) updated for MPI, changed output to ASCII
Cu   01 Jun 08 Rewritten by ATP
C ----------------------------------------------------------------------
C     use mpi
      implicit none
C Passed parameters
      integer nbas,nf,iter,ifi,ierr
      double precision bas(3,nbas),vel(3,nbas),eps,zeta,zacc,zsqacc,
     .                 veps,time, tube(6)
C Local parameters
      integer nbas0,iprint,i1mach, e(4), ipr, iomode
      real(8), parameter :: autime = 0.048377d0
      character(len=6), parameter :: fmt1='(3i10)'
      character(len=9), parameter :: fmt3='(6g24.16)', fmt4='(3g24.16)'
C For MPI ...
      integer procid,master,mpipid
      logical ioerr

      ierr = 0
      ioerr = .false.
      master = 0
C     comm = mpi_comm_world

      procid = mpipid(1)
C     call mpi_comm_rank(comm, procid, ierr)

      if (procid == master) iomode = ifi
C     call mpi_bcast(iomode, 1, mpi_integer, master, comm, ierr)
      call mpibc1(iomode,1,2,.false.,' ',' ')

      if (procid == master) then
      if (iomode > 0) then
        rewind ifi
        read (ifi,fmt1,iostat=e(1)) nbas0,nf,iter
        read (ifi,fmt3,iostat=e(2)) eps,zeta,zacc,zsqacc,veps,time
        read (ifi,fmt4,iostat=e(3)) bas
        read (ifi,fmt4,iostat=e(4)) vel
        ioerr = e(1) /= 0 .or. any(e > 0) .or. nbas /= nbas0
        ipr = iprint()
        if ( ipr >= 20 .and. (.not. ioerr)) then
          call awrit0('IOSTRT: read from STRT file',' ',120,i1mach(2))
          call awrit0('        ... continuing MD run',' ',120,i1mach(2))
          call awrit5('        eps=%d, zeta=%d, zacc=%d, veps=%d,'//
     .                ' time=%d',' ',120,i1mach(2),eps,zeta,zacc,veps,
     .                time*autime)
        else if (ipr >= 10 .and. ioerr) then
          call awrit0(
     &              ' IOSTRT: Error reading STRT file',' ',60,i1mach(2))
          rewind ifi
          iter = 0
        endif
      elseif (iomode < 0) then
        write (-ifi,fmt1) nbas,nf,iter
        write (-ifi,fmt3) eps,zeta,zacc,zsqacc,veps,time
        write (-ifi,fmt4) bas
        write (-ifi,fmt4) vel
      endif
      end if

C     call mpi_bcast(ioerr,1,mpi_logical,master,comm,ierr)
      call mpibc1(ioerr,1,1,.false.,' ',' ')

      if ((.not. ioerr) .and. (iomode > 0)) then
        if (procid == master) then
          e(1:2) = (/ nf, iter/)
          tube = (/eps,zeta,zacc,zsqacc,veps,time/)
        end if

C       call mpi_bcast(e, 2, mpi_integer,master,comm,ierr)
C       call mpi_bcast(tube,6,mpi_real8,master,comm,ierr)
        call mpibc1(e,2,2,.false.,' ',' ')
        call mpibc1(tube,6,4,.false.,' ',' ')

        if (procid /= master) then
          nf     = e(1)
          iter   = e(2)
          eps    = tube(1)
          zeta   = tube(2)
          zacc   = tube(3)
          zsqacc = tube(4)
          veps   = tube(5)
          time   = tube(6)
        end if

C       call mpi_bcast(bas, 3*nbas, mpi_real8, master, comm, ierr)
C       call mpi_bcast(vel, 3*nbas, mpi_real8, master, comm, ierr)
        call mpibc1(bas,3*nbas,4,.false.,' ',' ')
        call mpibc1(vel,3*nbas,4,.false.,' ',' ')

      end if

      if ((ierr == 0 .and. ioerr) .or. (ierr /= 0)) ierr = 1

      end subroutine iostrt
