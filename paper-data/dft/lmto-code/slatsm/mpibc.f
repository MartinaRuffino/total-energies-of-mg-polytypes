      subroutine mpibc1(vec,n,cast,mlog,funnam,label)
C- Broadcasts a vector from master node to the world (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   n     :length of vector
Ci   cast  :cast of vector:
Ci         : 1 logical
Ci         : 2 int
Ci         : 4 double
Ci         : 6 double complex
Ci   mlog  : T write message to mlog file
Ci   funnam:string used in writing message (function name)
Ci   label :string used in writing message (variable name)
Co Outputs
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   25 Jun 09 Broadcast large DP arrays in blocks, to avoid error:
Cu             'Message size exceeds P4s maximum message size'
Cu             See parameter MAX_MSG_SIZE
Cu   09 Jul 07 Can broadcast logical vectors
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
      use mpi
      implicit none
      integer numprocs, ierr
      integer MAX_PROCS, MAX_MSG_SIZE
      parameter (MAX_PROCS = 100, MAX_MSG_SIZE=16777216)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master
C ... Passed parameters
      logical mlog
      integer i,n,cast
      double precision vec(n)
      character funnam*(*), label*(*)
C ... Local parameters
      integer nmsg

      if (n <= 0) return
      master = 0
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if (numprocs == 1) return

      if (cast == 1) then
        call MPI_BCAST(vec,n,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
      elseif (cast == 2) then
        call MPI_BCAST(vec,n,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      elseif (cast == 4) then
        if (n > MAX_MSG_SIZE) then
          do  i = 1, n, MAX_MSG_SIZE
            nmsg = min(MAX_MSG_SIZE,n-i+1)
            call MPI_BCAST(vec(i),nmsg,mpi_real8,master,MPI_COMM_WORLD,ierr)
          enddo
        else
          call MPI_BCAST(vec,n,mpi_real8,master,MPI_COMM_WORLD,ierr)
        endif
      elseif (cast == 6) then
        call MPI_BCAST(vec,2*n,mpi_real8,master,MPI_COMM_WORLD,ierr)
      else
        call rxi('mpibc1: cast not implemented',cast)
      endif

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1

        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .    ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'

        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,cast)
      endif

      end

      subroutine mpibcc(vec,n,mlog,funnam,label)
C- Special purpose mpibc1 for broadcasting character strings (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   n     :length of vector
Ci   mlog  : T write message to mlog file
Ci   funnam:string used in writing message (function name)
Ci   label :string used in writing message (variable name)
Co Outputs
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   09 Jul 07 First created
C ----------------------------------------------------------------------\
      use mpi
      implicit none
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master
C ... Passed parameters
      logical mlog
      integer n
      character vec(n)*(*)
      character funnam*(*), label*(*)
C ... Local parameters

      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs == 1) return
      if (n <= 0) return
      master = 0
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_BCAST(vec,n,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)

      end

      subroutine mpibc3(vec,n,cast,pid,mlog,funnam,label)
C- Broadcasts a vector from specified node to the world (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   n     :length of vector
Ci   cast  :cast of vector:
Ci         : 2 int
Ci         : 4 double
Ci         : 6 double complex
Ci   mlog  : T write message to mlog file
Ci   funnam:string used in writing message (function name)
Ci   label :string used in writing message (variable name)
Co Outputs
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
      use mpi
      implicit none
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master
C ... Passed parameters
      logical mlog
      integer n,cast,pid
      double precision vec(n)
      character funnam*(*), label*(*)
C ... Local parameters

      if (n <= 0) return
      master = 0
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs == 1) return
      if (cast == 2) then
        call MPI_BCAST(vec,n,MPI_INTEGER,pid,MPI_COMM_WORLD,ierr)
      elseif (cast == 4) then
        call MPI_BCAST(vec,n,mpi_real8,pid,MPI_COMM_WORLD,ierr)
      elseif (cast == 6) then
        call MPI_BCAST(vec,2*n,mpi_real8,pid,MPI_COMM_WORLD,ierr)
      else
        call rxi('mpibc3: cast not implemented',cast)
      endif

C      if (mlog) then
C        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C        call strcop(shortname(procid),name,10,'.',ierr)
C        namelen(procid) = ierr-1
C
C        call gettime(datim)
C        strn = ' '//funnam//' '//datim//' Process %i of %i on '
C     .    //shortname(procid)(1:namelen(procid))//' bcast '//label//
C     .    ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
C
C        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
C     .    cast)
C      endif

      end

      subroutine mpibc2(vec,n,cast,op,mlog,funnam,label)
C- Performs MPI_ALLREDUCE on a vector (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   n     :length of vector
Ci   cast  :cast of vector:
Ci         : 2 int
Ci         : 4 double
Ci         : 6 double complex
Ci   op    : 0 do nothing
Ci         : 1 maximum                       (MPI_MAX)
Ci         : 2 minimum                       (MPI_MIN)
Ci         : 3 sum                           (MPI_SUM)
Ci         : 4 product                       (MPI_PROD)
Ci         : 5 logical and                   (MPI_LAND)
Ci         : 6 bit-wise and                  (MPI_BAND)
Ci         : 7 logical or                    (MPI_LOR)
Ci         : 8 bit-wise or                   (MPI_BOR)
Ci         : 9 logical xor                   (MPI_LXOR)
Ci         :10 bit-wise xor                  (MPI_BXOR)
Ci         :11 max value and location        (MPI_MAXLOC) NOT IMPLEMENTED
Ci         :12 min value and location        (MPI_MINLOC) NOT IMPLEMENTED
Ci   mlog  : T write message to mlog file
Ci   funnam:string used in writing message (function name)
Ci   label :string used in writing message (variable name)
Co Outputs
Cl Local variables
Cr Remarks
Cr   ALLREDUCE sums the contributions from all the individual threads
Cr   The ALLREDUCE operation is specified by op
Cu Updates
Cu   07 Nov 12 Add op to enable different MPI_ALLREDUCE operations
Cu   25 Jun 09 Broadcast large DP arrays in blocks, to avoid error:
Cu             'Message size exceeds P4s maximum message size'
Cu             See parameter MAX_MSG_SIZE
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
      use mpi
      implicit none
      integer numprocs, ierr
      integer MAX_PROCS, MAX_MSG_SIZE
      parameter (MAX_PROCS = 100, MAX_MSG_SIZE=16777216)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master
C ... Passed parameters
      logical mlog
      integer n,cast,op
      double precision vec(n)
      character funnam*(*), label*(*)

C ... Local parameters
      integer i,nmsg,lop
      integer, allocatable :: ibuf(:)
      real(8) ,allocatable :: dbuf(:)
C ... Heap (needed only if avoid allocate)
      integer w(1)
      common /w/ w
      integer obuf

      if (n <= 0) return
      master = 0
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs == 1) return

      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )

      lop = 0
      select case (op)
        case ( 1); lop = MPI_MAX
        case ( 2); lop = MPI_MIN
        case ( 3); lop = MPI_SUM
        case ( 4); lop = MPI_PROD
        case ( 5); lop = MPI_LAND
        case ( 6); lop = MPI_BAND
        case ( 7); lop = MPI_LOR
        case ( 8); lop = MPI_BOR
        case ( 9); lop = MPI_LXOR
        case (10); lop = MPI_BXOR
        case (11); lop = MPI_MAXLOC
        case (12); lop = MPI_MINLOC
      end select

C      if (lop < 1 .or. lop > 10) return

      if (cast == 2) then
        allocate(ibuf(n), stat=ierr)
        call MPI_ALLREDUCE(vec,ibuf,n,MPI_INTEGER,lop,MPI_COMM_WORLD,ierr)
        call icopy(n,ibuf,1,vec,1)
        deallocate(ibuf, stat=ierr)
      elseif (cast == 4) then
        if (n > MAX_MSG_SIZE) then
          allocate(dbuf(MAX_MSG_SIZE), stat=ierr)
          do  i = 1, n, MAX_MSG_SIZE
            nmsg = min(MAX_MSG_SIZE,n-i+1)
            call MPI_ALLREDUCE(vec(i),dbuf,nmsg,mpi_real8,lop,MPI_COMM_WORLD,ierr)
            call dcopy(nmsg,dbuf,1,vec(i),1)
          enddo
          deallocate(dbuf, stat=ierr)
        else
          allocate(dbuf(n), stat=ierr)
          call MPI_ALLREDUCE(vec,dbuf,n,mpi_real8,lop,MPI_COMM_WORLD,ierr)
          call dcopy(n,dbuf,1,vec,1)
          deallocate(dbuf, stat=ierr)
        endif
      elseif (cast == 6) then
        allocate(dbuf(2*n), stat=ierr)
        call MPI_ALLREDUCE(vec,dbuf,2*n,mpi_real8,lop,MPI_COMM_WORLD,ierr)
        call dcopy(2*n,dbuf,1,vec,1)
        deallocate(dbuf, stat=ierr)
C        call defrr(obuf,2*n)
C        call MPI_ALLREDUCE(vec,w(obuf),2*n,
C     .    mpi_real8,lop,MPI_COMM_WORLD,ierr)
C        call dcopy(2*n,w(obuf),1,vec,1)
C        call rlse(obuf)
      else
        call rxi('mpibc2: cast not implemented',cast)
      endif

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//' allreduce '//label
        call awrit2(strn,' ',-256,lgunit(3),procid,numprocs)
      endif

      end

      subroutine mpibc4(veci,vecd,kpproc,n,cast)
C- Performs MPI_ALLGATHERV on a vector of blocks (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   kpproc:processor i will contain data in a continguous block
Ci         :whose length is
Ci         :  length(i) = [kpproc(i) - kpproc(i-1)]*n
Ci         :Offsets are ordered by processor; thus
Ci         :  offset(1) = 0
Ci         :  offset(i+1) = offset(i) + length(i)
Ci   n     :enters into length of vector
Ci   cast  :cast of vector:
Ci         : 2 int
Ci         : 4 double
Ci         : 6 double complex
Co Outputs
Cl Local variables
Cr Remarks
Cu Updates
Cu    2 Jan 11 First created
C ----------------------------------------------------------------------
      use mpi
      implicit none
      integer numprocs, ierr
      integer MAX_PROCS, MAX_MSG_SIZE
      parameter (MAX_PROCS = 100, MAX_MSG_SIZE=16777216)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master
C ... Passed parameters
      integer n,cast
      integer kpproc(0:*)
      integer veci(*)
      double precision vecd(*)
C ... Local parameters
      integer i,nmsg,ista,iend,nelts
      integer, allocatable :: ibuf(:),offset(:),length(:)
      real(8) ,allocatable :: dbuf(:)
C ... Heap (needed only if avoid allocate)
      integer w(1)
      common /w/ w
      integer obuf

      if (n <= 0) return
      master = 0
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs == 1) return
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )

      allocate (offset(0:numprocs), stat=ierr)
      allocate (length(0:numprocs), stat=ierr)
      offset(0) = 0

      do  i = 0, numprocs-1
        ista = kpproc(i)
        iend = kpproc(i+1)-1
        length(i) = (iend - ista + 1)*n
        offset(i+1) = offset(i) + length(i)
C       print *, 'ista=',i,ista,iend,offset(i),length(i)
      enddo
C     print *, offset(0:numprocs-1)
C     print *, length(0:numprocs-1)
      nelts = offset(numprocs)

      if (cast == 2) then
        call rx('mpibc4 not ready')
      elseif (cast == 4) then
        allocate(dbuf(nelts), stat=ierr)
        ista = kpproc(procid)
        call MPI_ALLGATHERV(vecd(1+offset(procid)),
     ,    length(procid),mpi_real8,dbuf,length,
     ,    offset,mpi_real8,MPI_COMM_WORLD,ierr)
        call dcopy(nelts,dbuf,1,vecd,1)
        deallocate(dbuf)
      elseif (cast == 6) then
        call rx('mipbc4 not ready')
      else
        call rxi('mpibc2: cast not implemented',cast)
      endif

      end

      subroutine mpibc5(vec,n,cast,tag,pidfrom,pidto)
C- Sends/receives vector from one process to another
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   n     :length of vector
Ci   cast  :cast of vector:
Ci         : 2 int
Ci         : 4 double
Ci         : 6 double complex
Ci  tag    : identifier of particular process
Co Outputs
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
      use mpi
      implicit none
C ... Passed parameters
      integer n,cast,tag,pidfrom,pidto
      double precision vec(n)
C ... Local parameters
      integer numprocs, ierr, nsend, mpicast, procid

C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      print 777, 'enter mpibc5',procid,pidfrom,pidto

      if (n <= 0) return
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs == 1) return

      nsend = n
      select case (cast)
        case (2);   mpicast = mpi_integer
        case (4,6); mpicast = mpi_real8
      case default
        call rxi('mpibc5: illegal value for cast:',cast)
      end select
      if (cast == 6) nsend = 2*n

C     print 777, 'mpibc5:',procid,n,nsend

      if (procid == pidfrom) then
C       print 777, 'send b',procid,pidfrom,pidto,vec(1),sum(vec)
C 777   format(a,3i5,10f14.7)
        call MPI_SEND(vec, nsend, mpicast, pidto, tag, MPI_COMM_WORLD, ierr)
      elseif (procid == pidto) then
C       print 777, 'recv b',procid,pidfrom,pidto,vec(1),sum(vec)
        call MPI_RECV(vec, nsend, mpicast, pidfrom, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
C       print 777, 'recv a',procid,pidfrom,pidto,vec(1),sum(vec)
      endif

      end

      integer function mpipid(mode)
C- Returns MPI procid
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return number of processors
Ci         :1 return procid
Ci         :2 calls MPI_BARRIER; returns ierr
Ci         :Otherwise, return 0
Co Outputs
Co   mpipid:procid or number of processors (see mode)
Cr Remarks
Cu Updates
Cu   24 Nov 05 Added mode 2
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
      use mpi
      implicit none
      integer numprocs, ierr, procid
      integer mode

      if (mode == 0) then
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
        mpipid = numprocs
      else if (mode == 1) then
        call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
        mpipid = procid
      else if (mode == 2) then
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        mpipid = ierr
      else
        mpipid = 0
      endif

      end

      real(8) function mpiquery(mode)
      use mpi
      integer mode

      mpiquery = 0
      if (mode == 0) then
        mpiquery = mpi_wtime()
      endif
      end

      subroutine mpiprn(shortname,length)
C- Returns MPI short proc name and length
C ----------------------------------------------------------------------
Co Outputs
Co   process name shortened to fewer than 11 chars, and length
Cr Remarks
Cu Updates
Cu   06 Jun 09 First created
C ----------------------------------------------------------------------
       use mpi
       implicit none
      character*(MPI_MAX_PROCESSOR_NAME) name
      integer resultlen, ierr, numprocs
      character*10 shortname
      integer length,i

      call mpi_comm_size( mpi_comm_world, numprocs, ierr )
      if ( numprocs > 1 ) then
      call mpi_get_processor_name(name, resultlen, ierr)
      call strcop(shortname,name,10,'.',i)
      length = i - 1
      else
      shortname = ' '
      length = 0
      endif

      end
C      subroutine fmain
C      implicit none
C      integer n,omppid,mpipid
C
C      n = 0
C
CC#ifdefC MPI
CC      print *, 'processor number, id:', mpipid(0),mpipid(1)
CC      call MPI_FINALIZE(n)
CC#endif
C
CC#ifdefC OPENMP
CC      n = omppid(0)
CC      print *, 'omppid for number of processors outside brackets:', n
CCC$OMP PARALLEL PRIVATE(n)
CC      n = omppid(0)
CC      print *, 'omppid for number of processors:', n
CCC$OMP END PARALLEL
CC
CC      n = omppid(1)
CC      print *, 'omppid(1), outside openmp PARALLEL:', n
CC
CCC$OMP PARALLEL PRIVATE(n)
CC      n = omppid(1)
CC      print *, 'omppid(1) for processor id:', n
CC!$OMP do
CC      do  n = 1, 10
CC      print *, 'omppid(1) in do loop: i, procid =', n,omppid(1)
CC      enddo
CC!$OMP end do
CC
CCC$OMP END PARALLEL
CC#endif
C
C      end
