C#define LM
C -------------- Main program for LM and related programs ----------
C This file contains the main program for most of the collection
C of main programs in the Questaal package.
C
C A brief explanation for most of the codes in the package
C can be found at https://questaal.org/docs/package_overview
C
C Many of the main programs are created from this file using the ccomp utility,
C documented at https://questaal.org/docs/misc/ccomp
C
C Updates
C   The repository keeps a log of all changes.
C   See also doc/ChangeLog for significant changes.
C ----------------------------------------------------
C ... Operating system permitting, The true main is in slatsm.a (fmain.c)

      subroutine fmain
      use mpi
      use structures
      implicit none

      integer procid, master, nproc
C#ifdefC MPE
C      include "mpef.h"
C#endif
C     integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer ierr
      integer, parameter :: MAX_PROCS = 100
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
C     double precision starttime, endtime
      logical mlog
C#ifdefC MPE
CC Event numbers:
C      include "fp/events.ins"
C#endif MPE


C Heap allocation
       integer wksize
C#ifdefC LMMC
C      parameter(wksize= 80 000 000)
C      integer w(wksize)
CC     Next two lines guarantee w is aligned along a d.p. boundary
C      double precision ws
C      equivalence (ws,w(1))
C      common /w/ w
C#else
      parameter(wksize= 0)
C#endif

C ... Controls for IO
      integer lstrn
      parameter (lstrn=2000)

      character prgnam*8, vrsion(2)*6, ext*20

C ... For structures
!       include 'structures.h'
      integer mxspec
      parameter (mxspec=256)
      character*8 slabl(mxspec)
      type(str_bz):: s_bz
      type(str_ctrl):: s_ctrl
      type(str_gw):: s_gw
      type(str_dmft):: s_dmft
      type(str_ham):: s_ham
      type(str_lat):: s_lat
      type(str_mix):: s_mix
      type(str_move):: s_move
      type(str_optic):: s_optic
      type(str_pot):: s_pot
      type(str_str):: s_str
      type(str_tb):: s_tb
      type(str_spec),pointer:: s_spec(:),s_spec0(:)
      type(str_site),pointer:: s_site(:),s_site0(:)
      integer,parameter :: nstrn=10
      type(str_strn) :: s_strn(nstrn)
C     type(str_symops),allocatable::  s_sym(:)

C     Used for initial dimensions of s_site, s_spec
      integer, parameter :: nbmx=5000

C ... Miscellaneous local variables
      character strn*1000,outs*20
      integer i,j,k,lc,stdo
C     ivn describes version:
C     ivn(:,1) corresponds to input file system
C     ivn(:,2) corresponds to program version (usually identical to ivn(:,1))
C     ivn(1,:) main version number.  Different versions have some incompatibilities.
C     ivn(2,:) minor version number.
C     ivn(3,1) revision within minor version
      integer ivn(3,2)
      logical swtmp
      integer, parameter :: NULLI=-99999
      logical, parameter :: T=.true., F=.false.
      procedure(real(8)) :: dglob,dlength,ddot,dmach
      procedure(logical) :: cmdopt
      procedure(integer) :: idalloc,allocvb   ! For dynamic memory allocation
      procedure(integer) :: fxst,fadd,fopn,lgunit,i1mach,fextg,a2vec,str_pack

C#ifdefC LMCHK | LMXBS | LMGPOL | LMSHF
C      integer auxmod
C#endif

C#ifdef LM | TBE | LMGF | LMF
      logical ltet
C#endif


C#ifdefC LMCTL
C      integer fopnn
C#endif
      integer fext
C#ifdefC TBE
C      integer ltbe
C      double precision mdprm(7)
C#endif

C ... Program-dependent name and help
C#ifdefC LMPG
C#elseif LM
C#elseifC LMFA
C      data prgnam /'LMFA'/
C#elseifC LMFGW
C      data prgnam /'LMFGW'/
C#elseifC LMFGWD
C      data prgnam /'LMFGWD'/
C#elseifC LMFGWS
C      data prgnam /'LMFGWS'/
C#elseifC LMFDMFT
C      data prgnam /'LMFDMFT'/
C#elseifC LMF
C      data prgnam /'LMF'/
C#elseifC LMAQU
C      data prgnam /'LMAQU'/
C#elseifC  LMDOS
C      data prgnam /'LMDOS'/
C#elseifC  LMCTL
C      data prgnam /'LMCTL'/
C#elseifC  LMIMP
C      data prgnam /'LMIMP'/
C#elseifC  LMCOR
C      data prgnam /'LMCOR'/
C#elseifC  LMFIT
C      data prgnam /'LMFIT'/
C#elseifC  LMPLAN
C      integer nbas,nbasp,nbaspp,nl
C      double precision emad,trumad,vmtz(2)
C      data prgnam /'LMPLAN'/
C#elseifC  LMAVGM
C      data prgnam /'LMAVGM'/
C#elseifC  LMMIX
C      data prgnam /'LMMIX'/
C#elseifC  LMSHF
C      data prgnam /'LMSHF'/
C      data auxmod /8/
C#elseifC  TBPG
C      data prgnam /'TBPG'/
C#elseifC  TBE
C      data prgnam /'TBE'/
C#elseifC  TBBND
C      data prgnam /'TBBND'/
C#elseifC  TBFIT
C      data prgnam /'TBFIT'/
C#elseifC  TBDOS
C      data prgnam /'TBDOS'/
C#elseifC  MMAG
C      data prgnam /'MMAG'/
C#endif

C ... Program-dependent cagetories
C#ifdefC LMPG
C      data prgnam /'LMPG'/
C#elseifC LMGPOL
C      data prgnam /'LMGPOL'/
C      data auxmod /32/
C#elseifC BLM
C      data prgnam /'BLM'/
C#elseifC LM67
C      data prgnam /'LM67'/
C#elseif LM
      data prgnam /'LM'/
C#elseifC LMMC
C      data prgnam /'LMMC'/
C#elseifC LMAQU | LMMIX
C      call rx('missing setup for lmaqu, lmmix')  ! Why???TK
C#elseifC LMGF
C      data prgnam /'LMGF'/
C#elseifC LMSTR
C      data prgnam /'LMSTR'/
C#elseifC LMCHK
C      data prgnam /'LMCHK'/ auxmod /1/
C#elseifC LMSCELL
C      procedure(logical) :: rdstrn
C      procedure(integer) :: nglob
C      integer nkd,nbx,ix(9),iplx(3,3)
C      double precision plat(3,3),plx(3,3),xx
C      data prgnam /'LMSCELL'/
C#elseifC LMXBS
C      data prgnam /'LMXBS'/
C      data auxmod /4/
C#endif

      integer:: nfilin,nrecs,fopna
C#ifdefC TBE
C      integer, parameter :: mxrecs=50000
C#else
      integer, parameter :: mxrecs=10000
C#endif
      integer, parameter :: recln0=120
      character*8 alabl
      character,allocatable:: recrd(:)
      real(8),parameter:: NULLR =-99999

C -------------- First executable statement ---------------
C     Link in routines that cause potential library conflicts
      call nada
C --- Version ---
      vrsion(1) = 'LM'
      vrsion(2) = ' '
      ivn(1,1) = 7
      ivn(2,1) = 14
      ivn(3,1) = 1
      ivn(:,2) = ivn(:,1)
C#ifdef LM | LMPG | LMGF | LMAVGM | LMIMP | LMCTL | LMSTR
      vrsion(2) = 'ASA'
C#elseifC LMFA | LMF | LMFGWD | LMFGWS
C      vrsion(2) = 'FP'
C#elseifC LMMC
C      vrsion(2) = 'MOL'
C      ivn(1,2) = 3
C      ivn(2,2) = 2
C#elseifC TBE
C      vrsion(2) = 'TB'
C      ivn(1,2) = 10
C      ivn(2,2) = 1
C      ivn(3,2) = 0
C#elseifC MMAG
C      vrsion(2) = 'MM'
C      ivn(1,2) = 2
C      ivn(2,2) = 1
C#endif
      stdo = lgunit(1)
      i = dglob('stdo',dble(stdo),1)
      call lodsyv('null',0,nullr,i)
      master = 0
      call mpi_comm_rank( mpi_comm_world, procid, ierr )
      call mpi_comm_size( mpi_comm_world, nproc, ierr )
      call finits(22,0,0,i)
      call logwarn(-1,'') ! Initialize warning log
C --- Help ---
      swtmp = .false.
      if (swtmp .or. cmdopt('--h',3,0,outs)) call lmhelp(prgnam,ivn)
      if (cmdopt('--version',9,0,outs)) then
        outs = ' ' // vrsion(1)
        call lmvsn(ivn,outs)
        strn = trim(outs) // '  ' // vrsion(2)
        call lmvsn(ivn(1,2),strn)
        call info0(0,0,0,trim(strn))
        call cexit(0,1)
      endif

C --- Dynamic memory allocation and other initialization ---
C ... Initialize s_strn ... needed for gfortran compiler
      do  i = 1, nstrn
          nullify(s_strn(i)%strn)
      enddo

      if (nproc > 1) then
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      mlog = cmdopt('--mlog',6,0,strn)
C#ifdefC MPE
C      ierr = MPE_INIT_LOG()
C      EVENT_START_RDCTRL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_RDCTRL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_UGCOMP = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_UGCOMP   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_FSMBL  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_FSMBL    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_PZHEV  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_PZHEV    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_KLOOP  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_KLOOP    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_MIXRHO = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_MIXRHO   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_SMHSBL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_SMHSBL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_AUGMBL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_AUGMBL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_HSIBL  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_HSIBL    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_RSIBL  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_RSIBL    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_RLOCBL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_RLOCBL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_DFRCE  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_DFRCE    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_BCAST  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_BCAST    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_ALLRED = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_ALLRED   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_BARRIER= MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_BARRIER  = MPE_LOG_GET_EVENT_NUMBER()
C#endif
C --- Dynamic memory allocation and other initialization ---
      if (procid == master) call headl2(prgnam,wksize,stdo)
C#ifdefC LMMC
C      call pshpr(1)
C      call wkinit(wksize)
C      call wkfast(T)
C      call poppr
C#endif
      if (procid == master) then
        i = fextg(ext)
      endif
      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
      if (procid == master) then
C#ifdefC MPE
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_RDCTRL,EVENT_END_RDCTRL,"rdctrl","pink")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_UGCOMP,EVENT_END_UGCOMP,"ugcomp","maroon")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_FSMBL,EVENT_END_FSMBL,  "fsmbl","aquamarine")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_PZHEV,EVENT_END_PZHEV,  "pzhev","brown")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_KLOOP,EVENT_END_KLOOP,  "k-loop","brown")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_MIXRHO,EVENT_END_MIXRHO,"mixrho","orange")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_SMHSBL,EVENT_END_SMHSBL,"smhsbl","blue")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_AUGMBL,EVENT_END_AUGMBL,"augmbl","cyan")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_HSIBL,EVENT_END_HSIBL,  "hsibl","gray")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_RSIBL,EVENT_END_RSIBL,  "rsibl","red")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_RLOCBL,EVENT_END_RLOCBL,"rlocbl","green")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_DFRCE,EVENT_END_DFRCE,  "dfrce","magenta")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_BCAST,EVENT_END_BCAST,  "broadcast","coral")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_ALLRED,EVENT_END_ALLRED,"allreduce","purple")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_BARRIER,EVENT_END_BARRIER,"barrier","yellow")
C#endif
        call gettime(datim)
        if (mlog) i = fopn('MLOG')
        if (.not. cmdopt('--input',6,0,strn)) call poseof(fopn('LOG'))
        if (mlog) then
          call awrit2(' '//prgnam//' '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//' is master',' ',
     .      256,lgunit(3),procid,nproc)
        endif
      else
        call strcat(ext,20,' ','_',1,' ',i)
        call bin2a(' ',0,0,procid,2,0,20,ext,i)
        ierr = fext(ext(1:i+1))
        if (mlog) ierr = fopn('MLOG')
        ierr = fextg(ext)
        call gettime(datim)
        if (mlog) then
          call awrit2(prgnam//' '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' file extension is '//ext(2:i+1),' ',
     .      256,lgunit(3),procid,nproc)
        endif
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      else
      call headl2(prgnam,wksize,stdo)
C     If using earlier than SLATSM.52, first argument should be 2
C      call finits(22,0,0,i)
C#ifdefC LMMC
C      call pshpr(1)
C      call wkinit(wksize)
C      call wkfast(T)
C      call poppr
C#endif
      if (.not. cmdopt('--input',6,0,strn)) call poseof(fopn('LOG'))
      i = fextg(ext)
      call word(ext,1,i,j)
      if (ext(i:i) == '.') i=i+1
      if (ext(j:j) == '.') j=j-1
      if (j >= i) call ptenv('EXT='//ext(i:j))
      endif

C#ifdefC BLM
C      call suctrl(prgnam)
C      call fexit(0,0,' ',0)
C#else

C ... Abort with error message if ctrl file is missing (swtmp = .true.)
      swtmp = .false.
      if (procid == master) then
C       if (cmdopt('--input',6,0,strn).and.nproc > 1)
C      .  call rx('--input not allowed with MPI')
      if (.not. cmdopt('--input',6,0,strn)) then
        if (fxst('CTRL') /= 1) then
          call awrit0(' '//prgnam//'%a:%9pmissing ctrl file',' ',80,
     .      i1mach(2))
          swtmp = .true.
        endif
      endif
      else
        if (cmdopt('--input',6,0,strn)) call cexit(0,1)
      endif

      call mpibc1(swtmp,1,1,.false.,prgnam,'error')
      if (swtmp) call cexit(-1,1)

C ... Set special file directory for temporary files
C     User may which to customize the directory
C     Default is to use the standard directory
C     sttmpd is located at the bottom of this file.
      call sttmpd

C ... File logical units
C     i = fadd('TMP',-1,4)
C     i = fadd('BAND',-1,4)
C#ifdefC TBE
C      i = fadd('STRT',-1,0)
C      i = fadd('QMOM',-1,0)
C#endif

C --- Set the top-level verbosity if spec'd from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i == 1) call setpr(k)
      endif

C --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
      allocate( recrd( 0:mxrecs*recln0-1 ) )
      if (procid == master) then
      if (.not.cmdopt('--input',7,0,strn)) then
        nfilin = fopna('CTRL',-1,1)
        swtmp = cmdopt('--show',6,0,strn)
        swtmp = swtmp .and. scan(strn(7:7),' p=')>0
        alabl = '#{}% ct '; if (swtmp) alabl = '#{}% ctp'
        call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs)
        i = 60
        if (swtmp) then
          i = 1
          call info0(i,0,0,' ------------------------ End '//
     .                     'of input file ----------------------')
        endif
        call info2(i,0,1,' '//prgnam//'%a : %i lines read from'//
     .    ' input file',nrecs,0)
        if (cmdopt('--showp',7,0,strn)) call cexit(0,1)
      endif
      endif
      if (nproc > 1) then
!         call mpibc1(nrecs,1,2,mlog,prgnam,'nrecs')
!         call mpibcc(recrd,recln0*(nrecs+1),mlog,prgnam,'ctrl')
        call mpi_bcast(nrecs,1,mpi_integer,master,mpi_comm_world,ierr)
        call mpi_bcast(recrd,recln0*nrecs,mpi_character,master,mpi_comm_world,ierr)
      endif

C --- Read input file ---
C#ifndef LM67
C     Ridiculously large s_site and s_spec are initially allocated;
C     later they are reduced to actual size.  rdctrl checks to ensure
C     that the initial s_site0 and s_spec0 are sufficiently large.
      allocate(s_site0(nbmx),s_spec0(nbmx))

C     Initialize strux
      call bcast_strx(-(2**9-1),s_bz,s_ctrl,s_ham,s_pot,s_lat,
     .  s_mix,s_spec0,s_site0,s_str,nbmx,nbmx)

C     Read input
      call rdctrl(recrd,recln0,nrecs,prgnam,vrsion,ivn,1,slabl,
     .  s_bz,s_ctrl,s_ham,s_pot,s_lat,s_mix,s_spec0,s_site0,nbmx,
     .  s_str,s_move,s_tb,s_optic,s_gw,s_dmft,s_strn)

C     Resize s_site, s_spec
      allocate(s_site(s_ctrl%nsite))
      do j = 1, s_ctrl%nsite ! use explicit loop to avoid stack issues when there are many atoms
        s_site(j) = s_site0(j)
      end do

      allocate(s_spec(s_ctrl%nspec))
      do j = 1, s_ctrl%nspec
        s_spec(j) = s_spec0(j)
      end do

      deallocate(s_site0,s_spec0)

      allocate(s_ctrl%spid(j))
      do  j = 1, s_ctrl%nspec
        s_ctrl%spid(j) = slabl(j)
      enddo
C#elseC
C      call rdctrlchk(recrd,recln0,nrecs,'LM',vrsion,ivn,F,slabl)
C      deallocate(recrd)
C      call rx0('lm67')
C#endif

       if (nproc > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_RDCTRL,procid,"rdctrl")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endif
       call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C#endif
       endif

C#ifdefC LMSCELL
CC ... Scale alat by factor, and plat,pos by its inverse
C      i = 8
C      if (cmdopt('--scala=',i,0,outs)) then
C        i = a2vec(outs,len(outs),i,4,', ',2,2,1,ix,xx)
C        if (i /= 1) call rx('lmscell failed to parse '//trim(outs))
C        s_lat%alat = s_lat%alat*xx
C        call info(10,0,0,' ... scale alat, inverse plat by factor %d',xx,0)
C        call dscal(9,1/xx,s_lat%plat,1)
C        nbx = nglob('nbasp')
C        do  i = 1, nbx
C          call dscal(3,1/xx,s_site(i)%pos,1)
C        enddo
C      endif
C
C      call setcg(s_lat,8,12)
C      s_lat%as = 1d0
C      s_lat%tol = dmach(2)
C      s_lat%nkdmx = 150000
C      s_lat%nkqmx = 150000
C      call lattic(s_lat,s_ctrl,s_site)
C      nkd = s_lat%nkd
C      plat = s_lat%plat
C      plx = s_lat%slat
C
CC ... Set default values for species data
C      call defspc(s_spec)
C
CC  ...Stack mode
C      if (cmdopt('--stack',7,0,strn)) then
C        call stackcel(strn(8:),s_ctrl,s_lat,s_site,s_spec,slabl)
C        call cexit(0,1)
C      endif
C
CC  ...Superlattice mode
C      if (plx(1,1) == NULLR .or. ddot(9,plx,1,plx,1) == 0) then
C        if (cmdopt('--plx=no',8,0,strn) .or. cmdopt('--noplx',7,0,strn)) then
C          strn(1:1) = '/'
C        elseif (cmdopt('--plx',5,0,strn)) then
C          strn = strn(7:)
C        else
C          write(stdo,'(/a)') ' Enter lattice vectors of supercell, or ? for help:'
C          strn = ' '; i = i1mach(1)
C          swtmp = rdstrn(i,strn,len(strn),F)
C        endif
C        call wordg(strn,101,', ',1,i,j)
C        if (i > j .or. strn(i:i) == '/') then
C          call dpzero(plx,9)
C        else if (strn(i:i) == '?') then
C          write(stdo,356)
C  356     format(/' Enter one of the following:'/
C     .      3x,'supercell lattice vectors (9 numbers), or'/
C     .      3x,'<enter> (take existing lattice vectors), or'/
C     .      3x,'m or p followed 9 integers ',
C     .      '(multiples of existing lattice vectors)'/,
C     .      3x,'a (or q) to abort lmscell')
C          call cexit(0,1)
C        else if (strn(i:i) == 'a' .or. strn(i:i) == 'q') then
C          call rx('aborting lmscell')
C        else if (strn(i:i) == 'm' .or. strn(i:i) == 'p') then
C          i = i+1
C          if (a2vec(strn,len(strn),i,2,', ',2,-3,9,ix,iplx) /= 9)
C     .      call rx('lmscell failed to read 9 values from string')
CC          do  j = 1, 3
CC          do  i = 1, 3
CC            plx(i,j) = plat(i,1)*iplx(1,j) +
CC     .                 plat(i,2)*iplx(2,j) +
CC     .                 plat(i,3)*iplx(3,j)
CC          enddo
CC          enddo
CC          print *, plx
CC          do  j = 1, 3
CC            plx(:,j) =  matmul(plat,iplx(:,j))
CC          enddo
C          plx = matmul(plat,iplx)
C        else
C          i = 0
C          if (a2vec(strn,len(strn),i,4,', ',2,-3,9,ix,plx) /= 9)
C     .      call rx('lmscell failed to read 9 values from string')
C        endif
C        if (ddot(9,plx,1,plx,1) == 0) call dcopy(9,plat,1,plx,1)
C      endif
C
C      call supcel(1,s_ctrl,s_lat,s_site,s_spec,s_pot,s_bz,s_move,slabl,plx,nkd,s_lat%dlv,nbx)
C      goto 1000
C#endif

C ... Input to scale MT radii (define RESIZE)
C#ifdef LM | LMCHK | LMPOL | LMPLAN | LMXBS | LMCOR | LMSTR | LMDOS
C#define RESIZE
C#endif
C#ifdefC LMGF | LMPG | LMF | LMFGWD | LMFGWS | LMFA | LMMC | LMSCELL
C#define RESIZE
C#endif

C --- Lattice setup ---
C#ifndef LMCTL | LMMC
      call setcg(s_lat,8,12)
      call lattic(s_lat,s_ctrl,s_site)
C#endif

C --- Generate symmetry operations; split species into classes  ---
      i = str_pack('symg',-2,s_strn,strn)
      if (strn == ' ') strn = 'find'
      if (cmdopt('--nosym',7,0,outs)) strn = ' '
      if (cmdopt('--socsym',7,0,outs)) s_lat%lsym = 1
      lc = 20
C#ifdefC TBE
C      mdprm = s_ctrl%mdprm
C      if (nint(mdprm(1)) == 1 .or. nint(mdprm(1)) == 2) then
C        strn = ' '
C      endif
C      lc = 40
C#endif
C#ifdefC MMAG
C      lc = 0
C#endif
C#ifdefC LMFA
CC      lc = 0
C#endif
C#ifdefC LMF
C      lc = 10
C#endif
C#ifndef LMFGW | LMFGWD | LMFGWS | LMFDMFT
      if (IAND(s_ctrl%lqp,1) == 0) lc = lc+2
C#endif
C#ifdefC LMF | LMFGWD | LMFGWS
C      lc = lc + 100
C#endif
C#ifndef LMMC
      call mksym(lc,slabl,strn,s_ctrl,s_lat,s_spec,s_site)
C#elseC
C      lc = 0
C#endif

C --- Allocate permanent class arrays, maps and other initialization ---
      if (mod(lc,100) >= 20)
     .  call clsprm(1,s_ctrl,s_ham,s_pot,s_spec,s_lat,s_bz,s_strn)

C --- Read available class parameters from file ---
C#ifdefC LMCTL | LMCHK
C      call aiocls(F,0,s_ctrl,s_ham,s_pot,s_spec,s_lat,1,0)
C#elseif LM | TBE | LMGF | LMSHF | LMPLAN
C#ifndef TBE
      call aiocls(F,12,s_ctrl,s_ham,s_pot,s_spec,s_lat,1,0)
C#endif
      call rdctrl(recrd,recln0,nrecs,prgnam,vrsion,ivn,2,slabl,
     .  s_bz,s_ctrl,s_ham,s_pot,s_lat,s_mix,s_spec,s_site,nbmx,
     .  s_str,s_move,s_tb,s_optic,s_gw,s_dmft,s_strn)
      if (nproc > 1) call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C#ifndef TBE
      call aiocls(F,17,s_ctrl,s_ham,s_pot,s_spec,s_lat,1,0)
      if (lc >= 20)
     .  call clsprp(1,s_ctrl,s_ham,s_pot,s_spec,s_lat,s_bz)
C#endif
C#endif
C     Discard contents of input file
      deallocate(recrd)
      if (procid == 0) call fclr('CTRL',-1)

C --- Optionally resize spheres ---
C#ifdef RESIZE
      strn = ' '; swtmp = cmdopt('--sfill',7,0,strn)
      call sfill(strn(8:),slabl,s_ctrl,s_lat,s_spec,s_site)
C#endif

C ... Set default values for species data
      call defspc(s_spec)

C ... Patch for now ... maybe replace
C#ifdef LM | TBE | LMGF | LMF
C#ifndef LMFGW | LMFGWD | LMFGWS
      ltet = IAND(s_ctrl%lmet,3) == 3 .or.
     .       IAND(s_ctrl%ldos,4+2+1) /= 0
      i = -2; if (mod(s_bz%lio,2) == 1) i = 0
      call mkqp(s_ctrl,s_bz,s_lat,ltet,F,1,i)
C#endif
C#endif

C ... quit after SHOW
      if (s_ctrl%quit == 1) then
        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif

C#ifdefC LMFA
C      call freeat(s_ctrl,s_ham,s_pot,s_spec,s_lat%tolft)
C#endif

C#ifdefC LMCHK & FP
C      if (cmdopt('--fp',4,0,strn)) then
C        call fpchk(s_spec,s_site)
C        call cexit(0,1)
C      endif
C#endif
C#ifdefC LMCHK | LMXBS | LMGPOL | LMSHF
C      if (cmdopt('--findes',8,0,strn) .or. cmdopt('--omax',6,0,strn)) auxmod = 128
C      call lmaux(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_pot,
C     .  s_str,s_spec,s_site,s_strn,slabl,auxmod)
C#endif

C#ifdef LM
      call lmasa(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_pot,s_str,s_spec,s_site,s_tb,s_optic,s_gw,s_strn)
C#endif

C#ifdefC LMFGWS
C      call sugws(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
C     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_strn)
C      goto 1000
C#endif

C#ifdefC LMF
C      call lmfp(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_move,s_pot,
C     .  s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
C#endif

C#ifdefC LMMC
C      if (cmdopt('--fit',5,0,strn)) then
C       call lmcfit('MCFIT',
C     .    s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,s_spec,
C     .    s_site,s_str,s_move,s_strn)
C      elseif (cmdopt('--atom',6,0,strn)) then
C        call lmca('MCA',
C     .    s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,s_spec,
C     .    s_site,s_str,s_move,s_strn)
C      elseif (cmdopt('--xbs',5,0,strn)) then
C        call lmxbs('MCXBS',
C     .    s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,s_spec,
C     .    s_site,s_str,s_move,s_strn)
C      elseif (cmdopt('--rho',5,0,strn)) then
C        call lmrho('MCXBS',
C     .    s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,s_spec,
C     .    s_site,s_str,s_move,s_strn)
C      else
C        call lmce('MCE',
C     .    s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,s_spec,
C     .    s_site,s_str,s_move,s_strn)
C      endif
C#endif

C#ifdefC LMSTR
C      call asastr(prgnam,s_ctrl,s_ham,s_lat,s_str,s_spec,s_site)
C#endif

C#ifdefC LMCTL
C      i = fopnn('LOG')
C      j = 0
C      if (cmdopt('-mad',4,0,strn)) j = 1
C      if (cmdopt('-spin1',6,0,strn)) j = j+10
C      if (cmdopt('-spinf',6,0,strn)) j = j+20
C      if (cmdopt('-enu',4,0,strn)) j = j+100
C      call shoctl(s_ctrl,s_spec,s_pot,j,i)
C      call fexit(0,0,' ',0)
C#endif

C#ifdefC LMGF
C      i = s_ctrl%lcgf
CC ... Special branch for exchange coupling
C      if (i >= 10 .and. i < 30) then
C        call exasa(s_bz,s_ctrl,s_ham,s_lat,s_pot,s_spec,s_site,
C     .    s_str,s_strn,slabl)
C      else
C        call lmasa(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
C     .    s_pot,s_str,s_spec,s_site,s_tb,s_optic,s_gw,s_strn)
C      endif
C#endif

C#ifdefC LMDOS
C      call asados(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_spec,s_site)
C#endif

C#ifdefC TBE
C      ltbe = 1
C      if (cmdopt('--band',6,0,outs)) ltbe = 2
C      if (cmdopt('--fit',5,0,outs)) ltbe = 4
C      call tbzint(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
C     .  s_pot,s_str,s_spec,s_site,s_tb,s_strn,ltbe)
C#endif

C#ifdefC LMIMP
C      call lmaux(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_pot,
C     .  s_str,s_spec,s_site,s_strn,slabl,2**6)
C#endif

C#ifdefC MMAG
CC     Uncomment these lines for Monte-Carlo atomistic simulations
CC     call mcasim(s_ctrl,s_spec,s_lat,s_move,s_strn)
C      call mmag(s_ctrl,s_ham,s_lat,s_move,s_spec,s_strn)
C#endif

C#ifdefC LMPLAN
C      call lmaux(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_pot,
C     .  s_str,s_spec,s_site,s_strn,slabl,2**1)
C#endif
C#endif

C -------------- End of program -------------
 1000 continue
      if (nproc > 1) then
C#ifdefC MPE
C      if (procid == master) i = fextg(ext)
C      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
C      i = 0
C      call skp2bl(ext,20,i)
C      ierr = MPE_FINISH_LOG(ext(2:i))
C#endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if ( procid == master ) then
        call rx0(prgnam//' on '//shortname(procid)(1:namelen(procid)))
      else
        call fexit(0,0,' ',0)
      endif
      else
      if (allocvb() /= 0) then
        i = idalloc(' ',1,1,1)
        if (i > 0) call info(30,1,0,' '//prgnam//'%a : maximum dynamic memory logged = %i MB',i,0)
      endif
      call rx0(prgnam)
      endif

      end
      subroutine lmvsn(ivn,strn)
C- Converts version into a string
      integer ivn(3)
      character strn*(*)
C     character cvn3*1
      if (ivn(3) == 0) then
        call awrit2('%a %i.%i',strn,len(strn),0,ivn(1),ivn(2))
      else
        call awrit3('%a %i.%i.%i',strn,len(strn),0,ivn(1),ivn(2),ivn(3))
C        cvn3 = char(ichar('a')+ivn(3)-1)
C        call awrit2('%a %i.%i.'//cvn3,strn,len(strn),0,ivn(1),ivn(2))
      endif
      end

      subroutine lmhelp(prgnam,ivn)
C- Help printout
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of main program
Ci   ivn   :program main version
Co Outputs
Co   message written to stdout
Cr Remarks
Cu Updates
Cu   11 Apr 03
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character prgnam*8
      integer ivn(3,2)
C ... Local parameters
      integer i1,i2
      character outs*1000

      call locase(prgnam)
      call info0(0,0,0,' usage:  '//prgnam//
     .  '%a [--OPTION] [-var-assign] [ext]')

      print 343
      print 344
  343 format(/' --h'/' --help',t17,'Print this message, and quit'
     .  /' --input',t17,
     .  'List categories, tokens, and data program expects, and quit'
     .  /' --show',t17,
     .  'Print control file after parsing by preprocessor,'/t17,
     .  'and echo input data as read from the control file'
     .  /' --showp',t17,
     .  'Same as --show, but quit after input parsed'
     .  /' --iactiv',t17,'(--no-iactiv) ',
     .  'Turn on (off) interactive mode'/
     .  t17,'This switch overrides input file setting',
     .  /' --pr#1[,#2...]',t17,
     .  'Set the verbosity (stack) to values #1,#2, ...'
     .  /' --time=#1[,#2]',t17,
     .  'Print timing info to # levels (#1=summary; #2=on-the-fly)'
     .  /' --shomem',t17,
     .  'Print memory allocation of some arrays (equivalent to IO_SHOMEM=1)'/
     .  /' -vnam=expr',t17,
     .  'Define numerical variable "nam"; set to result of ''expr'''
     .  /' -cnam=strn',t17,
     .  'Define character variable "nam"; set to ''strn'''
     .  )


  344 format(
     .  /' --rpos=filnam',t17,
     .  'After reading input file, read site positions from "filnam"'
     .  /' --fixlat',t17,'Adjust lattice vectors and point group ops, attempt to'
     .  /t17,'render them internally consistent'
     .  /' --fixpos[:tol=#]', ' Adjust positions slightly, rendering them'
     .  /t17,'as consistent as possible with the symmetry group'
     .  /' --nosym',t17,'Suppress symmetry operations')

      if (.true.) then

        outs = '%N '//prgnam//'%a-specific options:'
        call strip(outs,i1,i2)
        call info0(0,0,0,'%N '//prgnam//'%a-specific options:')

C#ifdefC LMFA
C        if (prgnam == 'lmfa') then
C          call info0(0,0,0,
C     .      '%N%1f--noopt        Suppress optimization of s.m. Hankel basis'//
C     .      '%N%1f--norscnst     In optimization of s.m. Hankel basis, do not constrain rsm < rmt'//
C     .      '%N%1f--plotwf       Writes atomic radial wave functions to disk files'//
C     .      '%N%1f--dumprho      Writes out the density for each atom to out.ext'//
C     .      '%N%1f--basp[~opts]  Turns on autofind EH,RSMH (better to use HAM_AUTOBAS)'//
C     .      '%N%1f               Options are :  ctrl  eh=#  eh2=#  rsmmx  incrlmx'//
C     .      '%N%1f--basfile=fn   write the autogenerated basis set parameters to file fn instead of default basp0'//
C     .      '%N%1f--getallloc    Look for local orbitals (better to use HAM_AUTOBAS)')
C        endif
C#endif

C#ifdefC LMFGWD
C        if (prgnam == 'lmfgwd') then
C          call info0(0,0,0,
C     .      ' --jobgw=#'//
C     .      '%6f-2 check GWinput -1 create GWinput  0 init mode  1 standard mode'//
C     .      '%N --job=#%8fsame as --jobgw=#'//
C     .      '%N --vxcsig'//
C     .      '%6f Write qsgw sigma in place of LDA vxc into file vxc.ext'//
C     .      '%N --shorbz=no'//
C     .      '%3f Suppress shortening of q'//
C     .      '%N   ...%9f the following apply to GWinput creation mode (--jobgw=-1):'//
C     .      '%N --gwin0%7f Read parts of template from GWin0'//
C     .      '%N --sigw'//
C     .      '%8f Add lines to GWinput for creating Sigma(omega).%N'//
C     .      '%15f Note: does not affect hsfp0_sc, but hsfp0 must be run with job 4'//
C     .      '%N --ib=list'//
C     .      '%5f list of QP levels for which to make sigma (for 1-shot)'//
C     .      '%N --make-Q0P%4f Make the Q0P file'//
C     .      '%N%1f ')
C        endif
C#endif

C#ifdefC LMFGW
C        if (prgnam == 'lmfgw') then
C          call info0(0,0,0,
C     .      'Under construction')
C        endif
C#endif

C#ifdefC LMFGWS
C        if (prgnam == 'lmfgws') then
C          call info0(0,0,0,
C     .      ' --sfuned'//
C     .      '%N%1f ')
C        endif
C#endif

C#ifdefC LM67
C        if (prgnam == 'lm67') then
C          call info0(0,0,0,
C     .      '%N%1f--prog=[LM|LMF|...] ')
C        endif
C#endif

C#ifdefC LMF
C        if (prgnam == 'lmf') then
C          call info0(0,0,0,
C     .      '%N%1f--rs=#1,#2,#3,#4,#5'//
C     .      '%N%6f#1=0 start from atm file; 1 from rst file;'//
C     .      ' 2 from rsta file'//
C     .      '%N%11fadd 10 to shift sm-rho 1st iter'//
C     .      '%N%11fadd 100 to rotate local rho 1st iter'//
C     .      '%N%6f#2=1 save rst file'//
C     .      '%N%6f(#3,#4,#5)=0 read (pos,E_f,pnu) from rst file')
C          call info0(0,0,0,
C     .      ' --shorten=no --shorbz=no --no-fixef0 --rdbasp --symsig[=no]'//
C     .      ' --wratrho --rdatrho[~v0]'//
C     .      '%N --rhopos --wrhomt --wpotmt --window=#,# --wden --etot'//
C     .      ' --ef=# --efrnge --minmax --zhblock'//
C     .      '%N   Special purpose modes:'//
C     .      '%N --optbas[:wbas][:sort][:etol][:spec=..] | --band[~options] | --wden[~options]'//
C     .      '%N   Switches to generate extra information:'//
C     .      '%N --pdos[~options] | --cls[~options] | --mull[~options]'//
C     .      ' | --jdosw=lst [--jdosw2=lst2]'//
C     .      ' --opt:read --opt:write'//
C     .      ' --wrhoat[:l=#][:lx=#]'//
C     .      '%N   Switches that invoke editors:'//
C     .      '%N --chimedit[~..] | --rsedit[~..] | --popted[~..] |'//
C     .      ' --wsig[~edit~..]'
C     .      )
C        endif
C#endif

C#ifdefC LMMC
C        if (prgnam == 'lmmc') then
C          call info0(0,0,0,
C     .      '%N%1f--atom --fit --rs --st'//
C     .       '%N --atom invokes the free atom program'//
C     .       '%N --fit  invokes the two-center fit'
C     .      )
C        endif
C#endif

C#ifdef LM
        if (prgnam == 'lm') then
          call info0(0,0,0,
     .      '%N%1f--rs=#1,#2 --band[~option...] --pdos[~option...] --rsedit[~option...] --zerq[~option...]'//
     .      '%N%1f--efrnge -mix=#1[,#2] --onesp --weula'//
     .      '%N%1fSee https://www.questaal.org/docs/input/commandline/#switches-for-lm')
        endif
C#endif

C#ifdefC LMGF
C        if (prgnam == 'lmgf') then
C          call info0(0,0,0,
C     .      '%N%1f--rs=#1,#2 --band[~option...] --pdos[~option...] --zerq[~option...]'//
C     .      '%1f--ef=#'//
C     .      '%N%1fExchange-mode-specific options:'//
C     .      '%N  --sites[:pair]:site-list'//
C     .      '%N  --wrsj[:j00][:amom][:sscl][:[g]scl=#][:tol=#] --wmfj '//
C     .      '--rcut=# --tolJq=# --amom[s]=mom1,mom2,... --2xmsh'//
C     .      '%N See https://www.questaal.org/docs/input/commandline/#switches-for-lmgf')
C        endif
C#endif

C#ifdefC LMDOS
C        if (prgnam == 'lmdos') then
C          call info0(0,0,0,
C     .      '%N%1f--dos~options  modifies number and kinds of dos generated.'//
C     .      '%N%16fSee Questaal command-line arguments documentation')
C          call info0(0,0,0,
C     .      ' --mull~options tells lmdos to read channels from a Mulliken fpopulation analysis.'//
C     .      '%N%16fSee Questaal command-line arguments documentation')
C          call info0(0,0,0,
C     .      ' --pdos~options tells lmdos to read channels from a partial wave decomposition of DOS.'//
C     .      '%N%16fSee Questaal command-line arguments documentation')
C          call info0(0,0,0,
C     .      ' --cls          tells lmdos to use file cls as weighting function')
C        endif
C#endif

C#ifdefC LMCHK
C        if (prgnam == 'lmchk') then
C          call info0(0,0,0,
C     .      '%N%1f--shell[:v][:e][:r=#][:sites:site-list]'//
C     .      '[:pairs:pair-list]...'//
C     .      '%N%8f...[:tab[=#]][:disp=fnam][:nn][:fn=fnam]'//
C     .      '%N --mino[:dxmx=#][:xtol=#][:maxit=#][:style=#]:list'//
C     .      '%N --getwsr find augmentation sphere radii'//
C     .      '%N --findes [--nescut=# --mino --wsite|--wsitex] locate possible sites for empty spheres'//
C     .      '%N --wpos=fnam  write site positions to file'//
C     .      '%N --syml~args  write symmetry lines file'//
C     .      '%N --omax=#     Override given SPEC_OMAX1'//
C     .      '%N --angles[~r=rmax][~sites=site-list]'//
C     .      '%N --euler[~r=rmax[,rmin]][~sign][~sites=site-list]'//
C     .      '%N --wsite[x][~short]~fn=fnam'//
C     .      '%N --terse')
C        endif
C#endif

C#ifdefC LMSCELL
C        if (prgnam == 'lmscell') then
C          call info0(0,0,0,
C     .    '%N  --stack[options~] invokes the superlattice editor'
C     .  //'%N  --wsite[x][~map][~short][~quad1][~fn=fnam]'
C     .  //'%N  --sort:plat | --sort:"expr [expr2] [expr3]"'
C     .  //'%N  --plx[~m]~#1,..,#9 | noplx'
C     .  //'%N  --rsta[,amom]'
C     .  //'%N  --ring:i1,i2 | swap:i1,i2[,i3,i4]'
C     .  //'%N  --sites:site-list'
C     .  //'%N  --shorten'
C     .  //'%N  --pl:expr'
C     .  //'%N  --sqs[~seed=#][~r2max=#][~r3max=#][~r3mode=#]'
C     .  //'%N  --wrsj[:fn=name][:scl=#]'
C     .  //'%N  --disp:fname:site-list'
C     .  //'%N  --xshft[x]=#,#,# Shift site positions by a uniform translation')
C
C        endif
C#endif

C#ifdefC LMSTR
C        if (prgnam == 'lmstr') then
C          call info0(0,0,0,
C     .      '%N%1f--chk%11fcompares file str.ext with str1.ext'//
C     .       '%N --plot[:con|:line[,v1x..z,v2x..z]|onec] '//
C     .       'plots envelope'//
C     .       '%N --pltg[:con|:line[,v1x..z,v2x..z]|onec] '//
C     .       'plots envelope, val-lap')
C        endif
C#endif

C#ifdefC LMXBS
C        if (prgnam == 'lmxbs') then
C          call info0(0,0,0,
C     .      '%N%1f-shift=x1,x1,x3 -spec -sites=list '//
C     .      '-dup=d1,d2,d3[,expr] -bs=val -ss=val')
C        endif
C#endif

C#ifdefC LMCTL
C        if (prgnam == 'lmctl') then
C          call info0(0,0,0,
C     .      '%N%1f -spin1 -spinf -mad -enu')
C        endif
C#endif

C#ifdefC LMPG
C        if (prgnam == 'lmpg') then
C          call info0(0,0,0,
C     .      '%N%1f -map -onesp')
C        endif
C#endif

C#ifdefC LMPLAN
C        if (prgnam == 'lmplan') then
C          call info0(0,0,0,
C     .      '%N%1f--pledit[~batch-editor-instructions]'//
C     .      '%N%1f---')
C        endif
C#endif

C#ifdefC LMSHF
C        if (prgnam == 'lmshf') then
C          call info0(0,0,0,
C     .      '%N%1f-enu=expr   linearize pot. pars around ''expr''')
C        endif
C#endif

C#ifdefC LMIMP
C        if (prgnam == 'lmimp') then
C          call info0(0,0,0,
C     .      '%N%1f -rs -4 -5 -3s -4s -47u -5s (1 is required)')
C        endif
C#endif

C#ifdefC LMCOR
C        if (prgnam == 'lmcor') then
C          call info0(0,0,0,
C     .      '%N%1f -findr')
C        endif
C#endif

C#ifdefC LMAVGM
C        if (prgnam == 'lmavgm') then
C          call info0(0,0,0,
C     .      '%N%1f -spin1')
C        endif
C#endif

C#ifdefC LMMIX
C        if (prgnam == 'lmmix') then
C          call info0(0,0,0,
C     .      '%N%1f-fn=mix-file-name -bin2a or -a2bin')
C        endif
C#endif

C#ifdefC TBE
C        if (prgnam == 'tbe') then
C          call info0(0,0,0,
C     .      '%N%1f--band[~option...] --wpos=fnam -cont -dumph'//
C     .      ' --st --md=# --mv=# --xyz=#'
C     .    //'%N%3fThe last 4 switches apply to'
C     .    //' molecular dynamics simulations'
C     .    //'%N%N  Mixing options'
C     .    //'%N --mxq          Mix multipole moments and magnetic'
C     .    //' moments togther: beta'
C     .    //'%N --mxq2         Mix multipole moments and magnetic'
C     .    //' moments separately: beta, betav'
C     .    //'%N --mxr          Mix Hamiltonian'
C     .    //'%N%N  MPI process mapping options. In most cases the'
C     .    //' default mapping is the optimal choice'
C     .    //'%N --kblocks=#    Number of rectangular processor blocks'
C     .    //' over which to distribute the k-points'
C     .    //'%N --nprows=#     Number of rows of the ScaLAPACK/BLACS'
C     .    //' process arrays/blocks.')
C
C          call info0(0,0,0,
C     .    '%N  Linear scaling playground options'
C     .    //'%N --linscl {msi|ls++}'
C     .    //'%N                Switch to linear scaling branch:'
C     .    //'%N                  msi : matrix sign'
C     .    //' iterations/polynomial mode'
C     .    //'%N                  ls++: use the LS++ library (MPI C++),'
C     .    //'from T. Kuehne et. al.'
C     .    //'%N --bgc          Set expected band gap centre'
C     .    //'%N --lstol        Set tolerance/accuracy for iterations'
C     .    //'%N                and threshold for values stored in'
C     .    //' sparse matrices (second argument, if present)'
C     .    //'%N --lsxxp        the LS++ p parameter (ls++ branch only)'
C     .    //'%N --lsxxb        the LS++ beta parameter (ls++ branch'
C     .    //' only)')
C        endif
C#endif

C#ifdefC MMAG
C        if (prgnam == 'mmag') then
C          call info0(0,0,0,'%N%1f--cont --wrsj[:fn=name]')
C        endif
C#endif

C#ifdefC BLM
C        if (prgnam == 'blm') then
C          call info0(0,0,0,
C     .      '%N --express[=n]%17pExpress style input file, level n (default=6 if n missing, 3 if switch missing)'//
C     .      '%N --asa%10ftailor input file to ASA calculation'//
C     .      '%N --gw[~sw]%6ftailor input file to GW calculation'//
C     .      '%N --gf%11fadd tokens for lmgf input file'//
C     .      '%N --pgf%10fadd tokens for lmpg input file'//
C     .      '%N --fpandasa%5ftags for both ASA and FP')
C          call info0(0,0,0,
C     .      '%N --gmax=#%7fassign # to plane wave cutoff GMAX.  No default is given for this number'//
C     .      '%N --nk=#[,#,#]%3fk-mesh for BZ integration.  No default is given for this number'//
C     .      '%N --nkgw=#[,#,#]%1fGW k-mesh.  Same as --gw~nk=...'//
C     .      '%N --nl=#%9fSet global maximum (l+1). Fixes lmxa --gw is also present'//
C     .      '%N --gwemax=#%5fsigma cutoff SIG_GMAX.  Same as --gw~emax='//
C     .      '%N --mag%10fset nsp=2 for spin polarized calculation'//
C     .      '%N --nit=#%8fspecify max number of iterations to self-consistency'//
C     .      '%N --conv=#%7fspecify tolerance for energy convergence (ITER_CONV)'//
C     .      '%N --convc=#%6fspecify tolerance for charge convergence (ITER_CONVC)'//
C     .      '%N --dhftol=#%5fset tolerance correction to Harris forces (ITER_DHFTOL)'//
C     .      '%N --convp=#%6fconv, convc, dhftol are assigned variables, controllable from command line'//
C     .      '%N --pbe%10fset switches for PBE functional (libxc)'//
C     .      '%N --ldau%9fAdd tags for ldau.  Absent lst, tags are added for each species.'//
C     .      '%N   %12f Use --ldau:spec=lst to add tags for species named in lst'//
C     .      '%N --pbesol%7fset switches for PBESOL functional (libxc)')
C          call info0(0,0,0,
C     .        ' --loc=#%8fset AUTOBAS_LOC'//
C     .      '%N --eloc=#%7fset AUTOBAS_ELOC'//
C     .      '%N --mto=#%8fset AUTOBAS_MTO'//
C     .      '%N --pmt[:emax=#] adjustments to prepare for the PMT basis'//
C     .      '%N --optics       template for OPTICS category'//
C     .      '%N --dv=string    variable declarations, e.g. --dv=idp=0'//
C     .      '%N --molstat      Add tag for molecular statics (lmf)'//
C     .      '%N --lmfit        Add tag for Levenberg-Marquardt fitting (lm)'//
C     .      '%N --noshorten%4fsuppress shortening of site positions')
C          call info0(0,0,0,
C     .      '%N ... the following concern augmentation spheres'//
C     .      '%N --omax=#[,#,#] set maximum sphere overlap to # (# may be negative).  0 is touching.'//
C     .      '%N --addes%8fadd tags to prepare for later addition of empty spheres'//
C     .      '%N --ehmx=#%7fset upper bound to LMTO sm-H energy (HAM_AUTOBAS_EHMX)'//
C     .      '%N --rsmmx=#%6fset upper bound to LMTO smoothing radius (HAM_AUTOBAS_RSMMX)'//
C     .      '%N --findes[~opt]%1fsearch for empty sphere sites to improve lattice packing'//
C     .      '%N %15ffor options, see Questaal command-line switches page'//
C     .      '%N --clfloat%6ffloat Pnu with traditional algorithm (lmf only)'//
C     .      '%N --frzwf%8fadd switches that can keep the basis set fixed')
C          call info0(0,0,0,
C     .      '%N ... the following concern how structural data is (read from or) written to disk')
C          call info0(0,0,0,
C     .      '%N --ctrl=fn%6fname the input file fn (default = actrl).  blm supplies extension.'//
C     .      '%N --rdsite[=fn]  read structural data from site file fn (default=sitein).  No extension appended.'//
C     .      '%N --scala=#%6fmultiply ALAT by #, and PLAT and POS by 1/#'//
C     .      '%N --scalp=#%6flike scala, but # specifies PLAT volume'//
C     .      '%N --xshft=#,#,#  Shift site positions by a uniform translation, Cartesian coordinates'//
C     .      '%N --xshftx=#,#,# Shift site positions by a uniform translation, Crystal coordinates'//
C     .      '%N --xpos%9fwrite site positions as multiples of PLAT'//
C     .      '%N --wsrmax=#%5fwhen finding sphere radii, do not permit value to exceed #. Default=3.3.'//
C     .      '%N%1f--wsite[~opt]%2fwrite site file for structural data'//
C     .      '%N%1f--wsitex[~opt]%1fsame as --wsite, but write positions in crystal coordinates'//
C     .      '%N --wpos=fn%6fwrite site positions to file fn.ext'//
C     .      '%N --nfile%8fAdd tags to control site file name with variable ''file'''//
C     .      '%N --molstat%6fAdd template for molecular statics')
C        endif
C#endif

C#ifdefC LMFDMFT
C          call info0(0,0,0,
C     .      '%N --ldadc=#         Set double counting to #'//
C     .      '%N --udrs            Generate and save output density'//
C     .      '%N --makesigqp       makes quasiparticle DMFT sigma, in format of QSGW sigma file'//
C     .      '%N --gprt[~options]  Write gloc to disk, with the following possible options:'//
C     .      '%N     ~rdsigr=fn    Read a new mesh of (frequencies, DMFT sigma) from file fn.ext'//
C     .      '%N                   Frequencies assumed to be real'//
C     .      '%N     ~rdsig=fn     Read a new mesh of (frequencies, DMFT sigma) from file fn.ext'//
C     .      '%N                   Frequencies assumed to be imaginary'//
C     .      '%N     ~mode=#       #=1:  local (default)'//
C     .      '%N                   #=2:  k-resolved gkloc; #=3: combine #1 and #2'//
C     .      '%N                   #=4:  write hkloc'//
C     .      '%N                   #=8:  write G'//
C     .      '%N                   #=18: write diagonal G(k,omega) into se file'//
C     .      '%N                   #=19: write diagonal sigma(i,omega) into se file'//
C     .      '%N     ~band@flags   Generate property along user-specified list of k-points'//
C     .      '%N     ~nom=#        number of Matsubara frequency points (defaults to all available points)')
C
C#endif

      endif

!       outs = '%N version '//prgnam
      outs = '%N '//trim(prgnam)//' v'
      call lmvsn(ivn(1,2),outs)
      call info0(0,0,0,trim(outs))
      call fexit(0,0,' ',0)
      end
      subroutine sttmpd
C- Creates special tmpdir for saving temporary files
C  User may wish to customize this routine.
      implicit none
      character tmpdir*100
      integer fopnT,ifi
C ... for henry, lm-MPIK
C      integer i1,i2,nw
C      character*40 strn

C     return
C ... Set customization of temporary directory here, if desired
C     This is usual default (current working directory)
      tmpdir = ' '
C     call gtenv('HOME',tmpdir)
C     call gtenv('TMPDIR',tmpdir)
C     call getenv('HOME',tmpdir)
C     call getenv('TMPDIR',tmpdir)
C     tmpdir = '/home/tmp/'

C ... Marco's cluster, lm-MPIK specific
C     tmpdir = '/local/svan2/tmp'

C ... for henry, lm-MPIK specific
C      call getenv('HOME',strn)
C      call strip(strn,i1,i2)
C      call wrdsg(strn(i1:i2),0,'/',nw)
C      call wordg(strn,0,'/',nw,i1,i2)
C      strn = strn(i1:)
C      if (strn == 'markv') strn = 'svan2'
C      call word(strn,1,i1,i2)
C      tmpdir = '/home/' // strn(i1:i2) // '/tmp'

C ... Set the directory
      ifi = fopnT(tmpdir,0,0,11)

C     debugging check
C      ifi = fopnT('tmp' ,-1,0,0)
C      call fshow
C      print *, ifi
C      write(ifi,*) 'test'
C      call rx('done')
      end
