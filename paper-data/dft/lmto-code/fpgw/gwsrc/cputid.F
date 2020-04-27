      subroutine cputid(ifile)
      implicit none
c cpu time in seconds and minutes
c ifile = file number where the cpu time information is
c         to be printed
c ifile = 0 ==> screen
      integer :: ifile
      character*(*):: message

      integer :: u
      real(8) :: cpu, cpuetime, cpuwtime, cpuw, cpumin
      real(8), save :: cpu0=-1.0d0, cpulast=0.0d0, cpuw0=0.0d0, cpuwlast=0.0d0

      procedure(real(8)) :: cpusec, dwtime

      if (ifile < 0) call rx('cputid: negative file logical unit')

      cpuetime = cpusec(); cpuwtime = dwtime()
      if (cpu0 < 0.d0) then
        cpu0 = cpuetime
        cpuw0 = cpuwtime
        cpulast = 0
        cpuwlast = 0
      endif
      cpu   = cpuetime - cpu0
      cpuw  = cpuwtime - cpuw0
      cpumin   = cpu/60.0d0
      u = ifile
      if(u == 0) u = 6
      write(u,"(' CPU= ',f12.3,' secs',f12.3,' mins','.   Since last call:',
     &   f13.6,' sec',f13.5,' (wall)')") cpu,cpumin,cpu-cpulast,cpuw-cpuwlast
      cpulast = cpu
      cpuwlast = cpuw
      return

      entry cputid2(message,ifile)

      if (ifile < 0) call rx('cputid2: negative file logical unit')
      cpuetime = cpusec(); cpuwtime = dwtime()
      if (cpu0 < 0.d0) then
        cpu0 = cpuetime
        cpuw0 = cpuwtime
        cpulast = cpu0
        cpuwlast = 0
      endif
      cpu   = cpuetime - cpu0
      cpuw  = cpuwtime - cpuw0
      cpumin   = cpu/60.d0
      u = ifile
      if(u == 0) u = 6
      write(u,"(a,' CPU= ',f12.3,' secs',f12.3,' mins','.   Since last call:',
     & f13.6,' sec',f13.5,' (wall)')") trim(message),cpu,cpumin,cpu-cpulast,
     & cpuw-cpuwlast
      cpulast = cpu
      cpuwlast = cpuw

      return
      end

C      module m_realtimediff
C      implicit none
C      integer,parameter:: nmax=10000
C      integer:: i1(nmax)=nmax*0
C      character(30):: msg(nmax)=''
C      real(8):: elapse(nmax)=nmax*0.0d0
C      real(8):: start(nmax)=nmax*0.0d0
C      end module m_realtimediff
C
C
C      subroutine realtimediff(id,msg0)
C      use m_realtimediff
C      implicit none
C      integer,intent(in):: id
C      character(*),intent(in)::msg0
C
C      integer::mode
C      real(8):: diff
C      integer:: i2,imax
C      character(20):: errmsg
C#define USEMPIWTIME
C#ifdef USEMPIWTIME
C      real(8):: mpi_wtime
C
C      if (id<=0 .or. id>nmax) then
C          write(errmsg,"(2i5)") id,nmax
C          call rx('fatal error, id<=0 or id>nmax, id,nmax='//trim(errmsg))
C      endif
C
C      mode=0
C      if (len_trim(msg0)>0) mode=1
C
C      if (mode == 0) then
C        i1(id)=i1(id)+1
C        start(id)=mpi_wtime()
C      else
C        diff=mpi_wtime()
C        msg(id)=msg0
C        elapse(id)=elapse(id)+diff-start(id)
C      endif
C#else
C      if (mode == 0) then
C        call system_clock(i1(id))
C      else
C        msg(id)=msg0
C        call system_clock(i2,irate,imax)
C        diff = i2-i1(id)
C        if (i1(id)>i2) diff=imax-i1(id)+i2
C        diff= dble(i2-i1(id))/dble(irate)  ! real time
C       ! write(*,*)'mpi=',str_rank,trim(msg),diff
C        elapse(id)=elapse(id)+diff
C      endif
C#endif
C
C      end subroutine realtimediff
C
C      subroutine print_realtimediff()
C      use m_mpi
C      use m_realtimediff
C      implicit none
C      integer:: i,j,isize
C      character(30):: str
C         do i=1,nmax
C           !if (i1(i) /= 0 .and. len_trim(msg(i))>0 ) then
C           if (i1(i) /= 0 ) then
C           str=msg(i)
C           do j=1,30
C             if (str(j:j) == ' ') str(j:j)='_'
C           enddo
C!           write(6,'(a,i4,1x,a,1x,F10.2,1x,a,1x,i5,1x,a)')'timediff,id=',i,'('//trim(str)//')',elapse(i),'sec' ,i1(i),'times'
C            write(6,*)'timediff,id=',i,trim(str),elapse(i),'sec', i1(i),'times'
C           endif
C         enddo
C      end subroutine print_realtimediff
