
   module mod_proc_stat
! Get process stats from the Linux procfs.
! Needless to say this is only working in Linux and could be useful for
! debugging rather than anything else.

   use iso_c_binding, only : c_size_t

   implicit none

   private
   integer, parameter ::lp = 8


   type proc_stat
      integer(lp)          :: pid
      character(len=100)   :: comm
      character            :: state
      real(lp)          :: ppid
      real(lp)          :: pgrp
      real(lp)          :: session
      real(lp)          :: tty_nr
      real(lp)          :: tpgid
      real(lp)          :: flags
      real(lp)          :: minflt
      real(lp)          :: cminflt
      real(lp)          :: majflt
      real(lp)          :: cmajflt
      real(lp)          :: utime
      real(lp)          :: stime
      real(lp)          :: cutime
      real(lp)          :: cstime
      real(lp)          :: priority
      real(lp)          :: nice
      real(lp)          :: num_threads
      real(lp)          :: itrealvalue
      real(lp)          :: starttime
      real(lp)          :: vsize
      real(lp)          :: rss
      real(lp)          :: rsslim
      real(lp)          :: startcode
      real(lp)          :: endcode
      real(lp)          :: startstack
      real(lp)          :: kstkesp
      real(lp)          :: kstkeip
      real(lp)          :: signal
      real(lp)          :: blocked
      real(lp)          :: sigignore
      real(lp)          :: sigcatch
      real(lp)          :: wchan
      real(lp)          :: nswap
      real(lp)          :: cnswap
      real(lp)          :: exit_signal
      real(lp)          :: processor
      real(lp)          :: rt_priority
      real(lp)          :: policy
      real(lp)          :: delayacct_blkio_ticks
      real(lp)          :: guest_time
      real(lp)          :: cguest_time
   end type

   public :: proc_stat, read_proc_stat, print_proc_stat

   contains

   subroutine read_proc_stat(ps)

      type(proc_stat), intent(inout) :: ps

      integer :: u, ierr

      open(unit=u, file='/proc/self/stat', form='formatted', action='read')
!     open(newunit=u, file='/proc/self/stat', form='formatted', action='read')

!       print *, 'newunit: ', u

      read(u,*) ps%pid, ps%comm, ps%state, ps%ppid, ps%pgrp, ps%session, &
         & ps%tty_nr, ps%tpgid, ps%flags, ps%minflt, ps%cminflt, ps%majflt, &
         & ps%cmajflt, ps%utime, ps%stime, ps%cutime, ps%cstime, ps%priority, &
         & ps%nice, ps%num_threads, ps%itrealvalue, ps%starttime, ps%vsize, ps%rss, &
         & ps%rsslim, ps%startcode, ps%endcode, ps%startstack, ps%kstkesp, ps%kstkeip, &
         & ps%signal, ps%blocked, ps%sigignore, ps%sigcatch, ps%wchan, ps%nswap, ps%cnswap, &
         & ps%exit_signal, ps%processor, ps%rt_priority, ps%policy, ps%delayacct_blkio_ticks, &
         & ps%guest_time, ps%cguest_time

      close(u)

   end subroutine read_proc_stat


   subroutine print_proc_stat(u, ps)
      integer, intent(in) :: u
      type(proc_stat), intent(in) :: ps

      write(u,*) 'pid                   : ', ps%pid
      write(u,*) 'comm                  : ', ps%comm
      write(u,*) 'state                 : ', ps%state
      write(u,*) 'ppid                  : ', ps%ppid
      write(u,*) 'pgrp                  : ', ps%pgrp
      write(u,*) 'session               : ', ps%session
      write(u,*) 'tty_nr                : ', ps%tty_nr
      write(u,*) 'tpgid                 : ', ps%tpgid
      write(u,*) 'flags                 : ', ps%flags
      write(u,*) 'minflt                : ', ps%minflt
      write(u,*) 'cminflt               : ', ps%cminflt
      write(u,*) 'majflt                : ', ps%majflt
      write(u,*) 'cmajflt               : ', ps%cmajflt
      write(u,*) 'utime                 : ', ps%utime
      write(u,*) 'stime                 : ', ps%stime
      write(u,*) 'cutime                : ', ps%cutime
      write(u,*) 'cstime                : ', ps%cstime
      write(u,*) 'priority              : ', ps%priority
      write(u,*) 'nice                  : ', ps%nice
      write(u,*) 'num_threads           : ', ps%num_threads
      write(u,*) 'itrealvalue           : ', ps%itrealvalue
      write(u,*) 'starttime             : ', ps%starttime
      write(u,*) 'vsize                 : ', ps%vsize
      write(u,*) 'rss                   : ', ps%rss
      write(u,*) 'rsslim                : ', ps%rsslim
      write(u,*) 'startcode             : ', ps%startcode
      write(u,*) 'endcode               : ', ps%endcode
      write(u,*) 'startstack            : ', ps%startstack
      write(u,*) 'kstkesp               : ', ps%kstkesp
      write(u,*) 'kstkeip               : ', ps%kstkeip
      write(u,*) 'signal                : ', ps%signal
      write(u,*) 'blocked               : ', ps%blocked
      write(u,*) 'sigignore             : ', ps%sigignore
      write(u,*) 'sigcatch              : ', ps%sigcatch
      write(u,*) 'wchan                 : ', ps%wchan
      write(u,*) 'nswap                 : ', ps%nswap
      write(u,*) 'cnswap                : ', ps%cnswap
      write(u,*) 'exit_signal           : ', ps%exit_signal
      write(u,*) 'processor             : ', ps%processor
      write(u,*) 'rt_priority           : ', ps%rt_priority
      write(u,*) 'policy                : ', ps%policy
      write(u,*) 'delayacct_blkio_ticks : ', ps%delayacct_blkio_ticks
      write(u,*) 'guest_time            : ', ps%guest_time
      write(u,*) 'cguest_time           : ', ps%cguest_time

   end subroutine print_proc_stat

   end module mod_proc_stat






