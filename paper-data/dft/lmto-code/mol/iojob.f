      logical function iojob(jobs,ifi)
C- i/o jobs category
      integer jobs(4), ifi
      logical scat,sw,rdstrn
      integer nlmx,i,j,ifi2,partok
      parameter (nlmx=500)
      character*80 t(nlmx),strn

      sw = scat(iabs(ifi),'JOB ',' ',.true.)
      if (ifi > 0) then
        iojob = .false.
        if (.not. sw) return
        backspace ifi
        if (.not. rdstrn(ifi,strn,72,.false.)) return
        call stlibv(1,0,.false.)
        call stlibv(3,1,.false.)
        call stlibv(12,1,.false.)
        call getcat(strn,'JOB ',' ',.true.)
        i = partok(strn,'JOB ',' ',jobs,' ',-4,2,0,.true.)
        iojob = .true.
      else
        ifi2 = -ifi
C Save everything after JOB line
        if (.not. sw) call poseof(ifi2)
        j = 0
        do  1  i = 1, nlmx
          read(ifi2,100,end=91) t(i)
  100     format(a80)
          j = i
    1   continue
        call fexit(0,9,' iojob: too many output lines',0)
   91   continue
C Overwrite JOB line
        if (scat(ifi2,'JOB ',' ',.true.)) then
          backspace ifi2
        else
          call poseof(ifi2)
        endif
        write(ifi2,'(''JOB   '',4i4)') jobs
C Restore rest of file
        do  2  i = 1, j
    2   write(ifi2,100) t(i)
        iojob = .true.
      endif
      end
