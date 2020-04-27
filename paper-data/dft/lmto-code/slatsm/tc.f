      subroutine tcinit(ltrace,lprof)
C- Setup for timer and workspace routines
C ----------------------------------------------------------------------
Ci Inputs
Ci   ltrace : >0 turn on running account of cpu usage
Ci            The value of ltrace is the depth to which quantities shown
Ci            <0 also turn on running account of workspace usage
Ci   lprof  : >0 print out summary of cpu and workspace usage
Ci            The value of lprof is the depth to which quantities shown
Cr Remarks
Cr   This routine logs cpu and workspace usage.
Cu Updates
Cu   17 May 14 (DP) drop w and print from one mpi process only
Cu        2012 (DP) collect walltime only instead of cputime
Cu   30 Jun 06 ltrace <0 => running workspace printout
Cu   15 Feb 02 tcprt now writes to a logical unit, rather than stdout
C ----------------------------------------------------------------------
      implicit none
      integer ltrace,lprof
      integer i,it,jcb,jlev,jt0,l,levx,lok,ltop,job
      integer, parameter :: nx = 200
      character*20 name(nx),nam0
      character*(*) str1,str2
      character*3  tt(0:20)
      integer :: ncall(nx),lev(nx),it1(nx),icb(nx),icall(0:30)
      integer iwk0(nx),iwk1(nx),iwk2(nx),iwk3(nx),iwmax
      integer :: lvl,lv,last,j1,j2,ilev
      integer ifi,iprint,err,stdo,lgunit
      logical :: lrt
      integer :: omppid, mpipid
      real :: ttot(nx),tavg,tm
      double precision cpusec
      character(len=83), parameter :: hdr = "("//
     ."'     calls      === cpu time ===  depth',i2,/,"//
     ."'                per call   total')"
      character(len=23), parameter :: tblfmt = "(i9,2(1x,f10.2),3x,30a)"
      integer :: nn = 0, ii = 0, level = 0, init = 0, ltr=0, lpr=0
      save

      lrt = mpipid(1) == 0 .and. omppid(1) == 0

      stdo = lgunit(1)
      if (ltrace /= 0) ltr = ltrace
      if (lprof > 0) lpr = lprof
      ltop = max0(iabs(ltr),lpr)
      jt0 = 1000*cpusec()
C     jt1 = 1000*cpusec()
      init = 1
      return

      entry tcn(str1)
      if (.not. lrt) return
      levx = level
      level = level+1
      if (level > ltop) return
      nam0 = str1
      job = 1
      goto 95

      entry tclev(str1,ilev)
      ilev = level-1
      str1 = ' '
      if (ilev >= 0) str1 = name(icall(ilev))
      return

      entry tcx(str2)
      if (.not. lrt) return
      levx = level-1
      level = level-1
      if (level+1 > ltop) return
      nam0 = str2
      job = 2

  95  continue
      if (init == 0) jt0 = 1000*cpusec()
      init = 1

c ... look for entry with same call path
      ii = 0
      do i = 1, nn, 1
        if (nam0 == name(i) .and. levx == lev(i)) then
          jcb = icb(i)
          lok = 1
          do jlev = levx-1, 0, -1
            if (jcb /= icall(jlev)) lok = 0
            jcb = icb(jcb)
          end do
          if (lok == 1) then
            ii = i
            go to 91
          end if
        end if
      end do

c ... add a new entry
      if (job == 2) call rxs('tcx: missing call to tcn: ',str2)
      nn = nn+1
      if (nn > nx) call rx('tcn: overflow')
      ncall(nn) = 0
C      iwk0(nn)=mnow*0.001+0.5
C      iwk1(nn)=-10.0
C      iwk2(nn)=mmax*0.001+0.5
C      iwk3(nn)=-10.0
      ttot(nn) = 0
      name(nn) = nam0
      lev(nn) = level-1
      icb(nn) = icall(max(level-2,0))
      ii = nn

c ... do real operations here
  91  continue
      if (job == 1) then
        it = 1000*cpusec()
        it1(ii) = it
        icall(level-1) = ii
         if (iabs(ltr) > levx .and. iprint() > 0) then
           write(stdo,100) 0.001*(it-jt0),str1
         endif
  100    format(' >>',f10.2,'   enter ',a:t43)
        return
      else
        it = 1000*cpusec()
        tm = (it-it1(ii))/1000d0
        ncall(ii) = ncall(ii)+1
        ttot(ii) = ttot(ii)+tm
        if (iabs(ltr) > levx .and. iprint() > 0) then
          write(stdo,101) 0.001d0*(it-jt0),str2,tm
  101     format(' >>',f10.2,'   exit  ',a,t35,f8.2)
        endif
        return
      endif

      entry tcprt(ifi)
      if (lpr == 0 .or. (.not. lrt)) return
      if (iprint() > 0) write (ifi,hdr) lpr
      it = 1000*cpusec()
      tm = (it-jt0)/1000d0
      if (iprint() > 0)
     .  write (ifi,tblfmt) 1,tm,tm,'main'
      lev(nn+1)=-1

      lvl=-1
      do  20  i = 1, nn
        lv = lev(i)
        if (lv+1 <= lpr) then

          if (lv > lvl) then
            tt(lv) = '|--'
            if (lv > 0) then
              if (tt(lv-1) == '|--') tt(lv-1) = '|  '
              if (tt(lv-1) == '`--') tt(lv-1) = '   '
            endif
          endif
          if (lv < lvl) tt(lv) = '|--'
          last = 1
          do  22  ii = i+1,nn+1
            if (lev(ii) == lv) last = 0
            if (lev(ii) < lv) goto 99
   22     continue
   99     continue
          if (last == 1) tt(lv) = '`--'

          tavg = ttot(i)/max(ncall(i),1)
          call strip(name(i),j1,j2)
          if (iprint() > 0)
     .      write (ifi,tblfmt) ncall(i),tavg,ttot(i),(tt(l),l=0,lev(i)),
     .        name(i)(j1:j2)
          lvl = lv
        endif
   20 continue

      end
      subroutine tc(string)
C- Routine for tracing program flow
      implicit none
      integer i1mach,ltc
      character*(*) string
      double precision x
      integer, save :: ltr=0,ncall=0

      ncall = ncall+1
      if (string == 'on' .or. string == 'ON') then
        ltr = 1
        call cpudel(i1mach(2),'set tc ...',x)
      elseif (string == 'off' .or. string == 'OFF') then
        ltr = 0
      elseif (string == 'tog' .or. string == 'TOG') then
        ltr = 1-ltr
      else
        if (ltr == 1) call cpudel(i1mach(2),string,x)
      endif
      return

      entry tcget(ltc)
      ltc = ltr
      end
