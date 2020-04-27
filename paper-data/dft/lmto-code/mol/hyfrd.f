      subroutine hyfrd(tabs,otspec,otdata,ntbl,ndata)
C- Disk read and memory allocation for TCF tables
      implicit none
      character*80 tabs
      integer otspec,otdata,ntbl,ndata
      integer ier,ifi,ifile,nfile,irep,nalf,ncof,iprint,nw,j1,j2
C     character fn*4
      real w(1)
      common /w/ w

C For MPI ...
      integer mpipid,procid,master,numprocs
      logical MPI,mlog,cmdopt
      character*80 outs
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      call pshpr(iprint()-10)

C --- Determine number and size of tables ---
      call defrr(otspec, 100)
      call defrr(otdata, 100000)
      ntbl = 0
      ndata = 0
      call words(tabs,nw)
      nfile = nw
      do 10 ifile = 1, nfile
        call word(tabs,ifile,j1,j2)
        ifi = 60+ifile
        if (procid == master) then
          call mopen(ifi,tabs(j1:j2),'u')
          rewind ifi
        endif
        do 11 irep=1,99
          if (procid == master) then
            call hyfin0(ncof,nalf,w(otspec),ier,ifi)
          endif
          call mpibc1(ncof,1,2,mlog,'hyfrd','ncof')
          call mpibc1(nalf,1,2,mlog,'hyfrd','nalf')
          call mpibc1(ier,1,2,mlog,'hyfrd','ier')
          call mpibc1(w(otspec),100,4,mlog,'hyfrd','otspec')
          if (ier /= 0) goto 92
          ntbl = ntbl + 1
          if (procid == master) then
            call hyfin2(ncof,nalf,w(otdata),ier,ifi)
          endif
          call mpibc1(w(otdata),ncof*nalf,4,mlog,'hyfrd','otdata')
          if (iprint() > 40) then
            call hyfsho(w(otspec))
          endif
          if (ier /= 0) call rx('hyfrd: error reading tcdata')
          ndata = ndata + ncof*nalf
   11   continue
   92   continue
   10 continue

C --- Allocate required arrays, read tables with hyfinp ---
      call rlse(otspec)
      call defrr(otspec, 100*ntbl)
      call defrr(otdata, ndata)
      if (procid == master) then
        call hyfinp(nfile,w(otspec),w(otdata),ntbl,ndata)
      endif
      call mpibc1(w(otspec),100*ntbl,4,mlog,'hyfrd','otspec')
      call mpibc1(w(otdata),ndata,4,mlog,'hyfrd','otdata')

C --- Close files ---
      if (procid == master) then
        do  20  ifile = 1, nfile
   20   call mclose(60+ifile)
      endif
      call poppr

      end
      subroutine hyfsho(t)
C- Displays characteristics of a TCF
      implicit none
      double precision t(100)
      character*7 s,dig(-1:9)*1
      integer nx1,nx2,j
      data dig /'-','0','1','2','3','4','5','6','7','8','9'/

      nx1 = nint(t(19))
      nx2 = nint(t(20))
      print 333, t(3), nint(t(7)), nint(t(12)), nint(t(9))
      write(s,'(7a1)') '<',(dig(nint(t(40+j))),j=1,nx1),'>'
      print 334, '1', t(16), t(1), s, (t(50+j),j=1,nx1)
      write(s,'(7a1)') '<',(dig(nint(t(60+j))),j=1,nx2),'>'
      print 334, '2', t(18), t(2), s, (t(70+j),j=1,nx2)
     .
  333 format('hyfsho:  d0=',f7.3,'  nalf=',i3,'  nxi=',i5,
     .  '  ncof=',i5)
  334 format('site ',a,':  ep=',f7.3,'  r=',f10.6,'  lxi=',a7,
     .  '  exi=',6f7.2)

      end

