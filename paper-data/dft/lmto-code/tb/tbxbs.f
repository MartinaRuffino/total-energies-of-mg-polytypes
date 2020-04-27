      subroutine tbxbs(s_spec,nbas,nspec,alat,pos,ipc,dclabl)
C- Make bs.ext file for xbs - adapted from nfp
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nspec,ipc(1)
      double precision pos(3,1),alat,dclabl(*)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      character*8 atid,spid,spidi,spidj

      integer ifi,ib,is,js,j,m,fopn
      double precision pscale,rscale,rbot,ravg,dtop,dd,atcol(3),amp1,
     .                 colxbs(3),rad,radxbs,bondgray,bmin,bmax,brad,
     .                 z,radi,radj,rmt,rmti,rmtj,dist,x
      logical gotcol,cmdopt,a2bin
      character*72 outs

      ifi = fopn('BS')

      pscale = alat
      rscale = 1d0

      print *,'TBXBS: build xbs file..'

      write(ifi,900)
C  911 format('* ',a)
  900 format('  ')

c write coordinates
      rbot = 1d10
      ravg = 0d0
      dtop = 0d0
      do ib = 1, nbas
        is = ipc(ib)
        call r8tos8(dclabl(is),atid)
        write(ifi,821) atid, (pos(m, ib)*pscale, m = 1, 3)
  821   format(' atom  ',a,' ',3f10.4)
        rmt = s_spec(is)%rmt
        radxbs = s_spec(is)%radxbs
        if (radxbs < 0) then
          radxbs = rmt
        endif
        rbot = dmin1(rbot,radxbs)
C        ravg = ravg+rmt/nbas
        ravg = ravg + radxbs/nbas
        dd = pos(1,ib)**2+pos(2,ib)**2+pos(3,ib)**2
        dtop = dmax1(dtop,dsqrt(dd))
      enddo
      amp1 = 0d0

c write species data
      write(ifi,900)
      do is = 1, nspec
        z = s_spec(is)%z
        rmt = s_spec(is)%rmt
        colxbs = s_spec(is)%colxbs
        radxbs = s_spec(is)%radxbs
        if (radxbs < 0) then
          radxbs = rmt
        endif
        if (colxbs(1) < 0) then
          gotcol = .false.
          call dpzero(atcol,3)
c     some arbitrary function to set gray color.. heavier means darker
          x = (z - 4)/4
          x = dmax1(x,0d0)
          atcol(1) = x/(1d0+x)
          atcol(1) = 1 - atcol(1)
        else
          gotcol = .true.
          call dcopy(3,colxbs,1,atcol,1)
        endif
        rad = dmax1(radxbs*rscale,radxbs)
        call r8tos8(dclabl(is),spid)
        if (gotcol) then
          write(ifi,822)  spid,rad,atcol
        else
          write(ifi,822)  spid,rad,atcol(1)
        endif
  822   format(' spec  ',a,'  ',f9.4,3f8.2)
      enddo

c write bond data
      bmax = -1
      if (cmdopt('-bmax=',6,0,outs)) then
        j = 6
        if (.not. a2bin(outs,bmax,4,0,' ',j,len(outs)))
     .    call rx('tbxbs: bad value -bmax')
      endif
      bondgray = 1d0
      write(ifi,900)
      do is = 1, nspec
        rmti = s_spec(is)%rmt
        radi = s_spec(is)%radxbs
        if (radi < 0) then
          radi = rmti
        endif
        call r8tos8(dclabl(is),spidi)
         do js = is, nspec
           call r8tos8(dclabl(js),spidj)
           rmtj = s_spec(js)%rmt
           radj = s_spec(js)%radxbs
           if (radj < 0) then
             radj = rmtj
           endif
           bmin = 0d0
           if (bmax < 0) then
             bmax = (radi+radj)*rscale*1.5d0
           endif
           brad = ravg*rscale/7d0
           write(ifi,823) spidi,spidj,
     .           bmin,bmax,brad,bondgray
  823       format(' bonds ',a6,a6,3f9.4,f8.2)
         enddo
      enddo

c write extra data
      dist = dtop*6*pscale
      write(ifi,901) amp1,dist
  901 format(/' amp ',f9.4/' scale  20'/' dist ',f9.4/' inc 5')

      call fclose(ifi)
      end

