      subroutine mkbs(atid,spid,ips,alat,nbas,nspec,z,rmt,pos,
     .                colxbs,radxbs)
C- Make bs.ext file for xbs - adapted from nfp
      implicit none
      integer nbas,nspec,ips(1)
      double precision z(1),rmt(1),pos(3,1),colxbs(3,1),radxbs(1),alat
C      character id*120
      character*8 atid(1),spid(1)

      integer ifi,ib,is,js,m
      double precision pscale,rscale,rbot,ravg,dtop,dd,atcol(3),amp1,
     .                 col(3),rad,bondgray,bmin,bmax,brad,dist,x
      logical gotcol

      ifi=83
      call mopen(ifi,'bs','f')
c scale factors for positions and sphere radii:
      pscale=alat
      rscale=1d0

C      write(ifi,911) id(1:60)
      write(ifi,900)
  911 format('* ',a)
  900 format('  ')

c write coordinates
      rbot=1d10
      ravg=0d0
      dtop=0d0
      do ib=1,nbas
        write(ifi,821) atid(ib),(pos(m,ib)*pscale,m=1,3)
  821   format(' atom  ',a,' ',3f10.4,2x,3f9.4)
        is=ips(ib)
        rbot=dmin1(rbot,rmt(is))
        ravg=ravg+rmt(is)/nbas
        dd=pos(1,ib)**2+pos(2,ib)**2+pos(3,ib)**2
        dtop=dmax1(dtop,dsqrt(dd))
      enddo
      amp1=0d0

c write species data
      write(ifi,900)
      do is=1,nspec
        if (colxbs(1,is) < 0) then
          gotcol = .false.
          call dpzero(atcol,3)
c     some arbitrary function to set gray color.. heavier means darker
          x=(z(is)-4)/4
          x=dmax1(x,0d0)
          atcol(1)=x/(1d0+x)
          atcol(1)=1-atcol(1)
        else
          gotcol = .true.
          call dcopy(3,colxbs(1,is),1,atcol,1)
        endif
        rad = dmax1(rmt(is)*rscale,radxbs(is))
        if (gotcol) then
          write(ifi,822)  spid(is),rad,atcol
        else
          write(ifi,822)  spid(is),rad,atcol(1)
        endif
  822   format(' spec  ',a,'  ',f9.4,3f8.2)
      enddo

c write bond data
      bondgray=1.0
      write(ifi,900)
      do is=1,nspec
         do js=is,nspec
            bmin=0d0
            bmax=(rmt(is)+rmt(js))*rscale*1.5
            brad=ravg*rscale/7.0
            write(ifi,823) spid(is),spid(js),
     .           bmin,bmax,brad,bondgray
  823       format(' bonds ',a6,a6,3f9.4,f8.2)
         enddo
      enddo

c write extra data
      dist = dtop*6*pscale
      write(ifi,901) amp1,dist
  901 format(/' amp ',f9.4/' scale  20'/' dist ',f9.4/' inc 5')

      call mclose(ifi)
      end

