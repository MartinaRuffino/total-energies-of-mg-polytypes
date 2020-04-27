      subroutine pgvpl(npl,pgplp,vshft,vpl)
C- Make plane-averaged potential shifts
      implicit none
      integer npl,pgplp(6,-1:npl)
      double precision vshft(*), vpl(-1:npl)
      integer ipl,ib,ib1,ib2
      double precision vbar,rmsdv,ddot,dsum

      do  10  ipl = -1, npl
        call gtibpl(ipl,npl,pgplp,ib1,ib2)
        vpl(ipl) = 0
        do  20  ib = ib1, ib2
   20   vpl(ipl) = vpl(ipl) + vshft(ib)
        vpl(ipl) = vpl(ipl)/(ib2-ib1+1)
   10 continue
      vbar = dsum(npl,vpl(0),1) / npl
      rmsdv = sqrt(dmax1(ddot(npl,vpl(0),1,vpl(0),1)-vbar**2*npl,0d0))

C ... Printout
      if (rmsdv < 1d-6) return
C     if (iprint() < 40 .and. ifi >= 0) return
      print '(4x,''PL'',6x,''vpl'')'
      do  30  ipl = -1, npl
        print 334, ipl,vpl(ipl)
  334   format(i6,2f12.6)
   30 continue

      end
