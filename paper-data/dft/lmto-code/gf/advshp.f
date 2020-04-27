      subroutine advshp(npl,mode,pgplp,glist,vshftp,vshft)
C- Add PL constant potential to selected ib
C-----------------------------------------------------------------------
Ci  npl  :number of principal layers (pgfset.f)
Ci mode: 1s digit
Ci       : 0  add vshftp to vshft; 1 subtract vshftp from vshft
Ci       : 10s digit
Ci       : 0 vshftp is a constant=vshftp(1);  1  vshftp depends on PL
Ci       : 100s digit
Ci       : 1 exclude PL -1, npl from shifts
Ci  gplp :index and dimensioning information for crystal subblocks.
Ci       :The meaning of pgplp depends on the context; see subas.f
Ci pgplp :index and dimensioning information for each PL (pgfset.f)
Ci glist :
Ci vshftp:
Co Outputs
Co vshft :array of site potential shifts
C-----------------------------------------------------------------------
      implicit none
C Passed variables
      integer npl,mode,glist(-1:*),pgplp(6,-1:npl)
      double precision vshftp(npl),vshft(*)
C Local variables
      integer ipl,kpl,kplv,ib,ib1,ib2,scrwid,iprint,lgunit
      double precision scale
      parameter (scrwid=80)
      integer mpipid,procid,master

      procid = mpipid(1)
      master = 0

C ... Printout
      if (iprint() > 60 .and. procid == master) then
        ipl = glist(0)
        kpl = glist(glist(-1)-1)
        if (mod(mode/100,10) == 1) then
         if (ipl < 0) ipl = glist(1)
         if (kpl == npl) kpl = glist(glist(-1)-2)
       endif
       call awrit5('%N ADVSHP: %?#n#subtract#add#'//
     .   '%?#n# PL-dependent V%j# V=%;7d# in PL %i..%i',' ',
     .   scrwid,lgunit(1),mod(mode,10),mod(mode/10,10),vshftp,ipl,kpl)
       if (mod(mode/10,10) == 1) then
         call awrit2(' V =%n:1;6,6d',' ',scrwid,lgunit(1),glist,vshftp)
       endif
      endif

      scale = 1
      if (mod(mode,10) == 1) scale = -1
      do  10  kpl = 1, glist(-1)
C   ... ipl is PL index, kplv is 1 or kpl, depending on mode
        ipl = glist(kpl-1)
        call gtibpl(ipl,npl,pgplp,ib1,ib2)
        if (mod(mode/100,10) == 1 .and.
     .    (ipl < 0 .or. ipl == npl)) goto 10
        kplv = 1
        if (mod(mode/10,10) == 1) kplv = kpl
        do  12  ib = ib1, ib2
          if (iprint() >= 70 .and. procid == master) then
          write(lgunit(1),333) ib,vshft(ib),vshft(ib)+scale*vshftp(kplv)
  333     format(i4,2f12.6)
          endif
          vshft(ib) = vshft(ib) + scale*vshftp(kplv)
   12   continue
   10 continue

      end

