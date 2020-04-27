      subroutine shoctl(s_ctrl,s_spec,s_pot,iopt,ifi)
C- Writes P, Q to stdout and log in format suitable for ctrl file
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nl nspin nclasp nccomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ics dclabl initc nrc
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pnu qnu pp ves
Cio    Passed to:  *
Ci Inputs
Ci   iopt  :1s digit    0 -> don't print out ves(rmt)
Ci         :           >0 -> do    print out ves(rmt)
Ci         :10s digit   1 spin-average moments
Ci         :10s digit   2 spin-flip moments
Ci         :100s digit  0 -> don't print out enu
Ci         :           >0 -> do    print out enu
Ci         :1000s digit 1 print out sum of sphere charges
Ci   ifi   :file logical unit
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   29 Jul 16 New option to print total ASA charge
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   30 Sep 04 Optionally writes enu
Cu   23 Jan 03 Extended options to average spins
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ifi,iopt
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      integer, allocatable :: lmx(:)
C ... Local parameters
      logical swves,swenu,swq
      integer nclspd,nl,nsp
      double precision xx,q(2)

      swves = mod(iopt,10) /= 0
      swenu = mod(iopt/100,10) /= 0
      swq   = mod(iopt/1000,10) /= 0

      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      nclspd = s_ctrl%nclasp + s_ctrl%nccomp

      allocate(lmx(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'lmxa',1,lmx,xx)

      if (mod(iopt/10,10) /= 0) then
        call spinav(mod(iopt/10,10)-1,nclspd,nl,nsp,s_pot%pnu,
     .    s_pot%qnu)
      endif

      call pvshoc(nclspd,s_ctrl%dclabl,nl,nsp,lmx,s_pot%pnu,s_pot%qnu,
     .  s_pot%pp,swves,swenu,s_pot%ves,s_ctrl%initc,s_ctrl%nrc,ifi)

      if (swq) then
        call asatotq(1,nclspd,nl,nsp,lmx,s_ctrl%nrc,s_pot%qnu,q,ifi)
      endif

      deallocate(lmx)

      end
      subroutine pvshoc(nclspd,dclabl,nl,nsp,lmx,pl,ql,pp,swves,swenu,
     .  ves,initc,nrc,ifi)
C- Printout ASA moments
C ----------------------------------------------------------------
Ci Inputs
Ci   clabl,nl,nsp,lmx,pnu,emom
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      logical swves,swenu
      integer nclspd,nl,nsp,lmx(0:nclspd),ifi,initc(0:nclspd),nrc(0:nclspd)
      double precision pl(0:nl-1,nsp,0:nclspd),pp(6,0:nl-1,nsp,0:*),
     .  ql(3,0:nl-1,nsp,0:nclspd),ves(0:*),dclabl(0:*)
C Local parameters
      integer ic,i,l,isp
      character clabl*8,strn*40

      do  ic = 0, nclspd-1
        call r8tos8(dclabl(ic),clabl)
        if (mod(initc(ic),2) /= 1) call info0(10,0,0,
     .    ' shoctl (warning) no moments for class '//clabl)
        call awrit2('        ATOM='//clabl//' P=%n;11,7D',' ',80,ifi,
     .    lmx(ic)+1,pl(0,1,ic))
        if (nsp == 2)
     .    call awrit2('%24f%n;11,7D',' ',80,ifi,lmx(ic)+1,pl(0,2,ic))
   23   format('        ATOM=',A8,' P=',20(4f11.7:,/24x))
   22   format('        ATOM=',A8,' P=',20(3f11.7:,/24x))
   21   format('        ATOM=',A8,' P=',20(2f11.7:,/24x))
        write(ifi,20) (((ql(i,l,isp,ic), i=1,3), l=0,lmx(ic)),isp=1,nsp)
   20   format(20x,'  Q=',20(3f11.7:,/24x))
        call awrit1('(20x,''ENU='',20(%if11.7:/24x))',
     .    strn,len(strn),0,lmx(ic)+1)
        if (swenu) write(ifi,strn)
     .    ((pp(1,l,isp,ic), l=0,lmx(ic)),isp=1,nsp)
        if (swves) write(ifi,43) ves(ic)
      enddo
   43 format(20x,'  V=',f11.7)

      end
      subroutine asatotq(ipr,nclspd,nl,nsp,lmx,nrc,qnu,q,ifi)
C- Compute and optionally print total ASA charge
      implicit none
      integer ipr,nclspd,nl,nsp,ifi,lmx(0:nclspd-1),nrc(0:nclspd-1)
      double precision qnu(3,0:nl-1,nsp,0:nclspd-1),q(2)
      integer ic,l

      q = 0
      do  ic = 0, nclspd-1
        do  l = 0, lmx(ic)
          q(1) = q(1) + qnu(1,l,1,ic)*nrc(ic)
          q(2) = q(2) + qnu(1,l,nsp,ic)*nrc(ic)
        enddo
      enddo

      if (ipr == 0) return

      call awrit5(' Total sphere charge %;12,7D %?#n==2#   moment %;12,7D ##',
     .  ' ',80,ifi,q(1)+q(2),nsp,q(1)-q(2),4,5)

      end
