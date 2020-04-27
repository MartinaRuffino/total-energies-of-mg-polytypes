      subroutine atfold(isw,adnf,nl,nsp,vmtz,pp,s_spec)
C- Automatic down-folding, 2nd generation LMTO
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: idxdn rmt lmxb
Co     Stored:    idxdn
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   isw:  :1 intermediate waves mapped to lower
Ci         :2 higher waves mapped to 'neglected'
Ci         :3 combination of 1 and 2
Ci         :5 all waves mapped to lower
Ci   adnf  :T automatic downfolding turned ON; see Remarks.
Ci         :F automatic downfolding turned OFF.  Initially assign
Ci          any zero values of idxdn to lower set (see Remarks)
Ci   nl    :(global maximum l) + 1
Ci   nsp   :number of spins, needed for dimensioning ppar
Ci   vmtz  :muffin-tin zero (asamad.f)
Ci   pp    :potential parameters (atomsr.f)
Ci          NB: although pp's are class-based, atfold assumes the
Ci          first nspec classes correspond to the nspec species.
Cr Remarks
Cr   There is an index in idxdn for each species and each l-channel.
Cr   It indicates in which set the orbital belongs:
Cr       1 - lower
Cr       2 - intermediate
Cr       3 - higher
Cr   If any of these are 0 on input, it will be automatically
Cr   assigned to a set using the following rules, provided adnf=T.
Cr   If |V^0 - vmtz|s^2 > 1  AND  |tau - 1| < tolh1
Cr      AND  (C - enu)sqrt p^gamma  >  tolh2  ==>  higher
Cr   If |gamma - gamma^safe|  >  |P^-1 - gamma^safe|
Cr      AND  (C - enu)sqrt p^gamma  >  toli   ==>  intermediate
Cr   Otherwise ==>  lower
Cr   You can make up your own downfolding rules.
Cu Updates
Cu   06 Sep 11 Started migration to f90 structures
Cu   10 Apr 02 Redimensionsed idxdn to accomodate larger lmxb
Cu   14 Sep 99 atfold rewrote to species-based idxdn, with new arguments
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical adnf
      integer isw,nl,nsp
      double precision pp(6,0:nl-1,nsp,1),avw,vmtz
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer, parameter :: n0=10, nkap0=4
      integer is,ic,l,isp,lgunit,iprint,idxdn(0:n0*nkap0),nspec,
     .  lmxb,stdo
      double precision tolh1,tolh2,toli,Gamma,work(5),rmax,dglob
      parameter ( tolh1 = 0.05d0, tolh2 = 0.3d0, toli = 0.25d0 )
      logical hcond(0:3),icond(0:2)
      procedure(integer) :: nglob

      avw = dglob('avw',0d0,0)
      nspec = nglob('nspec')
      stdo = nglob('stdo')
C     choose arbitrarily :  spins, classes not distinguished
      isp = 1
      do  is = 1, nspec
      idxdn(0:n0-1) = s_spec(is)%idxdn(1:n0,1)
      rmax = s_spec(is)%rmt
      lmxb = s_spec(is)%lmxb
C     Assume representative ic = is
      ic = is
C ... Any idxdn for l>lmxb is automatically '4'
        do  l = lmxb+1, nl-1
          idxdn(l) = 4
        enddo
C ... If no automatic downfolding, initial zero's to 1
      if (.not. adnf) then
          do  l = 0, lmxb
            if (idxdn(l) == 0) idxdn(l) = 1
          enddo
          if (isw == 0) cycle
      endif

      if (iprint() >= 40 .and. is == 1) write(stdo,1)
    1 format(' Atfold: checks for automatic downfolding'/
     .  ' spc l    |g - g_s|  |P^-1 - g_s| (C-enu)s(p) |V0-vmtz|s^2',
     .  ' |tau-1|  Conditions')

        do  l = 0, lmxb

        if (pp(5,l,isp,ic) == 0) then
          if (idxdn(l) == 0) idxdn(l) = 1
          cycle
        endif

C   --- V^0 ---
        work(1) = pp(2,l,isp,ic) - pp(3,l,isp,ic)**2/pp(5,l,isp,ic)

C   --- tau ---
        Gamma = pp(3,l,isp,ic)**2/pp(5,l,isp,ic)**2
        work(2) = 2*(2*l+1)**2*(2*l+3)*(avw/rmax)**(2*l+1)
        work(2) = work(2)/Gamma/rmax**2

C   --- ( C - enu ) sqrt p^gamma ---
        work(3) = (pp(2,l,isp,ic) - pp(1,l,isp,ic))
     .                 *dsqrt(pp(4,l,isp,ic))

C   --- Gamma^safe ---
        work(4) = (rmax/avw)**(2*l+1)/(4*(l+1)*(2*l+1))

C   --- Inverse unscreened potential function at enu (2nd order) ---
        work(5) =  pp(3,l,isp,ic)**2 /
     .            (pp(1,l,isp,ic) - pp(2,l,isp,ic)) + pp(5,l,isp,ic)

        icond(1) = dabs(pp(5,l,isp,ic) - work(4)) >
     .             dabs(work(5) - work(4))
        icond(2) =  work(3) > toli
        hcond(1) =  abs( work(1) - vmtz )*rmax**2 > 1d0
        hcond(2) =  work(3) > tolh2
        hcond(3) =  dabs( work(2) -  1d0 ) < tolh1
        if (iprint() >= 40) print 334, ic, l,
     .    dabs(pp(5,l,isp,ic)-work(4)),dabs(work(5)-work(4)),work(3),
     .    dabs(work(1)-vmtz)*rmax**2,dabs(work(2)-1d0),
     .    icond(1), icond(2), idxdn(l) == 3 .or. idxdn(l) == 0,
     .    hcond(1), hcond(2), hcond(3),
     .    idxdn(l) == 3 .or. idxdn(l) == 0
  334   format(2i3,5f12.5,2x,3l1,2x,4l1)

C   ... shift to intermediate if icond(0) is true
        icond(0) = icond(1) .and. icond(2) .and. idxdn(l) == 0
        icond(0) = icond(0) .or.
     .    (hcond(1) .and. hcond(2) .and. hcond(3) .and. isw >= 2)
        icond(0) = icond(0) .or. idxdn(l) > 2 .and. isw >= 2
        icond(0) = icond(0) .and. mod(isw,2) == 0

C   ... shift to high if hcond(0) is true
        hcond(0) = hcond(1) .and. hcond(2) .and. hcond(3) .and.
     .             idxdn(l) == 0 .and. isw < 2

        if (hcond(0)) then
C --- higher set ---
          idxdn(l) = 3
            if (iprint() >= 10) write (lgunit(1),3) ic,l
            if (iprint() >= 10) write (lgunit(2),3) ic,l
            cycle
        endif
        if (icond(0)) then
C --- intermediate set ---
          idxdn(l) = 2
            if (iprint() >= 10) write (lgunit(1),4) ic,l
            if (iprint() >= 10) write (lgunit(2),4) ic,l
            cycle
        endif
C --- lower set ---
          if (idxdn(l) == 2 .and. mod(isw,2) == 1 .or. idxdn(l)
     . > 2 .and. isw == 5) then
            if (iprint() >= 10) write (lgunit(1),5) ic,l
            if (iprint() >= 10) write (lgunit(2),5) ic,l
        endif
        if ( idxdn(l) == 0 ) idxdn(l) = 1
        enddo
      s_spec(is)%idxdn = reshape(idxdn(0:n0*nkap0-1), (/n0,nkap0/))
      enddo
 3    format(1x,'Atfold:  species',i4,',  l = ',i1,' to higher set')
 4    format(1x,'Atfold:  species',i4,',  l = ',i1,
     .                                           ' to intermediate set')
 5    format(1x,'Atfold:  species',i4,',  l = ',i1,' to lower set')
      end
