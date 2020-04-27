      subroutine gfzerq(mode,npl,nsp,pgplp,ipc,znuc,qc,glist,ef0,dos,
     .  dosne,nzne,zval,vshftp,vhold,rmsdv,rmsdq)
C- Estimate potential shift in PL to ensure neutrality
C-----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 impose neutrality for each layer separately
Ci          1 impose neutrality with a single shift for all layers
Ci            0 < ipl < npl
Ci          2 impose neutrality with a single shift for all layers
Ci   npl   :number of principal layers (pgfset.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pgplp :index and dimensioning information for crystal subblocks.
Ci          The meaning of pgplp depends on the context; see subas.f
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   znuc  :nuclear charges
Ci   qc    :sphere core charge
Ci   glist :list of PL for which to compute vshftp
Ci   ef0   :Fermi level (used for printout only)
Ci   dos   :density of states
Ci   dos(1):dos(1,isp,ipl) DOS(efermi) for spin isp and layer ipl
Ci         :dos(2,isp,ipl) integral of DOS to efermi.
Ci   dosne :dos at endpoints of nonequilibrium contour
Ci   nzne  :number of energy points along nonequilibrium contour
Ci   zval  :target electron charge.
Ci         :If zval<0, find target charge internally
Cio Inputs/Outputs
Cio  vhold :vhold(*,kpl): information passed to rfalsi, used to
Cio                        bracket minimum
Cio          mode=0: vhold generated for each layer separately
Cio          mode=1: vhold generated for layers combined -> vhold(*,1)
Cio          On output, vhold is updatd
Cio  vshftp :estimate to shift PL to neutral point
Cio         :On input, vshftp is estimate for prior iteration
Cio         :(0 for first iteration)
Cio         :On output, vshftp is updated as a new estimate
Co Outputs
Co   rmsdv :rms change in vshftp
Co   rmsdq :rms change in estimated charge
Cr Remarks
Cr   The simplest estimate for dv is:
Cr     dv = (ztarg-zhave)/dos
Cr   If two points in the integrated DOS are known, regula falsi is used.
Cu Updates
Cu   27 Jan 04 (S.Faleev) Argument list changed to handle non-equil.mode
Cu   27 Feb 02 Adapted from old pgzerq
C-----------------------------------------------------------------------
      implicit none
C Passed variables
      integer mode,npl,nsp,glist(-1:*),ipc(*),pgplp(6,-1:npl),nzne
      double precision dos(5,nsp,*),znuc(*),qc(*),vhold(12,npl),
     .  vshftp(npl),ef0,rmsdv,rmsdq,dosne(2,nsp,*)
C Local variables
      integer ipl,kpl,ipr,ib1,ib2,nipl
      double precision dv,dosi,dose,zval,ztarg,dvcap,dznow,sumq,sumde,
     .  sumz,sumdi,vold
      parameter (dvcap=.2d0)

      if (glist(-1) <= 0) return
      call getpr(ipr)
      if (glist(-1) == 1) then
        dose = dos(1,1,1) + (nsp-1)*dos(1,nsp,1)
        call info5(40,0,0,
     .    ' gfzerq: estimate '//
     .    'pot shift for charge neutrality,'//
     .    '%?#n==0# using dos(Ef0=%;5d) = %;3g# Ef0=%;5d#',
     .    nint(vhold(12,1)),ef0,dose,0,0)
      else
        call info2(40,0,0,
     .    ' gfzerq: estimate '//
     .    'pot shift for charge neutrality,'//
     .    '%?#n==0# from dos,##'//
     .    ' Ef0 = %;5d',nint(vhold(12,1)),ef0)
      endif

C        if (ipr >= 40) print 333
C  333   format('  PL   ztarg   vshiftp    vold       N(Ef0)',
C     .    '      D(Ef0)    est dV')

      rmsdv = 0
      rmsdq = 0
      sumq = 0
      sumz = 0
      sumde  = 0
      nipl = 0
      do  10  kpl = 1, glist(-1)
        ipl = glist(kpl-1)
        if (mode == 1 .and. (ipl == -1 .or. ipl == npl)) goto 10
        nipl = nipl+1
        call gtibpl(ipl,npl,pgplp,ib1,ib2)
        ztarg = zval
        if (ztarg < 0) call gtqval(ib1,ib2,znuc,qc,ipc,ztarg)
        dosi = dos(2,1,kpl) + (nsp-1)*dos(2,nsp,kpl)
        dose = dos(1,1,kpl) + (nsp-1)*dos(1,nsp,kpl)
        if (nzne >= 2) dose = dose+dosne(2,1,kpl)-dosne(1,1,kpl)
     .     + (nsp-1)*( dosne(2,nsp,kpl)-dosne(1,nsp,kpl) )
C   ... dznow is difference between target and actual charge
        dznow = dosi - ztarg
C   ... sumz is the sum total charge in all layers
        sumz  = sumz + ztarg
C   ... sumq is the difference between target and actual charge
        sumq  = sumq + dznow
        rmsdq = rmsdq + dznow**2
        if (ipl >= 0 .and. ipl < npl .or. mode == 2)
     .    sumde = sumde + dose
        if (mode == 0) then
          vold = vshftp(kpl)
          call dvdos(vshftp(kpl),dosi,dose,vhold(1,kpl),ztarg,dvcap,dv)
          call info8(20,0,0,' gfzerq: ef=%;3g  dq=%;3g  vold=%;3g'//
     .      '  dv=%;3g  vnew= %;6d  layer %i',
     .      ef0,dznow,vold,abs(dv),vshftp(kpl),ipl,0,0)
          rmsdv = rmsdv + dv**2
        endif
   10 continue

      sumdi = sumq + sumz
      if (mode == 0) then
        rmsdv = dsqrt(rmsdv/glist(-1))
        rmsdq = dsqrt(rmsdq/glist(-1))
      else
        vold = vshftp(1)
        call dvdos(vshftp,sumdi,sumde,vhold(1,1),sumz,dvcap,dv)
        rmsdv = abs(dv)
        rmsdq = abs(sumq)
      endif
      if (mode > 0) call info5(20,0,0,
     .  ' gfzerq: ef=%;3g  dq=%;3g  vold=%;3g  dv=%;3g  vnew= %;6d',
     .  ef0,rmsdq,vold,rmsdv,vshftp)
      if (mode == 0) call info2(20,0,0,
     .  '%9faverage <dq>=%;3g  <dv>=%;3g ',rmsdq,rmsdv)
      end
