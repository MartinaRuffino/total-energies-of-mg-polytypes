      subroutine optshs(ikp,nkp,nchi2,iabc,optme,g,ngrp,
     .  esciss,nfilo,nfiup,nemlo,nemup,nbmax,eband)
C- Nonlinear optical matrix elements
C  under construction, adapted from Sergey's shg_sme.f
Cu  18 Jun 18 Optical matrix elements not generated, but passed in optdme
      implicit none
C ... Passed parameters
      integer ikp,nkp,nfilo,nfiup,nemlo,nemup
      integer nbmax,nchi2,iabc(3,6),ngrp
      double precision esciss,eband(nbmax)
      double precision g(3,3,ngrp)
      complex(8) :: optme(nfilo:nemup,nfilo:nemup,3)
C ... Local parameters
C  ...  Formerly contained in 'PAR_OPT'
        integer kkmax,nfilmax,nempmax
        integer kk2,nfil2,nemp2,npmax
        integer kk3,nfil3,nemp3
        integer kk4,nfil4,nemp4

C  Optics (linear)
        parameter(kkmax   =  250)    !!! max.num. of k-points  for optics
        parameter(nfilmax =   25)    !!! max.num. filled bands for optics
        parameter(nempmax =   15)    !!! max.num. empty  bands for optics
C  SHG- static
        parameter(kk2      = 800)    !!! max.num. of k-points  for SHG
        parameter(nfil2    =  31)    !!! max.num. filled bands for SHG
        parameter(nemp2    =  30)    !!! max.num. empty  bands for SHG
C  Frequency- dependent SHG and LEO:
        parameter(npmax    =   6)    !!! max.num. of polar. (<= 6)
C  SHG- frequency dependent
        parameter(kk3      = 100)    !!! max.num. of k-points  for SHG
        parameter(nfil3    =  10)    !!! max.num. filled bands for SHG
        parameter(nemp3    =  20)    !!! max.num. empty  bands for SHG
C  LEO- frequency dependent
        parameter(kk4      = 100)    !!! max.num. of k-points  for LEO 600
        parameter(nfil4    =  10)    !!! max.num. filled bands for LEO 25
        parameter(nemp4    =  10)    !!! max.num. empty  bands for LEO 30
C   ... end of PAR_OPT

      double precision shg_me(kk2,6),shg_me_ee(kk2,6),
     .  shg_me_ei2(kk2,6),shg_me_ei3(kk2,6,nfil2,nemp2)
      common/shg_s/ shg_me,shg_me_ee,shg_me_ei2,shg_me_ei3
      integer il,ir
CC !!! FIX Fermi distribution function
      double precision f(nbmax)
C !!! FIX
      double precision ppp(3,nfil2+nemp2,nfil2+nemp2,2)
      integer j1,j2,j3
      integer i1,i2,i3
      integer i,j,k,kk,igrp

      double precision optr_a(3),opti_a(3)
      double precision optr_b(3),opti_b(3)
      double precision optr_c(3),opti_c(3)
      double precision optr_ag(3),opti_ag(3)
      double precision optr_bg(3),opti_bg(3)
      double precision optr_cg(3),opti_cg(3)
      complex*16 opt_ag(3),opt_bg(3),opt_cg(3)
      double precision sum_gr_ee(6),sum_gr_ei2(6),sum_gr_ei3(6)
      double precision xx1,xx2,xx3,xx4,xx5,xx6
      double precision e12,e23,e31
      double precision e12_0,e23_0,e31_0
      double precision f12,f23,f31
      integer ia,ib,ic
      double precision anor_a,anor_b,anor_c
      double precision gam, gam1,rat1
      double precision sum_ee(6),sum_ei2(6),
     .  sum_ei3(6,nfil2+nemp2,nfil2+nemp2)

      double precision delta

      double precision optrr(3),optii(3)

      if (nkp > kk2) call rx('Increase kk2 at optsms.f')
      if ((nfiup-nfilo+1) > nfil2) call rx('Increase nfil2 at optshs.f')
      if ((nemup-nemlo+1) > nemp2) call rx('Increase nemp2 at optshs.f')
      if (nemlo <= nfiup) call rx('Check nemlo,nfiup; SHG program is built for sem. only')

C      do  k = 1, 3
C        call yprmi('optme k = %i',k,0,3,optme(nfilo,nfilo,k),0,nemup-nfilo+1,nemup-nfilo+1,nemup-nfilo+1)
C      enddo

C ... Fermi distribution: assume insulator for now
      do  10  j = 1, nemup
      f(j) = 0
   10 if (j <= nfiup) f(j) = 1

C --- Copy matrix elements between all pairs of levels ---
      il = 0
      do  i = nfilo, nemup
      il = il+1
      ir = 0
      do  j = nfilo, nemup
        ir = ir+1
        if (ir < il) cycle

        do  k = 1, 3
          optrr(k) = optme(i,j,k)
          optii(k) = dimag(optme(i,j,k))
        enddo

C   ... Save for chi2.  grad is anti-hermitian operator
        ppp(1,il,ir,1) = optrr(1)
        ppp(2,il,ir,1) = optrr(2)
        ppp(3,il,ir,1) = optrr(3)
        ppp(1,il,ir,2) = optii(1)
        ppp(2,il,ir,2) = optii(2)
        ppp(3,il,ir,2) = optii(3)
        if (ir /= il) then
          ppp(1,ir,il,1) = -optrr(1)
          ppp(2,ir,il,1) = -optrr(2)
          ppp(3,ir,il,1) = -optrr(3)
          ppp(1,ir,il,2) = optii(1)
          ppp(2,ir,il,2) = optii(2)
          ppp(3,ir,il,2) = optii(3)
        endif
      enddo
      enddo

C   ... Setup for SHG.  Cutoff parameters:
        gam = 0.001d0
        gam1 = 0.01d0
        rat1 = 10d0
C   ... gap "scissors"
        delta = esciss

C   ... Setup to make integrand for SHG
        do  30  k = 1, nchi2
          shg_me(ikp,k) = 0d0
          shg_me_ee(ikp,k) = 0d0
          shg_me_ei2(ikp,k) = 0d0
          sum_ee(k) = 0d0
          sum_ei2(k) = 0d0
          do  31  i1 = 1, nfiup-nfilo+1
          do  31  i2 = 1, nemup-nemlo+1
   31     shg_me_ei3(ikp,k,i1,i2) = 0d0
   30     continue
        do  32  k = 1, nchi2
        do  32  i1 = 1, nfil2+nemp2
        do  32  i2 = 1, nfil2+nemp2
   32   sum_ei3(k,i1,i2) = 0d0

C   --- Integrand for SHG: pure interband + three-band intraband part --
        i1 = 0
        do  40  j1 = nfilo, nemup
        i1 = i1+1

        i2 = 0
        do  40  j2 = nfilo, nemup
        i2 = i2+1

        i3 = 0
        do  40  j3 = nfilo, nemup
        i3 = i3+1

C   ... Each of the j1,j2,j3 must be different
        if (((j1-j2)*(j2-j3)*(j3-j1)) == 0) goto 40
C   ... No processes of the type (ccc) or (vvv)
        if (j1 <= nfiup.and.j2 <= nfiup.and.j3 <= nfiup) goto 40
        if (j1 >= nemlo.and.j2 >= nemlo.and.j3 >= nemlo) goto 40

C   ... Energy Factors
        e12_0 = eband(j1) - eband(j2)
        e23_0 = eband(j2) - eband(j3)
        e31_0 = eband(j3) - eband(j1)
        e12 = e12_0 + delta*(f(j2)-f(j1))
        e23 = e23_0 + delta*(f(j3)-f(j2))
        e31 = e31_0 + delta*(f(j1)-f(j3))
        if (abs(e12) < gam .or. abs(e23) < gam .or.
     .      abs(e31) < gam) goto 40

C   ... Fermi Factors
        f12 = f(j1)-f(j2)
        f23 = f(j2)-f(j3)
        f31 = f(j3)-f(j1)
        do  51  i = 1, 3
          optr_a(i) = ppp(i,i1,i2,1)
          opti_a(i) = ppp(i,i1,i2,2)
          optr_b(i) = ppp(i,i2,i3,1)
          opti_b(i) = ppp(i,i2,i3,2)
          optr_c(i) = ppp(i,i3,i1,1)
          opti_c(i) = ppp(i,i3,i1,2)
   51   continue
        anor_a = 0d0
        anor_b = 0d0
        anor_c = 0d0
        do  52  i = 1, 3
          anor_a = anor_a + (optr_a(i)**2+opti_a(i)**2)
          anor_b = anor_b + (optr_b(i)**2+opti_b(i)**2)
          anor_c = anor_c + (optr_c(i)**2+opti_c(i)**2)
   52   continue
        if (abs(e12) <= gam1.and.sqrt(anor_a)/abs(e12) > rat1) then
          do  53  i = 1, 3
            optr_a(i) = 0d0
            opti_a(i) = 0d0
   53     continue
        endif
        if (abs(e23) <= gam1.and.sqrt(anor_b)/abs(e23) > rat1) then
          do  54  i = 1, 3
            optr_b(i) = 0d0
            opti_b(i) = 0d0
   54     continue
        endif
        if (abs(e31) <= gam1.and.sqrt(anor_c)/abs(e31) > rat1) then
          do  55  i = 1, 3
            optr_c(i) = 0d0
            opti_c(i) = 0d0
   55     continue
        endif

C   ... ??? Make symmetrized chi^2_{abc}, static limit ???
        call dpzero(sum_gr_ee, 6)
        call dpzero(sum_gr_ei3,6)
        do  60  igrp = 1, ngrp
          call grpop(optr_a,optr_ag,g,igrp)
          call grpop(opti_a,opti_ag,g,igrp)
          do  61  kk = 1, 3
   61     opt_ag(kk) = dcmplx(optr_ag(kk),opti_ag(kk))
          call grpop(optr_b,optr_bg,g,igrp)
          call grpop(opti_b,opti_bg,g,igrp)
          do  62  kk = 1, 3
   62     opt_bg(kk) = dcmplx(optr_bg(kk),opti_bg(kk))
          call grpop(optr_c,optr_cg,g,igrp)
          call grpop(opti_c,opti_cg,g,igrp)
          do  63  kk = 1, 3
  63     opt_cg(kk) = dcmplx(optr_cg(kk),opti_cg(kk))

          do  65  k = 1, nchi2
            ia = iabc(1,k)
            ib = iabc(2,k)
            ic = iabc(3,k)
            xx1 = dble(opt_ag(ia)*opt_bg(ib)*opt_cg(ic))
            xx2 = dble(opt_ag(ia)*opt_bg(ic)*opt_cg(ib))
            xx3 = dble(opt_ag(ib)*opt_bg(ia)*opt_cg(ic))
            xx4 = dble(opt_ag(ib)*opt_bg(ic)*opt_cg(ia))
            xx5 = dble(opt_ag(ic)*opt_bg(ia)*opt_cg(ib))
            xx6 = dble(opt_ag(ic)*opt_bg(ib)*opt_cg(ia))
C       ... Interband
            sum_gr_ee(k) = sum_gr_ee(k) + xx1 + xx2
C       ... Three-Band Intraband
            sum_gr_ei3(k) = sum_gr_ei3(k) +
     .         xx1 + xx2 + xx3 + xx4 + xx5 + xx6
   65     continue
   60   continue
        do  999 k = 1, nchi2
          sum_gr_ee(k) = sum_gr_ee(k)/ngrp
          sum_gr_ei3(k) = sum_gr_ei3(k)/ngrp
  999   continue

        do  70  k = 1, nchi2
C     ... Interband
          sum_ee(k) = sum_ee(k) + sum_gr_ee(k)/(e12_0*e23_0*e31_0)
     .      * (f23/e23 + f31/e31 - 2d0*f12/e12)/(e31-e23)
C         if (k==1) print *, k,sum_ee(k)
C     ... Three-Band Intraband
          if (f12 > 0d0) then
            sum_ei3(k,i1,i2) = sum_ei3(k,i1,i2)
     .        + sum_gr_ei3(k)/(e12_0*e23_0*e31_0)/2d0
     .        *(e31-e23)*f12    !!!/e12**3  *2d0 !!!*2. iz-za vybora 1-2
          endif
   70   continue
   40 continue


C --- Integrand for SHG: 2-band intraband part --
      i1 = 0
      do  140  j1 = nfilo, nemup
      i1 = i1+1

      i2 = 0
      do  140  j2 = nfilo, nemup
      i2 = i2+1

        if ((j1-j2) == 0) goto 140
C   ... Energy Factors
        e12_0 = eband(j1) - eband(j2)
        e12 = e12_0 + delta*(f(j2)-f(j1))
        if (abs(e12) < gam) goto 140

C   ... Fermi Factors
        f12 = f(j1)-f(j2)
        do  151  i = 1, 3
          optr_a(i) = ppp(i,i1,i1,1)-ppp(i,i2,i2,1)
          opti_a(i) = ppp(i,i1,i1,2)-ppp(i,i2,i2,2)
          optr_b(i) = ppp(i,i1,i2,1)
          opti_b(i) = ppp(i,i1,i2,2)
          optr_c(i) = ppp(i,i2,i1,1)
          opti_c(i) = ppp(i,i2,i1,2)
  151   continue

C   ... Symmetrize
        do  156  k = 1, nchi2
  156   sum_gr_ei2(k) = 0d0
        do  160  igrp = 1, ngrp
          call grpop(optr_a,optr_ag,g,igrp)
          call grpop(opti_a,opti_ag,g,igrp)
          do  161  kk = 1, 3
  161     opt_ag(kk) = dcmplx(optr_ag(kk),opti_ag(kk))
          call grpop(optr_b,optr_bg,g,igrp)
          call grpop(opti_b,opti_bg,g,igrp)
          do  162  kk = 1, 3
  162     opt_bg(kk) = dcmplx(optr_bg(kk),opti_bg(kk))
          call grpop(optr_c,optr_cg,g,igrp)
          call grpop(opti_c,opti_cg,g,igrp)
          do  163  kk = 1, 3
  163     opt_cg(kk) = dcmplx(optr_cg(kk),opti_cg(kk))

          do  165  k = 1, nchi2
            ia = iabc(1,k)
            ib = iabc(2,k)
            ic = iabc(3,k)
             sum_gr_ei2(k) = sum_gr_ei2(k) +
     .             ( dble(opt_ag(ia)*opt_bg(ib)*opt_cg(ic))
     .             + dble(opt_ag(ia)*opt_bg(ic)*opt_cg(ib))
     .             + dble(opt_ag(ib)*opt_bg(ia)*opt_cg(ic))
     .             + dble(opt_ag(ib)*opt_bg(ic)*opt_cg(ia))
     .             + dble(opt_ag(ic)*opt_bg(ia)*opt_cg(ib))
     .             + dble(opt_ag(ic)*opt_bg(ib)*opt_cg(ia)) )

  165     continue
  160   continue
        do  170  k = 1, nchi2
  170   sum_ei2(k) = sum_ei2(k) + sum_gr_ei2(k)*(-1d0)*f12
     .          /e12_0**2/e12**3/ngrp
  140 continue
C ... End of assembling 2-band intraband SHG integrand

C --- Total integrand for SHG ---
      do  200  k = 1, nchi2
        shg_me_ee(ikp,k) = sum_ee(k)
        shg_me_ei2(ikp,k) = sum_ei2(k)
  200 continue
      i1 = 0
      do  210  j1 = nfilo, nfiup
        i1 = i1+1

        i2 = 0
        do  210  j2 = nemlo, nemup
        i2 = i2+1

        do  210  k = 1, nchi2
        shg_me_ei3(ikp,k,i1,i2) = sum_ei3(k,i1,i2+nemlo-nfilo)
  210 continue

      end

