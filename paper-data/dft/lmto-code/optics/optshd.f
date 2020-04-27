      subroutine optshd(nf1,nf2,nl,nbas,isp,nsp,ikp,nkp,ldim,ipc,
     .  nchi2,iabc,
     .  nclass,gradm,icheck,
     .  g,ngrp,qp,esciss,
     .  nfilo,nfiup,nemlo,nemup,nbmax,eband,ccp)
C-
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ikp   :k-point label
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   ipc:the jth atom belongs to class ipc(j)
Ci   nchi2 :number of (abc) for which to calculate xi^(abc) (not used)
Ci   iabc  :(abc) axes in xi^(abc) for nonlinear optics
Ci   nclass:number of inequivalent classes
Ci   aamt,abmt,bamt,bbmt: full matrix elements of (phi2 grad phi1)
Ci         :for pairs of (phi,phidot) products.
Ci   icheck:not used
Ci   g     :point group operations
Ci   ngrp  :number of point group operations
Ci   qp    :k-point
Ci   esciss:Energy shift (scissors operator)
Ci   nfilo :first occupied band
Ci   nfiup :last  occupied band
Ci   nemlo :first unoccupied band
Ci   nemup :last  unoccupied band
Ci   nbmax :maximum number of bands
Ci   eband :energy bands; alias eb (sec*.f)
Ci   ccp   :decomposition of norm (aka dnpp), in true spherical harm.
Co Outputs
Co   ?
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Dec 13 Matrix elements now contained in gradm
Cu   13 Jul 05 (wrl) updated from original NL optics optsms
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nf1,nf2,nl,nbas,isp,nsp,ikp,nkp,ldim,nchi2,nclass,
     .        icheck,ngrp,nfilo,nfiup,nemlo,nemup,nbmax
      integer iabc(3,6),ipc(nbas)
      double precision esciss
      double precision eband(nbmax,nsp,nkp),g(3,3,48),qp(3)
      complex*16 ccp(0:nl*nl-1,nbas,2,ldim,nsp)
      double precision gradm(3,nf1,nf2,nl**2,nl**2,nsp,nclass)

C ... Local parameters
c!!! Parameters: kk3,nfil3,nemp3,npmax from 'PAR_OPT' (See that file)
      include 'PAR_OPT'

      integer il,ir,jsp,if1,if2
      integer j1,j2,j3
      integer i1,i2,i3

      double precision qmt_ee (npmax,nfil3,nemp3,kk3) !w
      double precision qmt_ei3(npmax,nfil3,nemp3,kk3) !w
      double precision pmt_ee (npmax,nfil3,nemp3,kk3) !2w
      double precision pmt_ei3(npmax,nfil3,nemp3,kk3) !2w
      common/shgmt/ qmt_ee,qmt_ei3,pmt_ee,pmt_ei3

      double precision anorm
      double precision f(nbmax) !!! Fermi distribution function
      double precision ppp(3,nfil3+nemp3,nfil3+nemp3,2)


      double precision optr_a(3),opti_a(3)
      double precision optr_b(3),opti_b(3)
      double precision optr_c(3),opti_c(3)
      double precision optr_ag(3),opti_ag(3)
      double precision optr_bg(3),opti_bg(3)
      double precision optr_cg(3),opti_cg(3)
      complex*16 opt_ag(3),opt_bg(3),opt_cg(3)

      double precision sum_gr1(npmax),sum_gr2(npmax),sum_gr3(npmax)
      double precision sum_gr4(npmax)
      double precision sum_e1(npmax,nfil3+nemp3,nfil3+nemp3)
      double precision sum_e2(npmax,nfil3+nemp3,nfil3+nemp3)
      double precision sum_e3(npmax,nfil3+nemp3,nfil3+nemp3)
      double precision sum_i1(npmax,nfil3+nemp3,nfil3+nemp3)
      double precision sum_i2(npmax,nfil3+nemp3,nfil3+nemp3)
      double precision sum_i3(npmax,nfil3+nemp3,nfil3+nemp3)
      double precision sum_i4(npmax,nfil3+nemp3,nfil3+nemp3)

      double precision E1,E2,E3,E12,E23,E31
      double precision E12_0,E23_0,E31_0
      double precision F1,F2,F3,F12,F23,F31
      integer i_a(3),i_b(3),i_c(3)  !3- maximal number of polarizations
      integer ia,ib,ic
      double precision anor_a,anor_b,anor_c
      double precision gam, gam1,rat1
      double precision Delta

c polarizations

c zinc-blende
C       data i_a/1,2,3/
C       data i_b/2,3,1/
C       data i_c/3,1,2/

c hexagonal
c       data i_a/3,3,1/
c       data i_b/3,1,3/
c       data i_c/3,1,1/

c trigonal (some components, 15R)
c       data i_a/2,2,1/
c       data i_b/2,1,1/
c       data i_c/2,1,2/

c    LiGaO2-wurtzite  ALSO: ZnGeN2, CdGeN2
c       data i_a/3,3,3/
c       data i_b/3,1,2/
c       data i_c/3,1,2/

      double precision c
      complex*16 opttt

      double precision opt1(3,2)
      integer igrp,j,kk,i,k,l1,m1,l2,ibas,m2,li,lf
      double precision optrr(3),optii(3)

      c=274.072d0

        if(nkp > kk3)stop'Increase kk3 at shg_sme.f'
        if((nfiup-nfilo+1) > nfil3)stop'Increase nfil3 at
     .          shg_smed.f'
        if((nemup-nemlo+1) > nemp3)stop'Increase nemp3 at
     .          shg_smed.f'
        if(nemlo <= nfiup)stop'Check nemlo,nfiup; SHGD program is built
     .          for sem. only'

c       print*,' shg_sta started!!! '
c       print*,' nfilo,nfiup,nemlo,nemup= ',nfilo,nfiup,nemlo,nemup

c  Fermi distribution
        do 801 j=1,nemup
        f(j)=0.d0
        if(j <= nfiup) f(j)=1.d0
  801   continue

        anorm=0.d0
        do 800 kk=1,3
        anorm=anorm+qp(kk)**2
  800   continue
        anorm=sqrt(anorm)

        if(anorm < 10.d0**(-5))then
        print*,' Gamma-point, Energies '
        print*,'E(k)= ',(eband(j,isp,ikp),j=1,nemup)
        print*,'f = ',(f(j),j=1,nemup)
        endif

C  Mat. el. interband transitions

         il=0
         do 4 i=nfilo,nemup
         il=il+1
         ir=0
         do 4 j=nfilo,nemup
         ir=ir+1

         if(ir < il)goto 4

c   electric dipole

            do 6 k=1,3

            opttt=(0.d0,0.d0)

            do 3 l1=0,nl-1
            do 3 m1=-l1,l1
            do 3 l2=0,nl-1

               if(l2 /= (l1+1). and. l2 /= (l1-1) ) goto 234
               if(k == 1)then
                  m2=m1+1
                  if(m2 > l2) goto 234
               endif
               if(k == 2)then
                  m2=m1-1
                  if(m2 < -l2)goto 234
               endif
               if(k == 3)m2=m1

               li=l1*l1+l1+m1
               lf=l2*l2+l2+m2

            jsp = 1
            do  ibas = 1, nbas
C           do  jsp = 1, nspc
            do  if2 = 1, nf2
            do  if1 = 1, nf1
              opttt = opttt +
     .          dconjg(ccp(li,ibas,jsp,if1,i))*ccp(lf,ibas,jsp,if2,j)
     .         *gradm(k,if1,if2,lf+1,li+1,isp,ipc(ibas))
            enddo
            enddo
            enddo


 234           continue
  3         continue

            opt1(k,1)= dreal(opttt)
            opt1(k,2)= dimag(opttt)

  6         continue
c
c  Cyclic coordinates everywhere in this program!!!
c  Return to Cartesian coordinates...
c
            optrr(1)= -(opt1(1,1)+opt1(2,1))/dsqrt(2.d0)
            optrr(2)= +(opt1(1,2)-opt1(2,2))/dsqrt(2.d0)
c           optrr(3)=   opt1(3,1)
            optrr(3)=  -opt1(3,1)

            optii(1)= -(opt1(1,2)+opt1(2,2))/dsqrt(2.d0)
            optii(2)= -(opt1(1,1)-opt1(2,1))/dsqrt(2.d0)
c           optii(3)=   opt1(3,2)
            optii(3)=  -opt1(3,2)

c real and imaginary parts of momentum operator

         ppp(1,il,ir,1)= optrr(1)
         ppp(2,il,ir,1)= optrr(2)
         ppp(3,il,ir,1)= optrr(3)

c        ppp(1,il,ir,1)= 1.d0
c        ppp(2,il,ir,1)= 1.d0
c        ppp(3,il,ir,1)= 1.d0

         ppp(1,il,ir,2)= optii(1)
         ppp(2,il,ir,2)= optii(2)
         ppp(3,il,ir,2)= optii(3)

c        ppp(1,il,ir,2)= 0.d0
c        ppp(2,il,ir,2)= 0.d0
c        ppp(3,il,ir,2)= 0.d0

c  ir <-> il, other half of matrix (gradient is anti-hermitian operator)

        if(ir /= il)then

         ppp(1,ir,il,1)=-optrr(1)
         ppp(2,ir,il,1)=-optrr(2)
         ppp(3,ir,il,1)=-optrr(3)

c        ppp(1,ir,il,1)=-1.d0
c        ppp(2,ir,il,1)=-1.d0
c        ppp(3,ir,il,1)=-1.d0

         ppp(1,ir,il,2)= optii(1)
         ppp(2,ir,il,2)= optii(2)
         ppp(3,ir,il,2)= optii(3)

c        ppp(1,ir,il,2)= 0.d0
c        ppp(2,ir,il,2)= 0.d0
c        ppp(3,ir,il,2)= 0.d0

        endif

  4      continue

c  cutoff parameters:
        gam=0.001d0
        gam1=0.01d0
        rat1=10.d0

        Delta = esciss

C       Delta=0.07d0
c       Delta=0.0d0      !!!No Gap "scissors"
c GaAs, GaP
c       Delta=0.095d0    !!!Gap "scissors" (GaAs)
c       Delta=0.095d0    !!!Gap "scissors" (GaP)
c SiC, GaN, AlN, BN
c       Delta=0.07d0     !!!Gap "scissors" (SiC,GaN)
c       Delta=0.14526d0  !!!Gap "scissors" (AlN,BN)
c Chalcopyrites
c       Delta=0.035d0    !!!Gap "scissors" (ZnGeP2)
c       Delta=0.018d0    !!!Gap "scissors" (CdGeAs2- ideal)
c       Delta=0.1d0      !!!Gap "scissors" (CdGeP2 - real, from LDA)
c       Delta=0.035d0    !!!Gap "scissors" (CdGeP2 - real, from Ham.)
c       Delta=0.07d0     !!!Gap "scissors" (ZnSiAs2- real, from LDA)
c       Delta=0.035d0    !!!Gap "scissors" (ZnSiP2 - real, LDA indir.)
c       Delta=0.07d0     !!!Gap "scissors" (ZnSiP2 - real, LDA direct)
c LiGaO2
c       Delta=0.133d0    !!!Gap "scissors" (LiGaO2)
c ZnGeN2
c       Delta=0.05d0     !!!Gap "scissors" (ZnGeN2)

c======================================
c  Integrand for SHG at a given k-point
c======================================

         do 300 k=1,npmax
c
         do 300 i1=1,nfil3
         do 300 i2=1,nemp3
c
         pmt_ee (k,i1,i2,ikp)=0.d0
         pmt_ei3(k,i1,i2,ikp)=0.d0
         qmt_ee (k,i1,i2,ikp)=0.d0
         qmt_ei3(k,i1,i2,ikp)=0.d0
c
  300    continue
c


         do 301 k=1,npmax
c
         do 301 i1=1,nfil3+nemp3
         do 301 i2=1,nfil3+nemp3
c
         sum_e1(k,i1,i2) =0.d0
         sum_e2(k,i1,i2) =0.d0
         sum_e3(k,i1,i2) =0.d0
c
         sum_i1(k,i1,i2) =0.d0
         sum_i2(k,i1,i2) =0.d0
         sum_i3(k,i1,i2) =0.d0
         sum_i4(k,i1,i2) =0.d0
c
  301    continue
c

c ============================
c   Three-Band Terms:
c ============================
         i1=0
         do 5 j1=nfilo,nemup
         i1=i1+1

         i2=0
         do 5 j2=nfilo,nemup
         i2=i2+1

         i3=0
         do 5 j3=nfilo,nemup
         i3=i3+1

c  all j1,j2,j3 should be different

         if(((j1-j2)*(j2-j3)*(j3-j1)) == 0) goto 5
c
c  Energy Factors
c
        E1 =eband(j1,isp,ikp)
        E2 =eband(j2,isp,ikp)
        E3 =eband(j3,isp,ikp)

        E12_0=eband(j1,isp,ikp)-eband(j2,isp,ikp)
        E23_0=eband(j2,isp,ikp)-eband(j3,isp,ikp)
        E31_0=eband(j3,isp,ikp)-eband(j1,isp,ikp)

        E12=eband(j1,isp,ikp)-eband(j2,isp,ikp)+Delta*(f(j2)-f(j1))
        E23=eband(j2,isp,ikp)-eband(j3,isp,ikp)+Delta*(f(j3)-f(j2))
        E31=eband(j3,isp,ikp)-eband(j1,isp,ikp)+Delta*(f(j1)-f(j3))

        if(abs(E12) < gam.or.abs(E23) < gam.or.abs(E31) < gam)
     .  goto 5
c       if(abs(E31-E23) < gam)goto 5
        if(abs(E31-E23) < 0.02)goto 5

c  Fermi Factors

        F1 =f(j1)
        F2 =f(j2)
        F3 =f(j3)

        F12=f(j1)-f(j2)
        F23=f(j2)-f(j3)
        F31=f(j3)-f(j1)
c
        do 501 i=1,3
c
         optr_a(i)=ppp(i,i1,i2,1)
         opti_a(i)=ppp(i,i1,i2,2)
c
         optr_b(i)=ppp(i,i2,i3,1)
         opti_b(i)=ppp(i,i2,i3,2)
c
         optr_c(i)=ppp(i,i3,i1,1)
         opti_c(i)=ppp(i,i3,i1,2)
c
 501    continue
c
        anor_a=0.d0
        anor_b=0.d0
        anor_c=0.d0
c
        do 502 i=1,3
        anor_a=anor_a+(optr_a(i)**2+opti_a(i)**2)
        anor_b=anor_b+(optr_b(i)**2+opti_b(i)**2)
        anor_c=anor_c+(optr_c(i)**2+opti_c(i)**2)
 502    continue
c
        if(abs(E12) <= gam1.and.sqrt(anor_a)/abs(E12) > rat1)then
        do 503 i=1,3
        optr_a(i)=0.d0
        opti_a(i)=0.d0
 503    continue
        endif
c
        if(abs(E23) <= gam1.and.sqrt(anor_b)/abs(E23) > rat1)then
        do 504 i=1,3
        optr_b(i)=0.d0
        opti_b(i)=0.d0
 504    continue
        endif
c
        if(abs(E31) <= gam1.and.sqrt(anor_c)/abs(E31) > rat1)then
        do 505 i=1,3
        optr_c(i)=0.d0
        opti_c(i)=0.d0
 505    continue
        endif
c
        do 302 k=1,npmax
        sum_gr1(k)=0.d0
        sum_gr2(k)=0.d0
        sum_gr3(k)=0.d0
  302   continue

        do 31 igrp=1,ngrp

        call grpop(optr_a,optr_ag,g,igrp)
        call grpop(opti_a,opti_ag,g,igrp)
        do 201 kk=1,3
        opt_ag(kk)=cmplx(optr_ag(kk),opti_ag(kk))
  201   continue

        call grpop(optr_b,optr_bg,g,igrp)
        call grpop(opti_b,opti_bg,g,igrp)
        do 202 kk=1,3
        opt_bg(kk)=cmplx(optr_bg(kk),opti_bg(kk))
  202   continue

        call grpop(optr_c,optr_cg,g,igrp)
        call grpop(opti_c,opti_cg,g,igrp)
        do 203 kk=1,3
        opt_cg(kk)=cmplx(optr_cg(kk),opti_cg(kk))
  203   continue

c  chi^2_{abc} - Matrix Elements for Frequency-Dependent

        do 32 k=1,npmax
c
        ia=iabc(1,k)
        ib=iabc(2,k)
        ic=iabc(3,k)
c
        sum_gr1(k)=sum_gr1(k) +1.d0/ngrp*
     .           (  dreal(opt_ag(ia)*opt_bg(ib)*opt_cg(ic))
     .             +dreal(opt_ag(ia)*opt_bg(ic)*opt_cg(ib))  )

        sum_gr2(k)=sum_gr2(k) +1.d0/ngrp*
     .           (  dreal(opt_cg(ia)*opt_bg(ib)*opt_ag(ic))
     .             +dreal(opt_cg(ia)*opt_bg(ic)*opt_ag(ib))  )

        sum_gr3(k)=sum_gr3(k) +1.d0/ngrp*
     .           (  dreal(opt_bg(ia)*opt_cg(ib)*opt_ag(ic))
     .             +dreal(opt_bg(ia)*opt_cg(ic)*opt_ag(ib))  )

  32    continue
c
  31    continue

        do 303 k=1,npmax

c Interband:

c e1 (three band part, 2w-term):

        sum_e1(k,i1,i2)=sum_e1(k,i1,i2) !!! i1,i2<=>n,m -ini. and fin.
     .          +2.d0*F12/(E31-E23)
     .          *sum_gr1(k) /(E12_0*E23_0*E31_0)

c e2 (three band part,  w-term):

        sum_e2(k,i1,i3)=sum_e2(k,i1,i3) !!! i1,i3<=>n,l -ini. and fin.
     .          +F31/(E31-E23)
     .          *sum_gr1(k)/(E12_0*E23_0*E31_0)

c e3 (three band part,  w-term):

        sum_e3(k,i3,i2)=sum_e3(k,i3,i2) !!! i3,i2<=>l,m -ini. and fin.
     .          +F23/(E31-E23)
     .          *sum_gr1(k)/(E12_0*E23_0*E31_0)

c Intraband (three band part)

c i1 (three band part, 2w-term):

        sum_i1(k,i1,i2)=sum_i1(k,i1,i2) !!! i1,i2<=>n,m -ini. and fin.
     .          -2.d0*(E31-E23)
     .          *sum_gr1(k)*F12 /(E12_0*E23_0*E31_0)      !!!/E12**2

c i2 (three band part,  w-term):

        sum_i2(k,i1,i2)=sum_i2(k,i1,i2) !!! i1,i2<=>n,m -ini. and fin.
     .          +0.5d0*(E23*sum_gr2(k)-E31*sum_gr3(k))
     .          /(E12_0*E23_0*E31_0)*F12                 !!!/E12**2

c i3 (three band part,  w-term):

        sum_i3(k,i1,i3)=sum_i3(k,i1,i3) !!! i1,i3<=>n,l -ini. and fin.
     .          +sum_gr1(k)/(E12_0*E23_0*E31_0)*F31*E12  !!!/E31**2

c i4 (three band part,  w-term):

        sum_i4(k,i3,i2)=sum_i4(k,i3,i2) !!! i3,i2<=>l,m -ini. and fin.
     .          -sum_gr1(k)/(E12_0*E23_0*E31_0)*F23*E12  !!!/E23**2

 303    continue
c
 5      continue

c!!! Three Band Terms only!
        goto 1000

c ===================
c  Intraband - Two-Band Terms:
c ===================
         i1=0
         do 1005 j1=nfilo,nemup
         i1=i1+1

         i2=0
         do 1005 j2=nfilo,nemup
         i2=i2+1

c  j1,j2 should be different

         if((j1-j2) == 0) goto 1005
c
c  Energy Factors
c
        E1 =eband(j1,isp,ikp)
        E2 =eband(j2,isp,ikp)

        E12_0=eband(j1,isp,ikp)-eband(j2,isp,ikp)

        E12=eband(j1,isp,ikp)-eband(j2,isp,ikp)+Delta*(f(j2)-f(j1))

        if(abs(E12) < gam)goto 1005

c  Fermi Factors

        F1 =f(j1)
        F2 =f(j2)
        F12=f(j1)-f(j2)

        do  i=1,3
c
         optr_a(i)=0.d0 !!! ppp(i,i2,i2,1)-ppp(i,i1,i1,1)
         opti_a(i)=ppp(i,i2,i2,2)-ppp(i,i1,i1,2)
c
         optr_b(i)=ppp(i,i1,i2,1)
         opti_b(i)=ppp(i,i1,i2,2)
c
         optr_c(i)=ppp(i,i2,i1,1)
         opti_c(i)=ppp(i,i2,i1,2)
c
       enddo

        do 1302 k=1,npmax
        sum_gr4(k)=0.d0
1302    continue

        do 1031 igrp=1,ngrp

        call grpop(optr_a,optr_ag,g,igrp)
        call grpop(opti_a,opti_ag,g,igrp)
        do 1201 kk=1,3
        opt_ag(kk)=cmplx(optr_ag(kk),opti_ag(kk))
1201    continue

        call grpop(optr_b,optr_bg,g,igrp)
        call grpop(opti_b,opti_bg,g,igrp)
        do 1202 kk=1,3
        opt_bg(kk)=cmplx(optr_bg(kk),opti_bg(kk))
1202    continue

        call grpop(optr_c,optr_cg,g,igrp)
        call grpop(opti_c,opti_cg,g,igrp)
        do 1203 kk=1,3
        opt_cg(kk)=cmplx(optr_cg(kk),opti_cg(kk))
1203    continue

        do 1032 k=1,npmax
c
        ia=i_a(k)
        ib=i_b(k)
        ic=i_c(k)

        sum_gr4(k)=sum_gr4(k) +1.d0/ngrp*
     .           (  dreal(opt_bg(ia)*opt_cg(ib)*opt_ag(ic))
     .             +dreal(opt_bg(ia)*opt_cg(ic)*opt_ag(ib))  )

1032    continue
c
1031    continue

        do 1303 k=1,npmax

c Intraband     DODELAT"!!!!!!!!!!!!!!!!!

c i1 (two+three band part, 2w-term):

        sum_i1(k,i1,i2)= sum_i1(k,i1,i2)  !!! i1,i2<=>n,m -ini. and fin.
     .          -8.0d0*sum_gr4(k)*F12   !!!/E12**4

c i2 (two+three band part,  w-term):

        sum_i2(k,i1,i2)= sum_i2(k,i1,i2)  !!! i1,i2<=>n,m -ini. and fin.
     .          -0.5d0*sum_gr4(k)*F12   !!!/E12**4

1303    continue

1005    continue

c!!! Three Band Terms only!
1000    continue
c         print*, 'OK 5'
c ===================
c Total contribution
c ===================

        i1=0
        do 802 j1=nfilo,nfiup
        i1=i1+1

        i2=0
        do 802 j2=nemlo,nemup
        i2=i2+1

        do 802 k=1,npmax

c 2*w term
        pmt_ee (k,i1,i2,ikp)=
     .          +sum_e1(k,i1,i2+nemlo-nfilo)  !!!interband
c
        pmt_ei3(k,i1,i2,ikp)=
     .          +sum_i1(k,i1,i2+nemlo-nfilo)  !!!intraband

c 1*w term
        qmt_ee (k,i1,i2,ikp)=
     .          +sum_e2 (k,i1,i2+nemlo-nfilo) !!!interband
     .          +sum_e3 (k,i1,i2+nemlo-nfilo)
c
        qmt_ei3(k,i1,i2,ikp)=
     .          +sum_i2 (k,i1,i2+nemlo-nfilo) !!!intraband
     .          +sum_i3 (k,i1,i2+nemlo-nfilo)
     .          +sum_i4 (k,i1,i2+nemlo-nfilo)

 802    continue

        return
C  99   stop ' OPWF read error'
        end

