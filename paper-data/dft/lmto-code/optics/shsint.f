         subroutine shsint(n1,n2,n3,nkp,nsp,ntet,idtet,
     .          vol_cell,esciss,nfilo,nfiup,nemlo,nemup,
     .          nbmax,eband,nchi2,iabc)
C-  BZ integration of static SHG by linear method
C ----------------------------------------------------------------------
Ci Inputs:
Ci  ntet, No. of different tetrahedra
Ci  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
Ci  irreducible k-points:
Ci  idtet(0,i), no. of tetrahedra of the i'th kind
Ci  nkp, no. of irreducible k-points; nb, no. of bands; nsp, see BNDASA;
Co Outputs:
Co  shg_stat, Integrated SHG for different tensor components
Cr Remarks:
C   none
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,nkp,nsp,ntet,nfilo,nfiup,nemlo,nemup,nbmax,nchi2
      integer iabc(3,6),idtet(0:4,*)
      double precision vol_cell,esciss
      double precision eband(nbmax,nsp,nkp)
C ... Local parameters
      include 'PAR_OPT'

      double precision shg_me(kk2,6) !!! 1-{123} zincblende
      double precision shg_me_ee(kk2,6)
      double precision shg_me_ei2(kk2,6)
c     double precision shg_me_ei3(kk2,6)
      double precision shg_me_ei3(kk2,6,nfil2,nemp2)

      common/shg_s/ shg_me,shg_me_ee,shg_me_ei2,shg_me_ei3

      double precision shg_stat(6) !!! 1-{123} zincblende
      double precision shg_stat_ee(6)
      double precision shg_stat_ei2(6)
      double precision shg_stat_ei3(6)

      double precision ds2,ds3
      double precision eci(4),ecj(4)
      double precision chucha
      double precision edif(4),e_mat(4,4)
      double precision Delta

      integer iq1,iq2,iq3,iq4,itet,k,ibf,ib,ibm,jb,ii,jj,i
      double precision volwgt,wtm,wt
      double precision wtm_ee,wtm_ei2,wtm_ei3
      double precision wt_ee,wt_ei2,wt_ei3
      double precision shgsi
      double precision pi,shgesu

      pi = 4d0*datan(1d0)

      if (nsp == 2) call rx('shg_sin not spin polarized')
        if(nkp > kk2)stop'Increase kk2 at shsint.f'
        if((nfiup-nfilo+1) > nfil2)stop'Increase nfil2 at
     .          shsint.f'
        if((nemup-nemlo+1) > nemp2)stop'Increase nemp2 at
     .          shsint.f'

        volwgt = dble(3-nsp)/(n1*n2*n3*6)

        Delta = esciss


        do 101 k=1,6
        shg_stat(k)=0.d0
        shg_stat_ee(k)=0.d0
        shg_stat_ei2(k)=0.d0
        shg_stat_ei3(k)=0.d0
  101   continue

C --- Loop over tetrahedra ---
        do  11  itet = 1, ntet
         iq1 = idtet(1,itet)
         iq2 = idtet(2,itet)
         iq3 = idtet(3,itet)
         iq4 = idtet(4,itet)

        do 9 k=1,6  !!! Polarizations, different SHG components:1-{123} zb

        wtm=
     .  0.25d0*(shg_me(iq1,k)+shg_me(iq2,k)+shg_me(iq3,k)+shg_me(iq4,k))

        wtm_ee=
     .  0.25d0*(shg_me_ee(iq1,k)+shg_me_ee(iq2,k)
     .         +shg_me_ee(iq3,k)+shg_me_ee(iq4,k))

        wtm_ei2=
     .  0.25d0*(shg_me_ei2(iq1,k)+shg_me_ei2(iq2,k)
     .         +shg_me_ei2(iq3,k)+shg_me_ei2(iq4,k))

c       wtm_ei3=
c    .  0.25d0*(shg_me_ei3(iq1,k)+shg_me_ei3(iq2,k)
c    .         +shg_me_ei3(iq3,k)+shg_me_ei3(iq4,k))

        wt     = volwgt*idtet(0,itet)*wtm
        wt_ee  = volwgt*idtet(0,itet)*wtm_ee
        wt_ei2 = volwgt*idtet(0,itet)*wtm_ei2


c       wt_ei3 = volwgt*idtet(0,itet)*wtm_ei3

        shg_stat(k)    =shg_stat(k)     + wt
        shg_stat_ee(k) =shg_stat_ee(k)  + wt_ee
        shg_stat_ei2(k)=shg_stat_ei2(k) + wt_ei2
c       shg_stat_ei3(k)=shg_stat_ei3(k) + wt_ei3

   9    continue   !!! polarizations

   11   continue   !!! tetrahedra

c New test for ei3:

c       goto 1000

        ds2=dsqrt(2.d0)*10.d0**(-6)
        ds3=dsqrt(3.d0)*10.d0**(-6)

C --- Loop over initial bands
        ibf=0
        do  20  ib = nfilo,nfiup
        ibf=ibf+1

C --- Loop over final bands
        ibm=0
        do  20  jb = nemlo,nemup
        ibm=ibm+1

C --- Loop over tetrahedra ---
        do  211  itet = 1, ntet
        iq1 = idtet(1,itet)
        iq2 = idtet(2,itet)
        iq3 = idtet(3,itet)
        iq4 = idtet(4,itet)

C ...      Set up energies at 4 corners of tetrahedron for ib-ini
           eci(1) = eband(ib,1,iq1) + ds2*0.d0
           eci(2) = eband(ib,1,iq2) + ds2*1.d0
           eci(3) = eband(ib,1,iq3) + ds2*2.d0
           eci(4) = eband(ib,1,iq4) + ds2*3.d0

C ...        Set up energies at 4 corners of tetrahedron for jb-fin
             ecj(1) = eband(jb,1,iq1) + ds3*0.d0 + Delta !!!Gap "scissors"
             ecj(2) = eband(jb,1,iq2) + ds3*1.d0 + Delta
             ecj(3) = eband(jb,1,iq3) + ds3*2.d0 + Delta
             ecj(4) = eband(jb,1,iq4) + ds3*3.d0 + Delta

        do 123 k=1,6  !!!polarizations

        wtm_ei3=
     .  0.25d0*(shg_me_ei3(iq1,k,ibf,ibm)+shg_me_ei3(iq2,k,ibf,ibm)
     .         +shg_me_ei3(iq3,k,ibf,ibm)+shg_me_ei3(iq4,k,ibf,ibm))

c Sostavit' chuchu iz 4-h energii s Log
        do 220 ii=1,4
        edif(ii)=ecj(ii)-eci(ii)
 220    continue
        do 221 ii=1,4
        do 221 jj=1,4
        e_mat(ii,jj)=edif(ii)-edif(jj)
 221    continue

c
        chucha=
     .          +0.5d0/e_mat(1,2)/e_mat(1,3)/e_mat(1,4)*dlog(edif(1))
     .          +0.5d0/e_mat(2,1)/e_mat(2,3)/e_mat(2,4)*dlog(edif(2))
     .          +0.5d0/e_mat(3,1)/e_mat(3,2)/e_mat(3,4)*dlog(edif(3))
     .          +0.5d0/e_mat(4,1)/e_mat(4,2)/e_mat(4,3)*dlog(edif(4))

c
c Znak - iz-za togo, chto E12**3 vybroshena, pri F12 ona otricatel'na.
c 6 - iz-za tetrahedron volume
        wt_ei3 =-volwgt*idtet(0,itet)*wtm_ei3*chucha *6.d0* 2.d0

        shg_stat_ei3(k)=shg_stat_ei3(k) + wt_ei3

  123    continue   !!! polarizations

  211    continue   !!! tetrahedra
c
  20     continue   !!! initial and final bands

 1000   continue

c
c  output!!!
c

C       print*,' Total: ** scale by -1 with new def sharm'
C       Revert to original def of sharm ... remove scaling by -1
        do 113 k=1,nchi2                !!! Polarizations
           shgesu= (shg_stat_ee(k)+shg_stat_ei3(k))
     .          /vol_cell*0.5d0*32.d0*5.834d0
           shgsi=shgesu*4*pi/3
        print '(" abc=",3i2," shg_stat=",f10.5," E-8 esu"," = ",f10.5," pm/V")'
     . ,(iabc(i,k),i=1,3),shgesu*(1),shgsi*(1)
  113   continue

        print*,' ee: '
        do 114 k=1,nchi2                !!! Polarizations
        shgesu=shg_stat_ee(k)/vol_cell*0.5d0*32.d0*5.834d0
        shgsi=shgesu*4*pi/3
        print '(" abc=",3i2," shg_stat_ee=",f10.5," E-8 esu"," = ",f10.5," pm/V")'
     . ,(iabc(i,k),i=1,3),shgesu*(1),shgsi*(1)
  114   continue

        print*,' ei2: '
        do 115 k=1,nchi2                !!! Polarizations
        shgesu=shg_stat_ei2(k)/vol_cell*0.5d0*32.d0*5.834d0
        shgsi=shgesu*4*pi/3
        print '(" abc=",3i2," shg_stat_ei2=",f10.5," E-8 esu"," = ",f10.5," pm/V")'
     . ,(iabc(i,k),i=1,3),shgesu*(1),shgsi*(1)
  115   continue

        print*,' ei3: '
        do 116 k=1,nchi2                !!! Polarizations
        shgesu=shg_stat_ei3(k)/vol_cell*0.5d0*32.d0*5.834d0
        shgsi=shgesu*4*pi/3
        print '(" abc=",3i2," shg_stat_ei3=",f10.5," E-8 esu"," = ",f10.5," pm/V")'
     . ,(iabc(i,k),i=1,3),shgesu*(1),shgsi*(1)
  116   continue

        end
