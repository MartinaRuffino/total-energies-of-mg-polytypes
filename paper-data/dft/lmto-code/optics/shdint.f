         subroutine shdint(ifio,n1,n2,n3,nkp,nsp,
     .          emin,emax,efermi,npts,ntet,idtet,
     .          vol_cell,esciss,
     .          nfilo,nfiup,nemlo,nemup,nbmax,eband)
C- BZ integration of frequency dependent SHG
C ----------------------------------------------------------------------
Ci Inputs:
Ci  eband, energy bands;
Ci  ntet, No. of different tetrahedra
Ci  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
Ci  irreducible k-points:
Ci  idtet(0,i), no. of tetrahedra of the i'th kind
Ci  nkp, no. of irreducible k-points; nb, no. of bands; nsp, see BNDASA;
Ci  emin, emax, npts: energy window and npts of points (see remarks)
Co Outputs:
Co  dos, Integrated Joint Density of States (idos)
Cr Remarks
Cr   Adapted from bzints to make joint density of states.
Cr   All energy differences between states below ef and states
Cr   above ef+emin are summed, and integrated over the BZ
Cr   Treatment near the critical points (ef and ef+emin) handled crudely
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifio,n1,n2,n3,nkp,nsp,npts,ntet,vol_cell,
     .  nfilo,nfiup,nemlo,nemup,nbmax,idtet(0:4,*)
      double precision emin,emax,efermi(2)
C ... Local parameters
      include 'PAR_OPT'
      double precision qmt_ee (npmax,nfil3,nemp3,kk3) !w
      double precision qmt_ei3(npmax,nfil3,nemp3,kk3) !w
      double precision pmt_ee (npmax,nfil3,nemp3,kk3) !2w
      double precision pmt_ei3(npmax,nfil3,nemp3,kk3) !2w
      common/shgmt/ qmt_ee,qmt_ei3,pmt_ee,pmt_ei3

      integer ibf,ibm

      double precision eband(nbmax,nsp,nkp),
     .  qint_ee (npmax,npts),pint_ee (npmax,npts),
     .  qint_ei3(npmax,npts),pint_ei3(npmax,npts),
     .  chi_ee  (npmax,npts),chi_ei3 (npmax,npts),
     .  chi_2w  (npmax,npts),chi_1w  (npmax,npts)
      double precision pi
      double precision esciss

c partial
      double precision shg_par(npts),shg_par_ee(npts),
     .  shg_par_ei3(npts)
      double precision f(nbmax) !!! Fermi distribution function
      double precision Delta

      integer ib,jb,iq1,iq2,iq3,iq4,itet,j,k,ipts,nfilem4
      double precision eci(4),ecj(4),volwgt,ef,
     .  wtmq_ee,wtmp_ee,wtmq_ei3,wtmp_ei3,
     .  wtq_ee, wtp_ee, wtq_ei3, wtp_ei3
      double precision yq_ee (npts),yp_ee (npts),
     .  yq_ei3(npts),yp_ei3(npts)
      double precision de,ei
      double precision ds2,ds3

C        print*, '3rd check'
      if(nkp > kk3)stop'Increase kk3 at shg_smed.f'
      if((nfiup-nfilo+1) > nfil3)stop'Increase nfil3 at
     .  shg_smed.f and shg_sind.f'
      if((nemup-nemlo+1) > nemp3)stop'Increase nemp3 at
     .  shg_smed.f and shg_sind.f'

      pi = 4d0*datan(1d0)

      volwgt = dble(3-nsp)/(n1*n2*n3*6)

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

      ds2=dsqrt(2.d0)*10.d0**(-6)
      ds3=dsqrt(3.d0)*10.d0**(-6)

      ef=efermi(1) + 0.001d0

      print*,' emin,emax,npts= ',emin,emax,npts
      print*,' vol_cell(a.u.)= ',vol_cell
      print*,' volwgt,nbmax= ',volwgt,nbmax
      print*,' nfilo,nfiup= ',nfilo,nfiup
      print*,' nemlo,nemup= ',nemlo,nemup
      print*,' ef= ',ef
      print*,' 6*n1*n2*n3, ntet,nsp = ',6*n1*n2*n3,ntet,nsp

c  Fermi distribution
      do 801 j=1,nemup
        f(j)=0.d0
        if(j <= nfiup) f(j)=1.d0
  801 continue

      do 111 k=1,npmax
      do 111 ipts=1,npts
c
        qint_ee (k,ipts)=0.d0
        pint_ee (k,ipts)=0.d0
c
        qint_ei3(k,ipts)=0.d0
        pint_ei3(k,ipts)=0.d0
c
        chi_ee  (k,ipts)=0.d0
        chi_ei3 (k,ipts)=0.d0
        chi_2w  (k,ipts)=0.d0
        chi_1w  (k,ipts)=0.d0
c
  111   continue

C --- Loop over initial bands
        ibf=0
        do  20  ib = nfilo,nfiup
          ibf=ibf+1

C --- Loop over final bands
          ibm=0
          do  20  jb = nemlo,nemup
            ibm=ibm+1


c partial
            do 21 ipts=1,npts
              shg_par_ee(ipts)=0.d0
              shg_par_ei3(ipts)=0.d0
   21       shg_par(ipts)=0.d0
c partial


c --- Fermi factor
c       F12=f(ib)-f(jb)
c       if(abs(F12-1.d0) > 0.000001d0)goto 20

C --- Loop over tetrahedra ---
            do  11  itet = 1, ntet
              iq1 = idtet(1,itet)
              iq2 = idtet(2,itet)
              iq3 = idtet(3,itet)
              iq4 = idtet(4,itet)
c       print*, 'iq= ',idtet(0,itet),iq1,iq2,iq3,iq4

C ...      Set up energies at 4 corners of tetrahedron for ib
              eci(1) = eband(ib,1,iq1) + ds2*0.d0
              eci(2) = eband(ib,1,iq2) + ds2*1.d0
              eci(3) = eband(ib,1,iq3) + ds2*2.d0
              eci(4) = eband(ib,1,iq4) + ds2*3.d0

C     ...        Set up energies at 4 corners of tetrahedron for jb
              ecj(1) = eband(jb,1,iq1) + ds3*0.d0 + Delta !!!Gap "scissors"
              ecj(2) = eband(jb,1,iq2) + ds3*1.d0 + Delta
              ecj(3) = eband(jb,1,iq3) + ds3*2.d0 + Delta
              ecj(4) = eband(jb,1,iq4) + ds3*3.d0 + Delta

c!!!TEST
c       do 116 iii=1,4
c       if((ecj(iii)-eci(iii)) < 0.1d0)then
c       ecj(iii)=eci(iii)+0.1d0
c       endif
c 116   continue
c!!!TEST

              do 123 k=1,npmax  !!! polarizations

c Interband:

                wtmq_ee=
     .            qmt_ee(k,ibf,ibm,iq1)+qmt_ee(k,ibf,ibm,iq2)+
     .            qmt_ee(k,ibf,ibm,iq3)+qmt_ee(k,ibf,ibm,iq4)

                wtmp_ee=
     .            pmt_ee(k,ibf,ibm,iq1)+pmt_ee(k,ibf,ibm,iq2)+
     .            pmt_ee(k,ibf,ibm,iq3)+pmt_ee(k,ibf,ibm,iq4)

                wtmq_ee=wtmq_ee/4.d0
                wtmp_ee=wtmp_ee/4.d0
c         wtmq_ee=1.d0
c         wtmp_ee=1.d0

                wtq_ee = volwgt*idtet(0,itet)*wtmq_ee
                wtp_ee = volwgt*idtet(0,itet)*wtmp_ee

                call opt0(eci,ecj,     emin,     emax,npts,yq_ee,ef,
     .            wtq_ee)
                call opt0(eci,ecj,2.d0*emin,2.d0*emax,npts,yp_ee,ef,
     .            wtp_ee)

                do 112 ipts=1,npts
                  qint_ee(k,ipts)=qint_ee(k,ipts)+yq_ee(ipts)
                  pint_ee(k,ipts)=pint_ee(k,ipts)+yp_ee(ipts)
  112           continue

c partial
                goto 251
                if(ib > 23.and.jb < 29)then

                  if(k == 3)then
                    do 212 ipts=1,npts
                      shg_par_ee(ipts)=shg_par_ee(ipts)+yp_ee(ipts)
  212               continue
                  endif

                endif
  251           continue

c Intraband (3 band terms):

                wtmq_ei3=
     .            qmt_ei3(k,ibf,ibm,iq1)+qmt_ei3(k,ibf,ibm,iq2)+
     .            qmt_ei3(k,ibf,ibm,iq3)+qmt_ei3(k,ibf,ibm,iq4)

                wtmp_ei3=
     .            pmt_ei3(k,ibf,ibm,iq1)+pmt_ei3(k,ibf,ibm,iq2)+
     .            pmt_ei3(k,ibf,ibm,iq3)+pmt_ei3(k,ibf,ibm,iq4)

                wtmq_ei3=wtmq_ei3/4.d0
                wtmp_ei3=wtmp_ei3/4.d0
c     wtmq_ei3=1.d0
c     wtmp_ei3=1.d0

                wtq_ei3 = volwgt*idtet(0,itet)*wtmq_ei3
                wtp_ei3 = volwgt*idtet(0,itet)*wtmp_ei3

                call
     .            opt0(eci,ecj,     emin,     emax,npts,yq_ei3,ef,
     .            wtq_ei3)
                call
     .            opt0(eci,ecj,2.d0*emin,2.d0*emax,npts,yp_ei3,ef,
     .            wtp_ei3)

                do 114 ipts=1,npts
                  qint_ei3(k,ipts)=qint_ei3(k,ipts)+yq_ei3(ipts)
                  pint_ei3(k,ipts)=pint_ei3(k,ipts)+yp_ei3(ipts)
  114           continue
c

c     partial
                goto 252
                if(ib > 23.and.jb < 29)then

                  if(k == 3)then
                    do 214 ipts=1,npts
                      shg_par_ei3(ipts)=shg_par_ei3(ipts)+yp_ei3(ipts)
  214               continue
                  endif

                endif
  252           continue

  123         continue          !!! polarizations
c
   11       continue            !!! tetrahedra
c

c     partial
            goto 253
            if(ib > 23.and.jb < 29)then

              stop 'nfilem4 not set'
              de=(emax-emin)/(npts-1)
              write(nfilem4,*)ib,jb
              do 213 ipts=1,npts
                ei=emin+(ipts-1)*de

                if(ei > 10.d0**(-5))then
                  shg_par_ee(ipts)=shg_par_ee(ipts)
     .              /vol_cell*pi*0.5d0*32.d0*5.834d0
                  shg_par_ei3(ipts)=shg_par_ei3(ipts)
     .              /vol_cell*pi*0.5d0*32.d0*5.834d0/4.d0/ei**2
                  shg_par(ipts)=shg_par_ee(ipts)+shg_par_ei3(ipts)

                  write(nfilem4,'(2f10.4)')ei,shg_par(ipts)
                endif

  213         continue

            endif

  253       continue

   20     continue              !!! initial and final bands

c
c  output!!!
c

c partial
c       goto 240
c partial

          de=(emax-emin)/(npts-1)
          do 113 ipts=1,npts
            ei=emin+(ipts-1)*de

            if(ei > 10.d0**(-5))then

              do 130 k=1,npmax

c     Interband - 3 band terms:

                qint_ee(k,ipts)=qint_ee(k,ipts)
     .            /vol_cell*pi*0.5d0*32.d0*5.834d0
                pint_ee(k,ipts)=pint_ee(k,ipts)
     .            /vol_cell*pi*0.5d0*32.d0*5.834d0

c     Intraband - 3 band terms:

                qint_ei3(k,ipts)=qint_ei3(k,ipts)
     .            /vol_cell*pi*0.5d0*32.d0*5.834d0/ei**2
                pint_ei3(k,ipts)=pint_ei3(k,ipts)
     .            /vol_cell*pi*0.5d0*32.d0*5.834d0/4.d0/ei**2

c     Total: Interband          - w+2w
                chi_ee (k,ipts)=qint_ee (k,ipts)+pint_ee (k,ipts)
c
c     Total: Intraband (3 band) - w+2w
                chi_ei3(k,ipts)=qint_ei3(k,ipts)+pint_ei3(k,ipts)

c     Total: 2w: Interband + Intraband
                chi_2w (k,ipts)=pint_ee (k,ipts)+pint_ei3(k,ipts)

c     Total: 1w: Interband + Intraband
                chi_1w (k,ipts)=qint_ee (k,ipts)+qint_ei3(k,ipts)

  130         continue

              write(ifio,'(7f10.4)')ei,
c    .          pint_ee (1,ipts)/10000.,pint_ee (2,ipts)/10000.,
c    .          pint_ee (3,ipts)/10000.,
c    .          pint_ei3(1,ipts)/10000.,pint_ei3(2,ipts)/10000.,
c    .          pint_ei3(3,ipts)/10000.

c     .         chi_ee (1,ipts),chi_ee (2,ipts),chi_ee (3,ipts),
c     .         chi_ei3(1,ipts),chi_ei3(2,ipts),chi_ei3(3,ipts)

     .          chi_2w (1,ipts),chi_2w (2,ipts),chi_2w (3,ipts),
     .          chi_1w (1,ipts),chi_1w (2,ipts),chi_1w (3,ipts)

c Interband - 3 band terms:

c    .  qint(1,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0,
c    .  qint(2,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0,
c    .  qint(3,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0,
c    .  pint(1,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0,
c    .  pint(2,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0,
c    .  pint(3,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0

c Intraband - 3 band terms:

c    .  qint(1,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/ei**2,
c    .  qint(2,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/ei**2,
c    .  qint(3,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/ei**2,
c    .  pint(1,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/4.d0/ei**2,
c    .  pint(2,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/4.d0/ei**2,
c    .  pint(3,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/4.d0/ei**2

c Intraband - 2 band terms:

c    .  qint(1,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/ei**4,
c    .  qint(2,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/ei**4,
c    .  qint(3,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/ei**4,
c    .  pint(1,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/16.d0/ei**4,
c    .  pint(2,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/16.d0/ei**4,
c    .  pint(3,ipts)/vol_cell*pi*0.5d0*32.d0*5.834d0/16.d0/ei**4
c
            endif

  113     continue

c     partial
C 240     continue

        end

