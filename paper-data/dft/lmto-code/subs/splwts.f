      subroutine splwts(nqp,nband,nbmx,nsp,wgts,evl,n,w,efermi,
     .                  metal,sumev,bndwts,wtot,entrpy,dosef,cv)
C- Make sampling weights for integrals under the Fermi surface
C-----------------------------------------------------------------------
Ci  Input
Ci    nqp : number of q-points; nband : number of bands
Ci    wgts: band weights
Ci    evl : energy eigenvalues
Ci    n   : n>0 Methfessel-Paxton polynomial order
Ci        : n<0 sampling done with Fermi-Dirac statistics
Ci    w   : n>0 gaussian width in Methfessel-Paxton integration (Ry)
Ci        : n<0 Temperature for Fermi distribution (Ry)
Ci    nbmx : first dimension of evl ;
Ci    metal : if F, weights unity below E_f and zero above.
Ci    efermi : Fermi energy
Co  Output
Co    bndwts : band and E_F - dependent k-point weights for integration
Co    wtot   : sum of all weights (charge)
Co    entrpy : electron entropy
Co    dosef  : DOS at Fermi energy
Co    cv     : electronic specific heat
Co             (only evaluated with Fermi-Dirac statistics)
Cr  Remarks
Cr    sum of occupied eigenvalues = sum_n,k  w_nk E_nk
Cr    w_nk are generalised occupation numbers;
Cr    see Needs et al. Phys Rev B 33 (1986) 3778, eqs 1 & 2.
Cu Updates
Cu   16 Jul 08 returns entropy as TS for all n
Cu   04 Aug 07 Generates dos(efermi), cv(T=w) for F-D statistics
Cu   02 May 07 (MvS) prints entropy to stdout
Cu   21 Jun 06 (ATP) generates entrpy as output
Cu   17 Jan 05 Output wtot
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nqp,nband,nbmx,nsp,n
      logical metal
      double precision wgts(nqp),evl(nbmx,nsp,nqp),w,efermi,sumev,
     .                 bndwts(nband,nsp,nqp),wtot,entrpy,dosef,cv
C ... Local parameters
      integer iqp,iband,isp,stdo
      double precision e,s,d,wt,x,xx,dsdt,tdsdt
      procedure(integer) :: iprint,i1mach,nglob

      stdo = nglob('stdo')
      sumev=0d0; wtot=0d0; entrpy=0d0; tdsdt=0; dsdt=0; dosef=0
      do  iqp = 1, nqp
        do  iband = 1, nband
          do  isp = 1, nsp
            e = evl(iband,isp,iqp)
            if (metal) then

C             debugging: check derivative numerically
C              x = (efermi - e) / (w+1d-7)
C              call delstp(n,x,d,s,sp)
C              x = (efermi - e) / (w-1d-7)
C              call delstp(n,x,d,s,sm)
C              dsdt1 = (sp-sm)/2d-7

              x = (efermi - e) / w
C             call delstp(n,x,d,s,xx)
              call delstd(n,x,d,s,xx,dsdt)
              if (abs(x) < 36) then
                dsdt = -(efermi-e)/w**2 * dsdt
              else
                dsdt = 0
              endif

CC             debugging: compare analytical, numerical derivative
C              if (abs(x) < 30) then
C                print 222, x,dsdt1,dsdt,dsdt-dsdt1
C  222           format(3f14.8,1pe12.3)
C              endif
            else
              s = 1d0
              if (e <= efermi) s = 0d0
              xx = 0
              d = 0
            endif
            wt = abs(wgts(iqp)) * (1d0 - s) / nsp
            bndwts(iband,isp,iqp) = wt
            dosef = dosef + d*abs(wgts(iqp))/w/nsp
            wtot = wtot + wt
            sumev = sumev + e * wt
            entrpy = entrpy + xx  * abs(wgts(iqp)) / nsp
            tdsdt  = tdsdt + dsdt * abs(wgts(iqp)) / nsp
          enddo
        enddo
      enddo
      tdsdt = tdsdt*w
      entrpy = entrpy*w
      if (n < 0) cv = tdsdt

C ... Print out band weights, if only 1 kp
      if (iprint() > 30 .and. nqp == 1) then
        do  isp = 1, nsp
          call info2(30,0,0,' SPLWTS: band weights .. '//
     .      '%?;(n==2);Spin %i;;%N       eval      weight',nsp,isp)
          do  iband = 1, nband
            write (stdo,20) evl(iband,isp,1),bndwts(iband,isp,1)
   20       format (4x,2f10.6)
          enddo
        enddo
      endif
C ... Print out various k-point integrations
      if (iprint() >= 10) then
        if (n >= 0) then
          call awrit6(' N=%i, W=%d, E_F=%d, sumev=%d, entropy term:'
     .                //' %d, %d electrons',' ',256,i1mach(2),
     .                n,w,efermi,sumev,entrpy,wtot)
        else
          call awrit5(' T=%dK, E_F=%d, sumev=%d, TS=%;3g,'
     .                 //' %d electrons',' ',256,i1mach(2),
     .                 0.1579d6*w,efermi,sumev,entrpy,wtot)
          call info5(30,0,0,
     .      '%fEntropy S=%;6d k_B,  specific heat TdS/dT=%;6d k_B'//
     .      ' (%?;n<0;Fermi-Dirac;sampling;)',entrpy/w,tdsdt,n,0,0)
        endif
C      call info5(10,0,0,' SPLWTS: Fermi energy:%;6d;'//
C     .  '  band energy=%;6d;  %;6d electrons  DOS(E_f)=%;4g',
C     .    efermi,sumev,wtot,dosef,0)
      endif
      end
