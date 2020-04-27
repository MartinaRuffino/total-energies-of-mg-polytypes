      subroutine madpot(nbas,ic1,ic2,nrc,ipc,dclabl,qt,vconst,rhrmx,
     .  rmax,dmad,bas,plat,vrl,noves,vrmax,ves,emad,trumad,vmtz,
     .  nangl,idcc,wts,etrms)
C- Calculate Madelung potential, Madelung energy, vmtz
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ic1,ic2: calculate Madelung potential for classes ic1 .. ic2
Ci   nrc   :number of memembers of each class
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   dclabl:class name, packed as a real number
Ci   qt    :electronic charge, less nuclear charge
Ci         :(excluding homogeneous nuclear charge in rhrmx(2))
Ci   vconst:constant potential shift added to estat potentials
Ci   rhrmx :constant background density (for Muffin-Tin correction)
Ci          rhrmx(1) = constant background density
Ci          rhrmx(2) = density to be incorporated in MT correction
Ci   rmax  :radius at which potential is evaluated, in a.u.
Ci   dmad  :Madelung matrix
Ci   bas   :basis vectors, in units of alat
Ci   plat  :primitive lattice vectors, in units of alat
Ci   vrl   :true => linear potential drop across plat(3)
Ci   noves :if true, evaluate vmtz only (using input ves)
Ci   vrmax :non-estat (XC) potential at wsr, for calc. vmtz
Ci   nangl :number of DLM classes.  nangl=0 unless DLM is active
Ci   idcc  :array pointing from DLM classes to parent classes
Ci   wts  :CPA weights for all classes
Ci   etrms:(CPA only) Madelung energy for each class is stored
Co Outputs
Co   ves:  Hartree potential on sphere surfaces
Co   emad: Madelung energy using convention as discussed in remarks
Co   trumad:True interatomic Madelung energy
Co   vmtz :muffin-tin zero (vmtz(2) = vmtz(up)-vmtz(dn))
Cr Remarks
Cr   Potential at site i, due to point charges from all points in the
Cr   lattice except i itself is given by (see madmat):
Cr
Cr   (1)  V_i(0) = 2 \sum_j M_ij DQ_j
Cr   DQ_j are the total charges and j runs over all sites in the
Cr   basis.  Rydberg units, with electronic charge positive,
Cr   nuclear charge negative.  The average Hartree potential at
Cr   the sphere radius rmax is
Cr
Cr   (2) V_i(rmax) = V_i(0)  +  2 DQ_i / rmax
Cr
Cr   where V_intra,i = 2 DQ_i / rmax is the intra-atomic contribution
Cr   to the potential.  The atomic code calculates potentials and total
Cr   energies subject to the boundary condition that the Hartree
Cr   potential at rmax is zero, regardless of the total charge inside
Cr   the sphere.  Then there is an onsite contribution to the
Cr   sphere energy,
Cr
Cr   (3) E_intra = 1/2 \int rho V_intra d^3r = 1/2 DQ_i V_intra,i
Cr
Cr   and the usual interatomic terms, \sum_i V_i(0) DQ_i, so that
Cr   the net Madelung contribution to the total energy is
Cr
Cr   (4) E_mad = 1/2 \sum_i V_i(rmax) DQ_i
Cr
Cr   Alternatively, the charges outside sphere i shift the potential
Cr   inside sphere i by a constant V_i(rmax).  This contributes
Cr   an extra term to the total energy Eqn 4, and also the
Cr   potential parameters for sphere i are shifted by V_i(rmax).
Cr
Cr   ---- Addition of background charge to sphere charges ---
Cr
Cr   ? The sphere charges consist of the electronic-nuclear charge qt,
Cr   less the smeared-out nuclear charge rhrmx(1)*vol
Cr
Cr   ---- MUFFIN-TIN Correction and homogeneous background ---
Cr
Cr   If the net system charge is a sum of nuclear charges, sphere
Cr   charges and a constant background charge of density n, the
Cr   charge DQ remains the net charge inside each sphere.  The
Cr   Madelung potential is calculated with charges Ztwid_i =
Cr   -DQ_i - Qtwid, ie the total charge in each sphere excluding
Cr   the charge Qtwid in sphere i from the homogeneous
Cr   background.  Eqn (2) still describes the potential at rmax.
Cr   Then for each R the Hartree energy of (Ztwid + homogeneous
Cr   background) is subtracted, since we use the Hartree energy
Cr   of the true density there. This term is
Cr
Cr   (5)  -3 Ztwid Qtwid / R  +  6 Qtwid^2 / 5 R
Cr
Cu Updates
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
Cu   25 Apr 12 (K Belashchenko) additions for CPA
Cu   21 Mar 03 New vconst
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical noves
      integer nbas,ic1,ic2,nangl
      integer ipc(nbas),nrc(*),idcc(*)
      double precision vconst,vrl,emad,trumad,qt(*),vrmax(2,*),rmax(*),bas(3,nbas),
     .  plat(3,3),dmad(nbas,nbas),ves(*),vmtz(2),dclabl(*),rhrmx(2),wts(*),etrms(22,*)
C ... Local parameters
      character clabl*8
      integer i,ibas,jbas,ic,icp,jc,ifi,ix,scrwid
      double precision vadd,sumdq,fpi,wt
      double precision dq,qih,qtw,ebak,ezv,sezv,sebak,sqb,avves,rmsv(2)
      parameter (scrwid=100)
C     Parameters for non-equilibrium mode
      double precision qlat(3,3),vol,b3(3),Lz,wk(nbas),xx,scrpar(3,ic2),spar
      procedure(integer) :: fxst,fopn,iprint,lgunit,awrite,isw
      procedure(real(8)) :: dlength

C ... Net, inhomogeneous, homo, background electronic & nuclear charges
      dq(ix)  = qt(ix) - fpi/3*rmax(ix)**3*rhrmx(2)
      qih(ix) = qt(ix) - fpi/3*rmax(ix)**3*rhrmx(1)
      qtw(ix) = fpi/3*rmax(ix)**3*(rhrmx(1)-rhrmx(2))
C     qbk(ix) = fpi/3*rmax(ix)**3*rhrmx(1)
C     zbk(ix) = fpi/3*rmax(ix)**3*rhrmx(2)

      fpi = 16*datan(1d0)

C ... Setup for linear potential drop across plat(3)
      if (abs(vrl) > 1d-7) then
        call dinv33(plat,1,qlat,vol)
        do  i = 1, 3
          b3(i) = qlat(i,3)/dlength(3,qlat(1,3),1)
        enddo
        Lz = b3(1)*plat(1,3) + b3(2)*plat(2,3) + b3(3)*plat(3,3)
      endif

      if (noves) then
        sumdq = 0
        if (iprint() >= 30) then
          do  i = 1, 2
            write (lgunit(i),3)
    3       format(/' MADPOT: potential not calculated.  Input Q,V:'/
     .      ' Class        Qtot       Qbak       V0        V0+Vxc')
          enddo
        endif
        do  ic = ic1, ic2
          sumdq = sumdq + nrc(ic)*dq(ic)
          if (iprint() >= 30) call r8tos8(dclabl(ic),clabl)
          do  i = 1, 2
            if (iprint() >= 30) write (lgunit(i),8) clabl,dq(ic),
     .          qtw(ic),ves(ic),ves(ic)+vrmax(1,ic)
          enddo
        enddo
          emad = 0
          trumad = 0
        goto 10
      endif

      avves = 0
      do  ibas = 1, nbas
        ic = ipc(ibas)
        if (ic < ic1 .or. ic > ic2) cycle
        ves(ic) = vconst
        if (abs(vrl) > 1d-7) then
          xx = (bas(1,ibas)*b3(1)+bas(2,ibas)*b3(2)+bas(3,ibas)*b3(3))/Lz*vrl
          ves(ic) = ves(ic) + xx
          wk(ic)  = xx
        endif
        do  jbas = 1, nbas
          jc = ipc(jbas)
          if (jc < ic1 .or. jc > ic2) cycle
          ves(ic) = ves(ic) + 2*qih(jc)*dmad(ibas,jbas)
        enddo
        avves = avves + (ves(ic) + 2*dq(ic)/rmax(ic))/nbas
      enddo
      if (iprint() >= 30) then
        do  i = 1, 2
          if (abs(vrl) <= 1d-7) then
            write (lgunit(i),4)
    4       format(/' Class',8x,'Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)')
          else
            write (lgunit(i),5)
    5       format(/' Class',8x,'Qtot       Qbak      vrl*z/L     Vmad     Vh(Rmax)    V(Rmax)')
          endif
        enddo
      endif
      emad   = 0
      trumad = 0
      sumdq  = 0
      sebak  = 0
      sezv   = 0
      sqb    = 0
      do  ic = ic1, ic2
        if (nrc(ic) == 0) cycle
        if (iprint() >= 30) call r8tos8(dclabl(ic),clabl)
        wt = 1d0
        if (nangl > 0) wt = wts(ic)
C       print *,'class',ic,'weight',wt
        if (wt /= 1d0 .and. wt /= 0d0) then
          call rx1('MADPOT: weight not 0 or 1 for class %i',ic)
        endif
        vadd = 2*dq(ic)/rmax(ic)*wt
        ebak = (3*qih(ic)*qtw(ic) + 6*qtw(ic)**2/5)/rmax(ic)
        ezv = 0.5d0*qih(ic)*ves(ic)
        if (iprint() >= 30 .and. ebak == 0) then
          do  i = 1, 2
          if (abs(vrl) <= 1d-7) then
C                write (lgunit(i),8) clabl,dq(ic),qtw(ic),ves(ic),ves(ic)
C     .                              +vadd,ves(ic)+vadd+vrmax(1,ic)
                call awrit5(' '//clabl//'%;11,6D%;11,6D%;11,6D%;11,6D%;11,6D',' ',scrwid,lgunit(i),
     .            dq(ic),qtw(ic),ves(ic),ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic))
            else
C                write (lgunit(i),9) clabl,dq(ic),qtw(ic),wk(ic),ves(ic),
C     .                              ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic)
                call awrit6(' '//clabl//'%;11,6D%;11,6D%;11,6D%;11,6D%;11,6D%;11,6D',' ',scrwid,
     .            lgunit(i),dq(ic),qtw(ic),wk(ic),ves(ic),ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic))
            endif
          enddo
        elseif (iprint() >= 30 .and. ebak /= 0) then
          do  i = 1, 2
            if (abs(vrl) <= 1d-7) then
                write (lgunit(i),8) clabl,dq(ic),qtw(ic),ves(ic),ves(ic)
     .                              +vadd,ves(ic)+vadd+vrmax(1,ic),
     .                              ezv-ebak
            else
                write (lgunit(i),9) clabl,dq(ic),qtw(ic),wk(ic),ves(ic),
     .                              ves(ic)+vadd,ves(ic)
     .                              +vadd+vrmax(1,ic),ezv-ebak
            endif
          enddo
        endif

        ebak = nrc(ic)*ebak
        ezv =  nrc(ic)*ezv
        sumdq = sumdq + nrc(ic)*dq(ic)
        trumad = trumad + ezv - ebak
        emad  = emad + ezv - ebak + .5d0*dq(ic)*vadd*nrc(ic)
        sezv  = sezv + ezv
        sebak = sebak + ebak
        sqb   = sqb + qtw(ic)*nrc(ic)
        ves(ic) = ves(ic) + vadd
        if (rhrmx(1) /= 0) ves(ic) = ves(ic) - avves
      enddo
C --- Assign ves to DLM classes (correct for different charge)
C ... This will be problematic if ic2 /= nclass
      if (nangl > 0) then
C   ... Read the CPA screening parameters if supplied in the scrpar file
        if (fxst('scrpar') == 1) then
          ifi = fopn('scrpar') ; rewind ifi
C         call iomagf(2,ic2,1,scrpar,scrpar,1,ifi)
          call ioextf(2,'scrpar',ic2,1,scrpar,scrpar,3,1,ifi)
          call fclr('scrpar',ifi)
        else
          scrpar = 0
        endif

        do  ic = ic2+1, ic2+nangl
          if (iprint() >= 30) call r8tos8(dclabl(ic),clabl)
          etrms(13,ic) = dq(ic)
          icp = idcc(ic)
          ves(ic) = ves(icp)
          vadd = 2*dq(ic)/rmax(ic)
          etrms(20,ic) = dq(ic)*(ves(ic) + .5d0*vadd)
C     ... Apply the screening correction for CPA
          spar = scrpar(1,icp)
          if (spar == 0) then
            continue
          elseif (spar < 1d0) then
            call rx1('Unphysical screening parameter for class',icp)
          else
            if (iprint() >= 30) then
              do  i = 1, 2
                write(lgunit(i),901) clabl,spar
              enddo
            endif
            vadd = vadd - 2*(dq(ic)-dq(icp))/spar
            etrms(20,ic) = etrms(20,ic) - 2*dq(ic)*(dq(ic)-dq(icp))/spar
          endif

          emad = emad + .5d0*dq(ic)*vadd*nrc(icp)*wts(ic)
          if (iprint() >= 30) then
            do  i = 1, 2
              write(lgunit(i),8) clabl,dq(ic),qtw(ic),ves(ic),
     .          ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic)
            enddo
          endif
          ves(ic) = ves(ic) + vadd
        enddo
      endif

C --- Reset Madelung potential if total system charge > .5 ---
      if (dabs(sumdq) > .01d0 .and. iprint() >= 10) print
     .    '('' (warning)  system has net charge: q='',f10.6)', sumdq
C      if (dabs(sumdq) > .5d0) then
C        emad = 0
C        if (iprint() >= 10) print
C     .    '('' madpot:  net charge too big; reset ves to zero'')'
C        call dpzero(ves(ic1),ic2-ic1+1)
C      endif

C --- Make vmtz by averaging over surfaces of spheres ---
   10 continue
      vmtz(2) = 0
      if (vmtz(1) == 0) call pvmadp(ic1,ic2,nrc,1,vrmax,rmax,ves,vmtz,rmsv)

C --- Printout ---
      if (iprint() >= 20 .and. .not. noves) then
        do  i = 1, 2
C   22   call awrit6(' Sum Q=%,6d  Emad=%,6d(%,6d)  '//
C     .      'Vmtz=%,6d%?;n; (up-dn=%,6d);;',' ',scrwid,lgunit(i),
C     .      sumdq, emad, trumad,vmtz,isw(vmtz(2) /= 0),vmtz(2))
        ix = awrite(' Sum Q=%,6d  Emad=%,6d(%,6d)'//
     .    '%?#n# Vconst=%,6;6d#%j#  Vmtz=%,6d%?;n; (up-dn=%,6d);;',' ',
     .    scrwid,lgunit(i),sumdq, emad, trumad,
     .    isw(vconst /= 0),vconst,
     .    vmtz,isw(vmtz(2) /= 0),vmtz(2))
      enddo

C  440   format(' Sum Q=',f10.6,'  Emad=',f10.6,'(',f10.6,')'/
C     .         ' Vmtz=',f9.6,' (up-dn)=',f9.6)
        if (rhrmx(1) /= 0) then
          print 6,sqb,sezv,sebak,sezv-sebak
    6     format(' SUMQB=',f9.6,'  EMAD=',f10.6,'(RB)',f10.6,'(SPHR)',
     .            '   DIFF=',f9.6)
          print 7,avves
    7     format(' Shifting average VES by',f10.6)
        endif

      endif
      if (iprint() >= 20 .and. noves) then
        do  i = 1, 2
          call awrit4(' Sum Q=%,6d  Vmtz=%,6d%?;n; (up-dn=%,6d);;',' ',
     .      scrwid,lgunit(i),sumdq, vmtz,isw(vmtz(2) /= 0),vmtz(2))
        enddo
C        do  122  i = 1, 2
C  122   write(lgunit(i),540) sumdq, vmtz
C  540   format(' SUM Q=',f10.6,'  VMTZ=',f9.6,'   VMTZ(UP-DN)=',f9.6)
      endif
    8 format(1x,a,5F11.6:'  ezv-ebak=',f9.6)
    9 format(1x,a,6F11.6:'  ezv-ebak=',f9.6)
  901 format(1x,a,' CPA screening correction with R=',F6.3)

      end
      subroutine pvmadp(ic1,ic2,nrc,lvrmx,vrmax,rmax,ves,vmtz,rmsv)
C- Make vmtz
C ----------------------------------------------------------------------
Ci Inputs
Ci   ic1   :starting class index
Ci   ic2   :ending class index
Ci   nrc   :number of members of this calss
Ci   lvrmx :1 if to add vrmax; 0 if to just average ves
Ci   rmax  :augmentation radius, in a.u.,
Ci   ves   :estat potential, by class
Co Outputs
Co   vmtz  :muffin-tin zero
Co   rmsv  :sum ves(ic)**2 wgt(ic)
Cu Updates
Cu   22 Mar 03 New argument rmsv
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lvrmx,ic1,ic2,nrc(ic2)
      double precision vrmax(2,ic2),rmax(ic2),ves(ic2),vmtz(2),rmsv(2)
C ... Local parameters
      integer ic
      double precision swgt,wgt

      vmtz(1) = 0d0
      vmtz(2) = 0d0
      swgt = 0d0
      rmsv(1) = 0
      rmsv(2) = 0
      do  ic = ic1, ic2
        wgt = nrc(ic)*rmax(ic)**3
        if (lvrmx == 1) then
          vmtz(1) = vmtz(1) + wgt*(vrmax(1,ic) + ves(ic))
          vmtz(2) = vmtz(2) + wgt*vrmax(2,ic)
          rmsv(1) = rmsv(1) + wgt*ves(ic)**2
          rmsv(2) = rmsv(2) + wgt*(vrmax(1,ic) + ves(ic))**2
        else
          vmtz(1) = vmtz(1) + wgt*ves(ic)
          rmsv(1) = rmsv(1) + wgt*ves(ic)**2
        endif
        swgt = swgt + wgt
      enddo
      vmtz(1) = vmtz(1)/swgt
      vmtz(2) = vmtz(2)/swgt
      rmsv(1) = rmsv(1)/swgt
      rmsv(2) = rmsv(2)/swgt

      end

      subroutine dplmom(nbas,bas,alat,plat,ipc,qt,qbar,vh,delq,delv)
C- Prints out contribution to dipole moment in plat(3) from excess charge
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   qt    :electronic charge, less nuclear charge
Ci         :(excluding homogeneous nuclear charge in rhrmx(2))
Ci   qbar  :Average charge (to remove monopole contribution)
Co Outputs
Co   h     :height of atom, defined as (qlat(3).bas)/|qlat(3)|
Co   vh    :resolution by site to dipole contribution
Co   delq  :net charge
Co   delv  :dipole moment
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   15 Sep 04 (S Faleev) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,ipc(nbas)
      double precision bas(3,nbas),alat,plat(3,3),qt(*),qbar,vh(nbas)
C ... Local parameters
      integer ib,i,ic
      double precision qlat(3,3),vol,b3(3),pi,z,Lz,dplm,delq,delv,dq
C     procedure(integer) :: nglob,iprint
      procedure(real(8)) :: dlength

      pi = 4*datan(1d0)
      call dinv33(plat,1,qlat,vol)
      do  i = 1, 3
        b3(i) = qlat(i,3)/dlength(3,qlat(1,3),1)
      enddo
      Lz = b3(1)*plat(1,3) + b3(2)*plat(2,3) + b3(3)*plat(3,3)

      dplm = 0d0; delq = 0
      do  ib = 1, nbas
        ic = ipc(ib)
        if (ic <= 0) cycle
        z = b3(1)*bas(1,ib) + b3(2)*bas(2,ib) + b3(3)*bas(3,ib)
        dq = qt(ic)-qbar
        delq = delq + dq
        dplm = dplm + 2*z*dq
        vh(ib) = 2*z*dq * 4*pi*Lz/vol/alat
      enddo
      delv = 4*pi*dplm*Lz/vol/alat

C     print *, sum(vh), delv

      call info2(20,0,0,' Dipole in plat(3) from charge distribution = %d'//
     .  '  Net charge = %d',delv,delq)
      end

C      subroutine uefld(nbas,bas,alat,plat,dmad)
CC- Add contribution of uniform electric field E=4*pi/vol*sum_i Q_i*z_i
CC- to Madelung constants due dipole moment of supercell
CC ----------------------------------------------------------------------
CCi Inputs
CCi Inputs
CCi   nbas  :size of basis
CCi   bas   :basis vectors, in units of alat
CCi   alat  :length scale of lattice and basis vectors, a.u.
CCi   plat  :primitive lattice vectors, in units of alat
CCo Input/Outputs
CCo   dmad  :constant E-field contribution added to Madelung matrix
CCu Updates
CCu   9 Feb 04 (S.Faleev) First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nbas
C      double precision alat,plat(3,3),bas(3,nbas),dmad(nbas,nbas)
CC ... Local parameters
C      integer i,j,ibas,jbas
C      double precision pi,b3(3),s(3),qlat(3,3),vol,z1,z2
C
C      pi = 4*datan(1d0)
C      call dinv33(plat,1,qlat,vol)
C      vol = dabs(vol)
C      do i = 1,3
C        b3(i) = qlat(i,3)/sqrt(qlat(1,3)**2+qlat(2,3)**2+qlat(3,3)**2)
C      enddo
C
CC ... Checks
C      call dpzero(s,3)
C      do j = 1,3
C        do i = 1,3
C          s(j) = s(j) + b3(i)*plat(i,j)
C        enddo
C      enddo
C      if ( (abs(s(1))+abs(s(2))) > 1d-10 .or. s(3) <= 0d0 )
C     .  call rx('uefld: wrong qlat(*,3)')
C
C      do  jbas = 1,nbas
C        do  ibas = 1, nbas
C          z1 = 0d0
C          z2 = 0d0
C          do  i = 1, 3
C            z1 = z1 + bas(i,ibas)*b3(i)
C            z2 = z2 + bas(i,jbas)*b3(i)
C          enddo
C          dmad(ibas,jbas) = dmad(ibas,jbas) + 4d0*pi*z1*z2/(vol*alat)
C        enddo
C      enddo
C
CCCCCCCCCCCCCCCC TEMP CCCCCCCCCCCCC
Cc      write(*,'(a,i4,4f12.6)')'uefld:nbas,b3,vol=',nbas,b3,vol
Cc      write(*,'(a,9f8.4)')    'uefld:plat=',plat
Cc      do i = 1,nbas
Cc      write(*,'(a,i4,3f12.6)')'uefld: ibas,bas=',i,bas(1,i),bas(2,i),
Cc     .        bas(3,i)
Cc      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      end

C      subroutine snot(strn,nbas,fc,ipc)
C      implicit none
C      character *(*) strn
C      integer nbas,ipc(nbas)
C      double precision fc(*)
C      integer ib
C      double precision fb(nbas)
C
C      do  ib = 1, nbas
C        fb(ib) = fc(ipc(ib))
C      enddo
C      call prmx(strn,fb,nbas,nbas,1)
C      end
