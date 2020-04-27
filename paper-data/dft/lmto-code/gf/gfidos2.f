      subroutine gfidos2(nl,nbas,moddos,lsw,ib1,ib2,indxsh,lidim1,lidim2,
     .  hdim1,hdim2,ldg,ldh,ghh,wtkp,zp0,wz,isp,nsp,nspc,pp,vRLshf,
     .  vshft,ipc,nrclas,iz,nz,nzp,ipl,npl,kpl,offpd,offgd,dosi,pdos,
     .  pdosl,gd,qnu,qne,dosne,orbtm,dmat,rhos,ncomp,theta,dlmcl,wts,gc,
     .  gcorr,mxy,qcorr,norb,sop)
C- Accumulate DOS and moments DOS for one CPA site, energy and kp
C-----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   moddos:specifies what kinds of quantities to calculate
Ci          1s digit
Ci           1 compute dosi (integrated dos and related quantities)
Ci           2 same as 1, but also set dos(1) to zero before starting
Ci         10s digit
Ci           1 compute partial dos, (see pdos, Outputs)
Ci           2 accumulate diagonal part of g, (see gd, Outputs)
Ci             (used eg for Pade approximant to fix Fermi level)
Ci           3 both 1 and 2
Ci          4s bit if set, pdos decomposed into both l and m
Ci             By default, pdos is contracted over m
Ci        100s digit
Ci           1 compute ASA moments of the DOS
Ci           2 compute Lloyd corrections
Ci           3 both 1 and 2
Ci           4 orbital moment (Add 4 to 100s digit)
Ci       1000s digit
Ci           1 Accumulate the site-diagonal density matrix, array dmat
Ci           2 Accumulate the full density matrix, array dmat
Ci           4 compute spin density matrix, array rhos
Ci      10000s digit
Ci           2 only process CPA sites with ncomp>1 (required option)
Ci     100000s digit (not implemented)
Ci           1 compute qnum, which is the m-resolved qnu
Ci   lsw   :switches
Ci         :1s digit
Ci         :0 real harmonics, ordered m = -l:l
Ci         :1 spherical harmonics.  Call lmorder to get m ordering.
Ci         :10s digit (not implemented)
Ci         :1 initialize to zero those the following quantities
Ci         :  that are to be accumulated: dosi,dmat,rhos
Ci   ib1   :accumulate information for block ib1 ... ib2
Ci   ib2   :accumulate information for block ib1 ... ib2
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   lidim1:offset to subtract from the hamiltian index table indxsh,
Ci         :for the gll supplied.  (Only a part of g is given in the layer GF case)
Ci         :Not used here because gc is used in place of gll.
Ci   lidim2
Ci   hdim1 :Not used now
Ci   hdim2
Ci   ldg
Ci   ldh
Ci   ghh   :Not used now
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci          If k-sum already performed, wtkp should be 1.
Ci   zp0   :complex energy
Ci   wz    :weight for complex energy integration
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   pp    :potential parameters (atomsr.f) needed to make energy
Ci          moments relative to enu.
Ci   vRLshf:array of orbital- and spin-dependent potential shifts
Ci   vshft :array of site potential shifts
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nrclas:number of atoms in the icomp class
Ci   iz    :index to quantities that keep energy-dependent
Ci          information (pdos,gd).
Ci   nz    :leading dimension of quantities that keep energy-dependent
Ci          information (pdos,gd). (for non-equilibrium mode: nz>nzp)
Ci   nzp   :number of energy points on equilibrium contour
Ci         :points iz>nzp are for nonequilibrium contour
Ci   ipl   :index to current principal layer (nonequil mode only)
Ci   npl   :total number of principal layers (nonequil mode only)
Ci   kpl   :index specifing what channel dos should be accumulated
Ci   orbtm
Ci   ncomp :number of CPA components for this site
Ci   theta :angles for the CPA components (DLM site or future noncollinear)
Ci   dlmcl :array pointing to the first CPA class for given CPA site
Ci   wts   :CPA weights
Ci   gc    :conditionally-averaged GF for one site (to be called with ib1=ib2)
Ci   gcorr :Lloyd correction for single-electron energy
Ci   norb  :number of orbitals for gc
Ci   sop   :spin-orbit parameters (atomsr.f).  Not used now
Ci Inputs/Outputs
Cio  offpd :On input offset to pdos array corresponding to site ib1
Cio        :On output, offset to pdos for site ib2+1
Cio  offgd :On input, offset to gd array corresponding to site ib1
Cio        :On output, offset to gd array for site ib2+1
Co Outputs
Co   dosi  :integrated dos and moments of dos for sites ib1..ib2
Co          1: dos at this e : -wtkp/pi Im G
Co          2: nos at this e : -wtkp/pi Im (wz*G)
Co          3: 1st energy mom : -wtkp/pi Im (wz*G z)
Co          4: 2nd energy mom : -wtkp/pi Im (wz*G z**2)
Co          5: projection of dos to fermi level (not calculated here)
Co   pdos  :energy resolved partial dos averaged over CPA components
Co   pdosl :energy resolved partial dos resolved by CPA components
Co   gd    :RL- diagonal Green's function added from diag. of gll
Co   qnu   :ASA energy moments
Co   qcorr :Lloyd corrections for single-electron energy
Co   mxy   :off-diagonal local moments for one site
Co   dmat  :density-matrix for basis site ib
Co   qne   :addition to charge by integration along nonequilib. contour
Co   dosne :dos at first and last points of nonequilibrium contour
Co         :(calc. with G^< for ipl<npl and with G^r for ipl=npl)
Co   rhos  :on-site spin density-matrix, contracted over m (untested)
Cl Local variables
Cl   offgi0:offset to first orbital of site ib1, in indxsh table
Cl   lmRL  :entry in indxsh for current orbital
Cl   li    :current orbital in gll, downfolding order, relative to ib1
Cl   lmgd  :current orbital in gd array, rel to 1st orbital of this site
Cl   lpd   :current channel pdos array, rel to 1st dos channel this site
Cl   ljh   :current row dimension for hh block
Cl   jsp   :index to spin-diagonal part of gll
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   ksp   :isp in the collinear case, and jsp in the noncollinear
Cl         :Thus, for the spinor parts of gll and gd:
Cl         :(jsp,jsp) = indices to diag part of gll for current spin
Cl         :(ksp,jsp) = indices to diag part of gd  for current spin
Cl   oghi  :offset to imaginary part of ghh
Cl   lhh   :true if current orbital belongs to higher downfolding block.
Cl   greal :Re(diagonal g) for current orbital, from gll or ghh
Cl   gimag :Im(diagonal g) for current orbital, from gll or ghh
Cl   gr12  :(nspc=2 only) like greal, but for off-diagonal in spins
Cl   gi12  :(nspc=2 only) like gimag, but for off-diagonal in spins
Cu Updates
Cu  18 Jun 18 Synchronize with updated spherical harmonics
Cu  01 Nov 11 (Belashchenko) Cloned from gfidos, processes one CPA site
Cu            Unused remnants of gfidos need to be cleaned up!
Cu  16 Nov 07 Added orbital-dependent potential shifts (vRLshf)
Cu  10 Jul 04 (S.Faleev) changes to handle noneqilibrium mode
Cu            Argument list changed.
Cu  18 Jun 04 (A Chantis) orbital moments in fully relativstic case
Cu   9 Mar 00 redesigned the iwaves downfolding, added h-waves
Cu            Argument list changed.
Cb Bugs
Cb  When both downfolding and nlo=nl*nl, the following line is wrong:
Cb     if (.not. lall) lpd = l
Cb  Site-diagonal density-matrix for iwaves is incorrectly generated.
Cb  This routine is supposed to process CPA sites only, so moddlm is
Cb  always supposed to be 2. Also, ib1=ib2 assumed. Remnants of
Cb  standard gfidos can be cleaned up.
C-----------------------------------------------------------------------
      implicit none
C Passed variables
      integer nl,nbas,ldg,kpl,ib1,ib2,isp,nsp,nspc,lsw,iz,nz,offpd,
     .  offgd,ipc(ib2),nrclas(*),indxsh(ldg),ldh,ipl,npl,nzp
      integer lidim1,lidim2,hdim1,hdim2
      integer ncomp,dlmcl,norb
      double precision ghh(ldh,nspc,*),gd(2,nz,nsp,nspc,0:*),zp0(2),
     .  wtkp,zp(2),wz(2),pdos(nz,nsp,0:1),dosi(5,nsp,kpl),vshft(nbas),
     .  pp(6,nl,nsp,*),qnu(3,nl,nsp,*),rhos(2,3,nl,2,2,1),
     .  orbtm(nl,nspc,nbas),qne(nsp,kpl),dosne(2,nsp,kpl),
     .  vRLshf(lidim2+hdim2,nsp),qcorr(2,nl,nsp,*)
C     double precision qnum(3,nl*nl,nsp,*)
      double complex gc(norb,2,norb,2,ncomp),gcorr(norb,ncomp,2),
     .  dmat(norb,norb,nsp,nspc,4,ncomp)
      double precision wts(ncomp),theta(ncomp,2),mxy(nl,2,*)
      real(8) sop(0:nl-1,nsp,nsp,9,*)
C Local variables
      logical lall,lpdos,lpade,ldosi,lmom,lgcorr,lhh,lxy,lorbm,lmomm
      integer i,ib,ic,im,is1,is2,jsp,ksp,l,landm(norb,2),ldmat,lhdim,
     .  li,ljh,lmRL,lmgd,lpd,m,moddos,mxorb,offgd0,offgi0,offpd0,oghi
      double precision pi,w(2),wt,zmenu(2),xx(2,0:2),yy(2,0:1),bync,
     .  greal,gimag,grc,gic,gr12,gi12
      double complex wei,ci,gcwk(norb,norb),xw
      double complex zmenuc(norb,2)
C For CPA
      integer moddlm,ic1,ic2,lmRLd,lmgdd,offpdl,ic0,icomp,ilm,n1,n2
      double precision wtd1,wtd2,pdosl(nz,nsp,0:1)
C     double precision qq,qqt                             !TEMP

      procedure(integer) :: ll,nglob

C ... Setup : decode what is to be calculated
      ldosi = mod(moddos,10) == 1 .or. mod(moddos,10) == 2
      lpdos = mod(mod(moddos/10,10),2) == 1
      lpade = mod(mod(moddos/10,10),4) > 1
      lall  = mod(moddos/10,10) >= 4
      lmom  = mod(mod(moddos/100,10),2) == 1
      lorbm = mod(mod(moddos/100,10)/4,2) == 1
      lgcorr  = mod(mod(moddos/100,10)/2,2) == 1
      ldmat = mod(moddos/1000,10)
      moddlm = mod(moddos/10000,10)
      lmomm = mod(mod(moddos/100000,10),2) == 1
      mxorb = nglob('mxorb')
      i = 1; if (mod(lsw,10) == 1) call lmorder(0,i,landm,landm); if (i==0) i = 2
      call lmorder(i,ll(norb),landm,landm(1,2))

      oghi  = ldh*nspc*mxorb
      call sanrg(.true.,iz,1,nz,'gfidos','iz')

c     if (moddlm /= 2) call rx1('GFIDOS2: moddlm = %i',moddlm)
      if (ib1 /= ib2) call rx('CPA and ib1!=ib2')
c     if (ldmat /= 0) then
c       print *,' WARNING: gfidos2 not ready for ldmat'
c     endif

      i = 1
      ci = dcmplx(0,1)
      pi = 4*datan(1d0)
C ... wt is -1/pi * kpt weight
      wt = -wtkp/pi
      if (nsp == 1) wt = wt*2
C ... want to keep passed energy zp0(1:2) = zp(1:2;izp) unchanged
      zp(1) = zp0(1)
      zp(2) = zp0(2)
C ... Nonequilibrium mode
      if (iz > nzp) then
         if (abs(wz(2)) > 1d-10) call rx('gfidos: wz(2) ' //
     .        'nonzero for non-equilibrium energy point')
         wz(2) = 0d0
C      G^< constructed only for ipl<npl; for ipl=npl it is G^r.
C      Integration along nonequilibrium contour with G^<
C      should be performed with real energy in integrand
         if (ipl < npl) zp(2) = 0d0
      endif

C ... w is -1/pi * kpt weight * energy weight
      w(1) = wz(1)*wt
      w(2) = wz(2)*wt
      if (mod(moddos,10) == 2 .and. iz <= nzp) then
        dosi(1,isp,kpl) = 0
        dosi(1,max(isp,nspc),kpl) = 0
      endif

C     Keep a local copy of offpd,offgd in case loop over second spin
      offpd0 = offpd
      offgd0 = offgd
C     Offset to indxsh
      offgi0 = mxorb*(ib1-1)

      lhdim = lidim2+hdim2

C --- For each of of the coupled spins, do ---
C     pass1 = .true.

C      qq = 0d0                      !TEMP
C      qqt = 0d0                     !TEMP

c     if (nspc == 2) print *,'WARNING(gfidos2): CPA and coupled spins'

      do  10  jsp = 1, nspc
      ksp = max(jsp,isp)

      if (lpdos) offpd = offpd0
      if (lpdos) offpdl = 0
      if (lpade) offgd = offgd0

C --- For each site in ib1..ib2, do ---
      lmRL = offgi0
      do  11  ib = ib1, ib2
        lmgd = -1
        if (moddlm == 2 .and. ncomp < 2) then
          do  l = 0, nl-1
            do  m = -l, l
              lmRL = lmRL + 1
              if (indxsh(lmRL) > lhdim) cycle
              lmgd = lmgd + 1
              if (lall) then
                lpd  = lmgd
              else
                lpd = l
              endif
            enddo
          enddo
          goto 12
        endif
        ic0 = ipc(ib)
        bync = 1d0 / nrclas(ic0)
        ljh = 0
C   ... ib is a CPA site, loop over angles
        ic1 = dlmcl
        ic2 = dlmcl + ncomp - 1
        lmRLd = lmRL
        lmgdd = lmgd
        do  ic = ic1, ic2
        lmRL = lmRLd
        lmgd = lmgdd
        icomp = ic - ic1 + 1
        lxy = (abs(dsin(theta(icomp,1))) > 1d-8) .and. lmom
        lxy = lxy .or. (nspc == 2)
        wtd1 = (1 + dcos(theta(icomp,1)))*wts(icomp)/2
        wtd2 = (1 - dcos(theta(icomp,1)))*wts(icomp)/2
        ilm = 0
        do  l = 0, nl-1
        do  m = -l, l
          ilm = ilm + 1
          im = landm(ilm,2)
C         lm = ilm - 1
C         Skip orbitals outside l,i,h blocks
          lmRL = lmRL + 1
          if (indxsh(lmRL) > lhdim) goto 23

C         Determine whether ll, ii, or hh  applies
          lhh = indxsh(lmRL) > lidim2

          if (lhh) then
            call rx('h-waves not allowed on CPA sites')
          else
            li = indxsh(lmRL) - lidim1
            greal = dble(gc(ilm,ksp,ilm,ksp,icomp))
            gimag = dimag(gc(ilm,ksp,ilm,ksp,icomp))
            grc   = dble(gcorr(ilm,icomp,ksp))
            gic   = dimag(gcorr(ilm,icomp,ksp))
            gr12  = 0
            gi12  = 0
            if (lxy .and. ksp == 1) then
              xw = gc(ilm,1,ilm,2,icomp) + gc(ilm,2,ilm,1,icomp)
              xx(1,0)  = w(1)*dimag(xw) + w(2)*dble(xw)
              xw = (gc(ilm,1,ilm,2,icomp) - gc(ilm,2,ilm,1,icomp))*ci
              xx(2,0)  = w(1)*dimag(xw) + w(2)*dble(xw)
              mxy(l+1,1,ic) = mxy(l+1,1,ic) + xx(1,0) * bync
              mxy(l+1,2,ic) = mxy(l+1,2,ic) + xx(2,0) * bync
            endif
          endif

C         Index to gd --- does not care which block g comes from
          lmgd = lmgd + 1

C         l- or lm- index to lpdos
          if (lall) then
            lpd  = lmgd
          else
            lpd = l
          endif

C     ... xx is -1/pi * gf * energy contour weight * kpt weight
          xx(1,0) = w(1)*greal - w(2)*gimag
          xx(2,0) = w(1)*gimag + w(2)*greal
          yy(1,0) = w(1)*grc - w(2)*gic
          yy(2,0) = w(1)*gic + w(2)*grc

C     ... Partial DOS
          if (lpdos) then
            pdos(iz,ksp,offpd+lpd) = pdos(iz,ksp,offpd+lpd)
     .        + wt*gimag*wtd1
            pdos(iz,3-ksp,offpd+lpd) = pdos(iz,3-ksp,offpd+lpd)
     .        + wt*gimag*wtd2
            pdosl(iz,ksp,offpdl+lpd) = pdosl(iz,ksp,offpdl+lpd)
     .        + wt*gimag
          endif

C     ... Save gd = diagonal part of g for this energy, orbital
          if (lpade) then
            gd(1,iz,ksp,jsp,offgd+lmgd) = gd(1,iz,ksp,jsp,offgd+lmgd) +
     .        wtkp*greal*wtd1
            gd(2,iz,ksp,jsp,offgd+lmgd) = gd(2,iz,ksp,jsp,offgd+lmgd) +
     .        wtkp*gimag*wtd1
            gd(1,iz,3-ksp,min(3-ksp,nspc),offgd+lmgd) =
     .        gd(1,iz,3-ksp,min(3-ksp,nspc),offgd+lmgd) +wtkp*greal*wtd2
            gd(2,iz,3-ksp,min(3-ksp,nspc),offgd+lmgd) =
     .        gd(2,iz,3-ksp,min(3-ksp,nspc),offgd+lmgd) +wtkp*gimag*wtd2
C           Retain the off-diagonal parts in spinor space
            if (nspc == 2 .and. lhh) then
              gd(1,iz,jsp,3-jsp,offgd+lmgd) = 0
              gd(2,iz,jsp,3-jsp,offgd+lmgd) = 0
            elseif (nspc == 2) then
c           gd(1,iz,jsp,3-jsp,offgd+lmgd)=gd(1,iz,jsp,3-jsp,offgd+lmgd)
c    .          + wtkp*gr12
c           gd(2,iz,jsp,3-jsp,offgd+lmgd)=gd(2,iz,jsp,3-jsp,offgd+lmgd)
c    .        + wtkp*gi12
            endif
          endif

C     ... PL-integrated dos and moments of dos
          if (ldosi) then
C            if (lsw/10 == 1 .and. pass1) then
C              call dpzero(dosi(1,ksp,kpl),5*nspc)
C            endif
            if (iz <= nzp) then
              dosi(1,ksp,kpl) = dosi(1,ksp,kpl) + wt*gimag*wtd1
              dosi(1,3-ksp,kpl) = dosi(1,3-ksp,kpl) + wt*gimag*wtd2
            else
C             DOS at endpoints of nonequilibrium contour
              if (nz <= nzp) call rx('gfidos: iz>nzp but nz<=nzp')
              if (iz == nzp+1)
     .           dosne(1,ksp,kpl) = dosne(1,ksp,kpl) + wt*gimag
              if (iz == nz)
     .           dosne(2,ksp,kpl) = dosne(2,ksp,kpl) + wt*gimag
C             if ( ipl > -1 .and. ipl < npl) then !TEMP
C                qq = qq + wt*gimag                     !TEMP
C                qqt = qqt+ wz(1)*wt*gimag              !TEMP
C             endif                                     !TEMP
            endif
C           if (isp == 1) print 333,iz,ib,lmgd,dosi(1,1,1)
C  333      format(3i4,f15.10)
C
C           Contribution to charge from nonequilibrium contour
            if (iz > nzp) qne(ksp,kpl) = qne(ksp,kpl) + xx(2,0)

            xx(1,1) = zp(1)*xx(1,0) - zp(2)*xx(2,0)
            xx(2,1) = zp(2)*xx(1,0) + zp(1)*xx(2,0)
            xx(1,2) = zp(1)*xx(1,1) - zp(2)*xx(2,1)
            xx(2,2) = zp(2)*xx(1,1) + zp(1)*xx(2,1)
C           No integration along noneqilibrium contour for left lead
C           Other energy integration includes both equilibrium and
C           nonequilibrium points
            if (ipl == -1 .and. iz > nzp) then
            else
            do  i = 0, 2
              dosi(2+i,ksp,kpl) = dosi(2+i,ksp,kpl) + xx(2,i)*wtd1
              dosi(2+i,3-ksp,kpl) = dosi(2+i,3-ksp,kpl) + xx(2,i)*wtd2
            enddo
            endif
          endif

C     ... Energy moments (z-enu)
          if (lmom .or. lgcorr .or. mod(ldmat,4) == 1) then
C            if (lsw/10 == 1 .and. ib == ib1) then
C              call dpzero(dosi(1,isp,kpl),5*nspc)
C            endif
            zmenu(1) = zp(1) - ! sop(l,ksp,ksp,1,ic) -
     .        (pp(1,l+1,ksp,ic)+vshft(ib)+vrlshf(indxsh(lmRL),ksp))
            zmenu(2) = zp(2)
            do  is1 = 1, 2
              zmenuc(ilm,is1) = dcmplx(zp(1),zp(2)) - ! sop(l,is1,is1,1,ic) -
     .          (pp(1,l+1,is1,ic)+vshft(ib)+vrlshf(indxsh(lmRL),is1))
            enddo
          endif

C     ... ASA energy moments qnu, possibly orbm and lgcorr
          if (lmom .or. lgcorr) then

            xx(1,1) = zmenu(1)*xx(1,0) - zmenu(2)*xx(2,0)
            xx(2,1) = zmenu(2)*xx(1,0) + zmenu(1)*xx(2,0)
            xx(1,2) = zmenu(1)*xx(1,1) - zmenu(2)*xx(2,1)
            xx(2,2) = zmenu(2)*xx(1,1) + zmenu(1)*xx(2,1)
            yy(2,1) = zmenu(2)*yy(1,0) + zmenu(1)*yy(2,0)

C           No integration along noneqilibrium contour for left lead
C           Other energy integration includes both equilibrium and
C           nonequilibrium points
            if (ipl == -1 .and. iz > nzp) then
            else
              if (lmom) then
                do  i = 1, 3
                  qnu(i,l+1,ksp,ic) = qnu(i,l+1,ksp,ic) + xx(2,i-1)*bync
                enddo
                if (nspc == 2 .and. lorbm) then
                  orbtm(l+1,ksp,ic) = orbtm(l+1,ksp,ic) + im*xx(2,0)*bync
                endif
              endif
C              if (lmomm) then
C                do  i = 1, 3
C                  qnum(i,ilm,ksp,ib) = qnum(i,ilm,ksp,ib) + xx(2,i-1)*bync
C                enddo
C              endif
              if (lgcorr) then
                do  i = 1, 2
                  qcorr(i,l+1,ksp,ic)=qcorr(i,l+1,ksp,ic)+yy(2,i-1)*bync
                enddo
              endif
            endif
          endif

C     ... Spin density matrix.  Should be hermitian in spin space,
C         but there is some deviation owing to integration errors.
          if (ldmat >= 4) then
C     ...   Temporary skip (not implemented yet)
            goto 23
            xx(1,0) = w(1)*gr12 - w(2)*gi12
            xx(2,0) = w(1)*gi12 + w(2)*gr12
            xx(1,1) = zmenu(1)*xx(1,0) - zmenu(2)*xx(2,0)
            xx(2,1) = zmenu(2)*xx(1,0) + zmenu(1)*xx(2,0)
            xx(1,2) = zmenu(1)*xx(1,1) - zmenu(2)*xx(2,1)
            xx(2,2) = zmenu(2)*xx(1,1) + zmenu(1)*xx(2,1)
            do  i = 1, 3
              rhos(1,i,l+1,jsp,jsp,ic)   = qnu(i,l+1,jsp,ic)
C             Im rhos = Re g;  Re rhos = -Im g
              rhos(1,i,l+1,jsp,3-jsp,ic) = rhos(1,i,l+1,jsp,3-jsp,ic) +
     .                                     xx(2,i-1)*bync
              rhos(2,i,l+1,jsp,3-jsp,ic) = rhos(2,i,l+1,jsp,3-jsp,ic) -
     .                                     xx(1,i-1)*bync
            enddo
          endif

   23     continue
C         pass1 = .false.
        enddo   ! Loop over m
        enddo   ! Loop over l
        if (lpdos) offpdl = offpdl+lpd+1

C   ... Site-diagonal density matrix
        if (mod(ldmat,4) == 1) then
          wei = dcmplx(w(1),w(2)) * (-ci)
          gcwk(:,:) = gc(:,ksp,:,ksp,icomp)
!           dmat(:,:,ksp,jsp,1,icomp) = dmat(:,:,ksp,jsp,1,icomp) + wei * gcwk(:,:)
          dmat(:,:,ksp,ksp,1,icomp) = dmat(:,:,ksp,ksp,1,icomp) + wei * gcwk(:,:)
          if (nspc == 1) cycle

          do n2 = 1, norb
          do n1 = 1, norb
            dmat(n1,n2,ksp,ksp,2,icomp) = dmat(n1,n2,ksp,ksp,2,icomp) + wei * zmenuc(n1,ksp) * gcwk(n1,n2)
            dmat(n1,n2,ksp,ksp,3,icomp) = dmat(n1,n2,ksp,ksp,3,icomp) + wei * gcwk(n1,n2) * zmenuc(n2,ksp)
            dmat(n1,n2,ksp,ksp,4,icomp) = dmat(n1,n2,ksp,ksp,4,icomp) + wei * zmenuc(n1,ksp) * gcwk(n1,n2) * zmenuc(n2,ksp)
          enddo
          enddo
          if (nspc == 2 .and. jsp == 1) then
            do is1 = 1, 2
              is2 = 3 - is1
              gcwk(:,:) = gc(:,is1,:,is2,icomp)
              do n2 = 1, norb
              do n1 = 1, norb
                dmat(n1,n2,is1,is2,1,icomp) = dmat(n1,n2,is1,is2,1,icomp) + wei                  * gcwk(n1,n2)
                dmat(n1,n2,is1,is2,2,icomp) = dmat(n1,n2,is1,is2,2,icomp) + wei * zmenuc(n1,is1) * gcwk(n1,n2)
                dmat(n1,n2,is1,is2,3,icomp) = dmat(n1,n2,is1,is2,3,icomp) + wei                  * gcwk(n1,n2) * zmenuc(n2,is2)
                dmat(n1,n2,is1,is2,4,icomp) = dmat(n1,n2,is1,is2,4,icomp) + wei * zmenuc(n1,is1) * gcwk(n1,n2) * zmenuc(n2,is2)
              enddo
              enddo
            enddo
          endif
        endif

        enddo   ! Loop over classes

   12   if (lpdos) offpd = offpd+lpd+1
        if (lpade) offgd = offgd+lmgd+1
   11 continue
   10 continue

C ... Entire density-matrix: dmat dimensioned dmat(ldg,nspc,ldg,nsp),
      if (mod(ldmat,4) == 2) then
c        call yscal((ldg*nspc)**2,w(2),-w(1),gll,gll(1,1,1,1,2),1)
        if (nspc == 1) then
c         call dpsadd(dmat,gll,ldg**2,1+ldg**2*(isp-1),1,1d0)
        else
c         call dpsadd(dmat,gll,8*ldg**2,1,1,1d0)
        endif
      endif

C      if (iz == nz .and. nspc == 2) then
C        do ib = ib1, ib2
C          ic = ipc(ib)
C          do l= 0, nl-1
C            do ksp = 1, 2
C              print*, ic,l,ksp,orbtm(l+1,ksp,ic)
C            enddo
C          enddo
C        enddo
C      endif

C     call zprm('gd',2,gd,nz*nsp*nspc,nz*nsp*nspc,offgd-offgd0)
C     call yprm('dosi',1,dosi,0,5*nsp,5*nsp,7)
      end
