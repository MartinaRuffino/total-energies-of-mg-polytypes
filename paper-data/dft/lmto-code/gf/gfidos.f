      subroutine gfidos(nl,nbas,moddos,lsw,ib1,ib2,indxsh,lidim1,lidim2,
     .  hdim1,hdim2,ldg,gll,ldh,ghh,wtkp,zp0,wz,isp,nsp,nspc,pp,vRLshf,
     .  vshft,ipc,nrclas,ncomp,nrhos,iz,nz,nzp,ipl,npl,kpl,offpd,offgd,dosi,
     .  pdos,gd,qnu,qne,dosne,orbtm,dmat,rhos)
C- Accumulate DOS and moments DOS for one PL, from one energy and kp
C-----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   moddos:specifies what kinds of quantities to calculate
Ci          1s digit
Ci           1 compute dosi (integrated dos and related quantities)
Ci           2 same as 1, but also set dosi(1) to zero before starting
Ci         10s digit
Ci           1 compute partial dos, (see pdos, Outputs)
Ci           2 accumulate diagonal part of g, (see gd, Outputs)
Ci             (used eg for Pade approximant to fix Fermi level)
Ci           3 both 1 and 2
Ci          4s bit if set, pdos decomposed into both l and m
Ci             By default, pdos is contracted over l
Ci        100s digit
Ci           1 compute asa moments of the dos
Ci       1000s digit
Ci           1 Accumulate the site-diagonal density matrix, array dmat
Ci           2 Accumulate the full density matrix, array dmat
Ci           4 compute spin density matrix, array rhos
Ci      10000s digit
Ci           0 normal mode
Ci           1 normal mode, but skip DLM sites with ncomp(ib)>1
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
Ci   gll   :Green's function for sites ib1..ib2.
Ci          gll(1,1) corresponds to 1st orbital in ib1
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
Ci   nrclas:nrclas(i) = number of atoms in the ith class
Ci   ncomp :ncomp(i) = number of CPA components at the ith basis site
Ci         :ncomp(i) should be 1 for non-CPA sites
Ci   iz    :index to quantities that keep energy-dependent
Ci          information (pdos,gd).
Ci   nz    :leading dimension of quantities that keep energy-dependent
Ci          information (pdos,gd). (for non-equilibrium mode: nz>nzp)
Ci   nzp   :number of energy points on equilibrium contour
Ci         :points iz>nzp are for nonequilibrium contour
Ci   kpl   :index specifing what channel dos should be accumulated
Ci   ipl   :index to current principal layer (nonequil mode only)
Ci   npl   :total number of principal layers (nonequil mode only)
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
Co          See Remarks about the meaning of dosi in the noncollinear case.
Co   pdos  :energy resolved partial dos
Co   gd    :RL- diagonal Green's function added from diag. of gll
Co   qnu   :ASA energy moments
Co   dmat  :density-matrix
Co   qne   :addition to charge by integration along nonequilib. contour
Co   dosne :dos at first and last points of nonequilibrium contour
Co         :(calc. with G^< for ipl<npl and with G^r for ipl=npl)
Co   rhos  :on-site spin density-matrix, possibly contracted over m
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
Cu  01 Nov 13 (Belashchenko) adjustments in preparation for fully relativistic GF
Cu  01 Nov 11 (Belashchenko) Added new mode to skip DLM classes
Cu  16 Nov 07 Added orbital-dependent potential shifts (vRLshf)
Cu  10 Jul 04 (S.Faleev) changes to handle noneqilibrium mode
Cu            Argument list changed.
Cu  18 Jun 04 (A Chantis) orbital moments in fully relativstic case
Cu   9 Mar 00 redesigned the iwaves downfolding, added h-waves
Cu            Argument list changed.
Cr Remarks
Cr    Spin-resolution for noncollinear case:
Cr    Spin-resolution of quantities such as qnu, pdos, have meaning
Cr    relative to the local quantization axis.
Cr    dosi, which is a global quantity, must be accumulated for a
Cr    global quantization axis.  In this context the spin-resolution
Cr    of dosi corresponds to the z component on this axis.
Cb Bugs
Cb  When both downfolding and nlo=nl*nl, the following line is wrong:
Cb     if (.not. lall) lpd = l
Cb  Site-diagonal density-matrix for iwaves is incorrectly generated.
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nbas,moddos,ldg,kpl,ib1,ib2,isp,nsp,nspc,lsw,iz,nz,
     .  offpd,offgd,ipc(ib2),nrclas(*),indxsh(*),ldh,ipl,npl,nzp,nrhos
      integer lidim1,lidim2,hdim1,hdim2
      double precision gll(ldg,nspc,ldg,nspc,2),ghh(ldh,nspc,*),
     .  gd(2,nz,nsp,nspc,0:*),zp0(2),
     .  wtkp,zp(2),wz(2),pdos(nz,nsp,0:*),dosi(5,nsp,kpl),vshft(nbas),
     .  pp(6,nl,nsp,*),qnu(3,nl,nsp,*),rhos(2,3,nrhos,2,2,1),
     .  orbtm(nl,nspc,nbas),
     .  qne(nsp,kpl),dosne(2,nsp,kpl),vRLshf(lidim2+hdim2,nsp)
      double complex dmat(nl**2,nl**2,nbas,nsp,nspc)
      integer ncomp(nbas)
C ... Local parameters
      logical lall,lpdos,lpade,ldosi,lmom,lhh
      integer i,ib,ic,im,j,jsp,ksp,l,ldmat,li,lmRL,lmgd,lpd,m,
     .  offgd0,offgi0,offpd0,mxorb,ljh,oghi,lhdim,ilm,k2,
     .  landm(nl**2,2)
      double precision pi,w(2),wt,zmenu(2),xx(2,0:2),bync,greal,gimag,
     .  gr12,gi12,gcwk(2)
      double complex wei,gc,ci
C     For DLM
      integer moddlm
C     For call to orbl
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
      integer nlm,norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),offi

      procedure(integer) :: ll,nglob

C ... Setup : decode what is to be calculated
      ldosi = mod(moddos,10) == 1 .or. mod(moddos,10) == 2
      lpdos = mod(mod(moddos/10,10),2) == 1
      lpade = mod(mod(moddos/10,10),4) > 1
      lall  = mod(moddos/10,10) >= 4
      lmom  = mod(mod(moddos/100,10),2) == 1
      ldmat = mod(moddos/1000,10)
      moddlm = mod(moddos/10000,10)
      mxorb = nglob('mxorb')
      i = 1; if (mod(lsw,10) == 1) call lmorder(0,i,landm,landm); if (i==0) i = 2
      call lmorder(i,nl-1,landm,landm(1,2))
      oghi  = ldh*nspc*mxorb
      call sanrg(.true.,iz,1,nz,'gfidos','iz')

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

C --- For each of the coupled spins, do ---
C     In noncollinear case, isp=1 always => need internal jsp=1..2
C     ksp is the current spin index in both cases:
C     ksp = isp in the collinear case
C         = jsp in the noncollinear case
C     jsp = 1   for independent spins, and spin index when nspc=2
      do  jsp = 1, nspc
      ksp = max(jsp,isp)

      if (lpdos) offpd = offpd0
      if (lpade) offgd = offgd0

C --- For each site in ib1..ib2, do ---
      lmRL = offgi0
      do  ib = ib1, ib2
        lmgd = -1
        if (moddlm == 1 .and. ncomp(ib) > 1) then
          do  l = 0, nl-1
            do  m = -l, l
              lmRL = lmRL + 1
              if (indxsh(lmRL) > lhdim) cycle
              lmgd = lmgd + 1
              lpd = l
              if (lall) lpd = lmgd
            enddo
          enddo
          goto 12
        endif
        ic = ipc(ib)
        bync = 1d0 / nrclas(ic)
        ljh = 0
        ilm = 0
        do  l = 0, nl-1
        do  m = -l, l
          ilm = ilm+1; im = landm(ilm,2)
          k2 = l+1
          if (nrhos == nl*nl) k2 = ilm

C         lm = ilm - 1
C         Skip orbitals outside l,i,h blocks
          lmRL = lmRL + 1
          if (indxsh(lmRL) > lhdim) cycle

C         Determine whether ll, ii, or hh  applies
          lhh = indxsh(lmRL) > lidim2

          if (lhh) then
            li = indxsh(lmRL) - (lidim2+hdim1)
            ljh = ljh+1
            greal = ghh(li,jsp,ljh)
C           gimag = ghh(li+oghi,jsp,ljh)
            gimag = ghh(li,jsp,ljh+oghi/ldh/nspc)
            gr12  = 0
            gi12  = 0
          else
            li = indxsh(lmRL) - lidim1
            greal = gll(li,jsp,li,jsp,1)
            gimag = gll(li,jsp,li,jsp,2)
            if (nspc == 2) then
              gr12  = gll(li,jsp,li,3-jsp,1)
              gi12  = gll(li,jsp,li,3-jsp,2)
            endif
          endif

C         Index to gd --- doesn't care which block g comes from
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

C     ... Partial DOS
          if (lpdos) pdos(iz,ksp,offpd+lpd) = pdos(iz,ksp,offpd+lpd) + wt*gimag

C     ... Accumulate gd = k-summed diagonal part of g for this energy, orbital
          if (lpade) then
            gd(1,iz,ksp,jsp,offgd+lmgd) = gd(1,iz,ksp,jsp,offgd+lmgd) + wtkp*greal
            gd(2,iz,ksp,jsp,offgd+lmgd) = gd(2,iz,ksp,jsp,offgd+lmgd) + wtkp*gimag
C           Retain the off-diagonal parts in spinor space
            if (nspc == 2) then
              gd(1,iz,jsp,3-jsp,offgd+lmgd) = gd(1,iz,jsp,3-jsp,offgd+lmgd) + wtkp*gr12
              gd(2,iz,jsp,3-jsp,offgd+lmgd) = gd(2,iz,jsp,3-jsp,offgd+lmgd) + wtkp*gi12
            endif
          endif

C     ... PL-integrated DOS and moments of DOS
C         For noncollinear case, g must be in global quantization axis
          if (ldosi) then
            if (iz <= nzp) dosi(1,ksp,kpl) = dosi(1,ksp,kpl) + wt*gimag
            if (iz > nzp) then !   DOS at endpoints of nonequilibrium contour
              if (nz <= nzp) call rx('gfidos: iz>nzp but nz<=nzp')
              if (iz == nzp+1) dosne(1,ksp,kpl) = dosne(1,ksp,kpl) + wt*gimag
              if (iz == nz)    dosne(2,ksp,kpl) = dosne(2,ksp,kpl) + wt*gimag
            endif

C           Contribution to charge from nonequilibrium contour
            if (iz > nzp) qne(ksp,kpl) = qne(ksp,kpl) + xx(2,0)

            xx(1,1) = zp(1)*xx(1,0) - zp(2)*xx(2,0)  ! First energy moment
            xx(2,1) = zp(2)*xx(1,0) + zp(1)*xx(2,0)
            xx(1,2) = zp(1)*xx(1,1) - zp(2)*xx(2,1)  ! Second energy moment
            xx(2,2) = zp(2)*xx(1,1) + zp(1)*xx(2,1)
C           No integration along noneqilibrium contour for left lead
C           Other energy integration includes both equilibrium and nonequilibrium points
            if (ipl == -1 .and. iz > nzp) then
            else
              do  i = 0, 2
                dosi(2+i,ksp,kpl) = dosi(2+i,ksp,kpl) + xx(2,i)
              enddo
            endif
          endif

C     ... ASA energy moments
          if (lmom) then
            zmenu(1) = zp(1) - (pp(1,l+1,ksp,ic)+vshft(ib)+vrlshf(indxsh(lmRL),ksp))
            zmenu(2) = zp(2)
            xx(1,1) = zmenu(1)*xx(1,0) - zmenu(2)*xx(2,0)
            xx(2,1) = zmenu(2)*xx(1,0) + zmenu(1)*xx(2,0)
            xx(1,2) = zmenu(1)*xx(1,1) - zmenu(2)*xx(2,1)
            xx(2,2) = zmenu(2)*xx(1,1) + zmenu(1)*xx(2,1)
C           No integration along noneqilibrium contour for left lead
C           Other energy integration includes both equilibrium and nonequilibrium points
            if (ipl == -1 .and. iz > nzp) then
            else
            do  i = 1, 3
              qnu(i,l+1,ksp,ic) = qnu(i,l+1,ksp,ic) + xx(2,i-1)*bync
            enddo
            if (nspc == 2) then
              orbtm(l+1,ksp,ic) = orbtm(l+1,ksp,ic) + im*xx(2,0)*bync
            endif
            endif
          endif

C     ... Spin density matrix.  Should be hermitian in spin space,
C         but there is some deviation owing to integration errors.
          if (ldmat >= 4) then
            if (.not. lmom) call rx('gfidos: dmat without lmom')
            do  i = 1, 3
              rhos(1,i,k2,jsp,jsp,ic) = rhos(1,i,k2,jsp,jsp,ic) + xx(2,i-1)*bync
            enddo
            xx(1,0) = w(1)*gr12 - w(2)*gi12
            xx(2,0) = w(1)*gi12 + w(2)*gr12
            xx(1,1) = zmenu(1)*xx(1,0) - zmenu(2)*xx(2,0)
            xx(2,1) = zmenu(2)*xx(1,0) + zmenu(1)*xx(2,0)
            xx(1,2) = zmenu(1)*xx(1,1) - zmenu(2)*xx(2,1)
            xx(2,2) = zmenu(2)*xx(1,1) + zmenu(1)*xx(2,1)
            do  i = 1, 3
C             Im rhos = Re g;  Re rhos = -Im g
              rhos(1,i,k2,jsp,3-jsp,ic) = rhos(1,i,k2,jsp,3-jsp,ic) + xx(2,i-1)*bync
              rhos(2,i,k2,jsp,3-jsp,ic) = rhos(2,i,k2,jsp,3-jsp,ic) - xx(1,i-1)*bync
            enddo
          endif

        enddo
        enddo

C   ... Site-diagonal density matrix
        if (mod(ldmat,4) == 1) then
C         This should be ldim2, not lidim2
          call orbl(ib,0,lidim2,indxsh,norb,ltab,ktab,offi,offl,nlm)
          offi = offi - offgi0
          wei = dcmplx(w(1),w(2)) * (-ci)
          do  j = 1, nlm
          do  i = 1, nlm
            gcwk(:) = gll(i+offi,jsp,j+offi,jsp,:)
            gc = dcmplx(gcwk(1),gcwk(2))
c           gcwk(:) = gll(j+offi,jsp,i+offi,jsp,:)
c           gc = gc - dcmplx(gcwk(1),-gcwk(2))
c           gc = gc/2
C           xx(1,0) = w(1)*gll(i+offi,jsp,j+offi,jsp,1)
C    .              - w(2)*gll(i+offi,jsp,j+offi,jsp,2)
c           xx(2,0) = w(1)*gll(i+offi,jsp,j+offi,jsp,2)
c    .              + w(2)*gll(i+offi,jsp,j+offi,jsp,1)
            dmat(i,j,ib,ksp,jsp) = dmat(i,j,ib,ksp,jsp) + wei*gc
            if (nspc == 2) then
              gcwk(:) = gll(i+offi,jsp,j+offi,3-jsp,:)
              gc = dcmplx(gcwk(1),gcwk(2))
c             gcwk(:) = gll(j+offi,3-jsp,i+offi,jsp,:)
c             gc = gc - dcmplx(gcwk(1),-gcwk(2))
c             gc = gc/2
c             xx(2,0) = w(1)*gll(i+offi,jsp,j+offi,3-jsp,2)
c    .                + w(2)*gll(i+offi,jsp,j+offi,3-jsp,1)
              dmat(i,j,ib,jsp,3-jsp) = dmat(i,j,ib,jsp,3-jsp) + wei*gc
            endif
          enddo
          enddo
        endif

   12   if (lpdos) offpd = offpd+lpd+1
        if (lpade) offgd = offgd+lmgd+1
      enddo
      enddo

C ... Entire density-matrix: dmat dimensioned dmat(ldg,nspc,ldg,nsp),
      if (mod(ldmat,4) == 2) then
         call yscal((ldg*nspc)**2,w(2),-w(1),gll,gll(1,1,1,1,2),1)
        if (nspc == 1) then
          call dpsadd(dmat,gll,ldg**2,1+ldg**2*(isp-1),1,1d0)
        else
          call dpsadd(dmat,gll,8*ldg**2,1,1,1d0)
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
