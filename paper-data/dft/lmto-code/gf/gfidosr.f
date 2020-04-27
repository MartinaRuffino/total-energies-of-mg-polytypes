      subroutine gfidosr(nl,nbas,moddos,lsw,ib,indxsh,lidim1,lidim2,
     .  hdim1,hdim2,ldg,ldh,ghh,wtkp,zp0,wz,isp,nsp,nspc,pprel,vRLshf,
     .  vshft,ipc,nrclas,iz,nz,nzp,ipl,npl,kpl,offpd,offgd,dosinew,pdos,
     .  pdosl,gd,qnur,dosne,orbtm,dmat,rhos,ncomp,theta,dlmcl,wts,gloc,
     .  gcorr,mxy,qcorr,norb,GFr)
C- Accumulate DOS and moments DOS for one CPA site, energy and kp
C-----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  gc ncomp norb dlmcl cpawt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed: thet
Cio    Passed to:  *
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   moddos:specifies what kinds of quantities to calculate
Ci          1s digit
Ci           1 compute dosi (integrated dos and related quantities)
Ci           2 same as 1, but also set dos(1) to zero before starting
Ci         10s digit (not implemented)
Ci           1 compute partial dos, (see pdos, Outputs)
Ci           2 accumulate diagonal part of g, (see gd, Outputs)
Ci             (used eg for Pade approximant to fix Fermi level)
Ci           3 both 1 and 2
Ci          4s bit if set, pdos decomposed into both l and m
Ci             By default, pdos is contracted over m
Ci        100s digit
Ci           1 compute asa moments of the dos
Ci       1000s digit (not implemented)
Ci           1 Accumulate the site-diagonal density matrix, array dmat
Ci           2 Accumulate the full density matrix, array dmat
Ci           4 compute spin density matrix, array rhos
Ci   ib1   :accumulate information for block ib1 ... ib2
Ci   ib2   :accumulate information for block ib1 ... ib2
Ci   zp0   :complex energy
Ci   wz    :weight for complex energy integration
Ci   pprel :Relativistic potential parameters in kappa-mu representation
Ci   vshft :array of site potential shifts
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nrclas:number of atoms in the icomp class
Ci   iz    :index to quantities that keep energy-dependent
Ci          information (pdos,gd).
Ci   nz    :leading dimension of quantities that keep energy-dependent
Ci          information (pdos,gd). (for non-equilibrium mode: nz>nzp)
Ci   nzp   :number of energy points on equilibrium contour
Ci         :points iz>nzp are for nonequilibrium contour
Ci   kpl   :index specifing what channel dos should be accumulated
Ci   gloc  :conditionally-averaged GF for one site, in lambda-mu repsn
Ci   GFr   : overlap integrals (Turek 6.141) (k1,k2,alpha1,alpha2,l,imu,type:
Co Outputs
Co   dosi  :integrated dos and moments of dos for sites ib1..ib2
Co          1: dos at this e : -1/pi Im G
Co          2: nos at this e : -1/pi Im (wz*G)
Co          3: 1st energy mom : -1/pi Im (wz*G z)
Co          4: 2nd energy mom : -1/pi Im (wz*G z**2)
Co          5: projection of dos to fermi level (not calculated here)
Co   qnur  :relativistic ASA energy moments
Co   q12   :charge from the non-diagonal <gdot|gdot>
Cl Local variables
Cl   greal :Re(diagonal g) for current orbital, from gll or ghh
Cl   gimag :Im(diagonal g) for current orbital, from gll or ghh
Cl   gr12  :like greal, but for off-diagonal in spins
Cl   gi12  :like gimag, but for off-diagonal in spins
Cl   1) <g|g>,<f|f>, 2) <g|gdot>,<f|fdot>,3) <gdot|g> ... 4) <gdot|gdot>,...)
Cu Updates
Cu  07 Aug 16 (MvS) redesign for >1 atom/cell, compatiblity with gfidos2
Cu  01 May 15 (Vishina) some adjustments for Dirac
Cu  13 Jun 14 (Belashchenko) Cloned from gfidos2 for Dirac
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nbas,ldg,kpl,ib,isp,nsp,nspc,lsw,iz,nz,offpd,
     .  offgd,ipc(ib),nrclas(*),indxsh(ldg),ldh,ipl,npl,nzp
      integer lidim1,lidim2,hdim1,hdim2
      integer ncomp,dlmcl,norb
      double precision ghh(ldh,nspc,*),gd(2,nz,nsp,nspc,0:*),
     .  wtkp,wz(2),pdos(nz,nsp,0:1),dosinew(5,nsp,kpl),vshft(nbas),
     .  pprel(5,0:nl-1,2*nl,2,2,*),qnur(4,0:nl-1,2*nl,2,2),rhos(2,3,nl,2,2,1),
     .  orbtm(nl,nspc,nbas),dosne(2,nsp,kpl),
     .  vRLshf(lidim2+hdim2,nsp),qcorr(2,nl,nsp,*)
C     double precision qnum(3,nl*nl,nsp,*)
      double complex gloc(2*norb,2*norb,ncomp),gcorr(norb,ncomp,2),
     .  dmat(norb,norb,nsp,nspc,4,ncomp),zp0
      double precision wts(ncomp),theta(ncomp,2),mxy(nl,2,*)
      double precision GFr(4,4,2,2,0:(nl-1),2*nl,4)
C ... Local parameters
      logical lall,lpdos,lpade,ldosi,lmom,lgcorr,lorbm,lmomm
      integer i,ic,l,ldmat,moddos,
     .  offgd0,offpd0,mxorb,nglob,oghi,lhdim
      double precision pi,wt,bync
      double complex ci,zp,xx(2,2,0:2),w
C For CPA
      integer moddlm,ic1,ic2,ic0,icomp
      double precision wtd1,wtd2,pdosl(nz,nsp,0:1)
C     double precision qq,qqt                             !TEMP
C Dirac-specific
      integer i1,i2,k1,k2,imu,idx(2,nl,2*nl)
      double precision dosisc(5,2,kpl),g22i(2,2),qtot(2)
      double precision q12,pq2,Mrel,enu,ms,ms2,mu,dosiupdo(5),
     .  gg(2,2,2,2,0:(nl-1),2*nl,4),ggl0(2,2,4),ggl(2,2,4),
     .  ff(2,2,2,2,0:(nl-1),2*nl,4),ffl0(2,2,4),ffl(2,2,4)
      double complex g22(2,2),g22ah(2,2),zmenu,zmenusq,zp1
C      logical ldosi,lmom
C      integer i,ib,ic,l,imu,k1,k2
C      integer idx(2,nl,2*nl)
C      double precision ggl(2,2,4), ffl(2,2,4),g22i(2,2),
C     .  ms,mu,ms2,dosisc(5,2,kpl),dosinew(5,2,kpl)
C      double precision
C      double precision dosiupdo(5),Mrel,mu,pq2
C      real(8) pi,wt,bync,enu,bync1
C     complex(8) w,wzp,g22ah(2,2),
C     .  zp,,zp1    !

CC For CPA
C      integer ic1,ic2,ic0,icomp
C      real(8) wtd1,wtd2

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
      oghi  = ldh*nspc*mxorb
      Mrel = 0
      call sanrg(.true.,iz,1,nz,'gfidos','iz')
c     if (moddlm /= 2) call rx1('GFIDOS2: moddlm = %i',moddlm)
c     if (ldmat /= 0) then
c       print *,' WARNING: gfidos2 not ready for ldmat'
c     endif

      ci = dcmplx(0,1)
      pi = 4*datan(1d0)
C ... wt is -1/pi * kpt weight
      wt = -wtkp/pi
C ... Keep passed energy zp0 = zp(izp) unchanged
      zp = zp0
C ... Nonequilibrium mode
      if (iz > nzp) then
         if (abs(wz(2)) > 1d-10) call rx('gfidos: wz(2) ' //
     .        'nonzero for non-equilibrium energy point')
         wz(2) = 0d0
C      G^< constructed only for ipl<npl; for ipl=npl it is G^r.
C      Integration along nonequilibrium contour with G^<
C      should be performed with real energy in integrand
         if (ipl < npl) zp = dble(zp0)
      endif

C ... w is -1/pi * kpt weight * energy weight
      w = dcmplx(wz(1),wz(2)) * wt
      if (mod(moddos,10) == 2 .and. iz <= nzp) then
        dosisc(1,:,kpl) = 0
        dosinew(1,:,kpl) = 0
      endif

      dosisc = 0

      dosiupdo = 0
      Mrel = 0
      q12 = 0
      pq2 = 0

C     Keep a local copy of offpd,offgd in case loop over second spin
      offpd0 = offpd
      offgd0 = offgd
C     Offset to indxsh
C     offgi0 = mxorb*(ib-1)
      lhdim = lidim2+hdim2
C     lmRL = offgi0
      call mstokm(2,nl,1,nl*nl,xx,xx,idx)
      ic0 = ipc(ib)
      bync = 1d0 / nrclas(ic0)
      ic1 = ic0 ; ic2 = ic0  ! No CPA for now

      do i1 = 1,2
        do i2 = 1,2
          gg(1,1,i1,i2,:,:,:) = GFr(1,1,i1,i2,:,:,:)
          ff(1,1,i1,i2,:,:,:) = GFr(2,2,i1,i2,:,:,:)
          gg(2,2,i1,i2,:,:,:) = GFr(3,3,i1,i2,:,:,:)
          ff(2,2,i1,i2,:,:,:) = GFr(4,4,i1,i2,:,:,:)
          gg(1,2,i1,i2,:,:,:) = GFr(1,3,i1,i2,:,:,:)
          ff(1,2,i1,i2,:,:,:) = GFr(2,4,i1,i2,:,:,:)
          gg(2,1,i1,i2,:,:,:) = GFr(3,1,i1,i2,:,:,:)
          ff(2,1,i1,i2,:,:,:) = GFr(4,2,i1,i2,:,:,:)
        enddo
      enddo

      do  ic = ic1, ic2
        icomp = ic - ic1 + 1
        wtd1 = (1 + dcos(theta(icomp,1)))*wts(icomp)/2
        wtd2 = (1 - dcos(theta(icomp,1)))*wts(icomp)/2
        do  l = 0, nl-1
C       enu = (pp(1,l+1,1,ic) + pp(1,l+1,2,ic))/2
        do  imu = 1, 2*(l+1)
          enu = pprel(5,l,imu,1,1,ic)
          mu  = imu - l - 1.5d0
          dosiupdo = 0

          ggl(:,:,:) = 0
          ggl0(:,:,:) = 0
          ffl(:,:,:) = 0
          ffl0(:,:,:) = 0
          do  k1 = 1,2
            do  k2 = 1,2
              ggl(:,:,:) = ggl(:,:,:) + gg(k1,k2,:,:,l,imu,:)
              ffl(:,:,:) = ffl(:,:,:) + ff(k1,k2,:,:,l,imu,:)
            enddo
          enddo
          do k1 = 1,2
!           ggl0(:,:,:) = ggl0(:,:,:) + gg(:,:,la1,la1,l,imu,:)
!           ffl0(:,:,:) = ffl0(:,:,:) + ff(:,:,la1,la1,l,imu,:)
            ggl0(:,:,:) = ggl0(:,:,:) + gg(k1,k1,:,:,l,imu,:)
            ffl0(:,:,:) = ffl0(:,:,:) + ff(k1,k1,:,:,l,imu,:)
          enddo

          k1 = idx(1,l+1,imu) ; k2 = idx(2,l+1,imu)
          g22(2,2) = gloc(k2,k2,icomp)
          if (k1 /= 0) then
            g22(1,1) = gloc(k1,k1,icomp)
            g22(1,2) = gloc(k1,k2,icomp)
            g22(2,1) = gloc(k2,k1,icomp)
          else
            g22(1,1) = 0; g22(1,2) = 0; g22(2,1) = 0
          endif
          call antiherm(g22,g22ah)
          g22i(1,1) = dimag(g22ah(1,1))
          g22i(1,2) = dimag(g22ah(1,2))
          g22i(2,1) = dimag(g22ah(2,1))
          g22i(2,2) = dimag(g22ah(2,2))

C     ... xx is -1/pi * gf * energy contour weight * kpt weight
          xx(:,:,0) = w*g22

          zmenu = zp - (enu+vshft(ib))
          zmenusq = zmenu*zmenu
C     ... PL-integrated dos and moments of dos
          if (ldosi) then
C            if (lsw == 1 .and. pass1) then
C              call dpzero(dosi(1,1,kpl),5*2)
C            endif

C             Need to insert the formula with overlap integrals of gg and ff
C            (this should make the imaginary part of DOS)

C --- dosi(1) using overlap integrals ---

!             qs = 0
!             ms = 0
!             if (imu == 1 .or. imu == 2*l+2) then
!                 qs = !qs(1,1,kpl) +
!     .   wt*dimag(g22(2,2))*(ggl(2,2,1)+zp*(ggl(2,2,2)+ggl(2,2,3))
!     .    + zp**2*ggl(2,2,4))+
!     .   wt*dimag(g22(2,2))*(ffl(2,2,1)+zp*(ffl(2,2,2)+ffl(2,2,3))
!     .    + zp**2*ffl(2,2,4))

!                  ms = !ms(1,1,kpl) +
!     .   - 2d0*mu/(-2*l-1d0)*wt*dimag(g22(2,2))*(ggl(2,2,1)+
!     .    zp*(ggl(2,2,2)+ggl(2,2,3)) + zp**2*ggl(2,2,4))+
!     .   2d0*mu/(2*l+3d0)*wt*dimag(g22(2,2))*(ffl(2,2,1)+
!     .    zp*(ffl(2,2,2)+ffl(2,2,3)) + zp**2*ffl(2,2,4))
!             else
!               do i1 = 1,2
!                  do i2 = 1,2
!                  qs = qs +
!     .            wt*dimag(g22(i1,i2))*(ggl(i1,i2,1)+
!     .    zp*(ggl(i1,i2,2)+ggl(i1,i2,3)) + zp**2*ggl(i1,i2,4))+
!     .            wt*dimag(g22(i1,i2))*(ffl(i1,i2,1)+
!     .    zp*(ffl(i1,i2,2)+ffl(i1,i2,3)) + zp**2*ffl(i1,i2,4))         !!! /4 -?
!                 enddo
!               enddo
!                 ms = !ms(1,1,kpl) +
!     .   -2d0*mu/(2*l+1d0)*wt*dimag(g22(1,1))*(ggl(1,1,1)+
!     .    zp*(ggl(1,1,2)+ggl(1,1,3)) + zp**2*ggl(1,1,4))+
!     .   2d0*mu/(-2*l+1d0)*wt*dimag(g22(1,1))*(ffl(1,1,1)+
!     .    zp*(ffl(1,1,2)+ffl(1,1,3)) + zp**2*ffl(1,1,4))
!     .    -dsqrt((2*l+1d0)**2-4*mu**2)/(2*l+1d0)*wt*dimag(g22(1,2))*
!     .    (ggl(1,2,1)+zp*(ggl(1,2,2)+ggl(1,2,3))+zp**2*ggl(1,2,4))
!     .    -dsqrt((2*l+1d0)**2-4*mu**2)/(2*l+1d0)*wt*dimag(g22(2,1))*
!     .    (ggl(2,1,1)+zp*(ggl(2,1,2)+ggl(2,1,3))+zp**2*ggl(2,1,4))
!     .   -2d0*mu/(-2*l-1d0)*wt*dimag(g22(2,2))*(ggl(2,2,1)+
!     .    zp*(ggl(2,2,2)+ggl(2,2,3)) + zp**2*ggl(2,2,4))+
!     .   2d0*mu/(2*l+3d0)*wt*dimag(g22(2,2))*(ffl(2,2,1)+
!     .    zp*(ffl(2,2,2)+ffl(2,2,3)) + zp**2*ffl(2,2,4))
!             endif
!              dosi(1,1,kpl) = dosi(1,1,kpl) + wtd1*(ms + qs)/2
!              dosi(1,2,kpl) = dosi(1,2,kpl) + wtd1*(qs - ms)/2

C --- dosi(1) NEW ---
            do i1 = 1,2
              dosiupdo(1) = dosiupdo(1) + wt*dimag(g22(i1,i1))
            enddo
            ms = 0
            if (imu == 1 .or. imu == 2*l+2) then
              ms = -2d0*mu/(-2*l-1d0)*wt*dimag(g22(2,2))*(ggl(2,2,1)+
     .          zp*(ggl(2,2,2)+ggl(2,2,3)) + zp**2*ggl(2,2,4))+
     .          2d0*mu/(2*l+3d0)*wt*dimag(g22(2,2))*(ffl(2,2,1)+
     .          zp*(ffl(2,2,2)+ffl(2,2,3)) + zp**2*ffl(2,2,4))
            else
              ms = -2d0*mu/(2*l+1d0)*wt*dimag(g22(1,1))*(ggl(1,1,1)+
     .          zp*(ggl(1,1,2)+ggl(1,1,3)) + zp**2*ggl(1,1,4))+
     .          2d0*mu/(-2*l+1d0)*wt*dimag(g22(1,1))*(ffl(1,1,1)+
     .          zp*(ffl(1,1,2)+ffl(1,1,3)) + zp**2*ffl(1,1,4))
     .          -dsqrt((2*l+1d0)**2-4*mu**2)/(2*l+1d0)*wt*dimag(g22(1,2))*
     .          (ggl(1,2,1)+zp*(ggl(1,2,2)+ggl(1,2,3))+zp**2*ggl(1,2,4))
     .          -dsqrt((2*l+1d0)**2-4*mu**2)/(2*l+1d0)*wt*dimag(g22(2,1))*
     .          (ggl(2,1,1)+zp*(ggl(2,1,2)+ggl(2,1,3))+zp**2*ggl(2,1,4))
     .          -2d0*mu/(-2*l-1d0)*wt*dimag(g22(2,2))*(ggl(2,2,1)+
     .          zp*(ggl(2,2,2)+ggl(2,2,3)) + zp**2*ggl(2,2,4))+
     .          2d0*mu/(2*l+3d0)*wt*dimag(g22(2,2))*(ffl(2,2,1)+
     .          zp*(ffl(2,2,2)+ffl(2,2,3)) + zp**2*ffl(2,2,4))
!             n12dot = q12dot + wt*dimag(g22(1,2))*ggl(1,2,4)+
!     .                 wt*dimag(g22(2,1))*ggl(2,1,4)+ wt*dimag(g22(1,2))*ffl(1,2,4)+
!     .                 wt*dimag(g22(2,1))*ffl(2,1,4)
!             n12dot = n12dot + wt*dimag(g22(1,2))*ggl(1,2,2)+
!     .                 wt*dimag(g22(2,1))*ggl(2,1,2)+ wt*dimag(g22(1,2))*ggl(1,2,3)+
!     .                 wt*dimag(g22(2,1))*ggl(2,1,3)+ wt*dimag(g22(1,2))*ffl(1,2,2)+
!     .                 wt*dimag(g22(2,1))*ffl(2,1,2)+ wt*dimag(g22(1,2))*ffl(1,2,3)+
!     .                 wt*dimag(g22(2,1))*ffl(2,1,3)
            endif
            dosinew(1,1,kpl) = dosinew(1,1,kpl)+(dosiupdo(1) - ms)/2*bync
            dosinew(1,2,kpl) = dosinew(1,2,kpl)+(dosiupdo(1) + ms)/2*bync
C --- dosi(1) NEW end

C --- dosi(1) without overlap integrals (average up + down)
            do i1 = 1,2
              dosisc(1,1,kpl) = dosisc(1,1,kpl) + wt*wtd1*dimag(g22(i1,i1))/2*bync
              dosisc(1,2,kpl) = dosisc(1,2,kpl) + wt*wtd1*dimag(g22(i1,i1))/2*bync
            enddo

            xx(:,:,1) = zp*xx(:,:,0)
            xx(:,:,2) = zp*xx(:,:,1)
            zp1 = zmenu

C --- dosi(2-4) using overlap integrals ---
!             do i = 0, 2
!               qs1 = 0; qs2 = 0
!               ms1 = 0; ms2 = 0
!               if (imu == 1 .or. imu == 2*l+2) then
!                  qs1 = !qs(1,1,kpl) +
!     .   dimag(xx(2,2,i))*(ggl(2,2,1)+zp1*(ggl(2,2,2)+
!     .    ggl(2,2,3))+ zp1**2*ggl(2,2,4))+
!     .   dimag(xx(2,2,i))*(ffl(2,2,1)+zp1*(ffl(2,2,2)+ffl(2,2,3))
!     .    + zp1**2*ffl(2,2,4))
!                  qs2 = !qs(1,1,kpl) +
!     .   dimag(xx(2,2,i))*(ggl(2,2,1)+zp1*(ggl(2,2,2)+
!     .    ggl(2,2,3))+ zp1**2*ggl(2,2,4))+
!     .   dimag(xx(2,2,i))*(ffl(2,2,1)+zp1*(ffl(2,2,2)+ffl(2,2,3))
!     .    + zp1**2*ffl(2,2,4))

!                   ms1 = !ms(1,1,kpl) +
!     .   - 2d0*mu/(-2*l-1d0)*dimag(xx(2,2,i))*(ggl(2,2,1)+
!     .    zp1*(ggl(2,2,2)+ggl(2,2,3)) + zp1**2*ggl(2,2,4))+
!     .   2d0*mu/(2*l+3d0)*dimag(xx(2,2,i))*(ffl(2,2,1)+
!     .    zp1*(ffl(2,2,2)+ffl(2,2,3)) + zp1**2*ffl(2,2,4))
!                   ms2 = !ms(1,1,kpl) +
!     .   - 2d0*mu/(-2*l-1d0)*dimag(xx(2,2,i))*(ggl(2,2,1)+
!     .    zp1*(ggl(2,2,2)+ggl(2,2,3)) + zp1**2*ggl(2,2,4))+
!     .   2d0*mu/(2*l+3d0)*dimag(xx(2,2,i))*(ffl(2,2,1)+
!     .    zp1*(ffl(2,2,2)+ffl(2,2,3)) + zp1**2*ffl(2,2,4))
!               else
!                 do i1 = 1,2
!                    qs1 = qs1 +
!     .    dimag(xx(i1,i1,i))*(ggl(i1,i1,1)+
!     .    zp1*(ggl(i1,i1,2)+ggl(i1,i1,3)) + zp1**2*ggl(i1,i1,4))+
!     .   dimag(xx(i1,i1,i))*(ffl(i1,i1,1)+
!     .    zp1*(ffl(i1,i1,2)+ffl(i1,i1,3)) + zp1**2*ffl(i1,i1,4))               !!! /4 -?
!                    qs2 = qs2 +
!     .    dimag(xx(i1,i1,i))*(ggl(i1,i1,1)+
!     .    zp1*(ggl(i1,i1,2)+ggl(i1,i1,3)) + zp1**2*ggl(i1,i1,4))+
!     .   dimag(xx(i1,i1,i))*(ffl(i1,i1,1)+
!     .    zp1*(ffl(i1,i1,2)+ffl(i1,i1,3)) + zp1**2*ffl(i1,i1,4))               !!! /4 -?
!                 enddo
!                     ms1 = !ms(1,1,kpl) +
!     .   -2d0*mu/(2*l+1d0)*dimag(xx(1,1,i))*(ggl(1,1,1)+
!     .    zp1*(ggl(1,1,2)+ggl(1,1,3)) + zp1**2*ggl(1,1,4))+
!     .   2d0*mu/(-2*l+1d0)*dimag(xx(1,1,i))*(ffl(1,1,1)+
!     .    zp1*(ffl(1,1,2)+ffl(1,1,3)) + zp1**2*ffl(1,1,4))
!     .   -2d0*mu/(-2*l-1d0)*dimag(xx(2,2,i))*(ggl(2,2,1)+
!     .    zp1*(ggl(2,2,2)+ggl(2,2,3)) + zp1**2*ggl(2,2,4))+
!     .   2d0*mu/(2*l+3d0)*dimag(xx(2,2,i))*(ffl(2,2,1)+
!     .    zp1*(ffl(2,2,2)+ffl(2,2,3)) + zp1**2*ffl(2,2,4))
!                     ms2 = !ms(1,1,kpl) +
!     .   -2d0*mu/(2*l+1d0)*dimag(xx(1,1,i))*(ggl(1,1,1)+
!     .    zp1*(ggl(1,1,2)+ggl(1,1,3)) + zp1**2*ggl(1,1,4))+
!     .   2d0*mu/(-2*l+1d0)*dimag(xx(1,1,i))*(ffl(1,1,1)+
!     .    zp1*(ffl(1,1,2)+ffl(1,1,3)) + zp1**2*ffl(1,1,4))
!     .   -2d0*mu/(-2*l-1d0)*dimag(xx(2,2,i))*(ggl(2,2,1)+
!     .    zp1*(ggl(2,2,2)+ggl(2,2,3)) + zp1**2*ggl(2,2,4))+
!     .   2d0*mu/(2*l+3d0)*dimag(xx(2,2,i))*(ffl(2,2,1)+
!     .    zp1*(ffl(2,2,2)+ffl(2,2,3)) + zp1**2*ffl(2,2,4))
!               endif
!                dosi(2+i,1,kpl) = dosi(2+i,1,kpl) + wtd1*(qs1 + ms2)/2
!                dosi(2+i,2,kpl) = dosi(2+i,2,kpl) + wtd1*(qs1 - ms2)/2
!            enddo

C --- dosi(2-4) NEW
            do i = 0,2
              do i1 = 1,2
                dosiupdo(2+i) = dosiupdo(2+i) + dimag(xx(i1,i1,i))
              enddo
              ms2 = 0
              if (imu == 1 .or. imu == 2*l+2) then
                ms2 = -2d0*mu/(-2*l-1d0)*dimag(xx(2,2,i))*(ggl(2,2,1)+
     .            zp1*(ggl(2,2,2)+ggl(2,2,3)) + zp1**2*ggl(2,2,4))+
     .            2d0*mu/(2*l+3d0)*dimag(xx(2,2,i))*(ffl(2,2,1)+
     .            zp1*(ffl(2,2,2)+ffl(2,2,3)) + zp1**2*ffl(2,2,4))
!              dosiupdo(2+i) = dimag(xx(2,2,i))
              else
                ms2 = -2d0*mu/(2*l+1d0)*dimag(xx(1,1,i))*(ggl(1,1,1)+
     .            zp1*(ggl(1,1,2)+ggl(1,1,3)) + zp1**2*ggl(1,1,4))+
     .            2d0*mu/(-2*l+1d0)*dimag(xx(1,1,i))*(ffl(1,1,1)+
     .            zp1*(ffl(1,1,2)+ffl(1,1,3)) + zp1**2*ffl(1,1,4))
     .            -2d0*mu/(-2*l-1d0)*dimag(xx(2,2,i))*(ggl(2,2,1)+
     .            zp1*(ggl(2,2,2)+ggl(2,2,3)) + zp1**2*ggl(2,2,4))+
     .            2d0*mu/(2*l+3d0)*dimag(xx(2,2,i))*(ffl(2,2,1)+
     .            zp1*(ffl(2,2,2)+ffl(2,2,3)) + zp1**2*ffl(2,2,4))
!            do i1 = 1,2
!               dosiupdo(2+i) = dosiupdo(2+i) + dimag(xx(i1,i1,i))
!            enddo
              endif
              dosinew(2+i,1,kpl) = dosinew(2+i,1,kpl)+(dosiupdo(2+i) - ms2)/2*bync
              dosinew(2+i,2,kpl) = dosinew(2+i,2,kpl)+(dosiupdo(2+i) + ms2)/2*bync
            enddo

C --- dosi(2-4) NEW end
          endif


C --- dosi(2-4) without overlap integrals; the way it was in the scalar case
          do i = 0,2
            do i1 = 1,2
              dosisc(2+i,1,kpl) = dosisc(2+i,1,kpl) + wtd1*dimag(xx(i1,i1,i))/2*bync
              dosisc(2+i,2,kpl) = dosisc(2+i,2,kpl) + wtd1*dimag(xx(i1,i1,i))/2*bync
            enddo
          enddo

C     ... ASA energy moments
          if (lmom) then
            zmenu = zp - (enu+vshft(ib))

            xx(:,:,1) = zmenu*xx(:,:,0)
            xx(:,:,2) = zmenu*xx(:,:,1)
            do i = 0, 2
              call antiherm(xx(:,:,i),g22ah)
              xx(:,:,i) = g22ah(:,:)
            enddo

            do  i = 1, 3
              if (k1 == 0) then
                qnur(i,l,imu,2,2) = qnur(i,l,imu,2,2) + dimag(xx(2,2,i-1))*bync
                pq2 = pq2 + qnur(3,l,imu,2,2)*pprel(4,l,imu,2,2,1)
              else
                qnur(i,l,imu,:,:) = qnur(i,l,imu,:,:) + dimag(xx(:,:,i-1))*bync
                if (i == 3) then
                 q12 = q12 +
     .              qnur(i,l,imu,1,2)*ggl0(1,2,4) +
     .              qnur(i,l,imu,2,1)*ggl0(2,1,4) +
     .              qnur(i,l,imu,1,2)*ffl0(1,2,4) +
     .              qnur(i,l,imu,2,1)*ffl0(2,1,4)
                endif
                do i1 = 1,2
                  do i2 = 1,2
                    pq2 = pq2 + qnur(3,l,imu,i1,i2)*pprel(4,l,imu,i1,i2,1)/bync
                  enddo
                enddo

              endif
            enddo
          endif

C     ... Accumulate magnetic moment from qnur
C          if (iz == nzp) then
C          mu  = imu - l - 1.5d0
C          if (imu == 1 .or. imu == 2*l+2) then
C            Mrel = Mrel - 2*mu*qnur(1,l,imu,2,2)
C          else
C            do i1 = 1,2
C              do i2 = 1,2
C                Mrel = Mrel - 2*mu*qnur(1,l,imu,i1,i2)
C              enddo
C            enddo
C          endif
C          endif

        enddo                   ! Loop over imu
      enddo                     ! loop over l
      enddo                     ! Loop over CPA components

C --- Printing qnur
!        stdo = nglob('stdo')
!        write(stdo,23)
!   23   format(' GFIDOSR qnur'/'   l  mu   ms1 ms2',5x,'q0',11x,'q1',11x,'q2')
!        do  l = 0, nl-1
!          do  imu = 1, 2*(l+1)
!            fac = dble(imu-l) - 1.5d0
!            do  i = 1, 2
!              do  k = 1, 2
!                if (dlength(3,qnur(1,l,imu,i,k),1) > 1d-6)
!     .          write(stdo,24) l,fac,i,k,(qnur(ir,l,imu,i,k), ir=1,3)
!              enddo
!            enddo
!          enddo
!        enddo
!   24   format(i4,f5.1,2i4,4f13.7)

C       print*,'GFIDOSR pq2 = ',pq2
!        print*,'GFIDOSR dosinew = ',dosinew

C       if (iz == nzp .or. .true.) then
C         call qtotrel(nl,1,qnur,qtot)
C         print *, 'GFIDOSR q,M = ',ib,qtot
C       endif

      end

      subroutine antiherm(A,B)
C Antihermitian part of 2*2 matrix
      implicit none
      complex(8) A(2,2),B(2,2)

      B = 0
      B(1,1) = cmplx(0d0,aimag(A(1,1)))
      B(1,2) = cmplx(.5d0*(real(A(1,2))-real(A(2,1))),.5d0*(aimag(A(1,2))+aimag(A(2,1))))
      B(2,1) = cmplx(.5d0*(real(A(2,1))-real(A(1,2))),.5d0*(aimag(A(2,1))+aimag(A(1,2))))
      B(2,2) = cmplx(0d0,aimag(A(2,2)))

      end subroutine antiherm
