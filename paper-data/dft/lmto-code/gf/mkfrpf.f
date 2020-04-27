      subroutine mkfrpf(s_site,nl,nbas,ipc,ipdc,lhdim,pfdim,indxsh,iopt,
     .  pp,pprel,vshft,z,P)
C- Relativistic potential functions parameterized from pot pars
C-----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ddpfr dpfr pfr
Cio    Passed to:  *
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ipdc  :hash (class control) array for CPA
Ci   lhdim :number of lower+intermediate+higher orbitals
Ci   pfdim :leading dimension for the P arrays (different from lhdim in CPA)
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   iopt  :a compound of digits determining what is generated:
Ci       10000s digit controls whether and how matrix form P is recorded in s_site(:).
Ci             What pot fun is generated depends on modeP (modeP=10s digit iopt)
Ci             As of now output is always recorded in vector form in array P, for non CPA sites.
Ci             This may change in future.
Ci             For some modeP, P is also recorded in matrix form:
Ci           0 s_site(:)%pfr (modeP=0)  s_site(:)%ddpfr (modeP=2)  s_site(:)%dpfr (modeP=3)
Ci           1 P is recorded in s_site(:)%pfra (modeP=0)
Ci         100s digit controls which repsn' bet P is made in
Ci           1 make P for bare representation (alpha=0)
Ci           2 make P for gamma representation (alpha=gamma)
Ci           3 make P for gamma representation averaged over spins
Ci          10s digit determines whether to make P or some derivative:
Ci           0 P <- potential function P
Ci           1 P <- P-dot
Ci           2 P <- -1/2 P-dotdot/P-dot
Ci           3 P <- sqrt(P-dot) --- actually returns A^T where A^T A = P-dot
Ci           4 P <- P^alpha [P^gamma]^-1 to convert g^alp -> g^gam
Ci           5 P <- (gamma-alpha) P^alpha [P^gamma]^-1 to convert g^alp -> g^gam
Ci           6 not used
Ci           7 P <- beta-alpha, where beta = gam or gamrel
Ci          1s digit
Ci           0 second-order P
Ci           1 third-order P
Ci           2 first-order P
Ci   pp    :scalar rel potential parameters, real, by class (atomsr.f)
Ci         :needed only to extract parameter alpha or gamma
Ci         :pp(5) gamma
Ci         :pp(6) alpha
Ci   pprel :relativistic potential parameters, real, by class (atomsr.f)
Ci         :imu ranges from  1 .. 2*(l+1)  and mu = imu-1 - (l+1/2)
Ci         :pprel(1,l,imu,isp,jsp) cgamma(l,imu,isp,jsp)
Ci         :pprel(2,l,imu,isp,jsp) gamma(l,imu,isp,jsp)
Ci         :pprel(3,l,imu,isp,jsp) W+, where (W+) W = delta(l,imu,isp,jsp)
Ci         :pprel(4,l,imu,isp,jsp) pgamma(l,imu,isp,jsp)
Ci         :pprel(5,l,imu,isp,jsp) enu(l,imu,isp,jsp)
Ci   vshft :array of site potential shifts
Ci   z     :complex energy
Co Outputs
Co   P     :some potential function in lms repn, depending on iopt
Cr Remarks
Cr   This is a relativistic analog of mkptfp, which see
Cr   s_site(ib)%pfr holds P in (norb*2,norb*2) matrix, in imu,lambda,lambda' rep
Cr
Cl Local variables
Cl   bet  representation in which P is to be made
Cl   irep 2 rotate to gamma repsn
Cl        3 rotate to gamma repsn averaged over spins, gamma defined by scal rel pp(5)
Cu Updates
Cu   18 Jun 18 Some redesign; synchronize with updated spherical harmonics
Cu   05 Jun 16 10000s iopt=1 -> write P to s_site(:)%pfra
Cu   19 Mar 16 (MvS) redesigned to make all potential functions internally consistent
Cu   16 Jan 16 (MvS) new 10000s digit iopt
Cu   10 Jan 16 distinguish irep=2 and irep=3.  New mode 6.
Cu   24 Apr 15 Bug fix : corrected proper treatment of vshft
Cu   19 Aug 13 (Belashchenko) Adapted for CPA
Cu   18 Jun 04 (A Chantis) working version
Cu   17 Mar 03 (A Chantis) first created
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nbas,ipc(nbas),lhdim,indxsh(*),iopt
      integer pfdim,ipdc(*)
      integer, parameter :: nsp=2
      double precision pp(6,nl,nsp,*),pprel(5,nl,2*nl,2,2,*),vshft(nbas)
      double complex z,P(pfdim,2,2),Ploc(nl*nl,2,2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer lmr,ibas,l,lmax,ic,ipr,modeP,irep,stdo,nglob,norb,ncomp,icomp,lm,recs,im,ilm,m
      double precision alp,gam2,enu,enu1,enu2
      double complex zloc,zero,one
      logical order3,savpfp
C ... Dynamically allocated arrays
      complex(8), allocatable,target, dimension(:,:,:,:) :: Pfib
C ... Dirac-specific
      integer imu,m1,m2,ms1,ms2,morder
      integer idxkm(2,0:nl-1,2*(nl+1)),i0(10)
      double precision crel(2,2),gamrel(2,2),delrel(2,2),unit(2,2)
      double precision mu,u,clebsh(2,2),zmcir(2,2),zmcii(2,2),
     .  W(2,2),WT(2,2),pprel1(2,2),pprel2(2,2),bet(2,2),xx
      double complex zmat(2,2),ep(2,2),sqep(2,2),epp(2,2),zmc(2,2),zmci(2,2),
     .  pbeti(2,2),palpi(2,2),pgam(2,2),zenu(2,2),ztmp(2,2),ztmp2(2,2),
     .  reskmu(2,2),reslms(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2)
C ..  For CPA
      integer idx,ib,ib1,ib2,icm,lmrd,idm,ic0
      logical ldsite
      real(8) cpagam(nl,nsp)

      procedure(integer) :: ll

C --- Setup and printout  ---
      call getpr(ipr)
C     ipr    = 60
      stdo   = nglob('stdo')
      irep   = mod(iopt/100,10)
      modeP  = mod(iopt/10,10)
C     order1 = mod(iopt,10) == 2
      order3 = mod(iopt,10) == 1
C     governs where in s_site to store P
      recs = modeP+1; if (modeP == 0 .and. mod(iopt/10000,10) == 1) recs = -1
      savpfp =  .true. ! Always write to P for now
C     orderi = mod(iopt,10) == 3
      i0 = 0

C     call setmorderdefault(1)
      call lmorder(0,morder,i0,i0)

C ... z, z_dot and z_dot_dot (2nd order potential fun)
      zloc = z   ! for debugging : zloc can be changed without affecting z
C     xx = 0
C     print *, '!! eps', sngl(dble(xx)); zloc = zloc + (1d-5,0d0)*xx

      call info0(50,0,0,' site  l  m%11fP(1,1)%17fP(2,1)%17fP(1,2)%17fP(2,2)')

      zero = dcmplx(0d0,0d0) ; one  = dcmplx(1d0,0d0)
      unit = zero ; unit(1,1) = one ; unit(2,2) = one
      zmat = zloc * unit
      ep = unit
      epp = zero

C --- For each site, make P or related function ---
      lmr = 0
      idx = nbas
      call dpzero(P,2*size(P))
      do  ib = 1, nbas

        lmr = nl*nl*(ib-1)
        ncomp = s_site(ib)%ncomp
        ldsite = ncomp > 1
        norb = s_site(ib)%norb
        lmax = ll(norb)

        call dpzero(Ploc,2*size(Ploc)) !! ?? not needed
        allocate(Pfib(norb,2,norb,2)); call dpzero(Pfib,2*size(pfib))

        if (recs == 0+1) then
          call dpzero(s_site(ib)%pfr,4*norb*norb*ncomp*2)
        elseif (recs == -1) then
          call dpzero(s_site(ib)%pfra,4*norb*norb*ncomp*2)
        elseif (recs == 2+1) then
          call dpzero(s_site(ib)%ddpfr,4*norb*norb*ncomp*2)
        elseif (recs == 3+1) then
          call dpzero(s_site(ib)%dpfr,4*norb*norb*ncomp*2)
        endif

C       CPA business needs rethinking
C       Returns gamma ppar, or if a CPA site, CPA-averaged gamma
        ic0 = s_site(ib)%dlmcl         ! Parent class
!       call cpagamma(ncomp,ic0,s_pot%dlmwt,nl,nsp,pp,cpagam) ! Not used yet

!       call basidx(ib,ncomp,idx,ldsite,ib1,ib2) ! ?? not needed
!        ic = ipc(ib)
!        icm = ipdc(ib)
        call mstokm(2,nl,ncomp,norb,xx,xx,idxkm) ! Indices for k,mu -> lms
        do  icomp = 1, ncomp
          ic = ipc(ib) ; ic0 = ic     ! class index to this site
          if (ldsite) then
            ic = ipdc(ib) + icomp - 1 ! class index to this component
            ic0 = ipdc(ib)            ! class index to first component
          endif

          do  l = 0, nl-1

            if (indxsh(lmr+l**2+1) > lhdim) cycle  ! Skip higher waves

C           Select repsn
C           gam = pp(5,l+1,1,ic)
            gam2 = (pp(5,l+1,1,ic) + pp(5,l+1,nsp,ic))/2
            alp = pp(6,l+1,1,ic)
            bet = alp*unit

            do  imu = 1, 2*(l+1)

C    ...    enu used for 3d order potential functions
            if (order3) then
C             What about vshft? From old code:
c             enu1 = pp(1,l+1,1,ic) + vshft(ib)
c             enu2 = pp(1,l+1,2,ic) + vshft(ib)
              enu = pprel(5,l+1,imu,1,1,ic) + vshft(ib)
              enu1 = enu ; enu2 = enu ! what about vshft
            endif

C     ..  . Extract pp(alpha)
            crel(:,:)   = pprel(1,l+1,imu,:,:,ic) + vshft(ib)*unit
            gamrel(:,:) = pprel(2,l+1,imu,:,:,ic)
            delrel(:,:) = pprel(3,l+1,imu,:,:,ic)

            if (irep == 1) bet = 0
            if (irep == 2) bet = gamrel
            if (irep == 3) bet = gam2*unit

            if (order3) then
              pgam(:,:) = dcmplx(pprel(4,l+1,imu,:,:,ic),0d0)
              zenu = 0 ; zenu(1,1) = zloc - enu1 ; zenu(2,2) = zloc - enu2
            endif

            if (order3) then
              ztmp = matmul(zenu,pgam)
              epp = 6*ztmp
              ztmp2 = matmul(zenu,ztmp)
              ep(:,:)  = unit(:,:) + 3*ztmp2
              ztmp = matmul(zenu,ztmp2)
              zmat = zloc*unit + ztmp
            endif

C       ... Setup for rotation of kappa-mu to lms rep'sn
            mu = imu - l - 1.5d0
            u = mu/(l+0.5d0)
            clebsh(1,1) =  dsqrt((1+u)/2) ; clebsh(1,2) = -dsqrt((1-u)/2)
            clebsh(2,2) =  clebsh(1,1) ; clebsh(2,1) = - clebsh(1,2)

C           zmc <- z-C; and local copy of W+ and W, where (W+) W = delta
            WT = delrel ; W = transpose(delrel)
            zmc = zmat - crel

C       ... Real, imaginary parts of (P^gam)^-1 = [W+ (z-C)^-1 W]
C           Note this is consistent with Schick's W, Eq. 26, PRB 54, 1610
            call zinv22(zmc,zmci)                    ! zmci <- (z-C)^-1
            zmcir(:,:) = dble(zmci(:,:)) ; zmcii(:,:) = dimag(zmci(:,:))
            pprel1 = matmul(matmul(WT,zmcir),W) ! Re W+ (z-C)^-1 W
            pprel2 = matmul(matmul(WT,zmcii),W) ! Im W+ (z-C)^-1 W

C       ... modeP = 2 : palpi <- (P^gam)^-1 = W+ (z-C)^-1 W
            if (modeP == 2) then
              palpi = dcmplx(pprel1,pprel2)
            endif

C       ... modeP = 4,5 : palpi <- (P^alpha)^-1
C       ... palpi = P^alp)^-1 = W+ (z-C)^-1 W + gamrel-alp, where W+ W = delta
            if (modeP == 4 .or. modeP == 5) then
              palpi = dcmplx(pprel1 + gamrel-alp*unit,pprel2)
            endif

C       ... Re (P^bet)^-1 = W+ (z-C)^-1 W + gamrel-bet, where W+ W = delta
            pbeti = dcmplx(pprel1 + gamrel-bet,pprel2)
C           call zprm('(P^bet)^-1 = del (z-C)^-1 del + gamrel-bet',2,pbeti,2,2,2)

C       --- Make Pot fun or related object in kappa-mu rep ---

C       ... modeP = 0: P^bet = (W+ (z-C)^-1 W + gamrel-bet)^-1
            if (modeP == 0) then
              call zinv22(pbeti,reskmu)

C       ... modeP = 2: -1/2*(P_dot_dot)*(P_dot)^(-1)
            elseif (modeP == 2) then
C             if (l == 2 .and. imu == 3) then
              call sqr2x2(ep,sqep)     ! sqrt(1 + 3 zmenu^2 pgam)

C             if (l == 2 .and. imu == 3) then
C               call prmx('gamma-alpha',gamrel-alp*unit,2,2,2)
C               call sqr2x2(dcmplx(gamrel-alp*unit),srgma) ! sqrt(gamma-alpha)
C               call prmx('C',crel,2,2,2)
C               call zprmx('z''-C',2,zmc,2,2,2)
C               call prmx('W=W',W,2,2,2)
C               call zprmx('sqrt(gamalp)',2,srgma,2,2,2)
C               call zprmx('Pbet',2,pbeti,2,2,2)
C               call zprmx('ep',2,ep,2,2,2)
C               call sqr2x2(ep,sqep)
C               call zprmx('sqrt(ep)',2,sqep,2,2,2)
C               call zprmx('mu+T = sqrt(pdot)',2,reskmu,2,2,2)
C               call zprmx('mu+T*mu = pdot',2,matmul(reskmu,transpose(reskmu)),2,2,2)
C             endif

C             term proportional to (bet-gam) : Simpler form that uses (alpha-gamma)^-1
C             seqp [W^T]^-1 {W gmb^-1 W^T + (z'-C)}^-1 sqep
C             call zinv22(dcmplx(gamrel-bet*unit),ztmp) ! (gam-bet)^-1
C             ztmp2 = matmul(matmul(dcmplx(W),ztmp),dcmplx(WT)) ! W (gam-bet)^-1 W^T
C             ztmp = zmC + ztmp2
C             call zinv22(ztmp,ztmp2) ! [W (gam-bet)^-1 W+T + (z'-C)]^-1
C             ztmp = matmul(sqep,matmul(ztmp2,sqep))
C             call zprmx('xx',2,ztmp,2,2,2)

C             term proportional to (bet-gam) : Stable form that uses (gam-bet)
C             seqp [W^T]^-1 {1 + gmb W^-1 (z'-C) [W^T]^-1}^-1 gmb W^-1 sqep
C             ? doesn't work
C              gmb = dcmplx(gamrel-bet*unit)
C              call zinv22(dcmplx(W),ztmp2) ! W^-1
C              ztmp = unit + matmul(gmb,matmul(matmul(ztmp2,zmC),transpose(ztmp2)))
C              call zinv22(ztmp,ztmp) ! ztmp = {1 + gmb W^-1 (z'-C) (W^T)^-1}^-1
CC             call zprmx('xx',2,ztmp,2,2,2)
C              reskmu = matmul(matmul(sqep,matmul(transpose(ztmp2),ztmp)),matmul(gmb,matmul(ztmp2,sqep)))
CC             call zprmx('xx2',2,reskmu,2,2,2)

C             term found empirically : M (Pbet-Palp) M^T where M^T = [W+ (z'-C)^-1 sqrt(e')]
              call zinv22(pbeti,pbeti) ! Replace [P^bet]^-1 with P^bet
              call zinv22(palpi,palpi) ! Replace [P^gam]^-1 with P^gam
              ztmp = matmul(matmul(dcmplx(WT),zmci),sqep)
C             call zprm('WT_zmCi_sqep',2,ztmp,2,2,2)
              reskmu = matmul(transpose(ztmp),matmul(palpi-pbeti,ztmp))
C             call zprm('correction',2,reskmu,2,2,2)

C             beta-independent contribution for third order potential functions
              call zinv22(sqep,ztmp2) ! (ep)^-1/2
              ztmp = -0.5d0*matmul(ztmp2,matmul(epp,ztmp2)) ! -1/2 (ep)^-1/2 epp (ep)^-1/2
C             call zprmx('xx1',2,reskmu,2,2,2)

C             -1/2 pdotdot/pdot = beta-independent term - beta-dependent term
              reskmu = ztmp + reskmu

C              if (l == 2 .and. imu == 3) then
C              call zprmx('-1/2 pdotdot/pdot',2,reskmu,2,2,2)
C              endif

C       ... Make sqrt(Pdot) or Pdot
C           Note: dPdot/dz = P d/dz (1/Pdot) P
C           P^-1 = W+ (z'-C)^-1 W + gamma-alpha   where
C           z' = z + p(z-enu)^3  => dz'/dz = 1+3p(z-enu)^2 = e'
C           Pdot = P d/dz (1/Pdot) P = [P W+ (z'-C)^-1] e' [(z'-C)^-1 W P]
C                = A+ A  where A+ = P^bet [W+ (z'-C)^-1 sqrt(e')]
C           A is proportional to normalization function A = sqrt(2/w)*N
            elseif (modeP == 1 .or. modeP == 3) then
              call zinv22(pbeti,pbeti) ! Replace [P^bet]^-1 with P^bet
              call sqr2x2(ep,sqep)     ! sqrt(1 + 3 zmenu^2 pgam)
              reskmu = matmul(pbeti,matmul(matmul(dcmplx(WT),zmci),sqep))
C             Make Pdot
              if (modeP == 1) reskmu = matmul(reskmu,transpose(reskmu))

C         ... modeP = 4: P^alpha [P^bet]^-1
              elseif (modeP == 4) then
                call zinv22(palpi,ztmp) ; reskmu = matmul(ztmp,pbeti)

C         ... modeP = 5: (beta-alpha) P^alpha [P^bet]^-1
              elseif (modeP == 5) then
C               call zinv22(palpi,ztmp) ; reskmu = matmul(matmul(ztmp,pbeti),dcmplx(bet-alp*unit))
                call zinv22(palpi,ztmp) ; reskmu = matmul(dcmplx(bet-alp*unit),matmul(ztmp,pbeti))

C         ... modeP = 7: (bet-alpha)
              elseif (modeP == 7) then
                reskmu = bet - alp * unit
              endif

C         --- Rotate to lms rep'sn from kappa-mu rep'sn, in m=l,...,-l, ms=1/2,-1/2 order ---
              if (.true.) then
C                call dpzero(reslms,2*size(reslms))
C                do  ms1 = 1, 2
C                  do  ms2 = 1, 2
C                    m1 = int(mu - (ms1-1.5d0))
C                    m2 = int(mu - (ms2-1.5d0))
CC                    if (iabs(m1) <= l .and. iabs(m2) <= l)
CC     .              call info2(1,0,0,'mu %d  m1 m2 ms1 ms2 %4:1,2i',mu,[m1,m2,ms1,ms2])
C                    if (iabs(m1) <= l .and. iabs(m2) <= l) then
C                      reslms(m1,m2,ms1,ms2) =
C     .                  clebsh(1,ms1)*reskmu(1,1)*clebsh(1,ms2) +
C     .                  clebsh(1,ms1)*reskmu(1,2)*clebsh(2,ms2) +
C     .                  clebsh(2,ms1)*reskmu(2,1)*clebsh(1,ms2) +
C     .                  clebsh(2,ms1)*reskmu(2,2)*clebsh(2,ms2)
C                    endif
C                  enddo
C                enddo
                call kmu2lms1(1,l,imu,nl,reskmu,reslms)

C           ... Add contribution from current imu to Pfib
C               Note that Pfib(ms1,ms2) = Pfib(ms2,ms1)^transpose
                lm = l*l + l + 1    ! index to m=0 for this l
                do  ms1 = 1, 2
                  m1 = imu - l - ms1
                  if (iabs(m1) > l) cycle
                  do  ms2 = 1, 2
                    m2 = imu - l - ms2
                    if (iabs(m2) > l) cycle
                    Pfib(lm+m1,ms1,lm+m2,ms2) = reslms(m1,m2,ms1,ms2) ! only for printout
                  enddo
                enddo
              endif

C         ... For lambda matrices, record in mu,lambda,lambda' rep in s_site(ib)%pfr
              if (recs == 0+1) then
                call savpr(icomp,norb,idxkm(1,l,imu),reskmu,s_site(ib)%pfr)
              elseif (recs == -1) then
                call savpr(icomp,norb,idxkm(1,l,imu),reskmu,s_site(ib)%pfra)
              elseif (recs == 2+1) then
                call savpr(icomp,norb,idxkm(1,l,imu),reskmu,s_site(ib)%ddpfr)
              elseif (recs == 3+1) then
                call savpr(icomp,norb,idxkm(1,l,imu),reskmu,s_site(ib)%dpfr)
              endif

            enddo               ! imu loop
          enddo                 ! l loop

          im = lmax+1
          if (morder==1 .or. morder==2) call pmorderx(11,1,1,2,im,i0,0,im**2,im**2,im**2,im**2,im**2,0,Pfib)
          if (savpfp .and. icomp == 1) call pfmat2vec(1,nl,lmr,pfdim,norb,indxsh,Pfib,P)

C         call yprm0('(1p,9e18.10)')
C         call yprmi('Pf (kmu) for ib=%i',ib,0,3,s_site(ib)%pfr,0,norb*2,norb*2,norb*2)

C     ... Printout
C         call yprmi('mkfrfp: Pfib for ib=%i',ib,0,3,Pfib,0,norb*2,norb*2,norb*2)
          if (ipr >= 50) then
            im = -1; if (morder==1 .or. morder==2) im = 1
            ilm = 0
            do  l = 0, lmax
              m1 = l**2 + 1 ; m2 = (l+1)**2
              do  m = -l, l
                ilm = ilm+1
                write(stdo,334) ib,l,m*im,
     .            Pfib(ilm,1,ilm,1),Pfib(min(max(ilm+im,m1),m2),2,ilm,1),
     .            Pfib(min(max(ilm-im,m1),m2),1,ilm,2),Pfib(ilm,2,ilm,2)
  334           format(3i4,8f10.5)
              enddo
            enddo
          endif

        enddo                   ! icomp loop
        deallocate(Pfib)

C       debugging check
C       call mstokm(1,nl,1,norb,Pfib,s_site(ib)%pfr,idxkm)
C       call yprmi('Pfib to lms from kmu',ib,0,3,Pfib,0,norb*2,norb*2,norb*2)

      enddo ! ib loop

C     call zprm('P in mkfrpf',2,p,pfdim,pfdim,4)

      end

      subroutine sqr2x2(mir,mat0)
C- Takes the square root of a complex non-symetric2x2 matrix.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mir   :Input matrxi
Co Outputs
Co   mat0  :mat0*mat0 = mir
Cr Remarks
Cr    Uses diagonalization (sqrt(A)=U^(-1)*sqrt(diag(A))*U)
Cr    Comment: Only a positive definite matrix has one and only square root
C ----------------------------------------------------------------------
      implicit none
      double complex mir(2,2), lamba(2),u(2,2),uinv(2,2)
      double complex a, b, c, d,a1(2), b1(2)
      double precision norma1(2)
      double complex eigen(2,2), mat(2,2), mat0(2,2)
      double complex mat2(2,2), mat2inv(2,2), rez(2,2)
      integer k

C ... Handle diagonal matrix separately
      if(dabs(dble(mir(1,2))) <= 1d-10 .and. dabs(dimag(mir(1,2))) <= 1d-10) then
      if(dabs(dble(mir(2,1))) <= 1d-10 .and. dabs(dimag(mir(2,1))) <= 1d-10) then

        mat0(1,1)=sqrt(mir(1,1))
        mat0(1,2)=0d0
        mat0(2,1)=0d0
        mat0(2,2)=sqrt(mir(2,2))

        goto 401
      endif
      endif

C ... Case of non-diagonal matrix. Eigenvalues and eigenvectors of mir
      a = 1.0d0
      b = mir(2,2)+mir(1,1)
      c = (mir(1,1)*mir(2,2)-mir(2,1)*mir(1,2))
      d = b**2-4d0*a*c
      lamba(1) = (b-sqrt(d))/(2d0*a)
      lamba(2) = (b+sqrt(d))/(2d0*a)

C ... Check for degeneracy
      if(lamba(1) == lamba(2)) then
        b1(1) = dcmplx(1.d0,0.5d0)
        a1(1) = -(mir(1,2)/(mir(1,1)-lamba(1)))*b1(1)
        b1(2) = dcmplx(1.d0,0.5d0)
        a1(2) = ((lamba(2)-mir(2,2))/mir(2,1))*b1(2)

       else

         do  k=1, 2

C      ... Choose one of the infinite set of eigenvectors
           b1(k) = dcmplx(1.d0*k,0.5d0)
           a1(k)=-(mir(1,2)/(mir(1,1)-lamba(k)))*b1(k)

C     ...  Normalization of eigenvectors
           norma1(k)=a1(k)*dconjg(a1(k))+b1(k)*dconjg(b1(k))
           a1(k)=a1(k)/dsqrt(norma1(k))
           b1(k)=b1(k)/dsqrt(norma1(k))

C     ... Test
c         print*, mir(1,1)*a1(k)+mir(1,2)*b1(k)-lamba(k)*a1(k)
c         print*, mir(2,1)*a1(k)+mir(2,2)*b1(k)-lamba(k)*b1(k)
C         print*, a1(k)*dconjg(a1(k))+b1(k)*dconjg(b1(k))

         enddo
       endif

C ...  Orthogonality check
C        ort = dconjg(a1(1))*a1(2) + dconjg(b1(1))*b1(2)

C        print*, 'ort=', ort
C ...  Gram-Schmidt orthogonalization
C        a1(2) = -ort*a1(1) + a1(2)
C        b1(2) = -ort*b1(1) + b1(2)
C        print*, 'ort=', dconjg(a1(1))*a1(2) + dconjg(b1(1))*b1(2)

       u(1,1) = a1(1)
       u(1,2) = a1(2)
       u(2,1) = b1(1)
       u(2,2) = b1(2)

       call zinv22(u,uinv)

       eigen(1,1) = sqrt(lamba(1))
       eigen(1,2) = 0d0
       eigen(2,1) = 0d0
       eigen(2,2) = sqrt(lamba(2))

       call zmpy22(u,eigen,mat)
       call zmpy22(mat,uinv,mat0)

  401  continue
       call zmpy22(mat0,mat0,mat2)
       call zinv22(mat2,mat2inv)
       call zmpy22(mat2inv,mir,rez)

      if(dabs(dble(rez(1,1))-1d0) >= 1d-6 .or.
     .   dabs(dble(rez(2,2))-1d0) >= 1d-6) then
         call rx('square root of matrix cannot be found')
      endif

      end

C      subroutine sqr2x2_newton(mat0,mat)
CC  ...Takes the square root of a complex 2x2 matrix.
CC     Uses modified Newton's algorithm (Mathematics of
CC     Computation, Volume 46, Issue 174, (1986), 537-549
CC     Formulas at page 544 (4.1), (4.2)
C
C      implicit none
C      double complex q0(2,2),mat0(2,2),q(2,2),mat(2,2),qn(2,2),
C     .               pn(2,2),matinv(2,2),qinv(2,2),testc(2,2),
C     .               testc1(2,2),testc2(2,2)
C      double precision rep(2,2),imp(2,2),req(2,2),imq(2,2),
C     .                 repn(2,2),impn(2,2),reqn(2,2),imqn(2,2)
C      integer i,j,nit
C
CC      do i=1, 2
CC       do j=1, 2
CC        print*, mat0(i,j)
CC       enddo
CC      enddo
C
CC      mat0(1,1)=dcmplx(98.8006190037678,-6.85158876278541)
CC      mat0(1,2)=dcmplx(-36.2759661516989,3.57466063430072)
CC      mat0(2,1)=dcmplx(-36.2759661516989,3.57466063430072)
CC      mat0(2,2)=dcmplx(98.8006190037678,-6.85158876278541)
C
C      q0(1,1)=dcmplx(1.0d0,0d0)
C      q0(2,2)=dcmplx(1.0d0,0d0)
C      q0(1,2)=dcmplx(0d0,0d0)
C      q0(2,1)=dcmplx(0d0,0d0)
C
C      do i=1, 2
C         do j=1, 2
C            mat(i,j) = mat0(i,j)
C            q(i,j) = q0(i,j)
C         enddo
C      enddo
C
C      nit=0
C
C 5    do i=1, 2
C         do j=1, 2
C            call zinv22(mat,matinv)
C            call zinv22(q,qinv)
C            pn(i,j)=0.5d0*(mat(i,j)+qinv(i,j))
C            qn(i,j)=0.5d0*(q(i,j)+matinv(i,j))
C         enddo
C      enddo
C
C      do i=1, 2
C         do j=1, 2
C            repn(i,j)=dble(pn(i,j))
C            impn(i,j)=dimag(pn(i,j))
C            reqn(i,j)=dble(qn(i,j))
C            imqn(i,j)=dimag(qn(i,j))
C
C            rep(i,j)=dble(mat(i,j))
C            imp(i,j)=dimag(mat(i,j))
C            req(i,j)=dble(q(i,j))
C            imq(i,j)=dimag(q(i,j))
C         enddo
C      enddo
C
CC      print*, nit
C
C      if(dabs(repn(1,1)-rep(1,1)) - dabs(repn(1,1))*1d-05) 6,6,10
C 6    if(dabs(repn(2,2)-rep(2,2)) - dabs(repn(2,2))*1d-05) 7,7,10
C 7    if(dabs(repn(1,2)-rep(1,2)) - dabs(repn(1,2))*1d-05) 8,8,10
C 8    if(dabs(repn(2,1)-rep(2,1)) - dabs(repn(2,1))*1d-05) 9,9,10
C
C 9    if(dabs(impn(1,1)-imp(1,1)) - dabs(impn(1,1))*1d-05) 11,11,10
C 11   if(dabs(impn(2,2)-imp(2,2)) - dabs(impn(2,2))*1d-05) 12,12,10
C 12   if(dabs(impn(1,2)-imp(1,2)) - dabs(impn(1,2))*1d-05) 13,13,10
C 13   if(dabs(impn(2,1)-imp(2,1)) - dabs(impn(2,1))*1d-05) 14,14,10
C
C
C 14   if(dabs(reqn(1,1)-req(1,1)) - dabs(reqn(1,1))*1d-05) 15,15,10
C 15   if(dabs(reqn(2,2)-req(2,2)) - dabs(reqn(2,2))*1d-05) 16,16,10
C 16   if(dabs(reqn(1,2)-req(1,2)) - dabs(reqn(1,2))*1d-05) 17,17,10
C 17   if(dabs(reqn(2,1)-req(2,1)) - dabs(reqn(2,1))*1d-05) 18,18,10
C
C 18   if(dabs(imqn(1,1)-imq(1,1)) - dabs(imqn(1,1))*1d-05) 19,19,10
C 19   if(dabs(imqn(2,2)-imq(2,2)) - dabs(imqn(2,2))*1d-05) 20,20,10
CC     Temporary fix of numerical absurd instability
C 20   if(nit >= 20) goto 22
C      if(dabs(imqn(1,2)-imq(1,2)) - dabs(imqn(1,2))*1d-05) 21,21,10
C 21   if(dabs(imqn(2,1)-imq(2,1)) - dabs(imqn(2,1))*1d-05) 22,22,10
C
C
C 10   do i=1, 2
C         do j=1, 2
C            mat(i,j)=dcmplx(repn(i,j),impn(i,j))
C            q(i,j)=dcmplx(reqn(i,j),imqn(i,j))
C         enddo
C      enddo
C
C      nit = nit + 1
C
C      if(nit > 200) then
C         call rx('matrix square root cannot be found')
CC         call partofmkfrpf(mat0,mat)
C      endif
C
C      goto 5
C
C 22   continue
C
CC ...Test if it worked
C      call zmpy22(mat,mat,testc1)
C      call zinv22(testc1,testc2)
C      call zmpy22(testc2,mat0,testc)
C
C      if(dabs(dble(testc(1,1))-1.0d0) > 1d-06
C     .     .or. dabs(dble(testc(2,2))-1d0) > 1d-06) then
C      call rx('matrix square root cannot be found')
C      endif
C
CC      print*,'...Start Test'
CC      do i=1, 2
CC         do j=1, 2
CC         print*, testc(i,j)
CC         enddo
CC      enddo
CC      print*,'...End Test'
C
C       do i=1, 2
C         do j=1, 2
C            mat0(i,j)=dcmplx(0d0,0d0)
C            q0(i,j)=dcmplx(0d0,0d0)
C            pn(i,j)=dcmplx(0d0,0d0)
C            qn(i,j)=dcmplx(0d0,0d0)
C         enddo
C      enddo
C
C      end



      subroutine savpr(icomp,norb,idx,fun,pr)
C- Record 2x2 kappa-mu matrix into pr (routine is required to manage dimensioning)
Ci Inputs
Ci   idx   :indices of kappa1 and kappa2 in the kappa-mu rep
Cr Remarks
Cr   fun is a 2x2 matrix with either (kappa,lambda) or (lambda,lambda) indices.
Cr   idx gives the 2 indices where to record it.
      implicit none
      integer icomp,norb,idx(2)
      double complex fun(2,2),pr(norb*2,norb*2,*)

      integer i1,i2

      i1 = idx(1) ; i2 = idx(2)

      pr(i2,i2,icomp) = fun(2,2)
      if (i1 == 0) return     ! case |mu| = l+1/2

      pr(i1,i1,icomp) = fun(1,1)
      pr(i1,i2,icomp) = fun(1,2)
      pr(i2,i1,icomp) = fun(2,1)

      end
