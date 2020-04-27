      subroutine mksopf(s_site,s_pot,nl,nbas,ipc,ipdc,pfdim,indxsh,iopt,pp,sop,socscl,vshft,z,P)
C- Potential functions parameterized from pot pars with spin-orbit coupling
C-----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr ddpfr dpfr
Cio    Passed to:  *
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ipdc  :hash (class control) array for CPA
Ci   iopt  :a compound of digits determining what is generated:
Ci       10000s digit controls whether and how matrix form P is recorded in s_site(:).
Ci             What pot fun is generated depends on modeP (modeP=10s digit iopt)
Ci             As of now output is always recorded in vector form in array P, for non CPA sites.
Ci             This may change in future.
Ci             For some modeP, P is also recorded in matrix form:
Ci           0 s_site(:)%pfr (modeP=0)  s_site(:)%ddpfr (modeP=2)  s_site(:)%dpfr (modeP=3)
Ci           1 P is recorded in s_site(:)%pfra (modeP=0)
Ci         100s digit
Ci           0 make P for alpha representation (bet=alpha)
Ci           1 make P for bare representation (bet=0)
Ci           2 make P for gamma representation (bet=gamma)
Ci           3 make P for gamma representation averaged over spins
Ci          10s digit - modeP determines whether to make P or some derivative:
Ci           0 P <- potential function P
Ci           1 P <- P-dot                                     (not implemented)
Ci           2 P <- -1/2 P-dotdot/P-dot
Ci           3 P <- sqrt(P-dot), choosing abs(sqrt(delta))
Ci           4 P <- P^alpha/P^gamma                           (not implemented)
Ci           5 P <- (gamma-alpha)                             (not implemented)
Ci          1s digit
Ci           0 second-order P
Ci           1 third-order P
Ci           2 first-order P
Ci   pp    :potential parameters, real, by class (atomsr.f)
Ci   sop   :sop(l,is1,is2,i=1..3) : matrix elements between orbitals
Ci         :of spin is1 and is2 for quantum number l.
Ci   socscl:scaling parameters for spin-orbit coupling, by class
Ci   vshft :array of site potential shifts
Ci   z     :complex energy
Cl Local variables
Cl   bet  representation in which P is to be made
Cl   irep 2 rotate to gamma repsn (as coded now, irep=2 maps to irep=3)
Cl        3 rotate to gamma repsn averaged over spins, gamma defined by scal rel pp(5)
Cr Remarks
Cr   This routine adds L.S terms by modifying the potential functions
Cr   Potential functions are stored only by site in s_site(:)%pfr
Cu Updates
Cu   18 Jun 18 Synchronize with updated spherical harmonics
Cu   02 Feb 18 (Belashchenko) Added SO contribution to sqrt(d/dz(z-C))
Cu   05 Jun 16 10000s iopt=1 -> write P to s_site(:)%pfra
Cu   16 Jan 16 (MvS) bug fix, make 10000s digit work
Cu   19 Aug 13 (Belashchenko) Adapted for CPA
Cu   18 Jun 04 (A Chantis) working version
Cu   17 Mar 03 (A Chantis) first created
C-----------------------------------------------------------------------
      use scale_matrix
      use structures
      implicit none
C ... Passed parameters
      integer nl,nbas,ipc(nbas),iopt,pfdim,ipdc(*),indxsh(pfdim)
      integer, parameter :: nsp=2
      double precision pp(6,nl,nsp,*),vshft(nbas),sop(0:nl-1,nsp,nsp,9,*),socscl(*)
      double complex z
      double complex P(pfdim,2,2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_pot)::   s_pot
C ... Dynamically allocated arrays
      integer, allocatable :: ll(:)
      real(8), allocatable :: lz(:,:),lp(:,:),lm(:,:)
      complex(8), allocatable :: lxx(:,:)
      complex(8), allocatable,target, dimension(:,:,:,:) :: zmCi,dzmC,wk2,wk3,pres
      complex(8), allocatable, dimension(:,:) :: wk4
      real(8), allocatable :: wrk(:)
      complex(8), pointer, dimension(:,:,:,:) :: Pfib
C ... Local parameters
      integer l,ic,ic0,ipr,modeP,irep,stdo,norb,ncomp,icomp,recs,i,morder,im
      logical order1,order3,savpfp
      integer, parameter :: nitmax=20
      real(8), parameter:: dCtol=1d-12
      complex(8), parameter :: znul = 0, zone = 1
C ... Dirac-specific
      integer m,m1,m2,is,is2,lmax,ierr,n2,lmr,ilm
      complex(8) sopz(0:nl-1,nsp,nsp),dsopz(0:nl-1,nsp,nsp)
      real(8), dimension(nl,2) :: enu,Ca,Cg,dela,delg,sqdel,pa,gam,gam2,alp,bet,pgam,xwk,ovl !,ovlm
      real(8) cpagam(nl,nsp),xx
C     real(8) sopa(3,2,2)
      complex(8), dimension(2) :: zz,ep,epp
      complex(8), dimension(nl,2) :: ep2,sqep,zwk
C     complex(8) :: z22(2,2), srz22(0:nl-1,nsp,nsp)
      procedure(integer) :: nglob
      procedure(real(8)) :: oalpha,dlength
C ..  For CPA
      integer ib
      logical ldsite

C --- Setup and printout  ---
      call getpr(ipr)
C     ipr    = 60
      stdo   = nglob('stdo')
      irep   = mod(iopt/100,10)
      modeP  = mod(iopt/10,10)
      order1 = mod(iopt,10) == 2
      order3 = mod(iopt,10) == 1
      is = mod(iopt/10000,10)
C     governs where in s_site to store P
      recs = modeP+1; if (modeP == 0 .and. mod(iopt/10000,10) == 1) recs = -1
      savpfp =  .true. ! Always write to P for now

      call lmorder(0,morder,l,m)

C     orderi = mod(iopt,10) == 3
C ... z, z_dot and z_dot_dot (2nd order potential fun)

      if (ipr >= 50) then
        write(stdo,333)
  333   format(' site  l   m',7x,'P(1,1)',14x,'P(2,1)',14x,'P(1,2)',14x,'P(2,2)')
      endif

C --- For each site, make P or related function ---
      lmr = 0
      do  ib = 1, nbas
        lmr = nl*nl*(ib-1)
        ncomp = s_site(ib)%ncomp
        ldsite = ncomp > 1
        norb = s_site(ib)%norb
        allocate(ll(norb),lz(norb,norb),lp(norb,norb),lm(norb,norb),lxx(norb,norb))
        allocate(zmCi(norb,2,norb,2),dzmC(norb,2,norb,2))
        allocate(wk2(norb,2,norb,2))
        allocate(wk3(norb,2,norb,2),pres(norb,2,norb,2))
        call orbmop(norb,ll,lz,lp,lm)
C       call prmx('lp',lp,norb,norb,norb)
C       call prmx('lm',lm,norb,norb,norb)
        if (recs == 0+1) then
          call dpzero(s_site(ib)%pfr,4*norb*norb*ncomp*2)
        elseif (recs == -1) then
          call dpzero(s_site(ib)%pfra,4*norb*norb*ncomp*2)
        elseif (recs == 2+1) then
          call dpzero(s_site(ib)%ddpfr,4*norb*norb*ncomp*2)
        elseif (recs == 3+1) then
          call dpzero(s_site(ib)%dpfr,4*norb*norb*ncomp*2)
        endif

C       Returns gamma ppar, or if a CPA site, CPA-averaged gamma
        ic0 = s_site(ib)%dlmcl         ! Parent class
        call cpagamma(ncomp,ic0,s_pot%dlmwt,nl,nsp,pp,cpagam)

        do  icomp = 1, ncomp
          ic = ipc(ib) ; ic0 = ic     ! class index to this site
          if (ldsite) then
            ic = ipdc(ib) + icomp - 1 ! class index to this component
            ic0 = ipdc(ib)            ! class index to first component
          endif
          lmax = sqrt(dble(norb))-1
          if ((lmax+1)**2 /= norb) call rx('mksopf: improper norb')

C     ... Extract pp(alpha)
          enu = pp(1,:,:,ic) + vshft(ib)
          Ca  = pp(2,:,:,ic) + vshft(ib)
          dela= pp(3,:,:,ic)**2
          pa  = pp(4,:,:,ic)
          gam = pp(5,:,:,ic)
          gam2(:,1) = (pp(5,:,1,ic) + pp(5,:,2,ic))/2
          gam2(:,2) = gam2(:,1)
          alp = pp(6,:,:,ic)

C     ... Make ovl(alpha)
          do  is = 1, 2
            do  l = 1, nl
              ovl(l,is) = oalpha(enu(l,is),Ca(l,is),dela(l,is),alp(l,is),gam(l,is))
            enddo
          enddo

C     ... Make p(gamma)
          xwk = 1 + (Ca-enu)*(gam-alp)/dela
          delg = dela*xwk**2
          Cg = enu + (Ca-enu)*xwk
          bet = alp             ! Default screening repsn is alpha
          if (order3) then
            do  is = 1, 2
              do  l = 1, nl
                pgam(l,is) = pa(l,is) - ovl(l,is)**2
              enddo
            enddo
          elseif (order1) then
            bet = cpagam ! Shift default screening to gamma
            Cg = Ca
            delg = dela
          endif
          sqdel = sqrt(delg)

C         Select screening representation bet if other than alpha
          if (irep == 1) bet = 0
          if (irep == 2) bet = cpagam
          if (irep == 3) then
            bet(:,1) = (cpagam(:,1) + cpagam(:,2))/2
            bet(:,2) = (cpagam(:,1) + cpagam(:,2))/2
          endif

          zmCi = 0; ep2 = 0
          wk3 = 0;
          do  l = 0, lmax   ! In this loop, l designates the true l
            m1 = l**2 + 1 ; m2 = (l+1)**2

            if (order3) then
              zz(:) = z + pgam(l+1,:)*(z - enu(l+1,:))**3
              ep(:) = 1 + 3 * pgam(l+1,:) * (z - enu(l+1,:))**2
              epp(:)  = 6 * pgam(l+1,:) * (z - enu(l+1,:))
              ep2(l+1,:) = -0.5d0 * epp(:)/ep(:)
              sqep(l+1,:) = sqrt(ep)
            else
              zz = z ; ep = 1 ; epp = 0; sqep(l+1,:) = 1
            endif

            do  is = 1, 2
              do  m = m1, m2
                zmCi(m,is,m,is) = zz(is) - Cg(l+1,is)
C               dzmC(m,is,m,is) = 1d0
                wk3(m,is,m,is) = ep(is)
              enddo
            enddo
          enddo ! end l loop

C     ... Subtract the SO Hamiltonian from z-C
          do  is = 1, 2
            do  is2 = 1, 2
              sopz(:,is,is2) =  socscl(ic) * (
     .                           sop(:,is,is2,1,ic)
     .        + (z-enu(:,is2)) * sop(:,is,is2,2,ic)
     .        + (z-enu(:,is))  * sop(:,is2,is,2,ic)
     .        + (z-enu(:,is))  * sop(:,is,is2,3,ic) * (z-enu(:,is2)))
              dsopz(:,is,is2) = socscl(ic) * (
     .          sop(:,is,is2,2,ic) + sop(:,is2,is,2,ic)
     .        + sop(:,is,is2,3,ic) * (z-enu(:,is) + z-enu(:,is2)))

C         ... Include sqrts in Eq. (19), (21) of supplement if order3
c             if (order3) then
c               do  l = 0, lmax
c                 ep(:) = 1 + 3 * pgam(l+1,:) * (z - enu(l+1,:))**2
c                 sopz(l,is,is2) = sopz(l,is,is2)*sqrt(ep(is)*ep(is2))
c                 ep(:) = 1 + pgam(l+1,:) * (z - enu(l+1,:))**2
c                 sopz(l,is,is2) = sopz(l,is,is2)/sqrt(ep(is)*ep(is2))
c               enddo
c             endif
            enddo
          enddo

C     ... zmCi <- z-C - xi(z)*L*S ; dzmC <- [d(zmCi)/dz]^(1/2) (approx)
C         Intermediate A = wk2 = d(zmCi)/dz-1
          lxx = lz ; call zmscal(norb,ll,sopz(:,1,1),lxx)
          zmCi(:,1,:,1) = zmCi(:,1,:,1) - lxx/2 ! spin 1 means (+1/2)
          lxx = lz ; call zmscal(norb,ll,dsopz(:,1,1),lxx)
          wk2(:,1,:,1) =   lxx/2 ! spin 1 means (+1/2)

          lxx = lz ; call zmscal(norb,ll,sopz(:,2,2),lxx)
          zmCi(:,2,:,2) = zmCi(:,2,:,2) + lxx/2 ! spin 2 means (-1/2)
          lxx = lz ; call zmscal(norb,ll,dsopz(:,2,2),lxx)
          wk2(:,2,:,2) = - lxx/2 ! spin 2 means (-1/2)

          lxx = lp ; call zmscal(norb,ll,sopz(:,2,1),lxx)
          zmCi(:,2,:,1) = - lxx/2             ! (1,2) goes with l+
          zmCi(:,1,:,2) = - transpose(lxx)/2  ! (2,1) goes with l-
          lxx = lp ; call zmscal(norb,ll,dsopz(:,2,1),lxx)
          wk2(:,2,:,1) =   lxx/2            ! (1,2) goes with l+
          wk2(:,1,:,2) =   transpose(lxx)/2 ! (2,1) goes with l-

C         call zprm('A',2,wk2,norb*2,norb*2,norb*2)

C     ... Taylor series for [d(zmCi)/dz]^(1/2) : 1 - A/2 - A^2/8 - A^3/16
          n2 = norb*2
C          call dpzero(dzmC,2*n2**2)
C          forall (m = 1:norb, is=1:2) dzmC(m,is,m,is) = 1
C          call daxpy(2*n2**2,-1/2d0,wk2,1,dzmC,1)
C          call zgemm('N','N',n2,n2,n2,(1d0,0d0),wk2,n2,wk2,n2,(0d0,0d0),wk3,n2)
C          call daxpy(2*n2**2,-1/8d0,wk3,1,dzmC,1)
C          call zgemm('N','N',n2,n2,n2,(1d0,0d0),wk2,n2,wk3,n2,(0d0,0d0),pres,n2)
C          call daxpy(2*n2**2,-1/16d0,pres,1,dzmC,1)
C         call zprm('dzmC 3rd order',2,dzmC,n2,n2,n2)

C         Newton Raphson to iterate square root to machine precision.
C         Call B exact dzmC.  Then B^2=1-A.  Iterate [B + B^-1 * (1-A)]/2 -> B
C         If initial B = 1, first iteration yields 1 - A/2, same as Taylor series
c         call dpzero(wk3,2*n2**2)
c         forall (m = 1:norb, is=1:2) wk3(m,is,m,is) = 1
          call dcopy(2*n2**2,wk3,1,dzmC,1) ! Initial estimate for B is unit matrix
          call daxpy(2*n2**2,-1d0,wk2,1,wk3,1) ! 1-A
          allocate(wrk(66*n2))
          do  i = 1, nitmax
            call dcopy(2*n2**2,dzmC,1,wk2,1) ! Copy of estimate for B to wk2
            call zqinv('N',wk2,n2,-66*n2,n2,wrk,n2,ierr)  ! wk2 = B^-1
C           call zprm('B^-1',2,wk2,norb*2,norb*2,norb*2)
            if (ierr /= 0) call rx1('mksopf: zqinv 1,ierr=%i',ierr)
            call zgemm('N','N',n2,n2,n2,(1d0,0d0),wk2,n2,wk3,n2,(0d0,0d0),pres,n2) ! B^-1 (1-A)
C           call zprm('B^-1 (1-A)',2,pres,norb*2,norb*2,norb*2)
            call daxpy(2*n2**2,-1d0,pres,1,dzmC,1) ! delta B
            xx = dlength(2*n2**2,dzmC,1)           ! measure of mean change
            call daxpy(2*n2**2,2d0,pres,1,dzmC,1)  ! B + B^-1(1-A)
            call dscal(2*n2**2,0.5d0,dzmC,1)       ! [B + B^-1(1-A)]/2
C           call zprm('dzmC after NR',2,dzmC,n2,n2,n2)
            if (i > nitmax-3) call info5(10,0,0,' mksopf ib=%i iter=%,2i RMS DQ=%,3;3g',ib,i,xx,4,5)
            if (xx < dctol) exit
          enddo
          if (xx >= dctol) call rx('mksopf failed to converge [d(zmCi)/dz]^(1/2)')
C         Debugging check: Compare B^2 to 1-A
          if (ipr >= 50) then
            call zgemm('N','N',n2,n2,n2,(1d0,0d0),dzmC,n2,dzmC,n2,(0d0,0d0),pres,n2) ! B*B
            call daxpy(2*n2**2,-1d0,wk3,1,pres,1) ! 1-A - B*B
            xx = dlength(2*n2**2,pres,1) ! measure of mean change
            call info5(10,0,0,' mksopf converged |dzmC - sqrt(1-A)| to %;3g in %i iter, site %i',xx,i,ib,4,5)
          endif
C         End of block making [d(zmCi)/dz]^(1/2)

          call dcopy(2*n2**2,zmCi,1,wk2,1) ! wk2 = z-C
          call zqinv('N',wk2,n2,-66*n2,n2,wrk,n2,ierr)
          if (ierr /= 0) call rx1('mksopf: zqinv 1,ierr=%i',ierr)
          zmCi = wk2 ! zmCi = wk2 = (z-C)^-1

C    ...  Make (P^alp) (P^bet)^-1
          if (modeP == 4 .or. modeP == 5) then
C           (P^alp)^-1 = D+ (z-C)^-1 D + gam-alp, where D+ D = delta
            zwk = gam - alp
            call sclmat(7,norb,nl,ll,wk2,dl=sqdel,dr=sqdel,zsh=zwk)  ! (P^alp)^-1
            call zqinv('N',wk2,n2,-66*n2,n2,wrk,n2,ierr) ! (P^alp)
            if (ierr /= 0) call rx1('mksopf: zqinv 2a,ierr=%i',ierr)
C           (P^bet)^-1 = D+ (z-C)^-1 D + gam-bet, where D+ D = delta
            zwk = gam - bet
            call zcopy(size(zmCi),zmCI,1,wk3,1)
            call sclmat(7,norb,nl,ll,wk3,dl=sqdel,dr=sqdel,zsh=zwk)  ! (P^bet)^-1

            allocate(wk4(n2,n2))
!             wk4 = matmul(reshape(wk2,(/n2,n2/)),reshape(wk3,(/n2,n2/)))
!             call zcopy(n2*n2,wk4,1,wk3,1)
! matmul (especially nested) can be flacky
            call zcopy(n2*n2,wk3,1,wk4,1)
            call zgemm('n','n',n2,n2,n2,zone,wk2,n2,wk4,n2,znul,wk3,n2)
            deallocate(wk4,wrk)
C           call zprm('  (P^alp) (P^bet)^-1',2,wk3,n2,n2,n2)
            goto 100 ! Nothing else needed
          endif

C     ... wk2 <- (P^bet) where  (P^bet)^-1 = D (z-C)^-1 D + gamma-bet  and  D D = delta
          zwk = gam - bet
          call sclmat(7,norb,nl,ll,wk2,dl=sqdel,dr=sqdel,zsh=zwk)
          call zqinv('N',wk2,n2,-66*n2,n2,wrk,n2,ierr)
          if (ierr /= 0) call rx1('mksopf: zqinv 2,ierr=%i',ierr)
          deallocate(wrk) ! wk2 = P

C    ...  wk3 <- wk2 * sqdel * zmCi * sqr(e')
          wk3 = zmCi ; call sclmat(1,norb,nl,ll,wk3,dl=sqdel)
c         if (order3) call sclmat(2,norb,nl,ll,wk3,zr=sqep)
          allocate(wk4(n2,n2))
!           wk4 = matmul(matmul(reshape(wk2,(/n2,n2/)),reshape(wk3,(/n2,n2/))),reshape(dzmC,(/n2,n2/)))
!           call zcopy(n2*n2,wk4,1,wk3,1)
! avoid crazy amount of temporaries like above. It is very easy to overwhelm the stack.
          call zgemm('n','n',n2,n2,n2,zone,wk2,n2,wk3,n2,znul,wk4,n2)
          call zgemm('n','n',n2,n2,n2,zone,wk4,n2,dzmC,n2,znul,wk3,n2)
          deallocate(wk4) ! wk3 = mu_R

C     --- Make the potential function or related object ---
  100     continue
C     ... modeP = 0: P = (D+ (z-C)^-1 D + gamma-alpha)^-1
          if (modeP == 0) then
            if (recs == 0+1) call zcopy(norb*norb*4,wk2,1,s_site(ib)%pfr(1,icomp),1)
            if (recs == -1)  call zcopy(norb*norb*4,wk2,1,s_site(ib)%pfra(1,icomp),1)
            Pfib => wk2

C     ... modeP = 1: P_dot = P^bet sqdel (z-C)^-1 e' (z-C)^-1 sqdel P^bet
          elseif (modeP == 1) then
C           pres = [P^bet * sqdel * zmCi * sqr(e')] * [sqr(e') * zmCi * sqdel * P^bet]
            call rx('mksopf: modeP = 1 not implemented')

C     ... modeP = 2: -1/2*(P_dot_dot)*(P_dot)^(-1)
          elseif (modeP == 2) then
C           pres = sqr(e') * gamma /sqdel * [P^bet * sqdel * zmCi * sqr(e')] - e''/(2e')
            pres = wk3
            call sclmat(5,norb,nl,ll,pres,zl=sqep*zwk/sqdel,zsh=ep2)
            if (recs == 2+1) call zcopy(norb*norb*4,pres,1,s_site(ib)%ddpfr(1,icomp),1)
            Pfib => pres

C     ... modeP = 3: sqrt(P-dot)
          elseif (modeP == 3) then
C           pres = [P^bet * sqdel * zmCi * sqr(e')]
            if (recs == 3+1) call zcopy(norb*norb*4,wk3,1,s_site(ib)%dpfr(1,icomp),1)
            Pfib => wk3

C     ... modeP = 4: P^alpha/P^bet, where bet = e.g. gamma or gamma-bar
          elseif (modeP == 4) then
            Pfib => wk3

C     ... modeP = 5: P^alpha/P^bet * (bet-alpha), where bet = e.g. gamma or gamma-bar
C     ... modeP = 5: (bet-alpha) * P^alpha/P^bet seems to be correct with current conventions
          elseif (modeP == 5) then
C           call zprm('  (P^alp) (P^bet)^-1',2,wk3,n2,n2,n2)
            zwk = bet-alp
C           call sclmat(2,norb,nl,ll,wk3,zr=zwk)
            call sclmat(1,norb,nl,ll,wk3,zl=zwk)
C           call zprm('  (bet-alp) (P^alp) (P^bet)^-1',2,wk3,n2,n2,n2)
            Pfib => wk3

C     ... bet-alp
          elseif (modeP == 7) then
            zwk = bet-alp
            call dpzero(wk3,2*size(wk3))
            call sclmat(4,norb,nl,ll,wk3,zsh=zwk)
            Pfib => wk3
          endif

C         Copy to vector form (P) for non-CPA sites.
          if (savpfp .and. icomp == 1) call pfmat2vec(1,nl,lmr,pfdim,norb,indxsh,Pfib,P)

C         call zprm('Pv',2,P(indxsh(lmr+1),1,1),pfdim,2*norb,4)

C     ... Printout
C         call yprmi('Pfib for ib=%i',ib,0,3,Pfib,0,norb*2,norb*2,norb*2)
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
        deallocate(zmCi,dzmC,wk2,wk3,pres,ll,lz,lp,lm,lxx)
C       call zprm('pfr matrix form',2,s_site(ib)%pfr,norb*2,norb*2,norb*2)
      enddo ! ib loop

C     call pmorder(1,1,nbas,1,nl,indxsh,0,pfdim,pfdim,4,0,4,P)
C     call zprm('Pv in mksopf',2,P,pfdim,pfdim,4)

      end

      subroutine pfmat2vec(opt,nl,lmr,pfdim,norb,iprmb,Pfib,Pfv)
C- Copy potential function to/from vector form for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0 do nothing
Ci         :1 copy Pfib to Pv
Ci         :2 copy Pv to Pfib
Ci   nl    :(global maximum l) + 1
Ci   lmr   :offset to first element in iprmb for this site.
Ci         :Usually lmr = nl*nl*(ib-1)
Ci   pfdim :Leading dimension of Pfv
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   norb  :leading dimensoon of Pfib
Cio Inputs/Outputs
Cio  Pfib  :Potential function in matrix form.  Input for opt=1, output for opt=2
Cio  Pfv   :Potential function in vector form.  Input for opt=2, output for opt=1
Cr Remarks
Cr   Pfv has compact vector form.   Not suitable for CPA, which has matrix form
Cr   Relation between matrix elements in vector and matrix form
Cr                     If m=l...-l           If m=-l...l
Cr   Pfv(i,2,1) =   Pfib(ilm-1,2,ilm,1)   Pfib(ilm+1,2,ilm,1)
Cr   Pfv(i,1,2) =   Pfib(ilm+1,1,ilm,2)   Pfib(ilm-1,1,ilm,2)
Cr   Reversal of m ordering merely permutes the rows of Pv.
Cb Bugs
Cb   This should be merged with pokepf
Cu Updates
Cu   30 May 18
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nl,lmr,norb,pfdim,iprmb(pfdim)
      double complex Pfv(pfdim,2,2),Pfib(norb,2,norb,2)
C ... Local parameters
      integer im,morder,l,m1,m2,ilm

      if (opt < 1) return
      if (opt /= 1) call dpzero(Pfib,2*size(Pfib))

      call lmorder(0,morder,[0],[0]); im = -1; if (morder==1 .or. morder==2) im = 1
      do  l = 0, nl-1
        m1 = l**2 + 1 ; m2 = (l+1)**2
        if (opt == 1) then
          do  ilm = m1, m2
            Pfv(iprmb(lmr+ilm),1,1) = Pfib(ilm,1,ilm,1)
            Pfv(iprmb(lmr+ilm),2,2) = Pfib(ilm,2,ilm,2)
            Pfv(iprmb(lmr+ilm),2,1) = Pfib(min(max(ilm+im,m1),m2),2,ilm,1)
            Pfv(iprmb(lmr+ilm),1,2) = Pfib(min(max(ilm-im,m1),m2),1,ilm,2)
          enddo
        else
          do  ilm = m1, m2
            Pfib(ilm,1,ilm,1)                    = Pfv(iprmb(lmr+ilm),1,1)
            Pfib(ilm,2,ilm,2)                    = Pfv(iprmb(lmr+ilm),2,2)
            Pfib(min(max(ilm+im,m1),m2),2,ilm,1) = Pfv(iprmb(lmr+ilm),2,1)
            Pfib(min(max(ilm-im,m1),m2),1,ilm,2) = Pfv(iprmb(lmr+ilm),1,2)
          enddo
        endif
      enddo
      end
