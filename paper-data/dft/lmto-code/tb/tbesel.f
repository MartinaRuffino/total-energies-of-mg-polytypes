      subroutine tbesel(tbc, mixQ,force,pv,s_site,ltb,nbas,nl,nlmq1,
     &                  nlmq,nsp,nclass,lmxl,ipc,dclabl,idxdn,symgr,
     &                  nsgrp,istab,bas,alat,qpol,indxcg,jcg,cg,gaunt,
     &                  ak,strx,dstrx,rho,mmom,drhos,ldim,rhoc,rhon,
     &                  qnu,stni,idu,uh,jh,qmpol,dh,pot0,ecorr,f,fnou,
     &                  fnom,ppdip,tpvq)
C- Make L-expanded electrostatic potential, multipole moments; force
C ----------------------------------------------------------------------
Ci Inputs:   rhoc=c*_RL c_RL' + cc;
Ci           qpol, 10 polarisation parameters;
Ci           rhon=c*_RL S c_RL' + cc (TB+U only);
Ci           drhosl=c*_RL dS_RLR'L'/dR c_R'L' + cc (ovlp only)
Ci           drhos: work array to make drhosl summed over all LL'
Ci           gaunt: Gaunt coefficients
Ci    force: if F, skip calculation of forces (see tbzint.f)
Ci    pv   : if F, skip calculation of pressure (see tbzint.f)
Ci    ldip : 3 include 'spherical' dipole correction to Ewald
Ci         : 2 include 'slab' dipole correction to Ewald
Ci         : any other number - skip the dipole correction
Ci    nlmq1: leading dimension of strx and vm (see Remarks)
Ci    nlmq : global L-cutoff for multipoles, leading dimension of dstrx
Ci           and qmpol, second dimension of strx
Ci    lmxl : species-resolved l-cutoffs for multipoles
Ci    stni : Stoner parameter (l=2)
Ci    idu  : for TB+U, index to which channels get additional potential
Ci    uh,uj: Hubbard U's and J's for TB+U and spin pol TB-L
Ci           uh is the screened U for TB+U. Unscreened U is taken from
Ci           start parameters (qnu)
Ci    rho  : s, p, d Mulliken charges, dimension nl,nsp,nbas (UL only)
Ci    qmpol: multipole moments on each site (from tbmpole)
Ci    ak   : product of Gaunt integrals for TB+U (see makcg9)
Cio Input/Outputs:
Cio   strx, dstrx: structure constants and their radial derivatives (input)
Cio   fstrx, dfstrx: local copy of strx, dstrx. With the --sfly switch
Cio                  both arrays are calculated in the program (see remarks)
Ci    mmom : magnetic moments on each site (from tbmpole)
Co Outputs:
Co    rho  : s, p, d Mulliken charges, dimension nl,nsp,nbas
Co    dh   : increment to on-site hamiltonian, \Delta H_{RLRL'},
Co    pot0 : monopole potential at site R (eq. 7.81, Finnis)
Co           ---all these in Rydberg atomic units.
Co    f    : force from e'static terms
Co    fnou : force from overlap dependence of monopoles (Hubbard part)
Co    fnom : force from overlap dependence of monopoles (Madelung part)
Co    ecorr: second order correction to total energy
Co           E_2 of Finnis' book, eq (7.75); or E^U in TB+U
Co    ppdip: accumulated dipole moment, RELAX=0 in ctrl is misused
Co           to permit the dipole to be accumulated over a subset
Co           of the sites (see tbtote also)
Co    tpvq : contribution of multipoles to 3pV (pressure)
Cl Local variables:
Co    vm   : Madelung potential due to charge transfer
Co    vu   : Hartree U part of the potential
Co    vj   : Stoner J part of the potential
Cr Remarks
Cr    The multipoles Q_L and components of electrostatic potential V_L
Cr    are defined here in terms of the conventions of Stone:
Cr    "Theory of intermolecular forces." If the multipole moments of the
Cr    charge at sites i are Q_L then the components of the potential
Cr    at sites j are defined from
Cr    V(r) = \Sum_L sqrt(4pi/(2l+1)) V_L r^l Y_L(r)
Cr         = \Sum_L sqrt(4pi/(2l+1))(2l+1)!! V_L J_L(r) where J_L are
Cr    Bessel functions as defined by Michael (see subroutine BESSL).
Cr    V_L = \Sum_L' sqrt((2l+1)*(2l'+1))/((2l+1)!!(2l'+1)!!) B_L'L Q_L'
Cr    where B_L'L are Michael's structure constants. Here the Q are
Cr    Stone's, not Jackson's multipole moments. They are related by
Cr    Q_L=sqrt(4pi/(2l+1)) Q_L(Jackson)
Cr    Jackson's are the ones used in the FP-LMTO and NFP programs.
Cr    That means that the conventional multipole moments (Jackson) are
Cr    these moments multiplied by sqrt{(2l+1)/4pi}. The V_L here are
Cr    defined such that the potential V(r) is expanded as
Cr    V(r)=\sum_L \sqrt(4pi/(2l+1)) V_L r^l Y_L(r),
Cr    therefore they must be multiplied by \sqrt{4pi/(2l+1)}
Cr    to get the V_L used in FP-LMTO, in which
Cr    V(r)=\sum_L V_L r^l Y_L(r).
Cr    The extra force from the electrostatics is \sum_L Q_L grad V_L
Cr    Note that in the Finnis 2nd order theory, multipoles are those
Cr    of the difference in charge with respect to free atoms. The old
Cr    MRS theory referred to total electronic and core charge is no
Cr    longer implemented. These theories are essentially equivalent
Cr    anyway.
Cr
Cr    If NOUAVG=F is set in CTRL, then U is averaged over all the
Cr    open channels on that site. Hence the Hubbard potential is
Cr    UdN where dN is the total charge transfer on that site. Otherwise
Cr    the Hubbard potential is V_l = \sum_l U_{ll'} dn_l' and dn_l is
Cr    the charge tansfer in the l-channel. This latter approach is
Cr    closer to the TB+U formula for E^U. In that case non diagonal
Cr    U_{ll'} must be defined. For now, we use U_{ll'}=min(U_l, U_l').
Cr    Other possibilities would be (U_l+U_l')/2 or \sqrt(U_l*U_l')
Cr    The Hubbard energy is then (1/2)\sum_{ll'}U_{ll'}dn_l dn_l'
Cr
Cr    In TB+U we determine on-site Coulomb integrals indirectly using
Cr    I=(U+2lJ)/(2l+1) in the d channels (l=2) and U is taken from
Cr    qnu and I is taken from the I= token in SPEC. In the s and p
Cr    channels we take U from qnu and set J=0. We then use F^0=U
Cr    and J=(5I-U)/4=(F^2+F^4)/14 if l=2. This furnishes us with the
Cr    Slater parameters F^0=U (all l); F^2=14J/2.6, F^4=1.6F^2 if l=2
Cr    F^2=F^4=0 if l<2.
Cr
Cr    The structure constants strx and dstrx may be made within tbesel
Cr    on the fly for each connecting vector (switch --sfly) or they
Cr    may be passed in as a lookup table (default). The latter increases
Cr    speed at the expense of memory.
Cr
Cr    L-dimensions of structure constants strx, nlmq1 and nlmq,
Cr    correspond to (l+2)^2 and (l+1)^2, respectively, where L is the
Cr    l-cutoff for multipoles, or equivalently the largest l for which
Cr    Delta_l'l"l\=0. For Hamiltonian or total energy calculations
Cr    nlmq suffices, but for forces We need to go one angular momentum up,
Cr    hence nlmq1. No such complication arises for dstrx though.
Cr
Cr    Meaning of input switches force and pv:
Cr    parameters force and pv do not necesserily coincide with those
Cr    stored in variable ltb (in 2^4 and 2^7 registers, respectively).
Cr    The latter define whether the force or pressure calculations are
Cr    required in principle (and hence affect dimension nlmq1) whereas
Cr    the former tell whether these should be done during current call
Cr    of tbesel. Both force and pv are set to .False. in tbzint.f if
Cr    self-consistency has not been reached yet.
Cr
Cu Upgrades
Cu    As of version 8.0, tbe uses (\delta q)*U rather than q*U as the
Cu    'Hubbard' potential. \delta q = q - q_0 and q_0 is computed in
Cu    each channel from qnu. The Hubbard U may now be the average
Cu    over all channels. The older theory (as in the MRS paper) is
Cu    no longer implemented.
Cu    TB+U retains all terms except the Hubbard potential which is
Cu    modified to be spin (s) and orbital dependent:
Cu    V^{s}_{LL'}=\sum_{L''L'''} V_{LL''L'''L'} dn^{-s}_{L''L'''} +
Cu                    (V_{LL''L'''L'}-V_{LL''L'L'''}) dn^{s}_{L''L'''}
Cu    n is the spin polarised density matrix (see tbfrce) and dn is
Cu    the difference given by subtracting the diagonal elements q_0,
Cu    the V are given in terms of the Slater radial integrals as
Cu    V_{LL'L''L'''}=\sum_{k=0}^{2lmax} R^k(ll'l''l''') A^k(LL'L''L''')
Cu    (k an even number less that 2l) with
Cu    A^k(LL'L''L''')=4pi/2k+1 \sum_{p=-k}^{k} C_{LL'''K} C_{L'L''K}
Cu    where K={kp} as L={lm} and C are the Gaunt (CG) coefficients.
Cu    As a first approximation to the Slater integrals we will have
Cu    R^0(llll) = F^0 = U, in all channels
Cu    R^2(2222) = F^2, R^4(2222) = F^4 (ie l=2, d-channel)
Cu    all others zero.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu    3 Mar 10 (SL)  L-summation limits, parameters force and pv
Cu   19 Jan 10 (ATP) Option for strux from lookup table
Cu   04 Jun 08 (ATP) Cleanup, multipole contribution to pressure
Cu   17 Mar 06 (ATP) Cleanup, strip out old MRS theory, start on
Cu                   non orthogonal TB-L
Cu    6 Jun 05 (ATP) Added TB+U
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   03 Dec 01 (ATP) bug fix
C ----------------------------------------------------------------------
      use tbprl, only : tbc_t
      use mpi

      use structures
      implicit none

      type(tbc_t), intent(in) :: tbc
      logical mixQ,force,pv
      integer ltb,nbas,nl,nlmq1,nlmq,nsp,nclass,nsgrp,istab(*),
     .        idxdn(0:nl-1,*),ldim
      integer lmxl(nclass),ipc(*),indxcg(*),jcg(*),idu(4,nbas)
      double precision ecorr
      double precision bas(3,*),symgr(*),qpol(10,nclass), !z(*),
     .  rho(nl,2,nbas),rhoc(nl**2,nl**2,nbas),alat,mmom(nbas),
     .  qnu(3,0:nl-1,nsp,nclass),stni(nclass),uh(4,nbas),jh(4,nbas),
     .  dclabl(nclass),qmpol(nlmq,nbas),cg(*),gaunt(9,9,25),
     .  dh(nl**2,nl**2,nbas,nsp),
     .  pot0(nbas),rhon(nl**2,nl**2,nbas,2),tpvq,
     .  No(3,nbas,nsp),f(3,nbas),fnou(3,nbas),fnom(3,nbas),
     .  Ak(9,9,9,9,3),ppdip(3),drhos(nbas,nbas,3)
      real(8), intent(in), target :: strx(nlmq, nlmq1, nbas, *),
     &                              dstrx(nlmq, nlmq , nbas, *)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)

      logical lov,MOL,bittst,UL,Uav,TBU,diagn,point,field,sfly,cmdopt
      integer i,k,ic,jc,ib,jb,l,lp,il,ilm,jlm,ilmp,ilmpp,ip,
     .        i1mach,iprint,
     .        nlm,ilm1,ilm2,ilm3,ilm4,
     .        isp,l1,l2,l3,l4,ifrlx(3),opi,opj,nlmi,nlmj,
     .        ilmi,ilmj,it(3),li,lj,li1,nlmi1
      integer ll,getnlm,parg,ione0
      double precision M,dq(0:2),vprt(36),dhd(25),qmp(25),
     .                 sumV,sumU,sumJ,tau(3),taua(3),
     .                 hl(100),bl(100),U,J,F2,F4,
     .                 dQtot,drho1,drho2,drho3,q0,
     .                 VL(9,9,9,9),Rk(9,9,9,9,3),Efield(3)
      double precision pi,dsum,dsqrt,ddot,getavU,getavJ,d1mach,dmin1
C     double precision dasum,efg(5)
      character*8 clabl
      character*120 strn
C...  Automatic arrays
      integer ione(nlmq1)
      double precision Ratm(nbas,3)
      double precision A(nbas),B(nbas)
      double precision vu(0:nl-1,nbas),vj(2,nbas)
      real(8), allocatable :: vm(:,:) !, drhos(:,:,:)
      real(8), pointer :: fstrx(:,:) =>null() ,fdstrx(:,:) => null()

      integer :: comm, vm_mt, err, pid, nproc
      logical :: lroot


      call tcn('tbesel')

      pi = 4d0*datan(1d0)

      comm   = tbc % c3d % comm
      lroot  = tbc % c3d % lrt
      pid    = tbc % c3d % id
      nproc  = tbc % c3d % sz



      no(1:3,1:nbas,1:nsp) = 0d0
      ppdip(1:3) = 0d0
      a(1:nbas) = 0d0
      vj(1:2,1:nbas) = 0d0
      lov = bittst(ltb,1)
      efield(1:3) = 0d0

      call info0(30,1,0,' TBESEL (Stone''s definitions for Q and V)')
      field = .false.

C ... get external electric field from command line (could put in ctrl)
      if (cmdopt('-efield=',8,0,strn)) then
        ip = 8
        call skipbl(strn,len(strn),ip)
        i = parg(' ',4,strn,ip,len(strn),', ',2,3,it,efield)
        if (i < 0) then
          call rxs2('TBESEL: failed to parse "',strn(1:ip+5),' ..."')
        else
          field = .true.
          if (iprint() > 20) then
            call awrit1(' TBESEL external electric field E= %3:-2,4;4d',
     .        ' ',120,i1mach(2),efield)
          endif
        endif
      endif
C ... Molecule or cluster
      MOL =  bittst(ltb,2**18)
C ... 2nd order self consistent TB (TB-L)
      UL = bittst(ltb,2**15)
C ... use U averaged over each channel at a site
      Uav = (.not. bittst(ltb,2**14))
      call rxx(UL .and. .not. Uav,
     .  ' TBESEL: set NOUAVG=F in ctrl.'//
     .  ' For an orbital dependent potential, use TB+U')
C ... TB+U:
      TBU = bittst(ltb,2**13)
C ... ignore non diagonal density matrix element (mostly for debugging)
      diagn = cmdopt('--diagn',7,0,strn)
C ... keep point charges only ---
      point = bittst(ltb,2**9)
      if (.not. point) point = cmdopt('--point',7,0,strn)
      if (point) call info0(10,0,0,' TBESEL: using point charges only')

C ... make strux on the fly to save memory
      sfly = cmdopt('--sfly',6,0,strn)

C --- sanity test ---
      if (TBU .and. UL) call rx('TBESEL: cannot have UL and TB+U')

C --- loop to get dipole moments on selected atoms (misuse RELAX=) ---
      do  ib = 1, nbas
        ic = ipc(ib)
        nlm = (lmxl(ic)+1)**2
        ifrlx = s_site(ib)%relax
        if (ifrlx(1) == 1 .and. nlm >= 4) then
          do  ilm = 2, 4
            ppdip(ilm-1) = ppdip(ilm-1) + qmpol(ilm,ib)
          enddo
        endif
      enddo

      call tcn('Madelung potential')
C --- get components of the Madelung potential ---

!       allocate(vm(nlmq1, tbc%esamap(pid)+1:tbc%esamap(pid+1)))
      allocate(vm(nlmq1, nbas))
      vm = 0d0

c...  preset ione(L) = (-1)^l
      ione0 = -1
      do il = 0, ll(nlmq1)
        ione0 = -ione0
        ione(il*il+1:(il+1)**2) = ione0
      enddo

C --- loop to make Madelung potential ---
      tpvq = 0d0
      do  ib = tbc%esamap(pid)+1, tbc%esamap(pid+1)
        ic = ipc(ib)
        li = lmxl(ic)
        nlmi  = (li+1)*(li+1)
        if (force .or. pv) then
          li1 = li + 1
          nlmi1  = (li1+1)*(li1+1)
        else
          li1 = li
          nlmi1  = nlmi
        endif

        do  jb = 1, nbas
          jc = ipc(jb)
          lj = lmxl(jc)
          nlmj = (lj+1)*(lj+1)

          fstrx  => strx (1:nlmj, 1:nlmi1, jb, ib-tbc%esamap(pid))
          do  ilm = 1, nlmi1
            vm(ilm,ib) = vm(ilm,ib)
     &                    + 2d0*sum(fstrx(1:nlmj,ilm)*qmpol(1:nlmj,jb))
          end do

          if (pv) then
            fdstrx => dstrx(1:nlmj, 1:nlmi , jb, ib-tbc%esamap(pid))
            do  ilm = 1, nlmi
               tpvq = tpvq - qmpol(ilm,ib)
     &                    * sum(fdstrx(1:nlmj,ilm)*qmpol(1:nlmj,jb))
            enddo
          endif

          if (lov) A(ib) = A(ib) + 2d0*fstrx(1,1)*qmpol(1,jb)

        enddo
C .. add electric field
        if (field) then
          vm(1,ib) = vm(1,ib)
     .             + sqrt(2d0) * alat * sum(efield(1:3)*bas(1:3,ib))
          if (nlmi1 >= 4) vm(2:4,ib) = vm(2:4,ib)
     .                                     + sqrt(2d0) * efield(1:3)
        endif
      enddo ! ib loop
      if (pv) call mpi_allreduce(mpi_in_place,tpvq,1,mpi_real8,
     &                                    mpi_sum,comm,err)

!       call vblshare( nlmq1*nbas, vm)

! This is simpler than the following but 2-3x slower
!       call mpi_allreduce(mpi_in_place, vm, nlmq1*nbas, mpi_real8,
!    &                                    mpi_sum, comm, err)


      call mpi_type_contiguous(nlmq1, mpi_real8, vm_mt, err)
      call mpi_type_commit(vm_mt, err)
      call mpi_allgatherv(mpi_in_place, 0, vm_mt, vm, tbc%escount,
     &                 tbc%esamap(0:nproc-1), vm_mt, comm, err)

      call mpi_type_free(vm_mt, err)

      if (lov) then
        if (.not. Uav) call rx('TBESEL: not ready for OVLP and U_l')
!         call vblshare( nbas, a)
        call mpi_allgatherv(mpi_in_place,0,mpi_real8,a,tbc%escount,
     &             tbc%esamap(0:nproc-1), mpi_real8, comm, err)

        do  ib = 1, nbas
          ic = ipc(ib)
          B(ib) = getavU(nl, nsp, qnu, idxdn, ic) * qmpol(1, ib)
          pot0(ib) = A(ib) + B(ib)
        enddo
      endif


      call tcx('Madelung potential')


      call tcn('U')
      nlm = nl**2
      dh(1:nlm,1:nlm,1:nbas,1:nsp) = 0d0
      Ratm(1:nbas,1:3) = 0d0
C ... "U * dN^2" contribution to E_2 or E^U:
      sumU = 0d0
      sumJ = 0d0
C --- get components of the Hubbard potential and double counting ---
      if (TBU) then
        do  ib = 1, nbas
          ic = ipc(ib)
C --- assemble matrix of Slater integrals ---
          call dcopy(19683,0d0,0,Rk,1)
          do  ilm1 = 1, nlm
            do  ilm2 = 1, nlm
              do  ilm3 = 1, nlm
                do  ilm4 = 1, nlm
                  l1 = ll(ilm1)
                  l2 = ll(ilm2)
                  l3 = ll(ilm3)
                  l4 = ll(ilm4)
C --- first approximation: R^k(ll'l''l''')=R^k(llll)delta(ll'l''l''')
                  if (l1 == l2 .and. l1 == l3 .and. l1 == l4) then
                    l = l1
                    U = qnu(3,l,1,ic)
                    Rk(ilm1,ilm2,ilm3,ilm4,1) = U
                    if (l == 2) then
                      F2 = (14d0/(4d0*2.6d0))*(5d0*stni(ic) - U)
                      call rxx(F2 < 0,
     .                         ' TBESEL: U, I parameter error. F2<0')
                      F4 = 1.6d0*F2
                      Rk(ilm1,ilm2,ilm3,ilm4,2) = F2
                      Rk(ilm1,ilm2,ilm3,ilm4,3) = F4
                      do  i = 1, 3
                        Ratm(ib,i) =  Rk(ilm1,ilm2,ilm3,ilm4,i)
                      enddo
                    else
                      Ratm(ib,l+1) = U
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo
C --- assemble matrix of interaction parameters V = R^k * A^k ---
          if (iprint() >= 60) write (*,10)
          call dcopy(6561,0d0,0,VL,1)
          do  ilm1 = 1, nlm
            do  ilm2 = 1, nlm
              do  ilm3 = 1, nlm
                do  ilm4 = 1, nlm
                  if ((idxdn(ll(ilm1),ic) == 1) .and.
     .                (idxdn(ll(ilm2),ic) == 1) .and.
     .                (idxdn(ll(ilm3),ic) == 1) .and.
     .                (idxdn(ll(ilm4),ic) == 1)) then
                    do  i = 1, 3
                     VL(ilm1,ilm2,ilm3,ilm4) = VL(ilm1,ilm2,ilm3,ilm4) +
     .               Rk(ilm1,ilm2,ilm3,ilm4,i)*Ak(ilm1,ilm2,ilm3,ilm4,i)
                    enddo

C --- verbose output ---
                    if (iprint() >= 60 .and.
     .                Rk(ilm1,ilm2,ilm3,ilm4,1) > d1mach(3))
     .                write(*,20)
     .                ilm1,ilm2,ilm3,ilm4,ll(ilm1),ll(ilm2),ll(ilm3),
     .                ll(ilm4),(Rk(ilm1,ilm2,ilm3,ilm4,k),k=1,3),
     .                (Ak(ilm1,ilm2,ilm3,ilm4,k),k=1,3),
     .                VL(ilm1,ilm2,ilm3,ilm4)
C ---------------------

                  endif
                enddo
              enddo
            enddo
          enddo
C --- put the on-site potential directly into dh ---
          if (iprint() > 40) write (*,50)
          do  isp = 1, nsp
            do  ilm1 = 1, nlm
              do  ilm2 = 1, nlm
                do  ilm3 = 1, nlm
                  do  ilm4 = 1, nlm
                    if ((idxdn(ll(ilm1),ic) == 1) .and.
     .                  (idxdn(ll(ilm2),ic) == 1) .and.
     .                  (idxdn(ll(ilm3),ic) == 1) .and.
     .                  (idxdn(ll(ilm4),ic) == 1)) then
                      if (diagn) then
                        if (ilm1 /= ilm2 .or. ilm3 /= ilm4) cycle
                      endif
                      U = VL(ilm1,ilm3,ilm4,ilm2)
                      J = VL(ilm1,ilm3,ilm2,ilm4)
                      drho1 = rhon(ilm3,ilm4,ib,2-isp/2)
                      drho2 = rhon(ilm3,ilm4,ib,isp)
                      if (ilm3 == ilm4) then
                        l = ll(ilm3)
                        q0 = qnu(1,l,isp,ic) / (2*l + 1)
                        drho1 = drho1 - q0
                        drho2 = drho2 - q0
                      endif
                      drho3 = rhon(ilm1,ilm2,ib,isp)
                      if (ilm1 == ilm2) then
                        l = ll(ilm1)
                        q0 = qnu(1,l,isp,ic) / (2*l + 1)
                        drho3 = drho3 - q0
                      endif
                      dh(ilm1,ilm2,ib,isp) = dh(ilm1,ilm2,ib,isp) +
     .                                       U*drho1 + (U-J)*drho2
                      sumU = sumU + (U*drho1 + (U-J)*drho2) * drho3

C --- verbose output ---
                      if ((iprint() >= 50 .and.
     .                (dabs(U)+dabs(J)) > d1mach(3)) .or.
     .                (iprint() > 40
     .                .and. ilm3 == ilm4 .and. ilm1 == ilm2))
     .                write (*,75)
     .                ilm1,ilm2,ilm3,ilm4,ll(ilm1),ll(ilm2),ll(ilm3),
     .                ll(ilm4),isp,drho1,drho2,U,J,U-J
C ----------------------

                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      else
C --- TB-L branch ---
        vu(0:nl-1,1:nbas) = 0d0
        do  ib = 1, nbas
          ic = ipc(ib)
          if (Uav) then
            U = getavU(nl,nsp,qnu,idxdn,ic)
            dQtot = qmpol(1,ib)
            sumU = sumU +  U * dQtot*dQtot
          endif
          do  l = 0, nl-1
            if (Uav) then
C ... no real need for this:
C              if (idxdn(l,ic) == 1) then
                vu(l,ib) = U * dQtot
C              endif
            else
              do  lp = 0, nl-1
                if (lp == l) then
                  U = qnu(3,l,1,ic)
                else
                  U = dmin1(qnu(3,l,1,ic),qnu(3,lp,1,ic))
                endif
                vu(l,ib) = vu(l,ib)
     .                      + U * (rho(lp+1,1,ib) - qnu(1,lp,1,ic))
                sumU = sumU + U * ((rho(l+1,1,ib) - qnu(1,l,1,ic)))
     .                          * ((rho(lp+1,1,ib) - qnu(1,lp,1,ic)))
              enddo
            endif
          enddo
C --- Spin polarised TB-L ---
          if (nsp == 2) then
            J = getavJ(nl,jh(1,ib),idxdn,ic)
            vj(1,ib) = - 0.5d0 * J * (qmpol(1,ib) + mmom(ib))
            vj(2,ib) = - 0.5d0 * J * (qmpol(1,ib) - mmom(ib))
            sumJ = sumJ - J * (0.25d0 * (qmpol(1,ib) + mmom(ib))**2
     .                      +  0.25d0 * (qmpol(1,ib) - mmom(ib))**2)
          endif
        enddo
      endif
      call tcx('U')

C --- increments to hamiltonian ---
      call tcn('delta H')
      do   ib = 1, nbas
        ic = ipc(ib)
        li = lmxl(ic)
        nlm  = (li+1)**2
        do  isp = 1, nsp
          if (iprint() > 50 .and. isp == 1) then
            if (UL) then
              write (*,150)
            else
              write (*,175)
            endif
          endif
          do  ilmp = 1, nl**2
            do  ilmpp = 1, nl**2
              if ((idxdn(ll(ilmp) ,ic) == 1) .and.
     .            (idxdn(ll(ilmpp),ic) == 1)) then
                if (.not. TBU) then
                  if (ilmp == ilmpp) then
                    dh(ilmp,ilmpp,ib,isp) = dh(ilmp,ilmpp,ib,isp) +
     .                vu(ll(ilmp),ib) + vj(isp,ib)
                  endif
                endif
                do  ilm = 1, nlm
                  call getM(ilm,ilmp,ilmpp,qpol(1,ic),M)
                  dh(ilmp,ilmpp,ib,isp) = dh(ilmp,ilmpp,ib,isp) +
     .                            vm(ilm,ib) * M * gaunt(ilmp,ilmpp,ilm)

C --- verbose output ---
                  if (iprint() > 50 .and. M /= 0d0 .and. isp == 1
     .                .and. gaunt(ilmp,ilmpp,ilm) /= 0d0) then
                    if (ilmp == ilmpp .and. UL) then
                      U = getavU(nl,nsp,qnu,idxdn,ic)
                      dQtot = dsum(nl,rho(1,1,ib),1)
     .                      - dsum(nl,qnu(1,0,1,ic),3)
                      write (*,250)
     .                  ilmp,ilmpp,ilm,ll(ilmp),ll(ilmpp),ll(ilm),
     .                  M,gaunt(ilmp,ilmpp,ilm),vm(ilm,ib),
     .                  dQtot,U,vu(ll(ilmp),ib),dh(ilmp,ilmpp,ib,isp)
                    else
                      write (*,210)
     .                  ilmp,ilmpp,ilm,ll(ilmp),ll(ilmpp),ll(ilm),
     .                  M,gaunt(ilmp,ilmpp,ilm),vm(ilm,ib),
     .                  dh(ilmp,ilmpp,ib,isp)
                    endif
                  endif
C ---------------------

                enddo
              endif
            enddo
          enddo
        enddo
      enddo
      call tcx('delta H')

C --- electrostatic energy and force ---
      call tcn('es forces')
C ... dQ * V :
      sumV = 0d0
!       write(1000+pid,*) 'newrun', tbc%esamap
      do  ib = 1, nbas
        ic = ipc(ib)
        nlmi = (lmxl(ic)+1)**2
        do  ilm = 1, nlmi
!           write(1000+pid,*)'qmpol,vm:',ib,ilm,qmpol(ilm,ib),vm(ilm,ib)
          sumV = sumV + qmpol(ilm,ib) * vm(ilm,ib)
        enddo
      enddo
      ecorr = 0.5d0*(sumV + sumU + sumJ)

      if (force) then
C ... Calculate es forces
c       call tbfrc2(nbas,nlmq1,nlmq,ipc,lmxl,vm,qmpol,indxcg,jcg,cg,f)
        call tbfrc3(nbas,nlmq1,nlmq,ipc,lmxl,vm,qmpol,f)

C ---   make terms in dervatives of the overlap for forces ---
        if (lov) then
!           allocate(drhos(nbas,nbas,3))
          call dpzero(fnou,3*nbas)
          call dpzero(fnom,3*nbas)
!           call dcopy(3*nbas**2,0d0,0,drhos,1)
!           opi = 0
!           do  ib = 1, nbas
!             ic = ipc(ib)
!             nlmi = getnlm(nl,idxdn,ic)
!             opj = 0
!             do  jb = 1, nbas
!               jc = ipc(jb)
!               nlmj = getnlm(nl,idxdn,jc)
!               do  k = 1, 3
!                 drhos(ib,jb,k) = 0d0
!                 do  ilmi = 1, nlmi
!                   do  ilmj = 1, nlmj
!                     drhos(ib,jb,k) = drhos(ib,jb,k)
!      .                + drhosl(opi+ilmi,opj+ilmj,k)
!                   enddo
!                 enddo
!               enddo
!               opj = opj + nlmj
!             enddo
!             opi = opi + nlmi
!           enddo
!           call rxx(opi /= ldim,' bug in tbesel: opi /= ldim')
!           call rxx(opj /= ldim,' bug in tbesel: opj /= ldim')
          if (iprint() > 40) then
            print *,   ' TBESEL: drhos:'
            do  k = 1, 3
              if (k. eq. 1) print *, '         x :'
              if (k. eq. 2) print *, '         y :'
              if (k. eq. 3) print *, '         z :'
              do  ib = 1, nbas
                write(*,400) (drhos(ib,jb,k),jb=1,nbas)
              enddo
            enddo
          endif

C ---     forces from overlap ---
          do  ib = 1, nbas
            do  jb = 1, nbas
              do  k = 1, 3
                fnou(k,ib) = fnou(k,ib)
     .                     - (B(ib) + B(jb)) * drhos(ib,jb,k)
                fnom(k,ib) = fnom(k,ib) - (vm(1,ib) + vm(1,jb))
     .                                  * drhos(ib,jb,k)
              enddo
            enddo
          enddo
!           deallocate(drhos)
        endif
        if (.not. MOL) then
          call symfor(nbas,1,symgr,nsgrp,istab,f)
          if (lov) then
            call symfor(nbas,1,symgr,nsgrp,istab,fnou)
            call symfor(nbas,1,symgr,nsgrp,istab,fnom)
          endif
        endif

      endif                   ! end of force if-loop
      call tcx('es forces')

C --- printout ---
      if (iprint() < 30) then
        call tcx('tbesel')
        return
      endif
      if (iprint() == 30) goto 1000
      do ib = 1, nbas
        ic = ipc(ib)
        li = lmxl(ic)
        nlmi  = (li+1)**2
        nlmi1  = nlmi
        if (force .or. pv) nlmi1  = (li+2)**2
        call r8tos8(dclabl(ic),clabl)
        print *
        if (mixQ) then
          call awrit0('Atom '//clabl,' ',180,i1mach(2))
        else
          if (nsp == 2) then
            if (nl == 3) then
              call awrit8('Atom '//clabl//
     .        '%cN_s=%d+%d N_p=%d+%d N_d=%d+%d Q=%d mom=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(1,2,ib),rho(2,1,ib),
     .        rho(2,2,ib),rho(3,1,ib),rho(3,2,ib),
     .        dsum(nl,rho(1,1,ib),1)+dsum(nl,rho(1,2,ib),1),
     .        dsum(nl,rho(1,1,ib),1)-dsum(nl,rho(1,2,ib),1))
            elseif (nl == 2) then
              call awrit6('Atom '//clabl//
     .        '%cN_s=%d+%d N_p=%d+%d Q=%d mom=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(1,2,ib),rho(2,1,ib),
     .        rho(2,2,ib),
     .        dsum(nl,rho(1,1,ib),1)+dsum(nl,rho(1,2,ib),1),
     .        dsum(nl,rho(1,1,ib),1)-dsum(nl,rho(1,2,ib),1))
            elseif (nl == 1) then
              call awrit4('Atom '//clabl//
     .        '%cN_s=%d+%d Q=%d mom=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(1,2,ib),
     .        dsum(nl,rho(1,1,ib),1)+dsum(nl,rho(1,2,ib),1),
     .        dsum(nl,rho(1,1,ib),1)-dsum(nl,rho(1,2,ib),1))
            endif
          else
            if (nl == 3) then
              call awrit4('Atom '//clabl//
     .        '%cN_s=%d N_p=%d N_d=%d Q=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(2,1,ib),
     .        rho(3,1,ib),dsum(nl,rho(1,1,ib),1))
            elseif (nl == 2) then
              call awrit3('Atom '//clabl//
     .        '%cN_s=%d N_p=%d Q=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(2,1,ib),
     .        dsum(nl,rho(1,1,ib),1))
            elseif (nl == 1) then
              call awrit2('Atom '//clabl//
     .        '%cN_s=%d Q=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),
     .        dsum(nl,rho(1,1,ib),1))
            endif
          endif
        endif
        if (iprint() > 31) then
          if (TBU) then
            print *, 'c*_RL S c_RL'':'
            print *, ' spin up:'
            do  i = 1, nl**2
              write (*,300) (rhon(i,k,ib,1),k=1,nl**2)
            enddo
            print *, ' spin down:'
            do  i = 1, nl**2
              write (*,300) (rhon(i,k,ib,2),k=1,nl**2)
            enddo
          endif
          print *, 'c*_RL c_RL'':'
          do  i = 1, nl**2
            write (*,300) (rhoc(i,k,ib),k=1,nl**2)
          enddo
        endif
        if (TBU) then
          if (nl == 3) then
            call awrit5('        U=%d, F^2=%d, F^4=%d, J=%d, I=%d',' ',
     .                  180,i1mach(2),Ratm(ib,1),Ratm(ib,2),Ratm(ib,3),
     .                  (Ratm(ib,2)+Ratm(ib,3))/14d0,stni(ic))
          else
            call awrit2('        U_s=%d, U_p=%d, F^2=0',' ',180,
     .                  i1mach(2),Ratm(ib,1),Ratm(ib,2))
          endif
        endif
        call awrit7('        M=%d %d %d %d %d %d %d',' ',180,i1mach(2),
     .   sqrt(4*pi/3)*qpol(1,ic),sqrt(4*pi/5)*qpol(2,ic),
     .   sqrt(4*pi/5)*qpol(3,ic),sqrt(4*pi/3)*qpol(4,ic),
     .   sqrt(4*pi/5)*qpol(5,ic),sqrt(4*pi/7)*qpol(6,ic),
     .   sqrt(4*pi/9)*qpol(7,ic))
        qmp = 0d0
        qmp(1:nlmi) = qmpol(1:nlmi,ib)
        call awrit3('        Q^e/e=%d, %3:1d, %5:1d',' ',180,
     .                i1mach(2),qmp(1),qmp(2),qmp(5))
        qmp(1) = qmp(1)/dsqrt(4d0*pi)
        if (nlmi >= 2) then
          qmp(2:4) = qmp(2:4)*dsqrt(3d0/(4d0*pi))
          if (nlmi >= 5)
     .      qmp(5:9) = qmp(5:9)*dsqrt(5d0/(4d0*pi))
        endif
        call awrit3('        Qmpol=%d, %3:1d, %5:1d',' ',180,
     .                i1mach(2),qmp(1),qmp(2),qmp(5))
        if (nlmi >= 10) then
          call dpcopy(qmpol(10,ib),qmp(10),1,7,dsqrt(7d0/(4d0*pi)))
          call awrit1('              %7:1d',' ',180,i1mach(2),qmp(10))
          if (nlmi >= 17) then
            call dpcopy(qmpol(17,ib),qmp(17),1,9,dsqrt(9d0/(4d0*pi)))
            call awrit1('              %9:1d',' ',180,i1mach(2),qmp(17))
          endif
        endif
        vprt = 0d0
        call dcopy(nlmi1,vm(1,ib),1,vprt,1)
        call awrit3('        e dV=  %d, %3:1d, %5:1d',' ',180,
     .              i1mach(2),vprt,vprt(2),vprt(5))
        if (nlmi1 > 9) then
        call awrit1('              %7:1d',' ',180,i1mach(2),vprt(10))
        if (nlmi1 > 16) then
        call awrit1('              %9:1d',' ',180,i1mach(2),vprt(17))
        if (nlmi1 > 25) then
        call awrit1('              %11:1d',' ',180,i1mach(2),vprt(25))
        endif
        endif
        endif
        if (field) then
          call awrit1('               includes E=%3:-2,4;4d',
     .                ' ',120,i1mach(2),efield)
          vprt(1) = sqrt(2d0)*alat*ddot(3,efield,1,bas(1,ib),1)
          vprt(2:4) = sqrt(2d0)*efield(1:3)
          call awrit2('             from field, e dV=  %d, %3:1d',
     .                ' ',180,i1mach(2),vprt,vprt(2))
        endif
        if (.not. TBU) then
          U = getavU(nl,nsp,qnu,idxdn,ic)
          J = getavJ(nl,jh(1,ib),idxdn,ic)
          if (nl == 3) then
            dq(0) = rho(1,1,ib) - qnu(1,0,1,ic)
            dq(1) = rho(2,1,ib) - qnu(1,1,1,ic)
            dq(2) = rho(3,1,ib) - qnu(1,2,1,ic)
            if (Uav) then
              call awrit3
     .        ('        Hubbard potential U=%d: '//
     .          '%1:1d, dq:%1:1d',' ',180,i1mach(2),
     .        U,vu(0,ib),qmpol(1,ib))
              if (nsp == 2) then
                call awrit3
     .          ('        Stoner potentials J=%d: '//
     .            '%2:1d, moment: %d',' ',180,i1mach(2),
     .          J,vj(1,ib),mmom(ib))
              endif
            else
              call awrit6
     .        ('        Hubbard potentials U=%d %d %d: '//
     .          '%3:1d, dq:%3:1d (total=%d)',' ',180,i1mach(2),
     .        qnu(3,0,1,ic),qnu(3,1,1,ic),qnu(3,2,1,ic),
     .        vu(0,ib),dq,dsum(3,dq,1))
            endif
          endif
          if (nl == 2) then
            dq(0) = rho(1,1,ib) - qnu(1,0,1,ic)
            dq(1) = rho(2,1,ib) - qnu(1,1,1,ic)
            if (Uav) then
              call awrit3
     .        ('        Hubbard potential U=%d: '//
     .          '%1:1d, dq:%1:1d ',' ',120,i1mach(2),U,
     .        vu(0,ib),qmpol(1,ib))
              if (nsp == 2) then
                call awrit3
     .          ('        Stoner potentials J=%d: '//
     .            '%2:1d, moment: %d',' ',180,i1mach(2),
     .          J,vj(1,ib),mmom(ib))
              endif
            else
              call awrit5
     .        ('        Hubbard potentials U=%d %d: '//
     .          '%2:-2d, dq:%2:-2d (total=%d)',' ',180,i1mach(2),
     .        qnu(3,0,1,ic),qnu(3,1,1,ic),
     .        vu(0,ib),dq,dsum(2,dq,1))
            endif
          endif
          if (nl == 1) then
            dq(0) = rho(1,1,ib) - qnu(1,0,1,ic)
            call awrit3
     .        ('        Hubbard potential U=%d: '//
     .        '%d, dq: %d ',' ',120,i1mach(2),U,vu(0,ib),dq)
          endif
          do isp = 1, nsp
            call dcopy(nl**2,dh(1,1,ib,isp),nl**2+1,dhd,1)
            if (isp == 1) then
              if (nl > 2) then
                call awrit3('        DH^e_LL=%1:-2d, %3:-2d, %5:-2d',
     .            ' ',180,i1mach(2),dhd,dhd(2),dhd(5))
              elseif (nl > 1) then
                call awrit2('        DH^e_LL=%1:-2d, %3:-2d',' ',180,
     .            i1mach(2),dhd,dhd(2))
              endif
            else
              if (nl > 2) then
                call awrit3('                %1:-2d, %3:-2d, %5:-2d',
     .            ' ',180,i1mach(2),dhd,dhd(2),dhd(5))
              elseif (nl > 1) then
                call awrit2('                %1:-2d, %3:-2d',' ',180,
     .            i1mach(2),dhd,dhd(2))
              endif
            endif
          enddo
          if (lov .and. iprint() >= 40) then
            call awrit3('        Overlap terms: B=%d A=%d pot0=%d',' ',
     .                  120,i1mach(2),B(ib),A(ib),pot0(ib))
          endif
        endif
        if (iprint() > 40 .or. (TBU .and. iprint() > 31)) then
          print *, 'DH^e_LL'':'
          do  isp = 1, nsp
            if (nsp == 2) then
              if (isp == 1) print*, ' spin up:'
              if (isp == 2) print*, ' spin down:'
            endif
            do  i = 1, nl**2
              write (*,300) (dh(i,k,ib,isp),k=1,nl**2)
            enddo
          enddo
        endif
C --- Electric field gradient ---
C        if (dasum(5,v(5),1) > 10*dsqrt(d1mach(3))) then
C          efg(1) = v(5)*dsqrt(4d0*pi/5d0)
C          efg(2) = v(6)*dsqrt(4d0*pi/5d0)
C          efg(3) = v(9)*dsqrt(4d0*pi/5d0)
C          efg(4) = v(7)*dsqrt(4d0*pi/5d0)
C          efg(5) = v(8)*dsqrt(4d0*pi/5d0)
C          call efgrad(z(ic),efg)
C        endif
        if (force .and. iprint() >= 20)
     .    call awrit1('        Force=%3:1d',' ',120,i1mach(2),f(1,ib))
      enddo
 1000 continue
      print *
      if (iprint() >= 30) then
        if (TBU) then
          write(*,425) sumV/2d0,sumU/2d0
        else
          if (nsp == 1) then
            write(*,450) sumV/2d0,sumU/2d0
          else
            write(*,475) sumV/2d0,sumU/2d0,sumJ/2d0
          endif
        endif
        write (*,500) ecorr
      endif

   10 format ('  L    L''   L''''  L'''''' l    l''   l''''  l''''''',
     .14x,' R^k',25x,' A^k',10x,'  V (accumulated)')
   20 format (8(i3,2x),7f10.6)
   50 format ('  L    L''   L''''  L'''''' l    l''   l''''  l'''''' s
     .      dn(-s)    dn(s)        U        J       U-J')
   75 format (9(i3,2x),5f10.6)
C 100 format ('  L''   L''''  L    l''   l''''  l      M         CG
C    .      rho_L''L''''')
  150 format ('  L''   L''''  L    l''   l''''  l      M         CG
     .   V_m        dQ      U     V_u       DH')
  175 format ('  L''   L''''  L    l''   l''''  l      M         CG
     .   V_m       DH')
C 120 format(1028g12.4)
C 200 format (6(i3,2x),2(2x,f6.2,2x),2f10.6)
  210 format (6(i3,2x),2(2x,f6.2,2x),1x,f10.6,25x,f10.6)
  250 format (6(i3,2x),2(2x,f6.2,2x),2f11.6,f4.1,2f10.6)
  300 format (5x,9f10.6)
  400 format (1028f10.6)
  425 format ('   (1/2) dQ dV             : ',f12.6/
     .        '   E_U                     : ',f12.6)
  450 format ('   (1/2) dQ dV             : ',f12.6/
     .        '   (1/2) U dN^2            : ',f12.6)
  475 format ('   (1/2) dQ dV             : ',f12.6/
     .        '   (1/2) U dN^2            : ',f12.6/
     .        '  -(1/2) J dN_u^2+dN_d^2   : ',f12.6)
  500 format ('   E_2                     : ',f12.6)

      call tcx('tbesel')
      end subroutine tbesel

      subroutine getM(ilm,ilmp,ilmpp,qpol,M)
C-
C ----------------------------------------------------------------------
Ci Inputs:
Ci
Co Outputs:
Co
Cr Remarks
Cr  The tight binding parameters in qpol are as follows
Cr  qpol(1) = M(011) = M(101)
Cr  qpol(2) = M(112)
Cr  qpol(3) = M(022) = M(202)
Cr  qpol(4) = M(121) = M(211)
Cr  qpol(5) = M(222)
Cr  qpol(6) = M(123)
Cr  qpol(7) = M(224)
Cr  These are converted into Stone's definitions by multiplying the
Cr  values from ctrl by \sqrt{4\pi/(2l+1)}
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer ilm,ilmp,ilmpp
      double precision qpol(10),M
C Local Variables
      integer ll,l,lp,lpp
      double precision fourpi,sqr4pi,fac,dsqrt,datan
      fourpi = 16d0*datan(1d0)
      sqr4pi = dsqrt(fourpi)

      l = ll(ilm)
      lp = ll(ilmp)
      lpp = ll(ilmpp)
      fac = dsqrt(fourpi/(2*l + 1))
      M = 0d0
      if (lp == 1 .and. lpp == 0 .and. l == 1) M = qpol(1)
      if (lp == 0 .and. lpp == 1 .and. l == 1) M = qpol(1)
      if (lp == 1 .and. lpp == 1 .and. l == 2) M = qpol(2)
      if (lp == 2 .and. lpp == 0 .and. l == 2) M = qpol(3)
      if (lp == 0 .and. lpp == 2 .and. l == 2) M = qpol(3)
      if (lp == 2 .and. lpp == 1 .and. l == 1) M = qpol(4)
      if (lp == 1 .and. lpp == 2 .and. l == 1) M = qpol(4)
      if (lp == 2 .and. lpp == 2 .and. l == 2) M = qpol(5)
      if (lp == 1 .and. lpp == 2 .and. l == 3) M = qpol(6)
      if (lp == 2 .and. lpp == 1 .and. l == 3) M = qpol(6)
      if (lp == 2 .and. lpp == 2 .and. l == 4) M = qpol(7)
      M = fac*M
      if (l == 0) then
        if (lp == lpp) then
          M = sqr4pi
        endif
      endif
      end subroutine getM

      double precision function getavU(nl,nsp,qnu,idxdn,ic)
C- Return average U over open channels
      implicit none
C Passed parameters
      integer nl,nsp,ic,idxdn(nl,*)
      double precision qnu(3,nl,nsp,*)
C Local variables
      integer l,isp,nchan

      nchan = 0
      getavU = 0d0
      do  isp = 1, nsp
        do  l = 1, nl
          if (idxdn(l,ic) == 1) then
            nchan = nchan + 1
            getavU = getavU + qnu(3,l,isp,ic)
          endif
        enddo
      enddo
      if (nchan. gt. 0) getavU = getavU / nchan
      end function getavU

      double precision function getavJ(nl,jh,idxdn,ic)
C- Return average J over open channels
      implicit none
C Passed parameters
      integer nl,ic,idxdn(nl,*)
      double precision jh(4)
C Local variables
      integer l,nchan

      nchan = 0
      getavJ = 0d0
      do  l = 1, nl
        if (idxdn(l,ic) == 1) then
          nchan = nchan + 1
          getavJ = getavJ + jh(l)
        endif
      enddo
      if (nchan. gt. 0) getavJ = getavJ / nchan
      end function getavJ

      integer function getnlm(nl,idxdn,ic)
C- Return number of lm's at species ic
      implicit none
      integer nl,ic,idxdn(0:nl-1,*),nlm,l
      nlm = 0
      do  l = 0, nl-1
       if (idxdn(l,ic) == 1) nlm = nlm + (2*l + 1)
      enddo
      getnlm = nlm
      end function getnlm

C      subroutine xxxind(dm)
CC- Shuffle 25X9 matrix from MSM lm indices to TBE lm indices
CC ----------------------------------------------------------------------
CCi Inputs:
CCi   dm : 25X9 structure constant sub block
CCo Outputs:
CCo   dm reordered
CCr Remarks
CCr   Michael's structure constants are ordered according to the scheme
CCr         1     2     3     4     5     6     7     8     9
CCr         1     y     z     x    xy    yz 3z^2-r^2 zx   x^2-y^2
CCr   while the TBE programs use the scheme
CCr         1     2     3     4     5     6     7     8     9
CCr         1     x     y     z    xy    yz    zx  x^2-y^2  3z^2-r^2
CCr   This widget rearranges the matrix made by Michael's HBLSTR into
CCr   the TBE order.
CCr   The l>2 ordering is unchanged
CC ----------------------------------------------------------------------
CC     implicit none
CC Passed Parameters
C      double precision dm(25,9)
CC Local Variables
C      double precision wk(25,9)
C      integer i, j, ind(9)
C      data ind /1,4,2,3,5,6,8,9,7/
C
C      do  1  i = 1, 9
C      do  1  j = 1, 9
C        wk(i,j) = dm(ind(i),ind(j))
C    1 continue
C      do  2  i = 10, 25
C      do  2  j = 1, 9
C        wk(i,j) = dm(i,ind(j))
C    2 continue
C      call dcopy(25*9,wk,1,dm,1)
C      end
C
C      subroutine xxxfac(dm)
CC- Include factor in 25X9 structure constant block
CC ----------------------------------------------------------------------
CCi Inputs:
CCi   dm : 25X9 structure constant sub block
CCo Outputs:
CCo   dm * fac
CCr Remarks
CCr   To convert to Stone's convention the l_i,l_j element is multiplied
CCr   by sqrt((2l_j+1)*(2l_i+1))/((2l_j+1)!!(2l_i+1)!!)
CC ----------------------------------------------------------------------
CC     implicit none
CC Passed Parameters
C      double precision dm(25,9)
CC Local Variables
C      integer i,j,li,lj,df(0:4),ll
C      double precision arg,fac
C      external ll
CC ... df(l) = (2l+1)!!
C      data df /1,3,15,105,945/
C
C      do  1  i = 1, 25
C      do  1  j = 1, 9
C        li = ll(i)
C        lj = ll(j)
C        arg = dble((2*lj+1))*dble((2*li+1))
C        fac = dsqrt(arg)/(df(lj)*df(li))
C        dm(i,j) = dm(i,j)*fac
C    1 continue
C      end
