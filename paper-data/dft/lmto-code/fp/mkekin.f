      subroutine mkekin(nbas,ndimh,s_site,s_spec,s_lat,s_ham,lcplxp,k1,
     .  k2,k3,vconst,smpot,smrho,sumev,sumtv)
C- Evaluate the valence kinetic energy
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec qkkl qhkl qhhl tauhh tauhk taukk sighh sighk
Ci                 sigkk pihh pihk pikk sohh sohk sokk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa kmxt lmxb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lncol iprmb
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   ndimh :dimension of hamiltonian
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   osig  :augmentation overlap integrals
Ci   otau  :augmentation kinetic energy integrals
Ci   oppi  :augmentation kinetic + potential integrals
Ci   k1..3 :dimensions smpot,smrho
Ci   vconst:constant potential added to hamiltonian
Ci   smpot :smooth input potential on uniform mesh (mkpot.f)
Ci         :smpot = Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
Ci   smrho :smooth output density n0 (no gaussians n0~-n0)
Ci   sumev :sum of eigenvalues
Co Outputs
Co   sumtv :kinetic energy
Cl Local variables
Cl   sraugm:sum_ib q * (tau+ppi-tau) : corresponds to valftr in mkpot.f
Cl   smresh:sm rho * sm V ; corresponds to valfsm in mkpot.f
Cl   lso   :1 include L.S coupling; 2 include LzSz part only
Co   sumso :site contribution to SO soupling
Cr Remarks
Cr   The valence kinetic energy is evaluated in the usual way as
Cr        sumtv = sumev - srhov
Cr   where sumev is the sum-of-eigenvalues and srhov is the integral
Cr   of the output density and input potential.
Cr   Integrals of the density with the xc potential are folded into the
Cr   electrostatic parts, that is:
Cr     V0 = V0es + V0xc  V1 = V1es + V1xc  V2 = V2es + V2xc
Cr   and are not discussed here.
Cr
Cr   mkekin make the electrostatic integral
Cr     int (n0~ Ves~ + n1 Ves1 - n2 Ves2~)                         (40)
Cr   as described in
Cr      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
Cr      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
Cr      ed. (Springer-Verlag, Berlin) 2000.
Cr   for an output density defined through (smrho,qkkl) and an input
Cr   potential defined through smpot and the matrix elements (ppi-tau).
C
Cr   Consider only one atom/cell for simplicity. The output density is
Cr     n = sum_ij Dij F~i F~j with Dij = {sum_n w_n z*_in z_jn}   (32)
Cr   i.e. the contraction of the density matrix with partial densities
Cr     F~i F~j = Fi Fj +
Cr               sum_kLk'L' C^i_kL (P~kLP~k'L' - PkLPk'L') C^j_k'L' (33)
Cr             = n0_ij + n1_ij - n2_ij
Cr   Note that local parts of the partial densities have two `levels'
Cr   of decomposition, namely at the `ij' level as in Eq. 33, or
Cr   at a still finer level in which the (kLk'L') indices are not
Cr   summed over.  Thus
Cr     n{1,2} = sum_ij D_ij n{1,2}_ij
Cr     n{1,2}_ij = sum_kLk'L' C^i_kL n{1,2}_kL,k'L' C^j_k'L'
Cr     n{1,2}_kL,k'L' = PkL Pk'L'
Cr   Note also that the 'k index' may be a sum over polynomials, or when
Cr   function `heads' are dealt with, the function itself, as described
Cr   in augmat.f.  As in making the matrix elements, we have to deal
Cr   with three cases, HH; HP; PP, but this is a inessential detail
Cr   needed only because representing H with sums of polynomials tends
Cr   to be ill-conditioned, and in the description below we ignore it.
Cr
Cr   Densities n0 and n2 have corresponding n0~ and n2~ which include
Cr   the additional multipole terms that guarantee n1 and n2~ have
Cr   the same multipole moments.  Thus:
Cr     n0~ = n0 + sum_M q_M G_M
Cr     n2~ = n2 + sum_M q_M G_M
Cr   where q_M are the difference in multipole moments between n1 and n2
Cr     q_M = int dr Y_M r^m (n1 - n2)
Cr   We can define partial densities for multipole contributions as well
Cr     n2~-n2 = sum_ij D_ij (n2~-n2)_ij
Cr     (n2~-n2)_ij = sum_M Q_ijM G_M
Cr                 = sum_kLk'L'M C^i_kL Q_kkLL'M G_M C^j_k'L'
Cr   with the two forms decomposing q_M into two levels:
Cr     q_M = sum_ij D_ij Q_ijM
Cr     Q_ijM = sum_kLk'L' C^i_kL Q_kkLL'M C^j_k'L'
Cr     Q_kkLL'M = int dr Y_M r^m (P~kL P~k'L' - PkL Pk'L')         (27)
Cr
Cr   Using the identity
Cr     n2~ - n2 = n0~ - n0 = sum_M q_M G_M
Cr   Eq. 40 is evaluated as
Cr     int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + sum_M q_M G_M (Ves0~-Ves2~)
Cr   The first term is evaluated on the mesh and stored in srmesh
Cr   The remaining terms amount to products of the density-matrix
Cr   and the ppi matrix elements.  Thus:
Cr     int n1 Ves1 = sum_ij D_ij int n1_ij Ves1
Cr     int n1_ij Ves1 = sum_kLk'L' C^i_kL int n1_kL,k'L' Ves1 C^j_k'L'
Cr                    = sum_kLk'L' C^i_kL pi1_kk'LL' C^j_k'L'
Cr   where pi1 is the first term of the pi matrix element, Eq. 29:
Cr     pi1_kk'LL' = P~kL V1 P~k'L'
Cr   Similarly for the second term, substituting n2 for n1 and
Cr   Ves2~ for Ves1.
Cr     int n2 Ves2~ = sum_ij D_ij int n2_ij Ves2~
Cr     int n2_ij Ves2~ = sum_kLk'L' C^i_kL int n2_kL,k'L' Ves2~ C^j_k'L'
Cr                    = sum_kLk'L' C^i_kL pi2_kk'LL' C^j_k'L'
Cr     pi2_kk'LL' = P~kL V1 P~k'L'
Cr   The last term just amounts to products of the density-matrix and
Cr   the remaining parts of the ppi matrix element:
Cr     pi_kk'LL'  = pi1_kk'LL' - pi2_kk'LL' + pi3_kk'LL'
Cr     pi3_kk'LL' = sum_M Q_kkLL'M int G_M (Ves0~ - Ves2~)
Cr   Evaluating the last term in the electrostatic integral we have
Cr     rhoV_MP = int sum_M q_M G_M (Ves0~ - Ves2~)
Cr             = int sum_ij D_ij sum_M Q_ijM G_M (Ves0~ - Ves2~)
Cr             = sum_ij D_ij sum_kLk'L'M C^i_kL pi3_kk'LL' C^j_k'L'
Cr   which follows using the relationship between Q_kkLL'M and Q_ijM
Cr   Using the definition of the local density-matrix (see rlocbl.f)
Cr      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
Cr   the electrostatic integral then becomes
Cr     int rhoVes = int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + rhoV_MP
Cr                = int n0 Ves0~
Cr                + sum_ij D_ij sum_kLk'L' C^i_kL pi_kk'LL' C^j_k'L'
Cr                = int n0 Ves0~ + sum_kLk'L' qpp'LL' pi_kk'LL'
Cr
Cu Updates
Cu   05 Jul 13 tso had been held in ppiz(...,2), now a separate array
Cu             Print out site-resolved tso
Cu   22 Nov 12 Replace qkkl with structures
Cu   10 Nov 11 Begin migration to f90 structures
Cu   29 Jun 05 (MvS) SO hamiltonian not included in d.c. terms when
Cu             evaluating kinetic energy.
Cu   27 Aug 01 Extended to local orbitals.
Cu   18 Jun 00 spin polarized
Cu   20 Jun 00 adapted from nfp get_ekin
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k1,k2,k3,nbas,ndimh,lcplxp
C     integer oppi(3,nbas),osig(3,1),otau(3,1)
      double precision sumev,sumtv,vconst
      double complex smpot(k1,k2,k3),smrho(k1,k2,k3)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Dynamically allocated local arrays
      real(8), allocatable :: tsos(:,:,:)
C ... Local parameters
      character outs*128
      integer ib,ipr,is,kmax,lgunit,lmxa,lmxh,n0,n1,n2,n3,
     .  ngabc(3),nglob,nkap0,nlma,nlmh,norb,nsp,nspc,stdo,lso,isw
      integer k
C     logical lgors
      parameter (n0=10,nkap0=4)
      integer nkaph,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .  ntab(n0*nkap0),blks(n0*nkap0)
      double precision qum1,qum2,sraugm,srhov,srmesh,sum1,sum2,sumh,
     .  sumso(4,2),sumq,sumt,vol,xx
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer, parameter :: lSOf=4,lSzLz=32,lSzLzp=64,lSzLzp0=128

      real(8), pointer :: qkk(:,:),qhk(:,:),qhh(:,:)
      real(8), pointer :: ppikk(:,:),ppihk(:,:),ppihh(:,:)
      real(8), pointer :: sigkk(:,:),sighk(:,:),sighh(:,:)
      real(8), pointer :: taukk(:,:),tauhk(:,:),tauhh(:,:)
      real(8), pointer :: tsokk(:,:),tsohk(:,:),tsohh(:,:)

C#ifdefC DEBUG
C      integer fopna
C#endif

      stdo = lgunit(1)
      nsp  = nglob('nsp')
      nspc = nglob('nspc')
      nkaph = nglob('nkaph')
      call getpr(ipr)
      ngabc = s_lat%nabc
      vol = s_lat%vol
      lso =   isw(IAND(s_ham%lncol,4) /= 0)
     .    + 2*isw(IAND(s_ham%lncol,lSzLz) /= 0)
     .    + 3*isw(IAND(s_ham%lncol,lSzLzp) /= 0)
     .    + 4*isw(IAND(s_ham%lncol,lSzLzp0) /= 0)

C     call zprm3('sm rho-out',0,smrho,k1,k2,k3*nsp)
C#ifdefC DEBUG
C      call info0(0,0,0,' rewinding pi file ...')
C      is = fopna('ppi',-1,4)
C      rewind is
C#endif

C --- Integral n0(out) (Ves0~ + Vxc0), contribution from mesh ---
C     Note that it does not include the term (n0~-n0) Ves0~
      call mshdot(vol,nsp,n1,n2,n3,k1,k2,k3,smpot,smrho,sum1,sum2)
      call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,qum1,qum2)
      srmesh = sum1 + vconst*qum1

C --- Integral rhout*veff, part from augmentation ---
      sraugm = 0d0
      if (lso == 0) then
        allocate(tsohh(1,1),tsohk(1,1),tsokk(1,1),tsos(1,1,1))
      else
        allocate(tsos(nbas,0:4,2))
        call dpzero(tsos,5*2*nbas)
      endif
      do  ib = 1, nbas
        is = s_site(ib)%spec

        lmxa = s_spec(is)%lmxa
        kmax = s_spec(is)%kmxt
        lmxh = s_spec(is)%lmxb
        if (lmxa == -1) cycle

        call orbl(ib,0,ndimh,s_ham%iprmb,norb,ltab,ktab,xx,offl,xx)
C       Block into groups of consecutive l
        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2

        qkk => s_site(ib)%qkkl
        qhk => s_site(ib)%qhkl
        qhh => s_site(ib)%qhhl

C       Take average of qhk(2,1)+qhk(2,1)* since make (1,2) block only
        if (nspc == 2) then
          k = size(qhk)/6
          call dsumdf(k*2,0.5d0,qhk(1,3),0,1,qhk(1,5),0,1)
        endif

        tauhh => s_site(ib)%tauhh
        tauhk => s_site(ib)%tauhk
        taukk => s_site(ib)%taukk
        sighh => s_site(ib)%sighh
        sighk => s_site(ib)%sighk
        sigkk => s_site(ib)%sigkk
        ppihh => s_site(ib)%pihh
        ppihk => s_site(ib)%pihk
        ppikk => s_site(ib)%pikk
        if (lso /= 0) then
          tsohh => s_site(ib)%sohh
          tsohk => s_site(ib)%sohk
          tsokk => s_site(ib)%sokk
        endif

        call pvgtkn(kmax,lmxa,nlma,nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .    tauhh,sighh,ppihh,ppihh,tsohh,
     .    tauhk,sighk,ppihk,ppihk,tsohk,
     .    taukk,sigkk,ppikk,ppikk,tsokk,
     .    lcplxp,lso,qhh,qhk,qkk,nsp,nspc,sumt,sumq,sumh,sumso)
C       Add site augmentation contribution to rhout * (ham - ke)
        sraugm = sraugm + sumh - sumt
        if (lso /= 0) then
          tsos(ib,1:4,1) = sumso(1:4,1)
          tsos(ib,1:4,2) = sumso(1:4,2)
          tsos(ib,0,1) = sum(tsos(ib,1:4,1))
          tsos(ib,0,2) = sum(tsos(ib,1:4,2))
        endif

C       Restore original (1,2) and (2,1) blocks of qhk
        if (nspc == 2) then
          k = size(qhk)/6
          call dsumdf(k*2,1d0,qhk(1,3),0,1,qhk(1,5),0,1)
        endif

      enddo
      if (lso == 0) deallocate(tsohh,tsohk,tsokk)

      srhov = srmesh + sraugm
      sumtv = sumev - srhov
      if (ipr >= 30) write(stdo,340) srmesh,sraugm,srhov,sumev,sumtv
  340 format(/' srhov:',3f14.6,' sumev=',f12.6,'   sumtv=',f12.6)

      if (ipr >= 45 .and. lso == 1) then
        write(stdo,341)
  341   format(/' Site resolution of <hso rhout>')
        outs = ' Site    ++%10f+-%10f-+%10f--%10ftot%9fhead'
        call arrprt(outs,
     .    '%,4i%;12,7D%;12,7D%;12,7D%;12,7D%;12,7D%;12,7D','Idddddd',
     .    nbas,0,1,0,'  | ',xx,tsos(1,1,1),tsos(1,2,1),tsos(1,3,1),
     .    tsos(1,4,1),tsos(1,0,1),tsos(1,0,2),xx)

C        outs = ' Site   Tso         Tso(hh)'
C        call arrprt(outs,'%,4i%;12,7D%;12,7D','Idd',nbas,0,3,
C     .    0,'  | ',xx,tsos(1,0,1),tsos(1,0,2),xx,xx,xx,xx,xx)
C        stop
      endif

      deallocate(tsos)
      end

      subroutine pvgtkn(kmax,
     .  lmxa,nlma,nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .  tauhh,sighh,ppihh,ppihhz,tsohh,
     .  tauhp,sighp,ppihp,ppihpz,tsohp,
     .  taupp,sigpp,ppipp,ppippz,tsopp,
     .  lcplxp,lso,qhh,qhp,qpp,nsp,nspc,sumt,sumq,sumh,sumso)
C- Local contribution to kinetic energy for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :cutoff in PkL expansion
Ci   lmxa  :dimensions sigpp, taupp
Ci   nlma  :L cutoff in PkL expansion
Ci   nkaph :dimensions augmentation matrices
Ci   norb  :number of orbitals for this site
Ci   ltab  :table of l quantum numbers for the orbitals
Ci   ktab  :table of k numbers (orbital type) for the orbitals
Ci   blks  :block size for grouping orbitals into blocks (gtbls1)
Ci   lmxh  :dimensions sighh, sighp, tauhh, tauhp
Ci   nlmh  :dimensions heads ppi and qhh and qhp
Ci   tauhh :head-head kinetic energy integrals (augmat.f)
Ci   sighh :head-head overlap integrals (augmat.f)
Ci   ppihh :head-head kinetic + potential integrals (augmat.f)
Ci   tauhp :head-tail kinetic energy integrals (augmat.f)
Ci   sighp :head-tail overlap integrals (augmat.f)
Ci   ppihp :head-tail kinetic + potential integrals (augmat.f)
Ci   taupp :tail-tail kinetic energy integrals (augmat.f)
Ci   sigpp :tail-tail overlap integrals (augmat.f)
Ci   ppipp :tail-tail potential integrals (augmat.f)
Ci   ppihhz:Complex analog of ppihh, used in place of ppihh if lcmplx>0
Ci   ppihpz:Complex analog of ppihp, used in place of ppihp if lcmplx>0
Ci   ppippz:Complex analog of ppipp, used in place of ppipp if lcmplx>0
Ci   tsoxx :is SO part of K.E..  It is hermitian e.g.
Ci         :tsohh(ik,jk,ilm,jlm,1,2) = tsohh*(jk,ik,jlm,ilm,2,1)
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   lso   :1 include L.S coupling; 2 include LzSz part only
Ci   qhh   :head-head density matrix for this site
Ci   qhp   :head-tail density matrix for this site
Ci   qpp   :tail-tail density matrix for this site
Ci   nsp   :number of spin channels
Ci   nspc  :2 for coupled spins; otherwise 1
Co Outputs
Co   sumt  :site contribution to kinetic energy
Co   sumq  :site contribution to overlap (charge ?)
Co   sumh  :site contribution to kinetic energy + potential
Co   sumso :site contribution to SO soupling
Co         :sumso(1:4,1) = total site contribution
Co         :sumso(1:4,2) = head part of site contribution
Cr Remarks
Cr   Tso should not be calculated this way ... Re[Im Thh* qhh] = 0
Cu Updates
Cu   09 Jul 13 Separates out SO contribution and prints explicitly
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   28 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kmax,lmxa,nlma,lmxh,nlmh,nsp,nspc,lcplxp,lso
      integer nkaph,norb,ltab(norb),ktab(norb),blks(norb)
      double precision sumt,sumq,sumh,sumso(4,2)
      double precision
     .  tauhh(nkaph,nkaph,0:lmxh,*),sighh(nkaph,nkaph,0:lmxh,*),
     .  tauhp(nkaph,0:kmax,0:lmxh,*),sighp(nkaph,0:kmax,0:lmxh,*),
     .  taupp(0:kmax,0:kmax,0:lmxa,*),sigpp(0:kmax,0:kmax,0:lmxa,*),
     .  ppihh(nkaph,nkaph,nlmh,nlmh,*),qhh(nkaph,nkaph,nlmh,nlmh,*),
     .  ppihp(nkaph,0:kmax,nlmh,nlma,*),qhp(nkaph,0:kmax,nlmh,nlma,*),
     .  ppipp(0:kmax,0:kmax,nlma,nlma,*),qpp(0:kmax,0:kmax,nlma,nlma,*)
      double complex
     .  ppihhz(nkaph,nkaph,nlmh,nlmh,nsp*nspc),
     .  ppihpz(nkaph,0:kmax,nlmh,nlma,nsp*nspc),
     .  ppippz(0:kmax,0:kmax,nlma,nlma,nsp*nspc),
     .  tsohh(nkaph,nkaph,nlmh,nlmh,2,2),
     .  tsohp(nkaph,0:kmax,nlmh,nlma,2,2),
     .  tsopp(0:kmax,0:kmax,nlma,nlma,2,2)
C ... Local parameters
      integer ilm1,ilm2,k1,k2,ll,nlm11,nlm12,nlm21,nlm22,i
      integer io1,io2,l1,l2
C     double precision xx
      double complex u(2,2),sum1,sum2,sum3,zdotc,sumt1,sumt2,sumt3
      complex(8), allocatable :: qhh3(:,:,:,:,:),phh3(:,:,:,:,:)
      complex(8), allocatable :: qhp3(:,:,:,:,:),php3(:,:,:,:,:)
      complex(8), allocatable :: qpp3(:,:,:,:,:),ppp3(:,:,:,:,:)


C ... Remove SO part from potential
C      if (lso /= 0) then
C        i = nkaph*nkaph*nlmh*nlmh
CC       call zprm('ppiz+so',2,ppihhz,i,i,nsp*nspc)
C        call daxpy(2*i*nsp*nspc,-1d0,tsohh,1,ppihhz,1)
CC       call ppi2z(6,nsp,nspc,nkaph,nkaph,nlmh,nlmh,xx,ppihhz,ppihhz)
C        i = nkaph*(1+kmax)*nlmh*nlma
C        call daxpy(2*i*nsp*nspc,-1d0,tsohp,1,ppihpz,1)
CC       call ppi2z(6,nsp,nspc,nkaph,1+kmax,nlmh,nlma,xx,ppihpz,ppihpz)
C        i = (1+kmax)*(1+kmax)*nlma*nlma
C        call daxpy(2*i*nsp*nspc,-1d0,tsopp,1,ppippz,1)
CC       call ppi2z(6,nsp,nspc,1+kmax,1+kmax,nlma,nlma,xx,ppippz,ppippz)
C      endif

      sumt1 = 0; sumt2 = 0; sumt3 = 0; sumso = 0
      if (nspc == 2) then

        allocate(qpp3(0:kmax,0:kmax,nlma,nlma,4))
        allocate(ppp3(0:kmax,0:kmax,nlma,nlma,4))
        i = (kmax+1)**2*nlma**2
        call rotspv(214,u,i,i,i,qpp,qpp,qpp3,qpp3)
        call rotspv(213,u,i,i,i,ppippz,ppippz,ppp3,ppp3)
        sum3 = zdotc(4*i,qpp3,1,ppp3,1)
        call rotspv(216,u,i,i,i,ppippz,ppippz,ppp3,ppp3)
        sum3 = zdotc(4*i,qpp3,1,ppp3,1)
        if (lso /= 0) then
        call rotspv(216,u,i,i,i,tsopp,tsopp,ppp3,ppp3)
        sumt3 = zdotc(4*i,qpp3,1,ppp3,1)
        do  k1 = 1, 4
          sumso(k1,1) =
     .      zdotc(i,qpp3(0,0,1,1,k1),1,ppp3(0,0,1,1,k1),1)
        enddo

        endif
C       Debugging
C        call snot(0,1,u)
C        call rotspv(014,u,i,i,i,qpp,qpp,qpp3,qpp3)
C        call rotspv(016,u,i,i,i,ppippz,ppippz,ppp3,ppp3)
CC       call zprm('qpp3',2,qpp3,i,i,4)
CC       call zprm('ppp3',2,ppp3,i,i,4)
C        sum3 = zdotc(4*i,qpp3,1,ppp3,1)

        deallocate(qpp3,ppp3)

        i = nkaph*(kmax+1)*nlmh*nlma
        allocate(qhp3(nkaph,0:kmax,nlmh,nlma,4))
        allocate(php3(nkaph,0:kmax,nlmh,nlma,4))
        call rotspv(214,u,i,i,i,qhp,qhp,qhp3,qhp3)
        call rotspv(213,u,i,i,i,ppihpz,ppihpz,php3,php3)
        ! Caution: some versions of zdotc in the libraries packaged with gfortran are buggy
        sum2 = zdotc(4*i,qhp3,1,php3,1)
        call rotspv(216,u,i,i,i,ppihpz,ppihpz,php3,php3)
        sum2 = zdotc(4*i,qhp3,1,php3,1)
        if (lso /= 0) then
        call rotspv(216,u,i,i,i,tsohp,tsohp,php3,php3) ! Noncollinear tsohp
        sumt2 = zdotc(4*i,qhp3,1,php3,1)
        do  k1 = 1, 4
          sumso(k1,1) = sumso(k1,1) +
     .      zdotc(i,qhp3(1,0,1,1,k1),1,php3(1,0,1,1,k1),1)
        enddo
        endif

C       Debugging
C        call snot(0,1,u)
C        call rotspv(014,u,i,i,i,qhp,qhp,qhp3,qhp3)
C        call rotspv(016,u,i,i,i,ppihpz,ppihpz,php3,php3)
C        sum2 = zdotc(4*i,qhp3,1,php3,1)

C        call zprm('qhp3',2,qhp3,i,i,4)
C        call zprm('php3',2,php3,i,i,4)

        deallocate(qhp3,php3)

        i = nkaph**2*nlmh**2
        allocate(qhh3(nkaph,nkaph,nlmh,nlmh,4))
        allocate(phh3(nkaph,nkaph,nlmh,nlmh,4))
        call rotspv(214,u,i,i,i,qhh,qhh,qhh3,qhh3)
        call rotspv(213,u,i,i,i,ppihhz,ppihhz,phh3,phh3) ! Collinear ppi
        sum1 = zdotc(4*i,qhh3,1,phh3,1)
        call rotspv(216,u,i,i,i,ppihhz,ppihhz,phh3,phh3) ! Noncollinear ppi
        sum1 = zdotc(4*i,qhh3,1,phh3,1)
        if (lso /= 0) then
        call rotspv(216,u,i,i,i,tsohh,tsohh,phh3,phh3) ! Noncollinear tsohh
        sumt1 = zdotc(4*i,qhh3,1,phh3,1)
        do  k1 = 1, 4
          sumso(k1,2) = zdotc(i,qhh3(1,1,1,1,k1),1,phh3(1,1,1,1,k1),1)
          sumso(k1,1) = sumso(k1,1) + sumso(k1,2)
        enddo
        endif
C       Debugging
C        call snot(0,1,u)
C        call rotspv(016,u,i,i,i,ppihhz,ppihhz,phh3,phh3)
C        call rotspv(014,u,i,i,i,qhh,qhh,qhh3,qhh3)
C        sum1 = zdotc(4*i,qhh3,1,phh3,1)
        deallocate(qhh3,phh3)

      endif

C#ifdefC DEBUG
C      integer fopna
C      i = fopna('ppi',-1,4)
C      call info0(0,0,0,' PVGTKN reading ppi')
C      read (i) ppihh
C      read (i) ppipp
C      read (i) ppihp
C#endif

C      print *, '!!'
C      ppihh = 0
C      tauhh = 0
C      ppihp = 0
C      tauhp = 0
C      ppipp = 0
C      taupp = 0

      sumt = 0d0; sumq = 0d0; sumh = 0d0

C --- For each spin, do ---
      do  i = 1, nsp

C ... Pkl*Pkl
      do  k1 = 0, kmax
      do  k2 = 0, kmax

      if (lcplxp == 0) then   ! Real potential
        do  ilm1 = 1, nlma
          l1 = ll(ilm1)
          sumt = sumt + qpp(k1,k2,ilm1,ilm1,i)*taupp(k1,k2,l1,i)
          sumq = sumq + qpp(k1,k2,ilm1,ilm1,i)*sigpp(k1,k2,l1,i)
          do  ilm2 = 1, nlma
          sumh = sumh + qpp(k1,k2,ilm1,ilm2,i)*ppipp(k1,k2,ilm1,ilm2,i)
          enddo
        enddo

      else                      ! Complex potential
        do  ilm1 = 1, nlma
        l1 = ll(ilm1)
        sumt = sumt + qpp(k1,k2,ilm1,ilm1,i)*taupp(k1,k2,l1,i)
        sumq = sumq + qpp(k1,k2,ilm1,ilm1,i)*sigpp(k1,k2,l1,i)
        enddo
        if (nspc == 1 .or. .true.) then
        do  ilm1 = 1, nlma
        do  ilm2 = 1, nlma
          sumh = sumh + qpp(k1,k2,ilm1,ilm2,i)*ppippz(k1,k2,ilm1,ilm2,i)
        enddo
        enddo
        else
          do  ilm1 = 1, nlma
          do  ilm2 = 1, nlma
            sumh = sumh +qpp3(k1,k2,ilm1,ilm2,i)*ppp3(k1,k2,ilm1,ilm2,i)
          enddo
          enddo
        endif
      endif                     ! real or complex branches
      enddo
      enddo                     !tail-tail

C ... Hsm*Hsm
      do  io2 = 1, norb
      if (blks(io2) /= 0) then
C       k2,l2 = k and starting l index for this block
        l2 = ltab(io2)
        k2 = ktab(io2)
        nlm21 = l2**2+1
        nlm22 = nlm21 + blks(io2)-1
        do  io1 = 1, norb
        if (blks(io1) /= 0) then
C         k1,l1 = k and starting l index for this block
          l1 = ltab(io1)
          k1 = ktab(io1)
          nlm11 = l1**2+1
          nlm12 = nlm11 + blks(io1)-1
          if (lcplxp == 0) then
          do  ilm1 = nlm11, nlm12
          l1 = ll(ilm1)
C          sumt = sumt + qhh(k1,k2,ilm1,ilm1,i)*tauhh(k1,k2,l1,i)
C          sumq = sumq + qhh(k1,k2,ilm1,ilm1,i)*sighh(k1,k2,l1,i)
C          do  ilm2 = nlm21, nlm22
C            sumh = sumh+ qhh(k1,k2,ilm1,ilm2,i)*ppihh(k1,k2,ilm1,ilm2,i)
C          enddo
          do  ilm2 = nlm21, nlm22
            if (ilm1 == ilm2) then
              sumt = sumt + qhh(k1,k2,ilm1,ilm2,i)*tauhh(k1,k2,l1,i)
              sumq = sumq + qhh(k1,k2,ilm1,ilm2,i)*sighh(k1,k2,l1,i)
            endif
            sumh = sumh+ qhh(k1,k2,ilm1,ilm2,i)*ppihh(k1,k2,ilm1,ilm2,i)
          enddo
          enddo

          else
          do  ilm1 = nlm11, nlm12
          l1 = ll(ilm1)
          do  ilm2 = nlm21, nlm22
            if (ilm1 == ilm2) then
              sumt = sumt + qhh(k1,k2,ilm1,ilm2,i)*tauhh(k1,k2,l1,i)
              sumq = sumq + qhh(k1,k2,ilm1,ilm2,i)*sighh(k1,k2,l1,i)
            endif
            sumh = sumh+qhh(k1,k2,ilm1,ilm2,i)*ppihhz(k1,k2,ilm1,ilm2,i)
          enddo
          enddo
          endif                 ! real or complex branches
        endif
        enddo
      endif
      enddo                     ! End of head-head block

C ... Hsm*Pkl
      do  io1 = 1, norb
      if (blks(io1) /= 0) then
C       k1,l1 = k and starting l index for this block
        l1 = ltab(io1)
        k1 = ktab(io1)
        nlm11 = l1**2+1
        nlm12 = nlm11 + blks(io1)-1
        do  k2 = 0, kmax

        if (lcplxp == 0) then ! Real potential
          do  ilm1 = nlm11, nlm12
          l1 = ll(ilm1)
          sumt = sumt + qhp(k1,k2,ilm1,ilm1,i)*tauhp(k1,k2,l1,i)
          sumq = sumq + qhp(k1,k2,ilm1,ilm1,i)*sighp(k1,k2,l1,i)
          do  ilm2 = 1, nlma
            sumh = sumh+qhp(k1,k2,ilm1,ilm2,i)*ppihp(k1,k2,ilm1,ilm2,i)
          enddo
          enddo
        else                    ! Complex potential
          do  ilm1 = nlm11, nlm12
          l1 = ll(ilm1)
          sumt = sumt + qhp(k1,k2,ilm1,ilm1,i)*tauhp(k1,k2,l1,i)
          sumq = sumq + qhp(k1,k2,ilm1,ilm1,i)*sighp(k1,k2,l1,i)
          do  ilm2 = 1, nlma
            sumh = sumh+qhp(k1,k2,ilm1,ilm2,i)*ppihpz(k1,k2,ilm1,ilm2,i)
          enddo
          enddo
          endif                 ! real or complex branches
        enddo
      endif
      enddo                     !head-tail
      enddo                     !Loop over spins
      if (nspc == 2) then
        sumh = sum1+sum2+sum3
        sumh = sum1+sum2+sum3 - (sumt1+sumt2+sumt3) ! Remove S.O.
C        sumso(1) = sumt1+sumt2+sumt3
C        sumso(2) = sumt1
      endif

C ... Restore SO part from potential
C      if (lso /= 0) then
C        i = nkaph*nkaph*nlmh*nlmh
CC       call zprm('ppiz+so',2,ppihhz,i,i,nsp*nspc)
CC       call ppi2z(5,nsp,nspc,nkaph,nkaph,nlmh,nlmh,xx,ppihhz,ppihhz)
C        call daxpy(2*i*nsp*nspc,1d0,tsohh,1,ppihhz,1)
CC       call zprm('ppiz+so',2,ppihhz,i,i,nsp*nspc)
C        i = nkaph*(1+kmax)*nlmh*nlma
CC       call ppi2z(5,nsp,nspc,nkaph,1+kmax,nlmh,nlma,xx,ppihpz,ppihpz)
C        call daxpy(2*i*nsp*nspc,1d0,tsohp,1,ppihpz,1)
C        i = (1+kmax)*(1+kmax)*nlma*nlma
CC       call ppi2z(5,nsp,nspc,1+kmax,1+kmax,nlma,nlma,xx,ppippz,ppippz)
C        call daxpy(2*i*nsp*nspc,1d0,tsopp,1,ppippz,1)
C      endif

C     call info2(0,0,0,' pvgtkn sumh= %,6;6d sumt= %,6;6d ',sumh,sumt)

      end
