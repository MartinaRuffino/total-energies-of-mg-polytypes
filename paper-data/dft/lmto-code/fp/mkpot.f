      subroutine mkpot(s_site,s_spec,s_lat,s_ham,s_pot,s_optic,nbas,lfrce,
     .  k1,k2,k3,smrho,s_rhat,qbg,smpot,qmom,vconst,ppn,hab,vab,sab,
     .  qval,qsc,gpot0,vval,fes,job,vorb,nlibu,lmaxu,lldau,nbf,bsite)
C- Make the potential from the density
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos class pnu pz v0 v1
Co     Stored:     saxis
Co     Allocated:  *
Cio    Elts passed:v0
Cio    Passed to:  rhomom smves vesgcm mshvmt symvvl ugcomp smvxt
Cio                smvxcm smcorm smvxc4 elocp locpot msh21c
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl z qc a nr rmt rg lfoca rfoca ctail etail stc
Ci                 lmxb p pz name lmxa rs3 eh3 vmtz orbp rsma idu kmxt
Ci                 mxcst coreh coreq
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  rhomom corprm smves vesgcm mshvmt symvvl ugcomp
Cio                smvxt smvxcm smcorm smvxc4 elocp uspecb locpot
Cio                gtpcor msh21c
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc vol ng plat alat pos kv nsgrp gv qlat awald tol
Ci                 nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:nabc ng kv gv ips0 bgv cy symgr ag cg indxcg jcg
Cio                qlv dlv
Cio    Passed to:  ioden2 symsmr smvxte smves vesft vesgcm mshvmt
Cio                symvvl ugcomp ggugbl gfigbl fklbl gklbl hhugbl
Cio                hhigbl phhigb hklbl hsmbl hgugbl smvxt smvxc2 vxcnlm
Cio                smvxcm smcorm locpot msh21c
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  eterms pnudef eula lrsa
Co     Stored:     eterms eula
Co     Allocated:  *
Cio    Elts passed:lncol lrsa eula
Cio    Passed to:  locpot
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  lcplxp vesrmt smvextc smcv0 shfac bxcscali socscl
Co     Stored:     smcv0
Co     Allocated:  smcv0
Cio    Elts passed:smvextc bxcscali v0beta smcv0 nlml
Cio    Passed to:  smvxte smvxt locpot
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic nlg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:kmxax rgrad rgrade
Cio    Passed to:  locpot
Ci   nbas  :size of basis
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rho2 rhoc
Cio    Passed to:  rhomom locpot
Ci Inputs
Ci   nbas  :size of basis
Ci   lfrce :nonzero =>  contribution to forces
Ci   k1,k2,k3 dimensions of smrho for smooth crystal density
Ci   smrho :smooth crystal density, on a uniform mesh
Ci   qbg   :homogeneous background charge
Ci   job   :1s digit
Ci         : 0 stops after energy terms
Ci         : 1 makes potpars also
Ci         :10s digit
Ci         :1 suppress updating potential used to make potpars
Ci         :100s digit
Ci         :1 exclude exchange-correlation potential
Ci         :1000s digit
Ci         :1 Generate valence-only rvvxcv,rvepsv,rvexv,rvecv
Ci         :10000s digit
Ci         :1 write sphere density for each site ib to file rhoMT.ib
Ci         :100000s digit
Ci         : 0 use density as is.
Ci         : 1 reset points of negative density to positive density
Ci             for purposes of calculating vxc
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci         : 3 like 1, but modifed rho is retained.
Ci ... The following are LDA+U inputs
Ci   vorb  :orbital dependent potential
Ci   nlibu :number of U blocks  used to dimension vorb
Ci   lmaxu :max l for U blocks  used to dimension vorb
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat beginning at dmats(*,lldau(ib))
Co Outputs:
Co   smpot :smooth potential on a uniform mesh:
Co         :Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
Co         :... If vext is added, potential is added into smpot
Co   qmom  :multipole moments of valence sphere densities
Co   vconst:additional constant potential term
Co   qval  :total valence charge, including semicore states
Co   qsc   :total charge from semicore states (local orbitals)
Co   ppn   :nmto-like potential parameters
Co   hab,vab,sab: augmentation matrices for local parts of the ham.
Co         :See Remarks in augmat.f for their generation and doc.
Co         :... If vext is present, hab and vab are modified
Co   vval  :coffs to YL expansion of es potential at MT boundary
Co   gpot0 :integrals of gaussians times electrostatic potential
Co         :... If vext is present, its contribution is added to gpot0
Co   fes   :contribution to the force from electrostatic + xc potential
Co         :... If vext is present, its contribution is added to fes
Co   sham->eterms various integrals for the total energy are stored:
Co         :(1)  ehar   --- not touched here
Co         :(2)  eks    --- not touched here
Co         :(3)  utot   = total electrostatic energy
Co         :(4)  valves = valence rho * estat potential
Co         :(5)  cpnves = core+nuc * estat potential
Co         :              Information only; does not affect total energy
Co         :(6)  rhoexc = rho * exc
Co         :(7)  rhovxc = rho * vxc
Co         :(8)  sumec  = sum of foca=0 core eigenvalues
Co         :              Note: an external potential will shift sumec
Co         :              There will be a corresponding shift in rhcvef1
Co         :              so that the kinetic energy will be minimally affected.
Co         :(9)  sumtc  = sum of all core kinetic energies
Co         :(10) rhcvef1= int (rhoc*vef1) for sites w/ lfoca=0
Co         :(11) valvef = smrhov * vsm + sum_ib valvef_ib
Co                        valvef_ib = rhov * vtrue - smrhov * vsm)_ib
Co         :(12) sumt0  = sum of frozen core kinetic energies
Co         :(13) dq1    --- not touched here
Co         :(14) dq2    --- not touched here
Co         :(15) amom   = system magnetic moment
Co         :(16) sumev  = sphere sum-of-eigenvalues --- not touched here
Co         :(17) rhvext = integral rho * vext
Co         :(18) rouvxt --- not touched here
Co         :(19) rhosig --- not touched here
Co         :(20) bfam   --- external B * local moment / 2
Co         :... If vext is added, valvef gets extra term
Co              rhvxt + sum_ib (rhov-smrhov) * vext_ib
Co              No other term in eterms is affected.
Cs Command-line switches
Cs   --oldvc  : Reference potential defined so average cell potential is zero
Cs   --wsmpot : Write smooth potential to file
Cl Local variables
Cl   rvmusm:int rhosm * vxc(rhosm+smcor1) where smcor1 is portion
Cl         :of smooth core density treated nonperturbatively.
Cl   rvepsm:int rhosm * exc(rhosm+smcor1) where smcor1 is portion
Cl         :of smooth core density treated nonperturbatively.
Cl   rvepsv:integral of valence density times exc(valence density)
Cl   rvvxcv:integral of valence density times vxc(valence density)
Cl   rhvsm :integral n0~ phi0~
Cl   sgp0  :compensating gaussians * sm-Ves = int (n0~-n0) phi0~
Cl   sgpe  :compensating gaussians * vext = int (n0~-n0) phiext
Cl   valfsm:rhvsm - sgp0 + rvmusm + fcvxc0 = n0 (phi0~ + Vxc(n0))
Cl         :NB sgp0 associated with the local parts, as is done in
Cl         :making the ppi matrix elements.
Cl   valfat:local contribution to density * veff
Cl         := sgp0 + sum_ib (rhovv1-rhovv2)_ib - focvxa_ib
Cl   rhovv1:generated by locpot = sum_ib {int[rho1*(v1-2*Z/r)]}_ib
Cl         :Note rhovef = rhovv1 - rhovv2
Cl   rhovv2:generated by locpot = sum_ib {int[(rho2*v2) + (n0~-n0)*gpotb] + focvxc}_ib
Cl   lso   :0 no SO coupling
Cl         :1 L.S coupling
Cl         :2 LzSz only
Cl         :3 LzSz + perturbative L+S-
Cl         :4 perturbative L.S
Cl   lcplxp:0 if ppi is real; 1 if ppi is complex
Cl   ns4   :1 if local density is not spin polarized
Cl         :2 if local spin density is collinear along z
Cl         :4 if local spin density is to be rotated;
Cl  vesrmt :vesrmt(ib) electrostatic potential at rmt, site ib
Cr Remarks
Cr *The total density is a sum of three terms,
Cr
Cr    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
Cr
Cr  The first term is the smooth density on a mesh of points; the
Cr  second is the true density and is defined on a radial mesh for each
Cr  sphere; the last is the 1-center expansion of the smooth density on
Cr  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
Cr  identical to the one-center expansion of n0(mesh).  The sum of the
Cr  three terms converges rapidly with l because errors in n_R(r) are
Cr  mostly canceled by errors in n0_R(r).)
Cr
Cr *Computation of the electrostatic energy:
Cr  We add and subtract a set of compensating gaussian orbitals
Cr
Cr    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
Cr
Cr  which render the integral of the local part (the last 3 terms)
Cr  zero in each RL channel.  The g_RL must be localized enough that
Cr  their spillout beyond the MT radius is negligible.
Cr
Cr  We define
Cr
Cr    n0~ = n0 + compensating gaussians
Cr
Cr  In the interstitial, the electrostatic potential of n0~ is the true
Cr  estat potential.  The potential of n0 is called phi0 and the
Cr  potential of n0~ is called phi0~.  The total electrostatic energy
Cr  is computed as
Cr    the electrostatic energy of  n0~  +
Cr    the electrostatic energy of (neutral) local parts
Cr
Cr  The first term is computed in subroutine smves;
Cr  the second term is computed in subroutine locpot.
Cb Bugs
Cb    Not really a bug, but need to check.
Cb    See about understanding and changing foxexc both in locpot
Cb    and smvxc.  The total energy seems to come out just right,
Cb    but focexc isn't real correction to rhov*exc.  Does it matter?
Cb    Maybe can't use lfoca=2 is all.
Cu Updates
Cu   27 Jan 17 New --vext~pot0 option to read reference potential from smpot0
Cu   25 Jan 17 New analysis of response to external potential
Cu   25 Jun 16 New ability to add external potential, through --vext
Cu   17 Jan 16 option to render LDA Bxc = 0 (:100s digit job = 2)
Cu   03 Aug 15 Bug fix for low-lying local orbitals:
Cu             rsml,ehl now calculated from the spin-averaged potential.
Cu             To retain the old convention (spin-1 potential) add 40 to s_ham%pnudef
Cu   09 Dec 13 First cut at using libxc functional for XC potential
Cu   15 Nov 13 (Ben Kaube) first cut at optical matrix elements
Cu             (MvS) first cut at the the noncollinear case
Cu   08 Sep 13 Returns parameters for core total energy in eterms
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   24 Jun 13 Passes shfac to locpot
Cu   25 Oct 11 Started migration to f90 structures
Cu   11 Mar 11 --rhopos switch replaced by 100000s digit job
Cu             local densities can be reset to be positive definite
Cu   01 Apr 10 Add external B field to potential.  New argument list
Cu   21 Apr 09 Handles GGA functionals
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   27 Apr 05 Added LDA+U stuff (Lambrecht)
Cu   24 Dec 04 Changes for full L.S coupling
Cu    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
Cu   15 Jul 04 (Chantis) Matrix elements Lz.Sz for spin-orbit coupling
Cu   14 Jan 02 rvexv and rvecv (T. Miyake)
Cu   17 Sep 01 Returns qsc.  Altered argument list.
Cu   24 Aug 01 Extended to local orbitals.  Altered argument list.
Cu             Local potentials now have correct boundary conditions
Cu   15 Aug 01 Generates rvepsv and rvvxcv.  Changed call to locpot.
Cu   20 Apr 01 Generates vesrmt
Cu   18 Apr 01 Added ability to exclude exchange-correlation potential
Cu   20 Feb 01 Added ppn to potential parameters generated
Cu   15 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp mk_potential.f
C ------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,job,lfrce,lcplxp,k1,k2,k3,n0,nppn,nab
      integer nlibu,lmaxu,lldau(nbas),nbf
      parameter (n0=10,nppn=12,nab=9)
      double precision qbg,qval,qsc,vconst
      double precision bsite(nbas,nbf,3),qmom(*),fes(3,*),gpot0(*),vval(*),
     .  hab(nab,n0,nbas),sab(nab,n0,nbas),vab(nab,n0,nbas),ppn(nppn,n0,nbas)
      double complex smrho(k1,k2,k3,*),smpot(k1,k2,k3,*)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_optic):: s_optic
      type(str_rhat)::  s_rhat(*)
      type(str_atparms) ::  s_atparms
C ... Dynamically allocated local arrays
      real(8), pointer :: vesrmt(:)
      real(8), allocatable :: hpot0(:),fxc(:)
      complex(8), allocatable :: smvxc(:),smvx(:),smvc(:),smexc(:)
      complex(8), allocatable :: zsmr(:,:)
C ... Local parameters
      character*80 outs
      integer i,j,ipl,ipr,kkk,lbf,lso,lxcfun,n1,n2,n3,ngabc(3),ns4,nsp,nspc,nspxc,stdl,stdo
      integer procid
      logical lpot0
      integer, parameter :: master=0
      double precision amom,amoml(0:3),cpnves,dq,eh,eks,focexc,
     .  bfam,focvxc,qsmc,smq,smag(0:3),sum2,rhoexc,rhoex,rhoec,
     .  rhovxc,rhvsm,rhvext,rcvxt,sgp0,sgpe,sumec,sumtc,uat,usm,
     .  utot,valfsm,valfat,valvef,valves,vol,vsum,rhcvef1,zsum,zvnsm,
     .  rhobg,pi,xx,xv(4),eps0
      double precision rvepsv(2),rvexv(2),rvecv(2),rvvxcv(2),fcexc0(2),
     .  rveps(2),rvvxc(2),repsm(2),repsmx(2),repsmc(2),
     .  rvmusm(2),rvepsm(2),rmusm(2),vxcavg(2),fcvxc0(2),fcex0(2),fcec0(2)
C     double precision tmp(2)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer, parameter :: lSOf=4,lSzLz=32,lSzLzp=64,lSzLzp0=128
      integer, parameter :: lBextz=256, lSOEfield=1024
      procedure(logical) :: cmdopt
      procedure(integer) :: a2vec,isum,iprint,isw,lgunit,nglob,mpipid
      procedure(real(8)) :: dlength
      integer,parameter :: NULLI=-99999
      integer :: ieps0 = 0   ! Default : no analysis of eps or addition of screening
      save ieps0

C     This ordering must match sham->eterms; see uham
      double precision eterms(22)
      equivalence (eterms(1),eh)
      equivalence (eterms(2),eks)
      equivalence (eterms(3),utot)
      equivalence (eterms(4),valves)
      equivalence (eterms(5),cpnves)
      equivalence (eterms(6),rhoexc)
      equivalence (eterms(7),rhovxc)
      equivalence (eterms(8),sumec)
      equivalence (eterms(9),sumtc)
      equivalence (eterms(10),rhcvef1)
      equivalence (eterms(11),valvef)
      equivalence (eterms(15),amom)
      equivalence (eterms(17),rhvext)
      equivalence (eterms(20),bfam)

C      print *, '!!'
C      job = job + 30000

      call tcn('mkpot')
      procid = mpipid(1)
      ipr = iprint()
      ipl = ipr
      lso =   isw(IAND(s_ham%lncol,4) /= 0)
     .    + 2*isw(IAND(s_ham%lncol,lSzLz) /= 0)
     .    + 3*isw(IAND(s_ham%lncol,lSzLzp) /= 0)
     .    + 4*isw(IAND(s_ham%lncol,lSzLzp0) /= 0)
      lbf =   isw(IAND(s_ham%lncol,8) /= 0)
     .    + 2*isw(IAND(s_ham%lncol,lBextz) /= 0)
      if (lso /= 0 .and. IAND(s_ham%lncol,lSOEfield) /= 0) lso = lso+10
      lcplxp = s_pot%lcplxp
      if (lso+lbf /= 0 .and. lcplxp == 0)
     .  call rx('mkpot: incompatible lso,lcplxp')
      if (isum(nbas,lldau,1) /= 0 .and. lcplxp == 0)
     .  call rx('mkpot: incompatible ldau,lcplxp')
      stdo = lgunit(1)
      stdl = lgunit(2)
      nsp  = nglob('nsp')
      nspc = nglob('nspc')
      lxcfun = nglob('lxcf')
      pi   = 4d0*datan(1d0)
      ngabc = s_lat%nabc
      vol = s_lat%vol
      eterms = s_ham%eterms
      ns4 = nsp*nspc; if (mod(s_ham%lrsa,10) == 0) ns4 = nsp ! See "Local variables"
      vesrmt => s_pot%vesrmt
      rhvext = 0
      eps0 = -NULLI             ! eps0 will remain NULLI unless ioden2 parses ~eps; see below

C     s_rhat => s_pot%rhat
C     smrho = reshape(s_pot%smrho,(/i1,i2,i3,nsp/))
C     smrho => s_pot%smrho

C ... Rotate sm rho to quantization axis
C      print *, sum(smrho(:,:,:,1)),sum(smrho(:,:,:,2))
      if (ns4 == 4) then
C      print *, sum(smrho(:,:,:,1)),sum(smrho(:,:,:,2))
C      print *, sum(abs(smrho(:,:,:,3))),sum(abs(smrho(:,:,:,4)))
        allocate(zsmr(k1*k2*k3,4))
        call rotspa(5,k1*k2*k3,k1*k2*k3,n1*n2*n3,smrho,smrho,zsmr)
C      print *, sum(smrho(:,:,:,1)),sum(smrho(:,:,:,2))
C      pause
      endif
C      print *, sum(smrho(:,:,:,1)),sum(smrho(:,:,:,2))

C --- Ensure density is positive definite ---
      if (mod(job,1000000) >= 100000) then
        call smrpos(smrho,vol,k1,k2,k3,n1,n2,n3,nsp,xx)
      endif

C --- Printout for smooth background charge ---
      if (qbg /= 0) then
        rhobg = (3d0/4d0/pi*vol)**(1d0/3d0)
        call info5(20,1,0,' Energy for background charge'//
     .    ' q=%d, radius r=%;3d :  E = 9/5*q*q/r = %;4d',
     .    qbg,rhobg,1.8d0*qbg*qbg/rhobg,0,0)
      endif

C --- Smooth electrostatic potential ---
C     print *, '!!'; call pshpr(40)
      call rhomom(0,1,nbas,s_site,s_spec,s_rhat,qmom,vsum)
C     call poppr
      vconst = -vsum/vol
C      if (ipr >= 40) write (stdo,334) vsum,vconst
C  334 format(' vsum=',f14.6,'   vconst=',f14.6)
      allocate(hpot0(nbas)); call dpzero(hpot0,nbas)

C ... Read in external potential, use smpot as work array
      if (cmdopt('--vext',6,0,outs)) then
        allocate(s_pot%smvextc(s_lat%ng*nsp))
        call ioden2(10,s_lat,outs(7:),nsp,nbas,xv,
     .    s_lat%nabc(1),s_lat%nabc(2),s_lat%nabc(3),
     .    smpot,xv,1,[0],lpot0,xv,0,[0],eps0)
C       Store in reciprocal space, packed form
        call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,-1)
        call gvgetf(s_lat%ng,1,s_lat%kv,k1,k2,k3,smpot,s_pot%smvextc)
C       Flag mkpot to set up analysis of screening vext
        if (eps0 /= -NULLI .and. ieps0 == 0) ieps0 = -1
        if (ieps0 == -1) then  ! Allocate smcv0
          if (eps0 == NULLI) eps0 = 1 ! No specification => unit eps=1
          call smvxte(-1,s_lat,s_pot,j,j,j,eps0,smpot,smrho) ! Get j = # nonzero s_pot%smvextc
          call ptr_pot(s_pot,8+1,'smcv0',j,2*nsp,xv)           ! Allocate smcv0
        endif
      endif

      call dpzero(smpot,2*k1*k2*k3*nsp)
      i = 1; if (cmdopt('--oldvc',7,0,outs)) i = 0
      call smves(i,nbas,s_site,s_spec,s_lat,k1,k2,k3,qmom,gpot0,vval,
     .  hpot0,sgp0,smrho,smpot,vconst,smq,qsmc,fes,rhvsm,zvnsm,zsum,vesrmt,qbg)

C ... Initialization for screening analysis
      if (ieps0 == -1) then   ! Copy or read initial potential from file; possibly make screening charge
        j = -7; if (lpot0) j = -8
        call smvxte(j,s_lat,s_pot,n1,n2,n3,eps0,smpot,smrho)
C       Remake etat potential with screening charge
        if (eps0 /= 1) then
          call info0(20,1,0,' Remake ves with screening charge')
          call smves(i,nbas,s_site,s_spec,s_lat,k1,k2,k3,qmom,gpot0,vval,
     .      hpot0,sgp0,smrho,smpot,vconst,smq,qsmc,fes,rhvsm,zvnsm,zsum,vesrmt,qbg)
        endif
        ieps0 = 1               ! Flags smvxte to do screening analysis
      endif

      if (cmdopt('--vext',6,0,outs)) then
        i = ieps0 ; if (mod(job,10) == 0 .and. ieps0 == 1) i = 2
        call smvxt(s_site,s_spec,s_lat,s_pot,nbas,k1,k2,k3,qmom,gpot0,
     .    hpot0,smrho,smpot,fes,rhvext,rcvxt,sgpe,eps0,i)
C       call smvxte(1,s_lat,s_pot,n1,n2,n3,eps0,smpot,smrho)
      endif

      smag = 0
      if (nsp == 2) then
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smag,sum2)
        smag(0) = 2*smag(0) - smq
C       call zprm3('smrho in mkpot',0,smrho,k1,k2,k3)
      endif
      if (ns4 == 4) then
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho(1,1,1,3),smag(1),sum2)
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho(1,1,1,4),smag(2),sum2)
        smag(1) = 2*smag(1)
        smag(2) =-2*smag(2)
        smag(3) =   smag(0)
        smag(0) = dlength(3,smag(1),1)
      endif

      deallocate(hpot0)
C     print *, k1,k2,k3
C     call zprm('smves',2,smpot,k1,k1,k2)

C --- Add smooth exchange-correlation potential ---
      call dpzero(rveps,2)
      call dpzero(rvvxc,2)
      call dpzero(rvepsv,2)
      call dpzero(rvexv,2)
      call dpzero(rvecv,2)
      call dpzero(rvvxcv,2)

C      call dbgsmrho('mkpot smrho before vxc',k1,k2,k3,s_lat,smrho)

      if (mod(job/100,10) /= 1) then
        kkk = k1*k2*k3

C       Case make vxc for spin averaged density
C       Combine densities, treat as though nsp=1
        nspxc = nsp
        if (nsp == 2 .and. mod(job/100,10) == 2) then
          call daxpy(k1*k2*k3*2,1d0,smrho(1,1,1,2),1,smrho,1)
          nspxc = 1
        endif

        allocate(smvxc(kkk*nspxc)); call dpzero(smvxc,2*kkk*nspxc)
        allocate(smvx(kkk*nspxc)); call dpzero(smvx,2*kkk*nspxc)
        allocate(smvc(kkk*nspxc)); call dpzero(smvc,2*kkk*nspxc)
        allocate(smexc(kkk)); call dpzero(smexc,2*kkk)
        allocate(fxc(3*nbas)); call dpzero(fxc,3*nbas)
C   ... Compute int smval*exc(smval) and int smval*vxc(smval)
C       where smval is smoothed valence density
        call pshpr(max(ipr-30,1))
        call smvxc2(0,s_lat,nspxc,lxcfun,s_pot%bxcscali,vol,n1,n2,n3,
     .    k1,k2,k3,smrho,smvxc,smvx,smvc,smexc,xx,rvepsv,rvexv,rvecv,
     .    rvvxcv,vxcavg)
C       call mshdot(vol,nspxc,n1,n2,n3,k1,k2,k3,smrho,smvxc,tmp,tmp(2))
        call poppr
        rvepsv(1) = rvepsv(1) + rvepsv(2)
        rvvxcv(1) = rvvxcv(1) + rvvxcv(2)
C   ... Smooth exchange-correlation potential, incl. smoothed core
        call dpzero(smvxc, kkk*nspxc)
        call dpzero(smexc, kkk)
        call smvxcm(s_site,s_spec,s_lat,nbas,nspxc,s_pot%bxcscali,lfrce,k1,k2,
     .    k3,smrho,smpot,smvxc,smvx,smvc,smexc,repsm,repsmx,repsmc,
     .    rmusm,rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,fxc)
        if (lfrce /= 0) call dpadd(fes,fxc,1,3*nbas,1d0)

        deallocate(smvxc,smvx,smvc,smexc,fxc)

C       Case vxc for average spin: restore spin 1 density, copy potential to second spin channel
        if (nsp == 2 .and. mod(job/100,10) == 2) then
          call daxpy(k1*k2*k3*2,-1d0,smrho(1,1,1,2),1,smrho,1)
          call dcopy(k1*k2*k3*2,smpot,1,smpot(1,1,1,2),1)
          repsm(2) = repsm(1)/2; repsm(1) = repsm(2)
          repsmx(2) = repsmx(1)/2; repsmx(1) = repsmx(2)
          repsmc(2) = repsmc(1)/2; repsmc(1) = repsmc(2)
          rmusm(2) = rmusm(1)/2; rmusm(1) = rmusm(2)
          rvmusm(2) = rvmusm(1)/2; rvmusm(1) = rvmusm(2)
        endif

      else
        call dpzero(repsm,2)
        call dpzero(repsmx,2)
        call dpzero(repsmc,2)
        call dpzero(rmusm,2)
        call dpzero(rvmusm,2)
        call dpzero(fcexc0,2)
        call dpzero(fcex0,2)
        call dpzero(fcec0,2)
        call dpzero(fcvxc0,2)
      endif

C     call dbgsmrho('mkpot smrho after vxc',k1,k2,k3,s_lat,smrho)

C ... Restore sm rho to global axis; rotate potential to axis
      if (ns4 == 4) then
        call rotspa(15,k1*k2*k3,k1*k2*k3,n1*n2*n3,smrho,smrho,zsmr)
        call rotspa(15,k1*k2*k3,k1*k2*k3,n1*n2*n3,smrho,smpot,zsmr)
        deallocate(zsmr)
      endif

      if (procid == master .and. cmdopt('--wsmpot',8,0,outs)) then
        call zprm3('$%N Dumping smpot to file out%N',0,smpot,k1,k2,k3)
      endif

C --- Make parameters for extended local orbitals ---
      i = job ; if (mod(job,10) /= 0 .and. mod(s_ham%pnudef/10,10) >= 4) i = 2
      call elocp(nbas,nsp,s_site,s_spec,i)

C --- Make local potential at atomic sites and augmentation matrices ---
      rhobg=qbg/vol
C     job = job + 200000
      call locpot(nbas,nsp,nspc,lso,lbf,lcplxp,s_pot%v0beta,s_pot,s_ham,s_lat,
     .  s_optic,s_site,s_spec,s_rhat,s_atparms,qmom,vval,gpot0,rhobg,nlibu,lmaxu,
     .  vorb,lldau,nbf,bsite,ppn,hab,vab,sab,qval,qsc,eterms,job)

C ... Combine spin-up and spin-down integrals
      repsm(1)  = repsm(1) + repsm(2)
      repsmx(1) = repsmx(1)+ repsmx(2)
      repsmc(1) = repsmc(1)+ repsmc(2)
      rmusm(1)  = rmusm(1) + rmusm(2)
      rvmusm(1) = rvmusm(1) + rvmusm(2)
      fcexc0(1) = fcexc0(1) + fcexc0(2)
      fcex0(1)  = fcex0(1) + fcex0(2)
      fcec0(1)  = fcec0(1) + fcec0(2)
      fcvxc0(1) = fcvxc0(1) + fcvxc0(2)
      s_atparms%rhoexc(1)  = s_atparms%rhoexc(1) + s_atparms%rhoexc(2)
      s_atparms%rhoex(1) = s_atparms%rhoex(1)+ s_atparms%rhoex(2)
      s_atparms%rhoec(1) = s_atparms%rhoec(1)+ s_atparms%rhoec(2)
      s_atparms%rhovxc(1)  = s_atparms%rhovxc(1) + s_atparms%rhovxc(2)
      s_atparms%focexc(1) = s_atparms%focexc(1) + s_atparms%focexc(2)
      s_atparms%focex(1)  = s_atparms%focex(1) + s_atparms%focex(2)
      s_atparms%focec(1)  = s_atparms%focec(1) + s_atparms%focec(2)
      s_atparms%focvxc(1) = s_atparms%focvxc(1) + s_atparms%focvxc(2)
      if (lbf > 0) bfam = s_atparms%bfam

C --- Total energy terms ---
C ... Integral of valence density times xc potential(valence density)
      if (mod(job,1000) < 100) then
        rvepsv(1) = rvepsv(1) + s_atparms%rvepsv
        rvexv(1)  = rvexv(1)  + s_atparms%rvexv
        rvecv(1)  = rvecv(1)  + s_atparms%rvecv
        rvvxcv(1) = rvvxcv(1) + s_atparms%rvvxcv
      endif
      if (mod(mod(job/1000,10),2) == 1) then
        rveps(1) = rveps(1) + s_atparms%rveps
        rvvxc(1) = rvvxc(1) + s_atparms%rvvxc
      endif

C ... Integral of valence density times estatic potential
      valves = rhvsm + s_atparms%vvesat

C ... Valence density times veff.
C    *Associate term (n0~-n0) Ves(n0~) with local part
C     because of the ppi matrix elements
C    *Also add fcvxc0(1) to smooth part because
C     rvmusm+fcvxc0 is perturbative approximation for rvmusm
C     when cores are not treated perturbatively.
      valfsm = rhvsm + rvmusm(1) - sgp0 - vconst*qbg
      valfat = s_atparms%rhovv1 + sgp0 - s_atparms%rhovv2
      valfsm = valfsm + fcvxc0(1)
      valfat = valfat - s_atparms%focvxc(1)

C ... Add integrals of density * external potential
      if (cmdopt('--vext',6,0,outs)) then
        call info0(31,1,0,
     .    ' Integrals of density with external potential%N'//
     .    '%21flocal%10fsm,loc%10fistl%11fsum')
        call info5(31,0,0,
     .    ' valence rho %;15,6D%;15,6D%;15,6D%;15,6D',
     .    s_atparms%rvext(1,1),s_atparms%rvext(2,1),rhvext,
     .    rhvext+s_atparms%rvext(1,1)-s_atparms%rvext(2,1),0)
        call info2(31,0,0,' gaussian part%14f%3;15,6D',
     .    [s_atparms%sgpote,sgpe,sgpe-s_atparms%sgpote],0)
        call info5(31,0,0,
     .    ' core rho    %;15,6D%;15,6D%;15,6D%;15,6D',
     .    s_atparms%rvext(1,2),s_atparms%rvext(2,2),rcvxt,
     .    rcvxt+s_atparms%rvext(1,2)-s_atparms%rvext(2,2),0)
        rhvext = rhvext + s_atparms%rvext(1,1)-s_atparms%rvext(2,1)
      endif

      valvef = valfsm + valfat

C ... Integral of core+nucleus times estatic potential (for printout only)
      cpnves = zvnsm + s_atparms%cpnves

C ... Total xc energy and potential integral
      focexc = fcexc0(1) - s_atparms%focexc(1)
C     focex  = fcex0(1)  - s_atparms%focex(1)
C     focec  = fcec0(1)  - s_atparms%focec(1)
      focvxc = fcvxc0(1) - s_atparms%focvxc(1)
      if (ipr >= 30 .and. dabs(focexc) > 1d-6)
     .  write (stdo,850) focexc,focvxc
  850 format(' foca xc integrals for spillout charge:',2f12.6)

      repsm(1) = repsm(1) + fcexc0(1)
      repsmx(1)= repsmx(1)+ fcex0(1)
      repsmc(1)= repsmc(1)+ fcec0(1)
      s_atparms%rhoexc(1) = s_atparms%rhoexc(1) - s_atparms%focexc(1)
      s_atparms%rhoex(1) = s_atparms%rhoex(1) - s_atparms%focex(1)
      s_atparms%rhoec(1) = s_atparms%rhoec(1) - s_atparms%focec(1)
      rmusm(1) = rmusm(1) + fcvxc0(1) + fcexc0(1)
      s_atparms%rhovxc(1) = s_atparms%rhovxc(1) - s_atparms%focvxc(1) - s_atparms%focexc(1)

      rhoexc = repsm(1) + s_atparms%rhoexc(1)
      rhoex  = repsmx(1)+ s_atparms%rhoex(1)
      rhoec  = repsmc(1)+ s_atparms%rhoec(1)
      rhovxc = rmusm(1) + s_atparms%rhovxc(1)

C ... Total electrostatic energy
      usm = 0.5d0*(rhvsm+zvnsm)
      uat = 0.5d0*(s_atparms%vvesat+s_atparms%cpnves)
      utot = usm + uat
      dq = smq + s_atparms%qloc + qsmc + s_atparms%qlocc + qbg - zsum
      amoml = smag+s_atparms%aloc
      amom = amoml(0)

C --- Printout ---
      if (ipr >= 20) write(stdo,'(1x)')
      if (ipr >= 30) then
        write (stdo,681)
        write (stdo,680) 'rhoval*vef ',valfsm,valfat,valvef,
     .                   'rhoval*ves ',rhvsm,s_atparms%vvesat,valves,
     .                   'psnuc*ves  ',zvnsm,s_atparms%cpnves,cpnves,
     .                   'utot       ',usm,uat,utot,
     .                   'rho*exc    ',repsm(1),s_atparms%rhoexc(1),rhoexc,
     .                   'rho*vxc    ',rmusm(1),s_atparms%rhovxc(1),rhovxc,
     .                   'valence chg',smq,s_atparms%qloc,smq+s_atparms%qloc
        if (nsp == 2)
     .  write (stdo,680) 'valence mag',smag(0),s_atparms%aloc(0),amom
        if (ns4 == 4) then
          write (stdo,680) '         Mx',smag(1),s_atparms%aloc(1),amoml(1)
          write (stdo,680) '         My',smag(2),s_atparms%aloc(2),amoml(2)
          write (stdo,680) '         Mz',smag(3),s_atparms%aloc(3),amoml(3)
        endif
C       call info2(50,0,0,'   (n0~-n0) (Ves+Vext) %;12,5D',sgp0,0)
        call info5(30,0,0,'   core charge%;21,6D%;17,6D%;17,6D',
     .    qsmc,s_atparms%qlocc,qsmc+s_atparms%qlocc,0,0)
        call info5(30,1,0,
     .    ' Charges:  valence%;12,5D   cores%;12,5D   nucleii%;12,5D%N'//
     .    '    hom background%;12,5D   deviation from neutrality: %;12,5D',
     .    smq+s_atparms%qloc,qsmc+s_atparms%qlocc,-zsum,qbg,dq)

      endif
  680 format(3x,a,4x,3f17.6)
  681 format(' Energy terms:',13x,'smooth',11x,'local',12x,'total')
      if (ipl >= 1 .and. stdl >= 0) then
        write (stdl,710) smq+s_atparms%qloc,smq,s_atparms%qloc,qbg,dq
  710   format('fp qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .         '  qbg',f11.6,' dQ',f11.6)
        if (nsp == 2) write (stdl,711) smag+s_atparms%aloc,smag,s_atparms%aloc
  711   format('fp mag',f11.5,'  sm',f11.5,'  loc',f11.5)
        write (stdl,720) rhovxc,rhoexc,utot
  720   format('fp pot  rvxc',f18.7,'  rexc',f18.7,'  rves',f16.7)
        write (stdl,721) rhoex,rhoec
  721   format('fp pot  rex ',f18.7,'  rec',f19.7)
        if (mod(mod(job/1000,10),2) == 1) then
        write (stdl,722) rvvxcv(1),rvepsv(1)
  722   format('fp pot  rvvxcv',f16.7,'  rvexcv',f16.7)
        write (stdl,723) rvexv(1),rvecv(1)
  723   format('fp pot  rvexv',f17.7,'  rvecv',f17.7)
        endif
c

      endif
      if (dabs(dq) > 1d-3 .and. iprint() > 0) call awrit1(
     .  ' (warning) system not neutral, dq=%d',' ',80,stdo,dq)
      s_ham%eterms = eterms
      if (associated(s_pot%smvextc)) deallocate(s_pot%smvextc)

C#ifdefC DBGVXT
CC      call mkchksum(120+mod(job,10),nsp,nbas,k1,k2,k3,smpot,smrho,vconst,s_pot,hab,vab,sab,s_site,s_rhat)
C#endif

      call tcx('mkpot')

      end

C#ifdefC DBGVXT
C      subroutine mkchksum(job,nsp,nbas,k1,k2,k3,smpot,smrho,vconst,s_pot,hab,vab,sab,s_site,s_rhat)
CC- For debugging
C      implicit none
CC ... Passed parameters
C      integer job,nsp,nbas,k1,k2,k3
C      integer, parameter :: n0=10,nppn=12,nab=9
C      double complex smrho(k1,k2,k3,nsp),smpot(k1,k2,k3,nsp)
C      double precision vconst,
C     .  hab(nab,n0,nsp,nbas),vab(nab,n0,nsp,nbas),sab(nab,n0,nsp,nbas)
CC ... Dynamically allocated local arrays
CC ... For structures
C      include 'structures.h'
CC     type(str_lat)::   s_lat
CC     type(str_ham)::   s_ham
C      type(str_pot)::   s_pot
C      type(str_rhat)::  s_rhat(nbas)
C      type(str_site)::  s_site(nbas)
CC ... Local parameters
C      integer ib,n1,n2,n3
C
CC     real(8), pointer :: sig(:,:),tau(:,:),ppi(:,:),tso(:,:)
C
C      if (mod(job,10) == 0) return
C
C      call info0(1,1,0,' mkpot : cksum potential, density arrays')
C      n1 = k1; n2 = k2; n3 = k3
C
C      ib = mod(job/10,10)
C      if (ib /= 0) then
C        if (ib == 2) then
C        call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,-1)
C        call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,-1)
C        endif
C        call info2(1,0,0,'# cksum global arrays in %?#(n==2)#G#R# space',ib,2)
C        call info5(1,0,0,'# smrho %;10F  smpot %;10F  vconst %;9F',sum(smrho),sum(smpot),vconst,4,5)
C        if (ib == 2) then
C        call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)
C        call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1)
C        endif
C      endif
C
C
C      ib = mod(job/100,10)
C      if (ib /= 0) then
C
C      call info0(1,1,0,'# cksum site arrays')
C      call info0(1,0,0,'# ib   sigkk     taukk     pikk         rho1      rho2       vrmt')
C      do  ib = 1, nbas
C
C        call info5(1,0,0,' %,4i %3:1;9F   %2:1;9F  %;9F',ib,
C     .    [sum(s_site(ib)%sigkk), sum(s_site(ib)%taukk), sum(s_site(ib)%pikk)],
C     .    [sum(s_site(ib)%rho1),sum(s_site(ib)%rho2)],s_pot%vesrmt(ib),5)
C
CC    .    [sum(hab(:,:,:,ib)),sum(vab(:,:,:,ib)),sum(sab(:,:,:,ib))],
C
C      enddo
C      endif
C
C      end
C#endif
