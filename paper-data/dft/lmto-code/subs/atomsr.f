C#define ALENA
      subroutine atscpp(s_ctrl,s_spec,s_pot,is,ic,dclabl,rmax,mode,nl,
     .  nsp,initc,rhozbk,avw,pnu,qnu,ves,dv,eula,bscal,facso,neul,bxc,
     .  ekap,pp,pprel,etot,sumev,qtot,amgm,amag,rhrmx,vrmax,thrpv,sop,
     .  rgrad,pmpol,vintra,clabl,eterms)
C- Make one atom self-consistent and generate new potential parameters
C ----------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lsx loptc smalit nclass nclasp nccomp nspec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lasa lncol lrel lcd lscr lham
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  idmod lmxa z a nr coreh coreq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  gtpcor
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     GFr
Co     Allocated:  *
Cio    Elts passed:qnur
Cio    Passed to:  *
Ci Inputs
Ci   is    :species index
Ci   ic    :class index
Ci   dclabl:class names, packed as real numbers
Ci   rmax  :augmentation radius, in a.u.
Ci   mode  : controls program flow of atom program
Ci         :1s digit = imake : tells sphere program what to generate
Ci         :   Note: nonzero 100s digit may modify this digit.
Ci         :   Also see 100s digit for definition of "given potential"
Ci         :0  Read double counting terms from atom files.
Ci         :   No potential or pot pars calculated; no atom file written
Ci         :   This is a 'cheap' version of imake=4, where
Ci         :   instead of computing the d.c. terms from the given potential,
Ci         :   they are copied from atom files.
Ci         :1  Make self-consistent potential from given pnu,qnu
Ci         :   No potential parameters calculated.
Ci         :2  Like imake=0, but make ppars from given potential
Ci         :   Passed pnu,qnu written to disk.
Ci         :   NB: potential, ppars need not be related to pnu,qnu
Ci         :3  Make self-consistent potential from given pnu,qnu
Ci         :   and potential parameters from resulting potential
Ci         :4  Make sphere double-counting terms from given potential
Ci         :   and supplied moments P,Q.  No internal self-consistency
Ci         :   in the potential, nor are potential parameters generated;
Ci         :   no atom file written.  Using this mode to make terms
Ci         :   in Kohn-Sham energy, e.g. v_in, rho_out.
Ci         :10s digit
Ci         :1  For input v, use s_pot for site ic; do not attempt to read from disk
Ci         :2  For input v, use s_pot for site ic; attempt to read from disk if missing
Ci         :*  For output v, write potential to s_pot(:,ic) if s_pot(:,ic) allocated
Ci         :100s digit: deals with generation of ppars from potential
Ci         :0  normal operation: potential may generated internally by
Ci         :   sphere program (see 1s digit), and becomes the "given potential"
Ci         :1  Potential not generated self-consistently
Ci         :   potential parameters are made if given potential is available
Ci         :   Note: 1s modes 1 and 3 are nonsensical in this case.
Ci         :   If 1s mode= 1 -> aborts; if 1s digit 3, switched to mode 2
Ci         :1000s digit: deals with renormalization of qnu from potential
Ci         :   Potential at rmt
Ci   nl    :(global maximum l) + 1, for dimensioning pp
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   initc : record of what parameters are available.
Ci           1 P,Q   2 pp   4 sop   8 vintra  16 pmpol  32 rgrad
Ci   rhozbk: constant nuclear background density (jellium)
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   ves   :l=0 electrostatic potential at rmax (see Remarks)
Ci   dv    :constant potential shift, passed through to vxc0sp
Ci   ekap  :LMTO energy
Ci   rhozbk:constant background charge, excluded in electrostatics
Ci   bscal :scaling of LDA XC field.  Use bscal=1 for normal vxc.
Ci   facso :Scale SO coupling L.S by facso
Ci Inputs/Outputs
Cio  initc : record of what parameters are available.
Cio        : 1 P,Q   2 pp   4 sop   8 vintra  16 pmpol  32 rgrad
Cio        : may be modiifed on output
Co Outputs
Co   etot,sumev: total and band structure energy of atom, which
Co         :can be made within KKR viewpoint if moments P,Q
Co         :generate sphere potential
Co   qtot  :total charge within sphere
Co   amgm  :difference between spin up and spin down charge (mag mom)
Co   amag  :(noncollinear case) magnetic moment (vector)
Co   rhrmx :density at rmax (not generated for imake=4)
Co   vrmax :total l=0 potential at rmax boundary, by class
Co         :(not generated for imake=4)
Co   sop   :matrix elements of spin orbit coupling
Co   rgrad :radial matrix elements of gradient operator; see rgrme
Co   pmpol :integral (phi-or-phidot * phi-or-phidot * r**l) :
Co         :matrix elements of w.f. * wf * r**l for multipole moments
Co   pp    :(depending on imake) potential parameters
Co   bxc   :(noncollinear) orientation of XC field
Co         :NB: if imake=4, bxc is an input, not an output
Co   pprel :(depending on imake) relativistic potential parameters
Co   eterms:integrals for the total E are accumulated for this sphere
Co         :(1)  ehar   --- not touched here
Co         :(2)  eks    --- not touched here
Co         :(3)  utot   = total electrostatic energy
Co         :(4)  valves --- not used by ASA
Co         :(5)  cpnves --- not used by ASA
Co         :(6)  rhoexc = rho * exc
Co         :(7)  rhovxc = rho * vxc (not needed for total energy)
Co         :(8)  sumec  = sum-of-core eigenvalues
Co         :(9)  sumtc  = sum-of-core K.E (not needed for total energy)
Co         :(10) xcore  = rhoc * total potential
Co         :(11) valvef = rhov * total potential
Co         :(12) sumt0  --- not used by ASA
Co         :(13) dq1    --- not used by ASA
Co         :(14) dq2    --- not used by ASA
Co         :(15) amgm   = system magnetic moment
Co         :(16) sumev  = sphere sum-of-eigenvalues
Co         :(17) rinvxt --- not touched here
Co         :(18) rouvxt --- not touched here
Co         :(19) bmval  = M<B> : magnetic contribution to valvef (asadc)
Cs Command-line switches
Cs   --keepsignm
Cs   -elin=
Cs   --dumprho
Cl Local variables
Ci   imake :generally same as mode, but may be altered; see mode above
Ci   lgdd  :0 q2 is coefficient to phidot**2 - p phi**2
Ci         :1 q2 is coefficient to phidot**2 + phi*phidotdot
Ci         :Both produce the same integrated density; see Remarks.
Ci         :lgdd=0 follows Stuttgart conventions.
Ci         :4 Add 4 to make (phidot,phidotdot) by integrating radial SE outward only
Cl   lfree :true for free atom (program overrides passed rmax);
Cl   lfrz  :0 for soft core
Cl         :1 for frozen core
Cl         :2 for soft core, spin-averaged potential
Cl   lso   :true if to calc spin-orbit coupling parameters.
Cl   lintra:true to calculate intra-atomic dC_i/dq_j
Cl   nrmix :1, maximum number of iterations towards self-consistency
Cl          2, number of previous iterations Anderson mixing for charge
Cl             nrmix>10: no amix; set beta to nrmix/100
Cl    nr   :number of points on the radial mesh
Cl   havedc:T if double-counting terms returned; otherwise F
Cl   GFr(4,4,2,2,l,imu,4):overlap integrals for the calculation of averages [Turek (6.141)] (see rdeq)
Cr Remarks
Cr   Charges qtot and amgm are always calculated.  Sign
Cr   convention of charge is + for electrons, - for nucleii.
Cr
Cr   Potential parameters and total energy are calculated from atomsc
Cr   subject to boundary condition that the potential at rmax
Cr   is zero, regardless of the total net charge inside the sphere.
Cr   See subroutine madpot for discussion of how this choice affects
Cr   the Madelung energy.
Cr
Cr   The potential parameters, however may be shifted by the
Cr   electrostatic energy ves depending on the value of idmod.
Cr   In the default mode of operation (idmod 0 or 1), potential
Cr   parameters pp are calculated for the enu corresponding to the
Cr   pnu; the pp's enu and c are shifted however by the electrostatic
Cr   potential ves so that the pp's as returned to caller correspond
Cr   to electrostatic potential = ves at rmax.  If idmod is 2, the
Cr   potential parameters are generated around enu given by the input
Cr   enu (again assuming that enu corresponds to electrostatic
Cr   potential = ves at rmax).  There is no shifting of enu and c
Cr   unless the parameters are generated internally from potpar.
Cr
Cr   Regardless of the value of idmod, self-consistency is achieved
Cr   keeping pnu and moments fixed.  This means that if any idmod is 2
Cr   the potential and moments are self-consistent with respect
Cr   to the potential parameters.
Cr
Cb Bugs
Cb   Inconsistency in lmx(ic) always assumed to be nl-1 (see sop)
Cu Updates
Cu   28 Oct 17 Extended functionality of mode (formerly imake), flexibility
Cu             in input stream for potential.
Cu   25 Sep 17 New 10s digit imake: supersed atom file w/ potential s_pot%v0
Cu   23 Mar 17 New option (s_ctrl%lasa,4096) to compute phi,phidot
Cu             from outward radial integration only
Cu   01 Aug 16 Redsign of Dirac solver
Cu   21 Apr 15 Further adjustments for fully relativistic GF
Cu   10 Oct 14 Further adjustments in preparation for fully relativistic GF
Cu   14 Nov 13 Some adjustments in preparation for fully relativistic GF
Cu   30 Jun 13 Implement scaling of SO coupling
Cu   08 Apr 13 Complete conversion to f90 pointers
Cu   28 Jan 13 core can be calculated with spin-averaged potential
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   09 Apr 11 Different treatment for generating bxc, noncol case
Cu   26 Oct 08 Do not implement LDAU, IDU modes 4,5 here (move to suldau)
Cu   21 Dec 05 (wrl) allow potential shifts to mimic LDA+U
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   18 Jun 04 (A Chantis) working fully relativistic code
Cu   21 Apr 04 Changes for l- and m-dependent XC fields
Cu    4 Apr 04 Additions for magnetic field in the presence of
Cu             orbital-dependent XC field.
Cu   19 Sep 03 (ATP) Enabled partial core occupation (core holes)
Cu   18 Mar 03 (A Chantis) relativistic potential parameters.
Cu             Altered argument list.
Cu   28 Feb 03 imake=4 implemented.  New argument eterms
Cu   22 Feb 03 Make and Printout <v_xc> inside sphere
Cu   15 Feb 03 SO parameters now include matrix element for ASA Bfield
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   28 Apr 98 potpar can make rgrad
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer,parameter :: nrmx=5001, NULLI=-99999
      integer is,ic,nl,nsp,nr,mode,initc,neul
      double precision amgm,avw,bscal,facso,etot,rhozbk,thrpv
      double precision qtot(ic),amag(3),eterms(22),
     .  rhrmx(ic),vrmax(2,ic),ves(ic),dv(ic),eula(neul,3),
     .  pmpol(nl,nl,2*nl-1,3,nsp,ic),vintra(nl*nl,nsp*nsp,ic),
     .  sop(0:nl-1,nsp,nsp,9,ic),rgrad(4,2,nl,nsp,ic),
     .  pnu(*),qnu(3,nl,nsp,*),pp(6,nl,nsp,ic),bxc(3,ic),
     .  pprel(5,nl,2*nl,2,2,ic),rmax(ic),dclabl(ic)
      character clabl*8
C ... Dynamically allocated local arrays
      real(8), allocatable :: v(:,:),rho(:,:),rhoi(:,:),rhoc(:,:)
      real(8), allocatable :: rofi(:,:),g(:),gp(:)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot) ::  s_pot
      type(str_spec)::  s_spec(*)
C ... Local parameters
      logical lscf,sw,swc,havev,haveso,havegr,havemp,haveva,havepp
      logical lrell,lfree,lxcg,lso,lbf,lintra,loptc,lmpol,lneedv,havedc
      integer ifi,nitmax,ipr,lgdd,lmx,i,j,intopt,imake0,imake1,imake2,imake3,
     .  nlspic,nmix,nn,idmodl(10),nrmix(2),lrel,kcor,lcor,lfrz
      integer, parameter :: n0=10, ncmx=200, nvmx=20
      integer idmod(n0),idu(4)
      double precision a,sumec,sumtc,sumev,ekin,utot,rhoeps,rhomu,rhov,
     .  rmx,ekap,exc(2),thrpvl(10),xx(3),z,qc,ec(ncmx),ev(nvmx),avvxc,
     .  xcore,qval,pnuloc(nl*nsp),qnuloc(3*nl*nsp),bhat(3),bmval,
     .  qcor(2),amomnc,t(n0,2),scalsoc,pzloc(nl),qrel
      double precision uh(4),jh(4)
C     double precision GFr(4,4,2,2,0:(nl-1),2*nl,4)
      double precision qnurl(4,nl,2*nl,2,2)
C     double precision pz(n0,2),qz(3,n0,3)
C     integer idmoz(n0)
      character job*3, lbl*8, outs*100, outs1*100
      procedure(logical) :: aiogen,aiomom,aiopar,aiopot,aiosop,aiorme,aiocor,aiova,aiomp
      procedure(integer) :: nglob,isw,lgunit,fopna
      procedure(real(8)) :: ddot

C ... External calls
      external addzbk,amagnc,asadc,atomsc,awrit0,awrit3,awrit4,awrit5,
     .         awrit6,bessl2,clebsh,config,daxpy,dcopy,dfclos,dpadd,
     .         dpcopy,dpscop,dpzero,dvset,fclose,fclr,fdpp,getpr,getqvc,
     .         gintsr,gtpcor,icopy,iinit,info2,info5,mpint,newrho,phidx,
     .         poiss0,poppr,potpar,prrmsh,pshpr,query,r8tos8,radgra,
     .         radmsh,radmwt,radwgt,rdeq,rgrme,rhocor,rmesh,rseq,rsq1,
     .         rxi,savvxc,setcc,soprm,tcn,tcx,v0intr,vxc0sp

C --- Setup ---
C     lscf: T, make sphere self-consistent given moments
C           F, only calculate charges qc, qtot and amgm

      call tcn('atscpp')
      call getpr(ipr)
      call r8tos8(dclabl(ic),clabl)
      call dpzero(pzloc,nl)

      imake0 = mod(mode,10)
      imake1 = mod(mode/10,10)
      imake2 = mod(mode/100,10)
      imake3 = mod(mode/1000,10)
      if (imake2 > 0) then
        if (imake0 == 1)
     .    call rx('inconsistent switches, "Make sphere potential and Use given potential"')
        if (imake0 == 3) then
          call info0(31,0,0,' Suppress generating potential from P,Q')
          imake0 = 2
        endif
      endif

      lrel   = mod(nglob('lrel'),10)
      lscf   = imake0 == 1 .or. imake0 == 3
      lneedv = imake0 == 2 .or. imake0 == 4
      lfree  = IAND(s_ctrl%lasa,8) /= 0
      lgdd   = isw(.not. IAND(s_ctrl%lasa,128) /= 0)
      if (IAND(s_ctrl%lasa,4096) /= 0) lgdd = lgdd + 4
      lmpol  = IAND(s_ctrl%lasa,32) /= 0
      lso    = IAND(s_ctrl%lncol,4) /= 0
      lbf    = IAND(s_ctrl%lncol,8) /= 0
      lrell  = IAND(mod(s_ctrl%lrel,10),-1) /= 0
      lfrz   = isw(IAND(s_ctrl%lcd,1) /= 0)
      if (lfrz == 0) lfrz = 2*isw(IAND(s_ctrl%lcd,16) /= 0)
      i = mod(s_ctrl%lscr,100)
      lintra = i >= 10 .and. i < 90 .or. mod(s_ctrl%lsx/10,2) /= 0
      loptc  = s_ctrl%loptc > 0
      call dpzero(bhat,3)
      nrmix = s_ctrl%smalit
      idmod = s_spec(is)%idmod
      lmx = s_spec(is)%lmxa
      z = s_spec(is)%z
C     This is now handled by LDA+U
C      idu = s_spec(is)%idu
C      uh = s_spec(is)%uh
C      jh = s_spec(is)%jh
      call iinit(idu,4)
      a = s_spec(is)%a
      nr = s_spec(is)%nr
      etot = 0
      amgm = 0
      sumev = 0
      call dpzero(exc,2)
      havedc = .false.
      havepp = mod(s_ctrl%initc(ic)/2,2)
      nitmax = nrmix(1)
      nmix = nrmix(2)
      thrpv = 0
      rmx = rmax(ic)
      if (lfree) then
        rmx = 50d0
        if (z < 10) rmx = 25
        if (z <= 6) rmx = 20
        if (nsp == 2 .and. nmix <= 10) nmix = 0
      endif
      nlspic = nl*nsp*(ic-1)
      lxcg = nglob('lxcf')/100 /= 0
      intopt = 10*nglob('lrquad')
      allocate(rofi(nr,2))
      call rmesh(z,rmx,isw(lrell),isw(lxcg),nrmx,a,nr)  ! mesh spacing, number of points
      call radmsh(rmx,a,nr,rofi)
      call radwgt(intopt,rmx,a,nr,rofi(1,2))

      ifi = fopna(clabl,30,0)
! Allocations are split to make it easier to find uninitialised arrays with valgrind
      allocate(v(nr,nsp)); v(1,1) = 0
      allocate(rho(nr,nsp))
      allocate(rhoi(nr,nsp))
      allocate(rhoc(nr,nsp))
      allocate(g(nr*2))
      allocate(gp(nr*2*4))

! Temporary measure, see where rho and rhoc are to be initialised properly.
      call dpzero(rho,size(rho)); call dpzero(rhoc,size(rhoc))

      haveso = mod(initc/4,2) == 1
      haveva = mod(initc/8,2) == 1
      havemp = mod(initc/16,2) == 1
      havegr = mod(initc/32,2) == 1
      job = 'gue'
      havev = .false.
      if (mod(imake1,4) > 0) then
        if (associated(s_pot%v0)) then
          havev = s_pot%v0(1,ic) /= NULLI
        endif
        if (havev) then
          idmodl(1:2) = shape(s_pot%v0)
          if (idmodl(1) < nr*nsp) call rx('s_pot%%v0 is improperly dimensioned')
          bhat = 0
          call dcopy(nr*nsp,s_pot%v0(1,ic),1,v,1)
        endif
      endif
      if (.not. havev .and. mod(imake1,2) == 0) then
        havev = aiopot(nr,nsp,a,rmx,bhat,v,ifi)
      endif

      if (imake3 /= 0) call rx('atscpp not ready for imake3')

C     if no bhat read in, flag that it is not there
      if (ddot(3,bhat,1,bhat,1) == 0) bhat(1) = NULLI
      if (havev) job = 'pot'
      call awrit0('%x     ... available from input:',outs,80,0)
      if (mod(initc,2) == 1) call awrit0('%a  p,q',outs,80,0)
      if (mod(initc/2,2) == 1) call awrit0('%a  ppar',outs,80,0)
      if (haveso) call awrit0('%a  so-par',outs,80,0)
      if (havegr) call awrit0('%a  rgrad-me',outs,80,0)
      if (haveva) call awrit0('%a  vintra',outs,80,0)
      if (havemp) call awrit0('%a  mp-par',outs,80,0)
      if (havev) call awrit0('%a  V',outs,80,0)
      swc = aiocor(nr,nsp,a,rmx,rhoc,sumec,sumtc,ifi)
      if (swc) call awrit0('%a  core',outs,80,0)
C     sw=T flags that some required information is missing, for printout
      sw = (lneedv .and. .not. havev) .or. (lscf .and. mod(initc,2) /= 1)
      i = -lgunit(1)
      call awrit0('%x ATSCPP: class '//clabl,outs1,80,0)
      if (imake0 == 0) call awrit0('%a   reading double-counting terms from GEN',outs1,80,0)
      if (lscf) call awrit0('%a  making sphere potential from P,Q',outs1,80,0)
      if (lneedv) call awrit0('%a  reading sphere potential from file',outs1,80,0)
      if (ipr >= 50 .or. sw) call awrit0('%a',outs1,-80,i)
      if (ipr >= 50 .or. sw) call awrit0('%a',outs,-80,i)
      if (lscf .and. mod(initc,2) /= 1) call rx('atscpp: missing P,Q for class '//clabl)
      if (imake0 == 2 .and. .not. havev) call rx('atscpp: missing potential for class '//clabl)
      if (lfrz == 1 .and. .not. swc) call rx('atscpp: missing core for frozen core '//clabl)
      if (lfree .and. imake0 /= 1) call rx('atscpp: lfree=T compatible only with imake=1')
      if (imake0 == 4 .and. .not. havev) return

      call dcopy(3*nl*nsp,qnu(1,1,1,ic),1,qnuloc,1)
      call dcopy(1*nl*nsp,pnu(1+nlspic),1,pnuloc,1)
      qnurl(1,1,1,1,1) = NULLI
      if (lrel == 2 .and. s_ctrl%stpnrho) call dcopy(4*8*nl*nl,s_pot%qnur(1,ic),1,qnurl,1)
C     Align relativistic enu,C to ves=0 unless imake0=0 ()
      if (imake0 /= 0) then ! If ppars are not touched, do not modify
        call shftppr(lrel,1,nl,nsp,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),ves(ic),ves(ic),.true.,.false.)
      endif

C ... Choose best noncollinear quantization axis available from ASA moments.
C     The ASA does not contain full information of the density, but the
C     average projection of magnetization onto each l is accomplished from
C     the output density matrix (routine amagn2).  Thus the self-
C     consistent qnu contain an m-averaged orbital-dependent magnetization.
C     Note: qnuloc is not necessarily supplied (imake0=0,2) in which case
C     amag is nonsensical.  But it is only used for printout.
      amomnc = 0
      if (nsp == 2 .and. neul > 0 .and. imake0 <= 4) then
        call pshpr(1)
        call amagnc(1,nl,1,xx,1,qnuloc,eula,neul,3,amag,amomnc,bxc(1,ic))
        call poppr
C        if (ddot(3,bxc(1,ic),1,bxc(1,ic),1) == 0) then
C          call dcopy(3,bxcl,1,bxc(1,ic),1)
C        endif
      endif

C      j = 23
CC     Case bhat is input
C      if (imake0 == 4) j = 21
CC     Rigid spin: no need to project moments onto Bxc
C*     if (neul == 1) j = j-10
C      j = j-10
CC     Collinear : only copy qnu to qnuloc
C      if (neul <= 0) j = 10
C      if (cmdopt('--keepsignm',11,0,outs)) j = j+100
C      call asprjq(j,clabl,nl,nsp,eula,neul,pnu(1+nlspic),
C     .  qnu(1+3*nlspic),pnuloc,qnuloc,bxc(1,ic),amomnc)

C     Extract core hole parameters
      call gtpcor(s_spec,is,kcor,lcor,qcor)

      scalsoc = 1d0
c     if (ic == 2) scalsoc = 1d-5

C --- Create potential from P,Q --
      if (lscf) then  ! imake0 = 1 or 3

        if (lrel == 2 .and. .not. s_ctrl%stpnrho) call rx('(abort) relativistic moments not available')
        call getqvc(nsp,nl,lmx,z,pnuloc,qnuloc,pzloc,0,0,kcor,lcor,qcor,qc,qtot(ic),amgm,0d0,0d0)
        qrel = 0 ; if (lrel == 2) qrel = sum(qnurl(1,:,:,1,1)+qnurl(1,:,:,2,2)) - (z-qc)
        ec(1) = 0
        call awrit3('%xATOM='//clabl//'%10pZ=%d  %?;n==1;frz;Qc=%d;',outs,100,0,z,lfrz,qc)
        call awrit8('%a  R=%1,6;6d  Qv=%1;6d  %?!n==2!mom=%1;5d!%0d!%a  %?!n==2!Qrel=%1;5d!%j!%a  a=%d  nr=%i',
     .    outs,100,0,rmx,qtot(ic),nsp,amgm,lrel,qrel,a,nr)
        if (nsp == 2 .and. bscal /= 1) call awrit1('%a  bscal=%1;6d',outs,100,0,bscal)
        do  j = 1, 2
          if (ipr >= 20 .or. j == 2) call awrit0('%a',outs,-100,-lgunit(j))
        enddo
C       Self-consistent potential with boundary condition ves(rmt)=0
        call atomsc(lgdd,nl,nsp,lmx,z,rhozbk,kcor,lcor,qcor,rmx,a,nr,rofi,ec,ev,pnuloc,qnuloc,qnurl,
     .    pzloc,bscal,idmod,v,rhoi,rho,rhoc,nmix,qc,sumec,sumtc,sumev,ekin,utot,rhoeps,etot,amgm,
     .    rhrmx(ic),vrmax(1,ic),qtot(ic),exc,job,nitmax,lfrz,scalsoc)
        havev = .true.
        if (lrel == 2) call dcopy(4*8*nl*nl,qnurl,1,s_pot%qnur(1,ic),1) ! update since qnur(4) might have been copied from enu
C       Not here ... pass qnurl to potpar and use either qnurl or pprel depending on idmod ??
C       Or generalize qnu2enu to copy to pprel depending on idmod.
C        if (lrel == 2) then
C          call qnu2enu(3,nl,lmx,pp,0d0,qnurl,pprel(1,1,1,1,1,ic)) ! Copy missing qnur(4) from pprel(5)
C        endif

C       call prrmsh('potential',rofi,v,nr,nr,2)
        call asadc(ic,nr,nsp,z,rofi,a,rofi(1,2),v,rhoc,rho,bscal,rhov,bmval,xcore,rhoeps,rhomu,utot,qval)
        eterms(3)  =  utot
        eterms(6)  =  rhoeps
        eterms(7)  =  rhomu
C       Suppress qc * ves(ic) to xcore and sumec; else add to both
        eterms(8)  =  sumec + ves(ic)*qc*1
        eterms(10) =  xcore + ves(ic)*qc*1
        eterms(9)  =  sumtc
C       Add qval * ves(ic) to rhov, since not incl. in asadc
        eterms(11) =  rhov + ves(ic)*qval
C       NB: in noncollinear case, this is moment || B
        eterms(15) =  amgm
C       Add qval * ves(ic) to sumev, since not incl. in atomsc
        eterms(16) =  sumev + ves(ic)*qval
C       Magnetic contribution to rhov
        eterms(19) =  bmval
        havedc = .true.
        if (lfrz /= 1) swc = .true.

C --- Compute double-counting terms for HK total energy ---
      else if (imake0 == 4) then
        call pshpr(1)
        call getqvc(nsp,nl,lmx,z,pnuloc,qnuloc,pzloc,0,0,kcor,lcor,qcor,qc,qtot(ic),amgm,0d0,0d0)
        ec(1) = 0
        call atomsc(lgdd,nl,nsp,lmx,z,rhozbk,kcor,lcor,qcor,rmx,a,nr,rofi,ec,ev,pnuloc,qnuloc,qnurl,pzloc,
     .    bscal,idmod,v,rhoi,rho,rhoc,nmix,qc,sumec,sumtc,sumev,ekin,utot,rhoeps,etot,amgm,
     .    xx(1),xx(2),qtot(ic),exc,job,0,lfrz,scalsoc)
        call asadc(ic,nr,nsp,z,rofi,a,rofi(1,2),v,rhoc,rho,bscal,rhov,bmval,xcore,rhoeps,rhomu,utot,qval)
        call poppr
        eterms(3)  =  utot
        eterms(6)  =  rhoeps
        eterms(7)  =  rhomu
C       Suppress qc * ves(ic) to xcore and sumec; else add to both
        eterms(8)  =  sumec + ves(ic)*qc*1
        eterms(10) =  xcore + ves(ic)*qc*1
        eterms(9)  =  sumtc
C       Add qval * ves(ic) to rhov, since not incl. in asadc
        eterms(11) =  rhov + ves(ic)*qval
        eterms(15) =  amgm
C       Magnetic contribution to rhov
        eterms(19) =  bmval
        havedc = .true.
C   ... Shift enu and c by crystal electrostatic potential
        call shftppr(lrel,1,nl,nsp,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),ves(ic),ves(ic),.false.,.true.)
        goto 99
C --- Else, try and read atomic parameters by aiogen ---
      else  ! imake0 = 0 or 2
        call getqvc(nsp,nl,lmx,z,pnuloc,qnuloc,pzloc,0,0,kcor,lcor,qcor,qc,qtot(ic),amgm,0d0,0d0)
        havedc = aiogen(lbl,xx,xx,nn,nn,nn,nn,xx,xx,xx,vrmax(1,ic),
     .    sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
        if (havedc) then
          qval = z - qc + qtot(ic)
          eterms(3)  =  utot
          eterms(6)  =  rhoeps
          eterms(8)  =  sumec + ves(ic)*qc*1
C         Add qval * ves(ic) to sumev, since not incl. in atomsc
          eterms(16) =  sumev + ves(ic)*qval
C         Not enough info to resolve xcore, rhov; stuff into rhov
          eterms(10) =  0 + ves(ic)*qc*1
          eterms(11) =  eterms(16) + sumec - ekin
        endif

C       vrmax from v(r) takes precedence, if it's available
        if (havev) then
          call dpscop(v,xx,1,nr,1,1d0)
          if (nsp == 2) then
            call dpscop(v,xx,1,2*nr,2,1d0)
            vrmax(1,ic) = (xx(1)+xx(2))/2 - 2*Z/rmx
            vrmax(2,ic) =  xx(1)-xx(2)
          else
            vrmax(1,ic) = xx(1) - 2*Z/rmx
            vrmax(2,ic) = 0
          endif
        endif
      endif

C --- Make sphere potential parameters ---
      avvxc = 0
      if (imake0 == 2 .or. imake0 == 3) then

        call icopy(nl,idmod,1,idmodl,1)
        if (IAND(s_ctrl%lasa,4096) /= 0) then
          forall ( i=1:nl) idmodl(i) = idmodl(i)+100
        endif

C        if (cmdopt('-elin=',6,0,outs)) then
C          i = 6
C          call rxx(.not. a2bin(outs,elin,4,0,' ',i,len(outs)),
C     .      'atomsr: failed to parse'//outs)
C          do  15  i = 0, nl*nsp-1
C            if (mod(idmod(1+mod(i,nl)),10) == 0) then
C              pp(1,i+1,1,ic) = elin
C              idmodl(mod(i,nl)+1) = 2
C            endif
C   15     continue
C        endif

C   --- Average exchange-correlation field ---

        if (nsp == 2) then
          intopt = 10*nglob('lrquad')
          call savvxc(nr,rho,rhoc,v,rofi(1,2),avvxc)
        endif

C   --- Second-generation potential parameters ---
C       Calculated for ves=0 on sphere surface
        call radmsh(rmx,a,nr,rofi)
        call setcc(lrell,1d0)
        if (lrel == 2) pprel(1:4,:,:,:,:,ic) = 0

        call potpar(nl,nsp,lmx,z,rmx,avw,ekap,lso.or.lbf.or.neul > 1,loptc,lmpol,a,nr,rofi,v,pnuloc,
     .    idmodl,exc,qnuloc,qnurl,idu,uh,jh,facso,thrpv,thrpvl,g,gp,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),
     .    sop(0,1,1,1,ic),pmpol(1,1,1,1,1,ic),rgrad(1,1,1,1,ic),t,s_pot%gfr(1:,ic))

        call setcc(lrell,1d0)
        havepp = .true.
        haveso = haveso .or. (lso.or.lbf.or.neul > 1)
        havegr = havegr .or. loptc
        havemp = havemp .or. lmpol
C   ... Shift enu and c by crystal electrostatic potential
        call shftppr(lrel,1,nl,nsp,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),ves(ic),ves(ic),.false.,.true.)

C   ... Second-order hamiltonian: set p^gamma = pph(4) = 0
        if (IAND(s_ctrl%lham,3) /= 0) then
          pp(4,:,:,ic) = 0
        endif

C   --- Intraatomic Coulomb d^2E/dq_i dq_j ---
        if (lintra) then
          call v0intr(nl,nsp,lmx,z,rhozbk,rmx,a,nr,rofi,pnuloc,qnuloc,pzloc,
     .      bscal,v,rhoi,rho,rhoc,g,gp,nmix,nitmax,qc,lfrz,avw,ekap,2,vintra(1,1,ic))
          haveva = .true.
        endif
      endif

C --- Reopen atomic file and write atomic data ---
      if (imake0 == 1) then
        initc = 1
      elseif (imake0 /= 0) then
        call dfclos(ifi)
        ifi = fopna(clabl,30,0)
        if (havedc) sw = aiogen(clabl,z,rmax(ic),lmx,nsp,lrell,nr,a,qc,qtot(ic),
     .                          vrmax(1,ic),sumec,sumev,thrpv,ekin,utot,rhoeps,etot,-ifi)
        j = 1; if (lrel == 2) j = ic
        if (mod(initc,2) == 1) then
          sw = aiomom(clabl,pnu(1+nlspic),qnu(1,1,1,ic),s_pot%qnur(1,j),idmod,
     .      nl,lmx,nsp,z,rhrmx(ic),vrmax(1,ic),(/0d0,0d0,1d0/),-ifi)
        endif
        if (havepp) sw = aiopar(clabl,lrel,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),ves(ic),nl,lmx,nsp,-ifi)
        if (haveso) sw = aiosop(clabl,sop(0,1,1,1,ic),nl,lmx,nsp,-ifi)
        if (havegr) sw = aiorme(clabl,rgrad(1,1,1,1,ic),nl,nsp,-ifi)
        if (haveva) sw = aiova(clabl,vintra(1,1,ic),nl,lmx,nsp,-ifi)
        if (havemp) sw = aiomp(clabl,pmpol(1,1,1,1,1,ic),nl,2*nl-2,nsp,-ifi)
        if (havev)  sw = aiopot(nr,nsp,a,rmax(ic),bhat,v,-ifi)
C       Copy potential to s_pot%v0
        if (havev .and. associated(s_pot%v0)) then
          idmodl(1:2) = shape(s_pot%v0)
          if (idmodl(1) < nr*nsp) call resizev0(s_pot,nr*nsp)
          call dcopy(nr*nsp,v,1,s_pot%v0(1,ic),1)
        endif
        if (swc) sw = aiocor(nr,nsp,a,rmax(ic),rhoc,sumec,sumtc,-ifi)
        call fclose(ifi)
        initc = 1  ! Moments available
        if (havepp) initc = initc+2  ! pp's available
        if (haveso) initc = initc+4
        if (haveva) initc = initc+8
        if (havemp) initc = initc+16
        if (havegr) initc = initc+32
      endif

C --- Printout of atomic parameters ---
C  20 continue
      if (ipr > 30 .and. imake0 /= 1) then
        write(lgunit(1),'(1x)')
        do  j = 1, 2
         if (havedc) then
           call awrit4(' v_rmax= %,6d%23petot= %,6d%?!n!%4f<v_xc>= %,6d',
     .       ' ',80,lgunit(j),vrmax(1,ic),etot,isw(avvxc /= 0d0),avvxc)
           if (ipr > 31 .and. lscf) then
             call awrit3(' thrpv=  %,6d%23pby l:%n:1,6d',' ',80,lgunit(j),thrpv,lmx+1,thrpvl(1))
           endif
         endif
         write(lgunit(j),'(1x)')
         if (.not. lfree .and. havepp) sw = aiopar(clabl,lrel,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),ves(ic),nl,lmx,nsp,-lgunit(j))
         if (haveso) sw = aiosop(clabl,sop(0,1,1,1,ic),nl,lmx,nsp,-lgunit(j))
         if (haveva) sw = aiova(clabl,vintra(1,1,ic),nl,lmx,nsp,-lgunit(j))
         if (havemp) sw = aiomp(clabl,pmpol(1,1,1,1,1,ic),nl,2*nl-2,nsp,-lgunit(j))
        enddo
      endif

C --- Cleanup ---
   99 continue
      call fclr(clabl,ifi)
      deallocate (v,rofi,rho,rhoi,rhoc,g,gp)
      call tcx('atscpp')

      end
      subroutine atomsc(lgdd,nl,nsp,lmax,z,rhozbk,kcor,lcor,qcor,rmax,a,nr,rofi,ec,ev,pnu,qnu,qnur,pz,
     .  bscal,idmod,v,rhoin,rho,rhoc,nrmix,qc,sumec,sumtc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,
     .  vrmax,qtot,exrmax,job,niter,lfrz,scalsoc)
C- Makes an atomic sphere self-consistent and get atomic charges
C ----------------------------------------------------------------
Ci Inputs
Ci   lgdd  :0 Add q2 [phidot**2 - p phi**2] to rho
Ci         :1 Add q2 [phidot**2 + phi*phidotdot] to rho
Ci         :Both produce the same integrated density; see Remarks.
Ci         :lgdd=0 is standard and follows Stuttgart conventions.
Ci         :4 Add 4 to make (phidot,phidotdot) by integrating radial SE outward only
Ci   nl    :dimensions pnu,qnu
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmax  :maximum l for this sphere
Ci   z     :nuclear charge
Ci   rhozbk:constant nuclear background density (jellium) added to z
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   rmax  :potential, density calculated to rmax
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Ci   qnu   :energy moments of charge (see Remarks)
Ci   qnur  :relativistic energy-weighted moments of the sphere charges, kmu repsn
Ci   bscal :factor to scale the effective xc magnetic field (for constr.B)
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
Ci   pz    :pnu for semicore states
Ci   v     :spherical potential (job='pot', else input v not used)
Ci   rhoin:input spherical charge density times 4*pi*r*r (job='rho', else input rhoin not used).
Ci         :Used internally as a work array; DESTROYED on output
Ci   nrmix :number of prior densities to mix to accelerate convergence
Ci         :to self-consistency (Anderson mixing); see Remarks.
Ci         :nrmix<0, linear mixing, with mixing beta |nrmix/100|
Ci   job   :(see Remarks)
Ci         :job='pot': start with potential
Ci         :job='rho': start with rhoin.
Ci         :job='gue': start with guessed charge density.
Ci   niter :number of iterations to attempt convergence to
Ci          self-consistency (see Remarks)
Ci         :If V is given, use niter=0 to generate rho for that V
Ci   lfrz  :0 for soft core
Ci         :1 for frozen core
Ci         :2 for soft core, spin-averaged potential
Cio Input/Outputs
Cio  ec    :core eigenvalues.  On input, these are guessed values.
Cio        :if input ec(1) = 0, atomsc makes an internal choice fo ec.
Cio  ev    :valence eigenvalues  On input, these are guessed values.
Cio        :if input ec(1) = 0, atomsc makes an internal choice for ev.
Co Outputs
Co   rofi  :dimensioned (nr,2).
Co         :(*,1) radial mesh points
Co         :(*,2) weights
Co   rho   :spherical charge density = 4 pi r**2 rhotrue
Co   rhoc  :core charge density (unchanged if lfrz=1)
Co   qc:   :core electronic charge
Co   sumec :core single particle energy
Co   sumtc :core kinetic energy (unchanged if lfrz=1)
Co   ekin  :total kinetic energy
Co   utot  :electrostatic energy
Co   rhoeps:exchange-correlation energy
Co   etot  :sphere total energy
Co   amgm  :difference between spin up and spin down charge
Co   rhrmx :true density at rmax
Co   vrmax :true potential at rmax
Co   qtot  :net charge in sphere
Co   exrmax:exchange-correlation energy at rmax
Cr Remarks
Cr   Boundary conditions pnu, moments qnu, and the electrostatic
Cr   potential at rmax are all that is required to uniquely determine
Cr   a self-consistent spherical charge and potential.  'job'
Cr   determines how the starting potential is made, but the final
Cr   potential  and density should be independent of the initial choice.
Cr
Cr   atomsc uses the boundary condition that the potential at rmax
Cr   is zero, regardless of the total net charge inside the sphere.
Cr   See subroutine madpot for discussion of how this choice affects
Cr   the Madelung energy.
Cr
Cr   Sphere total energy is sum of K.E., Hartree energy, XC energy:
Cr      etot = ekin + utot + rhoeps
Cr   The kinetic energy is computed via double-counting terms
Cr     ekin = sumev + sumec + dsumec - rhovh - rhomu
Cb Bugs
Cb   Total energy terms need to be cleaned up and simplified.
Cu Updates
Cu   23 Mar 17 New option (s_ctrl%lasa,4096) to compute phi,phidot
Cu             from outward radial integration only
Cu   03 Jun 14 Added semilocal orbitals
Cu   28 Jan 13 core can be calculated with spin-averaged potential
Cu   25 Apr 12 (K Belashchenko) new bscal
Cu   10 Apr 12 Repackaged radial mesh quadrature (radwgt)
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character job*3
      integer nr,nsp,nl,nrmix,niter,kcor,lcor,lgdd,lmax,lfrz,idmod(0:nl-1)
      integer, parameter :: ncmx=200, nvmx=20
      double precision ec(ncmx),ev(nvmx),rofi(nr,2),v(nr,nsp),exc(nr),
     .  rho(nr,nsp),rhoc(nr,nsp),rhoin(nr,nsp),pnu(nl,nsp),pz(nl),
     .  qnu(3,nl,nsp),bscal,z,rmax,a,qc,vrmax(2),exrmax(2),rhrmx,
     .  rhozbk,qcor(2),amgm,ekin,etot,qtot,rhoeps,sumec,sumev,sumtc,
     .  utot,scalsoc,qnur(4,0:nl-1,2*nl,2,2)
C ... Local parameters
C     logical :: debug=.false.
      logical last,ltmp
      integer ir,isp,iter,jmix,lrel,lrel2,ipr1,l,i,k,imu,stdo,intopt,ninout(nvmx)
      character strn*10
      double precision b,decay,dl,dold,drdi,drho,dsumec,ea,fac,pi,
     .  rho0t,rhomu,rhovh,rmsdel,ro,sumq,tolrsq,vhrmax,vnucl,vrhoc,vsum,
     .  zvnucl,rvh(2),rho0(2),reps(2),rmu(2),sec(2),stc(2),sev(2),sfac(4),xv
      real(8), allocatable :: rhom(:,:,:)
C     for Anderson mixing:
      integer,parameter :: NULLI=-99999
      integer nmix,lxcf
      double precision norm(10,10),awk(10,2),beta,beta1
C     tolch is tolerance for change in the charge density, tolv for rseq
      double precision tolch,tl
      parameter (tolrsq=1d-12)
      procedure(logical) :: cmdopt
      procedure(integer) :: iprint,nglob,amix
      procedure(real(8)) :: dasum,dlength,ddot

      sfac = 0

      if (lmax >= nl) call rx('atomsc:  lmax too large')
      lrel = mod(nglob('lrel'),10)
      lrel2 = nglob('lrel')
      tolch = 5d-5 ; if (lrel == 2) tolch = 2d-7
      stdo = nglob('stdo')
      pi = 4d0*datan(1d0)
      vnucl = 0; vrhoc = 0;
      b = rmax/(dexp(a*(nr-1)) - 1)
      vhrmax = 0d0
      nmix = min(max(nrmix,0),10)
      sec = 0d0; stc = 0d0
      if (allocated(rhom)) deallocate(rhom)
      allocate(rhom(nr*nsp,0:nmix+1,2))
      intopt = 10*nglob('lrquad')
      lxcf = nglob('lxcf')
      do  i = 1, nsp
      do  l = 0, lmax
        if (int(pnu(l+1,i)) /= int(pnu(l+1,nsp)))
     .    call rx('atomsc: spin mismatch in pnu')
      enddo
      enddo

C --- Core charge, radial mesh points and weights ---
      if (kcor /= 0) then
      if (qcor(1) /= 0 .or. qcor(2) /= 0) then
        call info5(30,0,0,'%9fAdd core hole:  kcor=%i  lcor=%i'//
     .    '  qcor=%d  amom=%d',kcor,lcor,qcor,qcor(2),0)
      endif
      endif
      call getqvc(nsp,nl,lmax,z,pnu,qnu,pz,0,0,kcor,lcor,qcor,qc,qtot,amgm,ec,ev)
C ... Guesses for core and valence eigenvalues
      if (ec(1) == 0) call getqvc(nsp,nl,lmax,z,pnu,qnu,pz,ncmx,nvmx,
     .                              kcor,lcor,qcor,qc,qtot,amgm,ec,ev)
      call radmsh(rmax,a,nr,rofi)
      call radwgt(intopt,rmax,a,nr,rofi(1,2))

C --- Moments printout ---
      if (iprint() >= 20) then
        ltmp = dasum(lmax+1,qnu(2,1,1),3)+dasum(lmax+1,qnu(3,1,1),3) == 0
        if (ltmp) then
          call awrit5('  Pl= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
     .      ' ',120,stdo,lmax+1,pnu,nsp,lmax+1,pnu(1,nsp))
          call dcopy(lmax+1,qnu,3,awk,1)
          call dcopy(lmax+1,qnu(1,1,nsp),3,awk(1,2),1)
          call awrit5('  Ql= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
     .      ' ',120,stdo,lmax+1,awk,nsp,lmax+1,awk(1,2))
        elseif (iprint() >= 30) then
          write (stdo,1)
    1     format('   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
          do  isp = 1, nsp
            do  l = 0, lmax
              dl = dtan(pi*(.5d0 - pnu(l+1,isp)))
              if (dabs(dl) > 9999) dl = 0
              write (stdo,2) l,pnu(l+1,isp),(qnu(i,l+1,isp),i=1,3),idmod(l),dl
    2         format(i4,4F13.7,i4,f13.7)
            enddo
          enddo
        endif

        if (lrel == 2 .and. iprint() >= 30) then
        write(stdo,22) ! sum(qnur(1,:,:,1,1)+qnur(1,:,:,2,2))
   22   format(' Relativistic moments'/
     .    '   l  mu   ms1 ms2',5x,'q0',11x,'q1',11x,'q2',11x,'enu')
        do  l = 0, lmax
          do  imu = 1, 2*(l+1)
            fac = dble(imu-l) - 1.5d0
            do  i = 1, 2
              do  k = 1, 2
                if (dlength(3,qnur(1,l,imu,i,k),1) > 1d-6) then
                  jmix=4; if (qnur(4,l,imu,i,k) == NULLI) jmix=3
                  write(stdo,23) l,fac,i,k,(qnur(ir,l,imu,i,k), ir=1,jmix)
                endif
              enddo
            enddo
          enddo
        enddo
   23   format(i4,f5.1,2i4,4f13.7)
        endif
      endif

C --- Initial charge density ---
      if (job == 'pot') then
C       call prrmsh('v for newrhorel',rofi,v,nr,nr,nsp)
        call newrhorel(z,lrel2,lgdd,nl,1,lmax,a,b,nr,rofi,v,rhoin,rhoc,kcor,lcor,qcor,
     .    pnu,qnu,qnur,pz,scalsoc,sec,stc,sev,ec,ev,sfac,ninout,tolrsq,nsp,lfrz,000)
C       call prrmsh('starting rho-rhoc for atom ',rofi,rhoin-rhoc,nr,nr,nsp)
        if (niter == 0) then
          call dcopy(nr*nsp,rhoin,1,rho,1)
!! --- Removed 'return' to initialised properly vars other than rho (sumec and sumtc).
!           return
        endif
      else if (job == 'gue') then
        decay = 1d0+z/10d0
        decay = dmin1(decay,5d0)
        decay = 5
C         write(6,335) decay
C  335   format(/' initialize exponential density with decay',f7.3)
        sumq = 0d0
        do  ir = 1, nr
          ro = dexp(-decay*rofi(ir,1))*rofi(ir,1)**2
          rhoin(ir,1) = ro
          sumq = sumq + ro*a*(rofi(ir,1)+b)
        enddo
        fac = z/(sumq*nsp)
        forall (ir = 1:nr) rhoin(ir,1) = rhoin(ir,1)*fac
        if (nsp == 2) then
          forall (ir = 1:nr) rhoin(ir,2) = rhoin(ir,1)
        endif
      else if (job /= 'rho') then
        call rx('atomsc: job not pot|rho|gue')
      endif
      call dpscop(rhoin,rhom(1,0,2),nr*nsp,1,1,1d0)

C --- Start self-consistency loop ---
      drho = 100d0
      last = .false.
      beta = .3d0 ; if (lrel == 2) beta = .2d0

      if (iprint() >= 30) write (stdo,3)
 3    format(/'  iter     qint',9x,'drho',10x,'vh0',10x,'rho0',10x,'vsum',5x,'beta')
      jmix = 0
      beta1 = beta
      if (nrmix < 0) beta1 = dble(-nrmix)/100

      rvh = 0
      reps = 0
      rmu = 0

      do  iter = 1, niter
        dold = drho; if (iter == 1) dold = 1

C       call prrmsh('rho generating v',rofi,rhoin-rhoc,nr,nr,nsp)

        tl = tolrsq
C       Hartree potential
        call addzbk(rofi,nr,nsp,rhoin,rhozbk,-1d0)
        call poiss0(z,a,b,rofi,rhoin,nr,vhrmax,v,rvh,vsum,nsp)

CC       debugging
CC       Let n(r) = true density.  Then rhoin(r) = 4*pi*r**2*n(r)
C
CC       Special test for constant density
CC       Substitute constant density of unit charge. Q = 1 = int [4*pi*r**2 n(r)] = int [rhoin]
CCC      Let rhoin(r) = A*r**2.  A fixed by int rhoin = 1 => A=3/R^3
C        forall (ir = 1:nr) rhoin(ir,1) = rofi(ir,1)**2*3/rofi(nr,1)**3
C        fac = ddot(nr,rhoin,1,rofi(1,2),1)
C        print *, 'charge', fac
C        call poiss0(z*0,a,b,rofi,rhoin,nr,2*fac/rofi(nr,1),v,rvh,vsum,nsp)
CC       call prrmsh('vh atom ',rofi,v,nr,nr,nsp)
CC       lscoul expects rhoL Y_L = r**2*n(r)
CC       rhoin = 4*pi r**2 n(r) = (4*pi) rhoL Y_L = rhoL/Y_L => rhoL = Y_L rhoin
C        fac = 1/dsqrt(4*pi)     ! This is Y0
C        call dscal(nr,fac,rhoin,1)
C        call pshpr(51)
C        call lscoul(rofi,nr,nsp,1,z*0,-.0002d0,rhoin,v)
C        call dscal(nr,fac,v,1) ! convert from v_L to true V
C        call prrmsh('v from lscoul',rofi,v,nr,nr,nsp)

        call addzbk(rofi,nr,nsp,rhoin,rhozbk,1d0)
        vnucl = v(1,1)
C       Exchange-correlation potential
        if (last .and. iprint() >= 50) call pshpr(80)
        reps = 0
        call vxc0sp(lxcf,rofi,rofi(1,2),rhoin,nr,v,exc,rho0,bscal,reps,rmu,nsp)

C       if (debug) then
C         print *, iter
C         call prrmsh('rho atom ',rofi,rhoin,nr,nr,nsp)
C         call prrmsh('vh+vxc atom ',rofi,v,nr,nr,nsp)
C       endif
        exrmax(1) = exc(nr)
        if (last .and. iprint() >= 50) call poppr
C       Get rhrmx, exrmax
        fac = 4*pi*rofi(nr,1)**2
        rhrmx = rhoin(nr,1)/fac
        if (nsp == 2) then
          exrmax(2) = exrmax(1)
          rhrmx = rhrmx + rhoin(nr,2)/fac
        endif
        ipr1 = 0
        if (last .and. iprint() >= 40) ipr1 = 1
        if (last .and. iprint() > 40) ipr1 = 2

        call newrhorel(z,lrel2,lgdd,nl,1,lmax,a,b,nr,rofi,v,rho,rhoc,kcor,lcor,qcor,
     .    pnu,qnu,qnur,pz,scalsoc,sec,stc,sev,ec,ev,sfac,ninout,tl,nsp,lfrz,ipr1)
C        if (debug .and. iter >= 3) then
C          print *, iter
C          call prrmsh('potential used to generate rho',rofi,v,nr,nr,nsp)
C          call prrmsh('valence rho generated by v',rofi,rho-rhoc,nr,nr,nsp)
C          call newrhorel(z,lrel2,lgdd,nl,1,lmax,a,b,nr,rofi,v,rho,rhoc,kcor,lcor,qcor,
C     .      pnu,qnu,qnur,pz,scalsoc,sec,stc,sev,ec,ev,sfac,ninout,tl,nsp,lfrz,ipr1)
C          call prrmsh('valence rho generated by v',rofi,rho-rhoc,nr,nr,nsp)
C        endif

        drho = 0
        sumq = 0
        vrhoc = 0
        rho0t = 0
        do  isp = 1, nsp
          rho0t = rho0t + rho0(isp)
          sumq = sumq + ddot(nr,rofi(1,2),1,rho(1,isp),1)
          do  ir = 1, nr
          drdi = a*(rofi(ir,1) + b)
          drho = drho+rofi(ir,2)/drdi*dabs(rho(ir,isp)-rhoin(ir,isp))
          enddo
          do  ir = 2, nr
            ea = (v(ir,isp)-2*z/rofi(ir,1))
            vrhoc = vrhoc+rofi(ir,2)*ea*rhoc(ir,isp)
          enddo
        enddo
C       call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)
        call dcopy(nr*nsp,rho,1,rhom,1)
        jmix = amix(nr*nsp,min(jmix,nmix),nmix,0,beta1,iprint()-70,
     .              .9d0,norm,awk(1,2),rhom,awk,rmsdel)
        call dpscop(rhom,rhoin,nr*nsp,1+nr*nsp*(nmix+2),1,1d0)
C
        if (last) exit
        if (iprint() >= 41 .or. iprint() >= 30 .and.
     .     (drho < tolch .or. iter == niter-1 .or. iter == 1))
     .      write (stdo,4) iter,sumq,drho,vnucl,rho0t,vsum,beta1
    4   format(i5,f12.6,1p,e12.3,0p,f14.4,e14.4,f14.4,f7.2)
C       call prrmsh('rhol for atom ',rofi,rhom,nr,nr,nl*nsp)
        last = (drho < tolch .or. iter == niter-1)
        if (iprint() > 100)  call query(' ',-1,0)
        jmix = jmix+1
C       Beta for next iteration
        if (nrmix < 0) cycle
        if (lrel == 2) then
          ro = beta1
          beta1 = min(max(dold/drho*ro,ro-.2d0,beta),1d0,ro+.2d0)
          if (nmix > 0 .and. drho < 1)
     .    beta1 = min(max(dold/drho*ro,ro-.1d0,beta),1d0,ro+.2d0)
          if (lrel == 2) beta1 = min(beta1,0.7d0)
        else
          beta1 = min(max((1-drho/dold)/beta1,beta1-.2d0,beta),1d0,beta1+.2d0)
          if (nmix > 0 .and. drho < 1) beta1 = 1
        endif
      enddo                     ! End of iteration loop
      deallocate(rhom)

      if (drho > tolch .and. niter > 0) call logwarn(1,'%N RSEQ warning!  sphere density not converged')
      if (iprint() >= 30) write(stdo,'(1x)')

C --- Collect terms for total energy ---
      rhoeps = 0
      rhomu  = 0
      sumev  = 0
      if (lfrz /= 1) then
        sumec = 0
        sumtc = 0
      endif
      rhovh  = 0
      do  isp = 1, nsp
        if (nsp == 2 .and. iprint() > 30) write (stdo,5) isp,
     .      v(nr,isp)-2*z/rmax,sev(isp),sec(isp),rvh(isp),reps(isp),
     .      rmu(isp)
    5   format(' Spin',i2,':',
     .  /' vrmax=',f12.5,'    sumev= ',f12.5,'    sumec=',f12.5,
     .  /' rhovh=',f12.5,'    rhoeps=',f12.5,'    rhomu=',f12.5)


        rhoeps = rhoeps + reps(isp)
        rhomu = rhomu + rmu(isp)
        sumev = sumev + sev(isp)
        if (lfrz /= 1) then
          sumec = sumec + sec(isp)
          sumtc = sumtc + stc(isp)
        endif
        rhovh = rhovh + rvh(isp)
      enddo
      zvnucl = -z*vnucl
      utot = .5d0*(rhovh + zvnucl)
C ... Correction to core eigenvalues if sumec not obtained from this V
      dsumec = vrhoc - (sumec-sumtc)
      ekin = sumev + sumec + dsumec - rhovh - rhomu
      etot = ekin + utot + rhoeps
      if (iprint() >= 40) write(stdo,6) sumev,sumec,vnucl,rhovh,
     .  zvnucl,utot,rhomu,rhoeps,dsumec,ekin,sumtc,etot
    6 format(/' sumev=',f13.6,'    sumec =',f13.6,'   vnucl =',f13.6
     .       /' rhovh=',f13.6,'    zvnucl=',f13.6,'   utot  =',f13.6
     .       /' rhomu=',f13.6,'    rhoeps=',f13.6,'   dsumec=',f13.6
     .       /' ekin= ',f13.6,'    tcore =',f13.6,'   etot  =',f13.6)
      vrmax(1) = -2*z/rmax
      do  isp = 1, nsp
        vrmax(1) = vrmax(1) + v(nr,isp)/nsp
      enddo
      vrmax(2) = 0
      if (nsp == 2) vrmax(2) = v(nr,1)-v(nr,2)

C --- Write out rho if requested ---
      if (cmdopt('--dumprho',8,0,strn)) then
        call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)
        call prrmsh('pot for atom ',rofi,v,nr,nr,nsp)
        call prrmsh('rhoc for atom ',rofi,rhoc,nr,nr,nsp)
        allocate(rhom(nr*nl*nsp,1,1)); call dpzero(rhom,nr*nl*nsp)

        call setpr(110)
        call newrhorel(z,lrel2,lgdd,nl,nl,lmax,a,b,nr,rofi,v,rhom,rhoc,kcor,lcor,qcor,
     .    pnu,qnu,qnur,pz,scalsoc,sec,stc,sev,ec,ev,sfac,ninout,tl,nsp,lfrz,000)

        call prrmsh('rhol for atom ',rofi,rhom,nr,nr,nl*nsp)
        deallocate(rhom)
      endif

CC#ifdef ALENA
      if (lrel == 2 .and. z /= 0) then
        xv = 0
        if (sfac(4) /= 0.0_8) xv = 2*sfac(4)/sfac(2)
        call info5(40,1,0,' Qv=%;6d  M(S)=%;6d  M(L)=%;6d  M(L+S)=%;6d  g=%;6d',
     .    sfac(1),sfac(2),sfac(3),sfac(4),xv)
      endif
CC#endif

      end
      subroutine addzbk(rofi,nr,nsp,rho,rhozbk,scale)
      implicit none
      integer nr,nsp
      double precision rofi(*),rho(nr,1),rhozbk,scale
      integer ir,isp
      double precision s

      if (rhozbk == 0) return
      s = 16*datan(1d0)*scale*rhozbk
      do  isp = 1, nsp
        do  ir = 2, nr
          rho(ir,isp) = rho(ir,isp)+s*rofi(ir)*rofi(ir)
        enddo
      enddo
      end
      subroutine poiss0(z,a,b,rofi,rho,nr,vhrmax,v,rhovh,vsum,nsp)
C- Hartree potential for spherical rho.
C  ---------------------------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rho   :spherical charge density times 4*pi*r*r
Ci   nr    :number of mesh points
Ci   vhrmax:estat potential at rmax, including nuclear term
Ci         :vhrmax fixes arbitrary constant of integration
Ci         :v(ir) is integrated numerically, then
Ci         :shifted by a constant so that v(nr)-2d0*z/rmax = vhrmax
Ci   nsp   :=1 spin degenerate, =2 non-degenerate
Ci   rofi  :radial mesh points
Co Outputs:
Co   v     :spherical Hartree potential, excluding nuclear term -2d0*z/r
Co         :Thus ves[electron+nuc] = v(nr) - 2d0*z/rmax
Co   vsum  :integral over that potential which is zero at rmax.
Co   rhovh :integral of rho*vh, where vh is the electrostatic
Co         :potential from both electronic and nuclear charge
Cr Remarks:
Cr    Solves Poisson's equation for given spherical charge density
Cr    and a specified value vhrmax at rofi(nr).
Cr    rho =  4*pi*r*r*rhotrue  but  v = vtrue
Cu Updates
Cu   10 Apr 12 Repackaged radial mesh quadrature (radwgt)
Cu             poiss0 no longer makes rofi, but takes rofi as input
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp
      double precision a,b,z,rofi(nr),v(nr,nsp),rho(nr,nsp),rhovh(nsp)
C ... Local parameters
      double precision pi
      double precision rmax,r2,r3,r4,f2,f3,f4,x23,x34,cc,bb,dd,df,drdi,
     .  r,srdrdi,g,f,y2,y3,y4,ro,vnow,a2b4,vhrmax,vsum,vhat0,wt(nr)
      integer ir,isp,intopt,nglob

      intopt = 10*nglob('lrquad')
      pi = 4d0*datan(1d0)
      rmax = rofi(nr)
      call radwgt(intopt,rmax,a,nr,wt)
C     call prrmsh('atomsc: integration weights',rofi,wt,nr,nr,1)

C --- Approximate rho/r**2 by cc + bb*r + dd*r*r near zero  ---
      r2 = rofi(2)
      r3 = rofi(3)
      r4 = rofi(4)
      f2 = 0d0
      f3 = 0d0
      f4 = 0d0
      do  isp = 1, nsp
        f2 = f2 + rho(2,isp)/r2**2
        f3 = f3 + rho(3,isp)/r3**2
        f4 = f4 + rho(4,isp)/r4**2
      enddo
      x23 = (r3*r3*f2 - r2*r2*f3)/(r3 - r2)
      x34 = (r4*r4*f3 - r3*r3*f4)/(r4 - r3)
      cc = (r2*x34 - r4*x23) / (r3*(r2 - r4))
      bb = ((r2+r3)*x34 - (r3+r4)*x23) / (r3*r3*(r4-r2))
      dd = (f2 - bb*r2 - cc)/r2**2

C --- Numerov for inhomogeneous solution ---
      a2b4 = a*a/4d0
      v(1,1) = 1d0
      df = 0d0
      do  ir = 2, 3
        r = rofi(ir)
        drdi = a*(r + b)
        srdrdi = dsqrt(drdi)
        v(ir,1) = v(1,1) - r*r*(cc/3d0 + r*bb/6d0 + r*r*dd/10d0)
        g = v(ir,1)*r/srdrdi
        f = g*(1d0 - a2b4/12d0)
        if (ir == 2) y2 = -2d0*f2*r2*drdi*srdrdi
        if (ir == 3) y3 = -2d0*f3*r3*drdi*srdrdi
        df = f - df
      enddo
      do  ir = 4, nr
        r = rofi(ir)
        drdi = a*(r + b)
        srdrdi = dsqrt(drdi)
        ro = 0d0
        do  isp = 1, nsp
          ro = ro + rho(ir,isp)
        enddo
        y4 = -2d0*drdi*srdrdi*ro/r
        df = df + g*a2b4 + (y4 + 10d0*y3 + y2)/12d0
        f = f + df
        g = f/(1d0 - a2b4/12d0)
        v(ir,1) = g*srdrdi/r
        y2 = y3
        y3 = y4
      enddo

C --- Add constant to get v(nr) = vhrmax ---
      vnow = v(nr,1) - 2d0*z/rmax
      do  ir = 1, nr
        v(ir,1) = v(ir,1) + (vhrmax - vnow)
      enddo

C --- Integral vh and  rho * vh ---
      rhovh(1) = 0d0
      rhovh(nsp) = 0d0
      vsum = 0d0
      vhat0 = 0d0
      do  ir = 2, nr
        r = rofi(ir)
        ro = 0d0
        do  isp = 1, nsp
          rhovh(isp) = rhovh(isp) + wt(ir)*rho(ir,isp)*(v(ir,1)-2d0*z/r)
          ro = ro + rho(ir,isp)
        enddo
        vhat0 = vhat0 + wt(ir)*ro*(1d0/r - 1d0/rmax)
        vsum = vsum + wt(ir)*r*r*(v(ir,1) - vhrmax)
      enddo
      vsum = 4d0*pi*(vsum - z*rmax*rmax)
      vhat0 = 2d0*vhat0 + 2d0*z/rmax + vhrmax
      v(1,1) = vhat0

C --- Copy to second spin channel if spin polarized ---
      if (nsp == 1) return
      do  ir = 1, nr
        v(ir,2) = v(ir,1)
      enddo
      end

C      subroutine fctp0(l,rofi,v,z,nr,nctp0,xrim,xmin,nsave)
CC- Initialize things for FCTP, which finds classical turning point
CC  ---------------------------------------------------
CCi Inputs
CCi   l     :angular quantum number
CCi   rofi  :radial mesh points
CCi   v     :spherical potential (atomsr.f)
CCi   z     :nuclear charge
CCi   nr    :number of radial mesh points
CCo Outputs:
CCo   nctp0:
CCo   xrim:
CCo   xmin:
CCo   nsave: Estimate for nctp
CC  ---------------------------------------------------
C      implicit none
C      integer l,nctp0,nr,nsave
C      double precision xmin,xrim,z,v(nr),rofi(nr)
CC Local
C      integer ir
C      double precision fllp1,r,x,xlast
C      fllp1 = l*(l + 1)
C      r = rofi(10)
C      x = fllp1/r/r - 2d0*z/r + v(10)
C      ir = 10
C   80 ir = ir + 1
C      xlast = x
C      r = rofi(ir)
C      x = fllp1/r/r - 2d0*z/r + v(ir)
C      if (x <= xlast .and. ir < nr) goto 80
C      nctp0 = ir - 1
C      xmin = xlast
C      r = rofi(nr)
C      xrim = fllp1/r/r - 2d0*z/r + v(nr)
C      if (xmin >= xrim-3d0) nctp0 = nr
C      if (xmin >= xrim-3d0) xmin = xrim
C      nsave = (nctp0 + nr)/2
C      end
C      subroutine fctp(e,nctp,nctp0,xrim,xmin,nsave,l,rofi,v,z,nr,a,b)
CC- Finds classical turning point for wave function
CC ----------------------------------------------------------------------
CCi Inputs
CCi   e     :energy
CCi   nctp0 :
CCi   xrim  :
CCi   xmin  :
CCi   nsave: Estimate for nctp
CCi   l     :l quantum number
CCi   rofi  :radial mesh points
CCi   v     :spherical potential
CCi   z     :nuclear charge
CCi   nr    :number of radial mesh points
CCi   a     :mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
CCi   b     :mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
CCo Outputs:
CCo   nctp:  rofi(nctp) is classical turning point
CCo   nsave: New estimate for NCTP, next iteration
CCr Remarks
CC  ---------------------------------------------------
C      implicit none
C      integer l,nctp,nctp0,nr,nsave
C      double precision a,b,e,xmin,xrim,z,v(nr),rofi(nr)
CC Local
C      integer irep,k,n1,n2,nlast,ntry
C      double precision dfdr,dvdr,fllp1,fntry,fofr,r,rtry,vme
C
C      fllp1 = l*(l+1)
C      if (nctp0 == nr) then
C        nctp = nr
C        return
C      endif
C      if (e > xrim) then
C        nctp = nr
C        return
C        endif
C      if (e < xmin) then
C        nctp = 2
C        return
C      endif
C      n1 = nctp0
C      n2 = nr
C      nctp = nsave
C      nlast = -10
C      do  20  irep = 1, 20
C      if (nctp > n2 .or. nctp < n1) nctp=(n1+n2+1)/2
C        r = rofi(nctp)
C        vme = (v(nctp)-e)
C        if (nctp > nr) call rx('fctp0: nctp gt nr')
C        k = min(nctp,nr-1)
C        dvdr = (v(k+1)-v(k-1))/(2.d0*a*(r+b))
CC       dvdr = (v(nctp+1)-v(nctp-1))/(2.d0*a*(r+b))
C        fofr = fllp1/r/r - 2.d0*z/r + vme
C        dfdr = -2.d0*fllp1/r/r/r+2.d0*z/r/r + dvdr
C        rtry = r - fofr/dfdr
C        rtry = dmax1(rtry,rofi(2))
C        fntry = dlog(rtry/b+1.d0)/a + 1.d0
C        ntry = fntry + .5d0
CC|      write(6,810) irep,n1,nctp,n2,fntry,ntry
CC|810   format(i6,'   n1,nctp,n2=',3i5,'   ntry=',f8.3,i6)
C        if (nlast == nctp) goto 98
C        if (fofr > 0.d0) n2=nctp
C        if (fofr < 0.d0) n1=nctp
C        nlast = nctp
C        nctp = ntry
C   20 continue
C   98 if (nctp == nctp0+1) nctp = 2
C      nsave = nctp
C      end

      subroutine fctp0(l,nr,rofi,v,z,nctp0)
C- Initialize things for FCTP, which finds classical turning point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   l     :angular momentum
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci   nr    :number of mesh points
Co Outputs:
Co   nctp0 :minimum of effective potential
Co          or nr if v(nr) > vmin + 3 Ry
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer, intent(in) :: l,nr
      integer, intent(out) :: nctp0
      real(8), intent(in) :: z,v(nr),rofi(nr)
C Local variables:
      integer :: ir

      real(8) :: veff,zz,fllp1,vi,vim1
C Statement functions:
      veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)

      zz = 2*z
      fllp1 = l*l+l

      ir = 10
      vim1 = veff(ir)
      ir = 11
      vi = veff(ir)
      do while (vi <= vim1 .and. ir < nr)
        ir = ir+1
        vim1 = vi
        vi = veff(ir)
      end do
      nctp0 = ir - 1
      vi = veff(nctp0)
      vim1 = veff(nr)
      if (vi >= vim1-3) nctp0 = nr

      end
      subroutine fctp(a,b,e,l,nctp0,nr,rofi,v,z,nctp)
C- Finds classical turning point
C  ---------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   e     :energy
Ci   l     :angular momentum
Ci   nctp0 :minimum of effective potential (fctp0)
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci Inputs/Outputs:
Cio  nctp  :rofi(nctp) is classical turning point
Cu Updates
Cu   13 Jun 00 Added safety checks to guard against jumps in potential
C  ---------------------------------------------------------------------
      implicit none
C Passed variables:
      integer nctp,nctp0,l,nr
      double precision e,rofi(nr),v(nr),z,a,b
C Local variables:
      integer ir,irep,n1,n2,nlast,ntry
      double precision r,veff,fllp1,fofr,dfdr,rtry,zz
C Intrinsic functions:
      intrinsic dlog,dmax1,min0
C Statement functions:
      veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)

      zz=z+z
      fllp1 = l*(l+1)

      if (nctp0 == nr .or. e > veff(nr)) then
        nctp = nr
      elseif (e < veff(nctp0)) then
        nctp = 2
      else
        n1 = nctp0
        n2 = nr-1
        nlast = -10
        do  irep = 1, 20
          if (nctp > n2 .or. nctp < n1) nctp = (n1 + n2)/2
          r = rofi(nctp)
          fofr =  veff(nctp)-e
          dfdr =-(fllp1+fllp1)/r/r/r + zz/r/r +
     .           (v(nctp+1) - v(nctp-1))/(2.d0*a*(r+b))
          rtry = dmax1(r-fofr/dfdr,rofi(2))
          ntry = dlog(rtry/b+1.d0)/a + 1.5d0

C         If there was a large change, check for safety
          if (nlast == nctp .and. iabs(ntry-nctp) > 10) then
            if (fofr > 0) then
              do  ntry = nctp-1, nctp0+1, -1
                fofr =  veff(ntry)-e
                if (fofr < 0) exit
                nctp = ntry
                n2 = ntry
              enddo
            else
              do  ntry = nctp+1, nr-1
                fofr =  veff(ntry)-e
                if (fofr > 0) exit
                nctp = ntry
                n1 = ntry
              enddo
            endif
          endif

C         Exit point
          if (nlast == nctp) then
            if (nctp == nctp0+1) nctp = 2
            return
          endif
          if (fofr > 0.d0)  then
            n2 = nctp
          else
            n1 = nctp
          endif
          nlast = nctp
          nctp = min0(ntry,nr-1)
        enddo
        if (nctp == nctp0+1) nctp = 2
      endif
      end

      subroutine newrhorel(z,lrel2,lgdd,nl,nlr,lmax,a,b,nr,rofi,v,rho,rhoc,kcor,lcor,qcor,
     .  pnu,qnu,qnur,pz,scalsoc,sumec,sumtc,sumev,ec,ev,sfac,ninout,tol,nsp,lfrz,ipr)
C- Makes spherical charge density for a spherical potential, either fully or scalar relativistic
C  ---------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Ci   lrel2 :0 for nonrelativistic, 1 for scalar relativistic, 2 relativistic
Ci         :Add 10 to calculate core relativistically
Ci   lgdd  :0 q2 is coefficient to phidot**2 - p phi**2
Ci         :1 q2 is coefficient to phidot**2 + phi*phidotdot
Ci         :Both produce the same integrated density; see Remarks.
Ci         :lgdd=0 follows Stuttgart conventions.
Ci         :Add 4 to generate phidot,dotdot from outward radial integration only
Ci   nl    :(global maximum l) + 1
Ci   nlr   :second dimension of rho:
Cr         :1 if spherical part of rho is to be generated
Ci         :nl if generated rho is to be decomposed by l
Ci         :In the latter case, the core is not calculated
Ci   lmax  :maximum l for a given site
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Ci   qnu   :energy moments of charge (see Remarks)
Ci   qnur  :relativistic energy moments of charge (1:3,l,mu,2*2)
Ci         :qnur(4,l,mu,2*2) is linearization energy, subject to b.c. ves(rmt)=0
Ci         :if qnur(4) is NULLI, it is assigned to enu generated from the scalar relativistic case.
Ci   pz    :boundary conditions for semicore levels, if they exist
Ci         :State treated as semicore if int(pz)=int(pnu)-1
Ci  scalsoc:not used now; intended to scale the speed of light
Ci   tol   :precision to which wave functions are integrated (scalar Dirac only)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lfrz  :0 for soft core
Ci         :1 for frozen core (do not make core rho)
Ci         :2 for soft core, spin-averaged potential
Ci   ipr   :0 no printout
Ci         :1 summary printout of core
Ci         :2 detailed printout of core
Co Outputs:
Co   ninout:mesh points connecting (in,out) integrations
Co   rho   :spherical charge density times 4*pi*r*r
Co   rhoc  :core density times 4*pi*r*r
Co   sumec :sum of core eigenvalues
Co   sumtc :core kinetic energy = sumec - v*rhoc
Co   sumev :sum of valence eigenvalues
Co   ec    :core eigenvalues
Co   ev    :valence eigenvalues
Co   sfac  :spectroscopic factors (lrel=2 only)
Co         :sfac(1) = charge
Co         :sfac(2) = Spin moment
Co         :sfac(3) = Orbital moment
Co         :sfac(4) = Total moment
Cr Remarks:
Cr  See remarks in newrho.f
Cu Updates
Cu   03 Jun 14 Added semilocal orbitals
C  ---------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nlr,lmax,nr,nsp,lrel2,kcor,lcor,lgdd,ipr,lfrz,ninout(*)
      double precision z,a,b,tol,sumev(nsp),sumec(nsp),sumtc(nsp),ec(*),
     .  ev(*),qcor(2),v(nr,nsp),rofi(nr,2),rho(nr,nlr,nsp),rhoc(nr,nsp),
     .  qnu(3,nl,nsp),pnu(nl,nsp),pz(nl),scalsoc,sfac(4)
      double precision qnur(4,0:(nl-1),2*nl,2,2)
C ... Local parameters
C     logical free
      integer, parameter :: n0=10, NULLI=-99999
      integer lr,ival,lrel,isp,l,ipnu,npnu,imu,la1,la2,kcc,ir,kap1,kap2,kcl(0:n0,2)
      double precision mu,rho1(nr,nlr),mq(nr,nlr),evall(0:n0,2),r,enu
      double precision q0r(2,2),q1r(2,2),q2r(2,2)
      double precision gsmt(2,2),Rgg(2,2),Rff(2,2),mqg1,mqg2,mqf1,mqf2
      double precision GFr(4,4,2,2,0:lmax,2*(lmax+1),4)
      double precision psi(nr,2,2,2),psd(nr,2,2,2),smalpa(2,2),gfn(2,2,2,2)
C     For magnetic moment
      double precision Qqlm(0:lmax,(2*lmax+2))
!     double precision mqlm(0:lmax,(2*lmax+2))
      double precision Jzq, Sq, Lq, Qq, Mtotq, Msq, Mlq, Mscalar, Mrel, Srel,Qscalar
      double precision Rggk(2,2,2,2),Rffk(2,2,2,2)
C     double precision fac,Rggx(2,2,2,2),Rffx(2,2,2,2)

C --- Call scalar relativistic density maker
      lrel = mod(lrel2,10)
C     If lrel<2, call newrho and exit
C     If lrel==2, omit generation of rho but generate nintout, ev, sumev, rhoc from scalar case.
      lr = nlr
C#ifdef ALENA
      if (lrel == 2) lr = -1
C#elseC
C#endif
      call newrho(z,lrel2,lgdd,nl,lr,lmax,a,b,nr,rofi,v,rho,rhoc,kcor,lcor,qcor,
     .  pnu,qnu,pz,scalsoc,sumec,sumtc,sumev,ec,ev,ninout,tol,nsp,lfrz,ipr)
C#ifdef ALENA
C#elseC
C      return
C#endif
      if (lrel /= 2) return

C ... Extract l-resolved linearization energies from scalar relativistic case
      ival = 0
      do  isp = 1, nsp
        do  l = 0, lmax
          npnu = 1
          if (pz(l+1) > 0 .and. int(pnu(l+1,1)-1) == mod(int(pz(l+1)),10)) npnu = 2
          if (npnu /= 1) call rx('newrhorel not ready for semicore valence')
          do  ipnu = 1, npnu    ! ipnu=1 (normal case) or 1..2 (semilocal case)
            ival = ival+1
            if (ipnu < npnu) cycle ! a semilocal orbital
            evall(l,isp) = ev(ival)
            kcl(l,isp) = ninout(ival)
          enddo
        enddo
      enddo
C ... Extract missing l,mu -resolved linearization energies and (in,out) matching from scalar relativistic case
      do  l = 0, lmax
        do  imu = 1, 2*l+2
          if (qnur(4,l,imu,1,1) == NULLI) qnur(4,l,imu,:,:) = (evall(l,1) + evall(l,2))/2
        enddo
      enddo
      if (z == 0) return

C     call setcc(lrel,1d0) !In rseq (sets speed of light in a common block)

C     debugging
      Jzq = 0; Sq=0; Lq=0; Qq=0; Mtotq=0; Qqlm=0; Mscalar = 0; Mrel = 0; Srel = 0; Qscalar = 0

C     print *, '!! omit core'; rho=0
      call dcopy(nr*nsp,rhoc,1,rho,1)

C --- Checking scalar moment
C      do i = 1,nl
C       Mscalar = Mscalar + qnu(1,i,1)-qnu(1,i,2)
C       Qscalar = Qscalar + qnu(1,i,1)+qnu(1,i,2)
C      enddo
C      print*,'Scalar moment M = ',Mscalar
C      print*,'Scalar moment M = ',Mscalar
C      print *, '!!'; Rggx = 0 ; Rffx = 0

C --- Loop over valence states ---
      lr = 1
      do  l = 0, lmax
        if (nlr == nl) lr = l+1
        do  imu = 1, 2*l+2
          mu  = imu - l - 1.5d0

C     ... Accumulate magnetic moment from qnur
          if (imu == 1 .or. imu == 2*l+2) then
            Mrel = Mrel + mu*qnur(1,l,imu,2,2)
            Srel = Srel + dsign(1d0,mu)*qnur(1,l,imu,2,2)
          else
            do la1 = 1,2
              do la2 = 1,2
                Mrel = Mrel + mu*qnur(1,l,imu,la1,la2)
                Srel = Srel + dsign(1d0,mu)*qnur(1,l,imu,la1,la2)
              enddo
            enddo
          endif

C          nn = int(pnu(l+1,1)) - l - 1 !!!!!!!!!!
C          ival = ival+1
C          eval = ev(ival)
C          val(1) = rmax
C          dl = dtan(pi*(0.5d0 - pnu(l+1,1)))
C          slo(1) = dl + 1
C          if (free) val(1) = 1d-30
C          if (free) slo(1) = -val(1)
C          call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v(1,1),g,sum,a,b,rofi,nr,nre,kcc)

          kcc = (kcl(l,1) + kcl(l,2))/2; if (nr-kcc <= 4) kcc = nr-8 ! (In, out) matching point
          enu = qnur(4,l,imu,1,1)
          call rdeqmu(enu,gfn,0,z,v,rofi,nr,nsp,kcc,a,b,l,nl-1,imu,psi,psd,gsmt,GFr)

C     ... Make small parameter p
          do  la1 = 1, 2
            do  la2 = 1, 2
              call prodrw(nr,4,psd(1,1,1,la1),psd(1,1,1,la2),rofi(1,2),smalpa(la1,la2))
            enddo
          enddo
          do la1 = 1,2
            do la2 = 1,2
              q0r(la1,la2) = qnur(1,l,imu,la1,la2)
              q1r(la1,la2) = qnur(2,l,imu,la1,la2)
              q2r(la1,la2) = qnur(3,l,imu,la1,la2)
C             q0r(la1,la2) = q0r(la1,la2) - smalpa(la1,la2)*qnur(3,l,imu,la1,la2)
            enddo
            q0r(la1,la1) = q0r(la1,la1) - smalpa(la1,la1)*qnur(3,l,imu,la1,la1)
          enddo

          do  ir = 2, nr
            r = rofi(ir,1)
            rho1(ir,lr) = 0; mq(ir,lr) = 0

            if (imu == 1 .or. imu == 2*l+2) then ! for a max imu (1 solution)
              Rgg(2,2) = q0r(2,2)*psi(ir,1,2,2)**2 + 2*q1r(2,2)*psi(ir,1,2,2)*psd(ir,1,2,2) + q2r(2,2)*psd(ir,1,2,2)**2
              Rff(2,2) = q0r(2,2)*psi(ir,2,2,2)**2 + 2*q1r(2,2)*psi(ir,2,2,2)*psd(ir,2,2,2) + q2r(2,2)*psd(ir,2,2,2)**2

              rho1(ir,lr) = Rgg(2,2) + Rff(2,2)
              Qqlm(l,imu) = Qqlm(l,imu) + (Rgg(2,2) + Rff(2,2))*rofi(ir,2)

              mq(ir,lr) = - 2d0*mu/(-2*l-1d0)*Rgg(2,2) + 2d0*mu/(2*l+3d0)*Rff(2,2)
            else                ! All other imu
              do kap1 = 1, 2
                do kap2 = 1, 2
                  do la1 = 1, 2
                    do la2 = 1, 2
                      Rggk(kap1,kap2,la1,la2) =
     .                    q0r(la1,la2)*psi(ir,1,kap1,la1)*psi(ir,1,kap2,la2) +
     .                    q1r(la1,la2)*(psi(ir,1,kap1,la1)*psd(ir,1,kap2,la2)+psd(ir,1,kap1,la1)*psi(ir,1,kap2,la2)) +
     .                    q2r(la1,la2)*psd(ir,1,kap1,la1)*psd(ir,1,kap2,la2)
                      Rffk(kap1,kap2,la1,la2) =
     .                    q0r(la1,la2)*psi(ir,2,kap1,la1)*psi(ir,2,kap2,la2)+
     .                    q1r(la1,la2)*(psi(ir,2,kap1,la1)*psd(ir,2,kap2,la2)+psd(ir,2,kap1,la1)*psi(ir,2,kap2,la2)) +
     .                    q2r(la1,la2)*psd(ir,2,kap1,la1)*psd(ir,2,kap2,la2)
                    enddo ! la2
                  enddo !la1
                  Rgg(kap1,kap2) =  sum(Rggk(kap1,kap2,:,:))
                  Rff(kap1,kap2) =  sum(Rffk(kap1,kap2,:,:))
                enddo ! kap2
                rho1(ir,lr) =  rho1(ir,lr) + Rgg(kap1,kap1) + Rff(kap1,kap1)
                Qqlm(l,imu) = Qqlm(l,imu) + (Rgg(kap1,kap1) + Rff(kap1,kap1))*rofi(ir,2)
              enddo ! kap1
              mqg1 = -2d0*mu/(2*l+1d0)*Rgg(1,1) - 2d0*mu/(-2*l-1d0)*Rgg(2,2)
              mqg2 = -dsqrt((2*l+1d0)**2-4*mu**2)/(2*l+1d0)*(Rgg(1,2)+Rgg(2,1))
              mqf1 = -2d0*mu/(-2*l+1d0)*Rff(1,1)-2d0*mu/(2*l+3d0)*Rff(2,2)
              mqf2 = 0          !-dsqrt((2*l+1d0)**2-4*mu**2)/(2*l+1d0)*(Rff(1,2)+Rff(2,1))
              mq(ir,lr) = mq(ir,lr) + mqg1 + mqg2 - mqf1 - mqf2
            endif ! imu > 1 and imu < 2l+2
            rho(ir,lr,1) = rho(ir,lr,1) + (rho1(ir,lr)-mq(ir,lr))/2 !/(2*l+2d0)   ! /(2l+1) for the test purpose only
            rho(ir,lr,2) = rho(ir,lr,2) + (rho1(ir,lr)+mq(ir,lr))/2 !/(2*l+2d0)
            Sq = Sq + mq(ir,lr)*rofi(ir,2)/2 ! Total spin
          enddo                 !  loop over ir

          Jzq = Jzq + mu*Qqlm(l,imu)
          Qq = Qq + Qqlm(l,imu)

        enddo                   !  loop over imu
      enddo                     !  loop over l

C      print *, sngl(rggx(1,1,:,:) + rffx(1,1,:,:))
C      print *, sngl(rggx(2,2,:,:) + rffx(2,2,:,:))
C      print *, sngl(rggx(1,1,:,:)+rggx(2,2,:,:) + rffx(1,1,:,:)+rffx(2,2,:,:))
C      print *, (sum(rggx(1,1,:,:)+rggx(2,2,:,:) + rffx(1,1,:,:)+rffx(2,2,:,:)))
C      imu = 3
C      print *, sum(Qqlm(2,:)),Qqlm(2,imu)
C     call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)

C ... Total magnetic dipole moment; spin mag moment; orbital MM in Bohr magneton
      Mtotq = Jzq + Sq
      Msq = 2*Sq
      Lq = Jzq - Sq
      Mlq = Lq
      sfac(1) = Qq
      sfac(2) = Msq
      sfac(3) = Lq
      sfac(4) = Mtotq

C     call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)

C     call rx0('newrhorel 1948')
      end

      subroutine newrho(z,lrel2,lgdd,nl,nlr,lmax,a,b,nr,rofi,v,rho,rhoc,kcor,lcor,qcor,
     .  pnu,qnu,pz,scalsoc,sumec,sumtc,sumev,ec,ev,ninout,tol,nsp,lfrz,ipr)
C- Makes charge density for a spherical potential and related quantites
C  ---------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Ci   lrel2 :0 for nonrelativistic, 1 for relativistic
Ci         :10s digit nonzero => calculate core with Dirac equation
Ci   lgdd  :0 q2 is coefficient to phidot**2 - p phi**2
Ci            This is standard, and follows Stuttgart conventions
Ci         :1 q2 is coefficient to phidot**2 + phi*phidotdot
Ci         :  This is approximate but decouples potential parameters from charges
Ci         :Both 0,1 produce the same integrated density; see Remarks.
Ci         :Add 4 to generate phidot,dotdot from outward radial integration only
Ci   nl    :dimensions pnu,qnu
Ci   nlr   :second dimension of rho:
Cr         : 1 generate rho(r) from spherical V, including core part
Ci         :nl Same as 1, except:
Ci         :   rho is to be decomposed by l
Ci         :   Core rho is not calculated or added to rho
Cr         :0  suppress generation of valence and core rho
Cr         :-1 suppress generation of valence rho.  rhoc is calculated if also lfrz=0
Ci         :
Ci   lmax  :maximum l for this sphere
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Ci   qnu   :energy moments of charge (see Remarks).
Ci         :Not used if nlr = 0 or nlr = -1
Ci   pz    :boundary conditions for semicore levels, if they exist
Ci         :State treated as semicore if int(pz)=int(pnu)-1
Ci   tol   :precision to which wave functions are integrated.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   g     :normalized wave function times r
Ci   gp    :energy derivatives of g
Ci   lfrz  :0 for soft core
Ci         :1 for frozen core (do not make core rho)
Ci         :2 for soft core, spin-averaged potential
Ci   ipr   :0 no printout
Ci         :1 summary printout of core
Ci         :2 detailed printout of core
Co Outputs:
Co   rho   :spherical charge density times 4*pi*r*r
Co   rhoc  :core density times 4*pi*r*r
Co   sumec :sum of core eigenvalues
Co   sumtc :core kinetic energy = sumec - v*rhoc
Co   sumev :sum of valence eigenvalues
Co   ec    :core eigenvalues
Co   ev    :valence eigenvalues
Cr Remarks:
Cr   rho is determined by boundary conditions pnu and moments qnu.
Cr   For rmax>10 phidot, phidotdot are set to zero.
Cr
Cr   Recall that val,slo correspond to u = r*phi, so
Cr     rmax * slo / val = D_u = D_phi + 1.
Cr     For val=rmax  slo = D + 1
Cr
Cr   Switch lgdd (1s digit) concerns the contribution of <phidot|phidot> to
Cr   the sphere charge.  The zeroth moment may be defined as the amount of
Cr   charge in channel l, i.e.
Cr     q^0 = (amount of <phi|phi>) +  p* (amount of <phidot|phidot>)
Cr   where
Cr     p = <phidot|phidot>
Cr   Then q^0 is not the amount of phi*phi inside the sphere, but rather
Cr   but instead, it is
Cr     (q_0 - p q_2)
Cr   since q2 is the amount of <phidot|phidot>.  The charge density is
Cr     rho= (q_0 - p q_2) phi*phi + 2 q_1 phi*phidot + q_2 phidot*phidot
Cr   This is the Stuttgart convention (lgdd=0)
Cr
Cr   Methfessel convention:  to avoid explicit dependence of rho on
Cr   p, he approximated p*phi*phi with -phi*phidotdot (they have the
Cr   same integrated charge).  Then
Cr     rho= q_0 phi*phi + 2 q_1 phi*phidot +
Cr          q_2 (phidot*phidot + phi*phidotdot)
Cr  Semilocal orbitals:
Cr    This routine treats a semilocal orbital as a valence orbital
Cr    with q0=(4*l+2)/nsp, q1=q2=0.
Cr    This is an approximation since:
Cr    * Moments q betw/ semilocal and valence orbitals are assumed 0
Cr    * In this version the semilocal orbital is normalized to unit
Cr      charge: no charge is allowed to spill out of the sphere.
Cr      The approximation could be improved by attaching a tail
Cr      and allowing the state to spill outside of the sphere.
Cr    * Neither of these approximations matter for free atoms,
Cr      which is all this option is used for so far.
Cu Updates
Cu   26 Aug 17 Free atom case and Dirac core levels sought (lrel2>10),
Cu             then Dirac valence levels also generated
Cu   23 Mar 17 New option (s_ctrl%lasa,4096) to compute phi,phidot
Cu             from outward radial integration only
Cu   03 Jun 14 Added semilocal orbitals
C  ---------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nlr,lmax,nr,nsp,lrel2,kcor,lcor,ipr,lfrz,lgdd,ninout(*)
      double precision z,a,b,tol,sumev(nsp),sumec(nsp),sumtc(nsp),ec(*),
     .  ev(*),qcor(2),v(nr,nsp),rofi(nr),rho(nr,nlr,nsp),rhoc(nr,nsp),
     .  qnu(3,nl,nsp),pnu(nl,nsp),pz(nl),scalsoc
C ... Dynamically allocated arrays
      real(8),pointer :: gloc(:)
C ... Local parameters
C     logical :: debug=.false.
      logical free
      integer, parameter :: n0=10
      integer konfig(0:n0),l,isp,ir,ival,lrel,lphidot,nn,nre,jr,lmaxc,lr,k,ipnu,npnu,kc
      double precision c,dl,dphi,dphip,eb1,eb2,eval,fllp1,gfac,p,
     .  phi,phip,pi,pl,q0,q00,q1,q2,r,rmax,ro,rocrit,slo(5),sum,tmc,val(5)
      real(8),target :: g(2*nr),gp(2*nr*4)
      integer nreD(0:lmax)
      double precision E,Ekmu(2*(2*lmax+1),0:lmax),tcorD(0:lmax),xx,evl(0:lmax,2)
      character*1 ang(6)
C     procedure(integer) :: nglob
      procedure(real(8)) :: dsum
C ... Heap
      data ang /'s','p','d','f','g','h'/

C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      lrel = mod(lrel2,10)
      call setcc(lrel,1d0)
      lr = 1
      lphidot = 1; if (mod(lgdd,10)/4 /= 0) lphidot = 4
      rocrit = 0.002d0/4
      pi = 4d0*datan(1d0)
      eb1 = -50d0
      eb2 =  50d0
      rmax = rofi(nr)
      free = (rmax > 9.99d0)
      call dpzero(evl,size(evl))
      call config(pnu,lmax,z,konfig,lmaxc)
C     Partial core occupation
      if (kcor > 0) then
        lmaxc = max(lmaxc,lcor)
        konfig(lcor) = max(konfig(lcor),kcor+1)
      endif
C     Semicore state as valence: remove from core konfig
      do  l = 0, lmax
        if (int(pnu(l+1,1)-1) == mod(int(pz(l+1)),10)) then
          if (kcor > 0 .and. l == lcor)
     .      call rx('newrho not implemented for pz and lcor')
          konfig(l) = mod(pz(l+1),10d0)
        endif
      enddo

C --- Calculate core density ---
      if (nlr == 1 .or. nlr == -1) then
        if (lfrz /= 1) then
          call dpzero(rhoc,nr*nsp)
          k = 0; if (lfrz == 2) k = 10
          if (lrel2 >= 10) k = k+100*mod(lrel2/10,10)
          call rhocor(k,z,lmaxc,nsp,konfig,a,b,nr,rofi,v,g,
     .                kcor,lcor,qcor,tol,ec,sumec,sumtc,rhoc,0d0,ipr)
        endif
        if (nlr == 1) call dcopy(nr*nsp,rhoc,1,rho,1)
      endif
C     call prrmsh('starting rhoc for atom ',rofi,rho,nr,nr,nsp)

      call setcc(lrel,scalsoc)
C --- Loop over valence states ---
      ival = 0
      do  isp = 1, nsp
        sumev(isp) = 0d0
        do  l = 0, lmax
        if (pz(l+1) > 0 .and. int(pnu(l+1,1)-1) == mod(int(pz(l+1)),10)) then
          npnu = 2
        else
          npnu = 1
        endif
        do  ipnu = 1, npnu  ! ipnu=1 (normal case) or 1..2 (semilocal case)

        if (ipnu < npnu) then  ! a semilocal orbital
          lr = 1
          q0 = (4*l+2)/nsp
          q1 = 0
          q2 = 0
          pl = mod(pz(l+1),10d0)
        else
          if (nlr == nl) lr = l+1
          q0 = max(qnu(1,l+1,isp),0d0)
          q1 = qnu(2,l+1,isp)
          q2 = qnu(3,l+1,isp)
          pl = pnu(l+1,isp)
        endif

        if (q0 < 1d-6) cycle ! No contribution to density
        nn = int(pl) - l - 1
        ival = ival+1
        eval = ev(ival)

        val(1) = rmax
        dl = dtan(pi*(0.5d0 - pl))
        slo(1) = dl + 1
        if (free) val(1) = 1d-30
        if (free) slo(1) = -val(1)
        call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,rofi,nr,nre,kc)
        ev(ival) = eval
        evl(l,isp) = eval
        ninout(ival) = kc
        sumev(isp) = sumev(isp) + eval*q0 + q1

        if (lrel2 >= 10 .and. isp == 1 .and. free .and. ipr > 0) then
          E = eval; xx = 0
          call rdeqcore(1,E,z,v,rofi,nr,nsp,nn,kc,a,b,l,
     .      Ekmu(1,l),[xx],tcorD(l),nreD(l),xx)
        endif

        if (q0 < 1d-6) cycle ! No contribution to density
        if (nlr <= 0) cycle   ! Density not calculated
        ro = g(nr)**2
        if (.not.free .and. ro < rocrit) write(*,1) l,nn,nre,ro
    1   format(' NEWRHO (warning): PHP,PHPP set to zero,l,nn,nre,rho=',3i5,2f8.4)
        if (free .or. ro < rocrit .or. ipnu < npnu) then
          call dpzero(gp,8*nr)
          p = 0
        else
          val(1) = val(1)/dsqrt(sum)
          slo(1) = slo(1)/dsqrt(sum)
C          call phidot(z,l,v(1,isp),eval,a,b,rofi,nr,g,val,slo,tol,
C     .                nn,gp,phi,dphi,phip,dphip,p)
          if (lphidot /= 1) then
            allocate(gloc(2*nr))
            call dcopy(2*nr,g,1,gloc,1)
          else
            gloc => g
          endif
          call phidx(lphidot,z,l,v(1,isp),0d0,0d0,rofi,nr,2,tol,eval,val,
     .      slo,nn,gloc,gp,phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
          if (.not. associated(gloc,g)) deallocate(gloc)
        endif
        fllp1 = l*(l+1)
C  ...  Case add q2 phi phidd rho
        if (mod(mod(lgdd,10),4) /= 0) then
          k = 2*nr
          do  ir = 2, nre
            jr = ir + nr
            r = rofi(ir)
            tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
            gfac = 1d0 + fllp1/(tmc*r)**2
            rho(ir,lr,isp) =  rho(ir,lr,isp) + q0*(gfac*g(ir)**2 + g(jr)**2) +
     .                               2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) +
     .                               q2*(gfac*(gp(ir)**2 + g(ir)*gp(ir+k)) +
     .                               gp(jr)**2 + g(jr)*gp(jr+k))

C            if (debug .and. ir >= 1000 .and. l == 1) then
C              print 333, ir,l,isp,rofi(ir),g(ir),gp(ir),gp(ir+k)
C  333         format(3i4,6f20.10)
C            endif

          enddo
C  ...  Case add -p q2 phi phi + q2 phidot^2 into rho
        else
          q00 = q0-p*q2
          do  ir = 2, nre
            jr = ir + nr
            r = rofi(ir)
            tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
            gfac = 1d0 + fllp1/(tmc*r)**2
            rho(ir,lr,isp) = rho(ir,lr,isp) + q00*(gfac*g(ir)**2 + g(jr)**2) +
     .                               2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) +
     .                               q2*(gfac*gp(ir)**2 + gp(jr)**2)
          enddo
        endif
        enddo ! loop over (semicore + valence)
        enddo ! loop over l
      enddo   ! loop over spins

      if (lrel2 >= 10 .and. free .and. ipr > 1) then
        if (nsp == 1) call info0(1,1,0,' Dirac valence levels:%N nl  chg%6feval(S)%6feval(D)')
        if (nsp == 2) call info0(1,1,0,' Dirac valence levels:%N nl  chg%12feval(S)%11feval(D)')
        do  l = 0, lmax
          if (qnu(1,l+1,1)+qnu(1,l+1,2) == 0) cycle
          call info8(1,0,0,' %i'//ang(l+1)//'  %;5F %n;11,6D  %n;11,6D',
     .      int(pnu(l+1,1)),max(qnu(1,l+1,1)+qnu(1,l+1,2),0d0),
     .      nsp,evl(l,:),2*(2*l+1),ekmu(1,l),7,8)
        enddo
      endif

      end

      subroutine potpar(nl,nsp,lmx,z,rmax,avw,ekap,lso,loptc,lmpol,a,nr,rofi,v,pnu,idmod,
     .  exc,qnu,qnur,idu,uh,jh,facso,thrpv,thrpvl,g,gp,pp,pprel,sop,pmpol,rgrad,t,GFr)
C- Generates potential parameters for given potential
C  ---------------------------------------------------
Ci Inputs:
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   ekap  :muffin-tin zero plus kap2
Ci   lso   :if true, make spin-orbit coupling parms sop
Ci   loptc :if true, make matrix elements <phi grad phi> etc
Ci   lmpol :if true, make multipole moments of phi,phidot
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :rofi(:,1) = radial mesh points
Ci         :rofi(:,2) = radial mesh quadrature weights
Ci   v     :spherical potential
Ci   pnu   :Determines b.c. for particular l channel; see remarks
Ci   idmod :0 or 1 : linearization energy enu for l channel determined
Ci         :         by input log derivative spec'd by pnu
Ci         :2      : linearization energy enu for l channel determined
Ci         :         from pp(1); pnu is output; see Remarks
Ci         :See also Remarks
Ci   exc   :XC energy used to calc. thrpv
Ci   qnu   :energy-weighted moments of the sphere charges
Ci         :used only in the calculation of 3pV
Ci   qnur  :relativistic moments. qnur(4) is used for linearization energy
Ci   idu   :idu(l+1)=0 => this l has no nonlocal U matrix
Ci                   1..3 this l has no nonlocal U matrix
Ci                   4    majority and minority C,E shifted by U; see remarks
Ci                   5    majority and minority C,E shifted by U and J; see remarks
Ci         :Note: this routine doesn't implement LDA+U; see remarks
Ci   uh    :U in LDA+U
Ci   jh    :J in LDA+U
Ci   facso
Co Outputs:
Co   pp:   :2nd generation potential parameters (see Remarks)
Co         :c = center of band = enu + omega-
Co         :sqrt(delta) has proper sign (that of phi-).
Co         :pp(1) enu
Co         :pp(2) calpha
Co         :pp(3) srdel = sqrt(delta) with proper sign (that of phi-).
Co         :pp(4) palpha
Co         :pp(5) gamma, or Q in Varenna notes
Co         :pp(6) alpha, or qbar in Varenna notes
Co   pprel :Relativistic potential parameters in kappa-mu representation
Co         :pprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),isp,jsp)  and  mu = imu-1 - (l+1/2)
Co         :pprel(1,l,imu,isp,jsp) cgamma(l,imu,isp,jsp)
Co         :pprel(2,l,imu,isp,jsp) gamma(l,imu,isp,jsp)
Co         :pprel(3,l,imu,isp,jsp) D+, where (D+) D = delta(l,imu,isp,jsp)
Co         :pprel(4,l,imu,isp,jsp) pgamma(l,imu,isp,jsp)
Co         :pprel(5,l,imu,isp,jsp) enu
Co         :If NULL, it will be averaged from scalar rel (pp(1)+pp(2))/2
Co         :Note that enu and C are SHIFTED by ves, while scalar rel enu,C are NOT.
Co   thrpv : 3 PV for sphere
Co   thrpvl: 3 PV for sphere, decomposed by l-channel and spin
Co   g     : wave function; see rseq
Co   gp    : phidot, except for ...
Co   sop   : spin-orbit coupling parameters
Co         : sop(l+1,is1,is2,i) matrix elts between spins is1 and is2
Co         : for quantum number l. Nine types of integrals are stored.
Co         : (i=1) <phi so  phi> (i=2) <phi so  dot> (i=3) <dot  so dot>
Co         : (i=4) <phi ||  phi> (i=5) <phi ||  dot> (i=6) <dot  || dot>
Co         : (i=7) <phi Bxc phi> (i=8) <phi Bxc dot> (i=9) <dot Bxc dot>
Co   pmpol :integral (phi-or-phidot * phi-or-phidot * r**l) :
Co         :matrix elements of w.f. * wf * r**l for multipole moments
Co   rgrad :Matrix elements of radial part of gradient operator
Co   GFr(4,4,2,2,l,imu,4):overlap integrals for the calculation of averages [Turek (6.141)] (see rdeq)
Cr Remarks:
Cr   pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Cr   pp's are generated for a specified potential.
Cr   The linearization energy enu is calculated from the input pnu.
Cr   Alternatively, by setting idmod=2 for a specified enu is frozen and the corresponding pnu calculated.
Cr
Cr   Potential parameters generated from (ekap=0):
Cr     omega(-) = -(phi/phidot) (-l-1-dnu)/(-l-1-dnudot)
Cr     omega(+) = -(phi/phidot) (l-dnu)/(l-dnudot)
Cr     phi(+) = phi + omega(+) phidot
Cr     phi(-) = phi + omega(-) phidot
Cr     C  = e + omega(-)
Cr     Delta = phi(-)**2 sdivw rmax/2
Cr     Gamma = (phi(-)/phi(+)) sdivw / (2 (2l+1))
Cr     where sdivw = (rmax/avw)**(l+l+1)
Cr   These conventions differ by those of MSM by the scale factor
Cr   sdivw for delta and gamma.
Cr   For general ekap (see Kanpur notes)
Cr     C = e - W{K,phi} / W{K,phidot}
Cr       = e - phi/phidot * (D{phi}-D{K})/(D{phidot}-D{K})
Cr     sqrt(delta) = - sqrt(avw/2) / W{K,phidot}
Cr                 = - sqrt(avw/2) /wsr / K_l / phidot /(D{phidot}-D{K})
Cr     gamma = W{J,phidot} / W{K,phidot}
Cr           = J_l/K_l * (D{phidot}-D{J})/(D{phidot}-D{K})
Cr   Uses W{a,b} = r*r*(ab' - a'b) = r*a*b*[D{b}-D{a}]
Cr
Cr   LDA+U parameters:  No LDA+U is implemented in this code.
Cr   For now, parameters are just a device to shift pot pars by a constant.
Cr   You can shift (enu,C) of both spins by U (idu=4)
Cr        or shift (enu,C) of spin1 by U, spin2 by J (idu=5)
Cr
Cr   potpar uses convention ves(rmt) = 0.
Cr   Caller should add shift to pp(1:2,:,:) and pprel(1) and pprel(5) returned by potpar
Cl Local variables
Cl   lpzi  :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz (none implemented for now)
Cr Updates
Cu  10 Apr 12 Repackaged radial mesh quadrature (radwgt)
Cu  21 Dec 05 (wrl) potential shifts to mimic LDA+U
Cu  11 Jul 05 Updated soprm call to render compatible with fp
Cu  18 Jun 04 (A Chantis) relativistic potential parameters
Cu   4 Apr 04 New matrix elements of <phi Bxc phi>
Cu  05 Jan 04 bug fix making SO matrix elements when nl ne lmx+1
Cu  07 Feb 03 When calc. SO matrix elements, also make <phi|phi>
Cu            for inclusion of external magnetic field.
Cu  20 Apr 98 Optical matrix elements adapted from Sergey Rashkeev
Cu            MvS spin-orbit matrix elements adapted from V. Antropov
C  ---------------------------------------------------
      implicit none
C ... Passed parameters
      logical lso,lmpol,loptc
      integer nl,lmx,nr,nsp,idmod(0:nl-1),idu(4),n0
      parameter (n0=10)
      double precision uh(4),jh(4),t(0:n0-1,2)
      double precision z,ekap,rmax,avw,a,thrpv,facso
      double precision rofi(nr,2),v(nr,nsp),g(*),gp(nr,4),pp(6,0:nl-1,nsp),
     .                 pnu(0:nl-1,nsp),qnu(3,0:nl-1,nsp),qnur(4,0:nl-1,2*nl,2,2),exc(2),
     .                 sop(0:nl-1,nsp,nsp,9),pmpol(nl,nl,2*nl-1,3,nsp),
     .                 rgrad(4,2,nl,nsp),pprel(5,0:nl-1,2*nl,2,2),
     .                 thrpvl(0:nl-1,nsp),GFr(4,4,2,2,0:lmx,2*(lmx+1),4)
C ... Local parameters
      logical lsor
      integer i,j,l,nn,nre,stdo,intopt,kc,lphidot
      integer lpzi(0:n0),nmsave(0:lmx)
      double precision eb1,eb2,pi,b,e,val(5),slo(5),sum,sdivw,xx,
     .  phi,dphi,phip,dphip,dlphi,dlphip,phplus,phmins,
     .  el(2),ql(2),e1me2,e1pe2,e1e2,d,dlmins,dlplus,
     .  fi(0:9),gi(0:9),ptmp1(2,2),ptmp2(2,2),ptmp3(2,2),ptmp4(2,2),
     .  ptmp5(2,2),ptmp6(2,2),ptmp7(2,2),ptmp8(2,2)
C     double precision omegam,omegap
      double precision wk,wkdot,wjdot,scl,enu(0:8,2),fac1,fac2
      real(8), allocatable :: wgt(:),dv(:),wkv(:),phiv(:,:,:,:)
      real(8),parameter::    NULLR =-99999d0, tol=1d-12
      procedure(integer) :: nglob
C ... For Dirac equation
      integer mumax,imu,lrel,la1,la2
      double precision scalede,phisq(0:lmx,nsp),gmt(2,2),fmt(2,2),gmtde(2,2),fmtde(2,2),gsmt(2,2)
      double precision psi(nr,2,2,2),psd(nr,2,2,2),smalpa(2,2),gfn(2,2,2,2)
      double precision c
      common /cc/ c

      stdo = nglob('stdo')
      lrel = mod(nglob('lrel'),10)
      lsor = lso .or. lrel == 2  ! for now make L.S parameters in Dirac case
      intopt = 10*nglob('lrquad')
      do  l = 0, min(lmx,3)
C       LDA+U-like shifts in C
        if (idu(l+1) <= 3) then
        elseif (idu(l+1) <= 4) then
          call info2(20,0,0,'%10fpotpar l=%i: enu and C shifted by %d',l,uh(l+1))
C         if (idu(l+1) <= 3) call rx('POTPAR: LDA+U not implemented in ASA')
        elseif (idu(l+1) == 5) then
          call info5(20,0,0,'%10fpotpar l=%i: enu+ and C+ shifted by'//
     .      ' %d;  enu- and C- shifted by %d',l,uh(l+1),jh(l+1),0,0)
        endif
      enddo

      if (lsor .or. lmpol .or. loptc) then
        allocate(phiv(nr,nl,nsp,2)) ! phiv(:,:,:,2) = phidot
        allocate(wkv(nr*4))
        allocate(dv(nr))
        allocate(wgt(nr))
        call radwgt(intopt,rmax,a,nr,wgt)
        call iinit(lpzi,n0+1)
      endif

      thrpv = 0
      eb1 = -20d0
      eb2 =  20d0
      pi = 4d0*datan(1d0)
      b = rmax/(dexp(a*nr - a) - 1d0)
      call bessl2(ekap*rmax**2,0,lmx+1,fi,gi)

      do  i = 1, nsp
        do  l = 0, lmx
          if (mod(idmod(l),10) == 2) then
            e = pp(1,l,i)
            call rsq1(0,e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
            xx = 0.5 - datan(slo(1)*rmax/val(1)-1)/pi
            if (xx <= 0) call rx('potpar: something wrong with pnu')
            pnu(l,i) = int(pnu(l,i)) + xx
          else
            e = -0.5d0
            nn = int(pnu(l,i))-l-1
            val(1) = rmax
            slo(1) = 1 + dtan(pi*(0.5d0 - pnu(l,i)))
          endif
          call rseq(eb1,eb2,e,tol,z,l,nn,val,slo,v(1,i),g,sum,a,b,rofi,nr,nre,kc)
          nmsave(l) = kc
          val(1) = val(1)/dsqrt(sum)
          slo(1) = slo(1)/dsqrt(sum)
          lphidot = 1; if (mod(idmod(l)/100,10) /= 0) lphidot = 4
          call phidx(lphidot,z,l,v(1,i),0d0,0d0,rofi,nr,2,tol,e,val,slo,nn,
     .      g,gp,phi,dphi,phip,dphip,pp(4,l,i),0d0,0d0,0d0,0d0)
          t(l,i) = e - (v(nr,i)-2*z/rofi(nr,1))
          phisq(l,i) = phi*phi
C     ... Keep local copies of phi and phidot for SO coupling
          if (lsor .or. lmpol .or. loptc) then
            call dpscop(g,phiv(1,l+1,i,1),nr,1,1,1d0)
            call dpscop(gp,phiv(1,l+1,i,2),nr,1,1,1d0)
            enu(l,i) = e
          endif
          dlphi =  rmax*dphi/phi
          dlphip = rmax*dphip/phip
          sdivw = (rmax/avw)**(l+l+1)
C     ... Following valid for ekap=0 only:
C         omegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
C         omegap = -(phi/phip)*(l-dlphi)/(l-dlphip)
C         phplus = phi + omegap*phip
C         phmins = phi + omegam*phip
C         pp(1,l,i) = e
C         pp(2,l,i) = e + omegam
C         pp(3,l,i) = phmins*dsqrt(sdivw*rmax/2)
C         pp(5,l,i) = phmins/(2*(2*l+1)*phplus)*sdivw
C         pp(6,l,i) = pp(5,l,i)
C     ... The following for general ekap:
          wk     = (dlphi-l)*gi(l)  + (l+l+1)*gi(l+1)
          wkdot  = (dlphip-l)*gi(l) + (l+l+1)*gi(l+1)
          wjdot  = (dlphip-l)*fi(l) + ekap*rmax**2*fi(l+1)/(l+l+1)
C         omegam = -(phi/phip)*wk/wkdot
C         scl    = Wronskian(phi,dot) ... should be 1.  Scale srdel to make compatible w/ ekap=0
          scl    = phi*phip*(dlphi-dlphip)*rmax
C     ... Force enu to stay really fixed for idmod=2
          if (mod(idmod(l),10) /= 2) pp(1,l,i) = e
          pp(2,l,i) = e - (phi/phip)*wk/wkdot
          pp(3,l,i) = -dsqrt(sdivw/2/rmax)/phip/wkdot*scl
          pp(5,l,i) = sdivw*wjdot/wkdot
          pp(6,l,i) = pp(5,l,i)
C     ... Constant shift for potential parameters, mimicking LDA+U
          if (idu(l+1) == 4) then
            pp(1,l,i) = pp(1,l,i) + uh(l+1)
            pp(2,l,i) = pp(2,l,i) + uh(l+1)
          endif
C     ... Spin-dependent constant shift for enu,C, mimicking LDA+U
          if (idu(l+1) == 5) then
            if (i == 1) then
              pp(1,l,i) = pp(1,l,i) + uh(l+1)
              pp(2,l,i) = pp(2,l,i) + uh(l+1)
            elseif (i == 2) then
              pp(1,l,i) = pp(1,l,i) + jh(l+1)
              pp(2,l,i) = pp(2,l,i) + jh(l+1)
            endif
          endif

C --- Calculate 3PV ala MSM J. Phys. F 12, 141, Eqn 2.24 ---
C         NB: formula valid only for ekap=0
          thrpvl(l,i) = 0
          if (qnu(1,l,i) <= 0 .or. qnu(3,l,i) <= 0 .or. exc(i) == 0) cycle
          d     = qnu(2,l,i)**2 - qnu(1,l,i)*qnu(3,l,i)
          e1pe2 = qnu(2,l,i)*qnu(3,l,i)/d
          e1e2  = qnu(3,l,i)**2/d
C The following two apply if have also third moment ...
C         e1pe2 = (qnu(2,l,i)*qnu(3,l,i) - qnu(1,l,i)*qnu(4,l,i))/d
C         e1e2  = (qnu(3,l,i)**2 - qnu(2,l,i)*qnu(4,l,i))/d
          e1me2 = dsqrt(max(e1pe2**2 - 4*e1e2,0d0))
          if (e1me2 /= 0) then
            el(1) = (e1pe2 + e1me2)/2
            el(2) = (e1pe2 - e1me2)/2
            ql(1) = ( qnu(2,l,i) - el(2)*qnu(1,l,i))/e1me2
            ql(2) = (-qnu(2,l,i) + el(1)*qnu(1,l,i))/e1me2
          else
            el(1) = e1pe2
            el(2) = e1pe2
            ql(1) = 0
            ql(2) = 0
          endif

C         First order approximation to phi, dphi
          phplus = phi + el(1)*phip
          phmins = phi + el(2)*phip
          dlplus = dlphi - el(1)/phi**2/rmax
          dlmins = dlphi - el(2)/phi**2/rmax

C#ifndef FIRST_ORDER_PRESSURE
C         Get phi, dphi by reintegrating wave function
          call rsq1(0,el(1)+e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          dlplus = rmax*(slo(1)-val(1)/rmax)/val(1)
          call gintsr(g,g,a,nr,z,el(1)+e,l,v(1,i),rofi,sum)
          phplus = val(1)/rmax/dsqrt(sum)
          call rsq1(0,el(2)+e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          dlmins = rmax*(slo(1)-val(1)/rmax)/val(1)
          call gintsr(g,g,a,nr,z,el(2)+e,l,v(1,i),rofi,sum)
          phmins = val(1)/rmax/dsqrt(sum)
C#endif  FIRST_ORDER_PRESSURE

          thrpvl(l,i) =
     .      rmax*ql(1)*phplus**2*
     .      (dlplus*(dlplus+1) - l*(l+1) + rmax**2*(el(1)+e-exc(i))) +
     .      rmax*ql(2)*phmins**2*
     .      (dlmins*(dlmins+1) - l*(l+1) + rmax**2*(el(2)+e-exc(i)))
          thrpv = thrpv + thrpvl(l,i)

        enddo
      enddo

      call info2(40,1,0,' K.E., phi(Rmt) by l: %n;12,6D',lmx+1,t(0,1))
      if (nsp == 2) call info2(40,0,0,'%14pspin 2: %n;12,6D',lmx+1,t(0,2))

C --- Make spin-orbit matrix elements of phi, phidot ---
      if (lsor) then
C       call pshpr(51)

C   ... Gradient of average v
        call dpcopy(v,wkv,1,nr,1d0/nsp)
        if (nsp == 2) call dpadd(wkv,v(1,2),1,nr,.5d0)
        call radgra(a,b,nr,rofi,wkv,dv)

C   ... Calculate matrix elements for L.S spin orbit coupling
        call soprm(5,lpzi,phiv,phiv(1,1,1,2),xx,nr,nsp,nl-1,lmx,v,dv,
     .    enu,xx,z,rofi,wgt,facso,sop,xx,' ')

C   ... Matrix elements for constant magnetic field
        call dvset(dv,1,nr,1d0)
        call soprm(2,lpzi,phiv,phiv(1,1,1,2),xx,nr,nsp,nl-1,lmx,v,dv,
     .    enu,xx,z,rofi,wgt,facso,sop(0,1,1,4),xx,' ')
C       Correct <dot||dot> term from (already known) diagonal part
C       <phi|phi> should be unity; <phi|dot> should be zero
        do  l = 0, lmx
          fac1 = pp(4,l,1)/sop(l,1,1,6)
          fac2 = pp(4,l,2)/sop(l,2,2,6)
          sop(l,1,1,6) = sop(l,1,1,6)*sqrt(fac1*fac1)
          sop(l,1,2,6) = sop(l,1,2,6)*sqrt(fac1*fac2)
          sop(l,2,1,6) = sop(l,2,1,6)*sqrt(fac2*fac1)
          sop(l,2,2,6) = sop(l,2,2,6)*sqrt(fac2*fac2)
        enddo

C   ... Matrix elements of XC field
        call dpcopy(v,dv,1,nr,0.5d0)
        call dpadd(dv,v(1,nsp),1,nr,-0.5d0)
C       call prrmsh('vxc ',rofi,dv,nr,nr,1)
        call soprm(2,lpzi,phiv,phiv(1,1,1,2),xx,nr,nsp,nl-1,lmx,v,dv,
     .    enu,xx,z,rofi,wgt,facso,sop(0,1,1,7),xx,' ')

C       debugging: generate M<B>
C        print *, '!!'
C        if (nsp == 2) then
C          wk = 0
C          do  l = 0, lmx
C            ql(1) = qnu(1,l,1) - pp(4,l,1)*qnu(3,l,1)
C            ql(2) = qnu(1,l,2) - pp(4,l,2)*qnu(3,l,2)
C            fac1 =
C     .        ql(1)*sop(l,1,1,7) +
C     .        qnu(2,l,1)*sop(l,1,1,8) +
C     .        qnu(3,l,1)*sop(l,1,1,9)
C            fac2 =
C     .        ql(2)*sop(l,2,2,7) +
C     .        qnu(2,l,2)*sop(l,2,2,8) +
C     .        qnu(3,l,2)*sop(l,2,2,9)
C            wk = wk + fac1-fac2
C            write(stdo,333) l, fac1,fac2,fac1-fac2
C  333       format(' l=',i1,'   q+<B>=',f12.6,'   q-<B>=',f12.6,
C     .             '   M<B>=',f12.6)
C          enddo
C          write(stdo,'('' <Bxc*M>'',f12.6)') wk
C        endif
      endif

C --- Fully relativistic potential parameters ---
      if (lrel == 2) then

        call dpzero(GFr,size(GFr))
        do  l = 0, lmx
          mumax= 2*(l+1)
          do  imu = 1, mumax    ! mu = dble(imu-l) - 1.5d0
C#ifdef ALENA | ALENAP
            if (mod(idmod(l),10) == 2) then
              if (pprel(5,l,imu,1,1) == NULLR) then
                pprel(5,l,imu,:,:) = (pp(1,l,1) + pp(1,l,2))/2
              endif
              e = pprel(5,l,imu,1,1)
              qnur(4,l,imu,1,1) = e
            else
              e = qnur(4,l,imu,1,1)
              if (e == NULLR) call rxi('sorry, no enu available for l=',l)
              pprel(5,l,imu,1,2) = 0; pprel(5,l,imu,2,1) = 0
              pprel(5,l,imu,1,1) = e; pprel(5,l,imu,2,2) = e
            endif
C           Make g (scal rel) as starting point to make grel
            call rsq1(0,e,l,z,v(1,1),nr,g,val,slo,nn,a,b,rofi,nr)
            eb1 = e - .3d0; eb2 = e + .3d0; xx = e
            call rseq(eb1,eb2,xx,tol,z,l,nn,val,slo,v(1,1),g,sum,a,b,rofi,nr,nre,kc)
            call rdeqmu(e,gfn,0,z,v,rofi,nr,nsp,nmsave(l),a,b,l,nl-1,imu,psi,psd,gsmt,GFr)

C       ... Make small parameter p
            do  la1 = 1, 2
              do  la2 = 1, 2
                call prodrw(nr,4,psd(1,1,1,la1),psd(1,1,1,la2),rofi(1,2),smalpa(la1,la2))
              enddo
            enddo
            pprel(4,l,imu,:,:) = smalpa(:,:)

C       ... The rest of the potential parameters
            gmt(:,:) = psi(nr,1,:,:)/rofi(nr,1)
            fmt(:,:) = psi(nr,2,:,:)/rofi(nr,1)
            gmtde(:,:) = psd(nr,1,:,:)/rofi(nr,1)
            fmtde(:,:) = psd(nr,2,:,:)/rofi(nr,1)
C           Current convention: enu has no ves, while c returned with ves added
            e = pprel(5,l,imu,1,1)
            call fdpp(e,e,[NULLR],0d0,gmt,fmt,gmtde,z,rmax,avw,l,nl-1,imu,pprel)
C#elseC
CC           prior to May 2015 : [NULLR] -> sop
C            scalede = rmax*0.5d0*(phisq(l,1) + phisq(l,2))
C            if (pprel(5,l,imu,1,1) == NULLR) pprel(5,l,imu,:,:) = (pp(1,l,1) + pp(1,l,2))/2
C            e = pprel(5,l,imu,1,1)
C            call rdeq_old(e,e,[NULLR],z,v,rofi,nr,2,a,b,l,nl-1,imu,scalede,gmt,fmt,gmtde,fmtde,gsmt,pprel)
C            call fdpp(e,e,[NULLR],0d0,gmt,fmt,gmtde,z,rmax,avw,l,nl-1,imu,pprel)
CC       ... Substitute scalar relativistic small parameter for fully relativistic one
C            if (z /= 0) then
C              do i = 1, 2
C                do j = 1, 2
C                  ptmp4(i,j) = 0d0
C                enddo
C                ptmp4(i,i) = pp(4,l,i)
C              enddo
C              call clebsh(l,imu,ptmp4,ptmp8)
C              do i = 1, 2
C                do j = 1, 2
C                  pprel(4,l,imu,i,j) = ptmp8(i,j)
C                enddo
C              enddo
C            endif
C#endif

            if (z == 0) then ! potential parameters for empty spheres
              do i = 1, 2
                do j = 1, 2
                  ptmp1(i,j) = 0d0
                  ptmp2(i,j) = 0d0
                  ptmp3(i,j) = 0d0
                  ptmp4(i,j) = 0d0
                enddo
                ptmp1(i,i) = pp(2,l,i)
                ptmp2(i,i) = pp(5,l,i)
                ptmp3(i,i) = pp(3,l,i)
                ptmp4(i,i) = pp(4,l,i)
              enddo
              call clebsh(l,imu,ptmp1,ptmp5)
              call clebsh(l,imu,ptmp2,ptmp6)
              call clebsh(l,imu,ptmp3,ptmp7)
              call clebsh(l,imu,ptmp4,ptmp8)
              do i = 1, 2
                do j = 1, 2
                  pprel(1,l,imu,i,j) = ptmp5(i,j)
                  pprel(2,l,imu,i,j) = ptmp6(i,j)
                  pprel(3,l,imu,i,j) = ptmp7(i,j)
                  pprel(4,l,imu,i,j) = ptmp8(i,j)
                  pprel(5,l,imu,i,j) = (pp(1,l,1)+pp(1,l,2))/2
C                 p2 for the empty spheres
C                 pprel(4,l,imu,i,j) = 0d0
                enddo
              enddo
            endif

          enddo
        enddo
      endif

C --- Matrix elements of wave function and its gradient ---
      if (loptc) then
        call rgrme(0,2,2,2,nsp,1,nl-1,lmx,nr,rofi,wgt,phiv,phiv,rgrad)
      endif

C --- Multipole moments of phi, phidot ---
      if (lmpol) then
C   ... Orthonormalize phi,phidot
        call soprm(4,lpzi,phiv,phiv(1,1,1,2),xx,nr,nsp,nl-1,lmx,v,xx,
     .    enu,xx,z,rofi,wgt,facso,sop,xx,' ')
        call mpint(phiv,phiv(1,1,1,2),nl-1,lmx,2*lmx,nr,nsp,rofi,rofi(1,2),pmpol)
      endif

      if (lsor .or. lmpol .or. loptc) then
        deallocate(wgt,dv,wkv,phiv)
      endif
      end
C      subroutine phidot(z,l,v,e,a,b,rofi,nr,g,val,slo,tol,nn,
C     .  gp,phi,dphi,phip,dphip,p)
CC- Generate Phidot,Phidotdot for a prescribed energy
CC ----------------------------------------------------------------
CCi Inputs:
CCi   z:     nuclear charge
CCi   l:     l quantum number for this g
CCi   v:     potential
CCi   a,b:   defines shifted logarithmic mesh (rmesh.f)
CCi   rofi:  list of points; must be consistent with a,b
CCi   nr:    number of mesh points
CCi   e:     Energy
CCi   val:   Value of r * wave function at sphere radius (rseq)
CCi   slo:   Derivative of r * wave function at sphere radius (rseq)
CCi   g:     Wave function times r normalized so that int (g*g) dr = 1
CCi   tol:   precision to which wave function is integrated
CCi   nn:    number of nodes
CCo Outputs:
CCo   gp:    first four energy derivatives to G
CCo   phi:   wave function at rmax, i.e. g/rmax
CCo   dphi:  slope of wave function at rmax, i.e. d(g/rmax)/dr
CCo   phip:  energy derivative of wave function at rmax
CCo   dphip: energy derivative of slope of wave function at rmax
CCo   p:     <gp**2> (potential parameter)
CCr Remarks:
CCr   This version makes energy derivatives by numerical differentiation
CCr   of wave function phi, and has the same calling sequence as the
CCr   analytic version phidot.  The only difference is that this routine
CCr   returns four derivatives of phi, whereas the analytic version only
CCr   returns two and is applicable only to the nonrelativistic case.
CC ----------------------------------------------------------------
C      implicit none
CC passed parameters
C      integer l,nr,nn
C      double precision z,e,a,b,val,slo,phi,dphi,phip,dphip,p,tol
C      double precision v(*),rofi(*)
C      double precision g(2*nr),gp(2*nr,4)
CC local variables
C      integer nre,i,iprint,nptdif
C      double precision rmax,eb1,eb2,dele,ddde,sum1,
C     .                 vali(5),sloi(5),ei(4),de1,de2,del1,del2
C      parameter (nptdif = 2)
C
C      rmax = rofi(nr)
C      eb1 = -50d0
C      eb2 = 20d0
C
Cc      dele = tol**.2d0
C      dele = .002d0
C      if (tol > 1d-9 .and. iprint() >= 0)
C     .  print *, 'phidot: tol too high for reliable num. diff'
C
C      ddde = -rmax/g(nr)**2
C      ei(1) = 1
C      ei(2) = -1
C      ei(3) = 1.5d0
C      ei(4) = -1.5d0
C      do  10  i = 1, nptdif
C        sloi(i) = slo + dele*ei(i)*ddde*val/rmax
C        ei(i) = e + dele*ei(i)
C        call rseq(eb1,eb2,ei(i),tol,z,l,nn,val,sloi(i),v,gp(1,i),
C     .            sum1,a,b,rofi,nr,nre,kc)
C        vali(i) = val/dsqrt(sum1)
C        sloi(i) = sloi(i)/dsqrt(sum1)
C   10 continue
C      de1  = (ei(1) - ei(2))/2
C      del1 = (ei(1) + ei(2))/2 - e
C      de2  = (ei(3) - ei(4))/2
C      del2 = (ei(3) + ei(4))/2 - e
CC     Energy derivatives of value and slope
C      call dfphi(de1,del1,de2,del2,1,val,vali,nptdif == 4)
C      call dfphi(de1,del1,de2,del2,1,slo,sloi,nptdif == 4)
C
C      call dfphi(de1,del1,de2,del2,2*nr,g,gp,nptdif == 4)
C      call gintsr(gp,gp,a,nr,z,e,l,v,rofi,p)
C
CC ... Get phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
C      phi = val/rmax
C      dphi = (slo - phi)/rmax
C      phip = vali(1)/rmax
C      dphip = (sloi(1) - phip)/rmax
C
C      if (iprint() >= 111) print 749, phi,dphi,phip,dphip
C  749 format(' PHIDOT:  phi,phip,phip,dphip=',4f12.6)
C      end
      subroutine gintsr(g1,g2,a,nr,z,e,l,v,rofi,sum)
C- Integrate inner product of two wave equations
C ----------------------------------------------------------------
Ci   g1,g2 :First and second radial wave functions
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   nr    :number of radial mesh points
Ci   z     :nuclear charge
Ci   e     :energy
Ci   l     :l quantum number of g1,g2
Ci   v     :spherical potential
Ci   rofi  :radial mesh points
Co Outputs:
Co   sum   :inner product
Cr Remarks:
Cu   10 Apr 12 Repackaged radial mesh quadrature (radwgt)
Cr   wt should be passed into this routine
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,l
      double precision a,z,e,sum,g1(nr,2),g2(nr,2),v(nr),rofi(nr)
C ... Local parameters
      integer i,intopt,nglob
      double precision fllp1,c,r,tmc,wt(nr)
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c
      tmc(i,r) = c - (v(i) - 2d0*z/r - e)/c

      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rofi(nr),a,nr,wt)
      fllp1 = l*(l+1)
      sum = 0d0
      do  i = 2, nr
        r = rofi(i)
        sum = sum + wt(i)*
     .    (g1(i,1)*g2(i,1)*(1+fllp1/(tmc(i,r)*r)**2) + g1(i,2)*g2(i,2))
      enddo

      end
      subroutine gintsl(g1,g2,a,nr,rofi,sum)
C- Integrate inner product of two wave equations, large component only
C ----------------------------------------------------------------------
Ci Inputs
Ci   g1    :first wave function
Ci   g2    :second wave function
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Co Outputs
Co   sum:   inner product
Cr Remarks
Cr   Like gintsr, but uses large component only (corresponds to c->infty)
Cu Updates
Cu   05 Apr 12 Replaced Simpson rule by weights
Cu   20 Apr 01 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr
      double precision a,g1(nr),g2(nr),rofi(nr)
C ... Local parameters
      integer intopt,nglob
      double precision sum,dot3,wt(nr)

      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rofi(nr),a,nr,wt)
      sum = dot3(nr,wt,g1,g2)
      end
C      subroutine asprjq(mode,clabl,nl,nsp,eula,neul,pnu,qnu,
C     .  pnuloc,qnuloc,bhat,amom)
CC- Find average magnetization axis and project ASA moments onto it
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :1s digit affects bhat
CCi         :0 use bhat as input
CCi         :1 Print out mag. dir. from eula and qnu; make amom
CCi         :2 make bhat from eula and qnu; make amom
CCi         :3 combination of 1+2
CCi         :10s digit affects qnuloc,pnuloc
CCi         :0 do nothing to qnuloc,pnuloc
CCi         :1 copy pnu->pnuloc and qnu->qnuloc
CCi         :2 copy and rotate to bhat coordinates
CCi         :10s digit affects B,qnuloc when magnetization < 0
CCi         :1 scale B by -1 when
CCi   mode  :0 do nothing; just return
CCi   clabl :class name (for printout only)
CCi   nl    :(global maximum l) + 1
CCi   nsp   :2 for spin-polarized case, otherwise 1
CCi   eula  :Euler angles for noncollinear spins
CCi   neul  :1, nl, or nl**2 if Euler angles are: l-independent,
CCi         :l-dependent, or lm-dependent, respectively
CCi   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
CCi          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
CCi   qnu   :energy-weighted moments of the sphere charges
CCio Inputs/Outputs
CCio  bhat  :direction vector for B-field; see 1s digit mode
CCo Outputs
CCo   amom  :projection of magnetization onto bhat
CCo   pnuloc:projection of pnu onto bhat
CCo   qnuloc:projection of qnu onto bhat
CCr Remarks
CCr   If asprjq computes the B-field, it is taken to be parallel to
CCr   the average magnetization, which is computed by summing the
CCr   (vector) magnetization over all orbitals.
CCu Updates
CCu   06 Apr 04 First created
CC ----------------------------------------------------------------------
CC     implicit none
CC ... Passed parameters
C      integer mode,nl,nsp,neul
C      double precision pnu(nl,nsp),qnu(3,nl,nsp),eula(neul,3),
C     .  amom,bhat(3),pnuloc(nl,nsp),qnuloc(3,nl,nsp)
C      character clabl*8
CC ... Local parameters
C      logical lwrite
C      integer i,k,l,m,ilm,stdo,lgunit,ipr,PRT1,PRT2,mode0,mode1,mode2
C      double precision ql,al,alpha,beta,gamma,hatm(3),ploc,qloc(3),
C     .  dploc,aloc(3),dotp,ddot,dsqrt,rotm(3,3),sal
C      character strn*5
C      parameter (PRT1=40,PRT2=50)
C
C      stdo = lgunit(1)
C      call getpr(ipr)
C
CC --- Construct the average magnetization (direction, amplitude) ---
C      mode0 = mod(mode,10)
C      mode2 = mod(mode/100,10)
C      if (mode0 /= 0 .and. nsp == 2) then
C      lwrite = .true.
C      if (mod(mode0,2) == 0 .or. ipr < PRT2) lwrite = .false.
C      if (neul <= 1) lwrite = .false.
C      strn = '    l'
C      if (neul == nl*nl) strn = '  ilm'
C      if (lwrite .and. mode0 >= 2) write(stdo,345) strn
C      if (lwrite .and. mode0 < 2) write(stdo,345) strn,'mhat.bxc'
C  345 format(a,'     ql       mom',9x,'alpha     beta      gamma',
C     .  19x,'mhat':17x,a)
C
C      ilm = 0
C      call dpzero(aloc,3)
C      sal = 0
C      do  l = 0, nl-1
C        do   m = -l, l
C          ilm = ilm+1
C          if (neul == nl) then
C            alpha = eula(l+1,1)
C            beta  = eula(l+1,2)
C            gamma = eula(l+1,3)
C          elseif (neul == nl*nl) then
C            alpha = eula(ilm,1)
C            beta  = eula(ilm,2)
C            gamma = eula(ilm,3)
C          elseif (neul == 1) then
C            alpha = eula(1,1)
C            beta  = eula(1,2)
C            gamma = eula(1,3)
C          else
C            call rxi('atscpp: bad value neul=',neul)
C          endif
C
CC         Charge, magnetic moments, quantization axis for these angles
C          ql = qnu(1,l+1,1)+qnu(1,l+1,2)
C          al = qnu(1,l+1,1)-qnu(1,l+1,2)
C          sal = sal + al/(2*l+1)
CC         The local magnetization points along V = zhat(loc).
CC         In global coordinates V = rotm^-1 zhat(loc) because
CC         rotm*V = zhat(loc) when V points along M as
CC         described in eua2rm.  V is then
CC         hatm(1) = dcos(alpha)*dsin(beta)
CC         hatm(2) = dsin(alpha)*dsin(beta)
CC         hatm(3) = dcos(beta)
C          call eua2rm(alpha,beta,gamma,rotm)
C          hatm(1) = rotm(3,1)
C          hatm(2) = rotm(3,2)
C          hatm(3) = rotm(3,3)
C
CC         Add to total magnetization
C          call daxpy(3,al/(2*l+1),hatm,1,aloc,1)
C
CC         Printout
C          if (neul == nl**2 .and. lwrite .and. mode0 >= 2) then
C            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
C     .        ilm, ql/(2*l+1), al/(2*l+1), alpha, beta, gamma, hatm
C          elseif (neul == nl**2 .and. lwrite .and. mode0 < 2) then
C            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
C     .        ilm, ql/(2*l+1), al/(2*l+1), alpha, beta, gamma, hatm,
C     .        ddot(3,hatm,1,bhat,1)
C          elseif (neul == nl .and. lwrite .and. mode0 >= 2) then
C            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
C     .        l, ql, al, alpha, beta, gamma, hatm
C            lwrite = .false.
C          elseif (neul == nl .and. lwrite .and. mode0 < 2) then
C            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
C     .        l, ql, al, alpha, beta, gamma, hatm, ddot(3,hatm,1,bhat,1)
C            lwrite = .false.
C          endif
C
C        enddo
C        lwrite = .true.
C        if (mod(mode0,2) == 0 .or. ipr < PRT2) lwrite = .false.
C      enddo
C      amom = dsqrt(ddot(3,aloc,1,aloc,1))
C
C      if (amom /= 0) then
CC       Assign bhat to point along magnetization direction
C        if (mode0 >= 2) call dpcopy(aloc,bhat,1,3,1/amom)
C
CC       If sum moments < 0, optionally reverse B
C        if (sal < 0 .and. mode2 == 1 .and. mode0 >= 2) then
C          if (mode0 >= 2) call dpcopy(aloc,bhat,1,3,-1/amom)
C        endif
C
CC       Printout of bhat, moment, angle
C        call info5(PRT1,0,0,' ATOM='//clabl//
C     .    '%a  bhat=%3;10,6D  |M|=%;6,6d  bhat.M/|M|=%;6d',bhat,
C     .    amom,ddot(3,aloc,1,bhat,1)/amom,0,0)
C
C      else
C        bhat(1) = 0
C        bhat(2) = 0
C        bhat(3) = 1
C      endif
C
CC     End of block to contruct bhat
C      endif
C
CC --- Rotate the qnu along projection bhat ---
C      mode1 = mod(mode/10,10)
C      if (mode1 == 0) return
C
C      if (nsp == 1 .or. mode1 == 1 .or.
C     .    ddot(3,bhat,1,bhat,1) == 0d0) then
C        call dcopy(3*nl*nsp,qnu,1,qnuloc,1)
C        call dcopy(1*nl*nsp,pnu,1,pnuloc,1)
C        return
C      endif
C
C      call dpzero(qnuloc,3*nl*nsp)
C      ilm = 0
C      do  l = 0, nl-1
C        call dpzero(aloc,3)
C        dploc = 0
C        do   m = -l, l
C          ilm = ilm+1
C          if (neul == nl) then
C            alpha = eula(l+1,1)
C            beta  = eula(l+1,2)
CC           gamma = eula(l+1,3)
C          elseif (neul == nl*nl) then
C            alpha = eula(ilm,1)
C            beta  = eula(ilm,2)
CC           gamma = eula(ilm,3)
C          elseif (neul == 1) then
C            alpha = eula(1,1)
C            beta  = eula(1,2)
CC           gamma = eula(1,3)
C          else
C            call rxi('atscpp: bad value neul=',neul)
C          endif
C
CC         Charge, magnetic moments, quantization axis for these angles
C          hatm(1) = dcos(alpha)*dsin(beta)
C          hatm(2) = dsin(alpha)*dsin(beta)
C          hatm(3) = dcos(beta)
C          dotp = hatm(1)*bhat(1) + hatm(2)*bhat(2) + hatm(3)*bhat(3)
C
C          do  i = 1, 3
C            aloc(i) = aloc(i) + (qnu(i,l+1,1)-qnu(i,l+1,2))/(2*l+1)*dotp
C          enddo
C          dploc = dploc + (pnu(l+1,1)-pnu(l+1,2))/(2*l+1)*dotp
C
C        enddo
C
C        do  i = 1, 3
C          qloc(i) = qnu(i,l+1,1) + qnu(i,l+1,2)
CC         aloc(i) = qnu(i,l+1,1) - qnu(i,l+1,2)
C          qnuloc(i,l+1,1) = (qloc(i) + aloc(i))/2
C          qnuloc(i,l+1,2) = (qloc(i) - aloc(i))/2
C        enddo
C        ploc = pnu(l+1,1) + pnu(l+1,2)
C        pnuloc(l+1,1) = (ploc + dploc)/2
C        pnuloc(l+1,2) = (ploc - dploc)/2
C
C      enddo
C
C      if (ipr >= PRT2) then
C        write(stdo,'(''  l isp'',19x,''qnu'',37x,''qloc'')')
C        do  k = 1, 2
C        do  l = 0, nl-1
C          write(stdo,'(2i3,1x,3f13.7:1x,3f13.7)')
C     .      l,k,(qnu(i,l+1,k),i=1,3),(qnuloc(i,l+1,k),i=1,3)
C        enddo
C        enddo
C      endif
C
C      end

C      subroutine qnu2enu(mode,nl,lmax,pp,ves,qnur,pprel)
CC- Copy qnur(3,...) to pprel(5,..) or vise-versa
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 copy all qnur(4,:,:,:,:) to pprel(5,:,:,:,:)
CCi         :1 copy all pprel(5,:,:,:,:) to
CCi         :2 copy all pp(1,:,:,:,:) to pprel(5,:,:,:,:)
CCi         :3 conditionally copy qnur(4,:,:,:,:)  to pprel(5,:,:,:,:) if pprel(5,:,:,1,1) = NULL
CCi         :4 conditionally copy pprel(5,:,:,:,:) to qnur(4,:,:,:,:)  if qnur(4,:,:,:,:)  = NULL
CCi         :5 conditionally copy pp(5,:,:)        to pprel(5,:,:,:,:) if pprel(5,:,:,1,1) = NULL
CCi         :6 copy all pp(1,:,:,:,:) to qnur(4,:,:,:,:)
CCi         :7 conditionally copy pp(1,:,:,:,:)    to qnur(4,:,:,:,:) if qnur(4,:,:,:,:) = NULL
CCi   nl    :(global maximum l) + 1
CCi   lmax  :maximum l for a given site
CCi   qnur
CCi   pprel
CCo Outputs
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   09 Apr 15
CC ----------------------------------------------------------------------
C      implicit none
C      integer mode,nl,lmax
C      double precision ves,qnur(4,0:nl-1,2*nl,2,2),pprel(5,0:nl-1,2*nl,2,2),pp(6,0:nl-1,2)
C      integer l,imu
C      integer,parameter :: NULLI=-99999
C
C      do  l = 0, lmax
C        do  imu = 1, 2*l+2
CC          if (mode == 0 .or. mode == 2 .and. pprel(5,l,imu,1,1) == NULLI) pprel(5,l,imu,:,:) = qnur(4,l,imu,1,1)
CC          if (mode == 1 .or. mode == 4 .and.  qnur(4,l,imu,1,1) == NULLI) qnur(4,l,imu,:,:) = pprel(5,l,imu,1,1)
C
C          if (mode == 0 .or. mode == 3 .and. pprel(5,l,imu,1,1) == NULLI) pprel(5,l,imu,:,:) = qnur(4,l,imu,1,1)
C          if (mode == 1 .or. mode == 4 .and.  qnur(4,l,imu,1,1) == NULLI) qnur(4,l,imu,:,:)  = pprel(5,l,imu,1,1) - ves
C          if (mode == 2 .or. mode == 5 .and. pprel(5,l,imu,1,1) == NULLI) pprel(5,l,imu,:,:) = (pp(1,l,1) + pp(1,l,2))/2
C          if (mode == 6 .or. mode == 7 .and.  qnur(4,l,imu,1,1) == NULLI) qnur(4,l,imu,:,:)  = (pp(1,l,1) + pp(1,l,2))/2 - ves
C
C        enddo
C      enddo
C      end
