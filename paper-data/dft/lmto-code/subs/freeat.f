       subroutine freeat(s_ctrl,s_ham,s_pot,s_spec,tolft)
C- For each species, makes free atom self-consistent
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  smalit lxcf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lcd
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  basopt pnudef seref
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name rsmfa rfoca coreq coreh z rmt a nr p q idmod
Ci                 lmxa pz eref rs3 eh3 vmtz rcfa lmxb idxdn orbp norp
Ci                 ntorb
Co     Stored:     z rmt a nr norp ntorb orbp p pz idxdn lmxb
Co     Allocated:  *
Cio    Elts passed:pz name
Cio    Passed to:  gtpcor ioorbp
Ci Inputs
Ci       tolft   : ! HAM_TOL variable
Co Outputs
Cs Command-line switches
Cs   --basfile=fn: Write basp file to fn, instead of default basp0
Cs   --basp      : autofind EH,RSMH (Note: better to use HAM_AUTOBAS)
Cs   --dumprho   : Writes out the density for each atom to out.ext
Cs               : Waits for prompt after each file is writen
Cs   --getallloc : Look for local orbitals (better to use HAM_AUTOBAS)
Cs   --noopt     : Suppress optimization of s.m. Hankel basis
Cs   --norscnst  : In optimization of s.m. Hankel basis,
Cs               : do not constrain rsm < rmt
Cs   --plotwf    : Writes atomic radial wave functions to disk files
Cl Local variables
Cl   ccof  :coefficient to fit of core tail to smoothed Hankel
Cl   ceh   :energy of core tail to smoothed Hankel
Cl   sumtc :core kinetic energy
Cr Remarks
Cu Updates
Cu   29 Jun 17 --basp has new switches ~ctrl and ~rsmmx
Cu             routine fabasp can restrict autogen to undefined RSMH,EH
Cu   18 Dec 16 New options to --basp : --basp~eh1=#~eh2=#~incrlmx
Cu   03 Jun 14 Enable deep local and valence orbitals both to have charge
Cu   14 Nov 13 Some adjustments in preparation for fully relativistic GF
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   27 Sep 12 When optimizating s.m. Hankel basis, check parms stay in bounds
Cu   06 Sep 11 Started migration to f90 structures
Cu   10 May 09 New autogeneration of basis parameters
Cu   01 Feb 06 Enables renormalized free atom density
Cu   01 Jul 05 Skips spheres with Z=0 and R=0
Cu   21 Jun 04 Added fit of sm. Hankel tails to local orbitals
Cu   18 Sep 03 (ATP) Enabled partial core occupation
Cu   06 Sep 03 Constrain rsm in fit to FA wave function
Cu   18 Mar 03 Altered sign of magnetic moment to conform to std
Cu   19 Apr 02 Redesigned freats call to avoid the use of structures
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   22 Mar 01 Added printout of reference energy
Cu   10 Jun 00 spin polarized
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      double precision tolft ! HAM_TOL variable
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
C ... Local parameters
      logical ltmp
      integer, parameter :: n0=10, nkap0=4, nrmx=5001, nxi0=10
      integer autob2,i,ifi,is,kcor,lcor,incrl,lmxa,lxcfun,nmcore,nr,
     .  nrmt,nsp,nspec,nxi,stdo,pnudef,lmxb2(2),iqocc(5)
      integer nrmix(2),idmod(n0)
      character*8 spid,chole*8,strn*80,fn*80
      double precision qc,ccof,ceh,z,rmt,rfoca,rsmfa,qcor(2),a,sumec,
     .  sumtc,eref,seref,etot,xx
      double precision hfc(nxi0,2),exi(nxi0),hfct(nxi0,2)
      double precision v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),rofi(nrmx*2)
      double precision pnu(n0,2),pz(n0,2),qat(n0,2),rcfa(2),gtop(2)
      double precision rtab(n0,2),etab(n0,2),orbp(n0,2,nkap0)
      double precision rs3,eh3,vmtz,basopt(15),pnuopt(n0),pzopt(n0)
      procedure(integer) :: nglob,fopna,iprint,lgunit,isw,iofa
      procedure(logical) :: cmdopt,ioorbp

      ifi = fopna('atm',-1,0)
      rewind ifi
      nspec = nglob('nspec')
      nsp   = nglob('nsp')
      stdo = lgunit(1)
      exi(1) = -1
      exi(2) = -2
      exi(3) = -4
      exi(4) = -6
      exi(5) = -9
      exi(6) = -15
      nxi = 6
      gtop(1) = 0d0
      gtop(2) = 0d0

      call dpzero(hfct,2*nxi0)
      basopt = s_ham%basopt
      pnudef = mod(s_ham%pnudef,10)

      call info5(21,1,0,' FREEAT: G tolerance for envelope functions : %;5g'//
     .  '%?#(n==0)##%N Search criteria for local orbitals:  QLOC: qistl > %;3g  ELOC: eig < %;6d Ry#',
     .  tolft,mod(int(basopt(1))/10,10),basopt(9),basopt(3),5)

C --- For each species, do ---
      do  is = 1, nspec

C       Unpack s_ctrl data
        nrmix = s_ctrl%smalit
        lxcfun = s_ctrl%lxcf
        nmcore = isw(IAND(s_ctrl%lcd,16)/=0)

C       Unpack species data
        spid = s_spec(is)%name
        rsmfa = s_spec(is)%rsmfa
        rfoca = s_spec(is)%rfoca
        qcor = s_spec(is)%coreq
        chole = s_spec(is)%coreh
        call gtpcor(s_spec,is,kcor,lcor,qcor)
        z = s_spec(is)%z
        rmt = s_spec(is)%rmt
        a = s_spec(is)%a
        nrmt = s_spec(is)%nr
        if (z == 0 .and. rmt == 0) goto 10
        call dpzero(qat,n0*2)
        pnu = s_spec(is)%p
        qat = s_spec(is)%q
        idmod = s_spec(is)%idmod
        lmxa = s_spec(is)%lmxa
        pz = s_spec(is)%pz
        eref = s_spec(is)%eref
        rs3 = s_spec(is)%rs3
        eh3 = s_spec(is)%eh3
        vmtz = s_spec(is)%vmtz
        rcfa = s_spec(is)%rcfa

C       These are starting values.  Some elements may be used instead of autogenerated value
        orbp = s_spec(is)%orbp
        do  i = 1, 2
          call dcopy(n0,orbp(1,1,i),1,rtab(1,i),1)
          call dcopy(n0,orbp(1,2,i),1,etab(1,i),1)
        enddo

        call freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,kcor,lcor,qcor,
     .    nmcore,nrmix,1,lxcfun,z,rmt,a,nrmt,pnu,pz,qat,rs3,eh3,vmtz,
     .    rcfa,idmod,lmxa,eref,basopt,pnudef,rtab,etab,hfc,hfct,nr,rofi,
     .    rho,rhoc,qc,ccof,ceh,sumec,sumtc,v,etot,pnuopt,pzopt,gtop,tolft)

C       Repack mesh parms, in case a mesh was selected by freats
        s_spec(is)%z = z
        s_spec(is)%rmt = rmt
        s_spec(is)%a = a
        s_spec(is)%nr = nrmt
C       s_spec(is)%lmxa = lmxa

C   ... Pack the basis into sspec
        call dpzero(orbp,n0*2*nkap0)
C       orbp = s_spec(is)%orbp
C       lmxb = s_spec(is)%lmxb
        do  i  = 1, 2
          call dcopy(n0,rtab(1,i),1,orbp(1,1,i),1)
          call dcopy(n0,etab(1,i),1,orbp(1,2,i),1)
        enddo
        s_spec(is)%norp = 2
        s_spec(is)%ntorb = n0
        s_spec(is)%orbp = orbp

C   --- File write ---
        if (iprint() > 40) write(stdo,230) spid
  230   format(/' write free atom data for species  ',a)
C       Copy second spin channel of rho,v
        if (nsp == 2 .and. nr > nrmt) then
          call dcopy(nrmt,rho(1+nr),1,rho(1+nrmt),1)
          call dcopy(nrmt,rhoc(1+nr),1,rhoc(1+nrmt),1)
          call dcopy(nrmt,v(1+nr),1,v(1+nrmt),1)
        endif
        i = iofa(63-2,spid,nxi0,nxi,exi,hfc,hfct,rsmfa,z,rmt,
     .    a,nrmt,nsp,qc,ccof,ceh,sumtc,rho,rhoc,v,-ifi)

        if (mod(int(basopt(1))/100,10) /= 0) then
          call dcopy(n0,pnuopt,1,pnu(1,1),1)
          call dcopy(n0,pnuopt,1,pnu(1,nsp),1)
          s_spec(is)%p = pnu
        endif

        if (mod(int(basopt(1))/10,10) /= 0) then
          call dcopy(n0,pzopt,1,pz(1,1),1)
          call dcopy(n0,pzopt,1,pz(1,nsp),1)
          s_spec(is)%pz = pz
        endif

   10   continue

      enddo
      call fclose(ifi)

      if (basopt(1) /= 0) then
        if (cmdopt('--basfile=',10,0,fn)) then
        else
          if (cmdopt('--usebasp',9,0,fn)) then
            fn = '--basfile=basp'
          else
            fn = '--basfile=basp0'
          endif
        endif
        call info0(20,1,0,' FREEAT:  writing file '//trim(fn(11:)))
        incrl = 0
        if (cmdopt('--basp',6,0,strn)) call swbasp(strn(7:),0,xx,xx,xx,incrl)
        do  is = 1, nspec
          autob2 = basopt(10)
          call fadflb(mod(autob2,10),99,s_spec(is)%z,lmxb2,iqocc) ! rsmh made to lmxb2
          do  i = lmxb2(1)+1, max(lmxb2(1),s_spec(is)%lmxb)+incrl
            if (s_spec(is)%orbp(i+1,1,1) == 0) then
              s_spec(is)%orbp(i+1,1:2,1) = s_spec(is)%orbp(i,1:2,1)
            endif
          enddo
        enddo
        ifi = fopna(trim(fn(11:)),-1,0)
        rewind ifi
        i = 1
        if (mod(int(basopt(1))/100,10) /= 0) i = i+100000
        if (mod(int(basopt(1))/10,10) /= 0)  i = i+10000
        ltmp = ioorbp(i,2,1,nspec,s_spec,0,-ifi)

        if (gtop(1) /= 0 .and. gtop(2) /= 0) then
          call info2(20,1,0,' FREEAT:  estimate HAM_GMAX from RSMH:'//
     .      '  GMAX=%,1;d (valence)  %,1;d (local orbitals)',
     .      gtop,gtop(2))
        elseif (gtop(1) /= 0) then
          call info2(20,1,0,' FREEAT:  estimate HAM_GMAX from RSMH:'//
     .      '  GMAX=%,1;d',gtop,gtop(2))
        elseif (gtop(2) /= 0) then
          call info2(20,1,0,' FREEAT:  estimate HAM_GMAX from RSMH:'//
     .      '  GMAX=%,1;d (local orbitals)',gtop(2),gtop(2))
        endif

      endif

      if (iprint() > 30)  then
        seref = s_ham%seref
        call awrit1('%x%N Sum of reference energies: %1;6d',' ',80,stdo,seref)
      endif

      end

      subroutine freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,kcor,lcor,qcor,
     .  nmcore,nrmix,lwf,lxcfun,z,rmt,a,nrmt,pnu,pz,qat,rs3,eh3,vmtz,
     .  rcfa,idmod,lmxa,eref,basopt,pnudef,rtab,etab,hfc,hfct,nr,rofi,
     .  rho,rhoc,qc,ccof,ceh,sec,stc,v,etot,pnuopt,pzopt,gtop,tolgv)
C- Makes one free atom self-consistent, fits rho tails to smoothed Hankels
C ----------------------------------------------------------------------
Ci Inputs
Ci   spid  :species label (for printout)
Ci   is    :species index
Ci   nxi0  :nxi0: leading dimension of hfc,hfct
Ci   nxi   :number of hankel functions used in fitting of tails
Ci   exi   :hankel energies used in fitting of tails
Ci   rfoca :smoothing radius for hankel fit to core
Ci   rsmfa :smoothing radius for hankel fit to valence
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   nrmix :nrmix(1) = maximum number of interations in sphere
Ci         :           before giving up on self-consistency
Ci         :nrmix(2) = no prior iterations Anderson mixing in
Ci         :           self-consistency cycle.
Ci   lwf   :1 print information about wave functions
Ci   lxcfun:selects LDA+GGA exchange-correlation functional
Ci   z     :nuclear charge
Ci   rmt   :augmentation radius, in a.u.
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nrmt  :number of mesh points from origin to rmt
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pz    :boundary conditions for local orbitals
Ci   qat   :valence charges for each l channel
Ci   rs3   :minimum allowed smoothing radius in attaching Hankel tails
Ci         :to local orbitals
Ci   eh3   :Hankel energy when attaching Hankel tails to high-lying
Ci         :local orbitals
Ci   vmtz  :parameter used in attaching Hankel tails to local orbitals
Ci         :It is used as a constant shift to Hankel energies for the
Ci         :fitting of local orbitals to Hankel tails. Thus vmtz
Ci         :is an estimate for the potential at the MT radius.
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
Ci   lmxa  :augmentation l-cutoff
Ci   eref  :reference energy (used for printout)
Ci   basopt:Parameters used generation of LMTO basis parms
Ci         : 1 = autob
Ci         :     1s   digit 1 or 3 Autogenerate RSMH,EH
Ci         :                2 or 4 Autogenerate RSMH,EH, RSMH2,EH2
Ci         :     10s  digit Find and estimate PZ from free atom wf
Ci         :                0 do not modify PZ
Ci         :                1 or 2 Set PZ satisfying criteria
Ci         :                1 PZ specified by caller takes precedence
Ci         :                2 PZ satisfying criteria is always autogenerated
Ci         :     100s digit Estimate P from free atom wf
Ci         :                0 do not modify P
Ci         :                1 Autogenerate P for l<=lmxb
Ci         :                2 Autogenerate P for l<=lmxa
Ci         : 2 = global cutoff to lmxb, 1st kappa ... no longer used
Ci         : 3 = elocmx: Set local orbital PZ set when E>elocmx
Ci         : 4 = rsmmx : maximum rsm, in units of rmt
Ci         : 5 = ehmx  : maximum eh, Ry.  Not used if ehmx = NULLI
Ci         : 6 = esprd : default spread in EH,EH2
Ci         : 7 = modeV : specifies a mode to modify potential
Ci         : 8 = vbar  : parameter in V modification
Ci         : 9 = pqrmx : Set local orbital PZ set when q(r>rmt)>pqrmx
Ci         :10 = 1s  digit LMTO basis parameters
Ci         :           : 0 standard basis, typically spd
Ci         :           : 1 hyperminimal basis
Ci         :           : 2 hyperminimal basis + 1 higher l
Ci         :           : 3 same as 0
Ci         :           : 4 increment standard basis by adding 1 l
Ci         :     10s digit
Ci         :           : 1 Set defaults for GW calculation
Ci         :     100s digit
Ci         :           : 0 Traditional mode
Ci         :           : 1 Jerome's experimental mode
Ci         : when basopt(10)'s 100s digit set (function fabaspj)
Ci         : 11 = deep energy for V(r)/ctp basis scheme
Ci         : 12 = shallow energy for V(r)/ctp basis scheme
Ci         : 13 = Hankel energy for V(r)/ctp basis scheme
Co Outputs
Co   rtab  :smoothing radius for optimized wave function
Co   etab  :energy for optimized wave function
Co   hfc   :fit coeffs for valence density,
Co   hfct  :contains fit coeffs for full density (not calc. now)
Co   nr    :number of radial mesh points for spherical rho
Co   rofi  :rofi(1..nr)=radial mesh for points
Co         :rofi(nr+1..2*nr)=radial mesh weights
Co   rho   :free-atom valence density
Co   rhoc  :free-atom core density
Co   qc    :Sphere core charge
Co   ccof  :coefficient to fit of core tail to smoothed Hankel
Co   ceh   :energy of core tail to smoothed Hankel
Co   sec   :sum of core eigenvalues
Co   stc   :core kinetic energy
Co   v     :spherical potential
Cl Local variables
Cl   itab  :itab(l+1)=1  a wave function was optimzed for this l
Cl   pnul  :EITHER : pnu for valence state, OR
Cl         :local orbital if DEEPER than valence (pz<pnu)
Cr Remarks
Cu Updates
Cu   03 Apr 18 (J. Jackson) re-run atom solver when semicore LO reduce Qc
Cu             atm file now has Qc corresponding to SC-LO case
Cu   26 Mar 18 (J. Jackson) call function fabaspj for basis setup **experimental**
Cu             adds basopt(10:13) for specific variables
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 May 09 New call to fabasp, autogeneration of basis parameters
Cu   01 Feb 06 Enables renormalized free atom density
Cu   19 Apr 02 Redesigned input to avoid the use of structures
Cu   10 Apr 02 Redimensionsed etab,rtab to accomodate larger lmax
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nrmt,is,nxi0,nxi,nrmix(2),lwf,lxcfun,kcor,lcor,pnudef,nmcore
      integer, parameter :: nrmx=5001, n0=10, ncmx=200, nvmx=20, NULLI=-99999
      character*8 spid
      double precision rsmfa,rfoca,qc,ccof,ceh,sec,stc,z,rmt,a,eref,basopt(15),
     .  v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),hfc(nxi0,1),hfct(nxi0,1),qat(n0,2),
     .  exi(1),rtab(n0,2),etab(n0,2),rofi(nrmx*2),rs3,eh3,vmtz,rcfa(2),qcor(2)
      double precision pnu(n0,2),pz(n0,2),pnuopt(n0),pzopt(n0),gtop(2),tolgv
C ... Dynamically allocated arrays
      real(8), allocatable :: g(:),psi(:)
C ... Local parameters
      integer autob,autob2,autopz,i,ifi,incrl,ipr,isp,l,lfrz,lmxa,lplawv,lplfa,lrel,lxcg,
     .  modeV,nitmax,nmix,nr,nsp,stdl,stdo
      integer idmod(n0),irchan(n0),lmxb2(2),iqocc(5)
      character str*8,strn*64
      double precision rmax,b,etot,dq,sumev,ekin,utot,rhoeps,amgm,rhrmx,qvt,
     .  qtot,qct,qvin,qcin,r,wt0,wt1,qtt,qtin,pnul,pzl,vbar,elocmx,pqrmx,
     .  gtopl,xx,eh0,eh2,qcstart
      double precision ec(ncmx),ev(nvmx),etabl(n0,2),exrmax(2),pl(n0,2),pzloc(n0),
     .  qatl(n0,2),ql(3,n0,2),rhoin(nrmx*2),rhot(nrmx*2),rtabl(n0,2),vrmax(2)
      double precision qlm(0:2,n0*n0,2),qnur(0:3,0:n0-1,2*n0,2,2)
C     double precision sum(0:2)
C     Parameters for optimized free-atom wave functions
      double precision ehmx,rsmmx,esprd,bscal
      integer itab(n0,2)
      procedure(logical) :: cmdopt
      procedure(integer) :: nglob,lgunit,fopna,iprint,isw,cmdoptswx

      stdo   = lgunit(1)
      stdl   = lgunit(2)
      ipr    = iprint()
      lrel   = mod(nglob('lrel'),10)
      nsp    = nglob('nsp')
      lfrz   = 0 + 2*nmcore
      lxcg   = lxcfun/100
      autob  = basopt(1)
      elocmx = basopt(3)
      rsmmx  = basopt(4)
      ehmx   = basopt(5)
      esprd  = basopt(6)
      modeV  = basopt(7)
      vbar   = basopt(8)
      pqrmx  = basopt(9)
      autob2 = basopt(10)
      !tolgv  = 5d-6 ! HAM_TOL now passed as argument
      call dpzero(pzloc,n0)
      if (z == 0 .and. lmxa == -1) then
        call dcopy(n0,pz,1,pzopt,1)
        return
      endif

      if (ehmx == NULLI) then
        ehmx = -0.1d0
        if (mod(autob2/10,10) > 0) ehmx = -0.3d0
        if (mod(autob2/10,10) > 0 .and. z <= 10) ehmx = -0.4d0
      endif

C ... Increase lmxa if less than lmxb
      if (mod(autob,10) /= 0 .or. cmdopt('--basp',6,0,strn)) then
        call fadflb(mod(autob2,10),99,z,lmxb2,iqocc)
        i = max(lmxb2(1),lmxb2(2))
        if (lmxb2(1) > lmxa) then
          call info2(2,1,0,' lmfa species '//spid//
     .      '%a:  lmxa=%i < autogen lmxb (=%i).'//
     .      '  Increase file LMXA, or STRUC_NL',lmxa,lmxb2(1))
          call rx('lmfa')
        endif
      endif

C --- Get species data ----
      call dpzero(pl,2*n0)
      call dpzero(ql,2*3*n0)
      qatl = qat
      do  i = 1, nsp
        do  l = 0, lmxa
          pnul = pnu(l+1,i)
          pzl  = mod(pz(l+1,1),10d0)
          if (pzl /= 0) then
C           Make pz the valence state: hardwire qat
            if (int(pnul-1) == int(pzl)) then
              if (qat(l+1,1) /= 0) then ! require 2 valence p.q.n's
C     call info5(21,1,0,' freeat, l=%i:'//
C     .        '  PZ=%d(Q=%d) and P=%;3d(Q=%d) both valence',
C     .        l,pzl,dble(4*l+2),pnul,qat(l+1,i))
              else
                pnul = pzl
                qatl(l+1,1) = 4*l+2
                qatl(l+1,2) = 0
              endif
            elseif (int(pnul+1) /= int(pzl)) then
              call fexit3(-1,111,' Exit -1 freeat, l=%i:  '//
     .          'sc PZ=%d incompatible with valence P=%;3d',l,pzl,pnul)
            endif
          endif
          pl(l+1,i) = int(pnul) + .5d0
          ql(1,l+1,i) = qatl(l+1,1)/nsp
          if (nsp == 2) then
            if (pnul == 0) pl(l+1,i) = pl(l+1,1)
            ql(1,l+1,i) = qatl(l+1,1)/nsp - qatl(l+1,2)/2*dble(2*i-3)
          endif
          ql(2,l+1,i) = 0d0
          ql(3,l+1,i) = 0d0
        enddo
      enddo
      call getqvc(nsp,n0,lmxa,z,pl,ql,pz,0,0,kcor,lcor,qcor,
     .            qc,qtot,amgm,xx,xx)
      if (ipr >= 20) then
        call awrit6('%N Species '//spid//'%a:  Z=%d'//
     .    '  Qc=%d  R=%1,6;6d  Q=%1;6d%?#n==2#  mom=%1;5d#%0d#%a',
     .    ' ',80,stdo,z,qc,rmt,qtot,nsp,amgm)
        do  l = 0, lmxa
          pzl  = mod(pz(l+1,1),10d0)
          if (pzl>0 .and. int(pzl)==int(pl(l+1,1)-1)) then
              call info5(21,0,0,'%8fl=%i:'//
     .        '  PZ=%d(Q=%d) and P=%;3d(Q=%d) both in valence',
     .        l,pzl,dble(4*l+2),pl(l+1,1),qatl(l+1,1))
          endif
       enddo
      endif

C --- Set up radial mesh for free atom ---
      rmax = 50d0
      if (z <  10) rmax = 25
      if (z <=  6) rmax = 20
      call pshpr(1)
      call rmesh(z,rmt,lrel,lxcg,nrmx,a,nrmt)
      call poppr
      b = rmt/(dexp(a*nrmt-a)-1d0)
      nr = 1d0+dlog(1d0+rmax/b)/a
      if (mod(nr,2)==0) nr = nr-1
      rmax = b*(dexp(a*(nr-1))-1d0)
      call info5(21,0,0,' mesh:   rmt=%,6;6d  rmax=%,6;6d'//
     .  '  a=%d  nr=%i  nr(rmax)=%i',rmt,rmax,a,nrmt,nr)

C --- Make atom self-consistent ---
      nitmax = nrmix(1)
      nmix = nrmix(2)
      nmix = -30
C     call pshpr(min(iprint(),40))
      ec(1) = 0
      call dpzero(qnur,size(qnur))
      bscal = 1d0  ! We don't have a number for a species now

      if (lrel == 2) then
        call qnu2qnur(1,n0,nsp,ql,qlm,qnur) ! Only 0th moment for free atom
      endif

      call atomsc(0,n0,nsp,lmxa,z,0d0,kcor,lcor,qcor,rmax,a,nr,
     .  rofi,ec,ev,pl,ql,qnur,pz,bscal,idmod,v,rhoin,rho,rhoc,nmix,
     .  qc,sec,stc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,dq,
     .  exrmax,'gue',nitmax,lfrz,1d0)
C     call poppr

C --- setup SC-LO immediately because these affect Qcore
      if (mod(autob/10,10) > 0) then !search for LO
        autopz = mod(autob/10,10) + 10
        call getsclo(autopz,pqrmx,elocmx,z,a,nr,b,nrmt,lmxa,pl,ql,pz,nsp,v,rofi)

C   ... Repartition core-valence charge owing to change in SCLO
        if (cmdopt('--usebasp',9,0,strn)) then ! possibly update core density

        qcstart = qc
        call getqvc(nsp,n0,lmxa,z,pl,ql,pz,0,0,kcor,lcor,qcor,qc,qtot,amgm,xx,xx)

        if (qc /= qcstart) then
          call info2(21,1,0,' Core charge changed from %d to %d: remake (core,valence) density',qcstart,qc)
          call pshpr(iprint()-30)        !avoid printout from repeated atom solve
          call atomsc(0,n0,nsp,lmxa,z,0d0,kcor,lcor,qcor,rmax,a,nr,rofi,ec,ev,
     .      pl,ql,qnur,pz,bscal,idmod,v,rhoin,rho,rhoc,nmix,qc,sec,stc,sumev,
     .      ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,dq,exrmax,'gue',nitmax,lfrz,1d0)
          call poppr
        end if
        end if ! lusebasp
      end if ! Search for LO

      call info8(20,0,0,'%?#n>=30#%N## sumev=%,6;6d  etot=%,6;6d'//
     .  '  eref=%,6;6d%?#n#  diff= %,6;6d',ipr,sumev,etot,
     .  eref,isw(eref/=0),etot-eref,7,8)

      call dcopy(lmxa+1,ql,3,qatl,1)
      call dcopy(lmxa+1,ql(1,1,nsp),3,qatl(1,2),1)
      call awrit4('fa  Pl %n:-1d  Ql %n:-1d',' ',120,stdl,lmxa+1,pl,lmxa+1,qatl)
      if (nsp == 2) call awrit4('fa  Pl2%n:-1d  Ql2%n:-1d',
     .  ' ',120,stdl,lmxa+1,pl(1,nsp),lmxa+1,qatl(1,nsp))

      if (dabs(dq) > 1d-5) call info2(10,0,0,' freeat (warning) atom not neutral, Q=%d',dq,2)

C .. Subtract core from density to make valence density
      do  isp = 1, nsp
      do  i = 1, nr
        rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)
        rho(i+(isp-1)*nr)  = rho(i+(isp-1)*nr)-rhoc(i+(isp-1)*nr)
      enddo
      enddo

C --- Renormalize atom density or potential ---
      call ivset(irchan,1,n0,0)
      call rnatm(pl,qatl,n0,irchan,lmxa,z,a,b,rofi,ev,nr,rcfa,nsp,v,rho)
C     call prrmsh('starting total rho',rofi,rhot,nr,nr,nsp)
      do  isp = 1, nsp
      do  i = 1, nr
        rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)+rhoc(i+(isp-1)*nr)
      enddo
      enddo
C     call prrmsh('ending total rho',rofi,rhot,nr,nr,nsp)

C --- Print info about free-atom wavefunctions ---
      if (lwf /= 0) then
      if (ipr > 30) then
        allocate(g(nr*2),psi(nr*(lmxa+1)*nsp))
        lplawv = 0
        if (ipr >= 50) lplawv = 1
        call pratfs(spid,lplawv,z,a,nr,rmax,nrmt,lmxa,pl,nsp,v,rofi,
     .    g,psi)
        deallocate(g,psi)
      endif

C --- Optimise smooth-Hankel basis ---
      call dvset(rtabl,1,n0,-1d0)
      call dpzero(etabl,n0)
      i = 1
      if (.not. cmdopt('--noopt',7,0,strn)) then
      if (cmdopt('--norscnst',10,0,strn)) i = 0
      if (z > 0) then
        call optfab(i,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,nsp,v,
     .    rofi,spid,itab,rtabl,etabl)
C   ... Fit value and slope of local orbitals
        call ftfalo(i,z,a,nr,rmax,nrmt,rmt,lmxa,pnu,pz,rs3,eh3,vmtz,
     .    nsp,v,rofi,spid,tolgv,gtopl)
        gtop(2) = max(gtop(2),gtopl)
      endif
      endif
      endif

C --- Autogenerate basis parameters EH,RSMH,EH2,RSMH2 ---
      if (mod(autob,10) /= 0 .or. cmdopt('--basp',6,0,strn)) then
        call fadflb(mod(autob2,10),lmxa,z,lmxb2,iqocc)
C       Use true valence pnu (instead of pz when local orbital present)
        do  l = 0, lmxa
          pnul = pnu(l+1,1)
          pl(l+1,1) = int(pnul) + .5d0
        enddo
        call dpzero(pnuopt,n0)
        if (mod(autob/100,10) /= 0) call dcopy(lmxa+1,pnu,1,pnuopt,1)
        eh0 = NULLI; eh2 = NULLI
        if (cmdopt('--basp',6,0,strn)) call swbasp(strn(7:),autob,eh0,eh2,rsmmx,incrl)
        if (mod(autob2/100,10) == 0) then ! Traditional mode
          call fabasp(autob,pnudef,lmxb2,modeV,vbar,eh0,eh2,0.8d0,-2d0,rsmmx*rmt,
     .    ehmx,elocmx,pqrmx,esprd,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,pz,
     .    nsp,v,rofi,spid,itab,rtab,etab,pnuopt,pzopt,tolgv,gtopl)
        elseif (mod(autob2/100,10) == 1) then ! Jerome's experimental mode
C         provide also the HAM_AUTOBAS_GW flag for choosing a deeper
C         default EH
          call vrbasp(autob,z,a,nr,b,nrmt,lmxa,pl,ql,pz,
     .      nsp,v,rofi,itab,rtab,etab,pnuopt,pzopt,tolgv,gtop,basopt(11),basopt(12),basopt(13))
        else
          call rx('unknown autobas mode')
        endif
        gtop(1) = max(gtop(1),gtopl)
        if (mod(autob2,10) == 1) then
          do  l = 0, max(lmxb2(1),lmxb2(2))
            if (iqocc(l+1) == 0) then
              rtab(l+1,1) = 0
              rtab(l+1,2) = 0
            endif
          enddo
        endif
      endif

C --- Print charges within/outside MT sphere ---
      if (z > 0) then
      qvt = 0d0
      qct = 0d0
      qvin = 0d0
      qcin = 0d0
      do  isp = 1, nsp
      do  i = 1, nr
        r = rofi(i)
        wt0 = rofi(i+nr)
C       wt0 = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
        wt1 = wt0
        if (i > nrmt) wt1 = 0d0
        if (i==1 .or. i==nrmt) wt1 = a*(r+b)/3d0
        if (i==1 .or. i==nr)   wt0 = a*(r+b)/3d0
        qvin = qvin + wt1*rho(i+(isp-1)*nr)
        qcin = qcin + wt1*rhoc(i+(isp-1)*nr)
        qvt = qvt + wt0*rho(i+(isp-1)*nr)
        qct = qct + wt0*rhoc(i+(isp-1)*nr)
      enddo
      enddo
      qtt = qvt + qct
      qtin = qvin + qcin
      if (ipr >= 40) write (stdo,550) qvin,qvt-qvin,qvt,qcin,qct-qcin,
     .  qct,qtin,qtt-qtin,qtt
  550 format(/' Charges:     inside',7x,'outside',7x,'sum'
     .   /' valence',3f13.6/' core   ',3f13.6/' total  ',3f13.6)

      write (stdl,710) z,rmax,qc,qct-qcin,dq,etot
  710 format('fa Z',f6.1,'   rm',f7.2,'  qc',f6.2,'  qspl',f8.5,
     .   '  dq',f8.5,'  Etot',f15.6)
      else
        qvt = 0
      endif

C      call radsum(nr,nr,1,nsp,rofi(nr+1),rhoc,sum(0))
C      call radsum(nr,nr,1,1,rofi(nr+1),rhoc,sum(1))
C      sum(2) = sum(0) - sum(1)
C      print *, sum(1:nsp)
C      print *, qct
C     call rx('done')

C --- Attach smooth Hankel tails to valence density ---
C     lplfa = nglob('lplfa')
      lplfa = 0
      if (qvt > 1d-6) then
        if (lplfa == 1) then
          write(stdo,344)
  344     format(/' write plot file with valence density..')
          if (is <  10) write (str,'(''pl'',i1)') is
          if (is >= 10) write (str,'(''pl'',i2)') is
          ifi = fopna(str,-1,0)
          write (ifi,490) spid,rmt,rsmfa,nxi
  490     format('# fit to fa density: ',a/
     .       '# rmt=',f7.3,'   rsm=',f7.3,'   nxi=',i2)
          call fclose(ifi)
        endif
        call tailsm(0,nr,nrmt,nsp,a,b,rmt,rsmfa,nxi0,nxi,exi,rofi,
     .     rho,rhot,hfc,hfct)
C       call prrmsh('rho-fa',rofi,rho,nr,nr,1)
      else
        call dpzero(hfc, nxi0*nsp)
        call dpzero(hfct, nxi0*nsp)
      endif

C --- Fit analytical expression to tail of core density ---
      call fctail(nr,nrmt,a,b,rfoca,rofi,rhoc,ccof,ceh)

      end

      subroutine pratfs(spid,lplawv,z,a,nr,rmax,nrmt,lmaxa,pl,nsp,v,
     .  rofi,g,psi)
C- Prints out core and valence energy levels of free-atom
C ----------------------------------------------------------------------
Ci Inputs
Ci   spid  :species label
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rmax  :free-atom radius, in a.u.
Ci   nrmt  :mesh for MT radius
Ci   lmaxa :muffin-tin cutoff
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   g     :normalized wave function times r (work array)
Co Outputs
Co   psi   :normalized wave functions for each l
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n0,nr,lmaxa,nsp,nrmt,lplawv
      parameter (n0=10)
      double precision z,a,rmax,pl(n0,nsp),v(nr,nsp),rofi(nr),
     .  g(2*nr),psi(nr,0:lmaxa,nsp)
C ... Local parameters
      integer fopna,isp,l,konfig,nn,nre,i,konf,kc,
     .  konfg(0:8),ifi,stdo,lmaxc
      double precision ev(0:20),pi,b,tol,eb1,eb2,dl,val,slo,sum,pzero,
     .  pmax,ctp,ecor,rmt
      procedure(integer) :: nglob
      character*1 lsym(0:n0-1), cc, str*15, spid*8
      data lsym /'s','p','d','f','g','5','6','7','8','9'/

      stdo = nglob('stdo')
      pi   = 4d0*datan(1d0)
      if (lmaxa>n0-1) call rx('pratfs:  lmax too large')
      b = rmax/(dexp(a*nr-a)-1d0)
      rmt = b*(dexp(a*nrmt-a)-1d0)
      tol = 1d-8
      write(stdo,580)
  580 format(/' Free-atom wavefunctions:')

      do  isp = 1, nsp

C --- Valence states ---
      if (isp == 1) write(stdo,401)
      if (isp == 2) write(stdo,'(/'' spin 2:'')')
      eb1 = -50d0
      eb2 =  50d0
      do  l = 0, lmaxa
        konfig = pl(l+1,isp)
        dl = dtan(pi*(0.5d0-pl(l+1,isp)))
        nn = konfig-l-1
        ev(l) = -0.5d0
        val = rmax
        slo = dl+1
        if (rmax > 9.99d0) then
          val = 1d-30
          slo = -val
        endif
        call rseq(eb1,eb2,ev(l),tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,
     .    rofi,nr,nre,kc)
C       call gintsr(g,g,a,nr,z,ev(l),l,v,rofi,sum) ! Should be normalized
        call gintsl(g,g,a,nr,rofi,sum)    ! Small component not important
        call gintsl(g,g,a,nrmt,rofi,pmax) ! for large r
        sum = sum - pmax
        call ppratf(ev(l),z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
        cc = ' '
        if (dabs(ctp-rmax) < 1d-3) cc = '*'
        write(stdo,400) konfig,lsym(l),ev(l),pzero,pmax,ctp,cc,sum
  400   format(i4,a1,f14.5,2x,3f12.3,a,f12.6)
  401   format(' valence:',6x,'eval',7x,'node at',6x,'max at',7x,
     .    'c.t.p.   rho(r>rmt)')

C   ... Copy valence wavefunction to psi
        do  i = 1, nr
          psi(i,l,isp) = g(i)
        enddo
      enddo

C --- Core states ---
      write(stdo,403)
      eb1 = -2.5d0*z*z-5d0
      eb2 = 50d0
      call config(pl,lmaxa,z,konfg,lmaxc)

      do  konf = 1, 8
      do  l = 0, min(konf-1,lmaxc)
      konfig = konfg(l)
      if (konf >= konfig) cycle
      nn = konf-l-1
      ecor = -50d0
      val = 1d-30
      slo = -val
      call rseq(eb1,eb2,ecor,tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,rofi,
     .  nr,nre,kc)
C     call gintsr(g,g,a,nr,z,ecor,l,v,rofi,sum) ! Should be normalized
      call gintsl(g,g,a,nr,rofi,sum)    ! Small component not important
      call gintsl(g,g,a,nrmt,rofi,pmax) ! for large r
      sum = sum - pmax
      call ppratf(ecor,z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
C      write(stdo,400) konf,lsym(l),ecor,pzero,pmax,ctp,' ',sum
  403 format(/' core:        ecore',7x,'node at',6x,'max at',
     .   7x,'c.t.p.   rho(r>rmt)')
      call info8(10,0,0,'%,4i'//lsym(l)//'%;14,5D  %;12,3D%;12,3D%;12,3D %;12,6D',
     .  konf,ecor,pzero,pmax,ctp,sum,0,0)
      enddo
      enddo
      enddo

C --- Write file with valence wavefunctions
      if (lplawv == 1) then
        write (str,'(''wf_'',a)') spid
        write (stdo,344) str
  344   format(/' Write valence wavefunctions to plot file: ',a)
        ifi = fopna(str,-1,0)
        write (ifi,490) spid,rmax,rmt,nr,lmaxa,nr,1+nsp*(lmaxa+1)
  490   format('# Free-atom wavefunctions (divided by r) for species ',
     .     a/'# rmax=',f7.3,'   rmt=',f7.3,'   nr=',i5,'   lmax=',i3/
     .    '% rows ',i5,' cols ',i3)
        do  50  i=1,nr
          write (ifi,495) rofi(i),((psi(i,l,isp),l=0,lmaxa),isp=1,nsp)
  495     format(f9.5,1p,16d14.5)
   50   continue
        call fclr(str,ifi)
      endif

      end

      subroutine ppratf(e,z,nr,nre,rofi,a,b,v,g,pzero,pmax,ctp)
C- Find outermost node and maximum of wavefct
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :wave function eigenvalue
Ci   z     :nuclear charge
Ci   nr    :number of radial mesh points
Ci   nre   :last point for which wf is calculated
Ci   rofi  :radial mesh points
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   v     :spherical potential (atomsr.f)
Ci   g     :normalized wave function times r
Co Outputs
Co   pzero :outermost node
Co   pmax  :outermost maximum
Co   ctp   :classical turning point
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nre
      double precision a,b,ctp,e,pmax,pzero,z,rofi(nr),v(nr),g(nr)
C ... Local parameters
      integer i,ir
      double precision g1,g2,rho1,rho2,rho3,x

C ... Find the classical turning point
      do  20 i = nr-1, 5, -1
        ir = i
        if (e > v(i)-2d0*z/rofi(i)) goto 21
   20 continue
   21 g1 = e-v(ir) + 2d0*z/rofi(ir)
      g2 = e-v(ir+1) + 2d0*z/rofi(ir+1)
      ctp = rofi(nr)
      if (g1*g2 < 0d0) ctp = (rofi(ir)*g2-rofi(ir+1)*g1)/(g2-g1)

C ... Find the outermost node
      do  10  i = nre-1, 5, -1
        ir = i
        if (g(i)*g(i+1) < 0d0) goto 11
   10 continue
   11 continue
      pzero = 0d0
      g1 = g(ir)
      g2 = g(ir+1)
      if (ir > 5) pzero = (rofi(ir)*g2-rofi(ir+1)*g1)/(g2-g1)

C ... Find the outermost maximum
      do  30  i = nre-2, 5, -1
        ir = i
        rho1 = g(i)*g(i)
        rho2 = g(i+1)*g(i+1)
        rho3 = g(i+2)*g(i+2)
        if (rho1 < rho2) goto 31
   30 continue
   31 pmax = 0
      if (ir > 5) then
        x = -0.5d0*(rho3-rho1)/(rho1+rho3-2*rho2)
        pmax = b*(dexp(a*(ir+x))-1d0)
      endif
      end

      subroutine optfab(isw,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,nsp,
     .  v,rofi,spid,itab,rtab,etab)
C- Optimise a minimal smooth-Hankel basis for the free atom.
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :1 constrain rsm to be <= rmt
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rmax  :muffin-tin radius, in a.u.
Ci   nrmt  :number of points between 0..rmt
Ci   rmt   :muffin-tin radius, in a.u.
Ci   lmxa  :muffin-tin l-cutoff
Ci   pl    :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   ql    :sphere moments
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   spid  :species label
Co Outputs
Co   itab  :itab(l+1)=1  optimized wave function was found for this l
Co   rtab  :smoothing radius for optimized wave function
Co   etab  :energy for optimized wave function
Cl Local variables
Cr Remarks
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmxa,nr,nrmt,nsp,n0,isw
      parameter (n0=10)
      integer itab(n0,2)
      double precision a,rmax,rmt,z
      double precision rofi(*),v(nr,nsp),pl(n0,nsp),ql(3,n0,nsp)
      character spid*8
C ... Dynamically allocated arrays
      real(8), allocatable :: h(:),g(:),gp(:)
      real(8), allocatable :: psi(:)
C ... Local parameters
      character strn*80
      double precision rtab(n0,2),etab(n0,2)
      integer ipr,irep,isp,istife,istifr,jpr,konfig,l,
     .  lplawv,lrel,nn,nrep,stdo,stdl
      double precision b,deh,deh0,dphi,dphip,drsm,drsm0,e1,e2,e3,eadd,
     .  eaddx,eh,elim1,elim2,enew,enu,eval,p,phi,phip,pnu,qvl,radd,
     .  raddx,rlim1,rlim2,rnew,rsm,stife,stifr,sume1,sume2,qrmt
      procedure(integer) :: nglob,iprint,lgunit
      procedure(logical) :: cmdopt

      ipr = iprint()
      stdo = lgunit(1)
      stdl = lgunit(2)
      lrel = mod(nglob('lrel'),10)
      if (z < 0.99d0) lrel = 0
      if (lmxa > 8) call rx('optfab:  lmax too large')
      b = rmax/(dexp(a*nr-a)-1d0)
      allocate(h(nr),g(2*nr),gp(2*nr*4))

      do  isp = 1, nsp

      call info5(20,0,0,'%?#n>=30#%N##'//
     .    ' Optimise free-atom basis for species '//spid//
     .    '%a, %?#n>1#spin 2#rmt=%;7g#',ipr,isp,rmt,0,0)

C --- Parameters for minimisation ---
      drsm0 = 0.1d0
      rlim1 = 0.3d0
      rlim2 = 2*rmt
      raddx = 0.2d0

      deh0  = 0.05d0
      elim1 = -5.0d0
      elim2 = -0.10d0
C     elim2 = -0.20d0
      eaddx = 0.099d0
      jpr=0
      if (ipr >= 50) jpr=1

C --- Loop over bound valence states ---
      if (ipr > 20) write (stdo,261)
      sume1 = 0d0
      sume2 = 0d0
      do  l = 0, lmxa
        itab(l+1,isp) = 0
        konfig = pl(l+1,isp)
        nn = konfig-l-1
        qvl = ql(1,l+1,isp)
c   ... get exact fa wavefunction, eigval, pnu at rmt
        call popta3(0,l,z,nn,nr,nrmt,rofi,v(1,isp),a,b,eval,pnu,g)
        if (eval > 0d0) cycle
        sume1 = sume1 + qvl*eval
C   ... Potential parameters at MT sphere
        call popta4(l,z,rmt,nrmt,rofi,v(1,isp),g,gp,
     .    a,b,pnu,enu,p,phi,dphi,phip,dphip)
        rsm = rmt
        eh = -1
        if (jpr > 0) write (stdo,340)
  340   format('  L   parin    aux      E1       E2       E3',
     .     '       stiff    Eout     parout')
        do  irep = 1, 50
          nrep = irep
C     ... Get center energy
          call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e2,qrmt)
C     ... Vary rsm
          drsm = drsm0
          call popta1(rsm+drsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
          call popta1(rsm-drsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
          call popta2(l,rsm,eh,drsm,e1,e2,e3,rlim1,rlim2,raddx,rnew,
     .       stifr,jpr)
C     ... Vary eh
          deh = deh0
c         if (eh+deh>-0.01d0) deh=-eh-0.01d0
          call popta1(rsm,eh+deh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
          call popta1(rsm,eh-deh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
          call popta2(l,eh,rsm,deh,e1,e2,e3,elim1,elim2,eaddx,enew,
     .       stife,jpr)

          radd = rnew-rsm
          eadd = enew-eh
          rsm = rnew
          eh = enew
          if (dabs(radd)<5d-3 .and. dabs(eadd)< 5d-3) goto 90
        enddo
   90   continue
C   ... End of iteration loop

        sume2 = sume2 + qvl*e2
        if (ipr > 20)
     .  write (stdo,260) l,nrep,rsm,eh,stifr,stife,e2,eval,pnu,qvl
  260   format(i2,i4,2f8.3,1x,2f9.1,1x,2f10.5,f8.2,f7.2)
  261   format(' l  it    Rsm      Eh     stiffR   stiffE',
     .     '      Eval      Exact     Pnu    Ql')
        istifr = stifr+0.5d0
        istife = stife+0.5d0
        write (stdl,710) l,nrep,rsm,eh,istifr,istife,e2,eval,pnu,qvl
  710   format('fa op',i2,i4,2f7.3,'  stf',2i6,'  ev',2f9.5,
     .     '  pq',2f6.2)

C   ... Possibly constrain rsm
        if (mod(isw,10) == 1 .and. rsm > rmt) then
        if (ipr > 20)
     .  write(stdo,'('' ... rsm exceeded rmt .. repeat with rsm=rmt'')')
        rsm = rmt
        sume2 = sume2 - qvl*e2

        do  irep = 1, 50
          nrep = irep
C     ... Get center energy
          call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e2,qrmt)
C     ... Vary eh
          deh = deh0
c         if (eh+deh>-0.01d0) deh=-eh-0.01d0
          call popta1(rsm,eh+deh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
          call popta1(rsm,eh-deh,l,z,rmt,nr,nrmt,rofi,h,
     .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
          call popta2(l,eh,rsm,deh,e1,e2,e3,elim1,elim2,eaddx,enew,
     .       stife,jpr)

          eadd = enew-eh
          eh = enew
          if (dabs(eadd)< 5d-3) exit
        enddo
C 190   continue
C   ... End of iteration loop

        sume2 = sume2 + qvl*e2
        if (ipr > 20)
     .  write (stdo,260) l,nrep,rsm,eh,stifr,stife,e2,eval,pnu,qvl
        istife = stife+0.5d0
        write (stdl,710) l,nrep,rsm,eh,istifr,istife,e2,eval,pnu,qvl

        elseif (mod(isw,10) == 2) then
          call rx('opfab not ready for this mode')

        endif

        itab(l+1,isp) = 1
        rtab(l+1,isp) = rsm
        etab(l+1,isp) = eh

      enddo
      if (ipr > 20) write (stdo,320) sume1,sume2,sume2-sume1
  320 format(' eigenvalue sum:  exact',f10.5,'    opt basis',f10.5,
     .   '    error',f8.5)
      write (stdl,720) sume1,sume2,sume2-sume1
  720 format('fa op sumev',f11.5,'   opt basis',f11.5,'   err',f9.5)

      enddo

C --- Make plot file ---
C     lplawv=nglob('lplawv')
      lplawv = 0
      if (cmdopt('--plotwf',8,0,strn)) lplawv = 1
      if (lplawv == 1) then
        if (nsp == 2) call rx('plotwf is not spinpol yet')
        allocate(psi(nr*lmxa))
        call popta5(lmxa,rtab,etab,itab,z,pl,rmax,rmt,nr,nrmt,
     .     rofi,psi,v,g,a,b,spid)
        deallocate(psi)
      endif

      deallocate(h,g,gp)

      end

      subroutine fabasp(autob,pnudef,lmxb,modeV,vbar,eh0,eh02,rsmn,ehmn,
     .  rsmx,ehmx,elocmx,pqrmx,esprd,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,pz,
     .  nsp,v,rofi,spid,itab,rtab,etab,pnuopt,pzopt,tolgv,gtop)
C- Choose a smooth-Hankel basis from the free-atom potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   autob :1s digit
Ci         :1 or 3 1-kappa parameters RSMH, EH
Ci         :2 or 4 2-kappa parameters RSMH, EH, RSMH2, EH2
Ci         :3 or 4: use any of predefined (i.e. input) RSMH, EH, RSMH2, EH2
Ci         :        See Input/Outputs
Ci         :10s digit autogenerate PZ
Ci         :0 do not modify PZ
Ci         :1 or 2 Set PZ satisfying criteria
Ci         :1 PZ specified by caller takes precedence
Ci         :2 PZ satisfying criteria is always autogenerated
Ci         :100s digit
Ci         :0 do not modify P
Ci         :1 Autogenerate P for l<=lmxb
Ci         :2 Autogenerate P for l<=lmxa
Ci  pnudef :mode that controls default lower bound of pnu
Ci         :0 Use version 6 defaults
Ci         :1 Use defaults tailored for LDA
Ci         :2 Use defaults tailored for GW
Ci   lmxb  :lmxb(1) l-cutoff for RSMH; lmxb(2) l-cutoff for RSMH2
Ci   modeV :Mode by which potential may be altered
Ci   vbar  :parameter in altering potential
Ci   rsmn  :smallest allowed rsm
Ci   ehmh  :smallest allowed eh
Ci   rsmx  :largest allowed rsm
Ci   ehmx  :largest allowed eh
Ci   elocmx:Set local orbital PZ set when E>elocmx
Ci   pqrmx :Set local orbital PZ set when q(r>rmt)>pqrmx
Ci   esprd :spread between high,low eh
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rmax  :Largest r for free atom
Ci   nrmt  :number of points between 0..rmt
Ci   rmt   :augmentation radius, in a.u.
Ci   lmxa  :muffin-tin l-cutoff
Ci   pl    :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   ql    :sphere moments
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   spid  :species label
Cio Inputs/Outputs
Cio  rtab  :smoothing radius for optimized wave function
Cio        :rtab(l+1,:) < 0   : value not defined; always set it
Cio        :rtab(l+1,:) >= 0  : if autob>2, keep value unchanged; else set it
Cio  etab  :energy for optimized wave function
Cio        :etab(l+1,:) = NULL: value not defined; always set it
Cio        :etab(l+1,:)!= NULL: if autob>2, keep value unchanged; else set it
Co Outputs
Co   itab  :itab(l+1,1) = 1: parameters RSMH,EH were chosen for this l
Co         :itab(l+1,2) = 1: parameters RSMH2,EH2 were chosen for this l
Co         :itab(l+1,:) =-1: parameters RSMH unchanged from input
Co   pnuopt:
Co   pzopt :
Cl Local variables
Cr Remarks
Ci   "Automatic basis" finder.  Each l is fit independently
Ci   Overwrites each of RSMH, EH, RSMH2, EH2, for each l depending on
Ci   autob and whether the parameter is "defined" (see Inputs/Outputs)
Ci     1. value is defined and 1s digit autob < 1 or 2
Ci        Restore input parameter after fitting algorithm described below
Ci     2. Otherwise, use algorithm below to find parameter.
Ci        Fit parameter overwrites input value
Ci   Fit algorithm proceeds as follows:
Ci      Fit is supplied, according to itab
Ci      If eh0 is not NULL, use eh0 for eh.
Ci      Otherwise:
Ci      If eh<ehmx and rsm<rsmx, fit parameter is used
Ci      Otherwise:
Ci      If eh>ehmx and rsm>rsmx, assign ehmx=eh, rsm=rsmx
Ci      Otherwise:
Ci      If eh>ehmx assign ehmx=eh, and refit rsm
Ci      Otherwise:
Ci      If rsm>rsmx, assign rsm=rsmx
Ci      Orbitals not fit are assigned ehmx=eh, rsm=rsmx
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   18 Dec 10 autogen pnu keep above cutoff as specified by pnudef
Cu   10 May 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmxa,nr,nrmt,nsp,n0,autob,lmxb(2),modeV,pnudef
      parameter (n0=10)
      integer itab(n0,2)
      double precision a,rmax,rmt,z,rsmn,rsmx,ehmn,ehmx,esprd,vbar,eh0,eh02,
     .  elocmx,pqrmx
      double precision rofi(*),v(nr,nsp),pl(n0,nsp),ql(3,n0,nsp),pz(n0)
      double precision pnuopt(n0),pzopt(n0),tolgv,gtop
      double precision rtab(n0,2),etab(n0,2)
      character spid*8
C ... Dynamically allocated arrays
      real(8), allocatable :: h(:),g(:),gp(:)
C ... Local parameters
      logical ltmp
      character c1*1,c2*1
      integer ipr,ir,irep,istife,istifr,jpr,konfig,l,
     .  nulli,lrel,nn,nrep,stdo,stdl,autop,autopz
      double precision b,deh,deh0,dphi,dphip,drsm,drsm0,e1,e2,e3,eadd,
     .  elo,eaddx,eh,eh2,elim1,elim2,enew,enu,eval,p,phi,phip,pnu,qvl,
     .  radd,raddx,rlim1,rlim2,rnew,rsm,rsm2,stife,stifr,sume1,sume2,
     .  q,qrmt,vtrue,pi,pdef(n0,2),pfree(n0),gam,gmax
      double precision rtab0(n0,2),etab0(n0,2)
      double precision veff(nr)
      parameter (NULLI=-99999)
      procedure(integer) :: nglob,iprint,lgunit

      call iinit(itab,n0*2)
      if (autob == 0) return
      call dcopy(size(rtab0),rtab,1,rtab0,1)
      call dcopy(size(etab0),etab,1,etab0,1)
      call dpzero(rtab, n0*2)
      call dpzero(etab, n0*2)
      ipr = iprint()
      stdo = lgunit(1)
      stdl = lgunit(2)
      lrel = mod(nglob('lrel'),10)
      if (z < 0.99d0) lrel = 0
      if (lmxa > 8) call rx('fabasp:  lmax too large')
      b = rmax/(dexp(a*nr-a)-1d0)
      if (rsmx == 0 .or. rsmx == NULLI) rsmx = rmt
      pi = 4d0*datan(1d0)
      gtop = 0
      autopz = mod(autob/10,10)
      autop  = mod(autob/100,10)

C ... Set P to defaults for l>lmxb
      if (autop /= 0) then
        call dpzero(pdef,n0*2)
        call defpq(10+pnudef,z,lmxa,1,pdef,0)
        call defpq(20+pnudef,z,lmxa,1,pfree,0)
        do  l = max(lmxb(1)+1,2), lmxa
          pnuopt(l+1) = pdef(l+1,1)
        enddo
      endif

C ... Effective potential which generates RSMH,EH
      if (nsp == 1) then
        call dpscop(v,veff,nr,1,1,1d0)
      else
        call dpscop(v,veff,nr,1,1,0.5d0)
        call dpsadd(veff,v(1,2),nr,1,1,0.5d0)
      endif
      if (vbar /= NULLI .and. modeV == 1) then
        do  ir = 2, nr
          vtrue = veff(ir) - 2*z/rofi(ir)
          if (vtrue > vbar) then
            veff(ir) = vbar + 2*z/rofi(ir)
          endif
        enddo
      elseif (vbar /= NULLI .and. modeV == 2) then
        do  ir = nrmt+1, nr
          veff(ir) = vbar + 2*z/rofi(ir)
        enddo
      endif

      call info5(20,0,0,'%N Make LMTO basis parms for species '//spid
     .  //'%a to lmxb=%i, rmt=%;5g  vbar=%;5g',lmxb,rmt,vbar,0,0)
      allocate(h(nr))
      allocate(g(2*nr))
      allocate(gp(2*nr*4))

C --- Setup for minimisation ---
      drsm0 = 0.1d0
      rlim1 = 0.3d0
      rlim2 = 2*rmt
      raddx = 0.2d0

      deh0  = 0.05d0
      elim1 = -5.0d0
      elim2 = -0.10d0
C     elim2 = -0.20d0
      eaddx = 0.099d0
      jpr=0
      if (ipr >= 50) jpr=1

C ... Default value for Pnu (autop == 1 -> default set by input)
      if (autop == 2) then
        do  l = 0, lmxa
          if (autop == 1 .and. pnuopt(l+1) /= 0) cycle
          pnuopt(l+1) = pdef(l+1,1)
        enddo
      endif

C --- Find MTO parameters RSMH,EH ---
      if (ipr >= 20) write (stdo,261)
      sume1 = 0d0
      sume2 = 0d0
C      lrtab = 0
C      if (dasum(lmxb(1)+1,rtab,1) /= 0) lrtab = 1 ! rsmh already available
      do  l = 0, lmxb(1)
        gmax = 0
        konfig = (pl(l+1,1)+pl(l+1,nsp))/2
        nn = konfig-l-1
        qvl = (ql(1,l+1,1)+ql(1,l+1,nsp))/2
C   ... Get exact fa wavefunction, eigval, pnu at rmt
        call popta3(0,l,z,nn,nr,nrmt,rofi,veff,a,b,eval,pnu,g)
        if (autop/=0) then
          pnuopt(l+1) = pnu
C         Shouldn't be too low
C         if (l >= 2) pnuopt(l+1) = max(pnu,pdef(l+1,1))
          if (pdef(l+1,1) < 0.5d0) pnuopt(l+1) = max(pnu,pdef(l+1,1))
C         Fractional part of pnuopt should be at least pfree
          pfree(l+1) = pfree(l+1) - int(pfree(l+1)) + int(pnuopt(l+1))
          pnuopt(l+1) = max(pnuopt(l+1),pfree(l+1))
        endif

        nrep = 1
        if (eval <= 0d0) then
        sume1 = sume1 + qvl*eval
C   ... Potential parameters at MT sphere
        call popta4(l,z,rmt,nrmt,rofi,veff,g,gp,a,b,pnu,enu,p,phi,dphi,phip,dphip)
        rsm = rsmx
        eh = -1; if (eh0 /= NULLI) eh = eh0
        if (jpr > 0) write (stdo,340)
  340   format('  L   parin    aux      E1       E2       E3',
     .     '       stiff    Eout     parout')
        do  irep = 1, 50
          nrep = irep
C     ... Expectation value of sm. Hankel at this rsm,eh
          call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .       veff,a,b,enu,p,phi,dphi,phip,dphip,e2,qrmt)
C     ... Vary rsm
          drsm = drsm0
          call popta1(rsm+drsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .       veff,a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
          call popta1(rsm-drsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .       veff,a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
          call popta2(l,rsm,eh,drsm,e1,e2,e3,rlim1,rlim2,raddx,rnew,
     .       stifr,jpr)

C     ... Vary eh
          deh = deh0
c         if (eh+deh>-0.01d0) deh=-eh-0.01d0
          call popta1(rsm,eh+deh,l,z,rmt,nr,nrmt,rofi,h,
     .       veff,a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
          call popta1(rsm,eh-deh,l,z,rmt,nr,nrmt,rofi,h,
     .       veff,a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
          call popta2(l,eh,rsm,deh,e1,e2,e3,elim1,elim2,eaddx,enew,
     .       stife,jpr)

          rnew = min(rnew,rsmx)
          radd = rnew-rsm
          rsm = rnew
          enew = min(ehmx,enew); if (eh0 /= NULLI) enew = eh0
          eadd = enew-eh
          eh = enew

          if (abs(radd)<5d-3 .and. abs(eadd)< 5d-3) exit
        enddo

C   ... Estimate maximum FT G for
        gam = rsm*rsm/4d0; gmax = 1d0
        do  irep = 1, 10
          gmax = sqrt(-log(tolgv/gmax**l)/gam)
        enddo
        sume2 = sume2 + qvl*e2
        c1 = ' '
        if (rsm == rsmx) c1 = '*'
        c2 = ' '
        if (eh == ehmx .or. eh == eh0) c2 = '*'
        if (ipr >= 20)
     .  write (stdo,260) l,nrep,rsm,c1,eh,c2,e2,eval,pnu,qvl,gmax
  260   format(i2,i4,f8.3,a1,f8.3,a1,1x,2f10.5,f8.2,f7.2:f6.1)
  261   format(' l  it    Rsm       Eh      ',
     .     '  Eval      Exact     Pnu    Ql   Gmax')
        istifr = stifr+0.5d0
        istife = stife+0.5d0
        write (stdl,710) l,nrep,rsm,eh,istifr,istife,e2,eval,pnu,qvl,gmax
  710   format('fa op',i2,i4,2f7.3,'  stf',2i6,'  ev',2f9.5,
     .     '  pq',2f6.2,f6.1)

        if (rsm < rsmn .or. eh < ehmn) then
          call info5(20,0,0,' ... l=%i  fit rsm=%,4;4d out of '//
     .      'range ... revert to default rsmx=%,4;4d',l,rsm,rsmx,4,5)
        else
          itab(l+1,1) = 1
        endif
        endif

C   ... No parameter generated; make default
        if (itab(l+1,1) == 0) then
          rsm = rsmx; gam = rsm*rsm/4d0; gmax = 1d0
          do  irep = 1, 10
            gmax = sqrt(-log(tolgv/gmax**l)/gam)
          enddo
          eh = ehmx; if (eh0 /= NULLI) eh = eh0

          if (eval < 0.0_8) then
            call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .        veff,a,b,enu,p,phi,dphi,phip,dphip,e2,qrmt)
            if (ipr >= 20)
     .        write (stdo,260) l,nrep,rsm,'+',eh,'+',e2,eval,pnu,qvl
          end if
        endif

        gtop = max(gtop,gmax)

C       if (lrtab == 0 .or. mod(autob,10) > 2) then
        itab(l+1,1) = 1
        rtab(l+1,1) = rsm
        etab(l+1,1) = eh
C       endif

      enddo

C --- Second kappa parameters ---
      if (mod(autob,10) == 2 .or. mod(autob,10) == 4) then
      do  l = 0, min(lmxb(1),lmxb(2))

        elo = ehmx - esprd ; if (eh0 /= NULLI) elo = eh0
        if (itab(l+1,1) == 1) then
          rsm  = rtab(l+1,1)
          rsm2 = rtab(l+1,1)
          eh   = etab(l+1,1)
          eh2  = eh - esprd
          if (eh0 == NULLI) then
          if (eh2 < elo) then
            eh2 = eh
            eh = eh2 + esprd
          endif
          if (eh > ehmx) then
            eh  = ehmx
            eh2 = ehmx - esprd
          endif
          endif
          if (eh02 /= NULLI) eh2 = eh02
          rtab(l+1,1) = rsm
          rtab(l+1,2) = rsm2
          etab(l+1,1) = eh
          etab(l+1,2) = eh2
          if (z == 0) rtab(l+1,2) = 0d0
        endif
      enddo
      endif

C --- Optionally restore input parameters ---
      if (mod(autob,10) <= 2) then
        do  ir = 1, 2
          do  l = 0, lmxb(ir)
            if (rtab0(l+1,ir) >= 0) then
              rtab(l+1,ir) = rtab0(l+1,ir)
            endif
            if (etab0(l+1,ir) /= NULLI) then
              etab(l+1,ir) = etab0(l+1,ir)
            endif
          enddo
          if (mod(autob,10) == 2 .or. mod(autob,10) == 4) exit
        enddo
      endif

C --- Print autogenerated Pnu ---
      if (autop /= 0) then
        call info2(20,1,0,' Autogenerated Pnu: %n:1,3;3d',lmxa+1,pnuopt)
      endif

C --- Autofind PZ ---
      call dcopy(n0,pz,1,pzopt,1)
      if (autopz /= 0) then
        call info2(20,0,0,'%N Find local orbitals which satisfy '//
     .    'E > %;3g Ry  or  q(r>rmt) > %;4g',elocmx,pqrmx)

      do  l = 0, lmxa
        konfig = (pl(l+1,1)+pl(l+1,nsp))/2
        nn = konfig-l-1
        pzopt(l+1) = 0
        if (mod(int(pz(l+1)),10) > int(pl(l+1,1))) then
          pzopt(l+1) = pz(l+1)
          call info2(30,0,0,' l=%i  high-lying local orbital '//
     .    'specified.  Use: PZ=%,3;3d',l,pzopt(l+1))
          cycle
        endif
        if (autopz == 1) pzopt(l+1) = pz(l+1)
        if (nn == 0) cycle ! Candidates must have a state below this one
        qvl = (ql(1,l+1,1)+ql(1,l+1,nsp))/2
        ltmp = qvl > 0 .and. pz(l+1) == 0 ! states w/ charge have deep cores
        ! Exception: s states in column I metals
        ltmp = ltmp .and. (l>0 .or. (z/=11 .and. z/=19 .and. z/=37 .and. z/=55 .and. z/=87))
        if (ltmp) cycle

C   ... FA wavefunction of core, eigval, pnu at rmt
        pnu = konfig-1 + 0.5d0
        nn = nn-1
        call popta3(0,l,z,nn,nr,nrmt,rofi,veff,a,b,eval,pnu,g)
        call gintsl(g,g,a,nr,rofi,q)
        call gintsl(g,g,a,nrmt,rofi,qrmt)
        if ((autopz == 2 .or. pz(l+1) <= 0) .and. eval>elocmx .or. q-qrmt>pqrmx) pzopt(l+1) = 10+pnu
        call info5(30,0,0,' l=%i  eval=%,3;3d  Q(r>rmt)=%,4;4G  PZ=%,3;3d  Use: PZ=%,3;3d',
     .    l,eval,q-qrmt,pnu,pzopt(l+1))
      enddo
      endif

      deallocate(h,g,gp)
      end

      subroutine ftfalo(icst,z,a,nr,rmax,nrmt,rmt,lmxa,pnu,pz,rs3,eh3,
     .  vmtz,nsp,v,rofi,spid,tolgv,gtop)
C- Fit value and slope of local orbitals to smoothed Hankel
C ----------------------------------------------------------------------
Ci Inputs
Ci   icst  :1 constrain rsm to be <= rmt
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rmax  :muffin-tin radius, in a.u.
Ci   nrmt  :number of points between 0..rmt
Ci   rmt   :muffin-tin radius, in a.u.
Ci   lmxa  :muffin-tin l-cutoff
Ci   pl    :boundary conditions for valence wavefunctions.
Ci   pz    :boundary conditions for local orbital. pz=0 -> no loc. orb.
Ci         :10s digit controls how local orbital included in hamiltonian
Ci         :10s digit nonzero -> smooth Hankel tail is attached.
Ci   rs3   :minimum allowed smoothing radius in attaching Hankel tails
Ci         :to local orbitals
Ci   eh3   :Hankel energy when attaching Hankel tails to high-lying
Ci         :local orbitals
Ci   vmtz  :parameter used in attaching Hankel tails to local orbitals
Ci         :It is used as a constant shift to Hankel energies for the
Ci         :fitting of local orbitals to Hankel tails. Thus vmtz
Ci         :is an estimate for the potential at the MT radius.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   spid  :species label
Ci   tolgv :tolerance in FT representation of wave functions : see gtop
Ci         :
Ci          functions are less than tol.
Co Outputs
Co   gtop  :Fourier repsn of local wave function should include G vectors up
Co         :to gtop
Cl Local variables
Cr Remarks
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   16 Jun 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmxa,nr,nrmt,nsp,n0,icst
      double precision a,rmax,rmt,z,rs3,eh3,vmtz
      parameter (n0=10)
      double precision rofi(*),v(nr,nsp),pz(n0,nsp),pnu(n0,nsp),tolgv,
     .  gtop
C ... Dynamically allocated arrays
      real(8), allocatable :: h(:),g(:),gp(:),psi(:)
C ... Local parameters
      logical cmdopt
      character spid*8, strn*80, flg(2)*1
      integer ipr,iprint,i,konfig,l,info,lgunit,nn,stdo,
     .  lplawv,loclo,nfit,isw,irep
      double precision b,dasum,dphi,dphip,e2,eh,eval,p,phi,phip,
     .  pnul,rsm,rsmin,rsmax,ekin
C     emin and emax are the maximum allowed ranges in Hankel energies
C     for the fitting of local orbitals to Hankel tails.
      double precision emin,emax,tphi,gmax,gam
C     For plotting wave functions
      integer itab(n0,2)
      double precision rtab(n0,2),etab(n0,2),pl(n0,nsp),qrmt
      data flg/'*',' '/

C     return

      ipr = iprint()
      stdo = lgunit(1)
C     stdl = lgunit(2)
      if (lmxa > 8) call rx('ftfalo:  lmax too large')
      b = rmax/(dexp(a*nr-a)-1d0)
      nfit = 0
      gtop = 0d0
      if (dasum(lmxa+1,pz,1) == 0) return
      allocate(h(nr),g(2*nr),gp(2*nr*4))

      do  i = 1, nsp

C --- Loop over local orbitals ---
C      sume1 = 0d0
C      sume2 = 0d0
      do  l = 0, lmxa

        itab(l+1,i) = 0
        pnul = pnu(l+1,i)
        pl(l+1,i) = pnu(l+1,i)
        konfig = mod(pz(l+1,1),10d0)

C       Skip all but local orbitals with tails attached
        if (mod(pz(l+1,1),100d0) < 10) cycle

C       Case local orbital deeper than valence
        if (int(pnul-1) == int(mod(pz(l+1,1),10d0))) then
          loclo = 1
C         Not needed, actually, since overwritten by popta3
C         pnul = mod(pz(l+1,1),10d0)

C       Case local orbital higher than the valence state
        elseif (int(pnul+1) == int(mod(pz(l+1,1),10d0))) then
          pnul = mod(pz(l+1,1),10d0)
          loclo = 0

C       Local orbital neither one: error
        else
          call fexit3(-1,111,' Exit -1 freeat, l=%i:  sc '//
     .      'PZ=%d incompatible with valence P=%;3d',l,pz(l+1,1),pnul)
        endif

C       Skip high-lying local orbitals unless specifically sought
        if (loclo == 0 .and. .not. cmdopt('--getallloc',11,0,strn))
     .    cycle

        nfit = nfit + 1
        if (nfit == 1) then
          call info2(20,1,0,
     .    ' Fit local orbitals to sm hankels, species '//spid//
     .    '%a, rmt=%;7g',rmt,0)
          if (ipr >= 20) write (stdo,261)
        endif

C   ... Get exact fa wavefunction, eigval, pnu_l at rmt
        if (loclo == 1) then
          nn = konfig-l-1
          call popta3(0,l,z,nn,nr,nrmt,rofi,v(1,i),a,b,eval,
     .      pnul,g)
C       Finish if in future, need w.f. at r>rmt
C        else
C          call popta3(1,l,z,nn,rmt,nr,nrmt,rofi,v(1,i),a,b,eval,
C     .      pnul,g)
        endif
        pl(l+1,i) = pnul

C   ... Potential parameters at MT sphere
        call popta4(l,z,rmt,nrmt,rofi,v(1,i),g,gp,
     .    a,b,pnul,eval,p,phi,dphi,phip,dphip)

C   ... Set conditions on envelope functions ... For now
        rsmin = rs3
        rsmax = 5
        if (icst == 1) rsmax = rmt
C       Use r->infty value for energy
        eh = min(-.02d0,eval)

C   ... Match Hankel to phi,dphi
C        rsm = rsmin
C        emax = -.02d0
C        emin = -5d0
C        call mtchre(100,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,phi,
C     .    dphi,rsm,eh,ekin,info)

C   ... Match slope and K.E. of Hankel to phi,dphi
        tphi = eval - (v(nrmt,i)-2*z/rmt)
        rsm = 0
        eh = min(eval-vmtz,-.02d0)
        emax = -.02d0
        emin = -10d0
C       if (ipr >= 20) call pshpr(max(ipr,50))
        call mtchre(003,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,tphi,
     .    dphi,rsm,eh,ekin,info)
C       if (ipr >= 20) call poppr
C       Match failed ... turn up verbosity and repeat for info
        gmax = 0
        if (info == -1) then
          call info2(0,2,1,
     .      ' *** ftfalo (fatal) cannot fit smooth Hankel to w.f.'//
     .      ' class '//spid//
     .      '%N ... possibly reduce RS3 (current value = %,1d)',rs3,0)
          call pshpr(max(ipr,110))
          call mtchr2(1,l,emin,emax,(emin+emax)/2,
     .      rmt,phi,dphi,rsmin,eh,ekin,i)
          call poppr
C         call pshpr(max(ipr,110))
C         call mtchre(103,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,tphi,
C    .      dphi,rsm,eh,ekin,info)
          call fexit2(-1,111,
     .      ' Exit -1 : ftfalo : failed to match log der=%,1;3d'//
     .      ' to envelope, l=%i',dphi/phi,l)
C   ... get GMAX Fourier cutoff for this wave function
        else
          gam = rsm*rsm/4d0
          gmax = 1d0
          do  irep = 1, 10
            gmax = dsqrt(-dlog(tolgv/gmax**l)/gam)
C           write(stdo,895) irep,gmax
C 895       format('irep,gmax=',i5,f12.6)
            gtop = max(gtop,gmax)
          enddo
        endif

C  ... Get energy of this wave function
        call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,
     .    v(1,i),a,b,eval,p,phi,dphi,phip,dphip,e2,qrmt)

        if (ipr >= 20) then
          if (gmax > 0) then
            write (stdo,260) l,rsm,eh,qrmt,e2,eval,pnul,tphi,ekin,
     .        flg(2-isw(dabs(ekin-tphi)>1d-5)),gmax
          else
            write (stdo,260) l,rsm,eh,qrmt,e2,eval,pnul,tphi,ekin,
     .        flg(2-isw(dabs(ekin-tphi)>1d-5)),gmax
          endif
        endif

  260   format(i2,2f7.3,3f10.5,f8.3,2f9.4,a1:f6.1)
  261   format(' l   Rsm    Eh     Q(r>rmt)   Eval',
     .    '      Exact     Pnu     K.E.   fit K.E.  Gmax')
        itab(l+1,i) = 1
        rtab(l+1,i) = rsm
        etab(l+1,i) = eh

C  10 continue
      enddo
C  80 continue
      enddo

C --- Make plot file ---
C     lplawv=nglob('lplawv')
      lplawv = 0
      if (cmdopt('--plotwf',8,0,strn)) lplawv = 1
      if (lplawv == 1) then
        if (nsp == 2) call rx('optfab is not spinpol yet')
        allocate(psi(nr*lmxa*nsp))
        call popta5(lmxa,rtab,etab,itab,z,pl,rmax,rmt,nr,nrmt,
     .     rofi,psi,v,g,a,b,spid)
        deallocate(psi)
      endif

      deallocate(h,g,gp)

      end

      subroutine popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,v,a,b,enu,p,
     .  phi,dphi,phip,dphip,eval,qrmt)
C- Calculate expectation value for smooth Hankel
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsm   :smoothing radius of basis function
Ci   eh    :energy of basis function
Ci   l     :l quantum number
Ci   z     :nuclear charge
Ci   rmt   :muffin-tin radius
Ci   nr    :number of radial mesh points
Ci   nrmt  :number points to muffin-tin radius
Ci   rofi  :radial mesh points
Ci   h     :work array
Ci   v     :spherical potential (atomsr.f)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   enu   :enu's for making charge density
Ci   p     :<gp**2> (potential parameter)
Co   phi   :wave function at rmt
Co   dphi  :radial derivative of of phi at rmt
Co   phip  :energy derivative of phi
Co   dphip :radial derivative of dphi
Co Outputs
Co   eval  :expectation value
Co   qrmt  :fraction of (wave function)^2 for r>rmt
Cr Remarks
Cu Updates
Cu   24 Sep 04 return qrmt
Cu   16 Jun 04 Adapted to new hansmd, mtchae
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer l,nr,nrmt
      double precision a,b,dphi,dphip,eh,enu,eval,p,phi,phip,rmt,rsm,z,
     .  rofi(nr),h(nr),v(nr),qrmt
C ... Local parameters
      integer i
      double precision alfa,beta,det,drdi,hlap,
     .  hum,hum1,hum2,r,sum,sum1,sum2,tum2,vum2,wt
C     double precision xi(0:20)

      double precision hs(0:l),dhs(0:l),ddhs(0:l)

C     pi = 4d0*datan(1d0)
C     asm = 1d0/rsm
C     lp1 = l+1

C ... Integrals over smooth Hankel on mesh
C     gfac = (asm*asm/pi)**1.5d0 * dexp(eh*rsm*rsm/4d0)
C     ta2 = 2d0*asm*asm
      tum2 = 0d0
      sum2 = 0d0
      vum2 = 0d0

      do  i = nrmt, nr
        r = rofi(i)

C   ... Make r*h and r Laplacian h, including L^2
        call hansmd(2,r,eh,rsm,l,hs,dhs,ddhs,det,det,det)
        h(i) = hs(l)*r
        hlap = ddhs(l)*r
CC      Old : r*h and r Laplacian h, including L^2
CC      h = r*radial part of sm. Hankel
C       call hansmr(r,eh,asm,xi,l)
C       h(i) = xi(l)*(r**lp1)
CC      radial part of Gaussian
C       gl = gfac * dexp(-asm*asm*r*r) * ta2**l * (r**lp1)
CC      r * (nabla_r - l(l+1)/r^2) h_l
C       hlap = -4d0*pi*gl - eh*h(i)

C  ...  Accumulate <h h>, <h v h>, <h -nabla h>
        wt = 2*(mod(i+1,2)+1)/3d0
        if (i==nrmt .or. i==nr) wt = 1d0/3d0
        drdi = a*(r+b)
        sum2 = sum2 + wt*drdi*h(i)*h(i)
        vum2 = vum2 + wt*drdi*h(i)*h(i)*(v(i)-2d0*z/r)
        tum2 = tum2 + wt*drdi*h(i)*(-hlap)
      enddo
      hum2 = tum2+vum2

C --- BC's: match phi,phidot to envelope at RMT ---
      call mtchae(0,rsm,eh,l,rmt,phi,dphi,phip,dphip,alfa,beta)
CC    OLD matching
CC    Match value, slope fl,dfl to linear combination of phi,phidot
C     call hansmr(rmt,eh,asm,xi,l+1)
CC    Value and radial derivative of h (JMP 39, 3393, Eq. 4.7)
C     fl = xi(l)*rmt**l
C     flp1 = xi(l+1)*rmt**(l+1)
C     dfl = l*fl/rmt-flp1
CC    Match fl,dfl to linear combination of phi,phidot
CC    Use  phi=phi(R); phip=phidot(R) dphi=phi'(R); dphip=phidot'(R)
CC    (phi  phip ) (alpha)   (fl )    (alpha)    1  (dphip -phip) (fl )
CC    (          ) (     ) = (   ) -> (     ) = --- (           ) (   )
CC    (dphi dphip) (beta )   (dfl)    (beta )   det (-dphi  phi ) (dfl)
C     det = phi*dphip-dphi*phip
C     alfa = (fl*dphip-dfl*phip)/det
C     beta = (dfl*phi-fl*dphi)/det

C     O = alpha^2 <phi | phi> + beta^2 <phidot | phidot>
      sum1 = alfa*alfa + beta*beta*p
      hum1 = alfa*alfa*enu + alfa*beta + beta*beta*enu*p

      sum = sum1+sum2
      hum = hum1+hum2
      eval = hum/sum

      qrmt = sum2/sum
      end

      subroutine popta2(l,x0,y0,dx,e1,e2,e3,xmin,xmax,xshx,xnew,stiff,jpr)
C- Find minimum from three values
C ----------------------------------------------------------------------
Ci Inputs
Ci   l     :angular momentum
Ci   x0    :starting value
Ci   y0    :used for printout
Ci   dx    :excursion in x for numerical differentiation
Ci   e1    :function value at x0-dx
Ci   e2    :function value at x0
Ci   e3    :function value at x0+dx
Ci   xmin  :boundary: estimated minimum must be >= xmin
Ci   xmax  :boundary: estimated minimum must be <= xmax
Ci   xshx  :maximum step size
Ci   jpr   :printout verbosity
Co Outputs
Co   xnew  :new estimate for the minimum
Co   stiff :estimated curvature
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer jpr,l
      double precision dx,e1,e2,e3,stiff,x0,xmax,xmin,xnew,xshx,y0
C ... Local parameters
      integer ie0,lgunit,stdo
      double precision a,aa,b,c,ee1,ee2,ee3,een,enew,xadd

      stdo = lgunit(1)
      c = e2
      b = (e3-e1)/(2*dx)
      a = (e1+e3-2*e2)/(2*dx*dx)
      if (a <= 0d0) then
        xadd = -xshx
        enew = e1
        if (e3 < e1) xadd = xshx
        if (e3 < e1) enew = e3
      else
        xadd = -b/(2*a)
        enew = a*xadd*xadd + b*xadd + c
      endif
      aa = 2*1d3*a

      if (xadd > xshx)  xadd = xshx
      if (xadd < -xshx) xadd = -xshx
      xnew = x0+xadd
      if (xnew > xmax) xnew = xmax
      if (xnew < xmin) xnew = xmin

      ie0 = e2
      ee1 = 1d3*(e1-ie0)
      ee2 = 1d3*(e2-ie0)
      ee3 = 1d3*(e3-ie0)
      een = 1d3*(enew-ie0)
      stiff = aa

      if (jpr>0) write (stdo,810)l,x0,y0,ee1,ee2,ee3,aa,een,xnew
  810 format(i3,f8.3,f8.3,f10.3,2f9.3,f9.1,f10.3,f8.3,a)

      end
      subroutine popta3(mode,l,z,nn,nr,nrmt,rofi,v,a,b,evl,pnu,g)
C- Get exact fa wavefunction, eigval, pnu at Rmt.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 boundary condition is val,slo = 0 at nr
Ci         :1 boundary condition is that w.f. satisfy pnu at nrmt
Ci         :  (under development)
Ci   l     :angular momentum
Ci   z     :nuclear charge
Ci   nr    :number of radial mesh points
Ci   nrmt  :number of radial mesh points to rmt
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio Inputs/Outputs
Cio  nn    :number of nodes (input mode 0; output mode 1)
Cio  pnu   :boundary condition at rmt (output mode 0; input mode 1)
Co Outputs
Co   g     :normalized wave function times r
Co   evl   :eigenvalue
Cl Local variables
Cl   rmt   :muffin-tin radius, in a.u.
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,l,nn,nr,nrmt
      double precision a,b,evl,pnu,z,rofi(nr),v(nr),g(nr*2)
C ... Local parameters
      integer lgunit,nre,stdo,konfig,nri,nn2,iprint,kc
      double precision d0l,p0l,dphi,drdi,du,eb1,eb2,g1,g2,g3,g4,g5,pi,
     .  slo,slou,sum,tol,val,valu,dnu,rmt

      stdo = lgunit(1)
      pi = 4d0*datan(1d0)

      eb1 = -30
      eb2 = 20
      tol = 1d-10
      val = 1d-30
      slo = -val
      evl = -0.5d0
      nri = nr
      rmt = rofi(nrmt)
      if (mode == 1) then
        konfig = pnu
        nn = konfig-l-1
        dnu = dtan(pi*(0.5d0-pnu))
        val = rmt
        slo = dnu+1d0
        nri = nrmt
      endif
      call rseq(eb1,eb2,evl,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nri,nre,
     .  kc)
      if (mode == 1) then
C       integration becomes rather strange for r>>rmt.
C       Need to truncate radius.
        call rsq1(nri,evl,l,z,v,nr,g,val,slo,nn2,a,b,rofi,nr)
C       call prrmsh('g',rofi,g,nr,nr,1)
        call rx('not finished mode 1')
      endif
      g1 = g(nrmt-2)
      g2 = g(nrmt-1)
      g3 = g(nrmt)
      g4 = g(nrmt+1)
      g5 = g(nrmt+2)
      drdi = a*(rmt+b)
      valu = g3
      slou = (-2*g5+16*g4-16*g2+2*g1)/(24d0*drdi)
      du   = rmt*slou/valu
      dphi = du-1
      pnu  = nn+l+1 + (0.5d0-datan(dphi)/pi)

C ... Don't set too low..
      d0l = l
      p0l = nn+l+1 + 0.5d0-datan(d0l)/pi
      p0l = nn+l+1 + 0.1d0
      if (l == 0)  p0l = nn+l+1 + 0.3d0
      if (l == 2)  p0l = nn+l+1 + 0.2d0
      if (pnu < p0l) then
        if (iprint() > 40) write (stdo,145) l,pnu,p0l
  145   format(' ... l=',i1,'  increase Pnu=',f8.3,'  to ',f8.3)
        pnu = p0l
      endif

      end
      subroutine popta4(l,z,rmt,nrmt,rofi,v,g,gp,a,b,pnu,enu,p,phi,dphi,
     .  phip,dphip)
C- Potential parameters at MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   l     :angular momentum
Ci   z     :nuclear charge
Ci   rmt   :muffin-tin radius, in a.u.
Ci   nrmt  :number of radial mesh points to rmt
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   g     :normalized wave function times r
Ci   gp    :energy derivative(s) of g
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   enu   :enu's for making charge density
Co Outputs
Co   phi   :wave function at rmt
Co   dphi  :radial derivative of of phi at rmt
Co   phip  :energy derivative of phi
Co   dphip :radial derivative of dphi
Co   p     :<gp**2> (potential parameter)
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer l,nrmt
      double precision a,b,dphi,dphip,enu,p,phi,phip,pnu,rmt,z
      double precision rofi(nrmt),v(nrmt),g(nrmt),gp(nrmt,4)
C ... Local parameters
      integer konfig,nn,nre,kc,lwronsk
      double precision dnu,eb1,eb2,pi,slo(5),sum,tol,val(5)
      procedure(integer) :: nglob

      pi = 4d0*datan(1d0)
      eb1 = -30
      eb2 = 20
      tol = 1d-10

      lwronsk = nglob('wronsk')
      konfig = pnu
      nn = konfig-l-1
      dnu = dtan(pi*(0.5d0-pnu))
      val(1) = rmt
      slo(1) = dnu+1d0
      enu=-0.5d0

      call rseq(eb1,eb2,enu,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nrmt,nre,kc)
      val(1) = val(1)/dsqrt(sum)
      slo(1) = slo(1)/dsqrt(sum)

C      call phidot(z,l,v,enu,a,b,rofi,nrmt,g,val,slo,tol,nn,gp,phi,dphi,
C     .  phip,dphip,p)

      call phidx(10*lwronsk+1,z,l,v,0d0,0d0,rofi,nrmt,2,tol,enu,val,slo,nn,g,gp,
     .  phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
C     dphip = (slo(2)-phip)/rmt


c|    write(stdo,200) l,enu,p,phi,dphi,phip,dphip
c|200 format(' PP',i2,'  e',f10.5,'  p',f10.5,'  bc',4f10.5)

      end

      subroutine popta5(lmax,rtab,etab,itab,z,pl,rmax,rmt,nr,nrmt,
     .  rofi,psi,v,g,a,b,spid)
C- Write wave functions to plot file
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmax  :maximum l for a given site
Ci   rtab  :smoothing radii for wavefunction, each l
Ci   etab  :smoothed hankel energies for wavefunction, each l
Ci   itab  :1 if a wave function calculated, 0 if not
Ci   z     :nuclear charge
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   rmax  :muffin-tin radius, in a.u.
Ci   rmt   :muffin-tin radius
Ci   nr    :number of radial mesh points
Ci   nrmt  :number points to muffin-tin radius
Ci   rofi  :radial mesh points
Ci   psi   :wave function tabulated on the rofi mesh
Ci   v     :spherical potential (atomsr.f)
Ci   g     :normalized wave function times r
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   spid
Co Outputs
C    wave functions written to disk
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer itab(0:*),lmax,nr,nrmt,n0
      parameter (n0=10)
      double precision a,b,rmax,rmt,z,rtab(0:*),etab(0:*),
     .  rofi(nr),psi(nr,0:*),g(nr,2),v(nr),pl(0:n0-1)
      character spid*8
C ... Local parameters
      integer i,ifi,konfig,l,lgunit,lp1,m,n,nn,nre,fopna,stdo,kc
      integer ltab(n0)
      double precision asm,dfl,drdi,eb1,eb2,eh,evl,fac,fl,flp1,r,rsm,
     .  slo,sum1,sum2,tol,val,wt,xi(0:20)
      character str*32

      stdo = lgunit(1)
      eb1 = -20
      eb2 = 20
      tol = 1d-8
      n = 0

      do  l = 0, lmax
        if (itab(l) == 0) cycle
        n = n+1
        ltab(n) = l
        lp1 = l+1
        rsm = rtab(l)
        eh = etab(l)
        asm = 1d0/rsm
        konfig = pl(l)
        nn = konfig-l-1

C ...   Smooth hankel fct outside rmt
        sum2 = 0d0
        do  i = nrmt, nr
          r = rofi(i)
          call hansmr(r,eh,asm,xi,l)
          psi(i,n) = xi(l)*(r**lp1)
          wt = 2*(mod(i+1,2)+1)/3d0
          if (i==nrmt .or. i==nr) wt=1d0/3d0
          drdi = a*(r+b)
          sum2 = sum2 + wt*drdi*psi(i,n)**2
        enddo

C ...   Attach numerical solution inside MT sphere
        call hansmr(rmt,eh,asm,xi,l+1)
        fl   = xi(l)*rmt**l
        flp1 = xi(l+1)*rmt**(l+1)
        dfl  = l*fl/rmt-flp1
        val = rmt*fl
        slo = rmt*dfl+fl
        evl = -0.5d0
        call rseq(eb1,eb2,evl,tol,z,l,nn,val,slo,v,
     .     g,sum1,a,b,rofi,nrmt,nre,kc)
        fac = val/(g(nrmt,1)*dsqrt(sum1+sum2))
        do  i = 1, nrmt
          psi(i,n) = fac*g(i,1)
        enddo
        fac = 1d0/dsqrt(sum1+sum2)
        do  i = nrmt+1,nr
          psi(i,n) = psi(i,n)*fac
        enddo
      enddo

C ... Write the plot file
      write (str,'(''wfa_'',a)') spid
      write (stdo,344) str
  344 format(/' Write fit wavefunctions to plot file: ',a)
      ifi = fopna(str,-1,0)
      write (ifi,490) spid,rmax,rmt,(ltab(i),i=1,n)
  490 format('# Free-atom opt basis (divided by r) for species ',
     .   a/'# rmax=',f7.3,'   rmt=',f7.3,
     .   '   l=',8i2)
        write (ifi,'(''% rows '',i5,'' cols '',i3)') nr,n+1
      do  i = 1, nr
        write (ifi,495) rofi(i),(psi(i,m),m=1,n)
  495   format(f9.5,1p,8d14.5)
      enddo
      call fclose(ifi)

      end

      subroutine fctail(nr,nrmt,a,b,rsm,rofi,rhoc,c,eh)
C- Fit one Hankel to tail of core density.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nrmt  :number points to muffin-tin radius
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rsm   :smoothing radius
Ci   rofi  :radial mesh points
Ci   rhoc  :core density
Co Outputs
Co   c     :coefficient to fit of rhoc(spin+)+rhoc(spin-)
Co   eh    :energy
Cl Local variables
Cl   rmt   :muffin-tin radius
Cr Remarks
Cb Bugs
Cb   Should this be fit to smoothed function??
Cu Updates
Cu   19 Apr 02 Make rmt a local variable.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nrmt
      double precision a,b,c,eh,rmt,rsm,rofi(nr),rhoc(nr,2)
C ... Local parameters
      integer i,nsp,lgunit,stdo,ipr,nglob
      double precision ak1,akap,fit,q,q0,r,s,v0,wt
      character sout*80

      call getpr(ipr)
      stdo = lgunit(1)
      nsp =  nglob('nsp')
      rmt = rofi(nrmt)
      q0 = 0d0
      do  i = nrmt, nr
        r = rofi(i)
        wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
        if (i==nrmt .or. i==nr) wt = a*(r+b)/3d0
        q0 = q0 + wt*(rhoc(i,1)+rhoc(i,nsp))/(3-nsp)
      enddo
      v0 = (rhoc(nrmt,1)+rhoc(nrmt,nsp))/(3-nsp)/(rmt*rmt)
      sout = ' '
      call awrit3('%?#(n>=30)#%N## coretail: q=%;3g, rho(rmt)=%;3g.',
     .  sout,len(sout),0,ipr,v0,q0)
C      if (ipr >= 20) write (stdo,339) v0,q0
C  339 format(/' coretail:  rho(rmt)=',f12.8,'   charge=',f12.8)
      if (dabs(q0) < 1d-6) then
        c = 0d0
        eh = -1d0
        return
      endif

C ... Parameters of hankel fct
      s = dsqrt(rmt**4 * v0**2 + 4*rmt*q0*v0)
      ak1 = (rmt*rmt*v0+s)/(2d0*q0)
C     ak2 = (rmt*rmt*v0-s)/(2d0*q0)
c|      write(stdo,975) ak1,ak2
c|  975 format('ak1,ak2=',2f14.8)
      akap = ak1
      c = rmt*v0*dexp(akap*rmt)
      eh = -akap*akap

      if (ipr >= 20) then
        call awrit2('%a  Fit with Hankel e=%;5g  coeff=%;5g',sout,
     .    len(sout),-stdo,eh,c)
      endif

C ... Test
      if (ipr > 30) then
      write (stdo,501)
      q = 0d0
      do  i = nrmt, nr
        r = rofi(i)
        wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
        if (i==nrmt .or. i==nr) wt = a*(r+b)/3d0
        fit = c*dexp(-akap*r)*r
        q = q+wt*fit
        if ((rhoc(i,1)+rhoc(i,nsp))/(3-nsp) < 1d-8) goto 90
        if (mod(i,5)==0 .or. i==nrmt)
     .     write (stdo,500) r,(rhoc(i,1)+rhoc(i,nsp))/(3-nsp),fit
  500   format(f12.6,2f14.8)
  501   format(6x,'r',12x,'rhoc',10x,'fit')
      enddo
   90 continue
c|      v=c*dexp(-akap*rmt)/rmt
c|      write(stdo,885) q,q0,v,v0
c|  885 format('q,q0,v,v0=',4f14.8)
      endif

c ... look at smoothed core..
c|      rg=0.4
c|      qc=36
c|      sum0=-c*dexp(eh*rg*rg/4d0)/eh
c|      cg=qc-sum0
c|      write(stdo,888) qc,sum0,cg
c|  888 format(' qcore=',f10.4,'  sum0=',f12.6,'   cg=',f12.6)
c|      ag=1d0/rg
c|      fac=4d0*pi*(ag*ag/pi)**1.5d0
c|      q=0d0
c|      do i=1,nr
c|        r=rofi(i)
c|        wt=2*(mod(i+1,2)+1)*a*(r+b)/3d0
c|        if (i==1 .or. i==nr) wt=a*(r+b)/3d0
c|        call hansmr(r,eh,ag,xi,1)
c|        fit=c*xi(0)*r*r + cg*fac*dexp(-ag*ag*r*r)*r*r
c|        if (rhoc(i)>1d-10) write(49,490) r,rhoc(i),fit
c|  490   format(f12.6,2f16.8)
c|        q=q+wt*fit
c|      enddo
c|      write(stdo,965) q
c|  965 format(' integral over smoothed core:',f10.5)

      end

      subroutine swbasp(sopts,autob,eh1,eh2,rsmmx,incrl)
C- Parse switches to control or modify sm Hankel basis
      implicit none
      character sopts*(*)
      double precision eh1,eh2,rsmmx
      integer autob,incrl
C ... Local variables
      integer j1,j2,m,i,iv(10)
      character dc*1
      integer, parameter :: NULLI=-99999
      procedure(integer) :: parg

      eh1 = NULLI; eh2 = NULLI; incrl = 0

      dc = sopts(1:1)
      if (dc /= ' ') then
C   ... Return here to resume parsing for arguments
        j2 = 0
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) == dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then
          elseif (sopts(j1:j1+3) == 'ctrl')  then
            if (mod(autob,10) > 2) autob = autob - 2
          elseif (sopts(j1:j1+5) == 'rsmmx=')  then
            m = 0
            i = parg('rsmmx=',4,sopts(j1:),m,len(sopts(j1:)),', '//dc,1,1,iv,rsmmx)
            if (i < 1) goto 999
          elseif (sopts(j1:j1+2) == 'eh=')  then
            m = 0
            i = parg('eh=',4,sopts(j1:),m,len(sopts(j1:)),', '//dc,1,1,iv,eh1)
            if (i < 1) goto 999
          elseif (sopts(j1:j1+3) == 'eh2=')  then
            m = 0
            i = parg('eh2=',4,sopts(j1:),m,len(sopts(j1:)),', '//dc,1,1,iv,eh2)
            if (i < 1) goto 999
          elseif (sopts(j1:j1+6) == 'incrlmx')  then
            incrl = 1
          else
            goto 999
          endif
          goto 10
        endif
      endif
      return

  999 continue
      call rxs('swbasp: failed to parse options:',sopts)

      end
