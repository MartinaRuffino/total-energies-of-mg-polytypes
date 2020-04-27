      subroutine suctrl(prgnam)
C- Creates ctrl file from init file
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of calling program
Co Outputs
Co   an input file is created, named `actrl'
Cs Command-line switches
Cs   --addes[~sw] : Add tags to prepare for adding empty spheres, especially through
Cs                : lmchk --findes --wsitex
Cs                : Switches:
Co                : ~end          Tells blm to put the E species last in SPEC
Cs   --asa[~sw]   : Tailor input file for ASA calculation
Cs                : Switches:
Cs                : ~rmaxs=#      Range for structure constants
Cs                : ~pfree        Add HAM_PMIN=-1 : restricts pnu to stay above free electron limit
Cs   --dv=nam=expr ...
Cs                : Add declarations to a %const directive (express must be <= 4)
Cs                : Note: no check is made on the syntax!
Cs   --ehmx=#     : Upper bound to EH (sets HAM_AUTOBAS_EHMX, used by lmfa)
Cs   --express[=n]: Use express mode level n (default n is 6); if switch is absent, default n is 3
Cs                :  n  mode
Cs                :  0  standard  All input through standard categories
Cs                :               No supporting comments are given.
Cs                :               Writes site file with --wsite or --wsitex
Cs                :
Cs                :           ... For n>0, an EXPRESS category is created
Cs                :               and a site file is created.  Thus input
Cs                :               conditions are set up so lattice and site
Cs                :               information are read through the site file
Cs                :               As the express level increase, the
Cs                :               input files become simpler, but contain
Cs                :               less information.
Cs                :
Cs                :  1  Expert    Like mode 0, but EXPRESS category is added.
Cs                :               Tags duplicated by EXPRESS are retained to
Cs                :               facilitate editing by the user.
Cs                :               For duplicated tags, EXPRESS takes precedence.
Cs                :               Input is terse with no supporting comments
Cs                :
Cs                :  2  Verbose   Similar to mode 1, with comments
Cs                :  3  Large     Similar to mode 2, duplicate tags removed
Cs                :  4  Standard  Most tags covered by defaults are removed.
Cs                :  5  Simple    No variables or expressions are used
Cs                :  6  Light     Some nonessential tags are removed
Cs                :  7  Skeleton  Minimal input file
Cs
Cs   --findes[~sw]: Find and add sites for empty spheres
Cs                : switches:
Cs                : ~rmin=#       set cutoff for smallest sphere size
Cs                : ~float        Set orbitals to float
Cs                : ~nesmx=#      Limit the number of E sites to #
Cs                : ~nspmx=#      Limit the number of E species to #
Cs                : ~1spec        Limit the number of E species to 1
Cs                : ~shorten=#,#,# shorten positions.  See shorps.
Cs                : ~sclwsr=#     Use in conjunction with ~omax=#.  Sets SPEC_SCLWSR when
Cs                :               omax is restored after ES are found and spheres are resized
Cs                : ~omax=#       Temporarily sets maximum overlap of actual atoms
Cs                :               when finding empty spheres.  Setting the overlap small
Cs                : ~mino         Refine initial placements, minimizinng sphere overlaps
Cs   --fixlat     : Adjust lattice vectors and group operations to render
Cs                : them internally consistent
Cs   --fixpos     : Make small adjustments to given positions if
Cs                : they are compatiable with, or increase symmetry
Cs   --nfile      : Lines for conditional read of site files siten in express category
Cs   --fpandasa   : Tags for both ASA and FP
Cs   --gf         : Add lines for lmgf --- implies ASA
Cs   --gmax=#     : If gmax is known in advance, assign # to GMAX=
Cs   --gwemax=#   : Same as gw~emax=#
Cs   --gw[~sw]    : Tailor input file for GW calculations (implies fp)
Cs                : Options:
Cs                :   ~pbtol=strn     Assign GW_PBTOL to strn, overriding default
Cs                :   ~emax=#         Set GW energy cutoff SIGP_EMAX, overriding default
Cs                :   ~gcutb=#        Set G-vector cutoff for interstitial part of one-particle objects
Cs                :   ~gcutx=#        Set G-vector cutoff for interstitial part of two-particle objects
Cs                :   ~emax=#         Set GW energy cutoff SIGP_EMAX, overriding default
Cs                :   ~nk=#,[#2][,#3] specifies k mesh for GW_NKABC
Cs   --input      : show tokens program will attempt to read, without reading anything
Cs   --loc=#      : Set AUTOBAS_LOC
Cu   --eloc=#     : Set AUTOBAS_ELOC
Cs   --mto=#      : Set AUTOBAS_MTO
Cs   --mag        : Modify input file for a magnetic calculation
Cs   --molstat    : Add tag for molecular statics (lattice relaxation) (lmf)
Cs   --lmfit      : Add tag for Levenberg-Marquardt fitting (lm)
Cs   --nk=        : --nk=#1[,#2,#3] specifies k mesh
Cs   --nkgw=      : Same as gw~nk=...
Cs   --nl=        : set maximum (global l)+1 ... becomes lmxa+1 if gw is also set
Cs   --noshorten | --shorten=no
Cs                : write site positions as is: do not shorten them
Cs   --pmt[:emax=#] Adjustments to prepare for the PMT basis
Cs   --optics     : Add template optics category
Cs   --frzwf      : Extra switches to keep the basis set fixed
Cs   --ewald      : Sets EWALD_TOL
Cs   --clfloat    : Use classical algorithm to float Pnu
Cs   --pbe        : Set switches for PBE functional (libxc)
Cs   --pgf        : Add lines for lmpg --- implies ASA (not implemented)
Cs   --pr#1       : set verbosity to #1
Cs   --conv=#     : Set tolerance for energy convergence (ITER_CONV)
Cs   --convc=#    : Set tolerance for charge convergence (ITER_CONVC)
Cs   --dhftol=#   : Set tolerance correction to Harris forces (ITER_DHFTOL)
Cs   --convp      : conv, convc, dhftol assigned variables, controllable from command line
Cs   --omax=#1,#2,#3
Cs   --omax=#1    : set maximum sphere overlaps to #1 {#2 and #3 concern ES)
Cs   --rsmmx=#    : Upper bound to RSMH (sets HAM_AUTOBAS_RSMMX, used by lmfa)
Cs   --scala=#    : Multiply alat by #, divide lattice and site vectors by #
Cs   --scalp=#    : Multiply alat by factor, divide lattice and site vectors by same factor
Cs                : Factor set so that volume of plat = #
Cs   --show       : show input file after preprocessing
Cs   --wpos=      : --wpos=fn writes site positions to file fn
Cs   --ctrl=fn    : Write input file to fn.ext instead of actrl.ext
Cs   --rdsite=fn  : Read site data from file fn
Cs   --wsite      : Write site file with positions in Cartesian coordinates
Cs                : Site file is automatically created in express mode (express>0)
Cs   --wsitex[~sw]: Write site file with positions as multiples of plat
Cs                : Optional switches:
Cs                : ~quad1 shift site positions to first quadrant
Cs   --wsitex     : Write site file with positions as multiples of plat
Cs   --wsrmax=#   : Use WSRMAX as constraint when finding sphere radii
Cs   --xpos       : Save site positions in generated ctrl file as multiples of plat
Cs   --xshft=     : Shift site positions by a uniform translation, Cartesian coordinates
Cs   --xshftx=    : Ditto, but shift in lattice vector coordinates
Cs   Examples:
Cs    *Test different styles of input files.  init.sbsei can be found in testing
Cs     blm --nk=-100 --gmax=8.1 init.sbsei --express=#
Cs     You should find that a self-consistent density created from an input file
Cs     generated by express=0 is the same for express=0, 2, 4, 6, or 7
Cr    *Test scala, scalp, fixpos
Cs     blm --scala=8.2874732 --fixpos=1e-2 --gw --wsitex init.sbsei
Cs     blm --scalp=1 --fixpos=1e-2 --gw --wsitex init.sbsei
Cr    *Test other switches
Cs     blm --gw --addes --fixpos:tol=1e-2 --scalp=1 --xshftx=0,0,-0.0398106/2 --wsitex sbsei
Cl Local variables
Cl   lcart    : lcart(ib)  coordinates
Cl            :     0     Crystal coordinates
Cl            :     1     Cartesian coordinates
Cl            :     2     Conventional unit cell
Cr Remarks
Cr    Generates files actrl.ext and site.ext
Cu Updates
Cu   09 Apr 19 New LATTIC_ORIGIN
Cu   14 Jul 18 New --convp
Cu   15 Nov 17 asa has some optional switches
Cu   20 Sep 17 MIX has variable beta if express is 3 or less
Cu   20 Jul 17 New LDAU tags
Cu   11 Jul 17 New 1spec option to --findes, new option --addes,end
Cu   07 Jul 17 d local orbitals extended to In and Tl
Cu   05 Jul 17 Bug fix: initialize z,rmt,xpos,amom,pnu,pz,rsmh,pl,ql
Cu   01 Jul 17 New --gw~eloc=#, --gw~gcutb=#, --gw~gcutx=#, --gw~pbtol=strn --gw~nk=#
Cu   28 Jun 17 New --nosort or --sort=no
Cu   24 Jun 17 --findes has ~float option; default for gwemax now depends on volume/cell
Cu             New --gwemax=#; new PBTOL
Cu   19 May 17 wrsmax now has default 3.3
Cu   14 May 17 New --rsmmx --ehmx --pbesol -pbe
Cu   03 Mar 17 New --dv and user specified tokens in SPEC (SPEC_TOKEN)
Cu             IDMOD=.,.,4 for transition metals and lgw
Cu             small A for PBE or LREL>11
Cu   05 Jan 17 suctrl can read IDMOD,P,PZ from SPEC and REL from HAM and FSMOM from BZ
Cu   19 May 16 New --rdsite and --findes switches
Cu   22 Apr 16 suctrl redesigned for new 'express mode' input style
Cu   12 Apr 16 Bug fix in default values for gw gcutb, gcutx
Cu   05 Aug 15 NSPEC tag commented out when --wsite is used
Cu   15 Apr 14 New --mmag switch and rdtok 'MMOM'
Cu   25 Mar 14 Extra printout associated with --gw and new --nk switch
Cu   07 Feb 14 Some upgrades, including --xpos switch
Cu   29 Dec 13 More updates, including --wsite and -gw switch
Cu   19 Oct 12 Several updates, including --asa switch
Cu   08 Jul 07 Reads in effective command-line arguments from cat CMD
Cu   13 Dec 03 Uses plat(3) when scaling translation part of symgrp
Cu   06 Nov 01 Initially created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character prgnam*8
C ... Dynamically allocated arrays
      integer, parameter :: mxrec=5000,recl0=128,nspmx=200,nbmx=2000,ngmx=48,reclprs=500
      integer, parameter :: nrmx=5001, nxi0=10, n0=10
      integer,allocatable:: istab(:),iprm(:,:),lockc(:),ics(:),nrclas(:)
      real(8),allocatable:: wk(:),zc(:),rmtc(:)
      character,allocatable:: instr(:)*(recl0)
      character(len=8),allocatable :: slabll(:)
      character(len=72),allocatable :: slabl2(:)
      integer, allocatable :: ntab(:),iax(:)
C ... Local parameters
      logical lout,lasa,lfp,lgw,lgf,lpg,lpbe,lpbesol,lwtag,lplatcv,lfindes,lpmt,loptics
      integer ifin,ifout,swt(7),catlen,infonl,is,kount,j1,iv(10)
      integer mxcsiz,nttab
      double precision xx,xv(10)
      character*50 outs*128,ext,errmsg, cxx*10, alabl*8, prs*(reclprs)
      character(len=8) :: slabl(nspmx)
      real(8), target :: pos(3,nbmx)
      integer adloc(nspmx)
      double precision z(nspmx),rmt(nspmx),xpos(3,nbmx),amom(4,nbmx),pnu(4,nbmx),pz(4,nbmx),rsmh(4,nbmx),pl(n0),ql(n0,2)
      double precision uh(4,2,nspmx)
      character*(recl0) aa
      integer, target :: ips(nbmx)
      real(8) :: mxcst(nspmx) = 0d0
      integer nrc(nspmx),lcart(nbmx),lock(nspmx),lockes(nspmx),lamom(nspmx),lldau(nspmx),lidmod(nspmx)
      integer lpnu(nspmx),lpz(nspmx),lrsmh(nspmx),lmxbs(nspmx),idxdn(n0,nspmx),idmod(n0,nspmx),modep(3)
      integer helppr,infopr,stdo,i,k,nspec,ib,nbas,iz,ng,nggen,nl,nclass,mxclas,nspec0,nbas0
      integer iprc  ! verbosity of comments included in inputfile
      integer,parameter :: NULLI=-99999
      procedure(logical) cmdopt,a2bin
      procedure(integer) isw,a2vec,iosite,lgunit,iprint,nargf,fgtcat,rdtok,fopna,fxst,cmdoptswx,wordsw
      procedure(real(8)) dlength,avwsr,dglob
C ... For structures
!     include 'structures.h'
      type(str_ctrl) ::  s_ctrl
      type(str_lat) ::  s_lat

C ... Default input values
      integer :: nit, lmet=5, loc=1, mto=4 !, nldau=0
      integer :: nkabclm(0:3)=(/0,NULLI,0,0/), nkabcgw(0:3)=(/0,NULLI,0,0/)
      real(8) :: conv=1d-5, convc=3d-5, dhftol=NULLI, gcutb=NULLI, gcutx=NULLI, gwemax=NULLI
      logical :: ccor= .true., adnf=.false.
      integer :: gamma, sx=0, scr=4, lso=-1
      integer :: gwsig=12, fpandasa=0
      integer :: nz=16, gfmode=1, pwmode=0
      real(8) :: ef=0d0, pwemax=3d0, fsmom=NULLI, ehmx=NULLI, rsmmx=NULLI

C ... For the lattice and space group
      integer havep,quswt
      logical aunits,linfo,lshow,ltmp,addes,lfrzw,v6float
      logical :: lesend = .false. ! If T, species with atomic number 0 are written last
      logical :: lmino = .false.  ! If T, shift ES to minimize overlap after they are found
      logical :: lpfree = .false. ! If T, add HAM PMIN=-1
      integer :: lconvp = 0       ! If >0, tags EXPRESS conv, convc, dhftol are given variables

C     Following Stuttgart:
C     plat(1) = Lattice vectors actually used
C     plat(2) = Lattice vectors of the most compact unit cell.
C     plat(3) = Lattice vectors of the conventional unit cell.
C               i.e. that with maximum number of right angles.
C               The conventional cell may not be a primitive cell
C               But it is the set which the chemists use and
C               the set for which basis vectors are scaled
C               to convert to cartesian coordinates.
C              *If plat are specified, plat(3) assigned to plat(1)
C     plat(4) = Lattice vectors of the standard unit cell; see stplat
C     plat(5) = Lattice vectors of the standard unit cell, possibly
C               rotated (roplat.f)
      double precision a,b,c,alat(2),plat(3,3,5),qlat(3,3),qlatcv(9),rbohr,gmax,eloc,platm(3,3)
      double precision dorig(3),alpha,beta,gam,pltol,fptol,vol,sclwse
      double precision gen(9,ngmx),g(9,ngmx),ag(3,ngmx)
      integer express,lmxa,lorigin,lrel,nsp,ngen,lsmalla,mxesspec,mxaddes
      integer isym(7),ivec(10),it(10),lmxb(2),ishores(3)
      character*30 csym(7),grpnam*10,gens*128,strn*256,strn2*128,dc*1,varstrn*128,pbtol*64,fn*120
      character*72 spectok(nspmx),specincl(nspmx)  ! Use specified species tokens
      real(8) :: rmaxs = NULLI  ! (ASA) sets STR_RMAXS = range for structure constants
C     For determining sphere MT radii
      double precision omax1(3),omax2(3),omaxe(3),wsmax
C     for cnvsop
C     character*72 symops(ngmx)
      parameter(rbohr=0.529177d0)

C ... data initializations
      data csym/7*'undef'/

C ... Setup
      stdo = lgunit(1)
      varstrn = ' '
      forall (is = 1:nspmx) spectok(is) = ' '
      forall (is = 1:nspmx) specincl(is) = ' '
      forall (is = 1:nspmx) lmxbs(is) = NULLI
      lgf = cmdopt('--gf',4,0,outs)
      lpg = cmdopt('--pgf',5,0,outs)
      lpbe = cmdopt('--pbe',5,0,outs)
      lpbesol = cmdopt('--pbesol',8,0,outs)
      lpmt = cmdopt('--pmt',5,0,outs)
      lfrzw = cmdopt('--frzwf',7,0,outs)
      loptics = cmdopt('--optics',8,0,outs)
      v6float = cmdopt('--clfloat',9,0,outs)
      allocate(instr(mxrec))
      s_ctrl%rmines = NULLI; mxesspec = NULLI
      pbtol = '3e-4,3e-4,1e-3'
      call dpzero(z,nspmx)
      call dpzero(rmt,nspmx)
      call dpzero(xpos,3*nbmx)
      call dpzero(amom,4*nbmx)
      call dpzero(pnu,4*nbmx)
      call dpzero(pz,4*nbmx)
      call dpzero(rsmh,4*nbmx)
      call dpzero(pl,n0)
      call dpzero(ql,n0*2)


      lasa = cmdopt('--asa',5,0,outs) .or. lgf
      lfp = .not. lasa
      if (lasa) then
        dc = outs(6:6)
        if (wordsw(outs,dc,'pfree','',j1) > 0) lpfree = .true.
        k = wordsw(outs,dc,'rmaxs=','',j1) + 6
        if (k > 6) then
          if (a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,1,ivec,rmaxs) < 1)
     .      call rx('blm failed to parse '//trim(outs))
        endif
      endif

      if (lpmt) pwmode=11
      if (cmdopt('--pmt:emax=',11,0,outs)) then
        i = 11
        is = a2vec(outs,len(outs),i,4,', ',2,2,1,it,pwemax)
        if (is /= 1) call rx('blm failed to parse '//trim(outs))
      endif
      lgw = cmdopt('--gw',4,0,outs)
      if (lgw) then
        dc = outs(5:5)
        lfp = .true.
        i = cmdoptswx('--gw','pbtol=','') + 6
        if (i>6) pbtol = outs(i:)
        i = cmdoptswx('--gw','emax=','') + 4
        if (i>6) then
          is = a2vec(outs,len_trim(outs),i,4,', '//dc,3,2,1,ivec,gwemax)
          if (is /= 1) call rx('blm failed to parse '//trim(outs))
        endif
        i = cmdoptswx('--gw','gcutx=','') + 5
        if (i>6) then
          is = a2vec(outs,len_trim(outs),i,4,', '//dc,3,2,1,ivec,gcutx)
          if (is /= 1) call rx('blm failed to parse '//trim(outs))
        endif
        i = cmdoptswx('--gw','gcutb=','') + 5
        if (i>6) then
          is = a2vec(outs,len_trim(outs),i,4,', '//dc,3,2,1,ivec,gcutb)
          if (is /= 1) call rx('blm failed to parse '//trim(outs))
        endif
        i = cmdoptswx('--gw','nk=','') + 2
        if (i>6) then
          nkabcgw(0) = a2vec(outs,len(outs),i,2,', '//dc,3,2,3,ivec,nkabcgw(1))
          if (nkabcgw(0) < 0) call rx('blm failed to parse '//trim(outs))
        endif
      endif

C      nldau = isw(cmdopt('--ldau',4,0,outs))
C      ldauspec = ' '
C      if (nldau > 0) then
C        dc = outs(7:7)
C        print *, dc,' ',outs
C        k = index(outs,dc//'spec=')
C        if (k == 0) then
C          nldau = -1
C        else
C          ldauspec = outs(k+6:)
C        endif
C      endif
      lasa = .not. lfp
      call iinit(idxdn,size(idxdn))
      if (cmdopt('--fpandasa',10,0,outs) .or. cmdopt('--asaandfp',10,0,outs)) then
        lfp = .true.
        lasa = .true.
      endif
      lrel = -1
      if (cmdopt('--dv=',5,0,outs)) varstrn = outs(6:)

      lfindes = cmdopt('--findes',8,0,outs)
      dc = outs(9:9)
      k = wordsw(outs,dc,'rmin=','',j1) + 5
      if (k > 5) then
        if (a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,1,ivec,s_ctrl%rmines) < 1)
     .    call rx('blm failed to parse '//trim(outs))
      endif
      k = wordsw(outs,dc,'nspmx=','',j1) + 6
      if (k > 6) then
        if (a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,1,ivec,mxesspec) < 1)
     .    call rx('blm failed to parse '//trim(outs))
      endif
      mxaddes = 0
      k = wordsw(outs,dc,'nesmx=','',j1) + 6
      if (k > 6) then
        if (a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,1,ivec,mxaddes) < 1)
     .    call rx('blm failed to parse '//trim(outs))
      endif
      ishores = 0
      k = wordsw(outs,dc,'shorten=','',j1) + 8
      if (k > 8) then
        is = a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,3,iv,ishores)
        if (is /= 3) call rx('suctrl: failed to parse '//trim(outs))
      endif
      lmino = wordsw(outs,dc,'mino','',j1) > 0
      i = cmdoptswx('--findes','1spec','') + 5
      if (i > 5) mxesspec = 1
      omaxe = NULLI
      k = wordsw(outs,dc,'omax=','',j1) + 5
      if (k > 5) then
        is = a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,3,iv,omaxe)
        if (is < 0) call rx('suctrl: failed to parse '//trim(outs))
        if (is < 2) omaxe(2) = omaxe(1)
        if (is < 3) omaxe(3) = omaxe(2)
      endif
      sclwse = NULLI
      k = wordsw(outs,dc,'sclwsr=','',j1) + 7
      if (k > 7) then
        if (a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,1,ivec,sclwse) < 1)
     .    call rx('blm failed to parse '//trim(outs))
      endif

      addes = cmdopt('--addes',7,0,outs)
      if (addes) then
        dc = outs(8:8)
        lesend = wordsw(outs,dc,'end',dc//' ',k) > 0
      endif

      addes = cmdopt('--addes',7,0,outs) .or. lfindes

      ltmp = cmdopt('--express',9,0,outs)
      express = 6*isw(cmdopt('--express',9,0,outs))
      if (ltmp .and. outs(10:10) == '=') then
        i = 10
        if (a2vec(outs,len(outs),i,2,', ',2,2,1,it,express) < 0)
     .    call rx(' blm failed to parse ' // trim(outs))
      elseif (ltmp) then
        express = 6   ! For novices
      else
        express = 3   ! Large
      endif
      if (lpg) call rx('--pgf not implemented, sorry')
      nsp=1; if (cmdopt('--mag',5,0,outs)) nsp=2
      if (nsp == 2) lso = 0
      lamom = 0; lldau = 0; lpnu = 0 ; lpz = 0; lidmod = 0; lrsmh=0; uh = 0
      wsmax = 3.3d0
      fpandasa = isw(lasa .and. lfp)
      if (lasa .and. lfp .and. express > 4) then
        call fexit(-1,1,
     .    ' blm: for both ASA and FP, express mode must be 4 or smaller',0)
      endif
      iprc = 1 ! normal comments
      if (express >= 2 .and. express <= 7) iprc = 2 ! Lots of comments
      if (express < 2 .or. express >= 9) iprc = 0 ! No comments

      xx = dglob('lrquad',0d0,1)
      gens = ' '
      call stswt('init io,read lunit',stdo,swt)
      if (cmdopt('--pr',4,0,outs) .or. cmdopt('-pr',3,0,outs)) then
        i = 4
        if (cmdopt('-pr',3,0,outs)) i = 3
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,it,ivec)
        call pshpr(ivec)
      endif
      linfo = .false.
      infonl = -1
C     Turn on 'help' mode
      if (cmdopt('--input',7,0,outs)) then
        call stswt('io,help',0,swt)
        linfo = .true.
C     Print out as read in; infonl = number of 'info' newlines
      elseif (cmdopt('--show',6,0,outs)) then
        call stswt('io,rw',0,swt)
        infonl = 0
      endif
      lshow = quswt('io,w',swt) == 1
      if (cmdopt('--iactiv',7,0,outs)) call initqu(.true.)
      infopr = 30
      helppr = 999
      if (linfo) then
        helppr = 0
        infopr = 999
      endif

C ... Open init file for reading and ctrl file for output
      call fextg(ext)
      if (cmdopt('--rdsite',8,0,strn2)) then
        j1 = 9
        dc = strn2(j1:j1)
        outs = 'sitein'
        if (dc == '=') outs = strn2(j1+1:)
C       Read alat,plat,nbas,nspec, and labels
        nspec = 0
        j1 = iosite(135000,3d0,0,trim(outs),ifin,slabl,alat,plat,nbas,
     .    nspec,xx,xx,xx,xx,xx,xx,xx)
        if (j1 >= 0)
     .    j1 = iosite(8000,3d0,0,trim(outs),ifin,slabl,alat,plat,nbas,
     .    nspec,pos,xx,xx,xx,ips,xx,xx)
        if (j1 < 0) call rx('suctrl failed to read site file')
        ifin = 0 ! flag that nothing to be read from init file
        call info2(10,0,0,
     .      ' ... site file contains %i species and %i atoms',nspec,nbas)

C        call info(20,0,0,'     Find atomic numbers for each species ...',0,0)
C        do  i = 1, nspec
C          call zslabl(1,slabl(i),iz)
C          if (iz == -1) call rxs(
C     .      'species label does not correspond to any formula: ',slabl(i))
C          z(i) = iz
C        enddo
C        call dvset(rmt,1,nspec,0d0)

        plat(:,:,3) = 0
C       call dcopy(9,plat,1,plat(1,1,3),1)
        lcart(1:nbas) = 0
        goto 10  ! skip read init file

      elseif (fxst('INIT') /= 1) then
        errmsg = 'missing setup file, init'//ext
        goto 999
      endif
      ifin  = fopna('INIT',-1,1)

C ... get command-line arguments in CMD
      call stswt('cat,opt,mxrec',0,swt)
      instr = ' '
      k = fgtcat(ifin,'CMD ','command-line switches can be specified '
     .  //'in this category',swt,.true.,mxrec,recl0,instr)
      if (k > 0) then
        call acmdop(instr(1)(4:),len(instr)*mxrec-4,0)
        if (lshow) call acmdop(instr(1)(4:),len(instr)*mxrec-4,1)
      endif

C --- Lattice vectors and symmetry group ---
      call info(infopr,1,infonl,' ... Reading lattice data',0,0)

      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'LATTICE ',
     .  'Specifies the lattice vectors and symmetry group',
     .  swt,.true.,mxrec,recl0,instr)

      call info(helppr,1,0,'   The lattice vectors may be defined '
     .  //'indirectly through%N%5ftoken SPCGRP=, or explicitly using '//
     .  'PLAT=.',0,0)
      call info(helppr,0,0,'%3fIn the former case, the program'//
     .  ' completes the basis%N%5fto make it compatible with the'//
     .  ' space group.',0,0)
      call info(helppr,0,1,'%3fIn the latter case, the basis is '//
     .  'enlarged only if the%N%5fuser additionally '//
     .  'supplies group generators (token GENS=)',0,0)

C ... Look for space group
      call stswt('token,alt,fixlen',0,swt)
C     Try to read space group symbol
      havep = -1
      k = rdtok('SPCGRP=',instr,
     .  'The space group symbol or space group number (1..230)',
     .  ' ',', ',1,swt,1,1,xx,grpnam,prs)
      if (k == 1) havep = 0
C     Choice of origin
      if (k == 1 .or. linfo) then
        call stswt('token,opt',0,swt)
        lorigin = 1
        i = rdtok('ORIGIN=',instr,
     .    '%i,,4;;Choice of origin (1 or 2), if space group is also supplied',
     .    ' ',', ',2,swt,1,1,lorigin,cxx,prs)
        call sanrg(.true.,lorigin,1,2,'','choice of origin')
      endif

C     Or primitive lattice vectors
      call stswt('token,reqd,fixlen',0,swt)
      if (k == 0 .or. linfo) then
        k = rdtok('PLAT=',instr,
     .  ' %;7d,,Primitive lattice vectors (dimensionless)',
     .  ' ',' ,',4,swt,1,9,plat,cxx,prs)

        havep = 1
        if (.not. lshow)
     .    call info(infopr,0,infonl,' ... read plat',0,0)

        call stswt('token,opt',0,swt)
        k = rdtok('GENS=',instr,
     .  'Generators of the space group',' ',', ',1,swt,1,1,xx,gens,prs)
        if (k == 0) then
          call info(infopr,0,0,' ... no group generators',0,0)
        elseif (.not. lshow) then
          call info(infopr,0,0,' ... read generators',0,0)
        endif
      endif

C ... Check whether units are in AA or a.u.
      call stswt('token,opt',0,swt)
      aunits = (1 == rdtok('UNITS=A',instr,
     .  'Specifies that lattice dimensions are in Angstrom',' ',', ',1,
     .  swt,1,0,xx,cxx,prs))

      call info(helppr,1,0,'   ... If the lattice vectors are '
     .  //'defined through token SPCGRP=,%N%7fone or more of the '//
     .  'following is required.',0,0)
      a = 0
      b = 0
      c = 0
      alpha = 0
      beta = 0
      gam = 0
      if (havep /= 1 .or. linfo) then
      call stswt('token,opt,fixlen',0,swt)
      k = rdtok('A=',instr,
     .  '%;7d,,`A'' lattice parameter, in a.u. (Angstrom if UNITS=A)',
     .  ' ',', ',4,swt,1,1,a,cxx,prs)
      k = rdtok('B=',instr,'%;7d,,`B'' lattice parameter',
     .  ' ',', ',4,swt,1,1,b,cxx,prs)
      k = rdtok('C=',instr,'%;7d,,`C'' lattice parameter',
     .  ' ',', ',4,swt,1,1,c,cxx,prs)
      k = rdtok('ALPHA=',instr,'%;7d,,angle(b,c), in degrees',
     .  ' ',', ',4,swt,1,1,alpha,cxx,prs)
      k = rdtok('BETA=',instr,'%;7d,,angle(c,a), in degrees',
     .  ' ',', ',4,swt,1,1,beta,cxx,prs)
      k = rdtok('GAMMA=',instr,'%;7d,,angle(a,b), in degrees',
     .  ' ',', ',4,swt,1,1,gam,cxx,prs)
      endif

      call info(helppr,1,0,'   ... If the lattice vectors are '
     .  //'defined through token PLAT=,%N%7fthe '//
     .  'following is required.',0,0)
      alat(1) = 0
      call stswt('token,o,fixlen',0,swt)
      if (havep == 0) call stswt('token,i',0,swt)
      if (havep == 1) call stswt('token,r',0,swt)
      k = rdtok('ALAT=',instr,
     .  '%;8g,,scaling of lattice vectors PLAT, in a.u.'//
     .  ' (Angstrom if UNITS=A)',' ',', ',4,swt,1,1,alat,cxx,prs)

C ... Construct the lattice vectors from input data
C     Lattice vectors are supplied; construct alternate sets
      if (havep == 1 .and. .not.linfo) then
        alat(2) = alat(1)
        isym(3) = 1
        call mkplat(alat,csym,isym,.true.,plat)
C       If plat specified, use them to define pos = plat * X
        call dcopy(9,plat(1,1,1),1,plat(1,1,3),1)

C     Find space group; plat from isym, a,b,c and alpha.,beta,gamma
      elseif (havep == 0 .and. .not.linfo) then
        call info(infopr,1,0,' ... Generating space group from SPCGRP='//grpnam,0,0)
        gens = ' '
        call gengrp(csym,dorig,gens,grpnam,-1,lorigin,isym,ngen,plat,qlat)
        if (isym(7) == 0) call rxs('illegal space group',grpnam)
        call roplat(a,alat,alpha,b,beta,c,gam,isym,
     .    plat(1,1,3),plat(1,1,5))
C        if (cmdopt('--conventional',14,0,outs)) then
C          call dcopy(9,plat(1,1,3),1,plat(1,1,5),1)
C        endif
        call dcopy(9,plat(1,1,5),1,plat(1,1,1),1)
        call dcopy(9,plat(1,1,5),1,plat(1,1,4),1)
        call cpplat(plat(1,1,1),plat(1,1,2))
      endif
C ... Convert to atomic units, if not using them already
      if (aunits) alat(1) = alat(1)/rbohr

C --- HAM category ---
      call stswt('cat,opt,mxrec',0,swt)
      k = fgtcat(ifin,'HAM ',
     .  'Optional parameters in the HAM category',
     .  swt,.true.,mxrec,recl0,instr)
      if (k > 0 .or. linfo) then
        call info(infopr,0,0,' ... Found and reading from HAM category ...',0,0)
        call stswt('token,opt',0,swt)
        k = rdtok('REL=',instr,'%i,,4;;Specifies relativistic treatment',
     .    ' ',', ',2,swt,1,1,lrel,cxx,prs)
      endif

C --- BZ category ---
      call stswt('cat,opt,mxrec',0,swt)
      k = fgtcat(ifin,'BZ ',
     .  'Optional parameters in the BZ category',
     .  swt,.true.,mxrec,recl0,instr)
      if (k > 0 .or. linfo) then
        call info(infopr,0,0,' ... Found and reading from BZ category ...',0,0)
        call stswt('token,opt',0,swt)
        k = rdtok('FSMOM=',instr,'%i,,4;;Specifies relativistic treatment',
     .    ' ',', ',4,swt,1,1,fsmom,cxx,prs)
      endif

C --- SITE category ---
      call info(infopr,1,infonl,' ... Reading site data ...',0,0)

      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'SITE ','Specifies basis.  At least one site must be specified.',
     .  swt,.true.,mxrec,recl0,instr)

C ... Read in site positions
      nspec = 0
      ib = 1
      catlen = swt(2)
      swt(2) = 1
   40 continue
        swt(1) = swt(2)
        swt(2) = catlen
        call stswt('token,opt,fixlen,mult',0,swt)
        k = rdtok('ATOM=',instr,
     .    'Species label.  There is one token for each site.',
     .    ' ',', ',1,swt,1,1,xx,alabl,prs)
        call stswt('token,nomult',0,swt)
C   ... No further labels found -> terminates infinite loop
        if (k == 0) goto 42

C   ... Add new species, if this label is new
        if (.not.linfo) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',i,k,.false.)
          if (i < 0) then
             i = nspec
             nspec = nspec+1
             slabl(nspec) = alabl
           endif
          ips(ib) = i+1
        endif
        call info(helppr,0,0,'%5fFollowing each label, '//
     .    'information about the particular site is supplied:',0,0)
        call stswt('token,alt,fixlen',0,swt)
C       Basis vectors as multiples of plat
        lcart(ib) = 0
        k = rdtok('X=',instr,'4;;Site positions, as (fractional) '//
     .    'multiples of the lattice vectors',' ',' ,',4,swt,1,3,
     .    pos(1,ib),cxx,prs)
C       Basis vectors as multiples of the conventional unit cell
        if (k == 0 .or. linfo) then
          k = rdtok('C=',instr,'4;;Site positions, as (fractional) '//
     .      'multiples of the lattice vectors of the conventional unitcell',
     .      ' ',' ,',4,swt,1,3,pos(1,ib),cxx,prs)
          if (k /= 0) lcart(ib) = 2
        endif
C       Basis vectors in cartesian coordinates
        call stswt('token,reqd,fixlen',0,swt)
        if (k == 0 .or. linfo) then
          k = rdtok('POS=',instr,'4;;Site positions, cartesian '//
     .      'coordinates',' ',' ,',4,swt,1,3,pos(1,ib),cxx,prs)
          lcart(ib) = 1
        endif

C       print *, ib,slabl(ib),z(ib),rmt(ib)
        ib = ib+1
        if (ib > nbmx) call rxi(
     .    'too many sites. Increase nbmx in suctrl; now',nbmx)
        if (nspec > nspmx) call rxi(
     .    'too many species. Increase nspmx in suctrl; now',nspmx)
        if (.not.linfo) goto 40
C ... End of loop
   42 continue
      nbas = ib-1

      if (.not.linfo) then
      if (nbas == 0) then
        call rx('No sites found')
      endif
      if (iprint() >= infopr) call
     .  awrit3('  %i site%?;(n==1);;s;, %i species found',' ',80,stdo,nbas,nbas,nspec)

      call info2(infopr+10,0,0,'%5fScaling and shortening basis vectors',0,0)
C     qlat = (plat+)^-1
      do  ib = 1, nbas
        k = 1; if (lcart(ib) == 2) k = 3
        call dinv33(plat(1,1,k),1,qlat,xx)

C       posc+ = plat posp+
        if (lcart(ib) == 0 .or. lcart(ib) == 2) then
          call dcopy(3,pos(1,ib),1,xv,1)
          call dmpy(plat(1,1,k),3,1,xv,3,1,pos(1,ib),3,1,3,1,3)
        endif
        if (.not. cmdopt('--noshorten',11,0,outs) .and. .not. cmdoptswx('--shorten','=no','')>0) then
          call shorbz(pos(1,ib),pos(1,ib),plat(1,1,k),qlat)
        endif
      enddo
      endif

C ... Re-entry point for --rdsite
   10 continue

C     Initialize values z,rmt to mark they haven't been set
      call dvset(z,1,nspec,-1d0)
      call dvset(rmt,1,nspec,0d0)

C --- Override default convc, conv, dhftol ---
      i = 7
      if (cmdopt('--conv=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,conv)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

      i = 8
      if (cmdopt('--convc=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,convc)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

      i = 9
      if (cmdopt('--dhftol=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,dhftol)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

      lconvp = isw(cmdopt('--convp',7,0,strn)) + isw(cmdopt('--molstat',9,0,strn))

C --- Set constraining wsrmax ---
      i = 9
      if (cmdopt('--wsrmax=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,wsmax)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

C --- Set rsmmx ---
      i = 8
      if (cmdopt('--rsmmx=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,rsmmx)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

C --- Set ehmx ---
      i = 7
      if (cmdopt('--ehmx=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,ehmx)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

C --- Scale alat and plat by specified factor ---
      i = 8
      if (cmdopt('--scala=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,xx)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
        alat(1) = alat(1)*xx
        call info(infopr,0,0,'     scaling alat, inverse plat by factor %d',xx,0)
        call dscal(9,1/xx,plat,1)
        call dscal(3*nbas,1/xx,pos,1)
      endif
      i = 8
      if (cmdopt('--scalp=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,xx)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
        call dinv33(plat,1,qlat,xv)
        xx = (xv(1)/xx)**(1d0/3d0)
        call info(infopr,0,0,'     scaling alat, inverse plat by factor %d',xx,0)
        alat(1) = alat(1)*xx
        call dscal(9,1/xx,plat,1)
        call dscal(3*nbas,1/xx,pos,1)
      endif

      i = 8
      if (cmdopt('--xshft=',i,0,strn) .or. cmdopt('--xshftx=',i+1,0,strn)) then
        ltmp = cmdopt('--xshftx=',i+1,0,strn)
        if (ltmp) i=i+1
        is = a2vec(strn,len(strn),i,4,', ',2,2,3,it,xpos)
        if (is /= 3) call rx('blm failed to parse '//trim(strn))
        xv(1:3) = xpos(:,1)
        if (ltmp) call dgemm('T','N',3,1,3,1d0,plat,3,xpos,3,0d0,xv,3)
        do  ib = 1, nbas
          pos(1:3,ib) = pos(1:3,ib) + xv(1:3)
        enddo
      endif

C --- SPEC category ---
      if (ifin == 0) goto 34  ! data not read from ifin
      call info(infopr,1,infonl,' ... Reading species data ...',0,0)

      call stswt('cat,opt,mxrec',0,swt)
      k = fgtcat(ifin,'SPEC ',
     .  'Specifies additional species information '//
     .  '%N%3f(overriding default values or possible input from '//
     .  'SITE category)'//
     .  '%N%3fLabels should coincide with labels in SITE category%N',
     .  swt,.true.,mxrec,recl0,instr)

C ... Make a list of species
      catlen = swt(2)
      swt(2) = 1
      is = 0
      kount = 0
   30 continue
        swt(1) = swt(2)
        swt(2) = catlen
        call stswt('token,opt,fixlen,mult',0,swt)
        k = rdtok('ATOM=',instr,
     .    'Species label.  Enter one for each species.',
     .    ' ',', ',1,swt,1,1,xx,alabl,prs)
        call stswt('token,nomult',0,swt)
C   ... No further labels found -> terminates infinite loop
        if (k == 0) goto 32
        call info(helppr,0,0,'%5f... Following each label, the'//
     .    ' following may be supplied:',0,0)

C       If this species missing from list of species:
        if (.not.linfo) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',is,k,.false.)
          if (is < 0) then
            call info(infopr-20,0,0,'%N%5f(warning) ... label `'
     .        //alabl//'%a'' does not match any sites ... '//
     .        'species ignored',0,0)
            goto 30
          else
            is = is+1
          endif
        endif
        kount = kount+1

        k = rdtok('Z=',instr,'%;7d,,4;;Atomic number.  If not '//
     .    'supplied, Z is inferred from the label',' ',', ',4,swt,
     .    is,1,z,cxx,prs)
        k = rdtok('R=',instr,'%;7d,,4;;Sphere radius',' ',', ',
     .    4,swt,is,1,rmt,cxx,prs)

        call stswt('token,opt,varlen',0,swt)
        k = rdtok('MMOM=',instr,'%;7d,,4;;Magnetic moments by l',
     .    ' ',' ,',4,swt,1,4,amom(1,is),cxx,prs)
        lamom(is) = iabs(k)

!       if (nsp == 2) then
          k = rdtok('UH=',instr,'%;7d,,4;;Hubbard parameters U by l',
     .      ' ',' ,',4,swt,1,4,uh(1,1,is),cxx,prs)
          lldau(is) = iabs(k)
          k = rdtok('JH=',instr,'%;7d,,4;;Hubbard parameters J by l',
     .      ' ',' ,',4,swt,1,4,uh(1,2,is),cxx,prs)
          lldau(is) = max(lldau(is),iabs(k))
!       endif

        k = rdtok('P=',instr,'%;7d,,4;;Log derivative parameters by l',
     .    ' ',' ,',4,swt,1,4,pnu(1,is),cxx,prs)
        lpnu(is) = iabs(k)

        k = rdtok('PZ=',instr,'%;7d,,4;;Local orbitals by l',
     .    ' ',' ,',4,swt,1,4,pz(1,is),cxx,prs)
        lpz(is) = iabs(k)

        k = rdtok('RSMH=',instr,'%;7d,,4;;Envelope smoothing radii by l',
     .    ' ',' ,',4,swt,1,4,rsmh(1,is),cxx,prs)
        lrsmh(is) = iabs(k)

        k = rdtok('IDMOD=',instr,'%i,,4;;Freeze log derivative parameters',
     .    ' ',' ,',2,swt,1,4,idmod(1,is),cxx,prs)
        lidmod(is) = iabs(k)

        k = rdtok('LMX=',instr,'%i,,4;;l-cutoff for basis',' ',', ',2,swt,1,1,lmxbs(is),cxx,prs)

        k = rdtok('CSTRMX=',instr,'%;7d,,4;;Constrain MT radius when finding sphere radii',' ',', ',
     .    4,swt,is,1,mxcst,cxx,prs)

C       call stswt('token,fixlen',0,swt)
        k = rdtok('TOKEN=',instr,
     .  ',,4;;User-specified tokens for this species',' ',' ',1,swt,1,1,xx,spectok(is),prs)

        k = rdtok('INCLUDE',instr,
     .  ',,4;;Add contents of file to this species%N%7fCan include directives, but be sure to preserve indentation.',
     .    ' ',' ',1,swt,1,1,xx,specincl(is),prs)

        if (.not.linfo) goto 30
C ... End of loop
   32 continue
      if (iprint() >= infopr) call awrit2
     .  ('%?;n;%N%5f; ;species data read for %i species',' ',80,stdo,lshow,kount)

   34 continue
      call info(infopr+20,0,0,'     Find missing atomic numbers ...',0,0)
      if (.not.linfo) then
      do  i = 1, nspec
        if (z(i) == -1) then
          call zslabl(1,slabl(i),iz)
          if (iz == -1) call rxs(
     .    'species label does not correspond to any formula: ',slabl(i))
          z(i) = iz
        endif
        adloc(i) = 0
C       P.Q.N of valence orbital for d state, to be augmented by local orbital, included as PZ=..
        if (z(i) >= 21 .and. z(i) <= 30) adloc(i) = adloc(i) + 300  ! Sc to Cu, Zn. No Ga since 4d is valence
        if (z(i) >= 39 .and. z(i) <= 49) adloc(i) = adloc(i) + 400  ! Y  to Ag, Cd, In
        if (z(i) >= 71 .and. z(i) <= 81) adloc(i) = adloc(i) + 500  ! Lu to Au, Hg, Tl
C       P.Q.N of valence orbital for d state, to be included as P=n+.92
C        if (z(i) == 49) adloc(i) = adloc(i) + 4000000  ! In ... not needed as lmfa uses this as default anyway
C        if (z(i) == 81) adloc(i) = adloc(i) + 5000000  ! Tl ... not needed as lmfa uses this as default anyway
      enddo
      endif

      if (linfo) return

C --- Complete basis for compatibility with supplied symmetry group ---
      call info(infopr,1,0,' ... Complete basis for supplied symmetry',0,0)
      fptol = 1d-5
      if (cmdopt('--fixpos',8,0,strn)) then
        j1 = 8+1
        if (strn(9:13) == ':tol=') then
          j1 = 13
        endif
        if (strn(9:9) /= ':' .or.
     .    .not. a2bin(strn,fptol,4,0,' ',j1,len(strn))) fptol = 1d-5
      endif
      if (cmdopt('--fixlat',8,0,strn)) then
        j1 = 8+1
        if (strn(9:13) == ':tol=') then
          j1 = 13
        endif
        if (strn(9:9) /= ':' .or.
     .    .not. a2bin(strn,pltol,4,0,' ',j1,len(strn))) pltol = 1d-5
      else
        pltol = 0
      endif

      allocate(istab(ngmx*nbmx))
      if (gens == ' ') gens = 'i*i'
      call gensym(slabl,gens,0,.true.,.false.,pltol,fptol,.false.,nbas,nspec,
     .  ngmx,plat,plat(1,1,3),0,xx,pos,ips,nrc,ng,g,ag,ngen,gen,aa,
     .  nggen,isym,istab)
      deallocate(istab)
C     At this point we can determine the space group name number
C     but to do so we need the following routines taken from
C     Stuttgart code:
C       cnvsop to make a set of generators (Stuttgart conventions)
C       pntgrp fills out isym(4..6), in particular isym(6)
C     Then we can call gengrp to get space group number
C     But this takes too much stuff for something not essential.
C     if (havep == 1) then
C       call cnvsop(.true.,plat,symops,g,ag,ng,1)
C       call pntgrp(csym,isym,ng,symops)
C       call gengrp(csym,dorig,gens,grpnam,1,k,isym,ngen,plat,qlat)
C     endif

      if (.not. cmdopt('--noshorten',11,0,outs) .and. .not. cmdoptswx('--shorten','=no','')>0) then
        call info2(infopr+10,1,0,' ... Shortening basis vectors',0,0)
        call dinv33(plat,1,qlat,xx)
        do  ib = 1, nbas
          call shorbz(pos(1,ib),pos(1,ib),plat,qlat)
        enddo
      endif

C --- Sort basis according to species index ---
      if (.not. cmdopt('--nosort',8,0,outs) .and. .not. cmdoptswx('--sort','no','')>0) then
        allocate(iprm(nbas,1),wk(3*nbas))
        call ivheap(1,nbas,ips,iprm,101)
        call ivprm(1,nbas,ips,wk,iprm,1)
        call dvprm(3,nbas,pos,wk,iprm,1)
        deallocate(iprm,wk)
      endif

C --- Printout final lattice information ---
      if (iprint() >= infopr) then
        lplatcv = dlength(9,plat(1,1,3),1) /= 0
        if (lplatcv) call dinv33(plat(1,1,3),1,qlatcv,xx)
        call dinv33(plat,1,qlat,xx)
        call dgemm('T','N',3,3,3,1d0,plat(1,1,3),3,qlat,3,0d0,platm,3)
        vol = abs(xx*alat(1)**3)
        call awrit2('%N Lattice vectors : alat = %,6,6g a.u.  vol=%,6,6g a.u.',
     .    ' ',80,stdo,alat(1),vol)
        write(stdo,"(/t18,'Plat',t43,'Conventional unit cell',10x,'As multiples of Plat')")
        if (lplatcv) then
          do  k = 1, 3
            call info5(2,0,0,'%f%3;11,7D%2f%3;11,7D%2f%3;8,4D',plat(1,k,1),plat(1,k,3),platm(1,k),4,5)
          enddo
        else
          do  k = 1, 3
            call info2(2,0,0,'%f%3;11,7D',plat(1,k,1),0)
          enddo
        endif
        write(stdo,"(/1x,'Basis vectors after sorting and shortening:')")
        write(stdo,357)
  357   format(' site spec',8x,'pos (Cartesian coordinates)',9x,
     .  'pos (multiples of plat)')
        do  i = 1, nbas
C         posp = posc (plat+)^-1  and  posp+ = (plat)^-1 posc+
          call dgemm('T','N',3,1,3,1d0,qlat,3,pos(1,i),3,0d0,xv,3)
          call info5(2,0,0,'%,4i  '//slabl(ips(i))//'%;10,6D%2;11,6D%f%3;11,6D',
     .      i,pos(1,i),pos(2,i),xv,0)
C         print 345, i, slabl(ips(i)), (pos(k,i),k=1,3), (xv(k),k=1,3)
          xpos(1:3,i) = xv(1:3)
          if (lplatcv) then
          call dgemm('T','N',3,1,3,1d0,qlatcv,3,pos(1,i),3,0d0,xv,3)
C         if (iprint() >= infopr+10) print 346, (xv(k),k=1,3)
          if (iprint() >= infopr+10) call info2(2,0,0,
     .      '%13fin terms of conventional unit cell%3;11,6D',xv,0)
          endif
  345     format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
  346     format(13x,'in terms of conventional unit cell',3f11.6)
        enddo
      endif

      if (cmdopt('--wpos=',7,0,strn)) then
          call iopos(.true.,0,strn(8:),nbas,pos,xx)
C         call fclr(strn(8:),-1)
        endif

      if (linfo) return

C ----------------------- rearrange species, E at end ----------------
      if (lesend) then
        allocate(iprm(nspec,2),slabll(nspec),ics(nbas),wk(n0*nspec),slabl2(nspec))
        is = 0
        do                  ! Loop twice, 1st time for z>0, second for z=0
        do  i = 1, nspec
          if (z(i)>0.5d0 .and. lesend .or. z(i)<=0.5d0 .and. .not.lesend) then
            is = is+1
            iprm(i,1) = is
            iprm(is,2) = i  ! Need both forward and inverse permutation
          endif
        enddo
        if (.not. lesend) exit
        lesend = .false.
        enddo
        lesend = .true.
        call icopy(nbas,ips,1,ics,1)
        do  i = 1, nbas
          is = ips(i)
          ips(i) = iprm(is,1)
        enddo
        forall (is = 1:nspec) slabll(is) = slabl(is)
        forall (is = 1:nspec) slabl(iprm(is,1)) = slabll(is)
        forall (is = 1:nspec) slabl2(is) = spectok(is)
        forall (is = 1:nspec) spectok(iprm(is,1)) = slabl2(is)
        forall (is = 1:nspec) slabl2(is) = specincl(is)
        forall (is = 1:nspec) specincl(iprm(is,1)) = slabl2(is)
        call dvprm(1,nspec,z,wk,iprm(1,2),1)
        call dvprm(4,nspec,amom,wk,iprm(1,2),1)
        call dvprm(4,nspec,pnu,wk,iprm(1,2),1)
        call dvprm(4,nspec,pz,wk,iprm(1,2),1)
        call dvprm(4,nspec,rsmh,wk,iprm(1,2),1)
        call ivprm(n0,nspec,idmod,wk,iprm(1,2),1)
        deallocate(iprm,slabll,ics,wk,slabl2)
      endif

C ----------------------- Determine sphere radii ----------------
      call info(10,1,0,' ... Make sphere radii',0,0)
      xx = dglob('lrel',1d0,1)
      xx = dglob('nsp',1d0,1)
      xx = dglob('lxcf',2d0,1)
C     Determine initial guess for sphere radii
      call iinit(lock,nspec)
      do  is = 1, nspec
        if (mxcst(is) == 0) cycle
        lock(is) = 2
        rmt(is) = mxcst(is)
      enddo
      forall (i=1:3) modep(i) = 2
      call makrm0(101,nspec,nbas,alat,plat,pos,slabl,ips,modep,lock,z,rmt)
      call dpzero(omax1,3)  ! Scale sphere radii to touching
      if (lasa) omax1 = (/.16d0,.18d0,.20d0/) ! ASA overlaps
      if (lasa .and. lfp) omax1 = (/.10d0,.10d0,.10d0/)  ! In between

      if (cmdopt('--omax=',7,0,outs)) then
        i = 7
        is = a2vec(outs,len_trim(outs),i,4,', ',2,3,3,iv,omax1)
        if (is < 0) call rx('suctrl: failed to parse '//trim(outs))
        if (is < 2) omax1(2) = omax1(1)
        if (is < 3) omax1(3) = omax1(2)
      endif
      omax2=(/1d0,1d0,1d0/)
      call sclwsr(100,nbas,nbas,nspec,alat,plat,pos,ips,modep,slabl,z,
     .  lock,1d0,wsmax,omax1,omax2,rmt)

C ----------------------- Find empty spheres ----------------
      nspec0 = nspec
      nbas0 = nbas
      if (lfindes) then
        nclass = nspec
        mxclas = 5*nbas
        allocate(zc(mxclas),rmtc(mxclas),lockc(mxclas),ics(nclass),nrclas(mxclas))
        do  i = 1, nclass
          ics(i) = i
          nrclas(i) = 0
        enddo
        do  i = 1, nbas
          nrclas(ips(i)) = nrclas(ips(i))+1
        enddo

C       Initialize strux
        call ptr_ctrl(s_ctrl,0,' ',0,0,0,xx)
        call ptr_lat(s_lat,0,' ',0,0,0,xx)
C        call bcast_strx(2**9+2**6,s_bz,,s_ctrl,s_bz,s_bz,s_lat,
C     .    s_bz,s_bz,s_bz,s_bz,nbmx,nbmx)

C       Fill in elements of s_ctrl and s_lat needed by findes
        s_ctrl%nesabc(1:3) = 100
        if (s_ctrl%rmines == NULLI) s_ctrl%rmines = 1
        s_ctrl%rmaxes = 2.5d0
        s_ctrl%modep(1:3) = 2
        s_ctrl%sclwsr = 0
        s_ctrl%omax1 = omax1
        s_ctrl%omax2 = omax2
        s_ctrl%wsrmax = 0d0
        s_ctrl%nspec = nspec
        allocate(s_ctrl%clabl(nbmx),s_ctrl%spid(nbmx))
        do  i = 1, nclass
          s_ctrl%clabl(i) = slabl(i)
          s_ctrl%spid(i) = slabl(i)
        enddo
        call ptr_ctrl(s_ctrl,4+1,'ipc',nbas,0,0,ips)
        call ptr_ctrl(s_ctrl,4+1,'ips',nbas,0,0,ips)

        s_lat%avw = avwsr(plat,alat,vol,nbas)
        call ptr_lat(s_lat,1,'pos',3,nbmx,0,0)
        call dcopy(3*nbas,pos,1,s_lat%pos,1)
        s_lat%plat0(:,:) = plat(:,:,1)
        s_lat%qlat = qlat

        if (omaxe(1) /= NULLI) then
          call info2(2,0,0,' temporary scaling of spheres while finding ES: use omax =%s,%3:1;1d%%',100*omaxe,2)
          call sclwsr(100,nbas,nbas,nspec,alat,plat,pos,ips,modep,slabl,z,
     .      lock,1d0,wsmax,omaxe,omax2,rmt)
          s_ctrl%omax1 = omaxe
        endif

        lockes(1:nspec) = 2
        nl = 2

        call spec2c(nspec,nclass,ics,rmt,rmtc,z,zc,lockes,lockc)
        call findes(1,mxaddes,NULLI,s_ctrl,s_lat,alat,nbas,nclass,nl,
     .    nrclas,5*nbas,ng,plat,g,ag,ishores,lockc,rmtc,zc)
        if (nbas > nbas0) then
          call dcopy(3*(nbas-nbas0),s_lat%pos(1,nbas0+1),1,pos(1,nbas0+1),1)
        endif
        do  i = nspec+1, nclass
          slabl(i) = s_ctrl%clabl(i)
          z(i) = 0
        enddo

        call dcopy(nclass,rmtc,1,rmt,1)
        deallocate(lockc,zc,rmtc)

        call icopy(nbas,s_ctrl%ips,1,ips,1)
        nspec = nclass

        if (lesend .and. z(nspec0) <= 0.5d0) nspec0 = nspec0-1 ! Maybe merge last species with new ES
        if (mxesspec /= NULLI .and. nspec-nspec0 > mxesspec) then
          do  i = nbas0+1, nbas
            if (ips(i) > nspec0+mxesspec) then
              ips(i) = nspec0+mxesspec
            endif
          enddo
          rmt(nspec0+mxesspec) = min(rmt(nspec0+mxesspec),rmt(nspec)) ! choose the smallest size
          nspec = nspec0 + mxesspec
        endif

        if (omaxe(1) /= NULLI .or. (lasa .and. .not. lfp)) then
          s_ctrl%omax1 = omax1
          if (sclwse /= NULLI) s_ctrl%sclwsr = sclwse
          call info2(2,0,0,' rescale spheres use omax =%s,%3:1;1d%%  sclwsr = %d',100*omax1,s_ctrl%sclwsr)
          k = s_ctrl%sclwsr/10
          call sclwsr(10*k,nbas,nbas,nspec,alat,plat,pos,ips,modep,slabl,z,
     .      lock,1d0,wsmax,omax1,omax2,rmt)
        endif

        if (lmino) then
          mxcsiz = 0
          call pshpr(iprint()-20)
          allocate(ntab(nbas+1))
          call pairs(nbas,nbas,alat(1),plat,[5d0],pos,
     .      [-1],3,-1,[0],nttab,ntab,iax,mxcsiz)
          call poppr

          call ovmin('~z',nbas,nbas,alat,plat,rmt,rmt,
     .      slabl,ips,modep,z,ntab,iax,pos,1)
        endif

      endif

C ----------------------- Create the ctrl file ----------------
      if (.not. cmdopt('--ctrl=',7,0,outs)) then
        outs = '--ctrl=actrl'
      endif
      call info(10,1,0,' ... Create input file '//trim(outs(8:))//trim(ext)//
     .  ' (express mode %i)',express,0)
      call sanrg(.true.,express,0,7,'','express mode')
      ifout = fopna(trim(outs(8:)),-1,2)
      rewind ifout
C     print *, '!!'; ifout = stdo
      if (iprc >= 1 .or. .true.) then
        write(ifout,"('# Autogenerated from ',a/'#')",advance='no')
     .    'init'//trim(ext)//' using:'
        j1 = nargf()
        do  i = 0, j1-1
          call getarf(i,strn)
          write(ifout,"(1x,a)",advance='no') trim(strn)
        enddo
        write(ifout,*)
      endif
      lout = .true.
      call stswt('io,write lunit',ifout,swt)

C ... Default values
      nit = 10
      if (cmdopt('--nit=',6,0,outs)) then
        i = 6
        if (a2vec(outs,len(outs),i,2,', ',2,2,1,it,nit) < 0)
     .    call rx(' blm failed to parse ' // trim(outs))
      endif

      mto = 4
      if (cmdopt('--mto=',6,0,outs)) then
        i = 6
        if (a2vec(outs,len(outs),i,2,', ',2,2,1,it,mto) < 0)
     .    call rx(' blm failed to parse ' // trim(outs))
      endif

      loc = 1
      if (cmdopt('--loc=',6,0,outs)) then
        i = 6
        if (a2vec(outs,len(outs),i,2,', ',2,2,1,it,loc) < 0)
     .    call rx(' blm failed to parse ' // trim(outs))
      endif


      gmax = NULLI
      if (cmdopt('--gmax=',7,0,strn)) then
        i = 7; is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,gmax)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

      eloc = NULLI
      if (cmdopt('--eloc=',7,0,strn)) then
        i = 7; is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,eloc)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

      if (cmdopt('--gwemax=',9,0,strn)) then
        i = 9; is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,gwemax)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
      endif

      if (cmdopt('--nk=',5,0,strn)) then
        i = 5
        nkabclm(0) = a2vec(strn,len(strn),i,2,', ',2,2,3,it,nkabclm(1))
        if (nkabclm(0) < 0) call rx('blm failed to parse '//trim(strn))
      endif

      if (cmdopt('--nkgw=',7,0,strn)) then
        i = 7
        nkabcgw(0) = a2vec(strn,len(strn),i,2,', ',2,2,3,it,nkabcgw(1))
        if (nkabcgw(0) < 0) call rx('blm failed to parse '//trim(strn))
      endif
      if (nkabclm(0) > 0 .and. nkabcgw(0) == 0) nkabcgw = nkabclm

C ... Mode- or program-dependent defaults
      if (lasa .and. .not. lfp) lmet = 1
      gamma = sx
      if (lgw) then
        xx = (vol/nbas)/135d0  ! Volume relative to Si (as a reference)
        if (gcutb == NULLI) gcutb = 2.7d0/xx**(1d0/3d0)
        if (gcutx == NULLI) gcutx = 2.2d0/xx**(1d0/3d0)
        if (gwemax == NULLI) gwemax = max(2d0,1.9d0/xx**(2d0/3d0))
      endif

C --- Default values for algebraic variables used in later expressions ---
      if (express <= 4) then
        if (iprc >= 1) then
          write(ifout,"(/'# Variables entering into expressions parsed by input')")
        endif
        i = 1; if (cmdoptswx('--findes','float','') > 0) i = 11
        call awrit5('%% const nit=%i%?;n; asa=0;;%?;n; asa=t;;%?;n; les=%i;%j;',
     .    ' ',120,ifout,nit,fpandasa,isw(lasa)-fpandasa,isw(addes),i)
        call awrit3('%% const %?;n;%jmet=(asa?1:5);met=%i;'//
     .    '%?;n; loptic=0 mefac=1;;',
     .    ' ',120,ifout,isw(lasa),lmet,isw(loptics))
        if (nsp == 2) then
          call awrit1('%% const nsp=%i so=0',' ',120,ifout,nsp)
        elseif (express <= 3) then
          call awrit1('%% const so=0 nsp=so?2:1',' ',120,ifout,nsp)
        endif
        if (lpbesol) then
          call awrit1('%% const lxcf=0 lxcf1=116 lxcf2=133'//
     .      '%?;(n>1);%35p# PBESOL functional;;',' ',120,ifout,iprc)
        elseif (lpbe) then
          call awrit1('%% const lxcf=0 lxcf1=101 lxcf2=130'//
     .      '%?;(n>1);%35p# PBE functional;;',' ',120,ifout,iprc)
        else
          call awrit1('%% const lxcf=2 lxcf1=0 lxcf2=0'//
     .      '%?;(n>1);%35p# for PBE use: lxcf=0 lxcf1=101 lxcf2=130;;',' ',120,ifout,iprc)
        endif
        if (isw(lasa)-fpandasa>0 .and. lfrzw) then
          call awrit0('%% const frzwf=f%35p# Set to T to freeze log derivatives',' ',120,ifout)
        elseif (lfrzw) then
          call awrit2('%% const alat=%d dalat=0 frzwf=f'//
     .      '%?@(n>1)@ # Use dalat to change volume',' ',120,ifout,alat,iprc)
        endif
        if (express <= 3 .and. lfp) then
          call awrit3('%% const pwmode=%i pwemax=%d'//
     .    '%?;(n>1);%35p# Use pwmode=1 or 11 to add APWs',' ',120,ifout,pwmode,pwemax,iprc)
        endif
        if (lgw) then
          call awrit5('%% const sig=%i gwemax=%,1;2d gcutb=%,1;1d gcutx=%,1;1d'//
     .      '%?;n;  # GW-specific;;',' ',
     .      120,ifout,gwsig,gwemax,gcutb,gcutx,iprc)
        endif
        if (lasa .and. express <= 4) then
          call awrit4('%% const ccor=%l sx=%i gamma=sx scr=%i'//
     .      '%?;(n>1);%35p# ASA-specific;;',
     .      ' ',120,ifout,ccor,sx,scr,iprc)
        endif
        if (lasa .and. rmaxs /= NULLI) then
          call awrit1('%% const rmaxs=%d%35p# Range for structure constants',
     .      ' ',120,ifout,rmaxs)
        endif
        if (lgf) then
          call awrit3('%% const gfmode=1 nz=%i ef=%d c3=t'//
     .      '%?;(n>1);%35p# lmgf-specific variables;;',' ',120,ifout,nz,ef,iprc)
        endif

C   ... Miscellaneous variables
        if (express <= 4) then
          strn = '% const'

          if (nkabclm(0) <= 1) then
            call awrit2('%a nkabc=%?;(n==0);0;%i;',strn,len(strn),0,nkabclm(0),nkabclm(1))
          elseif (nkabclm(0) == 2 .and. nkabclm(1) == nkabclm(2)) then
            call awrit2('%a nk1=%i nk2=nk1',strn,len(strn),0,nkabclm(1),nkabclm(2))
          elseif (nkabclm(0) == 2) then
            call awrit2('%a nk1=%i nk2=%i',strn,len(strn),0,nkabclm(1),nkabclm(2))
          elseif (nkabclm(0) == 3 .and. nkabclm(1) == nkabclm(2) .and. nkabclm(1) == nkabclm(3)) then
            call awrit3('%a nk1=%i nk2=nk1 nk3=nk2',strn,len(strn),0,nkabclm(1),nkabclm(2),nkabclm(3))
          elseif (nkabclm(0) == 3 .and. nkabclm(1) == nkabclm(2)) then
            call awrit2('%a nk1=%i nk2=nk1 nk3=%i',strn,len(strn),0,nkabclm(1),nkabclm(3))
          elseif (nkabclm(0) == 3 .and. nkabclm(2) == nkabclm(3)) then
            call awrit3('%a nk1=%i nk2=%i nk3=nk2',strn,len(strn),0,nkabclm(1),nkabclm(2),nkabclm(3))
          elseif (nkabclm(0) == 3) then
            call awrit3('%a nk1=%i nk2=%i nk3=%i',strn,len(strn),0,nkabclm(1),nkabclm(2),nkabclm(3))
          endif

          if (lgw .and. nkabcgw(0) <= 1) then
            call awrit2('%a nkgw=%?;n;nkabc;%i;',strn,len(strn),0,isw(nkabcgw(0)==0.and.nkabcgw(0)<=1),nkabcgw(1))
          elseif (lgw .and. nkabcgw(0) == 2) then
            call awrit2('%a nkgw1=%i nkgw2=%i',strn,len(strn),0,nkabcgw(1),nkabcgw(2))
          elseif (lgw .and. nkabcgw(0) == 3 .and. nkabcgw(1) == nkabcgw(2) .and. nkabcgw(1) == nkabcgw(3)) then
            call awrit1('%a nkgw1=%i nkgw2=nkgw1 nkgw3=nkgw2',strn,len(strn),0,nkabcgw(1))
          elseif (lgw .and. nkabcgw(0) == 3 .and. nkabcgw(1) == nkabcgw(2)) then
            call awrit2('%a nkgw1=%i nkgw2=nkgw1 nkgw3=%i',strn,len(strn),0,nkabcgw(1),nkabcgw(3))
          elseif (lgw .and. nkabcgw(0) == 3 .and. nkabcgw(2) == nkabcgw(3)) then
            call awrit2('%a nkgw1=%i nkgw2=%i nkgw3=nkgw2',strn,len(strn),0,nkabcgw(1),nkabcgw(2))
          elseif (lgw .and. nkabcgw(0) == 3) then
            call awrit3('%a nkgw1=%i nkgw2=%i nkgw3=%i',strn,len(strn),0,nkabcgw(1),nkabcgw(2),nkabcgw(3))
          endif

          if (lfp) then
            call awrit2('%a gmax=%?;n;0;%d;',strn,len(strn),0,isw(gmax<=0),gmax)
          endif
          if (express <= 3) then
            call awrit0('%a beta=.3',strn,len(strn),0)
          endif
          if (lpmt) then
            call awrit0('%a '//'kmxa={pwemax>2?flor(2+pwemax):4}',strn,len(strn),0)
          endif
          if (varstrn /= ' ') then
            call awrit0('%a '//trim(varstrn),strn,len(strn),0)
          endif
          write(ifout,"(a)") trim(strn)
         endif

         if (lconvp>0) then
           xx = dhftol; if (xx == NULLI) xx = 0
           call awrit3('%% const conv=%g convc=%g dhftol=%g',' ',120,ifout,conv,convc,xx)
         endif

        write(ifout,"('')")
      endif

C --- VERS ---
      if (.not. lfp) then
        call awrit1('VERS%6pLM:7 ASA:7%?;n; # FP:7;;',' ',120,ifout,isw(express <= 4))
      else if (.not. lasa) then
        call awrit1('VERS%6pLM:7 FP:7%?;n; # ASA:7;;',' ',120,ifout,isw(express <= 4))
      else
        call awrit1('VERS%6pLM:7 FP:7%?;n; ASA:7;;',' ',120,ifout,isw(express <= 4))
      endif

      if (express <= 5) then
      call awrit1('IO%6pSHOW=f HELP=f IACTIV=f VERBOS=35%?;(n<=4);,35  OUTPUT=*;;',
     .  ' ',120,ifout,express)
      endif

C --- EXPRESS ---
      if (express > 0) then
        write(ifout,"('EXPRESS')")
        if (iprc > 0) write(ifout,"('# Lattice vectors and site positions')")
        if (cmdopt('--nfile',7,0,outs)) then
          write(ifout,358)
  358     format('%ifdef file'/'  file=   site{file}'/'%endif')
        endif

        write(ifout,"('  file=   site')")
        if (addes) write(ifout,"('# file= essite')")

C   ... Basis set
        if (lfp .and. express > 1) then
          if (iprc > 0) write(ifout,"(/'# Basis set')")
          call awrit4('%?;n;#; ; gmax=%10p'//
     .      '%?;(n<=4);{gmax}%j;%d;'//
     .      '%?;(n>1);%35p# PW cutoff for charge density;;',
     .      ' ',120,ifout,isw(.not.lfp .or. lfrzw),express,gmax,iprc)
          ivec(1:2) = [1,1]; if (lgw) ivec(1) = 2; if (v6float) ivec(2) = 0
          call awrit5('  autobas[pnu=1 loc=%i%?;n; eloc=%d;%j; mto=%i lmto=%?;n;5;4;%-1j%?;n; gw=%-1j%i',
     .      strn,len(strn),0,loc,isw(loc>0.and.eloc/=NULLI),eloc,mto,isw(lgw))
          call awrit6('%a%?;n; rsmmx=%d;%j;%?;n; ehmx=%d;%j;%?;n; pfloat=%s,%2i;%j;]',
     .      strn,len(strn),-ifout,
     .      isw(rsmmx /= NULLI),rsmmx,
     .      isw(ehmx /= NULLI),ehmx,
     .      isw(v6float.or.lgw),ivec)

        endif

C   ... Self-consistency
        if (iprc > 0) write(ifout,"(/'# Self-consistency')")
        if (lasa) then
        call awrit3(
     .    '  nit=    %?;(n<=4);{abs(nit)}%j;%i;'//
     .    '%?;(n>1);%35p# Maximum number of iterations;;',
     .    ' ',120,ifout,express,nit,iprc)
        else
        call awrit3(
     .    '  nit=    %?;(n<=4);{nit}%j;%i;'//
     .    '%?;(n>1);%35p# Maximum number of iterations;;',
     .    ' ',120,ifout,express,nit,iprc)
        endif
        if (express <= 3) then
          call awrit1('  mix=    B2,b={beta},k=7'//
     .      '%?;(n>1);%35p# Charge density mixing parameters;;',' ',120,ifout,iprc)
        else
          call awrit1('  mix=    B2,b=.3,k=7'//
     .      '%?;(n>1);%35p# Charge density mixing parameters;;',' ',120,ifout,iprc)
        endif
        if (express <= 6) then
          call awrit3('  conv=   %?#n#{conv}%j#%g#'//
     .      '%?;(n>1);%35p# Convergence tolerance (energy);;',' ',120,ifout,lconvp,conv,iprc)
        endif
        call awrit3('  convc=  %?#n#{convc}%j#%g#'//
     .    '%?;(n>1);%35p# tolerance in RMS (output-input) density',' ',120,ifout,lconvp,convc,iprc)
        if (dhftol /= NULLI) then
          call awrit3('  dhftol= %?#n#{dhftol}%j#%g#'//
     .      '%?;(n>1);%35p# tolerance in correction to Harris Forces',' ',120,ifout,lconvp,dhftol,iprc)
        endif

C   ... Brillouin zone
        if (iprc > 0) write(ifout,"(/'# Brillouin zone')")
        if (nkabclm(0) == 0) nkabclm(0) = 1
        if (express >= 5) then
          call awrit3('  nkabc= %n:1i%?;(n>1);%35p# 1 to 3 values;;',
     .    ' ',120,ifout,nkabclm,nkabclm(1),iprc)
        else
          if (nkabclm(0) <= 1) then
            call awrit1('  nkabc=  {nkabc}%?;(n>1);%35p# 1 to 3 values;;',' ',120,ifout,iprc)
          elseif (nkabclm(0) == 2) then
            call awrit1('  nkabc=  {nk1},{nk2}%?;(n>1);%35p# 1 to 3 values;;',' ',120,ifout,iprc)
          elseif (nkabclm(0) == 3) then
            call awrit1('  nkabc=  {nk1},{nk2},{nk3}%?;(n>1);%35p# 1 to 3 values;;',' ',120,ifout,iprc)
          endif
        endif
        if (express <= 6) then
        call awrit3('  metal=  %?;(n<=4);{met}%j;%i;'//
     .    '%?;(n>1);%35p# Management of k-point integration weights in metals;;',
     .    ' ',120,ifout,express,lmet,iprc)
        endif

C   ... Potential
        lwtag = nsp == 2 .or. express <= 6
        if (lwtag) then
        if (iprc > 0) write(ifout,"(/'# Potential')")
        if (nsp == 2 .or. express <= 3) then
          call awrit3('  nspin=  %?;(n<=4);{nsp}%j;%i;'//
     .    '%?;(n>1);%35p# 2 for spin polarized calculations;;',
     .      ' ',120,ifout,express,nsp,iprc)
          call awrit3('  so=     %?;(n<=4);{so}%j;%i;'//
     .    '%?;(n>1);%35p# 1 turns on spin-orbit coupling;;',
     .      ' ',120,ifout,express,max(lso,0),iprc)
        endif

        if (express <= 4) then
          call awrit1('%2pxcfun=  {lxcf},{lxcf1},{lxcf2}'//
     .      '%?;(n>1);%35p# set lxcf=0 for libxc functionals',' ',120,ifout,iprc)

        elseif (express <= 6 .and. lpbe) then
          call awrit1('%2pxcfun=  0,101,130'//
     .      '%?;(n>1);%35p# Use 0,x,x for libxc, e.g. 0,101,130 for PBE',
     .      ' ',120,ifout,iprc)
        elseif (express <= 6) then
          call awrit1('%2pxcfun=  2'//
     .      '%?;(n>1);%35p# Use 0,x,x for libxc, e.g. 0,101,130 for PBE',
     .      ' ',120,ifout,iprc)
        endif
        endif

C   ... End of EXPRESS
        write(ifout,"('')")

      endif


C --- SYMGRP ---
C     Not needed, but written for informational purposes
      if (express <= 5) then
      k = index(aa,'i*i'); if (k /= 0) aa(k:k+2) = ' '
      call skpblb(aa,len(aa),k)
      k = k+1
      if (k > 0) then
        i = 0
        call skipbl(aa,k,i)
        write(ifout,'(''#SYMGRP '',a)') aa(i+1:k)
      else
        write(ifout,'(''SYMGRP  find'')')
      endif
      endif

C --- HAM ---
      if (lasa .and. cmdopt('--lmfit',7,0,strn)) write(ifout,"('% const rdvext=0')")

      lwtag = (lfp .and. express <= 5) .or. (lasa .and. express <= 2) .or. lgw
     .        .or. lrel > 0 .or. (lfrzw .and. lpfree) .or. cmdopt('--lmfit',7,0,strn)
      if (lfp .and. lfrzw) then
        call awrit3('HAM   GMAXDA='//
     .    '%?;(n<=4);{gmax}%j;%d;'//
     .    '%?;(n>1);%35p# PW cutoff for charge density;;'//
     .    '%N%6fFRZWF={frzwf}%16f# Freeze phi,dot',
     .    ' ',120,ifout,express,gmax,iprc)
        lwtag = .false.
      elseif (express <= 2 .and. lfp) then
        call awrit4('HAM   GMAX%?;n;DA=;=;'//
     .    '%?;(n<=4);{gmax}%j;%d;'//
     .    '%?;(n>1);%35p# PW cutoff for charge density;;',
     .    ' ',120,ifout,isw(lfrzw),express,gmax,iprc)
        lwtag = .false.
      elseif (lasa .and. (lfrzw .and. .not. lpfree .or. .not. lfrzw .and. lpfree)) then
        if (lfrzw)  call awrit0('HAM   FRZWF={frzwf}%35p# Freeze all Pnu',' ',120,ifout)
        if (lpfree) call awrit0('HAM   PMIN=-1%35p# Constrain Pnu > P_free',' ',120,ifout)
        lwtag = .false.
      else if (lwtag) then
        write(ifout,"('HAM')")
      endif

      if (lwtag .and. lasa) then
        if (lfrzw)  call awrit0('      FRZWF={frzwf}%35p# Freeze all Pnu',' ',120,ifout)
        if (lpfree) call awrit0('      PMIN=-1%35p# Constrain Pnu > P_free',' ',120,ifout)
      endif

      if (lasa .and. cmdopt('--lmfit',7,0,strn)) write(ifout,"(6x,'RDVEXT={rdvext}')")

      if (lrel > 0) then
        call stswt('lunit',ifout,swt)
        k = rdtok('REL=',instr,'%i,,6;;Specifies relativistic treatment',
     .    ' ',', ',2,swt,is,1,lrel,cxx,prs)
      endif

      if (express <= 2 .and. lfp) then  ! tags duplicate those in EXPRESS
        ivec(1:2) = [1,1]; if (lgw) ivec(1) = 2; if (v6float) ivec(2) = 0
        call awrit5('%x%6fAUTOBAS[PNU=1 LOC=%i%?;n; ELOC=%d;%j; MTO=%i LMTO=%?;n;5;4;%-1j GW=%i',
     .      strn,len(strn),0,loc,isw(loc>0.and.eloc/=NULLI),eloc,mto,isw(lgw))
        call awrit6('%a%?;n; RSMMX=%d;%j;%?;n; EHMX=%d;%j;%?;n; PFLOAT=%s,%2i;%j;]',
     .      strn,len(strn),-ifout,
     .      isw(rsmmx /= NULLI),rsmmx,
     .      isw(ehmx /= NULLI),ehmx,
     .      isw(v6float.or.lgw),ivec)
      endif

      if (express <= 3 .and. lfp .and. lfrzw) then
        call awrit0('%6pPWMODE={pwmode} PWEMIN=0 PWEMAX={pwemax}/(1+{dalat/alat})^2 OVEPS=0',
     .    ' ',120,ifout)
      elseif (express <= 3 .and. lfp) then
        call awrit6('%?;(n<5); ;#;%-1j%6pPWMODE=%?;(n<=4);{pwmode}%j;%i; PWEMIN=0 '//
     .    'PWEMAX=%?;(n<=4);{pwemax}%j;%d; OVEPS=0'//
     .    '%?;(n>1); # For APW addition to basis;;',
     .    ' ',120,ifout,express,pwmode,express,pwemax,express,iprc)
      endif

      if (lpmt) then
        call awrit0('%6pTOL=1d-16',' ',120,ifout)
      endif

      if (express <= 2) then ! tags duplicate those in EXPRESS
        write(ifout,"(6x,'XCFUN={lxcf},{lxcf1},{lxcf2}':t35,a)")
      endif

      if (express <= 3 .and. lfp) then
        call awrit1('%6fFORCES={so==0} ELIND=-0.7 '//
     .  '%?;(n<=2);NSPIN={nsp} SO={so} ;',' ',120,ifout,express)
      elseif (express <= 5 .and. lfp) then
        call awrit1('%6fFORCES=%i ELIND=-0.7 ',' ',120,ifout,isw(lso<=0))
      endif

      if (lgw) then
        if (express <= 4) then
          call awrit1('%6fRDSIG={sig} SIGP[EMAX={gwemax}]'//
     .      '%?;(n>1);%2f# Add self-energy to LDA;;',
     .      ' ',120,ifout,iprc)
        else
          call awrit3('%6fRDSIG=%i SIGP[EMAX=%d]'//
     .      '%?;(n>1);%35p# Add self-energy to LDA;;',
     .      ' ',120,ifout,gwsig,gwemax,iprc)
        endif
      endif

C --- OPTICS ---
      if (loptics) then
        call awrit1('OPTICS%N%6f'//
     .    'MODE=%?;(n<=4);{loptic};0; NPTS=1001 WINDOW=0,1 MEFAC=%-1j%?;(n<=4);{mefac};1; LTET=3',
     .    ' ',120,ifout,express)
      endif

C --- OPTIONS ---
      if (lasa) then
        if (iprc >= 2) write(ifout,"('')")
        write(ifout,"('OPTIONS')")
        if (express <= 4) then
          write(ifout,106)
  106     format(6x,'ASA[ CCOR={ccor} TWOC=F ADNF=F GAMMA={gamma}]',' SCR={scr} SX={sx}')
        else
          call awrit5('%6fASA[ CCOR=%l SX=%i TWOC=F ADNF=%l GAMMA=%i SCR=%i ]',
     .      ' ',120,ifout,ccor,sx,adnf,gamma,scr)
        endif
      endif

C --- STR ---
      if (lasa .and. rmaxs /= NULLI) then
        if (rmaxs /= NULLI)  call awrit0('STR   RMAXS={rmaxs}%35p# Range of strux',' ',120,ifout)
      endif

C --- ITER ---
      if (cmdopt('--lmfit',7,0,strn)) then
        write(ifout,"('% const lmfit=0')")
      endif
      if (express <= 2) then
        write(ifout,"('ITER  ',a,'  NIT={nit}  CONVC=1e-5')") 'MIX=B2,b={beta},k=7'
      endif
      if (cmdopt('--lmfit',7,0,strn)) then
        if (express > 2) write(ifout,"('ITER')")
        write(ifout,"(6x,'FIT[MODE={lmfit} NBFIT=1 NBFITF=1 LAM=1 SCL=5 SHFT=1 WT=1,.2,.1]')")
      endif

C --- EWALD ---
      if (cmdopt('--ewald=',8,0,strn)) then
        i = 8; is = a2vec(strn,len(strn),i,4,', ',2,2,1,it,xx)
        if (is /= 1) call rx('blm failed to parse '//trim(strn))
        call awrit1('EWALD TOL=%g',' ',120,ifout,xx)
      elseif (lpmt) then
        write(ifout,"('EWALD TOL=1d-16')")
      endif

C --- GF ---
      if (lgf) then
        if (iprc >= 2) write(ifout,"('')")
        call awrit3('GF%6pMODE=%?;(n<=4);{gfmode}%j;%i; '//
     .               'GFOPTS=%?@(n<=4)@{?~c3~p3;~p2;}@p3@',
     .    ' ',120,ifout,express,gfmode,express)
      endif

C --- BZ ---
      lwtag = fsmom /= NULLI
      if (express <= 2) then
        strn = 'BZ    METAL={met}'
        if (nkabclm(0) <= 1) then
          call awrit0('%a  NKABC={nkabc}',strn,120,-ifout)
        elseif (nkabclm(0) == 2) then
          call awrit0('%a  NKABC={nk1},{nk2}',strn,120,-ifout)
        elseif (nkabclm(0) == 3) then
          call awrit0('%a  NKABC={nk1},{nk2},{nk3}',strn,120,-ifout)
        endif
        lwtag = .false.
      endif
      if (lgf) then
        call awrit5('%?;(n>2);BZ;;%-1j'//
     .    '%6pEMESH=%?;(n<=4);{nz}%j;%i;,10,-1,%?;(n<=4);{ef}%j;%d;,.5,.3'//
     .    ' %?@n@%35p# nz-pts;contour mode;emin;emax;ecc;bunching @@',
     .    ' ',120,ifout,express,nz,express,ef,iprc)
        lwtag = .false.
      endif
      if (lwtag) then
        write(ifout,"('BZ')")
      endif
      if (fsmom /= NULLI) then
        call awrit1('%6fFSMOM=%g',' ',120,ifout,fsmom)
      endif

C     write(ifout,"('#',7x,'INVIT=f',a)")
C    .  ' # slows execution slightly.  A bit slower but more reliable'

C ... Get basis lmx. If it has an f orbital possibly downfold it
      nl = 4
      if (cmdopt('--nl=',5,0,outs)) then
        i = 5
        is = a2vec(outs,len_trim(outs),i,2,', ',2,2,1,it,nl)
        if (is /= 1) call rx('blm failed to parse '//trim(outs))
      endif
      do  is = 1, nspec
        j1 = 0
        if (lgw) j1 = 5
        call deflx(3,j1,z(is),lmxb,lmxa)
        if (z(is) == 0 .or. (z(is) <=3 .and. lmxb(1) < rmt(is))) then
          lmxb(1) = rmt(is)*(2d0/1.8d0)
          lmxa    = min(lmxb(1)+1,2)
          if (z(is) > 0) lmxa = lmxb(1)+1
        endif
        if (lmxbs(is) > NULLI) lmxb(1) = lmxbs(is)
        if (lasa) then
          if (cmdopt('--nl=',5,0,outs)) then
            lmxb(1) = min(nl-1,lmxb(1))
          endif
          nl = max(nl,lmxb(1)+1)
          if (lmxb(1) == 3) then
            pl(:) = 0; ql(:,1) = 0
            call defpq(1,z(is)+0,lmxb,1,pl,ql)
            if (ql(4,1) == 0) idxdn(4,is) = -1
            pl(:) = 0; ql(:,1) = 0
            call defpq(1,z(is),lmxb,1,pl,ql)
            if (ql(4,1) == 0) idxdn(4,is) = idxdn(4,is)-1
            if (idxdn(4,is) == -2) idxdn(4,is) = 2
          endif
        else
          nl = max(nl,lmxa+1)
        endif
      enddo
      if (lgw) nl = max(nl,5)

C --- GW ---
      if (lgw) then
        strn = 'GW    NKABC='
        if (express >= 5) then
          call awrit2('%a%n:1i',strn,120,0,nkabcgw,nkabcgw(1))
        elseif (nkabcgw(0) <= 1) then
          call awrit0('%a{nkgw}',strn,120,0)
        elseif (nkabcgw(0) == 2) then
          call awrit0('%a{nkgw1},{nkgw2}',strn,120,0)
        elseif (nkabcgw(0) == 3) then
          call awrit0('%a{nkgw1},{nkgw2},{nkgw3}',strn,120,0)
        endif
        call awrit4(
     .    '%a GCUTB=%?@(n<=4)@{gcutb}%j@%,2;2d@ GCUTX=%?@(n<=4)@{gcutx}%j@%,2;2d@ DELRE=.01 .1'//
     .      '%N%6fGSMEAR=0.003 PBTOL='//trim(pbtol),
     .    strn,len(strn),-ifout,express,gcutb,express,gcutx)
      endif

C --- STRUC ---
      if (express >= 3 .and. lfp .and. lfrzw) then
        call awrit0('STRUC DALAT={dalat}',' ',120,ifout)
      endif
      if (express >= 3) goto 70                 ! Skip STRUC
      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'STRUC ','Specifies the lattice vectors basis size',
     .  swt,.false.,mxrec,recl0,instr)

      if (addes) then
        call awrit5('      NBAS=%i+{les?%i:0}  NL=%i  NSPEC=%i+{les?%i:0}',' ',120,ifout,
     .    nbas0,nbas-nbas0,nl,nspec0,nspec-nspec0)
      else
        call stswt('lunit',0,swt)
        k = rdtok('    NL=',instr,'default value of 1+l for basis and augmentation',
     .    ' ',' ',2,swt,1,1,nl,cxx,prs)
        k = rdtok('NBAS=',instr,'%i,,number of atoms',' ',' ',2,swt,1,1,nbas,cxx,prs)
        call stswt('lunit',ifout,swt)
        k = rdtok('NSPEC=',instr,'number of species',' ',' ',2,swt,1,1,nspec,cxx,prs)
      endif

      if (cmdopt('--wsitex',8,0,strn) .or. cmdopt('--wsite',7,0,strn)) then
        write(ifout,110)
  110   format('      FILE=site       # supersede NBAS, ALAT, and PLAT in ctrl file')
        if (addes) then
        write(ifout,113)
  113   format('#     FILE=essite     # use this file after running lmchk --findes --wsite')
        endif
      endif

      call stswt('lunit',0,swt)
      call awrit1('      ALAT= %;8,1d',' ',120,ifout,alat)
      if (lfrzw) call awrit0('      DALAT={dalat}',' ',120,ifout)

C      k = rdtok('ALAT=',instr,
C     .  '%;8g,,scaling of lattice vectors, in a.u.',
C     .  ' ',', ',4,swt,1,1,alat,cxx,prs)
C      k = rdtok('DALAT=',instr,
C     .  '%;8g,,is added to alat after input is read',
C     .  ' ',', ',4,swt,1,1,0d0,cxx,prs)
      call stswt('lunit',ifout,swt)
      k = rdtok('    PLAT=',instr,' %;7d,,Primitive lattice vectors (dimensionless)',
     .  ' ',' ,',4,swt,1,9,plat,cxx,prs)

   70 continue

C --- SPEC category ---
      lsmalla = 0; if (lpbe .or. lrel == 2) lsmalla = 1; if (lrel > 10) lsmalla = 2
      if (iprc >= 2) write(ifout,"('')")
      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'SPEC ',
     .  'Specifies additional species information '//
     .  '%N%3f(overriding default values or possible input from '//
     .  'SITE category)'//
     .  '%N%3fLabels should coincide with labels in SITE category%N',
     .  swt,.true.,mxrec,recl0,instr)

      call awrit3('%,2d %,2d %,2d',outs,30,0,omax1(1),omax1(2),omax1(3))
      if (lasa .and. addes) then
        write(ifout,"('# SCLWSR=11 WSRMAX=3.3 OMAX1=',a)") trim(outs)
      elseif (addes) then
        write(ifout,"('# SCLWSR=21 WSRMAX=3.3')")
      elseif (lasa) then
        write(ifout,"('# SCLWSR=11 WSRMAX=3.3 OMAX1=',a)") trim(outs)
      endif
      call stswt('token,opt,fixlen,mult',0,swt)
      catlen = swt(2)
      swt(2) = 1
      is = 0
      kount = 0
  130 continue
        if (lout) then
          if (kount >= nspec) goto 132
          alabl = slabl(is+1)
        else
          swt(1) = swt(2)
          swt(2) = catlen
        endif
        call stswt('lunit',0,swt)
        k = rdtok('ATOM=',instr,
     .    '%16p,,Species label.  Enter one for each species.',
     .    ' ',', ',1,swt,1,1,xx,alabl,prs)

C       No further labels found -> terminates infinite loop
        if (k == 0) goto 132
        call info(helppr,0,0,'%5f... Following each label, the'//
     .    ' following may be supplied:',0,0)

C       If this species missing from species list, do
        if (.not.linfo) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',is,k,.false.)
          if (is < 0) then
            call info(infopr-20,0,0,'%N%5f(warning) ... label `'
     .        //alabl//'%a'' does not match any sites ... '//
     .        'species ignored',0,0)
            goto 130
          else
            is = is+1
          endif
        endif
        kount = kount+1

C       Tokens for this species
        k = rdtok('Z=',instr,'%;4,0D%b %b,,0;;Atomic number.  If not '//
     .    'supplied, Z is inferred from the label',' ',', ',4,swt,
     .    is,1,z,cxx,prs)

        strn = '%;9,6D,,4'
        if (express <= 4 .and. lfp .and. lfindes) then
          strn = '%;9,6D                ,,4'
          if (z(is) == 0) strn = '%;9,6D*{(les==11)?0:1},,4'
        endif
        k = rdtok('R=',instr,strn,' ',', ',4,swt,is,1,rmt,cxx,prs)
        j1 = 0
        if (lgw) j1 = 5
        call deflx(3,j1,z(is),lmxb,lmxa)
        if (z(is) == 0 .or. (z(is) <=3 .and. lmxb(1) < rmt(is))) then
          lmxb(1) = rmt(is)*(2d0/1.8d0)
          lmxa    = min(lmxb(1)+1,2)
          if (z(is) > 0) lmxa = lmxb(1)+1
        endif
        if (lmxbs(is) > NULLI) lmxb(1) = lmxbs(is)
        if (cmdopt('--nl=',5,0,outs)) then
          lmxb(1) = min(nl-1,lmxb(1))
        endif

C       This string is complex ... print out in suctrl.f
C       llast = lamom(is) == 0 .and. idxdn(4,is) == 0 .and. .not. lpmt
        if (isw(lasa)-fpandasa>0) then
C         if (llast) call stswt('lunit',ifout,swt)
          k = rdtok('LMX=',instr,'%i,,4;;l-cutoff for basis',
     .      ' ',', ',4,swt,1,1,lmxb,cxx,prs)
        else
          i = lmxb(1); ! if (lpmt) i = min(lmxb(1),2)
          k = rdtok('LMX=',instr,'%i,,4;;l-cutoff for basis',
     .      ' ',', ',4,swt,1,1,i,cxx,prs)
C         if (llast) call stswt('lunit',ifout,swt)
          if (lgw) lmxa = nl-1
C          k = rdtok('LMXA=',instr,'%i,,4;;l-cutoff for augmentation',
C     .      ' ',', ',4,swt,1,1,lmxa,cxx,prs)
          call awrit1('%a LMXA=%i',prs,len(prs),0,lmxa)
          if (lpmt) then
C           llast = lamom(is) == 0 .and. idxdn(4,is) == 0
C           if (llast) call stswt('lunit',ifout,swt)
            if (express <= 4) then
              call awrit0('%a KMXA={kmxa}',prs,len(prs),0)
            else
              i = max(4,int(pwemax-2)+4)
              k = rdtok('KMXA=',instr,'%i,,4;;k-cutoff for augmentation',
     .          ' ',', ',4,swt,1,1,i,cxx,prs)
            endif
          endif
        endif

        if (idxdn(4,is) /= 0) then
C         llast = lamom(is) == 0
C         if (llast) call stswt('lunit',ifout,swt)
C         k = rdtok('IDXDN=',instr,'%i,,4;;',' ',' ,',2,swt,1,4,idxdn(1,is),cxx,prs)
          call awrit2('%a  %s,IDXDN=%ni',prs,len(prs),0,lmxb(1)+1,idxdn(1,is))
        endif

        if (lamom(is) /= 0) then
C         llast =  .not. lpmt
C         if (llast) call stswt('lunit',ifout,swt)
C          k = rdtok('MMOM=',instr,'%;6d,,0',' ',' ,',4,
C     .      swt,1,lamom(is),amom(1,is),cxx,prs)
          call awrit2('%a %s,MMOM=%nd',prs,len(prs),0,lamom(is),amom(1,is))
        endif

C   ... Automatically add PZ=...,n+1.4
        if (lgw .and. mod(adloc(is)/100,10) /= 0 .and. pz(3,is) == 0) then
          pz(3,is) = mod(adloc(is)/100,10) + 1.4d0
          lpz(is) = max(lpz(is),3)
        endif

C   ... Automatically add P=...,n+0.92
        if (lgw .and. mod(adloc(is)/1000000,10) /= 0 .and. pnu(3,is) == 0) then
          pnu(3,is) = mod(adloc(is)/1000000,10) + 0.92d0
          lpnu(is) = max(lpnu(is),3)
        endif

        if (lpnu(is) /= 0) then
          call awrit2('%a %s,P=%nd',prs,len(prs),0,lpnu(is),pnu(1,is))
        endif

        if (lpz(is) /= 0) then
          call awrit2('%a %s,PZ=%nd',prs,len(prs),0,lpz(is),pz(1,is))
        endif

        if (lrsmh(is) /= 0) then
          call awrit2('%a %s,RSMH=%nd',prs,len(prs),0,lrsmh(is),rsmh(1,is))
        endif

        if (lidmod(is) /= 0) then
          call awrit2('%a %s,IDMOD=%ni',prs,len(prs),0,lidmod(is),idmod(1,is))
        endif

        if (lsmalla /= 0) then
          xx = .015d0
          if (lsmalla > 1) xx = .01d0
          call awrit1('%a %s,A=%d',prs,len(prs),0,xx)
        endif

        if (spectok(is) /= ' ') then
          call awrit0('%a '//trim(spectok(is)),prs,len(prs),0)
        endif

C       Write the string to ifout
        k = rdtok(' ',' ',' ',' ',' ,',0,[0,0,0,0,4,ifout,0],0,0,[0],cxx,prs)

        if (lldau(is) /= 0) then
          write(ifout,"('% ifdef ldau')")
          call ivset(iv,1,4,0)
          do  i = 1, 4
            if (abs(uh(i,1,is))+abs(uh(i,2,is)) /= 0) iv(i) = 12
          enddo
          call awrit6('%16f%s,IDU=%ni  UH=%nd  JH=%nd',' ',80,ifout,
     .      lldau(is),iv,lldau(is),uh(1,1,is),lldau(is),uh(1,2,is))
          write(ifout,"('% endif')")
        endif

        if (specincl(is) /= ' ') then
          print *, 'need to write file here:', specincl(is)
        endif

        if (.not.linfo) goto 130

C ... End of loop
  132 continue

C     No ES species yet.. Add one but comment it out
      if (lesend .and. mxesspec == NULLI) then
        if (lgw) then
          call awrit1('# ATOM=E        Z=  0  R= 2*{(les==11)?0:1}  LMX=2 LMXA=%i',' ',80,ifout,nl-1)
        else
          call awrit0('# ATOM=E        Z=  0  R= 2*{(les==11)?0:1}  LMX=2',' ',80,ifout)
        endif
      endif

C --- SITE category ---
      i = 0;if (express > 0) i = 15001
      if (cmdopt('--wsite',7,0,strn)) then
         i = 15001
         dc = strn(8:8)
       endif
      if (cmdopt('--wsitex',8,0,strn)) then
        i = 15011
        dc = strn(9:9)
      endif

      if (i /= 0) then
        fn = 'site' ! Default file name

C   ... Shift sites to first quadrant
        if (wordsw(strn,dc,'quad1','',j1) > 0) then
          call info0(30,0,0,' ... shifting sites to first quadrant')
          do  k = 1, nbas
C           posp = posc (plat+)^-1  and  posp+ = (plat)^-1 posc+
            call dgemm('T','N',3,1,3,1d0,qlat,3,pos(1,k),3,0d0,xv,3)
            ltmp = .true.
            do while (ltmp)
              ltmp = xv(1)<-fptol .or. xv(2)<-fptol .or. xv(3)<-fptol
              if (xv(1) < -fptol) xv(1) = xv(1)+1
              if (xv(2) < -fptol) xv(2) = xv(2)+1
              if (xv(3) < -fptol) xv(3) = xv(3)+1
            enddo
C           posc = posp (plat+) and  posc+ = (plat) posp+
            call dgemm('N','N',3,1,3,1d0,plat,3,xv,3,0d0,pos(1,k),3)
          enddo
        endif

C   ... Set output site file name
        if (wordsw(strn,dc,'fn','= ',j1) /= 0) then
          fn = strn(j1+1:)
        endif

        j1 = iosite(i,3d0,0,trim(fn),i,slabl,alat,plat,nbas,
     .              nspec,pos,xx,xx,xx,ips,xx,xx)
      endif

      if (express >= 3) goto 80 ! Skip SITE
      if (cmdopt('--wsite',7,0,strn) .or. cmdopt('--wsitex',8,0,strn)) then
        write(ifout,108)
  108   format('SITE  FILE=site')
        if (addes) then
          write(ifout,112)
  112     format('#SITE FILE=essite')
        endif
      elseif (cmdopt('--wsite',7,0,strn)) then
        write(ifout,108)
      endif


      call stswt('cat,reqd,mxrec',0,swt); call stswt('lunit',ifout,swt)
      k = fgtcat(ifin,'SITE ',
     .  'Specifies basis.  '//
     .  'At least one must be specified.%N%3f'//
     .  'Labels must coincide with specifications'//
     .  ' in SPEC category.%N',
     .  swt,.true.,mxrec,recl0,instr)

      ib = 1
      catlen = swt(2)
      swt(2) = 1
      is = 0
  140 continue
        if (lout) then
          if (ib > nbas) goto 142
          is = ips(ib)
          alabl = slabl(is)
        else
          swt(1) = swt(2)
          swt(2) = catlen
        endif
        call stswt('lunit',0,swt)
        call stswt('token,opt,fixlen,mult',0,swt)
        k = rdtok('ATOM=',instr,'%16p,,Species label.',
     .    ' ',', ',1,swt,1,1,xx,alabl,prs)
C   ... No further labels found -> terminates infinite loop
        if (k == 0) goto 142
C   ... Add new species, if this label is new
        if (.not.linfo .and. .not.lout) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',i,k,.false.)
          if (i < 0) then
             i = nspec
             nspec = nspec+1
             slabl(nspec) = alabl
           endif
          ips(ib) = i+1
        endif
        is = ips(ib)

        call info(helppr,0,0,'%5f... Following each label, the'//
     .    ' following must be supplied:',0,0)
        call stswt('token,reqd,fixlen',0,swt)
        call stswt('lunit',ifout,swt)
        if (cmdopt('--xpos',6,0,strn)) then
        k = rdtok('XPOS=',instr,'%,;11,7D,,4;;Site positions, units of '
     .    //'lattice vectors',' ',' ,',4,swt,1,3,xpos(1,ib),cxx,prs)
        else
        k = rdtok('POS=',instr,'%,;11,7D,,4;;Site positions, Cartesian '
     .    //'coordinates',' ',' ,',4,swt,1,3,pos(1,ib),cxx,prs)
        endif

        ib = ib+1
        if (ib > nbmx) call rxi(
     .    'too many sites. Increase nbmx in suctrl; now',nbmx)
        if (nspec > nspmx) call rxi(
     .    'too many species. Increase nspmx in suctrl; now',nspmx)
        if (.not.linfo) goto 140
C ... End of loop
  142 continue
      nbas = ib-1
   80 continue
      deallocate(instr)

C --- Molecular statics ---
      if (lfp .and. cmdoptswx('--molstat','','') == 1) then
        write(ifout,355)
  355   format('% ifdef minx'/
     .    'DYN     MSTAT[MODE=6 HESS=t XTOL=.001 GTOL=0 STEP=.010 NKILL=0] NIT={minx}'/
     .    '% endif')
      endif

C --- START category ---
      if (lasa) then
        if (iprc >= 2) write(ifout,"('')")
        if (express <= 4) then
          write(ifout,"('START CNTROL={nit<=0} BEGMOM={nit<=0}')")
        else
          write(ifout,"('START CNTROL=T BEGMOM=T')")
        endif
      endif

      return
C --- Error exit ---
  999 continue
      outs = prgnam // '%a : ' // errmsg
      call fexit(1,9,outs,0)

      end
