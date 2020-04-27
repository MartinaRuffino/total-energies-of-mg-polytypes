C#ifdefC FP
C      integer function iorsf(mode,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
C     .  fid,nbas,nat,nspec,qnu,nit,lbin,ifi)
C#else
      integer function iorsa(mode,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .  fid,nbas,nat,nspec,qnu,nit,lbin,ifi)
C#endif
C- I/O for charge density
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lrs
Cio    Passed to:  bcast_strx iinit mpibc1
C    ... ASA
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos class spec force vel pnu pz
Co     Stored:     pos pos0 force vel spec clabel pnu pz
C    ... FP
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos class spec force vel pnu pz v0 v1
Co     Stored:     pos pos0 force vel spec clabel pnu pz
Co     Allocated:  v0 v1
Cio    Elts passed:rhoc rho1 rho2
Cio    Passed to:  *
Ci     Elts read:  pos class spec force vel pnu pz v0 v1 bxc cpawt omg
Ci                 omgn domg gc gcorr j0 pdos rho1 rho2 rhoc rho1x
Ci                 rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl sighh
Ci                 sighk sigkk tauhh tauhk taukk pihh pihk pikk sighhx
Ci                 sighkx sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx
Ci                 thet
Co     Stored:     pos pos0 force vel spec clabel pnu pz bxc cpawt omg
Co                 omgn domg gc gcorr j0 pdos rho1 rho2 rhoc rho1x
Co                 rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl sighh
Co                 sighk sigkk tauhh tauhk taukk pihh pihk pikk sighhx
Co                 sighkx sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx
Co                 thet v0 v1
Co     Allocated:  v0 v1
Cio    Elts passed:rhoc rho1 rho2
Cio    Passed to:  bcast_strx
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name a nr rmt z lmxa lmxl kmxt p pz lfoca qc idmod
Ci                 rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
Ci                 chfa rsmfa rhoc
Co     Stored:     name lmxa a nr rmt z lmxl kmxt p pz qc lfoca rfoca
Co                 coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa rhoc
Co     Allocated:  rhoc
Cio    Elts passed:rhoc
Cio    Passed to:  bcast_strx
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc ag bgv cg cy dlv gv gvq indxcg ips0
Ci                 istab jcg kv igv igv2 kv2 pos qlv symgr
Co     Stored:     alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                 jcg kv igv igv2 kv2 pos qlv symgr
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bcast_strx iinit mpibc1
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  (FP) smrho
Ci                 (ASA) qbyl vesrmt
Co     Stored:     (FP) smrho
Co                 (ASA) ves v0
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bcast_strx iinit mpibc1
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Stored:     ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bcast_strx iinit mpibc1
C    ... FP
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos class spec force vel pnu pz v0 v1 rho1 rho2 rhoc
Ci                 rho1x rho2x rhocx bxc cpawt omg omgn domg gc gcorr
Ci                 j0 pdos qhhl qhkl qkkl eqhhl eqhkl eqkkl sighh sighk
Ci                 sigkk tauhh tauhk taukk pihh pihk pikk sighhx sighkx
Ci                 sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx thet
Co     Stored:     pos pos0 force vel spec clabel pnu pz bxc cpawt omg
Co                 omgn domg gc gcorr j0 pdos rho1 rho2 rhoc rho1x
Co                 rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl sighh
Co                 sighk sigkk tauhh tauhk taukk pihh pihk pikk sighhx
Co                 sighkx sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx
Co                 thet v0 v1
Co     Allocated:  v0 v1 rho1 rho2 rhoc rho1x rho2x rhocx
Cio    Elts passed:rho1 rho2 rhoc rho1x rho2x rhocx
Cio    Passed to:  dfratm pvsms2 bcast_strx
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name a nr rmt z lmxa lmxl kmxt p pz lfoca qc idmod
Ci                 rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
Ci                 chfa rsmfa rhoc
Co     Stored:     name lmxa a nr rmt z lmxl kmxt p pz qc lfoca rfoca
Co                 coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa rhoc
Co     Allocated:  rhoc
Cio    Elts passed:rhoc
Cio    Passed to:  dfratm pvsms2 bcast_strx
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc dist ag bgv cg cy dlv gv gvq indxcg
Ci                 ips0 istab jcg kv igv igv2 kv2 pos qlv symgr
Co     Stored:     alat plat nabc ag bgv cg cy dlv gv gvq indxcg ips0
Co                 istab jcg kv igv igv2 kv2 pos qlv symgr
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bcast_strx iinit mpibc1
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  ves aamom bxc cp ddpf dddpf ddpfr dlmwt dpf dpfr
Ci                 gibbs gma gmar grrme mad mxy palp papg pf pfnc pfr
Ci                 pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt rhos
Ci                 rhrmx sop thetcl vdif vintr vrmax vshft smpot smrho
Ci                 smrout
Co     Stored:     ves aamom bxc cp ddpf dddpf ddpfr dlmwt dpf dpfr
Co                 gibbs gma gmar grrme mad mxy palp papg pf pfnc pfr
Co                 pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt rhos
Co                 rhrmx sop thetcl vdif vintr vrmax vshft smpot smrho
Co                 smrout
Co     Allocated:  *
Cio    Elts passed:smrho
Cio    Passed to:  bcast_strx iinit mpibc1
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Stored:     ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bcast_strx iinit mpibc1
Ci Inputs
Ci   mode  : ASA I/O
Ci         : 0 normal ASA I/O mode
Ci         :-1 (input only) read nbas,nat,nspec from file, return nsp0 in iors
Ci         : 1 Writing to ASA format from FP code
Ci         :-4 (input for ASA only) read the following from file:
Ci         :    s_lat%alat s_lat%plat
Ci         :    s_bz%ef s_bz%def s_bz%w s_bz%n
Ci         :    s_spec(is)%name s_spec(is)%a s_spec(is)%nr s_spec(is)%rmt s_spec(is)%z
Ci         :    s_spec(is)%lmxa s_spec(is)%lmxl s_spec(is)%kmxt
Ci         :    s_spec(is)%nxi s_spec(is)%exi s_spec(is)%chfa s_spec(is)%rsmfa
Ci         :    s_site(ib)%pos s_site(ib)%pos0 s_site(ib)%force s_site(ib)%vel
Ci         :    s_spec(is)%p s_spec(is)%pz
Ci         :    s_site(ib)%spec s_site(ib)%pnu s_site(ib)%pz
Ci         : 3 ASA-CPA mode: Reads/writes header and structural data only.
Ci         :   Data associated spheres (P,Q) are not touched by iorsa
Ci         : FP I/O, mesh density
Ci         : 1 normal FP I/O
Ci         :-1 (input only) read nbas,nat,nspec from file, return nsp0 in iors
Ci         :-2 (input only) FP species and site data; poke into strux:
Ci         :    lat->alat plat ngabc
Ci         :    bz
Ci         :    site->spec clabel pos pos0 force vel lmxa
Ci         :    spec->name a nr rmt z lmxa lmxl kmxt p pz qc ctail etail
Ci                    stc rhoc nxi exi chfa rsmfa lfoca rfoca
Ci         :-3 (input only) only site data returned:
Ci         :    site->spec clabel pos pos0 force vel lmxa
Ci   fid   :string containing identification
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lbin  :T file I/O in binary mode
Ci         :F file I/O in ascii mode
Cio Inputs/Outputs
Cio  qnu   :vector of moments (ASA only --- mode=0)
Cio  nit   :iteration number
Co  Outputs
Co   iorsa :-1 if file I/O failed
Co         :0  if file I/O successful
Cr Remarks
Cr   The density consists of a smooth part (smrho) plus
Cr   nbas atom-centered densities inside the MT spheres.
Cr   Their sum is the full charge density.
Cr   The local density is represented as the difference of the
Cr   two valence components in s_site%rho[12], the core in s_site%rhoc
Cr   Density in the MT spheres:
Cr      mesh parameters rmt,nr,a;
Cr      total density rho (times r**2) to lmxl;
Cr      a spherical potential v0 defining the wave functions within rmt
Cr      pnu and idmod to lmxa
Cr      NOTE: on input, arrays for rhoat and v0 are allocated here
Cr   Smooth density
Cr      real part of complex*16 array smrho contains the density
Cr      k1,k2,k3 are the physical dimensions of the array
Cr      n1,n2,n3 are the dimensions of the mesh.
Cr   Additional information stored:
Cr      fid: file identifier, a string of length 64 characters or less.
Cr      parameters relating to coordinates and molecular dynamics.
Cr   On input, iors tries to transform data format where needed:
Cr      lmxl incompatible: pad or truncate l-components
Cr      FT mesh changed: map onto new mesh
Cl Local variables
Cl   lrs switches:  0=>use rst file data; 1=>ignore rst file data
Cl   lrs(1) site positions
Cl   lrs(2) starting fermi level
Cl   lrs(3) starting pnu's
Cm MPI
Cm   Master process reads and broadcasts line by line. err= and end= are
Cm   troublesome since the slave processes may hang if the rst file is
Cm   corrupted or incompatible. For now if iors returns < 1 lmfp will
Cm   exit and hope the slave processes follow suit!
Cu Updates
Cu   17 Apr 18 reading mode -2 : better handling of missing floating sites
Cu   28 Oct 17 ASA can optionally read moments and potential
Cu   20 Jul 17 Enable translation of augmentation mesh to different mesh by interpolation
Cu   06 May 15 New mode 3
Cu   12 Nov 13 Preparation for noncollinear FP density
Cu   18 Jul 08 New mode -3
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   12 Jul 12 Bug fix, lshear=T, MPI case
Cu   25 Oct 11 Started migration to f90 structures
Cu   01 Jul 08 New mode -2
Cu   25 Aug 07 Bug fix, case floating orbitals not positioned at end
Cu   10 Jul 07 Cleaner error exit, MPI
Cu   20 Jun 06 Repackaged MPI
Cu   07 Jul 05 (version update 1.04)
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   04 Feb 05 Spin-split file non-spin-polarized density (FP)
Cu   26 Apr 03 Added MPI parallelization for ASA
Cu   11 Jan 03 Bug fix: corrected calc. vol (might have been<0)
Cu   10 Dec 02 File shears atom positions by shear transformation
Cu             (slat->plat) (file plat)^-1
Cu   19 Feb 02 (version update 1.03)
Cu             File now contains nspec
Cu             File contents backwardly compatible with prior versions.
Cu             New mode (-1)
Cu             Routine's argument list changed.
Cu   15 Feb 02 (ATP) Added MPI parallelization for fp
Cu   15 Jan 02 ascii version now labels site rho1,rho2,..
Cu   27 Aug 01 Extended to local orbitals.
Cu   17 May 01 Added ASA file I/O.  New argument list.
Cu   27 Apr 01 Added lbin switch
Cu   25 Jun 00 spin polarized
Cu   21 Apr 00 Adapted from nfp rw_rs
C  ----------------------------------------------------------------------
      use structures
      use mpi
      implicit none
C ... Passed parameters
      integer, parameter :: n0=10
      logical lbin
      integer mode,nbas,nat,nspec,nit,ifi
      double precision qnu(3,n0,*)
      double precision qnur(32*n0*n0,1)
      character fid*(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_bz)::    s_bz
C ... Dynamically allocated local arrays
      real(8), allocatable :: rwgt(:),oldmsh(:,:)
      complex(8), allocatable :: smrho(:),h(:)
      real(8), pointer :: v0(:,:),v1(:,:),rhoc(:,:),p_v0(:,:),p_v1(:,:)
      real(8), pointer :: rho1(:,:),rho2(:,:)
      real(8),allocatable :: rofi(:),rofi0(:)
C ... Local parameters
      real(8),parameter :: NULLR=-99999
      logical isanrg,lfail,ltmp1,ltmp2,lshear,mlog
      integer procid,master,mpipid,nproc,err
      integer i,i0,i1,i2,i3,i4,ib,ibaug,ic,intopt,iors,ipr,is,is0,isp,
     .  isw,jb,jfi,k,k1,k11,k2,k21,k3,k31,kmax,kmax0,kmxv,l,lfoc,lfoc0,
     .  lgunit,lmxa,lmxa0,lmxb,lmxb0,lmxl,lmxl0,lmxr,lmxv,lmxv0,lrel,
     .  lrel0,maxis,n1,n11,n2,n21,n3,n31,nat0,nbas0,ngabc(3),nlml,nlml0,
     .  nspc,nspc0,nr,nr0,nsp,nsp0,nspec0,nxi,stdo,ns4,iv(10)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer idmod(n0),idmoz(n0),lrs(10),nl,jx,npoly,ilm
      double precision a,a0,alat,alat0,cof,dr,drmx,dv,dvmx,eh,fac,qc,rfoc,rfoc0,rh,
     .  rmt,rmt0,rsma,rsma0,rsmfa,rsmr,rsmr0,rsmv,rsmv0,stc,sum,vfac,
     .  vol,vol0,vs,vs1,xx,z,z0,tolr
      double precision exi(n0),force(3),hfc(n0,2),plat(3,3),plat0(3,3),
     .  pnu(n0,2),pnus(n0,2),pnz(n0,2),pnzs(n0,2),pos(3),pos0(3),bxc(3),
     .  ql(n0,2*n0),qnul(3,n0,2),qlat(3,3),qlat0(3,3),vec0(3),vel(3),vrmax(2),wk(100)
      character spid*8,spid0*8,fid0*68,line*20,msg*23,use*80,ignore*80,
     .  msgw*17,datimp*32,usernm*32,hostnm*32,jobid*32,ffmt*32,ifmt*32,dc*1
      procedure(logical) aiomom,aiopot,cmdopt
      procedure(integer) nglob,iprint,wordsw,a2vec
      procedure(real(8)) dlength

C ... data statements
      data vec0 /0d0,0d0,0d0/

C ... MPI setup
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,ignore)
      intopt = 10*nglob('lrquad')
      nullify(rhoc)

      alat = s_lat%alat
      plat = s_lat%plat
      ngabc = s_lat%nabc
      if (mode >= 0) then
        call dinv33(plat,1,qlat,fac)
        vol = dabs(fac)*alat**3
      endif
      ipr    = iprint()
      stdo   = lgunit(1)
      if (mode == 3) then
        vs   =  1.06d0
      elseif (mode == 0 .or. mode == -4) then
        vs   = -1.05d0
      elseif (mode == -3 .or. mode == -2 .or. mode == -1 .or. mode == 1) then
C#ifdefC FP
C        vs   =  1.04d0
C#else
        vs   = -1.05d0
C#endif
      endif
      nsp    = nglob('nsp')
      nspc   = nglob('nspc')
      nspc = 1  ! for now
      ns4    = nsp*nspc
      lrel   = mod(nglob('lrel'),10)
      lrs(1) = isw(IAND(s_ctrl%lrs,32) /= 0)
      lrs(2) = isw(IAND(s_ctrl%lrs,64) /= 0)
      lrs(3) = isw(IAND(s_ctrl%lrs,128) /= 0)
      if (mode == -2 .or. mode == -4) then
        lrs(1) = 0
        lrs(3) = 0
        maxis = 0               ! highest species index read so far...to look for floating sites
      elseif (mode == -3) then
        lrs(1) = 0
        lrs(2) = 1
        lrs(3) = 1
      endif
      msg    = '         File mismatch:'
      msgw   = '         warning:'
      iors   = -1
      line = 'header'
C#ifdefC FP
C      if (mode == 1) call fftz30(n1,n2,n3,k1,k2,k3)
C#endif
      ffmt = '(5f15.10)'
      ifmt = '(20i5)'
C ... Hardwired for now
      nl = nglob('nl')

C ... Message
      if (ipr >= 10) then
        fid0 = 'read'
        i = 5
        if (ifi < 0) then
          fid0 = 'write'
          i = i+1
        endif
        fid0(i:) = ' restart file ('
        i = i+15
        fid0(i:) = 'ascii'
        if (lbin) fid0(i:) = 'binary'
        i = i+5
        if (lbin) i = i+1
        fid0(i:) = ', asa'
C#ifdefC FP
C        if (mode == -1) fid0(i:) = ' dimensioning data'
C        if (mode == 1) fid0(i:) = ', mesh density'
C#endif
        i = len_trim(fid0)+1
        fid0(i:i) = ')'
        write(stdo,700) fid0(1:i)
  700   format(/' iors  : ',a)
      endif

C --- Input ---
      if (ifi > 0) then
        jfi = ifi
        if (procid == master) then
          rewind jfi
        endif
        use    = '         use from  restart file:'
        ignore = '         ignore in restart file:'

C ...   Read header and general information
        line = 'header'
C       MPI check to see if at least 1st record can be read
C       Abort with error message if file is missing (lfail = .true.)
        lfail = .false.
        if (nproc > 1) then
          if (procid == master) then
            if (lbin) then
              lfail = .true.
              read (jfi,end=2,err=2) vs1
              lfail = .false.
              rewind jfi
            else
              lfail = .true.
              read (jfi,*,end=2,err=2)
              read (jfi,*,end=2,err=2) vs1
              lfail = .false.
              rewind jfi
            endif
          endif
    2     continue
          call mpibc1(lfail,1,1,mlog,'iors','read error')
          if (lfail) goto 998
        endif
        if (procid == master) then
        if (lbin) then
          read(jfi,end=998,err=998) vs1
          read(jfi) fid0
          read(jfi) datimp,usernm,hostnm
          if (abs(vs1) <= 1.021d0) then
            nspec0 = nspec
            read(jfi) nbas0,nsp0,nspc0,lrel0
            nat0 = nbas0
          else if (abs(vs1) <= 1.031d0) then
            read(jfi) nbas0,nsp0,nspc0,lrel0,nspec0
            nat0 = nbas0
          else
            read(jfi) nbas0,nat0,nsp0,nspc0,lrel0,nspec0
          endif
          read(jfi) nit
          read(jfi) alat0,vol0,plat0
        else
          read(jfi,*,end=998,err=998)
          read(jfi,*,end=998,err=998) vs1
          read(jfi,'(1x,a68)',end=998,err=998) fid0
          read(jfi,'(1x,a32)',end=998,err=998) datimp
          read(jfi,'(1x,a32)',end=998,err=998) usernm
          read(jfi,'(1x,a32)',end=998,err=998) hostnm
          if (abs(vs1) <= 1.021d0) then
            nspec0 = nspec
            read(jfi,*,end=998,err=998) nbas0,nsp0,nspc0,lrel0,nit
            nat0 = nbas0
          else if (abs(vs1) <= 1.031d0) then
            read(jfi,*,end=998,err=998) nbas0,nsp0,nspc0,lrel0,nit,nspec0
            nat0 = nbas0
          else
           read(jfi,*,end=998,err=998) nbas0,nat0,nsp0,nspc0,lrel0,nit,nspec0
          endif
          read(jfi,*) alat0,vol0
          read(jfi,*) plat0
        endif
        if (mode == -2 .or. mode == -4) then
          s_lat%alat = alat0; s_lat%plat = plat0
          alat = alat0; vol = vol0
          call dcopy(9,plat0,1,plat,1)
          call dinv33(plat,1,qlat,fac)
        endif
        fid = fid0
        call strip(fid,i,i1); call strip(datimp,i,i2);
        call strip(usernm,i,i3); call strip(hostnm,i,i4)
        if (ipr >= 40) write (stdo,12) fid(1:i1),usernm(1:i3),hostnm(1:i4),datimp(1:i2)

        if (mode == -1) then
          nbas = nbas0; nat  = nat0; nspec = nspec0; iors = nsp0
          goto 900
        endif
        endif  ! master node branch
        call mpibc1(vs1,1,4,mlog,'iors','vs1')
        call mpibc1(plat0,9,4,mlog,'iors','plat')
        call mpibc1(nbas0,1,2,mlog,'iors','nbas0')
        call mpibc1(nsp0,1,2,mlog,'iors','nsp0')
        call mpibc1(nat0,1,2,mlog,'iors','nat0')
        call mpibc1(nspc0,1,2,mlog,'iors','nspc0')
        call mpibc1(lrel0,1,2,mlog,'iors','lrel0')
        call mpibc1(plat,9,4,mlog,'iors','plat')
        call mpibc1(nit,1,2,mlog,'iors','nit')
        allocate(oldmsh(4,nspec)); call dpzero(oldmsh,size(oldmsh))

        call dinv33(plat0,1,qlat0,fac)
        lshear = dlength(9,plat-plat0,1) > 1d-6

C       if (isanrg(nspec0, nspec,nspec,msg,'nspec',.true.)) goto 999

C       Number of real atoms may not increase
        if (nat > nat0) then
          if (isanrg(nat0,  nat,nat,msg,'nat',.true.)) goto 999
        elseif (nat == nat0 .and. nbas /= nbas0) then
          call info2(10,0,0,
     .      '%9f(warning) mismatch in nbas ... skipping sites'//
     .      '%N%18f expected nbas=%i but rst file has nbas=%i',
     .      nbas,nbas0)
C       OK to reduce nat (e.g. E.S.); these sites will be skipped
        elseif (nat < nat0) then
          call info0(20,0,0,'%9f(warning) rst mismatch in nat ... skipping sites')
        endif
C       if (isanrg(nsp0,  nsp,nsp,  msg,'nsp', .true.)) goto 999
        if (nsp0 < nsp) call info0(20,0,0,
     .  '%9f(warning) rst file not spin pol .. splitting spins')
        if (nsp0 > nsp) call rx(
     .    'IORS : file is SPIN POL but calculation is not')
        if (isanrg(nspc0, nspc,nspc,msg,'nspc',.true.)) goto 999
        lfail = isanrg(lrel0, lrel,lrel,msgw,'lrel',.false.)
        call fsanrg(abs(vs1),1.00d0,abs(vs),0d0,' ','file''s version',.true.)

C   --- Read smooth charge density ---
        if (mode == 0 .or. mode==-4 .or. mode==3) then
C   ... FP mesh density
        elseif (mode==1 .or. mode==-2 .or. mode==-3) then
C#ifdefC FP
C        line = 'smoothed density'
C        if (procid == master) then
C        if (lbin) then
C          read(jfi,err=999,end=999) n11,n21,n31
C        else
C          read(jfi,*,err=999,end=999)
C          read(jfi,*,err=999,end=999) n11,n21,n31
C        endif
C        if (mode == -2 .or. mode == -3) then
C          n1 = n11
C          n2 = n21
C          n3 = n31
C          if (mode == -2) s_lat%nabc = ngabc
C          call fftz30(n1,n2,n3,k1,k2,k3)
C          i = k1*k2*k3
CC         Read into dummy; discard
C          allocate(smrho(i*nsp)); call dpzero(smrho,2*i*nsp)
C          call dpdftr(n1,n2,n3,k1,k2,k3,nsp0,smrho,lbin,jfi)
C          deallocate(smrho)
C        elseif (n11==n1 .and. n21==n2 .and. n31==n3) then
C          call dpdftr(n1,n2,n3,k1,k2,k3,nsp0*nspc0,s_pot%smrho,lbin,jfi)
C          if (nsp > nsp0) then
CC           call zprm3('smrho before split',0,s_pot%smrho,k1,k2,k3)
C            i = k1*k2*k3*2
C            call dscal(i,0.5d0,s_pot%smrho,1)
C            call dpscop(s_pot%smrho,s_pot%smrho,i,1,1+i,1d0)
CC           call zprm3('smrho after split',0,s_pot%smrho,k1,k2,k3*2)
C          endif
Cc   ... or read and remesh
C        else
C          if (ipr >= 10) write (stdo,3) n11,n21,n31,n1,n2,n3
C    3     format(9x,'remesh density from  ',i4,'  *',i4,'  *',i4,
C     .       '    to  ',i4,'  *',i4,'  *',i4)
C          call fftz30(n11,n21,n31,k11,k21,k31)
C          i = k11*k21*k31
C          allocate(h(i*nsp)); call dpzero(h,2*i*nsp)
C          call dpdftr(n11,n21,n31,k11,k21,k31,nsp0,h,lbin,jfi)
C          if (nsp > nsp0) then
CC           call zprm3('smrho before split',0,h,k1,k2,k3)
C            i = k11*k21*k31*2
C            call dscal(i,0.5d0,h,1)
C            call dpscop(h,h,i,1,1+i,1d0)
C          endif
C          call pshpr(1)
C          i = 0
C          if (n1 == 2*n11 .and. n2 == 2*n21 .and. n3 == 2*n31) i=3
C          call chgmsh(i,nsp,plat,n11,n21,n31,k11,k21,k31,h,
C     .                      plat,n1,n2,n3,k1,k2,k3,s_pot%smrho)
C          call poppr
C          deallocate(h)
C        endif
C        if (mode == -2 .or. mode == -3) goto 5
CC   ... If cell volume changed, scale smooth density to maintain charge
C        vfac = vol0/vol
C        if (dabs(vfac-1d0) > 1d-8) then
C          if (ipr >= 10) write (stdo,4) vol0,vfac
C    4     format(9x,'volume changed from',f8.2,
C     .      ' :  scale smooth density by',f8.4)
C          call dpcopy(s_pot%smrho,s_pot%smrho,1,2*k1*k2*k3*nsp,vfac)
C        endif
C        endif
C        call mpibc1(s_pot%smrho,2*k1*k2*k3*ns4,4,mlog,'iors','smrho')
C    5   continue
C#else
        call rx('smooth density not implemented')
C#endif
        endif

C   --- Read information related to dynamics ---
C       For compatibility with nfp, read record into wk
        if (procid == master) then
        if (lbin) then
          call dpdump(wk,100,jfi)
        else
          read(jfi,*,err=999,end=999)
          call dfdump(wk,100,jfi)
        endif
        endif
!         call mpibc1(wk,3,4,mlog,'iors','wk')
        call mpi_bcast(wk,4,mpi_real8,master,mpi_comm_world,err)

        if (lrs(2) /= 0) then
          call awrit0('%a ef window,',ignore,len(ignore),0)
        else
          call awrit0('%a ef window,',use,len(use),0)
          s_bz%ef = wk(1)
          s_bz%def = wk(2)
          s_bz%w = wk(3)
          s_bz%n = nint(wk(4))
        endif

C   --- Read atomic positions,forces,velocities ---
        line = 'site data'
        if (lrs(1) /= 0) then
          call awrit0('%a positions,',ignore,len(ignore),0)
        else
          call awrit0('%a positions,',use,len(use),0)
          if (lshear) call awrit0('%a%7b(sheared),',use,len(use),0)
C         Must be enough sites to read from rst file
          if (nbas0 < nbas) then
            call info0(2,0,0,'%9f oops ... cannot use '//
     .      'file site positions (site mismatch)')
            if (isanrg(nbas0,nbas,nbas,msg,'nbas',.true.)) goto 999
          endif
        endif

        do  ib = 1, nbas0
          if (procid == master) then
            if (lbin) then
              read(jfi,err=999,end=999) jb,pos,force,vel
            else
              read(jfi,'(6x,i10)',err=999,end=999) jb
              read(jfi,*,err=999,end=999) pos
              read(jfi,*,err=999,end=999) force
              read(jfi,*,err=999,end=999) vel
            endif
          endif
          call mpibc1(force,3,4,mlog,'iors','force')
          call mpibc1(pos,3,4,mlog,'iors','pos')
          call mpibc1(vel,3,4,mlog,'iors','vel')

          if (mode == -2 .or. mode == -3 .or. mode == -4) then
            s_site(ib)%pos = pos
            s_site(ib)%pos0 = pos
            s_site(ib)%force = force
            s_site(ib)%vel = vel
          endif

          if (ib > nbas) cycle

C         rst file positions in pos0
          if (lshear) then
            call dgemm('T','N',3,1,3,1d0,qlat0,3,pos,3,0d0,pos0,3)
            call dgemm('N','N',3,1,3,1d0,plat,3,pos0,3,0d0,pos,3)
          endif
          call dcopy(3,pos,1,pos0,1)

          if (lrs(1) /= 0) then
            pos = s_site(ib)%pos
            call dpzero(vel,3)
          endif
          s_site(ib)%pos = pos
          s_site(ib)%pos0 = pos0
          s_site(ib)%force = force
          s_site(ib)%vel = vel
        enddo

        if (mode == 3) iors = 0
        if (mode == -3 .or. mode == 3) goto 900

C   --- Read information for local densities ---
        if (lrs(3) /= 0) then
          call awrit0('%a pnu,',ignore,len(ignore),0)
        else
          call awrit0('%a pnu,',use,len(use),0)
        endif
        if (ipr >= 10) then
          call awrit1('%a%?#(p>32)#%b # *#',use,len(use),-stdo,0)
          call awrit1('%a%?#(p>32)#%b # *#',ignore,len(ignore),-stdo,0)
        endif
        if (procid == master) then
          if (.not. lbin) read(jfi,*,err=999,end=999)
        endif
        ibaug = 0
        do  ib = 1, nbas
          ic = s_site(ib)%class
          is = s_site(ib)%spec
C         mode = -2: get species is, lmxa, spid
          if (mode == -2 .and. lbin) then
C           Possibly a floating orbital.
C           If 1st entry > nspec => probably nr, otherwise assume floating
            read(jfi) is
            backspace jfi
            if (is > nspec) then
              call info2(10,0,0,
     .          '%9f(warning) species mismatch (%i) for ib=%i ... '//
     .          'ignore and assign to species 0',is,ib)
              s_site(ib)%spec = 0
              cycle
            endif
            read(jfi) is,spid,lmxa,lmxl,nr,rmt,a,z,qc
            backspace jfi
C           call spacks(1,'spec name',sspec,spid,is,is)
            s_spec(is)%name = spid
            s_site(ib)%spec = is
            maxis = max(maxis,is)
C           call spacks(1,'site clabel',ssite,spid,ib,ib)
            s_site(ib)%clabel = spid
            s_spec(is)%lmxa = lmxa
          elseif (mode == -2 .or. mode == -4) then
C           If ib > nat, assume floating orbital with spid = E
            if (ib > nat) then
              call info2(10,0,0,
     .          '%9f(warning) site %i exceeds nat=%i ... ignore and assign to species 0',ib,nat)
C             s_spec(is)%name = 'Null'
              s_site(ib)%spec = 0
C             s_spec(is)%lmxa = -1
              cycle
            endif
            read(jfi,'(19x,a8)') spid
            read(jfi,*)
            read(jfi,*) is,lmxa,lmxl,nr
C           call spacks(1,'spec name',sspec,spid,is,is)
            s_spec(is)%name = spid
            s_site(ib)%spec = is
            maxis = max(maxis,is)
            s_spec(is)%lmxa = lmxa
            backspace jfi
            backspace jfi
            backspace jfi
          endif
C         is = -1 -> spec struc does not have these parameters
C         lskip = .false.
          if (is == 0 .and. mode == 0) is = -1
          if (is /= -1) then
C           call spacks(0,'spec name',sspec,spid,is,is)
            spid = s_spec(is)%name
            a = s_spec(is)%a
            nr = s_spec(is)%nr
            rmt = s_spec(is)%rmt
            z = s_spec(is)%z
            lmxa = s_spec(is)%lmxa
            lmxl = s_spec(is)%lmxl
            kmax = s_spec(is)%kmxt
C           lmxa = -1 -> site has no augmentation
            if (lmxa == -1) then
C             if (ib > nat0) cycle
C             lskip = .true.
              cycle
            endif
          endif
          ibaug = ibaug+1
          if (procid == master) then
          if (lbin) then
            read(jfi) is0,spid0,lmxa0,lmxl0,nr0,rmt0,a0,z0,qc
            if (ib > nbas) cycle
          else
            read(jfi,'(19x,a8)') spid0
            read(jfi,*)
            read(jfi,*) is0,lmxa0,lmxl0,nr0
            read(jfi,*)
            read(jfi,*) rmt0,a0,z0,qc
            if (ib > nbas) cycle
          endif
          if (ipr >= 40) then
            if (ib <= nat) write (stdo,13) ib,is0,spid0
            if (ib > nat) write (stdo,13) ib,is0,spid0,' (skip)'
          endif
c     ... read(but don't use) extra info since record is present
          if (lbin) then
            read(jfi) rsma0,rsmr0,rsmv0,lmxv0,lmxr,lmxb0,kmax0
          else
            read(jfi,*)
            read(jfi,*) rsma0,rsmr0,rsmv0,lmxv0,lmxr,lmxb0,kmax0
          endif
          endif
          call mpibc1(lmxa0,1,2,mlog,'iors','lmxa0')
          call mpibc1(lmxl0,1,2,mlog,'iors','lmxl0')
          call mpibc1(nr0,1,2,mlog,'iors','nr0')
          call mpibc1(a0,1,4,mlog,'iors','a0')
          call mpibc1(qc,1,4,mlog,'iors','qc')
C         is = -1 -> pack spec struc with file parameters
          if (is == -1 .or. mode == -2 .or. mode == -4) then
            is = is0
C           call spacks(1,'spec name',sspec,spid0,is,is)
            s_spec(is)%name = spid0
            spid = spid0
            s_spec(is)%a = a0
            s_spec(is)%nr = nr0
            s_spec(is)%rmt = rmt0
            s_spec(is)%z = z0
            s_spec(is)%lmxa = lmxa0
            s_spec(is)%lmxl = lmxl0
            s_spec(is)%kmxt = kmax0
            a = s_spec(is)%a
            nr = s_spec(is)%nr
            rmt = s_spec(is)%rmt
            z = s_spec(is)%z
            lmxa = s_spec(is)%lmxa
            lmxl = s_spec(is)%lmxl
            kmax = s_spec(is)%kmxt
          endif

C     ... Read augmentation data
          call dpzero(pnu,2*n0)
          call dpzero(pnz,2*n0)
          call dpzero(ql,2*n0)
C         if (lmxa0 < lmxa) call defpq(0,z,lmxa,nsp,pnu,ql)
          pnu = s_spec(is)%p
          pnz = s_spec(is)%pz
          call iinit(idmod,n0)
          call iinit(idmoz,n0)
          if (procid == master) then
          do  isp = 1, nsp0
            if (lbin) then
              read(jfi) (pnu(l+1,isp), l=0,lmxa0)
              read(jfi) (pnz(l+1,isp), l=0,lmxa0)
            else
              if (isp == 1) read(jfi,*)
              read(jfi,*) (pnu(l+1,isp), l=0,lmxa0)
              read(jfi,*) (pnz(l+1,isp), l=0,lmxa0)
            endif
            if (nsp > nsp0) then
              do  l = 0, lmxa0
                pnu(l+1,2) = pnu(l+1,1)
                pnz(l+1,2) = pnz(l+1,1)
              enddo
            endif
          enddo
          endif
          do  isp = 1, nsp
            call mpibc1(pnu(1,isp),lmxa0+1,4,mlog,'iors','pnu')
            call mpibc1(pnz(1,isp),lmxa0+1,4,mlog,'iors','pnu')
C       ... For backwards compatibility: prior versions wrote pnu for pnz
            do  l = 0, lmxa0
              if (pnu(l+1,isp) == mod(pnz(l+1,isp),10d0)) pnz(l+1,isp) = 0
            enddo
          enddo
C         read(but don't use) site idmod for compatibility with nfp
          if (procid == master) then
          if (lbin) then
            read(jfi) (idmod(l+1), l=0,lmxa0)
            read(jfi) (idmoz(l+1), l=0,lmxa0)
          else
            read(jfi,*) (idmod(l+1), l=0,lmxa0)
            read(jfi,*) (idmoz(l+1), l=0,lmxa0)
          endif
          endif

C     ... Store pnu and pz in site and spec data
          pnus = s_spec(is)%p
          pnzs = s_spec(is)%pz
          if  (lrs(3) == 0 .or. mode == -2 .or. mode == -4) then
            s_site(ib)%pnu = pnu
            s_site(ib)%pz = pnz
            s_spec(is)%p = pnu
            s_spec(is)%pz = pnz
            if (mode == -2 .or. mode == -4) then
              pnus = s_spec(is)%p
              pnzs = s_spec(is)%pz
            endif

C     ... Verify lowest valence pnu compatible with file
          lfail = .false.
          ltmp1 = .false.
          ltmp2 = .false.
          do  i = 1, lmxa+1
            vec0(1) = mod(pnz(i,1),10d0)
            if (vec0(1) == 0) vec0(1) = pnu(i,1)
            vec0(2) = mod(pnzs(i,1),10d0)
            if (vec0(2) == 0) vec0(2) = pnus(i,1)
            ltmp1 = ltmp1 .or. int(pnu(i,1)) /= int(pnus(i,1))
            ltmp2 = ltmp2 .or. int(mod(pnz(i,1),10d0)) /= int(mod(pnzs(i,1),10d0))
            lfail = lfail .or. min(int(pnu(i,1)),int(vec0(1))) /=
     .                         min(int(pnus(i,1)),int(vec0(2)))
          enddo
          if (ltmp1 .and. ipr>=20) then
            write (stdo,14) ib,spid,'file pnu',(pnu(i,1),i = 1,lmxa+1)
            write (stdo,15) 'given pnu is',(pnus(i,1),i = 1,lmxa+1)
            call logwarn(2,'%26fwarning!  princ. q.n. mismatch')
          endif
          if (ltmp2 .and. ipr>=20) then
            write (stdo,14) ib,spid,'file pz ',(pnz(i,1),i = 1,lmxa+1)
            write (stdo,15) 'given pz  is',(pnzs(i,1),i = 1,lmxa+1)
            call logwarn(2,'%26fwarning!  local orbital mismatch')
          endif
          if (lfail .and. lrs(3) == 1) call rx('iors: file''s pnu is incompatible with input')

          nlml0 = (lmxl0+1)**2
          nlml = (lmxl+1)**2
          if (nr <= 0)   nr = nr0
          if (a <= 1d-6) a = a0

          if (procid == master) then
C     ... Sanity checks, or inform about changed parameters
C         if (is0 /= is) call xxerri('species pointer',is,is0)
          tolr = 0.03d0*rmt0
          call readrssw(fid0,i,[-1])
          if (i > 0) then
            dc = fid0(5:5)
            i = wordsw(fid0,dc,'dr','= ',k)
            if (i /= 0) then
              i = a2vec(fid0,len_trim(fid0),k,4,', '//dc,3,3,1,iv,tolr)
            endif
          endif

          call fsanrg(rmt0,rmt,rmt,1d-3,'file:         '//trim(spid)//':','rmt',.false.)
C#ifdefC FP
C          call fsanrg(rmt,rmt0-2*tolr,rmt0+2*tolr,1d-6,msg//'         '//trim(spid)//':','',.true.)
C#else
          call fsanrg(rmt,rmt0-2*tolr,rmt0+2*tolr,1d-6,msg//'         '//trim(spid)//':','',.false.)
C#endif
          call fsanrg(z0,z,z,1d-6,'file:         '//trim(spid)//':','z',.true.)
          call fsanrg(a0,a,a,0d-9,'file:         '//trim(spid)//':','a',.false.)
          lfail = isanrg(nr0,nr,nr,'file:         '//trim(spid)//':','nr',.false.)
          if (isanrg(lmxl,0,lmxa,'file:         '//trim(spid)//':','lmxl',.false.)) goto 999

C#ifdefC FP
C          if (kmax0 /= kmax .and. ipr>=10)
C     .      write(stdo,16) ib,spid,'kmax',kmax0,kmax
C#endif
          if (lmxa0 /= lmxa .and. ipr>=10)
     .      write(stdo,16) ib,spid,'lmax',lmxa0,lmxa

          endif

C         Case read but skip over this site data
          else
            lmxl = lmxl0
            nr = nr0
          endif

C     --- Allocate and read arrays for local density and potential ---
          nlml0 = (lmxl0+1)**2
          nlml = (lmxl+1)**2

          allocate(v0(nr,nsp),v1(nr,nsp))

C     ... FP local densities rho1,rho2,rhoc and potentials v0, v1
          if (mode == 1 .or. mode == -2) then
C#ifdefC FP
C            call dfratm(s_site,s_spec,1,ib,ib,xx)
C
CC       ... Read rho1,rho2,rhoc,v0,v1 from disk
C            if (procid == master) then ! Allocate temporary arrays if lattice mismatch
C              if (nr == nr0 .and. a == a0 .and. abs(rmt-rmt0) < 1d-8) then
C                rho1 => s_site(ib)%rho1
C                rho2 => s_site(ib)%rho2
C                rhoc => s_site(ib)%rhoc
C                p_v0 => v0
C                p_v1 => v1
C              else  ! Otherwise read into final destination
C                allocate(rho1(nr0,nlml*2*2),rho2(nr0,nlml*2*2))
C                allocate(rhoc(nr0,2),p_v0(nr0,2),p_v1(nr0,2))
C              endif
C
C              if (.not. lbin) read(jfi,*)
C              call dpdbyl(rho1,nr0,nlml0,nlml,nsp0*nspc0,ns4,lbin,jfi)
C              if (.not. lbin) read(jfi,*)
C              call dpdbyl(rho2,nr0,nlml0,nlml,nsp0*nspc0,ns4,lbin,jfi)
C              if (nlml0 > nlml.and.ipr>=10)
C     .          write(stdo,202) ib,spid,'truncate',nlml0,nlml
C              if (nlml0 < nlml.and.ipr>=10)
C     .          write(stdo,202) ib,spid,'inflate',nlml0,nlml
C  202         format(9x,'site',i4,', species ',a,
C     .               ': ',a,' local density from nlm=',i3,' to',i3)
C              if (.not. lbin) read(jfi,*)
C              call dpdbyl(rhoc,nr0,1,1,nsp0,ns4,lbin,jfi)
C              if (.not. lbin) read(jfi,*)
C              call dpdbyl(p_v0,nr0,1,1,nsp0,nsp,lbin,jfi)
C              if (.not. lbin) read(jfi,*)
C              call dpdbyl(p_v1,nr0,1,1,nsp0,nsp,lbin,jfi)
C              if (nsp0 < nsp) then
C                call dscal(nr0*2,2d0,v0,1)
C                call dscal(nr0*2,2d0,v1,1)
C              endif
CC             call prmx('iors: v0',v0,nr0,nr0,nsp)
C
CC             General polynomial interpolation.  Works, but slow.
C              if (nr /= nr0 .or. a /= a0 .or. abs(rmt-rmt0) > 1d-8) then
C                oldmsh(1,is) = 1; oldmsh(2,is) = nr0; oldmsh(3,is) = a0; oldmsh(4,is) = rmt0
C                allocate(rofi(nr),rofi0(nr0))
C                call radmsh(rmt,a,nr,rofi); call radmsh(rmt0,a0,nr0,rofi0)
C                dvmx = 0; drmx = 0
C                jx = 0
C                do  i = 1, nr
C                  npoly = 7 ; if (rofi(i) > rofi0(nr0)) npoly = 5
C                  do  isp = 1, nsp
C                    call polint(rofi0,p_v0(1,isp),nr0,npoly,rofi(i),0d0,0,jx,v0(i,isp),dv)
C                    dvmx = max(dvmx,abs(dv))
C!                    call polint(rofi0,p_v1(1,isp),nr0,npoly,rofi(i),0d0,0,jx,v1(i,isp),dv)
C                    dvmx = max(dvmx,abs(dv))
C                    call polint(rofi0,rhoc(1,isp),nr0,npoly,rofi(i),0d0,0,jx,s_site(ib)%rhoc(i,isp),dr)
C                    drmx = max(drmx,abs(dr))
C                    do  ilm = 1, nlml
C                      k = ilm+(isp-1)*nlml
C                      call polint(rofi0,rho1(1,k),nr0,npoly,rofi(i),0d0,0,jx,s_site(ib)%rho1(i,k),dr)
C                      drmx = max(drmx,abs(dr))
C                      call polint(rofi0,rho2(1,k),nr0,npoly,rofi(i),0d0,0,jx,s_site(ib)%rho2(i,k),dr)
C                      drmx = max(drmx,abs(dr))
C                    enddo
C                  enddo
C                enddo
C
C!          call snot
C!          print *, associated(v0),associated(v1)
C!          call ptr_site(s_site,4+1,'v1',ib,nr,nsp,v1)
C!          deallocate(v1)
CC                call prrmsh('v0(init) ',rofi0,p_v0,nr0,nr0,nsp)
CC                call prrmsh('v0(final)',rofi,v0,nr,nr,nsp)
CC                call prrmsh('v1(init) ',rofi0,p_v1,nr0,nr0,nsp)
CC                call prrmsh('v1(final)',rofi,v1,nr,nr,nsp)
CC                call prrmsh('rhoc(init) ',rofi0,rhoc,nr0,nr0,nsp)
CC                call prrmsh('rhoc(final)',rofi,s_site(ib)%rhoc,nr,nr,nsp)
CC                call prrmsh('rho1(init) ',rofi0,rho1,nr0,nr0,nlml*nsp)
CC                call prrmsh('rho1(final)',rofi,s_site(ib)%rho1,nr,nr,nlml*nsp)
CC                call prrmsh('rho2(init) ',rofi0,rho2,nr0,nr0,nlml*nsp)
CC                call prrmsh('rho2(final)',rofi,s_site(ib)%rho2,nr,nr,nlml*nsp)
C                 call info8(20,0,0,'%8f interpolate (a=%d,nr=%i,rmt=%,6;6d)->'//
C     .             '(a=%d,nr=%i,rmt=%,6;6d)  errr=%,3;3g  errV=%,3;3g',
C     .             a0,nr0,rmt0,a,nr,rmt,drmx,dvmx)
C                 deallocate(rofi,rofi0)
C                 if (nr == nr0 .and. a == a0) then
C                   nullify(rho1,rho2,rhoc,p_v0,p_v1)
C                 else
C                   deallocate(rho1,rho2,rhoc,p_v0,p_v1)
C                 endif
C              endif
C            endif
C            if (mode == -2) then
C              deallocate(s_site(ib)%rho1,s_site(ib)%rho2,s_site(ib)%rhoc)
C            else
C              i = nr*ns4
C              call mpibc1(s_site(ib)%rho1,i*nlml,4,mlog,'iors','rho1')
C              call mpibc1(s_site(ib)%rho2,i*nlml,4,mlog,'iors','rho2')
C              call mpibc1(s_site(ib)%rhoc,nr*nsp,4,mlog,'iors','rho3')
C              i = nr*nsp
C              call mpibc1(v0,i,4,mlog,'iors','v0')
C              call mpibc1(v1,i,4,mlog,'iors','v1')
C
C            endif
C#else
            call rxi('iors: not implemented, mode',mode)
C#endif
C     ... ASA sphere density
          else
            if (lbin) call rx('iors not ready for binary ASA input')
C           call dpdbyl(v0,nr0,1,1,nsp0,nsp,lbin,jfi)
            call dpzero(ql,3*n0*nsp)
            if (procid == master) then
              qnur(1,ib/ib) = NULLR
              read(jfi,'(a)') jobid
              backspace jfi
              if (jobid(1:7) == 'MOMNTS:') then
                ltmp1 = aiomom(spid0,pnu,qnu(1,1,1+(ib-1)*nsp),qnur(1,ib/ib),idmod,n0,lmxa,nsp,z,rh,vrmax,bxc,jfi)
                if (.not. ltmp1) call rxs('iors failed to read sphere moments for site ',spid0)
              else
                call info0(10,0,0,'%9f(warning) failed to read sphere moments for species '//spid0)
              endif
              if (vs1 < -1.041d0) then
                read(jfi,43) vrmax(1)
   43           format(12x,f11.7)
                if (ic == 0) ic = is
                s_pot%ves(ic) = vrmax(1)
              endif
              read(jfi,*,end=998,err=998) fid0; backspace jfi
              if (fid0(1:4) == 'POT:') then
                nr=0; isp=0; a=0; rmt=0
                ltmp1 = aiopot(nr,isp,a,rmt,[NULLR],[NULLR],jfi)
                if (.not. ltmp1) call rxs('iors failed to read potential for site ',spid0)
                deallocate(v0)
                allocate(v0(nr,isp)); v0(1,1) = NULLR
                if (.not. aiopot(nr,isp,a,rmt,[NULLR],v0,jfi))
     .          call rxs('iors failed to read potential for site ',spid0)
C               let a,nr be determined by file, require match to rmt, and spin for now
                call fsanrg(rmt,s_spec(is)%rmt,s_spec(is)%rmt,1d-6,'iors:','rmax',.false.)
                call sanrg(.true.,isp,nsp,nsp,'iors:','file nsp')
                if (associated(s_pot%v0)) then
                  iv(1:2) = shape(s_pot%v0)
                  if (iv(1) < nr*nsp) call resizev0(s_pot,nr*nsp)
                  ic = s_site(ib)%class
                  call dcopy(nr*nsp,v0,1,s_pot%v0(1,ic),1)
C                 For now, if nr mismatch, flag that v0 is missing
                  if (abs(rmt-s_spec(is)%rmt) < 1d-6) then
                    call info2(60,0,0,'%9fread v0(1:%i,1:nsp) into class %i',nr,ic)
                  else
                    s_pot%v0(1,ic) = NULLR
                    call info2(60,0,0,'%9fskip reading v0',nr,ic)
                  endif
                endif
!               s_pot%nrmxv = max(s_pot%nrmxv,nr)
              endif
              if (lrs(3) == 0) s_site(ib)%pnu = pnu
            endif  ! procid = master
          endif    ! asa sphere density parameters for this site

C     ... store data in strucs
          s_spec(is)%a = a
          s_spec(is)%nr = nr
          s_spec(is)%rmt = rmt
          s_spec(is)%z = z
          s_spec(is)%lmxa = lmxa
          s_spec(is)%lmxl = lmxl
          s_spec(is)%kmxt = kmax
          s_spec(is)%qc = qc

          call ptr_site(s_site,4+1,'v0',ib,nr,nsp,v0)
          call ptr_site(s_site,4+1,'v1',ib,nr,nsp,v1)
          deallocate(v0,v1)

          if (mode == 0) then
            call mpibc1(pnu,n0*nsp,4,mlog,'iors','pnu')
            call dcopy(n0*nsp,pnu,1,s_site(ib)%pnu,1)
          endif

        enddo  ! loop over sites
        if (isanrg(ibaug,nat,nat,msg,'nat',.false.)) goto 999
        if (mode == -2) nspec = maxis

C#ifdefC FP
CC       Rotate local density after reading
C        if (lshear .and. IAND(s_ctrl%lrs,256) /= 0) then
C          call dgemm('N','T',3,3,3,1d0,plat,3,qlat0,3,0d0,wk,3)
CC          print *, '!!'
CC          print "(9f12.6)", plat
CC          print "(9f12.6)", qlat0
CC          call dgemm('N','T',3,3,3,1d0,plat0,3,qlat0,3,0d0,wk,3)
C          call pvsms2(s_site,s_spec,wk,nbas,nsp)
CC         stop
C        endif
C#endif

C   --- Read data on free-atom core states and fit to fa density ---
        line = 'species data'
        if (procid == master) then
          if (.not. lbin) read(jfi,*,err=999,end=999)
        endif
        do  is = 1, nspec
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          lmxa = s_spec(is)%lmxa
          spid = s_spec(is)%name
          if (lmxa == -1) cycle
          if (procid == master) then
          if (lbin) then
            read(jfi,err=999,end=999) nr0,a0,qc,cof,eh,stc,lfoc0,rfoc0
          else
            read(jfi,*,err=999,end=999)
            read(jfi,*,err=999,end=999)
            read(jfi,*,err=999,end=999) nr0,lfoc0
            read(jfi,*,err=999,end=999)
            read(jfi,*,err=999,end=999) a0,qc,cof,eh
            read(jfi,*,err=999,end=999)
            read(jfi,*,err=999,end=999) stc,rfoc0
          endif
          lfail = isanrg(nr0,nr,nr,'file: species '//trim(spid)//':','nr',.false.)
          lfoc = s_spec(is)%lfoca
C#ifdefC FP
C          lfail = isanrg(lfoc0,lfoc,lfoc,'file: species '//trim(spid)//':','foca',.false.)
C#endif
          call fsanrg(a0,a,a,0d-9,msgw//' species '//trim(spid)//':','file mesh a',.false.)
          s_spec(is)%lfoca = lfoc0
          s_spec(is)%rfoca = rfoc0
          endif
          call mpibc1(s_spec(is)%lfoca,1,2,mlog,'iors','lfoca')
          call mpibc1(s_spec(is)%rfoca,1,4,mlog,'iors','rfoca')
          call mpibc1(qc,1,4,mlog,'iors','qc')
          call mpibc1(cof,1,4,mlog,'iors','cof')
          call mpibc1(eh,1,4,mlog,'iors','eh')
          call mpibc1(stc,1,4,mlog,'iors','stc')

C     ... FP core densities
          if (mode == 1 .or. mode == -2) then
          allocate(rhoc(nr,2))
          if (oldmsh(1,is) == 0) then
            rho1 => rhoc
          else                  ! Otherwise read into final destination
            allocate(rho1(nr0,nsp0))
          endif

          if (procid == master) then
            if (lbin) then
              call dpdump(rho1,nr0*nsp0,jfi)
            else
              call dfdump(rho1,nr0*nsp0,jfi)
            endif

C           General polynomial interpolation.  Works, but slow.
            if (oldmsh(1,is) /= 0) then
              if (nr0 /= oldmsh(2,is) .or. a0 /= oldmsh(3,is)) call rx('bug in iors')
              rmt0 = oldmsh(4,is);  rmt = s_spec(is)%rmt
              allocate(rofi(nr),rofi0(nr0))
              call radmsh(rmt,a,nr,rofi); call radmsh(rmt0,a0,nr0,rofi0)
              drmx = 0
              jx = 0
              do  i = 1, nr
                npoly = 7 ; if (rofi(i) > rofi0(nr0)) npoly = 5
                do  isp = 1, nsp
                  call polint(rofi0,rho1(1,isp),nr0,npoly,rofi(i),0d0,0,jx,rhoc(i,isp),dr)
                  drmx = max(drmx,abs(dr))
                enddo
              enddo
C              call prrmsh('rhoc(init) ',rofi0,rho1,nr0,nr0,nsp)
C              call prrmsh('rhoc(final)',rofi,rhoc,nr,nr,nsp)
              call info8(20,0,0,'%8f interpolate (a=%d,nr=%i,rmt=%,6;6d)->'//
     .          '(a=%d,nr=%i,rmt=%,6;6d)  errr=%,3;3g',
     .          a0,nr0,rmt0,a,nr,rmt,drmx,8)
              deallocate(rofi,rofi0,rho1)
            endif

C           Spin-split core density
            if (nsp > nsp0) then
              i = nr0
              call dscal(i,0.5d0,rhoc,1)
              call dpscop(rhoc,rho1,i,1,1+i,1d0)
            endif
          endif

          if (mode == 1) then
            call mpibc1(rhoc,nr*nsp,4,mlog,'iors','rhoca')
          else
            deallocate(rhoc)
            s_spec(is)%coreh = ' '
            s_spec(is)%pb1 = ' '
            s_spec(is)%pb2 = ' '
          endif
          call dpzero(exi,n0)
          call dpzero(hfc,n0*2)
          if (procid == master) then
          if (lbin) then
            read(jfi,err=999,end=999) rsmfa,nxi
            read(jfi,err=999,end=999)
     .        ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp0)
          else
            read(jfi,*,err=999,end=999)
            read(jfi,*,err=999,end=999) rsmfa,nxi
            read(jfi,*,err=999,end=999) (exi(i),i=1,nxi)
            call dfdump(hfc(1,1),nxi,jfi)
            if (nsp0 == 2) call dfdump(hfc(1,2),nxi,jfi)
          endif
          if (nsp > nsp0) then
            i = n0
            call dscal(i,0.5d0,hfc,1)
            call dpscop(hfc,hfc,i,1,1+i,1d0)
          endif
          endif
          call mpibc1(rsmfa,1,4,mlog,'iors','rsmfa')
          call mpibc1(nxi,1,2,mlog,'iors','nxi')
          call mpibc1(exi,nxi,4,mlog,'iors','exi')
          call mpibc1(hfc,nsp*nxi,4,mlog,'iors','hfc')

C     ... For now, mode 0 (ASA) stores no core data
          elseif (mode == 0 .or. mode == -4) then
            rsmfa = 0
            nxi = 0
            call dpzero(exi,n0)
            call dpzero(hfc,2*n0)
          endif
C         call p_spec_hfc(sspec,is,rsmfa,nxi,exi,hfc)
C         Pack species data
          s_spec(is)%ctail = cof
          s_spec(is)%etail = eh
          s_spec(is)%stc = stc
          if (associated(rhoc)) then
            call ptr_spec(s_spec,4+1,'rhoc',is,nr*nsp,0,rhoc)
            deallocate(rhoc)
          endif
          s_spec(is)%nxi = nxi
          s_spec(is)%exi = exi
          s_spec(is)%chfa = hfc
          s_spec(is)%rsmfa = rsmfa
        enddo

C   ... Copy or rescale cores, in case foca was switched on or off
        if (mode == 1) then
        do  ib = 1, nbas
          is = s_site(ib)%spec
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          lmxa = s_spec(is)%lmxa
          lfoc = s_spec(is)%lfoca
          qc = s_spec(is)%qc
          if (lmxa == -1) cycle
          if (lfoc > 0) then
            call dpcopy(s_spec(is)%rhoc,s_site(ib)%rhoc,1,nr*nsp,1d0)
          else
            allocate(rwgt(nr))
            call radwgt(intopt,rmt,a,nr,rwgt)
            call radsum(nr,nr,1,nsp,rwgt,s_site(ib)%rhoc,sum)
            fac = 1d0
            if (dabs(sum) > 1d-6) fac = qc/sum
            if (dabs(fac-1d0) > 1d-7 .and. ipr>=30)
     .         write(stdo,787) ib,qc,sum,fac
  787       format(' fix core chg: ib=',i4,'  qc,sum,fac=',3f12.6)
            call dpcopy(s_site(ib)%rhoc,s_site(ib)%rhoc,1,nr*nsp,fac)
            deallocate(rwgt)
          endif
        enddo
        endif

C   ... MPI broadcast all structures
        call bcast_strx(1+2+2**2+2**3+2**5+2**6,
     .    s_bz,s_ctrl,xx,s_pot,s_lat,xx,s_spec,s_site,xx,nspec,nbas)

        if (mode == 1) then
        elseif (mode == 0) then
          call mpibc1(qnu,3*n0*nsp*nbas,4,mlog,'iors','qnu')
        endif

C   ... debugging check
C        ib = 1; is = s_site(ib)%spec
C        a = s_spec(is)%a; nr = s_spec(is)%nr; rmt = s_spec(is)%rmt
C        allocate(rofi(nr))
C        lmxa = s_spec(is)%lmxa; lfoc = s_spec(is)%lfoca
C        print *, s_spec(is)%lfoca
C        call radmsh(rmt,a,nr,rofi)
C        call prrmsh('rhoc ',rofi,s_site(1)%rhoc,nr,nr,nsp)
C        call rx0('done')

C --- Output ---
      else
        if (procid /= master) then
          iors = 0
          goto 900
        endif
        jfi = -ifi
        rewind jfi
C ...   Write header and general information
        fid0 = fid
        call strip(fid0,i0,i1)
        call fextg(datimp)
        jobid = datimp(2:)
        call ftime(datimp)
        hostnm = ' '
        usernm = ' '
        call gtenv('HOST',hostnm)
        call gtenv('USER',usernm)
        call strip(datimp,i,i2)
        call strip(usernm,i,i3)
        call strip(hostnm,i,i4)
        if (ipr >= 40) write(stdo,12) fid(1:i1),
     .     usernm(1:i3),hostnm(1:i4),datimp(1:i2)

  721   format('----------------------- ',a,' -----------------------')

        if (lbin) then
          write(jfi) vs
          write(jfi) fid0
          write(jfi) datimp,usernm,hostnm,jobid
          if (abs(vs) <= 1.021d0) then
            write(jfi) nbas,nsp,nspc,lrel
          else if (abs(vs) <= 1.031d0) then
            write(jfi) nbas,nsp,nspc,lrel,nspec
          else
            write(jfi) nbas,nat,nsp,nspc,lrel,nspec
          endif
          write(jfi) nit
          write(jfi) alat,vol,plat
        else
          write (jfi,18) 'header'
          write(jfi,'(1f15.10,t50,a)') vs, '! version; id; date; name; host'
          write(jfi,'(1x,a68)') fid0
          write(jfi,'(1x,a32)') datimp
          write(jfi,'(1x,a32)') usernm
          write(jfi,'(1x,a32)') hostnm
C         write(jfi,'(''  nbas nsp  nspc lrel nit'')')
          if (abs(vs) <= 1.021d0) then
            write(jfi,'(5i5,t50,a)') nbas,nsp,nspc,lrel,nit,'! nbas nsp nspc lrel nit'
          else if (abs(vs) <= 1.031d0) then
            write(jfi,'(6i5,t50,a)') nbas,nsp,nspc,lrel,nit,nspec,'! nbas nsp nspc lrel nit nspec'
          else
            write(jfi,'(7i5,t50,a)') nbas,nat,nsp,nspc,lrel,nit,nspec,'! nbas nat nsp nspc lrel nit nspec'
          endif
          write(jfi,'(f15.10,1pe20.10,t50,a)') alat,vol,'! alat vol; plat'
          write(jfi,'(3f15.10)') plat
        endif

C   --- Write smooth charge density ---
C#ifdefC FP
C        if (lbin) then
C          write(jfi) n1,n2,n3
C        else
C          write (jfi,18) 'smooth density'
C          write(jfi,ifmt) n1,n2,n3
C        endif
C        call dpdftr(n1,n2,n3,k1,k2,k3,nsp*nspc,s_pot%smrho,lbin,-jfi)
C#endif

C   --- Write information related to dynamics ---
        call dpzero(wk,100)
        wk(1) = s_bz%ef
        wk(2) = s_bz%def
        wk(3) = s_bz%w
        wk(4) = s_bz%n
        if (lbin) then
          call dpdump(wk,100,-jfi)
        else
          write (jfi,18) 'geometric data: ef def w ... '
          call dfdump(wk,100,-jfi)
        endif
        do  ib = 1, nbas
          pos = s_site(ib)%pos
          force = s_site(ib)%force
          vel = s_site(ib)%vel
          if (lbin) then
            write(jfi) ib,pos,force,vel
          else
            write(jfi,'('' Site:'',i5)') ib
            write(jfi,345) pos,   'position'
            write(jfi,345) force, 'force'
            write(jfi,345) vel,   'velocity'
  345       format(3f15.10,t50,'! ',a)
          endif
        enddo

!!        print *, ifi,mode; stop '1345'

C   --- Write information for local densities ---
        if (mode == 3) then
          iors = 0
          goto 900
        endif
        if (ipr >= 50) write(stdo,8)
    8   format(/9x,'ib:spc la ll   rmt     nr   a     pnu')
        if (.not. lbin) write(jfi,18) 'site densities'
        do  ib = 1, nbas

          ic = s_site(ib)%class
          is = s_site(ib)%spec
          spid = s_spec(is)%name
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          z = s_spec(is)%z
          qc = s_spec(is)%qc
          idmod = s_spec(is)%idmod
          rsma = s_spec(is)%rsma
          lmxa = s_spec(is)%lmxa
          lmxl = s_spec(is)%lmxl
          lmxb = s_spec(is)%lmxb
          kmxv = s_spec(is)%kmxv
          rsmv = s_spec(is)%rsmv
          kmax = s_spec(is)%kmxt
          pnu = s_site(ib)%pnu
          pnz = s_site(ib)%pz
          v0 => s_site(ib)%v0
          v1 => s_site(ib)%v1
          if (lmxa == -1) cycle

          if (lbin) then
            write(jfi) is,spid,lmxa,lmxl,nr,rmt,a,z,qc
          else
            write(jfi,'('' site'',i5,''   spec: '',a8)') ib, spid
            write(jfi,'('' spec  lmxa lmxl  nr'')')
            write(jfi,ifmt) is,lmxa,lmxl,nr
            write(jfi,'(t7,'' rmt'',t21,''a'',t36,''z'',t51,''qc'')')
            write(jfi,ffmt) rmt,a,z,qc
          endif
C     ... Some extra info... lots of it useless or obsolete
          lmxr = 0
          lmxv = 0
          rsmr = 0
          if (lbin) then
            write(jfi) rsma,rsmr,rsmv,lmxv,lmxr,lmxb,kmax
          else
            write(jfi,'('' Next two lines for nfp compatibility'')')
            write(jfi,ffmt) rsma,rsmr,rsmv
            write(jfi,ifmt) lmxv,lmxr,lmxb,kmax
          endif

C     ... Write augmentation data
          do  isp = 1, nsp
            if (lbin) then
              write(jfi) (pnu(l+1,isp), l=0,lmxa)
              write(jfi) (pnz(l+1,isp), l=0,lmxa)
            else
              if (isp == 1) write(jfi,'('' pnu, pz, idmod, idmoz'')')
              write(jfi,ffmt) (pnu(l+1,isp), l=0,lmxa)
              write(jfi,ffmt) (pnz(l+1,isp), l=0,lmxa)
            endif
          enddo
C         Write for compatibility with nfp
          if (lbin) then
            write(jfi) (idmod(l+1), l=0,lmxa)
            write(jfi) (idmod(l+1), l=0,lmxa)
          else
            write(jfi,ifmt) (idmod(l+1), l=0,lmxa)
            write(jfi,ifmt) (idmod(l+1), l=0,lmxa)
          endif

C     ... Write arrays for local density and potential
          nlml = (lmxl+1)**2
C#ifdefC FP
C          if (.not. lbin) write(jfi,'('' rho1'')')
C          call dpdbyl(s_site(ib)%rho1,nr,nlml,nlml,ns4,ns4,lbin,-jfi)
C          if (.not. lbin) write(jfi,'('' rho2'')')
C          call dpdbyl(s_site(ib)%rho2,nr,nlml,nlml,ns4,ns4,lbin,-jfi)
C          if (.not. lbin) write(jfi,'('' rhoc'')')
C          call dpdbyl(s_site(ib)%rhoc,nr,1,1,nsp,nsp,lbin,-jfi)
C          if (.not. lbin) write(jfi,'('' v0'')')
C          call dpdbyl(v0,nr,1,1,nsp,nsp,lbin,-jfi)
C          if (.not. lbin) write(jfi,'('' v1'')')
C          call dpdbyl(v1,nr,1,1,nsp,nsp,lbin,-jfi)
C#else
          if (lbin) then
            call rx('iorsa not ready for binary mode')
C           call dpdump(pl,nl*nsp,-jfi)
C           call dpdump(s_site(ib)%rho1,nl*nsp*3,-jfi)
C           call dpdbyl(v0,nr,1,1,nsp,nsp,lbin,-jfi)
          else
            rh = 0
            vrmax(1) = 0
            vrmax(2) = 0

            if (mode == 0) then
              ltmp1 = aiomom(spid,pnu,qnu(1,1,1+(ib-1)*nsp),[NULLR],idmod,n0,lmxa,nsp,0d0,rh,vrmax,(/0d0,0d0,0d0/),-jfi)
              vrmax(1) = s_pot%ves(ic)
            elseif (mode == 1) then
              call dpzero(qnul,3*n0*nsp)
              forall (i=1:n0, is0=1:nsp) qnul(1,i,is0) = s_pot%qbyl(i+n0*(is0-1),ib)
              ltmp1 = aiomom(spid,pnu,qnul,[NULLR],idmod,n0,lmxa,nsp,0d0,rh,vrmax,(/0d0,0d0,0d0/),-jfi)
              vrmax(1) = s_pot%vesrmt(ib)
            endif
            write(jfi,143) vrmax(1)
  143       format('  Ves(rmax)=',f11.7)
            if (mode == 0 .and. associated(s_pot%v0)) then
              ic = s_site(ib)%class
              ltmp1 = aiopot(nr,nsp,a,rmt,[NULLR],s_pot%v0(1,ic),-jfi)
            elseif (mode == 1) then
              ltmp1 = aiopot(nr,nsp,a,rmt,[NULLR],v0,-jfi)
            endif
          endif
C#endif

          if (ipr >= 50) then
            write(stdo,349) ib,spid,lmxa,lmxl,rmt,nr,a,(pnu(l+1,1),l=0,lmxa)
            if (nsp == 2)  write(stdo,350) (pnu(l+1,2), l=0,lmxa)
C           if (nspc == 2) write(stdo,350) (pz(l+1,1), l=0,lmxa)
C            if (nspc == 2 .and. nsp == 2)
C     .        write(stdo,350) (pz(l+1,2), l=0,lmxa)
          endif
        enddo  ! Loop over ib


C   --- Write data on free-atom core states and fit to fa density ---
        if (.not. lbin) write(jfi,721) 'species data'

        do  is = 1, nspec
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          qc = s_spec(is)%qc
          lmxa = s_spec(is)%lmxa
          lfoc = s_spec(is)%lfoca
          rfoc = s_spec(is)%rfoca
          cof = s_spec(is)%ctail
          eh = s_spec(is)%etail
          stc = s_spec(is)%stc
          nxi = s_spec(is)%nxi
          exi = s_spec(is)%exi
          hfc = s_spec(is)%chfa
          rsmfa = s_spec(is)%rsmfa
          if (lmxa == -1) cycle

          if (lbin) then
            write(jfi) nr,a,qc,cof,eh,stc,lfoc,rfoc
          else
            spid0 = s_spec(is)%name
            write(jfi,'('' spec'',i5,2x,'' label '',a)') is, spid0
            write(jfi,'(''  nr   lfoc'')')
            write(jfi,ifmt) nr,lfoc
            write(jfi,'(t7,''  a '',t21,''qc'',t36,''cof'',t61,''eh'')')
            write(jfi,'(2f15.10,1p,2e25.15)') a,qc,cof,eh
            write(jfi,'(t17,'' stc'',t31,''rfoc'')')
            write(jfi,'(1p,2e25.15)') stc,rfoc
          endif
C     ... For now, FP only: ASA stores no core data
C#ifdefC FP
C          if (lbin) then
C            call dpdump(s_spec(is)%rhoc,nr*nsp,-jfi)
C            write(jfi) rsmfa,nxi
C            write(jfi) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp)
C          else
C            call dfdump(s_spec(is)%rhoc,nr*nsp,-jfi)
C            write(jfi,'('' core tail: rsmfa,n,e,cof=''/f15.10,i4)')
C     .        rsmfa,nxi
C            write(jfi,ffmt) (exi(i),i=1,nxi)
C            call dfdump(hfc(1,1),nxi,-jfi)
C            if (nsp == 2) call dfdump(hfc(1,2),nxi,-jfi)
C          endif
C#endif
        enddo

      endif

  349 format(i11,':',a4,2i2,f9.5,i5,f6.3,1x,8f6.3)
  350 format(41x,8f6.3)
  364 format(/9x,'ib:spc la ll   rmt     nr   a     pnu')

      iors = 0
      if (allocated(oldmsh)) deallocate(oldmsh)

C     Exit for both FP and ASA
  900 continue
C#ifdefC FP
C      iorsf = iors
C#else
      iorsa = iors
C#endif
      return

C ... Error handling
  998 continue
      if (allocated(oldmsh)) deallocate(oldmsh)
      if (ipr > 1) write(stdo,'('' iors  : empty file ... nothing read'')')
      return
  999 continue
      if (ipr > 1)
     .  write(stdo,'('' iors  : read failed in: '',a)') line

   12 format(9x,'id -  ',a/9x,'written by -  ',a,' on ',a,' at: ',a)
   13 format('   atom',i4,'    species',i4,':',a:a)
   14 format(9x,'site',i4,':',a,':',a,' is',8F6.2)
   15 format(26x,a,8F6.2)
   16 format(9x,'site',i4,', species ',a,': augmentation ',a,
     .        ' changed from',i2,' to',i2)
   17 format(9x,'site',i4,', species ',a,': ',a,
     .        ' local density from nlm=',i3,' to',i3)

   18 format('----------------------- ',a,' -----------------------')
   19 format(3F15.10,t50,'! ',a)

      end

C     Dummy entry to distinguish ASA and FP iors
C      integer function iorsfp(mode,s_ctrl,s_site,s_spec,s_lat,s_pot,
C     .  s_bz,fid,nbas,nat,nspec,qnu,nit,lbin,ifi)
C      implicit none
CC ... For structures
C      include 'structures.h'
C      type(str_ctrl)::  s_ctrl
C      type(str_site)::  s_site(*)
C      type(str_spec)::  s_spec(*)
C      type(str_lat)::   s_lat
C      type(str_pot)::   s_pot
C      type(str_bz)::    s_bz
CC ... Passed parameters
C      logical lbin
C      integer mode,nbas,nat,nspec,nit,ifi,n0
C      parameter (n0=10)
C      double precision qnu(3,n0,*)
C      character fid*(*)
C      integer iors
C
C      iorsfp = iors(mode,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
C     .  fid,nbas,nat,nspec,qnu,nit,lbin,ifi)
C      end
C      subroutine snot
C      end
