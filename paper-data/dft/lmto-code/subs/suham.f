      subroutine suham(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_str,s_strn)
C- Hamiltonian setup
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbasp nspec nl lham lcgf lgen3 lncol nspin lpgf ldlm
Ci                 nclasp ipc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lfp lgen3 ips ipc ncomp idcc
Cio    Passed to:  pp2alp
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  ng alat tolft plat qlat kv kv2 igv igv2 nabc gmax
Co     Stored:     igv igv2 ng gmax nabc
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:pos gv qlat igv igv2
Cio    Passed to:  sugvec sugvec0
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxl p pz idxdn rmt lmxb ncomp name ngcut orbp
Co     Stored:     idxdn ngcut orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  atfold makidx nscpa showbs sugcut uspecb
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec clabel
Co     Stored:     pnu pz norb
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  setnorb showbs pvioeu
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  pwmode pwemin pwemax npwpad lncol neula qss nbf
Co     Stored:     ndham ndofH ldham lmxax npwmin npwpad hord
Co     Allocated:  offH iprmb bdots
Cio    Elts passed:offH iprmb eula magf bdots
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz
Co     Stored:     nlma nlml
Co     Allocated:  pti
Cio    Elts passed:pp pti
Cio    Passed to:  *
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:kaps alph
Cio    Passed to:  rdstrx pp2alp
Cio Inputs
Ci   gfopts:string containing switches for GF program
Cl Local variables
Cl   ndim  :total number of lmto orbitals = nl**2 * nbas
Cl   npwmin:lower limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   npwmax:upper limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   nqdiv: loop over mesh of q-points to estimate npwmax
Cl        : nqdiv is fineness of q-mesh.
Cl  npwpad: a 'safety' padding to npwmax in case npwmax
Cl        : underestimates actual upper limit
Cr Remarks
Cr   This routine generates energy-independent hamiltonian setup.
Cr  *It generates and packs a table of hamiltonian offsets offH,
Cr   orbital permutation indices oindxo.
Cr
Cr  *For the ASA 2nd generation LMTO:
Cr   Extract order of potential function from gfopts
Cr   Transform pp's to alpha representation
Cu Updates
Cu   17 Jul 13 Print out bfield in FP case
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   18 Jan 12 (Belashchenko) Added call to setnorb
Cu             Modified call to pp2alp with new mode 10 for DLM
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Oct 10 Setup to keep 2nd gen ASA strux in memory (rdstrx)
Cu   07 Jul 08 Make sham->ndham = estimate for upper dimension of
Cu             hamiltonian, including possible PW part
Cu             Make sham->lmxax = largest lmax in basis
Cu   18 Apr 05 Force small parameter p -> 0 in 2-center turned on (ASA)
Cu   14 Feb 03 Makes and packs magnetic field B.sigma
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_strn) :: s_strn(*)
C ... Local parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer bitand,hord,i,i1,iprint,lcgf,ldham(16),ldim,lgen3,
     .  lham,lidim,lihdim,lncol,lpgf,nbasp,ndim,neul,nl,
     .  nsp,nspc,nspec,nspx,partok,lfp,nvi,nvl,ngabc(3),
     .  ib,is,lmxa,lmxl,stdo,stdl,nglob,nkaph,isw,nbf,lmxax,pwmode,
     .  j1,j2,j3,m,npw,npwmin,npwmax,ndham,ldlm,mode,str_pack
      double precision pwemin,pwemax,plat(3,3),qlat(3,3),q(3),Gmin,Gmax
      integer n0,nqdiv,npwpad
      parameter (n0=10,nqdiv=12)
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      integer ng
      logical ltmp,bittst,adnf
      character*80 outs,gfopts
      double precision qss(4),vmtz,pnu(n0,2),pz(n0,2),alat,tolgv
      integer, parameter :: lBextz=256
C ... Needed for iostr
      logical rdstrx
C     integer nttab
      double precision ckbas,cksumf,xx

C --- Setup ---
      nbasp = s_ctrl%nbasp
      nspec = s_ctrl%nspec
      nl = s_ctrl%nl
      lham = s_ctrl%lham
      lcgf = s_ctrl%lcgf
      lgen3 = s_ctrl%lgen3
      lfp = mod(s_ctrl%lfp,2)
      lncol = s_ctrl%lncol
      nsp = s_ctrl%nspin
      lpgf  = s_ctrl%lpgf(1)
      nkaph = nglob('nkaph')
      pwmode = s_ham%pwmode
      pwemin = s_ham%pwemin
      pwemax = s_ham%pwemax
      npwpad = s_ham%npwpad
      nspc  = 1
      nspx  = nsp
      if (bitand(lncol,1) /= 0) then
        nspc = 2
        nspx = 1
      endif
      ndim = nbasp * nl**2 * nkaph

C --- Load strux into memory ---
      if (lfp == 0 .and. lgen3 == 0) then
        ckbas = cksumf(s_lat%pos,3*nbasp)
        ltmp = rdstrx(00008,s_str,'STR', nl,nbasp,ckbas)
        ltmp = rdstrx(00008,s_str,'SDOT',nl,nbasp,ckbas)
      elseif (lgen3 /= 0) then
C       Strux not needed if calculated by Ewald
        if (IAND(s_ctrl%lgen3,2) == 0) then
        ckbas = cksumf(s_lat%pos,3*nbasp)
        ltmp = rdstrx(10008,s_str,'STR', nl,nbasp,ckbas)
        ltmp = rdstrx(10008,s_str,'SDOT',nl,nbasp,ckbas)
        endif
      endif

C     call defi(offH, -n0H*nkap0*(nbasp+1))
      call ptr_ham(s_ham,8+1,'offH',n0H*nkap0,(nbasp+1),xx)
      call ptr_ham(s_ham,8+1,'iprmb',ndim+3,0,xx)
C     call defi(oiprmb,ndim+3)
      stdo = nglob('stdo')
      stdl = nglob('stdl')

C --- Automatic downfolding (2nd generation only?) ---
      if (lgen3 == 0 .and. lfp == 0) then
        adnf = bittst(lham,4) .and. mod(s_ctrl%lrel,10) /= 2
C       Suppress if in gamma repsn (not needed for GF mode)
C       if (IAND(s_ctrl%lham,128) /= 0) adnf = 0
        vmtz = s_pot%vmtz
        call atfold(0,adnf,nl,nsp,vmtz,s_pot%pp,s_spec)
      endif

C --- Hamiltonian offsets, orbital permutation table ---
      call iinit(ldham,16)
      if (mod(pwmode,10) == 2) then
        call info0(20,0,0,' suham: LMTO basis will be excluded')
      else
        call makidx(nl,nkaph,1,nbasp,0,s_spec,s_ctrl%ips,s_ham%offH,
     .    s_ham%iprmb,ldham)
      endif
C     Dimension of LMTO portion of hamiltonian
      s_ham%nlmto = ldham(1)
      s_ham%ndham = ldham(1)

C ... Make rest of ldham
      ldham(4) = nspc
      ldham(5) = ldham(1) * nspc
      ldham(6) = ldham(2) * nspc
      ldham(7) = ldham(3) * nspc
      ldham(8) = nspx
      s_ham%ndofH = nbasp+1
      s_ham%ldham = ldham

C ... Assign norb for sites
      call setnorb(s_site,s_ham%offH,nbasp)

C ... Printout orbital positions in hamiltonian, resolved by l
      if (iprint() >= 50) then
        call showbs(s_site,s_spec,nkaph,0,s_ham%iprmb,ldham)
      endif

C ------- Potential- and implementation- specific setup -------

C --- FP setup ---
      if (lfp /= 0) then
        nvi = 0
        nvl = 0
        lmxax = -1
        do  ib = 1, nbasp
          is = s_site(ib)%spec
          lmxa = s_spec(is)%lmxa
          lmxl = s_spec(is)%lmxl
          pnu = s_spec(is)%p
          pz = s_spec(is)%pz
          call dcopy(n0,pz,1,pz(1,2),1)
C         Augmentation dimensioning parameters
          nvi = nvi + (lmxa+1)**2
          nvl = nvl + (lmxl+1)**2
C         Poke spec starting pnu to site
          s_site(ib)%pnu = pnu
          s_site(ib)%pz = pz
C         Find largest lmxa
          lmxax = max(lmxax,lmxa)
        enddo
        s_ham%lmxax = lmxax
        s_pot%nlma = nvi
        s_pot%nlml = nvl
        call info5(30,0,0,' suham :  %i augmentation'//
     .    ' channels, %i local potential channels  Maximum lmxa=%i',
     .    nvi,nvl,lmxax,0,0)

        ng = s_lat%ng
        alat = s_lat%alat
        tolgv = s_lat%tolft
        call sugcut(1,nspec,s_spec,alat,ng,s_lat%gv,tolgv)

C   ... PW setup : estimate upper bound to number of G vectors
C       to set up upper bound to hamiltonian dimension
        pwmode = s_ham%pwmode
        pwemax = s_ham%pwemax
        if (pwemax > 0 .and. mod(pwmode,10) > 0) then

          alat = s_lat%alat
          plat = s_lat%plat
          qlat = s_lat%qlat
          Gmin = dsqrt(pwemin)
          Gmax = dsqrt(pwemax)
          call iinit(ngabc,3)
          if (mod(pwmode/10,10) == 1) then
            call info0(70,1,0,' Estimate max size of PW basis from'//
     .        'combinations of recip. lattice vectors ...')
            npwmax = -1
            npwmin = 99999
            call pshpr(iprint()-40)
            do  j1 = 0, nqdiv
            do  j2 = 0, nqdiv
            do  j3 = 0, nqdiv
              do   m = 1, 3
                q(m) = (qlat(m,1)/nqdiv)*j1 +
     .                 (qlat(m,2)/nqdiv)*j2 +
     .                 (qlat(m,3)/nqdiv)*j3
              enddo
C             call shorbz(q,q,qlat,plat)
              call sugvec(s_lat,16,q,Gmin,Gmax,ngabc,0,npw)
              npwmin = min(npwmin,npw)
              npwmax = max(npwmax,npw)
            enddo
            enddo
            enddo
            call poppr
            if (npwpad < 0) then
              npwpad = max(nint((npwmax-npwmin)*0.2d0),3)
            endif
          else
            call dpzero(q,3)
            call pshpr(iprint()-40)
            call sugvec(s_lat,16,q,Gmin,Gmax,ngabc,0,npw)
C            call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
C     .        xx,xx,xx)
            call poppr
            npwmin = npw
            npwmax = npw
            npwpad = 0
          endif
          ndham = npwmax + npwpad
          if (mod(pwmode,10) /= 2) ndham = ldham(1) + npwmax + npwpad
          s_ham%npwmin = npwmin
          s_ham%npwpad = npwpad
          s_ham%ndham = ndham
          if (mod(pwmode/10,10) == 1) then
            call info2(20,1,0,' suham:  q-dependent PW basis with'//
     .        '  Emin = %d < E < %d.',pwemin,pwemax)
            call info5(30,0,0,'%9fEst. min,max PW dimension = %i,%i.'//
     .        '  Use npwpad = %i => ndham = %i',
     .        npwmin,npwmax,npwpad,ndham,0)
            if (stdl > 0)
     .        call awrit3('dm  nmto %i  npw(est) %i  tot %i',' ',80,
     .        stdl,ldham(1),npwmax,ndham)
          else
            call info5(20,0,0,' suham:  PW basis with  %d < E < '//
     .        '%d  =>  npw = %i,  ndham = %i',
     .        pwemin,pwemax,npw,ndham,0)
            if (stdl > 0)
     .        call awrit3('dm  nmto %i  npw %i  tot %i',' ',80,
     .        stdl,ldham(1),npwmax,ndham)
          endif

C    ...  Printout APW basis
          if (iprint() >= 40) then
            call dpzero(q,3)
            call info0(40,1,-1,' G vectors at the Gamma point:')
            call pshpr(iprint())
            if (iprint() >= 50) call setpr(100)
            call sugvec(s_lat,16+9,q,Gmin,Gmax,ngabc,0,npw)
            call setpr(10)
            call sugvec0(s_lat)
            call poppr
          endif
        else
          if (stdl > 0)
     .      call awrit3('dm  nmto %i  npw %i  tot %i',' ',80,
     .      stdl,ldham(1),0,ldham)
        endif
      endif

C --- Third-generation NMTO ---
      if (lgen3 /= 0) then

      endif

C --- Second-generation LMTO, ASA ---
C ... Green's function-specific initialization
      if (lcgf /= 0 .or. lpgf /= 0) then
C       gfopts = ' '
        i1 = str_pack('gfopt',-2,s_strn,gfopts)
        if (i1 > 0) then
C         gfopts = sstrn(i1:i2)
          if (gfopts /= ' ') then
          call partk0(0,len(gfopts),1,-1,0,len(gfopts),-1,31,.false.)
          hord = 2
          ltmp = .false.
          i = partok(gfopts,'pz',  ' ;',ltmp,' ',0,0,0,0)
          if (ltmp) hord = 4
          i = partok(gfopts,'p3',  ' ;',ltmp,' ',0,0,0,0)
          if (ltmp) hord = 3
          i = partok(gfopts,'p1',  ' ;',ltmp,' ',0,0,0,0)
          if (ltmp) hord = 1
          s_ham%hord = hord
          if (iprint() >= 30) call awrit1(' SUHAM : 2nd generation'//
     .      ' ASA; potential functions P%?#n<4#%-1j%i#(z)#',' ',
     .      80,stdo,hord)
        endif
        endif
      endif

C ... 2nd gen LMTO-specific initialization
C     Make pti here since pp in gamma rep
      if (lcgf == 0 .and. lfp == 0 .and. lgen3 == 0) then
        call ptr_pot(s_pot,8+1,'pti',ndim*nsp,0,xx)
        call makipt(nl,nbasp,nsp,s_ctrl%ipc,lihdim,s_ham%iprmb,
     .    s_pot%pp,ndim,s_pot%pti)
      endif

C ... Transform pp's to alpha representation
      if (lfp == 0 .and. lgen3 == 0) then
        mode = 0
        ldlm = s_ctrl%ldlm
        if (ldlm /= 0) mode = 10
        call pp2alp(mode,s_ctrl,s_str,lham,s_pot%pp)
      endif

C ... Printout Euler angles
      lncol = s_ham%lncol
      neul = s_ham%neula
      qss = s_ham%qss
      if (iprint() >= 30 .and. bitand(lncol,1) /= 0) then
        call awrit4('%x          %?#n#Qss%3:1;6d angle %;6d.  #%2j#'//
     .    '%?#n#Euler angles:#noncollinear hamiltonian',outs,80,
     .    -stdo,bittst(lncol,2),qss,qss(4),isw(bittst(lncol,1).and.
     .    iprint() > 30))
        if (iprint() > 30 .and. bittst(lncol,1))
     .    call pvioeu(101,s_site,s_ham%eula,nbasp,neul)
      endif

C ... Make B.sigma for magnetic field
      if (lfp == 0 .and. lgen3 == 0 .and. bittst(lncol,8)) then
        nbf = s_ham%nbf
        call ptr_ham(s_ham,8+1,'bdots',4*lihdim,0,xx)
        call mkbfld(nl,nbasp,lihdim,s_ham%iprmb,s_ham%magf,nbf,
     .    s_ham%eula,neul,s_ham%bdots)
      endif


C ... Make Bz.sigmaz for magnetic field
      if (lfp /= 0 .and. iand(lncol,lBextz) /= 0) then
        nbf = s_ham%nbf
        call prbfield(nl,nbasp,s_ham%iprmb,s_ham%magf,nbf)
      endif

C     s_ham%obdots = obdots

      end
      subroutine showbs(s_site,s_spec,nkaph,iprma,iprmb,ldham)
C- Makes hamiltonian offsets and permutation indices
C ----------------------------------------------------------------------
Ci Inputs
Ci   nkaph :number of types of one l-quantum number in the basis
Ci   iprma :if supplied, iprma holds a permutation table of site
Ci         :indices.  Orbital of the hamiltonian are ordered by
Ci         :sites iprma(1),iprma(2),iprma(3),...
Ci         :iprma(1)=0=> iprma is not supplied.  Then the orbitals
Ci         :are ordered by sites 1,2,3,...
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:   spid
Ci     Stored:
Co Outputs
Co   Routine only prints out information
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nkaph,iprma(*)
      integer iprmb(*),ldham(4)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer nkap0,n0
      parameter (nkap0=4,n0=10)
      integer l,ndim,ipr,nglob,nbas,ldim,off,offs,
     .  specw,fieldw,iorb,offsi,ib,is,norb,
     .  ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      logical lprmib
      character sout*140,spid*8,fmt*20,spdlab*7,lsym*1

C --- Setup ---
      nbas = nglob('nbas')
      ldim = ldham(1)
      if (nkaph > nkap0) call rx('showbs: increase nkap0')
C     stdo = lgunit(1)
C     stdl = lgunit(2)
      lprmib = iprma(1) > 0
      if (lprmib) call rx('showbs not ready for iprma')
      call getpr(ipr)
      if (ipr < 10) return
      spdlab = 'spdfghi'

C     Field width = number of digits each offset uses * 2 + padding
      fieldw = int(dlog(dble(ldim))/dlog(10d0)+1d-14)+1

C     Get longest spec name + site
      specw = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        spid = s_spec(is)%name
        specw = max(specw,len(trim(spid)))
      enddo
      call awrit2('(i4,3x,a%i,i%i)',fmt,len(fmt),0,specw+2,fieldw+1)
      offs = 4+2 + specw+2 + 2*fieldw + 6

      call info0(0,1,0,
     .  ' Orbital positions in hamiltonian, resolved by l:%N'//
     .  ' Site  Spec  Total    By l ...')

C     Loop over sites
      do  ib = 1, nbas
        is = s_site(ib)%spec
        spid = s_spec(is)%name
        call orbl(ib,0,ldim,iprmb,norb,ltab,ktab,off,offl,ndim)

        write(sout,fmt) ib, spid, offl(1)+1
        call awrit1('%a:%i',sout,len(sout),0,offl(1)+ndim)

C       Loop over orbitals
        do  iorb = 1, norb
          l   = ltab(iorb)
C           = ktab(iorb)  <- not needed for single-kappa hamiltonians
          off = offl(iorb)

          offsi = (iorb-1)*(2*fieldw + 5)
          lsym = spdlab(l+1:l+1)
          call awrit3('%np%i:%i('//lsym//')',sout,len(sout),0,
     .      offs+offsi,offl(iorb)+1,offl(iorb)+2*l+1)
        enddo
        call info0(0,0,0,trim(sout))

      enddo

      end
