C#define PRTNOCR
      subroutine rsedit(sopts,mode,s_ctrl,s_site,s_spec,s_lat,s_pot,
     .  s_ham,s_bz,nbas,nat,nspec)
C- Restart file editor
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  zbak lrs
Co     Stored:     lrs
Co     Allocated:  *
Cio    Elts passed:lrs
Cio    Passed to:  iorsf bcast_strx iinit mpibc1
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec rho1 rho2 rhoc rho1x rho2x rhocx
Co     Stored:     vel pnu pz
Co     Allocated:  rho1 rho2 rhoc rho1x rho2x rhocx
Cio    Elts passed:rho1 rho2 rhoc rho1x rho2x rhocx
Cio    Passed to:  dfratm
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name nr lmxl
Co     Stored:     lmxa
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  dfratm
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc vol ng qlat awald tol nkd nkq dist ag
Ci                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg kv igv
Ci                 igv2 kv2 pos qlv symgr
Co     Stored:     alat plat nabc ag bgv cg cy dlv gv gvq indxcg ips0
Co                 istab jcg kv igv igv2 kv2 pos qlv symgr
Co     Allocated:  *
Cio    Elts passed:gv kv cg indxcg jcg cy qlv dlv
Cio    Passed to:  rdovfa ovlpfa ovlocr hxpbl ghibl hklbl gklbl iorsf
Cio                bcast_strx iinit mpibc1 prsed2 prsed4
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  nlml nlma ves aamom bxc cp ddpf dddpf ddpfr dlmwt
Ci                 dpf dpfr gibbs gma gmar grrme mad mxy palp papg pf
Ci                 pfnc pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp
Ci                 qt rhos rhrmx sop thetcl vdif vintr vrmax vshft
Ci                 smpot smrho smrout
Co     Stored:     ves aamom bxc cp ddpf dddpf ddpfr dlmwt dpf dpfr
Co                 gibbs gma gmar grrme mad mxy palp papg pf pfnc pfr
Co                 pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt rhos
Co                 rhrmx sop thetcl vdif vintr vrmax vshft smpot smrho
Co                 smrout
Co     Allocated:  *
Cio    Elts passed:smrho rhat
Cio    Passed to:  rdovfa iorsf bcast_strx iinit mpibc1
C
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Stored:     ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  iorsf bcast_strx iinit mpibc1
Ci Inputs
Ci   sopts :command options performed automatically, before reading
Ci         :from standard input
Ci   mode  :0 ASA I/O
Ci         :1 FP  I/O, mesh density
Ci         :-1 (input only) read nbas,nat,nspec from file
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Co Outputs
Co   rsedit never returns.
Co   rst file can be written.
Cr Remarks
Cr   The density consists of a smooth part (smrho) plus
Cr   nbas atom-centered densities inside the MT spheres.
Cr   Their sum is the full charge density.
Cr   The local density is represented as the difference of the
Cr   two valence components in s_site%rho[12], plus the core density.
Cr   Density in the MT spheres:
Cr      mesh parameters rmt,nr,a;
Cr      total density rho (times r**2) to lmxl;
Cr      a spherical potential v0 defining the wave functions within rmt
Cr      pnu and idmod to lmxa
Cr   Smooth density
Cr      real part of complex*16 array smrho contains the density
Cr      k1,k2,k3 are the physical dimensions of the array
Cr      n1,n2,n3 are the dimensions of the mesh.
Cr   For mapping a file to cell with a different mesh, do e.g.
Cr     lmscell ctrl.fe --wsite~map~fn=site2
Cr     lmf '--rsedit~scell#mshdat#m123=16,16,16#n123=32,32,32#fn=in.fe site2'
Cu Updates
Cu   25 Aug 17 New scell#mshdat
Cu   07 Aug 17 New symrho and verb
Cu   04 Aug 15 Some redesign.  New 'set asite' and 'add asite'
Cu   22 Jul 13 new flsite
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu             New scell command
Cu   25 Oct 11 Started migration to f90 structures
Cu   29 Jun 11 New allnc option
Cu   10 Mar 11 new rmsite and ssave options
Cu   19 Dec 08 New exch option
Cu   01 Jul 08 First created
C  ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
      integer mode,nbas,nat,nspec,ifi,n0
      parameter (n0=10)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_pot),target::   s_pot
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
C ... Local Dynamical arrays
      type(str_pot), pointer::   s_pot1
      type(str_site), allocatable::  s_site1(:),s_site2(:)
      type(str_spec), allocatable::  s_spec1(:),s_spec2(:)
      integer, allocatable:: ilst(:,:)
      integer, allocatable :: imap(:),ipsx(:)
C     real(8), pointer  :: p_v0(:),p_v1(:),p_v01(:),p_v11(:)
C     real(8), pointer  :: p_v0(:),p_v1(:)
      character, allocatable:: slabl(:)*8
      real(8),pointer :: posx(:,:)
      complex(8),allocatable:: meshdat(:,:,:),meshnew(:,:,:)
      complex(8),allocatable:: smrho(:,:,:,:),smrho2(:,:,:,:),cg1(:)
C ... Local parameters
      type(str_lat)::   s_lat2
      type(str_pot)::   s_pot2
      logical lnsave,lbin,lsopts,have2s,havers,lfail,lrmsite,lchsite
      integer nglob,fopna,fopng,iorsf,stdo,iosite,ioibmap
      integer i,j,k,l,im,m,n,j1,j2,js1,js2,nit(2),nspeck(2),isync,nw,
     .  lexch,nspec10,nbas10,lmxa,specsize,ssizei,ng,isite,n123(3),m123(3)
      integer nsp,ngabc(3),n1,n2,n3,k1,k2,k3,lrssav,nvl,nchan,iosits
      integer ngabc2(3)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C     integer oorh1,oorh2
C     integer osmrhf,oqmom
      integer nlst(2),ix(10),nbask(2),natk(2),isw
C     integer idup(3)
      double precision qbg,fac(2),alat,plat(3,3),vol,xx,xxv(20)
      double precision Qi_c0
C     double precision plat2(3,3)
      character dc*1, dc2*1, fn*120, fileid*68, outs*150, spid*8
C     for supercell
      integer nbx,ngmx,k1x,k2x,k3x,ibx,ib,is
      double precision alatx,volx,plax(3,3),qlat(3,3),qlx(3,3),dlength
      double precision QL_c0,QLs,QLc,Qi,QL
      procedure(logical) :: cmdopt
      procedure(integer) :: wordsw,a2vec,rdm

C ... data statements
C     data vec0 /0d0,0d0,0d0/

C     Defaults
      stdo = nglob('stdo')
      fn = 'rst1'
      havers = .false.
      have2s = .false.
C     # spec and # sites in spec 2nd atom density
      nbask(1) = nbas
      natk(1) = nat
      nspeck(1) = nspec
      nbask(2) = nbas
      nspeck(2) = nspec
      natk(2) = nat
      nit(1) = 1
C     true if rs data available, but not saved
      lnsave = .false.
      nsp  = nglob('nsp')
      qbg  = s_ctrl%zbak(1)
      nvl = s_pot%nlml
      nchan = s_pot%nlma
      alat = s_lat%alat
      plat = s_lat%plat
      s_pot1 => s_pot
      lrssav = s_ctrl%lrs

C     Allocate local species structures; copy original to them.
C     Make space for extra species
      nspec10 = nspec+10
      allocate(s_spec1(nspec10))
      allocate(s_spec2(nspec10))
      specsize = ssizei(s_spec2(1),s_spec2(2))
      s_spec1(1:nspec) = s_spec(1:nspec)
      s_spec2(1:nspec) = s_spec(1:nspec)
C     Allocate local site structures; copy original to them
C     Make space for extra sites
      nbas10 = nbas+10
      allocate(s_site1(nbas10))
      allocate(s_site2(nbas10))
      s_site1(1:nbas) = s_site(1:nbas)
      s_site2(1:nbas) = s_site(1:nbas)
C     Allocate second lat structure; copy original to it
      s_lat2 = s_lat
C     Allocate second pot structure; copy original to it
      s_pot2 = s_pot

      dc = sopts(1:1)
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the restart file editor. ',
     .    'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the restart file editor.  ',
     .    'Enter ''?'' to see options.')
        lsopts = .false.
      endif
      if (mode /= 1) call rx('rsedit not ready for mode ne 1')

C ... Return here to resume parsing for arguments
   10 continue
      if (lsopts) then
        js2 = js2+1
        if (js2 > len(sopts)) then
          lsopts = .false.
          goto 10
        endif
        if (sopts(js2:js2) == dc) goto 10
        js1 = min(len(sopts),js2)
        call nwordg(sopts,0,dc,1,js1,js2)
        if (js2 < js1) lsopts = .false.
      endif

C 306 format(' Failed to parse string ',a,' ... try again.')
  100 continue
C#ifdef PRTNOCR
      print '(/'' Option : '',$)'
C#elseC
C      print '(/'' Option : '')'
C#endif
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif
      call locase(outs)

C ... Parse and execute the next command
      if (.false.) then

      elseif (outs == ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',
     .    ' ''?'' to see menu')
        goto 10

C ... Read 2nd density from restart file
      elseif (outs(1:6) == 'readb ' .or. outs(1:6) == 'reada ' .or.
     .        outs(1:5) == 'read ' .or. outs(1:7) == 'readfa ') then
        if (outs(1:7) == 'readfa ') then
C         Make space for smooth density
          ngabc2 = s_lat%nabc
          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
C         call defcc(osmrh2, -k1*k2*k3*nsp)
          allocate(s_pot2%smrho(k1*k2*k3,nsp))
C         s_pot%osmrho = osmrh2
          call dfratm(s_site2,s_spec2,1,1,nbask(2),s_pot2%rhat) ! allocate rho2 density
C         call defi(oorh2,3*nbask(2))    ! rho2 density in rho[12c]x
          call rdovfa(nbask(2),nspeck(2),s_site2,s_spec2,s_lat2,s_pot2,s_ham,qbg)
        else
          call word(outs,2,j1,j2)
          if (outs(1:6) == 'reada ') then
            if (j2 < j1) then
              fn = 'rsta'
            else
              fn = outs(j1:j2)
            endif
            ifi = fopna(fn,-1,1)
            lbin = .false.
          else
            if (j2 < j1) then
              fn = 'rst'
            else
              fn = outs(j1:j2)
            endif
            ifi = fopna(fn,-1,4+1)
            lbin = .true.
          endif
          s_ctrl%lrs = 32*0+64+128*0
C         Get the number of atoms and species
          k = iorsf(-1,s_ctrl,xx,xx,s_lat2,s_pot2,s_bz,
     .      fileid,nbask(2),natk(2),nspeck(2),xx,nit(2),lbin,ifi)
          call info5(10,0,0,' rsedit: '//
     .      'read 2nd density from file "'//trim(fn)//'" ...  '//
     .      'nbas = %i  nat = %i  nspec = %i',nbask(2),natk(1),nspeck(2),
     .      0,0)
C         Allocate species and site structures, fill with available info
          deallocate(s_spec2,s_site2)
          allocate(s_spec2(nspeck(2)))
          allocate(s_site2(nbask(2)))
          call bcast_strx(-1-2,xx,xx,xx,xx,xx,xx,s_spec2,s_site2,xx,
     .      nspeck(2),nbask(2))
          do  i = 1, nspeck(2)
            s_spec2(i)%p = 0
            s_spec2(i)%pz = 0
          enddo
          call pshpr(1)
          k = iorsf(-2,s_ctrl,s_site2,s_spec2,s_lat2,s_pot2,s_bz,
     .      fileid,nbask(2),natk(2),nspeck(2),xx,nit(2),lbin,ifi)
          call poppr
C         Check for unassigned species => reduce basis
          k = 0
          do  i = 1, nbask(2)
            j = s_site2(i)%spec
            if (j == 0) exit
            k = i
          enddo
          if (k /= nbask(2)) then
            call info2(10,0,0,'%9f(warning) only %i sites had '//
     .        'associated species ... reducing nbas(2)',k,0)
          endif
          nbask(2) = k
C         Make space for smooth density
          ngabc2 = s_lat2%nabc
          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
          allocate(s_pot2%smrho(k1*k2*k3,nsp))
          allocate(s_pot2%smpot(k1*k2*k3,nsp))
C         Pack volume
          call dfratm(s_site2,s_spec2,8+1,1,nbask(2),s_pot2%rhat) ! link rho2 to s_pot2%rhat
          k = iorsf(1,s_ctrl,s_site2,s_spec2,s_lat2,s_pot2,s_bz,
     .      fileid,nbask(2),natk(2),nspeck(2),xx,nit(2),lbin,ifi)
          call fclose(ifi)
        endif

        have2s = .true.
        s_ctrl%lrs = lrssav

        call info0(2,1,0,' Second Density:')
        call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
        call prsed2(nbask(2),nsp,s_site2,s_spec2,s_lat2,s_pot2%smrho,k1,k2,k3)

        goto 10

C --- Read restart file ---
      elseif (outs(1:4) == 'rsb ' .or. outs(1:4) == 'rsa ' .or.
     .        outs(1:3) == 'rs' .or. outs(1:5) == 'rsfa ') then
        if (outs(1:5) == 'rsfa ') then
          call rdovfa(nbask,nspeck,s_site1,s_spec1,s_lat,s_pot,s_ham,qbg)
          fileid = 'rsedit overlap fa densities'
        else
          call word(outs,2,j1,j2)
          if (outs(1:4) == 'rsa ') then
            if (j2 < j1) then
              fn = 'rsta'
            else
              fn = outs(j1:j2)
            endif
            ifi = fopna(fn,-1,1)
            lbin = .false.
          else
            if (j2 < j1) then
              fn = 'rst'
            else
              fn = outs(j1:j2)
            endif
            ifi = fopna(fn,-1,4+1)
            lbin = .true.
          endif
          call info0(10,0,-1,' rsedit: read restart file "'//trim(fn)//'"')
          k = iorsf(1,s_ctrl,s_site1,s_spec1,s_lat,s_pot,s_bz,fileid,
     .      nbask(1),natk(1),nspeck(1),xx,nit,lbin,ifi)
          call fclose(ifi)
          if (k < 0) call rx('abort : failed to read rst file')
        endif
        call dfratm(s_site1,s_spec1,8,1,nbask(1),s_pot1%rhat)
        call info5(10,0,0,'%9pfile has nbas = %i  nat = %i  '//
     .    'nspec = %i',nbask(1),natk(1),nspeck(1),0,0)
        havers = .true.
        lnsave = .true.
        call info2(2,1,0,
     .    ' Density from %?#n==0#rs file:#overlapping free atoms:#',
     .    isw(outs(1:5) == 'rsfa '),0)
        ngabc = s_lat%nabc
        call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
        call prsed2(nbask(1),nsp,s_site1,s_spec1,s_lat,s_pot%smrho,k1,k2,k3)

C --- Write atomic densities to file ---
      elseif (outs(1:7) == 'wrhoat ') then
        if (havers) then
          call prsed3x(nbask(1),nsp,s_site1,s_spec1)
        else
          call info0(0,0,0,' no rst file read ... nothing writtten')
        endif

C --- Symmetrize charge density ---
      elseif (outs(1:7) == 'symrho ') then
        if (havers) then
          call symrho(s_site,s_spec,s_lat,3,s_pot%smrho,s_pot%rhat,[xx],[xx],[xx])
        else
          call info0(0,0,0,' no rst file read ... nothing symmetrized')
        endif

C --- Merge spins ---
      elseif (outs(1:6) == 'onesp ') then

        if (.not. havers) then
          call info0(0,0,0,' no rst file read ... nothing changed')
          goto 98
        endif

        call prsed3x(nbask(1),nsp,s_site1,s_spec1)

C --- Make density for supercell  ---
      elseif (outs(1:5) == 'scell') then

        n123 = 0
        dc2 = outs(6:6)
        if (wordsw(outs,dc2,'mshdat','',j1) > 0) then ! Make supercell of file on regular mesh
          if (wordsw(outs,dc2,'fn','=',j1) <= 0) then ! dimensions of original mesh data
            call info0(0,0,0,' mshdat requires argument fn= ...'); goto 98
          endif
          call nwordg(outs,1,dc//dc2//' ',1,j1,j2); fn = outs(j1+1:j2); ifi = fopng(fn,-1,1)
          if (wordsw(outs,dc2,'n123','=',j1) <= 0) then ! dimensions of original mesh data
            call info0(0,0,0,' mshdat requires argument n123= ...'); goto 98
          endif
          if (a2vec(outs,len_trim(outs),j1,2,', '//dc2,3,3,3,ix,n123) /= 3) goto 98
          if (wordsw(outs,dc2,'m123','=',j1) <= 0) then  ! dimensions of new mesh data
            call info0(0,0,0,' mshdat requires argument m123= ...'); goto 98
          endif
          if (a2vec(outs,len_trim(outs),j1,2,', '//dc2,3,3,3,ix,m123) /= 3) goto 98
          allocate(meshnew(n123(1),n123(2),n123(3)*2),meshdat(m123(1),m123(2),m123(3)))
          call zrdm3(ifi,meshdat,m123(1),m123(2),m123(3),m123(1),m123(2),m123(3))
C         call zprm3('meshdat',0,meshdat,m123(1),m123(2),m123(3))
          call info2(10,0,0,' read mesh data with %sx%3i divisions, interpolating to %sx%3i',m123,n123)
        endif

C   ... read new site information from site file
        call words(outs,nw)
        if (nw /= 2) then
          call info0(0,0,0,' mshdat requires 2nd arg (site file) ...')
          goto 98
        endif
        call word(outs,2,j1,j2)
C       Read alat,plat,nbas of supercell
        j = iosite(7000,3d0,0,outs(j1:j2),ifi,s_ctrl%spid,alatx,plax,nbx,
     .    nspec,xx,xx,xx,xx,xx,xx,xx)
        call info5(10,0,0,' scell has %i sites with lattice parameter'//
     .    ' %d = %d*alat%N plax %9:1;6,6d%N plat %9:1;6,6d',
     .    nbx,alatx,alatx/alat,plax,plat)
        call dinv33(plat,1,qlat,vol)  ! Make vol
        call dinv33(plax,1,qlx,volx)  ! Make qlx,volx

        call info2(10,0,0,' Volume ratio %d',abs(volx/vol),0)
        if (abs(dble(nbx)/nbas-abs(volx/vol)) > 1d-5) then
          write (stdo,13)
   13     format(' supercell has incommensurate lattice ... skip scell')
          goto 98
        endif
C       Read pos,ips,imap of supercell
        allocate(posx(3,nbx),ipsx(nbx))
        j = iosite(8000,3d0,0,outs(j1:j2),ifi,s_ctrl%spid,alatx,plax,nbx,
     .    nspec,posx,xx,xx,xx,ipsx,xx,xx)
        allocate(imap(nbx)); call iinit(imap,nbx)
        if (ioibmap(ifi,nbx,imap) /= nbx) then
          call info0(1,0,0,' rsedit could not read map data from site file ... '//
     .      ' nothing done')
          goto 10
        endif

C   ... Fill out supercell-specific parts of s_lat2
        deallocate(s_site2)
        allocate(s_site2(nbx))
        s_lat2%alat = alatx
        s_lat2%plat = plax
        s_lat2%plat0 = plax
        deallocate(s_site2)
        allocate(s_site2(nbx))
        call lattic(s_lat2,s_ctrl,s_site2)

        if (n123(1) > 0) then

          call chgmsh(13,1,plat,m123(1),m123(2),m123(3),m123(1),m123(2),m123(3),meshdat,
     .                     plax,n123(1),n123(2),n123(3),n123(1),n123(2),n123(3),meshnew)

          call zprm3('$',0,meshnew,n123(1),n123(2),n123(3))
          call rx0(' data written to file out')
        endif


        if (.not. havers) then
          call info0(2,0,0,' no rst file read ... nothing changed')
          goto 98
        endif

        call info0(20,1,0,
     .    ' configure supercell lattice and mesh parameters')

C       Map smooth mesh density
        if (dlength(9,alatx*plax-alat*plat,1) > 1d-6) then

          if (s_lat%gmax <= 0) then
          call info0(0,0,0,' scell requires you specify gmax ... sorry')
          goto 98
          endif
          ngabc2 = 0
          call mshsiz(alatx,plax,0,s_lat%gmax,ngabc2,ngmx)
          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1x,k2x,k3x)
          s_lat2%nabc = ngabc2
          allocate(s_pot2%smrho(k1x*k2x*k3x,nsp))

          call chgmsh(10,nsp,plat,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .      s_pot%smrho,plax,ngabc2(1),ngabc2(2),ngabc2(3),k1x,k2x,k3x,
     .      s_pot2%smrho)

        endif

C   ... Copy local site data
        nbask(2) = nbx
        natk(2) = 0
        do  ibx = 1, nbx
          ib = imap(ibx)
          is = s_site(ib)%spec
          if (s_spec(is)%lmxa > -1) then
            natk(2) = natk(2) + 1
          endif
          s_site2(ibx)%pnu = s_site1(ib)%pnu
          s_site2(ibx)%pz  = s_site1(ib)%pz
          s_site2(ibx)%v0 => s_site1(ib)%v0
          s_site2(ibx)%v1 => s_site1(ib)%v1
          s_site2(ibx)%rho1 => s_site1(ib)%rho1
          s_site2(ibx)%rho2 => s_site1(ib)%rho2
          s_site2(ibx)%rhoc => s_site1(ib)%rhoc

          s_site2(ibx)%pos(1:3) = posx(1:3,ibx)
          s_site2(ibx)%pos0(1:3) = posx(1:3,ibx)
          s_site2(ibx)%force = s_site1(ib)%force
          s_site2(ibx)%vel = s_site1(ib)%vel
          s_site2(ibx)%spec = s_site1(ib)%spec
          s_site2(ibx)%class = s_site1(ib)%class
          s_site2(ibx)%clabel = s_site1(ib)%clabel

        enddo

        fn = 'rst2'
C       ifi = fopna(fn,-1,0); lbin = .false.
        ifi = fopna(fn,-1,4); lbin = .true.
        fileid = 'supercell autogenerated by lmf'
        k = iorsf(1,s_ctrl,s_site2,s_spec1,s_lat2,s_pot2,s_bz,
     .    fileid,nbx,natk(2),nspeck(2),xx,nit,lbin,-ifi)
        call rx0('wrote restart file rst2')

C --- Insert a floating orbital at a site ---
C     addflt site-index spec-index posx posy posz
      elseif (outs(1:7) == 'addflt ') then

        if (.not. havers) then
          call info0(0,0,0,' no rst file read ... nothing changed')
          goto 98
        endif

C       Read which entry to insert, and coordinates
        call words(outs,nw)
        if (nw /= 6) goto 98
        call word(outs,2,j1,j2)
        i = j1-1
        i = a2vec(outs,len(outs),i,2,', ',2,3,2,ix,nlst(1))
        if (i <= 0 .or. nlst(1) > nbask(1)+1) then
          call info0(0,0,0,'%6p... site index out of range: '//outs)
          goto 98
        endif
C       Must point to a valid species
        if (i <= 0 .or. nlst(2) > nspeck(1)) then
          call info0(0,0,0,'%6p... spec index out of range: '//outs)
          goto 98
        endif
        lmxa = s_spec1(nlst(2))%lmxa
        if (lmxa >= 0) then
          call info2(0,0,0,'%6p... spec %i not a floating orbital',
     .      nlst(2),0)
          goto 98
        endif

        i = j1-1
        i = a2vec(outs,len(outs),i,4,', ',2,3,5,ix,xxv)
        if (i /= 5) then
          call info0(0,0,0,'%6p... invalid site positions')
          goto 98
        endif

        call rx('rsedit line 420')
C       call prsed3(-nlst(1),nbask(1),nsp,s_site1,s_spec1)
        s_site1(nlst(1))%spec = nlst(2)
        s_site1(nlst(1))%pos =  xxv(3:5)
        s_site1(nlst(1))%pos0 = xxv(3:5)
        call dpzero(xxv,20)
        s_site(nlst(1))%vel = xxv(1:3)
        s_site(nlst(1))%pnu(1:n0,1:2) = reshape(xxv,[10,2])
        s_site(nlst(1))%pz(1:n0,1:2) = reshape(xxv,[10,2])

C --- Append a species from 2nd density to first density ---
      elseif (outs(1:8) == 'addspec ') then

        if (.not. havers) then
          call info0(0,0,0,' no rst file read ... nothing changed')
          goto 98
        endif

C       Read which species in 2nd density to append
        call words(outs,nw)
        if (nw /= 2) goto 98
        call word(outs,2,j1,j2)
        i = 0
        i = a2vec(outs(j1:j2),len(outs(j1:j2)),i,2,' ',1,1,1,ix,nlst(1))
        if (i <= 0 .or. nlst(1) < 0 .or. nlst(1) > nspeck(2)) then
        call info0(0,0,0,'%6p... invalid species index : '//outs(j1:j2))
        goto 98
        endif

        if (.not. have2s .and. nlst(1) > 0) then
          call info0(0,0,0,' no 2nd density ... nothing changed')
          goto 98
        endif
C       i = nint(sspec(1))
        if (nlst(1) > 0) then
          nspeck(1) = nspeck(1)+1
          s_spec1(nspeck(1)) = s_spec2(nlst(1))
        else
          nspeck(1) = nspeck(1)+1
          spid = 'FLT'
          s_spec1(nspeck(1))%name = spid
          s_spec1(nspeck(1))%z = 0
          s_spec(nspeck(1))%lmxa = -1

        endif

        spid = s_spec1(nspeck(1))%name
        call info2(10,0,0,' added species %i, id = '//spid,nspeck(1),0)

C --- Remove a site from the table and remake a charge density ---
C     Local densities sites not removed are preserved.
C     Interstitial density ni is modified to preserve neutrality
C     this is accomplished by adding the the difference delta ni
C     arising from a Mattheis construction, with the site removed,
C     and without.
C     Let Qi = int(ni) and QL = site charge
C     Qi_c0, QL_c0 = charges read from rst file, before substitution
C     Qi_m0, QL_m0 = charges from Mattheis constr, before substitution
C     Qi_mf, QL_mf = charges from Mattheis constr, after substitution
C     Qi_cf, QL_cf = charges after substitution
C     QL_cf and QL_c0 are the same, except one site is removed.
C     QL_cf is constructed to satisfy charge neutrality.
C     In general Qi + QL = 0  =>  Qi_cf + QL_cf = 0
C     Also Qi_cf = Qi_c0 + alfa*(Qi_mf - Qi_m0) = -QL_cf
C     => alfa = (-Qi_c0 - QL_cf)/(Qi_mf - Qi_m0)
C     Finally ni_mf = ni_m0 + alfa*(ni_mf - ni_m0)
C     QLs = net charge at site to be removed

      elseif (outs(1:7) == 'rmsite ' .or. outs(1:7) == 'chsite ') then

        if (.not. havers) then
          call info0(0,0,0,' no rst file read ... nothing changed')
          goto 98
        endif

C       Find indices to sites to remove or replace
        lrmsite = outs(1:7) == 'rmsite '
        lchsite = outs(1:7) == 'chsite '
        call word(outs,2,j1,j2)
        if (j2 < j1) goto 98
        call numsyv(n)
        call lodsyv('nbas',0,dble(nbask(1)),i)
        i = 0
        j = 2 ; if (lrmsite) j=1
        m = a2vec(outs(j1:),len(outs(j1:)),i,2,', '//dc,3,3,j,ix,nlst)
        if (m /= j) goto 98
        call clrsyv(n)
        do  n = 1, m
          if (nlst(n) <= 0 .or. nlst(n) > nbask(1)) then
          call info0(0,0,0,'%6p... invalid site index : '//outs(j1:j2))
          goto 98
          endif
        enddo
        j1 = j1+i-1
        ib = nlst(1); ibx = nlst(2)
        vol = s_lat%vol

C   ... Lump to add to sm rho from removal of site ib
C       Make a smoothed version of rho1-rho2, on a mesh
C       Integrated density QLs should be close to q1(ib)-q2(ib)
C       Charge too lumpy ... skip this branch
C        ng = s_lat%ng; allocate(cg1(ng)); call dpzero(cg1,2*ng)
C        n = 2
C        call rhgcmp(2,ib,ib,s_site1,s_spec1,s_lat,s_pot%rhat,n,ng,cg1)
C        allocate(smrho(k1,k2,k3,1))
C        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cg1,smrho)
C        call fftz3(smrho,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,1,0,1)
C        call mshint(vol,1,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
C     .    smrho,QLs,xx)
C

C       QLc = charge that will be actually lost
        call rhomom(2,ib,ib,s_site1,s_spec1,s_pot%rhat,QLc,xx)

        is = s_site1(ib)%spec
        if (s_spec1(is)%nxi <= 0) then
          call info0(0,0,0,' no coefficients for atom density')
          goto 98
        endif

C   ... Smooth interstitial density from site ib
        ng = s_lat%ng
        allocate(cg1(ng*nsp),smrho(k1,k2,k3,nsp))
        call ovlpfa(s_site,s_lat,s_ham,ib,ib,s_spec1(is)%nxi,n0,
     .    s_spec1(is)%exi,s_spec1(is)%chfa,s_spec1(is)%rsmfa,ng,ng,
     .    s_lat%gv,cg1)
        call gvputf(ng,nsp,s_lat%kv,k1,k2,k3,cg1,smrho)
        call fftz3(smrho,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,nsp,0,1)
        call mshint(vol,1,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .    smrho,QLs,xx)

C   ... Show charges before removing or changing site
        call mshint(vol,nsp,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .    s_pot%smrho,Qi,xx)
        call rhomom(2,1,nbask(1),s_site1,s_spec1,s_pot1%rhat,QL,xx)
        call info5(10,1,0,
     .    ' Starting charges:   Qi = %;6,6d '//
     .    ' Sum QL = %;6,6d  tot = %;6,6d '//
     .    ' QL(ib) = %;6,6d  Qfa(ib) = %;6,6d',
     .      Qi,QL,Qi+QL,QLc,QLs)

        if (lrmsite) then !Remove ib site density
          call prsed3(ib,nbask(1),s_site1)
          natk(1) = natk(1)-1 ! One fewer augmented atom
          nbask(1) = nbask(1)-1   ! One fewer augmented atom

        else                    !substitute ib site density with ibx density
          is = s_site1(ib)%spec
          if (s_spec1(is)%lmxa /= -1) natk(1) = natk(1)-1 ! One fewer augmented atom
          is = s_site1(ibx)%spec
          if (s_spec1(is)%lmxa /= -1) natk(1) = natk(1)+1 ! One extra augmented atom
          call prsed3(ib,-ibx,s_site1) ! Substitute the density
        endif
        call dfratm(s_site1,s_spec1,8,1,nbask(1),s_pot1%rhat) ! Reassociate s_rhat

C   ... Lump to add to sm rho from addition of site density ibx
C       Charge too lumpy ... skip this branch
C        if (.not. lrmsite) then
CC         Make a smoothed version of rho1-rho2, on a mesh
CC         Integrated density QLs should be close to q1(ib)-q2(ib)
C          call dpzero(cg1,2*ng)
C          n = 2
C          call rhgcmp(2,ib,ib,s_site1,s_spec1,s_lat,s_pot%rhat,n,ng,
C     .      cg1)
C          call gvputf(ng,1,s_lat%kv,k1,k2,k3,cg1,smrho)
C          call fftz3(smrho,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,1,0,1)
C          call mshint(vol,1,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
C     .      smrho,QLs,xx)
C
CC         QLc = charge that will be actually gained
C          call rhomom(2,ib,ib,s_site1,s_spec1,s_pot%rhat,QLc,xx)
C
CC     ... Add scaled lump to smrho to conserve total charge
C          call daxpy(2*k1*k2*k3*nsp,-QLc/QLs,smrho,1,s_pot%smrho,1)
C
C        endif

C   ... Lump to add to sm rho from addition of new site density ib
        if (lrmsite) then

        call mshint(vol,nsp,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .    s_pot%smrho,Qi,xx)
        call rhomom(2,1,nbask(1),s_site1,s_spec1,s_pot1%rhat,QL,xx)
        call info5(10,1,0,
     .    ' After substitution: Qi = %;6,6d '//
     .    ' Sum QL = %;6,6d  tot = %;6,6d '//
     .    ' Old Qfa(ib) = %;6,6d',
     .      Qi,QL,Qi+QL,QLs,0)

        else

C        if (s_spec1(is)%nxi <= 0) then
C          call info0(0,0,0,' no coefficients for atom density')
C          goto 98
C        endif

C   ... Smooth interstitial density from new site ib
        is = s_site1(ib)%spec
        ng = s_lat%ng
        allocate(smrho2(k1,k2,k3,nsp))
        call ovlpfa(s_site,s_lat,s_ham,ib,ib,s_spec1(is)%nxi,n0,
     .    s_spec1(is)%exi,s_spec1(is)%chfa,s_spec1(is)%rsmfa,ng,ng,
     .    s_lat%gv,cg1)
        call gvputf(ng,nsp,s_lat%kv,k1,k2,k3,cg1,smrho2)
        call fftz3(smrho2,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,nsp,0,1)
        call mshint(vol,1,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .    smrho2,QLs,xx)
C       Overwrite smrho with smrho(new ib) - smrho(old ib)
        call dscal(2*k1*k2*k3*nsp,-1d0,smrho,1)
        call daxpy(2*k1*k2*k3*nsp,1d0,smrho2,1,smrho,1)
        deallocate(smrho2)

C   ... Show charges after removing or changing site
        call mshint(vol,nsp,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .    s_pot%smrho,Qi,xx)
        call rhomom(2,1,nbask(1),s_site1,s_spec1,s_pot1%rhat,QL,xx)
        call info5(10,1,0,
     .    ' After substitution: Qi = %;6,6d '//
     .    ' Sum QL = %;6,6d  tot = %;6,6d '//
     .    ' Qfa(new ib) = %;6,6d',
     .      Qi,QL,Qi+QL,QLs,0)

C       Remake QLs, now shift in FA charge
        call mshint(vol,1,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .    smrho,QLs,xx)

        endif

C   ... Add sm FA site site density (or difference) to istl density
        call info2(10,0,0,' Add %;6,6d * rhos(ib) to smrho ',
     .    -(Qi+QL)/QLs,0)
        call daxpy(2*k1*k2*k3*nsp,-(Qi+QL)/QLs,smrho,1,s_pot%smrho,1)

C   ... Printout of final density
        call mshint(vol,nsp,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .    s_pot%smrho,Qi,xx)
        call rhomom(2,1,nbask(1),s_site1,s_spec1,s_pot1%rhat,QL,xx)
        call info5(10,0,0,
     .    ' Final charges: Qi = %;6,6d '//
     .    ' Sum QL = %;6,6d  tot = %;6,6d ',
     .      Qi,QL,Qi+QL,0,0)

       deallocate(cg1,smrho)

C --- Verbosity ---
      elseif (outs(1:5) == 'verb ') then

C       Read which entry to insert, and coordinates
        call word(outs,2,j1,j2)
        if (j2 < j1) goto 98
        i = j1-1
        i = a2vec(outs,len(outs),i,2,', ',2,3,2,ix,j)
        call pshpr(j)
        call shopr

C --- Show ---
      elseif (outs(1:5) == 'show ') then
        if (havers) then

          vol = s_lat%vol
          call mshint(vol,nsp,ngabc(1),ngabc(2),ngabc(3),k1,k2,k3,
     .      s_pot%smrho,Qi_c0,xx)
          call rhomom(2,1,nbask(1),s_site1,s_spec1,s_pot1%rhat,QL_c0,xx)
          call info5(10,1,0,' Interstitial charge: %;6,6d'//
     .      '  Local charge: %;6,6d  sum: %;6,6d',
     .      Qi_c0,QL_c0,Qi_c0+QL_c0,0,0)

          call info2(2,1,0,' Density from rs file:  '//
     .      '%i sites, %i species',nbask(1),nspeck(1))
          call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
          call prsed2(nbask(1),nsp,
     .      s_site1,s_spec1,s_lat,s_pot%smrho,k1,k2,k3)
        endif

        if (have2s) then
          call info2(2,1,0,' Second Density:  '//
     .      '%i sites, %i species',nbask(2),nspeck(2))

          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
          call prsed2(nbask(2),nsp,
     .      s_site2,s_spec2,s_lat2,s_pot2%smrho,k1,k2,k3)
        endif

C ... Spin averaging of densities, or  Exchange l=0 parts of
C     set object 1|2 list1 list2
C     exch 1|2 site1 site2  or exch 1|2 flip site1 site2
      elseif (outs(1:4) == 'set ' .or. outs(1:5) == 'exch ') then

      lexch = 0
      if (outs(1:5) == 'exch ') lexch = 1
      call words(outs,nw)
      if (nw < 4) goto 98

C     Switches for exchange
      if (lexch /= 0) then

C       Determine whether first or 2nd density
        call word(outs,2,j1,j2)
        if (outs(j1:j2) == '1') then
          lexch = 1
        elseif (outs(j1:j2) == '2') then
          lexch = 2
        else
          call info0(0,0,0,' 2rd argument must be 1 or 2 in '//
     .      trim(outs))
          goto 98
        endif

C       Check for flip
        call word(outs,3,j1,j2)
        if (outs(j1:j2) == 'flip') then
          if (nw /= 5) goto 98 ! Must be 5 arguments in this case
          isync = 4
          nw = 3
        else
          isync = 2
          nw = 2
        endif
        im = 0

C       Get 1st and 2nd site
        call numsyv(n)
        do  k = 1, 2
          call word(outs,nw+k,j1,j2)
          call lodsyv('nbas',0,dble(nbask(lexch)),i)
          i = 0
          i = a2vec(outs(j1:j2),len(outs(j1:j2)),i,2,' ',1,1,1,ix,
     .      nlst(k))
          if (i <= 0 .or. nlst(k) > nbask(lexch)) then
            call info0(0,0,0,
     .        '%6p... invalid site index : '//outs(j1:j2))
            goto 98
          endif
        enddo
        call clrsyv(n)

C       Carry out operation on 1st or second density
        if (lexch == 1) then
          ix(1) = nbask(1); ix(2) = nbask(1)
          if (havers) then
            call prsed1(im*10+isync,1,nlst,ix,nsp,
     .        s_site1,s_spec1,s_site1,s_spec1,
     .        s_pot%smrho,ngabc(1),ngabc(2),ngabc(3),
     .        s_pot%smrho,ngabc(1),ngabc(2),ngabc(3),
     .        fac,lfail)
            if (lfail) then
              call info0(0,0,0,' Densities failed to synchronize')
              goto 98
            endif
          else
            call info0(0,0,0,' no rst data ... cannot do: '//trim(outs))
            goto 10
          endif
        else
          ix(1) = nbask(2); ix(2) = nbask(2)
          if (have2s) then
            call prsed1(im*10+isync,1,nlst,ix,nsp,
     .        s_site2,s_spec2,s_site2,s_spec2,
     .        s_pot2%smrho,ngabc2(1),ngabc2(2),ngabc2(3),
     .        s_pot2%smrho,ngabc2(1),ngabc2(2),ngabc2(3),
     .        fac,lfail)
          else
            if (lfail) then
              call info0(0,0,0,' Densities failed to synchronize')
              goto 98
            endif
            goto 10
          endif
        endif

C     Switches for set
      else
        call word(outs,2,j1,j2)
        im = -1
        if (outs(j1:j2) == 'all') im = 30 + 4+2+1
        if (outs(j1:j2) == 'n') im = 31
        if (outs(j1:j2) == 'v') im = 4
        if (outs(j1:j2) == 'pnu') im = 2
        if (outs(j1:j2) == 'ntrue') im = 1
        if (outs(j1:j2) == 'asite') im = 7
        if (im == -1) then
          call info0(0,0,0,'     ... mode "'//outs(j1:j2)//'" not recognized')
          goto 98
        endif

C       asite implies site index as fifth argument
        isite = 0  ! If no specific site
        if (outs(j1:j2) == 'asite') then
          if (nw < 5) goto 98
          call word(outs,5,j1,j2)
          i = 0
          i = a2vec(outs(j1:j2),len(outs(j1:j2)),i,2,', ',2,3,1,ix,isite)
          if (i /= 1) goto 98
        endif

C       Fourth argument same for each density
        call word(outs,4,j1,j2)
        if (outs(j1:j2) == 'zers') then
          k = im*10+0
        elseif (outs(j1:j2) == 'zerq') then
          k = im*10+1
        elseif (outs(j1:j2) == 'flip') then
          k = im*10+2
        elseif (outs(j1:j2) == 'onesp') then
          k = im*10+4
        else
          call info0(0,0,0,' failed to recognize 4th argument in: '//trim(outs))
          goto 98
        endif

C       Carry out operation on 1st or second density
        call word(outs,3,j1,j2)
        if (outs(j1:j2) == '1') then
          if (havers) then
            call prsed4(k,isite,nbask(1),nspeck(1),nsp,
     .        s_site1,s_spec1,s_lat,s_pot%smrho,k1,k2,k3)
          else
            call info0(0,0,0,' no rst data ... cannot do: '//trim(outs))
            goto 10
          endif
        elseif (outs(j1:j2) == '2') then
          if (have2s) then
            call prsed4(k,isite,nbask(2),nspeck(2),nsp,s_site2,s_spec2,s_lat2,
     .        s_pot2%smrho,k1,k2,k3)
          else
            call info0(0,0,0,' no 2nd density... cannot do: '//
     .        trim(outs))
            goto 10
          endif
        else
          call info0(0,0,0,' failed to recognize 3rd argument in '//
     .      trim(outs))
          goto 98
        endif
      endif

      lnsave = .true.

C ... Linear combination of densities
      elseif (outs(1:4) == 'add ') then
        if (.not. (havers .and. have2s)) goto 97

        call words(outs,nw)
        if (nw < 3) goto 98
        call word(outs,2,j1,j2)
        im = -1
        if (outs(j1:j2) == 'asite') im = 8
        if (outs(j1:j2) == 'all') im = 38
        if (outs(j1:j2) == 'n') im = 31
        if (outs(j1:j2) == 'v') im = 4
        if (outs(j1:j2) == 'pnu') im = 2
        if (outs(j1:j2) == 'ntrue') im = 1
        if (outs(j1:j2) == 'rhoc') im = 9  ! Not fully implemented
        if (im == -1) then
          call info0(0,0,0,
     .      '     ... mode "'//outs(j1:j2)//'" not recognized')
          goto 98
        endif
        call word(outs,3,j1,j2)
        if (outs(j1:j2) == 'sync') then
          isync = 0
C         Require 3 words or 5 words
          if (im > 8 .and. nw /= 3 .or.
     .        im <= 8 .and. nw /= 5) goto 98
          nw = 3
        else
          isync = 1
C         Require 4 words or 6 words
          if (im > 8 .and. nw /= 4 .or.
     .        im <= 8 .and. (nw /= 4 .and. nw /= 6)) goto 98
          j = 0
          j = a2vec(outs(j1:),len(outs(j1:)),j,4,', ',2,-3,2,ix,fac)
          if (j /= 2) then
            call word(outs,4,j,j2)
            call info0(0,0,0,'     ... Could not parse arguments : '//outs(j1:j2))
            goto 98
          endif
        endif
        call numsyv(n)
        if (nw == 6) then ! get site lists
          do  k = 1, 2
            call word(outs,4+k,j1,j2)
            call lodsyv('nbas',0,dble(nbask(1)),i)
            if (k == 2) call lodsyv('nbas',0,dble(nbask(2)),i)
C           call shosyv(0,0,0,6)
            call mkils0(outs(j1:j2),nlst(k),j)
            if (nlst(k) <= 0) then
              call info0(0,0,0,'%6p... Bad or null list : '//outs(j1:j2))
              goto 98
            endif
            if (k == 1) then
              allocate(ilst(nlst(k),2))
              call mkilst(outs(j1:j2),nlst(k),ilst)
              call imxmn(nlst(k),ilst,1,l,m)
            elseif (k == 2) then
              if (nlst(2) == 1) then
                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
                call ivset(ilst(1,2),2,nlst(1),ilst(1,2))
              elseif (nlst(2) /= nlst(1)) then
                call info2(0,0,0,'%5p... 2nd list: '//outs(j1:j2)//
     .            ' has %i entries when %i are needed',nlst(2),
     .            nlst(1))
                goto 98
              else
                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
              endif
C             Check that elements in list are within range
              call imxmn(nlst(k),ilst(1,k),1,l,m)
            endif
            if (l < 1 .or.  k == 1 .and. m > nbask(1)
     .        .or.  k == 2 .and. m > nbask(2)) then
              call info0(0,0,0,
     .          '     ... Entry in list out of range: '//outs(j1:j2))
              goto 98
            endif
          enddo
        else
          if (nbask(1) /= nbask(2)) then
            call info0(0,0,0,
     .        '%10fadda requires 1-1 correspondence between sites'//
     .        ' in 1st and 2nd densities')
            goto 98
          endif
          nlst(1) = nbask(1)
          nlst(2) = nbask(1)
          allocate(ilst(nlst(1),2))
          do  i = 1, nlst(1)
            ilst(i,1:2) = i
          enddo
        endif
        call clrsyv(n)

        call prsed1(im*10+isync,nlst,ilst,nbask,nsp,
     .    s_site1,s_spec1,s_site2,s_spec2,
     .    s_pot%smrho,ngabc(1),ngabc(2),ngabc(3),
     .    s_pot2%smrho,ngabc2(1),ngabc2(2),ngabc2(3),
     .    fac,lfail)
        if (lfail) then
          call info0(0,0,0,' Densities failed to synchronize')
        elseif (isync == 0) then
          call info0(0,0,0,' Densities synchronized')
        elseif (isync == 1) then
          call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
          call info0(2,1,0,' New density from rs file:')
          call prsed2(nbask(1),nsp,
     .      s_site1,s_spec1,s_lat,s_pot%smrho,k1,k2,k3)
        endif

        deallocate(ilst)
        lnsave = .true.
        goto 10

C ... Save site file
      elseif (outs(1:6) == 'ssave ' .or. outs(1:7) == 'ssavex ') then

        call word(outs,2,j1,j2)
        if (j2 < j1) then
          print *, 'ssave requires file name as second argument'
          goto 98
        endif
        allocate(slabl(nspeck(1)))
        do  i = 1, nspeck(1)
          slabl(i) = s_spec(i)%name
        enddo
        if (iosits(1000*(2+4+8)+1,3d0,0,outs(j1:j2),ifi,slabl,alat,plat,
     .    nbask(1),nspeck(1),s_spec1,s_site1) < 0)
     .    call rx('failed to write site file')
        deallocate(slabl)

C ... Save rst file
      elseif (outs(1:5) == 'save ' .or. outs(1:6) == 'savea ') then
        lbin = outs(1:5) == 'save '
        call word(outs,2,j1,j2)
        if (j2 >= j1) fn = outs(j1:j2)
        if (lbin) then
          if (j2 < j1) fn = 'rst'
          ifi = fopna(fn,-1,4)
        else
          if (j2 < j1) fn = 'rsta'
          ifi = fopna(fn,-1,0)
        endif
        call info0(10,0,-1,
     .    ' rsedit: write restart file "'//trim(fn)//'"')
        k = iorsf(1,s_ctrl,s_site1,s_spec1,s_lat,s_pot,s_bz,fileid,
     .    nbask(1),natk(1),nspeck(1),xx,nit,lbin,-ifi)

        lnsave = .false.

C ... abort
      elseif (outs(1:2) == 'a ') then
        call rx0('aborting rs editor ... no file written')

C ... quit
      elseif (outs(1:2) == 'q '. or. outs(1:5) == 'quit ') then
        if (lnsave) then
          print '('' rs file not saved ... really quit?'')'
          read(*,'(a150)') outs
          call locase(outs)
          if (.not. (outs(1:1) == 'y' .or. outs(1:1) == 'q'))
     .      goto 10
        endif
        call rx0('exit rs editor')

C ... help
      elseif (outs == '?') then
        print 310
        print 311
        print 312
        print 313
        print 314
  310   format(
     .    ' Select one of these options:'/
     .  t4,'rs [fn]',t15,
     .    'read density and associated parameters from binary restart ',
     .    'file.'/t15,'Use "fn" for file name, ',
     .    'or "rst.ext" if "fn" is not supplied'/
     .  t4,'rsa [fn]',t15,
     .    'read density and associated parameters from ascii restart ',
     .    'file.'/t15,'Use "fn" for file name, ',
     .    'or "rst.ext" if "fn" is not supplied'/
     .  t4,'rsb [fn]',t15, 'same as "rs"'/
     .  t4,'rsfa [fn]',t15,
     .    'generate density and associated parameters from atm ',
     .    'file.'/t15,'Use "fn" for file name, ',
     .    'or "atm.ext" if "fn" is not supplied')
  311   format(/
     .  t4,'read [fn]',t15,
     .    'read 2nd density from binary restart file.  Use "fn" for ',
     .    'file name,'/t15,'or "rst.ext" if "fn" is not supplied'/
     .  t4,'reada [fn]',t15,
     .    'read 2nd density from ascii restart file.  Use "fn" for ',
     .    'file name,'/t15,'or "rsta.ext" if "fn" is not supplied'/
     .  t4,'readfa [fn]',t15,' same as rsfa, but for 2nd density'/
     .  t4,'readb [fn]',t15, 'same as "read"'//

     .  t4,'wrhoat',t15, 'Write atomic densities to rhoMT.{1..nbas}'//

     .  t4,'symrho',t15, 'Symmetrize site densities'//

     .  t4,'scell fn',t15, 'Create rst file for supercell'/t15,
     .    'Specify supercell lattice with site file "fn"'//

     .  t4,'rmsite #1',t15, 'remove site #1 from site data'/

     .  t4,'chsite #1,#2  copy site #2 density into site #1'//

     .  t4,'addspec #',t15,
     .    'append species # from 2nd density to 1st density.'/t15,
     .    'If # is zero, append a generic floating orbital species'/

     .  t4,'addflt #1 #2  posx posy posz'/
     .    t15,'insert a floating orbital at site #1.'/t15,
     .    '#2 = species index; it must point to a floating orbital.'//

     .  t4,'show',t15, 'Show summary information about densities'//
     .  t4,'...The following affect the charge or magnetic moment',
     .    ' of a density.'/
     .  t4,'set all|allnc|n 1|2  zers|zerq|flip|onesp'/
     .  t4,'set v|ntrue|pnu 1|2  zers|zerq|flip|onesp'/
     .  t4,'set asite       1|2  zers|zerq|flip|onesp site#'/
     .  t14,'* The first argument specifies which objects are affected.'
     .    /
     .  t16,'"all" applies to all quantities,'/
     .  t16,'"allnc" to all quantities excluding the core,'/
     .  t16,'"n" to the three density types,'/
     .  t16,'"v" to the spherical potential, ',
     .    '"ntrue" to the true local density,'/
     .  t16,'"pnu" to the linearization parameters'/
     .  t14,'* The second argument refers to the first or second ',
     .    'density'/
     .  t14,'* Third argument = zers:  zeros spin part of density:   ',
     .    'n+ - n- = 0'/
     .  t14,'  Third argument = zerq:  zeros charge part of density: ',
     .    'n+ + n- = 0'/
     .  t14,'  Third argument = flip:  exchange n+, n- '/
     .  t14,'  Third argument = onesp: condense to nsp=1')

  312   format(/
     .  t4,'...The following put a linear combination of the ',
     .    'two densities into the first, or'/
     .  t4,'checks whether the two densities are well enough ',
     .    'synchronized to make it possible'/
     .  t4,'add all|allnc|n',t26,'sync'/
     .  t4,'add v|ntrue|pnu',t26,'sync lst1 lst2'/
     .  t4,'add all|allnc|n',t26,'fac1 fac2'/
     .  t4,'add asite|v|ntrue|pnu',t26,'fac1 fac2 [lst1 lst2]'/
     .  t14,'* The first argument''s meaning is the same as the ',
     .    'first argument to "set"'/
     .  t14,'* Modifier "sync" checks whether the densities are ',
     .  'compatible, making no changes.'/
     .  t14,'* If "sync" is missing, the objects specified will be'/
     .  t16,'replaced by fac1*[1st type] + fac2*[second type].'/
     .  t14,'* For the site-only cases "v" and "ntrue," ',
     .  'site list "lst1" is required.'/
     .  t16,'"lst2" is a matching site list for the 2nd density, OR'/
     .  t16,'"lst2" is single site (same site added to all sites)'/)

  313   format(
     .  t4,'...The following exchanges the l=0 parts of two site',
     .    ' densities '/
     .  t4,'exch 1|2',t22,'[flip] site1 site2'/
     .    )

  314   format(
     .    t4,'save  [fn]',t15,'saves restart data in ',
     .    'binary restart file (name="rst" unless fn supplied)'/
     .    t4,'savea [fn]',t15,'saves restart data in ',
     .    'ascii restart file (name="rsta" unless fn supplied)'/
     .    t4,'ssave  fn',t15,'saves site data in file fn'/
     .    t4,'q',t15,'to quit the editor'/
     .    t4,'a',t15,'to abort')

      else
        print '(1x,''unrecognized option: '',a)', trim(outs)

      endif
      goto 10

   97 call info0(0,0,0,'%10pBoth densities must be supplied before'
     .  //'invoking this command ')
   98 call info0(0,0,0,' rsedit:  improper usage of '//trim(outs)//
     .  ' ... nothing done')
      goto 10

      end

      subroutine prsed1(mode,nlst,ilst,nbask,nsp,
     .  s_site1,s_spec1,s_site2,s_spec2,
     .  smrho,k11,k21,k31,smrh2,k12,k22,k32,fac,lfail)
C- Add fac*2nd density to density
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 check synchronization
Ci         :1 add the density
Ci         :2 swap l=0 part of 2 local densities.  nlst should be 1
Ci         :4 like 2, but exchange spin in the swap
Ci         :10s digit:
Ci         :0 work with true density only
Ci         :1 work with pnu
Ci         :2 work with v0,v1
Ci         :4 work with n1,n2,n0
Ci         :Any combination of the above is allowed
Ci         :8 same as 7, but exclude core density
Ci         :9 core densities only

Ci         :10s digit: true site quantities
Ci         :1 modify true density, incl core
Ci         :2 modify pnu
Ci         :4 modify v0,v1
Ci         :8 1+2+4, and also species core density
Ci         :9 species core densities only
Ci         :100s digit: interstitial quantities
Ci         :1 modify smooth density
Ci         :2 modify local smooth density
Ci         :Any combination is allowed

Ci   nlst  :number of site densities to which 2nd density is added
Ci   ilst  :list of indices to site densities
Ci         :ilst(:,1) is list of sites to be modified
Ci         :ilst(:,2) is list from which 2nd density is to be found
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci
Ci   site1 :struct for site-specific information; see routine usite
Ci   spec1 :struct for species-specific information; see routine uspec
Ci   orh1  :vector of offsets to arrays containing local site densities
Ci   site2 :Analog of site1 for 2nd density
Ci   spec2 :Analog of spec1 for 2nd density
Ci   orh2  :Analog of orh1 for 2nd density
Ci   fac   :Add fac * 2nd density into density
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Jun 08
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nlst,ilst(nlst,2),nbask(2),nsp,
     .  k11,k21,k31,k12,k22,k32
C     integer orh1(3,*),orh2(3,*)
      double precision fac(2)
      double complex smrho(k11,k21,k31,nsp),smrh2(k12,k22,k32,nsp)
      logical lfail
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site1(*),s_site2(*)
      type(str_spec)::  s_spec1(*),s_spec2(*)
C ... Local Dynamical arrays
      real(8), pointer  :: p_v01(:,:),p_v11(:,:),p_v02(:,:),p_v12(:,:)
C ... Local parameters
      logical listl,ln1,ln2,lc,lv0,lpnu,bittst,lcs
      integer il,ii,ib,jb,is,js,lmxa,lmxl,lmxa2,lmxl2,nr1,nr2,nrmx
      integer nglob,stdo,ipr,nlml1,nlml2,isync,n0,isw
      integer, allocatable:: flgspc(:)
      parameter (nrmx=5001,n0=10)
      double precision rmt1,rmt2,a1,a2,rofi(nrmx)
      double precision pnu(n0,2),pnu2(n0,2),pnz(n0,2),pnz2(n0,2)
      character spid1*8,spid2*8
      logical lfails

      stdo = nglob('stdo')
      call getpr(ipr)
      isync = mod(mode,10)

      il  = mod(mode/10,10)
      if (il == 8) il = 8+7
      ii  = mod(mode/100,10)

      listl = bittst(ii,1)      ! Modify smrho
      ln2   = bittst(ii,2)      ! Modify n2

      ln1   = bittst(il,1)      ! Modify n1
      lc    = ln1               ! Modify core
      lcs   = bittst(il,8)      ! Modify species core
      lpnu  = bittst(il,2)      ! Modify pnu
      lv0   = bittst(il,4)      ! Modify v0

      lfails = .false.
      if (listl) then
          call info0(0,1,-1,' Checking whether sm density meshes synchronize ...')
        if (k11 /= k12 .or. k21 /= k22 .or. k31 /= k32) then
          write(stdo,212) k11,k21,k31,k12,k22,k32
  212     format(' mesh mismatch:'/
     .      ' first  density n1 n2 n3 = ',3i5/
     .      ' second density n1 n2 n3 = ',3i5)
          lfails = .true.
          if (isync == 1) call rx('sm density mesh mismatch')
        else
          call info0(0,0,0,' ok')
        endif
        if (isync == 1) then
          call info2(0,1,-1,' ... Replace smooth density n0 by linear '
     .      //'combination %d*[n0] + %d*[2nd n0]',fac(1),fac(2))
          call lcrho(k11*k21*k31*2,nsp,1,1,fac(1),fac(2),smrho,smrh2)
        endif
      endif

      if (isync >= 2) then
        call info5(0,1,0,' ... Swap '//
     .    '%?#n==4#and spin flip ##'//
     .    '%?#n==0# true site# site #%-1j'//
     .    '%?#n==0# density n1##%-1j'//
     .    '%?#n>=4# densities##%-1j'//
     .    '%?#n>=8# (exclude core)##'//
     .    '%?#n>=1# pnu##',
     .    isync,il,mod(il,4),mod(il,8),mod(il,2))
        if (isync == 3 .or. isync > 4 .or. nlst > 1) then
          call info0(0,1,0,' ... illegal parameters (aborting)')
          return
        endif
      elseif (isync == 1) then
        call info8(0,1,0,' ... Replace site ('//
     .    '%?#n#true-rho ##'//
     .    '%?#n#sm-rho ##'//
     .    '%?#n#pnu ##'//
     .    '%?#n#pot ##'//
     .    '%?#n#spec-core ##'//
     .    '%a)  by linear combination %d*[1st] + %d*[2nd]',
     .    isw(ln1),isw(ln2),isw(lpnu),isw(lv0),isw(lcs),fac(1),fac(2),0)
      elseif (isync == 0) then
        call info0(0,1,0,' Checking whether radial meshes synchronize ...')
      endif
      call info2(0,0,0,' site 1 : %n:1,3i',nlst,ilst(1,1))
      call info2(0,0,0,' site 2 : %n:1,3i',nlst,ilst(1,2))

  200 format(/' ... Improper synchronization between radial meshes'/
     .  '  site1',8x,'site2',8x,'nr1  nr2',3x,'rmt1   rmt2',4x,' a1    a2')

      lfail = .false.
      allocate(flgspc(nbask(1)))
      call iinit(flgspc,nbask(1))
      do  il = 1, nlst

        ib = ilst(il,1)
        if (ib > nbask(1))
     .    call rxi('illegal site index for density:',ib)
        jb = ilst(il,2)
        if (jb > nbask(2))
     .    call rxi('illegal index for 2nd density:',jb)

        is = s_site1(ib)%spec
        pnu = s_site1(ib)%pnu
        pnz = s_site1(ib)%pz
        p_v01 => s_site1(ib)%v0
        p_v11 => s_site1(ib)%v1
        js = s_site2(jb)%spec
        pnu2 = s_site2(jb)%pnu
        pnz2 = s_site2(jb)%pz
        p_v02 => s_site2(jb)%v0
        p_v12 => s_site2(jb)%v1
        lmxa = s_spec1(is)%lmxa
        lmxl = s_spec1(is)%lmxl
        lmxa2 = s_spec2(js)%lmxa
        lmxl2 = s_spec2(js)%lmxl
        spid1 = s_spec1(is)%name
        spid2 = s_spec2(is)%name

        if (lmxa == -1) cycle

        a1 = s_spec1(is)%a
        nr1 = s_spec1(is)%nr
        rmt1 = s_spec1(is)%rmt
        a2 = s_spec2(js)%a
        nr2 = s_spec2(js)%nr
        rmt2 = s_spec2(js)%rmt

C       Radial meshes must match
        if (abs(a1-a2) > 1d-6 .or. nr1 /= nr2 .or.
     .      abs(rmt1-rmt2) > 1d-6) then

          if (.not. lfail) write(stdo,200)
          lfail = .true.
          write(stdo,202) ib,spid1,jb,spid2,
     .      nr1,nr2,rmt1,rmt2,a1,a2
  202     format(i4,1x,a,i4,1x,a,2i5,1x,2f7.4,1x,2f6.3)
          if (isync == 1) call rx('site density mismatch')

        endif
        if (isync == 0) cycle

        nlml1 = (lmxl+1)**2
        nlml2 = (lmxl2+1)**2
        call radmsh(rmt1,a1,nr1,rofi)

C       True local density
        if (ln1) then
          if (isync < 2) then
            call lcrho(nr1,nsp,nlml1,nlml2,fac(1),fac(2),
     .        s_site1(ib)%rho1,s_site2(jb)%rho1)
          else
            call swrho(isync/4,nr1,nsp,1,nlml1,nlml2,
     .      s_site1(ib)%rho1,s_site2(jb)%rho1)
          endif
        endif

C       Smooth local density and site core density
        if (ln2 .and. isync < 2) then
          call lcrho(nr1,nsp,nlml1,nlml2,fac(1),fac(2),
     .        s_site1(ib)%rho2,s_site2(jb)%rho2)
        elseif (ln2) then
          call rx('not ready for this swap')
        endif
        if (lc .and. isync < 2) then
          call lcrho(nr1,nsp,1,1,fac(1),fac(2),
     .      s_site1(ib)%rhoc,s_site2(jb)%rhoc)
        elseif (lc) then
          call rx('not ready for this swap')
        endif

C       Potentials v0,v1
        if (lv0) then
          if (isync < 2) then
            call lcrho(nr1,nsp,1,1,fac(1),fac(2),p_v01,p_v02)
            call lcrho(nr1,nsp,1,1,fac(1),fac(2),p_v11,p_v12)
          else
            call rx('not ready for this swap')
          endif
        endif

C       Pnu modified
        if (lpnu) then
          if (isync < 2) then
            call lcrho(n0,nsp,1,1,fac(1),fac(2),pnu,pnu2)
            s_site1(ib)%pnu = pnu
            call lcrho(n0,nsp,1,1,fac(1),fac(2),pnz,pnz2)
            s_site1(ib)%pz = pnz
          else
            call rx('not ready for this swap')
          endif
        endif

C       Species core density modified if im ge 4
        if (lcs .and. flgspc(is) == 0) then
          if (isync < 2) then
            call lcrho(nr1,nsp,1,1,fac(1),fac(2),
     .        s_spec1(is)%rhoc,s_spec2(js)%rhoc)
            if (lcs .and. fac(1) == 0.and.fac(2) == 1) then
              s_spec1(is)%ctail = s_spec2(js)%ctail
              s_spec1(is)%etail = s_spec2(js)%etail
              s_spec1(is)%stc   = s_spec2(js)%stc
              s_spec1(is)%lfoca = s_spec2(js)%lfoca
              s_spec1(is)%rfoca = s_spec2(js)%rfoca
              s_spec1(is)%nxi   = s_spec2(js)%nxi
              s_spec1(is)%exi   = s_spec2(js)%exi
              s_spec1(is)%chfa  = s_spec2(js)%chfa
              s_spec1(is)%rsmfa = s_spec2(js)%rsmfa
            endif
          else
            call rx('not ready for this swap')
          endif
          flgspc(is) = 1
        endif

      enddo

      if (lfail) then
        call info0(0,0,0,' ... Radial meshes failed to synchronize')
      elseif (isync == 0) then
        call info0(0,0,0,' ... Radial meshes synchronized')
C      elseif (isync == 1) then
C        call info0(0,0,0,' ... Linear combination successful')
      endif

      lfail = lfail .or. lfails

      end

      subroutine prsed2(nbas,nsp,s_site,s_spec,s_lat,smrho,k1,k2,k3)
C- Print out density information
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspec
Ci   sspec :struct for species-specific information; see routine uspec
Ci   nspec2:number of species for second density
Ci   sspec2
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   28 Jun 08 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,k1,k2,k3
      double complex smrho(k1,k2,k3,nsp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Local parameters
      integer nrmx,nlmx
      parameter (nrmx=5001, nlmx=64)
      integer ib,is,nr,lmxa,lmxl,lmxb,nlml,ipr,nglob,stdo,n1,n2,n3,
     .  intopt
      double precision z,qc,rg,a,rmt,rsma
      double precision pi,srfpi,y0,rofi(nrmx),rwgt(nrmx),qv1,qv2,am1,am2
      double precision ddot,smq,sum2,alat,plat(3,3),qlat(3,3),vol,smag
      character spid*8
C     double complex smpot(k1,k2,k3,2)

      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')
      call getpr(ipr)
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0 = 1d0/srfpi
      n1 = k1
      n2 = k2
      n3 = k3
      alat = s_lat%alat
      plat = s_lat%plat
      call dinv33(plat,1,qlat,vol)
      vol = dabs(vol)*(alat**3)
      if (nsp == 2) then
        call daxpy(k1*k2*k3*2,1d0,smrho(1,1,1,2),1,smrho,1)
      endif
      call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smq,sum2)
      smag = 0
      if (nsp == 2) then
        call daxpy(k1*k2*k3*2,-1d0,smrho(1,1,1,2),1,smrho,1)
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smag,sum2)
        smag = 2*smag - smq
      endif

C      if (nsp == 1) write(stdo,100)
C      if (nsp == 2) write(stdo,201)

      call info5(0,0,-1,' Smooth density: vol=%,2;2d  n123 = %i %i %i',
     .  vol,k1,k2,k3,0)
      call info5(0,0,0,'  qsm = %,6;6d%?#n==2#  mag. mom = %,6;6d##',
     .  smq,nsp,smag,0,0)

      if (nsp == 1) write(stdo,200)
      if (nsp == 2) write(stdo,201)
  200 format(/' Site parameters and charges:'/
     .  '   site',10x,'z     rmt     nr   a  lmxl   qtrue',6x,'q2')
  201 format(/' Site parameters and charges:'/
     .  '   site',10x,'z     rmt     nr   a  lmxl   qtrue',
     .  5x,'mtrue',5x,'q2        m2')

      do  ib = 1, nbas

        is = s_site(ib)%spec
        z = s_spec(is)%z
        qc = s_spec(is)%qc
        rg = s_spec(is)%rg
        spid = s_spec(is)%name
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        rsma = s_spec(is)%rsma
        lmxa = s_spec(is)%lmxa
        lmxl = s_spec(is)%lmxl
        lmxb = s_spec(is)%lmxb
        if (lmxa == -1) cycle
        nlml = (lmxl+1)**2

        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)

        call splrho(0,nsp,nr,nlml,s_site(ib)%rho1,s_site(ib)%rho2,s_site(ib)%rhoc)
        qv1   = srfpi*ddot(nr,rwgt,1,s_site(ib)%rho1,1)
        qv2   = srfpi*ddot(nr,rwgt,1,s_site(ib)%rho2,1)
        call splrho(1,nsp,nr,nlml,s_site(ib)%rho1,s_site(ib)%rho2,s_site(ib)%rhoc)
        am1   = 2*srfpi*ddot(nr,rwgt,1,s_site(ib)%rho1,1) - qv1
        am2   = 2*srfpi*ddot(nr,rwgt,1,s_site(ib)%rho2,1) - qv2

        if (nsp == 1) then
          write(stdo,202) ib,spid,z,rmt,nr,a,lmxl,qv1,qv2
        else
          write(stdo,202) ib,spid,z,rmt,nr,a,lmxl,qv1,am1,qv2,am2
        endif
  202   format(i4,1x,a,f6.1,f9.5,i5,f6.3,i3,4f10.5)

   10   continue
      enddo


      end

      subroutine prsed3(ibx,nbx,s_site)
C- Remove or substitute one element in site data
C ----------------------------------------------------------------------
Ci   ibx   : remove or substitute site ibx
Ci   nbx   : if nbx>0, shift s_site(ib+1) to s_site(ib) for ib>=ibx
Ci         : thus removing ibx from s_site
Ci         : nbx should be the size of the basis, nbas
Ci         : if nbx<0, replace s_site(ibx) with s_site(-nbx)
Cu Updates
Cu   30 Aug 13 Rewritten
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ibx,nbx
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ib2,kb

      ib2 = nbx-1  ! if removing a site from the table
      if (nbx < 0) ib2 = ibx

      do  ib = ibx, ib2

        kb = ib+1  ! if removing a site from the table
        if (nbx < 0) kb = -nbx ! if substituting one site

        if (nbx > 0) then
        s_site(ib)%pos = s_site(kb)%pos
        endif
C       s_site(ib)%class = s_site(ib)%class ! not meaningful
        s_site(ib)%spec = s_site(kb)%spec
        s_site(ib)%pnu = s_site(kb)%pnu
        s_site(ib)%pz = s_site(kb)%pz

        s_site(ib)%v0 => s_site(kb)%v0
        s_site(ib)%v1 => s_site(kb)%v1
        s_site(ib)%rho1 => s_site(kb)%rho1
        s_site(ib)%rho2 => s_site(kb)%rho2
        s_site(ib)%rhoc => s_site(kb)%rhoc

      enddo
      end

      subroutine prsed3x(nbas,nsp,s_site,s_spec)
C- Write atomic densities to disk
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer nrmx,nlmx
      parameter (nrmx=5001, nlmx=64)
      integer ib,is,nr,lmxa,lmxl,nlml
      double precision a,rmt,xx
      double precision rofi(nrmx)

      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        lmxl = s_spec(is)%lmxl

        if (lmxa == -1) cycle
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        nlml = (lmxl+1)**2
        call radmsh(rmt,a,nr,rofi)
        call wrhomt(1,'rhoMT.','density',ib,s_site(ib)%rho1,xx,rofi,
     .    nr,nlml,nsp)
      enddo

      end

      subroutine prsed4(mode,isite,nbas,nspec,nsp,
     .  s_site,s_spec,s_lat,smrho,k1,k2,k3)
C- Remove spin or charge part of atomic densities
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu v0 v1
Co     Stored:     pnu
Co     Allocated:  *
Cio    Elts passed:rho1 rho2 rhoc
Cio    Passed to:  prsed2
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z qc rg name a nr rmt rsma lmxa lmxl lmxb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rhoc
Cio    Passed to:  prsed2
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  prsed2
Ci Inputs
Ci   mode  :0 remove spin part
Ci         :1 remove charge part
Ci         :2 flip spins
Ci         :4 merge into one spin
Ci         :10s digit: true site quantities
Ci         :1 modify true density, incl core
Ci         :2 modify pnu
Ci         :4 modify v0,v1
Ci         :100s digit: interstitial quantities
Ci         :1 modify smooth density
Ci         :2 modify local smooth density
Ci         :Any combination is allowed
Ci   isite :If nonzero, modify only this site
Ci   nbas  :number of atoms
Ci   nspec :number of species
Ci   smrho :smooth density on uniform mesh
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Jun 08
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,isite,nbas,nspec,nsp,k1,k2,k3
      double complex smrho(k1,k2,k3,nsp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Local parameters
      logical listl,ln1,ln2,lc,lv0,lpnu,lspec,bittst
      integer nrmx,nlmx,n0
      parameter (nrmx=5001, nlmx=64, n0=10)
      integer ib,is,nr,lmxa,lmxl,lmxb,nlml,ipr,nglob,stdo,ii,il,mode0,intopt,isisite
      double precision z,qc,rg,a,rmt,rsma,dglob,xx
      double precision pi,srfpi,y0,rofi(nrmx),rwgt(nrmx),pnu(n0,2),pnz(n0,2)
      real(8), pointer  :: p_v0(:,:),p_v1(:,:)
C     double precision dval


      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')
      call getpr(ipr)
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0 = 1d0/srfpi
      mode0 = mod(mode,10)
      il  = mod(mode/10,10)
      ii  = mod(mode/100,10)

      listl = bittst(ii,1)      ! Modify smrho
      ln2   = bittst(ii,2)      ! Modify n2

      ln1   = bittst(il,1)      ! Modify n1
      lc    = ln1               ! Modify core
      lpnu  = bittst(il,2)      ! Modify pnu
      lv0   = bittst(il,4)      ! Modify v0
      lspec = ln1               ! modify core in species

      if (mode0 == 0 .and. nsp == 1) return
      if (mode0 == 2 .and. nsp == 1) return

      if (listl) then
C       Replace smrho-up,smrho-down with smrho, smrho(up-down)
        call splrho(30,nsp,k1*k2*k3*2,1,smrho,smrho,smrho)
        if (mode0 /= 4) then
        call prsed5(mod(mode0,4),smrho,k1*k2*k3*2,1,nsp)
C       Replace smrho, smrho(up-down) with smrho-up,smrho-down
        call splrho(31,nsp,k1*k2*k3*2,1,smrho,smrho,smrho)
        endif
      endif

C      if (nsp == 1) write(stdo,200)
C      if (nsp == 2) write(stdo,201)
C  200 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
C     .  6x,'q2')
C  201 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
C     .  5x,'mtrue',5x,'q2        m2')

      do  ib = 1, nbas

        if (isite /= 0 .and. ib /= isite) cycle
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
        pnz = s_site(ib)%pz
        p_v0 => s_site(ib)%v0
        p_v1 => s_site(ib)%v1
        z = s_spec(is)%z
        qc = s_spec(is)%qc
        rg = s_spec(is)%rg
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        rsma = s_spec(is)%rsma
        lmxa = s_spec(is)%lmxa
        lmxl = s_spec(is)%lmxl
        lmxb = s_spec(is)%lmxb
        if (lmxa == -1) cycle
        nlml = (lmxl+1)**2

        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)

C       Replace rho-up,rho-down with rho, up-down
        call splrho(0,nsp,nr,nlml,
     .    s_site(ib)%rho1,s_site(ib)%rho2,s_site(ib)%rhoc)

C       Modify local densities
        if (mode0 /= 4) then
          if (ln1) call prsed5(mod(mode0,4),s_site(ib)%rho1,nr,nlml,nsp)
          if (lc) call prsed5(mod(mode0,4),s_site(ib)%rhoc,nr,1,nsp)
          if (ln2) call prsed5(mod(mode0,4),s_site(ib)%rho2,nr,nlml,nsp)
        endif

C       Restore rho-up,rho-down from rho, up-down
        if (mode0 /= 4) then
          call splrho(1,nsp,nr,nlml,s_site(ib)%rho1,s_site(ib)%rho2,s_site(ib)%rhoc)
        endif

C       Modify potentials v0,v1: m
        if (lv0) then
          call splrho(30,nsp,nr,1,p_v0,xx,xx)
          call prsed5(mod(mode0,4),p_v0,nr,1,nsp)
          call splrho(31,nsp,nr,1,p_v0,xx,xx)
          call splrho(30,nsp,nr,1,p_v1,xx,xx)
          call prsed5(mod(mode0,4),p_v1,nr,1,nsp)
          call splrho(31,nsp,nr,1,p_v1,xx,xx)
        endif

C       Modify pnu
        if (lpnu) then
C         call prmx('pnu',pnu,n0,n0,2)
          call splrho(30,nsp,n0,1,pnu,xx,xx)
          call prsed5(mod(mode0,4),pnu,n0,1,nsp)
          call splrho(31,nsp,n0,1,pnu,xx,xx)
          s_site(ib)%pnu = pnu
          call splrho(30,nsp,n0,1,pnz,xx,xx)
          call prsed5(mod(mode0,4),pnz,n0,1,nsp)
          call splrho(31,nsp,n0,1,pnz,xx,xx)
          s_site(ib)%pz = pnz
C         call prmx('pnu',pnu,n0,n0,2)
        endif

      enddo

C     Modify species core density
      if (lspec) then
        isisite = 0; if (isite /= 0) isisite = s_site(isite)%spec
        do  is = 1, nspec
          if (s_spec(is)%lmxa == -1) cycle
          if (isisite /= 0 .and. isisite /= is) cycle
          nr = s_spec(is)%nr
          call splrho(30,nsp,nr,1,s_spec(is)%rhoc,xx,xx)
          if (mode0 /= 4) then
          call prsed5(mod(mode0,4),s_spec(is)%rhoc,nr,1,nsp)
          call splrho(31,nsp,nr,1,s_spec(is)%rhoc,xx,xx)
          endif
        enddo
      endif

C ... Reduce number of spins from 2 to 1
      if (mode0 == 4) then
        nsp = 1
        is = dglob('nsp',1d0,1)
      endif

      call info0(1,1,0,' New density:')
      call prsed2(nbas,nsp,s_site,s_spec,s_lat,smrho,k1,k2,k3)

      end
      subroutine prsed5(mode,rho,nr,nlml,nsp)
C- Zero out spin part (mode=0) or charge part (mode=1) of rho
C  Density rho is passed as (up+dn,up-dn)
Ci   mode  :0 remove spin part
Ci         :1 remove charge part
Ci         :2 flip spins

      implicit none
      integer mode,nr,nlml,nsp
      double precision rho(nr,nlml,nsp)

C     Flip sign out up-dn
      if (mode == 2) then
        call dscal(nr*nlml,-1d0,rho(1,1,2),1)
        return
      endif
C     Swap up+dn,up-dn
      if (mode == 0) then
        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
      endif
C     Zero out channel 1 (up-dn if mode=0, otherwise up+dn)
      call dpzero(rho,nr*nlml)
C     Swap up+dn,up-dn
      if (mode == 0) then
        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
      endif
      end
C      subroutine prsed6(idup,plat,plat2,
C     .  nbas,nbas2,nsp,site1,spec1,orh1,
C     .  site2,spec2,orh2,
C     .  smrho,k11,k21,k31,
C     .  smrh2,k12,k22,k32)
CC- Duplicate structure in a supercell
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 check synchronization
CCi         :1 add the density
CCi         :2 swap l=0 part of 2 local densities.  nlst should be 1
CCi         :4 like 2, but exchange spin in the swap
CCi         :10s digit:
CCi         :0 modify true density only
CCi         :1 modify pnu
CCi         :2 modify v0,v1
CCi         :4 modify n1,n2,n0
CCi         :Any combination is allowed
CCi   nlst  :number of site densities to which 2nd density is added
CCi   ilst  :list of indices to site densities
CCi         :ilst(:,1) is list of sites to be modified
CCi         :ilst(:,2) is list from which 2nd density is to be found
CCi   nsp   :2 for spin-polarized case, otherwise 1
CCi
CCi   site1 :struct for site-specific information; see routine usite
CCi   spec1 :struct for species-specific information; see routine uspec
CCi   orh1  :vector of offsets to arrays containing local site densities
CCi   site2 :Analog of site1 for 2nd density
CCi   spec2 :Analog of spec1 for 2nd density
CCi   orh2  :Analog of orh1 for 2nd density
CCo Outputs
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   28 Jun 08
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer idup(3),nbas,nbas2,nsp,
C     .  k11,k21,k31,k12,k22,k32
C      integer orh1(3,*),orh2(3,*)
C      double precision site1(1),spec1(1),site2(1),spec2(1)
C      double precision plat(3,3),plat2(3,3)
C      double complex smrho(k11,k21,k31,nsp),smrh2(k12,k22,k32,nsp)
CC ... Local parameters
C      integer i1,i2,i3,ii(3)
C      integer nglob,stdo,ipr,ssize
C      double precision trans(3)
CC      integer il,ib,jb,is,js,lmxa,lmxl,lmxa2,lmxl2,nr1,nr2,nrmx
CC      integer nglob,stdo,ipr,nlml1,nlml2,isync,im,n0
CC      integer orhoc1,orhoc2
CC      integer ,allocatable:: flgspc(:)
CC      parameter (nrmx=5001,n0=10)
CC      double precision rmt1,rmt2,a1,a2,rofi(nrmx),pnu(n0,2),pnu2(n0,2)
CC      character spid1*8,spid2*8
CC      logical lfails
CC ... Heap
C      integer w(1)
C      common /w/ w
C
C      ssize = site1(1)
C      stdo = nglob('stdo')
C      call getpr(ipr)
C
C      do i3 = 1, idup(3)
C      ii(3) = i3
C      do i2 = 1, idup(2)
C      ii(2) = i2
C      do i1 = 1, idup(1)
C      ii(1) = i1
C        print *, ii
C      enddo
C      enddo
C      enddo
C
C      stop
C
C      end
