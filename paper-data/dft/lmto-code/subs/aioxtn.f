      subroutine aioxtn(s_ctrl,s_spec,s_site,s_lat,s_bz,slabl,pnuc,qnuc)
C- I/O moments and pp's for all classes, other formats
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: nbas nbasp nclass nl nspin lrs
Co     Stored:    nbas nbasp lrs
Co     Allocated: *
Cio    Elts Passed:lrs
Cio    Passed to: iorsa bcast_strx iinit mpibc1
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa rmt z coreh coreq name a nr lmxl kmxt p pz lfoca
Ci                qc idmod rsma lmxb kmxv rsmv rfoca ctail etail stc
Ci                nxi exi chfa rsmfa rhoc
Co     Stored:    name lmxa a nr rmt z lmxl kmxt p pz qc lfoca rfoca
Co                coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa rhoc
Co     Allocated: rhoc
Cio    Elts Passed:rhoc
Cio    Passed to: getq gtpcor iosits iorsa bcast_strx
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: class pos spec vel eula pl relax vshft ndelta delta
Ci                bxc cpawt omg omgn domg gc gcorr j0 pdos rho1 rho2
Ci                rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl
Ci                eqkkl sighh sighk sigkk tauhh tauhk taukk pihh pihk
Ci                pikk sighhx sighkx sigkkx tauhhx tauhkx taukkx pihhx
Ci                pihkx pikkx thet v0 v1 force pnu pz
Co     Stored:    pos spec vel eula pl relax vshft bxc cpawt omg omgn
Co                domg gc gcorr j0 pdos rho1 rho2 rhoc rho1x rho2x
Co                rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl sighh sighk
Co                sigkk tauhh tauhk taukk pihh pihk pikk sighhx sighkx
Co                sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx thet v0
Co                v1 pos0 force clabel pnu pz
Co     Allocated: v0 v1
Cio    Elts Passed:rhoc rho1 rho2
Cio    Passed to: iosits bcast_strx iorsa aioxt1 asars1
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat plat vol awald nkd nkq nabc ag bgv cg cy dlv gv
Ci                gvq indxcg ips0 istab jcg kv igv igv2 kv2 pos qlv
Ci                symgr
Co     Stored:    alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                jcg kv igv igv2 kv2 pos qlv symgr
Co     Allocated: *
Cio    Elts Passed:pos dlv qlv
Cio    Passed to: iorsa bcast_strx iinit mpibc1
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Stored:    ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: iorsa bcast_strx iinit mpibc1
Ci Inputs
Ci   slabl :vector of species labels
Ci   pnuc  :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci          This array is class-based
Ci   qnuc  :energy-weighted moments of the sphere charges
Ci          This array is class-based
Co Outputs
Co   *When importing from atom files, output are updated
Co    atom files.
Co   *When importing from rst files, output is pnuc,qnuc
Cl Local variables
Cl         :
Cb Bugs
Cb   When importing from atom files,
Cb   file importing for layer case (nbasp>nbas) fails to generate
Cb   proper Madelung potential.
Cr Remarks
Cr
Cu Updates
Cu   17 Jun 13 Completed f90 replacement of f77 pointers
Cu   25 Oct 11 Started migration to f90 structures
Cu   04 Nov 04 Upgrade of functionality
Cu   05 Mar 02 added to the restart file editor
Cu   01 Mar 02 First cut at import mode from rst files (ASA mode)
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      double precision pnuc(*),qnuc(*)
      character*8 slabl(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_lat)::   s_lat
      type(str_bz)::    s_bz
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      real(8), allocatable :: pos(:,:)
      integer, allocatable :: ipc(:),ipa(:)
      real(8), pointer :: pnus(:),qnus(:),wkl(:)
      real(8), allocatable :: qcs(:),qts(:)
      real(8), allocatable :: amoms(:)
      real(8), allocatable :: pos1(:)
      real(8), allocatable :: pos2(:)
      real(8), allocatable :: pos3(:)
      real(8), allocatable :: qnui(:)
      real(8), allocatable :: v(:)
      real(8), allocatable :: qtl(:)
      real(8), allocatable :: rmaxl(:)
      real(8), allocatable :: vrmaxl(:)
      real(8), allocatable :: vesl(:)
      real(8), allocatable :: mad(:)
C ... Local parameters
      real(8),parameter :: NULLR=-99999d0
      integer ic,ifi,fopna,fopnx,nr,lmx,nsp,idmod(9),nl,nbas,nbasp,
     .  nbasn,nclass,j1,j2,nkd,nkq,i,j,k,mode,stdo,lio,nspec,nati
      procedure(integer) :: nglob,iorsa,lgunit,iosits
      double precision a,rmax,vrmax(2),rhrmx(2),ves,dval,
     .  z,qc,dq,sumec,sumev,thrpv,ekin,utot,rhoeps,etot(26),
     .  pnl(2*9),qnl(3*2*9),pp(3*2*9),evrl(9*2),dsum,wk(10),xx
      double precision awald,alat,vol,plat(3,3),tau(3),emad,
     .  trumad,vmtz(2),dglob,eterms(22)
      integer ov,is,nbasi,a2vec,
     .  nspeci,nit,irs,n0,nlsts,nw,
     .  nlstd,ip,lsts(10000),lstd(10000),oqcs
      parameter (n0=10)
      character outs*80, clabl*8, alabel*4, fnam*80
      logical sw,cmdopt,lrel,lgen,lmom,lpp,lpot,lqtot,lves,a2bin
      logical aiogen,aiomom,aiopar,aiopot
      logical io4g,io4mm,io4pp,io4pt,io3smm,lstutt,lrs
      logical io4sg,io4smm,io4spp,io4spt
      logical io5sg,io5smm,io5spt
      logical io47um,io47up
      character*32 jobid

      nbas = s_ctrl%nbas
      nbasp = s_ctrl%nbasp
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      alat = s_lat%alat
      plat = s_lat%plat
      nspec = nglob('nspec')
      allocate(ipc(nbasp))
      call sitepack(s_site,1,nbas,'class',1,ipc,xx)
      alat = s_lat%alat
      vol = s_lat%vol
      stdo = lgunit(1)

      call dpzero(vrmax,2)
c     if charges at all sites are found, lqtot=.true.
      lqtot = .true.
      lves = .false.
      allocate(qtl(nclass)); call dpzero(qtl,nclass)
      allocate(rmaxl(nclass)); call dpzero(rmaxl,nclass)
      allocate(vrmaxl(nclass*2)); call dpzero(vrmaxl,nclass*2)
      allocate(vesl(nclass)); call dpzero(vesl,nclass)

      lrs = cmdopt('-rs',3,0,outs)

c   ... import from rst files
        if (lrs) then

          irs = s_ctrl%lrs

          call rx('aioxtn update this code')
C          integer, pointer :: ips(:)
C          allocate(ips(nbasp))
C          call spackv(10,'site spec',ssite,1,nbasp,ips)
C          real(8), pointer :: spot(:)
C          call defdr(ospot,1000)
C          allocate(spot(1000))
C          call upot(spot,-1,-000,0,0,0)
C          call redfrr(ospot,nint(dval(spot,1)))
C
C          real(8), pointer :: sbz(:)
C          call defdr(osbz,1000)
C          allocate(sbz(1000))
C          call ubz(sbz,-1,-000,0,0,0)
C          call redfrr(osbz,nint(dval(sbz,1)))
C
C          real(8), pointer :: slat(:)
C          call defdr(oslat,1000)
C          allocate(slat(1000))
C          call ulat(slat,-1,0,0,0,0)
C          call redfrr(oslat,nint(dval(slat,1)))
C
C          real(8), pointer :: pnus(:)
C          real(8), pointer :: qnus(:)
C          integer, pointer :: ipa(:)
C          call defi(oipa,nbasp)
C          allocate(pnus(n0*nsp*nbasp)); call dpzero(pnus,n0*nsp*nbasp)
C          allocate(qnus(3*n0*nsp*nbasp)); call dpzero(qnus,3*n0*nsp*nbasp)
C          allocate(ipa(nbasp))
C          call asars1(30,s_site,nbasp,nsp,nl,pnus,qnus,
C     .      ipc,ipa,pnuc,qnuc)
C          call rlse(oipa)


   99     print 301
  301     format(//
     .      ' Welcome to the restart file editor.')
  100     continue
          print '('' Option (enter ? for menu) : '',$)'
          outs = ' '
          read(*,'(a80)') outs
          call locase(outs)
          if (.false.) then
          elseif (outs == ' ') then
            print '('' Enter ''''q'''' to exit, ''''a'''' to abort'')'
            goto 100
          elseif (outs(1:1) == 'a') then
            call rx0('exiting aioxtn without rewriting rst file')
          elseif (outs(1:1) == 'q') then
            goto 199
          elseif (outs(1:1) == 'i') then
            goto 101
          elseif (outs(1:3) == 'dup') then
            call words(outs,nw)
            if (nw == 1) then
              print 315
  315         format(/' Enter the lower and upper range of sites ',
     .          'to duplicate: ',$)
              read(*,*) i,j
            elseif (nw == 2 .or. nw == 3) then
              ip = 3
              k = a2vec(outs,len(outs),ip,2,', ',2,-3,-2,lsts,lstd)
              if (k < 0) then
                print 306, outs(1:ip+1)
                goto 100
              elseif (k == 1) then
                i = lstd(1)
                j = i
              elseif (k == 2) then
                i = lstd(1)
                j = lstd(2)
              endif
            endif
            if (j < i) then
              print *, 'upper range < lower range ... skipping'
              goto 100
            endif

            nbasn = nbasp + j-i+1
C           duplicate ssite
            call rx('aioxtn update this code')


C           duplicate pnu
            k = n0*nsp
            wkl => pnus
            allocate(pnus(k*nbasn))
            call dcopy(k*nbasp,wkl,1,pnus,1)
            call dmscop(pnus,k,wkl,k,1,k,i,nbasp,1,j+1,1d0)
C           call aioxt3(k,wkl,nbasp,pnus,nbasn)
C           duplicate qnu
            k = 3*n0*nsp
            allocate(qnus(k*nbasn))
            wkl => qnus
            call dcopy(k*nbasp,wkl,1,qnus,1)
            call dmscop(qnus,k,wkl,k,1,k,i,nbasp,1,j+1,1d0)
C           call aioxt3(k,wkl,nbasp,qnus,nbasn)
C           Update nbas,nbasp
            nbas  = nbas  + j-i+1
            nbasp = nbasn
            i = dglob('nbasp',dble(nbasp),1)
            i = dglob('nbas',dble(nbas),1)
            s_ctrl%nbas = nbas
            s_ctrl%nbasp = nbasp
C           Update ipc,ips,z,lmxa
            call rx('aioxtn update this code')
C           call defi(oipc,nbasp) -> ipc
C            call spackv(10,'site class',ssite,1,nbasp,ipc)
C            integer, pointer :: ips(:)
C            allocate(ips(nbasp))
C            call spackv(10,'site spec',ssite,1,nbasp,ips)

            print '('' basis enlarged to'',i5,'' atoms'')', nbasp
            goto 100

          elseif (outs(1:5) == 'shift') then
            call words(outs,nw)
            if (nw == 1) then
              print 311
  311         format(/' Enter the lower and upper range of sites ',
     .          'to shift: ',$)
              read(*,*) i,j
            elseif (nw == 2 .or. nw == 3) then
              ip = 5
              k = a2vec(outs,len(outs),ip,2,', ',2,-3,-2,lsts,lstd)
              if (k < 0) then
                print 306, outs(1:ip+1)
                goto 100
              elseif (k == 1) then
                i = lstd(1)
                j = i
              elseif (k == 2) then
                i = lstd(1)
                j = lstd(2)
              endif
            endif
            if (j < i) then
              print *, 'upper range < lower range ... skipping'
              goto 100
            endif
            print 314
  314       format(/' Enter dx dy dz: ',$)
            read(*,*) tau
            allocate(pos(3,nbasp))
            call sitepack(s_site,1,nbas,'pos',3,xx,pos)
            call aioxt2(i,j,tau,pos)
            call sitepack(s_site,1,nbas,'-pos',3,xx,pos)
            deallocate(pos)
            goto 100

          elseif (outs(1:1) == 's') then


            call rx('need allocate lmx,z')
            allocate(qcs(nbasp),qts(nbasp),amoms(nbasp))
            call getq(nsp,n0,lmx,nbasp,z,pnus,qnus,
     .        s_ctrl%ics,s_spec,qcs,qts,amoms)

            allocate(pos(3,nbasp))
            call sitepack(s_site,1,nbas,'pos',3,xx,pos)
            allocate(pos1(nbasp),pos2(nbasp),pos3(nbasp))
            call dmscop(pos1,1,pos,3,1,1,1,nbasp,1,1,1d0)
            call dmscop(pos2,1,pos,3,2,2,1,nbasp,1,1,1d0)
            call dmscop(pos3,1,pos,3,3,3,1,nbasp,1,1,1d0)

            outs = ' Site   Z      x         y         z'//
     .        '         dq     amom'
            print *
C            call arrprt(outs,'%,4i%:-3,4;1d%:-3,4;4d%:-2,4;4d',
            call arrprt(outs,
     .        '%,4i%;6,1D%;10,5D%;10,5D%;10,5D%:-3,4;4d%:-2,4;4d',
     .        'Idddddd',nbasp,0,2,0,'  | ',xx,
     .        z,pos1,pos2,pos3,qts,amoms,xx)
            wk(1) = dsum(nbasp,qts,1)
            wk(2) = dsum(nbasp,amoms,1)
            call awrit2(' sum-of-sphere charges=%:,6;6d'//
     .        '  sum-of-sphere moments=%:,6;6d',' ',80,
     .        stdo,wk,wk(2))
            call rlse(oqcs)
            print *, ' '
            goto 100

          elseif (outs == 'wsite') then
            lio = 1000*(2+4+8+16+32) + 1
            if (iosits(lio,3d0,0,'site',ifi,slabl,alat,plat,nbasp,nspec,
     .        s_spec,s_site) < 0)
     .        call rx('failed to write ssite')
            goto 100

          elseif (outs == '?') then
            print 310
  310       format(/' Select one of these options:'/
     .        t4,'i',t10,'to import a portion of a restart file'/
     .        t4,'dup',t10,
     .                   'to duplicate some site data, enlarging basis'/
     .        t4,'s',t10,'to show some site data'/
     .        t4,'shift',t10,'to shift some site positions'/
     .        t4,'wsite',t10,'to write site file'/
     .        t4,'q',t10,'to quit the editor'/
     .        )
            goto 100

          else
            print '(1x,''unrecognized option: '',a10)', outs
            goto 100

          endif

C     ... Import some pnu,qnu from a rst file
  101     continue
          fnam = ' '
          print 302
  302     format(/' Enter the full name of a restart file ',
     .            '(enter nothing to exit): ',$)
          read(*,'(a80)') fnam
          if (fnam == ' ') goto 100
          ifi = fopnx(fnam,72,-1,-1)
          if (ifi /= 1) then
            print *, 'sorry, file doesn''t seem to exist.  try again.'
            goto 101
          endif
C         call pshpr(110)
          ifi = fopnx(fnam,72,0,-1)
C         call poppr

          i = iorsa(-1,s_ctrl,xx,xx,s_lat,xx,s_bz,
     .      jobid,nbasi,nati,nspeci,xx,nit,.false.,ifi)

          call word(fnam,1,j1,j2)
          print 303, fnam(j1:j2), nbasi, nspeci
  303     format(' file ',a,' contains nbas =',i4,' and nspec =',i4)

C         setup for file read ... much of the data will be discarded
          s_ctrl%lrs = 0

          call rx('aioxtn update this code: allocate ssite and sspec')

          do  i = 1, nspeci
            s_spec(i)%name = 'zz'
          enddo
          allocate(qnui(3*n0*nsp*nbasi))
          call dpzero(qnui,3*n0*nsp*nbasi)

          i = iorsa(0,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,jobid,nbasi,
     .      nati,nspeci,qnui,nit,.false.,ifi)

          s_ctrl%lrs = irs
          call fclose(ifi)


  110     print 304, fnam(j1:j2)
  304     format(/' Enter list of sites from input file ',a,
     .      ' (? for help) : ',$)
          read(*,'(a80)') outs

          if (outs == ' ') then
            goto 101

          elseif (outs == '?') then
          print 305
  305     format(/
     .      ' If you enter a single number, you can poke data from'/
     .      ' this site into a group of sites.'/
     .      ' If you enter a list, the list of destination sites'/
     .      ' you enter following this will have to be the same length'/
     .      ' as the list you enter now.'/
     .      ' See documentation for site-list syntax.')

          goto 110
          endif

          call mkils0(outs,nlsts,ip)
          if (nlsts < 0) then
            print 306, outs(1:ip+1)
  306       format(' Failed to parse string ',a,' ... try again.')
            goto 101
          endif
          call mkilst(outs,nlsts,lsts)
          print 307, outs(1:ip+1), nlsts
  307     format(' List ',a,' contains',i4,' sites.')


  120     print 308
  308     format(/' Enter list of destination sites : ',$)
          read(*,'(a80)') outs

          if (.false.) then
          elseif (outs == ' ' .or.  outs == '?') then
            print 312
  312       format(/' Enter ''q'' to quit'/
     .        ' Enter ''?'' to see this message'/
     .        ' Otherwise, enter a site-list ',
     .        ' (see documentation for syntax)')
            if (nlsts == 1)  print 313
  313       format(2x,' Because you imported only one site from the',
     .        ' input file,'/
     .        2x,' you can enter an implicit list with z=#')
            goto 120
          elseif (outs == 'q') then
            goto 99
          elseif (outs(1:2) == 'z=') then
            if (nlsts /= 1) then
              print *, 'Sorry, z= option allowed only when you import',
     .          ' from one site'
              goto 120
            endif
            i = 2
            if (.not. a2bin(outs,wk,4,0,' ',i,len(outs)-1)) then
              print *, 'Sorry, couldn''t parse string ... try again'
              goto 120
            endif
            nlstd = 0
            do  i = 1, nbasp
              if (dval(z,i) == wk(1)) then
                nlstd = nlstd + 1
                lstd(nlstd) = i
              endif
            enddo
            call word(outs,1,j1,j2)
            print 307, outs(j1:j2), nlstd
            if (nlstd == 0) then
              print *, 'Sorry, no sites with this z ... try again'
              goto 120
            endif

          else

            call mkils0(outs,nlstd,ip)
            if (nlstd < 0) then
              print 306, outs(1:ip+1)
              goto 101
            endif
            call mkilst(outs,nlstd,lstd)
            print 307, outs(1:ip+1), nlstd
          endif

          if (nlsts > 1 .and. nlsts /= nlstd) then
            print 309
  309       format(' Size of lists do not match ... try again.')
            goto 101
          endif


C         Poke imported site->pnu and qnui to pnus and qnus
          mode = 1
          if (nlsts == 1) mode = 0
          call aioxt1(mode,nlstd,lsts,lstd,nsp,nbasi,s_site,
     .      qnui,nbasp,pnus,qnus)

          goto 101


C     ... Cleanup and exit rst editor
  199     continue
C         call defi(oipa,nbasp)
          allocate(ipa(nbasp))
          call asars1(31,s_site,nbasp,nsp,nl,pnus,qnus,ipc,ipa,pnuc,
     .      qnuc)
          deallocate(ipa)
          call rx('no more ipc')

          return

        endif

c     call pshpr(110)
    5 continue
      print 357
  357 format(/' aioxtn: read gen   mom   ppar  pot')
      do  10  ic = 1, nclass
        lgen = .false.
        lmom = .false.
        lpp  = .false.
        lpot = .false.
        rhrmx(1) = 0
        rhrmx(2) = 0
        call r8tos8(s_ctrl%dclabl(ic),clabl)
        is = s_ctrl%ics(ic)
        lmx = s_spec(is)%lmxa
        rmax = s_spec(is)%rmt
        z = s_spec(is)%z
        lstutt = cmdopt('-3s ',4,0,outs) .or.
     .           cmdopt('-4s ',4,0,outs) .or.
     .           cmdopt('-47u',4,0,outs) .or.
     .           cmdopt('-5s ',4,0,outs)

c   ... stuttgart lmtoxx don't lower case, and no extension
        if (lstutt) then
          ifi = fopnx(clabl,2,0,-1)
        else
          ifi = fopna(clabl,-1,0)
        endif
        rewind ifi

        if (cmdopt('-4 ',4,0,outs) .or. cmdopt('-3s ',4,0,outs)) then
          if (cmdopt('-4 ',4,0,outs)) then
          lgen = io4g(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .      qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          lmom = io4mm(alabel,pnl,qnl,idmod,nl,lmx,nsp,rhrmx,vrmax,ifi)
          lpp  = io4pp(alabel,pp,ves,nl,lmx,nsp,ifi)
          else
          lgen = io4sg(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .      qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          lmom = io3smm(alabel,pnl,qnl,idmod,nl,lmx,nsp,rhrmx,vrmax,ifi)
          lpp  = io4spp(alabel,pp,ves,nl,lmx,nsp,ifi)
          endif
          allocate(v(nr*nsp))
          lpot = io4pt(nr,nsp,a,rmax,v,ifi)
          clabl = alabel
        elseif (cmdopt('-5 ',3,0,outs)) then
          call rx('aioxtn not ready for v5')
        elseif (cmdopt('-4s',3,0,outs)) then
          alabel = clabl
          lgen = io4sg(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .      qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          lmom = io4smm(alabel,pnl,qnl,evrl,idmod,nl,lmx,nsp,vrmax,ifi)
          lpp  = io4spp(alabel,pp,ves,nl,lmx,nsp,ifi)
          allocate(v(nr*nsp))
          lpot = io4spt(nr,nsp,a,rmax,v,ifi)
          call fclose(ifi)
          ifi = fopna(clabl,-1,0)
        elseif (cmdopt('-47u',4,0,outs)) then
          nr = 0
          alabel = clabl
c         aiogen same as 5s
          lgen = io5sg(a,etot,lmx,lrel,nr,nsp,qc,dq,rmax,z,ifi)
          ekin  = etot(2)
          utot  = etot(3)
          sumec = etot(5)
          sumev = etot(6)
          rhoeps = etot(10)
          thrpv = 0
          lmom = io47um(alabel,lmx,nl,nsp,pnl,idmod,qnl,vrmax,ifi)
c     ... for now this doesn't mean anything, so we leave it out.
          lpp = io47up(alabel,lmx,nl,nsp,pp,ves,ifi)
          allocate(v(nr*nsp))
c         aiopot same as 5s
          lpot = io5spt(nr,nsp,a,rmax,v,ifi)
          call fclr(clabl,ifi)
          ifi = fopna(clabl,-1,0)
        elseif (cmdopt('-5s',3,0,outs)) then
          nr = 0
c     ... lgen picks up for us a,nr
          alabel = clabl
          lgen = io5sg(a,etot,lmx,lrel,nr,nsp,qc,dq,rmax,z,ifi)
          ekin  = etot(2)
          utot  = etot(3)
          sumec = etot(5)
          sumev = etot(6)
          rhoeps = etot(10)
          thrpv = 0
          lmom = io5smm(alabel,idmod,lmx,nl,nsp,pnl,qnl,vrmax,ifi)
          lpp  = .false.
c     ... for now this doesn't mean anything, so we leave it out.
c          lpp  = io5spp(alabel,lmx,nl,nsp,pp,ves,ifi)
          allocate(v(nr*nsp))
          lpot = io5spt(nr,nsp,a,rmax,v,ifi)
          if (lpot .and. lves) then
            call dpscop(vesl,ves,1,ic,1,1d0)
            call daxpy(nr*nsp,-1d0,ves,0,v,1)
          endif
          call fclr(clabl,ifi)
c          call fclose(ifi)
c          ifi = fopna(clabl,-1,0)
        else
          call rx('aioxtn: need conversion -4, -5 or -3s, -4s or -5s')
        endif
        if (lgen) then
          call dvset(qtl,ic,ic,dq)
          call dvset(rmaxl,ic,ic,rmax)
        else
          lqtot = .false.
        endif
        print 358, clabl, lgen,lmom,lpp,lpot
  358   format(2x,a8,5(l6))

        if (lstutt) then
          call fclose(ifi)
          ifi = fopna(clabl,-1,0)
          rewind ifi
        endif
        rewind ifi
        if (lgen) sw = aiogen(clabl,z,rmax,lmx,nsp,lrel,nr,a,qc,
     .    dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,-ifi)
        if (lmom) sw = aiomom(clabl,pnl,qnl,[NULLR],idmod,nl,lmx,nsp,z,rhrmx,vrmax,(/0d0,0d0,1d0/),-ifi)
        if (lpp) sw = aiopar(clabl,0,pp,xx,ves,nl,lmx,nsp,-ifi)
        if (lpot) sw = aiopot(nr,nsp,a,rmax,-99d0,v,-ifi)
        if (lpot) call rlse(ov)
        call fclr(clabl,ifi)
   10 continue

c     sphere charges are known ... get madelung potential to shift pot
      if (lstutt .and. .not. lves) then
        clabl = 'no'
        if (lqtot) clabl = 'yes'

        print '(/'' ... were all sphere charges found? '',a)', clabl
        if (lqtot) then

        print '('' ... undo madelung shift in site potentials'')'
        alat = s_lat%alat
        vol = s_lat%vol
        awald = s_lat%awald
        nkd = s_lat%nkd
        nkq = s_lat%nkq
        allocate(mad(nbasp**2))
        call madmat(nbasp,s_lat%pos,awald,alat,vol,s_lat%dlv,nkd,
     .    s_lat%qlv,nkq,mad)
        vmtz(1) = 0
        call madpot(nbasp,1,nclass,s_ctrl%nrc,ipc,s_ctrl%dclabl,
     .    qtl,0d0,rhrmx,rmaxl,mad,xx,xx,0d0,.false.,vrmaxl,
     .    vesl,emad,trumad,vmtz,0,xx,xx,eterms)
        deallocate(mad)
        print '(/'' ... re-read files, shifting potentials by ves'')'
        lves = .true.
        goto 5
        endif
      endif


      end

      subroutine aioxt1(mode,nlst,lsts,lstd,nsp,nbasi,s_site,
     .  qnui,nbas,pnu,qnu)
C- Poke P,Q parameters from source list into destination list
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 import data is for one site only (lsts(1))
Ci         :1 import Pi,Qi(lsts(1..nlst) to Pd,Qd(lstd(1..nlst)
Ci   nlst  :number of sites for which to distribute imported P,Q
Ci   lsts  :list of sites from which to distribute imported P,Q
Ci   lstd  :list of sites to which to distribute imported P,Q
Ci   nbasi :number of sites for which imported data is available
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ssitei:struct for site-specific information; see routine usite
Ci         :It contains pnu for imported data
Ci   qnui  :energy-weighted moments of the sphere charges to import
Co Outputs
Co   pnu   :boundary conditions.  In this routine array is site-based
Co   qnu   :energy-weighted moments of the sphere charges
Co         :In this routine array is site-based
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Feb 02 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nlst,nsp,lsts(nlst),lstd(nlst),n0,nbasi,nbas
      parameter (n0=10)
      double precision qnui(3,n0,nsp,*)
      double precision pnu(n0,nsp,nbas),qnu(3,n0,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ibs,ibd
      double precision ploc(n0,2)

      print '(''      source  dest'')'
      do  ib = 1, nlst
        if (mode == 0) then
          ibs = lsts(1)
        else
          ibs = lsts(ib)
        endif
        ibd = lstd(ib)
        print 333, ib, ibs,ibd
  333   format(i4,2i6,$)
        if (ibs > nbasi) then
          print 334, 'ib(s)'
  334     format(' *** ', a,' exceeds max val ... copy suppressed')
        elseif (ibd > nbas) then
          print 334, 'dest'
        else
          print *, ' '
          ploc = s_site(ibs)%pnu
          call dcopy(n0*nsp,ploc,1,pnu(1,1,ibd),1)
          call dcopy(3*n0*nsp,qnui(1,1,1,ibs),1,qnu(1,1,1,ibd),1)
        endif
      enddo
      end

c ---- file i/o atomic data, v4 ---
      logical function io4g(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     . qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
c- file i/o for atomic general data.  ifi>0 for read, <0 for write
c ----------------------------------------------------------------
ci inputs
ci
co outputs
co
cr remarks
cr
c ----------------------------------------------------------------
      implicit none
c passed parameters
      integer lmx,nsp,nr,ifi
      double precision z,rmax,a,
     . qc,dq,vrmax(2),sumec,sumev,thrpv,ekin,utot,rhoeps,etot
      character*4 alabel
      logical lrel
c local parameters
      logical scat

      io4g = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'gen:',':',.true.)) return
        read(ifi,200) alabel,lmx,nsp,rmax,lrel,nr,a
        read(ifi,201) z,qc,dq,vrmax,rhoeps
        read(ifi,202) sumec,utot,ekin,sumev,etot,thrpv
        io4g = .true.
      else
        write(-ifi,'(''gen:'')')
        write(-ifi,100) alabel,lmx,nsp,rmax,lrel,nr,a
        write(-ifi,101) int(z),int(qc),dq,vrmax,rhoeps
        write(-ifi,102) sumec,utot,ekin,sumev,etot,thrpv
      endif
  100 format(3x,a4,'  lmx=',i1,'  nspin=',i1,'  rmax=',f9.6,'  rel=',l1,
     .       '  nr=',i4,'  a=',f5.3)
  101 format(3x,'z=',i2,  '  qc=',i2,'  qtot=',f9.6,'  vrmax=',2f9.6,
     .          '  rhoeps=',f12.6)
  102 format('   sumec=',f15.7,'  utot=',f15.7,'  ekin= ',f15.7/
     .       '   sumev=',f15.7,'  etot=',f15.7,'  thrpv=',f15.7)

  200 format(3x,a4,6x,i1,8x,i1,7x,f9.6,6x,l1,5x,i4,4x,f4.3)
  201 format(3x,2x,f2.0,5x,f2.0,7x,f9.6,8x,2f9.6,9x,f12.6)
  202 format(2(9x,f15.7,7x,f15.7,8x,f15.7/))

      end

      logical function io3smm(albl,pl,ql,idmod,nl,lmax,nsp,rh,vrmax,ifi)
c- file i/o for moments.  ifi>0 for read, <0 for write
c ----------------------------------------------------------------
ci inputs
ci   albl,nl,lmax,nsp
ci   ifi:  logical unit: positive for read, negative for write
ci   pl:   related to log derivatives (see atom)
ci   ql:   moments q (see atom) (ifi > 0)
co   rh,vrmax density and vxc at rmax (ifi>0)
co outputs
co   pl,ql:   moments q (see atom) (ifi < 0)
co   aiomom:true unless read error or category not found
co   rh,vrmax density and vxc at rmax (ifi<0)
cr remarks
cr   input for second spin defaults to that for first spin if not
cr   available for second spin.
cr   bug in reading vrmax for nsp=2
c ----------------------------------------------------------------
      implicit none
c passed parameters
      integer nl,lmax,nsp,ifi
      character*4 albl
      integer idmod(0:nl-1)
      double precision pl(0:nl-1,nsp),ql(3,0:nl-1,nsp),rh,vrmax(2)
c local parameters
      double precision pi,dl,xx(4)
      integer i,l,ii,j,k,iprint
      logical scat
C#ifdefC AIX
C      character*72 s
C#endif
      io3smm = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        io3smm = .true.
        backspace ifi
        rh = 0
        vrmax(2) = 0
        read(ifi,'(21x,f8.6)') vrmax(1)
        read(ifi,*)
        do  10  i = 1, nsp
        do  16  l = 0, lmax
          read(ifi,100,err=12) ii, (xx(j), j=1,4), k
          if (ii /= l) call rx('IO3SMM: bad l quantum number')
          idmod(l) = k
          pl(l,i) = xx(1)
          pl(l,nsp) = xx(1)
          do  15  ii = 1, 3
            ql(ii,l,nsp) = xx(1+ii)/nsp
            ql(ii,l,i) = xx(1+ii)
   15     continue
   16   continue
C --- Patch for bug in AIX err= ---
C#ifdefC AIX
C        if (i == 1 .and. nsp == 2) then
C          if (.not. rdstrn(ifi,s,72,.false.)) goto 18
C          if (s(1:5) /= '   0 ') goto 18
C          backspace ifi
C        endif
C#endif
   10   continue
      else
        pi = 4*datan(1.d0)
        write(-ifi,'(''MOMNTS:  '',a4,''  rho,vrmax='',3f10.6)')
     .    albl, rh, vrmax
        write(-ifi,891)
  891   format(
     .  '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          dl = dtan(pi*(.5d0 - pl(l,i)))
          if (dabs(dl) > 9999) dl = 0
          write(-ifi,100) l,pl(l,i),(ql(ii,l,i),ii=1,3),idmod(l),dl
   11   continue
      endif
      return

  100 format(i2,4f11.7,i3,f13.7,f11.7)

C --- handle read exception ---
C#ifdefC AIX
C   18 l = 0
C      i = 2
C#endif
   12 continue
      if (l == 0 .and. i == 2) then
        if (iprint() >= 20)
     .  print *, 'IO3SMM, ATOM=',albl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  40  l = 0, lmax
          do  40  ii = 1, 3
          ql(ii,l,1) = ql(ii,l,2)
   40   continue
      else
        io3smm = .false.
      endif
      end

      logical function io4mm(albl,pl,ql,idmod,nl,lmax,nsp,rh,vrmax,ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs
Ci   albl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   pl:   related to log derivatives (see atom)
Ci   ql:   moments q (see atom) (ifi > 0)
Co   rh,vrmax density and Vxc at rmax (ifi>0)
Co Outputs
Co   pl,ql:   moments q (see atom) (ifi < 0)
Co   aiomom:true unless read error or category not found
Co   rh,vrmax density and Vxc at rmax (ifi<0)
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   Bug in reading vrmax for nsp=2
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nl,lmax,nsp,ifi
      character*4 albl
      integer idmod(0:nl-1)
      double precision pl(0:nl-1,nsp),ql(3,0:nl-1,nsp),rh,vrmax(2)
C Local parameters
      double precision pi,dl,xx(4)
      integer i,l,ii,j,k,iprint
      logical scat
C#ifdefC AIX
C      character*72 s
C#endif

      io4mm = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        io4mm = .true.
        backspace ifi
        read(ifi,'(25x,3f10.6)') rh, vrmax
        read(ifi,*)
        do  10  i = 1, nsp
        do  16  l = 0, lmax
          read(ifi,100,err=12) ii, (xx(j), j=1,4), k
          if (ii /= l) call rx('IO4MM: bad l quantum number')
          idmod(l) = k
          pl(l,i) = xx(1)
          pl(l,nsp) = xx(1)
          do  15  ii = 1, 3
            ql(ii,l,nsp) = xx(1+ii)/nsp
            ql(ii,l,i) = xx(1+ii)
   15     continue
   16   continue
C --- Patch for bug in AIX err= ---
C#ifdefC AIX
C        if (i == 1 .and. nsp == 2) then
C          if (.not. rdstrn(ifi,s,72,.false.)) goto 18
C          if (s(1:5) /= '   0 ') goto 18
C          backspace ifi
C        endif
C#endif
   10   continue
      else
        pi = 4*datan(1.d0)
        write(-ifi,'(''MOMNTS:  '',a4,''  rho,vrmax='',3f10.6)')
     .    albl, rh, vrmax
        write(-ifi,891)
  891   format(
     .  '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          dl = dtan(pi*(.5d0 - pl(l,i)))
          if (dabs(dl) > 9999) dl = 0
          write(-ifi,100) l,pl(l,i),(ql(ii,l,i),ii=1,3),idmod(l),dl
   11   continue
      endif
      return

  100 format(i4,4f13.7,i4,f13.7)

C --- handle read exception ---
C#ifdefC AIX
C   18 l = 0
C      i = 2
C#endif
   12 continue
      if (l == 0 .and. i == 2) then
        if (iprint() >= 20)
     .  print *, 'IO4MM, ATOM=',albl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  40  l = 0, lmax
          do  40  ii = 1, 3
          ql(ii,l,1) = ql(ii,l,2)
   40   continue
      else
        io4mm = .false.
      endif
      end

      logical function io4pp(albl,pp,ves,nl,lmax,nsp,ifi)
C- File I/O for potential parameters.
C ----------------------------------------------------------------
Ci Inputs
Ci   albl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   ves:  electrostatic potential at rmax
Ci   pp:   potential parameters (ifi > 0) read from file ifi
Co Outputs
Co   pp:   potential parameters (ifi < 0) written to file -ifi
Co   aiopar:true unless read error or category not found
Cr Remarks
Cr   11 Apr 94 Added convention: reading 1/p=0 => set p to 0
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*4 albl
      integer nl,lmax,nsp,ifi
      double precision pp(6,0:nl-1,nsp),ves
C Local parameters
      double precision xx,x2(6)
      equivalence (xx,x2)
      integer i,l,ll,k,iprint,i1mach
      logical scat,sw
C#ifdefC AIX
C      logical rdstrn
C      character*72 s
C#endif

      sw = .false.
      io4pp = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        io4pp = .true.
        backspace ifi
        read(ifi,'(17x,f12.8)') ves
        read(ifi,*)
        do  10  i = 1, nsp
        do  16  l = 0, lmax
          read(ifi,20,err=12) ll, (x2(k), k=1,6)
          if (ll /= l) call rx('IO4PP: bad l quantum number')
C Map delta into sqrt(delta), preserving sign and 1/sqrt(p) into p
          x2(3) = dsign(1.d0,x2(3))*dsqrt(dabs(x2(3)))
          if (x2(4) == 0) then
            if (iprint() > 40)
     .        print *, 'io4pp: encountered 1/p=0 ... set p to 0'
          else
            x2(4) = 1/x2(4)**2
          endif
          do  15  k = 1, 6
          pp(k,l,i) = x2(k)
   15     pp(k,l,nsp) = x2(k)
   16     continue
C --- Patch for bug in AIX err= ---
C#ifdefC AIX
C        if (i == 1 .and. nsp == 2) then
C          if (.not. rdstrn(ifi,s,72,.false.)) goto 18
C          if (s(1:3) /= ' 0 ') goto 18
C          backspace ifi
C        endif
C#endif
   10   continue
      else
        write(-ifi,21) albl, ves
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          xx = pp(3,l,i)
          if (pp(4,l,i) == 0) then
            write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .        0d0, (pp(k,l,i), k=5,6)
            sw = .true.
          else
            write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .        1/dsqrt(pp(4,l,i)), (pp(k,l,i), k=5,6)
        endif
   11   continue
        if (iprint() >= 30  .and. sw)
     .    write(i1mach(2),*) 'io4pp: encountered p=0 ... wrote 1/p=0'
        io4pp = .true.
      endif
      return

C --- handle read exception ---
C#ifdefC AIX
C   18 l = 0
C      i = 2
C#endif
   12 continue
      if (l == 0 .and. i == 2) then
        if (iprint() >= 20)
     .  print *, 'IO4PP, ATOM=',albl,
     .  ':  spin 2 input missing; taken from spin 1 data'
      else
        io4pp = .false.
      endif

   20 format(i2,3f12.8,f12.7,3f12.8)
   21 format('PPAR:  ',a4,'  ves=',f12.8/ ' l',5x,'e_nu',10x,'C',8x,
     .       '+/-del',5x,'1/sqrt(p)',6x,'gam',9x,'alp')

      end

      logical function io4pt(nr,nsp,a,rmax,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci   nr,nsp,a,rmax
Ci   v, if file write
Co Outputs
Co   v, if file read
Cr Remarks
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, aiopot returns true only if the category is found,
Cr    the file's value of a and nr match input and rmax is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,v(nr,nsp)
C Local parameters
      integer i,isp,nr2,nsp2
      double precision a2,rmax2
      logical scat

      io4pt = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'POT:',':',.true.)) return
        read(ifi,102,err=15) nr2,nsp2,a2,rmax2
        if (nr == 0) nr=nr2
        if (nsp2 == 0) nsp=nsp2
        if (a2 == 0) a=a2
        if (rmax == 0) rmax=rmax2
        if (a2 /= a .or. nr /= nr2 .or.
     .      dabs(rmax2-rmax) > .03) goto 15
        do  10  isp = 1, min0(nsp2,nsp)
c          read(ifi,101,err=15) (v(i,isp),i = 1,nr)
          read(ifi,101) (v(i,isp),i = 1,nr)
          do  10  i = 1, nr
          v(i,nsp) = v(i,isp)
   10   continue
        io4pt = .true.
   15   continue
      else
        write(-ifi,'(''POT:'')')
        write(-ifi,102) nr,nsp,a,rmax
        do  11  isp = 1, nsp
          write(-ifi,101) (v(i,isp),i = 1,nr)
   11   continue
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f12.5)
      return
      end

C ---- File I/O atomic data, Stuttgart v4 ---
      logical function io4sg(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .  qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
      implicit none
C Passed parameters
      integer lmx,nsp,nr,ifi
      double precision z,rmax,a,
     . qc,dq,vrmax(2),sumec,sumev,thrpv,ekin,utot,rhoeps,etot
      character*8 alabel
      logical lrel
C Local parameters
      logical scat

      io4sg = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'GEN:',':',.true.)) return
        read(ifi,200) alabel,lmx,nsp,rmax,lrel,nr,a
        read(ifi,201) z,qc,dq,vrmax,rhoeps
        read(ifi,202) sumec,utot,ekin,sumev,etot,thrpv
        io4sg = .true.
      else
        write(-ifi,'(''GEN:'')')
        write(-ifi,100) alabel,lmx,nsp,rmax,lrel,nr,a
        write(-ifi,101) int(z),int(qc),dq,vrmax,rhoeps
        write(-ifi,102) sumec,utot,ekin,sumev,etot,thrpv
      endif
  100 format(3X,A4,'  LMX=',I1,'  NSPIN=',I1,'  RMAX=',F9.6,'  REL=',L1,
     .       '  NR=',I3,'  A=',F4.3)
  101 format(3X,'Z=',I2,  '  QC=',I2,'  QTOT=',F9.6,'  VRMAX=',2F9.6,
     .          '  RHOEPS=',F12.6)
  102 format(
     .  '   SUMEC=',F15.6,'  UTOT=',F15.6,'  EKIN= ',F15.6/
     .  '   SUMEV=',F15.6,'  ETOT=',F15.6,'  THRPV=',F15.6)

  200 format(3X,A4,6X,I1,8X,I1,7X,F9.6,6X,L1,5X,I3,4X,F4.3)
  201 format(3X,2X,F2.0,5X,F2.0,7X,F9.6,8X,2F9.6,9X,F12.6)
  202 format(2(9X,F15.6,7X,F15.6,8X,F15.6/))

      end
      logical function io4smm(alabl,pl,ql,evrl,idmod,nl,lmax,nsp,vrmax,
     .        ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   pl:   related to log derivatives (see atom)
Ci   ql:   moments q (see atom) (ifi > 0)
Ci   evrl: eshift (if idmod /= 0, ifi > 0)
Co Outputs
Co   pl,ql:   moments q (see atom) (ifi < 0)
Ci   evrl: eshift (if idmod /= 0, ifi < 0)
Co   io4smm:true unless read error or category not found
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nl,lmax,nsp,ifi
      character*8 alabl
      integer idmod(0:nl-1)
      double precision pl(0:nl-1,nsp),ql(3,0:nl-1,nsp),evrl(0:nl-1,nsp),
     .                 vrmax
C Local parameters
      double precision pi,dl,xx(6)
      integer i,l,ii,j,k,iprint
      logical scat
      external scat,iprint

      io4smm = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        io4smm = .true.
        backspace ifi
        read(ifi,'(21x,f10.6)') vrmax
        read(ifi,*)
        do  10  i = 1, nsp
          do  10  l = 0, lmax
          read(ifi,100,err=12) ii, (xx(j), j=1,4), k, (xx(j), j=5,6)
          if (ii /= l) stop 'IO4SMM: bad l quantum number'
          idmod(l) = k
          evrl(l,i) = xx(6)
          evrl(l,nsp) = xx(6)
          pl(l,i) = xx(1)
          pl(l,nsp) = xx(1)
          do  15  ii = 1, 3
            ql(ii,l,nsp) = xx(1+ii)/nsp
            ql(ii,l,i) = xx(1+ii)
   15     continue
   10   continue
      else
        pi = 4*datan(1.d0)
        write(-ifi,'(''MOMNTS:  '',A4,''  VRMAX='',F10.6)') alabl, vrmax
        write(-ifi,891)
  891   format(
     .  ' L',6X,'PL',9X,'Q0',9X,'Q1',9X,'Q2',4X,' ID ',5X,'DL',9X,'ESH')
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          dl = dtan(pi*(.5d0 - pl(l,i)))
          write(-ifi,100) l,pl(l,i),(ql(ii,l,i),ii=1,3),idmod(l),dl,evrl(l,i)
   11   continue
      endif
      return

  100 format(i2,4f11.7,i3,f12.7,f11.7)

C --- handle read exception ---
   12 continue
      if (l == 0 .and. i == 2) then
        if (iprint() >= 20)
     .  print *, 'IO4SMM, ATOM=',alabl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  40  l = 0, lmax
          do  40  ii = 1, 3
          ql(ii,l,1) = ql(ii,l,2)
   40   continue
      else
        io4smm = .false.
      endif
      end
      logical function io4spp(alabl,pp,ves,nl,lmax,nsp,ifi)
C- File I/O for potential parameters.
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   ves:  electrostatic potential at rmax
Ci   pp:   potential parameters (ifi > 0) read from file ifi
Co Outputs
Co   pp:   potential parameters (ifi < 0) written to file -ifi
Co   io4spp:true unless read error or category not found
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*8 alabl
      integer nl,lmax,nsp,ifi
      double precision pp(6,0:nl-1,nsp),ves
C Local parameters
      double precision xx,x2(6)
      equivalence (xx,x2)
      integer i,l,ll,k,iprint
      logical scat
      external scat,iprint

      io4spp = .false.

      if (ifi > 0) then
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        io4spp = .true.
        backspace ifi
        read(ifi,'(17x,f10.6)') ves
        read(ifi,*)
        do  10  i = 1, nsp
          do  10  l = 0, lmax
          read(ifi,20,err=12) ll, (x2(k), k=1,6)
          if (ll /= l) stop 'IO4SPP: bad l quantum number'
C Map delta into sqrt(delta), preserving sign and 1/sqrt(p) into p
C#ifndef NOMAP_PP
          x2(3) = dsign(1.d0,x2(3))*dsqrt(dabs(x2(3)))
          x2(4) = 1/x2(4)**2
C#endif
          do  15  k = 1, 6
          pp(k,l,i) = x2(k)
   15     pp(k,l,nsp) = x2(k)
   10   continue
      else
        write(-ifi,21) alabl, ves
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          xx = pp(3,l,i)
c aek avoid 1/sqrt(p) to be infinity
          if ( pp(4,l,i) < 1.d-99 ) then
          write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .                   0.0d0, (pp(k,l,i), k=5,6)
          else
          write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .                   1/dsqrt(pp(4,l,i)), (pp(k,l,i), k=5,6)
          endif
   11   continue
      endif
      return

C --- handle read exception ---
   12 continue
      if (l == 0 .and. i == 2) then
        if (iprint() >= 20)
     .  print *, 'IO4SPP, ATOM=',alabl,
     .  ':  spin 2 input missing; taken from spin 1 data'
      else
        io4spp = .false.
      endif

   20 format(i2,6f12.7)
   21 format('PPAR:  ',A4,'  VES=',F10.6/ ' L',8X,'ENU',9X,'C',8X,
     .       '+/-DEL',4X,'1/SQRT(P)',6X,'GAM',9X,'ALP')

      end
      logical function io4spt(nr,nsp,a,rmax,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci   nr,nsp,a,rmax
Ci   v, if file write
Co Outputs
Co   v, if file read
Cr Remarks
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, io4spt returns true only if the category is found,
Cr    the file's value of a and nr match input and rmax is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,v(nr,nsp)
C Local parameters
      integer i,isp,nr2,nsp2
      double precision a2,rmax2
      logical scat

      io4spt = .false.
      if (ifi > 0) then
        if (.not. scat(ifi,'POT:',':',.true.)) return
        read(ifi,102,err=15) nr2,nsp2,a2,rmax2
        if ( dabs(a2-a) > 1.d-4 .or. nr /= nr2 .or.
     .      dabs(rmax2-rmax) > .03) goto 15
        do  10  isp = 1, min0(nsp2,nsp)
c          read(ifi,101,err=15) (v(i,isp),i = 1,nr)
          read(ifi,101) (v(i,isp),i = 1,nr)
          do  10  i = 1, nr
          v(i,nsp) = v(i,isp)
   10   continue
        io4spt = .true.
   15   continue
      else
        write(-ifi,'(''POT:'')')
        write(-ifi,102) nr,nsp,a,rmax
        do  11  isp = 1, nsp
          write(-ifi,101) (v(i,isp),i = 1,nr)
   11   continue
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f12.5)
      return
      end

C ---- File I/O atomic data, Stuttgart v47+u ---
C ... aiogen is same as 5g
      logical function io47um(clabl,lmx,nl,nsp,pnl,idmod,qnl,vrmax,ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------------
Ci Inputs:
Ci   clabl :name of the different inequivalent atom
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Co Outputs ( Inputs if file write: ifi < 0)
Cio  pnl   :principal quantum number
Cio  qnl   :moments
Cio  vrmax :potential at rmax=wsr
Cr Remarks:
Cr   aiomom:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   zeroth moment in atomic file is taken by convention to be
Cr   electronic charge q0 = m0 + p^gamma m2,  rather than m0.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer nl,lmx,nsp,ifi
C-ise
      integer idmod(0:nl-1)
C-ise
      double precision pnl(0:nl-1,*),qnl(0:2,0:nl-1,*),vrmax
      character*4 clabl
C Local variables:
      integer iost,isp,l,l2,j,k,idm
      double precision pi,dl,xx(6)
      logical scat
      character chr
      character*72 messg
      parameter(pi=3.14159265358979324d0)
C External calls:
      external dcopy,dscal,errmsg,scat
C Intrinsic functions:
      intrinsic dtan

      io47um = .false.
      if (ifi > 0) then
C ----  Read moments and p's from atomic file
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        backspace ifi
        read(ifi,200,err=15,end=15,iostat=iost)vrmax
        read(ifi,*  ,err=15,end=15,iostat=iost)
        do isp = 1, 2
          do l = 0, lmx
           read(ifi,201,err=15,end=15,iostat=iost)chr
           if (chr >= '0'.and.chr <= '9') then
             backspace ifi
             read(ifi,301,err=15,end=15,iostat=iost)
     .            l2,(xx(j),j=1,4),idm
             if (isp <= nsp) then
               idmod(l) = idm
               pnl(l,isp) = xx(1)
               call dcopy(3,xx(2),1,qnl(0,l,isp),1)
             elseif (l == l2) then
               pnl(l,1) = (pnl(l,1)+xx(1))*0.5d0
               do k=0,2
                 qnl(k,l,1)= qnl(k,l,1)+xx(k+2)
               enddo
             endif
           else
             if (isp <= nsp) then
               write(messg,400)clabl,isp,l
               call errmsg(messg,1)
               goto 15
             else
               io47um = .true.
               return
             endif
           endif
          enddo
        enddo
        io47um = .true.
      else
C ----  write atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,vrmax
        do isp = 1, nsp
          do l = 0, lmx
            dl = dtan(pi*(0.5d0 - pnl(l,isp)))
            write(-ifi,301,err=15,iostat=iost)
C-ise     .       l,pnl(l,isp),(qnl(k,l,isp),k=0,2),0,dl
     .       l,pnl(l,isp),(qnl(k,l,isp),k=0,2),idmod(l),dl
          enddo
        enddo
      endif
      return

   15 continue
      if (iost /= 0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi > 0 .and. l == 0 .and. isp == 2) then
        write(messg,401) clabl
        call errmsg(messg,1)
        call dcopy(  nl,pnl (0,1)  ,1,pnl (0,2),1)
        call dcopy(3*nl,qnl (0,0,1),1,qnl (0,0,2),1)
        call dscal(6*nl,0.5d0,qnl,1)
        io47um = .true.
      endif

200   format(21x,f10.6)
201   format(1x,a1)
300   format('MOMNTS:  ',a4,'  VRMAX=',f10.6,
     .  /' L',6x,'PL',9x,'Q0',9x,'Q1',9x,'Q2   IDMOD',5x,'DL')
301   format(i2,4f11.7,i3,f12.7,f11.7)
400   format(' IO47UM: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing moments.$')
401   format(' IO47UM: ATOM ',a4,':spin 2 input missing;',
     .       ' moments split.$')
402   format(' IO47UM: I/O error, iostat=',i3,'$')
      end

      logical function io47up(clabl,lmx,nl,nsp,pp,ves,ifi)
C- File I/O for moments.  ifi > 0 for read, < 0 for write
C ----------------------------------------------------------------------
Ci Inputs: (pp: if file write : ifi < 0)
Ci   ifi   :logical unit: positive for read, negative for write
Ci   clabl :name of the different inequivalent atom
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Cio Inputs/Outputs:
Cio  pp    :potential parameters
Cio  ves   :electrostatic potential at wsr
Cr Remarks:
Cr   aiopar:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer npp
      parameter(npp=6)
      integer nl,lmx,nsp,ifi
      double precision pp(npp,0:nl-1,nsp),ves
      character*4 clabl
C Local variables:
      integer isp,l,l2,k,iost
      double precision x2(6),tiny,vmax
      logical scat
      parameter(tiny=1.d-6,vmax=.999999d9)
      character  chr
      character*4 cdum
      character*80 messg
C External calls:
      external errmsg,scat
C Intrinsic functions:
      intrinsic dabs,dsign,dsqrt

      io47up = .false.

      if (ifi > 0) then
C ----- Read potential parameters from atomic file
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        backspace (ifi)
        read(ifi,500,err=15,end=15,iostat=iost)cdum,ves
        do isp = 1, nsp
          do l = 0, lmx
            read(ifi,200,err=15,end=15,iostat=iost)chr
            if ( chr >= '0'.and. chr <= '9') then
              backspace ifi
              read(ifi,*  ,err=15,end=15,iostat=iost)l2,(x2(k),k=1,6)
            else
              write(messg,401)clabl,isp,l
              call errmsg(messg,1)
              goto 15
            endif
C --------- Map delta into sqrt(delta) preserving sign, 1/sqrt(p) into p
            x2(3) = dsqrt(dabs(x2(3)))*dsign(1.d0,x2(3))
C-ise            x2(4) = dmax1(x2(4),tiny)
C-ise            x2(4) = 1/x2(4)/x2(4)
            call dcopy(6,x2,1,pp(1,l,isp),1)
          enddo
        enddo
        io47up = .true.
      else
C ----- Write potential parameters to atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,ves
        do isp = 1, nsp
          do l = 0, lmx
            do k=1,6
              x2(k) = pp(k,l,isp)
            enddo
C-ise            if (dabs(x2(4)) <= 1.d0/vmax/vmax)  x2(4) = 1.d0/vmax/vmax
            x2(3) = x2(3)*x2(3)*dsign(1.d0,x2(3))
C-ise            x2(4) = 1/dsqrt(x2(4))
            write(-ifi,301,err=15,iostat=iost)l
            do k=1,6
              if(x2(k) > -9.9999999d0.and.x2(k) < 99.9999999d0) then
                write(-ifi,302,err=15,iostat=iost)x2(k)
              else
                write(-ifi,303,err=15,iostat=iost)x2(k)
              endif
            enddo
            write(-ifi,304,err=15,iostat=iost)
          enddo
        enddo
      endif
      return

15    continue
      if (iost /= 0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi > 0 .and. l == 0 .and. isp == 2) then
        write(messg,400) clabl
        call errmsg(messg,1)
        call dcopy(6*(lmx+1),pp(1,0,1),1,pp(1,0,2),1)
        io47up = .true.
      endif

200   format(1x,a1)
300   format('PPAR:  ',a4,'  VES=',f10.6/ ' L',5x,'ENU',9x,'C',8x,
     .       '+/-DEL',7x,'P',8x,'GAMMA',6x,'ALPHA')
301   format(1x,i1,$)
302   format(f11.7,$)
303   format(e11.5e1,$)
304   format()
400   format(' IO47UP: ATOM ',a4,':spin 2 input missing; ',
     .          'taken from spin 1.$')
401   format(' IO47UP: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing potential parameters$')
402   format(' IO47UP: I/O error, iostat=',i3,'$')
500   format(7x,a4,6x,f10.6/,1x)
      end

C ---- File I/O atomic data, Stuttgart v56 ---
      logical function io5sg(a,etot,lmx,lrel,nr,nsp,qc,qtot,
     .                       wsr,z,ifi)
C- File I/O for atomic general data.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs:
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci Inputs/Outputs:
Cio  etot  :total energy of the atom
Cio  lmx   :maximum l-quantum number inside the sphere
Cio  lrel  :true, for relativistic calculation, false otherwise
Cio  nr    :number of mesh points
Cio  nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Cio  qc    :core electronic charge
Cio  qtot  :difference between spin up and spin down charge
Cio  wsr   :Wigner-Seitz radius
Cio  z     :nuclear charge
C ----------------------------------------------------------------
      implicit none
C Passed variables:
      integer lmx,nsp,nr,ifi
      double precision a,etot(13,0:2),qc,qtot,wsr,z
      logical lrel
C Local variables:
      integer i,iqc,iost,isp,iz,iz2,lmx2,lgunit,nr2,nsp2
      double precision a2,dsum,wsr2,tiny
      logical scat
      character*500 messg
      parameter(tiny=1.d-5)
C External calls:
      external  dsum,errmsg,scat,lgunit
C Intrinsic functions:
      intrinsic idnint

      io5sg = .false.

      iz=idnint(z)

      if (ifi > 0) then
        if (.not. scat(ifi,'GENERAL:',':',.true.)) return
        backspace ifi
        read(ifi,300,err=15,end=15,iostat=iost)lmx2,nsp2,wsr2,
     .                                         lrel,nr2,a2
        read(ifi,301,err=15,end=15,iostat=iost)iz2,iqc,qtot
        qc = iqc
        if (nsp == 2) then
          do  10  isp = 1, nsp2
            read(ifi,302,err=15,end=15,iostat=iost)i
            read(ifi,303,err=15,end=15,iostat=iost)(etot(i,isp),i=1,13)
   10     continue
        endif
        read(ifi,303,err=15,end=15,iostat=iost)(etot(i,0),i=1,13)
C        if (dabs(a-a2) > tiny .or.
C     .      lmx /= lmx2 .or. nr /= nr2 .or. nsp /= nsp2 .or.
C     .      dabs(wsr-wsr2) > tiny .or. iz /= iz2) then
C          write(messg,400)clabl,a,a2,lmx,lmx2,nr,nr2,nsp,
C     .                    nsp2,wsr,wsr2,iz,iz2
C          call errmsg(messg,-1)
C   ... Require match for these only wsr,z only
        if (lmx /= lmx2 .or. nsp /= nsp2 .or. dabs(wsr-wsr2) > tiny
     .      .or. iz /= iz2) then
          print *, 'io5sg: file mismatch'
          io5sg = .false.
          return
        else
          io5sg = .true.
        endif
        a = a2
        nr = nr2
      else
       if (-ifi /= lgunit(1)) then
         write(-ifi,300,err=15,iostat=iost)lmx,nsp,wsr,lrel,nr,a
         write(-ifi,301,err=15,iostat=iost)iz,idnint(qc),qtot
       endif
       if (nsp == 2) then
         do  20  isp = 1, nsp
           write(-ifi,302,iostat=iost)isp
           write(-ifi,303,iostat=iost)(etot(i,isp),i=1,13)
   20    continue
       endif
       write(-ifi,303,iostat=iost)(etot(i,0),i=1,13)
       io5sg = .true.
      endif
      return

C --- handle read exception ---
15    continue
      write(messg,401)iost
      call errmsg(messg,1)

300   format('GENERAL: LMX=',i1,'  NSPIN=',i1,'   WSR=',f9.6,'  REL=',
     .       l1,'  NR=',i3,'  A=',f4.3)
301   format(3x,'Z=',i2,  '  QC=',i2,'  QTOT=',f9.6)
302   format(/' SPIN',i2,':')
303   format(/' ETOT  =',f15.7,2x,'EKIN = ',f15.7,2x,'UTOT =',f15.7,
     .       /' EMAD  =',f15.7,2x,'SUMEC= ',f15.7,2x,'SUMEV=',f15.7,
     .       /' RHOV  =',f15.7,2x,'UNUCL= ',f15.7,2x,'UELEC=',f15.7,
     .       /' RHOEPS=',f15.7,2x,'REPSNL=',f15.7,2x,'RHOMU=',f15.7,
     .       /' RMUNL =',f15.7,2x)

400   format(' AIOGEN: difference between data in atomic file ',a4,
     .       '|        and actual values:',
     .       '|        A    ',': actual value:',f8.4 ,'   read:',f8.4,
     .       '|        LMX  ',': actual value:',i8   ,'   read:',i8,
     .       '|        NR   ',': actual value:',i8   ,'   read:',i8,
     .       '|        NSP  ',': actual value:',i8   ,'   read:',i8,
     .       '|        WSR  ',': actual value:',f8.4 ,'   read:',f8.4,
     .       '|        Z    ',': actual value:',i8   ,'   read:',i8,'$')
401   format(' AIOGEN: I/O error, iostat=',i3,'$')
      end
      logical function io5smm(clabl,idmenu,lmx,nl,nsp,pnl,qnl,vrmax,ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------------
Ci Inputs:
Ci   clabl :name of the different inequivalent atom
Ci   idmenu:if lmul=t and idmenu=1 => orbital collected to a single enu
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Co Outputs ( Inputs if file write: ifi < 0)
Cio  pnl   :principal quantum number
Cio  qnl   :moments
Cio  vrmax :potential at rmax=wsr
Cr Remarks:
Cr   aiomom:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   zeroth moment in atomic file is taken by convention to be
Cr   electronic charge q0 = m0 + p^gamma m2,  rather than m0.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer idmenu(0:*),nl,lmx,nsp,ifi
      double precision pnl(0:nl-1,*),qnl(0:2,0:nl-1,*),vrmax
      character*4 clabl
C Local variables:
      integer iost,isp,l,l2,j,k,ii
      double precision pi,dl,xx(6)
      logical scat
      character chr
      character*72 messg
      parameter(pi=3.14159265358979324d0)
C External calls:
      external dcopy,dscal,errmsg,scat
C Intrinsic functions:
      intrinsic dtan

      io5smm = .false.
      if (ifi > 0) then
C ----  Read moments and p's from atomic file
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        backspace ifi
        read(ifi,200,err=15,end=15,iostat=iost)vrmax
        read(ifi,*  ,err=15,end=15,iostat=iost)
        do  10  isp = 1, 2
          do  12  l = 0, lmx
           read(ifi,201,err=15,end=15,iostat=iost)chr
           if (chr >= '0'.and.chr <= '9') then
             backspace ifi
             read(ifi,301,err=15,end=15,iostat=iost)
     .            l2,(xx(j),j=1,4),ii
             if (isp <= nsp) then
               pnl(l,isp) = xx(1)
               call dcopy(3,xx(2),1,qnl(0,l,isp),1)
               idmenu(l) = ii
             elseif (l == l2) then
               pnl(l,1) = (pnl(l,1)+xx(1))*0.5d0
               do  14  k=0,2
                 qnl(k,l,1)= qnl(k,l,1)+xx(k+2)
   14          continue
               idmenu(l) = ii
             endif
           else
             if (isp <= nsp) then
               write(messg,400)clabl,isp,l
               call errmsg(messg,1)
               goto 15
             else
               io5smm = .true.
               return
             endif
           endif
   12    continue
   10   continue
        io5smm = .true.
      else
C ----  write atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,vrmax
        do  20  isp = 1, nsp
          do  22  l = 0, lmx
            dl = dtan(pi*(0.5d0 - pnl(l,isp)))
            write(-ifi,301,err=15,iostat=iost)
     .       l,pnl(l,isp),(qnl(k,l,isp),k=0,2),idmenu(l),dl
   22     continue
   20   continue
      endif
      return

   15 continue
      if (iost /= 0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi > 0 .and. l == 0 .and. isp == 2) then
        write(messg,401) clabl
        call errmsg(messg,1)
        call dcopy(  nl,pnl (0,1)  ,1,pnl (0,2),1)
        call dcopy(3*nl,qnl (0,0,1),1,qnl (0,0,2),1)
        call dscal(6*nl,0.5d0,qnl,1)
        io5smm = .true.
      endif

200   format(21x,f10.6)
201   format(1x,a1)
300   format('MOMNTS:  ',a4,'  VRMAX=',f10.6,
     .  /' L',6x,'PL',9x,'Q0',9x,'Q1',9x,'Q2    IDMENU',5x,'DL')
301   format(i2,4f11.7,2x,i3,1x,f12.7,f11.7)
400   format(' AIOMOM: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing moments.$')
401   format(' AIOMOM: ATOM ',a4,':spin 2 input missing;',
     .       ' moments split.$')
402   format(' AIOMOM: I/O error, iostat=',i3,'$')
      end
      logical function io5spp(clabl,lmx,nl,nsp,pp,ves,ifi)
C- File I/O for moments.  ifi > 0 for read, < 0 for write
C ----------------------------------------------------------------------
Ci Inputs: (pp: if file write : ifi < 0)
Ci   ifi   :logical unit: positive for read, negative for write
Ci   clabl :name of the different inequivalent atom
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Ci Inputs/Outputs: (if file read : ifi > 0)
Cio  pp    :potential parameters
Cio  ves   :electrostatic potential at wsr
Cr Remarks:
Cr   aiopar:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer npp
      parameter(npp=11)
      integer nl,lmx,nsp,ifi
      double precision pp(npp,0:nl-1,nsp),ves
      character*4 clabl
C Local variables:
      integer isp,l,l2,k,iost
      double precision x2(7),tiny,vmax
      logical scat
      parameter(tiny=1.d-6,vmax=.999999d9)
      character  chr
      character*4 cdum
      character*80 messg
C External calls:
      external errmsg,scat
C Intrinsic functions:
      intrinsic dabs

      io5spp = .false.

      if (ifi > 0) then
C ----- Read potential parameters from atomic file
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        backspace (ifi)
        read(ifi,300,err=15,end=15,iostat=iost)cdum,ves
        do  10  isp = 1, nsp
          do  12  l = 0, lmx
            read(ifi,200,err=15,end=15,iostat=iost)chr
            if ( chr >= '0'.and. chr <= '9') then
              backspace ifi
              read(ifi,*  ,err=15,end=15,iostat=iost)l2,(x2(k),k=1,7)
            else
              write(messg,401)clabl,isp,l
              call errmsg(messg,1)
              goto 15
            endif
            call dcopy(4,x2,1,pp(1,l,isp),1)
            pp( 9,l,isp)=x2(5)
            pp(10,l,isp)=x2(6)
            pp(11,l,isp)=x2(7)
   12     continue
   10   continue
        io5spp = .true.
      else
C ----- Write potential parameters to atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,ves
        do  20  isp = 1, nsp
          do  22  l = 0, lmx
            do  24  k=1,4
              x2(k) = pp(k,l,isp)
   24       continue
            x2(5) = pp( 9,l,isp)
            x2(6) = pp(10,l,isp)
            x2(7) = pp(11,l,isp)
            write(-ifi,301,err=15,iostat=iost)l
            do  25  k=1,7
              if(x2(k) > -9.9999999d0.and.x2(k) < 99.9999999d0) then
                write(-ifi,302,err=15,iostat=iost)x2(k)
              else
                write(-ifi,303,err=15,iostat=iost)x2(k)
              endif
   25       continue
            write(-ifi,304,err=15,iostat=iost)
   22     continue
   20   continue
      endif
      return

15    continue
      if (iost /= 0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi > 0 .and. l == 0 .and. isp == 2) then
        write(messg,400) clabl
        call errmsg(messg,1)
        call dcopy(6*(lmx+1),pp(1,0,1),1,pp(1,0,2),1)
        io5spp = .true.
      endif

200   format(1x,a1)
300   format('PPAR:  ',a4,'  VES=',f10.6/ ' L',7x,'ENU',8x,'S',7x,
     .       '  O   ',3x,'   P     ',5x,' PHI ',5x,'PHIDOT',
     .       4x,' PGAMMA')
301   format(1x,i1,$)
302   format(f11.7,$)
303   format(e11.5e1,$)
304   format()
400   format(' AIOPAR: ATOM ',a4,':spin 2 input missing; ',
     .          'taken from spin 1.$')
401   format(' AIOPAR: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing potential parameters$')
402   format(' AIOPAR: I/O error, iostat=',i3,'$')
      end
      logical function io5spt(nr,nsp,a,wsr,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   nr    :number of mesh points
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Ci Inputs:/Outputs (i if file write (ifi < 0), o if file read (ifi > 0))
Cio  v     :spherical potential
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio  wsr   :Wigner-Seitz radius
Co   aiopot:see remarks
Cr Remarks:
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, wsr,
Cr    followed by the potential.
Cr    On reading, aiopot returns true only if the category is found,
Cr    the file's value of a and nr match input and wsr is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer ifi,nr,nsp
      double precision a,wsr,v(nr,nsp)
C Local variables:
      integer i,iost,isp,nr2,nsp2
      double precision a2,wsr2,tiny
      logical scat
      character*300 messg
      parameter(tiny=1.d-5)
C External call
      external errmsg,scat
C Intrinsic functions:
      intrinsic dabs,min0

      io5spt = .false.
      if (ifi > 0) then
C ----- read file
        if (.not. scat(ifi,'POT:',':',.true.)) return
        read(ifi,301,err=15,end=15,iostat=iost)nr2,nsp2,a2,wsr2
        if (dabs(a2-a) > tiny .or. nr /= nr2 .or.
     .      dabs(wsr2-wsr) > tiny) then
          write(messg,400)nr,nr2,a,a2,wsr,wsr2
          call errmsg(messg,-1)
        else
          do  10  isp = 1, min0(nsp2,nsp)
            read(ifi,300,err=15,end=15,iostat=iost)(v(i,isp),i = 1,nr)
            do  12  i = 1, nr
   12       v(i,nsp) = v(i,isp)
   10     continue
          io5spt = .true.
        endif
      else
C ----- write file
        write(-ifi,'(''POT:'')',err=15,iostat=iost)
        write(-ifi,301,err=15,iostat=iost)nr,nsp,a,wsr
        do  20  isp = 1, nsp
          write(-ifi,300,err=15,iostat=iost)(v(i,isp),i=1,nr)
   20   continue
        io5spt = .true.
      endif
      return

C --- handle read exception ---
15    continue
      write(messg,401)iost
      call errmsg(messg,1)

300   format(1p,5d16.9)
301   format(2i5,2f12.5)
400   format(' AIOPOT: difference between data in atomic file ',
     .       '|        and actual values:',
     .       '|        NR ',': actual value:',i8  ,'   read:',i8,
     .       '|        A  ',': actual value:',f8.4,'   read:',f8.4,
     .       '|        WSR',': actual value:',f8.4,'   read:',f8.4,'$')
401   format(' AIOPOT: I/O error, iostat=',i3,'$')
      end
      subroutine aioxt2(i1,i2,tau,pos)
C- Shifts pos(1..3,i1..i2) by a const tau
      implicit none
      integer i1,i2
      double precision tau(3),pos(3,*)
      integer i,j
      do  i = i1, i2
        do  j = 1, 3
          pos(j,i) = pos(j,i) + tau(j)
        enddo
      enddo
      end
C      subroutine aioxt3(k,s_site1,n1,s_site2,n2)
CC ... For structures
C      include 'structures.h'
CC      type(str_site)::  s_site1(*),s_site2(*)
C      integer n1,n2,k
C      double precision ssite1(k,n1),ssite2(k,n2)
C
C      end
