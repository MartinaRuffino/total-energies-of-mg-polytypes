      subroutine supcel(mode,s_ctrl,s_lat,s_site,s_spec,s_pot,s_bz,s_move,
     .  slabl,plx,nkd,dlat,nbx)
C- Generates basis vectors in a supercell
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nspec nclass nl nspin nccomp cllst clp clssl
Ci                 group initc ics idcc ipc ipcp ips mxcst nrc ncomp
Ci                 pgfsl pgfvl pgord pgplp spid dclabl clabl pos rmax
Ci                 nclasp
Co     Stored:     cllst clp clssl group initc ics idcc ipc ipcp ips
Co                 mxcst nrc ncomp pgfsl pgfvl pgord pgplp spid clabl
Co                 dclabl pos rmax
Co     Allocated:  *
Cio    Elts passed:lrs dclabl
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 asars2
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat avw nabc ag bgv cg cy dlv gv gvq indxcg
Ci                 ips0 istab jcg kv igv igv2 kv2 pos qlv symgr s_sym
Co     Stored:     alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                 jcg kv igv igv2 kv2 pos qlv symgr s_sym
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos class force vel eula pl relax vshft pnu pz
Ci                 v0 v1 bxc cpawt omg omgn domg dmat gc gcu gcorr gii
Ci                 sfvrtx j0 pdos rho1 rho2 rhoc rho1x rho2x rhocx qhhl
Ci                 qhkl qkkl eqhhl eqhkl eqkkl sighh sighk sigkk tauhh
Ci                 tauhk taukk pihh pihk pikk sohh sohk sokk sighhx
Ci                 sighkx sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx
Ci                 thet clabel
Co     Stored:     class pos pos0 force vel spec clabel pnu pz bxc
Co                 cpawt omg omgn domg dmat gc gcu gcorr gii sfvrtx j0
Co                 pdos rho1 rho2 rhoc rho1x rho2x rhocx qhhl qhkl qkkl
Co                 eqhhl eqhkl eqkkl sighh sighk sigkk tauhh tauhk
Co                 taukk pihh pihk pikk sohh sohk sokk sighhx sighkx
Co                 sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx thet
Co                 v0 v1
Co     Allocated:  v0 v1
Cio    Elts passed:rhoc rho1 rho2
Cio    Passed to:  asars iorsa bcast_strx asars1 iopos
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z name a nr rmt lmxa lmxl kmxt p pz lfoca qc idmod
Ci                 rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
Ci                 chfa rsmfa rhoc
Co     Stored:     name lmxa a nr rmt z lmxl kmxt p pz qc lfoca rfoca
Co                 coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa rhoc
Co     Allocated:  rhoc
Cio    Elts passed:rhoc rmt nr a idmod z
Cio    Passed to:  asars iorsa bcast_strx asars2 iosits
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  ves aamom bxc cp ddpf dddpf ddpfr dlmwt dmatk dpf
Ci                 dpfr gibbs gma gmar grrme mad mxy palp papg pf pfnc
Ci                 pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qnur qpp
Ci                 qt rhat rnew rhos rhrmx socscl sop shfac thetcl vdif
Ci                 vintr vrmax vshft smpot smrho smrout GFr
Co     Stored:     ves aamom bxc cp ddpf dddpf ddpfr dlmwt dmatk dpf
Co                 dpfr gibbs gma gmar grrme mad mxy palp papg pf pfnc
Co                 pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qnur qpp
Co                 qt rhat rnew rhos rhrmx socscl sop shfac thetcl vdif
Co                 vintr vrmax vshft smpot smrho smrout GFr
Co     Allocated:  ves
Cio    Elts passed:pnu qnu rhrmx vrmax
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 asars2
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  ef def w n dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Stored:     ef def w n dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  kt ts tstot tsequ
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  sqssim
Ci Inputs
Ci   mode  :1s digit
Ci         :0 Generate nbx only and return; no other options exercised.
Ci         :100s digit: return mapping of original sites to superlattice sites
Ci   slabl :vector of species labels
Ci   plx   :supercell primitive lattice vectors
Ci   nkd   :number of direct lattice vectors (lattic.f)
Ci   dlat  :direct lattice vectors
Co Outputs
Co   nbx   :number of sites in superlattice
Cl Local variables
Cl   qlx   :supercell reciprocal lattice
Cl   ilstt :A subset of original sites to be generated into supercell.
Cl         :List may be null (nlst=0) in which case all sites are to be used.
Cl   ... The following apply to the SQS maker
Cl   nbs   :number of elements in ibs
Cl   ibs   :list of sites participating in SQS and corresponding species
Cl         :ibs(i) = site participating in SQS.
Cl   ipss  :(SQS) permuted species table
Cl   r3max :Constraint; triplets excluded for which (di+dj)/2 > r3max
Cl   alf2  :exponential decay of 3-center interactions: e2 = scal2*exp(-alf2*r2c)
Cl   scal2 :prefactor for scaling of pair correlation energy
Cl   alf3  :exponential decay of 3-center interactions: e3 = scal3*exp(-alf3*r3max)
Cl   scal3 :prefactor for scaling of triplet correlation energy
Cl   iaxs  :pair table for SQS.  Same form as standard iax table (no on-site terms)
Cl   iax3  :triplets table for SQS.  Refers to pair table iaxs
Cl         :iax3(1,i3c) = kb
Cl         :iax3(2,i3c) = ib
Cl         :iax3(3,i3c) = jb
Cl         :iax3(4,i3c) = ip : points to pair index iaxs(:,ip) connecting ib and kb
Cl         :iax3(5,i3c) = jp : points to pair index iaxs(:,jp) connecting jb and kb
Cl         :iax3(6,i3c) = kp : points to pair index iaxs(:,kp) connecting ib and jb
Cl   imap  :imap(i,1) indicates which site site i in the SL originated from
Cl         :imap(i,2) not used now
Cs Command-line switches
Cs   --rsta       :Makes supercell ASA restart file from original rsta file
Cs   --rsta,amom  :Same as --rsta, but in noncolinear case, if moment is < 0
Cs   --sqs        :
Cs   --seed=#     : Used with SQS
Cs   --wsite...   :Writes a site-file to disk.  ... refers to arguments:
Cs                :... = [~map]~fnam
Cs                :fnam is the name of the site file name
Cs                :Optional ~map appends the correspondence of the original
Cs                :sites to the supercell sites.  Maps are used elsewhere,
Cs                :e.g. for supercell FP charge densities and self-energies
Cs   --wsitex...  :Same as --wsite, but write positions as multiples of plat
Cs                :See --wsite for args "..." following --wsitex.
Cs   --sort~plat  :order by (integer) multiples of original lattice vectors.
Cs                :sites a common original unit cell are kept together.
Cs   --sort:"expr1 [expr2] [expr3]"  :orders site table according to expressions
Cs                :involving Cartesian components of positions, e.g. --sort:'x3 x2'
Cs   --wpos=fn    :Write positions to file fn
Cs   --ring:      :ring:i1,i2 shifts i1..i2-1 one register higher, and site
Cs                :i2 becomes site it
Cs   --swap:      :swap:i1,i2 swaps pairs i1 and i2
Cs   --first      :not documented
Cs   --sites~lst  :Make supercell of subset of sites in original basis
Cs   --keep       :Do not search for sites in supercell, but preserve sites as given.
Cs                :Neither site positions nor nbas will be altered.
Cs   --disp       :Displace a set of atoms in the neighborhood of a given one
Cs                :e.g. --disp:tab2:style=3:21x
Cs   --shorten=   :Suppress all shortening of basis vectors (?)
Cs   --refsite... :compute shifts relative to reference site file, scale them
Cs                :and add scaled shifts to reference
Cs                :tags:
Cs                :~scl=# (required) scaling of shifts.  #=1 > positions are unchanged
Cs                :~fn=ref-site-file (required) name of reference site file
Cs                :~plx (optional) measure shifts in units of respective lattice vectors
Cs                :~swap (optional) tells supcel to exchange the role of 'current' and 'reference'
Cs   --dlat:      :Use these multiples of plat as basis vectors
Cs   --wrsj       :Write pairwise exchange interactions
Cs   --xshft[x]=..:translate site positions by constant vector (Cartesian or multiples of plx)
Cs   --scala=#    :Multiply alat by #, divide lattice and site vectors by #
Cr Remarks
Cr   The first occurence of any site in the supercell corresponds to
Cr   the corresponding atom in the original cell.
Cr   Guaranteed because dlat(:,1) = (0,0,0)
Cb Bugs
Cb   No check is made to see whether ther are enough dlat to encompass
Cb   all lattice vectors
Cu Updates
Cu   31 Oct 18 Fix --sort:expr
Cu   05 Oct 18 New --sort~plat; fix --keep
Cu   11 Jul 18 New --refsite and --keep
Cu   04 Sep 16 SQS maker first SQS atom to start at any site. Still must be contiguous
Cu   24 Aug 16 Finds q-points in original cell that map to q=0 in superlattice
Cu   22 Jun 16 Extra arguments to --wsite.  Note slight change in syntax :
Cu             use   wsite[:...]:fn=flnam  rather than wsite[:...]:filnam
Cu   27 Jul 15 First cut at SQS maker
Cu   21 Jul 13 New --ring and --swap options
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   30 Apr 13 new --wsite~map option
Cu             completed (mostly) migration to f90 structures
Cu   25 Oct 11 Started migration to f90 structures
Cu   04 Feb 11 New --rsta
Cu   01 Aug 08 New --wsitex
Cu   19 Apr 03 Changed distortion to read from displacements file
Cu   12 Apr 03 Some bug fixes when partial list of sites saved
Cu   09 Apr 03 Added pairwise distortion; also writes original cell
Cu             site sdata first
Cu   24 Nov 02 R.S. pairwise hamiltonian may be created for supercell
Cu             from file rsj
Cu   26 Nov 02 supcel can restrict supercell sites to a specified list
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*8 slabl(1)
      integer nkd,mode,nbx
      double precision dlat(3,nkd)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_pot)::   s_pot
      type(str_bz)::    s_bz
      type(str_move)::  s_move
C ... Dynamically allocated local arrays
      integer, allocatable :: iax(:),iaxs(:,:),ntabs(:)
      integer, allocatable :: iax3(:,:),iax3c(:,:),ibs(:),ipss(:),ntab3s(:),sigms(:)
      integer,pointer :: ips(:)
      real(8),pointer :: z(:,:)
      integer,allocatable :: ntab(:),imap(:,:),ilstt(:),lst1(:),iwk(:),iprm(:),ipc(:)
      real(8),allocatable :: pos(:,:),pos2(:,:),dpos(:,:,:),rtab(:,:),rtabs(:,:),hpair(:),
     .  disp(:),aamom(:),amom2(:),pnu(:,:,:),qnu(:,:,:,:),e2(:),e3(:),qsup(:,:),dl(:,:)
C ... Local parameters
      integer,parameter:: NULLI = -99999
      type(str_site),pointer:: s_ssite(:),s_sitel(:)
      type(str_spec),pointer:: s_specl(:)
      character*72 alabl*8,outs,sortex(5),out2,jfilnm,dc*1
      logical lsort,lfirst,lwmap,lswap,lcycle,ltmp,lkeep,lsplat
      integer i,j,k,m,ib,jb,kb,nbas,is,js,lok,lio,ifi,ibx,
     .  nspec,ixv(36),imark(36),it(3),ip,jp,parg,j1,j2,nexpr,nbxsav,a2vec,nclus3,n3c,
     .  nmark,nfilj,nfilm,nw,stdo,nn,nc,sqspid(2),iv(12)
      integer mxcsiz,nbasp,ifirsj,iorsj,nttab,nlst,ilst,nlst1,nspsqs,r3mode,ipass,ib0,nttab_save
      integer nclass,nl,nsp,n0,niax,mode0,nbs
      parameter (n0=10,niax=10)
      integer mxcsz
      double precision sclrsj,pi,alf2,alf3,esqs(0:3)
      double precision plat(3,3),plx(3,3),qlx(3,3),qlat(3,3),platm(3,3)
      double precision dd(3),tau(3),x(3),vol,volx,tol,alat,rcut,xx,scal2,scal3
C     For SQS
      integer nttabs,nttab3     ! number of pairs, triples with involving an SQS alloy element
      double precision r2max,r3max
      procedure(integer) :: isw,lgunit,fopnx,fopna,fopng,fxst,fopn,rdm,iprint,wordsw,iosite,iosits
      procedure(logical) :: cmdopt,latvec
      procedure(real(8)) :: dlength,dglob,dsum,ddot,avwsr
C     For reference cell
      integer,parameter :: mxspec=256
      integer nbasl,nspecl
      double precision alatl,platl(3,3)
      character*8, pointer :: slabll(:)

C ... Initialization
      alat = s_lat%alat
      plat = s_lat%plat
      nbas = s_ctrl%nbas
      nspec = s_ctrl%nspec
      tol = 1d-6
      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      mode0 = mod(mode,10)
      lkeep = cmdopt('--keep',6,0,outs)

C --- Each class by itself ---
      allocate(ipc(nbas))
      forall (i=1:nbas) ipc(i) = i
      call sitepack(s_site,1,nbas,'-class',1,ipc,xx)
      deallocate(ipc)

C --- Get list of lattice vectors to include ---
      nmark = 0
      if (cmdopt('--dlat:',7,0,outs)) then
        m = 7
        k = a2vec(outs,72,m,2,',; ',3,-3,36,ixv,imark)
        nmark = k/3
        if (3*nmark /= k) call rxi('SUPCEL: ''--dlat:'' takes'//
     .    ' multiples of 3 arguments, but read',k)

        do  ib = 1, nmark
C         Assemble vector as specified multiples of prim. latt. vec.
          call dpzero(x,3)
          call daxpy(3,dble(imark(3*ib-2)),plat(1,1),1,x,1)
          call daxpy(3,dble(imark(3*ib-1)),plat(1,2),1,x,1)
          call daxpy(3,dble(imark(3*ib-0)),plat(1,3),1,x,1)
C         Mark which of dlat it belongs to
          do  k = 1, nkd
            if (tol > abs(x(1)-dlat(1,k))+abs(x(2)-dlat(2,k))+abs(x(3)-dlat(3,k))) then
              imark(ib) = k
              goto 14
            endif
          enddo
          call rx('supcel: missing dlat')
   14     continue
        enddo
      endif

C --- Optionally get a list (subset) of sites in original basis ---
      nlst = 0
      allocate(ilstt(nbas)); call iinit(ilstt,nbas)
      if (cmdopt('--sites',7,0,outs)) then
        i = 8
        allocate(ips(nbas),z(nspec,1))
        ips(:) = s_site(1:nbas)%spec
        z(:,1) = s_spec(1:nspec)%z
C       call awrit2('%N Z= %n:1,3;3d',' ',100,6,nspec,z)
        dc = outs(i:i)
        call baslst(0,11,outs(i:),j2,ips,nbas,slabl,z,0,' ',x,nlst,ilstt)
        if (nlst == 0) ilstt(1) = 0
        deallocate(ips,z)
      endif

C --- Supercell reciprocal lattice vectors ---
      call info0(10,1,0,' Basis vectors of supercell:')
C     write(stdo,'(/a)') ' Basis vectors of supercell:'
      call dinv33(plx,1,qlx,volx)
      call dinv33(plat,1,qlat,vol)
      if (dabs(volx) < 1d-6)
     .  call rx('supercell vectors are linearly dependent')
      if (iprint() > 20) then
        write(stdo,351)
  351   format(/15x,'plx',28x,'plx/plat',27x,'qlx')
        call dgemm('T','N',3,3,3,1d0,plx,3,qlat,3,0d0,platm,3)
        do  k = 1, 3
          call info5(2,0,0,'%3;10,5D   %3;10,5D   %3;10,5D',
     .      plx(1,k),platm(k,1:3),qlx(1,k),0,0)
        enddo
      endif

C --- Printout: angles between plx ---
      if (iprint() > 20) then
      do  i = 1, 3
        j = mod(i,3)+1
        tau(i) = ddot(3,plx(1,i),1,plx(1,j),1)/dlength(3,plx(1,i),1)/dlength(3,plx(1,j),1)
      enddo
      call info5(10,1,0, ' Angles between SL vectors : %;11,6D (1,2)  %;11,6D (2,3)  %;11,6D (3,1)',
     .  180/pi*dacos(tau(1)),180/pi*dacos(tau(2)),180/pi*dacos(tau(3)),4,5)
      call info5(10,0,0, ' cosine of angles%12f%;11,6D (1,2)  %;11,6D (2,3)  %;11,6D (3,1)',
     .  tau(1),tau(2),tau(3),4,5)
      endif

      if (.not. latvec(3,tol,qlat,plx)) then
        call info0(10,1,0, ' *** Warning!  supercell lattice '//
     .    'vectors are not multiples of plat!')
      endif

      call ivset(ixv,1,3,2)
      if (cmdopt('-shorten=',9,0,outs) .or.
     .  cmdopt('--shorten=',10,0,outs)) then
        ip = 9
        if (cmdopt('--shorten=',10,0,outs)) ip = 10
        call skipbl(outs,len(outs),ip)
        if (parg(' ',2,outs,ip,len(outs),', ',2,3,it,ixv) < 0)
     .    call rxs2('supcel: failed to parse "',outs(1:ip+5),' ..."')
        print '(/''  supcel: shorten vectors with mode ='',3i2)',ixv
      endif

C --- Count supercell basis vectors (nbx) ---
C     and make table mapping supercell sites to original sites
      write (stdo,*) ' '
      nbx = 0
      ilst = 0
      do  ib = 1, nbas
        if (nlst > 0) then
    6     continue
          if (ilst >= nlst) cycle
          if (ilstt(1+ilst) < ib) then
            ilst = ilst+1
            goto 6
          endif
          if (ilstt(1+ilst) /= ib) cycle
        endif

        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        alabl = s_spec(is)%name
        if (lkeep) then
          nbx = nbx+1
        else
          do  k = 1, nkd
            call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
          enddo
        endif
      enddo
C     Make table imap mapping supercell sites to original sites
      allocate(imap(nbx,2),dl(3,1)); call iinit(imap,nbx*2)
      nbx = 0
      ilst = 0
      call info0(20,0,0,'    ib spec%21ftau%18f+ plat  n(ib)')
      do  ib = 1, nbas
        if (nlst > 0) then
    8     continue
          if (ilst >= nlst) cycle
          if (ilstt(1+ilst) < ib) then
            ilst = ilst+1
            goto 8
          endif
          if (ilstt(1+ilst) /= ib) cycle
        endif
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        alabl = s_spec(is)%name
        i = 0
        do  k = 1, nkd
          if (lkeep) then
            nbx = nbx+1
            lok = 1
            dd = tau
          else
            call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
          endif
          if (lok == 1) then
            xx = (dd(1)-tau(1))**2+(dd(2)-tau(2))**2+(dd(3)-tau(3))**2
            x(1) = (dd(1)-tau(1))
            x(2) = (dd(2)-tau(2))
            x(3) = (dd(3)-tau(3))
C           dl = translation as multiple of original lattice vector
            call dgemm('T','N',3,1,3,1d0,qlat,3,x,3,0d0,dl,3)
            if (dlength(3,dl(:,1)-nint(dl(:,1)),1) > tol) call rx('supcel: noninteger lattice translation')
            imap(nbx,1) = ib
            i = i+1
            if (iprint() >= 20) write (stdo,3) ib,alabl,(dd(m),m=1,3),(nint(dl(m,1)),m=1,3),i
    3       format(i6,2x,a,3F11.6,2x,3i3,i5)
          endif
          if (lkeep) exit
        enddo
      enddo
      deallocate(dl)

C --- Printout ---
      if (iprint() >= 20) then
        write (stdo,4) abs(vol),abs(volx),abs(volx/vol)
    4   format(/' cell vol=',f10.6,'     scell vol=',f12.6,'     ratio=',f7.2)
        write (stdo,7) nbas,nbx,dble(nbx)/nbas
    7   format(' cell nbas=',i9,'     scell nbas=',i11,'     ratio=',f7.2)
        if (nlst == 0) then
          if (abs(dble(nbx)/nbas-abs(volx/vol)) > 1d-5) write (stdo,13)
        else
          if (abs(dble(nbx)/nlst-abs(volx/vol)) > 1d-5) write (stdo,13)
        endif
      endif

      nc = nint(abs(volx/vol))
      if (abs(nc - abs(volx/vol)) < tol) then
        k = maxval(platm)+.9
        allocate(qsup(3,nc))
        call psup6(qlx,plat,k,nc,m,qsup)
        call info2(20,1,0,
     .    ' Found %i q point%-1j%?#(n>1)#s## in original cell that map to q=0:'//
     .    '%N%13fq (2*pi/a)%20fq (multiples of qlat)'
     .    ,m,2)
        do  i = 1, m
          call dgemm('T','N',3,1,3,1d0,plat,3,qsup(1,i),3,0d0,x,3)
C         x = matmul(plat,qsup(1:3,i))
          call info2(2,0,0,'%3;11,6D   %3;11,6D',qsup(1,i),x)
        enddo
      endif

      if (mode0 == 0) then
        deallocate(ilstt)
        return
      endif

C --- Set up expressions to sort ---
      lsplat = .false.
      if (cmdopt('--sort',6,0,outs)) then
        lsplat = outs(8:11) == 'plat'
      endif
      lsort  = cmdopt('--sort:',7,0,outs) .and. .not. lsplat
      lfirst = cmdopt('--first',7,0,outs)
      lswap  = cmdopt('--swap:',7,0,outs)
      lcycle = cmdopt('--ring:',7,0,outs)
      call rxx(isw(lsort)+isw(lsplat)+isw(lfirst)+isw(lcycle)+isw(lswap)>1,
     .  'supcel: switches --sort --sort:plat --first --swap --ring are mutually exclusive')

      nexpr = 0
      if (lsort) then
        i = 7
        if (cmdopt('--sort:',7,0,outs)) i = 8
        do  k = 1, 3
          call word(outs(i:),k,j1,j2)
          if (j2 < j1) exit
          nexpr = nexpr+1
          sortex(nexpr) = outs(i-1+j1:i-1+j2)
        enddo
        outs = ' '
        do  k = 1, nexpr
          outs = trim(outs) // '   ' // sortex(k)
        enddo
        call info2(20,1,0,' Sort positions with %i expressions: '//trim(outs),nexpr,0)
      endif

C --- Get list of positions; sort them ---
C     Output is permutation table iprm
      if (lswap) then
        allocate(iprm(nbx))
        do  i = 1, nbx
          iprm(i) = i
        enddo
        lswap = cmdopt('--swap:',7,0,outs)
        m = 7
        k = a2vec(outs,72,m,2,', ',2,-2,36,ixv,imark)
        if (k <= 0 .or. mod(k,2) /= 0)
     .    call rx('supcel: invalid use of --swap')
        if (lsort .or. lfirst .or. lcycle) call info0(20,1,0,
     .    ' supcel (warning) --sort, --first, --ring ignored')
        do  i = 1, k, 2
          j = imark(i)
          m = imark(i+1)
          if (j>nbx .or. m>nbx)
     .      call rx('supcel: invalid use of --swap')
          iprm(j) = m
          iprm(m) = j
        enddo
      elseif (lcycle) then
        allocate(iprm(nbx))
        do  i = 1, nbx
          iprm(i) = i
        enddo
        if (lsort .or. lfirst) call info0(20,1,0,
     .    ' supcel (warning) --sort, --first ignored')
        m = 7
        lcycle = cmdopt('--ring:',m,0,outs)
        k = a2vec(outs,72,m,2,', ',2,-2,2,ixv,imark)
        if (k /= 2) call rx('supcel: invalid use of --ring')
        do  i = min(imark(1),imark(2))+1, max(imark(1),imark(2))
          iprm(i) = i-1
        enddo
        j = min(imark(1),imark(2))
        m = max(imark(1),imark(2))
        iprm(j) = m
      elseif  (lsort .or. lfirst .or. lsplat) then
        allocate(pos(3,nbx),dl(3,nbx))
        allocate(iwk(nbx)); call iinit(iwk,nbx)
        allocate(iprm(nbx))
        nbxsav = nbx
        nbx = 0
        ilst = 0
        do  ib = 1, nbas
          if (nlst > 0) then
    9       continue
            if (ilst >= nlst) cycle
            if (ilstt(1+ilst) < ib) then
              ilst = ilst+1
              goto 9
            endif
            if (ilstt(1+ilst) /= ib) cycle
          endif

          is = s_site(ib)%spec
          tau = s_site(ib)%pos
C         alabl = s_spec(is)%name
          do  k = 1, nkd
            call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
            if (lok == 1) then
              call dcopy(3,dd,1,pos(1,nbx),1)
              iwk(nbx) = nbx
C             case lfirst: if dd = tau, this was the original vector
              if (lfirst .or. lsplat) then
                xx = (dd(1)-tau(1))**2+(dd(2)-tau(2))**2+(dd(3)-tau(3))**2
                x(1) = (dd(1)-tau(1))
                x(2) = (dd(2)-tau(2))
                x(3) = (dd(3)-tau(3))
C               dl = translation as multiple of original lattice vector
                if (lsplat) then
                  call dgemm('T','N',3,1,3,1d0,qlat,3,x,3,0d0,dl(1,nbx),3)
                elseif (latvec(1,tol,qlx,x)) then
                  iwk(nbx) = ib-nbas
                endif
              endif
            endif
          enddo
        enddo ! complete supercell list.  Hang onto dl for subsequent printout

C     call prmx('unsorted pos',pos,nbx,nbx,3)

        if (lsort) then
          allocate(pos2(nbx,3),z(nbx,2))
          call dmcpy(pos,3,1,pos2,1,nbx,3,nbx)
          do  ibx = 1, nbx
            is = s_site(imap(ibx,1))%spec
            z(ibx,1) = s_spec(is)%z; z(ibx,2) = is
          enddo
          call mapdatx(nexpr,sortex,' ',iprm,nbx,3,z,pos2,[xx])
          deallocate(z); allocate(z(nexpr,nbx))
          call dmcpy(pos2,nbx,1,z,1,nexpr,nbx,nexpr)
          call dvheap(nexpr,nbx,z,iprm,1d-7,101)
          deallocate(pos2,z)
        elseif (lsplat) then
          call dvheap(3,nbx,dl,iprm,1d-8,101)
C         call yprm('iprm',0,iprm,0,nbx,nbx,1)
        else
          call ivheap(1,nbx,iwk,iprm,1)
C     call yprm('ordering',0,iwk,0,nbx,nbx,1)
C     call yprm('iprm',0,iprm,0,nbx,nbx,1)
        endif

        if (iprint() > 40) then
        call info0(20,1,0,'   Superlattice positions in permuted order')
        if (lsplat) then
          call info0(20,0,0,'   ib spec%21fpos%30fpos/plx%16f+ plat')
        else
          call info0(20,0,0,'   ib spec%21fpos')
        endif
        do  i = 1, nbx
          ibx = iprm(i)
          ib = imap(ibx,1)
          is = s_site(ib)%spec
          tau = s_site(ib)%pos
          alabl = s_spec(is)%name
          if (lsplat) then
C           dl = translation as multiple of original lattice vector
            call dgemm('N','N',3,1,3,1d0,plat,3,dl(1,ibx),3,0d0,dd,3)
            x(1) = dd(1)+tau(1)
            x(2) = dd(2)+tau(2)
            x(3) = dd(3)+tau(3)
            call dgemm('T','N',3,1,3,1d0,qlx,3,x,3,0d0,tau,3)
            write (stdo,2) ib,alabl,(x(m),m=1,3),(tau(m),m=1,3),(nint(dl(m,ibx)),m=1,3)
          else
            write (stdo,2) ib,alabl,pos(:,ibx)
          endif
    2     format(i5,2x,a,3F11.6,2x,3F11.6,2x,3i3,i5)
        enddo
        endif

        deallocate(iwk,pos)
      endif

C --- Read R.S. hamiltonian from file rsj ---
C     ifirsj is the file logical unit for real-space hamiltonian
C     ifirsj = 0 flags that no file will be written
C     Switch syntax:  --wrsj[:fn=name][:scl=#]
      ifirsj = 0
      if (cmdopt('--wrsj',6,0,jfilnm)) then

        out2 = jfilnm(7:)
        if (jfilnm(1:5) == '-wrsj') out2 = jfilnm(6:)
C     This is default name for supercell rsj file
        jfilnm = 'rsj2'

        dc = out2(1:1)
        sclrsj = 1
        if (dc /= ' ') then
C     ... Return here to resume parsing for arguments
          j2 = 0
   10     continue
          j2 = j2+1
          if (out2(j2:j2) == dc) goto 10
          j1 = min(len(out2),j2)
          call nwordg(out2,0,dc//' ',1,j1,j2)
          if (j2 >= j1) then
            if (.false.) then
            elseif (out2(j1:j1+3) == 'scl=')  then
              j = 0
              i = parg('scl=',4,out2(j1:),j,len(out2(j1:)),dc//' ',1,1,
     .          k,sclrsj)
              if (i <= 0) goto 11

            elseif (out2(j1:j1+2) == 'fn=')  then
              jfilnm = out2(j1+5:j2)
            else
              goto 11
            endif
            goto 10
   11       continue
            call rxs('supcel: failed to parse wrsj option', out2)
          endif
        endif

        call info0(20,1,0,
     .    ' ... Writing supercell r.s. hamiltonian to file '//jfilnm)

C     ... Setup for rsj file read
        ifirsj = 0
        mxcsiz = 0
        nbasp = nbas
        if (fxst('rsj') /= 1) then
          call rx('supcel: specified --wrsj but no file rsj')
        endif
        ifirsj = fopna('rsj',-1,0)
        allocate(aamom(nbas))
        nfilj = iorsj(0,ifirsj,nbas,alat,plat,aamom,xx,xx,xx,xx,x)

C       This makes iax table and rtab suitable for iorsj
        allocate(pos(3,nbas))
        call sitepack(s_site,1,nbas,'pos',3,xx,pos)

        allocate(ntab(nbasp+1))
        call pairs(nbas,nbasp,1d0,plat,x(1:1)/2*(1+1d-6),pos,
     .    [-1],3,-1,[0],nttab,ntab,iax,mxcsiz)
        allocate(rtab(3,nttab))
        call mkrtab(000,1d0,plat,pos,iax,nttab,pos,rtab)

C       Read pair hamiltonian for entries matching iax table
        allocate(hpair(nttab))
        nfilm = iorsj(03,ifirsj,nbas,xx,plat,aamom,iax,rtab,ntab,hpair,x)

        if (nfilj/=nfilm) call info2(20,0,0,'%8f(warning) only'//
     .    ' matched %i out of %i file entries',nfilm,nfilj)
        if (nfilm/=nttab) then
          call info2(20,0,0,'%8f(warning) iax table has %i pairs'//
     .      ' but matched only %i file entries',nttab,nfilm)

C         Purge iax table of unused entries
          call psup3(nlst,ilstt,iax,nttab,nfilm)
          call symiax(1,plat,nbas,xx,xx,xx,0,ntab,iax,nttab,mxcsz)
          if (nttab /= nfilm) then
            call fexit2(-1,111,' Exit -1 : iorsj: bug in symiax:'//
     .        ' nttab=%i, nfilm=%i',nttab,nfilm)
          endif

C         Remake rtab
          call mkrtab(000,1d0,plat,pos,iax,nttab,pos,rtab)

C         Re-read pair hamiltonian with new iax table
          nfilm = iorsj(03,ifirsj,nbas,xx,plat,aamom,iax,
     .      rtab,ntab,hpair,x)
          if (nttab /= nfilm) call rx('supcel: bug in iorsj')

        endif

C       Close input file rsj
        call fclose(ifirsj)

C       Open supercell file; write header
        allocate(amom2(nbx))
        call cpvprm(0,1,nbx,imap,aamom,amom2)
        ifirsj = fopna(jfilnm,-1,0)
        i = iorsj(0,-ifirsj,nbx,alat,plx,amom2,xx,xx,xx,xx,xx)

        deallocate(aamom,pos)
      endif

      if (.not. (cmdopt('--wsite',7,0,outs)) .and.
     .    .not. (cmdopt('--wpos',5,0,outs))  .and.
     .    .not. (cmdopt('--sqs',5,0,outs)) ) return

C --- Create supercell s_ssite ---
      allocate(s_ssite(nbx))

      nspec = 0
      nbx = 0
      ilst = 0
      do  ib = 1, nbas
        if (nlst > 0) then
   12     continue
          if (ilst >= nlst) cycle
          if (ilstt(1+ilst) < ib) then
            ilst = ilst+1
            goto 12
          endif
          if (ilstt(1+ilst) /= ib) cycle
        endif

        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        alabl = s_spec(is)%name
        nspec = max(nspec,is)
        do  k = 1, nkd

          if (lkeep) then
            nbx = nbx+1
            lok = 1
            dd = tau
          else
            call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
          endif
          if (lok == 1) then

C       ... Copy parameters from site ib to site nbx
C           12 Apr 03 ssite(1) may not be struc size if 1st ib<>1
            s_ssite(nbx)%spec = is
            s_ssite(nbx)%pos = dd
            s_ssite(nbx)%class = s_site(ib)%class
            s_ssite(nbx)%force = s_site(ib)%force
            s_ssite(nbx)%vel = s_site(ib)%vel
            s_ssite(nbx)%eula = s_site(ib)%eula
            s_ssite(nbx)%pl = s_site(ib)%pl
            s_ssite(nbx)%relax = s_site(ib)%relax
            s_ssite(nbx)%vshft = s_site(ib)%vshft
            s_ssite(nbx)%pnu = s_site(ib)%pnu
            s_ssite(nbx)%pz = s_site(ib)%pz
            s_ssite(nbx)%v0 => s_site(ib)%v0
            s_ssite(nbx)%v1 => s_site(ib)%v1

C       ... Write rsj to supercell file
            if (ifirsj /= 0) then
              call psup2(ifirsj,ib,nbx,ntab,rtab,sclrsj,hpair)
              deallocate(ntab,rtab,hpair)
            endif

C       ... PL index
            if (cmdopt('--pl:',5,0,out2)) then
              call prjpos(2,ixv,xx,plat,1,dd,x)
              call mapdat(1,out2(6:),' ',[m],1,3,x,[xx])
C     call spack(2,'site pl',ssite,1,m,j1,j2,k)
C     call dvset(sssite,m+(nbx-1)*ndsite,m+(nbx-1)*ndsite,x)
              s_ssite(nbx)%pl = nint(x(1))
            endif
          endif
          if (lkeep) exit
        enddo
      enddo

C --- Permute site strux if sorted ---
      if (lsort .or. lsplat .or. lfirst .or. lswap .or. lcycle) then
        call siteperm(nbx,iprm,s_ssite)
        allocate(iwk(nbx))
        call ivprm(1,nbx,imap,iwk,iprm,1)
        deallocate(iwk,iprm)
      endif

      allocate(pos(3,nbx))
      call sitepack(s_ssite,1,nbx,'pos',3,xx,pos)
C     call prmx('pos',pos,3,3,nbx)

C --- Displacement around an atom ---
      if (cmdopt('-disp',5,0,out2) .or. cmdopt('--disp',6,0,out2)) then
        allocate(dpos(3,nbx,1)); call dpzero(dpos,3*nbx)
        dc = out2(6:6)
        if (out2(1:6) == '--disp') dc = out2(7:7)
        call wrdsg(out2,1,dc//' ',nw)
        if (nw < 3) call rx('supcel: --disp requires 2 arguments')

C   ... Read displacements file
        call wordg(out2,0,dc//' ',2,j1,j2)
        ifi = fopn(out2(j1:j2))
        rewind ifi
C   ... Determine size of file
        nn = 0
        nc = 0
        j = rdm(ifi,0,0,' ',xx,nn,nc)
        if (j < 0) call rxs
     .    ('supcel: failed to read displacements file ',out2(j1:j2))
        if (nc /= 6) call rxs('supcel: displacements file '//
     .    'apparently does not have 6 columns: ',out2(j1:j2))
        allocate(disp(nn*nc))
        rewind ifi
        j = rdm(ifi,0,nn*nc,' ',disp,nn,nc)
        call psup5(nn,disp,rcut)
        call wordg(out2,0,dc//' ',3,i,j)
        call info(30,1,0,' ... Displacements from file '//
     .    out2(j1:j2)//' for site list='//out2(i:j)//'. rcut=%d',rcut,0)

        nlst1 = 0
        allocate(lst1(nbx)); call iinit(lst1,nbx)
        call wordg(out2,0,dc//' ',3,j1,j2)
        if (j2 < j1) call rx('supcel: --disp requires 2 arguments')
        allocate(ips(nbx),z(nspec,1))
        ips(:) = s_site(1:nbx)%spec
        z(:,1) = s_spec(1:nspec)%z
        call baslst(1,11,out2(j1-1:),j2,ips,nbx,slabl,z,0,' ',x,nlst1,lst1)
        deallocate(ips,z)

C   ... Make iax table and rtab
        mxcsiz = 0
C       call setpr(50)
        call pshpr(iprint()-10)
        allocate(ntab(nbx+1))
        call pairs(nbx,nbx,1d0,plx,[rcut/2*(1+1d-6)],pos,
     .    [-1],3,-1,[0],nttab,ntab,iax,mxcsiz)
        allocate(rtab(3,nttab))
        call mkrtab(000,1d0,plx,pos,iax,nttab,pos,rtab)
        call poppr
        call ivset(ixv,1,3,2)
        call shorps(nttab,plx,ixv,rtab,rtab)

C       Get displacements
        call psup4(iax,nttab,nlst1,lst1,rtab,disp,nn,dpos)

C       Repack into to s_ssite
        call daxpy(3*nbx,1d0,dpos,1,pos,1)
        call sitepack(s_ssite,1,nbx,'-pos',3,xx,pos)
        deallocate(ntab,dpos,rtab,disp,lst1)
        deallocate(iax)
      endif

C --- Rotation of lattice ... no use for it, since can be done on original lattice ---
C      i = 6
C      if (cmdopt('--rot=',i,0,out2)) then
C        call a2rotm(out2(i+1:),.false.,iprint()-10,platm) ! platm used as a dummy here
C        call lattdf(2,platm,plx,nbx,pos,0,[0d0],[0d0])
C      endif

C --- Uniform translation of positions within lattice ---
      i = 8
      if (cmdopt('--xshft=',i,0,out2) .or. cmdopt('--xshftx=',i+1,0,out2)) then
        ltmp = out2(1:9) == '--xshftx='
        if (ltmp) i=i+1
        is = a2vec(out2,len(out2),i,4,', ',2,2,3,it,tau)
        if (is /= 3) call rx('blm failed to parse '//trim(out2))
        dd = tau
        if (ltmp) call dgemm('N','N',3,1,3,1d0,plx,3,tau,3,0d0,dd,3)
        forall (ib = 1:nbx) pos(1:3,ib) = pos(1:3,ib) + dd(1:3)
        call sitepack(s_ssite,1,nbx,'-pos',3,xx,pos)
      endif

C --- SQS ---
      if (cmdopt('--sqs',5,0,out2)) then

C       Default values for SQS
        r2max = 1.3d0           ! range of pair table
        r3max = 0.75d0          ! maximum perimeter for triples
        alf2 = 1                ! decay parameter for pair weighting
        alf3 = 1                ! decay parameter for triplet weighting
        scal2 = 1               ! prefactor for pair weighting
        scal3 = 1               ! prefactor for triplet weighting
        nspsqs = 2              ! Number of species to mix.  Must be 2 for now
        sqspid(1) = 1           ! First alloy species
        sqspid(2) = 2           ! Second alloy species
        r3mode = 1              ! specifies how clusters are included

        dc = out2(6:6); outs = ' '
        if (dc /= ' ') then
C     ... Return here to resume parsing for arguments
          j2 = 6
          do
            j2 = j2+1
            if (out2(j2:j2) == dc) cycle
            j1 = min(len(out2),j2)
            call nwordg(out2,0,dc//' ',1,j1,j2)
            if (j2 < j1) exit
            if (.false.) then
            elseif (out2(j1:j1+4) == 'seed=')  then
              if (j1+5 > j2) call rx('sqs: bad argument to seed')
              i = j1+4
              j = a2vec(out2,j2,i,2,dc//' ',2,3,1,iv,k)
              if (j /= 1) call rxs('sqs: cannot parse ', out2(j1:j2))
              call ran1in(k)
              call awrit1('%a  seed=%i',outs,len(outs),0,k)
              cycle
            elseif (out2(j1:j1+5) == 'r2max=')  then
              if (j1+6 > j2) call rx('sqs: bad argument to r2max')
              i = j1+5
              j = a2vec(out2,j2,i,4,dc//' ',2,3,1,iv,r2max)
              if (j /= 1) call rxs('sqs: cannot parse ', out2(j1:j2))
              cycle
            elseif (out2(j1:j1+5) == 'r3max=')  then
              if (j1+6 > j2) call rx('sqs: bad argument to r3max')
              i = j1+5
              j = a2vec(out2,j2,i,4,dc//' ',2,3,1,iv,r3max)
              if (j /= 1) call rxs('sqs: cannot parse ', out2(j1:j2))
              cycle
            elseif (out2(j1:j1+4) == 'alf2=')  then
              if (j1+5 > j2) call rx('sqs: bad argument to r3max')
              i = j1+4
              j = a2vec(out2,j2,i,4,dc//' ',2,3,1,iv,alf2)
              if (j /= 1) call rxs('sqs: cannot parse ', out2(j1:j2))
              cycle
            elseif (out2(j1:j1+6) == 'r3mode=')  then
              if (j1+7 > j2) call rx('sqs: bad argument to seed')
              i = j1+6
              j = a2vec(out2,j2,i,2,dc//' ',2,3,1,iv,r3mode)
              if (j /= 1) call rxs('sqs: cannot parse ', out2(j1:j2))
              cycle
            else
              call rxs('sqs: failed to parse argument, ',out2(j1:j2))
            endif
          enddo
        endif

C   ... Count number of sites in the species and identify first site
C       Note : algorithm below expects contiguous SQS sites, starting at ib1
C        ib1 = 0; nsqs = 0
C        do  ib  = 1, nbx
C          is = s_ssite(ib)%spec
C          if (is /= sqspid(1) .and. is /= sqspid(2)) cycle
C          nsqs = nsqs + 1
C          if (ib1 == 0) ib1 = ib
C        enddo
C        call info5(20,1,1,
C     .    ' Initiate SQS algorithm for %i sites starting at %i.  r2max=%;4d  r3max=%;4d  r3mode=%i'//
C     .    trim(outs),nsqs,ib1,r2max,r3max,r3mode)

        call info5(20,1,1,' Initiate SQS algorithm : r2max=%;4d  r3max=%;4d  r3mode=%i'//
     .    trim(outs),r2max,r3max,r3mode,0,0)

C   ... Make iaxs, ntabs, rtab for all pairs with species 1 or 2 at the origin
C       Lattice vectors are plx and positions pos
        mxcsiz = 0
C       call setpr(50)
        call pshpr(iprint()-10)
        allocate(ntab(nbx+1))
        if (allocated(iax)) deallocate(iax)
        call pairs(nbx,nbx,alat,plx,[4*r2max*s_lat%avw/2],pos,[-1],3,-1,[0],nttab,ntab,iax,mxcsiz)
        allocate(rtab(3,nttab))
        call mkrtab(000,1d0,plx,pos,iax,nttab,pos,rtab)
        call poppr
        call ivset(ixv,1,3,2)

C       Pass 1: Count number of elements in each region
C       Pass 2: fill out ntabs, iaxs, e2 = length of connecting vectot
C       Bug: loop may assume SQS atoms are contiguous in Table (not sure)
        do  ipass = 1, 2
          nttabs = 0; ib0 = 0
          do  i = 1, nttab
            iax(7+niax*(i-1)) = 0
            iax(8+niax*(i-1)) = 0
            ib = iax(1+niax*(i-1)); is = s_ssite(ib)%spec
            jb = iax(2+niax*(i-1)); js = s_ssite(jb)%spec
C           Both elements in pair must belong to species 1 or 2
            if (is /= sqspid(1) .and. is /= sqspid(2)) cycle
            if (js /= sqspid(1) .and. js /= sqspid(2)) cycle
            xx = dlength(3,rtab(1,i),1)
            if (xx == 0 .or. xx > r2max+tol) cycle
            nttabs = nttabs+1
            if (ib > nbx) call rx('something wrong with iax table')
C            if (ipass == 1 .and. iprint() >= 50) then
C              write(stdo,'(i4,i4,3f11.6,f9.4)') ib,jb,rtab(1:3,i),xx
C            endif
            if (ipass == 2) then
              if (((ib/=ib0 .and. ib0/=0) .or. nttabs == nttab_save) .and. iprint()>=40) then
                write(stdo,345) ib0,ntab(ib0+1)-ntab(ib0),ntabs(ib0)-ntabs(ib0-1)+isw(ib==ib0)
  345           format(' supcel, ib=',i3,':  pair table pared from',i4,' to',i4)
              endif
              e2(nttabs) = xx
              ntabs(ib) = nttabs
              call icopy(niax,iax(1+niax*(i-1)),1,iaxs(1,nttabs),1)
              ib0 = ib          ! last participating site
            endif
          enddo                 ! Loop over pairs
          nttab_save = nttabs

C     ... Pass 1 : allocate pair tables.  Pass 2: tidy up ntabs
          if (ipass == 1) then
            allocate(e2(nttabs))
            allocate(iaxs(niax,nttabs),ntabs(0:nbx)); ntabs = 0
          elseif (ipass == 2) then
            call ivset(ntabs(1),ib0+1,nbx,nttabs) ! cumulative pairs for sites above final
          endif
        enddo                   ! ipass

C   ... Make initial ipss table, and ibs table
        allocate(ipss(nbx))
        ipss(:) = s_ssite(:)%spec
        nbs = 0
        do  ib = 1, nbx
          if (ntabs(ib)-ntabs(ib-1) <= 0) cycle
          nbs = nbs+1
        enddo
C        if (iprint()>=40) then
C          write(stdo,345) nbs,ntab(nbs+1)-ntab(nbs),ntabs(nbs)-ntabs(nbs-1)
C        endif

C   ... Number and list of all participating sites
        allocate(sigms(nbx))
        call mcasms(0,1,nbx,ipss,sqspid,sigms) ! Make sigms
        allocate(ibs(nbs))
        nbs = 0
        do  ib = 1, nbx
          if (ntabs(ib)-ntabs(ib-1) <= 0) cycle
          nbs = nbs+1
          ibs(nbs) = ib
        enddo

C   ... Count number of triplets nttab3
        call ppair5(1,nbx,plx,pos,-1d0,ntabs,iaxs) ! Set up iax(7,:) in preparation for pair3c
        allocate(rtabs(3,nttabs))
        call mkrtab(000,1d0,plx,pos,iaxs,nttabs,pos,rtabs)
        call pshpr(1)
        nttab3 = 0
        do  ib = 1, nbx
          nclus3 = ntabs(ib)-ntabs(ib-1) ! Consider all triplets in full cluster
          if (nclus3 <= 0) cycle
          call pair3c(10*r3mode,ib,iaxs,ntabs,nclus3,rtabs,r3max,n3c,it,it)
          nttab3 = nttab3+n3c
        enddo
        call poppr
        allocate(iax3c(3,nttab3),iax3(niax,nttab3))
C       Identify all triplets entering into 3-body correlation energy, as entries into iaxs table
C       List is sorted now ... maybe not needed
        allocate(ntab3s(0:nbx)); ntab3s = 0
        nttab3 = 0
C       call pshpr(iprint()-10)
        do  ib = 1, nbx
          nclus3 = ntabs(ib)-ntabs(ib-1) ! Consider all triplets within cluster defining pairs
          ntab3s(ib) = nttab3
          if (nclus3 <= 0) cycle
          call pair3c(201+10*r3mode,ib,iaxs,ntabs,nclus3,rtabs,r3max,n3c,iax3c(1,1+nttab3),iax3(1,1+nttab3))
          nttab3 = nttab3 + n3c
          ntab3s(ib) = nttab3
        enddo
C       call poppr
        deallocate(iax3c)
        call info5(30,1,0,' SQS:  %i pairs, %i triplets among %i SQS sites',nttabs,nttab3,nbs,0,0)

C   ... Convert pair length into pair correlation energy table e2
C       For now, expresssion is e2 = scal2*exp(alf2*d)
        call idmxmn(nttabs,e2,1,j1,j2); dd(1) = e2(j1); dd(2) = e2(j2)
        do  i = 1, nttabs
          xx = e2(i)
          e2(i) = scal2*dexp(-alf2*xx)
        enddo
        call info5(10,0,0,'%7fsmallest pair length %,4;4d  energy %,4;4d%N%7flargest  pair length %,4;4d  energy %,4;4d',
     .    dd(1),e2(j1),dd(2),e2(j2),0)
        call info2(30,0,0,'%7fpair energy e2 scales as %d*exp(-%d*d)',scal2,alf2)
        if (iprint() >= 50) then
          write(stdo,"(/' Pairs retained for SQS'/,'  ib  jb',9x,'--- r(jb)-r(ib) ---',10x,'d',10x,'e2')")
          do  i = 1, nttabs
            ib = iaxs(1,i); jb = iaxs(2,i)
            write(stdo,'(i4,i4,3f11.6,f9.4,f12.6)') ib,jb,rtabs(1:3,i),dlength(3,rtab(1,i),1),e2(i)
          enddo
        endif

C   ... Triplet correlation energy table e3
        allocate(e3(max(nttab3,1))); e3(1) = 0
        do  i = 1, nttab3
          kb = iax3(1,i)
          ib = iax3(2,i)
          jb = iax3(3,i)
          ip = iax3(4,i)
          jp = iax3(5,i)
          dd(1) = dlength(3,rtabs(1,ip),1)
          dd(2) = dlength(3,rtabs(1,jp),1)
          dd(3) = dlength(3,rtabs(1:3,ip)-rtabs(1:3,jp),1)
          e3(i) = dsum(3,dd,1)
        enddo
        if (nttab3 == 0) then
          j1 = 0 ; j2 = 0
        else
          call idmxmn(nttab3,e3,1,j1,j2)
        endif
        dd(1) = e3(j1); dd(2) = e3(j2)
        forall (i = 1:nttab3) e3(i) = scal3*dexp(-alf3*e3(i))
        call info5(10,0,0,'%7fsmallest triplet length %,4;4d  energy %,4;4d%N%7flargest  triplet length %,4;4d  energy %,4;4d',
     .    dd(1),e3(j1),dd(2),e3(j2),0)

C   ... Set up site permutation table; populated with unpermuted lattice
        allocate(iprm(nbx))
        forall (i=1:nbx) iprm(i) = i

C   ... Permute lattice into SQS
C#ifdefC DEBUG
C        iprm(2) = 13
C        iprm(13) = 2
C        iprm(4) = 14
C        iprm(14) = 4
C        iprm(1) = 10
C        iprm(10) = 1
C#endif
        if (s_move%tstot == NULLI)
     .    call rx('supcel: no annealing time specified : annealing requires DYN_TTOT')
        call sqssim(s_move,nbx,ipss,iprm,sigms,ntabs,iaxs,e2,ntab3s,iax3,e3,esqs)

C       Update imap
        allocate(iwk(nbx))
        call icopy(nbx,imap,1,iwk,1)
        forall (i = 1:nbx) imap(i,1) = iwk(iprm(i))
        call icopy(nbx,imap(1,2),1,iwk,1)
        forall (i = 1:nbx) imap(i,2) = iwk(iprm(i))
        deallocate(iwk)

        call info0(20,2,0,' SQS completed ... permute species and write site permutations to file siteperm')
        ifi = fopna('siteperm',-1,0)
        rewind (ifi)
        call awrit4('%% rows %i cols 2    %i SQS sites out of %i  e=%d',' ',80,ifi,nbx,nbs,nbx,esqs)
C       call info2(0,0,0,'%n:1i',nbx,iprm)
        do  i = 1, nbx
          write(ifi,"(2i5)") i, iprm(i)
          s_ssite(i)%spec = ipss(iprm(i))
        enddo
        if (.not. cmdopt('--wsite',7,0,outs)) return
      endif

C --- Read ASA rst file, write for supercell ---
      if (cmdopt('--rsta',5,0,out2)) then
        call info0(30,1,0,' Creating supercell restart file rstax ...')
        nclass = s_ctrl%nclass
        nl = s_ctrl%nl
        nsp = s_ctrl%nspin
        allocate(pnu(nl,nsp,nbas),qnu(3,nl,nsp,nbas))
        call ptr_pot(s_pot,8+1,'ves',nbas,0,xx)
        ifi = fopna('rsta',-1,0)
        call asars(0,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,pnu,qnu,.false.,ifi)
        if (cmdopt('--rsta,amom',10,0,out2)) then
          do  ib  = 1, nbas
            x = s_ssite(ib)%eula
            if (dsum(nl,qnu(1,1,1,ib),3)-dsum(nl,qnu(1,1,2,ib),3) < 0) then
              call dswap(nl,pnu(1,1,ib),1,pnu(1,2,ib),1)
              call dswap(3*nl,qnu(1,1,1,ib),1,qnu(1,1,2,ib),1)
              if (x(1) > 0) then
                x(1) = x(1)-pi
              else
                x(1) = x(1)+pi
              endif
              x(2) = pi-x(2)
            endif
            s_ssite(ib)%eula = x
C     111       format(i4,2f12.6,3x,3f12.6)
C     print 111, ib,
C     .        dsum(nl,qnu(1,1,1,ib),3)+dsum(nl,qnu(1,1,2,ib),3),
C     .        dsum(nl,qnu(1,1,1,ib),3)-dsum(nl,qnu(1,1,2,ib),3),x
          enddo
        endif
        ifi = fopna('rstax',-1,0)
        x(1) = dglob('nbasp',dble(nbx),1)
        call asars(0,s_ctrl,s_ssite,s_spec,s_lat,s_pot,s_bz,pnu,qnu,.false.,-ifi)
      endif

C --- Modify positions as excursions around a reference site file ---
      allocate(slabll(mxspec))
      if (cmdopt('--refsite',9,0,outs)) then
        dc = outs(10:10)

        ltmp = wordsw(outs,dc,'plx','',j1) > 0 ! Shifts in terms of multiples of latvec
        lswap = wordsw(outs,dc,'swap','',j1) > 0 ! swap role of site and reference files

        if (wordsw(outs,dc,'scl=','',j1) == 0) call rx('--refsite requires tag "scl="')! Scale for shifts
        call nwordg(outs,0,dc//' ',1,j1,j2)
        i = 0
        if (a2vec(outs(j1:j2),len_trim(outs(j1:j2)),i,4,dc//' ',2,3,1,iv,sclrsj) /= 1)
     .    call rx('supcel failed to parse scl=')

        i = wordsw(outs,dc,'fn=','',j1) + 4
        if (i < 5 .or. dc == ' ') call rx('supcel: --refsite requires file name')
        j1 = i
        call nwordg(outs,0,dc//' ',1,j1,j2)
        ifi = fopnx(outs(j1:j2),70,-1,-1)
        if (ifi /= 1) then
          call info0(1,0,0,' supcel: no file "'//outs(j1:j2)//'.ext" ... try "'//outs(j1:j2)//'"')
          ifi = fopng(outs(j1:j2),-1,1)
        else
          ifi = fopna(outs(j1:j2),-1,0)
        endif

C   --- Scale positions by fac*(change in pos reference site file)
        i = 135002     ! (128+7)*1000 + 2 : read alatl,nbasl,nspecl,slabbl file ifi
        nspecl = 0
        rewind ifi
        j = iosite(i,3d0,0,'',ifi,slabll,alatl,platl,nbasl,nspecl,xx,xx,xx,xx,xx,xx,xx)
        call info5(1,1,0,' Scaled shifts relative to ref "'//outs(j1:j2)//
     .    '": read %i species, %i sites  scale=%g%?#n# (latvec units)##'//
     .    '%N compare lattice vectors alat*plat:'//
     .    '%N%15fplat%29fref%30fdiff',
     .    nspecl,nbasl,sclrsj,isw(ltmp),5)
        do  k = 1, 3
          call info5(2,0,0,'%3;10,5D   %3;10,5D   %3;10,5D',
     .      alat*plx(1:3,k),alatl*platl(1:3,k),alat*plx(1:3,k)-alatl*platl(1:3,k),4,5)
        enddo
        xx = avwsr(platl,alatl,x(1),1); xx = avwsr(plx,alat,x(2),1)
        call info5(2,0,0,' supercell volume %g  reference volume %g  ratio %;10F',
     .    x(1),x(2),x(1)/x(2),4,5)
        call dinv33(platl,1,platm,xx)

C   ... Read data from site file
        allocate(s_sitel(nbasl),s_specl(nspecl))
        forall (is=1:nspecl) s_specl(is)%name = slabll(is)
        rewind ifi
        j = iosits(8072,3d0,0,'',ifi,slabll,alatl,platl,nbasl,nspecl,s_specl,s_sitel)

        call info2(1,0,0,'%N Site positions and shifts in %?#n#multiples of latvec#Cartesian coordinates#:'//
     .    '%N%14fscell%27fref%25fuscaled shift%22fnew pos',1,2)
        allocate(dpos(3,nbx,4))
        do  ib = 1, nbasl
          if (ltmp) then
            call dgemm('T','N',3,1,3,1d0,qlx,3,s_ssite(ib)%pos,3,0d0,dpos(1:3,ib,1),3)
            call dgemm('T','N',3,1,3,1d0,platm,3,s_sitel(ib)%pos,3,0d0,dpos(1:3,ib,2),3)
          else
            dpos(1:3,ib,1) = alat*s_ssite(ib)%pos
            dpos(1:3,ib,2) = alatl*s_sitel(ib)%pos
          endif
          if (lswap) then
            dd = dpos(1:3,ib,1); dpos(1:3,ib,1) = dpos(1:3,ib,2); dpos(1:3,ib,2) = dd
          endif
          dpos(1:3,ib,3) = dpos(1:3,ib,1) - dpos(1:3,ib,2)
          dpos(1:3,ib,4) = sclrsj*dpos(1:3,ib,3) + dpos(1:3,ib,2)

          call info5(2,0,0,'%3;10,6D  %3;10,6D  %3;10,6D  %3;10,6D',dpos(1:3,ib,1),dpos(1:3,ib,2),dpos(1:3,ib,3),dpos(1:3,ib,4),5)

          if (ltmp) then
            call dgemm('N','N',3,1,3,1d0,plx,3,dpos(1:3,ib,4),3,0d0,pos(1,ib),3)
          else
            pos(1:3,ib) = dpos(1:3,ib,4)/alat
          endif
          s_ssite(ib)%pos(1:3) = pos(1:3,ib)

        enddo
      endif

C --- Write site file ---
      lio = 0; lfirst = cmdopt('--wsite',7,0,out2)
      if (out2(1:8) == '--wsitex') then
        lio = 1000*(1+2+4+8+16+32) + 1 + 10
        j2 = 9
      elseif (out2(1:7) == '--wsite') then
        lio = 1000*(1+2+4+8+16+32) + 1
        j2 = 8
      endif

      if (lio /= 0) then
        outs = 'site'           ! Default file name
        dc = out2(j2:j2)
        lwmap = .false.         ! do not append mapping data
        if (dc == '=') then
          outs = out2(j2+1:)
C  ... Loop until arguments are exhausted
        elseif (dc /= ' ') then
          do  while (.true.)
            j2 = j2+1
            if (out2(j2:j2) == dc) cycle
            j1 = min(len(out2),j2)
            call nwordg(out2,0,dc//' ',1,j1,j2)
            if (j1 > j2) exit
            if (.false.) then
            elseif (out2(j1:j1+4) == 'short')  then
              lio = lio - 1000*iand(lio/1000,32)
            elseif (out2(j1:j1+2) == 'map')  then
              lwmap = .true.
            elseif (out2(j1:j1+2) == 'fn=')  then
              outs = out2(j1+3:)
            elseif (out2(j1:j1+4) == 'quad1')  then
              call info0(30,0,0,' ... shifting sites to first quadrant')
              do  ib = 1, nbx
C               posp = posc (plat+)^-1  and  posp+ = (plat)^-1 posc+
                call dgemm('T','N',3,1,3,1d0,qlx,3,s_ssite(ib)%pos,3,0d0,tau,3)
                ltmp = .true.
                do  while (ltmp)
                  ltmp = tau(1) < -tol .or. tau(2) < -tol .or. tau(3) < -tol
                  if (tau(1) < -tol) tau(1) = tau(1)+1
                  if (tau(2) < -tol) tau(2) = tau(2)+1
                  if (tau(3) < -tol) tau(3) = tau(3)+1
                enddo
C               posc = posp (plat+) and  posc+ = (plat) posp+
                call dgemm('N','N',3,1,3,1d0,plx,3,tau,3,0d0,s_ssite(ib)%pos,3)
              enddo
            else
              call rxs('supcel: failed to parse --wsite option', out2)
            endif
          enddo
        endif

        if (iosits(lio,3d0,0,outs,ifi,slabl,alat,plx,nbx,nspec, s_spec,s_ssite) < 0)
     .    call rx('failed to write ssite')
        if (lwmap) then
          write(ifi,"('#map',24i4)") imap(:,1)
        endif
      endif

C --- Write positions to file ---
      if (cmdopt('--wpos=',7,0,outs)) then
        call iopos(.true.,0,outs(8:),nbx,pos,s_site)
      endif

      deallocate(ilstt,pos)
      return
   13 format(' (warning) supercell has incommensurate lattice')
C     999 continue

      end
      subroutine psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
C- Kernel called by supcel: to decide whether a vector should be added
C ----------------------------------------------------------------------
Ci Inputs
Ci   nmark :number of specified multiples of lattice vectors, if any
Ci         :0 => multiples of lattice vectors were not specified
Ci   imark :points to entry of specified lattice vectors in dlat
Ci   tau   :position vector for this site
Ci   k     :which of the lattice vectors in dlat to check
Ci   dlat  :a list of lattice vectors, of which dlat(k) is to be used here.
Ci   plx   :(passed to shorps; used if vectors are to be shortened)
Ci   qlx   :reciprocal lattice vectors to supercell
Ci   ixv   :(passed to shorps; used if vectors are to be shortened)
Ci         :0 suppresses shifts along plat(j)
Ci         :1 shifts to unit cell at origin (pos in 1st quadrant)
Ci         :2 shifts to minimize length of pos
Co Outputs
Co   lok   :0 do not add this vector to the list
Co         :1 add this vector to the list
Co   dd    :connecting vector
Co   nbx   :incremented by one if this vector should be included
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nmark,imark(nmark),lok,nbx,ixv(3),k
      double precision dlat(3,k),tau(3),plx(3,3),qlx(3,3),dd(3)
C ... Local parameters
      double precision x(3)
      integer m
      logical cmdopt
      character*72 outs

      dd(1) = dlat(1,k)+tau(1)
      dd(2) = dlat(2,k)+tau(2)
      dd(3) = dlat(3,k)+tau(3)
      if (nmark > 0) then
        lok = 0
        do  m = 1, nmark
          if (imark(m) == k) lok = 1
        enddo
      else
        lok = 1
        do  m = 1, 3
          x(m) = dd(1)*qlx(1,m)+dd(2)*qlx(2,m)+dd(3)*qlx(3,m)
          if (x(m) < -0.001d0 .or. x(m) > 0.999d0) lok = 0
        enddo
      endif
      if (lok == 1) then
        nbx = nbx+1
        if (cmdopt('-shorten',8,0,outs)) call shorps(1,plx,ixv,dd,dd)
        if (cmdopt('--shorten',9,0,outs)) call shorps(1,plx,ixv,dd,dd)
      endif
      end

      subroutine psup2(ifi,ib,ib2,ntab,rtab,sclrsj,rsj)
C- Kernel called by supcel: write r.s. pairwise ham for this site
C  Note that jb2 of supercell is not known, so we just write zero
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   ib    :site index of original cell
Ci   ib2   :site index of supercell
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   rsj
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Nov 02
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,ib,ib2,ntab(ib+1)
      double precision sclrsj
C     integer niax
C     parameter (niax=10)
C     integer iax(niax,*)
      double precision rsJ(*),rtab(3,*)
C ... Local parameters
      integer it

      do  it = ntab(ib)+1, ntab(ib+1)
C        call awrit6('%,4i%,4i%;12,7D%;12,7D%;12,7D%;12,4D',' ',80,ifi,
C     .    ib2,0,rtab(1,it),rtab(2,it),rtab(3,it),sclrsj*rsj(it))
        call awrit6('%,6i%,6i%;12,7D%;12,7D%;12,7D %;6g',' ',80,ifi,
     .    ib2,0,rtab(1,it),rtab(2,it),rtab(3,it),sclrsj*rsj(it))
      enddo
      end

      subroutine psup3(nlst,ilst,iax,nttab,n)
C- Kernel called by supcel to purge pairs outside of site list
C  used when a subset of all sites is specified
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlst  :number of sites in site list
Ci   ilst  :site list
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Co Outputs
Co   n     :number of pairs after purging
Co         :NB: n is not set if nlst=0
Cr Remarks
Cr
Cu Updates
Cu   26 Nov 02 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax,nlst,ilst(nlst),nttab
      parameter (niax=10)
      integer iax(niax,nttab)
C ... Local parameters
      integer i,ib,jb,low,n

      if (nlst == 0) return
      n = 0
      do  i = 1, nttab
        ib = iax(1,i)
        low = 1
        call hunti(ilst,nlst,ib,0,low)
        if (low >= nlst) then
C         print *, 'i, iax',i,iax(1,i)
          iax(1,i) = 0
        elseif (ilst(low+1) /= ib) then
C         print *, 'i, iax',i,iax(1,i)
          iax(1,i) = 0
        else
          jb = iax(2,i)
          call hunti(ilst,nlst,jb,0,low)
          if (low >= nlst) then
C           print *, 'i, iax',i,iax(1,i),iax(2,i)
            iax(1,i) = 0
          elseif (ilst(low+1) /= jb) then
C           print *, 'i, iax',i,iax(1,i),iax(2,i)
            iax(1,i) = 0
          else
            n = n+1
          endif
        endif
      enddo

      call info2(20,0,0,' supcel: site-list further reduced pair'//
     .  ' table to %i pairs',n,0)

      end

      subroutine mapdat(nexpr,expr,sincl,iwk,nr,nc,dat,dat2)
C- Replace columns of dat with algebraic expressions of them
      implicit none
C     Passed parameters
      integer nexpr,nr,nc,iwk(nr)
      character*(*) expr(1), sincl
C     dat2 is not used unless poins are excluded
      double precision dat(nr,max(nexpr,nc)),dat2(nr,nexpr)
C     Local variables
      integer ir,iv0,ival,i,j,ii,jr
      character*4 xn, outs*80
      logical a2bin,logi

      call numsyv(iv0)
      jr = 0
      do  ir = 1, nr
        call clrsyv(iv0)
        call addsyv('i',dble(ir),ival)

C   ... Load data table
        do  j = 1, nc
          ii = 1
          xn = 'x   '
          call bin2a('',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ir,j),ival)
        enddo

C   ... Exclude points if not satisfy sincl
        logi = .true.
        if (sincl /= ' ') then
          j = 0
          if (.not. a2bin(sincl,logi,0,0,' ',j,-1))
     .      call rx('mapdat:  error parsing sincl')
        endif

C   ... Put expressions of these vars into dat
        do  i = 1, nexpr
          j = 0
          if (.not. a2bin(expr(i),dat(ir,i),4,0,' ',j,-1)) then
            outs = expr(i)
            call skpblb(expr(i),len(expr(i)),ii)
            call fexit(-1,1,'MAPDAT: failed to parse expr: '
     .        //outs(1:j+1)//'<-->'//outs(j+2:ii+1),0)
          endif
        enddo
        if (logi) then
          jr = jr+1
          iwk(jr) = ir
        endif
      enddo
C     If some points have been excluded, shrink dat
      if (jr /= nr) then
        call dcopy(nr*nexpr,dat,1,dat2,1)
        call rowpmt(nr,jr,nexpr,iwk,dat2,dat)
        nr = jr
      endif
      call clrsyv(iv0)

      end

      subroutine mapdatx(nexpr,expr,sincl,iwk,nr,nc,z,dat,dat2)
C- Replace columns of dat with algebraic expressions of them
C  This is a variation of mapdat
C  where z is an additional array separate from data
C  The first column is loaded as 'z', the second as 'is'
      implicit none
C     Passed parameters
      integer nexpr,nr,nc,iwk(nr)
      character*(*) expr(nexpr), sincl
C     dat2 is not used unless points are excluded
      double precision dat(nr,max(nexpr,nc)),dat2(nr,nexpr),z(nr,2)
C     Local variables
      integer ir,iv0,ival,i,j,ii,jr
      character*4 xn, outs*80
      logical a2bin,logi

      call numsyv(iv0)
      jr = 0
      do  ir = 1, nr
        call clrsyv(iv0)
        call addsyv('i',dble(ir),ival)

C   ... Load data table
        do  j = 1, nc
          ii = 1
          xn = 'x   '
          call bin2a('',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ir,j),ival)
        enddo
C       Extra variables associated with z
        call addsyv('z',z(ir,1),ival)
        call addsyv('is',z(ir,2),ival)

C   ... Exclude points if not satisfy sincl
        logi = .true.
        if (sincl /= ' ') then
          j = 0
          if (.not. a2bin(sincl,logi,0,0,' ',j,-1))
     .      call rx('mapdat:  error parsing sincl')
        endif

C   ... Put expressions of these vars into dat
        do  i = 1, nexpr
          j = 0
          if (.not. a2bin(expr(i),dat(ir,i),4,0,' ',j,-1)) then
            outs = expr(i)
            call skpblb(expr(i),len(expr(i)),ii)
            call fexit(-1,1,'MAPDAT: failed to parse expr: '
     .        //outs(1:j+1)//'<-->'//outs(j+2:ii+1),0)
          endif
        enddo
        if (logi) then
          jr = jr+1
          iwk(jr) = ir
        endif
      enddo
C     If some points have been excluded, shrink dat
      if (jr /= nr) then
        call dcopy(nr*nexpr,dat,1,dat2,1)
        call rowpmt(nr,jr,nexpr,iwk,dat2,dat)
        nr = jr
      endif
      call clrsyv(iv0)

      end

      subroutine msort(nr,ncs,ncw,iprm,wk,s,ssort)
C- Sort matrix according to increasing wk.
      implicit none
      integer nr,ncs,ncw,iprm(nr)
      double precision wk(ncw,nr),s(nr,ncs),ssort(nr,ncs)
      integer ir,ic,ix

      call dvheap(ncw,nr,wk,iprm,0d0,1)
      do  ir = 1, nr
        ix = iprm(ir)
        do  ic = 1, ncs
          ssort(ir,ic) = s(ix,ic)
        enddo
      enddo
      end

      subroutine rowpmt(nf,nt,nc,ipmt,afrom,ato)
      implicit none
      integer nf,nt,nc,ipmt(nt),i
      double precision afrom(nf,nc), ato(nt,nc)

      do  i = 1, nt
        call dcopy(nc,afrom(ipmt(i),1),nf,ato(i,1),nt)
      enddo
      end

      subroutine psup4(iax,nttab,nlst1,lst1,rtab,disp,nn,dpos)
C- Kernel called by supcel to displace pairs radially
C ----------------------------------------------------------------------
Ci Inputs
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   nlst1 :number of sites of "center" type
Ci   lst1  :list of sites of "center" type
Ci   nlst2 :number of sites of "neighbor" type
Ci   lst2  :list of sites of "neighbor" type
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Co Inputs/Outputs
Cio   pos  :On input, basis vectors
Cio        :On output, displaced basis vectors
Cu Updates
Cu   07 Apr 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlst1,lst1(nlst1),niax,nttab,nn
      double precision rtab(3,*),dpos(3,*),disp(nn,6)
      parameter (niax=10)
      integer iax(niax,*)
C ... Local parameters
      integer iat,low1,i1,i2,ipr,count,stdo,lgunit,i,j
      double precision rhat(3),tol

      count = 0
      low1 = 0
      call getpr(ipr)
      ipr = 100
      stdo = lgunit(1)
      tol = 1d-6

C      if (ipr >= 50) write(stdo,330) nlst1
C  330 format(' supcel: add displacement around',i4,' site(s)'/
C     .  '  ia  ib',15x,'rb-ra',27x,'delta pos(ib)')
      if (ipr >= 50) write (stdo,1)
    1 format('  ia  ib',15x,'rb-ra',27x,'delta pos(ib)')

      do  iat = 1, nttab

C   ... If the center atom isn't in the list, skip this pair
        i1 = iax(1,iat)
        call hunti(lst1,nlst1,i1,0,low1)
        if (low1 >= nlst1) cycle
        if (i1 /= lst1(low1+1)) cycle

C   ... See if this pair is in any of in the disp list
        call dcopy(3,rtab(1,iat),1,rhat,1)
        do  i = 1, nn
          if (abs(rhat(1)-disp(i,1)) < tol) then
            if (abs(rhat(2)-disp(i,2)) < tol) then
C             We have a match!  Add to dpos
              if (abs(rhat(3)-disp(i,3)) < tol) then
                i2 = iax(2,iat)
                dpos(1,i2) = dpos(1,i2) + disp(i,4)
                dpos(2,i2) = dpos(2,i2) + disp(i,5)
                dpos(3,i2) = dpos(3,i2) + disp(i,6)
                if (ipr >= 50) write(stdo,333)
     .            i1,i2,(rtab(j,iat),j=1,3),(disp(i,j),j=4,6)
  333           format(2i4,3f11.6,2x,3f11.6)
                count = count+1
              endif
            endif
          endif
        enddo

      enddo
      call info(30,0,0,' ... Added displacements for %i site(s),'//
     .  ' %i pair(s)',nlst1,count)
      end
      subroutine psup5(nn,disp,rmax)
C- Determine largest length of position vectors
      implicit none
      integer nn
      double precision disp(nn,3),rmax,ddot
      integer i

      rmax = 0
      do  i = 1, nn
        rmax = max(rmax,ddot(3,disp(i,1),nn,disp(i,1),nn))
      enddo
      rmax = dsqrt(rmax)
      end
      subroutine psup6(qlx,plat,nm,nqmx,nq,q)
C- Finds q vectors that are G vectors in a supercell
      implicit none
      integer nq,nm,nqmx
      double precision qlx(3,3),plat(3,3),q(3,nqmx)
      integer j1,j2,j3,i,iq
      double precision G(3),qlat(3,3)
      real(8), parameter :: tol = 1d-6
      procedure(logical) :: latvec

C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
      call mkqlat(plat,qlat,G)  ! G is a dummy argument

C --- Loop over range of G vectors ---
      nq = 0    ! Running count of number of inequivalent q-points
      do  j1 = -nm, nm
      do  j2 = -nm, nm
      do  j3 = -nm, nm

C       Make a G vector
        do  i = 1, 3
          G(i) = qlx(i,1)*j1+qlx(i,2)*j2+qlx(i,3)*j3
        enddo
        call shorps(1,qlat,[1,1,1],G,G)

        do  iq = 1, nq
          if (latvec(1,tol,plat,G(1:3)-q(1:3,iq))) goto 10  ! Already in list
        enddo
        nq = nq+1
        if (nq > nqmx) call rx('bug in psup6')
        q(1:3,nq) = G(1:3)
   10   continue

C       if (nq == nqmx) return

      enddo
      enddo
      enddo

      end
