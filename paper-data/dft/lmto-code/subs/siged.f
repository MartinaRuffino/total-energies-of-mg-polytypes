      subroutine siged(mode,sopts,s_ctrl,s_spec,s_lat,s_ham,s_gw,nsp,
     .  ndhrs,ndimh,iprmb,hreal,ntab,iax,hrs,lrsw)
C- Sigma editor
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:spid
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  idxdn lmxb lmxa ncomp
Co     Stored:     idxdn
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  makidx nscpa
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat pos nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:symgr
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  rsrnge offH rsstol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  eseavr nkabc lshft
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :compound of these digits:
Ci         :1s digit 1 : merge up and down spins, set nsp=1, and exit
Ci         :1s digit 2 : average up and down spins, and exit
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ndhrs :dimensions hrs
Ci   ndimh :hamiltonian dimension
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   hreal :1 if hrs is real
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Cio Inputs/Outputs
Cio  hrs   :real-space hamiltonian, which may be modified on exit
Cl Local variables
Cl   plat  :primitive lattice vectors, in units of alat
Cl   pos   :basis vectors
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Oct 18 Added spin flip or zers
Cu   14 Jul 17 Added spin flip
Cu   09 Dec 15 New mode 2
Cu   14 Oct 13 showp and zero have new conditions for pairs
Cu   20 Sep 13 showp redesigned
Cu   29 Aug 13 major redesign
Cu   01 May 06 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
      logical lrsw
      integer mode,nsp,ndhrs,ndimh,iprmb(*),hreal
      integer n0,nkap0,n0H,niax
      parameter (n0=10,nkap0=4,n0H=5,niax=10)
      integer iax(niax,*),ntab(*)
      double precision hrs(2,ndhrs,ndhrs,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_gw)::    s_gw
C ... Local Dynamical arrays
      integer, allocatable :: iax0(:)
      logical, allocatable :: lzero(:)
      integer, allocatable :: imap(:),ipsx(:),iprmbx(:),ilst(:)
      integer, allocatable :: ntabx0(:),ntabx(:),iaxx(:,:)
      real(8), allocatable :: rtab(:,:),hrsx(:,:,:,:,:),qpgw(:),delT(:,:),wk(:),rtabx0(:,:)
      integer, pointer:: offH(:,:),offHx(:,:)
      real(8),pointer :: pos(:,:),posx(:,:)
C     type(str_site), allocatable::  s_sitex(:)
C     Allows calls to mkilssa, which allocates and assigns list
      interface
         subroutine mkilssa(iopt,slst,nlist,list,irnge)
         integer iopt,nlist,irnge(2)
         integer, allocatable, target :: list(:)
         character*(*) slst
        end subroutine mkilssa
      end interface
C ... Fixed parameters
      integer(4),parameter:: NULLI =-99999
C ... Local parameters
      logical lin,lsopts,ltrans,lhavex,ltmp,lfbz,lshorten
      character outs*150, strn*80, dc*1, dc2*1, fn*120
      integer i,i1,ib,ib1,ib2,ibx,ifi,ifis2,ik1,
     .  ioibmap,iorb,ip,jp,ipx,jpx,isp,isub(4),iq,it,it1,it2,j,j1,
     .  j2,jb,jbx,js1,js2,kcplx,l1,lascii,lrot,mxcsiz,mxkp,nbas,nbx,
     .  ndimhx,ndimx,nkaph,nl,nlmto2,norb,nqp,nspec,nttab,nttabx0,nttabx,
     .  nw,stdo,nsgrp,ibx1,ipr,ifis,nspx,ndhrsx,n,k,hrealx,nscnts,isw,
     .  nsites,lprint,nsp2,nspc2,ndimh2,nqp2
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),blks(n0*nkap0),ntabh(n0*nkap0)
      integer ix(100),iwk(100),nkxyz(3),nkxyz2(3),lshft(3),lshft2(3),ldhamx(16)
      double precision alat,alatx,rr,drr2,rcut,vol,volx,rsrangx
      double precision xx(10),plat(3,3),plx(3,3),vlat(3),rtabx(3),
     .  qlx(3,3),qb(3,3),eseavr(2),eseavr2(2),rsstol,rsrnge
      character*1 :: star(0:1)
      procedure(logical) :: isanrg
      procedure(integer) :: nglob,wordsw,mkilsd,a2vec,fopna,iosite
      procedure(real(8)) :: dlength,ddot,dglob
      procedure(real(8)) :: dsum

C ... Default range: all pairs, ham subblock is (1:ndhrs,1:ndhrs)
      nbas = nglob('nbas')
      nspec = nglob('nspec')
      nl = nglob('nl')
      nkaph = nglob('nkaph')
      it1 = 1
      it2 = ntab(nbas+1)
      stdo = nglob('stdo')
      rcut = 2
      call getpr(ipr)
      if (hreal /= 0) call rx('siged not ready for hreal')
      isub(1) = 1
      isub(2) = ndhrs
      isub(3) = 1
      isub(4) = ndhrs
      alat = s_lat%alat
      plat = s_lat%plat
      pos => s_lat%pos
      nttab = ntab(nbas+1)
      rsrnge = s_ham%rsrnge
      rsrangx= s_ham%rsrnge
      offH => s_ham%offH
      eseavr = s_gw%eseavr
      ltrans = .true.
      lhavex = .false.
      lascii = 0
      kcplx = 1
      lrot = 0
      rsstol = s_ham%rsstol
      lfbz = .false.
      star(0) = '*'
      star(1) = ' '

C ... Merge up, down spin potentials into average
C     print *, dsum(4*ndhrs*ndhrs*it2,hrs,1)
      if (mod(mode,10) == 1 .or. mod(mode,10) == 2) then
        do  i = it1, it2
C         print *, i
C         call zprm('h1,i',2,hrs(1,1,1,1,i),ndhrs,ndhrs,ndhrs)
C         call zprm('h2,i',2,hrs(1,1,1,2,i),ndhrs,ndhrs,ndhrs)
          call daxpy(2*ndhrs*ndhrs,1d0,hrs(1,1,1,2,i),1,hrs(1,1,1,1,i),1)
          call dcopy(2*ndhrs*ndhrs,hrs(1,1,1,1,i),1,hrs(1,1,1,2,i),1)
C         call zprm('hi',2,hrs(1,1,1,1,i),ndhrs,ndhrs,ndhrs)
        enddo
C       print *, dsum(4*ndhrs*ndhrs*it2,hrs,1)
        call dscal(4*ndhrs*ndhrs*it2,0.5d0,hrs,1)
        s_gw%eseavr(1) = (eseavr(1) + eseavr(2))/2
        s_gw%eseavr(2) = (eseavr(1) + eseavr(2))/2

C       print *, dsum(4*ndhrs*ndhrs*it2,hrs,1)
        if (mod(mode,10) == 2) return
        nsp = 1
        rr = dglob('nsp',1d0,1)
        return
      endif

C      print *, '!!'
C      call dvset(hrs,1,2*ndhrs*ndhrs*nsp*it2,-1d0)

      dc = sopts(1:1)
      if (sopts(2:5) == 'edit') then
        outs = sopts(6:)
        sopts = outs
        dc = sopts(1:1)
      endif
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the static self-energy file editor. ',
     .    'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the static self-energy file editor.  ',
     .    'Enter ''?'' to see options.')
        lsopts = .false.
      endif

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

      write(stdo,"(/' Option : ')",advance='no')
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif
      call locase(outs)

C ... Parse and execute the next command
      j1 = 1
      if (outs(1:4) == 'edit') then
        goto 10

C ... Flag outfile to be R.S.
      elseif (outs == 'rs') then
        lrsw = .true.
        goto 10

C --- Null input ---
      elseif (outs == ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',' ''?'' to see menu')
        goto 10

C --- Specify full BZ ---
      elseif (outs(1:4) == 'fbz') then
        lfbz = .true.
        goto 10

C --- Specify pair ---
C      elseif (outs(1:4) == 'pair') then
C        call words(outs,nw)
C        if (nw == 1) then
C          it1 = 1
C          it2 = ntab(nbas+1)
C        else
C          i = 5
C          nw = a2vec(outs,len(outs),i,2,', ',2,-3,2,ix,iwk)
C          if (nw <= 0) then
C            goto 99
C          else if (nw == 1) then
C            iwk(2) = iwk(1)
C          endif
C          it1 = min(iwk(1),iwk(2))
C          it2 = max(iwk(1),iwk(2))
C        endif
C
C        ib = iax(1,it1)
C        jb = iax(2,it1)
C
C        if (it1 == it2) then
C          print '('' edit pair'',i6,''   ib ='',i4,''   jb ='',i4)',
C     .      it1,ib,jb
C
C          ib = iax(1,it1)
C          isub(1) = 1
C          call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(2))
C          jb = iax(2,it1)
C          isub(3) = 1
C          call orbl(jb,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(4))
C        else
C          print '('' edit range'',i5,''  to'',i5)', it1,it2
C        endif
C        goto 10

C --- Summarize hamiltonian or pair table ---
      elseif (outs(1:3) == 'sum' .or. outs(1:4) == 'sumx') then

        ltmp = outs(1:4) == 'sumx'
        if (ltmp .and. .not. lhavex) then
          call info0(2,0,0,' read auxilliary hrsx before using sumx')
          goto 98
        endif

        k = 2                   ! next argument in outs
C       Re-entry point for optional arguments
   15   continue
        call word(outs,k,j1,j2)
        if (j2 < j1) then ! No more words

C   ... Mode "sum ib=..."
        elseif (outs(j1:j1+1) == 'ib') then
          if (ltmp) then
            call info0(2,1,1,' ib option not available with sumx')
            goto 98
          endif
          ib1 = 1; ib2 = nbas
          if (outs(j1:j1+2) == 'ib=') then
          j = 3
          j = a2vec(outs(j1:),len(outs(j1:)),j,2,', ',2,-2,2,ix,iwk)
          if (j <= 0) goto 99
          if (j == 1) iwk(2) = iwk(1)
          ib1 = min(iwk(1),iwk(2))
          ib2 = max(iwk(1),iwk(2))
          endif

          do  ib = ib1, ib2
          call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,xx)
          call gtbsl1(4,norb,ltab,ktab,xx,xx,ntabh,blks)
          write(stdo,359) ib
  359     format(/
     .      ' Summary of hamiltonian subblocks for site',i4/
     .      '     kap  lstart start   end   size')
          do  iorb = 1, norb
          if (blks(iorb) /= 0) then
            l1  = ltab(iorb) ! l index for this subblock of sigma
            ik1 = ktab(iorb) ! kappa index
            i1 = offl(iorb) ! orbital index in iprmb order
            outs = ' '
C            call awrit2('%i,%i',outs,len(outs),0,1+i1-offl(1),
C     .                  i1-offl(1)+blks(iorb))
C            print '(3i7,5x,a10,i3)', ik1,l1,1+i1, outs(1:10), blks(iorb)
            print '(5i7)', ik1,l1,1+i1, i1+blks(iorb), blks(iorb)

          endif
          enddo
          enddo
          goto 10

        else
          j = 0
          j = a2vec(outs(j1:),len(outs(j1:)),j,2,', ',2,-2,2,ix,iwk)
          if (j <= 0) goto 99
          if (j == 1) iwk(2) = iwk(1)
          it1 = min(iwk(1),iwk(2))
          it2 = max(iwk(1),iwk(2))
          if ((it1 < 1.or.it2 > nttab) .and. .not. ltmp) then
           call info2(2,1,0,' pair range (%i,%i) out of bounds',it1,it2)
           goto 98
          endif

          print 367
  367     format(5x,'pair   ib    jb       iax',9x,'dist  norbi norbj')
          do  it = it1, it2
            ib = iax(1,it)
            jb = iax(2,it)
            if (ltmp) then
              rr = dsqrt(drr2(plx,posx(1,ib),posx(1,jb),
     .            iaxx(3,it),iaxx(4,it),iaxx(5,it),vlat))
              print "(i8, 3i6, 2i3, f12.5)", it,iax(1:5,it),rr
            else
              rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .          iax(3,it),iax(4,it),iax(5,it),vlat))
              call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,i)
              call orbl(jb,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,j)
              print "(i8, 3i6, 2i3, f12.5, 2i5)", it,iax(1:5,it),rr,i,j
            endif
          enddo
          goto 10

        endif

C   ... Default mode (no arguments to sum)
        call info2(2,1,0,' Summary of pair table for'//
     .    '%?#n# auxilliary## R.S. sigma%N'//
     .    ' site  start   end  size',isw(ltmp),0)
        if (ltmp) then
          do  ib = 1, nbx
            write(stdo,357)
     .        ib,ntabx(ib)+1,ntabx(ib+1),ntabx(ib+1)-ntabx(ib)
          enddo
        else
          do  ib = 1, nbas
            write(stdo,357) ib,ntab(ib)+1,ntab(ib+1),ntab(ib+1)-ntab(ib)
          enddo
        endif
  357   format(i4,2i7,i5)

        goto 10

C --- Site-local sigma ---
        write(stdo,358)
  358   format(/' Zeroing out off-diagonal parts of sigma')

        do  i = it1, it2
          ib = iax(1,i)
          if (ntab(ib)+1 == i) then
            call info2(10,0,0,' retaining pair %i, ib=%i',i,ib)
          else
            do  isp = 1, nsp
              call dpzero(hrs(1,1,1,isp,i),2*ndhrs*ndhrs)
            enddo
          endif
        enddo

        goto 10

C --- Summarize hamiltonian or pair table ---
      elseif (outs(1:4) == 'spin') then

        dc2 = outs(5:5)
        j = 1 ! marks tag as 'flip', if present
        ltmp = wordsw(outs,dc2,'flip',' ='//dc2,j1) > 0
        if (.not. ltmp) then
          j = 2 ! marks tag as 'zers', if present
          ltmp = wordsw(outs,dc2,'zers',' ='//dc2,j1) > 0
        endif

        if (ltmp) then
          if (nsp == 1) then
            call info0(1,0,0,'only one spin ... nothing to do')
            goto 10
          endif

C     ... Get range, if supplied
          if (outs(j1:j1) == '=') then ! case flag is 'flip=list'
            call nwordg(outs,0,dc2//' ',1,j1,j2); j1 = j1+1
            call mkilssa(11,outs(j1:j2),nsites,ilst,[1,ntab(1+nbas)])
          else
            nsites = ntab(1+nbas); allocate(ilst(nsites))
            forall (i = 1:nsites) ilst(i) = i
          endif

          call dpzero(xx,4)
          call info2(1,0,0,' %?#n==1#swapping#zeroing# spins, %i pairs ...',j,nsites)
          do  it = 1, nsites
            i = ilst(it)
C            xx(1) = xx(1) + dsum(2*ndhrs*ndhrs,hrs(1,1,1,1,i),1)
C            xx(2) = xx(2) + dsum(2*ndhrs*ndhrs,hrs(1,1,1,2,i),1)
            if (j == 1) then
              call dswap(2*ndhrs*ndhrs,hrs(1,1,1,1,i),1,hrs(1,1,1,2,i),1)
            else
              call daxpy(2*ndhrs*ndhrs,1d0,hrs(1,1,1,2,i),1,hrs(1,1,1,1,i),1)
              call dscal(2*ndhrs*ndhrs,0.5d0,hrs(1,1,1,1,i),1)
              call dcopy(2*ndhrs*ndhrs,hrs(1,1,1,1,i),1,hrs(1,1,1,2,i),1)
            endif
C            xx(3) = xx(3) + dsum(2*ndhrs*ndhrs,hrs(1,1,1,1,i),1)
C            xx(4) = xx(4) + dsum(2*ndhrs*ndhrs,hrs(1,1,1,2,i),1)
          enddo
C         call info2(30,0,0,' sum check hrs before flip %2:1;9,9F  after flip %2:1;9,9F',xx,xx(3))
C         print *, dsum(4*ndhrs*ndhrs*ntab(1+nbas),hrs,1)

        endif

        if (j == 1) then
          call dswap(1,eseavr(1),1,eseavr(2),1)
          s_gw%eseavr = eseavr
        else if (nsites == ntab(1+nbas)) then ! average if all sites were averaged
          s_gw%eseavr(1) = (eseavr(1) + eseavr(2))/2
          s_gw%eseavr(2) = (eseavr(1) + eseavr(2))/2
          eseavr = s_gw%eseavr
        endif
        deallocate(ilst)

        goto 10

C --- Change average sigma for high energy ---
      elseif (outs(1:6) == 'sigbar') then

        k = 2                   ! next argument in outs
        call word(outs,k,j1,j2)
        if (j2 <= j1) goto 98
        j = 0
        j = a2vec(outs(j1:),len(outs(j1:)),j,4,', ',2,-2,2,ix,xx)
        if (j /= 1 .and. j /= 2) goto 98
        eseavr(1:j) = xx(1:j)
        s_gw%eseavr = eseavr

C --- Dump hrs to file ---
      elseif (outs(1:6) == 'rssave') then

        k = 2                   ! next argument in outs
        call word(outs,k,j1,j2)
        fn = 'rssigm'
        if (j2 >= j1) fn = outs(j1:j2)
        ifis = fopna(trim(fn),-1,4)
        nkxyz = s_gw%nkabc
        lshft = s_gw%lshft
C       nqp is for full BZ but shouldn't matter for R.S. sigm
        nqp = nkxyz(1)*nkxyz(2)*nkxyz(3)
        rewind ifis
        lascii = 0  ! Only binary files for now
        call iosigh(lascii*10+0,100,nsp,1,ndimh,ndimh,nkxyz(1),nkxyz(2),
     .    nkxyz(3),nqp,iq,lshft(1),lshft(2),lshft(3),-ifis,eseavr)
        call iosigrs(lascii*10+0,100,1d-6,hreal,nsp,nbas,ndhrs,nttab,
     .    s_lat%alat,plat,s_lat%pos,s_ham%nprs,s_ham%iaxs,
     .    s_ham%hrs,-ifis)
        call info2(10,0,0,' wrote hrs to file '//trim(fn)//
     .    ', %i sites, %i pairs',nbas,nttab)

        goto 10

C --- Read hrsx from file ---
      elseif (outs(1:6) == 'readx') then

        k = 2                   ! next argument in outs
        call word(outs,k,j1,j2)
        fn = 'rssigm'
        if (j2 >= j1) fn = outs(j1:j2)
        ifis = fopna(trim(fn),-1,4)
        rewind ifis
        lascii = 0  ! Only binary files for now
C       Read in all dimensioning parameters
        call iosigh(lascii*10+0,i,nsp2,nspc2,ndimh2,nlmto2,
     .    nkxyz2(1),nkxyz2(2),nkxyz2(3),nqp2,iq,
     .    lshft2(1),lshft2(2),lshft2(3),ifis,eseavr2)
C       Re-read header, requiring match for  mode, nsp
        rewind ifis
        call iosigh(lascii*10+3,100,nsp,nspc2,ndimh2,nlmto2,
     .    nkxyz2(1),nkxyz2(2),nkxyz2(3),nqp2,iq,
     .    lshft2(1),lshft2(2),lshft2(3),ifis,eseavr2)

C       Read mode,hrealx,nttabx,nspx,nbx,ndhrsx
        call iosigrs(lascii*10+3,i,1d-6,hrealx,nspx,nbx,ndhrsx,nttabx,
     .    xx,xx,xx,xx,xx,xx,ifis)
C        if (.not. isanrg(i,100,100,'iosigrs:',"file's mode",.false.)
C     .    goto 98 ! should never happen, since 1st header has this already
        ltmp = isanrg(nspx,nsp,nsp,'iosigrs:',"file's nsp",.false.).or.
     .  isanrg(hrealx,hreal,hreal,'iosigrs:',"file's hreal",.false.).or.
     .  isanrg(ndhrsx,ndhrs,ndhrs,'iosigrs:',"file's ndhrs",.false.)
        if (ltmp) goto 98
        if (allocated(ntabx)) then
          deallocate(ntabx,iaxx,posx,hrsx)
        endif
        allocate(ntabx(nbx+1),iaxx(niax,nttabx),posx(3,nbx))
        allocate(hrsx(2,ndhrs,ndhrs,nspx,nttabx))

C       Read alatx,plx,posx
        call iosigrs(lascii*10+4,i,1d-6,hrealx,nspx,nbx,ndhrsx,nttabx,
     .    alatx,plx,posx,ntabx,iaxx,hrsx,ifis)
C       Read ntabx,iaxx,hrsx
        call iosigrs(lascii*10+8,i,1d-6,hrealx,nspx,nbx,ndhrsx,nttabx,
     .    alatx,plx,posx,ntabx,iaxx,hrsx,ifis)

        call fclose(ifis)
        call info5(10,0,0,' read hrsx from file '//trim(fn)//', '//
     .    '%i sites, %i pairs, ham dimension %i',nbx,nttabx,ndimh2,0,0)
        lhavex = .true.

        goto 10

C --- Insert portions of hrsx into hrs ---
      elseif (outs(1:5) == 'merge') then

        if (.not. lhavex) then
          print *, 'siged: readx required before merge'
          goto 98
        endif

        k = 2                   ! next argument in outs
        call word(outs,k,j1,j2)
        if (j2 < j1) goto 98
C       Switch -sigbar
        if (index(outs(j1:j2),'-sigbar') == 1) then ! Copy
          eseavr = eseavr2
          s_gw%eseavr = eseavr
          k = k+1
          call word(outs,k,j1,j2)
        endif
        nsites = mkilsd(outs(j1:j2),-1,j)
        if (nsites <= 0) then
          call info0(0,0,0,'%6p... Bad or null list : '//outs(j1:j2))
          goto 98
        endif
        allocate(ilst(nsites))
        if (mkilsd(outs(j1:j2),nsites,ilst) /= nsites)
     .    call rx('bug in mkilsd')
        if (nsites == 1 .and. ilst(1) == 0) then
          deallocate(ilst)
          nsites = nbas
          allocate(ilst(nbas))
          do  i = 1, nbas
            ilst(i) = i
          enddo
        endif
        call imxmn(nsites,ilst,1,i,j)
        if (i < 0 .or. j > nbas) goto 97

C        if (nttabx /= nttab) then
C          write (stdo,"(' not ready for incommensurate rssig')")
C          goto 99
C        endif

        call info0(10,0,0,' Copying the following pairs from hrsx'//
     .    ' and their conjugates'//
     .    '%N    site  start   end  number matched')

        nw = ndhrs*ndhrs*nsp
        do  i = 1, nsites
          ib = ilst(i)
C         Loop over all pairs connected to point equivalent to ib
          j1 = ntab(ib)
          n = 0
          do  ip = ntab(ib)+1, ntab(ib+1)
C           Find equivalent pair in hrsx table
            ibx = ib
C           lok = 0
            do  ipx = ntabx(ibx)+1, ntabx(ibx+1)
              if (iaxx(1,ipx) /= iax(1,ip) .or.
     .            iaxx(2,ipx) /= iax(2,ip) .or.
     .            iaxx(3,ipx) /= iax(3,ip) .or.
     .            iaxx(4,ipx) /= iax(4,ip) .or.
     .            iaxx(5,ipx) /= iax(5,ip)) cycle
              n = n+1
              call dcopy(2*nw,hrsx(1,1,1,1,ipx),1,hrs(1,1,1,1,ip),1)
              jp = iax(6,ip)
              jpx = iaxx(6,ipx)
              call dcopy(2*nw,hrsx(1,1,1,1,jpx),1,hrs(1,1,1,1,jp),1)
C             lok = 1
              exit
            enddo
C            if (lok == 0) then
C              print *, 'hi'
C            endif
          enddo
          write(stdo,355) ib,ntab(ib)+1,ntab(ib+1),ntab(ib+1)-ntab(ib),
     .      n,star(isw(ntab(ib+1)-ntab(ib) == n))
  355     format(5i7,1x,a1)
        enddo

        deallocate(ilst)
        goto 10

C --- Zero hrs_RR' outside given range ---
      elseif (outs(1:5) == 'zero ') then

        nw = ndhrs*ndhrs*nsp
        k = 2                   ! next argument in outs
        xx(1) = 0
        call word(outs,k,j1,j2)
        if (j2 <= j1) goto 98
C       Switch shorten
        if (index(outs(j1:j2),'-short') == 1) then ! Shorten or short
          lshorten = .true.
          k = k+1
          call word(outs,k,j1,j2)
        endif

C   ... Translate conditions into pair list
        allocate(lzero(nttab))
        call sig_pair_list(outs,k,nbas,plat,pos,nttab,ntab,iax,lzero,lprint,j1)
        if (j1 == 0) goto 10  ! illegal radius; nothing to do
        if (j1 < 0) then
          j1 = -j1; goto 97
        endif

C   ... Zero pairs satisfying conditions
        do  ib = 1, nbas
          do  ip = ntab(ib)+1, ntab(ib+1)
            if (.not. lzero(ip)) cycle
            call dpzero(hrs(1,1,1,1,ip),2*nw)
          enddo
        enddo

C   ... Collapse ntab,iax,hrs
        if (lshorten) then
          call pshorten(nbas,nttab,ndhrs,nsp,lzero,ntab,iax,hrs)
C         Debugging check: remake iax(6) through pairc call.
C         Should not be necessary
C         call pshpr(41)
C         call ppair1(40,1,nbas,nbas,-1,alat,plat,pos,xx,
C     .     nttab,ntab,iax,1000)
C         call ppair1(20,1,nbas,nbas,-1,alat,plat,pos,xx,
C     .     nttab,ntab,iax,1000)
C         call poppr
C         stop
        endif

        deallocate(lzero)
        goto 10

C --- Supercell of hrs ---
      elseif (outs(1:5) == 'scell') then

        lrsw = .true.           ! Write R.S. sigma by default
        k = 2                   ! next argument in outs
        nkxyz2 = 0              ! Flag nkabc for scell not read
        eseavr2 = eseavr
        hrealx = hreal
        dc2 = outs(6:6)

        if (wordsw(outs,dc2,'sigq','',j1) > 0) then ! write sigma(k)
          lrsw = .false.
        endif
        if (wordsw(outs,dc2,'sigrs','',j1) > 0) then ! write sigma(T). Imples nkgw=1
          lrsw = .true.
          nkxyz2 = [1,1,1]
        endif
        if (wordsw(outs,dc2,'range','=',j1) > 0) then ! dimensions of original mesh data
          if (a2vec(outs,len_trim(outs),j1,4,' '//dc2,2,1,1,ix,rsrangx) < 1) goto 98
          print *, '(warning) range has no effect yet ...'
        endif
        if (wordsw(outs,dc2,'nq','=',j1) > 0) then ! dimensions of original mesh data
          if (a2vec(outs,len_trim(outs),j1,2,', '//dc2,3,3,3,ix,nkxyz2) /= 3) goto 98
        endif
        call word(outs,2,j1,j2)
        if (j2 <= j1) goto 98

C       Read alat,plat,nbas of supercell
        j = iosite(7000,3d0,0,outs(j1:j2),ifi,s_ctrl%spid,alatx,plx,nbx,
     .    nspec,xx,xx,xx,xx,xx,xx,xx)
        call info5(10,0,0,' scell has %i sites with lattice parameter'//
     .    ' %d = %d*alat%N plx %9:1;6,6d%N plat%9:1;6,6d',
     .    nbx,alatx,alatx/alat,plx,plat)
        call dinv33(plx,1,qlx,volx)  ! Make qlx
        call dinv33(plat,0,xx,vol)
        call info2(10,0,0,' Volume ratio %d',abs(volx/vol),0)
        if (abs(dble(nbx)/nbas-abs(volx/vol)) > 1d-5) then
          write (stdo,13)
   13     format(' supercell has incommensurate lattice ... skip scell')
          goto 99
        endif
C       Read pos,ips,imap of supercell
        allocate(posx(3,nbx),ipsx(nbx))
        j = iosite(8000,3d0,0,outs(j1:j2),ifi,s_ctrl%spid,alatx,plx,nbx,
     .    nspec,posx,xx,xx,xx,ipsx,xx,xx)
        allocate(imap(nbx)); call iinit(imap,nbx)
        if (ioibmap(ifi,nbx,imap) /= nbx) then
          call info0(1,0,0,' siged could not read map data from site file ... nothing done')
          goto 10
        endif

C   ... Raw pair table for rs hamiltonian
        call info(30,0,0,' siged: make neighbor table for SL r.s.'//
     .    'hamiltonian using range = %;4d * alat',s_ham%rsrnge,0)
C       allocate(ntabx0(nbx+1),offHx(n0H*nkap0,(nbx+1)))
        allocate(ntabx0(nbx+1))
        call pshpr(min(ipr-10,30))
        mxcsiz = 0
        call pairs(nbx,nbx,1d0,plx,[s_ham%rsrnge/2],posx,
     .    [-1],3,-1,[0],nttabx0,ntabx0,iax0,mxcsiz)
        call poppr

C   ... Make ntabx for supercell
        allocate(ntabx(nbx+1),offHx(n0H*nkap0,(nbx+1)))
        ntabx(1) = 0
        offHx(1,1) = 0
        do  ibx = 1, nbx
          ib = imap(ibx)
          ntabx(ibx+1) = ntabx(ibx) + ntab(ib+1) - ntab(ib)
        enddo
        nttabx = ntabx(nbx+1)

C   ... Generate SL pair table by matching corresponding rtab
        if (allocated(rtab)) deallocate(rtab)
        allocate(rtab(3,nttab),rtabx0(3,nttabx0),iaxx(niax,nttabx))
        call mkrtab(000,1d0,plat,pos,iax,nttab,pos,rtab)
        call mkrtab(000,1d0,plx,posx,iax0,nttabx0,posx,rtabx0)
        do  ibx = 1, nbx
          ib = imap(ibx)
          ipx = ntabx(ibx)
C         Loop over all pairs connected to point equivalent to ibx
          do  ip = ntab(ib)+1, ntab(ib+1)
            ipx = ipx+1
            iaxx(1,ipx) = ibx
C           Find equivalent pair in SL table
            do  jp = ntabx0(ibx)+1, ntabx0(ibx+1)
              rtabx(:) = rtabx0(:,jp) - rtab(:,ip)
              if (dlength(3,rtabx,1) > 1d-6) cycle
C             Sanity check
              if (ibx /= iax0(1+niax*(jp-1)))
     .          call rx('bug in siged')
              jbx = iax0(2+niax*(jp-1))
              iaxx(2,ipx) = jbx
              rtabx(:) = rtabx0(:,jp) - (posx(:,jbx)-posx(:,ibx))
              call dgemm('T','N',3,1,3,1d0,qlx,3,rtabx,3,0d0,xx,3)
              iaxx(3,ipx) = nint(xx(1))
              iaxx(4,ipx) = nint(xx(2))
              iaxx(5,ipx) = nint(xx(3))
C             Sanity check
              if (dlength(3,xx(1:3)-iaxx(3:5,ipx),1) > 1d-6)
     .          call rx('bug in siged')
              iaxx(8,ipx) = ip
              goto 25
            enddo
C           Should never happen
            print *, ibx,ip
            call rx('bug in siged')
   25       continue
          enddo
        enddo
        if (ipx /= nttabx) call rx('bug in siged')
        deallocate(rtabx0,iax0,ntabx0)

C   ... Determine translations for phase shifts
        ibx1 = 0
        allocate(delT(3,nbx),wk(3*nbx)); delT = 0
C       allocate(s_sitex(nbx))
C        do  ib = 1, nbas
C          call hunti(imap,nbx,ib,0,ibx1); ibx1 = ibx1+1 !first occurence of ibx
C          if (imap(ibx1) /= ib) call rx('bug in siged')
C          delT(:,ibx1) = posx(:,ibx1) - pos(:,ib)
C          do  ibx = ibx1, nbx
C            if (imap(ibx) /= ib) exit
C            delT(:,ibx) = posx(:,ibx) - pos(:,ib)
CC            s_sitex(ibx)%clabel = ' ' ; i = 0
CC            call bin2a(' ',0,0,ib,2,0,8,s_sitex(ibx)%clabel,i)
C          enddo
C        enddo
CC       Transform delT to plx units for to set up chksgr
C        call dcopy(3*nbx,delT,1,wk,1)
C        call dgemm('T','N',3,nbx,3,1d0,qlx,3,wk,3,0d0,delT,3)
C        deallocate(wk)
CC       lrot = 2
CC       Not needed, but save for record
CC       call iopos(.true.,1,'phase',nbx,delT,s_sitex)
CC       deallocate(s_sitex)

C   ... Fill out colums 6,7,10 of the iaxx table
        call ppair1(0,1,nbx,nbx,-1,alatx,plx,posx,rsrnge,nttabx,ntabx,iaxx,mxcsiz)

C   ... Copy hrs to hrsx
        allocate(hrsx(2,ndhrs,ndhrs,nsp,nttabx))
        nw = ndhrs*ndhrs*nsp
        do  ipx = 1, nttabx
          ip = iaxx(8,ipx)
          call dcopy(2*nw,hrs(1,1,1,1,ip),1,hrsx(1,1,1,1,ipx),1)
C         delT = matmul(plx,iaxx(3:5,ipx))
C         print *, ip,ipx,ip-ipx
          iaxx(8,ipx) = 0  ! Proper table has no mapping
        enddo

C   ... Set up k mesh
C       Use full BZ since symops not known ??
        if (nkxyz2(1) == 0) then
          nkxyz2 = s_gw%nkabc
        endif
        lshft = s_gw%lshft
        mxkp = nkxyz2(1)*nkxyz2(2)*nkxyz2(3)
        allocate(qpgw(3*mxkp))
        nsgrp = s_lat%nsgrp
C       If lattice unchanged, preserve symmetry
C        if (dlength(9,alatx*plx-alat*plat,1) > 1d-6) then
C          nsgrp = 1
C          lfbz = .true.
C        endif
        if (lfbz) nsgrp = 1
        call info(20,1,0,' k-mesh for supercell ...',0,0)
        call bzmshp(' SIGED',1,nkxyz2,lshft,plx,s_lat%symgr,nsgrp,nsgrp,
     .    lfbz,0,0,qb,nqp,qpgw,xx,xx,xx)

C   ... Hamiltonian dimensioning arrays for supercell
        ndimx = nbx * nl**2 * nkaph
        allocate(iprmbx(ndimx+2))
        call iinit(ldhamx,16)
        call makidx(nl,nkaph,1,nbx,0,s_spec,ipsx,offHx,iprmbx,ldhamx)
        ndimhx = offHx(1,nbx+1)

        if (lrsw) then
          fn = 'sigm2rs'
        else
          fn = 'sigm2'
        endif
        ifis2 = fopna(trim(fn),-1,4)
        nscnts = 0
        if (lrsw) nscnts = 100
        call iosigh(0,nscnts,nsp,1,ndimhx,ndimhx,nkxyz2(1),nkxyz2(2),nkxyz2(3),
     .    nqp,iq,lshft(1),lshft(2),lshft(3),-ifis2,eseavr2)
        if (lrsw) then
          call iosigrs(lascii*10+0,nscnts,1d-6,hrealx,nsp,nbx,ndhrs,
     .      nttabx,alatx,plx,posx,ntabx,iaxx,hrsx,-ifis2)
          call rx0('done writing rs sigma, file '//fn)
        else
        call chksgr(1,ltrans,lascii,kcplx,plx,nsp,ndimhx,ifis2,hreal,
     .    iprmbx,nttabx,iaxx,hrsx,ndhrs,rsstol,nqp,qpgw,nbx,lrot,xx,
     .    delT)
          call rx0('done writing sigma(k), file '//fn)
        endif

C --- Get subblock of hrs ---
C      elseif (outs(1:3) == 'sub') then
C
CC   ... Subblock is entire hrs in the absence of explicit arguments
C        ib = iax(1,it1)
C        isub(1) = 1
C        call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(2))
C        jb = iax(2,it1)
C        isub(3) = 1
C        call orbl(jb,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(4))
C        call words(outs,nw)
C        if (nw > 1) then
C          call word(outs,2,j1,j2)
C          i = j1-1
C          nw = a2vec(outs,len(outs),i,2,', ',2,-3,4,ix,iwk)
C          if (nw /= 4) goto 99
C          call icopy(4,iwk,1,isub,1)
C          goto 10
CC        else
CC          j1 = 1
CC          goto 99
C        endif
C
C --- Mask subblock of hrs ---
      elseif (outs(1:5) == 'masks') then

        call rx('masks branch assumes it1 == it2 ... fix')
        call rx(' need isub(1..4)')
        call words(outs,nw)
        if (nw == 1) then
          print *, 'no subblocks specified ... nothing to mask'
          goto 10
        endif

        do  i = isub(1), isub(2)
        do  j = isub(3), isub(4)

C         See whether (i,j) is inside any mask
          lin = .false.
          do  i1 = 2, nw
            call word(outs,i1,j1,j2)
            l1 = j1-1
            j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,4,ix,iwk)
            if (j1 /= 4) goto 99
            lin = lin .or. (i >= iwk(1) .and. i <= iwk(2) .and.
     .                      j >= iwk(3) .and. j <= iwk(4))
          enddo

C         If (i,j) is outside all masks, zero it out
          if (.not. lin) then
            hrs(1,i,j,1,it1) = 0
            hrs(2,i,j,1,it1) = 0
            hrs(1,i,j,nsp,it1) = 0
            hrs(2,i,j,nsp,it1) = 0
          endif


        enddo
        enddo

        goto 10

C --- Mask hrs not in range ---
      elseif (outs(1:5) == 'maskn') then

        call rx('fix maskn branch')

        call words(outs,nw)
        if (nw == 1) then
          print *, 'no ranges specified ... nothing to mask'
          goto 10
        endif

        do  it = it1, it2

C         See whether (i,j) is inside any mask
          lin = .false.
          do  i1 = 2, nw
            call word(outs,i1,j1,j2)
            l1 = j1-1
            j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,2,ix,iwk)
            if (j1 /= 2) goto 99
            lin = lin .or. (it >= iwk(1) .and. it <= iwk(2))
          enddo

C         If (it) is outside all masks, zero it out
          if (.not. lin) then
            call dpzero(hrs(1,1,1,1,it),2*ndhrs*ndhrs*nsp)
          endif

        enddo

        goto 10

C --- Show sigma for range of pairs ---
      elseif (outs(1:5) == 'shows' .or. outs(1:6) == 'showsd') then

        call words(outs,nw)

        do  i = 1, max(nw-1,1)
          if (nw == 1) then
            call icopy(4,isub,1,iwk,1)
          else
            call word(outs,i+1,j1,j2)
            l1 = j1-1
            if (outs(1:6) == 'showsd') then
              j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,2,ix,iwk)
              if (j1 /= 2) goto 99
            else
              j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,4,ix,iwk)
              if (j1 /= 4) goto 99
            endif
          endif

          do  it = it1, it2
          ib = iax(1,it)
          jb = iax(2,it)
          rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .      iax(3,it),iax(4,it),iax(5,it),vlat))
          write(stdo,332)
          write(stdo,334) iax(1,it),iax(2,it),
     .      (vlat(j),j=1,3), rr, (iax(j,it), j=3,5), it

          do  isp = 1, nsp
            if (outs(1:6) == 'showsd') then
              write(stdo,'(19f12.6)')
     .          (hrs(1,j1,j1,isp,it), j1=iwk(1),iwk(2))
              write(stdo,'(19f12.6)')
     .          (hrs(2,j1,j1,isp,it), j1=iwk(1),iwk(2))
            else
            call awrit5('rows (%i,%i) cols (%i,%i) spin %i',strn,80,0,
     .        iwk(1),iwk(2),iwk(3),iwk(4),isp)
            call ywrm(0,strn,3,stdo,'(19f12.6)',
     .        hrs(1,iwk(1),iwk(3),isp,it),1,ndhrs,
     .        iwk(2)-iwk(1)+1,iwk(4)-iwk(3)+1)
            endif
          enddo
          enddo
        enddo

        goto 10

C --- Show pairs within table ---
      elseif (outs(1:5) == 'showp' ) then

C   ... Translate conditions into pair list
        allocate(lzero(nttab))
        call sig_pair_list(outs,2,nbas,plat,pos,nttab,ntab,
     .    iax,lzero,lprint,j1)
        if (j1 == 0) goto 10  ! illegal radius; nothing to do
        if (j1 < 0) then
          j1 = -j1; goto 97
        endif

        i = 0
        nw = ndhrs*ndhrs*nsp
        do  it = it1, it2
          ib = iax(1,it)
          jb = iax(2,it)
          rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .      iax(3,it),iax(4,it),iax(5,it),vlat))
          if (ib /= i) then
            write(stdo,332)
  332       format(/'  ib  jb',9x,'--- r(jb)-r(ib) ---',10x,
     .        'd       -x-plat-   pair     <s>')
C            write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut
C  345       format(' pairc, ib=',i3,':',i4,' neighbors in range',f7.3)
          endif
          i = ib
          if (lzero(it)) then
            if (lprint == 1) then
              n = 2*nw
              xx(1) = dsqrt(ddot(n,hrs(1,1,1,1,it),1,hrs(1,1,1,1,it),1))
            endif
            if (lprint == 0) then
              write(stdo,334) iax(1,it),iax(2,it),
     .          (vlat(j),j=1,3), rr, (iax(j,it), j=3,5), it
            else
              write(stdo,334) iax(1,it),iax(2,it),
     .          (vlat(j),j=1,3), rr, (iax(j,it), j=3,5), it, xx(1)
            endif
  334       format(i4,i4,3f11.6,f9.4,3x,3i3,i6,2f12.7)
          endif
        enddo

        goto 299

C --- abort ---
      elseif (outs(1:1) == 'a') then
        call rx0('aborting sigma editor ... no file written')

C --- quit ---
      elseif (outs(1:1) == 'q') then
        goto 199

C --- help ---
      elseif (outs == '?') then
        print 310
  310   format(/
C     .  ' Editor works with a range of pairs,',
C     .    ' or with a particular pair.'/
C     .    ' In the latter case, you can edit a subblock',
C     .    ' of h(rs) for that pair'//
     .    ' Select one of these options:'//
C    .  t4,'i',t11,'to import a portion of a sigma file'/
C     .  t4,'pair',t11, '[#1 [#2]] specifies a pair or range of pairs',
C     .    ' used by other options.'/
C     .           ,t11, 'Omit #2 to select a specific pair.'/
C     .           ,t11, 'No arguments => select all pairs.'//
     .  t4,'sum', t11, '[#[,#]] : summarizes site pair table'/
     .            t11, 'Optional [#,#] gives detailed information ',
     .                  'about a specific range of pairs'/
     .  t4,'sum', t11, 'ib=#1,#2 : shows hamiltonian subblocks ',
     .                  'for sites in (#1,#2)'/
     .  t4,'sum', t11, 'ib : summarizes hamiltonian subblocks '
     .                  'for all sites'/)
        print 311
  311   format(
C     .  t4,'sub',  t11, 'r1,r2,c1,c2 : delimit subblock of hrs ',
C     .              '(used by other options)'/
C     .  t4,'sub',  t11, '              work with entire block of hrs'//
     .  t4,'sloc'  t11, 'retain only block diagonal parts of sigma'//
     .  t4,'scell[options] snam : map hrs into a supercell specified by site file ''snam'''/
     .             t11, 'options:'/
     .             t11, '@sigq writes sigma(k)'/
     .             t11, '@sigrs writes sigma(T)'/
     .             t11, '@range=# limits the range of sigma(T) when matching'/
     .             t11, '@nq=#,#,# sets the supercell k mesh'//
     .  t4,'sigbar',t11,'#1 [#2]: set avg sigma for high-lying states'/
     .              t11,'Use argument  #2  in spin polarized case'//

     .  t4,'rssave',t11,'[fn] : write rs sigma hrs to file'/
     .              t11,"If not specified, file name is 'rssigm'."//

     .  t4,'readx',t11,'[fn] : read rs sigma and associated files',
     .                  ' to auxilliary storage'/
     .             t11,'for subsequent pasting into hrs.'/
     .              t11,"If not specified, file name is 'rssigm'."//

     .  t4,'merge',t11,'[sigbar] site-list : merge pairs in hrsx to '
     .                    'corrresponding pairs in hrs '/
     .             t11,'for all pairs connected to sites in site-list'/
     .  t4,'merge',t11,'[sigbar] 0 :  merge hrsx into hrs for all pairs',
     .                  ' in common'//
     .  t11,"optional '-sigbar' imports avg sigma from auxilliary file",
     .    "reducing size of rs sigma"//


     .  t4,'shows' t11, '[r1,r2,c1,c2] shows subblock of sigma, for '/
     .             t11, 'list of pairs specified by ''pair'' command')

        print 312
  312   format(/

     .  t4,'spin[options]'/
     .             t11, 'options:'/
     .             t11, '@flip[=list] exchanges (up,down) spins'/
     .             t11, '@zers[=list] zero out spin-up - spin-dn'/
     .             t11, 'optional list performs function for pairs in list'//

     .  t4,'zero',t11,'[-shorten] CONDITION '
     .  'zero hrs for pairs satisfying CONDITION described below'//

     .  t4,'showp',t11,'CONDITION [rms]: ',
     .  'show pairs satisfying CONDITION described below'/
     .  t11,"Optional 'rms' prints out rms sum of sigma elements",
     .    " for each pair"//

     .  t11,'CONDITION defining pairs for zero or showp is one of: '/
     .  t13,'ib=site-list',t56,'pairs connected to site ib'/
     .  t13,'ip=pair-list',t56,'pairs in list and conjugate'/
     .  t13,'r<#1 [ib=site-list] [p<#2 p1,p2,p3]',t56,
     .    'pairs separated by distance r<#1'/
     .  t13,'r>#1 [ib=site-list] [p<#2 p1,p2,p3]',t56,
     .    'pairs separated by distance r>#1'//
     .  t11,'The r< and r> cases optional permit ',
     .    'further restrictions:'/
     .  t13,'ib=site-list',t30,
     .    'Restrict search to sites in site-list'/
     .  t13,'p<#2 p1,p2,p3',t30,
     .    'Both sites in the pair must be within neighborhood'/
     .  t30,'of (p1,p2,p3), i.e. for each site,'/
     .  t30,'d = |(px-p1,py-p2,pz-p3)| where (px,py,pz) = site posn'/
     .  t13,'Example:  r<.45 ib=2,3 p<.45 -.25,-.25,-.25 rms'/

     . /t4,'q',t11,'quit the editor, writing sigma file'/
     .  t4,'a',t11,'abort without writing sigma file'/
     .  )
        goto 10

      else
        print '(1x,''unrecognized option: '',a10)', outs
        goto 10

      endif

C --- Cleanup and exit rst editor ---
  199 continue
      return

C --- Parsing failure ---
   97 call info0(0,0,0,' siged:  improper arguments to '//
     .  trim(outs)//' ... nothing done')
      goto 299

   98 call info0(0,0,0,' siged:  improper arguments: '//
     .  trim(outs(j1:))//' ... nothing done')
      goto 299

   99 call info0(0,0,0,' siged:  failed to parse arg starting at: '//
     .  trim(outs(j1:))//' ... nothing done')
      goto 299

  299 continue
      if (allocated(ilst)) deallocate(ilst)
      if (allocated(rtab)) deallocate(rtab)
      if (allocated(lzero)) deallocate(lzero)
      goto 10

      end
      subroutine pshorten(nbas,nttab,ndhrs,nsp,lzero,ntab,iax,hrsz)
C- Eliminate pairs from potential represented in standard real-space format
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ndhrs :dimensions hrsz
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lzero :lzero(i) = T => retain pair i in iax,hrs
Ci   lzero :lzero(i) = F => remove pair i from iax,hrs
Cio Inputs/Outputs
Cio  ... each of these must be supplied on input and
Cio      may be shortened on output
Cio  nttab :total number of pairs in neighbor and iax (pairc.f)
Cio  iax   :neighbor table containing pair information (pairc.f)
Cio  hrsz  :R.S. hamiltonian
Co Outputs
Co  ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Co        :
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   Shortens
Cu Updates
Cu   23 Aug 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nttab,ndhrs,nsp
      integer niax
      parameter (niax=10)
      logical lzero(nttab)
      integer ntab(nbas+1),iax(niax,nttab)
      double complex hrsz(ndhrs,ndhrs,nsp,nttab)
C ... Local Dynamical arrays
      integer, allocatable :: imap(:,:)
C ... Local parameters
      integer i,ib,nelt,nttabn,oldi,newi,oldci,newci

      allocate(imap(nttab,2)) ! Mapping of old to new table

C ... Condense ntab, iax, hrs
      nttabn = 1
      nelt = ndhrs*ndhrs*nsp
      do  i = 1, nttab
        imap(i,1) = nttabn      ! New location in array
        imap(i,2) = iax(6,i)    ! Original location of conjugate pair
        ib = iax(1,i)
        ntab(ib+1) = nttabn-1 ! possibly no occurences in new pair table for this site
        if (lzero(i)) cycle ! Eliminate this pair
        if (nttabn /= i) then
          call icopy(niax,iax(1,i),1,iax(1,nttabn),1)
          call dcopy(2*nelt,hrsz(1,1,1,i),1,hrsz(1,1,1,nttabn),1)
        endif
        nttabn = nttabn+1
        ntab(ib+1) = nttabn-1
      enddo

C ... Update iax(6,:)
      do  oldi = 1, nttab
        if (lzero(oldi)) cycle ! Pair no longer in table
        newi = imap(oldi,1)
        oldci = imap(oldi,2)
        newci = imap(oldci,1)
        iax(6,newi) = newci
        iax(6,newci) = newi
      enddo

      call info2(30,0,0,' pshorten:  R.S. array reduced '//
     .  'from %i to %i pairs',nttab,nttabn-1)
      nttab = nttabn-1

      deallocate(imap)

      end

      subroutine sig_pair_list(sopts,kw,nbas,plat,pos,nttab,ntab,iax,
     .  lincl,lprint,j1)
C- Extract a list of pairs from one of a possible set of specifications
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :arguments specifying which pair to include
Ci         :sopts should take one of the following forms:
Ci         : [no arguments]
Ci           returns all pairs
Ci         : r<#           [trailing args**]
Ci         : returns pairs separated by distance r<#
Ci         : # is a length (dimensionless; alat*# is true length)
Ci         : r>#           [trailing args]
Ci         : p<# #1,#2,#3  [trailing args]
Ci         : ib=site-list* [trailing args]
Ci         : ip=pair-list* [trailing args]
Ci         :*  see slatsm/mkilst for syntax of integer lists
Ci         :** trailing args affects lprint.  At present can be only rms.
Ci         : Example (k=2):   showp r<.5 rms
Ci     kw  :Parse sopts starting at the kth word
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors, in units of alat
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Outputs
Co   lincl :One value for each element in pair table:
Co         :T include this pair
Co         :F exclude this pair
Co   lprint:0 by default
Co         :1 if final argument 'rms'  encountered
Co   j1    :0   illegal input value; do nothing
Co         :<0  improper arguments; starting at -j1
Co         :>0  normal return; sopts parsed to j1
Cu Updates
Cu   14 Oct 13  Adapted from siged
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nttab,niax,lprint,j1,kw
      parameter (niax=10)
      character sopts*(*)
      double precision plat(3,3),pos(3,nbas)
      integer iax(niax,nttab),ntab(nbas+1)
      logical lincl(nttab)
C ... Local Dynamical arrays
      real(8), allocatable :: rtab(:,:),posx(:,:)
      integer, allocatable :: ilst(:)
C ... Local parameters
      logical lin,ltmp
      integer i,ib,jb,ip,isw,j,j2,jp,k,k2,n_sitel,nsites,stdo
      integer ix(100)
      double precision xx(10)
      procedure(real(8)) :: dlength
      procedure(integer) :: nglob,a2vec,mkilsd,iinear

C ... Setup for all variations of grouping pairs
      stdo = nglob('stdo')
      nsites = 0                ! No site list by default
      allocate(rtab(3,nttab))
C     call mkrtab(000,alat,plat,pos,iax,nttab,pos,rtab)
      call mkrtab(000,1d0,plat,pos,iax,nttab,pos,rtab) ! alat=1 => rtab,pos same dim
      lprint = 0                ! No extra arguments by default
      k = kw                    ! First argument in sopts
      call word(sopts,k,j1,j2)
      if (j2 <= j1) then      ! Nothing more to do
        lincl = .true.          ! All pairs included in the absence of specifications
        goto 90
      endif
      lincl = .false.           ! Initially all pairs excluded; populate included ones

C ... Conditions 'r<', 'r>'
      if (index(sopts(j1:j2),'r<') == 1 .or.
     .    index(sopts(j1:j2),'r>') == 1) then
        k = k+1
        ltmp = index(sopts(j1:j2),'r<') == 1
        lincl = .false.         ! Initially all pairs excluded; populate included ones
        j1 = j1+1
        xx(1) = 0
        if (a2vec(sopts,j2,j1,4,', ',2,-3,1,ix,xx) /= 1) then
          j1 = -j1; goto 99
        endif
        call info2(10,0,0,' Pairs for range %?#n#<#># %d'//
     .    '%N site  nclus  incl   excl',isw(ltmp),xx)
        if (xx(1) <= 0) then
          call info0(0,0,0,' range <=0 ... nothing done')
          j1 = 0; goto 99
        endif

C   ... Look for restrictions on sites
        nsites = n_sitel(sopts,k,'ib=',j1)
        if (nsites < 0) goto 99
        if (nsites > 0) then
          allocate(ilst(nsites))
          if (mkilsd(sopts(j1+3:),nsites,ilst) /= nsites)
     .      call rx('bug in mkilsd')
        endif

C   ... For each site, do
        do  ib = 1, nbas
          k2 = iinear(nsites,ib,ilst,1)
          if (k2 /= 0) then   ! a list: cycle if ib not in it
            if (ilst(k2) /= ib) cycle
          endif
          k2 = 0
          do  ip = ntab(ib)+1, ntab(ib+1)
            lin = dlength(3,rtab(1,ip),1) < xx(1)
C           Bug in gfortran ... do the long way
C           lincl(ip) = lin .xor. (.not. ltmp)
            if (lin) then
              lincl(ip) = .false. ! Pair excluded
              if (ltmp) lincl(ip) = .true. ! Pair included
            else
              lincl(ip) = .true. ! Pair included
              if (ltmp) lincl(ip) = .false. ! Pair included
            endif
            if (lincl(ip)) then
              k2 = k2+1
              jp = iax(6,ip)
              lincl(jp) = .true. ! Conjugate pair included
            endif
          enddo
          i = ntab(ib+1)-ntab(ib)
          write(stdo,354) ib,i,k2,i-k2
  354     format(i4,3i7)
        enddo


C   ... Check for additional restrictions on absolute position of sites
        call word(sopts,k,j1,j2)
        if (j2 <= j1) goto 90
        if (index(sopts(j1:j2),'p<') /= 1) goto 90
        k = k+1
        j1 = j1+1
        xx(1) = 0
        if (a2vec(sopts,j2,j1,4,', ',2,-3,1,ix,xx) /= 1) then
          j1 = -j1; goto 99
        endif
        call word(sopts,k,j1,j2)
        k = k+1
        j = j1-1
        if (a2vec(sopts,j2,j,4,', ',2,-3,-3,ix,xx(2)) /= 3) then
          j1 = -j1; goto 99
        endif
        call info5(10,0,0,' Pairs for |pos-(%d,%d,%d)| < %d'//
     .    '%N site  nclus  incl   excl',xx(2),xx(3),xx(4),xx,0)

C   ... Shift positions by (p1,p2,p3) and shorten
        allocate(posx(3,nbas))
        call dcopy(3*nbas,pos,1,posx,1)
        do  ib = 1, nbas
          posx(1:3,ib) = posx(1:3,ib) - xx(2:4)
        enddo
        call shorps(nbas,plat,(/52,2,2/),posx,posx)

C       Exclude pairs for which either ib or jb not within (p1,p2,p3)
        do  ib = 1, nbas
          k2 = 0
          do  ip = ntab(ib)+1, ntab(ib+1)
            if (.not. lincl(ip)) cycle
            k2 = k2+1
            jb = iax(2,ip)
            if (dlength(3,posx(1,ib),1) < xx(1) .and.
     .          dlength(3,posx(1,jb),1) < xx(1)) cycle

            k2 = k2-1
            lincl(ip) = .false. ! Pair excluded
            jp = iax(6,ip)
            lincl(jp) = .false. ! Conjugate pair excluded
          enddo
          i = ntab(ib+1)-ntab(ib)
          write(stdo,354) ib,i,k2,i-k2
        enddo

C ... Pairs connected to specified site or pair indices
      elseif (index(sopts(j1:j2),'ib=') == 1 .or. ! Site
     .        index(sopts(j1:j2),'ip=') == 1) then ! Pair

C   ... Make list
        nsites = n_sitel(sopts,k,sopts(j1:j1+2),j1)
        if (nsites <= 0) goto 99
        allocate(ilst(nsites))
        if (mkilsd(sopts(j1+3:),nsites,ilst) /= nsites)
     .    call rx('bug in mkilsd')

        if (index(sopts(j1:j2),'ib=') == 1) then ! Site case
          call info2(10,0,0,' Show pairs connected to '//
     .      '%?#n>1#%-1j%i sites: %-1j%n:1i%-2j#%-1j#'//
     .      '%?#n==1#%-1jsite %i##',nsites,ilst)
          call imxmn(nsites,ilst,1,i,j)
          if (i < 0 .or. j > nbas) then
            j1 = -1; goto 99
          endif

          do  i = 1, nsites
            ib = ilst(i)
C           Loop over all pairs connected to ib
            j1 = ntab(ib)
            do  ip = ntab(ib)+1, ntab(ib+1)
              lincl(ip) = .true. ! Pair included
              jp = iax(6,ip)
              lincl(jp) = .true. ! Pair included
            enddo
            write(stdo,354) ib,ntab(ib+1)-ntab(ib)
          enddo

        else                    ! pair case
          call info2(10,0,0,' Pair table for '//
     .      '%?#n>1#%-1j%i pairs: %-1j%n:1i%-2j#%-1j#'//
     .      '%?#n==1#%-1jpair %i##',nsites,ilst)
          call imxmn(nsites,ilst,1,i,j)
          if (i < 0 .or. j > nttab) then
            j1 = -1; goto 99
          endif

          do  i = 1, nsites
            ip = ilst(i)
            lincl(ip) = .true.  ! Pair included
            jp = iax(6,ip)
            lincl(jp) = .true.  ! Pair included
          enddo

        endif
      endif


C ... Switches specifying what property to print
   90 continue
      call word(sopts,k,j1,j2)
      if (j2 <= j1) goto 99
      if (sopts(j1:j2) == 'rms') then
        lprint = 1
      else
        j1 = -j1; goto 99
      endif

   99 continue
      if (allocated(ilst)) deallocate(ilst)
      if (allocated(rtab)) deallocate(rtab)

      end

      integer function n_sitel(sopts,k,token,j1)
      implicit none
      integer k,j1
      character sopts*150,token*3
      integer mkilsd,j,j2


      call word(sopts,k,j1,j2)
      n_sitel = 0
      if (index(sopts(j1:j2),token) /= 1) return
      k = k+1
      n_sitel = mkilsd(sopts(j1+3:j2),-1,j)
      if (n_sitel <= 0) then
        call info0(0,0,0,'%6p... Bad or null list : '//sopts(j1:j2))
        j1 = -j1
      endif
      end
