      subroutine stackcel(sopts,s_ctrl,s_lat,s_site,s_spec,slabl)
C- Stacks cells together to make a supercell
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nspec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos
Co     Stored:     pos
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  ssiteadd
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  sspecadd
Ci Inputs
Ci   sopts :string defining joined cell
Ci   slabl :vector of species labels
Cr Remarks
Cr Silly example (site site3.nb from nb-fe tutorial).
Cr Add one plane, shift it, and insert another plane in the middle
Cr lmscell -vfile=3 ctrl.nb --stack~file@site3.nb~addpos@dx=0,0,1@targ=4:6~file@site3.nb@insert=4@dpos=0,0,1~show@planes
Cb Bugs
Cu Updates
Cu   30 Jan 18  Add autogeneration of principal layer assignments for layer code
Cu   09 Dec 17  Added insert option to file instruction
Cu   10 Jul 17  Added sort instruction
Cu   04 Jul 17  First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*(*) sopts
      character*8 slabl(1)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C     Dummy structures ... used as dummy arguments
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
      type(str_mix)::   s_mix
      type(str_str)::   s_str
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      character*8, pointer :: slabll(:),slabl0(:)
      integer,allocatable :: ilsts(:),iblst(:),iplan(:),natpl(:),ips(:),pl(:)
      real(8),allocatable :: posl(:,:),pos(:,:),pos2(:,:),plan(:),wk(:,:)
      type(str_site), allocatable::  s_site0(:),s_sitel(:),s_ssite(:)
      type(str_spec), allocatable::  s_spec0(:),s_specl(:),s_sspec(:)
      integer, allocatable :: iax(:),ntab(:)
      real(8),allocatable:: qnus(:,:,:,:)
C ... Local parameters
      integer,parameter:: NULLI = -99999
      integer,parameter :: mxspec=256,n0=10,niax=10
      real(8),parameter :: tiny = 1d-4

      character*256 outs,dc*1,dc2*1,fn*120,sortex(3),strn,errmsg

      logical ltmp,lsopts,lall
      logical :: newspec = .false.
      integer i,j,ib,is,lio,ifi,jfi,nbas,nbas0,nspec,nsp,j1,j2,js1,js2,stdo,iv(12),ishores(3),iany(1)
      integer nspecl,nbasl,nspecx,nsitex,nplane,ndup,nexpr,n,iv0,nblst,nlst,nlsts,insert,lrs
      integer nbasp,nbaspp,npadl,npadr,mxcsiz,nttab
      double precision pi,alatl,platl(3,3),platr(3,3),xv(10),normal(3,3),any(1)
      double precision plat(3,3),plx(3,3),plati(3,3),qlat(3,3)
      double precision tol,fptol,alat,xx,h(3),vol(3),z
      double precision scale,stretch,dpos(3),pad(3),shft(3),range
      procedure(logical) :: cmdopt,latvec,a2bin,rdstrn
      procedure(integer) :: fxst,fopn,fhndl,fopna,fopnx,fopng,parg2,nglob
      procedure(integer) :: isw,lgunit,rdm,iprint,iosite,iorsa,a2vec,wordsw,iosits,mkilsd
      procedure(real(8)) :: dotprd,dlength,avwsr,ddot,dglob

      alat = s_lat%alat
      plat = s_lat%plat
      call dinv33(plat,1,qlat,xx)
      nbas = s_ctrl%nbas
      npadl = 0; npadr = 0
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      tol = 1d-6
      fptol = 1d-5
      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
!      print *, fhndl('site')
!      call fclr('site',-1)

      call ptr_pot(s_pot,0,' ',0,0,[0])
      s_bz%ef=0; s_bz%def=0; s_bz%w=0; s_bz%n=0

C ... Make s_spec0 and s_ssite0 initial s_spec and s_site
      allocate(s_site0(nbas),s_spec0(nspec),slabl0(nspec))
      nspecx = 0; nsitex = 0
      call sspecadd(nspecx,s_spec0,nspec,0,[0],s_spec) ! copy existing s_spec to s_spec0
      forall (ib=1:nbas) s_site(ib)%qnu(1,1,1) = NULLI
      call ssiteadd(nsitex,s_site0,nbas,0,[0],s_site) ! copy existing structure
      forall (is=1:nspec) slabl0(is) = slabl(is)

      dc = sopts(1:1)
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the superlattice editor.  Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the superlattice editor.  Enter ''?'' to see options.')
        lsopts = .false.
      endif

C --- Return here to resume parsing for arguments ---
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

C     Entry point for instruction given in outs
   20 continue
      call locase(outs)
      nbas0 = nbas   ! Number of sites prior to instruction
      nbasp = nbas + npadl + npadr

C ... Parse and execute the next command

C --- Null input ---
      if (outs == ' ') then
        print 304
  304   format(' enter ''q'' to exit, ''a'' to abort',' ''?'' to see menu')
        goto 10

C --- Read parameters new insertion site file ---
      elseif (outs(1:5) == 'rsasa')  then

        call rdasapq(outs(6:),dc,s_ctrl,s_site0,s_spec0,s_pot,1,nbas,0,errmsg)
        if (errmsg /= ' ') goto 97
        goto 10

C --- Read parameters new insertion site file ---
      elseif (outs(1:4) == 'file')  then

        dc2 = outs(5:5)
        j1 = 6
        call nwordg(outs,1,dc//dc2//' ',1,j1,j2)
        if (j1 > j2) then
          call info0(10,0,0," no file specified ... nothing done")
          goto 10
        endif
        fn = outs(j1:j2)
        ifi = fopnx(fn,72,-1,-1)
        if (ifi /= 1) then
          call info0(10,0,0," missing file '"//trim(fn)//"' ... nothing done")
          goto 10
        endif
        ifi = fopng(fn,-1,0)
        scale = nulli; stretch = nulli; dpos = nulli; pad = nulli; ndup = 1

C   ... Scale insertion lattice vectors
        i = wordsw(outs,dc2,'scale','= ',j1);
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,1,iv,scale) /= 1) goto 98
        endif

C   ... Stretch insertion lattice vectors along 3rd axis
        i = wordsw(outs,dc2,'stretch','= ',j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,1,iv,stretch) /= 1) goto 98
        endif

C   ... Shift insertion basis vectors
        i = wordsw(outs,dc2,'dpos','= ',j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,3,iv,dpos) /= 3) goto 98
        endif

C   ... Read pad(1), which adds plat(:,3)*pad(1) to plat(:,3) before insertion
        i = wordsw(outs,dc2,'pad','= ',j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,1,iv,pad) < 1) goto 98
        endif

        insert = NULLI
        i = parg2(outs,dc2,'insert=',',',' '//dc2,2,1,1,iany,insert,any)
        if (insert /= NULLI .and. (insert <= 0 .or. insert >= nbas)) then
          call info0(10,0,0," illegal insert ... nothing done")
          goto 10
        endif

C        i = wordsw(outs,dc2,'insert','= ',j1)
C        if (i /= 0) then
C          if (a2vec(outs,len_trim(outs),j1,2,', '//dc,3,3,1,iv,insert) /= 1) goto 98
C          if (insert <= 0 .or. insert >= nbas) call rx('bad insert')
C        endif

C   ... Read ndup, which does the insertion ndup times
        i = wordsw(outs,dc2,'dup','= ',j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,2,' '//dc2,2,2,1,iv,ndup) < 1) goto 98
        endif

C   ... Read pad(2), which shifts new basis vectors by plat(:,3)*pad(2) before insertion
        i = wordsw(outs,dc2,'shft','= ',j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,1,iv,pad(2)) < 1) goto 98
        endif

        if (wordsw(outs,dc2,'rsasa','',lrs) == 0) lrs = 0

        goto 30

C --- Removes sites from list ---
C     rmsite@targ=list
      elseif (outs(1:6) == 'rmsite')  then

        j1 = 8
        dc2 = outs(j1-1:j1-1)

C       Required site list
        if (wordsw(outs,dc2,'targ=','',j1) /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          nlsts = mkilsd(outs(j1:j2),-1,iv)
          allocate(ilsts(nlsts))
          if (nlsts <= 0) call rxs(' bad or null site list : ',outs(j1:j2))
          nlsts = mkilsd(outs(j1:j2),nlsts,ilsts)  ! Do not sort
          if (maxval(ilsts(1:nlsts)) > nbas) call rxs(' bad or null site list : ',outs(j1:j2))
          if (nlsts <= 0 .or. nlsts > nbas) goto 98
        else
          call info0(10,0,0," no list specified ... nothing done")
          goto 10
        endif
        call ilst2a(ilsts,nlsts,strn)
        call info2(10,1,0,' remove %?;n>1;%-1j%i sites :;site; '//trim(strn),nlsts,2)

        allocate(iblst(nbas))
        forall (ib = 1:nbas) iblst(ib) = ib
        forall (i = 1:nlsts) iblst(ilsts(i)) = 0
        call siteperm(nbas,iblst,s_site0)
        s_ctrl%nbas = nbas
        deallocate(ilsts,iblst)
        outs = 'show pos'
        goto 20

C --- Compare current site positions to file ---
      elseif (outs(1:6) == 'cmppos')  then

        j1 = 8
        dc2 = outs(j1-1:j1-1)

C   ... Read file data into pos
        i = wordsw(outs,dc2,'fn=','',j1)
        if (i == 0) then
          call info0(10,0,0," no file specified ... nothing done")
          goto 10
        endif
        call nwordg(outs,0,dc//dc2,1,j1,j2)
        fn = outs(j1:j2)
        ifi = fopnx(fn,70,-1,-1)
        if (ifi /= 1) then
          call info0(10,0,0," missing file '"//trim(fn)//"' ... nothing done")
          goto 10
        endif
        ifi = fopna(fn,-1,0)

C   ... Read file data into pos; load current positions into posl
        allocate(posl(3,nbas),pos(3,nbas))
        call sitepack(s_site,1,nbas,'pos',3,xx,posl)
        call ioposlu(.false.,-1,ifi,nbas,pos,s_site)
        call fclr(' ',ifi)

C       Optional site list
        if (wordsw(outs,dc2,'targ=','',j1) /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          nlsts = mkilsd(outs(j1:j2),-1,iv)
          allocate(ilsts(nlsts))
          if (nlsts <= 0) call rxs(' bad or null site list : ',outs(j1:j2))
C         call mkilssr(11,outs(j1:j2),nlsts,ilsts,[1,nbas])
          nlsts = mkilsd(outs(j1:j2),nlsts,ilsts)  ! Do not sort
          if (maxval(ilsts(1:nlsts)) > nbas) call rxs(' bad or null site list : ',outs(j1:j2))
          if (nlsts <= 0 .or. nlsts > nbas) goto 98
        else
          allocate(ilsts(nbas))
          forall (ib=1:nbas) ilsts(ib) = ib
          nlsts = nbas
        endif
        call ilst2a(ilsts,nlsts,strn)
        call info2(10,1,0,' compare positions for %?;n>1;%-1j%i sites :;site; '//trim(strn),nlsts,2)

C       Set ltmp if to shorten change in position
        ltmp = wordsw(outs,dc2,'shorten',dc//dc2//' ',j1) > 0

C   ... Printout
        call info0(2,1,0,'%15fGiven%31fFile%32fshift')
        xv(1:3) = 0
        do  i = 1, nlsts
          ib = ilsts(i)
          dpos(:) = pos(:,ib)-posl(:,ib)
          if (ltmp) call shorps(1,plat,(/72,2,2/),dpos,dpos)
          call info5(2,0,0,'%3;11,6D   %3;11,6D   %3;11,6D',
     .      posl(1,ib),pos(1,ib),dpos,4,5)
          pos(:,ib) = dpos(:)
          xv(1:3) = xv(1:3) + dpos(:)/nlsts
        enddo
        call info2(2,0,0,'%57fAverage shift =%3;11,6D',xv,2)

C   ... Write shifts into file
        i = wordsw(outs,dc2,'wdx=','',j1)
        if (i > 0) then
          call nwordg(outs,0,dc2,1,j1,j2)
          call iopos(.true.,0,outs(j1:j2),nbas,pos,s_site)
        endif

        deallocate(posl,pos,ilsts)
        goto 10

C --- Replace or shift a subset of current site positions ---
C     newpos|addpos[@shorten=#,#,#][@targ=lst][@src=list]  @dp=#,#,# | @dx=#,#,# | @fn=filename
      elseif (outs(1:6) == 'newpos' .or. outs(1:6) == 'addpos') then

        j1 = 8
        dc2 = outs(j1-1:j1-1)
        allocate(pos(3,nbas),ips(nbas))
        call sitepack(s_site0,1,nbas,'spec',1,ips,xv)
        call sitepack(s_site0,1,nbas,'pos',3,xx,pos)

        ishores = 0
        i = wordsw(outs,dc2,'shorten=','',j1) + 8
        if (i > 8) then
          is = a2vec(outs,len_trim(outs),i,2,', '//dc2,3,2,3,iv,ishores)
          if (is /= 3) call rx('stackel: failed to parse '//trim(outs))
        endif

C       Optional site list
        if (wordsw(outs,dc2,'targ=','',j1) /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          nlsts = mkilsd(outs(j1:j2),-1,iv)
          allocate(ilsts(nlsts))
          if (nlsts <= 0) call rxs(' bad or null site list : ',outs(j1:j2))
C         call mkilssr(11,outs(j1:j2),nlsts,ilsts,[1,nbas])
          nlsts = mkilsd(outs(j1:j2),nlsts,ilsts)  ! Do not sort
          if (maxval(ilsts(1:nlsts)) > nbas) call rxs(' bad or null site list : ',outs(j1:j2))
          if (nlsts <= 0 .or. nlsts > nbas) goto 98
        else
          allocate(ilsts(nbas))
          forall (ib=1:nbas) ilsts(ib) = ib
          nlsts = nbas
        endif
        call ilst2a(ilsts,nlsts,strn)
        call info2(10,1,0,' shift positions for %?;n>1;%-1j%i sites :;site; '//trim(strn),nlsts,2)
        allocate(posl(3,nlsts)); call dpzero(posl,size(posl))

C   ... Read the shift according to mode
        ltmp = .false.
        i = wordsw(outs,dc2,'dx','= ',j1) ! Fixed shift, cartesian coordinates
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,3,iv,xv(4)) /= 3) goto 98
C         posp+ = (plat)^-1 pos+
          call dgemm('T','N',3,1,3,1d0,qlat,3,xv(4),3,0d0,xv(7),3)
          ltmp = .true.
        endif
        i = wordsw(outs,dc2,'dp','= ',j1) ! Fixed shift, units of plat
        if (i /= 0 .and. .not. ltmp) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,3,iv,xv(7)) /= 3) goto 98
C         Convert pos to Cartesian coordinates
          call dgemm('N','N',3,1,3,1d0,plat,3,xv(7),3,0d0,xv(4),3)
          ltmp = .true.
        endif
        if (ltmp) then
          call info2(10,0,0,' uniform shift :%3;10,6D (Cart) =%3;10,6D (plat)',xv(4),xv(7))
          forall (i=1:nlsts) posl(:,i) = xv(4:7)
        else ! Read shifts from file
          i = wordsw(outs,dc2,'fn=','',j1)
          if (i == 0) then
            call info0(10,0,0," no shift specified ... nothing done")
            goto 10
          endif
          call nwordg(outs,0,dc2,1,j1,j2)
          fn = outs(j1:j2)
          ifi = fopnx(fn,72,-1,-1)
          if (ifi /= 1) then
            call info0(10,0,0," missing file '"//trim(fn)//"' ... nothing done")
            goto 10
          endif
C         Extract size of positions file; read into pos2
          call info0(2,0,-1,' read positions from file '//trim(fn)//'.')
          ifi = fopng(fn,-1,0)
          n = 0
          call ioposlu(.false.,0,ifi,n,xv,s_site)
          allocate(pos2(3,n)); call dpzero(pos2,size(pos2))
          rewind ifi
          call ioposlu(.false.,0,ifi,n,pos2,s_site)

C     ... Copy subset of pos2 into posl
C         Optional src list
          if (wordsw(outs,dc2,'src=','',j1) /= 0) then
            call nwordg(outs,0,dc2//' ',1,j1,j2)
            nblst = mkilsd(outs(j1:j2),-1,iv)
            allocate(iblst(nblst))
            if (nblst <= 0) call rxs(' band or null band list,',outs)
C           call mkilssr(11,outs(j1:j2),nblst,iblst,[1,nbas]) ! do not sort
            nblst = mkilsd(outs(j1:j2),nblst,iblst)
            if (nblst <= 0 .or. nblst > nbas) goto 98
            call ilst2a(iblst,nblst,strn)
            call info2(10,0,0,' src %?;n>1;%-1j%i sites;site; from file : '//trim(strn),nblst,2)
            if (nblst < nlsts) then
              call info0(10,0,0," fewer source elements than destination elements ... nothing done")
              goto 10
            endif
            if (maxval(iblst(1:nblst)) > nbas) call rxs(' bad or null site list : ',outs(j1:j2))
          else
            allocate(iblst(nbas))
            forall (i=1:nbas) iblst(i) = i
            nblst = nbas
          endif
C         Make a list of positions
          if (nblst > nlsts) then
            call info2(10,0,0," more source sites than destination sites ... use first %i",nlsts,2)
          endif
          forall (i=1:nlsts) posl(:,i) = pos2(:,iblst(i))
          deallocate(iblst,pos2)
          call fclr(' ',ifi)
        endif
C       At end of this block, nlsts and posl(1:nlsts) shall be avalailable

C   ... For each element in list, replace pos or add posl to pos
        ltmp = outs(1:6) == 'newpos'  ! Copy, rather than add
        iv(1:3) = ishores; iv(1) = iv(1)+60
        write(stdo,357)
        do  i = 1, nlsts
          ib = ilsts(i)
          is = ips(ib)
          xv(1:3) = pos(1:3,ib) + posl(:,i)
          if (ltmp) xv(1:3) = posl(:,i)
          if (sum(ishores) /= 0) then
            call shorps(1,plat,iv,xv,xv)
          endif
C         posp+ = (plat)^-1 pos+
          call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
          pos(1:3,ib) = xv(1:3)
          print 345, ib, slabl(is), (xv(j),j=1,6)
        enddo

        call sitepack(s_site0,1,nbas,'-pos',3,xx,pos)
        call dcopy(3*nbas,pos,1,s_lat%pos,1)

C   ... Write shifts into file
        i = wordsw(outs,dc2,'wdx=','',j1)
        if (i > 0) then
          call nwordg(outs,0,dc2,1,j1,j2)
          call iopos(.true.,0,outs(j1:j2),nbas,pos,s_site0)
        endif

        deallocate(posl,ilsts,pos,ips)
        goto 10

C --- scale lattice
      elseif (outs(1:5) == 'scale') then
        j1 = 6
        if (a2vec(outs,len_trim(outs),j1,4,', '//dc,3,3,1,iv,stretch) /= 1) goto 98
        alat = alat*stretch

        call info2(10,1,0,'   Lattice constant scaled to %,6;6d',alat,2)

        outs = 'show size'
        goto 20

C --- Stretch lattice
      elseif (outs(1:7) == 'stretch') then
        dc2 = outs(8:8)
        i = wordsw(outs,'','stretch',' '//dc2,j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc,3,3,1,iv,stretch) /= 1) goto 98
        endif

C       Convert positions to multiples of plat
        call dinv33(plat,0,plati,xx)
        allocate(posl(3,nbas),pos(3,nbas))
        call sitepack(s_site0,1,nbas,'pos',3,xx,posl)
        call dgemm('N','N',3,nbas,3,1d0,plati,3,posl,3,0d0,pos,3)
C       Stretch plat
        call dscal(3,stretch,plat(1,3),1)
C       Convert positions back to Cartesian coordinates
        call dgemm('N','N',3,nbas,3,1d0,plat,3,pos,3,0d0,posl,3)
        call sitepack(s_site0,1,nbas,'-pos',3,xx,posl)
        call dcopy(3*nbas,posl,1,s_lat%pos,1)
        deallocate(posl,pos)

        outs = 'show size pos'
        goto 20

C --- Pad third lattice vector ---
      elseif (outs(1:3) == 'pad') then
        dc2 = outs(4:4)
        i = wordsw(outs,'','pad',' '//dc2,j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc,3,3,3,iv,pad) /= 3) goto 98
        endif
        call daxpy(3,1d0,pad,1,plat(1,3),1)

        outs = 'show size pos'
        goto 20

C --- Reorder site positions ---
      elseif (outs(1:7) == 'newspec') then

        dc2 = outs(8:8)
        if (dc2 == ' ') then
          newspec = .not. newspec
        else
          i = wordsw(outs,'','newspec',' '//dc2,j1)
          if (i /= 0) then
            if (a2vec(outs,len_trim(outs),j1,0,', '//dc,3,3,1,iv,newspec) /= 1) goto 98
          endif
        endif
        goto 10

C --- Set PL ---
      elseif (outs(1:5) == 'setpl') then

        j1 = 7
        dc2 = outs(j1-1:j1-1)

C       Read range (required)
        if (wordsw(outs,dc2,'range=','',j1) /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          i = 0
          if (a2vec(outs(j1:j2),len_trim(outs(j1:j2)),i,4,', '//dc2,3,3,1,iv,range) /= 1) goto 98
        else
          call info0(10,0,0," no range specified ... nothing done")
          goto 10
        endif

        if (wordsw(outs,dc2,'rsasa','',lrs) == 0) lrs = 0

        call dinv33(plat,1,qlat,vol)
        call dpcopy(qlat(1,3),normal,1,3,1/dlength(3,qlat(1,3),1))

        allocate(plan(nbas),iplan(nbas),natpl(nbas),ips(nbas),pl(nbas*5),pos(3,nbas),slabll(nbas))
        call sitepack(s_site0,1,nbas,'spec',1,ips,xx)
        call sitepack(s_site0,1,nbas,'pos',3,xx,pos)
        call ivset(pl,1,nbas*5,NULLI) ! initialize pl to NULL to flag not set
        forall (is=1:nspec) slabll(is) = s_spec0(is)%name
        platl = s_lat%platl
        platr = s_lat%platr

        if (wordsw(outs,dc2,'shorps=','',j1) /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          i = 0
          if (a2vec(outs(j1:j2),len_trim(outs(j1:j2)),i,2,', '//dc2,3,3,3,iv,ishores) /= 3) goto 98
          plx = plat
          plx(:,3) = plat(:,3) - 2*platl(:,3) - 2*platr(:,3)
          call shorps(nbas,plx,ishores,pos,pos)
          call sitepack(s_site0,1,nbas,'-pos',3,xx,pos)
        endif

C       Make array plan = height of each atom projected onto to normal
        ltmp = wordsw(outs,dc2,'sort','',j1) > 0  ! Flag to sort by height
   15   continue
        do  ib = 1, nbas
          z = ddot(3,pos(1,ib),1,normal,1) / dlength(3,normal,1)
          plan(ib) = z
          if (ib == 1) cycle
          if (z-plan(ib-1) < -tiny .and. .not. ltmp) then
            call info0(10,0,0," sites not ordered by height ... nothing done")
            goto 10
          endif
        enddo

C       Sort by height h
        if (ltmp) then
          call dvheap(1,nbas,plan,iplan,0d0,101)
          call siteperm(nbas,iplan,s_site0)
          call sitepack(s_site0,1,nbas,'pos',3,xx,pos)
          ltmp = .false.
          goto 15
        endif

C       Update s_lat%pos
        call dcopy(3*nbas,pos,1,s_lat%pos,1)

C       Find all sites in first PL: defined as sites whose normal is < PLATL . normal
C       and  all sites in last  PL: defined as sites whose normal is > (PLAT-PLATR) . normal
C       For now all sites in active layer are PL 1
        xv(1) = dabs(ddot(3,platl(1,3),1,normal,1) / dlength(3,normal,1)) ! highest projection for first PL
        xv(2) = dabs(ddot(3,platr(1,3),1,normal,1) / dlength(3,normal,1)) ! lowest projection for last PL
        xv(3) = dabs(ddot(3,plat(1,3),1,normal,1) / dlength(3,normal,1)) ! lowest projection for last PL
        xv(3) = xv(3) - 2*xv(2) - 2*xv(1)  ! Remove double padding of plat
        if (xv(3) < 0) call rx('faulty PLATL or PLATR')
        do  ib = 1, nbas
          if (plan(ib)-plan(1) < xv(1)-tiny) pl(ib) = 1  ! This will be the first element in the active layer
        enddo
        do  ib = 1, nbas
          if (plan(ib)-plan(1) > xv(3)-xv(2)-tiny) pl(ib) = 0  ! Temporarily 0 ... later becomes 1
        enddo
        npadl = 0; npadr = 0
        do  ib = 1, nbas
          if (pl(ib) == 1) npadl = npadl+1
          if (pl(ib) == 0) npadr = npadr+1
        enddo
        if (npadl*npadr == 0) then
          call info0(10,0,0," no sites in L or R region ... nothing done")
          goto 10
        endif
        s_ctrl%npadl = npadl
        s_ctrl%npadr = npadr
        nbasp = nbas + npadl + npadr
        nbaspp = nbas + 2*(npadl + npadr)
        xv(4) = dglob('nbasp',dble(nbasp),1)

C       Make padded basis
        allocate(posl(3,nbaspp))
        call pgbasp(nbas,npadl,npadr,s_lat%pos,plat,platl,platr,posl)
        call dscal(9,2d0,platl,1)
        call pgbasp(nbasp,npadl,npadr,posl,plat,platl,platr,s_lat%pos)
        call dscal(9,.5d0,platl,1)
        deallocate(posl)

C       Enlarge s_pot%ves
        if (associated(s_pot%ves)) call ptr_pot(s_pot,2,'ves',nbasp,0,[0d0])

C       Make s_sitel ... block assumes all the L PL sites precede the right PL sites
        allocate(s_sitel(nbasp-nbas))
        npadl = 0; npadr = 0
        do  ib = 1, nbas
          if (pl(ib) == 1) then
            npadl = npadl+1
            j = npadl
            pl(nbas+j) = 0
            pl(nbasp+j) = -1
          elseif (pl(ib) == 0) then
            npadr = npadr+1
            j = npadl+npadr
            pl(nbas+j) = 2
            pl(nbasp+j) = 3

          else
            pl(ib) = 1
            cycle
          endif
          pl(ib) = 1
          s_sitel(j)%spec = s_site0(ib)%spec
          s_sitel(j)%pos(1:3) = s_lat%pos(1:3,nbas+j)
          s_sitel(j)%pnu = s_site0(ib)%pnu
          s_sitel(j)%pz = s_site0(ib)%pz
          s_sitel(j)%qnu = s_site0(ib)%qnu
          s_sitel(j)%force = s_site0(ib)%force
          s_sitel(j)%vel = s_site0(ib)%vel
        enddo

C        if (lrs > 0) then
C          print *, nbas
C          stop '674 s_sitel'
C        endif

        if (associated(s_pot%ves)) then
          forall (ib = 1:npadl) s_pot%ves(nbas+ib) = s_pot%ves(ib)
          forall (ib = nbas-npadr+1:nbas) s_pot%ves(ib+npadl+npadr) = s_pot%ves(ib)
        endif

C       Create s_site0 for padded basis
        allocate(s_ssite(nbas))
        nsitex = 0
        call ssiteadd(nsitex,s_ssite,nbas,0,[0],s_site0) ! Copy existing structure to s_ssite
        deallocate(s_site0); allocate(s_site0(nbasp))
        nsitex = 0
        call ssiteadd(nsitex,s_site0,nbas,0,[0],s_ssite) ! Copy existing structure to enlarged s_site0
        call ssiteadd(nsitex,s_site0,nbasp-nbas,0,[0],s_sitel) ! Add padded sites to s_site0
        if (nsitex /= nbasp) call rx('bug in stackcel')
        deallocate(s_ssite,s_sitel)

C       Create neighbor table for padded basis
        mxcsiz = 0; allocate(ntab(nbasp+1))
        call pairs(nbas,nbasp,alat,plat,[range/2],s_lat%pos,
     .    [-1],2,1,pl,nttab,ntab,iax,mxcsiz)

C       Make the PL indices consistent with neighbor table
        call makepl(nbas,npadl,npadr,ntab,iax,pl)

        deallocate(s_ctrl%pgfsl); allocate(s_ctrl%pgfsl(nbaspp))
        call icopy(nbaspp,pl,1,s_ctrl%pgfsl,1)
        forall (ib=1:nbasp) s_site0(ib)%pl = pl(ib)
        forall (ib=1:nbasp) s_site0(ib)%class = ib

        outs = 'show planes'
        deallocate(plan,iplan,natpl,ips,pl,pos,slabll)
        goto 20

C --- Display a position in both Cartesian coordinates and lattice vector multiples  ---
      elseif (outs(1:8) == 'cart2lat' .or. outs(1:8) == 'lat2cart') then

        j1 = 10
        dc2 = outs(j1-1:j1-1)

        ltmp = .false.
        if (wordsw(outs,dc2,'q','',j1) /= 0) ltmp = .true.

        i = wordsw(outs,'','x=',' '//dc2,j1)
        if (i /= 0) then
          j1 = j1-1
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc2,3,3,3,iv,xv) /= 3) goto 98
        else
          call info0(10,0,0,'no coordinate specifed (x=) ... nothing to do')
          goto 10
        endif

        call dcopy(9,plat,1,platl,1)
        if (ltmp) call dcopy(9,qlat,1,platl,1)
        call dinv33(platl,0,plati,xx)

C       Convert positions to multiples of plat
        if (outs(1:8) == 'cart2lat') then
          call dgemm('N','N',3,1,3,1d0,plati,3,xv,3,0d0,xv(4),3)
        else
          call dcopy(3,xv,1,xv(4),1)
        endif

C       Convert positions back to Cartesian coordinates
        call dgemm('N','N',3,1,3,1d0,platl,3,xv(4),3,0d0,xv,3)

        call info0(2,0,0,'     x  (Cartesian coordinates)           x  (multiples of lat)')
        call info2(2,0,0,'%3;11,6D   %3;11,6D',xv(1),xv(4))

        goto 10

C --- Reorder site positions ---
      elseif (outs(1:4) == 'sort') then

        j1 = 6
        dc2 = outs(j1-1:j1-1)

        nexpr = 0
        do  i = 1, 3
          call nwordg(outs,1,dc//dc2,1,j1,j2)
          if (j2 < j1) exit
          nexpr = nexpr+1
          sortex(nexpr) = outs(j1:j2)
          j1 = j2+2
        enddo
        if (nexpr == 0) then
          call info0(10,0,0,'no expressions ... nothing to do')
          goto 10
        endif

        allocate(ips(nbas),pos(3,nbas),posl(3,nbas),wk(nexpr,nbas),ilsts(nbas))
        call sitepack(s_site0,1,nbas,'spec',1,ips,xx)
        call sitepack(s_site0,1,nbas,'pos',3,xx,pos)
        call dcopy(3*nbas,pos,1,posl,1)

C   ... Specify site list directly
        if (wordsw(outs,dc2,'targ=','',j1) /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          nlsts = mkilsd(outs(j1:j2),-1,iv)
          if (nlsts /= nbas) call rxs(' bad or null site list : ',outs(j1:j2))
          call mkilssr(11,outs(j1:j2),nlsts,ilsts,[1,nbas]) ! Make sure no duplicates
          if (nlsts /= nbas) call rxs(' bad or null site list : ',outs(j1:j2))
          nlsts = mkilsd(outs(j1:j2),nlsts,ilsts)  ! Do not sort
          if (nlsts <= 0 .or. nlsts > nbas) goto 98

C   ... Site list through algebraic expressions
        else

          call numsyv(iv0)
          call dinv33(plat,1,qlat,vol)
          do  ib = 1, nbas
C           call shosyv(0,0,0,6)

C           projection of bas along qlat(i)
            do  i = 1, 3
              h(i) = ddot(3,pos(1,ib),1,qlat(1,i),1)
              dpos(i) = ddot(3,pos(1,ib),1,plat(1,i),1)
            enddo

            call lodsyv('ib',1,dble(ib),n)
            i = ips(ib)
            call lodsyv('is',1,dble(i),n)
            call lodsyv('z',1,s_spec(i)%z,n)
            call lodsyv('x1',1,pos(1,ib),n)
            call lodsyv('x2',1,pos(2,ib),n)
            call lodsyv('x3',1,pos(3,ib),n)
            call lodsyv('h1',1,h(1),n)
            call lodsyv('h2',1,h(2),n)
            call lodsyv('h3',1,h(3),n)
C           call lodsyv('h',1,h(3),n)
            call lodsyv('p1',1,dpos(1),n)
            call lodsyv('p2',1,dpos(2),n)
            call lodsyv('p3',1,dpos(3),n)

C           call shosyv(0,0,0,6)

            do  i = 1, nexpr
              j = 0
              if (.not. a2bin(sortex(i),wk(i,ib),4,0,' ',j,-1)) then
                call info0(1,0,0,'... failed to parse expr:'//trim(sortex(i))//
     .            ' ... nothing done')
                call clrsyv(iv0)
                goto 10
              endif
            enddo
            call clrsyv(iv0)
          enddo
          call dvheap(nexpr,nbas,wk,ilsts,0d0,101)
        endif

C       Make permutation table ilsts, permute site structure
        call siteperm(nbas,ilsts,s_site0)
        call sitepack(s_site0,1,nbas,'pos',3,xx,pos)
        call dcopy(3*nbas,pos,1,s_lat%pos,1)

        deallocate(ips,pos,posl,wk,ilsts)
        outs = 'show pos'
        goto 20

C --- Assemble basp file ---
      elseif (outs(1:4) == 'basp') then

C       Get file name and write header
        j1 = 6
        dc2 = outs(j1-1:j1-1)
        i = wordsw(outs,dc2,'fn=','',j1)
        call nwordg(outs,1,dc//dc2//' ',1,j1,j2)
        if (i == 0 .or. j1 > j2) then
          call info0(10,0,0," no file specified ... nothing done")
          goto 10
        endif
        fn = outs(j1:j2)
        ifi = fopng(fn,-1,0)
        rewind ifi; write(ifi,'(''BASIS:'')')

        j1 = 5
        do while (outs(j1:j1) == dc2)
          call nwordg(outs,1,dc//dc2//' ',2,j1,j2)
          if (outs(j1:j1+2) == 'fn=') goto 80  ! this is the output file

          jfi = fopnx(outs(j1:j2),72,-1,-1)
          if (jfi /= 1) then
            call info0(10,0,0," missing file '"//trim(outs(j1:j2))//"' ... nothing done")
            call fclr(' ',ifi)
            goto 10
          endif
          jfi = fopng(outs(j1:j2),-1,0)
          call info0(10,0,0,' concatenating file '//outs(j1:j2))
          do while (rdstrn(jfi,strn,len(strn),.false.))
            if (strn(1:6) == 'BASIS:' .or. strn == ' ') cycle
            write(ifi,'(a)') trim(strn)
          enddo
   80     j1 = j2+1
        enddo

        call info0(10,0,0,' wrote basis file '//trim(fn))
        goto 10

C --- Printout structural information ---
      elseif (outs(1:4) == 'show')  then

        j1 = 6
        dc2 = outs(j1-1:j1-1)

        call nwordg(outs,1,dc//dc2//' ',1,j1,j2)
        lall = j1>j2
        i = wordsw(outs,dc2,'all',dc//dc2//' ',j1)
        lall = lall .or. i>0

        call dinv33(plat,1,qlat,vol)
        call dpcopy(qlat(1,3),normal,1,3,1/dlength(3,qlat(1,3),1))

C   ... Size of existing lattice
        i = wordsw(outs,dc2,'size',dc//dc2//' ',j1)
        if (i>0 .or. lall) then
          xx = avwsr(plat,alat,vol(1),1)
          h(1) = ddot(3,plat(1,3),1,normal,1)
          write(stdo,351)
  351     format(/t17,'Plat',t55,'Qlat')
          do  i = 1, 3
            call info2(2,0,0,'%3;11,6D     %3;11,6D',plat(1,i),qlat(1,i))
          enddo
          call info5(10,1,0,'   Cell volume %;9F  nspec %i  nbas %i  height %;9F  normal to plane %s,%3,6d',
     .      vol,nspec,nbas,h(1),normal)
        endif

        allocate(plan(nbasp),iplan(nbasp),natpl(nbasp),ips(nbasp),pl(nbasp),pos(3,nbasp),slabll(nbasp))
        call ivset(pl,1,nbasp,0)
        if (size(s_ctrl%pgfsl) > nbasp) call icopy(nbasp,s_ctrl%pgfsl,1,pl,1)
        call sitepack(s_site0,1,nbasp,'spec',1,ips,xx)
        call sitepack(s_site0,1,nbasp,'pos',3,xx,pos)
        forall (is=1:nspec) slabll(is) = s_spec0(is)%name

C   ... Site positions
        i = wordsw(outs,dc2,'pos',dc//dc2//' ',j1)
        if (i>0 .or. lall) then
          write(stdo,357)
  357     format(/' site spec',8x,'pos (Cartesian coordinates)',9x,'pos (multiples of plat)')
C         qlat = (plat+)^-1
          call dinv33(plat,1,qlat,xx)
          do  i = 1, nbasp
            is = ips(i)
            call dpscop(pos,xv,3,3*i-2,1,1d0)
C           posp+ = (plat)^-1 pos+
            call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
            print 345, i, slabll(is), (xv(j),j=1,3), (xv(3+j),j=1,3)
  345       format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
          enddo
        endif

C   ... Projection onto normal
        i = wordsw(outs,dc2,'planes',dc//dc2//' ',j1)
        if (i>0 .or. lall) then
          call gtplan(plat,npadl,npadr,nbasp,pos,normal,ips,pl,slabll,nplane,plan,iplan,natpl)
        endif

C   ... Sphere overlaps
        i = wordsw(outs,dc2,'ovl',dc//dc2//' ',j1)
        if (i>0 .or. lall) then
          s_ctrl%nclass = nspec; s_ctrl%nclasp = nspec
          call ptr_ctrl(s_ctrl,8+1,'rmax',nspec,0,0,xx)
          call ptr_ctrl(s_ctrl,8+1,'ics',nspec,0,0,xx)
          call ptr_ctrl(s_ctrl,8+1,'ipc',nbas,0,0,xx)
          if (associated(s_ctrl%clabl)) deallocate(s_ctrl%clabl)
          allocate(s_ctrl%clabl(nspec))
          do  j = 1, nspec
            s_ctrl%rmax(j) = s_spec0(j)%rmt
            s_ctrl%ics(j) = j
            s_ctrl%clabl(j) = slabll(j)
          enddo
          do  ib = 1, nbas
            is = s_site0(ib)%spec
            s_ctrl%ipc(ib) = is
          enddo
          s_str%rmax = 0; s_str%mxnbr = 0
          call lmaux('stack',s_bz,s_ctrl,s_ham,s_lat,s_mix,s_pot,
     .      s_str,s_spec0,s_site0,' ',slabll,1)
        endif

        deallocate(plan,iplan,natpl,ips,pl,pos,slabll)

        goto 10

C --- Write to positions file ---
      elseif (outs(1:4) == 'wpos') then

C   ... Write site positions into positions file
        j1 = 6
        dc2 = outs(j1-1:j1-1)
        i = wordsw(outs,dc2,'fn=','',j1)
        if (i > 0) then
          call nwordg(outs,0,dc2,1,j1,j2)
          fn = outs(j1:j2)
        else
          fn = 'pos'
        endif
        allocate(pos(3,nbas))
        call sitepack(s_site0,1,nbas,'pos',3,xx,pos)
        call iopos(.true.,0,trim(fn),nbas,pos,s_site0)

        goto 10

C --- Put sites in first quadrant ---
      elseif (outs(1:13) == 'shorten:quad1') then

        call info0(30,0,0,' shifting sites to first quadrant')

        call dinv33(plat,1,qlat,xx)
        do  ib = 1, nbas
C         posp = posc (plat+)^-1  and  posp+ = (plat)^-1 posc+
          call dgemm('T','N',3,1,3,1d0,qlat,3,s_site0(ib)%pos,3,0d0,xv,3)
          ltmp = .true.
          do  while (ltmp)
            ltmp = xv(1)<-fptol .or. xv(2)<-fptol .or. xv(3)<-fptol
            if (xv(1) < -fptol) xv(1) = xv(1)+1
            if (xv(2) < -fptol) xv(2) = xv(2)+1
            if (xv(3) < -fptol) xv(3) = xv(3)+1
          enddo
C         posc = posp (plat+) and  posc+ = (plat) posp+
          call dgemm('N','N',3,1,3,1d0,plat,3,xv,3,0d0,s_site0(ib)%pos,3)
        enddo

        goto 10

C --- Write superlattice file ---
      elseif (outs(1:5) == 'wsite')  then

        fn = 'sites' ! Default file name
        lio = 1000*(1+2+4+8+16*0+32) + 4+1
        dc2 = outs(6:6)
        j1 = 7
        if (dc2 == 'x') then
          lio = 1000*(1+2+4+8+16*0+32) + 4+1 + 10
          dc2 = outs(7:7)
          j1 = 8
        endif

C   ... Shift sites to first quadrant ... move to independent instruction
C        i = wordsw(outs,dc2,'quad1','',j1)
C        if (i /= 0) then
C
C          call info0(30,0,0,' ... shifting sites to first quadrant')
C
C          call dinv33(plx,1,qlat,xx)
C          do  ib = 1, nbas
CC           posp = posc (plat+)^-1  and  posp+ = (plat)^-1 posc+
C            call dgemm('T','N',3,1,3,1d0,qlat,3,s_site0(ib)%pos,3,0d0,xv,3)
C            ltmp = .true.
C            do  while (ltmp)
C              ltmp = xv(1)<-fptol .or. xv(2)<-fptol .or. xv(3)<-fptol
C              if (xv(1) < -fptol) xv(1) = xv(1)+1
C              if (xv(2) < -fptol) xv(2) = xv(2)+1
C              if (xv(3) < -fptol) xv(3) = xv(3)+1
C            enddo
CC           posc = posp (plat+) and  posc+ = (plat) posp+
C            call dgemm('N','N',3,1,3,1d0,plx,3,xv,3,0d0,s_site0(ib)%pos,3)
C          enddo
C        endif

C   ... Short output
        i = wordsw(outs,dc2,'short','',j1)
        if (i /= 0) then
          lio = lio - 1000*iand(lio/1000,32)
        endif

C   ... Set output site file name
        i = wordsw(outs,dc2,'fn','= ',j1)
        if (i /= 0) then
          j1 = j1+1
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          fn = outs(j1:j2)
        endif

C   ... Undo padding in layer case
        if (s_ctrl%lpgf(1) > 0) then
          plat(:,3) = plat(:,3) - 2*platl(:,3) - 2*platr(:,3)
        endif

        if (iosits(lio,3d0,0,fn,ifi,slabl0,alat,plat,nbas,nspec,
     .    s_spec0,s_site0) < 0) call rx('failed to write ssite')

C   ... Write restart file
        i = wordsw(outs,dc2,'rsasa','',j1)
        if (i /= 0) then
C     ... Supplied file name
          fn = 'rsta'
          if (outs(j1:j1) == '=') then
            j1 = j1+1
            call nwordg(outs,0,dc//dc2//' ',1,j1,j2)
            fn = outs(j1:j2)
          endif
          jfi = fopng(trim(fn),-1,0)
          nbasp = nbas+npadl+npadr
          allocate(qnus(3,n0,nsp,nbasp))
          do  ib = 1, nbasp
            call dcopy(3*n0*nsp,s_site0(ib)%qnu,1,qnus(1,1,1,ib),1)
            is = s_site0(ib)%spec
            s_spec0(is)%lmxb = 0; s_spec0(is)%kmxt = 0
          enddo
          nullify(s_pot%v0)
          i = iorsa(0,s_ctrl,s_site0,s_spec0,s_lat,s_pot,s_bz,'made by lmscell',nbasp,nbasp,nspec,qnus,0,.false.,-jfi)

        endif

        call rx0('site file '//trim(fn)//' written')

C --- Quit ---
      elseif (outs(1:5) == 'q')  then

        call rx0('quit stackcel editor')

C --- Help ---
      elseif (outs == '?') then
        print 310
        print 311
        print 312
        print 313
        print 314
        print 315
        print 316
        print 317
        print 318
        print 319
        print 320
        print 321
        print 322

  310   format(/
     .    ' Each instruction can accept options that modify its function.'/
     .    ' Options are separated by a delimiter (first char after keyword),'/
     .    " e.g a space ' ', or '@' as assumed below."/
     .    ' Select one of the following instructions:'//
     .    t4,'file@filnam[@options]'/
     .    t6,"read site file 'filnam' and add to existing structure"/
     .    t6,'filnam is full file name, including extension'/
     .    t6,'Options'/
     .    t8,'@insert=#',t25,'insert new sites before site #.  Sites below # remain as is'/
     .    t8,'@scale=#',t25,'scale the lattice constant of the insertion lattice by #'/
     .    t8,'@stretch=#',t25,'stretch 3rd lat vec P(3) of the insertion lattice by #'/
     .    t25,'Basis vectors are kept as fixed multiples of P(3)'/
     .    t8,'@dpos=x,y,z',t25,'rigidly shift the basis vectors of the insertion lattice'/
     .    t25,'by vector (x,y,z), expressed in Cartesian coordinates'/
     .    t8,'@dup=#',t25,'repeat insertion # times')

  311   format(/
     .    t4,'show@[options]'/
     .    t6,'shows information about the current structure'/
     .    t6,'Options'/
     .    t8,'@all',t25,'show all of the information below (same as no arguments)'/
     .    t8,'@size',t25,'Size of the structure'/
     .    t8,'@pos',t25,'positions within the structure'/
     .    t8,'@ovl',t25,'sphere overlaps'/
     .    t8,'@planes',t25,'project positions onto axis normal to basal plane')

  312   format(/
     .    t4,'scale@#',t25,'Scales lattice constant by #'/
     .    t4,'pad@#1,#2,#3',t25,'Adds a vector to p3.'/
     .    t4,'stretch@#',t25,'Stretches third lattice vector p3 by #.'/
     .    t25,'Basis vectors are kept as fixed multiples of p3.')

  313   format(/
     .    t4,'newpos|addpos[@options]@dpos-specification'/
     .    t6,'shifts or replaces a subset of site positions by dpos'/
     .    t6,'By default the operation is performed on all sites, but the target list can be restricted'/
     .    t6,'The source for dpos is specified in one of three ways:'/
     .    t8,'@dx=#,#,#',t20,'uniform dpos for all sites, Cartesian coordinates'/
     .    t8,'@dp=#,#,#',t20,'uniform dpos for all sites, multiples of plat'/
     .    t8,'@fn=filnam',t20, 'dpos is read from positions file ''filnam'' (source)'/
     .    t20,'The source list cannot be smaller than the target list'/
     .    t20,'By default elements 1,2,3 ..., nbas of the source are applied to the target'/
     .    t6,'Options:'/
     .    t8,'@targ=lst',t20,'Restrict the operation to subset of target sites given by lst.'/
     .    t8,'@src=slst',t20,'(use with @fn=) cull a list of source sites to be copied to target.'/
     .    t20,'The source becomes the list of positions defined by slst'/
     .    t8,'@wdx=nwfile',t20,'writes shifts to a positions file nwfile.ext')

  314   format(/
     .    t4,'cmppos[@options]@fn=file'/
     .    t6,'reads positions from file, and displays shifts relative to current positions'/
     .    t6,'Options:'/
     .    t8,'@shorten',t20,'shortens each shifted site by adding lattice vectors'/
     .    t8,'@targ=lst',t20,'Restrict the operation to subset of target sites given by lst.'/
     .    t8,'@wdx=nwfile',t20,'writes shifts to positions file nwfile.ext')

  315   format(/
     .    t4,'sort@expr1[@expr2][@expr3] | sort@targ=list'/
     .    t6,'reorders sites by sorting them with expressions of site data, or by specifying the permutations directly'/
     .    t6,'Variables that can be used in expressions:'/
     .    t8,'is,ib',t25,'site and species index'/
     .    t8,'x1,x2,x3',t25,'site positions, Cartesian coordinates'/
     .    t8,'p1,p2,p3',t25,'site positions, projections onto plat'/
     .    t8,'h1,h2,h3',t25,'site positions, projections onto qlat'/
     .    t8,'Example',t25,
     .    'order sites by distance from the (p1,p2) plane:  sort h3'//
     .    t4,'shorten:quad1'/
     .    t6,'shifts all site positions to lie in the first quadrant')

  316   format(/
     .    t4,'newspec|newspec=t|newspec=f'/
     .    t6,'If newspec is T, any sites added will be assigned new species names'/
     .    t6,'newspec toggles the current status'/
     .    t6,'newspec=t sets newspec to T'/
     .    t6,'newspec=f sets newspec to F')

  317   format(/
     .    t4,'rsasa@options'/
     .    t6,'reads ASA-style rsta file (for ASA only). Contents must match current structure.'/
     .    t6,'Options:'/
     .    t6,'fn=filename :  name of restart file'/
     .    t6,'rdim[=a,nr,z,lma,class,ves,rmt] : Read parameters from rsta file.'/
     .    t6,'If parameters are not specified, it reads a,nr,z,lmxa,class.'//
     .    t4,'rmsite@targ=list'/
     .    t6,'Removes a list of sites from structure')

  318   format(/
     .    t4,'cart2lat[@q]@x=#1,#2,#3 | lat2cart[@q]@x=#1,#2,#3'/
     .    t6,'Convert Cartesian coordinate <b>(#1,#2,#3)</b> into multiples of lattice vectors or vice-versa.'/
     .    t6,'Optional @q convertq a q point into multiples of reciprocal lattice vectors.'//
     .    t4,'setpl@range=#'/
     .    t6,'Assigns to each site a principal layer index, determined by the range of the pairs table,'/
     .    t6,'assembled from pairs within # of each other (# is in atomic units)'/
     .    t6,'Principal layers are used by lmpg and this switch works only if PGF_MODE is nonzero.')

  319   format(/
     .    t4,'basp@f1@f2...@fn=filename'/
     .    t6,'(for lmf) assembles a basp file from files f1,f2,...'/
     .    t6,'Result is written to file ''filename'' (no extension is appended).'/
     .    t6,'f1, f2, ... should have the standard format for basp files')

  320   format(/
     .    t4,'wsite@[options]'/
     .    t4,'wsitex@[options]'/
     .    t6,'writes a site file corresponding to the current structure to disk'/
     .    t6,'wsitex writes positions as fractional multiples of lattice vectors'/
     .    t6,'Options'/
     .    t8,'@fn=file',t25,'names site file as file.ext (default is sites.ext)'/
     .    t8,'@short',t25,'write site file in short format'/
     .    t8,'@quad1',t25,'translate sites to first quadrant')

  321   format(/
     .    t4,'wpos[@fn=nwfile]',t25,'writes current positions to positions file nwfile.ext')

  322   format(/
     .    t4,'q',t6,'quit the editor')

        goto 10

C --- Option not recognized ---
      else
        call info0(0,0,0,' stackcel:  "'//trim(outs)//
     .    '" not recognized ... enter "?" for a list of options')
        goto 10

      endif

      lrs = 0  ! No associated rsta file

C --- Add new struture to existing lattice --
   30 continue

C --- Read dimensioning data from the next site file ---
      allocate(slabll(mxspec))
      i = 135002               ! Load species labels; file already open
      nspecl = 0
      rewind ifi
      j = iosite(i,3d0,0,trim(fn),ifi,slabll,alatl,platl,nbasl,nspecl,xx,xx,xx,xx,xx,xx,xx)
      xx = avwsr(platl,alatl,vol(2),1)

C --- Printout ---
C      call info8(20,1,0,' read file '//trim(fn)//
C     .    ' %i sites.  alat=%d%?;n;  scaled to %d;%j;%?;n;  stretch=%d;%j;%?;n;  %s,dpos=%3d;%j;',
C     .    nbasl,alatl,isw(scale /= NULLI),scale*alatl,isw(stretch /= NULLI),stretch,
C     .    isw(dpos(1) /= NULLI),dpos)
      call info8(20,1,0,' read file '//trim(fn)//
     .  ' %i sites.  alat=%d%?;n;  scaled to %d;%j;%?;n;  stretch=%d;%j;%?;n;  pad=%d;%j;',
     .  nbasl,alatl,isw(scale /= NULLI),scale*alatl,isw(stretch /= NULLI),stretch,
     .  isw(pad(1) /= NULLI),pad(1))
      if (scale /= NULLI) alatl = alatl*scale
      do  i = 1, 3
        xv(i) = 180/pi*dotprd(3,platl(1,i),1,plat(1,i),1)
        xv(3+i) = (alatl*dlength(3,platl(1,i),1))/(alat*dlength(3,plat(1,i),1))
        call info5(1,0,0,' plat(%i) :  angle to reference %;5F  scale %;9F',i,xv(i),xv(3+i),0,0)
        if (i<3 .and. abs(xv(i))>tol) call logwarn(2,'%12fwarning!  angle mismatch')
        if (i<3 .and. abs(xv(3+i)-1)>tol) call logwarn(2,'%12fwarning!  scale mismatch')
      enddo

      if (dabs(xv(1))+dabs(xv(2)) > tol .or. abs(xv(3+1)-1)>tol .or. abs(xv(3+2)-1)>tol) then
        call info0(2,0,0,' ... nothing done')
        call fclose(ifi)
        goto 10
      endif

C --- Read data from site file ---
      allocate(s_sitel(nbasl),s_specl(nspecl))
      forall (is=1:nspecl) s_specl(is)%name = slabll(is)
      j = iosits(40072,3d0,0,fn,ifi,slabll,alatl,platl,nbasl,nspecl,s_specl,s_sitel)

C --- Stretch lattice ---
      if (stretch /= NULLI) then
C       Convert positions to multiples of platl
        call dinv33(platl,0,plati,xx)
        allocate(posl(3,nbasl),pos(3,nbasl))
        call sitepack(s_sitel,1,nbasl,'pos',3,xx,posl)
        call dgemm('N','N',3,nbasl,3,1d0,plati,3,posl,3,0d0,pos,3)
C       Stretch platl
        call dscal(3,stretch,platl(1,3),1)
C       Convert positions back to Cartesian coordinates
        call dgemm('N','N',3,nbasl,3,1d0,platl,3,pos,3,0d0,posl,3)
        call sitepack(s_sitel,1,nbasl,'-pos',3,xx,posl)
        deallocate(posl,pos)
      endif

C ... Printout of insertion volume, height
      call dinv33(platl,1,qlat,vol(3))
      call dpcopy(qlat(1,3),normal,1,3,1/dlength(3,qlat(1,3),1))
      xx = avwsr(platl,alatl,vol(3),1)
      h(2) = ddot(3,platl(1,3),1,normal,1)
      call info5(20,0,0,' Insertion volume %;9F (initially %;9F)  height %;9F  normal %s,%3,6d',
     .  vol(3),vol(2),h(2),normal,5)

C --- Make superlattice plx ---
      if (insert == NULLI .or. .true.) then
        if (pad(1) /= NULLI) plat(:,3) = plat(:,3)*(1+pad(1))
        call dcopy(9,plat,1,plx,1)
        plx(:,3) = plat(:,3) + platl(:,3)
        call dinv33(plx,1,qlat,vol(1))
        call dpcopy(qlat(1,3),normal,1,3,1/dlength(3,qlat(1,3),1))
        xx = avwsr(plx,alat,vol,1)
        h(3) = ddot(3,plx(1,3),1,normal,1)
        call info5(20,0,0,' Combined  volume=%;9F  height %;9F  normal %s,%3d',vol,h(3),normal,4,5)
        call info0(20,0,0,'%16pPlat%54pQlat')
        do  i = 1, 3
          call info2(2,0,0,'%3;11,6D     %3;11,6D',plx(1,i),qlat(1,i))
        enddo
      else
        call dcopy(9,plat,1,plx,1)
        call info2(2,0,0,' insert => no change in volume',insert,2)
      endif

C --- Create superlattice s_ssite ---
C ... Count number of species in joint lattice
      nlst = 0; allocate(ilsts(nspecl)); call ivset(ilsts,1,nspecl,0)
      allocate(ips(nbas0)); call sitepack(s_site0,1,nbas0,'spec',1,ips,xx)
      do  i = 1, nspecl
        if (newspec) call uniqclabl(1,1,nbas0,ips,2,slabl0,slabll(i))
C        if (newspec) then
C          do  n = 1, 999
C            call clabel0(slabll(i),1,n,clabl) ! Make a new label
C            do  ib = 1, nbas0  ! Loop over sites prior to insert instruction
C              j = s_site0(ib)%spec
C              if (slabl0(j) == clabl) exit ! Match found
C            enddo
C            if (slabl0(j) /= clabl) exit ! If no match, we are done
C          enddo
C          slabll(i) = clabl
C       endif
C       If new species, append to existing list
        do  j = 1, nspec
          if (slabl0(j) == slabll(i)) exit
          if (j == nspec) then
            nlst = nlst+1
            ilsts(nlst) = i
          endif
        enddo
      enddo
      forall (is=1:nspecl) s_specl(is)%name = slabll(is)

! rmt (and probably other properties) are not initialised for new species. fixes test lmscell_6 on avx512
      do i = 1, nlst
        s_specl(ilsts(i)) % rmt = 0
        s_specl(ilsts(i)) % z = 0
      end do

      deallocate(ips)

      call info2(20,1,0,' Insertion requires %i new species:',nlst,2)
      do  i = 1, nlst
        call info0(20,0,0,' '//trim(slabll(ilsts(i))))
      enddo

C ... Create species table of merged lattice
      allocate(s_ssite(nbas+nbasl),s_sspec(nspec+nlst))
      nspecx = 0
      call sspecadd(nspecx,s_sspec,nspec,0,[0],s_spec0) ! Copy existing structure
      if (nlst > 0) then
      call sspecadd(nspecx,s_sspec,nspecl,nlst,ilsts,s_specl) ! Add new species to list
      endif
      if (nspecx /= nspec+nlst) call rx('bug in stackcel')

C ... Associate species in s_sitel with merged species table
      do  ib = 1, nbasl
        is = s_sitel(ib)%spec
        do  i = 1, nspecx
          if (s_sspec(i)%name == slabll(is)) then
            s_sitel(ib)%spec = i
            exit
          endif
          if (i == nspecx) call rx('bug in stackcel')
        enddo
      enddo
      deallocate(slabll)

C ... Shift positions in new sites by plat + dpos
      call info0(20,1,0,' Merging basis ... add uniform shift to inserted sites')
      if (dpos(1) == NULLI) dpos = 0  ! dpos from dpos=
      call info2(20,0,0,' from plat  %3:1;9F',plat(1,3),2)
      call info2(20,0,0,' from dpos  %3:1;9F',dpos,2)
      shft(1:3) = dpos(1:3)
      if (insert == NULLI) then
        shft(1:3) = dpos(1:3) + plat(:,3) ! Add basic shift : lattice vector of existing cell
      endif
      if (pad(2) /= NULLI) then
        call info2(20,0,0,' from shft  %3:1;9F',plat(:,3)*pad(2),2)
        shft(1:3) = shft(1:3) + plat(:,3)*pad(2)
      endif
      call info2(20,0,0,' net shift  %3:1;9F',shft,2)
      do  ib = 1, nbasl
        s_sitel(ib)%pos(1:3) = s_sitel(ib)%pos(1:3) + shft(1:3)
        s_sitel(ib)%qnu(1,1,1) = NULLI
      enddo

C ... Read info from rsta file
      if (lrs > 0) then
        call rdasapq(outs(lrs:),dc,s_ctrl,s_sitel,s_specl,s_pot,1,nbasl,nbas,errmsg)
        if (errmsg /= ' ') goto 97
      endif

C ... Create s_ssite for merged lattice
      nsitex = 0
      call ssiteadd(nsitex,s_ssite,nbas,0,[0],s_site0) ! Copy existing structure to s_ssite
      call ssiteadd(nsitex,s_ssite,nbasl,0,[0],s_sitel) ! Add new sites to list
      if (nsitex /= nbas+nbasl) call rx('bug in stackcel')

C ... Copy superlattice s_sspec, s_ssite into s_spec0, s_site0; make new slabl0
      deallocate(s_site0,s_spec0,slabl0); allocate(s_site0(nsitex),s_spec0(nspecx),slabl0(nspecx))
      nbas = nsitex; nspec = nspecx; plat = plx
      nspecx = 0; nsitex = 0
      call sspecadd(nspecx,s_spec0,nspec,0,[0],s_sspec)
      call ssiteadd(nsitex,s_site0,nbas,0,[0],s_ssite)
      forall (is=1:nspec) slabl0(is) = s_spec0(is)%name

C ... Make permutation table if an insertion
      if (insert /= NULLI) then
        deallocate(ilsts); allocate(ilsts(nbas),iblst(nbas))
        forall (ib=1:nbas) iblst(ib) = ib
        forall (ib=insert:nbas) iblst(ib) = ib+nbasl
        forall (ib=1:nbasl) iblst(ib+nbas-nbasl) = ib+insert-1
        forall (ib=1:nbas) ilsts(iblst(ib)) = ib  ! siteperm expects reverse order
        call siteperm(nbas,ilsts,s_site0)
        deallocate(iblst)
      endif

      s_ctrl%nbas = nbas
      s_ctrl%nspec = nspec
      call ptr_lat(s_lat,1,'pos',3,nbas,0,0)
      call sitepack(s_site0,1,nbas,'pos',3,xx,s_lat%pos)
      s_lat%plat = plat

C ... Set up for the next site file
      deallocate(s_sitel,s_specl,ilsts,s_ssite,s_sspec)
      ndup = ndup - 1
      if (ndup > 0) goto 30
      call fclose(ifi)
      goto 10

   97 continue
      call info0(2,0,0,' '//trim(errmsg)//' ... nothing done')
      goto 10

   98 continue
      call rx('stackcel : failed to parse '//trim(outs))

      end

      subroutine sspecadd(nspec,s_spec,nspecadd,nlst,ilsts,s_specadd)
C- Append species to existing structure
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nspec    :number of species currently in the table
Ci   nspecadd :number of species to add
Ci   nlst     :Size of ilsts.  If 0, include all species
Ci   ilsts    :sorted list of species in s_specadd to add.
Ci            :Not used if nlst is zero.
Cio  s_specadd
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   03 Jul 17
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nspec,nspecadd,nlst,ilsts(nlst)
C ... For structures
!       include 'structures.h'
      type(str_spec)::  s_spec(nspec+nspecadd),s_specadd(nspecadd)
C ... Local parameters
      integer nspecx,is,ilst

      nspecx = 0
      ilst = 0
      do  is = 1, nspecadd

C       Skip this species if list is present and species is not in list
        if (nlst > 0) then
          do while (ilst < nlst)
            if (ilsts(1+ilst) == is) exit
            ilst = ilst+1
          end do
          if (ilst == nlst) cycle
        endif

C   ... Append this species to s_spec
        nspecx = nspecx+1
        s_spec(nspec+nspecx)%name  = s_specadd(is)%name
        s_spec(nspec+nspecx)%a     = s_specadd(is)%a
        s_spec(nspec+nspecx)%lmxa  = s_specadd(is)%lmxa
        s_spec(nspec+nspecx)%lmxl  = s_specadd(is)%lmxl
        s_spec(nspec+nspecx)%lfoca = s_specadd(is)%lfoca
        s_spec(nspec+nspecx)%nr    = s_specadd(is)%nr
        s_spec(nspec+nspecx)%z     = s_specadd(is)%z
        s_spec(nspec+nspecx)%rmt   = s_specadd(is)%rmt
        s_spec(nspec+nspecx)%hcr   = s_specadd(is)%hcr
        s_spec(nspec+nspecx)%idmod = s_specadd(is)%idmod

        s_spec(nspecx)%rsma       = 0
        s_spec(nspecx)%rsmv       = 0
        s_spec(is)%qc             = 0
        s_spec(is)%pz             = 0
        s_spec(is)%rfoca          = 0
        s_spec(is)%ctail          = 0
        s_spec(is)%etail          = 0
        s_spec(is)%stc            = 0

!       Nothing else needed for now
!       s_spec(nspecx)%rfoca       = s_specadd(is)%rfoca

      enddo
      nspec = nspec+nspecx

      end

      subroutine ssiteadd(nsite,s_site,nsiteadd,nlst,ilstb,s_siteadd)
C- Append new site structure to existing site structure
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsite     :Number of sites already in s_site
Ci   nsiteadd  :Number of sites to add
Ci   nlst      :If > 0, list of sites which must be ... ?
Ci   ilstb     :list of sites for nlst
Cio Structures
Cio  s_site    :site struct to be enlarge
Cio  s_siteadd : enlarge s_site with s_siteadd
Cr Remarks
Cr   s_site must be allocated at least as large as nsite+nsiteadd
Cu Updates
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nsite,nsiteadd,nlst,ilstb(nlst)
C ... For structures
!       include 'structures.h'
      type(str_site)::  s_site(nsite+nsiteadd),s_siteadd(nsiteadd)
C ... Local parameters
      integer nsitex,ib,ilst

      nsitex = 0
      ilst = 0
      do  ib = 1, nsiteadd

C       Skip this site if list is present and site is not in list
        if (nlst > 0) then
   12     continue
          if (ilst >= nlst) cycle
          if (ilstb(1+ilst) < ib) then
            ilst = ilst+1
            goto 12
          endif
          if (ilstb(1+ilst) /= ib) cycle
        endif

C   ... Append contents of this site to s_site
        nsitex = nsitex+1

        s_site(nsite+nsitex)%spec = s_siteadd(ib)%spec
        s_site(nsite+nsitex)%pos = s_siteadd(ib)%pos
        s_site(nsite+nsitex)%class = s_siteadd(ib)%class
        s_site(nsite+nsitex)%force = s_siteadd(ib)%force
        s_site(nsite+nsitex)%vel = s_siteadd(ib)%vel
        s_site(nsite+nsitex)%eula = s_siteadd(ib)%eula
        s_site(nsite+nsitex)%pl = s_siteadd(ib)%pl
        s_site(nsite+nsitex)%relax = s_siteadd(ib)%relax
        s_site(nsite+nsitex)%vshft = s_siteadd(ib)%vshft
        s_site(nsite+nsitex)%pnu = s_siteadd(ib)%pnu
        s_site(nsite+nsitex)%pz  = s_siteadd(ib)%pz
        s_site(nsite+nsitex)%qnu = s_siteadd(ib)%qnu
        s_site(nsite+nsitex)%pz = s_siteadd(ib)%pz
        s_site(nsite+nsitex)%v0 => s_siteadd(ib)%v0
        s_site(nsite+nsitex)%v1 => s_siteadd(ib)%v1

      enddo
      nsite = nsite+nsitex

      end

      subroutine makepl(nbas,npadl,npadr,ntab,iax,pl)
C- Build PL from given neighbor table
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   npadl :number of sites in first PL
Ci   npadr :number of sites in last PL
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Co Outputs
Co   pl    :Principal layer indices
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Jan 18
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, parameter :: niax=10, NULLI = -99999
C ... Local parameters
      integer nbas,npadl,npadr,iax(niax,*),pl(nbas+2*npadl+2*npadr),ntab(nbas+npadl+npadr+1)
      integer npl,ib,jb,ipr,i12(2),nbasp

C ... Build up PL from the left
C     Initial prior layer is left padded layer
      i12(1) = nbas+1; i12(2) = nbas+npadl; npl = 1
      do  ib = 1, nbas
        pl(ib) = NULLI
        do  ipr = ntab(i12(1))+1, ntab(i12(2)+1)
          if (iax(2,ipr) /= ib) cycle
          pl(ib) = npl
          exit
        enddo
        if (pl(ib) == NULLI .and. npl == 1 .and. ib <= npadl) then ! Stick left PL until npadl
          pl(ib) = 1
        elseif (pl(ib) /= NULLI .and. npl == 1 .and. ib > npadl) then ! Stick left PL until npadl
          print *, ib
          call rx('makepl problem with left layer')
        elseif (pl(ib) == NULLI) then ! New PL
          i12(1) = i12(2)+1
          if (npl == 1) i12(1) = 1 ! Special handling if prior PL is L-padded layer
          i12(2) = ib-1
          npl = npl+1
          pl(ib) = npl
        endif
      enddo



C ... Check whether npl should be merged with npl-1
      if (npl > 1) then

C       Propagate all sites ib >= nbas-npadr to PL of nbas-npadr
        npl = pl(nbas-npadr+1)
        forall (jb = nbas-npadr+2:nbas) pl(jb) = npl

        i12(1) = nbas+npadl+1; i12(2) = nbas+npadl+npadr
        do  ib = nbas, 1, -1
          do  ipr = ntab(i12(1))+1, ntab(i12(2)+1) ! Pairs in padded layer
            if (pl(ib) == npl) cycle   ! Already the last layer
            if (iax(2,ipr) /= ib) cycle ! Connective vector does not touch this site
C           Right layer spills into this site ... merge npl with npl-1
            npl = npl-1
            forall (jb = ib+1:nbas) pl(jb) = npl
            cycle
          enddo
        enddo
      endif

C ... Set left right padded and doubly padded pl index to 0 and -1
      nbasp = nbas+npadl+npadr
      forall (ib = nbas+1:nbas+npadl) pl(ib) = 0
      forall (ib = nbasp+1:nbasp+npadl) pl(ib) = -1

C ... Set right padded and doubly padded pl index to npl+1 and npl+2
      forall (ib = nbas+npadl+1:nbasp) pl(ib) = npl+1
      forall (ib = nbasp+npadl+1:nbasp+npadl+npadr) pl(ib) = npl+2

      end
