      subroutine plana(npadl,npadr,nbaspp,slabl,
     .  s_lat,s_spec,s_site,ves,vshft,vconst,vrl,pnu,qnu)
C- Plane analysis for lmto programs
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat plat plat2 avw vol platl platr
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: z qc lmxa rmt coreh coreq
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: getq gtpcor iosits
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: class spec clabel pos vel eula pl relax vshft ndelta
Ci                delta bxc cpawt omg omgn domg gc gcorr j0 pdos rho1
Ci                rho2 rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl
Ci                eqhkl eqkkl sighh sighk sigkk tauhh tauhk taukk pihh
Ci                pihk pikk sighhx sighkx sigkkx tauhhx tauhkx taukkx
Ci                pihhx pihkx pikkx thet v0 v1
Co     Stored:    spec pos vel eula pl relax vshft bxc cpawt omg omgn
Co                domg gc gcorr j0 pdos rho1 rho2 rhoc rho1x rho2x
Co                rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl sighh sighk
Co                sigkk tauhh tauhk taukk pihh pihk pikk sighhx sighkx
Co                sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx thet v0
Co                v1
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: iosits bcast_strx
Ci Inputs
Ci   nbaspp:size of doubly padded basis (layer programs)
Ci          nbaspp = nbas + 2*(nbas(left bulk) + nbas(right bulk))
Ci   slabl :list of species labels.
Ci   ves   :electrostatic potential at MT boundaries
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Co Outputs
Cu Updates
Cu   08 Jul 17 Started reworking stretch; new show and wsite
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   03 Nov 04 update last arg in partok call (need integer)
Cu   31 Oct 03 plana orders planes L-C-R for layer geometry
Cu   01 Mar 02 plana can write to site file
Cu             Altered argument list
Cu   12 Feb 02 Adapted for ASA v6
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer npadl,npadr,nbaspp
      character*8 slabl(*)
      double precision ves(*),pnu(*),qnu(*),vshft(nbaspp),vconst(3),vrl
c     double precision pnu(nl,nsp,1),qnu(3,nl,nsp,1),pp(*),bas(3,1)
C ... For structures
!       include 'structures.h'
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: iplan(:),ilsts(:),iblst(:)
      integer, allocatable :: natpl(:)
      real(8), allocatable :: plan(:)
      integer, allocatable :: iwk(:)
      real(8), allocatable :: bas2(:),pos(:,:),posl(:,:),pos2(:,:)
      real(8), allocatable :: wk(:)
C      real(8), allocatable :: vshfp(:)
C      real(8), allocatable :: vdipl(:)
      real(8), allocatable :: potqz(:)
C ... Local parameters
      logical ltmp
      integer ipc(nbaspp),ips(nbaspp),lmx(nbaspp),pl(nbaspp)
      double precision qc(nbaspp),qt(nbaspp),dq(nbaspp)
      double precision z(nbaspp),bas(3,nbaspp),wsr(nbaspp)
      integer nplane,mcont2,recln0,i,iv(10),iarg(10),j,j1,j2,js1,js2,
     .  k,m,n,ifi,ib,ic,is,l,lio,stdo,ip,nlsts,nblst
      integer nbas,nclspp,nl,nsp,mxint,nbasp,nclasp,nspec
      double precision xv(10),darg(10),cellen,ddot,avw,carea,vol,hhat(3),hvc
      double precision alat,plat(3,3),plat0(3,3),qlat(3,3),normal(3),dpos(3)
C     double precision hdotb
      double precision cphi,sphi,ctheta,stheta,phi,xx,stretch,
     .  tmp(3,3),tmp2(3,3),tmp3(3,3),pi,vol0,glat(3,3),r(5),r2(5),
     .  platl(3,3),plati(3,3),platr(3,3)
      character*8 nam, dc*1, dc2*1, sopts*256, outs*256, fn*120, strn*120
C     character dc*1, fn*120, fn2*120, outs*256, strn*80
      parameter (mcont2=17+2-1, recln0=72)
      logical sw,lsopts
      character*8 clabl(nbaspp)
      character*80 carg(10)
      procedure(logical) :: a2bin,cmdopt
      procedure(integer) :: nglob,iosits,i1mach,a2vec,wordsw,iosite,lgunit,iclbas,mkilsd
      procedure(integer) :: fext,fopn,fopnx,fopna
      procedure(real(8)) :: dlength,avwsr

C     data ncont2 /-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/

      xv = 0

      stdo = nglob('stdo')
      sopts = ' '
      if (cmdopt('--pledit',8,0,outs)) then
        sopts = outs(9:)
      endif

      nl = nglob('nl')
      nsp = nglob('nsp')
      nbas = nglob('nbas')
      nspec = nglob('nspec')
      nbasp  = (nbas+nbaspp)/2
      call dpzero(qc,nbaspp); call dpzero(qt,nbaspp); call dpzero(dq,nbaspp)

C     Make ipc, ips, clabl
      do  ib = 1, nbaspp
        ic = s_site(ib)%class
        is = s_site(ib)%spec
        ipc(ib) = ic
        ips(ib) = is
        pl(ib) = s_site(ib)%pl
        clabl(ic) = s_site(ib)%clabel
        z(ic) = s_spec(is)%z
        qc(ic) = s_spec(is)%qc
        lmx(ic) = s_spec(is)%lmxa
        wsr(ic) = s_spec(is)%rmt
      enddo
      nclasp = mxint(nbasp,ipc)
      nclspp = mxint(nbaspp,ipc)
C     Unpack basis
      call sitepack(s_site,1,nbaspp,'pos',3,xx,bas)
      alat = s_lat%alat
      plat = s_lat%plat

C     Unit vector normal to the plane
      call dinv33(s_lat%plat,1,qlat,vol)
      do  i = 1, 3
        normal(i) = qlat(i,3)/dlength(3,qlat(1,3),1)
      enddo

      avw = s_lat%avw
      vol = s_lat%vol
      hvc = dsqrt(ddot(3,normal,1,normal,1)) ! should be 1 now
      if (pnu(1) /= 0) then
      call getq(nsp,nl,lmx,nclspp,z,pnu,qnu,0,s_spec,qc,qt,dq)
      endif

      call dcopy(9,plat,1,plat0,1)
      rewind (fopn('PLAN'))
      allocate(plan(nbasp),iplan(nbasp),natpl(nbasp))
      call gtplan(plat,npadl,npadr,nbasp,bas,normal,ipc,pl,clabl,
     .  nplane,plan,iplan,natpl)

C --- Input from from stdin ---
  120 continue

      dc = sopts(1:1)
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the layer file editor.  Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the layer file editor.  Enter ''?'' to see options.')
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

C ... re-entry if command already given
   19 continue

C --- Parse and execute the next command ---
      if (.false.) then

C --- Null input ---
      elseif (outs == ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',' ''?'' to see menu')
        goto 10

C --- Read self-energy from file ---
C      elseif (outs(1:5) == 'rplan ') then
C        call rx('need pnu qnu z wsr idxdn for ioplan')
C        ifi = fopn('PLAN')
C        call ioplan(nbas,nl,nsp,nclspp,plat,
C     .    clabl,z,wsr,pnu,qnu,idxdn,ipc,bas,iplan,npl,ifi)
C
CC     ... Begin loop that makes various adjustments
C          do  i = 1, 3
C          hhat(i) = plat(i,3) / hvc
C          enddo
C          do  i = 1, npl
C          ib = i+nbasp
C          ic = i+nclspp
C          ipc(ib) = ic
C
CC      ... Check for replication of class names
Cc          do  43  j = 1, nclspp
C            do  j = 1, ic-1
C            if (clabl(j) /= clabl(ic)) cycle
C              print 2, clabl(ic), j
C    2         format('   Class ',a4,' already named in class',i3)
C              write(*,3,advance='no')
C    3         format('   Confirm name, or type new name: ')
C              read(*,'(a4)') nam
C            if (nam /= ' ')  clabl(ic) = nam
C          enddo
C
CC     ... Adjust bas to proper projection onto h
C          hdotb = ddot(3,bas(1,ib),1,hhat,1)
C          bas(1,ib) = bas(1,ib) + hhat(1)*(darg(1) - hdotb)
C          bas(2,ib) = bas(2,ib) + hhat(2)*(darg(1) - hdotb)
C          bas(3,ib) = bas(3,ib) + hhat(3)*(darg(1) - hdotb)
C            print 4,clabl(ic),(bas(j,ib),j=1,3)
C    4       format(' Atom ',a4,' inserted at:',3F12.5)
C          enddo
C
C        nbas = nbas+npl
C        nbasp = nbasp+npl
C        nclasp = nclasp+npl
C        nclspp = nclspp+npl
C
CC   ... Remap class table, finding equivalent names
C        allocate(iwk(nbasp))
C        call fixcl(nbasp,nclspp,ipc,clabl,iwk)
C        deallocate(iwk)
C
CC   ... Make new plane table
C        deallocate(plan,iplan,natpl)
C        allocate(plan(nbasp),iplan(nbasp),natpl(nbasp))
C        call gtplan(npadl,npadr,nbasp,bas,normal,ipc,pl,clabl,
C     .    nplane,plan,iplan,natpl)
C        goto 10

C --- Write a single plane to plan file  ---
C      elseif (outs(1:6) == 'wplan ') then
C        call rx('need pnu qnu z wsr idxdn for ioplan')
C        call ioplan(nbas,nl,nsp,nclspp,plat,
C     .    clabl,z,wsr,
C     .    pnu,qnu,
C     .    idxdn,
C     .    ipc,bas,iplan,iarg,-fopn('PLAN'))
C        goto 10

C --- Skip over planes ---
C      elseif (outs(1:6) == 'wplan ') then
C        call rx('need pnu qnu z wsr for ioplan')
C        ifi = fopn('PLAN')
C        rewind ifi
C        do  i = 1, iarg(1)-1
C          print *, '... skipping past plane ', i
C          call ioplan(nbas,nl,nsp,nclspp,plat,
C     .      clabl,z,wsr,
C     .      pnu,qnu,
C     .      idxdn,
C     .      ipc,bas,iplan,npl,ifi)
C        enddo
C        goto 10

C --- Change extension to name  ---
C      elseif (outs(1:4) == 'ext ') then
C        call fclose(fopn('PLAN'))
C        ifi = fext(carg)
C        ifi = fopn('PLAN')
C        goto 10

C --- Show ---
      elseif (outs(1:4) == 'show')  then
        dc2 = outs(5:5)

        j1 = 6
        call nwordg(outs,1,dc//dc2//' ',1,j1,j2)
        ltmp = j1>j2
        i = wordsw(outs,dc2,'all',dc//dc2//' ',j1)
        ltmp = ltmp .or. i>0

C   ... Size of existing lattice
        i = wordsw(outs,dc2,'size',dc//dc2//' ',j1)
        if (i>0 .or. ltmp) then
          call dinv33(plat,1,qlat,vol)
C         call dpcopy(qlat(1,3),normal,1,3,1/dlength(3,qlat(1,3),1))
          xx = avwsr(plat,alat,vol,1)
          xv(1) = ddot(3,plat(1,3),1,normal,1)
          write(stdo,351)
  351     format(/t17,'Plat',t55,'Qlat')
          do  i = 1, 3
            call info2(2,0,0,'%3;11,6D     %3;11,6D',plat(1,i),qlat(1,i))
          enddo
          call info5(10,1,0,'   Cell volume %;9F  nspec %i  nbas %i  height %;9F  normal to plane %s,%3,6d',
     .      vol,nspec,nbas,xv(1),normal)
        endif

C   ... Site positions
        i = wordsw(outs,dc2,'pos',dc//dc2//' ',j1)
        if (i>0 .or. ltmp) then
          write(stdo,357)
  357     format(/' site spec',8x,'pos (Cartesian coordinates)',9x,'pos (multiples of plat)')
C         qlat = (plat+)^-1
          call dinv33(plat,1,qlat,xx)
          do  i = 1, nbas
            is = ips(i)
            call dpscop(bas,xv,3,3*i-2,1,1d0)
C           posp+ = (plat)^-1 pos+
            call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
            print 345, i, slabl(is), (xv(j),j=1,3), (xv(3+j),j=1,3)
  345       format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
          enddo
        endif

C   ... Projection onto normal
        i = wordsw(outs,dc2,'planes',dc//dc2//' ',j1)
        if (i>0 .or. ltmp) then
          call gtplan(plat,0,0,nbas,bas,normal,ips,pl,slabl,nplane,plan,iplan,natpl)
        endif

        goto 10

C --- Replace or shift a subset of current site positions ---
C     newpos|addpos[@targ=lst][@src=list]  @dp=#,#,# | @dx=#,#,# | @fn=filename
      elseif (outs(1:6) == 'newpos' .or. outs(1:6) == 'addpos') then

        j1 = 8
        dc2 = outs(j1-1:j1-1)
        call dinv33(plat,1,qlat,xx)

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
          call info2(10,1,0,' shift :%3;10,6D (Cart) =%3;10,6D (plat) :',xv(4),xv(7))
          allocate(posl(3,nlsts)); call dpzero(posl,size(posl))
          forall (i=1:nlsts) posl(:,i) = xv(4:7)
        else ! Read shifts from file
          i = wordsw(outs,dc2,'fn=','',j1)
          if (i == 0) then
            call info0(10,0,0," no shift specified ... nothing done")
            goto 10
          endif
          call nwordg(outs,0,dc2,1,j1,j2)
          fn = outs(j1:j2)
          ifi = fopnx(fn,70,-1,-1)
          if (ifi /= 1) then
            call info0(10,0,0," missing file '"//trim(fn)//"' ... nothing done")
            goto 10
          endif
C         Extract size of positions file; read into pos2
          call info0(2,0,-1,' read positions from file '//trim(fn)//'.')
          ifi = fopna(fn,-1,0)
          n = 0
          call ioposlu(.false.,0,ifi,n,bas,s_site)
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
          allocate(posl(3,nlsts)); call dpzero(posl,size(posl))
          forall (i=1:nlsts) posl(:,i) = pos2(:,iblst(i))
          deallocate(iblst,pos2)
          call fclr(' ',ifi)
        endif
C       At end of this block, nlsts and posl(1:nlsts) shall be avalailable

C   ... For each element in list, replace base or add posl to bas
        ltmp = outs(1:6) == 'newpos'  ! Copy, rather than add
        write(stdo,357)
        do  i = 1, nlsts
          ib = ilsts(i)
          is = ips(ib)
          xv(1:3) = bas(1:3,ib) + posl(:,i)
          if (ltmp) xv(1:3) = posl(:,i)
C         posp+ = (plat)^-1 pos+
          call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
          print 345, ib, slabl(is), (xv(j),j=1,6)
          bas(1:3,ib) = xv(1:3)
        enddo

        deallocate(posl,ilsts)
        goto 10

C --- Compare current site positions to file ---
C     cmppos[@shorten][@wdx=filename]@fn=filename
      elseif (outs(1:6) == 'cmppos')  then

        j1 = 8
        dc2 = outs(j1-1:j1-1)

C   ... Read file data into pos
        i = wordsw(outs,dc2,'fn=','',j1)
        if (i == 0) then
          call info0(10,0,0," no file specified ... nothing done")
          goto 10
        endif
        call nwordg(outs,0,dc2,1,j1,j2)
        fn = outs(j1:j2)
        ifi = fopnx(fn,70,-1,-1)
        if (ifi /= 1) then
          call info0(10,0,0," missing file '"//trim(fn)//"' ... nothing done")
          goto 10
        endif
        ifi = fopna(fn,-1,0)
        allocate(pos(3,nbas))
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
          dpos(:) = pos(:,ib)-bas(:,ib)
          if (ltmp) call shorps(1,plat,(/72,2,2/),dpos,dpos)
          call info5(2,0,0,'%3;11,6D   %3;11,6D   %3;11,6D',
     .      bas(1,ib),pos(1,ib),dpos,4,5)
C         bas(:,ib) = bas(:,ib) + dpos(:)
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

        deallocate(pos,ilsts)
        goto 10

C --- Stretch plat ---
      elseif (outs(1:6) == 'normal') then
        i = wordsw(outs,'','normal','= ',j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc,3,3,3,iv,normal) /= 3) goto 98
        endif
        call dscal(3,1/dlength(3,normal,1),normal,1)
        goto 10

C --- Stretch plat ---
      elseif (outs(1:7) == 'stretch') then
        i = wordsw(outs,'','stretch','= ',j1)
        if (i /= 0) then
          if (a2vec(outs,len_trim(outs),j1,4,', '//dc,3,3,1,iv,stretch) /= 1) goto 98
        endif

C       Convert positions to multiples of plat
        call dinv33(plat,0,plati,xx)
        allocate(pos(3,nbas))
        call dgemm('N','N',3,nbas,3,1d0,plati,3,bas,3,0d0,pos,3)
C       Stretch plat
        call dscal(3,stretch,plat(1,3),1)
C       Convert positions back to Cartesian coordinates
        call dgemm('N','N',3,nbas,3,1d0,plat,3,pos,3,0d0,bas,3)
        deallocate(pos)

C   ... Make new normal vector
        call dinv33(plat,1,qlat,vol)
        do  i = 1, 3
          normal(i) = qlat(i,3)/dlength(3,qlat(1,3),1)
        enddo

        outs = 'show all'
        goto 19

C --- Sum charge in plane ---
      elseif (outs(1:5) == 'qsum ') then
        ip = 5
        xv(1:4) = 0
        ip = a2vec(outs,len(outs),ip,4,' ',1,-2,-4,iarg,xv)
        iarg(1:2) = xv(1:2)
        darg(1:2) = xv(3:4)
        if (iarg(1) /= xv(1) .or. iarg(2) /= xv(2))
     .    call fexit(-1,119,'qsum expects first 2 args to be integer',0)
        allocate(bas2(nbaspp),wk(nbaspp))
        allocate(potqz(nplane))
        cellen = alat*hvc
        carea =  avw
        carea =  vol / cellen
C       call dplmom(nbaspp,s_lat%pos,alat,plat,ipc,qt,0d0,wk,xv(1),xv(2))
        call xxpl(nbasp,npadl,npadr,alat,plat,s_lat%pos,ipc,nplane,plan,iplan,ves,vshft,vrl,iarg,qt,
     .    clabl,cellen,carea,natpl,wsr,darg,potqz)
        deallocate(potqz,wk)
        goto 10

C --- Write superlattice file ---
      elseif (outs(1:5) == 'wsite')  then

        fn = 'sites' ! Default file name
        lio = 15001
        dc2 = outs(6:6)
        j1 = 7
        if (dc2 == 'x') then
          lio = 15001 + 10
          dc2 = outs(7:7)
          j1 = 8
        endif

C   ... Short output
        i = wordsw(outs,dc2,'short','',j1)
        if (i /= 0) then
          lio = lio - 1000*iand(lio/1000,32)
        endif

C   ... Set output site file name
        i = wordsw(outs,dc2,'fn','= ',j1)
        if (i /= 0) then
          fn = outs(j1+1:)
        endif

        j1 = iosite(lio,3d0,0,fn,i,slabl,alat,plat,nbas,
     .    nspec,bas,xx,xx,xx,ips,xx,xx)

        call rx0('site file '//trim(fn)//' written')

C --- Slide planes ---
      elseif (outs(1:6) == 'slide ') then
        call rx('slide must parse darg')
        do  i = 1, 3
          hhat(i) = plat(i,3) / hvc
        enddo
        do  ib = 1, nbasp
        j = iplan(ib)
        if (j >= iarg(1) .and. j <= iarg(2)) then
          call daxpy(3,darg(1),hhat,1,bas(1,ib),1)
            print 5,ib,(bas(j,ib),j=1,3)
    5       format(' Sliding site',i3,' to',3F10.5)
        endif
        enddo

C   ... Make new plane table
        deallocate(plan,iplan,natpl)
        allocate(plan(nbasp),iplan(nbasp),natpl(nbasp))
        call gtplan(npadl,npadr,nbasp,bas,normal,ipc,pl,clabl,
     .    nplane,plan,iplan,natpl)
        goto 10

C --- Rename classes ---
      elseif (outs(1:5) == 'cnam ') then
        allocate(iwk(nbasp))
        call fixcl(nbaspp,nclspp,ipc,clabl,iwk)
        deallocate(iwk)
        do  ic = 1, nclspp
          print 6,clabl(ic),ic,nclspp
    6     format(' Set name of class ',a4,', number',i3,' of',i3,
     .      ' to:')
          read(*,'(a4)') nam
          if (nam /= ' ')  clabl(ic) = nam
        enddo
        goto 10

C --- Shift ves in planes ---
      elseif (outs(1:7) == 'vshift ') then
        call rx('vshfit must parse iarg, darg')
        do  ic = 1, nclspp
        ib = iclbas(ic,ipc,nbaspp)
        if (ib == 0) cycle
        j = iplan(ib)
        if (j >= iarg(1) .and. j <= iarg(2)) then
            print 7,clabl(ic),ves(ic),ves(ic)+darg(1)
    7       format(' Shift ves of class ',a4,' from',f10.5,' to',f10.5)
          ves(ic) = ves(ic)+darg(1)
        endif
        enddo
      goto 10

C --- Rotate plat ---
      elseif (outs(1:4) == 'rot ') then
        call rx('rot must parse iarg, darg')
      pi = 4*datan(1d0)
      darg(1) = darg(1)*pi
C        do  i = 1, 3
C        hhat(i) = planvc(i,3)/hvc
C        enddo
c      hhat(1) = .6123724
c      hhat(2) = .6123724
c      hhat(3) = .5
      if (dabs(hhat(2))+dabs(hhat(2)) < 1d-8) then
        phi = 0
      else
        phi = datan2(hhat(2),hhat(1))
      endif
      cphi = dcos(phi)
      sphi = dsin(phi)
      ctheta = hhat(3)
      stheta = dsqrt(1-ctheta**2)
      tmp(1,1) = cphi*ctheta
      tmp(2,1) = sphi*ctheta
      tmp(3,1) = -stheta
      tmp(1,2) = -sphi
      tmp(2,2) = cphi
      tmp(3,2) = 0
      tmp(1,3) = cphi*stheta
      tmp(2,3) = sphi*stheta
      tmp(3,3) = ctheta
        print 16,tmp
      ctheta = dcos(darg(1))
      stheta = dsin(darg(1))
      call dpzero(tmp2,9)
      tmp2(1,1) = ctheta
      tmp2(2,1) = stheta
      tmp2(1,2) = -stheta
      tmp2(2,2) = ctheta
      tmp2(3,3) = 1
      call dmpy(tmp,3,1,tmp2,3,1,tmp3,3,1,3,3,3)
      call dmpy(tmp3,3,1,tmp,1,3,tmp2,3,1,3,3,3)
        print 8
    8   format(/9x,'ROTATION MATRIX')
        print 16,tmp2
      call dmpy(tmp2,3,1,plat,3,1,tmp,3,1,3,3,3)
        print 14
        print 15,((plat(m,k),m=1,3),(tmp(m,k),m=1,3),k=1,3)
      call dcopy(9,tmp,1,plat,1)
      allocate(bas2(3*nbasp))
      call dmpy(tmp2,3,1,bas,3,1,bas2,3,1,3,nbasp,3)
      call dcopy(3*nbasp,bas2,1,bas,1)
      deallocate(bas2)
      goto 10

C --- Display projection of basis vectors onto planvc ---
      elseif (outs(1:5) == 'proj ') then

        print *, 'proj has not been redesigned yet, sorry'
        goto 10

        call dinv33(plat,0,tmp,vol0)
        print *, ' ... Basis in multiples of padded plat'
        print 17
        do  ib = 1, nbasp
          call dmpy(tmp,3,1,bas(1,ib),3,1,tmp2,3,1,3,1,3)
          print 18,ib,(bas(j,ib),j=1,3),(tmp2(j,1),j=1,3)
        enddo

        call rx('redo the rest of proj')

C        print *, ' ... Repeat projection after subtracting padded Plat(3)'
C        darg(1) = 0
CC   ... Adjust bas to proper projection onto h
C        do  i = 1, 3
C          hhat(i) = planvc(i,3)/hvc
C        enddo
C        print 17
C        do  ib = 1, nbasp
C          hdotb = ddot(3,bas(1,ib),1,hhat,1)
C          tmp(1,1) = bas(1,ib) + hhat(1)*(darg(1) - hdotb)
C          tmp(2,1) = bas(2,ib) + hhat(2)*(darg(1) - hdotb)
C          tmp(3,1) = bas(3,ib) + hhat(3)*(darg(1) - hdotb)
C          print 18,ib,(bas(j,ib),j=1,3),(tmp(j,1),j=1,3)
C        enddo
        goto 10

C --- Atomic positions within range of plane ---
C To circumscribe 3 planes, first 4 classes only, use eg
C pos 6 6 2 const p1<=4*sqrt(2)&p2<=2.01&p3>1&p3<1.5&ic<=4
C and extract with grep pos log.ext | pextract ssssssssss hijb
C Variables are p1,p2,p3,x,y,z,ic,n, ib(first occurrence)
      elseif (outs(1:4) == 'pos ') then
        call rx('pos must parse iarg, darg')
C   ... Shift all basis vectors to first cell
        call dinv33(plat,0,glat,vol0)
        print *, ' ... shifting basis vectors to first cell'
        do  ib = 1, nbasp
          call dmpy(glat,3,1,bas(1,ib),3,1,r,3,1,3,1,3)
          call dpzero(r2,3)
          do  i = 1, 3
            r(i) = -nint(r(i)-.49999d0)
            call daxpy(3,r(i),plat(1,i),1,r2,1)
          enddo
          call daxpy(3,1d0,r2,1,bas(1,ib),1)
          call awrit4('ib=%i %6psh=%3:1,5;5d (%3;0d)%42pnew bas='//
     .      '%3:1,6;6d',carg,80,i1mach(2),ib,r2,r,bas(1,ib))
        enddo
C   ... Generate table of all positions
        l = (iarg(1)+1)*(iarg(2)+1)*(iarg(3)+1)*nbasp
        allocate(bas2(5*l),wk(l))
        call dcopy(3,normal,1,hhat,1)
        call dscal(3,1/hvc,hhat,1)
        l = 0
        do  ib = 1, nbasp
        do  i = 0, iarg(1)
        do  j = 0, iarg(2)
        do  k = 0, iarg(3)
          call dpzero(r(3),3)
          call daxpy(3,dble(i),plat(1,1),1,r(3),1)
          call daxpy(3,dble(j),plat(1,2),1,r(3),1)
          call daxpy(3,dble(k),plat(1,3),1,r(3),1)
          call daxpy(3,1d0,bas(1,ib),1,r(3),1)
          r(1) = ddot(3,hhat,1,r(3),1)
          r(2) = ib
          call dpscop(r,bas2,5,1,5*l+1,1d0)
          if (carg(1) == ' ') then
            l = l+1
          else
            call lodsyv('ib',1,dble(ib),n)
            call lodsyv('ic',1,dble(ipc(ib)),n)
            call lodsyv('x',1,r(3),n)
            call lodsyv('y',1,r(4),n)
            call lodsyv('z',1,r(5),n)
            call lodsyv('n',1,r(1),n)
            call rx('revisit planvc')
C            do  n = 1, 3
C              r2(n) = ddot(3,planvc(1,n),1,r(3),1)/
C     .          dsqrt(ddot(3,planvc(1,n),1,planvc(1,n),1))
C            enddo
            call lodsyv('p1',1,r2(1),n)
            call lodsyv('p2',1,r2(2),n)
            call lodsyv('p3',1,r2(3),n)
            n = 0
            if (.not. a2bin(carg,sw,0,0,' ',n,-1)) cycle
            if (sw) l = l+1
          endif
        enddo
        enddo
        enddo
        enddo
        call dvshel(5,l,bas2,wk,0)
C   ... Write out table of positions within (lo,hi)
        print 9
    9   format(' ic       ib',13x,'projection',26x,'position')
        write (lgunit(2),11)
   11   format(5x,'ic',7x,'ib',13x,'projection',22x,'position')
        m = 0
        do  i = 0, l-1
          call dpscop(bas2,r,5,5*i+1,1,1d0)
          j = nint(r(2))
          k = ipc(j)
          if (darg(4) == -99 .and. darg(5) == -99 .or.
     .       (r(1) >= darg(4) .and. r(1) <= darg(5))) then
            m = m+1
            call rx('revisit planvc')
C            do  n = 1, 3
C              r2(n) = ddot(3,planvc(1,n),1,r(3),1)/
C     .              dsqrt(ddot(3,planvc(1,n),1,planvc(1,n),1))
C            enddo
              print 12,k,clabl(k),j,(r2(n),n=1,3),r(3),r(4),r(5)
   12         format(i3,1x,a4,i4,2(f12.6,2F11.6))
              write (lgunit(2),13) k,clabl(k),j,(r2(n),n=1,3),(r(n),n=3,5)
   13         format('pos',i4,1x,a4,i4,3F11.6,3F10.5)
          endif
        enddo
        k = nbasp*(iarg(1)+1)*(iarg(2)+1)*(iarg(3)+1)
        call awrit3(' pos: %i points written; %i satisfied constraints'//
     .    ' out of %i',carg,80,i1mach(2),m,l,k)
        deallocate(bas2,wk)
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
        call iopos(.true.,0,trim(fn),nbas,bas,s_site)

        goto 10

C --- Write to site file ---
      elseif (outs(1:6) == 'wsite ') then
        call rx('wsite must parse carg')
C     k will be the number of sites to write to site file (normally nbas)
      k = nbas
C     tmp will be the plat written to site file (normally plat)
      call dcopy(9,plat,1,tmp,1)

C     Some padding sites .. possibly alter k, tmp
      if (nbasp > nbas) then
        platl = s_lat%platl
        platr = s_lat%platr
        if (carg(1) == '-ppad') then
          k = nbaspp
          carg(1) = carg(2)
        elseif (carg(1) == '-pad') then
          k = nbasp
          call daxpy(3,-1d0,platl(1,3),1,tmp(1,3),1)
          call daxpy(3,-1d0,platr(1,3),1,tmp(1,3),1)
          carg(1) = carg(2)
        else
          call daxpy(3,-2d0,platl(1,3),1,tmp(1,3),1)
          call daxpy(3,-2d0,platr(1,3),1,tmp(1,3),1)
        endif
      endif
      if (carg(1) == ' ') carg(1) = 'site'

      lio = 1000*(2+4+8+16+32) + 1
      j = iosits(lio,3d0,0,carg,ifi,slabl,alat,tmp,k,nspec,
     .  s_spec,s_site)
      call fclose(ifi)
      goto 10

      elseif (outs(1:2) == '? ') then
        print 310
  310   format(
     .    ' Select one of the following instructions.'/
     .    " Delimit options with a space or another character, e.g. '@'"/
     .   )

        print 311
  311   format(/
     .    t4,'show@[options]'/
     .    t6,'shows information about the current structure'/
     .    t6,'Options:'/
     .    t8,'@all',t20,'show all of the information below (same as no arguments)'/
     .    t8,'@size',t20,'Size of the structure'/
     .    t8,'@pos',t20,'positions within the structure'/
     .    t8,'@plane',t20,'project positions onto axis normal to basal plane')

        print 312
  312   format(/
     .    t4,'newpos|addpos[@options]@dpos-specification'/
     .    t6,'shifts or replaces a subset of site positions by dpos'/
     .    t6,'By default the operation is peformed on all sites, but the target list can be restricted'/
     .    t6,'dpos is specified in one of three ways:'/
     .    t8,'@dx=#,#,#',t20,'constant dpos for all sites, Cartesian coordinates'/
     .    t8,'@dp=#,#,#',t20,'constant dpos for all sites, multiples of plat'/
     .    t8,'@fn=file',t20, 'dpos is read from positions file (source)'/
     .    t20,'Elements 1,2,3 ... of the source are applied to elements 1,2,3 ... of the target'/
     .    t20,'The source list cannot be smaller than the target list'/
     .    t20,'See optional @src=lst to extract a subset of elements from the file'/
     .    t6,'Options:'/
     .    t8,'@targ=lst',t20,'Restrict the operation to subset of target sites given by lst.'/
     .    t8,'@src=slst',t20,'(use with @fn=) cull a list of source sites to be copied to target.'/
     .    t20,'The source becomes the list of positions defined by slst'/
     .    t8,'@wdx=nwfile',t20,'writes shifts to a positions file nwfile.ext')

        print 313
  313   format(/
     .    t4,'cmppos[@options]@fn=file'/
     .    t6,'reads positions from file.ext, displays shifts relative to current positions'/
     .    t6,'Options:'/
     .    t8,'@shorten',t20,'shortens each shifts by adding possible lattice vectors'/
     .    t8,'@targ=lst',t20,'Restrict the operation to subset of target sites given by lst.'/
     .    t8,'@wdx=nwfile',t20,'writes shifts to positions file nwfile.ext')

        print 314
  314   format(/
     .    t4,'stretch h',t20,'stretches plat(3) by h'//
     .    t4,'qsum ip1 ip2 [efcell efblk]'/t20,
     .    'calculates dipole between planes ip2 and ip2, work function and Ef')

        print 315
  315   format(/
     .    t4,'wsite@[options]'/
     .    t4,'wsitex@[options]'/
     .    t6,'writes a site file corresponding to the current structure to disk'/
     .    t6,'wsitex writes positions as fractional multiples of lattice vectors'/
     .    t6,'Options:'/
     .    t8,'@fn=file',t20,'names site file as file.ext (default is sites.ext)'/
     .    t8,'@short',t20,'write site file in short format')

        print 316
  316   format(/
     .    t4,'wpos@fn=file'/
     .    t6,'write site positions to file.ext')

        print 317
  317   format(/
     .    t4,'q',t20,'quit the editor')


C       print 310
C  310   format(
C     .    ' Select one of these options:'/
C     .    t4,'proj',t15,'displays projection of basis vectors into normal plane'/
C     .    t4,'rot th',t15,'rotates plat by th about vc'/
C     .    t4,'qsum ip1 ip2 [efcell efblk]'/t15,
C     .    'calculates dipole between planes ip2 and ip2, work function and Ef'/
C     .    t4,'vshift n1 n2 v'/t15,'shifts potential in planes n1..n2 by v'/
C     .    t4,'name ext',t15,'changes file extension to ext'/
C     .    t4,'cnam',t15,'renames classes'/
C     .    t4,'write',t15,'writes CLASS, SITE and START to log file'/
C     .    t4,'wsite [-pad | -ppad] [filename]'/t15,
C     .    'writes site data to site file.  If not specified, filename is ''site'''/
C     .    t4,'pos [l m n lo hi] [const expr]'/t15,
C     .    'logs atom positions for lattice vectors (0..l,0..m,0..n)'/
C     .    t15,'projecting between (lo,hi) onto normal.'/
C     .    t15,'Logical expr uses vars x,y,z for position'/
C     .    t15,'(eg x<2&y>1) and px,py,pz or n for projection onto normal'/
C     .    t4,'rplan h',t15,'reads next plane from plan file and puts at h'/
C     .    t4,'splan n',t15,'skips to plane n to plan file'/
C     .    t4,'wplan n',t15,'writes plane n to plan file'/
C     .    t4,'rplan (x,y,z) h',t15,'same but puts at h / h.(x,y,z)'/
C     .    t4,'stretch h',t15,'stretches plat(3) by h'/
C     .    t4,'slide n1 n2 (x,y,z) h'/t15,'slides planes n1..n2 by h'/
C     .    t4,'stretch (x,y,z) h'/t15,'stretches plat(3) by h / h.(x,y,z)'/
C     .    t4,'q',t15,'quit editor')
        goto 10

C --- Quit ---
      elseif (outs(1:1) == 'q') then
        call rx0('plana')

      else
        print '(1x,''unrecognized option: '',a10)', outs
        goto 10

C --- Write file (?) ---
C      else
CC     print *, 'implementing write'
C      do  ic = 1, nclspp
C        mapc(ic) = ic
C      enddo
C      call dscal(nclspp,1/avw,wsr,1)
C      call prtbas(nbasp,ipc,clabl,1,nclasp,mapc,-1,1,bas,z,avw,wsr,
C     .            ves,switch,nl,pnu,qnu)
C      call dscal(nclspp,avw,wsr,1)
C      goto 10
      endif
   14 format(/13x,'ORG PLAT',26x,'NEW PLAT')
   15 format(3F10.5,5x,3F10.5)
   16 format(3F10.5)
   17 format(' IB',14x,'BAS',33x,'PROJ')
   18 format(i3,3F10.5,5x,3F10.5)

   98 continue
      call rx('plana : failed to parse '//trim(outs))


      end
      subroutine ioplan(nbas,nl,nsp,nclspp,plat,
     .  clabl,z,wsr,pnu,qnu,idxdn,ipc,bas,iplane,ipl,ifi)
C- I/o of site and class data for a plane of atoms
Ci   ipl write plane i to file (write only)
Co   ipl number of atoms in plane read (read only)

      implicit none

      integer ifi

C Parameters for plane
      integer iplane(2),ipl

C Parameters and for options, structure, basis
      integer nbas,nl,nsp,nclspp
      double precision plat(9)

C Class parameters
      double precision z(*),wsr(*)
      character*8 clabl(1)
      integer idxdn(nl,*)
      double precision pnu(nl,nsp,1),qnu(3,nl,nsp,1)

C Site parameters
      double precision bas(3,nbas)
      integer ipc(nbas)

C local variables
      integer np, jfi, i, ic, j, k, m, ib
      double precision plat0(9)

      if (ifi < 0) goto 1
C --- Read branch ---
      ic = nclspp
      ib = nbas
      jfi = ifi
      read(jfi,*,end=120,err=120) np, plat0
      do  101  i = 1, 6
        if (dabs(plat(i)-plat0(i)) < 1d-4) cycle
        print *, 'ioplan:  plat mismatch in plane file'
  101 continue
      do  110  i = 1, np
        ic = ic+1
        ib = ib+1
        read(jfi,335) m, clabl(ic)
        read(jfi,*) z(ic), wsr(ic), (bas(m,ib), m=1,3)
        read(jfi,*) (idxdn(j,ic), j=1,nl)
        read(jfi,*) ((pnu(j,k,ic), j=1,nl), k=1, nsp)
        read(jfi,*) (((qnu(m,j,k,ic), m=1,3), j=1,nl), k=1, nsp)
        print *, 'read site', i, ', class ', clabl(ic)
  110 continue
      ipl = np
      return
  120 ipl=0
      return

C --- Write branch ---
    1 continue
      jfi = -ifi
C Count number of atoms to write
      np = 0
      do  5  i = 1, nbas
        if (iplane(i) == ipl) np = np+1
    5 continue

      write(jfi,*) np, plat
      do  10  i = 1, nbas
        if (iplane(i) == ipl) then
          ic = ipc(i)
          print *, 'writing site', i, ', class ', clabl(ic)
          write(jfi,335) ic, clabl(ic)
  335     format(i2,1x,a4)
          write(jfi,*) z(ic), wsr(ic), (bas(m,i), m=1,3)
          write(jfi,*) (idxdn(j,ic), j=1,nl)
          write(jfi,*) ((pnu(j,k,ic), j=1,nl), k=1, nsp)
          write(jfi,*) (((qnu(m,j,k,ic), m=1,3), j=1,nl), k=1, nsp)
        endif
   10 continue

      end
      subroutine xxpl(nbasp,npadl,npadr,alat,plat,pos,ipc,nplane,plane,iplane,ves,vshft,vrl,iarg,
     .  qt,clabl,cellen,carea,natpln,rmax,darg,potqz)
      implicit none
      integer nbasp,npadl,npadr,nplane,iarg(2),ipc(nbasp),iplane(nbasp),natpln(nbasp)
      double precision cellen,carea,alat,vrl,darg(2),plat(3,3),
     .  vshft(nbasp),pos(3,nbasp),
     .  plane(nplane),potqz(nplane),
     .  qt(*),ves(*),rmax(*)
      character*8 clabl(nbasp), cc*4
C ... Dynamically allocated arrays
      integer, allocatable :: ibpl(:)
      real(8), allocatable :: vdipib(:)
C Local variables
      integer i,j,ib,ic,imnarg,imxarg,imnpln,stdo
      double precision sumq,field0,pi8,dipole,vintra,rytoev,qdipol,dipmad,
     .  diphar,wfbpln,wfbmad,wfbhar,wfsmad,wfshar,efpln,efmad,efhar,xx,delq,delv
      double precision vmadpl(nplane),vharpl(nplane),vshfpl(nplane),
     .  qplan(nplane),qzpl(nplane),qlat(3,3),vol,b3(3),Lz
!      double precision delv(nplane)
      procedure(integer) :: nglob,iprint,lgunit
      procedure(real(8)) :: dlength
      parameter (rytoev = 13.605826d0)

      if (iarg(1) > 0 .and. iarg(2) > 0) then
        continue
      elseif (iarg(1) > 0 .or. iarg(2) > 0) then
        call rx('need specify two planes, e.g. "qsum 2 4"')
      endif

      pi8 = 32 * datan(1.d0)
      stdo = nglob('stdo')

C --- Calculate planar average of Madelung and Hartree potentials ---
      call dpzero(qplan,nplane)
      call dpzero(vmadpl,nplane)
      call dpzero(vharpl,nplane)
      call dpzero(qzpl,nplane)
      call dpzero(vshfpl,nplane)
      do  i = 1, nbasp
        qplan(iplane(i)) = qplan(iplane(i)) + qt(ipc(i))
        vintra = 2*qt(ipc(i)) / rmax(ipc(i))
        vmadpl(iplane(i)) = vmadpl(iplane(i))
     .        + (ves(ipc(i)) - vintra) / natpln(iplane(i))
        vharpl(iplane(i)) = vharpl(iplane(i))
     .        + ves(ipc(i)) / natpln(iplane(i))
        vshfpl(iplane(i)) = vshfpl(iplane(i)) + vshft(i) / natpln(iplane(i))
      enddo

C ... Setup for linear potential drop across plat(3)
      if (vrl /= 0) then
        call dinv33(plat,1,qlat,vol)
        do  i = 1, 3
          b3(i) = qlat(i,3)/dlength(3,qlat(1,3),1)
        enddo
        Lz = b3(1)*plat(1,3) + b3(2)*plat(2,3) + b3(3)*plat(3,3)
      endif

      allocate(vdipib(nbasp))
      call dplmom(nbasp,pos,alat,plat,ipc,qt,0d0,vdipib,delq,delv)

C --- Dipole moment from plane (1,i), eliminating monopole term ---
C      allocate(ipcl(nbasp),vdipib(nbasp))
C!     forall (j=1:nbasp) pos(3,j) = pos(3,j)+1 ! Debugging: confirm dpl indep of shift
C      do  i = 1, nplane
C        call icopy(nbasp,ipc,1,ipcl,1)
C        forall (j=1:nbasp, iplane(j) > i) ipcl(j) = 0
C!        call pshpr(1)
C        call dplmom(nbasp,pos,alat,plat,ipcl,qt,0d0,vdipib,xx,delv(i))
C!        call poppr
CC       call dplmom(nbasp,pos,alat,plat,ipcl,qt,xx/sum(min(ipcl,1)),vdipib,delq,delv(i))
C      enddo
C      deallocate(ipcl)

      write(stdo,335)
  335 format(/'   ib      Class',7x,'Plane',6x,'Qtot      vrl*z/L   Vh(Rmax)    vshft    dipl cntr')

      allocate(ibpl(nbasp))
      call ivheap(1,nbasp,iplane,ibpl,1)
      xx = 0
      do  j = 1, nbasp
        ib = ibpl(j)
        ic = ipc(ib)
        if (ib <= nbasp-npadl-npadr) then
          cc = ' '
        elseif (ib <= nbasp-npadr) then
          cc = '(L) '
        else
          cc = '(R) '
        endif
        if (vrl /= 0) then
          xx = (pos(1,ib)*b3(1)+pos(2,ib)*b3(2)+pos(3,ib)*b3(3))/Lz*vrl
        endif
        write(stdo,334) ib, cc, ipc(ib), clabl(ic), iplane(ib), qt(ic),xx,ves(ic),vshft(ib),vdipib(ib)
  334   format(i5,a4,i4,':',a8,i5,2x,6f11.6)
      enddo

C --- Find constant-field term to make potential periodic ---
      sumq = 0d0
      qdipol = 0d0
      field0 = 0d0
      imnarg = min0(iarg(1),iarg(2))
      imxarg = max0(iarg(1),iarg(2))
      imnpln = 1

      do  20  i = 1, nplane
        sumq = sumq + qplan(i)
        if ( (i >= imnarg) .and. (i <= imxarg) )
     .          qdipol = qdipol + qplan(i)
        field0 = field0 - alat*plane(i)*qplan(i)
        if (plane(i) < plane(imnpln)) imnpln = i

        do  15  j = 1, nplane
          if (plane(j) <= plane(i))
     .      qzpl(i) = qzpl(i) + alat*( plane(i) - plane(j) )*qplan(j)
   15   continue

   20 continue

      field0 = field0 + (cellen + alat*plane(imnpln))*sumq
      field0 = field0 / cellen

C --- Printout ---
      if (iprint() < 20) return

      do  25  i = 1, 2
   25 write(lgunit(i),332) pi8*field0 / carea

C  332 format(//,' Constant field that makes V periodic ',f12.8,//,
C     .          2x,'Plane       qt         z*qt       potqz   ',
C     .          '   vmadpl      vharpl      moment',/)
  332 format(//,' Constant field to make V periodic = ',f12.8,//
     .          2x,'Plane       qt         z*qt       potqz   ',
     .          '   vmadpl      vharpl      vshfpl        z   ',/)


C --- Calculate potential from planar charges ---
      do  30  i = 1, nplane
        potqz(i) = (pi8 / carea)*( -qzpl(i)
     .    + alat*( plane(i) - plane(imnpln) )*field0 )

        write(stdo,333) i, qplan(i),qzpl(i),potqz(i),vmadpl(i),vharpl(i),vshfpl(i),
     .    plane(i)

  333   format(i5,1x,10f12.6)

   30 continue

C --- Calculate dipole (3 ways) ---
      if (iarg(1) > 0 .and. iarg(2) > 0) then
        dipole = potqz(iarg(2)) - potqz(iarg(1))
        dipmad = vmadpl(iarg(2)) - vmadpl(iarg(1))
        diphar = vharpl(iarg(2)) - vharpl(iarg(1))
      else
        dipole = 0
        dipmad = 0
        diphar = 0
      endif

      do  40  i = 1, 2
   40 write(lgunit(i),337) sumq, iarg(1), iarg(2), qdipol, dipole,
     .  rytoev*dipole, dipmad, rytoev*dipmad, diphar, rytoev*diphar

  337 format(//,' Total Charge = ',f12.7,//,' Dipoles for planes ',
     .  i3,' through ',i3,/,' Charge in these planes = ',f12.7,//,
     .  ' Dipole(planar charges) =     ',f12.8,' = ',f12.8,' eV',
     .  ' (from potqz)'/,
     .  ' Dipole(Madelung potential) = ',f12.8,' = ',f12.8,' eV',
     .  ' (from vmadpl)'/,
     .  ' Dipole(Hartree potential) =  ',f12.8,' = ',f12.8,' eV',
     .  ' (from vharpl)')

      if ( dabs(darg(1)) + dabs(darg(2)) < 1.d-5 ) return

C --- Calculate difference between the vacuum and the
C     bulk Fermi level (or VBM) for the surface (3 ways) ----
      wfbpln = rytoev*( dipole - darg(2) )
      wfbmad = rytoev*( dipmad - darg(2) )
      wfbhar = rytoev*( diphar - darg(2) )

C --- Calculate surface work function - difference
C     between vacuum and the supercell Fermi level (2 ways) ---
      wfsmad = rytoev*( vmadpl(iarg(2)) - darg(1) )
      wfshar = rytoev*( vharpl(iarg(2)) - darg(1) )

C --- Calculate Fermi level relative to
C     the bulk Fermi level (or VBM) (3 ways) ---

      efpln = darg(1) - potqz(iarg(1)) + potqz(iarg(2))
     .  - 0.5d0*vmadpl(iarg(2)) - 0.5d0*vharpl(iarg(2)) - darg(2)
      efmad = darg(1) - vmadpl(iarg(1)) - darg(2)
      efhar = darg(1) - vharpl(iarg(1)) - darg(2)

      do  50  i = 1, 2
   50 write(lgunit(i),338) darg(1), rytoev*darg(1), darg(2),
     .    rytoev*darg(2), wfbpln, wfbmad, wfbhar, wfsmad,
     .    wfshar, efpln, rytoev*efpln, efmad, rytoev*efmad,
     .    efhar, rytoev*efhar

  338 format(' Fermi level of supercell = ',f10.6,' = ',f12.8,' eV',/,
     .  ' Fermi level of bulk =      ',f10.6,' = ',f12.8,' eV',//,
     .  ' Separation between bulk Fermi level (or VBM) and vacuum:',//,
     .  '    PE threshold(planar charges) =     ',f12.8,' eV',/,
     .  '    PE threshold(Madelung potential) = ',f12.8,' eV',/,
     .  '    PE threshold(Hartree potential) =  ',f12.8,' eV',//,
     .  ' Separation between supercell Fermi level and vacuum:',//,
     .  '    Work Function(Madelung potential) = ',f12.8,' eV',/,
     .  '    Work Function(Hartree potential) =  ',f12.8,' eV',//,
     .  ' Supercell Fermi level relative to bulk Fermi level ',
     .  '(or VBM):',//,
     .  '    Fermi level(planar charges) =     ',f12.8,
     .  ' = ',f12.8,' eV',/,
     .  '    Fermi level(Madelung potential) = ',f12.8,
     .  ' = ',f12.8,' eV',/,
     .  '    Fermi level(Hartree potential) =  ',f12.8,
     .  ' = ',f12.8,' eV',/)

      end
C      subroutine lmplio(recln0,unit,offset,iosw,ido,
C     .  plat,planvc,iout,iarg,darg,carg)
CC- Control for plane subroutines
C
C      implicit none
C      double precision plat(3,3)
C
CC control parameters
C      integer recln
C      integer unit,offset,ido,iosw,recln0
C      integer mxchr,loop0,nlin,nlist,ctlen
C      parameter (mxchr=20,ctlen=120)
C      character ctbl(mxchr,2)*(ctlen), cxx*1
C
CC output
C      character*80 carg(10)
C      integer iout, iarg(10)
C      double precision darg(10)
C
CC local variables
C      integer recl0
C      parameter (recl0=72)
C      character*1 recrd(0:1)
C      character*(recl0) a
C      integer i,i1,iosw2,i0,partok,i1mach
C      logical F,sw
C      double precision cost,hvc,xx
C
CC --- Common blocks ---
C      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
C     .        iend,ichoos,nchoos,optio
C      logical noerr
C
C      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
C     .               iend,ichoos,nchoos,noerr,optio
C      common /ww/ recrd
CC     data T /.true./, F /.false./
C
C      recoff = offset
C      iout = 0
C      hvc = 1  ! should be 1 now
C
CC --- Do according to switch iosw ---
C      do  1000  iosw2 = iosw, iosw + iosw/3
C      optio = iosw2 - 2*(iosw/3)
C
C      goto (1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
C     .  170), ido+1
C      i = recln(recln0)
C      return
C
C    1 if (unit /= i1mach(1)) rewind unit
C      nrecs = 0
C      call rdfiln(unit,'#{}%',0,loop0,nlin,xx,100,
C     .  xx,nlist,cxx,ctbl,mxchr,a,recrd(recoff),reclen,nrecs)
Cc      print *,nrecs, ' lines'
C      return
C
CC --- IO ---
C   10 call getcat(recrd,'io ',' ',f)
C      if (.not. noerr) goto 1000
C      i0 = partok(recrd(catbeg),'show=','=',sw,' ',-1,0,0,0)
C      if (noerr .and. sw) iosw = 3
C      i0 = partok(recrd(catbeg),'verbos=','=',i1,' ',-1,2,0,0)
C      if (noerr) call pshprt(i1)
C      i0 = partok(recrd(catbeg),'iactiv=','=',sw,' ',-1,0,0,0)
C      if (noerr) call initqu(sw)
C      goto 1000
C
C   20 call getcat(recrd,'help ',' ',f)
C      if (noerr) then
C        iosw = 0
C        return
C      endif
C      goto 1000
C
C   30 call getcat(recrd,'write ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) print *,
C     .  '  ... writes CLASS, SITE and START to log file'
C      if (optio == 1 .and. noerr) iout = 3
C      goto 1000
C
C   40 call getcat(recrd,'rplan ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) then
C        print *,
C     .    '  ... rplan h  reads next plane from plan file and puts at h'
C        print *,
C     .    '  ... rplan (x,y,z) h  same but puts at h / h.(x,y,z)'
C      endif
C      call getxyz('rplan ',4,darg,i)
C      if (optio == 1 .and. noerr) then
C        iout = 4
CC --- evaluate h from projection ---
C        if (i > 0) then
C          cost = plat(i,3)/hvc
C          darg(1) = darg(1)/cost
C        endif
C      endif
C      goto 1000
C
C   50 call getcat(recrd,'wplan ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) print *,
C     .  '  ... wplan n  writes plane n to plan file'
C      i0 = partok(recrd(catbeg),'wplan ',' ',iarg,' ',-1,2,0,1)
C      if (optio == 1 .and. noerr) iout = 5
C      goto 1000
C
C   60 call getcat(recrd,'splan ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) print *,
C     .  '  ... splan n  skips to plane n in plan file'
C      i0 = partok(recrd(catbeg),'splan ',' ',iarg,' ',-1,2,0,1)
C      if (optio == 1 .and. noerr) iout = 6
C      goto 1000
C
C   70 call getcat(recrd,'name ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0)
C     .  print *, '  ... name ext changes file extension to ext'
C      i = partok(recrd(catbeg),'name ',' ',4,carg,1,1,0,1)
C      if (optio == 1 .and. noerr) iout = 7
C      goto 1000
C
C   80 call getcat(recrd,'stretch ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) then
C        print *,
C     .  '  ... stretch h  stretches plat(3) by h, or:'
C        print *,
C     .  '  ... stretch (x,y,z) h  stretches plat(3) by h / h.(x,y,z)'
C      endif
C      call getxyz('stretch ',4,darg,i)
C      if (optio == 1 .and. noerr) then
C        iout = 8
CC --- evaluate h from projection ---
C        if (i > 0) then
C          cost = plat(i,3) / hvc
C          darg(1) = darg(1)/cost
C        endif
C      endif
C      goto 1000
C
C   90 call getcat(recrd,'qsum ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) then
C        print *, '  ... qsum ip1 ip2 [efcell efblk]'
C        print *, '      ',
C     .  'calculates dipole betw/ planes ip2 and ip2, work func, and Ef'
C      endif
C      if (optio == 1 .and. noerr) iout = 9
C      i0 = partok(recrd(catbeg),'qsum ',' ',darg,' ',4,4,0,1)
C      iarg(1) = idnint(darg(1))
C      iarg(2) = idnint(darg(2))
C      if (dabs(iarg(1)-darg(1)) > 1d-5 .or.
C     .    dabs(iarg(2)-darg(2)) > 1d-5 .and. optio == 1)
C     .  call fexit(-1,119,'qsum expects first 2 args to be integer',0)
C      darg(1) = darg(3)
C      darg(2) = darg(4)
C      goto 1000
C
C  100 call getcat(recrd,'slide ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) then
C        print *,
C     .  '  ... slide n1 n2 (x,y,z) h  slides planes n1..n2 by h'
C      endif
C      i0 = partok(recrd(catbeg),'slide ',' ',iarg,' ',-2,2,0,1)
C      call getxyz('slide ',3,darg,i)
C      if (optio == 1 .and. noerr) then
C        iout = 10
CC --- evaluate h from projection ---
C        if (i > 0) then
C          cost = plat(i,3) / hvc
C          darg(1) = darg(1)/cost
C        endif
C      endif
C      goto 1000
C
C  110 call getcat(recrd,'cnam ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0)
C     .  print *, '  ... cnam renames classes'
C      if (optio == 1 .and. noerr) iout = 11
C      goto 1000
C
C  120 call getcat(recrd,'vshift ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) then
C        print *,
C     .  '  ... vshift n1 n2 v  shifts potential in planes n1..n2 by v'
C      endif
C      i0 = partok(recrd(catbeg),'vshift ',' ',darg,' ',-3,4,0,1)
C      if (optio /= 1) goto 1000
C      iarg(1) = nint(darg(1))
C      iarg(2) = nint(darg(2))
C      if (dabs(iarg(1)-darg(1)) > 1d-5 .or.
C     .    dabs(iarg(2)-darg(2)) > 1d-5 .and. optio == 1)
C     .  call fexit(-1,119,'vshift expects first 2 args to be integer',0)
C      darg(1) = darg(3)
C      if (optio == 1 .and. noerr) iout = 12
C      goto 1000
C
C  130 call getcat(recrd,'rot ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) print *,
C     .  '  ... rot x rotates plat by x (fractions of pi) about vc'
C      i0 = partok(recrd(catbeg),'rot ',' ',darg,' ',-1,4,0,1)
C      if (optio == 1 .and. noerr) iout = 13
C      goto 1000
C
C  140 call getcat(recrd,'proj ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0)
C     .  print *,
C     .  '  ... displays projection of basis vectors into normal plane'
C      if (optio == 1 .and. noerr) iout = 14
C      goto 1000
C
C  150 call getcat(recrd,'pos ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) then
C        print 345
C  345   format(
C     .'  ... pos [l m n lo hi]  [const expr] logs atom positions for'/
C     .'      lattice vectors (0..l,0..m,0..n) projecting betw. (lo,hi)'/
C     .'      onto normal.  Logical expr uses vars x,y,z for position'/
C     .'      (eg x<2&y>1) and px,py,pz or n for projection onto normal')
C        goto 100
C      endif
C      if (optio == 1 .and. noerr) iout = 15
C      if (optio == 1) darg(4)=-99
C      if (optio == 1) darg(5)=-99
C      i0 = partok(recrd(catbeg),'pos ',' ',darg,' ',5,4,0,1)
C      iarg(1) = idnint(darg(1))
C      iarg(2) = idnint(darg(2))
C      iarg(3) = idnint(darg(3))
C      if (i0 == 4) darg(5) = darg(4)
C      i0 = partok(recrd(catbeg),'const ',' ',60,carg,1,1,0,0)
C      goto 1000
C
C  160 call getcat(recrd,'wsite ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) then
C        print *, '  ... writes site data to site file'
C        print *, '      Syntax is : wsite [-pad | -ppad] [filename]'
C        print *, '        -pad and -ppad add padded'
C        print *, '          (or doubly padded) sites to site file.'
C        print *, '      If [filename] is not included, site filename'
C        print *, '      defaults to ''site''.'
C      endif
CC     There seems to be a bug in partok.  Do this way
C      i = partok(recrd(catbeg),'wsite ',' ',8,carg,1,1,0,0)
C      if (carg(1)(1:1) == '-') then
C        i = partok(recrd(catbeg),'wsite ',' ',8,carg(2),2,1,0,0)
C      endif
C
C      if (optio == 1 .and. noerr) iout = 16
C      goto 1000
C
C  170 call getcat(recrd,'q ',' ',f)
C      if (.not. noerr) goto 1000
C      if (iosw == 0) print *, '  ... quits program'
C      if (optio == 1 .and. noerr) call fexit(-1,119,' ',0)
C      goto 1000
C
C 1000 continue
C      end
      subroutine getxyz(ccat,nch,h,i)
C- Get distance along a line given projection along component
      use iolib_cb
      implicit none
      integer i,nch
      character*(*) ccat
      double precision h

C --- iolib functions ---
      integer partok, i0

C --- Common blocks ---
      character*1 recrd(0:1)
!       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .        iend,ichoos,nchoos,optio
!       logical noerr
!
!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio
      common /ww/ recrd

C local variables
C      logical T,F
C      data T /.true./, F /.false./


      nchoos = nch
      i = -1
      i0 = partok(recrd(catbeg),'x ',' ',h,' ',-1,4,0,1)
      if (noerr) then
        i=1
        goto 85
      endif
      i0 = partok(recrd(catbeg),'y ',' ',h,' ',-1,4,0,1)
      if (noerr) then
        i=2
        goto 85
      endif
      i0 = partok(recrd(catbeg),'z ',' ',h,' ',-1,4,0,1)
      if (noerr) then
        i=3
        goto 85
      endif
      if (nchoos >= 4)
     .  i0 = partok(recrd(catbeg),ccat,' ',h,' ',-1,4,0,1)
   85 continue

      end
      subroutine gtplan(plat,npadl,npadr,nbasp,bas,r,ipc,pl,clabl,nplane,planes,
     .  iplane,natpln)
C- Enumerate all planes in basis, make table of planes and pointers to table
C ----------------------------------------------------------------
Ci Inputs
Ci   nbasp,bas
Ci   r:  vector defining plane normal
Ci   clabl, ipc (for printout only)
Co Outputs
Co   nplane,planes:  number of planes and component along r
Co   iplane:         maps basis index to plane index
Co   natpln:         number of atoms in a plane
Cr Remarks
Cr   planes and iplane must be dimensioned sufficiently large
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer npadl,npadr,nbasp,nplane,iplane(*),natpln(*),ipc(*),pl(*)
      character*8 clabl(1)
      double precision bas(3,nbasp),r(3),planes(*),plat(3,3)
C ... Local parameters
      integer i,j,ii,iprint,lgunit,nbas,stdo
      double precision ddot,z,tiny,dsqrt,hold,v(3)
      character*4 cc
      parameter (tiny = 1d-4)

      nbas = nbasp - npadl - npadr
      stdo = lgunit(1)

C --- Find all planes, generate pointer table ---
C     and determine number of atoms per plane
      nplane = 0
      do  ii = 1, nbasp
        call pplan1(1,ii,nbas,npadl,npadr,i)
        z = ddot(3,bas(1,i),1,r,1) / dsqrt(ddot(3,r,1,r,1))
        do  j = 1, nplane
          if (dabs(z-planes(j)) < tiny) then
            iplane(i) = j
            natpln(j) = natpln(j) + 1
            goto 10
          endif
        enddo
        nplane = nplane+1
        planes(nplane) = z
        iplane(i) = nplane
        natpln(nplane) = 1
   10 continue
      enddo

C --- Printout ---
      if (iprint() < 30) return
      write(stdo,333) nplane, r
      write(stdo,336) (ddot(3,plat(1,i),1,r,1) / dsqrt(ddot(3,r,1,r,1)), i=1,3)
  333 format(/' Gtplan:',i3,' different planes normal to (', 3f12.6,')')
  336 format(' Projection of normal onto P1, P2, P3:  ',3f12.6,
     .  /,'   ib      Class',7x,'Plane  PL ',6x,
     .  'x-x.h     y-y.h     z-z.h      h       del h',/)

C     First offset is [last plane - projection of plat onto normal]
      hold = planes(iplane(nbasp)) - ddot(3,plat(1,3),1,r,1) / dsqrt(ddot(3,r,1,r,1))
      do  ii = 1, nbasp
        call pplan1(1,ii,nbas,npadl,npadr,i)
        z = ddot(3,bas(1,i),1,r,1) / dsqrt(ddot(3,r,1,r,1))
        call dcopy(3,bas(1,i),1,v,1)
        call daxpy(3,-z/dsqrt(ddot(3,r,1,r,1)),r,1,v,1)
        if (i <= nbas) then
          cc = ' '
        elseif (i <= nbas+npadl) then
          cc = '(L) '
        else
          cc = '(R) '
        endif
        write(stdo,334) i, cc, ipc(i), clabl(ipc(i)), iplane(i), pl(i),
     .  (v(j), j=1,3), planes(iplane(i)), planes(iplane(i))-hold
  334   format(i5,a4,i4,':',a8,i5,i5,3x,5f10.5)
        hold = planes(iplane(i))
      enddo
      end
      subroutine fixcl(nbas,nclspp,ipc,clabl,iwk)
      implicit none
      integer nbas,nclspp,ipc(16),iwk(16)
      character*8 clabl(8)

C Local variables
      integer ib,ic,i,j,ndel

c      nclspp = 20
c      clabl(17)='a'
c      clabl(18)='b'
c      clabl(19)='c'
c      clabl(20)='d'
c      IPC(2) = 18
c      IPC(3) = 20
c      IPC(14) = 18
c      IPC(15) = 20
c      PRINT *, IPC

C --- Eliminate replication of class names ---
      do  10  ic = 1, nclspp
        do  13  j = 1, ic-1
          if (clabl(j) /= clabl(ic)) cycle
          print 432, clabl(ic), ic, j
  432     format(' fixcl:  reset class ',a4,', number',i3,
     .           ' to earlier class',i3)
          do  20  ib = 1, nbas
            if (ipc(ib) == ic) ipc(ib) = j
   20     continue
          goto 10
   13   continue
   10 continue

C --- Eliminate unused classes ---
      call icopy(nbas,ipc,1,iwk,1)
      call ishell(nbas,iwk)
      ib = 1
   40 ib = ib+1
        if (iwk(ib) == iwk(ib-1) .or. iwk(ib) == iwk(ib-1)+1)
     .    goto 49
        print *, ' fixcl:  discarding unused class(es)',
     .    iwk(ib-1)+1, ' through ', iwk(ib)-1

C --- Collapse tables clabl and ipc --
        ndel =  iwk(ib) - iwk(ib-1) - 1

c        PRINT *, CLABL
        do  45  ic = iwk(ib), nclspp
   45   clabl(ic-ndel) = clabl(ic)
        nclspp = nclspp - ndel
c        PRINT *, CLABL

        do  48  i = 1, nbas
   48   if (ipc(i) >= iwk(ib)) ipc(i) = ipc(i)-ndel
        call icopy(nbas,ipc,1,iwk,1)
        call ishell(nbas,iwk)

c        PRINT *, IPC
c        PRINT *, IWK

   49 if (ib < nbas) goto 40
C Highest class is no larger than iwk(nb) ...
      nclspp = min(nclspp,iwk(nbas))

C For debugging, check ipc
      do  50  ib = 1, nbas
   50 if (ipc(ib) > nclspp)
     .    call fexit(-1,119, 'fixcl: bad ipc',0)

c      PRINT *, IPC
      end
      subroutine pplan1(mode,ib,nbas,npadl,npadr,ibp)
C- Return a permuted site index that orders sites L,1:nbas,R
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 for forward permutation, 1 for reverse
Ci   ib    :unpermuted site index
Ci   nbas  :size of basis without padding
Ci   npadl :number of L padding sites
Ci   npadr :number of R padding sites
Co Outputs
Co   ibp   :forward permutation reorders sites as follows:
Co         :ib  =  1:nbas,1+nbas:npadl+nbas,1+nbas+npadl:npadr+nbas+npadl
Co         :ibp =  npadl+1:npadl+nbas,1:npadl,1:nbas,1:npadr
Co         :reverse permutation not implemented
Co         :ib  =  1:nbas,1+nbas:npadl+nbas,1+nbas+npadl:npadr+nbas+npadl
Co         :ibp =  npadl+1:npadl+nbas,1:npadl,1:nbas,1:npadr
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   31 Oct 03 First created
C ----------------------------------------------------------------------
      implicit none
      integer mode,ib,nbas,npadl,npadr,ibp

      if (mode == 0) then
        if (ib <= nbas) then
          ibp = ib+npadl
        elseif (ib <= nbas+npadl) then
          ibp = ib - nbas
        else
          ibp = ib
        endif
      else
        if (ib <= npadl) then
          ibp = ib+nbas
        elseif (ib <= nbas+npadl) then
          ibp = ib - npadl
        else
          ibp = ib
        endif
      endif
      end
