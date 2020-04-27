      subroutine sfill(sfargs,slabl,s_ctrl,s_lat,s_spec,s_site)
C- Adjust sphere sizes to fill to specified fraction of cell volume
C ----------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nspec modep npadl npadr nclasp ics
Ci                 sclwsr omax1 omax2 wsrmax rmax
Co     Stored:     rmax
Co     Allocated:  *
Cio    Elts passed:rmax dclabl ipc
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw alat plat vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rmt z mxcst
Co     Stored:     rmt
Co     Allocated:  *
Cio    Elts passed:mxcst
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   sfargs string specifying how rmax to be scaled.
Ci          First character is the delimiter separating switches.
Ci          Last argument specifies which classes are to be included.
Ci          For now, there are no special switches.
Ci   slabl :vector of species labels
Cr Remarks
Cr   Using automatic scaling routine sclwsr, program tries to
Cr   rescale sphere radii so that sum-of-sphere volumes = specified
Cr   volume fraction of cell volume.
Cr   Constraints are imposed through the following (see routine sclwsr)
Cr     1 maximum overlaps may not exceed omax1 and omax2
Cr     2 specific radii maybe locked (2's digit of spec->mxcst)
Cr     3 If wsrmax is nonzero, no sphere may exceed wsrmax.
Cr
Cr   OLD method: user supplies class list and program scales those
Cr   classes to satisfy volume requirements.  This method is kept
Cr   for now, but will be discarded in future.
Cr   A subset of the nclass rmax are scaled to fulfill some criterion.
Cr   For now rmax are scaled fill the cell volume.
Cu Updates
Cu   12 Nov 17 Removed OLD style; new style has CL options
Cu   24 Jul 15 Replace dclabl with clabl
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   06 Sep 11 Started migration to f90 structures
Cu   17 May 02 Modified MT radii scaling to lower priority for E.S.
Cu   23 Apr 02 Added automatic rescaling using sclwsr.  New arg list.
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*(*) sfargs
      character*8 slabl(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_lat):: s_lat
      type(str_spec):: s_spec(*)
      type(str_site):: s_site(*)
C ... Dynamically allocated arrays
      integer,allocatable :: ips(:),locks(:),ivx(:)
      real(8),allocatable :: zc(:),rmts(:)
C ... Local parameters
      character dc*1
      integer nbas,nbasp,nclass,nspec,npadl,npadr,nclasp
      integer ls,i,iv(10),j1,ipr,ic,k,is,opt
      integer modep(3)
      double precision plat(3,3),alat,avw,vol,wsrmax,omax1(3),omax2(3),volfac,xx
      procedure(integer) :: wordsw,a2vec,lgunit,parg,bitand,isw

C ... Do nothing if no specifications
      call getpr(ipr)

      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nspec = s_ctrl%nspec
      modep = s_ctrl%modep
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      nclasp = s_ctrl%nclasp
      avw = s_lat%avw
      alat = s_lat%alat
      plat = s_lat%plat
      vol = s_lat%vol
      nbasp = nbas + npadl + npadr
      allocate(zc(nclasp))
      allocate(ips(nbasp))
      allocate(rmts(nspec))
      do  is = 1, nspec
        rmts(is) = s_spec(is)%rmt
      enddo
      do  ic = 1, nclasp
        is = s_ctrl%ics(ic)
        zc(ic) = s_spec(is)%z
      enddo
      do  i = 1, nbas
        ips(i) = s_site(i)%spec
      enddo

      ls = len(sfargs)
C     if (sfargs(1:4) /= 'auto') goto 100

C ... Use sclwsr
C  10 continue

C ... Override default sclwsr
      dc = sfargs(1:1)
      if (dc /= ' ') then
        k = wordsw(sfargs,dc,'sclwsr=','',j1) + 7
        if (k > 7) then
          i = a2vec(sfargs,len_trim(sfargs),k,4,', '//dc,3,2,1,iv,s_ctrl%sclwsr)
          if (i < 1) call rx('lmchk failed to parse '//trim(sfargs))
        endif

      endif

C ... Setup defaults
      volfac = s_ctrl%sclwsr
      if (volfac == 0) return
      omax1 = s_ctrl%omax1
      omax2 = s_ctrl%omax2
      wsrmax = s_ctrl%wsrmax
      opt = volfac/10
      opt = 10*opt
      volfac = volfac-opt
      allocate(locks(nspec))
      call spec2class(s_spec,nspec,-1,'mxcst',1,locks,xx)
      do  i = 1, nspec
        locks(i) = bitand(s_spec(i)%mxcst,2)
      enddo
      i = 3; if (modep(3) /= 2) i = 2
      if (i == 2) call rx('sfill needs update for 2D')

C ... Override default
      if (dc /= ' ') then
        k = wordsw(sfargs,dc,'wsrmax=','',j1) + 7
        if (k > 7) then
          if (a2vec(sfargs,len_trim(sfargs),k,4,', '//dc,3,2,1,iv,wsrmax) < 1)
     .      call rx('sfill failed to parse '//trim(sfargs))
          s_ctrl%wsrmax = wsrmax
        endif

        k = wordsw(sfargs,dc,'lock=','',j1) + 5
        if (k > 5) then
          allocate(ivx(nspec))
          is = a2vec(sfargs,len_trim(sfargs),k,2,', '//dc,3,2,nspec,ivx,locks)
          if (is < 0) call rx('sfill: failed to parse '//trim(sfargs))
          deallocate(ivx)
          do  i  = 1, is
            if (locks(i) == 1) locks(i) = 2
          enddo
        endif

      endif

      call info5(20,1,0,'%N SFILL:  automatic scaling of sphere radii, omax =%s,%3:1;1d%%'//
     .  '%?#n#  Constrain rmax<=%d##',100*omax1,isw(wsrmax /= 0),wsrmax,4,5)
      if (sum(locks(1:nspec)) /= 0) then
        call info2(20,0,0,'         lock:%n:1i',nspec,locks)
      endif

      k = wordsw(sfargs,dc,'omax=','',j1) + 5
      if (k > 5) then
        is = a2vec(sfargs,len_trim(sfargs),k,4,', '//dc,3,2,3,iv,omax1)
        if (is < 0) call rx('sfill: failed to parse '//trim(sfargs))
        if (is < 2) omax1(2) = omax1(1)
        if (is < 3) omax1(3) = omax1(2)
        s_ctrl%omax1 = omax1
      endif

      call sclwsr(opt,nbas,nbasp,nspec,alat,plat,s_lat%pos,ips,
     .  modep,slabl,zc,locks,volfac,wsrmax,omax1,omax2,rmts)

C     Poke back into spec->rmt
      do  is = 1, nspec
        s_spec(is)%rmt = rmts(is)
      enddo
C     Poke back into array->rmax, if it has been allocated
      if (associated(s_ctrl%rmax)) then
        do  ic = 1, nclasp
          is = s_ctrl%ics(ic)
          s_ctrl%rmax(ic) = rmts(is)
        enddo
      endif

C ... MPI broadcast s_ctrl%rmax and s_spec
C     (not needed since all processors go through this setp)
C     call bcast_strx(1,w,w,w,w,w,w,s_spec,w,w,nspec,nbas)

      deallocate(locks)
      return

C --- OLD style ---
C  100 continue
C      j1 = 0
C      call skipbl(sfargs,ls,j1)
C      j1 = j1+1
C      dc = sfargs(j1:j1)
C      j1 = j1+1
C      lstyle = 1
C
CC ... Return here to resume parsing for arguments
C  101 continue
C      call nwordg(sfargs,0,dc//' ',1,j1,j2)
C      j2 = j2+1
C
CC ... Parse special arguments
C      if (sfargs(j2:j2) /= ' ')  then
C        m = j1-1
C        i = parg('style=',2,sfargs,m,ls,dc,1,1,iv,lstyle)
C        j1 = j2+1
C        goto 101
C      endif
C
CC --- Change subset of rmax to fill cell volume ---
CC ... Make the list of classes to include in the scaling, backup rmax
C      allocate(clst(nclass))
C      allocate(rmx2(nclass))
C      call clist(lstyle,sfargs(j1:j2),s_ctrl%clabl,zc,nclass,nlist,clst)
C      call dcopy(nclass,s_ctrl%rmax,1,rmx2,1)
C
CC ... Cell volume
CC ... Sum of sphere volumes
C      call pshpr(1)
C      call ovlchk(nbas,nbasp,s_lat%pos,alat,s_ctrl%rmax,s_ctrl%rmax,
C     .  s_ctrl%clabl,s_ctrl%ipc,modep,plat,fovl,vols)
C
CC ... vol0 <- volume of all spheres not in list
C      call psfil2(nclass,nlist,clst,1d0,rmx2,0d0,rmx2,
C     .  s_ctrl%rmax)
C      call ovlchk(nbas,nbasp,s_lat%pos,alat,s_ctrl%rmax,s_ctrl%rmax,
C     .  s_ctrl%clabl,s_ctrl%ipc,modep,plat,fovl,vol0)
C
CC ... scale <- relative change in sphere sizes to make vols = vol
C      scale = ((vol-vol0)/(vols-vol0))**(1d0/3d0)
C      call psfil2(nclass,nlist,clst,1d0,rmx2,scale,rmx2,
C     .  s_ctrl%rmax)
C
CC ... Poke back into sspec
C      do  ic = 1, nclasp
C        rmax = s_ctrl%rmax(ic)
C        is = s_ctrl%ics(ic)
C        s_spec(is)%rmt = rmax
C      enddo
C      if (ipr >= 10) then
C        call awrit6('%N SFILL:  Cell vol %;8g   Sum of sphere vol '//
C     .    '%;8g (%1;5d)%N%9fRelative vol of %i classes to resize: %;4d'
C     .    //'  scale=%;5d',' ',160,lgunit(1),vol,vols,vols/vol,
C     .    nlist,(vols-vol0)/vols,scale)
CC        call awrit2(' change rmax in classes:  %n:1i',' ',80,
CC     .    lgunit(1),nlist,clst)
C        if (ipr >= 30 .and. dabs(scale-1d0) > 1d-5) then
C          call awrit0('  Class%8fOld rmax    New rmax',' ',80,
C     .      lgunit(1))
C          do  i = 1, nlist
C            ic = clst(i)
C            print 351, ic,s_ctrl%clabl(ic),rmx2(ic),s_ctrl%rmax(ic)
C  351       format(i3,2x,a,f10.6,f12.6)
C          enddo
C        endif
C      endif
C
CC ... Debugging to make sure it worked
C      call ovlchk(nbas,nbasp,s_lat%pos,alat,s_ctrl%rmax,s_ctrl%rmax,
C     .  s_ctrl%clabl,s_ctrl%ipc,modep,plat,fovl,vols)
C
C      call poppr
C      stop
      end
      subroutine psfil2(nclass,nlist,list,s1,rmax1,s2,rmax2,rmax)
C- Copy s1*rmax1 to rmax, except for elements in list copy rmax2*scale
      implicit none
      integer nclass,nlist,list(nlist)
      double precision s1,s2,rmax1(nclass),rmax2(nclass),rmax(nclass)
      integer i

      call dpcopy(rmax1,rmax,1,nclass,s1)
      do  i = 1, nlist
        rmax(list(i)) = rmax2(list(i))*s2
      enddo
      end
