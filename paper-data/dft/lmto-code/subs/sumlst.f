      subroutine sumlst(lopt,nchmx,nbas,ng,s_site,s_spec,
     .  sopts,mode,nsites,lsites,lmxch,nchan,lmdim,lchan,nll)
C- Set up list of sites and mode for Mulliken or partial DOS projection
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxb lmxa
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs:
Ci   lopt  :1s digit: choice of l channels
Ci         :  0 for partial DOS (lcutoff to lmxa; see Remarks)
Ci         :  1 for Mulliken analysis (lcutoff to lmxb; see Remarks)
Ci         :10s digit what to return
Ci         :  0 return all outputs
Ci         :  1 return nsites and lmxch only
Ci              (NB: if sopts contains site list, lsites also returned)
Ci              (NB: if sopts contains nl=, value of nll is altered)
Ci   nchmx :dimensions array lchan.
Ci         :nchmx=0 -> lchan is not generated
Ci   nbas  :size of basis
Ci   ng    :# of space group ops --- checks the m-resolved case
Ci   sopts  from command line (see Remarks) with
Ci           the --mull or --pdos stripped
Co Outputs:
Co   mode  :Whether to resolve by site, by l, or by lm (see remarks)
Co   nsites:number of sites for which partial dos is generated
Co   lsites:list of sites  for which partial dos is generated
Co   nchan :total number of DOS channels
Co   lmxch :maximum basis l for all of the atoms in the list
Co   lmdim :leading dimension of lchan
Co   lchan :array of channel numbers for each site,l,m; see Remarks
Co         :(lopt=1 only; not used for lopt=0)
Co   nll   :if nonzero, use lmax = nll-1 for all sites in making
Co         :dos channels.  NB: input value of nll passed by caller
Co         :is its default value.
Cr Remarks
Cr   sumlst:options passed in sopts
Cr   [~mode=#][~sites=list][~group=lst1,...][~nl=#][~lcut=l,...][~lmin=l,...]]
Cr   modes are 0: all sites atom-resolved (default --mull)
Cr             1: all sites l-resolved
Cr             2: all sites lm-resolved (default --pdos)
Cr                Note 3 is internally added to mode if a site list is
Cr                specified
Cr
Cr   The --mull (lopt=1) and --pdos (lopt=0) modes are similar.  They
Cr   differ in the default l-cutoff (lmxa for lopt=0, lmxb for lopt=1).
Cr
Cr   The norm of the overlap matrix can be decomposed in the following
Cr   ways:
Cr
Cr     Mulliken analysis: orbitals that have (R,Ylm) character have a
Cr     projection onto the overlap (see mullmf.f).  The sum of all
Cr     components adds to 1, though the decomposition depends on choice of
Cr     basis
Cr
Cr     Partial DOS: the norm can projected onto the augmentation sphere
Cr     at R, with a decomposition by l and m.
Cr
Cr   In either case each piece of the decomposition can be mapped onto a
Cr   particular (R,lm) channel.  In the Mulliken case, multiple-kappa basis
Cr   sets will have more than one orbital mapped onto the same channel.
Cr
Cr   sumlst returns a list of channels for which the decomposition is to be
Cr   mapped.  Channels are either of (R,lm) (mode=2 or 5), contracted over
Cr   l (mode=1 or 4), or contracted over both l and m (mode=0 or 3) For
Cr   given mode and triplet (R,l,m) sumlst returns a table fo channel
Cr   indices specifying which channel a component of the norm corresponding
Cr   to (R,l,m) gets mapped.
Cr
Cr   Thus sumlst returns lchan(1,ib), lchan(1..l,ib), or lchan(1..lm,ib),
Cr   depending on mode.
Cr   This is used by routines that take information from eigenvectors
Cr   and project partial dos into channels, e.g. mullmf.f and mkpdos.f
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   05 Mar 10 (MvS) bugs fix for site-dependent lmax.
Cu             New qualifier lcut=. lmdim added to argument list
Cu   30 Aug 04 Changed first argument to integer lopt
Cu   18 Jan 02 (MvS) redesigned switches for compatibility with site
Cu                   lists in other context
Cu   27 Mar 01 (MvS) extended to handle both --mull and --pdos
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*(*) sopts
      integer lopt,nchmx,nbas,ng,nsites,mode,lsites(*),nchan,lmxch,
     .  lmdim,lchan(*),nll
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer mxmode
      parameter (mxmode=6)
      character dc*1,modstr(0:mxmode-1)*128
      integer lmax,isite,lbot,ltop,ilm,nlmin
      integer i,j,ipr,iprmin,ib,is,iprint,lgunit,stdo,j2,j1,
     .  m,iv(10),parg,lstyle,lopt0,lopt1,nlcut,ngrp,jb,iinear,ksite
      integer nlist,igrp,lgrp(nbas)
      integer,allocatable:: lcut(:),lmin(:),grptbl(:,:),iwk(:)
      double precision z
      data modstr /'all sites atom-projected',
     .             'all sites l-projected','all sites lm-projected',
     .             'site list atom-projected',
     .             'site list l-projected','site list lm-projected' /

      lopt0 = mod(lopt,10)
      lopt1 = mod(lopt/10,10)
      ipr = iprint()
      iprmin = 30
      mode = -1
      dc = sopts(1:1)
      stdo = lgunit(1)
      mode = 2
      nsites = 0
      if (lopt0 == 1) mode = 0
      if (dc /= ' ') then
C   ... Return here to resume parsing for arguments
        j2 = 0
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) == dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then
C     ... DOS mode; see Remarks
          elseif (sopts(j1:j1+4) == 'mode=')  then
            m = 0
            i = parg('mode=',2,sopts(j1:),m,len(sopts(j1:)),
     .        dc//' ',1,1,iv,mode)
            if (i <= 0) goto 999
            call sanrg(.true.,mode,0,5,' sumlst:','mode')
C            if (ipr >= 10) then
C              if ((mode == 2 .or. mode == 5) .and. ng > 1)
C     .          call awrit0(' **WARNING** sumlst: '//
C     .          'for lm-decomposition suppress symops',' ',128,stdo)
C            endif
C     ... Specify upper limit to l cutoff + 1
          elseif (sopts(j1:j1+2) == 'nl=')  then
            m = 0
            i = parg('nl=',2,sopts(j1:),m,len(sopts(j1:)),
     .        dc//' ',1,1,iv,nll)
            if (i <= 0) goto 999
C     ... Specify individual l-cutoffs
          elseif (sopts(j1:j1+4) == 'lcut=')  then
            call mkils0(sopts(j1+5:j2),nlcut,j)
            if (nlcut < 0) goto 999
            if (allocated(lcut)) deallocate(lcut)
            allocate(lcut(max(nbas,nlcut)))
            call ivset(lcut,1,max(nbas,nlcut),-2)
            call mkilst(sopts(j1+5:j2),nlcut,lcut)
C     ... Specify individual l-minima
          elseif (sopts(j1:j1+4) == 'lmin=')  then
            call mkils0(sopts(j1+5:j2),nlmin,j)
            if (nlmin < 0) goto 999
            allocate(lmin(max(nbas,nlmin)))
            call ivset(lmin,1,max(nbas,nlmin),-2)
            call mkilst(sopts(j1+5:j2),nlmin,lmin)
C     ... Site list
          elseif (sopts(j1:j1+5) == 'sites=')  then
            if (mode < 3) mode = mode+3
C           We need to make slabl,z before using other lstyle
            lstyle = 1
            call slist(lstyle,sopts(j1+6:j2),' ',z,nbas,nsites,lsites)
C     ... Group list
          elseif (sopts(j1:j1+5) == 'group=')  then
            is = j1
            if (mode < 3) mode = mode+3
            j1 = j1+6
            call wrdsg(sopts(j1:j2),1,'; '//dc,ngrp)
            allocate(grptbl(nbas,ngrp))
            call ivset(grptbl,1,nbas*ngrp,0)
            i = 0
            do  while (j1 <= j2)
              i = i+1
              call nwordg(sopts,1,'; '//dc,1,j1,j)
              call mkilss(1,sopts(j1:j),nlist,grptbl(1,i))
              if (nlist > nbas) call rxs('sumlst: bad list ',sopts(j1:))
              j1 = j+2
            enddo
            if (ngrp /= i)
     .        call rxs('sumlst: failed to parse ',sopts(is:))
C           Get nsites,lsites
            nsites = 0
            do  i = 1, ngrp
            do  isite = 1, nbas
              if (grptbl(isite,i) == 0) cycle
              nsites = nsites+1
              ib = grptbl(isite,i)
              if (ib > nbas) call rxs('sumlst:  ib>nbas in list ',sopts(is:))
              if (nsites > nbas) call rxs('sumlst:  nsites>nbas in list ',sopts(is:))
              lsites(nsites) = ib
              lgrp(nsites) = i
            enddo
            enddo

C           Sanity check: no duplication of basis allowed
            if (allocated(iwk)) deallocate(iwk)
            allocate(iwk(nsites))
            call ivheap(1,nsites,lsites,iwk,1)
            do  i = 2, nsites
              if (lsites(iwk(i)) /= lsites(iwk(i-1))) cycle
              call rx2('sumlst: duplicate %i site in group %i',lsites(iwk(i)),lgrp(iwk(i)))
            enddo
          else
            goto 999
          endif
          goto 10
        endif
      endif
      if (ipr >= iprmin+10) write(stdo,1)
    1 format (' sumlst:  Site  lmin   lmax   group')
      if (mode < 3) nsites = nbas

      if (nsites <= 0) call rxi
     .  ('sumlst: --sites=list is required input for mode =',mode)

C --- Generate lsites,lmxch ---
      lmxch = -1
      do  isite = 1, nsites
        if (mode < 3) then
          ib = isite
          if (lopt1 == 0) then
            lsites(isite) = ib
          endif
        else
          ib = lsites(isite)
        endif
        call sanrg(.true.,ib,1,nbas,' sumlst:','site')
        is = s_site(ib)%spec
        if (lopt0 == 1) lmax = s_spec(is)%lmxb
        if (lopt0 == 0) lmax = s_spec(is)%lmxa
        if (nll > 0) lmax = min(nll-1,lmax)
        lmxch = max(lmxch,lmax)
      enddo
      if (lopt1 == 1) return

C --- Generate nchan,lchan ---
C     if (mode < 3) mode = mode+3
C     ksite = loop counter running nsite sites.
C     isite = current channel index: assemble lchan(1..lmdim,isite)
C             Note isite = ksite unless arranged by group
C     ib    = site index for this isite
C     i  = index to next new channel, poked into lchan
      i = 0
      do  isite = 1, nsites
        ib = lsites(isite)
C       Group ordering: isite -> group index; loop over sites in group
        if (allocated(grptbl)) then
          igrp = lgrp(isite)
          jb = grptbl(1,igrp)
        else
          jb = ib
        endif

        is = s_site(ib)%spec
        if (lopt0 == 1) lmax = s_spec(is)%lmxb
        if (lopt0 == 0) lmax = s_spec(is)%lmxa
        if (nll > 0) lmax = min(nll-1,lmax)
        if (allocated(lcut)) then
          j = lcut(isite)
          if (allocated(grptbl)) j = lcut(igrp)
          if (j /= -2) lmax = min(lmax,j)
        endif
        lbot = 0
        if (allocated(lmin)) then
          j = lmin(isite)
          if (allocated(grptbl)) j = lmin(igrp)
          if (j /= -2) lbot = max(0,j)
        endif
        if (ipr >= iprmin+10) then
          if (allocated(grptbl)) then
            write(stdo,2) ib, lbot, lmax, lgrp(isite)
          else
            write(stdo,2) ib, lbot, lmax
          endif
        endif
    2   format (9x,2i5,i6:4x,i3)
        if (mode == 0 .or. mode == 3) then
          lbot = 1
          ltop = 1
          lmdim = 1
        elseif (mode == 1 .or. mode == 4) then
          lbot = lbot+1
          ltop = lmax+1
          lmdim = lmxch+1
        elseif (mode == 2 .or. mode == 5) then
          lbot = (lbot+1)**2
          ltop = (lmax+1)**2
          lmdim = (lmxch+1)**2
        endif
C       Create channels for this site
        if (ib == jb) then
        if (nchmx > 0) then
          do  ilm = 1, lmdim
            if (ilm >= lbot .and. ilm <= ltop) then
              i = i+1
              if (i > nchmx) call rxi('sumlst: bndfp needs nchmx>',i)
              call mchan(lmdim,0d0,0d0,0,nsites,0,isite,ilm,-1,i,lchan)
            else
              call mchan(lmdim,0d0,0d0,0,nsites,0,isite,ilm,-1,0,lchan)
            endif
          enddo
        else
          i = i + ltop-lbot+1
        endif
C       Use existing channels for this site
        elseif (nchmx > 0) then
C         Site index for first member of group
          ksite = iinear(nsites,jb,lsites,1)
          call icopy(lmdim,lchan(1+lmdim*(ksite-1)),1,
     .                     lchan(1+lmdim*(isite-1)),1)
        endif
      enddo
      nchan = i

C --- Printout --
      call info5(iprmin,1,0,
     .  ' sumlst:  %?#n#Mulliken#Partial# DOS mode %i ('//
     .  trim(modstr(mode))//')  %i sites %i channels',
     .  lopt0,mode,nsites,nchan,0)

      if (allocated(grptbl)) deallocate(grptbl)
      return

  999 continue
      call rxs('sumlst: failed to parse options in ',sopts(j1:))

      end
C      subroutine snot(lmdim,nsites,lchan)
C      implicit none
CC Passed Parameters
C      integer lmdim,nsites,lchan(lmdim,nsites)
CC Local Variables
C      integer i,isite
C
C      do  isite = 1, nsites
C        print 333, isite, (lchan(i,isite), i=1,lmdim)
C  333   format(i4,2x,12i4)
C      enddo
C      end
