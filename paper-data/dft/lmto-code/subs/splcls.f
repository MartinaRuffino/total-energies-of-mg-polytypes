      subroutine splcls(nosplt,bas,nbas,ng,istab,nspec,slabl,nclass,ipc,ics,nrc)
C- Splits species into classes
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nosplt:   T copy class and species
Ci   bas,nbas: dimensionless basis vectors, and number
Ci   nspec:    number of species
Ci   ipc:
Ci   slabl:    on input, slabl is species label
Ci   ng:       number of group operations
Ci   istab:    site ib is transformed into istab(ib,ig) by grp op ig
Ci   slabl:    species labels
Cio Inputs/Outputs:
Cio   ipc:     on input, site j belongs to species ipc(j). ipc(j=1:nbas)
Cio            on output site j belongs to class ipc(j)
Co  Outputs:
Ci   ics:      ics(i=1:nclass) class i belongs to species ics(i)
Co   nclass:   number of classes
Co   nrc:      number of classes per each species
Cu Updates
Cu   05 Apr 18 splcls now calls pvsym2 to get class labels, for consistency
Cu   05 Apr 16 (DMT) reindex ipc more reliably and quickly
Cu   11 Mar 13 species which have no sites are assigned no class
Cu   04 Apr 03 Search for equivalent classes more thorough
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed parameters:
      logical nosplt
      integer nbas,nspec,nclass,ng,istab(nbas,ng),ipc(nbas),
     .  ics(*),nrc(nspec)
      double precision bas(3,*)
      character*8 slabl(nspec)
C Local parameters:
!       include 'structures.h'
      type(str_site):: s_site(1)  ! Not used
      real(8), allocatable :: dclabl(:)
      integer ib,ic,icn,iclbsj,ig,jb,m,i,is,ipr,idx,stdo,ns
      logical lyetno,lsplit
      character*80 outs,clabl*8
      procedure(integer) :: nglob

      stdo = nglob('stdo')
      call getpr(ipr)

      nrc = 0
      do  i = 1, nbas
        nrc(ipc(i)) = 1
      enddo

      ns = 0   ! Number of species used
      do  i = 1, nspec
        if (nrc(i) /= 0) then
          ns = ns+1
          ics(ns) = i
        end if
      end do
      nclass = ns ! at this stage only, nclass and nrc could change with symmetry (or rather lack thereof) later

! --- Prepare to use nrc temporarily as function mapping species -> class.
!     Same thing can be achieved without ics, by cumulative sum over nrc,
!     then quickly reversed by step differentiation
      do i = 1, ns
        nrc(ics(i)) = i
      end do

! --- Reindex the classes skipping unused species, before this step ipc assumed
!     class == specie. This is the only step walking over all atoms besides init.
      do  i = 1, nbas
        ipc(i) = nrc(ipc(i))
      enddo

! --- Reverse nrc to its proper meaning. Multiplicity of classes for each species,
!     set to 1 at this stage, to be expanded later.
      where (nrc /= 0) nrc = 1

      lsplit = .false.  ! Turns .true. if any species is split
      if ( .not. nosplt) then
C --- For each species, split to make equivalent classes ---
      ic = 1
      do while (ic <= nclass)
        is = ics(ic)
C       if (nrs(is) == 0) cycle
        ib = iclbsj(ic,ipc,-nbas,1)
C   ... No sites of this class ... skip
        if (ib < 0) goto 11
        lyetno = .true.
C   ... For each basis atom in this class, do
        do  jb = 1, nbas
C         print *, 'jb=',ib,jb,nclass
          if (ipc(jb) == ic) then
C      ... If there is a g mapping ib->jb, sites are equivalent
            do  ig = 1, ng
              if (istab(ib,ig) == jb) goto 6
            enddo
C      ... If there is a g mapping jb->ib, sites are equivalent
            do  ig = 1, ng
              if (istab(jb,ig) == ib) goto 6
            enddo
C          04 Apr 03 consider these other possibilities
C      ... If there is a g mapping ib->kb,jb, sites are equivalent
            do  ig = 1, ng
              if (istab(istab(ib,ig),ig) == jb) goto 6
            enddo
C      ... If there is a g mapping jb->kb,ib, sites are equivalent
            do  ig = 1, ng
              if (istab(istab(jb,ig),ig) == ib) goto 6
            enddo
C      ... There wasn't an equivalent site
            if (ipr >= 70) then
              write(stdo,400) slabl(is),ib,(bas(m,ib),m = 1,3),
     .                                  jb,(bas(m,jb),m = 1,3)
            endif
C      ... If the classes haven't been split yet, do so
            if (lyetno) then
              nclass = nclass+1
              icn  =  nclass
              ics(icn) = is
              nrc(is) = nrc(is)+1
              lyetno = .false.
              lsplit = .true.
            endif
            if (nclass > nbas) then
              call rx('splcls:  problem with istab')
            endif
            icn  =  nclass
            ipc(jb)=  icn
          endif
    6   enddo
   11   ic = ic + 1
      enddo
      endif

C --- Printout ---
      if (.not. lsplit .or. ipr < 20) return
      call info5(20,1,0,' SPLCLS:  %i species '//
     .  '%?;n;(%i participating) ;%j;split into %i classes',
     .  nspec,nspec-ns,ns,nclass,0)
      if (ipr <= 30) return

      allocate(dclabl(nclass))
      call pvsym2(2,nbas,nclass,ics,ipc,nspec,slabl,s_site,dclabl,[0])

      call info0(30,0,0,' Species  Class      Sites...')
      do  is = 1, nspec
        if (nrc(is) == 1 .and. ipr < 40) cycle
        outs = ' '//slabl(is)
C       Loop should be big enough to encompass all occurences of species
C       idx = idx'th repetition of this spec
C       ic =  index to class pointing to idx'th repetition of spec
        do  idx = 1, max(nbas,nspec)
          ic = iclbsj(is,ics,-nclass,idx)
          if (ic <= 0) exit ! no such occurrence
          call r8tos8(dclabl(ic),clabl)
          outs = ' '
          if (idx == 1) outs = ' '//slabl(is)
          call awrit1('%(n>9?9:10)p%i:'//clabl,outs,80,0,ic)
          do  ib = 1, nbas
            if (ipc(ib) == ic)
     .        call awrit1('%a%(p>20?p:20)p %i',outs,80,0,ib)
          enddo
          call awrit0(outs,' ',-80,stdo)
        enddo
      enddo
  400 format(' SPLCLS: species: ',a,'has inequivalent positions:'/
     .       '  IB: ',i3,',  POS=',3f10.5/
     .       '  JB: ',i3,',  POS=',3f10.5)

      deallocate(dclabl)
      end
