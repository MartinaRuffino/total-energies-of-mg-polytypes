      subroutine iopos(lio,mode,filnam,nbas,bas,s_site)
C- File I/O of site positions
C ----------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  clabel
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci  lio:    true for write, false for read
Ci  mode:   0 read all sites as a matrix.
Ci
Ci          1 read each line in the form
Ci            ib x y z
Ci          Here ib must be between 1 and nbas.
Ci          Not all sites need be specified; sites not
Ci          specified are not changed.
Ci
Ci         -1 mode unknown; must be specified by input file;
Ci          see Remarks.
Ci filnam:  file name containing data
Ci  nbas:   size of basis
Cio bas:    is read in (lio=F), or written out (lio=T)
Cr Remarks
Cr   File mode can be specified by input file, which it does
Cr   by first line beginning with % and containing mode=#, e.g.
Cr     % mode=1
Cr
Cr   See also iosits.
Cm MPI
Cm   Master does I/O. If read then bas is broadcast.
Cu Updates
Cu   04 Nov 17 New ioposlu, which accepts file pointer instead of name,
Cu             and can also read the size of the file.
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   20 Feb 12 mode 1, write species name after bas
Cu   07 Aug 10 Repackaged to eliminate explicit MPI calls
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lio
      integer mode,nbas,ifil
      character*(*) filnam
      double precision bas(3,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site),intent(inout)::  s_site(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk(:)
C ... Local parameters
      logical usefilename
      integer j1,j2,ifi,fopna,j,rdm,ipr,i,lmode,ix(4)
      double precision xv(4)
      character fnam*100, clabl*8
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),
     .  vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)
C ... MPI
      logical mlog
      integer procid,master,mpipid
      procedure(logical) :: parstr,a2bin,cmdopt
      procedure(integer) :: lgunit,a2vec,isw

      data rdarg /'#{}% c'/

      mlog = cmdopt('--mlog',6,0,fnam)
      ifi = 0; usefilename = .true.

C     Entry point for ioposlu
   10 continue
      procid = mpipid(1)
      master = 0

C     Only master node does file I/O
      if (procid == master) then

C ... Open file
      j1 = 1; j2 = 1
      if (usefilename) then
        fnam = filnam
        call nword(fnam,1,j1,j2)
        i = 1
        if (lio) i = 2
        ifi = fopna(fnam(j1:j2),-1,i)
      endif
      rewind ifi
      call getpr(ipr)

C --- File write ---
      if (lio) then
C   ... Write, mode=1
        if (mode == 1) then
          write(ifi,'(''% mode=1'')')
          do  i = 1, nbas
            clabl = s_site(i)%clabel
            write(ifi,"(i4,3f13.7,3x,'# ',a8)")
     .         i,bas(1,i),bas(2,i),bas(3,i),clabl
          enddo
          call fclose(ifi)

C   ... Write, mode=0
        elseif (mode == 0) then
          allocate(wk(3*nbas))
          call dmcpy(bas,3,1,wk,1,nbas,3,nbas)
          call ywrm(0,' ',1,ifi,'(3f14.8)',wk,1,nbas,nbas,3)
          call fclose(ifi)
          deallocate(wk)
        else
          call rxi('IOPOS: not implemented, mode=',mode)
        endif

        call info5(10,1,1,' iopos : write %i sites into file '//
     .    '%?;n;'//fnam(j1:j2)//'%j;%i;, mode %i',nbas,isw(usefilename),ifi,mode,5)

C --- File read ---
      else
C   ... Look for mode on first line
        nr = 0
        call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
        backspace ifi
        lmode = 0
        if (recrd(1:1) == '%') then
          i = 0
          j = 0
          if (parstr(recrd,'mode=',len(recrd),5,'=',i,j)) then
            i = j
            if (.not. a2bin(recrd,lmode,2,0,' ',i,-1))
     .        call rxs('IOPOS failed to parse ...',recrd(j-5:))
          endif
        endif
        if (mode >= 0 .and. lmode /= mode)
     .    call rxi('IOPOS: mode mismatch: need mode',mode)
C   ... Read, mode=0
        if (lmode == 0) then
          allocate(wk(3*nbas*2))
          j = rdm(ifi,0,3*nbas,' ',wk,nbas,3)
          if (j < 0) call rxi('IOPOS: file incompatible with nbas=',nbas)
          call dmcpy(wk,nbas,1,bas,1,3,nbas,3)
          deallocate(wk)
          nr = nbas
C   ... Read, mode=1
        elseif (lmode == 1) then
          nr = 0
   41     call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .      ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
          if (nr >= 0) then
            if (recrd(1:1) == '%') goto 41
            i = 0
            j = a2vec(recrd,len(recrd),i,4,', ',2,3,4,ix,xv)
            if (j /= 4) call rxs(
     .        'IOPOS: failed to parse line:',recrd)
            i = nint(xv(1))
            call dcopy(3,xv(2),1,bas(1,i),1)
            goto 41
          endif
          nr = -nr-1
        else
          call rxi('IOPOS: unkown mode',lmode)
        endif
C     ... Show input
        if (ipr >= 0) then
          call awrit5(' iopos : read %i sites from file '//
     .      '%?;n;'//fnam(j1:j2)//'%j;%i;'
     .      //'%a, mode=%i%?#n>=40#:'//
     .      '%N  ib       x           y           z',
     .      ' ',100,lgunit(1),nr,isw(usefilename),ifi,lmode,ipr)
          if (ipr >= 40) then
            do  i = 1, nbas
              print 310, i, bas(1,i), bas(2,i), bas(3,i)
            enddo
  310       format(i4,3f12.6)
          endif
        endif
      endif
      call fclose(ifi)

      endif

      if (.not. lio) then
        call mpibc1(bas,3*nbas,4,mlog,'iopos','pos')
      endif
      return

      entry ioposlu(lio,mode,ifil,nbas,bas,s_site)
C- Same as iopos, but specify file through already-open logical unit
C  Additional feature: on file read, if input nbas is zero,
C  ioposlu reads nbas without reading bas, and returns nbas as output

      procid = mpipid(1)
      master = 0

      ifi = ifil; usefilename = .false.
      if (.not. lio .and. nbas == 0) then

C       Only master node does file I/O
        if (procid == master) then
          rewind ifi
          j = rdm(ifi,0,0,' ',xv,nbas,3)
        endif
        call mpibc1(nbas,1,2,0,'','')
        return
      endif

      goto 10  ! Join iopos

      end
      subroutine ioposs(lio,mode,filnam,nbas,s_site)
C- File i/o of site positions from site structure
Cu   05 Jul 13 Replace f77 pointers with f90 ones
      use structures
      implicit none
      logical lio
      integer mode,nbas
      character*(*) filnam
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: pos(:)

      call rx('update ioposs')

C ... Unpack site positions
      allocate(pos(3*nbas))
C ... File I/O
      call iopos(lio,mode,filnam,nbas,pos,s_site)
C ... Repack site positions
      deallocate(pos)
      end
