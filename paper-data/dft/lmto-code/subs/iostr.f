      logical function rdstrx(lio,s_str,name,nl,nbas,ckbas)
C- Return pointers to real-space strux; read from disk first call
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio   :specifies conditions for file i/o (condensed into digits)
Ci         :digit  specification
Ci         :1-10s   bits in 1-10s digit control program flow
Ci                  0 for file read
Ci                  1 for file write
Ci                  2 (read only) No header matching requirement;
Ci                    header parms returned into passed arguments
Ci                  4 read/write file header only and exit.
Ci                  8 Allocate arrays before reading
Ci                 16 file is formatted (not implemented)
Ci         :100-1000 bits in these digits contain info about structure of s
Ci                  1 s is complex (otherwise s is real)
Ci                  2 s is stored in slices according to site.
Ci                    row ordering that of the cluster determined
Ci                    by iax table.
Ci                    Default storage (2's bit zero):
Ci                    s is stored by (nds,nds) blocks
Ci                  4 s is the structure matrix; otherwise,
Ci                    s is the matrix B of 1-c expansion coffs
Ci                  8 file contains energy derivative of s
Ci         :10000-100000 info about the functions s corresponds to
Ci                  0 2nd generation Hankel functions with Andersen's
Ci                    conventions
Ci                  1 NMTO Hankel functions with Tank and Andersen
Ci                    conventions: alpha=alpha(2,2,nl**2,nbas,nkap)
Ci                    and s=s(nl**2,nl**2,nsite,nkap)
Ci                  2 Conventional Hankel functions, Methfessel defs
Ci   name  :file name holding strux
Ci   nl    :1 + global-maximum-l for (nl**2,nl**2) format
Ci          total dimension of s for vector format
Ci   nbas  :number of atoms in the basis
Ci   nkap  :number of energies for which strux are calculated
Ci   kap2  :Hankel energies
Ci   itral :characterizes structure matrix transformations,
Ci          (Used by NMTO only; see mktra2.f)
Ci   ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Ci   lmaxw :maximum l for Watson-sphere
Co Outputs:
Co   nitab :total number of inequivalent pairs in iax table (strscr.f)
Co          Unused for vector format
Co   oalpha:offset for tight-binding screening parameters
Co   oiax  :offset for iax array
Co   onpr  :offset for number of sites in each cluster
Co   os    :offset for screened structure constant matrix
Cr Remarks
Cr   A high-level reader for strux that can store strux in RAM
Cr   to avoid excessive disk reads.
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   13 Oct 10 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character name*(*)
      integer lio,nl,nbas
      double precision ckbas
C ... For structures
!      include 'structures.h'
      type(str_str):: s_str
      type(str_str0):: s_sta
C ... Local parameters
      logical iostr,iostr1,lok,ldum,isanrg
      integer k,lio01,lio45,getdig,ii,ifi,i,j
      integer lior,lmaxwr,itralr,nitab,nlr,nbr
      integer niax,nds,nttab,ncplx,nkapr
      integer mpipid,master,procid
      parameter (niax=10)
      integer nkapsr,nkapnr !haves(2),
      double precision ckbasr
c     integer oalpha,onpr,os
!       save haves
!       data haves /0,0/
      integer, save :: haves(2) = [0,0]

      if (name == 'str' .or. name == 'STR') then
        k = 1
      elseif (name == 'sdot' .or. name == 'SDOT') then
        k = 2
      else
        call rx('rdstrx: file name not recognized')
      endif

C ... File read
C     numprocs = mpipid(0)
      procid = mpipid(1)
      rdstrx = .true.

      if (haves(k) == 0) then
        lio01 = mod(lio/1,100)
        lio45 = mod(lio/10000,100)
        master = 0
        if (mod(lio01,4) /= 0) then
          call rxi('rdstrx: illegal option lio = ',lio)
        endif
        if (procid == master) then
          ii = 0

C         Check to see whether header can be read
          lior = 0
          lok = iostr1(lior,name,itralr,nttab,nitab,nlr,nbr,
     .      nkapr,lmaxwr,ckbasr,ifi)
          if (lok) then
          if (getdig(lio01,2,2) == 1) then
C           itral = itralr
            nl = nlr
            nbas = nbr
C           nkap = nkapr
C           lmaxw = lmaxwr
            ckbas = ckbasr
          else
            ldum = isanrg(mod(lior/10000,100),lio45,lio45,
     .        'file: RDSTRX:','strux type',.true.)
          endif

C         Use nds=nl**2 for now, until nds is passed
          nds = nlr*nlr
          call iostr5(lior,nkapr,nkapsr,nkapnr,ncplx)

C         lmaxwr, itral, nkap, kap2 are not given.
C         Use file values from iostr1 or read them now
          lok = iostr(lio,name,nl,nbas,nkapr,s_str%kaps,itralr,ckbas,
     .                lmaxwr,nitab,s_sta)
          nttab = s_sta%n(nbas+1)
          endif
        endif  ! procid = master

        call mpibc1(lok,1,1,.false.,'rdstrx','error')
        if (.not. lok)
     .    call rxs('RDSTRX: failed to read strux from file ',name)

C  ...  Load dimensioning parameters into s_str
        call mpibc1(lior,1,2,.false.,'rdstrx','liorr')
        call mpibc1(nkapr,1,2,.false.,'rdstrx','nkap')
        call mpibc1(nttab,1,2,.false.,'rdstrx','nttab')
        call mpibc1(nitab,1,2,.false.,'rdstrx','nitab')
        call mpibc1(nds,1,2,.false.,'rdstrx','nds')
        call mpibc1(lmaxwr,1,2,.false.,'rdstrx','lmaxw')
        call iostr5(lior,nkapr,nkapsr,nkapnr,ncplx)
        s_str%nkaps = nkapr
        s_str%nitab = nitab
        s_str%nttab = nttab
        s_str%lmaxw = lmaxwr
        s_str%nds   = nds

        if (getdig(lio01,3,2) == 0 .or. .not. lok) then
          deallocate(s_sta%a,s_sta%i,s_sta%n,s_sta%s)
          rdstrx = lok
          return
        endif

C  ...  Load array elements into s_str
        i = nds*nbas*nkapsr**2*nkapnr
        j = nds**2*nitab*nkapsr**2*ncplx*nkapnr

C       Other MPI processes must allocate for alpha,iax,ntab,s
        if (procid /= master) then
          allocate(s_sta%a(i))
          allocate(s_sta%i(niax*nttab))
          allocate(s_sta%n(nbas+1))
          allocate(s_sta%s(j))
        endif

        call mpibc1(s_sta%a,i,4,.false.,'rdstrx','alpha')
        call mpibc1(s_sta%i,niax*nttab,2,.false.,'rdstrx','iax')
        call mpibc1(s_sta%n,nbas+1,2,.false.,'rdstrx','npr')
        call mpibc1(s_sta%s,j,4,.false.,'rdstrx','s')

! C       Copy arrays alpha,iax,ntab,s into s_str
!         call ptr_str(s_str,4+1,'iax',niax*nttab,0,s_sta%i)
!         call ptr_str(s_str,4+1,'npr',nbas+1,0,s_sta%n)
!         if (k == 1) then
!           call ptr_str(s_str,4+1,'alp',i,0,s_sta%a)
!           call ptr_str(s_str,4+1,'s',j,0,s_sta%s)
!         else
!           call ptr_str(s_str,4+1,'adot',i,0,s_sta%a)
!           call ptr_str(s_str,4+1,'sdot',j,0,s_sta%s)
!         endif

        if (associated(s_str%iax)) deallocate(s_str%iax)
        allocate(s_str%iax(niax*nttab))
        s_str%iax(1:niax*nttab) = s_sta%i(1:niax*nttab)

        if (associated(s_str%npr)) deallocate(s_str%npr)
        allocate(s_str%npr(nbas+1))
        s_str%npr(1:nbas+1) = s_sta%n(1:nbas+1)

        if (k == 1) then
         if (associated(s_str%alph)) deallocate(s_str%alph)
         allocate(s_str%alph(i))
         s_str%alph(1:i) = s_sta%a(1:i)

         if (associated(s_str%s)) deallocate(s_str%s  )
         allocate(s_str%s(j))
!          s_str%s = s_sta%s(1:j)
         call dcopy(j, s_sta%s, 1, s_str%s, 1)
        else
         if (associated(s_str%adot)) deallocate(s_str%adot)
         allocate(s_str%adot(i))
         s_str%adot(1:i) = s_sta%a(1:i)

         if (associated(s_str%sdot)) deallocate(s_str%sdot  )
         allocate(s_str%sdot(j))
!          s_str%sdot(1:j) = s_sta%s(1:j)
         call dcopy(j, s_sta%s, 1, s_str%sdot, 1)
        end if

        deallocate(s_sta%a,s_sta%i,s_sta%n,s_sta%s)
        haves(k) = 1
      endif

      end

      logical function iostr(lio,name,nl,nbas,nkap,kap2,itral,ckbas,
     .  lmaxw,nitab,s_sta)
C- Real space structure constants file read or write
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio   :specifies conditions for file i/o (condensed into digits)
Ci         :digit  specification
Ci         :1-10s   bits in 1-10s digit control program flow
Ci                  0 for file read
Ci                  1 for file write
Ci                  2 (read only) No header matching requirement;
Ci                    header parms returned into passed arguments
Ci                  4 read/write file header only and exit.
Ci                  8 Add this number to allocate arrays before reading
Ci                    or to release arrays after writing
Ci                    Bits may be take in combination
Ci                 16 If read error, print warning but do not abort
Ci         :100-1000 bits in these digits contain info about structure of s
Ci                  1 s is complex (otherwise s is real)
Ci                  2 s is stored in slices according to site.
Ci                    row ordering that of the cluster determined
Ci                    by iax table.
Ci                    Default storage (2's bit zero):
Ci                    s is stored by (nds,nds) blocks
Ci                  4 s is the structure matrix
Ci                    Default: s is the matrix of expansion coffs
Ci                  8 file contains energy derivative of s
Ci         :10000-100000 info about the functions s corresponds to
Ci                  0 2nd generation Hankel functions with Andersen's
Ci                    conventions
Ci                  1 NMTO Hankel functions with Tank and Andersen
Ci                    conventions: alpha=alpha(2,2,nl**2,nbas,nkap)
Ci                    and s=s(nl**2,nl**2,nsite,nkap)
Ci                  2 Conventional Hankel functions, Methfessel defs
Ci   name  :file name holding strux
Ci   nl    :1 + global-maximum-l for (nl**2,nl**2) format
Ci          total dimension of s for vector format
Ci   nbas  :number of atoms in the basis
Ci   nkap  :number of energies for which strux are calculated
Ci   kap2  :Hankel energies
Ci   itral :characterizes structure matrix transformations,
Ci          (Used by NMTO only; see mktra2.f)
Ci   ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Ci   lmaxw :maximum l for Watson-sphere
Cio Inputs/Outputs:
Cio  nitab :total number of inequivalent pairs in iax table (strscr.f)
Cio         Unused for vector format
Cio  oalpha:offset for tight-binding screening parameters
Cio   iax  :offset for iax array
Cio  onpr  :offset for number of sites in each cluster
Cio  os    :offset for screened structure constant matrix
Cb Bugs
Cb   Only lio 10s digit only 0 for MPI
Cr Remarks
Cr   If file read, allocates memory for and reads alpha,iax,npr,s
Cr   Input parameters nl,nbas,bas,nkap,lmaxw,itral must match file
Cr   for iostr to complete file read.
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   06 Aug 06 Redesigned to work with 2-kappa strux
Cu   8 Jun 04  (MPI) read from master node only and broadcast
Cu             (ccomp with MASTERIO.) (Warning: not checked.)
Cu   1 Aug 98  revised for 3rd generation lmto
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character name*(*)
      integer itral,lio,nitab,nl,nbas,nkap,lmaxw
      double precision kap2(nkap),ckbas
C ... For structures
!      include 'structures.h'
      type(str_str0), intent(inout) :: s_sta
C ... Local parameters
      character locnam*8
      logical iostr1,iostr2,iostr4,lerr,ldum,isanrg
      integer nlr,nbr,ifi,ipr,nkr,nttab,lgunit,itralr,lmaxwr,niax,ncplx,
     .  lio01,lio23,lio45,getdig,nds,lscoff,nmto,lior,nblk,nkaps,nkapn
      integer mpipid,procid,vsn
      double precision tiny,ckbasr
C     double precision sttime,entime,MPI_WTIME
      parameter (tiny=1d-8,niax=10,vsn=1)

C --- Setup information ---
      call getpr(ipr)
      locnam = name
      lio01 = mod(lio/1,100)
      lio23 = mod(lio/100,100)
      lio45 = mod(lio/10000,100)
      ncplx = 1+getdig(lio23,0,2)

      procid = mpipid(1)
C     Leading dimension of 'alpha' array depends on whether it
C     corresponds to 2nd gen alpha or tral matrix

C --- File write ---
      if (mod(lio,2) /= 0) then

C       Use nds=nl**2 for now, until nds is passed
        nds = nl*nl
        call iostr5(lio,nkap,nkaps,nkapn,ncplx)
        lscoff = nds*nbas*nkaps**2*nkapn

C   ... Write header
        nttab = s_sta%n(nbas+1)
        if (.not. iostr1(lio,name,itral,nttab,nitab,nl,nbas,
     .    nkap,lmaxw,ckbas,ifi)) goto 99
C   ... Write strux, block format
        if (getdig(lio23,1,2) == 0) then
        nblk = nitab*nkaps**2*nkapn
          call rxx(getdig(lio23,0,2) == 1,'iostr not ready for complex')
          if (.not. iostr2(-ifi,nbas,nttab,nblk,lscoff,nds,nkap,
     .        s_sta%n,s_sta%i,kap2,s_sta%a,s_sta%s)) goto 99
C   ... Write strux, vector format
        else
          if (.not. iostr4(-ifi,ncplx,nl,nkap,kap2,s_sta%s))
     .      goto 99
        endif

        if (ipr >= 40) call awrit4(' IOSTR: wrote to file '
     .    //locnam//'%a %i sites'//
     .    '%?#n#, %i inequivalent sites#%j#'//
     .    '%?#n>1#, %-1j%i energies##',
     .    ' ',80,lgunit(1),nttab,nttab-nitab,nitab,nkap)
        call fclose(ifi)

C --- File read ---
      else
C   ... Open file, and read header
C#ifdefC MASTERIO
C        if (procid == 0) then
C#endif
        lior = 0
        lerr = .not. iostr1(lior,name,itralr,nttab,nitab,nlr,nbr,
     .    nkr,lmaxwr,ckbasr,ifi)
        if (lerr) goto 99
C       Copy file parameters to those in argument list
        if (getdig(lio01,1,2) == 1) then
          nl = nlr
          nbas = nbr
          nkap = nkr
          itral = itralr
          ckbas = ckbasr
          lmaxw = lmaxwr
        endif

        ldum = isanrg(nlr,nl,nl,'file: IOSTR:','nl',.true.)
        ldum = isanrg(nbr,nbas,nbas,'file: IOSTR:','nbas',.true.)
        ldum = isanrg(nkr,nkap,nkap,'file: IOSTR:','nkap',.true.)
        call fsanrg(ckbasr,ckbas,ckbas,tiny,'file: IOSTR:','cksum',
     .    .true.)
C   ... This is only a warning, not fatal error
        ldum = isanrg(lmaxwr,lmaxw,lmaxw,'file: IOSTR:','lmaxw',.false.)

C#ifdefC MASTERIO
C        endif
C
C        call mpibc1(lerr,1,0,.false.,'iostr','lerr')
C        call mpibc1(nl,1,2,.false.,'iostr','nl')
C        call mpibc1(nbas,1,2,.false.,'iostr','nbas')
C        call mpibc1(nkap,1,2,.false.,'iostr','nkap')
C        call mpibc1(itral,1,2,.false.,'iostr','itral')
C        call mpibc1(lmaxw,1,2,.false.,'iostr','lmaxw')
C        call mpibc1(nttab,1,2,.false.,'iostr','nttab')
C        call mpibc1(nitab,1,2,.false.,'iostr','nitab')
C#endif

C   ... If failed to read header, give up
        if (lerr) goto 99

C       Use nds=nl**2 for now, until nds is passed
        nds = nl*nl
        call iostr5(lio,nkap,nkaps,nkapn,ncplx)
        nmto = getdig(lio45,0,2)

C   ... Allocate memory for arrays
        if (getdig(lio01,3,2) == 1) then
C          if (getdig(lio01,1,2) == 1)
C     .      call rx('iostr: illegal combination of bits in lio')
          lscoff = nds*nbas*nkaps**2*nkapn

!           call d1alloc(s_sta%a,lscoff,0)
!           call i1alloc(s_sta%i,niax*nttab,0)
!           call i1alloc(s_sta%n,nbas+1,0)
!           call d1alloc(s_sta%s,nds**2*nitab*nkaps**2*ncplx*nkapn,0)

          allocate(s_sta%a(lscoff))
          allocate(s_sta%i(niax*nttab))
          allocate(s_sta%n(nbas+1))
          allocate(s_sta%s(nds**2*nitab*nkaps**2*ncplx*nkapn))
        endif

C   ... If only header is to be read, exit
        if (getdig(lio01,2,2) == 1) then
          iostr = .true.
          return
        endif

C        if (getdig(lio01,1,2) == 1)
C     .    call rx('iostr: illegal combination of bits in lio')

C   ... Read strux, block format
        if (getdig(lio23,1,2) == 0) then
          call rxx(getdig(lio23,0,2) == 1,'iostr not ready for complex')
C#ifdefC MASTERIO
C          if (procid == 0) then
C#endif
          nblk = nitab*nkap**2
          if (nmto == 1) nblk = nitab*nkap
          lerr = iostr2(ifi,nbas,nttab,nblk,lscoff,nds,nkap,
     .                  s_sta%n,s_sta%i,kap2,s_sta%a,s_sta%s)
C#ifdefC MASTERIO
C          endif
CC         sttime = MPI_WTIME()
CC         call info0(20,0,-1,' iostr: MPI broadcast strx data ...')
C          call mpibc1(lerr,1,0,.false.,'iostr','lerr')
C          call mpibc1(s_sta%n,nbas+1,2,.false.,'iostr','ntab')
C          call mpibc1(s_sta%i,niax*nttab,2,.false.,'iostr','iax')
C          call mpibc1(kap2,nkap,4,.false.,'iostr','kap2')
C          call mpibc1(s_sta%a,lscoff,4,.false.,'iostr','alp')
C          call mpibc1(s_sta%s,nds**2*nitab*nkap,4,.false.,'iostr','s')
CC         entime = MPI_WTIME()
CC         call info2(20,0,0,'  took %;1d sec',(entime-sttime),0)
C#endif
C   ... Read strux, vector format
        else
C#ifdefC MASTERIO
C          if (procid == 0) then
C#endif
          lerr = iostr4(ifi,ncplx,nl,nkap,kap2,s_sta%s)
        endif
C#ifdefC MASTERIO
C          call mpibc1(lerr,1,0,.false.,'iostr','lerr')
C          call mpibc1(s_sta%n,nbas+1,2,.false.,'iostr','ntab')
C          call mpibc1(s_sta%i,niax*nttab,2,.false.,'iostr','iax')
C          call mpibc1(kap2,nkap,4,.false.,'iostr','kap2')
C          call mpibc1(s_sta%s,nds**2*nitab*nkap*ncplx,4,.false.,'iostr',
C     .      's')
C        endif
C#endif
        if (.not. lerr) goto 99
        if (ipr >= 100) call awrit3(' IOSTR: read %i sites'//
     .    '%?#n#, %i inequivalent sites##',' ',80,
     .    lgunit(1),nttab,nttab-nitab,nitab)
      endif
      iostr = .true.
      return

C ... Error handling
   99 continue
      if ((lio/10) == 0) call rxs('IOSTR: file mismatch, file ',name)
      if (ipr >= 10)
     .  print *, 'IOSTR (warning): failed to read file ',name
      iostr = .false.

      end
      logical function iostr1(lio,name,itral,nttab,nitab,nl,nbas,
     .  nkap,lmaxw,ckbas,ifi)
C- Low level structure constants file open and header I/O
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Cio  lio   :1s digit 0 for read, otherwise write
Cio        :On file READ, returns with file contents of lio
Cio  name  :file name
Cio  itral :characterizes structure matrix transformations (mktral.f)
Cio  nttab :total number of pairs in neighbor and iax (pairs.f)
Cio  nitab :total number of inequivalent pairs in iax table (strscr.f)
Cio  nl    :global maximum l + 1
Cio  nbas  :number of atoms in the basis
Cio  nkap  :number of kinetic energies for which strux are calculated
Cio  lmaxw :maximum l for Watson-sphere
Cio  ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Co Outputs
Co   ifi   :file logical unit number
Co   iostr1:T if read successful, false otherwise
Cr Remarks
Cr   Reads screened strux header data from file, versions 0 and 1.
Cr
Cr   Structure of file(including header) version 2
Cr     record    content
Cr       1         -99  vsn  lio
Cr       2         itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
Cr   Structure of file(including header) version 1
Cr     record    content
Cr       1         -99  vsn  ihold
Cr       2         itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
Cr       ... end of header
Cr       3         ntab,iax,kap2,X,s
Cr                 where X=alpha or tral, depending on ihold
Cr
Cr   Structure of file(including header) version 0
Cr     record    content
Cr       1         itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
Cr       ... end of header
Cr       2         ntab,iax,kap2,alpha,s
Cu Updates
Cu   8 Jun 04 Implement version 1.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character name*(*)
      integer lio,itral,nttab,nitab,nl,nbas,nkap,lmaxw
      double precision ckbas
C ... Local parameters
      integer ifi,fopnx
      integer i,vsn,rdvsn
      character strn*120
      parameter (vsn=2)

      ifi = fopnx(name,100,16+8+4+0,-1)
      rewind ifi

C ... File READ
      if (mod(lio,2) == 0) then

C ... Check for version number.  If nonexistent, vsn=0
        read(ifi,end=99,err=99) i,rdvsn,lio
        if (i == -99 .and. rdvsn == 2) then
          read(ifi,end=99,err=99) itral,nttab,nitab,nl,nbas,nkap,lmaxw,
     .    ckbas
        elseif (i == -99 .and. rdvsn == 1) then
          read(ifi,end=99,err=99) itral,nttab,nitab,nl,nbas,nkap,lmaxw,
     .    ckbas
        else
          rdvsn = 0
          lio = 0
          backspace ifi
          read(ifi,end=99,err=99) itral,nttab,nitab,nl,nbas,nkap,lmaxw,
     .    ckbas
        endif

        strn = ' IOSTR : read file='//name//
     .    '%a lio=%i nkap=%i nl=%i nttab=%i'
        i = lio

C ... File WRITE
      else
C       This is version 0
C       write(ifi) itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
C       This is version 1
C        write(ifi) -99,vsn,ihold
C        write(ifi) itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
C       This is version 2
        write(ifi) -99,vsn,lio-1
        write(ifi) itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas

        strn = ' IOSTR : write file='//name//
     .    '%a lio=%i nkap=%i nl=%i nttab=%i'
        i = lio-1


      endif

      call info5(110,1,0,strn,i,nkap,nl,nttab,0)

      iostr1 = .true.
      return
   99 iostr1 = .false.
      end
      logical function iostr2(ifi,nbas,nttab,nblk,lscoff,nds,nkap,
     .  ntab,iax,kap2,scoffs,s)
C- Low level structure constants I/O
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   nbas  :number of atoms in the basis
Ci   nttab :total number of pairs in neighbor and iax (pairs.f)
Ci   nblk  :number of blocks of s to write
Ci   lscoff:size of scoffs matrix
Ci   nds   :leading dimensions of s
Ci   nkap  :number of kinetic energies for which strux are calculated
Cio Inputs/Outputs
Cio  ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairs.f)
Cio  iax   :array of parameters containing info about each pair
Cio  kap2  :interstitial kinetic energies
Cio  scoffs :tight-binding screening parameters
Cio  s     :real-space structure constant matrix
C ----------------------------------------------------------------------
      implicit none
      integer ifi,nbas,nttab,nblk,nds,nkap,niax,ntab(nbas+1)
      parameter (niax=10)
      integer iax(niax,nttab),lscoff
      double precision kap2(nkap),scoffs(lscoff),s(nds,nds,nblk)

      iostr2 = .true.
      if (ifi > 0) then
        read(ifi,end=99,err=99) ntab,iax,kap2,scoffs,s
      else
        write(-ifi) ntab,iax,kap2,scoffs,s
      endif
      return
   99 continue
      iostr2 = .false.
      end
      logical function iostr4(ifi,ncplx,ndim,nkap,kap2,s)
C- Low level structure constants I/O, vector format
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   ncplx :1 if real, 2 if complex
Ci   nl2   :leading dimensions of s
Ci   nkap  :number of kinetic energies for which strux are calculated
Cio Inputs/Outputs
Cio  kap2  :interstitial kinetic energies
Cio  s     :real-space structure constant matrix
C ----------------------------------------------------------------------
      implicit none
      integer ifi,ncplx,ndim,nkap,niax
      parameter (niax=10)
      double precision kap2(ncplx,nkap),s(ncplx,ndim,nkap)
      iostr4 = .true.
      if (ifi > 0) then
        read(ifi,end=99,err=99) kap2,s
      else
        write(-ifi) kap2,s
      endif
      return
   99 continue
      iostr4 = .false.
      end
      subroutine iostr5(lio,nkap,nkaps,nkapn,ncplx)
C- Return nkaps, nkapn, ncplx
      implicit none
      integer lio,nkap,nkaps,nkapn,ncplx
      integer lio23,lio45,nmto,getdig

      lio23 = mod(lio/100,100)
      lio45 = mod(lio/10000,100)
      nmto = getdig(lio45,0,2)
      nkaps = nkap
      nkapn = 1
      if (nmto == 1) then
        nkaps = 1
        nkapn = nkap
      endif
      ncplx = 1+getdig(lio23,0,2)
      end

      logical function iostr6(lio,name,nbas,nttab,nds,nbalp,nbstr)
C- Returns dimensioning parameters for structure constants arrays
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio   :specifies conditions for file i/o (condensed into digits)
Ci         :See routine iostr for description
Ci         :digits 100-1000 and 10000-100000 are used.
Ci   name  :file name holding strux
Co Outputs
Co   nbas  :number of sites used to generate pairs
Co         :npr should be dimensioned npr(nbas+1)
Co   nttab :Number of pairs in iax table
Co         :iax should be dimensioned iax(niax,nttab)
Co   nds   :leading dimension of alpha and s
Co   nbalp :number of blocks of alpha (blocks arranged by by site)
Co         :alpha is dimensioned alpha(nds,nks**2,nbas,nkn) where
Co         :nks = number of 2-kappa strux and
Co         :nkn = number of 1-kappa strux
Co         :iostr6 returns nbalp = nds*nks**2*nbas*nkn
Co   nbstr :number of blocks of strx s (blocks arranged by inequivalent pairs)
Co         :s is dimensioned s(nds,nds,nks,nks,nitab,nkn)
Co         :iostr6 returns nbstr = nds**2*nks**2*nitab*nkn*ncplx
Cb Bugs
Cb   Only lio 10s digit only 0 for MPI
Cr Remarks
Cr   Returns memory requirements for alpha,iax,s
Cu Updates
Cu   02 Nov 12 Adapted from iostr
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character name*(*)
      integer lio,nitab,nttab,nbalp,nbstr
C ... Local parameters
      character locnam*8
      logical iostr1,lerr
      integer nl,nbas,ifi,ipr,nkap,itral,lmaxw,ncplx,
     .  lio01,lio23,lio45,getdig,nds,lior,nks,nkn
      double precision ckbasr
C#ifdefC MASTERIO
C     integer procid
C#endif
C     double precision tiny
C     double precision sttime,entime,MPI_WTIME
C     parameter (tiny=1d-8,niax=10)

C --- Setup information ---
      call getpr(ipr)
      locnam = name
      lio01 = mod(lio/1,100)
      lio23 = mod(lio/100,100)
      lio45 = mod(lio/10000,100)
      ncplx = 1+getdig(lio23,0,2)

C   ... Open file, and read header
C#ifdefC MASTERIO
C        if (procid == 0) then
C#endif
      lior = 0
      lerr = .not. iostr1(lior,name,itral,nttab,nitab,nl,nbas,
     .  nkap,lmaxw,ckbasr,ifi)
      if (lerr) goto 99

C#ifdefC MASTERIO
C        endif
C
C#endif

C ... If failed to read header, give up
      if (lerr) goto 99

C     Use nds=nl**2 for now, until nds is passed
      nds = nl*nl
      call iostr5(lio,nkap,nks,nkn,ncplx)

C ... Return memory required for alpha, iax, s
C     nbalp = nds*nks**2*nbas*nkn
      nbalp = nks**2*nbas*nkn
C     nbstr = nds**2*nks**2*nitab*nkn*ncplx
      nbstr = nks**2*nitab*nkn*ncplx

      return

C ... Error handling
   99 continue
      if ((lio/10) == 0) call rxs('IOSTR: file mismatch, file ',name)
      if (ipr >= 10)
     .  print *, 'IOSTR (warning): failed to read file ',name
      iostr6 = .false.

      end
