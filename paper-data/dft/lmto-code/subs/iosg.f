      logical function iosg(lio,name,ivl,nl,nbas,ckbas,nttabg,s_stag)
C- Real space value-laplacian structure constants file read or write
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
Ci                 16 file is formatted (not implemented)
Ci   name  :file name holding strux (eg STRG)
Ci   ivl   :specifies functions used to built the value-laplacian strux
Ci   nl    :1 + global-maximum-l
Ci          specifies the leading dimensions of sg
Ci   nbas  :number of atoms in the basis
Ci   ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Cio Inputs/Outputs:
Cio  nttabg:number of strux for the value-laplacian basis
Cio  oiax  :offset for iax array
Cio  onpr  :offset for the total number of sites in each cluster ib, ntab(ib+1)
Cio  onprg :offset for the number of sites in each value-laplacian cluster, ntabg(ib)
Cio  osg   :offset for the value-laplacian structure constant matrix
Cb Bugs
Cb   MPI is not implemented
Cb   All sites are assumed non equivalent
Cb   Program does not check if ehvl used to build struxs are same as those in ctrl.*
Cb   in the reading mode
Cr Remarks
Cr   Structure constant matrix sg(nl,nl,2,2,nttabg) is treated as as 1D
Cr   array sg(ndim) throughout.
Cr
Cr   If file read and 8th bit of lio is on, allocates memory for and reads
Cr   iax, ntab, ntabg, and sg
Cr
Cr   Input parameters nl,nbas, and bas must match file
Cr   for iosg to complete file read.
Cr
Cr   Content of file STRG (including header), version 2
Cr     record  content
Cr       1     -77 vsn ihold
Cr       2     ivl,nttab,nttabg,nl,nbas,ckbas
Cr     .... end of header
Cr       3     ntab,ntabg,iax,sg,-999
Cr
Cr   In the reading mode program stops if either of the following happens:
Cr       first number is not -77
Cr       version number (vsn) is not 2
Cr       if ivl, nl, nbas or chbas do not match those in the input file
Cr       last number is not -999 (to verify if size of sg is correct)
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   24 Jan 08 (S. Lozovoi) adapted from iostr.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character name*(*)
      integer lio,nttabg,ivl,nl,nbas
      double precision ckbas
C ... For structures
!      include 'structures.h'
      type(str_str0):: s_stag
C ... Local parameters
      character locnam*8
      logical iosg1,iosg2,lerr
      integer nlr,nbr,ivlr,ifi,ipr,nttab,lgunit,niax,
     .  lio01,getdig,nds,lior,ndim,ng,ig
      integer vsn
      double precision tiny,ckbasr
      parameter (tiny=1d-8, niax=10, vsn=1)

C --- Setup information ---
      call getpr(ipr)
      locnam = name
      lio01 = mod(lio/1,100)

C --- File write ---
      if (mod(lio,2) /= 0) then

C       Use nds=nl**2 for now, until nds is passed
        nds = nl*nl
        ndim = (2*nds)**2*nttabg

C   ... Write header
        nttab = s_stag%n(nbas+1)
        if (.not. iosg1(lio,name,ivl,nttab,nl,nbas,
     .    nttabg,ckbas,ifi)) goto 99
C   ... Write strux
        if (.not. iosg2(-ifi,nbas,nttab,ndim,
     .    s_stag%n,s_stag%ng,s_stag%i,s_stag%s)) goto 99

        if (ipr >= 40) call awrit3(' IOSG: wrote to file '
     .    //locnam//'%a %i sites'//
     .    ' out of %i total',
     .    ' ',80,lgunit(1),nttabg,nttab,0)
        call fclose(ifi)

C --- File read ---
      else
C   ... Open file, and read header
        lior = 0
        lerr = .not. iosg1(lior,name,ivlr,nttab,nlr,nbr,
     .    nttabg,ckbasr,ifi)

C   ... If failed to read header, give up
        if (lerr) goto 99

C       Copy file parameters to those in argument list
        if (getdig(lio01,1,2) == 1) then
          nl = nlr
          nbas = nbr
          ckbas = ckbasr
        endif

        call sanrg(.true.,ivlr,ivl,ivl,'file: IOSG:','ivl')
        call sanrg(.true.,nlr,nl,nl,'file: IOSG:','nl')
        call sanrg(.true.,nbr,nbas,nbas,'file: IOSG:','nbas')
        call fsanrg(ckbasr,ckbas,ckbas,tiny,'file: IOSG:','cksum',
     .    .true.)

C   ... If only header is to be read, exit
        if (getdig(lio01,2,2) == 1) then
          iosg = .true.
          return
        endif

C   ... Allocate memory for arrays
        nds = nl*nl
        ndim = (2*nds)**2*nttabg
        if (getdig(lio01,3,2) == 1) then
          call i1alloc(s_stag%i,niax*nttab,0)
          call i1alloc(s_stag%n,nbas+1,0)
          call i1alloc(s_stag%ng,nbas,0)
          call d1alloc(s_stag%s,ndim,0)
        endif

C   ... Read strux
        lerr = iosg2(ifi,nbas,nttab,ndim,
     .               s_stag%n,s_stag%ng,s_stag%i,s_stag%s)
        if (.not. lerr) goto 99

C ... Sanity check
        ng = 0
        do  ig = 1, nbas
          ng = ng + s_stag%ng(ig)
        enddo
        call sanrg(.true.,ng,nttabg,nttabg,'file: IOSG:','nttabg')

C ... Normal exit
        if (ipr >= 10) call awrit3(' IOSG: read %i sites'//
     .    ' out of %i total',' ',80,
     .    lgunit(1),nttabg,nttab,0)
      endif
      iosg = .true.
      return

C ... Error handling
   99 continue
      if ((lio/10) == 0) call rxs('IOSG: file mismatch, file ',name)
      if (ipr >= 10)
     .  print *, 'IOSG (warning): failed to read file ',name
      iosg = .false.

      end

C      logical function iosgx(lio,name,ivl,nl,nbas,ckbas,nttabg,
C     .  oiax,onpr,onprg,osg)
CC- Real space value-laplacian structure constants file read or write
CC ----------------------------------------------------------------------
CCi Inputs
CCi   lio   :specifies conditions for file i/o (condensed into digits)
CCi         :digit  specification
CCi         :1-10s   bits in 1-10s digit control program flow
CCi                  0 for file read
CCi                  1 for file write
CCi                  2 (read only) No header matching requirement;
CCi                    header parms returned into passed arguments
CCi                  4 read/write file header only and exit.
CCi                  8 Add this number to allocate arrays before reading
CCi                    or to release arrays after writing
CCi                 16 file is formatted (not implemented)
CCi   name  :file name holding strux (eg STRG)
CCi   ivl   :specifies functions used to built the value-laplacian strux
CCi   nl    :1 + global-maximum-l
CCi          specifies the leading dimensions of sg
CCi   nbas  :number of atoms in the basis
CCi   ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
CCio Inputs/Outputs:
CCio  nttabg:number of strux for the value-laplacian basis
CCio  oiax  :offset for iax array
CCio  onpr  :offset for the total number of sites in each cluster ib, ntab(ib+1)
CCio  onprg :offset for the number of sites in each value-laplacian cluster, ntabg(ib)
CCio  osg   :offset for the value-laplacian structure constant matrix
CCb Bugs
CCb   MPI is not implemented
CCb   All sites are assumed non equivalent
CCb   Program does not check if ehvl used to build struxs are same as those in ctrl.*
CCb   in the reading mode
CCr Remarks
CCr   Structure constant matrix sg(nl,nl,2,2,nttabg) is treated as as 1D
CCr   array sg(ndim) throughout.
CCr
CCr   If file read and 8th bit of lio is on, allocates memory for and reads
CCr   iax, ntab, ntabg, and sg
CCr
CCr   Input parameters nl,nbas, and bas must match file
CCr   for iosg to complete file read.
CCr
CCr   Content of file STRG (including header), version 2
CCr     record  content
CCr       1     -77 vsn ihold
CCr       2     ivl,nttab,nttabg,nl,nbas,ckbas
CCr     .... end of header
CCr       3     ntab,ntabg,iax,sg,-999
CCr
CCr   In the reading mode program stops if either of the following happens:
CCr       first number is not -77
CCr       version number (vsn) is not 2
CCr       if ivl, nl, nbas or chbas do not match those in the input file
CCr       last number is not -999 (to verify if size of sg is correct)
CCu Updates
CCu  24 Jan 08 (S. Lozovoi) adapted from iostr.f
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      character name*(*)
C      integer lio,nttabg,ivl,nl,nbas,oiax,onpr,onprg,osg
C      double precision ckbas
CC ... Local parameters
C      character locnam*8
C      logical iosg1,iosg2,lerr
C      integer nlr,nbr,ivlr,ifi,ipr,nttab,lgunit,niax,
C     .  lio01,getdig,nds,lior,ndim,ng,ig
C      integer vsn
C      double precision tiny,ckbasr
C      parameter (tiny=1d-8, niax=10, vsn=1)
C
CC --- Setup information ---
C      call getpr(ipr)
C      locnam = name
C      lio01 = mod(lio/1,100)
C
CC --- File write ---
C      if (mod(lio,2) /= 0) then
C
CC       Use nds=nl**2 for now, until nds is passed
C        nds = nl*nl
C        ndim = (2*nds)**2*nttabg
C
CC   ... Write header
C        nttab = w(onpr+nbas)
C        if (.not. iosg1(lio,name,ivl,nttab,nl,nbas,
C     .    nttabg,ckbas,ifi)) goto 99
CC   ... Write strux
C          if (.not. iosg2(-ifi,nbas,nttab,ndim,
C     .        w(onpr),w(onprg),w(oiax),w(osg))) goto 99
C
C        if (ipr >= 40) call awrit3(' IOSG: wrote to file '
C     .    //locnam//'%a %i sites'//
C     .    ' out of %i total',
C     .    ' ',80,lgunit(1),nttabg,nttab,0)
C        call fclose(ifi)
C
CC --- File read ---
C      else
CC   ... Open file, and read header
C        lior = 0
C        lerr = .not. iosg1(lior,name,ivlr,nttab,nlr,nbr,
C     .    nttabg,ckbasr,ifi)
C
CC   ... If failed to read header, give up
C        if (lerr) goto 99
C
CC       Copy file parameters to those in argument list
C        if (getdig(lio01,1,2) == 1) then
C          nl = nlr
C          nbas = nbr
C          ckbas = ckbasr
C        endif
C
C        call sanrg(.true.,ivlr,ivl,ivl,'file: IOSG:','ivl')
C        call sanrg(.true.,nlr,nl,nl,'file: IOSG:','nl')
C        call sanrg(.true.,nbr,nbas,nbas,'file: IOSG:','nbas')
C        call fsanrg(ckbasr,ckbas,ckbas,tiny,'file: IOSG:','cksum',
C     .    .true.)
C
CC   ... If only header is to be read, exit
C        if (getdig(lio01,2,2) == 1) then
C          iosgx = .true.
C          return
C        endif
C
CC   ... Allocate memory for arrays
C        nds = nl*nl
C        ndim = (2*nds)**2*nttabg
C        if (getdig(lio01,3,2) == 1) then
C          call defi(oiax,niax*nttab)
C          call defi(onpr,nbas+1)
C          call defi(onprg,nbas)
C          call defdr(osg,ndim)
C        endif
C
CC   ... Read strux
C        lerr = iosg2(ifi,nbas,nttab,ndim,
C     .                 w(onpr),w(onprg),w(oiax),w(osg))
C        if (.not. lerr) goto 99
C
CC ... Sanity check
C        ng = 0
C        do  ig = 1, nbas
C          ng = ng + w(onprg+ig-1)
C        enddo
C        call sanrg(.true.,ng,nttabg,nttabg,'file: IOSG:','nttabg')
C
CC ... Normal exit
C        if (ipr >= 10) call awrit3(' IOSG: read %i sites'//
C     .    ' out of %i total',' ',80,
C     .    lgunit(1),nttabg,nttab,0)
C      endif
C      iosgx = .true.
C      return
C
CC ... Error handling
C   99 continue
C      if ((lio/10) == 0) call rxs('IOSG: file mismatch, file ',name)
C      if (ipr >= 10)
C     .  print *, 'IOSG (warning): failed to read file ',name
C      iosgx = .false.
C
C      end

      logical function iosg1(lio,name,ivl,nttab,nl,nbas,nttabg,
     .  ckbas,ifi)
C- Low level structure constants file open and header I/O
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Cio  lio   :1s digit 0 for read, otherwise write
Cio        :On file READ, returns with file contents of lio
Cio  name  :file name
Ci   ivl   :specifies functions used to built the value-laplacian strux
Cio  nttab :total number of pairs in neighbor table iax (pairc.f)
Cio  nl    :global maximum l + 1
Cio  nbas  :number of atoms in the basis
Cio  nttabg:total number of pairs in the value-laplacian cluster (pairg.f)
Cio  ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Co Outputs
Co   ifi   :file logical unit number
Co   iosg1:T if read successful, false otherwise
Cr Remarks
Cr   Reads strux header data for the value-laplacian basis from file,
Cr   corresponds to iostr version 2.
Cu Updates
Cu   24 Jan 08 Adapted from iostr1 (iostr.f)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character name*(*)
      integer lio,ivl,nttab,nttabg,nl,nbas
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

C ... Check for version number.  If nonexistent, stop
        read(ifi,end=99,err=99) i,rdvsn,lio
        if (i /= -77) call rxs('IOSG: file mismatch, file ',name)
        if (rdvsn == 2) then
          read(ifi,end=99,err=99) ivl,nttab,nttabg,nl,nbas,ckbas
        else
C ... Only ver=2 is currently permitted. Stop
          call rxi('iosg: version is not supported. ver = ',rdvsn)
        endif

        strn = ' IOSG : read file='//name//
     .    '%a lio=%i nl=%i nbas=%i nttab=%i nttabg=%i'
        i = lio

C ... File WRITE
      else
        i = lio-1
        write(ifi) -77,vsn,i
        write(ifi) ivl,nttab,nttabg,nl,nbas,ckbas

        strn = ' IOSG : write file='//name//
     .    '%a lio=%i nl=%i nbas=%i nttab=%i nttabg=%i'

      endif

      call info5(110,1,0,strn,i,nl,nbas,nttab,nttabg)

      iosg1 = .true.
      return

   99 iosg1 = .false.
      end

      logical function iosg2(ifi,nbas,nttab,ndim,ntab,ntabg,
     .  iax,sg)
C- Low level structure constants I/O
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   nbas  :number of atoms in the basis
Ci   nttab :total number of pairs in neighbor table iax (pairc.f)
Ci   ndim  :size of sg
Cio Inputs/Outputs
Cio  ntab  :ntab(ib) no. pairs in neighbor table iax preceding ib (pairc.f)
Cio  ntabg :ntabg(ib) no. pairs for the value-laplacian basis in cluster ib (pairg.f)
Cio  iax   :array of parameters containing info about each pair
Cio  sg    :real-space structure constant matrix for the value-laplacian basis
Co Outputs
Co   iosg2:T if read successful, false otherwise
Cr Remarks
Cr   nctrl is placed at the end of the file to ensure that array dimensions match
Cu Updates
Cu   24 Jan 08 Adapted from iostr4 (iostr.f)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,nbas,nttab,ndim
      integer niax
      parameter (niax=10)
      integer iax(niax,nttab),ntab(nbas+1),ntabg(nbas)
      double precision sg(ndim)
C ... Local parameters
      integer nctrl,ncr
      parameter (nctrl=-999)

      if (ifi > 0) then
        read(ifi,end=99,err=99) ntab,ntabg,iax,sg,ncr
        if (ncr /= nctrl)
     .    call rx('iosg2(r): array dimension mismatch')
      else
        write(-ifi) ntab,ntabg,iax,sg,nctrl
      endif

      iosg2 = .true.
      return

   99 continue
      iosg2 = .false.
      end
