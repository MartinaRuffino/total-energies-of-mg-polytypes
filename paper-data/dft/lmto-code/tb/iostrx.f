      logical function iostrx(lio,name,nl,nbas,nkap,kap2,itral,ckbas,
     .  lmaxw,nitab,oalpha,oiax,onpr,os)
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
Cio  oiax  :offset for iax array
Cio  onpr  :offset for number of sites in each cluster
Cio  os    :offset for screened structure constant matrix
Cb Bugs
Cb   Only lio 10s digit only 0 for MPI
Cr Remarks
Cr   If file read, allocates memory for and reads alpha,iax,npr,s
Cr   Input parameters nl,nbas,bas,nkap,lmaxw,itral must match file
Cr   for iostr to complete file read.
Cu Updates
Cu   06 Aug 06 Redesigned to work with 2-kappa strux
Cu   8 Jun 04  (MPI) read from master node only and broadcast
Cu             (ccomp with MASTERIO.) (Warning: not checked.)
Cu   1 Aug 98  revised for 3rd generation lmto
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character name*(*)
      integer itral,lio,nitab,nl,nbas,nkap,lmaxw,oalpha,oiax,onpr,os
      double precision kap2(nkap),ckbas
C ... Local parameters
      character locnam*8
      logical iostr1,iostr2,iostr4,lerr,ltmp
      integer nlr,nbr,ifi,ipr,nkr,nttab,lgunit,itralr,lmaxwr,niax,ncplx,
     .  lio01,lio23,lio45,getdig,nds,lscoff,nmto,lior,nblk,nkaps,nkapn
      integer mpipid,procid,vsn
      double precision tiny,ckbasr
C     double precision sttime,entime,MPI_WTIME
      parameter (tiny=1d-8,niax=10,vsn=1)
      procedure(logical) :: isanrg
C ... Heap
      integer w(1)
      common /w/ w

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
        nttab = w(onpr+nbas)
        if (.not. iostr1(lio,name,itral,nttab,nitab,nl,nbas,
     .    nkap,lmaxw,ckbas,ifi)) goto 99
C   ... Write strux, block format
        if (getdig(lio23,1,2) == 0) then
        nblk = nitab*nkaps**2*nkapn
          call rxx(getdig(lio23,0,2) == 1,'iostr not ready for complex')
          if (.not. iostr2(-ifi,nbas,nttab,nblk,lscoff,nds,nkap,
     .        w(onpr),w(oiax),kap2,w(oalpha),w(os))) goto 99
C   ... Write strux, vector format
        else
          if (.not. iostr4(-ifi,ncplx,nl,nkap,kap2,w(os)))
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

        call sanrg(.true.,nlr,nl,nl,'file: IOSTR:','nl')
        call sanrg(.true.,nbr,nbas,nbas,'file: IOSTR:','nbas')
        call sanrg(.true.,nkr,nkap,nkap,'file: IOSTR:','nkap')
        call fsanrg(ckbasr,ckbas,ckbas,tiny,'file: IOSTR:','cksum',.true.)
C   ... This is only a warning, not fatal error
        ltmp = isanrg(lmaxwr,lmaxw,lmaxw,'file: IOSTR:','lmaxw',.false.)

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
          call defdr(oalpha,lscoff)
          call defi(oiax,niax*nttab)
          call defi(onpr,nbas+1)
          call defdr(os,nds**2*nitab*nkaps**2*ncplx*nkapn)
        endif

C   ... If only header is to be read, exit
        if (getdig(lio01,2,2) == 1) then
          iostrx = .true.
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
     .                  w(onpr),w(oiax),kap2,w(oalpha),w(os))
C#ifdefC MASTERIO
C          endif
CC         sttime = MPI_WTIME()
CC         call info0(20,0,-1,' iostr: MPI broadcast strx data ...')
C          call mpibc1(lerr,1,0,.false.,'iostr','lerr')
C          call mpibc1(w(onpr),nbas+1,2,.false.,'iostr','ntab')
C          call mpibc1(w(oiax),niax*nttab,2,.false.,'iostr','iax')
C          call mpibc1(kap2,nkap,4,.false.,'iostr','kap2')
C          call mpibc1(w(oalpha),lscoff,4,.false.,'iostr','alp')
C          call mpibc1(w(os),nds**2*nitab*nkap,4,.false.,'iostr','s')
CC         entime = MPI_WTIME()
CC         call info2(20,0,0,'  took %;1d sec',(entime-sttime),0)
C#endif
C   ... Read strux, vector format
        else
C#ifdefC MASTERIO
C          if (procid == 0) then
C#endif
          lerr = iostr4(ifi,ncplx,nl,nkap,kap2,w(os))
        endif
C#ifdefC MASTERIO
C          call mpibc1(lerr,1,0,.false.,'iostr','lerr')
C          call mpibc1(w(onpr),nbas+1,2,.false.,'iostr','ntab')
C          call mpibc1(w(oiax),niax*nttab,2,.false.,'iostr','iax')
C          call mpibc1(kap2,nkap,4,.false.,'iostr','kap2')
C          call mpibc1(w(os),nds**2*nitab*nkap*ncplx,4,.false.,'iostr',
C     .      's')
C        endif
C#endif
        if (.not. lerr) goto 99
        if (ipr >= 100) call awrit3(' IOSTR: read %i sites'//
     .    '%?#n#, %i inequivalent sites##',' ',80,
     .    lgunit(1),nttab,nttab-nitab,nitab)
      endif
      iostrx = .true.
      return

C ... Error handling
   99 continue
      if ((lio/10) == 0) call rxs('IOSTR: file mismatch, file ',name)
      if (ipr >= 10)
     .  print *, 'IOSTR (warning): failed to read file ',name
      iostrx = .false.

      end
