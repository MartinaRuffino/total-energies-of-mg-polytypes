      subroutine iosigh(optio,mode,nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,lshft1,lshft2,lshft3,ifi,eseavr)
C- Read/write header information in sigma file
C ----------------------------------------------------------------------
Ci Inputs
Ci   optio :1s digit (read only)
Ci         :0 read header data into
Ci            nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3, [eseavr]
Ci            eseavr is set to NULLI if not present
Ci         :1 File data must passed matched value for
Ci            nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
Ci            eseavr is also read, if it is present
Ci         :2 file mode must passed matched value
Ci         :These digits may be taken in any combination
Ci         :10s digit
Ci         :2 file is I/O in ascii format
Ci         :otherwise, file in binary format
Ci   ifi   :File logical unit: ifi>0 for read; ifi<0 for write
Cio Inputs/Outputs
Cio  mode  : switches with information about the contents of sigma
Cio        :... 1-10s digit, file contains:
Cio        :0    sigma in the orbital basis
Cio        :1    sigma in the LDA basis, and LDA eigenvalues
Cio        :2    sigma in the lda basis (high energy parts replaced),
Cio              and LDA eigenvalues
Cio        :3    sigma in the orbital basis (high energy parts replaced)
Cio        :LW4  LDA eigenvalues and eigenvectors
Cio        :LW5  LDA eigenvalues and eigenvectors
Cio        :LW6  LDA eigenvalues and eigenvectors
Cl         :LW7      eigenvalues and eigenvectors
Cio        :LW8      eigenvalues and eigenvectors
Cio        :LW9  LDA vxc in the orbital basis
Cio        :10   exchange sigma_x in the orbital basis
Cio        :20   U matrix in the orbital basis
Cio        :... 100s digit, file contains:
Cio        :0    sigma (or object) in k-space
Cio        :1    sigma (or object) in R.S., as inverse Bloch sum
Cio  nsp   :2 for spin-polarized case, otherwise 1
Cio  nspc  :2 for noncollinear case, otherwise 1
Cio  ndham :maximum dimension of hamiltonian
Cio  nlmto :dimension of LMTO part of hamiltonian
Cio  nk1..nk3 : number of divisions in 1st,2nd,3rd Qlat
Cio  nqp   :number of k-points in the file.
Cio        :Usually, not not necessarily, the number of irreducible k-points.
Cio  nqpx  :number of k-points in the file.
Cio        :The total number of qp written to file (read/written only if ivsn>3).
Cio  lshft1:1 if offset BZ in Qlat(1), otherwise 0
Cio  lshft2:1 if offset BZ in Qlat(2), otherwise 0
Cio  lshft3:1 if offset BZ in Qlat(3), otherwise 0
Cio  eseavr:average value of self-energy for high-lying states
Cio        :on file READ (optio 0 only)
Cio        :if it is present in the header, it is returned.
Cio        :If it is not present, it is returned as NULLI
Cl Local variables
Cr Remarks
Cr   Reads or writes header file for sigma, with sanity checks
Cr   header contents:
Cr    1.  header, 2 records
Cr        record 1: 0, version-number, modef  (missing in early version)
Cr        record 2: nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
Cr    2.  qp offset (version 2 and later)
Cr   Body contents depend on mode.
Cr    ... 3..6 are repeated for nqp points (nqp is contained in header)
Cr    foreach iq, do
Cr    foreach isp, do
Cr       recrd   mode    mode-dependent record
Cr    3.  qp     all      5    qp(1:3) ndimhx nev
Cr    4.  eval   1-8    1-4,8  eval(nev); else no eval written
Cr    5.  occnum  6       6    dmatk(ndimhx)
Cr    6.  sigm   all      0    sigm(nlmto,nlmto)
Cr                       1-6   sigma or evec(ndimhx,ndimhx)
Cr                       7,8   sigma or evec(ndimhx,nev)
Cr    enddo
Cr    enddo
Cu Updates
Cu   16 Mar 18 iosigh v3: reads/writes nspc,ndham
Cu   14 Nov 12 Workaround to using backspace, to work on pgf compiler
Cu   27 Jun 08 Revised mode, and input checking.  Old optio=1 -> optio=3
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,lshft1,lshft2,lshft3,ifi,optio,mode
      double precision eseavr(2)
C ... Local parameters
      logical ltmp,isanrg
      integer jfi,i1,i1x,i2x,i2,i3,i4,i5,i6,i6x,i7,i8,i9,ivsn,modef,optio0,lascii
C     double precision xx
      integer NULLI
      procedure(logical) :: iosigh2
      parameter (NULLI=-99999)

      optio0 = mod(optio,10)
      lascii = mod(optio/10,10)

      if (ifi < 0) then
        ivsn = 3 ; if (nqp < nqpx) ivsn = 4
        jfi = -ifi
        rewind jfi
C       Write 0 => file contains version number (versions 1 and later)
        if (lascii == 2) then
          write(jfi,"('# header version mode')")
          write(jfi,"(i6,2i7)") 0,ivsn,mode
          if (ivsn > 3) then
            call rx('iosigh not ready for ascii ivsn>3')
          endif
          write(jfi,"('# nsp nspc ndham nlmto   -nkabc-    nqp  -lshft-     eseavr ...')")
          write(jfi,"(i4,i5,2i6,3i4,i6,3i3,4x,1p,2e20.13)")
     .      nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,eseavr(1:nsp)
          write(jfi,"('# header end')")
        else
          if (ivsn /= 0) write(jfi) 0, ivsn, mode
          if (ivsn > 3) then
            write(jfi) nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,lshft1,lshft2,lshft3,eseavr(1:nsp)
          else
          write(jfi) nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,eseavr(1:nsp)
        endif
        endif
      else
        rewind ifi
C       Read version, checking whether header is pre- version 1
        if (lascii == 2) then
          read(ifi,*)
          read(ifi,*,err=99,end=99) i1,ivsn,modef
        else
          read(ifi,err=99,end=99) i1,ivsn,modef
        endif
C       File header is pre-version 1; set version to zero
        if (i1 /= 0) then
          ivsn = 0; modef = 0; rewind ifi
        endif

C       Assign mode or check with file correspondence
        if (mod(optio0/2,2) /= 0) then
C         ltmp=isanrg(modef,mode,mode,'iosigh:','file''s mode',.true.)
        else
          mode = modef
        endif

C   ... Read header parameters or check with file correspondence
        if (mod(optio0,2) /= 0) then
          if (.not. iosigh2(ivsn,lascii,i1,i1x,i2x,i2,i3,i4,i5,i6,i6x,i7,i8,i9,eseavr,ifi)) goto 99
C         ltmp=isanrg(i1,nsp,nsp,'iosigh:','file''s nsp',.true.)
C         ltmp=isanrg(i1x,nspc,nspc,'iosigh:','file''s nspc',.true.)
C         ltmp=isanrg(i2x,ndham,ndham,'iosigh:','file''s ndham',.true.)
C         ltmp=isanrg(i2,nlmto,nlmto,'iosigh:','file''s nlmto',.true.)
C         ltmp=isanrg(i3,nk1,nk1,'iosigh:','file''s nk1',.true.)
C         ltmp=isanrg(i4,nk2,nk2,'iosigh:','file''s nk2',.true.)
C         ltmp=isanrg(i5,nk3,nk3,'iosigh:','file''s nk3',.true.)
C         ltmp=isanrg(i6,nqp,nqp,'iosigh:','file''s nqp',.true.)
C         ltmp=isanrg(i7,lshft1,lshft1,'iosigh:','file''s lshft1',.true.)
C         ltmp=isanrg(i8,lshft2,lshft2,'iosigh:','file''s lshft2',.true.)
C         ltmp=isanrg(i9,lshft3,lshft3,'iosigh:','file''s lshft3',.true.)
        else
          if (.not. iosigh2(ivsn,lascii,nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,
     .      lshft1,lshft2,lshft3,eseavr,ifi)) goto 99
        endif
      endif
      return

   99 continue
      call rx('iosigh: failed to read header from file')
      end

      logical function iosigh2(ivsn,lascii,nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,
     .  lshft1,lshft2,lshft3,eseavr,ifi)
C- Read 2nd line of header file
C ----------------------------------------------------------------------
Ci Inputs
Ci  lascii :2, File I/O in ASCII mode
Ci   ivsn  :version number
Ci   ifi   :File logical unit: ifi>0 for read; ifi<0 for write
Co Outputs
Co  nsp   :2 for spin-polarized case, otherwise 1
Co  nspc  :2 for noncollinear case, otherwise 1
Co  ndham :maximum dimension of hamiltonian
Co  nlmto :dimension of LMTO hamiltonian
Co  nk1..nk3 : number of divisions in 1st,2nd,3rd Qlat
Co  nqp   :number of irr. k-points in file
Co  nqpx  :total number of k-points in file (read only if ivsn>3)
Co  lshft1:1 if offset BZ in Qlat(1), otherwise 0
Co  lshft2:1 if offset BZ in Qlat(2), otherwise 0
Co  lshft3:1 if offset BZ in Qlat(3), otherwise 0
Co  eseavr:average value of self-energy for high-lying states
Co        :if it is present in the header, it is returned.
Co        :If it is not present, it is returned as NULLI
Cl Local variables
Cr Remarks
Cr   This is a kernel called by iosigh
Cr   Versions prior to version 2 do not contain nspc or ndham
C    In this case, nspc is set to 1, ndham to nlmto
Cu Updates
Cu   14 Nov 12  Passes ivsn to work around using backspace after file read
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ivsn,lascii,nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,ifi,
     .  lshft1,lshft2,lshft3
      double precision eseavr(2)
C ... Local parameters
      integer k
      integer NULLI
      parameter (NULLI=-99999)

C     get spin index to determine the number of eseavr to read
      if (lascii == 2) then
        read(ifi,*,err=99,end=99)
        read(ifi,*,err=99,end=99) nsp
      else
        read(ifi,err=99,end=99) nsp
      endif

C     Read 2nd line of header
      iosigh2 = .true.
      backspace ifi
      eseavr(1:nsp) = NULLI
      if (ivsn > 4) then
        call rx('iosigh: file version too new for this code to read')
      elseif (ivsn > 3 .and.  lascii == 2) then
        read(ifi,*,err=99,end=99)
     .  k,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,lshft1,lshft2,lshft3,eseavr(1:nsp)
        read(ifi,*,err=99,end=99) ! Line marking end of header
      elseif (ivsn > 2 .and.  lascii == 2) then
        read(ifi,*,err=99,end=99)
     .  k,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,eseavr(1:nsp)
        read(ifi,*,err=99,end=99) ! Line marking end of header
      elseif (lascii == 2) then
        read(ifi,*,err=99,end=99)
     .  k,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,eseavr(1:nsp)
        read(ifi,*,err=99,end=99) ! Line marking end of header
        ndham = nlmto; nspc = 1
      elseif (ivsn > 3) then
        read(ifi,err=11,end=11)
     .    k,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,nqpx,lshft1,lshft2,lshft3,eseavr(1:nsp)
      elseif (ivsn > 2) then
        read(ifi,err=11,end=11)
     .    k,nspc,ndham,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,eseavr(1:nsp)
      else
        read(ifi,err=11,end=11)
     .  k,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,eseavr(1:nsp)
        ndham = nlmto; nspc = 1
      endif
      return

C     Attempt to read data without ndham,nspc or eseavr, for backward compatibility
   11 continue
      eseavr = NULLI
C     Avoid using backspace ... doesn't work on pgf compiler
      rewind ifi
      if (ivsn > 0) read(ifi) k
C     backspace ifi
      read(ifi,err=99,end=99) nsp,nlmto,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
      ndham = nlmto; nspc = 1
      return

   99 continue
      iosigh2 = .false.
      end
