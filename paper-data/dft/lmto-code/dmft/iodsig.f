      integer function iodsig(mode,nomg,beta,ndsig,omega,sigma,ifi)
C- I/O of DMFT sigma or delta
C ----------------------------------------------------------------------
Cio Structures
Ci Inputs
Ci   mode  :1s digit specifies DFMT code
Ci         :0 no code; nothing to do
Ci         :1 K Haule's CTQMC code
Ci         :10s digit (file READ) specifies how to handle non existent file
Ci         :0 Always return, whether file read is successful or not
Ci         :1 missing file: abort
Ci         :2 abort if passed nomg or ndsig does not match file
Ci         :3 combination of 1+2
Ci         :4 return dimensioning parameters nomg and ndsig only (file read only)
Ci         :5 combination of 1+4
Ci         :100s digit applies to omega or sigma after reading or before writing
Ci         :0 Do not scale omega or sigma
Ci         :1 scale omega or sigma from eV to Ry
Ci         :2 scale omg,sigma from Ry to eV
Ci         :Add 4 if to treat input sigma as nonmagnetic, and split spins
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/outputs
Cio  nomg  :number of frequencies.
Cio        :Input for file write, or if 10s digit mode = 2 or 3
Cio        :If 10s digit mode >=4, nomg is output
Cio  beta  :inverse temperature, eV^-1.  Used only on file read to check omega(1)
Cio        :0 => no check is made.
Cio  ndsig :number of DMFT channels
Cio        :Required for file write, or if 10s digit mode = 2 or 3
Cio        :If 10s digit mode >=4, ndsig is output
Cio  omega :List of frequencies.
Cio        :Required for file write.
Cio        :On input if beta>0, a the first fequency is checked
Cio        :against the lowest Matsubara frequency.
Cio  sigma :self-energy
Cio        :Read or written unless 10s digit mode >=4
Co Outputs
Co  iodsig :(file READ)
Co         :0  normal return
Co         :-1 could not read elements from file
Co         :-2 file mismatch
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   31 Mar 17 (MvS) Extra switches in mode, e.g. to return dimensioning parameters
Cu   27 Feb 16 (MvS) additions to work in MPI mode
Cu   19 Oct 15  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ifi,nomg,ndsig
      double precision beta,omega(nomg),sigma(nomg,2*ndsig)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk(:,:)
C ... Local parameters
      logical ltmp,lsplit
      integer mode0,mode1,i,l,rdm,jfi,nomgl,lx
      double precision eVscaleRy
      real(8), parameter :: ry2ev = 13.60569193d0, pi = 4*datan(1d0)
      integer procid,master
      procedure(logical) :: isanrg
      procedure(integer) :: mpipid,isw

C ... Setup
      procid = mpipid(1)
      master = 0
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      if (mode0 == 0) return
      if (mode0 == 0) call rx('no DMFT code specified')
      if (mode0 /= 1) call rx('DMFT code must be 1 for now')
      select case (mod(mod(mode/100,10),4))
        case (0); eVscaleRy = 1
        case (1); eVscaleRy = 1/ry2ev
        case (2); eVscaleRy = ry2ev
        case default; call rxi('iodsig: bad 100s digit mode',mode)
      end select
      lsplit = mod(mode/100,10) >= 4

C --- File READ ---
      if (ifi > 0) then
        allocate(wk(nomg,2*ndsig+1))
        if (mod(mode1/4,2) == 1) then
          call info0(30,0,-1,' Read dimensioning parameters for DMFT sigma ...')
        else
          call info2(20,0,0,' Read DMFT sigma from file%?#(n)#%j ...#.  Scale omega, sigma by %;4d',
     .      isw(eVscaleRy == 1),eVscaleRy)
        endif

C   ... Count the number of rows and columns in the file
        nomgl=0; l=0
        If (procid == master) then
          rewind ifi
          i = rdm(ifi,0,0,' ',wk,nomgl,l)
        endif
        call mpibc1(i,1,2,.false.,'','')
        call mpibc1(nomgl,1,2,.false.,'','')
        call mpibc1(l,1,2,.false.,'','')
        if (i < 0) then
          call info0(10,0,0,' iodsig (warning) missing or improper file ...')
          if (mod(mode1,2) == 1)call rx('failed to read self-energy file')
          iodsig = -1
          return
        endif
        l = (l-1)/2  ! Number of correlated channels (subtract one for first column; account for sigma being complex)

        iodsig = 0

C   ... Return if seeking dimensioning parameters only
        if (mod(mode1/4,2) == 1) then
          nomg = nomgl
          ndsig = l
          call info2(30,0,0,' found %i frequencies and %i channels',nomg,ndsig)
          return
        endif

        if (lsplit) call info0(30,0,-1,' spin-splitting nonmagnetic sigma ...')

C   ... Check match in  number of frequencies, number of channels in sigma
        ltmp = isanrg(nomgl,nomg,nomg,'file: IODSIG:','number of frequencies',mod(mode1/2,2) == 1)
        if (ltmp) iodsig = -1
        lx = l ; if (lsplit) lx = 2*l
        ltmp = isanrg(lx,ndsig,ndsig,'file: IODSIG:','number of correlated channels',mod(mode1/2,2) == 1)
        if (ltmp) iodsig = -1
        if (iodsig == -1) return

C   ... Read omega,sigma into single wk array; unpack wk into omega,sigma
        if (procid == master) then
          rewind ifi
          i = rdm(ifi,0,nomg*(2*l+1),' ',wk,nomg,2*l+1)
        endif
        call mpibc1(i,1,2,.false.,'','')
        if (i /= 1) call rx('failed to read self-energy file')
        call mpibc1(wk,size(wk),4,.false.,'','')
        omega(1:nomg) = wk(:,1)*eVscaleRy
        sigma(1:nomg,1:2*l) = wk(1:nomg,2:2*l+1)*eVscaleRy
        if (lsplit) sigma(1:nomg,2*l+1:4*l) = sigma(1:nomg,1:2*l)
        deallocate(wk)
C       Check first frequency against ctrl file's beta
        if (beta > 0) then
          call fsanrg(omega(1)/eVscaleRy,pi/beta,pi/beta,1d-5,
     .      ' IODSIG:','pi/beta',mod(mode1/2,2) == 1)
        endif

C --- File WRITE ---
      else
        if (lsplit) call rx('iodsig not ready for lsplit and write')
        jfi = -ifi
        if (procid == master) then
          rewind jfi
          allocate(wk(2*ndsig,1))
          do  i = 1, nomg
            wk(:,1) = sigma(i,1:2*ndsig) * eVscaleRy
C           write(jfi,333) omega(i) * eVscaleRy, wk
            call awrit1('%;14,8D%$',' ',10000,jfi,omega(i)*eVscaleRy)
            do  l = 1, ndsig
              if (sigma(i,2*l-1)*eVscaleRy > 999d0  .or.  sigma(i,2*l)*eVscaleRy > 999d0 .or.
     .           -sigma(i,2*l-1)*eVscaleRy > 99.9d0 .or. -sigma(i,2*l)*eVscaleRy > 99.9d0) then
C                print *, i,l,sigma(i,2*l-1:2*l)*eVscaleRy
C                call awrit1(' %2;12,7D',' ',10000,6,sigma(i,2*l-1:2*l)*eVscaleRy)
                call awrit1(' %2:1;11F%$',' ',10000,jfi,sigma(i,2*l-1:2*l)*eVscaleRy)
              else
                call awrit1(' %2;12,7D%$',' ',10000,jfi,sigma(i,2*l-1:2*l)*eVscaleRy)
              endif
            enddo
            write(jfi,"('')")
  333       format(f14.8,100(1x,2f12.7))
          enddo
          deallocate(wk)
        endif
        iodsig = 0
      endif

      end
