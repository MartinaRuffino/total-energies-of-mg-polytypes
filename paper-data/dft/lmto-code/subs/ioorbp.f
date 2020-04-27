      logical function ioorbp(mode,nkapi,is1,is2,s_spec,nkmax,ifi)
C- File I/O of parameters defining orbitals in basis envelope functions
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: name lmxb idxdn orbp lmxa p pz norp ntorb
Co     Stored:    norp ntorb p pz orbp idxdn lmxb
Co           :Case mode = 1:
Co           :sspec->rtab table of sm-H smoothing radii is I/O
Co           :sspec->etab table of sm-H energies is I/O
Co     Allocated: *
Cio    Elts Passed:pz name
Cio    Passed to: *
Ci Inputs
Ci   mode  :A compound of the following digits:
Ci          1    parameters consist of rsmh,eh in sspec->orbp
Ci               and rsmh2,eh2, if nkapi=2
Ci          10   if to determine idxdn locally
Ci          100  if basis is allowed to be only partially complete:
Ci                  basis for missing species are unchanged
Ci          1000s digit (file read)
Ci          : 1 or 3 Only 1-kappa RSMH,EH are read
Ci          : 2 or 4 Allow 2-kappa RSMH,EH, RSMH2,EH2
Ci          : 1 or 2 Parameters are left untouched for any species
Ci          :        in which RSMH,EH are already supplied
Ci          10000s digit read or write PZ
Ci          : 1 or 2 read/write PZ
Ci          : 1      Supplied nonzero values of PZ take precedence
Ci          :        over contents basis file
Ci          : 2      Supplied PZ are all overwritten
Ci          100000s digit read or write P
Ci          : 1 or 2 read/write P
Ci          : 1      Supplied values of P for l>lmxb take precedence
Ci          :        over contents basis file
Ci          : 2      Fractional part of P for l>last read
Ci          :        are generated from free-electron like formula
Ci   nkapi :number of envelope fns of one l-quantum number in the basis
Ci   is1..2:range of species to read/write basis parameters
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Co Inputs/Outputs
Co  nkmax   :(used in file read only): largest nkapi encountered in is1..is2
Co          :10's digit set to 1 if any PZ is read
Co Outputs
Co  ioorbp  :false: failed to read basis from file
Co          :true: successfully read basis from file
Cl Local variables
Cl   orbp   :orbp(:,1,1) = RSMH1  orbp(:,2,1) = EH1
Cl          :orbp(:,1,2) = RSMH2  orbp(:,2,2) = EH2
Cl          :orbp(:,1,3) = PZ     orbp(:,2,3) = PZ (read from this file)
Cl          :orbp(:,1,3) = RSMZ   orbp(:,2,3) = EHZ (set in elocp for extended LO)
Cr Remarks
Cb Bugs
Cb   This routine doesn't check for local orbitals
Cu Updates
Cu   05 Jul 17 species data may now occur in any order
Cu   04 Jan 17 When writing P in the spin pol case, write spin-averaged P
Cu             Also attempt to read RSMH,EH to one higher l than indicated by lmxb
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Dec 10 Keep pnu from sspec for l > values of file pnu
Cu   24 Aug 09 Some bug fixes
Cu   09 May 09 Reads/writes combinations of P, PZ, (RSMH,EH), (RSMH2,EH2)
Cu   28 Apr 09 Small bug fixes: writes RMSH,EH for all nonzero values
Cu   02 Sep 03 ioorbp can read local orbitals.
Cu             Miscellaneous other additions to read mode
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   24 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nkapi,is1,is2,ifi,nkmax
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      character spid*8,sout*256,prtstr*80,unused*256
      integer jfi,n0,parg,norp,ntorb,nkap0,nkape
      integer procid,master,mpipid
      parameter (n0=10,nkap0=4,nkape=2)
      integer is,ks,lp1,jp,ikap,lmxa,lmxb,lmxi,j1,j2,ls,i,j,lpz,ipr,
     .  it(n0),idxdn(n0,nkap0),nli,mode0,mode1,mode2,mode3,mode4,mode5,konfig
      double precision rtab(n0,nkap0),etab(n0,nkap0),orbp(n0,2,nkap0),
     .  dasum,pz(n0,2),pi
      logical,allocatable:: lread(:)
      procedure(logical) :: scat,rdstrn
      procedure(integer) :: nglob

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      mode4 = mod(mode/10000,10)
      mode5 = mod(mode/100000,10)

      if (mode0 > 1 .or. mode1 > 1 .or.
     .  mode2 > 1 .or. mode4 > 2)
     .  call rxi('ioorbp: unexpected value for mode, mode=',mode)

      lpz = 0
      ioorbp = .true.
      ls = len(sout)
      if (is2 < is1) return
      call getpr(ipr)
      pi = 4d0*datan(1d0)
      procid = mpipid(1)
      master = 0
      unused = ' '

C --- File READ ---
      if (ifi > 0) then
        allocate(lread(is1:is2)); lread = .false.
        nkmax = 1
        ioorbp = .false.
        if (mode0 == 1 .or. mode4 == 1 .or. mode5 == 1) then
          if (.not. scat(ifi,'BASIS:',':',.true.)) goto 15
          if (mode2 == 1) ioorbp = .true.
C     ... For each line in file, do
          do
   18     if (.not. rdstrn(ifi,sout,ls,.false.)) exit
C         Blank any text after '#'
          i = scan(sout,'#')
          if (i>0) sout(i:) = ' '
          call word(sout,1,j1,j2)
          if (j1 > j2) cycle

C     ... Identify species ks
          do  is = is1, is2
            ks = is
            spid = s_spec(ks)%name
            if (spid == sout(j1:j2)) exit
            if (is == is2) then
              ks = len_trim(unused)
              unused(ks+2:) = sout(j1:j2)
              goto 18
            endif
          enddo
          lread(ks) = .true.

          prtstr = ' ioorbp: read species '//spid//'x'
          call dpzero(orbp,n0*2*nkap0)
          s_spec(ks)%norp = 2
          s_spec(ks)%ntorb = n0*nkape
          call iinit(idxdn,n0*nkap0)
          lmxb = s_spec(ks)%lmxb
          idxdn = s_spec(ks)%idxdn
C         Case determine idxdn locally
          if (mode1 /= 0) call iinit(idxdn,n0*nkape)
          if (lmxb+1 > n0) call rx('increase size(idxdn) in s_spec')

C     --- Unpack orbp; if RSMH present, skip file read ---
          if (mode3 <= 2) then
            orbp = s_spec(ks)%orbp
            if (dasum(n0,orbp(1,1,1),1) /= 0) then
C              call info0(20,0,0,' ioorbp: RSMH+EH already supplied '//
C     .          'for species '//spid//' ... skipping file input')
              if (dasum(n0,orbp(1,1,2),1) /= 0) nkmax = 2
              goto 12
            endif
            if (mode3 == 0) goto 12 ! mode3=0 => always skip reading RSMH
          endif

C     --- Read first set of parameters RSMH,EH ---
          if (mode0 == 1) then
C         Try to read RSMH with one extra l
            jp = j2
            i = parg(' RSMH=',4,sout,jp,ls,', ',-3,lmxb+2,it,orbp(1,1,1))
C         Read RSMH with required l
          if (i < 0) then
          jp = j2
          i = parg(' RSMH=',4,sout,jp,ls,', ',-3,lmxb+1,it,orbp(1,1,1))
          endif
C         No smoothing radii read ... skip this species
          if (i == 0 .or. i == -1) goto 12
          prtstr = trim(prtstr) // ' RSMH,EH'
C         After next block, lmxi = true lmxb for this rsmh,eh
          lmxi = lmxb
C         Fewer rsm read than lmxb : l's beyond last rsm -> high
          if (i < 0) then
            do  j  = -i, lmxb+1
C             print *, j, orbp(j,1,1)
              orbp(j,1,1) = 0
              idxdn(j,1) = 4
            enddo
            lmxi=-i-2
          endif
C         For l's > lmxb, set idxdn as follows:
          do  lp1 = 1, lmxi+1
C           if rsm > 0 and idxdn "undetermined", idxdn => low
            if (orbp(lp1,1,1) > 0 .and. idxdn(lp1,1) == 0) idxdn(lp1,1) = 1
C           if rsm < 0, idxdn -> high
            if (orbp(lp1,1,1) < 0) idxdn(lp1,1) = 4
C           if rsm = 0, idxdn -> intermediate.  ? Why is this here?
C           if (orbp(lp1,1,1) == 0) idxdn(lp1,1) = max(idxdn(lp1,1),2)
C           if rsm = 0, idxdn -> high
            if (orbp(lp1,1,2) == 0) idxdn(lp1,2) = 4
          enddo
C         Try to read EH with one extra l
          jp = j2
          i = parg(' EH=',4,sout,jp,ls,', ',-3,lmxi+2,it,orbp(1,2,1))
C         Read EH with required l
          if (i < 0) then
            jp = j2
            i = parg(' EH=',4,sout,jp,ls,', ',-3,lmxi+1,it,orbp(1,2,1))
          endif
          if (i == 0 .or. i == -1) goto 12
C         if (i < lmxi+1) goto 12
C         Fewer eh read than lmxb : l's beyond last eh -> high
          if (i < 0) then
            do  j  = -i, lmxi+1
c             print *, j, orbp(j,1,1)
              orbp(j,1,1) = 0
              idxdn(j,1) = 4
            enddo
            lmxi = -i-2
          endif
          endif

C     --- Read second set of parameters RSMH2,EH2 ---
          if (mode0 == 1 .and. nkapi == 2) then
            jp = j2
            call dvset(orbp(1,1,2),1,n0,0d0)
C           Try to read RSMH2 with one extra l
            jp = j2
            i=parg(' RSMH2=',4,sout,jp,ls,', ',-3,lmxb+2,it,orbp(1,1,2))
C           Read RSMH2 with required l
            if (i < 0) then
              jp = j2
              i=parg(' RSMH2=',4,sout,jp,ls,', ',-3,lmxb+1,it,orbp(1,1,2))
            endif
            if (mode3 == 1 .or. mode3 == 3) then
              call dvset(orbp(1,1,2),1,n0,0d0) ! Only use 1-kappa
              i = 0
            endif
            if (i /= 0) then
              prtstr = trim(prtstr) // ' RSMH2,EH2'
C             lmxi = true lmxb for this rsmh,eh
              lmxi = lmxb
              if (i < 0) lmxi=-i-2
              call ivset(idxdn(1,2),lmxi+2,n0,4)
              do  lp1 = 1, lmxi+1
C               if rsm > 0, idxdn => low
                if (orbp(lp1,1,2) > 0) idxdn(lp1,2) = 1
C               if rsm < 0, idxdn -> high
                if (orbp(lp1,1,2) < 0) idxdn(lp1,2) = 4
C               if rsm = 0, idxdn -> high
C               if (orbp(lp1,1,2) == 0) idxdn(lp1,2) = max(idxdn(lp1,2),2)
                if (orbp(lp1,1,2) == 0) idxdn(lp1,2) = 4
              enddo
              jp = j2
C             Try to read EH2 with one extra l
              jp = j2
              i=parg(' EH2=',4,sout,jp,ls,', ',-3,lmxi+2,it,orbp(1,2,2))
C             Read EH2 with required l
              if (i < 0) then
                jp = j2
                i=parg(' EH2=',4,sout,jp,ls,', ',-3,lmxi+1,it,orbp(1,2,2))
              endif
              if (i < lmxi+1) goto 12
              if (dasum(n0,orbp(1,1,2),1) /= 0) nkmax = 2
            else
              call ivset(idxdn(1,2),1,n0,4)
            endif
          endif
   12     continue

C     --- Read P ---
          call dpzero(orbp(1,1,3),2*n0)
          if (mode5 == 1 .or. mode5 == 2) then
          lmxa = s_spec(ks)%lmxa
          orbp(:,:,3) = s_spec(ks)%p
C     ... This branch sets default values
          if (mode5 == 2) then
            do  lp1 = 1, lmxa+1
              konfig = orbp(lp1,1,3)
              orbp(lp1,1,3) = konfig + 0.5d0-datan(lp1-1-.5d0)/pi
            enddo
          endif
          jp = j2
          i = parg(' P=',4,sout,jp,ls,', ',-3,lmxa+1,it,pz)
C         If file has P- values, use them
          if (i /= 0) then
            prtstr = trim(prtstr) // ' P'
            lmxi = lmxa
            if (i < 0) lmxi=-i-2
            if (i < lmxa+1) lmxi = i-1
            do  lp1 = 1, lmxi+1
              orbp(lp1,1,3) = pz(lp1,1)
            enddo
            call dcopy(n0,orbp(1,1,3),1,orbp(1,2,3),1)
            s_spec(ks)%p(1:n0,1:2) = orbp(1:n0,1:2,3)
          endif
          endif
          call dpzero(orbp(1,1,3),n0)

C     --- Read and pack PZ ---
C  15     continue
          if (mode4 == 1 .or. mode4 == 2) then
            lmxa = s_spec(ks)%lmxa
            orbp(:,:,3) = s_spec(ks)%pz
            if (mode4 == 2) then
              call dpzero(orbp(1,1,3),n0)
              s_spec(ks)%pz(1:n0,1:2) = orbp(1:n0,1:2,3)
            endif
            jp = j2
            i = parg(' PZ=',4,sout,jp,ls,', ',-3,lmxb+1,it,pz)
            if (i /= 0) then
              prtstr = trim(prtstr) // ' PZ'
              lmxi = min(lmxb,i-1)
              if (i < 0) lmxi=-i-2
              do  lp1 = 1, lmxi+1
                if (mode4 == 1 .and. orbp(lp1,1,3) /= 0) cycle
                orbp(lp1,1,3) = pz(lp1,1)
              enddo
              call dcopy(n0,orbp(1,1,3),1,orbp(1,2,3),1)
              s_spec(ks)%pz(1:n0,1:2) = orbp(1:n0,1:2,3)
            endif
C            j1 = 0
C            j2 = 0
C            do  lp1 = 1, lmxb+1
C              if (idxdn(lp1,1) == 1) j1 = 1
C              if (idxdn(lp1,2) == 1) j2 = 1
C            enddo
          endif

C         Set lpz if a local orbital specified
          orbp(:,:,3) = s_spec(ks)%pz
          if (dasum(n0,orbp(1,1,3),1) /= 0) lpz = 1

C     ... Reduce lmxb
          do  lp1 = lmxb+1, 1, -1
            if (idxdn(lp1,1) /= 4 .or. idxdn(lp1,nkapi) /= 4) goto 14
            if (orbp(lp1,1,3) >= 10) goto 14
            lmxb = lmxb-1
          enddo
   14     continue

C     ... Repack (Note: idxdn will be remade in suidx)
          s_spec(ks)%orbp = orbp
          s_spec(ks)%idxdn = idxdn
          s_spec(ks)%lmxb = lmxb

C         Printout what was read for this species
          prtstr(31:31) = ' '
          if (prtstr(33:33) == ' ') prtstr(33:33) = '-'
          call info0(20,0,0,trim(prtstr))
        enddo  ! terminates when all lines have been read from file
        endif  ! if (mode0 == 1 .or. mode4 == 1 .or. mode5 == 1)
        ioorbp = .true.

C       Check for orbitals in species not read
   15   continue
        do  ks = is1, is2
          if (.not. lread(ks)) then
            orbp = s_spec(ks)%orbp
            if (dasum(n0,orbp(1,1,2),1) /= 0) nkmax = 2
            if (dasum(n0,s_spec(ks)%pz,1) /= 0) lpz = 1
            call info0(20,0,0,'%9fmissing entry for species '//trim(s_spec(ks)%name))
          endif
        enddo
        if (unused /= ' ') call info0(20,0,0,' ioorbp: species present but not used:'//trim(unused))
        if (lpz == 1) nkmax = nkmax + 10
        deallocate(lread)

C --- File WRITE ---
      else
        if (procid == master) then
        jfi = -ifi
        write(jfi,'(''BASIS:'')')
        if (mode0 == 1 .or. mode4 /= 0) then
          do  is = is1, is2
          spid = s_spec(is)%name
          s_spec(is)%norp = 2
          s_spec(is)%ntorb = n0*nkape
          norp = s_spec(is)%norp
          ntorb = s_spec(is)%ntorb
          orbp = s_spec(is)%orbp
          lmxb = s_spec(is)%lmxb
C         Sanity checks
          call sanrg(.true.,norp,2,2,'ioorbp:','norp')
          call sanrg(.true.,ntorb,n0*nkapi,n0*nkape,'ioorbp:','ntorb')
          sout = ' ' // spid

          if (mode0 == 1) then
          call dpzero(etab,n0*2)
          call dpzero(rtab,n0*2)
          do  ikap = 1, nkapi
            call dcopy(n0,orbp(1,1,ikap),1,rtab(1,ikap),1)
            call dcopy(n0,orbp(1,2,ikap),1,etab(1,ikap),1)
            nli = 0
            do  i = n0, 1, -1
              if (orbp(i,1,ikap) <= 0) cycle
              nli = i
              exit
            enddo
            if (nli > 0) then
              call awrit6(
     .        '%a RSMH%?#n==2#2##=%n:1;3d EH%?#n==2#2##=%n:1;3d',
     .        sout,ls,0,
     .        ikap,nli,rtab(1,ikap),
     .        ikap,nli,etab(1,ikap))
              lmxb = max(lmxb,nli-1)
            endif
          enddo
          endif

C     ... Write P to lmxa+1
          if (mode5 >= 1) then
            lmxa = s_spec(is)%lmxa
            pz = s_spec(is)%p
            lpz = 0
            do  lp1 = 1, lmxa+1
              if (pz(lp1,1) > 0) lpz = lp1
            enddo
C           Take the spin average
            if (nglob('nsp') == 2) then
              pz(1:lpz,1) = (pz(1:lpz,1)+pz(1:lpz,2))/2
            endif
            if (lpz > 0)
     .        call awrit2('%a P=%n:1;3d',sout,ls,0,lpz,pz)
          endif

C     ... Write PZ
          if (mod(mode4,2) == 1) then
            pz = s_spec(is)%pz
            if (nglob('nsp') == 2) then
              pz(1:lpz,1) = (pz(1:lpz,1)+pz(1:lpz,2))/2
            endif
            lpz = 0
            do  lp1 = 1, lmxb+1
              if (pz(lp1,1) > 0) lpz = lp1
            enddo
            if (lpz > 0)
     .        call awrit2('%a PZ=%n:1;4d',sout,ls,0,lpz,pz)
          endif

          call awrit0('%a',sout,-ls,-jfi)

          enddo
        endif
        endif
      endif
      end
