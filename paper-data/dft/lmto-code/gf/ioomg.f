      subroutine rdomg(mode,s_site,nbas,nzp)
C- Read Omega files for all CPA sites and energy points
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb omg domg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode   : 0 read from file 'omega'//ib (site ib) into s_site(ib)%omg
Ci          : 1 read from file 'domega'//ib (site ib) into s_site(ib)%domg
Ci          : Add 10 if to require omega be in alpha-repsn
Ci          : Add 20 if to require omega be in gamma-repsn
Ci          : Add 30 if to require omega be in spin-avgd gamma-repsn
Ci          : Add 100 if 2x2 matrix (relativistic)
Ci   nbas  :size of basis
Ci   nzp   :number of energy points
Co Outputs
Co   omega for CPA site ib read from file 'om-hr'//ib or 'dom-hr'//ib
Co   and stored in s_site(ib)%omg
Cl Local variables
Cu Updates
Cu   19 Aug 13 (Belashchenko) Added relativistic case (2 x 2)
Cu   05 Jan 13 Routine cleaned up
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nzp
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      logical lok
      integer ncomp,norb,ib,ioomg,nread,nok,i,mod0,mod1,mod2,nspc
      integer mpipid,procid,master
      complex(8), pointer :: omega(:,:)

      procid = mpipid(1)
      master = 0
      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)
      mod2 = mod(mode/100,10)
      nspc = 1
      if (mod2 == 1) nspc = 2

C ... Count number of CPA sites
      nread = 0
      do  ib = 1, nbas
        ncomp = s_site(ib)%ncomp
        if (ncomp < 2) cycle
        nread = nread+1
      enddo
      call info2(20,0,0,' RDOMG : read omega file%?#(n==1)##s#%-1j, '//
     .  '%i%-1j site%?#(n==1)##s#, %i energy points',nread,nzp)

C --- Read Omega for site ib from file.  Unsuccessful read = Omega=0 ---
      nok = 0
      do  ib = 1, nbas
        lok = .true.
        ncomp = s_site(ib)%ncomp
        if (ncomp < 2) cycle
        norb = s_site(ib)%norb
        if (mod0 == 0) then
          omega => s_site(ib)%omg
        elseif (mod0 == 1) then
          omega => s_site(ib)%domg
        endif
c       do  izp = 1, nzp
        if (procid == master) then
          i = ioomg(100*mod1+10*mod0+1,nzp,norb,nspc,omega,ib)
          if (i /= nzp) lok = .false.
        endif
        call mpibc1(omega,norb*norb*2*nspc*nzp,6,.false.,'rdomg',
     .              'omega')
c       enddo
        if (lok) nok = nok+1
      enddo
      if (nok /= nread) then
        call info2(20,0,0,'         read unsuccessful for '//
     .    '%i site%-1j%?#(n==1)##s#; set Omega=0',nread-nok,0)
      endif

      end

      subroutine womg(mode,s_site,ib,nzp)
C- Write omega file to disk for site ib
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  omg domg norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode   : 0 write to file 'omega'//ib (site ib) s_site(ib)%omg
Ci          :   'om-hr' is also written (human-readable format)
Ci          : 1 write to file 'domega'//ib (site ib) s_site(ib)%domg
Ci          :   'dom-hr' is also written (human-readable format)
Ci          : Add 10 if to flag omega is in alpha-repsn
Ci          : Add 20 if to flag omega is in gamma-repsn
Ci          : Add 30 if to flag omega is in spin-avgd gamma-repsn
Ci          : Add 100 if 2x2 matrix (relativistic)
Ci   nbas  :size of basis
Ci   nzp   :number of energy points
Co Outputs
Co   omega for CPA site ib written to files:
Co      'omega'//ib or 'domega'//ib and also:
Co      'om-hr'//ib or 'dom-hr'//ib
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 Jan 13 Routine cleaned up
C ----------------------------------------------------------------------
C ... Passed parameters
      use structures
      implicit none
      integer mode,ib,nzp
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer norb,iomod,i,ioomg,mod1,mod2,nspc
      complex(8), pointer :: omega(:,:)

      mod1 = mod(mode/10,10)
      mod2 = mod(mode/100,10)
      nspc = 1
      if (mod2 == 1) nspc = 2
      if (mod(mode,10) == 0) then
        call info2(20,0,0,' WOMG : write omega files '//
     .    'for site %i, %i energy points',ib,nzp)
        iomod = 2
        omega => s_site(ib)%omg
      elseif (mod(mode,10) == 1) then
        call info2(20,0,0,' WOMG : write domega files '//
     .    'for site %i, %i energy points',ib,nzp)
        iomod = 12
        omega => s_site(ib)%domg
      else
        call rx1('WOMG: invalid mode = %i',mode)
      endif
      iomod = iomod + 100*mod1

      norb = s_site(ib)%norb

      i = ioomg(iomod,nzp,norb,nspc,omega,ib)
      end

      integer function ioomg(mode,izp,norb,nspc,omega,ib)
C- Read/write CPA coherent interactor Omega for site ib
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode   : 1s digit
Ci          : 0 read Omega for 1 z-point
Ci          : 1 read izp z-points (not implemented)
Ci          : 2 write izp z-points.  Both omega and om-hr are written
Ci          : 10s digit
Ci          : 0 filename is  'omega'//ib for site ib
Ci          : 1 filename is 'domega'//ib for site ib
Ci          : 100s digit
Ci          : 2 omega generated for gamma respsn
Ci          : 3 omega generated for spin-avereged gamma respsn
Ci          : 1000s digit
Ci          : 0 two separate spin components
Ci          : 1 2x2 matrix (relativistic)
Ci   ib     : Read omega for site ib
Ci   izp    : read/write omega for energy points iz=1...izp
Ci   norb   : number of orbitals at site ib
Ci   nspc   : number of coupled spins
Cio Inputs/Outputs
Ci   omega  : read from, or written to disk
Co  Outputs
Co   ioomg  :>0  number of energies read or written; read successful
Co          :<=0 number of energies successfully read or written
Cl Local variables
Cl   omfile : name of file to read or write (depends on 10s digit mode)
Cl          : File READ:  'omega'//ib or 'domega'//ib
Cl          : File WRITE: 'om-hr'//ib or 'dom-hr'//ib
Cu Updates
Cu   05 Jan 13 Routine cleaned up
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,izp,norb,ib,nspc
      double complex omega(norb,nspc,norb,2,*)
C ... Local parameters
      integer ifi,ifi2,fopn,i,j,k,n1,n2,iz,iz1,mod0,mod1,mod2,nglob,stdo
      integer fm,rdm
      double precision wk(2,norb,norb,2)
      logical success,lrew,known,lold,l22
      character *10 categ,cat2,omfile,omfil2
      character *4  azp,azp1
      logical scat

      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)
      mod2 = mod(mode/100,10)
      call sanrg(.true.,mod0,0,2,'ioomg:','mod0')
      stdo = nglob('stdo')
      ioomg = 0
      success = .false.

      i = 0
      call bin2a('',0,0,ib,2,0,4,azp1,i)

C ... High-resolution file name ('omega'//ib), logical unit
      if (mod1 == 0) then
        omfile = 'omega'//azp1(1:i)
      elseif (mod1 == 1) then
        omfile = 'domega'//azp1(1:i)
      else
        call rx1('IOOMG: invalid mode: %i',mode)
      endif
      ifi = fopn(omfile)
      rewind ifi
      if (izp == 1) then
        if (mod2 /= 0) then
          success = scat(ifi,'REPSN:',':',.false.)
          if (.not. success)
     .      call rx('ioomg: file missing REPSN category')
          backspace ifi
          read(ifi,*) categ,i
          call sanrg(.true.,i,mod2,mod2,'ioomg','REPSN')
        endif
      endif

C ... Case WRITE: Human-readable name ('om-hr'//ib), logical unit
      if (mod0 == 2) then
        if (mod1 == 0) then
          omfil2 = 'om-hr'//azp1(1:i)
        elseif (mod1 == 1) then
          omfil2 = 'dom-hr'//azp1(1:i)
        endif
        ifi2 = fopn(omfil2)
        rewind ifi2
      endif

C --- File READ omega(:,:,:,:,izp)
      if (mod0 == 0 .or. mod0 == 1) then
        known = .false. ; lrew = .true. ; l22 = .true.
        iz1 = izp ; if (mod0 == 1) iz1 = 1
        do  iz = iz1, izp
          i = 0
          call bin2a('',0,0,iz,2,0,4,azp,i)
          categ = 'ZP'//azp(1:i)//':'
          cat2  = ' ZP'//azp(1:i)//':'  ! OLD style for backward compatibility

          if (.not. known) then         ! If first point, establish new or old category format
            success = scat(ifi,categ,':',.true.)
            if (success) then
              lold = .false. ; known = .true.
            else
              success = scat(ifi,cat2,':',.true.)
              if (success) then
                lold = .true. ; known = .true.
              endif
            endif
          else                          ! Format known, read without rewinding
            if (lold) then
              success = scat(ifi,cat2 ,':',.false.)
            else
              success = scat(ifi,categ,':',.false.)
            endif
          endif
          if (success) then
            i = rdm(ifi,10001,2*norb**2,' ',wk,2*norb,norb)  ! Omega(:,1,:,1,iz)
            if (i == 1) then
              i = rdm(ifi,10001,2*norb**2,' ',wk(1,1,1,2),2*norb,norb)  ! Omega(:,nspc,:,2,iz)
            endif
            if (i == 1) then
              do  i = 1, norb
                do  j = 1, norb
                  omega(j,1,i,1,iz) = dcmplx(wk(1,i,j,1),wk(2,i,j,1))
                  omega(j,nspc,i,2,iz) = dcmplx(wk(1,i,j,2),wk(2,i,j,2))
                enddo
              enddo
              success = .true.
            else
              success = .false.
            endif
          endif
          if (success) then
            ioomg = ioomg + 1
          else
            call dpzero(omega(1,1,1,1,iz),norb*nspc*norb*2*2)
          endif
          if (success .and. nspc == 2 .and. l22) then
            i = rdm(ifi,10001,2*norb**2,' ',wk,2*norb,norb)  ! Omega(:,1,:,2,iz)
            if (i == 1) then
              i = rdm(ifi,10001,2*norb**2,' ',wk(1,1,1,2),2*norb,norb)  ! Omega(:,2,:,1,iz)
            endif
            if (i == 1) then
              do  i = 1, norb
                do  j = 1, norb
                  omega(j,1,i,2,iz) = dcmplx(wk(1,i,j,1),wk(2,i,j,1))
                  omega(j,2,i,1,iz) = dcmplx(wk(1,i,j,2),wk(2,i,j,2))
                enddo
              enddo
            else
              omega(:,1,:,2,iz) = dcmplx(0d0,0d0)
              omega(:,2,:,1,iz) = dcmplx(0d0,0d0)
              l22 = .false. ; known = .false.
            endif
          endif
        enddo ! iz
      endif

C --- File WRITE Omega(:,:,:,1:izp)
      if (mod0 == 2) then
        if (mod2 /= 0) then
          write(ifi,"('REPSN:',i3)") mod2
        endif
        do  iz = 1, izp
          i = 0
          call bin2a('',0,0,iz,2,0,4,azp,i)
          cat2 = 'ZP'//azp(1:i)//':'
          write(ifi,201) cat2
          write(ifi2,201) cat2
          assign 205 to fm
          if (norb > 9) assign 204 to fm

c         if (norb > 9) then
          write(ifi2,fm) (((omega(n1,min(k,nspc),n2,k,iz),n2=1,norb),
     .      n1=1,norb),k=1,2)
          if (nspc == 2) then
            write(ifi2,fm) (((omega(n1,k,n2,3-k,iz),n2=1,norb),n1=1,
     .        norb),k=1,2)
          endif
c         else
c           write(ifi2,205) (((omega(n1,min(k,nspc),n2,k,iz),n2=1,norb),
c    .        n1=1,norb),k=1,2)
c           if (nspc == 2) then
c             write(ifi2,205) (((omega(n1,k,n2,3-k,iz),n2=1,norb),n1=1,
c    .          norb),k=1,2)
c           endif
c         endif
          do  k = 1, 2
            write(ifi,202) ((omega(n1,min(k,nspc),n2,k,iz),n2=1,norb),
     .        n1=1,norb)
            write(ifi,201)
          enddo
          if (nspc == 2) then
            do  k = 1, 2
              write(ifi,202) ((omega(n1,k,n2,3-k,iz),n2=1,norb),
     .          n1=1,norb)
              write(ifi,201)
            enddo
          endif
        enddo
        success = .true.
        call fclr(omfil2,ifi2)
        ioomg = izp
      endif
      call fclr(omfile,ifi)
      return

C --- Error exit ---
    8 call rx('Error reading omega file')

C 200 format(1000(20f34.28/))
  201 format(a)
  202 format(18f18.12)
C  203 format(32(32f18.12/))
  204 format(2(16(32f8.4/)/))
  205 format(2(9(18f10.6/)/))
  801 format(' Omega file has no record for izp=',i3,
     .       ': set to zero')
      end

      subroutine cpomg(nbas,nspc,nzp,s_site)
C- Copy Omega to domg
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb omg
Co     Stored:     domg
Co     Allocated:  *
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   nspc  :number of coupled spins
Ci   nzp
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nspc,nzp,nthet
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,norb,len

      do  ib = 1, nbas
        nthet = s_site(ib)%ncomp
        if (nthet < 2) cycle
        norb = s_site(ib)%norb
        len = norb*norb*nspc*2*nzp
        call dpcopy(s_site(ib)%omg,s_site(ib)%domg,1,len*2,1d0)
      enddo

      end
