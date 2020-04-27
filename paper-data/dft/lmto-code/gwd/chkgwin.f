      subroutine chkgwin(s_site,s_spec,s_lat,s_gw,s_bz,nbas)
C- Performs some checks on the GWinput file
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu pz class
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa z p pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_gw  : not used now
Cio  s_bz  : not used now
Ci Inputs
Ci   nbas  :size of basis
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   02 May 13 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,n0,nkap0
      parameter (n0=10,nkap0=4)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(nbas)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_gw)::    s_gw
      type(str_bz)::    s_bz
C ... Local parameters
      logical ltmp,rdstrn,lerr,a2bin,lwarn,lerri,lwarni
      integer a2vec,fopng,i,iat,ibc,ibs,ic,iclbsj,ifii,ipr,is,j,j1,j2,k,
     .  konf,l,lmaxa,nbas0,ncore,nglob,nl,nncx,nnvv(n0,nbas),npqn,nsp
      integer iwk(300),ix(10),konfig(0:n0),lcutmx(nbas)
      integer,allocatable:: iocc(:,:,:,:),ips(:),ipc(:)
      double precision pnu(n0,2),pnz(n0,2),xx,vol,vol0,qsafe,QpGcut_cou,QpGcut_psi
      integer ntok,strnsiz,NULLI
      parameter (ntok=26,strnsiz=16,NULLI=-99999)
      character*(ntok*strnsiz) toklst
      character strn*512

C --- Setup ---
      call getpr(ipr)
      vol = s_lat%vol
      nsp  = nglob('nsp')
      nl   = nglob('nl')
      QpGcut_cou = NULLI; QpGcut_psi = NULLI
      toklst =
     .  '<PRODUCT_BASIS> '//
     .  'QpGcut_cou      '//
     .  'QpGcut_psi      '//
     .  'Verbose         '//
     .  'dw              '//
     .  'esmr            '//
     .  'iSigMode        '//
     .  'n1n2n3          '//
     .  'omg_c           '

C ... Count the number of atoms, excluding floating orbitals; check lmxa
      call info0(2,1,-1,' ... Checking for uniformity in lmxa ...')
      nbas0 = 0
      j = 0
      do  i = 1, nbas
        is = s_site(i)%spec
        lcutmx(i) = s_spec(is)%lmxa
        if (lcutmx(i) == -1) cycle  ! floating orbital
        if (lcutmx(i) /= lcutmx(1)) j = i
        if (s_spec(is)%z == 0) cycle
        nbas0 = nbas0 + 1
      enddo
      if (j /= 0) then
      call info5(2,1,1,' Encountered lmxa=%i for site %i but lmxa=%i'//
     .  ' for site 1.%N GW code does not yet run properly unless all '
     .  //'lmxa are identical.',lcutmx(j),i,lcutmx(1),0,0)
      call rx('modify ctrl file before proceeding')
      endif
      call info0(2,0,0,' ok')

      call info0(2,1,1,' ... Checking header section of GWinput file')

C ... Read first part of GWinput
      ltmp = .false.
      ltmp = .true.
      ifii = fopng('GWinput',-1,0)
      rewind ifii
      do
      if (.not. rdstrn(ifii,strn,len(strn),.false.)) exit
      call word(strn,1,j1,j2)
      if (j1 > j2) cycle
      if (strn(j1:j1) == '!') cycle
C     if (strn(j1:j1) >= '0' .and. strn(j1:j1) <= '9') cycle
      k = index(toklst,strn(j1:j2))
      if (k == 0) cycle
      if (mod(k-1,16) /= 0) cycle
      i = (k-1)/strnsiz

      select case (i)

        case (0 )  ! <PRODUCT_BASIS>
          call info0(2,1,1,
     .      ' ... Checking product basis section of GWinput file')
          if (.not. rdstrn(ifii,strn,len(strn),.false.)) goto 999
          if (.not. rdstrn(ifii,strn,len(strn),.false.)) goto 999 !Tolerance
          lerr = .true.
          call word(strn,1,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' PB tolerance has invalid argument')
          else
            qsafe = 3e-4
            if (xx < qsafe) then
              call info2(2,0,0,' PB tol=%d ... should be safe',
     .          xx,qsafe)
            else
              call info2(2,0,0,'*PB tol=%d ... monitor convergence '//
     .        '(sp parts especially affected.  %d is safer.)',xx,qsafe)
            endif
          endif
          if (.not. rdstrn(ifii,strn,len(strn),.false.)) goto 999
          if (.not. rdstrn(ifii,strn,len(strn),.false.)) goto 999 !lcutmx
          iwk(1:nbas) = lcutmx(1:nbas)
          j = 0
          j = a2vec(strn,len(strn),j,2,', ',2,-3,nbas0,ix,lcutmx)
          if (j /= nbas0) then
            call info2(2,0,0,'*Warning: expected %i values for lcutmx '
     .        //'but read only %i',nbas0,j)
            if (j < 0) call rx('fix GWinput file')
          else
            call info0(2,1,-1,
     .        ' ... Checking consistency between lcutmx,lmxa')
            iat = 0; lerr = .false.
            do  i = 1, nbas
              is = s_site(i)%spec
              lmaxa = s_spec(is)%lmxa
              if (lmaxa > -1) then
                iat = iat + 1
                if (lcutmx(iat) > lmaxa) then
                  call errnl(lerr)
                  call info5(2,0,0,' site %,3i: lcutmx=%i but lmxa=%i',
     .              i,lcutmx(iat),lmaxa,0,0)
                endif
              endif
            enddo
          endif
          if (lerr) then
          call info0(2,0,0,'*apparent inconsistency in lcutmx')
          else
            call info0(2,0,0,' ... ok')
          endif

          call info0(2,1,-1,' ... Checking radial wave functions table')
          if (.not. rdstrn(ifii,strn,len(strn),.false.)) goto 999

          do  i = 1, nbas
            is = s_site(i)%spec
            pnu = s_site(i)%pnu
            pnz = s_site(i)%pz
            lmaxa = s_spec(is)%lmxa
            if (lmaxa == -1) cycle
            call atwf(0,xx,lmaxa,xx,nsp,pnu,pnz,xx,xx,xx,xx,xx,k,
     .        ncore,konfig,xx,xx,xx)
            do  l = 0, lmaxa
              read(ifii,*) iwk(1:4)
              nnvv(l+1,i) = 2
              if (pnz(l+1,1) /= 0) nnvv(l+1,i) = 3
              konf = pnu(l+1,1)
              if (pnz(l+1,1) /= 0)
     .          konf = min(int(pnu(l+1,1)),int(mod(pnz(l+1,1),10d0)))
              nncx = konf - l - 1
              if (iwk(1) /= i .or. iwk(2) /= l .or.
     .          iwk(3) /= nnvv(l+1,i) .or. iwk(4) /= nncx) then
                print *, nnvv(l+1,i),nncx
                call errnl(lerr)
                call info5(2,0,0,
     .            '*Error: Expected ib,l,nnvv,nncc= %i %i %i %i '//
     .            'but encountered %4:1,1i',i,l,nnvv(l+1,i),nncx,iwk)
              endif
            enddo
            do  l = lmaxa+1, nl-1
              read(ifii,*) iwk(1:4)
              call rx('not finished')
C     write(jfi,'(4i5)') i,l,0,0
            enddo
          enddo
          if (lerr) then
            call rx('GWinput is not consistent with ctrl file')
          else
            call info0(2,0,0,' ... ok')
          endif

          call info0(2,1,-1,' ... Checking valence product basis table')
          if (.not. rdstrn(ifii,strn,len(strn),.false.)) goto 999
          call word(strn,1,j1,j2)
          if (strn(j1:j2) /= 'atom') call rx(
     .      ' Line should begin with "atom" but reads '//strn(j1:j2))
          lerr = .false.
          allocate(iocc(2,nkap0,0:nl-1,nbas)); iocc = NULLI
          allocate(ips(nbas),ipc(nbas))
          do  i = 1, nbas
            is = s_site(i)%spec
            ips(i) = is
            ipc(i) = s_site(i)%class
            lmaxa = s_spec(is)%lmxa
            if (lmaxa == -1) cycle
            pnu = s_spec(is)%p
            pnz = s_spec(is)%pz
C           read ib,l,radial pqn k, and iocc,uocc. 1st 3 must match ctrl
            do  l = 0, lmaxa
              npqn = 2
              if (pnz(l+1,1) /= 0) npqn = 3
              do  k = 1, npqn

C               Read table from disk; save occ list for further checks
                read(ifii,*) iwk(1:5)
                iocc(1,k,l,i) = iwk(4)
                iocc(2,k,l,i) = iwk(5)

C               basis, l, radial qn must match
                if (iwk(1) /= i .or. iwk(2) /= l .or. iwk(3) /= k) then
                  call errnl(lerr)
                  call info5(2,0,0,
     .              '*Error: Expected ib,l,n= %i %i %i '//
     .              'but encountered %3:1,1i',i,l,k,iwk,0)
                endif

              enddo
            enddo

          enddo
          if (lerr) then
            call rx('GWinput is not consistent with ctrl file')
          else
            call info0(2,0,0,' ... ordering ok')
          endif

          call info0(2,1,0,' ... Checking whether sites with common '//
     .      'species have equivalent lcut'//
     .      '%N   site  1st in species')
          lerr = .false.; lwarn = .false.
          do  i = 1, nbas
            lerri = .false.; lwarni = .false.
            is = s_site(i)%spec
            ibs = iclbsj(is,ips,-nbas,1)
            if (ibs == i) cycle
            lmaxa = s_spec(is)%lmxa
            if (lmaxa == -1) cycle

            call info5(2,0,-1,'%,6i%,9i',i,ibs,0,0,0)
            if (lcutmx(i) /= lcutmx(ibs)) then
C             call errnl(lwarni)
              lwarni = .true.
              call info5(2,0,0,
     .          '%4p warning: ib= %i expected lcut = %i but '//
     .          ' encountered %i',i,lcutmx(i),lcutmx(ibs),0,0)
            endif
            if (.not. (lerri .or. lwarni))
     .        call info0(2,0,0,'%2p ... match')
            lerr = lerr .or. lerri; lwarn = lwarn .or. lwarni
          enddo
          if (lwarn) then
            call info0(2,0,0,'%5p(warning)'//
     .        ' Sites within a species are not equivalent')
          else
            call info0(2,0,0,'%5pSites within a species are equivalent')
          endif

          call info0(2,1,0,' ... Checking whether sites in common '//
     .      'classes or species have equivalent occupation tables'//
     .      '%N   site  1st in class   1st in species')
          lerr = .false.; lwarn = .false.
          do  i = 1, nbas
            lerri = .false.; lwarni = .false.
            is = s_site(i)%spec
            ic = s_site(i)%class
            ibc = iclbsj(ic,ipc,-nbas,1)
            ibs = iclbsj(is,ips,-nbas,1)
            if (ibc == i .and. ibs == i) cycle
            lmaxa = s_spec(is)%lmxa
            if (lmaxa == -1) cycle
            call info5(2,0,-1,'%,6i%,9i%,15i',i,ibc,ibs,0,0)
            pnz = s_spec(is)%pz
            do  l = 0, lmaxa
              npqn = 2
              if (pnz(l+1,1) /= 0) npqn = 3
              do  k = 1, npqn

                if (iocc(1,k,l,ibc) /= iocc(1,k,l,i) .or.
     .              iocc(2,k,l,ibc) /= iocc(2,k,l,i)) then
                  call errnl(lerri)
                  call info8(2,0,0,
     .              '%3p*Error:   ib,l,n= %i %i %i '//
     .              'expected iocc,uocc = %i,%i '//
     .              'but encountered %i,%i',i,l,k,
     .              iocc(1,k,l,i),iocc(2,k,l,i),
     .              iocc(1,k,l,ibc),iocc(2,k,l,ibc),0)
                endif

                if (ibc /= ibs) then
                if (iocc(1,k,l,ibs) /= iocc(1,k,l,i) .or.
     .              iocc(2,k,l,ibs) /= iocc(2,k,l,i)) then
                  call errnl(lwarni)
                  call info8(2,0,0,
     .              '%4p warning: ib,l,n= %i %i %i '//
     .              'expected iocc,uocc = %i,%i '//
     .              'but encountered %i,%i',i,l,k,
     .              iocc(1,k,l,i),iocc(2,k,l,i),
     .              iocc(1,k,l,ibs),iocc(2,k,l,ibs),0)
                endif
                endif
              enddo
            enddo
            if (.not. (lerri .or. lwarni))
     .        call info0(2,0,0,'%2p ... match')
            lerr = lerr .or. lerri; lwarn = lwarn .or. lwarni
          enddo
          deallocate(iocc)

          if (lerr) then
            call rx(
     .        'Sites within a class are not equivalent')
          else
            call info0(2,0,0,'%5pSites within a class are equivalent')
          endif

          if (lwarn) then
            call info0(2,0,0,'%5p(warning)'//
     .        ' Sites within a species are not equivalent')
          else
            call info0(2,0,0,'%5pSites within a species are equivalent')
          endif

          call info0(2,1,-1,' ... Checking core product basis table')
          if (.not. rdstrn(ifii,strn,len(strn),.false.)) goto 999
          call word(strn,1,j1,j2)
          if (strn(j1:j2) /= 'atom') call rx(
     .      ' Line should begin with "atom" but reads '//strn(j1:j2))
          lerr = .false.
          allocate(iocc(4,10,0:nl-1,nbas)); iocc = NULLI
          do  i = 1, nbas
            is = s_site(i)%spec
            lmaxa = s_spec(is)%lmxa
            if (lmaxa == -1) cycle

            pnu = s_spec(is)%p
            pnz = s_spec(is)%pz
            call atwf(0,xx,lmaxa,xx,nsp,pnu,pnz,xx,xx,xx,xx,xx,k,ncore,
     .        konfig,xx,xx,xx)
C           read ib,l,radial pqn k, and iocc,uocc. 1st 3 must match ctrl
            do  l = 0, lmaxa
              do  k = l+1, mod(konfig(l),10)-1
C               Read table from disk; save occ list for further checks
                read(ifii,*) iwk(1:7)
                iocc(1,k,l,i) = iwk(4)
                iocc(2,k,l,i) = iwk(5)
                iocc(3,k,l,i) = iwk(6)
                iocc(4,k,l,i) = iwk(7)

C               basis, l, radial qn must match
                if (iwk(1) /= i .or. iwk(2) /= l .or. iwk(3) /= k-l)then
                  call errnl(lerr)
                  call info5(2,0,0,
     .              '*Error: Expected ib,l,n= %i %i %i '//
     .              'but encountered %3:1,1i',i,l,k-l,iwk,0)
                endif
              enddo
            enddo
          enddo

          if (lerr) then
            call rx('GWinput is not consistent with ctrl file')
          else
            call info0(2,0,0,' ... ordering ok')
          endif

          call info0(2,1,0,' ... Checking whether sites in common '//
     .      'classes or species have equivalent tables'//
     .      '%N   site  1st in class   1st in species')
          lerr = .false.; lwarn = .false.
          do  i = 1, nbas
            lerri = .false.; lwarni = .false.
            is = s_site(i)%spec
            ic = s_site(i)%class
            ibc = iclbsj(ic,ipc,-nbas,1)
            ibs = iclbsj(is,ips,-nbas,1)
            if (ibc == i .and. ibs == i) cycle
            lmaxa = s_spec(is)%lmxa
            if (lmaxa == -1) cycle
            call info5(2,0,-1,'%,6i%,9i%,15i',i,ibc,ibs,0,0)

            pnu = s_spec(is)%p
            pnz = s_spec(is)%pz
            call atwf(0,xx,lmaxa,xx,nsp,pnu,pnz,xx,xx,xx,xx,xx,k,ncore,
     .        konfig,xx,xx,xx)
C           read ib,l,radial pqn k, and iocc,uocc. 1st 3 must match ctrl
            do  l = 0, lmaxa
              do  k = l+1, mod(konfig(l),10)-1

                if (iocc(1,k,l,ibc) /= iocc(1,k,l,i) .or.
     .              iocc(2,k,l,ibc) /= iocc(2,k,l,i) .or.
     .              iocc(3,k,l,ibc) /= iocc(3,k,l,i) .or.
     .              iocc(4,k,l,ibc) /= iocc(4,k,l,i)) then
                  call errnl(lerri)
                  call info5(2,0,0,
     .              '%3p*Error:   ib,l,n= %i %i %i '//
     .              'expected iocc,uocc,ForX0,ForSxc = %4:1i '//
     .              'but encountered %4:1i',i,l,k-1,
     .              iocc(1,k,l,i),iocc(1,k,l,ibc))
                endif

                if (ibc /= ibs) then
                if (iocc(1,k,l,ibs) /= iocc(1,k,l,i) .or.
     .              iocc(2,k,l,ibs) /= iocc(2,k,l,i) .or.
     .              iocc(3,k,l,ibs) /= iocc(3,k,l,i) .or.
     .              iocc(4,k,l,ibs) /= iocc(4,k,l,i)) then
                  call errnl(lerri)
                  call info5(2,0,0,
     .              '%4pwarning:   ib,l,n= %i %i %i '//
     .              'expected iocc,uocc,ForX0,ForSxc = %4:1i '//
     .              'but encountered %4:1i',i,l,k-1,
     .              iocc(1,k,l,i),iocc(1,k,l,ibs))
                endif
                endif
              enddo
            enddo
            if (.not. (lerri .or. lwarni))
     .        call info0(2,0,0,'%2p ... match')
            lerr = lerr .or. lerri; lwarn = lwarn .or. lwarni
          enddo
          deallocate(iocc)

          if (lerr) then
            call rx('Sites within a class are not equivalent')
          else
            call info0(2,0,0,'%5pSites within a class are equivalent')
          endif

          if (lwarn) then
            call info0(2,0,0,'%5p(warning)'//
     .        ' Sites within a species are not equivalent')
          else
            call info0(2,0,0,'%5pSites within a species are equivalent')
          endif

          print *

          if (QpGcut_cou == NULLI) call rx('QpGcut_cou was not read')
          if (QpGcut_psi == NULLI) call rx('QpGcut_psi was not read')
          if (QpGcut_psi < QpGcut_cou)
     .      call rx('Error : Increase QpGcut_psi to be QpGcut_psi > QpGcut_cou')

          call rx0('done checking GWinput')

        case (1)   ! QpGcut_cou
          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' QpGcut_cou has invalid argument')
          else
            vol0 = vol/nbas0
            qsafe = 180d0/vol0
            if (xx+.1d0 < qsafe) then
              call info2(2,0,0,'*QpGcut_cou=%,1d '//
     .        'may be too low: %;1d is safer',xx,qsafe)
            else
              call info2(2,0,0,' QpGcut_cou=%d ... ok',xx,0)
            endif
            QpGcut_cou = xx
          endif
        case (2)   ! QpGcut_psi
          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' QpGcut_psi has invalid argument')
          else
            vol0 = vol/nbas0
            qsafe = 220d0/vol0
            if (xx+.1d0 < qsafe) then
              call info2(2,0,0,'*QpGcut_psi=%,1d '//
     .        'may be too low: %;1d is safer',xx,qsafe)
            else
              call info2(2,0,0,' QpGcut_psi=%d ... ok',xx,0)
            endif
          endif
          QpGcut_psi = xx
        case (3)   ! Verbose
          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' Verbose has invalid argument')
          else
            if (xx >= 40) then
              call info0(2,0,0,' Verbose ... ok')
            else
              call info0(2,0,0,' Use verbose>=40 to get timings ')
            endif
          endif

        case (4)   ! dw
          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' dw has invalid argument')
          else
            if (xx <= .01000001d0) then
              call info2(2,0,0,' dw=%d ... ok for QSGW '//
     .          '(but should be small to calc. susceptibility)',
     .          xx,0)
            elseif (xx <= .02000001d0) then
              call info2(2,0,0,' dw=%,2d '//
     .        'is a bit high: .01 is safe for QSGW',xx,0)
            else
              call info2(2,0,0,'*dw=%,2d '//
     .        'may be too high: .01 is safe for QSGW',xx,0)
            endif
          endif

        case (5)   ! esmr

          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' esmr has invalid argument')
          else
            call info2(2,0,0,' esmr=%,2g. '//
     .        '0.003 is typical for insulators.  '//
     .        'For some metals need esmr ~0.01',
     .        xx,0)
          endif

        case (6)   ! iSigMode

          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' iSigMode has invalid argument')
          else
            if (xx == 3) then
              call info0(2,0,0,' iSigMode = 3 ... ok for QSGW ')
            else
              call info2(2,0,0,'*iSigMode=%d. '//
     .        'Usually use iSigMode=3 for QSGW.',xx,0)
            endif
          endif

        case (7)   ! n1n2n3

          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            j = a2vec(strn(j1:),len(strn(j1:)),j,2,', ',2,-3,3,ix,iwk)
            if (j /= 3) then
              call info0(2,0,0,' n1n2n3 has invalid argument')
            else
              call info2(2,0,0,'*n1n2n3=%3:,2i ... '//
     .          'monitor convergence in n1n2n3',iwk,0)
            endif
          endif

        case (8)   ! omg_c

          lerr = .true.
          call word(strn,2,j1,j2)
          if (j1 <= j2) then
            j = 0
            lerr = .not. a2bin(strn(j1:j2),xx,4,0,' ',j,-1)
          endif
          if (lerr) then
            call info0(2,0,0,' omg_c has invalid argument')
          else
            if (xx >= .04d0) then
              call info2(2,0,0,' omg_c=%,2d should be ok for QSGW',xx,0)
            else
              call info2(2,0,0,' omg_c=%,2d : '//
     .        'check convergence in this parameter',xx,0)
            endif
          endif

      end select
      enddo


      call rx0('completed check')
  999 continue
      call rx('missing part of GWinput')

      end
      subroutine errnl(lerr)
C- Prints a newline if lerr is not set; then sets lerr to T.
      implicit none
      logical lerr

      if (.not. lerr) print *
      lerr = .true.
      end
