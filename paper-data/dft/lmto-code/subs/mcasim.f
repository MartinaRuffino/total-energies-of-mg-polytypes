      subroutine sqssim(s_move,nbas,ips,iprm,sigb,ntab,iax,e2,ntab3,iax3,e3,esqs)
C- Generate SQS via simulated annealing
C ----------------------------------------------------------------------
Ci Inputs
Ci     Elts read:  kt ts tstot tsequ
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci   nbas  :size of basis
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   sigb  :1 for species of first kind, -1 for species of 2nd kind
Ci   ntab  :ntab(ib)=offset to pair neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   e2    :configuration energy for each pair in pair table
Ci   ntab3 :ntab(ib)=offset to triples neighbor table for cluster ib (pair3c.f)
Cl   iax3  :triplets table for SQS.  Refers to pair table iax
Cl         :iax3(1,i3c) = kb
Cl         :iax3(2,i3c) = ib
Cl         :iax3(3,i3c) = jb
Cl         :iax3(4,i3c) = ip : points to pair index iaxs(:,ip) connecting ib and kb
Cl         :iax3(5,i3c) = jp : points to pair index iaxs(:,jp) connecting jb and kb
Cl         :iax3(6,i3c) = kp : points to pair index iaxs(:,kp) connecting ib and jb
Ci   e3    :configuration energy for each triplet in triplets table
Cio Inputs/Outputs
Cio  iprm  :site permutation table: {ib} are permuted to {iprm(ib)}
Co Outputs
Co   esqs  :configuration energy
Co         :esqs(0) = configuration energy e2**2 + e3**2
Co         :esqs(1) = not used
Co         :esqs(2) = pair energy e2
Co         :esqs(3) = triplet energy e3
Cl Local variables
Cl   iprmi :inverse of iprm
Cr Remarks
Cr Two-center Hamiltonian is  1/2 sum_ib,jb e_ib,jb sigma_ib sigma_jb
Cu Updates
Cu   04 Sep 16 sign of s_move%tstot : tstot<0 => flag to end on minimum energy
Cu   28 Jul 15 Redesigned
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer, parameter ::  niax=10
      integer nbas,ips(nbas),iprm(nbas),sigb(nbas)
      integer ntab(nbas+1),iax(niax,*)
      integer ntab3(nbas+1),iax3(niax,*)
C     double precision e2(nttab),e3(nttab3)
      double precision e2(*),e3(*)
C ... For structures
!      include 'structures.h'
      type(str_move)::  s_move
C ... Dynamically allocated local arrays
      integer, allocatable :: iprmi(:)
C ... Local parameters
      logical lundo,metropolis,endonmin
      integer ntry,i1,i2,jb,nclus,ip,ipr,nswap,stdo,nttab,nttab3,nochange
      integer ntot,nequ,i,ib
      double precision kt,ts,ttot,tequ
      double precision e2c,edc,esqs(0:3),eold(0:3),de,ebar,esqr,etol,emin,estart
      character*80 outs
      integer, parameter :: mxnochange = 100 000
      procedure(logical) :: cmdopt
      procedure(integer) :: nglob,lgunit,iprint
      procedure(real(8)) :: ddot
      procedure(real) :: ran1

      call tcn('sqssim')

C --- Setup and initial printout ---
      stdo = nglob('stdo')
      call getpr(ipr)
      allocate(iprmi(nbas))
      nttab = ntab(nbas+1)
      nttab3 = ntab3(nbas+1)
      nochange = 0

C ... Initial printout
      kt = s_move%kt            ! temperature
      ts = s_move%ts            ! time step
      ttot = dabs(s_move%tstot)
      tequ = s_move%tsequ
      ntot = ttot/ts
      nequ = tequ/ts
      endonmin = s_move%tstot < 0

      call awrit5(' mcasim:   %i mc steps (ttot=%;4d, ts=%;4d)'//
     .  ' incl. %i equil. steps,  kT=%;4d',
     .  outs,80,-lgunit(2),ntot,ttot,ts,nequ,kt)
      if (iprint() >= 20) then
        print *
        call awrit0('%a',outs,80,-lgunit(1))
      endif
C     Create inverse of iprm and initial sigma
      do  i = 1, nbas
        ib = iprm(i)
        iprmi(ib) = i
      enddo
      etol = (ddot(nttab,e2,1,e2,1) + ddot(nttab3,e3,1,e3,1))*1d-12
      emin = 9d9

      if (ran1() == 0) call rx('sqssim: random number generator not seeded')

C     Get starting esqs
      call mcasan(1,nbas,iprm,iprmi)  ! Sanity checks
      call mcasho(1,nbas,iprm,iprmi,sigb,ntab,iax,e2,ntab3,iax3,e3,eold)
C     call mcasan(1,nbas,iprm,iprmi,ntab,iax,ips,esite,e2,sigb,sigl,sigma,eold)
C     call mcasho(1,nbas,iprm,iprmi,ntab,iax,esite,e2,sigma,sigb,eold)
      esqs = eold
      estart = eold(0)
      ebar = 0
      esqr = 0
      nswap = 0

      if (cmdopt('--shopr',7,0,outs)) then
        call mcasme2s(1,nbas,iprm,iprmi,sigb,ntab,iax,e2,esqs(2),edc)
        call rx0('done')
      endif

      if (ntot == 0) call rx0('no iterations')


C --- For ntot time steps, do ---
      ntry = 0
!     do  while (nswap < ntot .or. (endonmin .and. abs(emin-esqs(0)) > etol) .and. nochange < mxnochange)
      do  while (nswap < ntot .or. (endonmin .and. abs(emin-esqs(0)) > etol))
        ntry = ntry+1

C   ... Decide on first of pair of atoms to swap
   12   continue
        ib = nbas*ran1() + 1
        ib = min(ib,nbas)
C   ... If no sigma, no contribution to hamiltonian ... skip
C       if (sigb(ib) == 0) goto 12
        if (ib > nbas) stop 'oops'
        i1 = iprmi(ib)

C   ... Decide on second of pair of atoms to swap
C       Require that it belong to cluster i1
        nclus = ntab(i1+1)-ntab(i1)
        if (nclus == 0) goto 12
        ip = nclus*ran1() + 1 + ntab(i1)
        i2 = iax(2,ip)
        jb = iprm(i2)
C   ... If no sigma, no contribution to hamiltonian ... skip
C        if (sigb(jb) == 0) goto 12
C   ... If equivalent sigma, no change in  hamiltonian ... skip
        if (sigb(jb) == sigb(ib)) goto 12
C   ... special constraints on movement
        if (ips(ib) > 2 .or. ips(jb) > 2) goto 12

C       Faster to compute change in 2-body part rather than recompute all of it
C       Subtract pair energy connected with the two sites before swap
        call mcasme2(ib,ib,iprm,iprmi,sigb,ntab,iax,e2,e2c,edc)
        esqs(2) = esqs(2) - (e2c-edc)*2 - edc
        call mcasme2(jb,jb,iprm,iprmi,sigb,ntab,iax,e2,e2c,edc)
        esqs(2) = esqs(2) - (e2c-edc)*2 - edc

C       Swap the atoms
        call mcasw(ib,jb,iprm,iprmi)

C       Add pair energy connected with the two sites after swap
        call mcasme2(ib,ib,iprm,iprmi,sigb,ntab,iax,e2,e2c,edc)
        esqs(2) = esqs(2) + (e2c-edc)*2 + edc
        call mcasme2(jb,jb,iprm,iprmi,sigb,ntab,iax,e2,e2c,edc)
        esqs(2) = esqs(2) + (e2c-edc)*2 + edc

C       No corresponding trick for three-body energy ... sum over all triplets
        call mcasme3(1,nbas,iprm,iprmi,sigb,ntab3,iax3,e3,esqs(3),edc)

        de = esqs(2)**2 - eold(2)**2 + esqs(3)**2 - eold(3)**2
        esqs(0) = esqs(0) + de

C       This should generate the same answer a bit slower)
C#ifdefC DEBUG
C        edc = esqs(0)
C        call pshpr(1)
C        call mcasho(1,nbas,iprm,iprmi,sigb,ntab,iax,e2,ntab3,iax3,e3,esqs)
C        call poppr(1)
C        if (abs(edc-esqs(0)) > 1d-12) stop 'oops'
C#endif

C  ...  Metropolis step : retain swap with probability exp(-de/kt)
        lundo = .false.; metropolis = .false.
        if (de > 0) then
          if (ran1() > dexp(-de/kt)) then
C           Undo atom swap
            call mcasw(jb,ib,iprm,iprmi)
            esqs = eold
            lundo = .true.
          else
            metropolis = .true.
          endif
        endif
        if (.not. lundo) then
          nswap = nswap+1
          nochange = 0
        else
          nochange = nochange + 1
        endif

C   ... Keep statistics
        if (.not. lundo .and. nswap > nequ) then
          emin = min(emin,esqs(0))
          ebar = ebar + esqs(0)
          esqr = esqr + esqs(0)**2
        endif

        if (ipr <= 50 .and. ipr > 30) then
          if (mod(ntry,1000) == 0)
     .      call awrit4(' sqssim completed %i+%i swaps in %i attempts'
     .      //' e=%d',' ',80,stdo,nswap-nequ,nequ,ntry,eold)
        endif
        if (ipr > 50 .or. ipr >= 40 .and. .not. lundo)  then
C          print 357, ib,jb,esqs(0),de,lundo,metropolis
C  357     format(' swap',i4,i4,' enow = ',f12.6,' de = ',f12.6,'  undo=',L1,'  metrop=',L1)
          print 357, ib,jb,esqs(0),de,metropolis
  357     format(' swap',i4,i4,' enow = ',f12.6,' de = ',f12.6,'  metrop=',L1)
        endif

        eold = esqs

C        print *, 'lundo=',lundo
C        call mcasho(1,nbas,iprm,iprmi,sigb,ntab,iax,e2,ntab3,iax3,e3,esqs)

        if (esqs(0) < etol) exit  ! Configuration energy is zero ... quit

        if (nochange >= mxnochange) then
          call info2(20,0,0,'no swaps in %i iterations ... quitting',nochange,2)
          exit
        endif
      enddo ! nswap

      ebar = ebar/(nswap-nequ)
      esqr = dsqrt(esqr/(nswap-nequ) - ebar**2)

      if (ipr >= 10) then
        do  ip = 1, 2
          call awrit3('%N sqssim completed %i+%i swaps in %i attempts:',
     .      ' ',120,lgunit(ip),nswap-nequ,nequ,ntry)
          call awrit6(' kt = %;4d  estart = %;4d  emin = %;4d  enow = %;4d  ebar = %;4d  rms(e-ebar) = %;4d ',
     .      ' ',120,lgunit(ip),kt,estart,emin,esqs(0),ebar,esqr)
        enddo
      endif

      call mcasan(1,nbas,iprm,iprmi)
      call mcasho(1,nbas,iprm,iprmi,sigb,ntab,iax,e2,ntab3,iax3,e3,esqs)
      call tcx('sqssim')
      end
      subroutine mcasw(ib1,ib2,iprm,iprmi)
C- Swap two entries in site permutation table
C ----------------------------------------------------------------------
Ci Inputs
Ci  ib1,ib2:pairs of sites to swap
Ci   iprm  :site permutation table: atom ib=iprm(i) is moved to site i
Ci   iprmi :inverse of iprm: atom ib is moved to site iprmi(ib)
Co Outputs
Co   iprm  :entries iprm(i1) and iprm(i2) are swapped, with
Co          ib1 = iprm(i1), ib2=iprm(ib2)
Co   iprmi :entries iprmi(ib1) and iprmi(ib2) are swapped
Co   sigma :sigma(ib1) = sigb(ib1)*sigl(i1) and
Co         :sigma(ib2) = sigb(ib2)*sigl(i2) are updated.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ib1,ib2,iprm(*),iprmi(*)
C ... Local parameters
      integer i1,i2

      i1 = iprmi(ib1)
      i2 = iprmi(ib2)
      iprmi(ib1) = i2
      iprmi(ib2) = i1
C#ifdefC DEBUG
C      call awrit4('  ... swap i1(ib1)=%i(%i) with i2(ib2)=%i(%i)',
C     .  ' ',80,6,i1,ib1,i2,ib2)
C      if (iprm(i1) /= ib1) stop 'bug in mcasw'
C      if (iprm(i2) /= ib2) stop 'bug in mcasw'
C#endif
      iprm(i1) = ib2
      iprm(i2) = ib1

CC     sigma(ib1) = sigb(ib1)*sigl(iprmi(ib1))
C      sigma(ib1) = sigb(ib1)*sigl(i2)
CC     sigma(ib2) = sigb(ib2)*sigl(iprmi(ib2))
C      sigma(ib2) = sigb(ib2)*sigl(i1)

      end

      subroutine mcasms(mode,i1,i2,ips,sqspid,sigma)
C- Generate sigma
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :determines how sigma is made
Ci  i1,i2  :range of sites for which sigma is made
Ci   ips   :species table: site i belongs to species ips(i)
Ci  sqspid :sqspid(1) = ID for first alloy species
Ci         :sqspid(2) = ID for 2nd   alloy species
Co Outputs
Co   sigma :sigma is generated for lattice sites i=i1..i2
Cr Remarks
Cr   mode=0:  sigma(i) = 1 if ips(i) = 1, sigma = -1 if ips(i) = 2
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,i1,i2,ips(i2),sqspid(2),sigma(i2)
C ... Local parameters
      integer i

      if (mode /= 0) call rxi('mcasms: mode not recognized, ',mode)
      do  i = i1, i2
        if (ips(i) == sqspid(1)) then
          sigma(i) = 1
        elseif (ips(i) == sqspid(2)) then
          sigma(i) = -1
        else
          sigma(i) = 0
        endif
      enddo

      end
      subroutine mcasme2(ib1,ib2,iprm,iprmi,sig2,ntab,iax,e2,e2c,edc)
C- Pair configuration energy for pairs associated with a range of sites
      implicit none
C ... Passed parameters
      integer, parameter ::  niax=10
      integer ib1,ib2
      integer iax(niax,*),ntab(*),iprm(*),iprmi(*),sig2(ib2)
      double precision e2(*),e2c,edc
C ... Local parameters
      integer ib,ibp,ip,jbp,jb
      logical lonsite


      call tcn('mcasme2')

C --- Pair configuration energy ---
      e2c = 0
      edc = 0
C     For each site do
      do  ib = ib1, ib2
        ibp = iprmi(ib)
C       Loop over each pair associated with site ib
        do  ip = ntab(ibp)+1, ntab(ibp+1)
          jbp = iax(2,ip)
          jb = iprm(jbp)
          e2c = e2c + sig2(ib)*sig2(jb)*(e2(ip)/2)
          lonsite = ib == jb .and. iax(3,ip) == 0 .and. iax(4,ip) == 0 .and. iax(5,ip) == 0
          if (lonsite) edc = edc + sig2(ib)*sig2(jb)*(e2(ip)/2)
C         if (lonsite) stop 'hi'
        enddo
      enddo

      call tcx('mcasme2')

      end

      subroutine mcasme2s(ib1,ib2,iprm,iprmi,sig2,ntab,iax,e2,e2c,edc)
C- Resolve pair configuration energy for display
      implicit none
C ... Passed parameters
      integer, parameter ::  niax=10
      integer ib1,ib2
      integer iax(niax,*),ntab(*),iprm(*),iprmi(*),sig2(ib2)
      double precision e2(*),e2c,edc
C ... Dynamically allocated local arrays
      integer, allocatable :: ssigpr(:)
C ... Local parameters
      integer ib,ibp,ip,jbp,jb,npr
      logical lonsite
      real(8), parameter ::  tol = 1d-6
      character*256 strn

      ip = maxval(ntab(1:ib2))
      allocate(ssigpr(0:ip))

      call info0(10,1,0,' site  sig')

C --- Pair configuration energy ---
      e2c = 0
      edc = 0
C     For each site do
      do  ib = ib1, ib2
        ibp = iprmi(ib)
C       Loop over each pair associated with site ib
        npr = 0; ssigpr(0) = 0; strn = ' '
        do  ip = ntab(ibp)+1, ntab(ibp+1)
          jbp = iax(2,ip)
          jb = iprm(jbp)
          e2c = e2c + sig2(ib)*sig2(jb)*(e2(ip)/2)
          lonsite = ib == jb .and. iax(3,ip) == 0 .and. iax(4,ip) == 0 .and. iax(5,ip) == 0
          if (lonsite) edc = edc + sig2(ib)*sig2(jb)*(e2(ip)/2)
          if (.not. lonsite) then
            if (ip > ntab(ibp)+1) then
              if (abs(e2(ip)-e2(ip-1)) > tol) write(strn(4*npr:4*npr),"('|')")
            endif
            npr = npr+1
            ssigpr(npr) = ssigpr(npr-1) + sig2(ib)*sig2(jb)
            write(strn(-3+4*npr:-1+4*npr),'(i3)') sig2(ib)*sig2(jb)
          endif
        enddo
        if (npr == 0) cycle
        write(*,"(i4,2x,a)") ib, trim(strn)
        call info2(10,0,0,'     %n,4i',npr,ssigpr(1))
      enddo

      end


      subroutine mcasme3(kb1,kb2,iprm,iprmi,sig3,ntab3,iax3,e3,e3c,edc)
C- Three-body configuration energy for triplets associated with a range of sites
      implicit none
C ... Passed parameters
      integer, parameter :: niax=10
      integer kb1,kb2,iprm(*),iprmi(*),sig3(*)
      integer iax3(niax,*),ntab3(*)
      double precision e3(*),e3c,edc
C ... Local parameters
      integer ib,jb,kb,ibp,jbp,kbp,ip
      double precision e3ci

      call tcn('mcasme3')

C --- Three-body configuration energy ---
      e3c = 0
      edc = 0  ! Should be zero
C     For each site do
      do  kb = kb1, kb2
        kbp = iprmi(kb)
C       Loop over each pair in the corresponding pair table
        e3ci = e3c
        do  ip = ntab3(kbp)+1, ntab3(kbp+1)
          ibp = iax3(2,ip)
          ib = iprm(ibp)
          jbp = iax3(3,ip)
          jb = iprm(jbp)
          e3c = e3c + sig3(ib)*sig3(jb)*sig3(kb)*(e3(ip)/3)
C         Check for possible triplets occupying only two sites. This should never happen
          if (ib == jb .or. ib == kb .or. jb == kb) then
            if (ib == jb .and. iax3(4,ip) == iax3(5,ip)) then
              call rx('oops! problem with dc in mcasme3')
            endif
C           edc = edc + sig3(ib)*sig3(jb)*(e3(ip)/3)
          endif
        enddo
C       print *, kb,ip,e3c-e3ci
      enddo

C     if (edc /= 0) call rx('oops! problem with dc in mcasme3')

      call tcx('mcasme3')

      end

      integer function iosiga(lread,fnam,nbas,ips,sigl,sigb,iprm,pos)
Cr Remarks
Cr   Loaded into variables table for each of the ib=1..nbas lines read:
Cr     is = ips(ib)
      implicit none
      integer lread,nbas,sigl(nbas),sigb(nbas),iprm(nbas),ips(nbas)
      character*(*) fnam
      double precision pos(3,nbas)
      integer ifi,ib,fopna,j1,j2,is,ival
      double precision px,py,pz
      character*1 ch
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),
     .  vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev),a2bin
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)

      iosiga = 0
      ifi = fopna(fnam,-1,0)
      rewind ifi
      if (lread == 0) then
        read(ifi,'(a)',end=99,err=99) ch
        nr = 0
        do  ib = 1, nbas
          call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .                ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
          if (nr < 0) goto 99
C         read(a,*,end=99,err=99) sigl(ib),sigb(ib),iprm(ib),px,py,pz
          is = ips(ib)
          call lodsyv('is',1,dble(is),ival)
          j1 = 0
          if  (.not. a2bin(a,sigl(ib),2,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,sigb(ib),2,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,iprm(ib),2,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,px,4,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,py,4,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,pz,4,0,' ',j1,-1)) goto 99
          if (abs(px-pos(1,ib))+abs(py-pos(2,ib))+abs(pz-pos(3,ib))
     . > 1d-6) goto 99
        enddo
        print *, 'iosiga:  read sigl,sigb,iprm from file ',fnam
      else
        write(ifi,
     .    '(''#   sigl sigb  iprm              pos'',21x,''spec'')')
        do  ib = 1, nbas
          write(ifi,333) sigl(ib),sigb(ib),iprm(ib),
     .      pos(1,ib),pos(2,ib),pos(3,ib),ips(iprm(ib))
  333     format(3i6,3f12.7,i6)
        enddo
        call word(fnam,1,j1,j2)
        print *, 'iosiga:  wrote sigl,sigb,iprm to file ',fnam(j1:j2)
      endif
      call fclose(ifi)

      return
C --- Error exit ---
   99 continue
      print *, 'iosiga: file read failed or file mismatch'
      iosiga = -1
      call fclose(ifi)
      end

      subroutine mcasho(ib1,ib2,iprm,iprmi,sigb,ntab,iax,e2,ntab3,iax3,e3,esqs)
C- Printout
C ----------------------------------------------------------------------
Ci Inputs
Ci  ib1,ib2:range of sites for which sig1 is calculated
Ci   iprm  :site permutation table: atom ib=iprm(i) is moved to site i
Ci   iprmi :inverse of iprm: atom ib is moved to site iprmi(ib)
Ci   sigb  :sigma for 2-center part of hamiltonian, sigb_i sigb_j e2_ij
Co Outputs
Co   esqs  :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ib1,ib2
      integer, parameter ::  niax=10
      integer sigb(ib2),iprm(ib2),iprmi(ib2)
      double precision e2(*),e3(*),esqs(0:3)
      integer ntab(ib2),iax(niax,*),ntab3(ib2+1),iax3(niax,*)
C ... Local parameters
      integer i,ncol
      double precision edc
      procedure(integer) :: iprint

      esqs = 0
      call mcasme2(ib1,ib2,iprm,iprmi,sigb,ntab,iax,e2,esqs(2),edc)
      call mcasme3(ib1,ib2,iprm,iprmi,sigb,ntab3,iax3,e3,esqs(3),edc)
      esqs(0) = esqs(2)**2 + esqs(3)**2
      if (iprint() < 30) return
      call info5(10,1,0,'  mcasho :  sites %i:%i  e2=%;4d  e3=%;4d  esqs=%;4d',
     .  ib1,ib2,esqs(2),esqs(3),esqs)
      ncol = min(4,ib2-ib1+1)
      call arrprt('  ib site sig','%,4i%,4i%,4i','Iii',
     .  ib2-ib1+1,0,ncol,0,'  | ',i,iprmi,sigb,i,i,i,i,i)

      end
      subroutine mcasan(ib1,ib2,iprm,iprmi)
C- Sanity checking
      implicit none
C ... Passed parameters
      integer ib1,ib2,iprm(ib2),iprmi(ib2)
C      parameter (niax=10)
C      integer sigb(ib2),sigl(ib2),sig1(ib2),ips(ib2),iax(niax,*)
C      double precision esite(ib2),epair(*),esqs
C ... Local parameters
      integer i,ib

C     Check that iprm and iprmi are consistent
      do  i = ib1, ib2
        ib = iprm(i)
        if (iprmi(ib) /= i) call rx('bug in mcasim')
C        if (sig1(ib) /= sigb(ib)*sigl(i)) then
C          call rx('bug in mcasim')
C        endif
      enddo

      end
