      subroutine aiocls(lio,mode,s_ctrl,s_ham,s_pot,s_spec,s_lat,ic1,ic2)
C- File I/O atomic data for classes ic1..ic2
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nl nspin nclass nccomp nclasp ics nrc zbak nbas
Ci                 initc dclabl
Co     Stored:     initc zbak
Co     Allocated:  *
Cio    Elts passed:ics dclabl initc lasa
Cio    Passed to:  *
Cio  s_ham  : (not used)
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  rhrmx
Co     Stored:     rhrmx
Co     Allocated:  *
Cio    Elts passed: pnu qnu vrmax pp pprel ves vintr pmpol sop grrme
Cio                rhrmx
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa idmod p q
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
Ci Inputs
Ci   lio    F for read, T for write
Ci   mode   ones digit
Ci          1 use data from first class corresponding to same species,
Ci            if data from own class is missing
Ci          2 like 1, but attempt to read data from disk anyway
Ci          3 make no attempt to read from disk but copy info from
Ci            first class corresponding to same species, if data missing
Ci          4 Add 4 if to use default P,Q when not otherwise supplied
Ci          10s digit
Ci          1 assemble background rho
Ci          100s digit
Ci          1 also attempt to read potential (ccomp file with READPOT)
Ci   ic1,ic2: range of classes to read data
Cr Remarks
Cr   Right now, read always takes data from file if available
Cu Updates
Cu   07 Aug 16 If Dirac qnur is available on disk, retain it even if lrel=1
Cu   01 Aug 16 Some changes for redsign of Dirac solver
Cu   10 May 15 Copies scalar relativistic enu to qnur if fully rel. is absent
Cu   22 Feb 15 Extension to fully relativistic GF
Cu   14 Nov 13 Some adjustments in preparation for fully relativistic GF
Cu   17 Jun 13 Completed f90 replacement of f77 pointers
Cu   08 May 13 Eliminate s_array
Cu   08 Sep 11 Begin migration to f90 structures
Cu   28 May 08 (Kirill) extensions for disordered local moments
Cu   09 Nov 07 Corrected sign of default moment (paioc2)
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   26 Apr 03 Added MPI calls
Cu   07 Feb 03 adjusted for redimensioned sop
Cu   30 May 02 Assign better default P
Cu   28 Apr 98 I/O of radial matrix elements of grad
Cu   28 Sep 00 Added setting default P,Q
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lio
      integer mode,ic1,ic2
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      real(8), allocatable :: lpp(:),lppr(:),lsop(:),lgrme(:),lva(:),lmp(:),pot(:),cor(:)
C ... Local parameters
      character*8 clabl,alabel, outs1*20, outs2*20, outs3*20, outs*80
      logical sw,aiomom,aiopar,aiopot,aiova,lpot,lcor,scat,haverel,
     .  aiogen,aiosop,aiorme,aiocor,aiomp,lrell,lgen
      logical havepq,havepp,haveso,haveop,haveva,havemp,havev,
     .        readpq,readpp,readso,readop,readva,readmp,readv
      integer,parameter :: n0=10, NULLI=-99999, nrmx=5001
      integer ic,is,lmx,k,kr,nl,nsp,ifi,jfi,fopn,lmxx,nspx,nrx,nr,
     .  idmod(n0),isw,bitand,i2,nclasp,iclbsj,icmap,jc,nbas,nclass,
     .  nclspp,iprint,lgunit,mode0,mode00,lrel,nangl,nclspd,nclsppd
      integer mpipid,procid
      double precision rhrmx,vrmax(2),bxc(3),bhat(3),ves,z,rmxx,ax,qc,dq,
     .  vrmxx(2),sumec,sumtc,sumev,thrpv,ekin,utot,rhoeps,etot,a,rmax,
     .  zbak(2),pdf(n0,2),qdf(n0,2),pnuloc(100),qnuloc(100),qnurloc(800),qlm(800)
      integer nr2,nsp2
      double precision a2,rmax2
      procedure(integer) asars2,nglob
C#ifdefC READPOT
C      double precision v0(nrmx*2)
C#endif

      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      nclass = s_ctrl%nclass
      nangl = s_ctrl%nccomp
      nclasp = s_ctrl%nclasp
      nclspd = nclasp + nangl
!     nlspcr = 8*nl*nl*max(nclsppd,nspec)
      lrel = mod(nglob('lrel'),10)
      haverel = lrel == 2     ! If all fully relativistic q's are read or copied, haverel remains true

C     call prmx('enter aiocls bxc',s_pot%bxc,3,3,2)

      mode0 = mod(mode,10)
      mode00 = mod(mode0,4)

C     call prmx('aiocls1 qnu',s_pot%qnu,3*nl*nsp,3*nl*nsp,nclass)

C#ifdefC READPOT
CC     Better to initially allocate with size one, reallocate as needed (see asars)
CC     Allocate global array for potential, if if not yet allocated
C      if (mod(mode/100,10) > 0) then
C        ic = 0
C        if (.not. associated(s_pot%v0)) then
C          ic = asars2(2,s_ctrl,s_spec,s_pot,1) ! ic = max value of nr if all classes available from disk
C         if (ic > 0) call ptr_pot(s_pot,8+1,'v0',ic*nsp,s_ctrl%nclasp+s_ctrl%nccomp,z)
C        endif
C      endif
C#endif

C ... MPI: only master does sphere program
      procid = mpipid(1)
      if (procid == 0) then

      allocate(lpp(100),lppr(3000),lsop(500),lgrme(500),lva(100),lmp(3000),pot(3000),cor(3000))
      i2 = ic2
      if (i2 == 0) i2 = nclspd
c     call awrit2('%n:1i',' ',100,6,nclasp,s_ctrl%ics)
      do  ic = ic1, i2

        is = s_ctrl%ics(ic)
        icmap = iclbsj(is,s_ctrl%ics,-nclasp,1)
        if (icmap == ic .or. mode00 == 0) icmap = 0
        lmx = s_spec(is)%lmxa
        z = s_spec(is)%z
        idmod = s_spec(is)%idmod
        pdf = s_spec(is)%p
        qdf = s_spec(is)%q
        if (s_ctrl%nrc(ic) == 0) call dpzero(qdf,2*n0)

        if (lmx == -1) then
          s_ctrl%initc(ic) = 0
          cycle
        endif

        call dpzero(pnuloc,100)
        call dpzero(qnuloc,100)
        call dpzero(qnurloc,800); qnurloc(1) = NULLI+1  ! => not yet read
!       if (lrel /= 2) qnurloc(1) = NULLI ! => not sought

C   --- Open the atom file ---
        outs1 = ' '; outs2 = ' '; outs3 = ' '
        call r8tos8(s_ctrl%dclabl(ic),clabl)
        zbak = s_ctrl%zbak; nbas = s_ctrl%nbas; nclass = s_ctrl%nclass
        if (mode00 /= 3) ifi = fopn(clabl)

C  --- Copy what is passed through to holding arrays ---
        k = nl*nsp
        kr = 8*nl*nl
        call pvaioc(s_ctrl%initc,1,ic,icmap,havepq,jc)
        readpq = .not. havepq .or. mode00 == 2 .and. jc /= ic
        readpq = readpq .and. mode00 < 3
        if (havepq) then
          call dpscop(s_pot%pnu,pnuloc,k,1+(jc-1)*k,1,1d0)
          call dpscop(s_pot%qnu,qnuloc,3*k,1+(jc-1)*3*k,1,1d0)
          if (s_ctrl%stpnrho) call dpscop(s_pot%qnur,qnurloc,4*kr,1+(jc-1)*4*kr,1,1d0)
          rhrmx = s_pot%rhrmx(jc)
C         call dpscop(s_pot%rhrmx,rhrmx,1,jc,1,1d0)
          call dpscop(s_pot%vrmax,vrmax,2,2*jc-1,1,1d0)
          if (ic /= jc) call awrit0('%a pq,',outs1,len(outs1),0)
C         call dmscop(bxc,3,s_pot%bxc,3,1,3,jc,jc,1,1,1d0)
        endif
        k = 6*nl*nsp
        call pvaioc(s_ctrl%initc,2,ic,icmap,havepp,jc)
        readpp = .not. havepp .or. mode00 == 2 .and. jc /= ic
        readpp = readpp .and. mode00 < 3
        if (havepp) then
          call dpscop(s_pot%pp,lpp,k,1+(jc-1)*k,1,1d0)
          if (lrel == 2) then
            k = 5*nl*2*nl*2*2
            call dpscop(s_pot%pprel,lppr,k,1+(jc-1)*k,1,1d0)
          endif
          call dpscop(s_pot%ves,ves,1,jc,1,1d0)
          if (ic /= jc) call awrit0('%a pp,',outs1,len(outs1),0)
        endif

        k = (nl*nsp)**2
        call pvaioc(s_ctrl%initc,8,ic,icmap,haveva,jc)
        haveva = haveva .and. associated(s_pot%vintr)
        readva = .not. haveva .or. mode00 == 2 .and. jc /= ic
        readva = readva .and. associated(s_pot%vintr) .and. mode00 < 3
        if (haveva) then
          call dpscop(s_pot%vintr,lva,k,1+(jc-1)*k,1,1d0)
          if (ic /= jc) call awrit0('%a va,',outs1,len(outs1),0)
        endif

        havev = .false.; readv = .false.
C#ifdefC READPOT
C        if (mod(mode/100,10) > 0) then
C        k = 0
C        if (associated(s_pot%v0)) k = size(s_pot%v0(:,ic))
C        call pvaioc(s_ctrl%initc,64,ic,icmap,havev,jc)
C        havev = havev .and. associated(s_pot%v0)
C        readv = .not. havev .or. mode00 == 2 .and. jc /= ic
C        readv = readv .and. associated(s_pot%v0) .and. mode00 < 3
C        if (havev) then
C          call dpscop(s_pot%v0,lva,k,1+(jc-1)*k,1,1d0)
C          if (ic /= jc) call awrit0('%a va,',outs1,len(outs1),0)
C        endif
C        endif
C#endif
        k = nl**2*(2*nl-1)*3*nsp
        call pvaioc(s_ctrl%initc,16,ic,icmap,havemp,jc)
        havemp = havemp .and. associated(s_pot%pmpol)
        readmp = .not. havemp .or. mode00 == 2 .and. jc /= ic
        readmp = readmp .and. associated(s_pot%pmpol) .and. mode00 < 3
        if (havemp) then
          call dpscop(s_pot%pmpol,lmp,k,1+(jc-1)*k,1,1d0)
          if (ic /= jc) call awrit0('%a mp,',outs1,len(outs1),0)
        endif
        call pvaioc(s_ctrl%initc,4,ic,icmap,haveso,jc)
        readso = .not. haveso .or. mode00 == 2 .and. jc /= ic
C       readso = readso .and. associated(s_pot%sop) .and. mode00 < 3
        readso = readso .and. mode00 < 3
        if (associated(s_pot%sop)) then
          k = size(s_pot%sop)
        else
          k = 0
        endif
        if (k <= 1) then; haveso=.false.; readso=.false.; endif
        k = nl*nsp*nsp*9
        if (haveso) then
          call dpscop(s_pot%sop,lsop,k,1+(jc-1)*k,1,1d0)
          if (ic /= jc) call awrit0('%a sop,',outs1,len(outs1),0)
        endif
        k = 8*nl*nsp
        call pvaioc(s_ctrl%initc,32,ic,icmap,haveop,jc)
        haveop = haveop .and. associated(s_pot%grrme)
        readop = .not. haveop .or. mode00 == 2 .and. jc /= ic
        readop = readop .and. associated(s_pot%grrme) .and. mode00 < 3
        if (haveop) then
          call dpscop(s_pot%grrme,lgrme,k,1+(jc-1)*k,1,1d0)
          if (ic /= jc) call awrit0('%a pp,',outs1,len(outs1),0)
        endif

C   --- File WRITE ---
        if (lio) then
          lgen = .false.
          lpot = .false.
          lcor = .false.
C     ... Pick up GEN and POT, if available, to save again
          if (scat(iabs(ifi),'GEN:',':',.true.)) then
            lgen = aiogen(alabel,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,dq,
     .        vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          endif
          if (scat(iabs(ifi),'POT:',':',.true.)) then
            read(ifi,102) nr,nsp,a,rmax
  102       format(2i5,2f12.5)
            lpot = aiopot(nr,nsp,a,rmax,-99d0,pot,ifi)
          endif
          lcor = aiocor(nr,nsp,a,rmxx,cor,sumec,sumtc,ifi)

          rewind ifi
          jfi = -ifi
          if (lgen) sw = aiogen(clabl,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,
     .      dq,vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,jfi)
          if (havepq) sw = aiomom(clabl,pnuloc,qnuloc,qnurloc,idmod,nl,lmx,nsp,z,rhrmx,vrmax,bxc,jfi)
          if (havepp) sw = aiopar(clabl,lrel,lpp,lppr,ves,nl,lmx,nsp,jfi)
          if (haveva) sw = aiova(clabl,lva,nl,lmx,nsp,jfi)
          if (havemp) sw = aiomp(clabl,lmp,nl,2*nl-2,nsp,jfi)
          if (haveso) sw = aiosop(clabl,lsop,nl,lmx,nsp,jfi)
          if (haveop) sw = aiorme(clabl,lgrme,nl,nsp,jfi)
          if (lpot)   sw = aiopot(nr,nsp,a,rmax,-99d0,pot,jfi)
          if (lcor) lcor = aiocor(nr,nsp,a,rmxx,cor,sumec,sumtc,jfi)

C   --- File READ ---
        else

C     ... Copy whatever is available on disk to holding arrays
          if (readpq .or. mode0 >= 4) then
            if (readpq) then
              rewind ifi
              readpq = aiomom(clabl,pnuloc,qnuloc,qnurloc,idmod,nl,lmx,nsp,z,rhrmx,vrmax,bxc,ifi)
              if (qnurloc(1) > NULLI+1) haverel = .true.
            endif
C           Couldn't read from atom file ; take default values
            if (readpq) call awrit0('%a pq,',outs2,len(outs2),0)
            if (mode0 >= 4 .and. .not. (readpq .or. havepq)) then
C             call dmcpy(pdf,n0,1,pnuloc,nl,1,nl,nsp)
              call paioc2(nsp,nl,n0,pdf,qdf,pnuloc,qnuloc)
              call awrit0('%a pq,',outs3,len(outs2),0)
              call dvset(vrmax,1,2,-.7d0)
              call dvset(bxc,1,3,0d0)
              rhrmx = .1d0
              readpq = .true.
            endif
          endif

          if (havepq.or.readpq) then
          if (lrel == 2) then
            if (qnurloc(1) == NULLI+1) then ! => not yet read ... estimate from qnu
            call info0(2,0,0,' AIOCLS: relativistic charges missing for class '//trim(clabl)//' ... estimate from scalar q')
C           sw = aiomom(clabl,pnuloc,qnuloc,qnurloc,idmod,nl,lmx,nsp,z,rhrmx,vrmax,bxc,-6)
            call qnu2qnur(1,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
C           Equivalent operation, 2 stage process, for testing
C           call qnu2qnur(20,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
C           call qnu2qnur(11,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge

C           Debugging: confirm that qnu -> qnur -> qnu conserves qnu
C            call qnu2qnur(0,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
CC           Equivalent operation, 2 stage process, for testing
CC           call qnu2qnur(20,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
CC           call qnu2qnur(11,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
C            qlm = NULLI; qnuloc=NULLI
C!           call qnu2qnurx(11,nl,nsp,qnuloc,qnurloc) ! Copy back to test
C!           call qnu2qnur(40,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
C!           Equivalent operation, 2 stage process
C            call qnu2qnur(50,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
C            call qnu2qnur(60,nl,nsp,qnuloc,qlm,qnurloc) ! Map, zeroing out off-digonal parts Q2 to conserve charge
C            sw = aiomom(clabl,pnuloc,qnuloc,qnurloc,idmod,nl,lmx,nsp,z,rhrmx,vrmax,bxc,-6)
C            stop
            endif

            call awrit0('%a qlr,',outs3,len(outs3),0)
          elseif (qnurloc(1) == NULLI+1) then
            qnurloc(1) = NULLI
          endif
          endif

          if (readpp) then
            readpp = aiopar(clabl,lrel,lpp,lppr,ves,nl,lmx,nsp,ifi)
            if (readpp) call awrit0('%a pp,',outs2,len(outs2),0)
            if (lrel == 2) then
              if (readpp) then
                call qnu2enu(4,nl,lmx,lpp,ves,qnurloc,lppr) ! Copy missing enu into qnur available from pprel
              elseif (aiopar(clabl,1,lpp,lppr,ves,nl,lmx,nsp,ifi)) then
                call qnu2enu(7,nl,lmx,lpp,ves,qnurloc,lppr) ! Copy missing enu into qnur available from pp
              endif
            endif
          endif
          if (readso) then
            readso = aiosop(clabl,lsop,nl,lmx,nsp,ifi)
            if (readso) call awrit0('%a so,',outs2,len(outs2),0)
          endif
          if (readop .and. associated(s_pot%grrme)) then
            readop = aiorme(clabl,lgrme,nl,nsp,ifi)
            if (readop) call awrit0('%a op,',outs2,len(outs2),0)
          endif
          if (readva .and. associated(s_pot%vintr)) then
            readva = aiova(clabl,lva,nl,lmx,nsp,ifi)
            if (readva) call awrit0('%a va,',outs2,len(outs2),0)
          endif
C#ifdefC READPOT
C          if (readv .and. associated(s_pot%v0)) then
C            nr2=0; nsp2=0; a2=0; rmax2=0
C            readv = aiopot(nr2,nsp2,a2,rmax2,bhat,v0,ifi)
C            if (readv) then
C              call awrit0('%a v,',outs2,len(outs2),0)
C              if (nr2*nsp > size(s_pot%v0(:,ic))) call rx('asars: v0 size mismatch')
C              readv = aiopot(nr2,nsp2,a2,rmax2,bhat,v0,ifi)
C            endif
C          endif
C#endif
          if (readmp .and. associated(s_pot%pmpol)) then
            readmp = aiomp(clabl,lmp,nl,2*nl-2,nsp,ifi)
            if (readmp) call awrit0('%a mp,',outs2,len(outs2),0)
          endif

C     ... Update what parameters are available
          s_ctrl%initc(ic) = isw(haveop.or.readop)*32+
     .                    isw(havev .or.readv )*32+
     .                    isw(havemp.or.readmp)*16+
     .                    isw(haveva.or.readva)*8 +
     .                    isw(haveso.or.readso)*4 +
     .                    isw(havepp.or.readpp)*2 +
     .                    isw(havepq.or.readpq)*1 +
     .                    s_ctrl%initc(ic) - bitand(s_ctrl%initc(ic),63)

          k = nl*nsp
          if (havepq .or. readpq) then
            call dpscop(pnuloc,s_pot%pnu,k,1,1+(ic-1)*k,1d0)
            call dpscop(qnuloc,s_pot%qnu,3*k,1,1+(ic-1)*3*k,1d0)
            call dpscop(qnurloc,s_pot%qnur(1,ic),4*kr,1,1,1d0)
            s_pot%rhrmx(jc) = rhrmx
C           call dpscop(rhrmx,s_pot%rhrmx,1,1,ic,1d0)
            call dpscop(vrmax,s_pot%vrmax,2,1,2*ic-1,1d0)
C           call dmscop(s_pot%bxc,3,bxc,3,1,3,1,1,1,ic,1d0)
          else
            haverel = .false.
          endif

          k = 6*nl*nsp
          if (havepp .or. readpp) then
            call dpscop(lpp,s_pot%pp,k,1,1+(ic-1)*k,1d0)
            if (lrel == 2) then
              k = 5*nl*2*nl*2*2
              call dpscop(lppr,s_pot%pprel,k,1,1+(ic-1)*k,1d0)
            endif
            call dpscop(ves,s_pot%ves,1,1,ic,1d0)
          endif
          k = (nl*nsp)**2
          if ((haveva .or. readva) .and. associated(s_pot%vintr))
     .      call dpscop(lva,s_pot%vintr,k,1,1+(ic-1)*k,1d0)
C#ifdefC READPOT
C          if ((havev .or. readv) .and. associated(s_pot%v0)) then
C            call dcopy(nr2*nsp,v0,1,s_pot%v0(1,ic),1)
C          endif
C#endif
          k = nl**2*(2*nl-1)*3*nsp
          if ((havemp .or. readmp) .and. associated(s_pot%pmpol))
     .      call dpscop(lmp,s_pot%pmpol,k,1,1+(ic-1)*k,1d0)
          k = nl*nsp*nsp*9
          if (haveso .or. readso)
     .      call dpscop(lsop,s_pot%sop,k,1,1+(ic-1)*k,1d0)
          k = 8*nl*nsp
          if ((haveop .or. readop) .and. associated(s_pot%grrme))
     .      call dpscop(lgrme,s_pot%grrme,k,1,1+(ic-1)*k,1d0)
        endif

        if (mode00 < 3) call fclr(clabl,ifi)
        if (iprint() > 40) then
          outs = ' '
          if (outs1 /= ' ') call awrit1('%x '//clabl//'%a: copied '//outs1//'%a%b from '//'class %i',outs,len(outs),0,jc)
          if (outs2 /= ' ') call awrit0('%x '//clabl//'%a: read '//outs2//'%a%b from '//'disk',outs,len(outs),0)
          if (outs3 /= ' ') then
            if (outs2 == ' ') then
              call awrit0('%x '//clabl//'%a: use defaults for: '//outs3//'%a%b',outs,len(outs),0)
            else
              call awrit0('%a; use defaults for: '//outs3//'%a%b',outs,len(outs),0)
            endif
          endif
          if (outs == ' ') call awrit0(' '//clabl//'%a: nothing read',outs,len(outs),0)
          call awrit0(' aiocls class'//outs,' ',-len(outs),lgunit(1))
        endif
C       Set flag if all relativistic q's were read or copied
        if (.not. lio .and. haverel) s_ctrl%stpnrho = .true.
      enddo
C     call awrit2('%n:1i',' ',100,6,nclasp,s_ctrl%initc)
      deallocate(lpp,lppr,lsop,lgrme,lva,lmp,pot,cor)
C     End of MPI master-only branch
      endif

C --- MPI broadcast everything passed out of aiocls ---
      call mpibc1(s_ctrl%initc,nclspd,2,.false.,'aiocls','initc')
      call mpibc1(zbak,2,4,.false.,'aiocls','zbak')
      call mpibc1(s_pot%pp,6*nl*nsp*nclspd,4,.false.,'aiocls','pp')
      if (lrel == 2) then
        k = 5*nl*2*nl*2*2
        call mpibc1(s_pot%pprel,k*nclspd,4,.false.,'aiocls','ppr')
        k = 4*nl*2*nl*2*2
        call mpibc1(s_pot%qnur,k*nclspd,4,.false.,'aiocls','qnur')
      endif
      call mpibc1(haveso,1,1,.false.,'aiocls','haveso')
      call mpibc1(readso,1,1,.false.,'aiocls','haveso')
      if (haveso .or. readso) then
        k = 9*nl*nsp*nsp*nclspd
        call mpibc1(s_pot%sop,k,4,.false.,'aiocls','sop')
      endif
      if (associated(s_pot%grrme)) then
        call mpibc1(s_pot%grrme,8*nl*nsp*nclspd,4,.false.,'aiocls','grrme')
      endif
C#ifdefC READPOT
C      if (associated(s_pot%v0)) then
C        call mpibc1(s_pot%v0,size(s_pot%v0),4,.false.,'aiocls','v0')
C      endif
C#endif
      call mpibc1(s_pot%pnu,nl*nsp*nclspd,4,.false.,'aiocls','pnu')
      call mpibc1(s_pot%qnu,3*nl*nsp*nclspd,4,.false.,'aiocls','qnu')
      call mpibc1(s_pot%rhrmx,nclspd,4,.false.,'aiocls','rhrmx')
      nclspp = 2*nclasp-nclass
      nclsppd = nclspp + nangl
      call mpibc1(s_pot%vrmax,2*nclsppd,4,.false.,'aiocls','vrmax')
C     call mpibc1(s_pot%bxc,3*nclasp,4,.false.,'aiocls','bxc')
      if (associated(s_pot%vintr)) then
        k = nclasp*(nl*nsp)**2
        call mpibc1(s_pot%vintr,k,4,.false.,'aiocls','vintr')
      endif
C     if (lgors('ctrl lasa,32',sctrl)) then
      if (associated(s_pot%pmpol)) then
        k = (2*nl-1)*nl**2*3*nsp*nclasp
        call mpibc1(s_pot%pmpol,k,4,.false.,'aiocls','pmpol')
      endif
      call mpibc1(s_pot%ves,nclsppd,4,.false.,'aiocls','ves')

C --- Assemble background rho (rmt) ---
C     if (mod(mode/10,10) /= 0 .and. lgors('ctrl lasa,64',sctrl)) then
      if (mod(mode/10,10) /= 0 .and. IAND(s_ctrl%lasa,64) /= 0) then
        zbak = s_ctrl%zbak
        nbas = s_ctrl%nbas
        nclass = s_ctrl%nclass
C        call awrit2('%n:1d',' ',100,6,nclass,s_pot%rhrmx)
C        call awrit2('%n:1i',' ',100,6,nclasp,s_ctrl%nrc)
        if (zbak(2) == 0) then
          do  ic = 1, nclass
            z = s_ctrl%dclabl(ic)
            call r8tos8(z,clabl)
            z = s_pot%rhrmx(ic)
            if (iprint() >= 30 .and. z == 0) print *,
     .        ' (AIOCLS, mtcor): rho(rmax) for class ',clabl,' is 0'
            zbak(2) = zbak(2) + (s_ctrl%nrc(ic)*z)/nbas
          enddo
          zbak(2) = zbak(2) * s_lat%vol
        endif
        s_ctrl%zbak = zbak
      endif

C     call prmx('aiocls2 qnu',s_pot%qnu,3*nl*nsp,3*nl*nsp,nclass)

      end

      subroutine prtcls(mode,s_ctrl,s_pot,s_spec,ic1,ic2,ifi)
C- Prints portions of ASA class data
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  Governs what to write:
Ci         :1   write gen data
Ci         :2   write P,Q
Ci         :4   write PPARs
Ci         :8   write intra-atomic density-density response matrix
Ci         :16  write multipole moments
Ci         :32  write spin-orbit coupling parameters.
Ci         :64  write optical matrix elements
Ci         :128 write spherical potential
Ci         :256 write core density
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp dclabl ics
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: pnu qnu rhrmx pp pprel osop grrme ovintr pmpol
Ci                ves
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idmod
Ci   ic1,ic2:range of classes to read data
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Co Outputs
Co   Some parameters are written to ifi
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Aug 09 Adapted from aiocls
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ifi,mode,ic1,ic2
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer,parameter :: n0=10, NULLI=-99999
      character*8 clabl
      logical lrel2,lgen,lpq,lpp,lva,lmp,lso,lpot,lcor,lop,sw,aiomom,aiopar,aiova,aiomp,aiosop,aiorme
      integer nglob,nclasp,nl,nsp,k,ic,is,lmx,idmod(n0),lrel,kr
      double precision rhrmx,vrmax(2),bxc(3),wk(100),wk2(1000),wkn(1000),ves
      real(8), pointer :: llmp(:)

C --- Setup ---
      nl = nglob('nl')
      nsp = nglob('nsp')
      lrel = mod(nglob('lrel'),10)
      nclasp = s_ctrl%nclasp

      lgen = mod(mode,2) /= 0
      lpq  = mod(mode/2,2) /= 0
      lpp  = mod(mode/4,2) /= 0
      lva  = mod(mode/8,2) /= 0
      lmp  = mod(mode/16,2) /= 0
      lso  = mod(mode/32,2) /= 0
      lop  = mod(mode/64,2) /= 0
      lpot = mod(mode/128,2) /= 0
      lcor = mod(mode/256,2) /= 0
      lrel2 = s_pot%qnur(1,1) == NULLI

      do   ic = ic1, ic2
        call r8tos8(s_ctrl%dclabl(ic),clabl)
        is = s_ctrl%ics(ic)
        lmx = s_spec(is)%lmxa
        idmod = s_spec(is)%idmod
        if (lgen) then
          call rx('prtcls not ready to write gen data')
        endif
        if (lpq) then
          k = nl*nsp
          kr = 8*nl*nl
          call dpscop(s_pot%pnu,wk,k,1+(ic-1)*k,1,1d0)
          call dpscop(s_pot%qnu,wk2,3*k,1+(ic-1)*3*k,1,1d0)
          if (lrel2) stop 'aiocls'
          if (lrel2) call dpscop(s_pot%qnur,wkn,4*kr,1+(ic-1)*4*kr,1,1d0)
C         call dpscop(s_pot%rhrmx,rhrmx,1,ic,1,1d0)
          rhrmx = s_pot%rhrmx(ic)
          call dpscop(s_pot%vrmax,vrmax,2,2*ic-1,1,1d0)
C         call dmscop(bxc,3,s_pot%bxc,3,1,3,ic,ic,1,1,1d0)
          stop 'aiomom aiocls'
          sw = aiomom(clabl,wk,wk2,wkn,idmod,nl,lmx,nsp,0d0,rhrmx,vrmax,bxc,-ifi)
        endif
        if (lpp) then
          k = 6*nl*nsp
          call dpscop(s_pot%pp,wk,k,1+(ic-1)*k,1,1d0)
          if (lrel == 2) then
            k = 5*nl*2*nl*2*2
            call dpscop(s_pot%pprel,wk2,k,1+(ic-1)*k,1,1d0)
          endif
          call dpscop(s_pot%ves,ves,1,ic,1,1d0)
          sw = aiopar(clabl,lrel,wk,wk2,ves,nl,lmx,nsp,-ifi)
        endif
        if (lva) then
          k = (nl*nsp)**2
          call dpscop(s_pot%vintr,wk,k,1+(ic-1)*k,1,1d0)
          sw = aiova(clabl,wk,nl,lmx,nsp,-ifi)
        endif
        if (lmp) then
          allocate(llmp(3000))
          k = nl**2*(2*nl-1)*3*nsp
          call dpscop(s_pot%pmpol,llmp,k,1+(ic-1)*k,1,1d0)
          sw = aiomp(clabl,llmp,nl,2*nl-2,nsp,-ifi)
          deallocate(llmp)
        endif
        if (lso) then
          k = nl*nsp*nsp*9
          call dpscop(s_pot%sop,wk2,k,1+(ic-1)*k,1,1d0)
          sw = aiosop(clabl,wk2,nl,lmx,nsp,-ifi)
        endif
        if (lop) then
          k = 8*nl*nsp
          call dpscop(s_pot%grrme,wk2,k,1+(ic-1)*k,1,1d0)
          sw = aiorme(clabl,wk2,nl,nsp,-ifi)
        endif
        if (lpot) then
          call rx('prtcls not ready to write pot data')
C         sw = aiopot(nr,nsp,a,rmax,-99d0,pot,-ifi)
        endif
        if (lcor) then
          call rx('prtcls not ready to write cor data')
C         lcor = aiocor(nr,nsp,a,rmxx,cor,sumec,sumtc,-ifi)
        endif

      enddo

      end


      subroutine pvaioc(initc,mask,ic0,icmap,lhave,ic)
C- Find whether data available either in class or mapped class
      implicit none
      logical lhave
      integer initc(*),mask,ic0,icmap
      integer ic
      ic = ic0
      lhave = mod(initc(ic)/mask,2) == 1
      if (.not. lhave .and. icmap /= 0) then
        lhave = mod(initc(icmap)/mask,2) == 1
        ic = icmap
      endif
      end
      subroutine paioc2(nsp,nl,n0,pat,qat,pnu,qnu)
C- Widget to copy pat,qat to pnu,qnu
      implicit none
      integer n0,nl,nsp
      double precision pat(n0,2),qat(n0,2),pnu(nl,nsp),qnu(3,nl,nsp)
      integer i,il

      do  i = 1, nsp
      do  il = 1, nl
C       pnu(il,i) = int(pat(il,i)) + .5d0
        pnu(il,i) = pat(il,i)
        qnu(1,il,i) = qat(il,1)/nsp
        if (nsp == 2) then
          if (pat(il,i) == 0) pnu(il,i) = pnu(il,1)
          qnu(1,il,i) = qat(il,1)/nsp + qat(il,2)/2*dble(3-2*i)
        endif
        qnu(2,il,i) = 0d0
        qnu(3,il,i) = 0d0
      enddo
      enddo
      end
      subroutine qnu2qnur(job,nl,nsp,qnu,qlm,qnur)
C- Estimate relativistic charges qnur from scalar relativistic qnu
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :-1 or -4 : zero out qnur(0:2) and return
Ci         :-1 leaves qnur(3) untouched while -1 sets qnur(3) to NULLI
Ci         :If job>=0:
Ci         :1s digit:
Ci         :1 set qnur(2,:,1,2) = qnur(2,:,2,1) = 0 (forward direction)
Ci         :  Omit qnur(2,:,1,2) = qnur(2,:,2,1) in conversion (reverse direction)
Ci         :  This conserves sphere charge but destroys exact mapping qnu <-> qnur
Ci         :2 set qnur(:,:,1,2) = qnur(:,:,2,1) = 0 (forward direction)
Ci         :  (not implemented)
Ci         :4 set qnur(3,:,:,:,:) = NULLI
Ci         :10s digit:
Ci         :0 convert qnu to qnur
Ci         :1 convert qlm to qnur
Ci         :2 convert qnu to qlm
Ci         :4 Add 4 to reverse the sense of conversion
Ci   qnu   :Scalar relativistic energy-weighted moments of the sphere charges
Ci   qlm   :Scalar relativistic energy-weighted moments of the sphere charges, resolved by m
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Outputs
Co   qnur  :Relativistic energy-weighted moments of the sphere charges
Cu Updates
Cu   09 Apr 15 Adapted from Kirill, with new job switch
Cu   09 Jun 14 (Kirill) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,nl,nsp
      double precision qnu(0:2,nl,nsp),qlm(0:2,nl*nl,nsp),qnur(0:3,0:nl-1,2*nl,2,2)
C ... Local parameters
      real(8) ql
      complex(8) qmat(nl*nl,2,nl*nl,2),qmrel(nl*nl*2,nl*nl*2),xx
C     complex(8) zsum
      integer idum(1),norb,idx(2,nl,2*nl),l,imu,ilm,ll,m,k1,k2,iprint,job1
      integer,parameter :: NULLI=-99999

      job1 = mod(job/10,10)
      if (job == -1 .or. job == -4 .or. job1 == 0 .or. job1 == 1) qnur(0:2,:,:,:,:) = 0
      if (job == -4 .or. mod(job,10) >= 4) qnur(3,:,:,:,:) = NULLI
      if (job == -1 .or. job == -4) return

      norb = nl*nl
      call dpzero(qmat,2*size(qmat)); call dpzero(qmrel,2*size(qmrel))
      if (job1 == 2 .or. job1 == 5) call dpzero(qlm,size(qlm)) ! Not necessary ?
      if (job1 == 4 .or. job1 == 6) call dpzero(qnu,size(qnu)) ! Not necessary ?

      call mstokm(2,nl,1,norb,xx,xx,idx) ! kap_1 = idx(1,l,imu) and kap_2 = idx(2,l,imu)

C ... Loop over 0th, 1st and 2nd moments
      do  m = 0, 2
C       Convert qnu to qlm or qmat
        if (job1 <= 2 .or. job1 == 6) then  ! qnu -> qmat or qlm -> qmat
          do  ilm = 1, norb
            l = ll(ilm)
            if (job1 == 0) then
              qmat(ilm,1,ilm,1) = qnu(m,l+1,1)/(2*l+1)
              qmat(ilm,2,ilm,2) = qnu(m,l+1,2)/(2*l+1)
            else if (job1 == 1) then
              qmat(ilm,1,ilm,1) = qlm(m,ilm,1)
              qmat(ilm,2,ilm,2) = qlm(m,ilm,2)
            else if (job1 == 2) then
              qlm(m,ilm,1) = qnu(m,l+1,1)/(2*l+1)
              qlm(m,ilm,2) = qnu(m,l+1,2)/(2*l+1)
            else
              qnu(m,l+1,1) = qnu(m,l+1,1) + qlm(m,ilm,1)
              qnu(m,l+1,2) = qnu(m,l+1,2) + qlm(m,ilm,2)
            endif
          enddo
          if (job1 == 2 .or. job1 == 6) cycle
          call mstokm(0,nl,1,norb,qmat,qmrel,idum) ! Convert q from (l,m) to (kappa,mu)
        endif
!       call prtarr(qmat,norb) !in gf/cpadlm.f

        do  l = 0, nl-1
          do  imu = 1, 2*(l+1)
            k1 = idx(1,l+1,imu) ; k2 = idx(2,l+1,imu)
            if (job1 < 4) then
              if (k1 == 0) then
                qnur(m,l,imu,2,2) = qmrel(k2,k2)
              else
                qnur(m,l,imu,1,1) = qmrel(k1,k1)
                qnur(m,l,imu,2,2) = qmrel(k2,k2)
C               Sphere charge is not conserved because:
C               p(la1,la2) = int dr phidot(r,la1)*phidot(r,la2) is nonzero
C                     while  int dr phi(r,la1)*phi(r,la2) is zero
                if ((m /= 2 .or. mod(job,10) /= 1) .and. mod(job,10) /= 2) then
                  qnur(m,l,imu,2,1) = qmrel(k2,k1)
                  qnur(m,l,imu,1,2) = qmrel(k1,k2)
                endif
              endif
            else
              if (k1 == 0) then
                qmrel(k2,k2) = qnur(m,l,imu,2,2)
              else
                qmrel(k1,k1) = qnur(m,l,imu,1,1)
                qmrel(k2,k2) = qnur(m,l,imu,2,2)
                if ((m /= 2 .or. mod(job,10) /= 1) .and. mod(job,10) /= 2) then
                qmrel(k1,k2) = qnur(m,l,imu,1,2)
                qmrel(k2,k1) = qnur(m,l,imu,2,1)
                endif
              endif
            endif
          enddo

          if (job1 < 4 .and. m == 0 .and. iprint() > 40) then
            ql = 0
            do  imu = 1, 2*(l+1)
              ql = ql + sum(qnur(m,l,imu,:,:))
            enddo
            if (ql /= 0) then
            call info2(40,0,0,' Estimate relativistic charges for l=%i  q=%,6d',l,ql)
            do  imu = 1, 2*(l+1)
              print "(f5.1,4f12.6)", dble(imu-l)-1.5d0, qnur(m,l,imu,:,:)
            enddo
          endif
          endif
        enddo

        if (job1 >= 4) then
          call mstokm(1,nl,1,norb,qmat,qmrel,idum) ! qmrel -> qmat
          do  ilm = 1, norb
            if (job1 == 4) then
              l = ll(ilm)
              qnu(m,l+1,1) = qmat(ilm,1,ilm,1)*(2*l+1)
              qnu(m,l+1,2) = qmat(ilm,2,ilm,2)*(2*l+1)
            elseif (job1 == 5) then
              qlm(m,ilm,1) = qmat(ilm,1,ilm,1)
              qlm(m,ilm,2) = qmat(ilm,2,ilm,2)
            endif
          enddo
        endif

      enddo

      if (job1 /= 2 .and. job1 /= 5 .and. job1 /= 6 .and. iprint() > 40) then
        call info2(40,0,0,' Scalar Dirac sphere charge=%d  Dirac sphere charge=%d',
     .    sum(qnu(0,:,:)),sum(qnur(0,:,:,1,1)+qnur(0,:,:,2,2)))
      endif

      end

      subroutine qnu2enu(mode,nl,lmax,pp,ves,qnur,pprel)
C- Copy qnur(3,...) to pprel(5,..) or vise-versa
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 copy all qnur(4,:,:,:,:) to pprel(5,:,:,:,:)
Ci         :1 copy all pprel(5,:,:,:,:) to qnur(4,:,:,:,:)
Ci         :2 copy all pp(1,:,:,:,:) to pprel(5,:,:,:,:)
Ci         :3 conditionally copy qnur(4,:,:,:,:)  to pprel(5,:,:,:,:) if pprel(5,:,:,1,1) = NULL
Ci         :4 conditionally copy pprel(5,:,:,:,:) to qnur(4,:,:,:,:)  if qnur(4,:,:,:,:)  = NULL
Ci         :5 conditionally copy pp(5,:,:)        to pprel(5,:,:,:,:) if pprel(5,:,:,1,1) = NULL
Ci         :6 copy all pp(1,:,:,:,:) to qnur(4,:,:,:,:)
Ci         :7 conditionally copy pp(1,:,:,:,:)    to qnur(4,:,:,:,:) if qnur(4,:,:,:,:) = NULL
Ci   nl    :(global maximum l) + 1
Ci   lmax  :maximum l for a given site
Ci   qnur
Ci   pprel
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   09 Apr 15
C ----------------------------------------------------------------------
      implicit none
      integer mode,nl,lmax
      double precision ves,qnur(4,0:nl-1,2*nl,2,2),pprel(5,0:nl-1,2*nl,2,2),pp(6,0:nl-1,2)
      integer l,imu
      integer,parameter :: NULLI=-99999

      do  l = 0, lmax
        do  imu = 1, 2*l+2
C          if (mode == 0 .or. mode == 2 .and. pprel(5,l,imu,1,1) == NULLI) pprel(5,l,imu,:,:) = qnur(4,l,imu,1,1)
C          if (mode == 1 .or. mode == 4 .and.  qnur(4,l,imu,1,1) == NULLI) qnur(4,l,imu,:,:) = pprel(5,l,imu,1,1)

          if (mode == 0 .or. mode == 3 .and. pprel(5,l,imu,1,1) == NULLI) pprel(5,l,imu,:,:) = qnur(4,l,imu,1,1)
          if (mode == 1 .or. mode == 4 .and.  qnur(4,l,imu,1,1) == NULLI) qnur(4,l,imu,:,:)  = pprel(5,l,imu,1,1) - ves
          if (mode == 2 .or. mode == 5 .and. pprel(5,l,imu,1,1) == NULLI) pprel(5,l,imu,:,:) = (pp(1,l,1) + pp(1,l,2))/2
          if (mode == 6 .or. mode == 7 .and.  qnur(4,l,imu,1,1) == NULLI) qnur(4,l,imu,:,:)  = (pp(1,l,1) + pp(1,l,2))/2 - ves

        enddo
      enddo
      end
