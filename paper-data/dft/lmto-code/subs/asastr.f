      subroutine asastr(prgnam,s_ctrl,s_ham,s_lat,s_str,s_spec,s_site)
C- Screened, real-space structure constants
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nspec nl nspin lpgf nclasp npadl npadr nbasp ips
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pgfsl ips
Cio    Passed to:  chkstr chksg
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat avw pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cy cg indxcg jcg pos
Cio    Passed to:  chkstr chksg
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  nkaps kaps rmax rfit lmaxw drwats iinv loka ivl
Ci                 rmaxg mxnbr amode nitab npr
Co     Stored:     kaps nttab nitab
Co     Allocated:  alp adot s npr iax sdot
Cio    Elts passed:kaps amode lequiv lshow nkaps npr s iax alph nitab
Cio    Passed to:  chkstr chksg
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxb hcr orbp kmxt rsma ehvl lmxa alpha rmt hsfitk
Co     Stored:     hcr
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hcr2a maadot chkstr e0parm chksg
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   prgnam:name of main program
Co Outputs
Co   strx are written to disk
Cl Local variables
Cl   amode :determines kind of strux to make:
Cl         :  0 2nd generation kappa=0
Cl         :    OKA conventions
Cl         :  1 SSSW, 1 kappa val, or 2 kappa val-slope
Cl         :    MSM conventions
Cl         :  2 NMTO, 1 kappa val, nkapn of them
Cl         :    OKA conventions
Cl         : 11 unit-zero val-slo Hankel SSWs, 2 (linked) energies
Cl         : 13 Similar to 11, but also Gaussian val-lap strux
Cl   rmaxs :sphere radius for structure constants, a.u.
Cl   nkaps :Number of functions per screened function
Cl         :with corresponding number of conditions
Cl   nkapn :Number of families of functions to make
Cl         (1 for now, nkap for NMTO)
Cl   ncplx :1 if strux are real; 2 if complex
Cr Remarks
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   01 Sep 11 Begin migration to f90 structures
Cu   04 Jul 08 (S. Lozovoi) val-lap structure constants
Cu   06 Aug 06 New 2-kappa strux
Cu   19 May 04 make strux work with Methfessel conventions for H,J
Cu   31 Aug 00 some adaptation to fit with NMTO
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*(*) prgnam*8
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_ham):: s_ham
      type(str_lat):: s_lat
      type(str_str):: s_str
      type(str_spec):: s_spec(*)
      type(str_site):: s_site(*)
      type(str_str0):: s_sta,s_stad,s_stag
C ... Dynamically allocated local arrays
      integer,allocatable :: lmxb(:)
C#ifndef LINUXF
      integer, allocatable, target :: iax(:)
C#elseC
C      integer, pointer :: iax(:) => null()
C#endif
      real(8),allocatable :: hcr(:,:)
      real(8),pointer :: pos(:,:)
      real(8),pointer :: alpha1(:),adot1(:)
      real(8),target,allocatable :: hcrl(:)
      real(8), allocatable :: rsmh(:),rsma(:)
      integer, allocatable :: kmx(:)
      real(8), allocatable :: tral(:),trad(:),rwats(:),ehvl(:)
C ... Local parameters
      character*120 outs
      logical cmdopt,strdif,ltmp,lsdot,iostr
      integer amode,awrite,fopn,fopno,i,iprint,isw,is,
     .  itral,j,lmaxw,lpbc,lpgf,mordrn,mxcsiz,n0,nbas,
     .  nbasp,nbaspp,nclasp,niax,nitab,nkaps,nl,nla,npadl,npadr,
     .  nttabg,liog,ivl,nsp,nspec,nttab,stdo,streqv,loka,nkapn,nds,lio,
     .  ncplx,nmto
      parameter (niax=10,n0=10)
      double precision avw,kap2(n0),rmaxs,rfit,plat(3,3),alat,
     .  drwats,xx,ckbas,cksumf,siinv(5),rmxg0,rmxg
      procedure(integer) :: nglob

C --- Unpack global variables ---
      stdo = nglob('stdo')
      nbas = s_ctrl%nbas
      nspec = s_ctrl%nspec
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      lpgf  = s_ctrl%lpgf(1)
C     lpbc  = 0 for pbc in 3 dimensions, 11 pgf padded geometry
      lpbc  = 0
      if (lpgf > 0) lpbc = 11
      nclasp = s_ctrl%nclasp
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr

      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas
      mordrn = 0
      ncplx = 1

C ... Global strux data
      nkaps = s_str%nkaps
      i = size(s_str%kaps)
      kap2(1:i) = s_str%kaps(1:i)
      rmaxs = s_str%rmax
      rfit = s_str%rfit
      amode = mod(s_str%amode,100)
      lmaxw = s_str%lmaxw
      drwats = s_str%drwats
      siinv(1:3) = s_str%iinv(1:3)
      loka = s_str%loka
      alat = s_lat%alat
      plat = s_lat%plat
      avw = s_lat%avw

      allocate(lmxb(nspec))
      allocate(hcr(nl,nspec))
      do  is = 1, nspec
        lmxb(is) = s_spec(is)%lmxb
        hcr(1:nl,is) = s_spec(is)%hcr(1:nl)
      enddo
      pos => s_lat%pos

C --- Parameters determining kind and number of strux ---
      call sanrg(.true.,amode,0,13,'ASASTR:','strux mode')
C     NMTO uses Tank's tral conventions instead of alpha
      if (amode == 2) then
        itral = 4
        lsdot = .true.
        nkapn = nkaps
        nkaps = 1
        nmto = 1
C     2rd generation only for kappa=0, screening spec'd by alpha
      elseif (amode == 0) then
        itral = 0
        nkaps = 1
        nkapn = 1
        kap2(1) = 0
        s_str%kaps(1) = kap2(1)
        lsdot = .true.
        nmto = 0
C       Find hcr in terms of alpha, for printout
        call hcr2a(s_spec)
C     All remaining modes
      else
        itral = 0
        lsdot = .false.   ! For now
        if (amode == 11 .or. amode == 13) then
          nkapn = 1
        else
          nkapn = nkaps
          nkaps = 1
        endif
        nmto = 0
        if (amode == 13) then
          ivl = s_str%ivl
          rmxg0 = s_str%rmaxg
          call sanrg(.true.,ivl,0,2,'ASASTR:','vl mode')
          call fsanrg(rmxg0,0d0,1d0,0d0,'ASASTR:','rmaxg/rmax',.true.)
c            print *,' asastr: ivl, rmxg0 =',ivl, rmxg0
        endif
      endif

C --- Compare structure constants ---
      if (cmdopt('--chk',5,0,outs)) then
        print *, prgnam // 'comparing file strx STR and STR1 ...'
        i = 16+8+4+1
        ltmp = strdif('STR','STR1')
        return
      endif

C --- Printout ---
      if (amode == 0) then
        call info2(10,1,0,' ASASTR: strux for 2nd gen LMTOs'//
     .    '%?#n#, and Sdot##'//
     .    '%?#n<0##%-1j%N%9fWatson sphere lmax=%i#',
     .    isw(lsdot),lmaxw)
      elseif (iprint() >= 10 .and. amode == 2) then
        call awrit4('%N ASASTR: NMTO S'//
     .    '%?#n#, Sdot## for'//
     .    '%?#n>1#%-1j %i energies%j# kap^2=%d#.'//
     .    '  Watson sphere%?#n<0#: none#%-1j lmax=%i#.',
     .    ' ',80,stdo,isw(lsdot),nkapn,kap2,lmaxw)
        i = awrite('%x%9fitral = %i:',outs,80,0,itral,xx,xx,xx,xx,
     .    xx,xx,xx)
        call mktrli(itral,outs(i+3:))
        call awrit0('%a',outs,80,-stdo)
      elseif (amode == 1 .or. amode == 11 .or. amode == 13) then
        call info5(10,1,0,' ASASTR: strux for Hankel-based '//
     .    '%?#n==1#unit-zero val#val-slo#%-1j'//
     .    '%?#n==13# and Gaussian-based val-lap##'//
     .    ' envelopes'//
     .    '%?#n#, and Sdot##'//
     .    '%?#n<0##%-1j%N%9fWatson sphere lmax=%i#',
     .    amode,isw(lsdot),lmaxw,0,0)
      endif
      if (amode > 0) then
        call info5(10,0,0,'%9fEnergies: %n:1d a.u.  avw = %,6;6d',
     .    nkaps*nkapn,kap2,avw,0,0)
      else
        call awrit1('         avw = %,6;6d',' ',80,stdo,avw)
      endif

C --- Default rmaxs, if input is not positive ---
      if (rmaxs <= 0d0) then
        rmaxs = 2.7d0*avw
        call info5(10,0,0,'%9f'//
     .    'Use default rmaxs = %;3d a.u. = %;3d*avw = %;3d*alat',
     .    rmaxs,rmaxs/avw,rmaxs/alat,0,0)
      endif

      if (siinv(1) /= 0)
     .  call awrit3('%9fIterative inversion: ncut = %d  nit = %d'//
     .  '  tol = %g',' ',80,stdo,siinv(1),siinv(2),siinv(3))

C ... Sanity checks on number of energies
      if (amode == 2) then
        call sanrg(.true.,nkapn,2,5,' ','number of energies NEL')
      elseif (amode == 1) then
        call sanrg(.true.,nkapn,1,1,' ','number of energies NEL') ! for now
      elseif (amode == 11 .or. amode == 13) then
        call sanrg(.true.,nkaps,2,2,' ','number of energies NEL')
      else
        call sanrg(.true.,nkaps,1,1,' ','number of energies NEL')
      endif

C --- Get neighbor table iax for each atom in the cluster ---
      if (lpbc == 0) then
        i = 3
        j = -1
C       if (.not. allocated(s_ctrl%pgfsl)) allocate(s_ctrl%pgfsl(1))
      elseif (lpbc == 1 .or. lpbc == 11) then
        i = 2
        j = 1
      else
        call rx('ASASTR: not implemented for lpbc>1')
      endif
C ... Make nttab,ntab,iax
      mxcsiz = s_str%mxnbr
      allocate(s_sta%n(nbasp+1))
      call pairs(nbas,nbasp,alat,plat,[rmaxs/2],pos,
     .  [-1],i,j,s_ctrl%pgfsl,nttab,s_sta%n,iax,mxcsiz)
      s_sta%i  => iax
      s_stad%n => s_sta%n
      s_stad%i => s_sta%i
      s_stag%n => s_sta%n
      s_stag%i => s_sta%i

C     Doctor the table
C     call snot(nttab,s_sta%n,s_sta%i,nbas)
C ... Patch iax(6) for the padded basis
      if (nbasp > nbas) call pairp6(nbas,npadl,npadr,s_sta%i,s_sta%n)
C ... Added pairs for the embedded cluster
      if (mordrn == 1) then
        call rx('pairec commented out for mordn')
C        mxnbr = 2*rmaxs**3*nbasp
C        call redfi(oiax, niax*mxnbr)
C        call pairec(nbas,0,0,s_ctrl%ipc,alat/alat,plate,w(obas),
C     .    w(orham),nttab,s_sta%n,s_sta%i,w(owk),mxcsiz)
C        call redfi(oiax, niax*nttab)
      endif
C     Make iax(9,:)
      call mkiaxd(nttab,lmxb,s_ctrl%ips,s_sta%i)
C     call rlse(olmx)

C     print *, '!! for now, assign nds to nl**2'
      nds = nl**2

C --- Screening parameters for all kappa's ---
      i = nds*nbaspp*nkaps**2*nkapn
      allocate(s_sta%a(i)); call dpzero(s_sta%a,i)
      allocate(s_stad%a(i)); call dpzero(s_stad%a,i)
      if (amode == 2 .or. amode == 0) then
        i = 4*nds*nbaspp*nkapn
        allocate(tral(i)); call dpzero(tral,i)
        allocate(trad(i)); call dpzero(trad,i)
      elseif (amode == 1 .or. amode == 11) then
        allocate(tral(1))
        allocate(trad(1))
      elseif (amode == 13) then
        allocate(tral(1))
        allocate(trad(1))
        allocate(s_stag%ng(nbas))
C       call defi(ontabg, nbas)
c       call defi(okmx,   nspec)
c ... read rsmh
c       do is = 1, nspec
c         call dpzero(orbp,2*n0*n0)
c         call dcopy(n0,orbp(1,1,1),1,rsmh(1,is),1)
c           print *, 'is, lmxb =',is,lmxb
c           print *, 'rsmh =',(rsmh(il,is),il=1,lmxb+1)
c       enddo
      else
        call sanrg(.true.,nkaps,1,1,' ','number of energies NEL')
      endif
      allocate(rwats(nbaspp))

C ... Some temporary arrays needed to generate alpha
      allocate(hcrl(nl*nspec))
      call spec2class(s_spec,nspec,0,'hcr',nl,xx,hcrl)
C     Unpack RSMH
      if (amode == 13) then
        allocate(rsmh(nl*nspec))
        allocate(kmx(nspec))
        allocate(rsma(nspec))
        call spec2class(s_spec,nspec,0,'orbp',nl,xx,rsmh)
        call spec2class(s_spec,nclasp,0,'kmxt',1,kmx,xx)
        call spec2class(s_spec,nspec,0,'rsma',1,xx,rsma)
        if (ivl /= 0) then
          allocate(ehvl(nl*nspec))
          call spec2class(s_spec,nspec,0,'ehvl',nl,xx,ehvl)
        else
          allocate(ehvl(1))
        endif
      endif

C ... Second generation: screening defined through spec->alpha
      if (amode == 0) then
        call rxx(loka /= 1,'2nd gen strux require OKA conventions')
C       Note: maadot and makalp should be merged
C       Need to extract alp from sspec, and pass to makalp
        call maadot(s_spec,nl,nbaspp,amode+10*mod(s_str%amode/100,10),
     .    s_ctrl%ips,avw,s_sta%a,s_stad%a,tral,trad)

C ... 1- and 2- kappa strux: screening defined through spec->hcr
      elseif (amode == 1 .or. amode == 11 .or. amode == 13) then
        call rxx(loka /= 0,'strux require MSM conventions')
        call makalp(nl,nds,nbas,nkaps,kap2,hcr,100*loka+1,
     .    s_ctrl%ips,lmxb,xx,s_sta%a)

c       Make ntab for Gaussians
        if (amode == 13) then
          if (nbas /= nbasp)
     .      call rx('ASASTR: padded basis is not implemented in mode 3')

          rmxg = rmxg0*rmaxs
          call pairg(nbas,0,-1,alat,plat,pos,rmxg/2,
     .      s_sta%n,s_sta%i,s_stag%ng,nttabg)
        endif

C ... NMTO: screening defined through spec->hcr.  Uses OKA now
      elseif (amode == 2) then
        call rxx(loka /= 1,'NMTO strux require OKA conventions')
        call dscal(nl*nspec,1/avw,hcr,1)
        call dscal(nkapn,avw**2,kap2,1)
        call mktra2(11,loka,nbaspp,s_ctrl%ips,nl,lmxb,avw,itral,kap2,
     .    nkapn,hcr,1,tral,trad,s_sta%a,s_stad%a,xx,xx)
      endif

C     For now keep work with w(..) and duplicate in strux
      i = nds*nbaspp*nkaps**2*nkapn
      call ptr_str(s_str,4+1,'alp',i,0,s_sta%a)
      call ptr_str(s_str,4+1,'adot',i,0,s_stad%a)

C --- Watson sphere radii ---
      if (lmaxw >= 0) then
        call mkrwat(alat/avw,pos,nbaspp,s_sta%i,s_sta%n,
     .    s_ctrl%ips,nl,nttab,plat,hcr,drwats,rwats)
C       call awrit2('%n:1d',' ',80,stdo,nbaspp,rwats)
      endif

C ... Count the number of inequivalent strx
      nitab = nttab
C     if (lgors('str lequiv,1',sstr)) then
      if (IAND(s_str%lequiv,1) == 1) then
        nitab = 0
        do  i = 1, nbas
          if (nmto == 0) then
            alpha1 => s_sta%a
            adot1  => s_stad%a
            ltmp = .true.
            nla = nl**2
          else
            alpha1 => hcrl
            adot1  => hcrl
            ltmp = .false.
            nla = nl
          endif
          j = streqv(s_sta%n,ltmp,nla,alpha1,adot1,i,0,plat,
     .      pos,s_sta%i,nitab)
        enddo
      endif
      deallocate(hcrl)

C --- Real-space screened structure constants ---
      i = nds**2*nkaps**2*nitab*nkapn
      allocate(s_sta%s(i)); call dpzero(s_sta%s,i)
      call ptr_str(s_str,8+1,'s',i,0,xx)
      if (lsdot) then
        allocate(s_stad%s(i)); call dpzero(s_stad%s,i)
      else
        allocate(s_stad%s(1))
      endif

C ... Make the strux
      if (amode == 1 .or. amode == 11 .or. amode == 13) then
        call strrs(nbas,npadl,npadr,alat,plat,pos,rwats,
     .    nl,kap2,nkaps,lmaxw,siinv,s_lat%cy,s_lat%cg,
     .    s_lat%indxcg,s_lat%jcg,lpgf,lsdot,nitab /= nttab,loka,
     .    s_sta%n,s_sta%i,s_sta%a,s_stad%a,s_sta%s,s_stad%s)
        if (amode == 13) then
c ... also make the value-laplacian unit basis out of G0 and G1/Hs/Hs-dot
          allocate(s_stag%s(4*nds**2*nttabg))
          call strg(ivl,nbas,nl,kmx,rsma,alat,plat,pos,
     .      s_ctrl%ips,rsmh,hcr,ehvl,
     .      s_lat%cy,s_lat%cg,s_lat%indxcg,s_lat%jcg,
     .      s_sta%n,s_sta%i,s_stag%ng,s_stag%s)
          deallocate(ehvl)
        endif
      else
        call strscr(loka,nbas,npadl,npadr,alat/avw,plat,pos,
     .    rwats,nl,kap2,nkapn,lmaxw,siinv,
     .    s_lat%cy,s_lat%cg,s_lat%indxcg,
     .    s_lat%jcg,lpbc == 11,lsdot,nitab /= nttab,s_sta%n,s_sta%i,
     .    s_sta%a,s_stad%a,tral,trad,s_sta%s,s_stad%s)
      endif
C ... Sort iax(7) by iax(1..5)
      call ppair5(1,nbasp,plat,pos,-1d0,s_sta%n,s_sta%i)
      if (amode == 13) then
        deallocate(rsmh,kmx,rsma)
      endif

C --- Write strux to disk ---
      ckbas = cksumf(pos,3*nbasp)
C ... for now
      lio = 1 + 100*(ncplx-1)
      if (amode == 1 .or. amode == 11) lio = lio+20000
      if (amode == 13) then
        lio = lio+20000
        liog = lio
      endif
      if (amode == 2) lio = lio+10000

C ... Copy data into s_str elements
      call ptr_str(s_str,4+1,'npr',nbasp+1,0,s_sta%n)
      call ptr_str(s_str,4+1,'iax',niax*nttab,0,s_sta%i)
      i = nds**2*nkaps**2*nitab*nkapn
      call ptr_str(s_str,4+1,'s',i,0,s_sta%s)
      if (.not. lsdot) i = 1
      call ptr_str(s_str,4+1,'sdot',i,0,s_stad%s)
      s_str%nttab = nttab
      s_str%nitab = nitab

C ... Write to disk
      ltmp = iostr(lio,'STR',nl,nbasp,max(nkaps,nkapn),kap2,itral,ckbas,
     .  lmaxw,nitab,s_sta)
      call fclose(fopno('STR'))
      if (lsdot) then
        ltmp = iostr(lio+800,'SDOT',nl,nbasp,max(nkaps,nkapn),kap2,
     .    itral,ckbas,lmaxw,nitab,s_stad)
      call fclose(fopno('SDOT'))
      else
        call dfclos(fopn('SDOT'))
      endif
      if (amode == 13) then
        call iosg(liog,'STRG',ivl,nl,nbas,ckbas,nttabg,s_stag)
        call fclose(fopno('STRG'))
        deallocate(s_stag%s,s_stag%ng)
      endif

      deallocate(s_sta%a,iax,s_sta%n,s_sta%s)
      s_sta%i => null()
      deallocate(tral,trad,rwats)

C --- Printout ---
C     if (lgors('str lshow,1',sstr)) then
      if (IAND(s_str%lshow,1) == 1) then
        call query('V>=30 to display structure constants',-1,0)
        lio = lio - 1 + 8
        i = nkaps
        if (amode == 2) i = nkapn
        if (amode == 13) liog = lio

        ltmp = iostr(lio,'STR',nl,nbasp,max(nkaps,nkapn),kap2,itral,
     .    ckbas,lmaxw,nitab,s_sta)
        call shostr(nds,nitab,nbasp,plat,pos,lio,s_sta%a,s_sta%i,
     .             s_sta%n,s_sta%s,nkaps,nkapn,kap2,1,1d0)
        if (lsdot) then
        call query('V>=40 to display sdot',-1,0)
        ltmp = iostr(lio+800,'SDOT',nl,nbasp,i,kap2,itral,ckbas,lmaxw,
     .             nitab,s_sta)
        call shostr(nds,nitab,nbasp,plat,pos,lio+800,s_sta%a,
     .    s_sta%i,s_sta%n,s_sta%s,nkaps,nkapn,kap2,1,1d0)
        endif
        if (amode == 13) then
          call iosg(liog,'STRG',ivl,nl,nbas,ckbas,nttabg,s_stag)
          call fclose(fopno('STRG'))
          call shosg(ivl,nl*nl,nttabg,nbas,plat,pos,s_stag,
     .      s_stag%n,s_stag%ng,s_stag%s,1,0.1d0)
        endif
      endif

C --- Plot screened orbital structure constants ---
      if (cmdopt('--plot',6,0,outs)) then
        call chkstr(s_ctrl,s_lat,s_spec,s_str)
      elseif (cmdopt('--pltg',6,0,outs)) then
        if (amode /= 13) then
          call info2(10,1,0,' ASASTR (warning) checking'//
     .      ' gaussian strux not generated with mode %i',amode,0)
        endif
        call chksg(s_ctrl,s_lat,s_spec,s_str)
      endif

      call rx0(prgnam)
      end

c...deb
C       subroutine xxn(s,nds)
C       integer nds,s(nds),i
C       print *, (s(i),i=1,nds)
C       end
C       subroutine xxr(s,nds)
C       integer nds,i
C       double precision s(nds)
C       print *, (s(i),i=1,nds)
C       end
c...deb

C      subroutine snot(nttab,ntab,iax,nbas)
C      implicit none
C      integer niax
C      parameter (niax=10)
C      integer nbas,nttab,ntab(nbas+1),iax(niax,nttab)
C
C
C      iax(:,2) = iax(:,5)
C      iax(:,3) = iax(:,6)
C      iax(:,4) = iax(:,10)
C      iax(6,2) = 4
C      iax(6,3) = 3
C      iax(6,4) = 2
C      ntab(2) = 2
C      ntab(3) = 4
C      nttab = 4
C
C      print *, 'modify nttab to', nttab
C      print *, 'iax is now'
C      print *, iax(1:5,1)
C      print *, iax(1:5,2)
C      print *, iax(1:5,3)
C      print *, iax(1:5,4)
C
C      end
