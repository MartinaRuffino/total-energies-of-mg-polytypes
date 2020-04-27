      subroutine mmag(s_ctrl,s_ham,s_lat,s_move,s_spec,s_strn)
C- Micromagnetics simulations
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nbasp nspec sdxsi nvario lpgf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  qss
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:eula
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat avw nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos istab symgr
Cio    Passed to:  *
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  nmodt modt ct prmint tstot kt ts tsequ gyro
Co     Stored:     ts
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mmdyn mmstp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rham z name
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  suemph
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  strn
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:strn
Cio    Passed to:  str_pack
Ci Inputs
Cl   nmodt :number of thermostat modes
Cl   modt  :modt(1..nmodt) : thermostat modes
Cl   ct    :ct(1..nmodt)   : thermostat mode coefficients
Cl   neq   :total number of integration variables
Cl   intgmd:integration mode
Cl         :1 Bulirsch-Stoer integration
Cl   amag  :average magnetization
Co Outputs
Cs Command-line switches
Cs   --cont   : continuation of MD run
Cs   --wcvmj  : Write r.s. exchange interactions to disk, formatted for CVM program
Cs   --wrsj   : Write r.s. exchange interactions to disk, for use by other programs
Cs   -tstart= :
Cr Remarks
Cr   mmag integrates the Landau-Lifshitz equations of motion
Cr   with global daemons.  It operates in one of two modes:
Cr
Cr 1.  An empirical (heisenberg) hamiltonian H is supplied.
Cr
Cr     It is supplied EITHER from the contents of a file, OR in the
Cr     form of a set of empirical rules.  In the former case,
Cr     structure element 'strn mmham' contains a string
Cr     FILE=filename; the contents of H are supplied by file
Cr     filename.  The syntax of that file are desribed in (and data
Cr     are read by) subroutine iorsj.  Alternatively, structure
Cr     element 'strn mmham' contains a set of rules that define the
Cr     Heisenberg coefficients.  In this latter case, H is assembled
Cr     by suemph.f, which also documents the syntax for the rules.
Cr
Cr     In this mode mmag calls mmdyn, which computes the forces
Cr     internally and integrates the equations of motion over the
Cr     specified interval (see struc smove).
Cr
Cr 2.  The magnetic forces are supplied by the calling program.
Cr
Cr     In this mode mmag calls mmdyn for a single time, supplying
Cr     the Eular angles and forces for that time.  mmdyn returns
Cr     with a new time for which it requires updated Eular angles
Cr     and forces.
Cr
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   06 Sep 11 Started migration to f90 structures
Cu   02 Oct 03 Implement --wcvmj to write CVM file
Cu   03 Jun 03 Writes to simpler sv file, rather than save file
Cu   24 Nov 02 Micromagnetics hamiltonian may be read from a file
Cu   17 Aug 01 updated from ASA package version 5.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mxcsiz
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_move)::  s_move
      type(str_spec)::  s_spec(*)
      type(str_strn) :: s_strn(*)
C ... Dynamically allocated local arrays
C#ifndef LINUXF
      integer, allocatable :: iax(:)
C#elseC
C      integer, pointer :: iax(:) => null()
C#endif
      integer, allocatable :: istb2(:)
      integer, allocatable :: ipair(:)
      integer, allocatable :: isite(:)
      real(8), allocatable :: epair(:)
      real(8), allocatable :: aamom(:)
      real(8), allocatable :: esite(:)
      real(8), allocatable :: rtab(:)
      real(8), allocatable :: hpair(:)
      real(8), allocatable :: huns(:)
      real(8), allocatable :: wk(:)
      real(8), allocatable :: intwk(:)
      real(8), allocatable :: frc(:)
C ... Local parameters
      character outs*256,jfilnm*72,hamrul*(2048)
      logical cmdopt,a2bin,parstr,rdstrn
      character htype(2)*5,dc*1
      double precision alat,plat(3,3),avw,xsi(3),qspirl(4),ct(3),
     .  amag(3),etot,tnow,xx,ttot,sclrsj,tolrsj
      integer fopna,i,ifi,jfi,j,j1,j2,k,lpbc,lpgf,ltype(2),
     .  modt(3),nbas,nbasp,neq,nmodt,nspec,nstep,nttab,ntype(2),nvarms,
     .  iorsj,iprint,parg,nfilm,nfilj,istyle,str_pack
      integer,allocatable :: ntab(:)
      integer niax
      parameter (niax=10)
C     Symmetry operations
      integer nsgrp
C     Specific to integration mode
      integer intgmd,mx
      real(8), allocatable :: rham(:)

      data htype /'h1:','h2:'/

      nbas = s_ctrl%nbas
      nbasp = s_ctrl%nbasp
      nspec = s_ctrl%nspec
      xsi = s_ctrl%sdxsi(1:3)
      nvarms = s_ctrl%nvario
      alat = s_lat%alat
      plat = s_lat%plat
      avw = s_lat%avw
      allocate(rham(nspec))
C     call spackv(10,'spec rham',sspec,1,nspec,rham)
      do  i = 1, nspec
        rham(i) = s_spec(i)%rham
      enddo

C     call spec2class(s_spec,nspec,-1,'rham',1,xx,rham)
      qspirl = s_ham%qss
      lpgf = s_ctrl%lpgf(1)
      lpbc = 0
      if (lpgf > 0) lpbc = 11

C --- Make or read the Hamiltonian ---
      ntype(1) = 1
      ntype(2) = 1
      ltype(1) = 0
      ltype(2) = -1
      j = str_pack('mmham',-2,s_strn,hamrul)
      allocate(aamom(nbasp))
      call dvset(aamom,1,nbasp,1d0)
      allocate(esite(nbasp))
      allocate(isite(nbasp))

C --- Read R.S. hamiltonian from file ---
      i = 0
      if (parstr(hamrul,'FILE=',len_trim(hamrul),5,'=',i,j)) then
        call word(hamrul(1+j:),1,j1,j2)
        jfilnm = hamrul(j+j1:j+j2)
        ifi = fopna(jfilnm,-1,0)
        call info0(20,1,0,' ... Read r.s. hamiltonian from file '//
     .    jfilnm(1:j2-j1+1))

C   ... Setup call for file read
        nfilj = iorsj(0,ifi,nbasp,alat,plat,aamom,xx,xx,xx,xx,xx)
        call rxx(nfilj < 0,'MMAG: missing or incompatible file'//
     .    ' lattice vectors')

C   ... Make a neighbor table that encompasses all file entries
        if (lpbc == 0) then
          i = 3
          j = -1
        elseif (lpbc == 1 .or. lpbc == 11) then
          i = 2
          j = 1
        else
          call rx('MMHAM: not implemented for lpbc>1')
        endif
        call info2(35,0,0,' ... Make neighbor table: rcut = %d * alat',
     .    xx,0)
        call pshpr(iprint()-10)
        mxcsiz = 0
        allocate(ntab(nbasp+1))
        call pairs(nbas,nbasp,alat,plat,[alat*xx/2*(1+1d-6)],s_lat%pos,
     .    [-1],i,j,[0],nttab,ntab,iax,mxcsiz)

C       For now make a copy
        allocate(rtab(3*nttab))
        call mkrtab(0,alat,plat,s_lat%pos,iax,nttab,s_lat%pos,rtab)
        call poppr

        allocate(epair(nttab)); call dpzero(epair,nttab)
        allocate(ipair(nttab)); call iinit(ipair,nttab)

C   ... Read Heisenberg coefficients from file
C       NB: note file's rtab does not scale by alat
        call dscal(3*nttab,1/alat,rtab,1)
        allocate(hpair(nttab)); call dpzero(hpair,nttab)
        nfilm = iorsj(3,ifi,nbasp,alat,plat,aamom,iax,rtab,
     .    ntab,hpair,xx)
        if (nfilm <= 0) call rx('mmag: no hamiltonian read')

        if (nfilj /= nfilm) call info2(20,0,0,'%8f(warning) only'//
     .    ' matched %i out of %i file entries',nfilm,nfilj)

C   ... Purge iax table of unused entries
        if (nfilm /= nttab) then
          call info0(20,0,0,' ... Reduce iax table and'//
     .      ' re-read pair hamiltonian')

C         Purge iax table of unused entries
          call symiax(1,plat,nbas,xx,xx,xx,0,ntab,iax,nttab,i)
          if (nttab /= nfilm) then
            call fexit2(-1,111,' Exit -1 : mmag: bug in symiax:'//
     .        ' nttab=%i, nfilm=%i',nttab,nfilm)
          endif
C         Remake rtab
          call mkrtab(000,1d0,plat,s_lat%pos,iax,nttab,s_lat%pos,rtab)

C         Re-read pair hamiltonian with new iax table
          nfilm = iorsj(2,ifi,nbasp,alat,plat,aamom,iax,
     .      rtab,ntab,hpair,xx)
          if (nttab /= nfilm) call rx('bug in iorsj')

C         Close input file rsj
          call fclose(ifi)

          nsgrp = s_lat%nsgrp
          allocate(istb2(nsgrp*nbas))
          call istbpm(s_lat%istab,nbas,nsgrp,0,istb2)
          allocate(huns(nttab)); call dpzero(huns,nttab)
          call dcopy(nttab,hpair,1,huns,1)
          i = 0 + 10*2*0 + 1000
C          call ppair1(40,1,nbas,nbas,s_ctrl%ips,alat,plat,s_lat%pos,xx,
C     .      nttab,ntab,iax,mxcsiz)
          call rsmsym(i,plat,1,xx,0,nbas,s_lat%pos,1,1,1,nttab,ntab,
     .      iax,s_lat%symgr,istb2,nsgrp,0,1,huns,hpair)
          deallocate(huns,istb2)

        endif

        call dscal(3*nttab,alat,rtab,1)
C       Distribute into epair and esite
        call pmmag1(nbas,ntab,rtab,hpair,esite,epair)
        deallocate(hpair)

C --- Use empirical rules to make hamiltonian ---
      else

C   ... Make iax and rtab
        if (lpbc == 0) then
          i = 3
          j = -1
        elseif (lpbc == 1 .or. lpbc == 11) then
          i = 2
          j = 1
        else
          call rx('MMHAM: not implemented for lpbc>1')
        endif

        call pshpr(iprint()-10)
        allocate(ntab(nbas+1))
        mxcsiz = 0
        call pairs(nbas,nbasp,alat,plat,rham,s_lat%pos,
     .    s_ctrl%ips,i,j,s_ctrl%pgfsl,nttab,ntab,iax,mxcsiz)

        allocate(rtab(3*nttab))
        call mkrtab(0,alat,plat,s_lat%pos,iax,nttab,s_lat%pos,rtab)
        call poppr

        allocate(epair(nttab)); call dpzero(epair,nttab)
        allocate(ipair(nttab)); call iinit(ipair,nttab)

        call suemph(s_spec,hamrul,ntab,iax,rtab,nbasp,
     .    nttab,s_ctrl%ips,ntype,ltype,htype,esite,isite,
     .    epair,ipair)
      endif

C --- Write R.S. Jij to file.  Switch syntax: --wrsj[:fn=name][:scl=#] ---
      if (cmdopt('--wrsj',6,0,outs)) then
        i = 6
        sclrsj = 1
        tolrsj = 1d-6
        jfilnm = 'rsj2'
        dc = outs(i+1:i+1)
        if (dc /= ' ') then
C     ... Return here to resume parsing for arguments
          j2 = i
   50     continue
          j2 = j2+1
          if (outs(j2:j2) == dc) goto 50
          j1 = min(len(outs),j2)
          call nwordg(outs,0,dc//' ',1,j1,j2)
          if (j2 >= j1) then
            if (.false.) then
            elseif (outs(j1:j1+3) == 'scl=')  then
              j = 0
              i = parg('scl=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,sclrsj)
              if (i <= 0) goto 59
            elseif (outs(j1:j1+3) == 'tol=')  then
              j = 0
              i = parg('tol=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,tolrsj)
              if (i <= 0) goto 59
            elseif (outs(j1:j1+2) == 'fn=')  then
              jfilnm = outs(j1+3:j2)
            else
              goto 59
            endif
            goto 50
   59       call rxs('mmag: failed to parse --wrsj options:',outs)
          endif
        endif

        call info0(20,0,0,' ... writing r.s. hamiltonian to file '//
     .    jfilnm)
        ifi = fopna(jfilnm,-1,0)
        call dscal(3*nttab,1/alat,rtab,1)
        i = iorsj(1,-ifi,nbasp,alat,plat,aamom,iax,rtab,ntab,epair,xx)
        call dscal(3*nttab,alat,rtab,1)
        call fclose(ifi)
        call rx0('mmag after writing file '//jfilnm)

C --- Write Jij to file suitable to CVM input ---
      elseif (cmdopt('--wcvmj',7,0,outs)) then
        i = 7
        istyle = 6
        sclrsj = 1
        tolrsj = 1d-6
        jfilnm = 'cvmj'
        dc = outs(i+1:i+1)
        if (dc /= ' ') then
C     ... Return here to resume parsing for arguments
          j2 = i
  150     continue
          j2 = j2+1
          if (outs(j2:j2) == dc) goto 150
          j1 = min(len(outs),j2)
          call nwordg(outs,0,dc//' ',1,j1,j2)
          if (j2 >= j1) then
            if (.false.) then
            elseif (outs(j1:j1+4) == 'mode2')  then
              istyle = 7
            elseif (outs(j1:j1+3) == 'scl=')  then
              j = 0
              i = parg('scl=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,sclrsj)
              if (i <= 0) goto 159
            elseif (outs(j1:j1+3) == 'tol=')  then
              j = 0
              i = parg('tol=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,tolrsj)
              if (i <= 0) goto 159
            elseif (outs(j1:j1+2) == 'fn=')  then
              jfilnm = outs(j1+3:j2)
            else
              goto 159
            endif
            goto 150
  159       call rxs('mmag: failed to parse --wcvmj options: ',outs)
          endif
        endif

        call info0(20,0,0,' ... writing r.s. hamiltonian to file '//
     .    jfilnm)
        ifi = fopna(jfilnm,-1,0)
        call dscal(3*nttab,1/alat,rtab,1)
        call dscal(nttab,sclrsj,epair,1)
        i = iorsj(istyle,-ifi,nbasp,alat,plat,aamom,iax,
     .      rtab,ntab,epair,xx)
        call dscal(3*nttab,alat,rtab,1)
        call fclose(ifi)
        call rx0('mmag after writing file '//jfilnm)

      endif

C ... Printout of initial conditions
      call pmmag2(nbas,aamom,ntab,iax,epair)
      allocate(wk(3*nbas+3)); call dpzero(wk,3*nbas+3)
      nmodt = s_move%nmodt
      modt = s_move%modt
      ct = s_move%ct

      neq = nbas*2 + nmodt
      intgmd = s_move%prmint(1)
      if (intgmd == 1) then
        mx = s_move%prmint(5)
        allocate(intwk(neq*(4+mx)+5))
      else
        call rxi('mmag: unknown integration method',intgmd)
      endif
      nstep = -1
      allocate(frc(nbas*3+3)); call dpzero(frc,nbas*3+3)

C --- Switches for restarting ---
      ifi = fopna('SV',-1,0)
      call poseof(ifi)
      tnow = 0
      j = 8
      if (cmdopt('-tstart=',j,0,outs))
     .  call rxx(.not.a2bin(outs,tnow,4,0,' ',j,len(outs)),
     .  'mmag: failed to parse -tstart=')
      j = 4

C     Switch --cont tells mmag to do the following:
C     1. Set tnow to last entry in save file and increment ttot
C     2. Read euler angles from eula-sv file
      if (cmdopt('--cont',j,0,outs)) then
        call rx('mmag not ready for new save')
        backspace ifi
        if (rdstrn(ifi,outs,len(outs),.false.)) then
        i = 0
        if (parstr(outs,'time=',len(outs)-5,5,'=',i,j)) then
        if (parg('time=',4,outs,i,len(outs),' ',1,1,k,tnow) == 1) then
          ttot = s_move%tstot
          ttot = ttot + tnow
          call rx('pack s_move tstot')
C         s_move%tstot = ttot
          call info2(20,1,0,
     .      ' ... read Euler angles from file eula-sv;'//
     .      ' set start time to %d',tnow,0)
          goto 30
        endif
        endif
        endif
        call info0(20,0,0,
     .    ' ... read Euler angles from file eula-sv;'//
     .    ' unable to read start time from save file')
   30   continue

        jfi = fopna('eula-sv',-1,0)
        rewind jfi
        j = 0
        call ioeula(nbas,3,s_ham%eula,j,xsi,jfi)
        call fclose(jfi)
        call rxx(j /= 1,'file eula-sv has neul>1')

      endif

C --- Micromagnetics dynamics ---
C  62 continue
      call mmdyn(s_move,nbas,xsi,s_ham%eula,qspirl,wk,epair,
     .  ipair,nstep,iax,nttab,intwk,frc,amag,aamom,
     .  etot,tnow,nvarms,ifi)

C --- Calculate the pairwise forces and continue ---
C      if (nstep > 0) then
C        etot = 0
C        call dpzero(frc,nbas*3+3)
C        call rx('mmag: need update mmpair')
C        call mmpair(nbas,nttab,iax,w(oebil),w(oibil),aamom,wk,
C     .   etot,amag,frc)
C        goto 62
C      endif

      end
      subroutine pmmag1(nbas,ntab,rtab,hpair,esite,ebil)
C- Distribute hpair into one-body and bilinear terms
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   rtab  :table of connecting vectors
Ci   hpair :pair hamiltonian.  On-site terms are those for which
Ci         :connecting vectors are zero
Co Outputs
Ci   esite :on-site terms
Ci   ebil  :bilinear terms
Cl Local variables
Cl         :
Cr Remarks
Cr   Assumes first entry in iax table is on-site term, if it exists
Cu Updates
Cu   12 Apr 03 Allow there to be no on-site term
Cu   24 Nov 02 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,ntab(nbas+1)
      double precision hpair(*),esite(nbas),ebil(*),rtab(3,*)
C ... Local parameters
      integer ib,it

      do  ib = 1, nbas
        it = ntab(ib)+1
        ebil(it) = hpair(it)
        if (rtab(1,it)**2+rtab(2,it)**2+rtab(3,it)**2 == 0) ebil(it)=0
        esite(ib) = hpair(ntab(ib)+1)
        do  it = ntab(ib)+2, ntab(ib+1)
          ebil(it) = hpair(it)
        enddo
      enddo

      end

      subroutine pmmag2(nbas,amom,ntab,iax,hpair)
C- Printout of energy for collinear spins, also decp'n by site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   amom  :vector of moments
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   hpair :pair hamiltonian.  On-site terms, if present, are
Ci         :assumed to have zero amplitude.
Ci         :connecting vectors are zero
Cr Remarks
Cr   Assumes first entry in iax table is on-site term
Cu Updates
Cu   24 Nov 02 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,ntab(nbas+1),niax
      parameter (niax=10)
      integer iax(niax,*)
      double precision hpair(*),amom(nbas)
C ... Local parameters
      integer ib,jb,it,nttab,iprint
      double precision dsum,ei,kboltz,facJ0,e0,amag0
C     kboltz is in units eV K^-1
      parameter (kboltz=8.617d-5)

      nttab = ntab(nbas+1)

C      x = dsum(nttab,hpair,1) / 2
C      call info2(20,1,0,' Collinear energy 1/2*sum_ij J_ij = %,8;8g',x,
C     .  0)

      e0 = 0
      amag0 = dsum(nbas,amom,1)
      do  it = 1, nttab
        ib = iax(1,it)
        jb = iax(2,it)
        e0 = e0 + dsign(1d0,amom(ib)*amom(jb))*hpair(it)/2
C        print 333, it, e0
C  333   format(i5,f17.7)
      enddo
      call info2(20,1,0,
     .  ' Collinear energy 1/2*sum_ij sgn(s_i s_j) J_ij = %,8;8g'//
     .  '  sum s_i = %;6d',e0,amag0)

      if (iprint() < 50) return

      facJ0 = 2d0/3d0 * 13.6d0 / kboltz
      call info0(0,0,0,
     .  '   site      J_00(mRy)     J_0(mRy)    2/3 J_0 (K)')

      do  ib = 1, nbas
C        e0 = hpair(ntab(ib)+1)
        e0 = 0
        ei = 0
        do  it = ntab(ib)+1, ntab(ib+1)
          ei = ei + hpair(it)
        enddo
        e0 = -0.5d0*e0
        ei = -0.5d0*ei
        call info5(0,0,0,' %,6i %;12,3D %;12,3D %;12,1D ',
     .    ib,e0*1000,ei*1000,facJ0*ei,0)
      enddo

      end
