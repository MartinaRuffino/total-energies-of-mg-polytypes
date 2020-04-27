      subroutine supgbc(s_ctrl,s_lat,s_ham,s_pot,s_str,s_spec,s_site,
     .  vRLshf,vshft,zp,nz,nzp,qp,nq,igf,kcplx,pgplp)
C- Set up boundary conditions for layer Green's function
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nl nbas nbasp nspin npl lpgf lham lrel lasa ldlm
Ci                 nccomp nclass nspec lncol
Co     Stored:     lham lasa
Co     Allocated:  *
Cio    Elts passed:lasa
Cio    Passed to:  supghz suhamz mkpotf makpfz pgsif
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos
Cio    Passed to:  pgsif
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lncol ldham hord nlibu udiag neula lgen3 lham offH
Co     Stored:     lncol
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  supghz suhamz mkpotf makpfz pgsif plhamnc plham
Cio                plhamso pgbevl
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz ves
Co     Stored:     *
Co     Allocated:  cp pf dpf ddpf dddpf pfr dpfr ddpfr papg gmar palp
Co                 gma
Cio    Elts passed: pf gma palp dpf ddpf dddpf pfr dpfr ddpfr papg gmar
Cio                pp sop socscl pprel dlmwt pfnc
Cio    Passed to:  supghz suhamz mkpotf makpfz mkptfp mksopf gvbma
Cio                pgsif plhamnc plham pgbevl
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  pgsif plhamnc plham pgbevl
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr z rmt lmxa hcr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  supghz suhamz mkpotf makpfz pgsif
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb class dlmcl pnu spec clabel v0
Co     Stored:     *
Co     Allocated:  pfr dpfr ddpfr
Cio    Elts passed:pfr dpfr ddpfr
Cio    Passed to:  supghz suhamz mkpotf makpfz mkptfp mkfrpf mksopf
Cio                gvbma pgsif plhamnc plham plhamso pgbevl
Ci Inputs
Ci   vRLshf:array of orbital- and spin-dependent potential shifts
Ci   vshft :array of potential shifts by site
Ci   zp,nz :list of complex energies, and number of points
Ci         :including points (if any) along non-equilibrium contour
Ci   nzp   :number of points on equilibrium contour
Ci   qp,nq :list of k-points, and number
Ci   igf   :kind of boundary condition to make.  The allowed values are
Ci          1 semi-infinite GF gLL and gRR
Ci          3 S_L+1,L g_LL S_L,L+1 and S_R-1,R g_RR S_R,R-1
Ci   kcplx :specifies form of complex rep'sn of gf generated
Ci          0: real, imaginary separated: gf = gf(ldrg,ldrg,1..2)
Ci          1: gf in double complex format
Ci          2: real, imaginary in columns : gf = gf(ldrg/2,1..2,ldrg/2)
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Co Outputs
Co   supgbc attempts to read GF for Left and Right PL from disk for each
Co          energy, spin, qp from files 'lbc' and 'rbc' (and for
Co          spin-pol case 'lbc2' and 'rbc2').  If supgbc fails to read
Co          the GF, they are made and stored in those files.
Cl Local variables
Cl  modsif :argument 'mode' passed to pgsif
Cl    lspc :0 if spins are to be treated as uncoupled
Cl         :1 if spins are to be treated as coupled
Cr Remarks
Cu Updates
Cu   13 Jan 18 Iterative refinement of surface gf
Cu   30 Apr 16 Surface GF can be noncollinear.
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   20 Nov 07 Additions for LDA+U implementation
Cu   12 Jan 04 (S.Faleev) Argument list changed for non-equilibrium mode
Cu   22 Mar 03 vshft changed for revised electrostatics
Cu   10 Jan 02 Adapted to new generation of potential functions
Cu   21 Dec 01 Altered argument list
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer igf,nz,nzp,nq,pgplp(6,-1:*)
      double precision zp(2,nz),qp(3,nq),vRLshf(*),vshft(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: g(:)
C ... Local parameters
      logical lrefineonly
      integer lncol,lspc
C     Bits for matching conditions (10-1000s digit of lio)
      integer, parameter :: MNC=1,MNB=2,MLB=8,MPL=16,MZP=64,MSC=128
C     Bits for selecting which data to I/O, iogfrs
      integer, parameter  :: ISH=1,IPOS=16,IPLMOD=64
      integer ifi(2,2),ipl,ipr,iq,isp,iz,k,kcplx,ld0,
     .  lgii,lrd,modsif,nbas,nl,npl,nsp,nbasp
C     integer mpipid,procid
      procedure(integer) isw,fopnx,iogfrs,frename
      double precision plat(3,3),xv(1),maxdel,max1it
      character*8 fnam(2), gfstrn(4)*8, src*16, dest*16
      data gfstrn / '  GF ', ' nothing', ' S*G*S+', ' S*G'/

C ... Setup
      call tcn('supgbc')
      call getpr(ipr)
      lgii = -1
      plat = s_lat%plat
      nl = s_ctrl%nl
      nbas = s_ctrl%nbas
      nbasp = s_ctrl%nbasp
      nsp = s_ctrl%nspin
      npl = s_ctrl%npl
      lncol = s_ham%lncol
      lspc = isw(s_ctrl%lpgf(2) >= 10 .and. lncol /= 0 .and. nsp == 2)
      if (nq < 1) call rx('supgbc requires at least 1 qp')

C     procid = mpipid(1)  ! Just a stub ... pgfasa not set up for MPI yet

C --- for Left, Right PL, do ---
      do  ipl = -1, npl, npl+1

      lrefineonly = s_ctrl%refinegstol(1) < 0

C --- Check boundary condition connecting to SIF GF, layer ipl ---
      modsif = 100*kcplx
      if (igf == 3) modsif = modsif+4
      if (ipl == -1) then
        fnam(1) = 'lbc'
        fnam(2) = 'lbc2'
C       vshf = vshfl
      else
        fnam(1) = 'rbc'
        fnam(2) = 'rbc2'
C       vshf = vshfr
        modsif = modsif+1
      endif

      ld0 = pgplp(4,ipl)
      do  isp = 1, nsp
        ifi(isp,1) = fopnx(fnam(isp),100,16+8+4+0,-1)
        ifi(isp,2) = ifi(isp,1)
        if (lspc == 1) exit
      enddo
      if (lspc == 0) then
        allocate(g(ld0*ld0))
      else
        allocate(g(ld0*2*ld0*2))
      endif

C ... Read header; return lrd=0 if read is successful
      call info0(20,0,-1,' SUPGBC: read '//trim(fnam(1))//' for'//trim(gfstrn(igf))//' ... ')
      call pshpr(1)
      lrd = 0
      do  isp = 1, nsp
      lrd = min(iogfrs(10000*IPLMOD+10*(MNC+MNB+MLB+MPL+MZP+MSC*lspc)+2,0,30,fnam,ifi(isp,1),
     .    1,1,nbas,k,zp,qp,plat,s_lat%pos,pgplp(1,ipl),xv,xv,0,0,kcplx,g),lrd)
      if (lrd == -1 .or. lspc == 1) exit  ! Skip reading header spin 2
      enddo
      call poppr
      call mpibc2(lrd,1,2,2,.false.,'','')  ! If any procid cannot read, set lrd<0 for all

C ... Read GF for each energy, spin, qp until error, or all are read
      call pshpr(ipr-10)
      if (lrd >= 0) then
        do  iz  = 1, nz
        do  isp = 1, nsp
          do  iq = 1, nq
C          zmv(1) = zp(1,iz)-vshf
C          zmv(2) = zp(2,iz)
            lrd = iogfrs(10000*(ISH+IPLMOD)+10*MZP+6,igf,30,fnam,ifi(isp,1),1,1,nbas,k,
     .        zp(1,iz),qp(1,iq),plat,s_lat%pos,pgplp(1,ipl),xv,xv,0,0,kcplx,g)
            if (lrd < 0) goto 20
          enddo
          if (lspc == 1) exit
        enddo
        enddo
      endif
   20 continue
      call mpibc2(lrd,1,2,2,.false.,'','')  ! If any procid cannot read, remake bc for all
      call poppr
      call info5(20,0,0,
     .  '%?;(n>=0);ok%j; missing -> generate for ipl=%i;'//
     .  '%?~n~ (collinear mode)~~'//
     .  '%?~n~ (noncollinear mode)~~',
     .  lrd,ipl,isw(lspc == 0.and.lncol /= 0),lspc,0)
      if (lrd < 0 .and. lrefineonly) then
        lrefineonly = .false.   ! Cannot refine if no gs available
        call info0(20,0,0,'%9fsuppress refine-only mode')
        s_ctrl%refinegstol(1) = abs(s_ctrl%refinegstol(1))
      elseif (lrefineonly) then
        call info0(20,0,0,'%9fRefine gs with iterative embedding')
      endif

C --- Create and store SIF GF, for each energy, spin, qp, layer ipl ---
      if (lrd < 0 .or. lrefineonly) then

C ... Assign output file to temporary name
      if (lrefineonly) then
        do  isp = 1, nsp
          ifi(isp,2) = fopnx(trim(fnam(isp))//'x',100,16+8+4+0,-1)
        enddo
      endif

C ... Write header
      call pshpr(1)
      do  isp = 1, nsp
        if (lrefineonly) then ! Read header of original file
        rewind ifi(isp,1)
        lrd = min(iogfrs(10000*IPLMOD+10*(MNC+MLB+MNB+MPL+MZP+MSC*lspc)+2,0,30,fnam,
     .    ifi(isp,1),1,1,nbas,k,zp,qp,plat,s_lat%pos,pgplp(1,ipl),xv,xv,0,0,kcplx,g),lrd)
        if (lrd < 0) call rx('failed to read header in file '//fnam)
        endif
        rewind ifi(isp,2)
        lrd = iogfrs(10000*(IPOS+IPLMOD)+10*(MNC+MLB+MNB+MPL+MZP+MSC*lspc)+3,0,0,fnam,
     .    ifi(isp,2),1,1,nbas,k,zp,qp,plat,s_lat%pos,pgplp(1,ipl),xv,xv,0,0,kcplx,g)
        if (lspc == 1) exit
      enddo
      call poppr
C     call mpibc2(lrd,1,2,2,.false.,'','')  ! If any procid cannot write
C     call snot('supgbc 235 start iz loop')

C ... Make GF and write for each energy, spin, qp
C     Read 1st gf with : mc -f9f12.6 -r:br:s=3:open lbc.pg
      maxdel = 0; max1it = 0
      do  iz  = 1, nz
      do  isp = 1, nsp

C   ... Hamiltonian setup for this energy and spin
C       Trick suham into making collinear p.f. make nspc=ldham(4)->1
        if (lspc == 0) then
          lncol = s_ham%lncol
          s_ham%lncol = 0
        endif
        call supghz(s_ctrl,s_site,s_ham,s_pot,s_spec,vRLshf,vshft,isp,zp(1,iz),isw(iz > nzp))
C       Restore proper lncol
        s_ham%lncol = lncol
        call pshpr(ipr-10)
        do  iq = 1, nq
          if (lrd < 0) goto 30
          s_ctrl%refinegstol(2) = 0
          if (lrefineonly) then ! Read existing g from file
            s_ctrl%refinegstol(1) = abs(s_ctrl%refinegstol(1))
            if (mod(modsif/10,10) /= 0) call rx('supgbc: illegal mode')
            modsif = modsif+20
            lrd = iogfrs(10000*(ISH+IPLMOD)+10*(MZP+MSC*lspc)+6,igf,30,fnam,ifi(isp,1),1,1,nbas,k,
     .        zp(1,iz),qp(1,iq),plat,s_lat%pos,pgplp(1,ipl),xv,xv,0,0,kcplx,g)
            if (lrd < 0) call rx('failed to read contents of file '//fnam)
C           k = ld0*s_ham%ldham(4)
C           call yprmi('file gs for iq=%i iz=%i',iq,iz,2+kcplx,g,k**2,k,k,k)
          endif
C          call info5(10,0,0,' supgbc: make gs for iz=%i isp=%i iq=%i q = %s,(%3;6d)',
C     .      iz,isp,iq,qp(1:3,iq),5)
          call pgsif(modsif,s_ctrl,s_lat,s_ham,s_pot,s_str,s_spec,s_site,
     .      ipl,pgplp,qp(1,iq),isp,ld0,ld0,lgii,g)
          if (lrefineonly) then
            s_ctrl%refinegstol(1) = -s_ctrl%refinegstol(1)
            modsif = modsif-20
          endif
C         zmv(1) = zp(1,iz)-vshf; zmv(2) = zp(2,iz)
C         Write to disk
C          k = ld0*s_ham%ldham(4)
C          call yprmi('write gs for iq=%i iz=%i',iq,iz,2+kcplx,g,k**2,k,k,k)
          lrd = iogfrs(10000*(ISH+IPLMOD)+10*(MSC*lspc)+7,igf,30,fnam,ifi(isp,2),1,1,nbas,k,zp(1,iz),
     .          qp(1,iq),plat,s_lat%pos,pgplp(1,ipl),xv,xv,0,0,kcplx,g)
          maxdel = max(maxdel,s_ctrl%refinegstol(2))
          max1it = max(max1it,s_ctrl%refinegstol(3))
C          print 444, ' procid,iz,isp,iq maxdel q',procid,iz,isp,iq,maxdel,qp(:,iq)
C  444     format(a,1x,4i4,1pg12.4,3x,0p,3f12.6)
        enddo                   ! iq
        call poppr
        if (lspc == 1) exit
      enddo                     ! isp
   30 continue
      enddo                     ! iz

C     call snot('supgbc 289 exit iz loop')

      call mpibc2(max1it,1,4,1,.false.,'','')
      call mpibc2(maxdel,1,4,1,.false.,'','')
      if (maxdel > 0) call info2(20,0,0,
     .  ' max change in gs from iterative refinement: 1st it %;4g  max %;4g ',max1it,maxdel)

      call clhamz(s_pot%pf)
C     call clpgv(s_pot)
      if (lrd < 0) call rxs('SUPGBC: failed to write file ',fnam)
      call fclose(ifi(1,2))  ! Close the file to flush the buffer
      if (nsp == 2) call fclose(ifi(2,2))

C     Rename temporary ?bcx and ?bc2x files to ?bc and ?bc2
      if (ifi(1,1) /= ifi(1,2)) then
        do  isp = 1, nsp
          call fclose(ifi(isp,1)) ! Close file where gs is read
          iq = fopnx(src,192,0,ifi(isp,2)) ! Full file name
          iq = fopnx(dest,192,0,ifi(isp,1)) ! Full file name
          iq = frename(trim(src),trim(dest))
          if (iq /= 0) call rx('supgbc: posix rename command failed')
          call fclr(' ',ifi(isp,2))
          if (lspc == 1) exit
        enddo
      endif

      endif

C ... Position file pointer past header
      call pshpr(1)
      do  isp = 1, nsp
        ifi(isp,1) = fopnx(fnam(isp),100,16+8+4+0,-1)
        lrd = iogfrs(10000*IPLMOD+10*(MNC+MNB+MLB+MPL+MZP+MSC*lspc)+2,0,30,fnam,ifi(isp,1),1,1,
     .    nbas,k,zp,qp,plat,s_lat%pos,pgplp(1,ipl),xv,xv,0,0,kcplx,g)
        if (lspc == 1) exit
      enddo
      call poppr

      deallocate(g)
      enddo ! ipl, ipr

C     In future iterations use disk data only
      s_ctrl%refinegstol(1) = abs(s_ctrl%refinegstol(1))

      call info0(20,0,0,'')

      call tcx('supgbc')
      end
C      subroutine snot(strn)
C      character*(*) strn
C
C      integer mpipid,procid
C      procid = mpipid(1)  ! Just a stub ... pgfasa not set up for MPI yet
C
C      print 333, procid,trim(strn)
C  333 format(1x,i3,1x,a)
C
C      end
