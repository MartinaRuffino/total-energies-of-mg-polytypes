      subroutine ephed(sopts,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
C- Electron phonon self-energy editor and generator
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: *
Co     Stored:    *
Cio    Passed to: *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat plat qlat npgrp nsgrp osymgr
Co     Stored:    *
Cio    Passed to: *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: *
Co     Stored:    *
Cio    Passed to: *
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read: nkabc
Co     Stored:    *
Cio    Passed to: *
Ci Inputs
Ci   sopts :command options performed automatically, before reading
Ci         :from standard input
Cr Remarks
Cr   ephed never returns.
Cl Local variables
Cl   agb   :holds spectral functions for multiple bands (dos branch)
Cl   lunits:0 eig,sigma,omega in Ry units
Cl         :1 eig,sigma,omega in eV units
Cl   limesh:0 => regular mesh is the original regular mesh containing SE
Cl         :1 => New regular mesh is to be constructed, interpolating sigma and maybe QP levels
Cl linputef:T when user supplies Ef
Cl  iocclst:list of intermediate bands to sum over when calculating e-ph sigma
Cl  nocclst:size of iocclst
Cl iocclstf:list of intermediate bands read from file (file analog of iocclst)
Cl nband_eph:size of iocclstf
Cl   iselst:list of bands for which to calculate diagonal e-ph self energy
Cl   nselst:size of iselst
Cl  iselstf:list of bands,  file analog of iselst
Cl nrenband_eph :size of iselstf
Cl  nphmode:Number of modes = natom * 3
Cl   haveg :0 Electron-phonon matrix elements file not accessed
Cl         :1 Electron-phonon matrix elements read or ready to read
Cl   irr   :Only used in post-processing of self-energy
Cl         :1  => qp match a regular mesh of points
Cl         :0  => qp do not match a regular mesh
Cl         :-1 => qp do not match a regular mesh; restrict calc to file q (no interpolation).
Cl         :      irr<0 specifies the same thing as limesh=3
Cl   ipqbnd:maps list of qp where spectra are calculated to qp supplied from input file.
Cl         :used when limesh is 2 or 3
Cl         :qp(iq) = qpse(ipqbnd(iq))
Cl   ef0   :Fermi level.  ef0=NULLI if not available.
Cl         :Note: ef0 is always in Ry units
Cl  chempot:chemical potential read from SE file.  Not used so far.  chempot=NULLI if not available.
Cl   nband :Number of one-particle states read from disk, or pared from disk read
Cl   nqg   :number of q-points on which electron-phonon g is given
Cl   qg    :q-points on which electron-phonon g is given
Cl   nqgi  :number of q-points on which electron-phonon g is interpolated to
Cl   qgi   :q-points on which electron-phonon g is interpolated to
Cl   nqgi
Cl   nkse  :number of points on which to calculate e-ph self-energy
Cl   qpse  :list of points on which to calculate e-ph self-energy
Cl   nqf   :number of qp in se file
Cl   qpg   :qp electron-phonon coupling matrix are available
Cl   qpf   :qp read from, written to se file
Cl   eigf  :file QP levels, read with self-energies
Cl         :eigf are relative to the Fermi level.  Units depend on lunits; see above
Cl   sexf  :file Fock exchange potential, if read
Cl   vxcf  :file LDA exchange-correlation potential, if read
Cl   nqgabc,nqgibz,nqgfbz,ipqg,wgibz:
Cl         :analogs of nkabc,nqibz,nqfbz,ipq,qibz,wibz
Cl         :for electron-phonon g interpolated to internal mesh to improve integration
Cl   lgetev:Calculate qp levels and use to improve interpolation
Cl   lqsgws:self-energy file has QSGW sigma
Cu Updates
Cu   04 Dec 18 First created.
C  ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_mix)::   s_mix
      type(str_move)::  s_move
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_optic):: s_optic
      type(str_gw)::    s_gw
      type(str_dmft)::  s_dmft
      type(str_strn) :: s_strn(*)
C ... Dynamically allocated arrays
      integer,pointer :: iocclstf(:),iocclst(:),iselstf(:),iselst(:)
      real(8),pointer:: qpf(:,:),qpse(:,:)
      real(8),pointer:: wibz(:),wgibz(:)
      real(8),pointer:: qfbz(:,:),qibz(:,:),qgibz(:,:)
      real(8),pointer:: qgfbz(:,:),eigg(:,:,:),sige(:,:,:)
      integer,pointer :: ipq(:,:,:),ipqg(:,:,:)
C     For Pade interpolation
      integer, allocatable :: iwk(:)
C ... Local parameters
      logical lnsave,lsopts,lse,ltet,linputef,lhf,llxc,lqsgws
      integer i,j,k,j1,j2,js1,js2,lunits,nbas,ndham,stdo,irr,limesh
      integer iq,iocc,ib,ix(10),nsp,nspc,nocclst,nselst
      integer nqgibz,nqgfbz,nkse
      integer nband,nbandx,lso,ifi
      integer npgrp,nsgrp,nqibz,isw,nqfbz,isp,ispx
      integer npoly,ioff
      integer haveg,nqf
      integer n1,n2,n3,nkabc(3),lshft(3)
      integer n1g,n2g,n3g,nqgabc(3)
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      equivalence (n1g,nqgabc(1)),(n2g,nqgabc(2)),(n3g,nqgabc(3))
      double precision qpi(3),TeV
      double precision qb(3,3),alat,vol,plat(3,3),qlat(3,3),kbT
      double precision rytoev,dymx,pi,ef0,evshft,sigdc
      integer NULLI,LW5
      integer nband_eph,nrenband_eph,nphmode,iqi,iqe,imode
      parameter (rytoev=13.60569193d0,npoly=6,dymx=1d-5,NULLI=-99999,LW5=5)
      real(8),parameter:: NULLR =-99999, a0=0.529177d0
      real(8):: eps=5d-3/rytoev
      complex(8),allocatable,target :: gmatrixf(:,:,:,:,:,:),gephi(:,:,:,:)
      complex(8),pointer :: gmatrix(:,:,:,:,:,:)
C     Parameters needed to call bndfp
      logical ltmp
      logical, parameter :: T=.true., F=.false.
      character dc*1, dc2*1, outs*256, tag*6, fn*120
C     Parameters for final state
C      integer ifinal,nfinal
C     double precision temp_1
C     double precision Epho, Phi_s, V_0
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))

      procedure(logical) :: latvec,isanrg,cmdopt
      procedure(integer) :: fopnx,fopna,fopng,fxst,nglob,a2vec,iprint,glitrp,wordsw,mkilsd
      procedure(real(8)) :: ddot,dglob,dlength

C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
c      double precision qk
c      integer jj1,jj2,jj3,k
c      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
c     .                    (jj2*ifac(2)-1)*qb(k,2) +
c     .                    (jj3*ifac(3)-1)*qb(k,3)

C     cyclic(i,j) = i-floor(dble(i)/j)*j  ! Cyclic modulus: returns n>=0 and n<j


      pi = 4*datan(1d0)
      kbT = NULLI
      nbas = nglob('nbas')
      nsp = nglob('nsp')
      stdo = nglob('stdo')
      haveg = 0  ! 0 electron-phonon matrix elements not read
                 ! i (i>0) electron-phonon matrix elements checked and ready to read k-point i
      lunits = 0 ! Assume Ry units
      lnsave = .false.          ! .true. when local SE generated
      linputef = .false.        ! .true. when user supplies Ef
      ef0 = NULLI               ! If not NULLI, ef0 is Fermi level got by running evsync
      sigdc = NULLI             ! If not NULLI, double counting to be added to sigma
      ioff = 0                  ! Offset to 1st band in QP levels corresponding to those in se
      nullify(qfbz,qibz,wibz,ipq,qpf,iocclstf,iocclst,iselstf,iselst)
      nocclst = NULLI; nselst  = NULLI
      nband = -1
      nspc = nglob('nspc')    ! 2 for noncollinear case
      lso = isw(IAND(s_ctrl%lncol,4) /= 0)
      lqsgws = .false.
      lhf = .false.; llxc = .false.

C ... Lattice information and symmetry
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      npgrp = s_lat%npgrp
      nsgrp = s_lat%nsgrp
      vol = s_lat%vol
      qpi = 0

C ... qp in full and irreducible BZ
      call info0(2,1,0,' EPHED: qp for irreducible and full BZ from gw nkabc')
      nkabc = s_gw%nkabc
      lshft = s_gw%lshft
      s_bz%nkabc = nkabc
      s_bz%lshft = lshft
      ltet = IAND(s_ctrl%lmet,3) == 3 .or. IAND(s_ctrl%ldos,4+2+1) /= 0
      nqfbz = n1*n2*n3
      if (nqfbz /= 0) then
        if (associated(qfbz)) deallocate(qfbz,wibz,ipq,qibz)
        allocate(qfbz(3,nqfbz),wibz(nqfbz),ipq(n1,n2,n3))
        wibz(1) = 0
        call bzmesh(plat,qb,n1,n2,n3,lshft,s_lat%symgr,1,ipq,
     .    qfbz,wibz,nqibz,nqfbz,0,0)
        deallocate(wibz)
        call mkqp(s_ctrl,s_bz,s_lat,ltet,F,1,-2)
        nqibz = s_bz%nkp
        allocate(qibz(3,nqibz),wibz(nqibz))
        call icopy(nqfbz,s_bz%ipq,1,ipq,1)
        call dcopy(nqibz,s_bz%wtkp,1,wibz,1)
        call dcopy(nqibz*3,s_bz%qp,1,qibz,1)
        irr = 1
        nqgfbz = nqfbz
        nqgibz = nqibz
        qgfbz => qfbz
        qgibz => qibz
        wgibz => wibz
        ipqg => ipq
      endif

      call pshpr(1)
      call supot(0,s_ctrl,s_lat,s_pot)
      call suham(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_str,s_strn)
      ndham = s_ham%ndham
      call poppr

      dc = sopts(1:1)
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the electron-phonon editor. ',
     .    'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the electron phonon editor.  ',
     .    'Enter ''?'' to see options.')
        lsopts = .false.
      endif

C ... Return here to resume parsing for arguments
   10 continue
      if (lsopts) then
        js2 = js2+1
        if (js2 > len(sopts)) then
          lsopts = .false.
          goto 10
        endif
        if (sopts(js2:js2) == dc) goto 10
        js1 = min(len(sopts),js2)
        call nwordg(sopts,0,dc,1,js1,js2)
        if (js2 < js1) lsopts = .false.
      endif

      write(stdo,"(/' Option : ')",advance='no')
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif

      call locase(outs)

C ... Parse and execute the next command
      if (.false.) then

C --- Null input ---
      elseif (outs == ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',' ''?'' to see menu')
        goto 10

C --- Read electron-phonon matrix elements from file ---
      elseif (outs(1:5)=='readg') then

        k = 5
        dc2 = outs(k+1:k+1)

C       Sanity check: nqibz must match contents of QIBZ, if file is present
        if (fopnx('QIBZ',2,-1,-1) == 1) then
          ifi = fopnx('QIBZ',2,1,-1)
          read(ifi,*) i
          call sanrg(.true.,i,nqibz,nqibz,'file:','qibz in file QIBZ')
          if (nqibz /= i) call rx('nqibz mismatch with file QIBZ')
          call fclose(ifi)
        endif

C       Parse optional switch fn
        fn = 'GMAT'
        i = wordsw(outs,dc2,'fn=','',j1)
        if (i /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          fn = outs(j1:j2)
        endif
        if (fopnx(fn,72,-1,-1) /= 1) then
          call info0(10,0,0,"missing file '"//trim(fn)//"' ... nothing done")
          goto 10
        endif
        ifi = fopng(fn,-1,4); rewind ifi

C       --------------- All of this data should be read from file, similar to ioseh
C       Input g has structure gmatrix(nband_eph,nrenband_eph,nphmode,nqibz,nqfbz,nsp))
        nband_eph = 12      ! Max number of bands in file for summing over intermediate states
        nrenband_eph = 12   ! Max number of bands in file at which to calculate diagonal sigma
        nphmode = 6         ! Number of modes = natom * 3

C       File list of intermediate states for pert theory.  Need add constraint: iselstf must be ordered.
        allocate(iocclstf(2*nband_eph)) !Use 2nd half of array as index when mapping iocclstf -> iocclst
        forall (ib=1:nband_eph) iocclstf(ib) = ib

C       File list of bands for which to calculated SE.  Need add constraint: iselstf must be ordered.
        allocate(iselstf(2*nrenband_eph)) !Use 2nd half of array as index when mapping iselstf -> iselst
        forall (ib=1:nrenband_eph) iselstf(ib) = ib

        allocate(gmatrixf(nband_eph,nrenband_eph,nphmode,nqibz,nqfbz,nsp))
        read(ifi) gmatrixf(1:nband_eph,1:nrenband_eph,1:nphmode,1:nqibz,1:nqfbz,1:nsp)
C       --------------- end of temporary file read

C   ... Parse optional occ list ("occ" means list of intermediate QP be summed over to make sigma)
        if (associated(iocclst) .and. .not. associated(iocclst,iocclstf)) deallocate(iocclst)
        iocclst => iocclstf  ! Default
        if (wordsw(outs,dc2,'occ=','',j1) > 0) then
          j2 = scan(outs(j1:),dc//dc2//' ') + j1-2
          nocclst = mkilsd(outs(j1:j2),-1,ix)
          if (nocclst <= 0) goto 99
          allocate(iocclst(nocclst))  ! Last half used temporarilty as a pointer?
          if (mkilsd(outs(j1:j2),nocclst,iocclst) < 0) call rx('sugws bug in mkilsd')
          if (.not. associated(iocclstf)) call rx('iocclstf not allocated')
          i = 0
          do  ib = 1, nocclst
            call hunti(iocclstf,nband_eph,iocclst(ib),0,i)
            if (iocclstf(i+1) /= iocclst(ib)) then
              call info2(0,0,0,' oops ... elements in occ= must belong list in g',1,2)
              goto 98
            endif
            iocclstf(ib+nband_eph) = i+1
          enddo
        endif

C   ... Parse optional ib list (list of QP at which diagonal sigma is to be calculated)
        if (associated(iselst) .and. .not. associated(iselst,iselstf)) deallocate(iselst)
        iselst => iselstf  ! Default
        if (wordsw(outs,dc2,'ib=','',j1) > 0) then
          j2 = scan(outs(j1:),dc//dc2//' ') + j1-2
          nselst = mkilsd(outs(j1:j2),-1,ix)
          if (nselst <= 0) goto 99
          allocate(iselst(nselst))
          if (mkilsd(outs(j1:j2),nselst,iselst) < 0) call rx('sugws bug in mkilsd')
          if (.not. associated(iselstf)) call rx('iselstf not allocated')
          i = 0
          do  ib = 1, nselst
            call hunti(iselstf,nrenband_eph,iselst(ib),0,i)
            if (iselstf(i+1) /= iselst(ib)) then
              call info2(0,0,0,' oops ... elements in ib= must belong list in g',1,2)
              goto 98
            endif
            iselstf(ib+nrenband_eph) = i+1
          enddo
        endif

        if (associated(iselst,iselstf) .and. associated(iocclst,iocclstf)) then
          gmatrix => gmatrixf
        else
          allocate(gmatrix(nocclst,nselst,nphmode,nqibz,nqfbz,nsp))
          forall (iocc=1:nocclst,ib=1:nselst,imode=1:nphmode,iqe=1:nqfbz,iqi=1:nqibz,isp=1:nsp)
     .      gmatrix(iocc,ib,imode,iqi,iqe,isp) =
     .      gmatrixf(iocclstf(nband_eph+iocc),iselstf(nrenband_eph+ib),imode,iqi,iqe,isp)
          deallocate(iocclstf,iselstf,gmatrixf)

        endif

        call info5(10,0,0,
     .    ' g matrix from file '//trim(fn)//
     .    ' to make SE for %i bands, %i qp   pert from %i bands, nq=%i   nsp=%i',nselst,nqfbz,nocclst,nqibz,nsp)
        call info2(40,0,0,' sum check %2:1;6g',sum(gmatrix),2)

        haveg = 1


C --- Set smearing ---
      elseif (outs(1:4) == 'eps ' .or. outs(1:4) == 'eps=') then

        if (outs(1:4) == 'eps ') then
          call word(outs,2,j1,j2)
        else
          call word(outs,1,j1,j2)
          j1 = j1+4
        endif

        if (j2 < j1) goto 99
        j = 0; j = a2vec(outs(j1:j2),len(outs(j1:j2)),j,4,' ',1,1,1,ix,eps)
        if (j <= 0) goto 99
        call info2(10,0,0,'%1fSmearing set to %g %?#n#eV#Ry#',eps,lunits)
        goto 10

C --- Set Fermi ---
      elseif (outs(1:3) == 'ef ' .or. outs(1:3) == 'ef=') then

        if (outs(1:3) == 'ef ') then
          call word(outs,2,j1,j2)
        else
          call word(outs,1,j1,j2)
          j1 = j1+3
        endif

        if (j2 <= j1) goto 99
        j = 0
        j = a2vec(outs(j1:j2),len(outs(j1:j2)),j,4,' ',1,1,1,ix,ef0)
        if (j <= 0) then
          ef0 = NULLI
          goto 99
        endif
        call info2(10,0,0,"%1fFermi level set to %g Ry",ef0,0)
        linputef = .true.
        goto 10

C --- Set units ---
      elseif (outs(1:6) == 'units' .or. outs(1:6) == 'units=') then

        if (outs(1:6) == 'units') then
          call word(outs,2,j1,j2)
        else
          call word(outs,1,j1,j2)
          j1 = j1+6
        endif

        if (j2 <= j1) then
          goto 99
        elseif (outs(j1:j2) == 'ry' .or. outs(j1:j2) == 'Ry') then
          lunits = 0
          call info0(10,0,0,' Units set to Ry')
        elseif (outs(j1:j2) == 'ev' .or. outs(j1:j2) == 'eV') then
          lunits = 1
          call info0(10,0,0,' Units set to eV')
        else
          goto 99
        endif

C --- Make self-energy, spectral functions to arbitrary q, or make DOS ---
      elseif (outs(1:6) == 'makese' .or. .false.) then

        dc2 = outs(6+1:6+1)
        lse  = outs(1:6) == 'makese'   ! T if to construct SE
        TeV = .0259d0  ! Default to room temperature
        lse = .true.
        if (haveg == 0) goto 95

C   ... Construct list of k-points where to make self-energy
        if (lse .and. wordsw(outs,dc2,'allk',dc//dc2//' ',j1) > 0) then
          nkse = nqibz
          qpse => qibz
        elseif (wordsw(outs,dc2,'ik=','',j1) > 0) then
          j = j1-1
          if (a2vec(outs,len(outs),j,2,dc//dc2//' ',3,1,1,ix,iq) /= 1) goto 99
          if (iq > nqf) goto 98
          allocate(qpse(3,1))
          ltmp = isanrg(iq,1,nqf,' ephed: ','iq',.true.)
          qpi(1:3) = qpf(1:3,iq) ! keep since qpse may be reallocated
          qpse(1:3,1) = qpf(1:3,iq)
        elseif (wordsw(outs,dc2,'q=','',j1) > 0) then
          j = j1-1
          allocate(qpse(3,1))
          if (a2vec(outs,len(outs),j,4,', '//dc//dc2,4,2,3,ix,qpse) /= 3) goto 99
        else
          goto 99
        endif

C   ... Optional switches
        isp = 1
        nqgabc = 0              ! flag that no interpolating mesh params nqgabc read

        if (wordsw(outs,dc2,'kt=','',j1)>0) then
          j = j1-1
          if (a2vec(outs,len(outs),j,4,' '//dc//dc2,3,1,1,ix,kbT) /= 1) goto 99
        endif

C        if (wordsw(outs,dc2,'domg=','',j1) > 0) then
C          print *, 'SE not ready for domg= option ...'; goto 99
C        endif

        if (wordsw(outs,dc2,'isp=','',j1) > 0) then
          if (nsp < 2) then
            call info0(1,0,0,' ignoring isp= ...')
          else
            j = j1-1
            if (a2vec(outs,len(outs),j,2,' '//dc//dc2,3,1,1,ix,isp) /= 1) goto 99
            if (isanrg(isp,1,nsp,' ephed: ','isp',.false.)) goto 99
            if (isp < 1 .or. isp > 2 .or. isp == 2 .and. nspc == 2) then
              call info0(1,0,0,'isp=2 but noncollinear ...')
              goto 98
            endif
          endif
        endif

        if (lse .and. wordsw(outs,dc2,'nq=','',j1) > 0) then
          j = j1-1
          j = a2vec(outs,len(outs),j,2,', '//dc//dc2,4,2,3,ix,nqgabc)
          if (j < 1) goto 99
          call fill3in(j,nqgabc)
        endif

C   ... This ends parsing of optional switches.  Sanity checks
        if (irr == 0) then
          call info0(10,0,0,"%1fNo connection established with BZ."//
     .      "%N Modify GW_NKABC and restart program.")
          goto 10
        endif

C   ... Check if require new interpolating mesh
        if (nqgabc(1) == 0) then ! points on the original regular mesh
          limesh = 0
          nqgabc = nkabc
        else ! A new regular mesh is to be constructed, interpolating sigma
          limesh = 1
        endif
C       End of inputs for making electron-phonon sigma

C   ... Make list of interpolating qp and corresponding evals.  Need
C       (1) points 1:nqgibz interpolated q mesh
C       (2) evals for all these points
C       This block makes:
C       nband  = highest band for which evals are needed
C       nqgibz = number of points in the interpolated mesh of the irreducible BZ
C       qgibz  = corresponding qp
C       nqgfbz = number of points in the interpolated mesh of the full BZ
C       qgfbz  = corresponding qp
C       wgibz  = weights for irr zone, interpolated mesh.
C       ipqg   = map nqgibz points in irr zone to the full zone

        call info0(10,0,0,' ')

        if (limesh == 0) stop 'handle limesh=0'

C       Make qgibz, nqgibz = irr qp and number on new regular mesh, and ipqg
        call pshpr(1)
        i = nqgabc(1)*nqgabc(2)*nqgabc(3)
        allocate(qgibz(3,i),wgibz(i),ipqg(nqgabc(1),nqgabc(2),nqgabc(3)))
        call bzmshp('sugws',3,nqgabc,lshft,plat,s_lat%symgr,nsgrp,
     .    nsgrp,.false.,0,0,qb,nqgibz,qgibz,wgibz,ipqg,ix)
        deallocate(wgibz)
        call poppr

C       Make qgfbz = (qp on new mesh in full BZ)
        nqgfbz = n1g*n2g*n3g
        allocate(qgfbz(3,nqgfbz),wgibz(nqgfbz),iwk(nqgfbz))
        wgibz(1) = 0
        call pshpr(1)
        call bzmesh(plat,qb,n1g,n2g,n3g,lshft,s_lat%symgr,1,iwk,qgfbz(1,1),wgibz,i,n1g*n2g*n3g,0,0)
        call poppr
        deallocate(wgibz,iwk)
        if (i /= n1g*n2g*n3g) call rx('bug in bzmesh')

        if (limesh == 0) then
          call info2(10,0,0,' Use as-given %sx%3i k-mesh',nqgabc,2)
        else if (limesh == 1) then
          call info2(10,0,0,' Interpolate bands and g to %sx%3i qp mesh  %i irr pts ',nqgabc,nqgibz)
        else
          call rx('bug in ephed')
        endif
        deallocate(ipqg)        ! It will be remade, this time with nqgibz points
        nband = max(maxval(iocclst),maxval(iselst))
        nbandx = nband*nspc

C       QP levels for the nqgibz points
        allocate(eigg(nbandx,nqgibz,nsp),sige(nbandx,nqgibz,nsp),qgibz(3,nqgibz),wgibz(nqgibz))
        allocate(ipqg(nqgfbz,1,1))
        i = 71                  ! Make qgibz(:,:), wgibz(:), ipqg(:,1,1))
        i = i+1000              ! Make eigg also
        if (linputef) i = i+2   ! use given ef0
        i = i+4000              ! Make spin weights if noncollinear
        call evalqp(i,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_move,s_pot,
     .    s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nbandx,ioff,nsp,nqgibz,
     .    nqgabc,nqgibz,0,lunits,ef0,eigg,qgibz,wgibz,ipqg,evshft)

C   --- Compute self-energy for each k in nkse ---
        ispx = isp; if (nspc == 2) ispx = 1
        allocate(gephi(nocclst,nselst,nphmode,nqgfbz))
        do  k = 1, nkse

C     ... Interpolate g to qgfbz
C     ... Interpolate g to each point in qgfbz
          do  iq = 1, nqgfbz

            call info5(10,1,0,' SUGWS:  interpolate g, '//
     .        '%i band(s), %i qp%?#(n==2)#, spin %i#%j#',
     .        nocclst,nqgfbz,nsp,ispx,5)

          if (glitrp(2,qpi,plat,nqfbz,qfbz,ipq,nqibz,nocclst,nselst,nphmode,
     .      gmatrix(1,1,1,1,k,isp),gephi(1,1,1,iq)) < 0) then
            call rx('failed to interpolate g')
          endif

        enddo  ! loop over q
        enddo

        stop

        goto 10

C ... show
      elseif (outs(1:5) == 'show ') then

C --- Abort ---
      elseif (outs(1:2) == 'a ') then
        call rx0('aborting sigfun editor ... no file written')

C --- quit ---
      elseif (outs(1:2) == 'q '. or. outs(1:5) == 'quit ') then
        if (lnsave) then
          print '('' se file not saved ... really quit?'')'
          read(*,'(a150)') outs
          call locase(outs)
          if (.not. (outs(1:1) == 'y' .or. outs(1:1) == 'q'))
     .      goto 10
        endif
        call rx0('exit sigfun editor')

C --- Help ---
      elseif (outs == '?') then
        print 310
        print 312
C       print 315
C       print 316
  310   format(
     .    ' Select one of these options:'/
     .    t4,'readg[Options]'/t12,'read electron-phonon matrix elements from file'/
     .    t12,"File name is 'GMAT' unless specified."/
     .    t12,"Optionals are separated by a delimiter, which is the first character, e.g. '@'"/
     .    t14,"@fn=nam",t25,"use nam for file name"/
     .    t14,"@occ=lst",t25,"include contribution to sigma from QP states lst"/
     .    t14,"@ib=lst",t25,"include contribution to sigma from QP states in list"/
     .    t14,"@useef",t25,"specify that the file chemical potential be the Fermi level"/
     .    t14,"@ib=list",t25,"Make e-ph sigma for QP states in this lst"//
     .
     .    t4,'units Ry|eV'/t12,
     .    'select units data for which data is to be returned.  Default is Ry.'//
     .    t4,'eps #',t12,'select smearing for sampling integration.'//
     .    t4,'ef #', t12,'Use # for Fermi level, overriding internally calculated value.')

  312   format(/t4,
     .    'makese ik=#|k=#1,#2,#3|allk|band[~args] [nq=#1[,#2,#3]] [isp=#]'/t12,
     .    'Make electron-phonon sigma for given k and band.  Writes to ??'/t12
     .    'ik|k  single k : index in QIBZ, or k in Cartesian coordinates ','(units of 2pi/alat)'/t12,
     .    'allk  sigma is made for all k in the irreducible BZ and written to disk'/t12,
     .    'band  sigma is made for k points along symmetry lines and written to disk'/t12,
     .    '      Optional ~args are parsed like options of the --band switch'/t12,
     .    '      In the multiple-k (latter) cases, data is written to disk as it is generated'/t12,
     .    'Options:'/t12
     .    'nq    interpolates g to new mesh for q-integration'/t12,
     .    'isp   functions for given spin (default=1)'
     .    )

  315   format(/
     .    t4,'qsgwh',t12,'Construct quasiparticle hamiltonian from lmf evals')

  316   format(/
     .    t4,'merges12'/t12,'merge spin1 and spin2 self-energies.'/
     .    t12,'Requires seq1.ext and seq2.ext reside on disk.'//
     .    t4,'savesea [fn]'/
     .    t4,'savese  [fn]'/t12,'saves q-interp self-energy, ',
     .    'ASCII or binary format.'/t12,
     .    "File name is 'seia.ext' or 'seib.ext' unless specified."//

     .    t4,'q',t12,'to quit the editor'/
     .    t4,'a',t12,'to abort')
      else
        call word(outs,1,j1,j2)
        call info0(0,0,0,' ephed:  "'//trim(outs(j1:j2))//
     .    '" not recognized ... enter "?" for a list of options')
      endif
      goto 10

   95 call info0(0,0,0,'%10pelectron-phonon file must be opened before invoking this command ... nothing done')
      goto 10
   96 call info0(0,0,0,' ephed:  missing tag '//trim(tag)//' in '//trim(outs)//' ... nothing done')
      goto 10
   97 call info0(0,0,0,'%10pself-energy must be read before invoking this command ')
   98 call info0(0,0,0,' ephed:  improper usage of '//trim(outs)//' ... nothing done')
      goto 10
   99 call info0(0,0,0,' ephed:  failed to parse arg starting at: '//trim(outs(j1:))//' ... nothing done')
      goto 10

      end

      integer function glitrp(mode,qpi,plat,nqfbz,qfbz,ipq,nq,
     .  nocclst,nselst,nphmode,gephq,gephi)
C- Linearly interpolate self-energy to specified q-point from mesh points
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :bit 1  Interpolate eigi (not implemented!)
Ci         :bit 2  Interpolate gi
Ci   qpi   :interpolate gephq to this q-point
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nqfbz :number of k-points in full BZ
Ci   qfbz  :k-points in full BZ
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci         :mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci         :ipq(1)=0 => not used
Ci   nq    :number of k-points in the irreducible BZ.
Ci         :NB: if ipq(1) is zero, nq should be nqfbz!
Ci  nocclst:number of intermediate bands, leading dimension of gephq and gephi
Ci  nselst :number of bands where sigma is calculated 2nd dimension of gephq and gephi
Cl  nphmode:Number of phonon modes
Ci  gephq :e-ph for each element in nocclst, nselst, nphmode, qp
Co Outputs
Co   gephi :gephq interpolated to qpi
Co   glitrp: 0 if interpolation is normal
Co         :-1 if interpolation failed
Cl Local variables
Cr Remarks
Cu Updates
Cu   04 Dec 18  Adapted from seqitp
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nq,nqfbz,mode,nphmode,ipq(nqfbz)
      integer nocclst,nselst
      real(8), intent(in)  :: qpi(3),qfbz(3,nqfbz),plat(3,3)
      real(8), intent(in)  :: gephq(nocclst,nselst,nphmode,nq)
      real(8), intent(out) :: gephi(nocclst,nselst,nphmode)
C ... Local parameters
      integer i,j,k,iqfbz(4,2),mode1,ib,iocc,imode
      double precision wt(4),xx,qp(3),qlat(3,3)
      integer, parameter :: NULLI=-99999
      procedure(integer) :: iprint,qplin
      procedure(real(8)) :: ddot

      mode1 = mod(mode/2,2)
      if (mode1 == 0) return

      call pshpr(1)
C     print *, '!!'; qfbz(1:3,1) = (/.4d0,.5d0,.6d0/)
      glitrp = qplin(qpi,plat,nqfbz,qfbz,ipq,iqfbz,wt)
C     If no mapping to IBZ, substitute full BZ for irr
      if (ipq(1) == 0) then
        call icopy(4,iqfbz(1,1),1,iqfbz(1,2),1)
        if (nq /= nqfbz)
     .    call rx('glitrp: inconsistent use of mapping to irr mesh')
      endif
      call poppr

      call info5(41,0,0,' Interp g for q=%3:1,5;5d, wt=%4:1,5;5d',qpi,wt,0,0,0)
      if (iprint() >= 45) then
        call mkqlat(plat,qlat,xx)
        do  j = 1, 4
          k = iqfbz(j,1)
          i = iqfbz(j,2)
          call shorbz(qfbz(:,k)-qpi(:),qp,qlat,plat)
          xx = dsqrt(ddot(3,qp,1,qp,1))
          call info8(40,0,0,' qi(%i) =%3;11,6D  |qi-qpi|=%,6;6d  iq=%i irr=%i wt=%,6;6d',
     .      j,qfbz(1,k),xx,iqfbz(j,1),i,wt(j),0,0)
        enddo
      endif
      if (glitrp < 0) then
        call info2(20,0,0,' wt gen failed for qp=%3:1,5;5d, wt =%4:1;9F',qpi,wt)
        return
      endif

      do  iocc = 1, nocclst
        do  ib = 1, nselst
          do  imode = 1, nphmode
            gephi(iocc,ib,imode) = 0
            do  j = 1, 4        ! Sum over interpolating qp
              gephi(iocc,ib,imode) = gephi(iocc,ib,imode) + wt(j)*gephq(iocc,ib,imode,iqfbz(j,2))
            enddo
          enddo
        enddo
      enddo

      end
