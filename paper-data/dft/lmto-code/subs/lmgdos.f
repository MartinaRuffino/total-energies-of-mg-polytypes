      subroutine lmgdos(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,s_pot,
     .  s_str,s_optic,vorb,dmatu,fmode0,nvfit,nfcof,ivcof,dosprm,wg,
     .  nchan,nspc,eband,nbmax,qnu,rhos,gdos)
C- Gradient of partial DOS wrt selected potential parameters
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspec nspin lncol lsx lham lasa lgen3
Ci                 loptc nbasp nclasp ipc
Co     Stored:     lham lasa
Co     Allocated:  *
Cio    Elts passed: lscr lasa ldos lham ipc dclabl ics lstonr ips rmax
Cio                lgen3 lqp ncomp idcc nrc
Cio    Passed to:  bndasa nmpot secmat secmtn secmt2 pp2alp asaddq
Cio                asaopm optint optin2 optinq intdrv dellmp
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw alat plat nsgrp nkd nkq awald vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:vol pos dlv qlv cg jcg indxcg
Cio    Passed to:  bndasa secmat secmtn sstrxq asaddq asaopm
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  idu lmxb lmxa a nr z rmt hcr name
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bndasa getidu sumlst nmpot secmat secmtn mullmf
Cio                iorbtm
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class v0
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bndasa getidu sumlst nmpot lmfitmr6 mullmf
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp ntet lmet n range w nevmx efmax lio ndos
Ci                 dosw zval fsmom qp
Co     Stored:     numq qp ef
Co     Allocated:  qp wtkb swtk
Cio    Elts passed:qp wtkp swtk wtkb idtet ipq n w def
Cio    Passed to:  bndasa subzi optint optin2 optinq intdrv gtpdss
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham qss neula eterms nbf nmto nlibu lmaxu udiag
Ci                 lsig kmto ndhrs
Co     Stored:     lham eterms
Co     Allocated:  *
Cio    Elts passed:iprmb eula bdots magf nprs iaxs hrs
Cio    Passed to:  bndasa nmpot secmat secmtn sstrxq asaddq
C
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  nrhos vmtz ves
Co     Stored:     *
Co     Allocated:  ppn
Cio    Elts passed:pp qpp sop bxc ppn pti
Cio    Passed to:  bndasa nmpot secmat secmtn secmt2 asaddq asaopm
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:alph nkaps iax kaps s sdot nitab adot
Cio    Passed to:  bndasa secmat secmtn secmt2 pp2alp dellmp
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ltet ocrng unrng dw window
Co     Stored:     ocrng unrng
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bndasa asaopm optint optin2 optinq intdrv
Ci Inputs
Ci   vorb  :orbital dependent potential matrices
Ci   dmatu :density matrix for LDA+U
Ci   nvfit :number of (independent) parameters to vary; see Remarks
Ci   dosprm: 1 = lower bound of ref dos
Ci         : 2 = upper bound of ref dos
Ci         : 3 = number of reference DOS points (integer)
Ci         : 4 = lower bound of fit dos
Ci         : 5 = upper bound of fit dos
Ci         : 6 = number of fit DOS points (integer)
Ci   wg    :DOS is broadened by a gaussian of width wg
Ci   nchan :number of DOS channels
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   eband :energy bands; alias eb (sec*.f)
Ci   nbmax :maximum number of bands
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   rhos  :spin density-matrix
Co Outputs
Co   gdos  :gradient of DOS
Cl Local variables
Cl         :
Cr Remarks
Cr   There are ncof coefficients all told.  They are one of three types:
Cr   - coefficient is free to vary
Cr     nvfit is the number of such coefficients.
Cr   - coefficient is frozen
Cr   - change in coefficient is linked to change in another coefficient.
Cu Updates
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   17 Nov 10 L-M fitting passes through fmode0
Cu   16 Apr 10 Allow gaussian broadening of DOS
Cu   25 Mar 10 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer fmode0,nbmax,nvfit,nfcof,ivcof(nfcof),nchan,nspc
      double precision gdos(nvfit,*)
      double precision wg,qnu(*),rhos(2,*),dosprm(6),eband(nbmax,*)
      double complex vorb(*),dmatu(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_optic):: s_optic
C ... Dynamically allocated arrays
      real(8),allocatable:: aamom(:),dosp(:,:),dosm(:,:)
C ... Local parameters
      integer nevmx,nbas,ivar,nfilem,fopna,i1,nsp,nspx,ldham(16),
     .  nl,nlink,mode2
      double precision sumev,amag(3),efermi,eps,xx
      integer,allocatable:: ivfit(:,:)
      procedure(integer) :: nglob
C ... MPI
      integer mpipid,procid,master

      procid = mpipid(1)
      mode2 = 100*fmode0
      master = 0
      nbas = nglob('nbas')
      nsp = nglob('nsp')
      nl = nglob('nl')
      nspx = nsp / nspc
      allocate(aamom(nbas))
      call lmfitmr2(-1,1,1,1,1,1,1,eps)

      call bndasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,s_pot,s_str,
     .  s_optic,vorb,dmatu,efermi,-1,0,0,0,xx,xx,eband,nbmax,nevmx,qnu,
     .  sumev,rhos,amag,aamom)

      ldham = s_ham%ldham
C     nl = s_ctrl%nl

      i1 = int(dosprm(6))
      allocate(dosp(i1,nchan*nspx))
      allocate(dosm(i1,nchan*nspx))
      allocate(ivfit(nfcof,2))

      do  ivar = 1, nvfit

        call lmfitmr6(10,s_site,nl,nbas,nsp,s_ham%iprmb,
     .    ldham,ivar,nfcof,ivcof,nlink,ivfit)

C   ... DOS at pp + eps
        call dellmp(1+mode2,s_ctrl,s_str,nlink,ivfit,nfcof,ivcof,
     .    s_pot%pp)
        call bndasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,s_pot,s_str,
     .    s_optic,vorb,dmatu,efermi,ivar,0,0,0,xx,xx,eband,nbmax,nevmx,
     .    qnu,sumev,rhos,amag,aamom)
        if (procid == master) nfilem = fopna('MOMS',-1,4)
        i1 = int(dosprm(6))
        call gtpdss(1,nfilem,s_bz,.false.,nchan,nsp,
     .    i1,dosprm(4),dosprm(5),wg,dosp)
        call fclose(nfilem)
C        call iodos(3,-fopna('dosp',-1,0),dosp,i1,nspx,i1,nchan,
C     .    dosprm(4),dosprm(5),nspx,0d0,0)
C        call fclose(fopna('dosp',-1,0))

C   ... DOS at pp - eps
        call dellmp(4+mode2,s_ctrl,s_str,nlink,ivfit,nfcof,ivcof,
     .    s_pot%pp)
        call bndasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,s_pot,s_str,
     .    s_optic,vorb,dmatu,efermi,ivar,0,0,0,xx,xx,eband,nbmax,nevmx,
     .    qnu,sumev,rhos,amag,aamom)
        if (procid == master) nfilem = fopna('MOMS',-1,4)
        i1 = int(dosprm(6))
        call gtpdss(1,nfilem,s_bz,.false.,nchan,nsp,
     .    i1,dosprm(4),dosprm(5),wg,dosm)
        call fclose(nfilem)
C        call iodos(3,-fopna('dosm',-1,0),dosm,i1,nspx,i1,nchan,
C     .    dosprm(4),dosprm(5),nspx,0d0,0)
C        call fclose(fopna('dosm',-1,0))

C   ... Restore pp to original value
        call dellmp(1+mode2,s_ctrl,s_str,nlink,ivfit,nfcof,ivcof,
     .    s_pot%pp)

C   ... Numerical derivative; transfer to gdos
        call daxpy(i1*nchan*nspx,-1d0,dosm,1,dosp,1)
        call dscal(i1*nchan*nspx,1/(2*eps),dosp,1)
        call dcopy(i1*nchan*nspx,dosp,1,gdos(ivar,1),nvfit)

      enddo

C     Debugging printout
C      call info0(30,0,0,' LMASA: saving gdos, file gdos')
C      call iodos(3,-fopna('gdos',-1,0),dosp,i1,nspx,i1,nchan,
C     .  dosprm(4),dosprm(5),nspx,0d0,0)
C      call fclose(fopna('gdos',-1,0))
C      call rx0('done')

      deallocate(dosp,dosm,aamom,ivfit)

      end

      subroutine dellmp(mode,s_ctrl,s_str,nlink,ivfit,nfcof,ivcof,pp)
C- Increment potential parameter for numerical differentiation
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :1 increment parameter ivar by epsilon
Ci         :2 decrement parameter ivar by epsilon
Ci         :3 increment parameter ivar by 2*epsilon
Ci         :4 decrement parameter ivar by epsilon
Ci         :10s digit
Ci         :0 2nd gen ASA pp's are in gamma repsn
Ci         :1 2nd gen ASA pp's are in alpha repsn
Ci         :100s digit
Ci         :1 ASA C and Delta
Ci         :2 ASA C+enu and Delta
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read:
Ci     Stored:
Ci     Passed to: pp2alp
Ci   nlink :total number of coefficients linked to given
Ci         :coefficients --- number of values to change
Ci   ivfit :vector of information about the ivar'th coefficient that is
Ci         :varied, including all coefficients linked to teh ivar'th one.
Ci         :nlink elements are returned
Ci         :ivfit(1:nlink,1) flags type of coefficient to vary
Ci         :ivfit(1:nlink,2) flags which coefficient within a type
Ci         :So far, only ASA variations C and Delta are implemented (mode=10)
Ci         :In that case: pp depends on l,isp, and class
Ci         :ivfit specifies which pp as follows:
Ci         : ivfit(:,1)   ivfit(:,2)        type
Ci         :    2         l+10*isp+100*ic   C parameter
Ci         :    3         l+10*isp+100*ic   Delta
Cio Inputs/Outputs
Cio  pp    :potential parameters in alpha or gamma repsn, depending
Cio        :on 10s digit mode (atomsr.f).  One parameter is incremented
Cio        :depending on ivfit.
Cio        :Manner of increment depends on mode.
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   23 Mar 10
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nlink,nfcof,ivcof(nfcof),ivfit(nfcof,2)
      double precision pp(6,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_str)::   s_str
C ... Local parameters
      integer nl,nsp,mode0,mode1,mode2
      procedure(integer) :: nglob
      double precision xx

      if (mod(mode,10) == 0) return
      nl = nglob('nl')
      nsp = nglob('nsp')
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)

C --- Backup copy of parameters ---
C      if (mode == -1) then
C      endif

C --- Increment potential parameter, depending on 1s digit mode ---
      if (mode0 >= 1 .and. mode0 <= 4) then

C       pp shifts must be in orthogonal representation
        if (ivfit(1,1) <= 6 .and. mode1 == 1) then
          call pp2alp(1,s_ctrl,s_str,0,pp)
C         allocate(o(nlspc))
C         call pptrns(1,nl,w(oipc),nclass,nsp,dum,nbas,pp,o)
C         deallocate(o)
        endif
C       call prmx('before pp',pp,6,6,nl*nsp*1)
        call lmfitmr2(2*mode0+10*mode2,nlink,xx,nfcof,nl,nsp,ivfit,pp)
C       call prmx('after pp',pp,6,6,nl*nsp*1)
C       Restore pp's to tb repsn
        if (ivfit(1,1) <= 6 .and. mode1 == 1) then
          call pp2alp(0,s_ctrl,s_str,0,pp)
C         allocate(o(nlspc))
C         call pptrns(0,nl,w(oipc),nclass,nsp,w(oalph),nbas,pp,o)
C         deallocate(o)
        endif
      endif

      end
