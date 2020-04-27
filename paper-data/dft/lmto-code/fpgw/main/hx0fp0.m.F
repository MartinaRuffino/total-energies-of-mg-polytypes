      program hx0fp0
C-  Calculate x0, W-V, \epsilon
C-------------------------------------------------
Cr Calculate chi0 by the following steps.
Cr  tetwt5: weights for tetrahedron integration
Cr  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
Cr  dpsion5: calculate real part by the Hilbert transformation from the Im part
Cr
Cr Be careful for the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
Cr      The indexing starts with core first and then valence on top of core
Cr      So n-index in "in" for valence electron is different from "inv".
Cr      Search in and inv.
C --- subroutine variables. but not perfect yet. In future, we need to make it automatic
C    i input
C    o output
C    m input and output
C
Cl Local variables
Cl  nmb     : total number of mixed product basis functions I,J in which <I|v|J> is calculated
Cl          : (formerly ngb)
Cl  nmbas   : number of mixed product basis functions used for susceptibility
Cl          : It may be nmb, but may be smaller in some contexts
Cl          : For magnetic susceptibility without local fields, nmbas = # magnetic atoms
Cl          : For charge susceptibility without local fields, nmbas = 1
Cl          : (formerly nmbas_in, nmbas, nmbas1, nmbas2)
Cl  nmbas1  : In this code, same as nmbas.  But it may be different in other contexts,
Cl          : so it is kept distinct for clarity
Cr comment rule
Cr C --- Chapter
Cr C ... section
Cr ! note for what it is
Cr c Comment out a line (this can be memo for development).

Cr eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
Cr (only epsPP_lmf_chipm mode works).
C-------------------------------------------------
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: readeval,init_readeigen,init_readeigen2
      use m_read_bzdata,only: read_bzdata,
     &  ngrp2=>ngrp,nqbz,nqibz,nqbzw,ntetf,n1,n2,n3,qbas,ginv,
     &  qbz,wbz,qibz,wibz,qbzw,
     &  idtetf,ib1bz,
     &  nstbz
      use m_genallcf_v3,only: genallcf_v3,
     &  nclass,natom,nspin,nl,nn,nnc, ngrp,
     &  nlmto,nlnmx,nctot,niw,nw,
     &  alat,ef,dw,delta,iclass,
     &  il,in,im,nlnm,nlnmv,
     &  nlnmc,icore,ncore,
     &  pos,ecore,symgg ! symgg=w(igrp)
      use keyvalue,only: getkeyvalue
      use m_pbindex,only: PBindex
      use m_readqgcou,only: readqgcou

      use mpi_mod, q_comm => comm_1, q_inter => inter_1, rank => mpi_rank

      implicit none

C ... For MPI
      integer num_q_groups
      integer q_start, q_step
      logical last_q,iqinit,MPIinitexists
      integer :: writer1, writer2, i_group, first_q_group, last_q_group

      real(8):: q(3),qgbin(3)
      real(8):: ua=1d0 ! this is a dummy.
!     integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2)
      integer:: ndble=8
      procedure(integer):: nword_in_dble
      real(8),allocatable::  wqt(:),q0i(:,:)
      integer,allocatable :: ngveccB(:,:), !,ngpn(:),ngcni(:),iqib(:),
     &   nx(:,:),nblocha(:),lx(:)
      complex(8),allocatable:: vcoul(:,:),zw(:,:),zw0(:,:),zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),
     &  symope(:,:),
     &  ppb(:,:),
     &  qbze(:,:),qibze(:,:)  !,ecore(:,:)  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: rcxq(:,:,:,:)
c     & ,rcxqmean(:,:,:,:),rcxqmeanc(:,:,:,:) !now rcxqmean is treated as a case of rcxq(nmbas,nmbas)

      logical :: tetra          !  tetrahedron method
      ! so as to reduce the memory usage.
      complex(8), parameter :: img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer,allocatable :: n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:)
      real(8):: ebmx
      integer:: nbmx,mtet(3)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)!,ekxx3(:,:) !,qbzw(:,: )
      real(8) :: anfvec(3)
      logical :: anfexist
      logical :: debug
      integer,allocatable:: ibasf(:)
      real(8),allocatable :: freq_r(:),freq_i(:),freqx(:),wx(:),expa(:),frhis(:)
!     real(8),allocatable :: transaf(:,:)
      integer::iverbose,idalloc
      integer nthreads,omppid
      real(8):: cpusec,walltime0,dwtime
      logical :: realomega, imagomega
c      real(8) wT0L0(2)
      complex(8),allocatable:: epsi(:,:),gbvec(:),zzr(:,:),x0mean(:,:,:),zzr0(:)
      complex(8) :: vcmean
      character*11 fileps
      character*16 filepsnolfc
      logical :: mpi__root=.true. !, paralellx0=.true. !, hist
!     integer :: mpi__size=1, mpi__rank=0
      character(5) :: charnum5

      real(8) :: Emin,Emax,omg2max
      real(8) :: ebasis = 0d0   ! Nonzero => make rcxq in E-basis.  ebasis = truncation
!     real(8) :: emin2,emax2,omg1max
      real(8), allocatable :: freqr2(:)

!     for tetwt5
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      real(8),allocatable :: whw(:)
      integer,allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)

      logical::imagonly=.false.,realonly=.false.
      integer::iopen,maxocc2,iclose,
     &  ixc,iqxini,iqxend,iqxendx,
     &  nw_input,
     &  ifhbe,
     &  nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &  nq0i,i,nq0ix,neps,ngcmx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,is,nwp,
     &   nxx,ifgb0vec
     &   ,iw,nwhis,nw2,noccxv,noccx,iqv
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ibas
     &   ,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0dummy
     &   ,kx,isf,job,nbnbx,nhwtot,noccxvx(2)=-9999
     &   ,jhwtot,ik,ibib,ib1,ib2,imode
     &   ,ngpmx

!     real(8):: epsrng,dwry,dwh
!     complex(8),allocatable:: wgt(:,:,:)

      real(8):: dnorm,omg_c,omg2
      integer::nwin, incwfin,  verbose
      real(8)::efin

      integer:: ngc,mrecg
      real(8):: quu(3),deltaq(3)
      logical:: omitqbz=.false.
      logical,allocatable :: iwgt(:,:,:,:)
      logical :: qbzreg
      integer:: nbcut,nbcut2

      integer,allocatable:: nstibz(:) !Nov2004 Miyake's tote
      real(8),allocatable:: wiw(:)
      real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree
     &   ,pi,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
      integer:: iqbz,iqindx,iflegas
     &   ,ifcor,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1
      real(8) :: eclda_bh,eclda_pz,wk4ec,faca
      logical ::legas
      real(8),allocatable:: ecut(:),ecuts(:) ,totexc(:), trpv(:),trlog(:)
      integer:: necut,iecut

! mode23
      integer:: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
      integer,allocatable:: nxx_r(:)
      real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
      character*3:: charnum3
      character*4:: charnum4
      complex(8),allocatable:: jcoup(:,:)
      real(8)::chg1,chg2,spinmom,schi=1d0
#ifdef COMMONLL
      integer::ll(51**2)
      common/llblock/ll
#else
      integer :: ll
      external ll
#endif
      integer:: ii,ixx

      complex(8),allocatable :: ppovl(:,:),x0meanx(:,:),ppovlz(:,:),ppovlzinv(:,:)

      logical ::  sergeyv=.true., chipm=.false.,nolfco=.false. !sergeyv only ngczero=.false.,
     &  ,epsmode=.false.,normalm=.false., eiqr=.false., crpa=.false., lqall
      logical:: localfieldcorrectionllw,localfieldcorrectionG
      real(8):: ddq(3)
      character*11 ::ttt

      logical ::smbasis

! Feb2006 time-reversal=off case
      logical :: timereversal
      integer:: npm ,jpm,ncc,nw_i
      real(8):: frr

      integer:: istat,imb,imb1,imb2
      integer,allocatable:: imbas(:), imbas_s(:)
!...
      complex(8),allocatable:: eiqrm(:)
      integer:: ifchipmn_mat, ifchipm_fmat !,ifchipm_mat
      real(8):: zz1,zz2,zz3,zz4,dumm1,dumm2
      real(8),allocatable:: mmnorm(:),
     &     momsite(:)
      real(8):: qvv(3)
!!
C     logical :: newaniso,newaniso2
      integer :: ifvcoud,idummy,igb1,igb2,nmbas1,nmbas,iq0
      complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:)
      real(8),allocatable:: vcousq(:)
      real(8):: fourpi,sqfourpi,tpioa,vcou1

!! Eq.(40) in PRB81 125102
c      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
!     complex(8),allocatable:: s_vc(:),vw_k(:),vw_ks(:)
      complex(8),allocatable:: llw(:,:),llwI(:,:)
      integer:: ifrcwx,ircw,nini,nend
C     Outputs generated by w0w0i
      complex(8),allocatable :: w0(:),w0i(:)
      complex(8) :: llmat(3,3)

!     complex(8):: vc1vc2
      integer,allocatable:: nwgt(:,:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
      integer:: l2nl,GWversion
      logical:: eibz4x0,tiii,iprintx,eibzmode
      character(128):: vcoudfile

      integer,parameter :: NULLI=-99999

!! -------------------------------------------------------------------
C      call MPI__Initialize()
C      call MPI__consoleout('hx0fp0')

      walltime0 = dwtime()
C$OMP parallel
      nthreads = omppid(0)
C$OMP end parallel


      hartree=2d0*rydberg()
      pi = 4d0*datan(1d0)
      fourpi=4d0*pi
      sqfourpi=sqrt(fourpi)
      localfieldcorrectionG = mod(GWversion(),10) /= 0 .or. localfieldcorrectionllw()
      call mpi_start
      if (mpi_size > 1) then
         call mpi_end()
         stop "MPI not yet implemented"
      endif
      writer1 = 0
      writer2 = 0
      if (rank == 0) then
         inquire(file='MPI.init',EXIST=MPIinitexists)
         if (MPIinitexists) then
            open(100,file="MPI.init", err=1234 )
            read(100,*) mpi_N1, mpi_N2
            close(100)
         else
            mpi_N1 = mpi_size
            mpi_N2 = 1
         endif
      endif
 1234 continue
#ifdef USE_MPI
      call MPI_Bcast(mpi_N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(mpi_N2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      if (rank == 0) write(6,*) "MPI initialization done"
#endif

C --- computational mode select
c takao keeps only the Sergey mode.
      write(6,"(a)")  'Choose computational mode:'
C     write(6,"(a)")  '    1   : normal, Im axis'
C      write(6,"(a)") '    2   : epsNolfc'
C      write(6,"(a)") '    3   : eps'
C      write(6,"(a)") '    22  : chi^+- NoLFC'
C      write(6,"(a)") '    23  : chi^+- '
      write(6,"(a)") '    11  : normal    Sergey'
C     write(6,"(a)") '    101 : Ecor'
      write(6,"(a)") '    202 : epsNoLFC  Sergey'
      write(6,"(a)") '    203 : eps       Sergey'
      write(6,"(a)") '    222 : chi^+- NoLFC Sergey'
      write(6,"(a)") '    223 : chi^+- Sergey'
      write(6,"(a)") '    12  : total energy Miyake Sergey'
      write(6,"(a)") '    -9999: just show version num'
c      write(6,"(a)") ' #2=iqxini   #3=iqxend' '
c      write(6,"(a)") '  10222 : <e^{iqr}|chi^+-|e^{iqr}> NoLFC'

      call getjobmode(ixc)
      if (ixc == NULLI) then
        if (rank == 0) read(5,*) ixc
      endif
#ifdef USE_MPI
      call MPI_BCAST(ixc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
      if (ixc==NULLI) call rx('Choose computational mode, job=#')

      call headver('hx0fp0',ixc)

      lqall=.true.
      if (ixc<=6.or.ixc==22.or.ixc==23.or.ixc==101) then
        write(*,*)'these modes are removed now'
        call rx( 'these modes are not supported')
! Sergey (Hilbert-transformation) modes
      elseif(ixc==11) then; write(6,*) " OK ixc=11 s.f.'s normal mode "
        normalm=.true.
        lqall=.false.
      elseif(ixc==111) then; write(6,*) " OK ixc=111 normal mode. fullband"
        normalm=.true.
      elseif(ixc==10011) then; write(6,*) " OK ixc=10011 crpa mode "
        normalm=.true.
        crpa=.true.
        call rx('crpa not implemented yet')
!     -- eps mode NoLFC
      elseif(ixc==202) then
        write(6,*) " OK ixc=202  sergey's eps mode Only NoLFC "
        realonly=.true.
        omitqbz=.true.
!     -- eps mode with LFC
      elseif(ixc==203) then
        write(6,*) " OK ixc=203 sergey's eps mode with LFC "
        realonly=.true.
c        iepsmode=203
        omitqbz=.true.
! Total energy modes
      elseif(ixc==12) then
        write(6,*) " ixc=12 Miyake's total energy Sergey--->need to fix this mode"
        call rx( " ixc=12 Miyake's total energy Sergey--->need to fix this mode")
        imagonly=.true.
        ecorr_on=901
!     -- chipm mode NoLFC
      elseif(ixc==222) then
        write(6,*) " OK ixc=222    chipm sergey's "
        realonly=.true.
        omitqbz=.true.
        eiqr =.false. ! .true. aug2012
!     -- chipm mode NoLFC
c      elseif(ixc==10222) then
c        write(6,*) " OK ixc=10222  <q|chipm_0|q> sergey"
c        sergeyv=.true.
c        realonly=.true.
c        omitqbz=.true.
c        eiqr =.true.
!     -- eps mode with LFC
      elseif(ixc==223) then
        write(6,*) " ixc=223  chipm with LFC sergey's -->commented out not. need to fix this mode if necessary."
        call rx( " ixc=223  chipm with LFC sergey's -->commented out not. need to fix this mode if necessary.")
        realonly=.true.
        omitqbz=.true.
        eiqr =.true.
      else
        call rx( ' hx0fp0: mode ixc is not appropriate')
      endif

      if (ixc==202.or.ixc==203.or.ixc==222.or.ixc==223) then
        epsmode=.true.
        if (mod(ixc,200)==22.or.mod(ixc,200)==23) chipm =.true.
        if (mod(ixc,10)==2)                       nolfco=.true.
      endif

C ... files for RPA correlation energy mode.
      if (ecorr_on > 0) then
        ieclog = 8155
        if (ecorr_on==901) then
          ieceig=8156
          open(ieceig,file='rpa_eigen.chk')
          close(ieceig,status='delete')
        endif
        open(ieclog, file='ecorr.chk')
      endif

C      if (newaniso()) then  ! only newaniso is supported now
C        newaniso2=.true.
C      else
C        call rx(' hx0fp0_sc: ony newaniso supported  now')
C      endif

C     if (rank == 0) print *,'Timereversal=',timereversal()
      call cputid (0)

C --- Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
!--------readin data set when you call read_BZDATA ---------------
!       integer::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
!c     &   n_index_qbz
!       integer:: n1,n2,n3
!       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
!       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
!     &    ,wibz(:),qbzw(:,:)
!       integer,allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
!     &    ,nstar(:),irk(:,:),nstbz(:)
!c  ,index_qbz(:,:,:)
!-----------------------------------------------------------------
      call read_BZDATA()

!! Use off-regular mesh (not including Gamma) for qbzreg()=F See hx0fp0.m.sc.F also.
!! This must be consistent with qg4gw.F-mkqg.F
      if (.not.qbzreg()) then
        deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
        do i=1,nqbz
          qbz(:,i) = qbz(:,i) - deltaq/2d0
          write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
         enddo
      endif
      if (rank == 0) then
      do i=1,nqbz
        if (i<10.or.i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
        if (i==10.and.nqbz>18) write(6,"('... ')")
      enddo
      endif

C     oct2005 cases not implemented
      if (smbasis().and.chipm) then
        print *,' smbasis=T & chipm=T is not implemented yet.'//
     &       ' Supply consistent MixSpin for smbasis!'//
     &       ' MixSpin should be converted at the end of hvccfp0.'
        call rx( ' smbasis=T & chipm=T is not implemented yet.')
      endif

c      call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
c      write(6,"(' ScaledGapX0=',f5.3)") sciss

C --- Readin by genallcf ---
      if (epsmode) then
        nwin = -999
      else
        nwin = 0      !Readin nw from NW file
      endif
      incwfin= 0  !use ForX0 for core in GWIN
      efin =  0d0 !readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if (ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWinput')
      tpioa=2d0*pi/alat

!...  These are allocated and setted by genallcf_v3
!      integer::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
!     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
!      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
!      character(120):: symgrp
!      character(6),allocatable :: clabl(:)
!      integer,allocatable:: iclass(:)
!     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
!     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
!     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
!     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
!     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
!     &    occv(:,:,:),unoccv(:,:,:)
!     &   ,occc(:,:,:),unoccc(:,:,:),
!     o    nocc(:,:,:),nunocc(:,:,:)
!      real(8), allocatable::
!     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)

      print *,' ncore=',ncore
      print *,' iclass=',iclass
! ...
      if (chipm.and.nspin==1) call rx( 'chipm mode is for nspin=2')
      nw_input = nw ; print *, 'nw delta=',nw_input,delta
      debug=.false.; if (verbose()>=100) debug=.true.
      iverbose = (verbose()/30)*10

ccccccccccccccccccccc
c      debug=.true.
cccccccccccccccccccc
      if (debug) print *,' end of genallc'
!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom) call rx( ' nclass /= natom ')

!! --- tetra or not
      if (delta <= 0d0) then
        tetra =  .true.
        delta = -delta
        print *,' hx0fp0: tetrahedron mode delta=',delta
      else
        tetra = .false. ! switch for tetrahedron method for dielectric functions
      endif

C ... read dimensions of h,hb
      ifhbe  = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      if (nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
      if (nqbz /=nqbzt) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')

!! --- Readin Offset Gamma --------
      if (debug) print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      write(6,*) ' ### nqibz nq0i=', nqibz,nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      do i=1,nq0i
        read (101, * ) wqt(i),q0i(1:3,i)
      enddo
      nq0ix = nq0i
      do i=1,nq0i
        if (wqt(i)==0d0) then
          nq0ix = i-1
          exit
        endif
      enddo
      neps = nq0i - nq0ix  ! number of zero weight q0p which are used for ixc=2 or 3 mode.
      write( 6,*) ' num of zero weight q0p=',neps
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      close(101)
C      if (.not.newaniso2) then
C        wqtsum = sum(abs(wqt(1:nq0i)))
C        call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
C        if (normalm.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p)) call rx( ' wqtsum of Q0P /=1 ')
C      endif

C --- readin by rdpp ; Radial integrals ppbrd and plane wave part
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx
! qibze(3,nqbze) qbze(3,nqibze)
      nqbze  = nqbz *(1 + nq0i)
      nqibze = nqibz + nq0i
      allocate( qbze(3, nqbze), qibze(3, nqibze))
      call dcopy(3*nqbz, qbz,  1, qbze,1)
      call dcopy(3*nqibz,qibz, 1, qibze,1)
      do i = 1,nq0i
        qibze(:,nqibz+i)  = q0i(:,i)
        ini = nqbz*(1 + i -1)
        do ix=1,nqbz
          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
        enddo
      enddo
      ngrpx = 1
      l2nl=2*(nl-1)
      allocate( nx(0:l2nl,nclass),  nblocha(nclass) ,lx(nclass),
     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:l2nl,nxx, nspin*nclass),
     &   cgr(nl**2,nl**2,(2*nl-1)**2,ngrpx), symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
      call rdpp_v3(nxx, nl,ngrpx, nn, nclass, nspin,symope,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
      if (.not.smbasis()) nblochpmx = nbloch + ngcmx
      allocate(ngveccB(3,ngcmx)) ! work arry

C ... for legas test
      legas = .false.
      INQUIRE (FILE = 'LEGAS', EXIST = legas)
      iqxend = nqibz + nq0i
      print *, ' nqibz nqibze=',nqibz,nqibze

!!... initialization of readEigen !readin m_hamindex
      call init_readeigen(ginv,nspin,nband,mrece) ! EVU EVD are readin in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)

C --- set realomega, imagomega tetra nw niw nwp ifgb0vec ------------
      if (ixc==101) then
        realomega =.false.
        imagomega =.true.
        nw = 0
        nwp= 0
      elseif( ixc==1 ) then !old imagw = 2 case
        realomega =.true.
        imagomega =.true.
      elseif (epsmode.and.(.not.sergeyv)) then !ixc==2.or.ixc==3.or.ixc==22.or.ixc==23 ) then
        call rx('not ready for old eps mode')
c$$$        print *, 'mode=',ixc
c$$$        realomega =.true.
c$$$        imagomega =.false.
c$$$        niw = 0
c$$$        call getkeyvalue("GWinput","EPSrange",epsrng )
c$$$        call getkeyvalue("GWinput","EPSdw",  dwry )
c$$$        epsrng=epsrng/2d0
c$$$        dw = dwry/2d0
c$$$        if (dw==0d0) then
c$$$          nw = 1
c$$$        else
c$$$          nw = (epsrng/2d0 - 1d-10)/(dw/2d0) + 2 !epsrng/2d0 corresponds to in a.u.
c$$$        endif

c$$$      if (ixc==6) then
c$$$! test modes. ixc=6 Spectrum function (Img part) along the Real axis with tetwt5. Histgram method.
c$$$        realomega = .true.
c$$$        imagomega = .false.
c$$$        tetra     = .true.
c$$$        niw = 0
c$$$! ---  For tetwt5 ---  the tetrahedron weight for spectrum function (imaginary part)
c$$$!   Histogram bins are specified by freq_r(1:nwp)
c$$$!     nwp=nw+1; frhis(1)=0
c$$$!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
c$$$!     The last  bin  is     [frhis(nw), frhis(nwp)].
c$$$
c$$$! ... These parameters specifies a test histogram; Sergey's mesh just for test modes.
c$$$        nw0 = 200    !100    800
c$$$        dwh = 0.01d0 !0.02d0 0.0025d0 !in hartree
c$$$        call findemaxmin(nband,qbze,nqbze,nspin,emax,emin)
c$$$c$$$ Even if we neglect this section, result is not so bad---anyway we may need to
c$$$c$$$ fix this block if necessary. Oct2005
c$$$c$$$        if (bzcase()==2) then
c$$$c$$$          allocate(qbz2(3,nqbz))
c$$$c$$$          do iq=1,nqbz
c$$$c$$$           qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
c$$$c$$$          enddo
c$$$c$$$          call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
c$$$c$$$          emax=max(emax,emax2)
c$$$c$$$          emin=min(emin,emin2)
c$$$c$$$          deallocate(qbz2)
c$$$c$$$        endif
c$$$        if (nctot > 0) Emin = minval(ecore(:,1:nspin))
c$$$        omg2max = (Emax-Emin)*.5+.2d0    !(in Hartree) covers  all relevant omega, +.2 for margin
c$$$        omg1max = dwh*(nw0-1)
c$$$        nwp = int(sqrt(omg2max*(2*nw0-1d0)/dwh-(nw0**2-3*nw0+1d0)))+1 ! + 1 for margin
c$$$        nw  = nwp-1
c$$$        write(*,*) Emax,Emin,nw0,nw  ! nwp is new max number in frequency array
c$$$        write(6,'(a32,2i7,2d15.3)')'hx0fp0: nw0,nw,omg1max,omg2max='
c$$$     &              , nw0,nw,  omg1max,omg2max
c$$$        if (nw <= nw0) stop 'hx0fp0:ixc==[456] nw2 <= nw'
c$$$        allocate(freq_r(0:nw))
c$$$        do iw=0,nw  !This is a test mesh by Sergey.Faleev
c$$$          if (iw+1<=nw0) then;  freq_r(iw)=dwh*iw
c$$$          else;  freq_r(iw)=dwh*((iw+1)**2+nw0**2-3*nw0+1)/(2*nw0-1d0)
c$$$          endif
c$$$        enddo
c$$$! freq_r(iw) is linear for iw<=nw and quadratic for nw<iw<=nw2
c$$$! freq_r(iw) chosen in such a way that it is continues with
c$$$!   nw nwp=nw+1 freq_r(1:nwp) are used after here.
c$$$        allocate(frhis(nwp))
c$$$        frhis=freq_r(0:nw)
c$$$        nwhis=nw
c$$$

C ... Sergey Faleev's normal mode.  First accumulate Im chi . Im chi is K-K transformed to obtain Re chi
      elseif(sergeyv) then
      realomega = .true.
      imagomega = .true.
      tetra     = .true.
      call getkeyvalue("GWinput","omg_c",omg_c )
      if (rank == 0) then
        write(6,"(/' Frequency parameters:')")
        write(6,"('   frequencies on Im omega axis    =',i5)") niw
        write(6,"('   mesh spacing dw along Real axis =',f12.6,' Ha')") dw
        write(6,"('   becomes 2*dw at:                 ',f12.6)") omg_c
        write(6,"('   Fermi energy                    =',f12.6)") ef
      endif

C      if (bzcase()==2) then
C        write(6,'(" dq_bzcase2=",3f9.4)')dq_bzcase2
C      endif
      call findemaxmin(nband,qbze,nqbze,nspin,emax,emin)

c$$$ Even if we neglect this section, result is not so bad---anyway we may need to
c$$$ fix this block if necessary. Oct2005
c$$$        if (bzcase()==2) then
c$$$          allocate(qbz2(3,nqbz))
c$$$          do iq=1,nqbz
c$$$           qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
c$$$          enddo
c$$$          call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
c$$$          emax=max(emax,emax2)
c$$$          emin=min(emin,emin2)
c$$$          deallocate(qbz2)
c$$$        endif
      if (nctot > 0) emin = minval(ecore(:,1:nspin))
      omg2max = (emax-emin)*.5d0+.2d0 !(in Hartree) covers  all relevant omega, +.2 for margin
      if (rank == 0) write(*,"('   emin emax',2f12.6,'  range',f12.6)") emin, emax, omg2max*2
      nw2 = int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
      allocate(freqr2(nw2))     !+1 b/c (iw-1)
      do  iw = 1, nw2           !linear + quadratic term
        freqr2(iw)=dw*(iw-1)+dw**2/2./omg_c*(iw-1)**2
      enddo
      if (nw2 < 2) call rx( "hx0fp0: nw2 < 2")
      if (dw*(nw-2) > freqr2(nw2-1)) call rx("hx0fp0: dw*(nw-2) > freqr2(nw2-1)")
      nwhis = nw2-1
      allocate(frhis(1:nwhis+1))
      frhis = freqr2(1:nwhis+1)
      if (rank == 0) write(*,"('   number of histogram bins         =',i4/)") nwhis

      nw = nw2-1                ! nw+1 is how many points of real omega we use
                                ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
                                ! maximum nw=nw2-1 because nwhis=nw2-1
      do iw = 3, nw2-1
                                !nw is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
                                !Here frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
                                !nw was constructed such that omg=dw*(nw-2)> all relevant frequencies needed
                                !for correlation Coulomb Wc(omg),
                                !and one more point omg=dw*(nw-1) needed for extrapolation.
                                !Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
                                !and one more point omg=frhis_m(nw) needed for extropolation
                                !used in subroutine alagr3z in  sxcf.f.
        omg2 = (freqr2(iw-2)+freqr2(iw-1))/2d0
        if (omg2 > dw*(nw_input-3)) then
          nw = iw
          exit                  ! 'nw_input' is only used to get maximum frequency for dressed coulomb line
        endif
      enddo
      if (epsmode.and.sergeyv) then !epsmode==202.or.iepsmode==203) then
        nw  = nwhis-1
        niw = 0
      endif
      allocate(freq_r(0:nw))
      freq_r(0) = 0d0
      do  iw = 1, nw
        freq_r(iw) = (frhis(iw)+frhis(iw+1))/2d0
      enddo
      nwp = nw+1
      delta = 0d0                 !only for writeemesh

!!    Tabulate frhis
C      if (onceww(1) .and. rank == 0) then
C        write(6,*) ' Frequency bins to accumulate Im chi (a.u.)'
C        do ihis= 1, min(10,nwhis)
C          write(6,"(' ihis Init  End=', i5,2f13.6)") ihis,frhis(ihis),frhis(ihis+1)
C        enddo
C        write(6,*) 'ihis  ...'
C        do ihis= max(min(10,nwhis),nwhis-10), nwhis
C          write(6,"(' ihis Init  End=', i5,2f13.6)") ihis,frhis(ihis),frhis(ihis+1)
C        enddo
C      endif

C     timereversal=F is implemented only for tetra=T and sergeyv=T
      npm = 1; nw_i = 0
C     iprintx = Timereversal(); print *,'Timereversal=',iprintx
      if (.not.timereversal()) then
        npm = 2
        nw_i = -nw
        if (.not.tetra)     call rx( ' Require tetra=T with timereversal=off')
        if (chipm)          call rx( ' Require timereversal with chipm.'//
     &    ' Takao think npm=1 works OK even for chipm=T and TimeReversal=off. But need to check' )
      endif
!! For chipm=T for TimeReversal=off(psi and psi* are not degenerated), Takao think npm=1 works OK. But need to check.
!! In addition, we have to consider Spin-orbit interaction correctly...

      if (realomega.and.(.not.sergeyv)) then !     elseif (realomega.and.(ixc==22.or.ixc==23)) then
        nwp = nw+1 !Feb2006 I want to unify the meaning of nwp as nwp=nw+1
        ! Now x0kf calculates one more freq_r along real axis.
        ! at freq_r(nw). --->It WAS not calculated. It wat up to freq_r(nw-1).
        ! But not tested well.
        allocate(freq_r(0:nw))
        do iw = 0,nw
          freq_r(iw) = dw*iw ! freq_r
        enddo
      endif

C ... get eigenvector corresponds to exp(iqr) (q is almost zero).
      if (epsmode) then !iepsmode/=0) then ;  print *, ' read in Mix0vec'
        allocate(epsi(nw_i:nw,neps)) !5July2005 nwp should be used after it is defined!
        ifgb0vec = iopen("Mix0vec",1,3,0)
      endif

C     write frequencies to file freq_r, used by hsfp0
      if (realomega .and. rank == 0) then
        open(UNIT=3111,file='freq_r')
        write(3111,"(2i8,'  !(a.u.=2Ry)')") nwp, nw_i
        do iw = nw_i, -1
          write(3111,"(d23.15,2x,i6)") -freq_r(-iw),iw
        enddo
        do iw= 0, nw
          write(3111,"(d23.15,2x,i6)") freq_r(iw),iw
        enddo
        close(3111)
      endif

C --- WVI only for imagonly for ixc==12
      if (imagonly.and.sergeyv) then
        realomega =.false.
        imagomega =.true.
      endif

C --- epsPP noLFC mode for ixc==13
      if (realonly.and.sergeyv) then
        realomega =.true.
        imagomega =.false.
      endif

C --- Imaginary frequencies
C     Points x for Gaussian quadrature between (0,1) and w=(1-x)/x
      if (imagomega) then
        allocate(freq_i(niw),freqx(niw),wx(niw),expa(niw) )
        call freq01 (niw,ua,freqx,freq_i,wx,expa)
        allocate(wiw(niw))
        do iw = 1, niw
          wiw(iw)=wx(iw)/(2d0*pi*freqx(iw)*freqx(iw))
        enddo
c        deallocate(freqx,wx,expa)! freqx,wx is used for Miyake's tote mode !Nov2004
      else
        allocate(wiw(1))
      endif

c ... for Miyake tote mode Nov2004
      if (ecorr_on>0) then       !it was bzcase()==2 Was it bug?
        allocate(nstibz(nqibz))
        do iq=1,nqibz
          iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
          nstibz(iq) = nstbz(iqbz)
        enddo
      endif

C --- Initialize tetrahedron integration
      call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband )
      call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
      mtet=(/1,1,1/)
      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
      if (tetra) then
        allocate(nbnb(nqbz,npm),nbnbtt(nqbz,npm),ekxx1(nband,nqbz),ekxx2(nband,nqbz))
      endif
C ... maximum no. occupied valence states
      noccxv = maxocc2(nspin,ef,nband,qbze,nqbze)
      if (noccxv>nband) call rx( 'hx0fp0: all the bands filled! too large Ef')
      noccx = noccxv + nctot
      nprecx = ndble            !We use double precision arrays only.

C --- write energy mesh ---
      ifemesh = iopen('emesh.hx0fp0',1,-1,0)
      if (.not.realomega) then; allocate(freq_r(nwp)); freq_r=0.0d0; endif
      if (.not.imagomega) then; allocate(freq_i(1)); freq_i=0.0d0; endif
      call writeemesh(ifemesh,freq_i,niw,freq_r,nwp,delta)

C ... ppb
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass,nspin) )

C ... This is just to get nblochpmx
      if (smbasis()) then
        call getngbpomat(nqibz+nq0i,nblochpmx,i)
      endif
      mrecl = nprecx*2*nblochpmx*nblochpmx/nword_in_dble()

C     Open WV.d and allocate zw
      if (rank == 0) then
        ifwd = iopen('WV.d',1,-1,0)
        write (ifwd,"(1x,i3,i10,2i6,i5,3i5)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz+nq0i-1,nw_i
        ifwd = iclose('WV.d')
      endif
      allocate(zw(nblochpmx,nblochpmx))
      nspinmx = nspin

!! --- READIN ANFcond. we don't use this path so much ---
      anfexist = .false.
      if (rank == 0) then
        inquire(file='ANFcond',EXIST=anfexist)
        if (anfexist) print *,'turn on ANFcond'
      endif
#ifdef USE_MPI
      call MPI_Bcast(anfexist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
#endif
      if (anfexist) then
        call rx('ANF mode no longer supported ... sorry')
      endif

C      if (anfexist) then
C        print *,' WARN: this mode is not maintained recently... this may cause a problem'
C        allocate(ibasf(natom),transaf(1:3,natom))
C        ifianf = 211
C        open(ifianf,file='ANFcond')
C        read(ifianf,*)
C        read(ifianf,*)
C        read(ifianf,*) anfvec(1:3)
C        write(6,'(" ANFcond=",3f13.6)') anfvec
C        do ibas=1,natom
C          read(ifianf,*) ibas1, ibasf(ibas)
C          if (ibas1/=ibas) call rx( ' readin ANFcond error:ibas1/=ibas')
C          write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
C        enddo
C        close(ifianf)
C        nspinmx =1
C      endif
C      if (smbasis() .and. anfexist) then
CC        print *,' smbasis=T & anfexist=T is not yet.'//
CC     &       'At leaset anfx0k shold be changed'
C        call rx( ' smbasis=T & anfexist=T is not yet.')
C      endif
C      if ((.not.timereversal()) .and. anfexist) then
C        print *,' timereversal=F & anfexist=T is not yet examined.'//
C     &       'At least you have to consider about anfx0k and so.'
C        call rx( ' timereversal=F & anfexist=T is not yet.')
C      endif

C ... Read here only for printout
      if (rank == 0) then
        call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
        call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
        if (nbcut+nbcut2 > 0) write(6,"(' nbcut nbcut2 =',2i5)") nbcut,nbcut2
      endif

C --- ppb ---
      do  is = 1,nspin          ! <Phi(R,L,n) Phi(R,L',n') B(R,i)> for each spin
        irot = 1
!       Parallelized with OMP and can run on all nodes: no I/O
        call ppbafp_v2(irot,ngrpx,is,nspin,
     .    il,in,im, nlnm,
     .    nl,nn,nclass,nlnmx,
     .    mdimx,lx,nx,nxx,  !Bloch wave
     .    cgr,nl-1,         !rotated CG
     .    ppbrd,            !radial integrals
     o    ppb(1,is))
      enddo
      if (rank == 0 .and. debug) print *, ' end of ppbafp_v2'

C ... Set iqxini
      if (omitqbz) then
        iqxini= nqibz + 1
C      elseif(newaniso2) then
C        iqxini=1
C      elseif(bzcase()==1) then
C        iqxini = 2
      else
        iqxini = 1
      endif

C... check write 1st part for Ec mode to ecorr.chk Nov2004
      if (ecorr_on>0) then
         call rx( ' need to check subroutine getwk and so on')
        if (ecorr_on >0.and. (.not.imagomega))
     &    call rx( ' hx0fp0: ecorr_on .and. (.not.imagomega)')
        write(ieclog,
     &          "('   iq                q                       wk')")
        do iqq = iqxini,iqxend
          call getwk(iqq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)
          write(ieclog,"(i5,3x,3f12.8, f15.5)") iqq, q, wk4ec
        enddo
        write(ieclog,*)
        write(ieclog,"('   iw omega(Ry)       wiw')")
        do iww=1,niw
          write(ieclog,"(i5, f10.5, f10.5)") iww,2d0*freq_i(iww),wiw(iww)
        enddo
        write(ieclog,*)
        write(ieclog,"(' Note:IntWgt=wk*wiw.',
     &     '  Ec =\sum_{k,iw) IntWgt(k,iw)*ecqw(k,iw)')")
        close(ieclog)
        open(ieclog,file="ecorr.chk",access='append')

        call getkeyvalue("GWinput","necut_p",necut, default=1 )
        allocate(totexc(necut),trpv(necut),trlog(necut))
        totexc = 0d0
        trpv   = 0d0
        trlog  = 0d0
        if (.not.sergeyv) necut=1
      else
        necut=1
      endif
      allocate(ecut(necut),ecuts(necut))
      call getkeyvalue("GWinput","ecut_p" ,ecut, necut,default=(/1d10/) )
      call getkeyvalue("GWinput","ecuts_p",ecuts,necut,default=(/1d10/) )
      call getkeyvalue("GWinput","ebasis",ebasis,default=0d0)
      if (necut/=1) call rx( 'Sorry! this version is not for necut=1')

      if (ebasis <= -1d0) ebasis = -nint(-ebasis)
!      if (rank == 0) then
        write(*,"(1x)")
        if (ebasis <= -1d0) then
          write(*,"( a,i4,a)") ' rcxq generated in E basis, remove',nint(-ebasis),' states'
        elseif (ebasis /= 0d0) then
          write(*,"( a,1p,E11.3)") ' rcxq generated in E basis, Ecut =', ebasis
        endif
        if (.not. timereversal()) write (*,"(' Time reversal suppressed')")
!      endif


      if (chipm.and.(.not.legas)) then
        if (smbasis()) then ! sanity check
          print *,'This version is not for spin sus for smbasis=T.'//
     &         'You need to modify MixSpin.'
          call rx( 'This version is not for spin sus for smbasis=T')
        endif

!Sep2006
c        nmbas=1
c        allocate(imbas(nmbas),imbas_s(nmbas))
c        imbas(1)=1
c        imbas_s(1)=1

c        allocate(imbas(natom),imbas_s(natom),iibas(natom))
c        print *,' goto MagAtom read natom=',natom
c        nmbas=natom
c        istat=-9999 ! istat=-9999 means noumber of readin arguments is returened in istat.
c        call getkeyvalue("GWinput","MagAtom",
c     &        imbas,nmbas,status=istat)
c        print *,' end of MagAtom read istat=',istat
c       nmbas = natom
c        do i=1, natom
c          iibas(i) = i
c        enddo
c        do i=1, istat
c          ibas = imbas(i)
c          if (ibas<0) iibas( abs(ibas) ) = ibas
c        enddo
c        write(6,"('Readin MagAtom nmbas =',i3,' imbas_s= ',10i3)")
c     &  nmbas, imbas(1:nmbas)
c        imbas_s(1:natom) =     iibas(1:natom)
c        imbas(1:natom)   = abs(iibas(1:natom))
c        write(6,"('   Use MagAtom nmbas =',i3,' imbas_s= ',10i3)")
c     &  nmbas, imbas_s(1:nmbas)

        nmbas=natom
        allocate(imbas(nmbas),imbas_s(nmbas))
        istat=-9999 ! istat=-9999 means noumber of readin arguments is returened in istat.
        call getkeyvalue("GWinput","MagAtom",
     &        imbas,nmbas,status=istat)
        nmbas = istat
        print *
        write(6,"('Readin MagAtom nmbas =',i3,' imbas= ',10i3)") nmbas,imbas(1:nmbas)
        imbas_s(1:nmbas) = imbas(1:nmbas)
        imbas(1:nmbas)   = abs(imbas(1:nmbas))
c
        allocate(jcoup(nw_i:nw,neps) )
        allocate( svec(nbloch,nmbas) )  !sep2006
        svec=0d0
        allocate(cvec(nbloch,nmbas),momsite(nmbas),mmnorm(nmbas))              !May2007
        cvec=0d0
        do imb = 1, nmbas
          ibas= imbas(imb)
c          open(ifv,file='MixSpin.'//charnum3(ibas))
          ifv = iopen ('MixSpin.'//charnum3(ibas),1,3,0)
          read(ifv,*) ibasx,lxx
          allocate(nxx_r(0:lxx))
          do i=0,lxx
            read(ifv,*) nxx_r(i)  !   write(6,"(2i5,d13.6)") nxx_r(i)
          enddo
          allocate(spinvec((lxx+1)**2,maxval(nxx_r)))
          allocate(consvec((lxx+1)**2,maxval(nxx_r)))
          spinvec=0d0
          do ilmx = 1, (lxx+1)**2
            lb = ll(ilmx )       !  print *,' lb=',lb,lxx,ilmx
            do ixx = 1, nxx_r(lb)  !  print *,' nn=',nn,nxx_r(lb)
              if (ilmx==1) then
                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),chg1,chg2
     &       ,consvec(ilmx,ixx)
              else
                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),dumm1,dumm2
     &       ,consvec(ilmx,ixx)
              endif
!           write(6,"(2i5,d13.6)") ilmx, ixx, spinvec(ilmx,ixx)
            enddo
          enddo

! Calculate ChiPM. So sign of omega should be correct.
          if (imb==1) then !determine spin direction with respect to ibas=imbas(imb=1)
            spinmom=(chg1-chg2)
            schi=1d0
            if (spinmom<0d0) then
              if (.not.sergeyv)  then
                freq_r= -freq_r
              else
                schi  = -1d0   ! This affects to dpsion. Obtained results
                ! should be the same in both mode.
              endif
            endif
          endif

C ... ReOrdering of spinvec in natom ordering...
          i=0
          if (ibas>1) i= sum(nblocha(1:ibas-1))
          do lb  = 0, lx (ibas)
            do nb  = 1, nx (lb,ibas)
              do mb  = -lb, lb
                i = i+1
                ilmx = lb**2+ lb+ mb +1
                svec(i,imb) = spinvec(ilmx,nb)
                cvec(i,imb) = consvec(ilmx,nb)
                write(6,"(' i lb mb svec svec**2=',3i4,2d13.5)")
     &      i,lb,mb,svec(i,imb),svec(i,imb)**2
              enddo
            enddo
          enddo
          deallocate(nxx_r,spinvec,consvec)
          close(ifv)
!       enddo

c         write(ifchipmlog,"(a)")
c     &   ' --- Table of mcm v.s. ss cutoff.'//
c     &   '  : mcm means <m|chi^{-1}|m>. '//
c     &   '  ss is the SVD values of chi^{-1}.'
c         if (.not.nolfco) then
c          write(ifchipmlog,"(13x,a)")'q              omega'//
c     &    '   mcm(ss<ssmin+2Ry)    mcm(ss<ssmin+10Ry)'//
c     &    '    mcm(ss<ssmin+5Ry)   mcm(ss<ssmin+100Ry)'//
c     &    '    mcm(all)'
c         endif
          mmnorm (imb) = sqrt(sum(svec(:,imb)**2))
          momsite(imb) = chg1-chg2
c        write(6,"(' svecsum=',d23.15)") sum(svec(:,imb)**2
c          write(ifchipmlog,"(2d23.15,' ! mmom mmnorm')")momsite(imb),mmnorm(imb)
          write(6,"( 'mmom mmnorm= ',2f14.10)")  momsite(imb),mmnorm(imb)
        enddo
c        ifx= iclose('ChiPM.log')

      endif ! (chipm.and.(.not.legas))

! I assume 1 is for majority for eiqr case.
c      if (ix==10222) then
cc        schi=1d0 !1d0 means Majority is isp=1. If Majority is isp=2, use schi=-1d0.
cc        allocate(jcoup(nw_i:nw,neps))
c        mmnorm=1d0
c      endif

c! nmbas is for rcxqmean
c      if (chipm .and. nolfco) then !.and. ix/=10222) then
cc        nmbas = nmbas
c      else
c        nmbas = 1
c      endif
c      if (epsmode.and.sergeyv.and.nolfco) then
c        allocate( rcxqmean(nwhis,npm,nmbas,nmbas))
c        if (debug) write(6,"('fff:',3i5)") nwhis,npm,nmbas
c      else
c        allocate( rcxqmean(1,1,1,1)) !dummy
c      endif


       if (chipm) allocate(eiqrm(nmbas))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cctakao
c$$$      allocate( x0meanx(nmbas,nmbas) )
c$$$      allocate( x0mat(nmbas,nmbas),x0matinv(nmbas,nmbas) )
c$$$      do 1101 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
c$$$         if (iq==iqxini+2) exit
c$$$         q = qibze(:,iq)
c$$$         print *,'aaaaaaaaaa q=',q
c$$$          read(ifgb0vec,*) qgbin(1:3),ngbin,igc0,dnorm
c$$$          if (sum(abs(q))==0d0)then
c$$$            if (sum(qgbin**2) >1d-7)stop'qgbin=0 xxx See hvccfp0'
c$$$          elseif(sum(abs(qgbin(1:3)-q)) >1d-8)then
c$$$            stop'qgbin inconsistent'
c$$$          endif
c$$$          write(6,"(' --- Readin Mix0vec: ',3d13.6,2i5,d18.8)")
c$$$     &      qgbin(1:3),ngbin,igc0,dnorm
c$$$c          if (ngb/=ngbin) stop 'hx0fp0: ngb/=ngbgin'
c$$$          ngb=ngbin
c$$$          write(6,"(' ngb nwp niw=',3i8)")ngb,nwp,niw
c$$$          nmbas=1
c$$$          allocate( gbvec(ngb),zzr(ngb,1),x0mean(nw_i:nw,1,1))
c$$$          x0mean=0d0
c$$$          do i=1,ngb
c$$$            read(ifgb0vec,"(4d24.15)") zz1,zz2,zz3,zz4
c$$$            gbvec(i)= dcmplx(zz1,zz2)
c$$$            zzr(i,1)= dcmplx(zz3,zz4)
c$$$          enddo
c$$$          write(6,"(' normchk=',255d23.15)") sum( dconjg(gbvec)*zzr(:,1) )
c$$$     &        ,sum(abs(gbvec(:))), sum(abs(zzr(:,1)))
c$$$            allocate(eiqrm(nmbas))
c$$$            do imb=1,nmbas
c$$$            eiqrm(imb)= sum( dconjg(gbvec(1:nbloch))*svec(1:nbloch,imb) )
c$$$            write(6,"(' <eiqr|m> =',255d23.15)") eiqrm(imb)
c$$$            if (imbas_s(imb)<-1) eiqrm(imb)= -eiqrm(imb)
c$$$            enddo
c$$$            write(6,"('<eiqr|m>:Set \pm in GWinput(for stuggard chi)')")
c$$$              iqixc2 = iq- (nqibz+nq0ix)
c$$$              ifx = iopen ('StonerNLFC.dat',1,3,0)
c$$$              read(ifx,*) jzero2
c$$$              ifx= iclose('StonerNLFC.dat')
c$$$
c$$$              ifchipm2=iopen(
c$$$     &                 'ChiPM'//charnum4(iqixc2)//'.nolfc.mat',1,3,0)
c$$$              do iw=1,10
c$$$              read(ifchipm2,
c$$$     &        '(36x,2x,20x,2x,255d23.15)') x0meanx(:,:)
c$$$              write(6,'("xxx x0mat=",255d13.5)') x0meanx
c$$$              x0matinv=x0meanx
c$$$              call matcinv(nmbas,x0matinv)
c$$$              do i=1,nmbas
c$$$                x0matinv(i,i)= x0matinv(i,i) - jzero2 ! (chipm_0^+-)^-1 - I
c$$$              enddo
c$$$              x0mat = x0matinv
c$$$              do i=1,nmbas
c$$$              x0mat(i,i) = x0mat(i,i)+ img*1d-30 ! to avoid inversion error.
c$$$              enddo
c$$$
c$$$              call matcinv(nmbas,x0mat) !this is full x0_+-
c$$$              trr = sum( eiqrm*matmul(x0mat,dconjg(eiqrm)) ) !*mmnorm
c$$$              write(6,
c$$$     &        '("ttt",3f12.8,2x,f10.5,2x,2d23.15,2x,2d23.15)') q, 2*schi*frr, trr,1d0/trr
c$$$              enddo
c$$$            ifx=iclose(ifchipm2)
c$$$ 1101  continue
c$$$          stop 'xxxxxxxxxxxxxxxxxxxxxxxx'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      eibzmode=eibz4x0()
!! -- EIBZ mode for nolfco -----------------------
!! If eibzmode=T, it is efficient but can slightly break crystal symmetry.
!! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
c      if (nolfco) then
c          eibzmode = .false.
c      endif
c!! ----------------------------------------------
!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
!!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
      if (eibzmode) then
!! commentout block inversion Use iqxendx=iqxend because of full inversion
        if (verbose()>=30) call cputid2(' start eibzmode block',0)
        iqxendx=iqxend
        if (epsmode) iqxendx=iqxend
        allocate(nwgt(nqbz,iqxini:iqxendx),igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
     .    eibzsym(ngrp,-1:1,iqxini:iqxendx))
!! Check timereversal is required for symmetrization operation or not. If tiii=timereversal=F is enforced,
!! the symmetrization procedure in x0kf_v4h becomes a little time-consuming.
        if (mpi__root) iprintx=.true.
        call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,
     i    timereversal(),ginv,iprintx,nwgt,igx,igxt,eibzsym,tiii)
        call PBindex(natom,lx,l2nl,nx) !all input. Returns required index stored in arrays in m_pbindex.
! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
        if (verbose()>=30) call cputid2(' complete eibzmode block',0)
        call readqgcou()        !no input. Read QGcou and store data into variables.
c     call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.
      else                      !dummy allocation to avoid undefined pointers
         iqxendx = iqxend
         allocate(nwgt(ngrp,iqxini:iqxendx),igx(1,1,iqxini:iqxendx),igxt(1,1,iqxini:iqxendx),eibzsym(1,1,iqxini:iqxendx))
         nwgt = 1
       endif
       if (rank == 0) print *,'Timereversal=',Timereversal()

!! == Calculate x0(q,iw) and W == main loop 1001 for iq.
!! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
!! Thus it is necessary to do iq=1 in advance to perform iq >nqibz, modify do 1001 loop.
       allocate(llw(nw_i:nw,nq0i),llwI(niw,nq0i))
       llw=1d99; llwI=1d99      ! for debugging

       if (rank == 0) then
         write(*,321) iqxend-iqxini+1, nq0i
  321    format(/' --- Start Q - loop ...',i5,' mesh points including',i5,' offset Gamma points')
      endif

      num_q_groups = mpi_N1
!     write(0,'(" [",I4,"] Ready for q-loop:",2F15.5)') rank, (dwtime() - walltime0)/60, cpusec()/60
      call mpi_make_groups(comm_world,iqxend-iqxini,num_q_groups,q_comm,q_inter)

#ifdef USE_MPI
      write(*,'("MPI: rank = ",I4,", q_comm%ID = ",I4,", q_comm%group_ID = ",I4)')
     &          rank, q_comm%ID, q_comm%group_ID

      write(*,'("MPI: rank = ",I4,", q_inter%ID = ",I4,", q_inter%group_ID = ",I4)')
     &          rank, q_inter%ID, q_inter%group_ID
#endif

!      write(0,'(" [",I4,"](",I3,",",I3,") Ready for q-loop:",2F15.5)') rank, q_comm%group_ID, q_comm%ID,
!     &     (dwtime() - walltime0)/60, cpusec()/60

      if (q_comm%size > 1) writer2 = 1
      q_start = q_comm%group_ID
      q_step =  num_q_groups
      last_q = .false.
#ifdef USE_MPI
      if (rank == 0) write (*,*) "Parallelizing q-loop over ",num_q_groups,"Processor groups."
#endif

!ML --------------------------------------------------------------------------------------
!ML MPI remark:
!ML
!ML   We need to figure out on which group the iq will run !
!ML
!ML   Then, the processor-group of iq=1 (group_ID == 0) can send an non-blocking
!ML   message to the group of (iqxend) to indicate that the files WV*.00001.*
!ML   are written.
      do i_group = 0, num_q_groups-1
        do iq = iqxini + i_group, iqxend, q_step
          if (iq==1)      first_q_group = i_group
          if (iq==iqxend) last_q_group  = i_group
        enddo
      enddo
#ifdef USE_MPI
      if (rank == 0) write(0,'(" [",I4,"] first and last q groups: ",2I10)') rank, first_q_group, last_q_group
#endif
!ML   Now we should have the group_ID which will work on last_q
!ML   last_q is known on all processes.
!ML
!ML   the head of group 0 can communicate with the head of last_q_group through q_inter
!ML
!ML   we can already initialise the receive buffer on last_q:
#ifdef USE_MPI
      if (first_q_group /= last_q_group) then
         if (q_inter%ID == last_q_group) then
            call MPI_Irecv(i_ready, 1, MPI_INTEGER, first_q_group, 1, q_inter%communicator, i_request, ierror)
            write(0,'(" [",I4,"] I will be running the last q vector:",2I10)') rank,q_inter%ID
            write(0,'(" [",I4,"] Waiting for message from :",I10)') rank,first_q_group
         endif
      else
         if (q_inter%ID == last_q_group) then
            write(0,'(" [",I4,"] I will be running the last q vector:",2I10)') rank,q_inter%ID
            write(0,'(" [",I4,"] No need to communicate!")') rank
         endif
      endif
#endif
      iqinit=.true.
      do 1001 iq = iqxini + q_start, iqxend, q_step ! q=(0,0,0) is omitted!  ! for each q, calculate chi(q)

#ifdef USE_MPI
C        if (q_comm%ID == 0) write (*,*) "q point ",iq," running on",q_comm%group_ID,"(",rank,")"
        write(*,*) "Beginning of q-loop:", rank, q_comm%group_ID, q_comm%ID, " iq = ",iq
!       write (0,'(" [",I4,"](",I3,",",I3,") Running iq = ",I4)') rank, q_comm%group_ID, q_comm%ID, iq
!       if (q_comm%ID == writer1 .and. (ixc==101.or.normalm)) ifrcw  = iopen('WVR.'//charnum5(iq),0,-1,-mrecl)
        if (q_comm%ID == writer2 .and. (ixc==101.or.normalm)) ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
#else
        if (normalm) ifrcw = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
        if (q_comm%ID == writer2 .and. (ixc==101.or.normalm)) ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
#endif

        call cputid (0)
        q = qibze(:,iq)
        if (q_comm%ID == 0) then
          write(6,"(' q vec =',3f14.6)") q
          call readqg('QGcou',q,ginv,quu,ngc,ngveccB) ! q was qq
        endif
#ifdef USE_MPI
        call MPI_BCast(ngc, 1, MPI_INTEGER, 0, q_comm%communicator, ierror)
        call MPI_BCast(quu, 9, MPI_REAL8, 0, q_comm%communicator, ierror)
        call MPI_BCast(ngveccB, 3*ngc, MPI_INTEGER, 0, q_comm%communicator, ierror)
#endif
!! Caution : confusing point
!!  ngc by QGcou is shown at the bottom of lqg4gw.
!!  ngc read from PPOVL are given by rdata4gw ---> ngc(iq>nqibz )=ngc for q=0
!!  Note that  ngc(iq>nqibz )=ngc (q=0), because when it is generated in mkqg.F

C       if (newaniso2.and.iq==1) then ! *sanity check
        if (iq==1) then ! *sanity check
          if (sum(q**2)>1d-10) then
            call rx( ' hx0fp0: sanity check. |q(1)| /= 0')
          endif
        endif

C   --- Read diagonalized Coulomb interaction ---
        ngb = nbloch + ngc !ngb is readin from vcoul 25jan2006
#ifdef USE_MPI
         write(0,"('[',I4,'](',I3,',',I3,')  do 1001: iq=',i6,' of',i6,'  q=',3f12.6)")
     &     rank,q_comm%group_ID, q_comm%ID, iq, iqxend, q
#else
         write(6,"(' do 1001: iq=',i6,' of',i6,'  q=',3f12.6)") iq, iqxend, q
#endif

!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
C       if (newaniso2.and.(.not.chipm)) then
        if ((.not.chipm)) then
C         Each q-point has its own file.  BUT, only the head of each q-group should do the reading!!
          iqv = iq
          if (iqv>nqibz .and. (ixc==1 .or. ixc==11)) iqv = 1 ! Use vcoul at q=(0,0,0) for offset gamma points
          if (q_comm%ID == 0) then
            vcoudfile='Vcoud.'//charnum5(iqv) ! file is closed at the end of do 1001.
            ifvcoud = iopen(trim(vcoudfile),0,-1,0)
C           rewind ifvcoud
            read(ifvcoud) i
            read(ifvcoud) qvv
          endif                 ! q_comm%ID == 0
#ifdef USE_MPI
          call MPI_Bcast(i,1,MPI_INTEGER,0,q_comm%communicator,ierror)
          call MPI_Bcast(qvv,3,MPI_REAL8,0,q_comm%communicator,ierror)
#endif
          if (i /= ngb) call rx('hx0fp0: ngb is not consistent with file '//vcoudfile)
          if (sum(abs(qvv-q))>1d-10 .and. iq<=nqibz) call rx( 'hx0fp0: q is not consistent with file '//vcoudfile)

C         Allocation needs to be done on all nodes
          if (allocated(zcousq)) deallocate( zcousq,vcousq )
          allocate(zcousq(ngb,ngb),vcousq(ngb))
          if (q_comm%ID == 0) then
            read(ifvcoud) vcousq
            read(ifvcoud) zcousq
          endif                 ! qcomm%ID == 0
#ifdef USE_MPI
!         if (q_comm%ID == 0) write(0,'(" [",I4,"] ngb = ",I6,"l ngb = ",I6)') rank ,ngb,ngb
          call MPI_Bcast(vcousq,ngb,MPI_REAL8,0,q_comm%communicator,ierror)
          call MPI_Bcast(zcousq,ngb*ngb,MPI_Double_Complex,0,q_comm%communicator,ierror)
#endif
          idummy = iclose(trim(vcoudfile))
          vcousq = sqrt(vcousq)
          if (allocated(zzr)) deallocate(zzr)
          allocate(zzr(1,1))    !dummy
          zzr = 0d0
        endif ! (newaniso2.and.(.not.chipm))

        if (localfieldcorrectionG .and. nolfco) call rx('inconsistent 1s digit GWversion and no local fields')
        if (chipm .and. nolfco) then ! nmbas has already been allocated as the No. of magnetic sites
!         nmbas = nmbas
        elseif (nolfco) then ! Just one basis function
          nmbas = 1
        else
          nmbas = ngb
        endif
        nmbas1 = nmbas

        if (chipm.and.nolfco) then

C         Even if these arrays are not accessed, we need to allocate
C         them to prevent the debug mode from failing.

          if (allocated(ppovlz)) deallocate(ppovlz)
          if (allocated(ppovlzinv)) deallocate(ppovlzinv)
          if (allocated(ppovl)) deallocate(ppovl)
          allocate(ppovl(0,0),ppovlz(0,0),   ppovlzinv(0,0))
C       elseif(newaniso2) then
!       Takao Apr 2019 has this
C        elseif(nolfco .and. nmbas1==1) then !for <e^iqr|x0|e^iqr>
C          if(allocated(ppovlz)) deallocate(ppovlz)
C          if(allocated(ppovl)) deallocate(ppovl)
C          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb)) !,   ppovlzinv(ngb,ngb))
C          call readppovl0(q,ngc,ppovl)
C          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
C          ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
C          write(6,*)'nnnnn',nbloch+ngc,ngb
        else
          if (allocated(ppovlz)) deallocate(ppovlz)
          if (allocated(ppovlzinv)) deallocate(ppovlzinv)
          if (allocated(ppovl)) deallocate(ppovl)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb),   ppovlzinv(ngb,ngb))
          call readppovl0(q,ngc,ppovl)
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
          ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
C         call zprm('ppovlz',ppovlz,ngb,ngb,ngb)
C         call zprm('zcousq',zcousq,ngb,ngb,ngb)
C         mc -qr zcousq -cc -t -qr ppovlz -x -p -v2dia -s-1 -v2dia -+ -px
C         write(*,*)'nbloch+ngc,ngb',nbloch+ngc,ngb
          write(*,*) 'sum(ppovl),sum(ppovlz)',sum(abs(ppovl)),sum(abs(ppovlz))
        endif

C   --- Set up for epsilon mode --
        if (epsmode) then
          iqixc2 = iq - (nqibz+nq0ix)
          read(ifgb0vec,*) qgbin(1:3),ngbin,igc0dummy,dnorm
          if (sum(abs(q)) == 0d0) then
            if (sum(qgbin**2) > 1d-7) call rx('qgbin=0 xxx See hvccfp0')
          elseif (sum(abs(qgbin(1:3)-q)) >1d-8) then
            call rx('qgbin inconsistent')
          endif
          write(6,"(' --- Readin Mix0vec: ',3d13.6,2i5,d18.8)") qgbin(1:3),ngbin
          print *,'ngb ngbin=',ngb,ngbin
          if (ngb/=ngbin) call rx('hx0fp0: ngb/=ngbgin')
          write(6,"(' ngb nwp niw=',3i8)")ngb,nwp,niw
          if (allocated(zzr0)) deallocate(zzr0)
          allocate(gbvec(ngb),zzr0(ngb))
          if ((.not.chipm).and.nolfco) then
            allocate(x0mean(nw_i:nw,1,1) )
            x0mean = 0d0
          endif
          do i = 1, ngb
            read(ifgb0vec,"(4d24.15)") zz1,zz2,zz3,zz4
! /dnorm at 14June2008
! gbvec and zzr0 do not change because we avoid /dnorm in hvccfp0.
! This means that we changed the definition of gbvec and zzr. See hvccfp0
! dcmplx(zz3,zz4) is zzr defined in hvccfp0; this equals <eiqr|B_I>.
! NOTE: dnorm = \sum_I <eiqr|I><I|eiqr> = volume of unit cell (or sum of MT volume if no IPW).
            gbvec(i)= dcmplx(zz1,zz2)/dnorm
            zzr0(i) = dcmplx(zz3,zz4)/dnorm
          enddo
          write(6,"(' normchk=',1p,255d23.15)") sum( dconjg(gbvec)*zzr0(:) )

          if (chipm) then
            ii = 0
            do imb = 1,nmbas
              ibas =  imbas(imb) !dec26 2006 bug: imb--->ibas
              if (ibas>1) i =  sum(nblocha(1:ibas-1))
              eiqrm(imb) =  sum( dconjg(gbvec(1:nbloch))*svec(1:nbloch,imb) )/mmnorm(imb) *dnorm
              write(6,"(' <eiqr|m>/|m| =',255d23.15)") eiqrm(imb)
              if (imbas_s(imb)<-1) eiqrm(imb) =  -eiqrm(imb)
           enddo
          endif

!! zzr is only for chipm.and.nolfco mode
          if (chipm .and.nolfco) then
            allocate(zzr(ngb,nmbas),x0mean(nw_i:nw,nmbas,nmbas))
            x0mean=0d0
            zzr   =0d0
            zzr(1:nbloch,1:nmbas) = svec(1:nbloch,1:nmbas)
          endif
!! ... Open ChiPM*.nolfc_mat
          if (wqt(iq-nqibz)==0d0.and.chipm) then
            ifchipmn_mat=iopen('ChiPM'//charnum4(iqixc2)//'.nlfc.mat',1,3,0)
            write(ifchipmn_mat,"(255i5)") nmbas
            write(ifchipmn_mat,"(255i5)") imbas(1:nmbas)
            write(ifchipmn_mat,"(255e23.15)") momsite(1:nmbas)
            write(ifchipmn_mat,"(255e23.15)")  mmnorm(1:nmbas)
            write(ifchipmn_mat,"(255e23.15)") eiqrm(1:nmbas)!if necessary, fix code to give eiqrm. takaoAug2012
c            write(ifchipmn_mat,"( ' Here was eiqrm: If needed, need to fix hx0fp0')")
            if (.not.nolfco) then
              ifchipm_fmat=iopen('ChiPM'//charnum4(iqixc2)//'.fmat',0,3,0)
              write(ifchipm_fmat) nbloch, natom,nmbas, iqxini,iqxend, nw_i,nw
              write(ifchipm_fmat) imbas(1:nmbas),momsite(1:nmbas),mmnorm(1:nmbas)
              write(ifchipm_fmat) nblocha(1:natom),svec(1:nbloch,1:nmbas)
              write(ifchipm_fmat) zzr0(1:nbloch) !zzr(1:nbloch,1)
            endif
          elseif(wqt(iq-nqibz)==0d0.and.(.not.chipm)) then
!! ... Open EPS* file
            filepsnolfc ='EPS'//charnum4(iqixc2)//'.nlfc.dat'
            ifepsdatnolfc = iopen ( filepsnolfc,1,3,0)
            write(ifepsdatnolfc,"(a)")' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
            if (.not.nolfco) then
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
              ifepsdat = iopen ( fileps,1,3,0)
              write(ifepsdat,"(a)") ' q(1:3)   w(Ry)   eps  epsi --- LFC included. '
            endif
          endif
        endif

C        if (epsmode.and.nolfco) then !iepsmode==202) then
Cc          rcxqmean=0d0
C        else
C          write(6,*) "rcxq alloc ngb nwhis npm ---",ngb,nwhis,npm
C          allocate(rcxq(ngb,ngb,nwhis,npm) )
C        endif

!! takao apr2012
        if (nolfco) then
          if (allocated(rcxq)) deallocate(rcxq)
          if (allocated(zxq)) deallocate(zxq)
          if (allocated(zxqi)) deallocate(zxqi)
          allocate(zxq(nmbas1,nmbas,nw_i:nw),zxqi(nmbas1,nmbas,niw))
        else
          allocate(zw0(ngb,ngb), zxq(ngb,ngb,nw_i:nw), zxqi(ngb,ngb,niw) )
        endif
        allocate(rcxq(nmbas1,nmbas,nwhis,npm))
        zxq=0d0; zxqi=0d0; if (sergeyv) rcxq = 0d0

C   --- For each spin, do ---
        do 1003 is = 1, nspinmx
          write(6,"(' ##### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ##### ')")iq, is, nqibz + nq0i,nspin
          if (debug) print *,' niw nw=',niw,nw
!! ==== spin chi_charge or chi_+- ====
          isf=is
          if (chipm) then
            write(6,*)" chi_+- mode ixc=",ixc
            if (is==1) isf=2
            if (is==2) isf=1
            if (sergeyv) then
              rcxq=0d0
            endif
          endif

!! ---  tetraini block tetra==.true.===============================1ini
!! --- tetrahedron method --- preparation
!!  nbnb    : total number of pairs (ib,jb) with non-zero weight.
!!  ib=n1b  : band index for occ.   1\ge n1b \ge nband+nctot.  !for chipm it can be unocc
!!            Valence-core order(Core index follows valence index).
!!  jb=n2b  : band index for unocc. 1\ge n2b \ge nband
!!  wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
          if (tetra) then
            if (q_comm%ID == 0)
     .        write(6,"(/' Frequency integration by tetrahedron method on real axis : niw nw nwhis=',3i5)") niw,nw,nwhis
!!          ekxx1  for rk, ekxx2 for q+rk See tetwt4
            do  kx = 1, nqbz
              if (q_comm%ID == 0) then
                call readeval(qbz(:,kx),is,ekxx1(1:nband, kx) )
                call readeval(q+qbz(:,kx),isf,ekxx2(1:nband, kx) )
              endif
            enddo
#ifdef USE_MPI
            call MPI_Bcast(ekxx1, nband*nqbz, MPI_REAL8, 0, q_comm%communicator, ierror)
            call MPI_Bcast(ekxx2, nband*nqbz, MPI_REAL8, 0, q_comm%communicator, ierror)
#endif

C       ... Get pairs (n1b n2b) with non-zero tetrahedron weights.
C           The pairs are not dependent on omega in the denominator of the dielectric function.
C           if (q_comm%ID == 0) write(6,"(' ... Initial call to tetwt5 for array dimensions')")
            job = 0
            if (npm==1) then
              ncc=0
            else
              ncc=nctot
            endif
            allocate(demin(nband+nctot,nband+ncc,nqbz,npm),
     .               demax(nband+nctot,nband+ncc,nqbz,npm),
     .               iwgt (nband+nctot,nband+ncc,nqbz,npm) )
            if (nctot==0) then
              deallocate(ecore)
              allocate(ecore(1,2)) ! Dummy to avoid pointer problems
            endif

            ! Dummy allocations to avoid problems with unallocated array pointers
            allocate(ibjb(1,1,1,1),ihw(1,1,1),nhw(1,1,1),jhw(1,1,1),whw(1))
            ! Until tetw5x has been parallelized, restrict its execution to the head of q_comm
            eibzmode = eibz4x0() .and. nmbas1 == nmbas
            if (q_comm%ID == 0) then
              call tetwt5x_dtet4(npm,ncc,
     i          q, ekxx1, ekxx2, qbas,ginv,ef,
     d          ntetf,nqbzw,nband,nqbz,
     i          nctot,ecore(1,is),idtetf,qbzw,ib1bz,
     i          job,
     o          iwgt,nbnb,      !job=0
     o          demin,demax,    !job=0
     i          frhis, nwhis,   ! job=1    not-used
     i          nbnbx,ibjb,nhwtot, ! job=1    not-used
     i          ihw,nhw,jhw,    ! job=1    not-used
     o          whw,            ! job=1    not-used
     i          iq,is,isf,nqibz,eibzmode,nwgt(:,iq),
     i          nbmx,ebmx,mtet)
!             write(0,'(" [",I4,"](",I3,",",I3,") tetwt5x_dtet4 initialization done.")') rank, q_comm%group_ID, q_comm%ID
              deallocate(ibjb,ihw,jhw,nhw,whw)

!             write(0,'(" [",I4,"](",I3,",",I3,") nband = ",I5)') rank, q_comm%group_ID, q_comm%ID, nband
!             write(0,'(" [",I4,"](",I3,",",I3,") nctot = ",I5)') rank, q_comm%group_ID, q_comm%ID, nctot
!             write(0,'(" [",I4,"](",I3,",",I3,") ncc   = ",I5)') rank, q_comm%group_ID, q_comm%ID, ncc
!             write(0,'(" [",I4,"](",I3,",",I3,") nqbz  = ",I5)') rank, q_comm%group_ID, q_comm%ID, nqbz
!             write(0,'(" [",I4,"](",I3,",",I3,") npm   = ",I5)') rank, q_comm%group_ID, q_comm%ID, npm

            endif ! q_comm%ID == 0
#ifdef USE_MPI
            call MPI_Bcast(nbnb,nqbz*npm,MPI_INTEGER,0,q_comm%communicator, ierror)
            call MPI_Bcast(iwgt,(nband+nctot)*(nband+ncc)*nqbz,MPI_LOGICAL,0,q_comm%communicator, ierror)
            call MPI_Bcast(demin,(nband+nctot)*(nband+ncc)*nqbz,MPI_INTEGER,0,q_comm%communicator, ierror)
            call MPI_Bcast(demax,(nband+nctot)*(nband+ncc)*nqbz,MPI_INTEGER,0,q_comm%communicator, ierror)
#endif

            nbnbx = maxval(nbnb(1:nqbz,1:npm)) !nbnbx = nbnbxx
            allocate(n1b(nbnbx,nqbz,npm),n2b(nbnbx,nqbz,npm))
            do jpm=1,npm
              call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx,
     o          n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
            enddo

C            n1b = -999; n2b = -999
C            do jpm=1,npm
C              call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx,
C     o          n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
C            enddo
C
C             call alignn12(0,nbnbx,nqbz,npm,nband,nctot,ncc,nbnb,n1b,n2b)
C             stop
C            call alignn12(1,nbnbx,nqbz,npm,nband,nctot,ncc,nbnb,n1b,n2b)
C
C            print *, 'n1b'
C            print *, n1b
C            print *, 'n2b'
C            print *, n2b


C       ... Sanity check
            if (sum(abs(nbnb-nbnbtt))/=0)then
              do ik=1,nqbz
                write(6,*)
                write(6,*) "nbnb  =",nbnb(ik,:)
                write(6,*) "nbnbtt=",nbnbtt(ik,:)
              enddo
              call rx( 'hx0fp0:sum(nbnb-nbnbtt)/=0')
            endif

            noccxv = maxval(noccxvx)
            noccx  = nctot + noccxv
            if (q_comm%ID == 0) write(*,322) nctot,noccxv
  322       format(' number of core levels:',i5/
     .            ' max number of occupied states:',i5)
            deallocate(iwgt)
          endif ! tetra initialization

c sergeyv=T only now feb2012
c$$$      call getkeyvalue("GWinput","TmpWWK",tmpwwk,default=.false.)
c$$$      if (tmpwwk)  ifwwk = iopen('tmp.wwk',0,-1,0)
c$$$C --- TetrahedronWeight_4 block. tettwt4  ixc<=5 ==============2ini
c$$$          if (tetra.and.(.not.sergeyv)) then !(ixc==101.or.ixc<=5.or.ixc==22.or.ixc==23)) then
c$$$             allocate( wgt(nband+nctot,nband,nqbz) )
c$$$             if (.not.tmpwwk) allocate( wwk(nbnbx,nqbz,niw+nw+1)) ! +1 is for nwp case
c$$$             if (tmpwwk) rewind(ifwwk)
c$$$             nwmax = niw +nw +1
c$$$             if (ixc==4.or.ixc==5) nwmax = nw+1
c$$$             do iw = 1, nwmax
c$$$                if (ixc==4.or.ixc==5.or.ixc==6) then
c$$$                   fff = freq_r(iw-1)+img*delta
c$$$                   write(6,"(' ### iw on real= ',i3,' omega(Ry)=',2d15.6)")iw-1,fff*2d0
c$$$                elseif(iw<=niw) then
c$$$                   fff = img * freq_i(iw) ! along img axis
c$$$                   write(6,"(' ### iw on imag= ',i3,' omega(Ry)=',2d15.6)")iw,fff*2d0
c$$$                endif
c$$$             enddo
c$$$
c$$$             do 1010 iw = 1,nwmax
c$$$                if (ixc==4.or.ixc==5.or.ixc==6) then
c$$$                   fff = freq_r(iw-1)+img*delta
c$$$!     write(6,"(' ### iw on real= ',i3,' omega(Ry)=',2d15.6)")
c$$$!     &          iw-1,fff*2d0
c$$$                else
c$$$                   if (iw<=niw) then
c$$$                      fff = img * freq_i(iw) ! along img axis
c$$$!     write(6,"(' ### iw on imag= ',i3,' omega(Ry)=',2d15.6)")
c$$$!     &           iw,fff*2d0
c$$$                   else
c$$$!     nw iw-niw ---> iw-niw-1
c$$$                      fff = freq_r(iw-niw-1) + img*delta ! delta is in a.u.  ! along real axis
c$$$!     write(6,"(' ### iw on real= ',i3,' omega(Ry)=',2d15.6)") iw,fff*2d0
c$$$                   endif
c$$$                endif
c$$$
c$$$                call tetwt4x_dev(fff,
c$$$     i               q, ekxx1, ekxx2, qbas,ginv,ef,
c$$$     d               ntetf, nqbzw, nband,nqbz,
c$$$     i               nctot,ecore(1,is),idtetf,qbzw,ib1bz, chipm, !ixc,
c$$$     o               wgt, nbnbtt,
c$$$     i       iq,is,isf,nqibz) ! new input for devided-tetrahedron method.
c$$$
c$$$!-------------------
c$$$
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$C ... Be careful about sign of spin direction.
c$$$c This may be confusing.
c$$$c Maybe better to make it clean with clearer notation.
c$$$              if (chipm.and.spinmom<0) wgt= dconjg(wgt)
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$              if (dreal(fff)==0d0) wgt = dreal(wgt) ! imag part zero clear for omg=0 !this is not necessary.
c$$$
c$$$              if (tmpwwk) then
c$$$                print *,' keyx rewind ifwwk nwmax=',iw,nbnbx,nqbz
c$$$                call rsvwwk4_w( wgt,nqbz,nband,nctot,nbnbx,
c$$$     i          n1b, n2b, nbnb, ifwwk)
c$$$              else
c$$$                call rsvwwk4( wgt,nqbz,nband,nctot,nbnbx,
c$$$     i          n1b, n2b, nbnb,
c$$$     o          wwk(1,1,iw) )
c$$$              endif
c$$$ccccccccccccccccccccccccccccc
c$$$c             write(996,"(' --- iq iw=',2i3)" ),iq,iw
c$$$c             do ik  = 1,nqbz
c$$$c             do ibib= 1,nbnb(ik)
c$$$c               write(996,"(4i4,2d16.8)" )
c$$$c     &         ik,iw,n1b(ibib,ik),n2b(ibib,ik),wwk(ibib,ik,iw)
c$$$c             enddo
c$$$c             enddo
c$$$cccccccccccccccccccccccccccc
c$$$c--------------------------
c$$$c              if (maxval(nbnb(1:nqbz,iw)) > nbnbx) then ! security check
c$$$c               print *,' hx0fp0: maxval(nbnb(1:nqbz,iw)) > nbnbx ',
c$$$c     &          nbnbx,nbnb(1:nqbz,iw)
c$$$c                stop ' hx0fp0: maxval(nbnb(1:nqbz,iw)) > nbnbx '
c$$$c              endif
c$$$c--------------------------
c$$$ccccccccccccccccccccccccccccccc
c$$$c            if (iq==2.and.iw==nwmax) write(996,*)' --- end of iq=2 ---'
c$$$c            if (iq==3.and.iw==1) then
c$$$c             write(996,"(' --- iq iw=',2i3)" ),iq,iw
c$$$c             do ik  = 1,nqbz
c$$$c             do ibib= 1,nbnb(ik)
c$$$c               write(996,"(' k ib jb=',3i4,' wwk=',2d16.8)" )
c$$$c     &         ik,n1b(ibib,ik),n2b(ibib,ik),wwk(ibib,ik,iw)
c$$$c             enddo
c$$$c             enddo
c$$$c             stop '----- test end ------------'
c$$$c           endif
c$$$ccccccccccccccccccccccccccccccc
c$$$ 1010       continue
c$$$            print *,'end of 1010 loop'
c$$$            deallocate(wgt)
c$$$            if (tmpwwk) then
c$$$              allocate( wwk(nbnbx,nqbz,nwmax) )
c$$$              print *,' wwk allocated'
c$$$              ifwwk=iclose('tmp.wwk')
c$$$c              rewind(ifwwk)
c$$$              ifwwk = iopen('tmp.wwk',0,-1,0)
c$$$              print *,' keyx rewind ifwwk nwmax=',nwmax,nbnbx,nqbz,niw,nw
c$$$              do iw=1,nwmax
c$$$                read(ifwwk) wwk(1:nbnbx,1:nqbz,iw)
c$$$              enddo
c$$$              ifwwk=iclose('tmp.wwk')
c$$$              print *,' readend of wwk'
c$$$            endif
c$$$          endif
c$$$C --- TetrahedronWeight_4 block end.========================2end
c$$$
c$$$
c$$$c==============x0kf_v2 block original(.not.sergeyv) case =====================3ini
c$$$          if (tetra.and.((.not.sergeyv).and.ixc/=4))then !ixc<=3.or.ixc==5.or.ixc==101.or.ixc==22.or.ixc==23)) then
c$$$            if (ixc==5) then
c$$$              wwk = dimag(wwk)
c$$$            endif
c$$$            print *,' --- goto x0kf_v2z'
c$$$            if (npm==2) stop 'x0kf_v2z may not work for npm==2'
c$$$            call x0kf_v2z(   ! lessl,
c$$$     i            wwk,n1b,n2b,nbnbx,nbnb, tetra, ! for tetra=T
c$$$     i              q, is,isf, nspin, !is,  !,ifcphi ifev1,ifev2 ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
c$$$     i              qbas,ginv,ef,qbz,wbz,   !iindxk,
c$$$cc     i            index_qbz,n_index_qbz, jan2004
c$$$     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, !mdim=nblocha
c$$$     i              ppb(1,is), !,pdb(1,is),dpb(1,is),ddb(1,is),
c$$$     i              freq_r,freq_i,delta,realomega,imagomega,
c$$$     i              icore,ncore,
c$$$     &              ecore(1,is),
c$$$     d            nlmto,nqbz,noccx,noccxv,
c$$$     d            nl,nclass,natom,nnc,
c$$$     d            nlnmx,mdimx,nbloch,nwp,  niw,
c$$$c     i    iq,ngpn,ngcni(iqx),ngpmx,ngcmx  , geigB(1,1,1,is), Jan2004
c$$$c     i    ngvecpB,ngveccB, ngveccB(1,1,iqx),
c$$$     i     iq,ngc,ngpmx,ngcmx  ,
c$$$c     i    ngveccB,  !ngveccB(1,1,iqx),
c$$$     i     nqbze,nband,nqibz,
c$$$     o     zxq,zxqi,            ! zxq... are accumulated variable for spins
c$$$     i     nolfco,zzr,    !for iepsmode
c$$$     o     x0mean, !   )  !for iepsmode
c$$$     i     vcoul,ngb ) ! for te
c$$$          endif
c$$$          if (epsmode.and.nolfco)
c$$$     &      print *,' nolfco --- end of x0kf_v2 sum x0mean',sum(abs(x0mean))
c$$$c=============x0kf_v2 block end===============================3end

C     --- Tetrahedron weights ---
C        sergeyv=T only now feb2012
C        if (ixc==4.or.sergeyv) then
!         Histogram secstions are specified by frhis(1:nwp)
!         The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!         The last  bin  is     [frhis(nw), frhis(nwp)].
!         nwp=nw+1; frhis(1)=0
          if (frhis(1)/=0d0) call rx('hx0fp0: require frhis(1)=0d0')

C     ... Indices for tetrahedron weights
          allocate(ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm))
          jhwtot = 1
          do jpm = 1, npm
            do ik = 1,nqbz
              do ibib = 1, nbnb(ik,jpm)
                call hisrange( frhis, nwhis,
     i               demin(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     i               demax(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     o               ihw(ibib,ik,jpm),nhw(ibib,ik,jpm))
!               write(0,'(" [",I4,"] iq = ",I3,", jhw(",I3,",",I3,",",I3,") = ",I6)') rank, iq, ibib, ik, jpm, jhwtot
                jhw(ibib,ik,jpm)= jhwtot
                jhwtot = jhwtot + nhw(ibib,ik,jpm)
              enddo ! ibib
            enddo ! ik
          enddo ! jpm
          nhwtot = jhwtot-1
          if (q_comm%ID == 0) write(*,323) nqbz,nbnbx,nhwtot
  323     format(' number of points in full BZ:',i5/
     .           ' max number of nonzero (ib,jb) pairs at any k:',i7/
     .           ' total number of histogram bins for all k:',i11)
          deallocate(demin,demax)
          allocate(whw(nhwtot),ibjb(nctot+nband,nband+ncc,nqbz,npm))

C     ... Populate ibjb for call to tetwt5x_dtet4
          whw = 0d0; ibjb = 0
          do jpm = 1, npm
            do ik  = 1, nqbz
              do ibib = 1,nbnb(ik,jpm)
                ib1 = n1b(ibib,ik,jpm)
                ib2 = n2b(ibib,ik,jpm)
                ibjb(ib1,ib2,ik,jpm) = ibib
              enddo
            enddo
          enddo

C     ... Generate the histogram weights whw
          job = 1
C         if (q_comm%ID == 0) write(6,*) 'begin making histogram weights'
          allocate(iwgt(1,1,1,1),demin(1,1,1,1),demax(1,1,1,1)) ! Dummies to avoid undefined pointers
          call tetwt5x_dtet4(npm,ncc,
     i            q,ekxx1,ekxx2,qbas,ginv,ef,
     d            ntetf,nqbzw,nband,nqbz,
     i            nctot,ecore(1,is),idtetf,qbzw,ib1bz,
     i            job,
     o            iwgt,nbnb,         ! job=0
     o            demin,demax,       ! job=0
     i            frhis,nwhis,       ! job=1
     i            nbnbx,ibjb,nhwtot, ! job=1
     i            ihw,nhw,jhw,       ! job=1
     o            whw,               ! job=1
     i            iq,is,isf,nqibz,eibzmode,nwgt(:,iq),
     i          nbmx,ebmx,mtet)
          if (q_comm%ID == 0) write(*,"(' completed histogram weights'$)")
          call cputid(0)
          print *, ' '
          deallocate(iwgt,demin,demax)
C         TetrahedronWeight_5 block end

C  --- x0k is the main routine to accumulate Im x0 ---
C      This routine is parallelized over k (further parallelization in progress)
C          symmetrize = (is==nspinmx)
C          symmetrize=.false.
C          call x0kf_v4hz(npm,ncc,
C     i           ihw,nhw,jhw,whw,nhwtot,
C     i           n1b,n2b,nbnbx,nbnb,
C     i           q,
C     i           nspin,is,isf,symmetrize,
C     i           qbas,ginv,qbz,wbz,
C     i           nblocha,nlnm,nlnmv,nlnmc,iclass,
C     i           ppb(1,is),
C     i           icore,ncore,
C     d           nlmto,nqbz,noccx,noccxv,
C     d           nl,nclass,natom,nnc,
C     d           nlnmx,mdimx,nbloch,nwhis,
C     i           iq,ngb,ngc,ngpmx,ngcmx,
C     i           nqbze,nband,nqibz,
C     o           rcxq,
C     i           nolfco,zzr,nmbas,ppovlz,ppovl,zcousq,
C     i           chipm,eibzmode,
C     i           nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq),ngrp,eibzsym(:,:,iq))

C         Reorder n1b, n2b so that core states precede valence states and is consistent
C         with ordering of the matrix elements of product and eigenfunction basis
          call alignn12(0,nbnbx,nqbz,0,npm,nband,nctot,ncc,nbnb,n1b,n2b)

C         Calculate the bare susceptibility rcxq
          call x0kr(npm,ihw,nhw,jhw,whw,nhwtot,n1b,n2b,nbnbx,nbnb,q,is,isf,qbas,qbz,
     i      nblocha,nlnmv,nlnmc,iclass,ppb(1,is),icore,ncore,nlmto,nqbz,noccx,noccxv,
     i      nl,nclass,natom,nnc,nlnmx,mdimx,nbloch,nwhis,iq,ngb,ngc,nband,
     i      nolfco,zzr,nmbas1,nmbas,ppovlz,ppovl,chipm,eibzmode,nwgt(1,iq),ebasis,
     o      rcxq)

C         Restore n1b,n2b to its original (valence,core) precedence
C         call alignn12(0,nbnbx,nqbz,0,npm,nband,nctot,ncc,nbnb,n1b,n2b)

C         Charge susceptibility in spin pol case will have is=2 and nspinmx=2 on 2nd spin pass
C         Spin susceptibility always has is=1 and nspinmx=2. Condition never satisfied.
          if (is==nspinmx .and. eibzmode .and. nmbas1 == ngb .and. nmbas == ngb) then
            if (ebasis /= 0) then
              call rcxq2Ebas(1,nmbas,nmbas,nwhis,npm,ngb,ppovlz,rcxq)
            endif
            call x0k_sym(npm,q,ginv,nbloch,nwhis,ngb,ngc,nolfco,zzr,nmbas,zcousq,chipm,ngrp,eibzsym(:,:,iq),rcxq)
          elseif (is==nspinmx .and. nmbas1 == ngb .and. nmbas == ngb) then
            call x0k_sym(npm,q,ginv,nbloch,nwhis,ngb,ngc,nolfco,zzr,nmbas,zcousq,chipm,1,[0d0],rcxq)
          endif

          write(*,'(" I am ",2I4,":  end of x0kr sum rcxq(",I4,",",I4,")=",2F15.5)')
     &         q_comm%ID, rank, iq, is, sum(abs(rcxq))

          if (sergeyv) deallocate(ihw,nhw,jhw,whw,ibjb)
          if (tetra.and.(.not.sergeyv)) deallocate(wwk) ! !ixc<=5.or.ixc==101.or.ixc==22.or.ixc==23))
          if (tetra) deallocate(n1b,n2b)
          iecut = 1

C     --- Hilbert transform.  Generate Real part from Imaginary part ---
          if (sergeyv .and. is==nspinmx.or.chipm) then
            if (q_comm%ID == 0)
     .        write(6,'(" ... begin dpsion5: generate Re chi from Im chi. nwhis nw_i niw nw_w nmbas1 nmbas=",6i5)')
     .        nwhis,nw_i,nw,niw,nmbas1,nmbas
!           Note: dpsion5 overwrites rcxq
            call dpsion5(frhis,nwhis,freq_r,nw,freq_i,niw,realomega,imagomega,rcxq,npm,nw_i,nmbas1,nmbas,
     o        zxq,zxqi,
     i        chipm,schi,is,ecut(iecut),ecuts(iecut))
            if (nolfco .and. epsmode) then
              do iw = nw_i, nw
                x0mean(iw,:,:)=zxq(:,:,iw)
              enddo
            endif
            if (q_comm%ID == 0) print *,'... end dpsion5 checksum W-v',sum(abs(zxq)),sum(abs(zxqi))
          endif                 ! sergeyv
 1003   continue                !end of spin loop
        if (allocated(rcxq)) deallocate(rcxq)

c$$$- 13July2005
c$$$    Takao needed to remove do 2005 loop so as to put dpsion2 into do 1003
c$$$    loop so as to make chipm mode work ok (do dpsion2 for eachs spin).
c$$$    If necessary we may need to recover do 2005 mode---this was for e_rpa test...
c$$$C... this loop is to calculate Ec for different ecut, ecuts.
c$$$        do 2005 iecut=1,necut !In usual mode, necut=1
c$$$c-----------------
c$$$        if (sergeyv) then !This was in 1003 loop as if (ixc==11.and.is==nspinmx) then
c$$$          if (debug) print *, ' --- goto dpsion2'
c$$$          call dpsion2(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega,
c$$$     i      ngb, rcxq,
c$$$     o      zxq, zxqi, ! zxq... are Not the accumulating variable for spins
c$$$     i      nolfco,chipm,is,  rcxqmean,ecut(iecut),ecuts(iecut), !iepsmode,
c$$$     o      x0mean)
c$$$          print *,' --- end of dpsion2'
c$$$          if (epsmode) print *,' end of dpsion sum x0mean',sum(abs(x0mean))
c$$$c          do iw=1,nwp
c$$$c            write(6,"(a,i5,2d13.5)")' sumcheck zxq=',iw,sum(abs(zxq(:,:,iw)))
c$$$c          enddo
c$$$        endif
c       write(ifdpin)  nwhis,nw,niw, ngb, frhis,freq_r, freq_i,rcxq

C  ---  RealOmega ---
        if (realomega) then
          if (chipm) then !ixc==22.or.ixc==23) then
            if (nspin==1) call rx( 'chipm modes are for nspin==2')
            if (anfexist) then  ! Not tested yet
!          if (anfexist) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
              write(6,*) 'I have not tested anfexist=t and ixc==22 23 '//
     &      '--- Maybe x0mean = 2d0* x0mean*2 will OK.'
              call rx( 'hx0fp0: tail output! Not yet for anfexist&ixc==22 23')
            endif
          elseif (epsmode.and.nolfco.and.(.not.chipm)) then !ixc==2.or.iepsmode==202) then
            if (nspin==1) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
            if (nspin==1) zxq = 2d0*zxq       !if paramagnetic, multiply x0 by 2
            if (anfexist) then
! Not tested yet
!          if (anfexist) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
              write(6,*) 'I have not tested anfexist=t and ixc==2 '//
     &     '--- Maybe x0mean = 2d0* x0mean*2 will OK.'
              call rx( 'hx0fp0: tail output! Not yet for anfexist&ixc==2')
            endif
          else
            if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
            if (anfexist) then  !Antiferro acceleration.
              do iw  = nw_i,nw  !max
                call anfx0k(natom,nclass,nblocha,iclass,pos,nbloch,ngc,
     i            q, ngveccB,qbas, ! for q+G
     i            anfvec,ibasf, ! these are antiferro informations.
     .            zxq(1,1,iw))  ! i/o
              enddo
            endif
          endif

c         write (ifxd,"(1x,3f10.4)") q(1),q(2),q(3)
c         write (ifrx) rxq,cxq
          if (epsmode) then
            if (nolfco) then
              ttt='without LFC'
            else
              ttt='with LFC'
            endif
            if (chipm) then
              write(6,*) '--- chi0_{+-}}^{-1}      --- '//ttt
            else
              write(6,*) '--- dielectric constant --- '//ttt
              write(6, *)" trace check for W-V"
            endif
          endif

C     ... prepare for iq0.
          iq0 = iq - nqibz
C          if (newaniso2 .and. ixc==11) then
c            if ( iq==1) then
c              ifepstinv = iopen('EPS0inv',0,-1,0)
c              write(ifepstinv) ngb
c            elseif (iq0>0) then
c              write(ifisk) ngb,nw_i,nw,niw
c              allocate( sk( ngb),  sks(ngb), w_k(ngb), w_ks(ngb) )
c              allocate( skI( ngb),  sksI(ngb), w_kI(ngb), w_ksI(ngb) )
cc            allocate( s_vc(ngb-1))
c            endif
cccc           if (iq<=nqibz) then
ccccc           endif
c            if (iq0 >= 1) then
c              ifepstinv = iopen('EPS0inv',0,0,0)
c              read(ifepstinv) ngb
c            endif
c             allocate(sk(ngb),skI(ngb))
C          endif
C          if (newaniso2) then
          if (allocated(epstilde)) deallocate(epstilde,epstinv)
          allocate(epstilde(ngb,ngb),epstinv(ngb,ngb))
C          endif

C         if (q_comm%ID == 0 .and. newaniso2.and.iq<=nqibz) then
          if (q_comm%ID == 0 .and. iq<=nqibz) then
            write(6,"(' === trace check for W-V, q =',3f12.6,' iq =',i3,' ===')") q,iq
            write(6,"(4x,' freq    iw',14x,'W-v(site)',25x,'W-v(tot)',24x,'diff')")
          endif
C         if (q_comm%ID == 0 .and. newaniso2.and.iq>nqibz) then
          if (q_comm%ID == 0 .and. iq>nqibz) then
            write(6,"(' === Contribution to eps(q=0), q =',3f12.6,' iq =',i3,' ===')") q,iq
            write(6,"(4x,' freq    iw',14x,'eps(LFC)',25x,'eps(noLFC)',23x,'diff')")
          endif

C     --- iw loop for real axis ---
          do  1015  iw = nw_i, nw
#ifdef USE_MPI
            writer1 = mod(iw,q_comm%size)
            recd = iw-nw_i+1
            if (q_comm%ID == writer1) ifrcw  = iopen('WVR.'//charnum5(iq)//'.'//charnum5(recd),0,-1,mrecl)
#endif
            frr = dsign(freq_r(abs(iw)),dble(iw))
            if (.not.epsmode) then
C             if (newaniso2.and.iq<=nqibz) then
              if (iq<=nqibz) then
              if (iq==1) then
                ix=1
                zw0(:,1)=0d0
                zw0(1,:)=0d0
              else
                ix=0
              endif

!!            Eqs.(37),(38) in PRB81 125102
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                  if (igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                enddo
              enddo
              epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
              if (ix == 0) then
                call matcinv(ngb,epstinv)
              else
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
              endif

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cmmm direct inversion vs. block inversion
c$$$                if (iq>nqibz) then
c$$$c direct inversion
c$$$                ix=0
c$$$                do igb1=ix+1,ngb
c$$$                do igb2=ix+1,ngb
c$$$                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$                  if (igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$                enddo
c$$$                enddo
c$$$                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$                do igb1=1+ix,ngb
c$$$                do igb2=1+ix,ngb
c$$$                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
c$$$                  if (igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
c$$$                enddo
c$$$                enddo
c$$$c              write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$              write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,1d0/epstinv(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$c             write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(1,2:10:3),zw0(1,63:70:3)
c$$$c block inversion
c$$$                ix=1
c$$$                do igb1=ix+1,ngb
c$$$                do igb2=ix+1,ngb
c$$$                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$                  if (igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$                enddo
c$$$                enddo
c$$$                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$                absq=sqrt(sum(q**2*tpioa**2))
c$$$                sk(  1:ngb)= zxq(1,1:ngb,iw)
c$$$                sks( 1:ngb)= zxq(1:ngb,1,iw)
c$$$                w_k(1) =0d0
c$$$                w_ks(1)=0d0
c$$$                w_k( 2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(vcousq(1)*sk(2:ngb)*vcousq(2:ngb),epstinv(2:ngb,2:ngb))
c$$$                w_ks(2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(epstinv(2:ngb,2:ngb),vcousq(1)*sks(2:ngb)*vcousq(2:ngb))
c$$$                llw(iw,iq0)=
c$$$     &             1d0
c$$$     &            -vcousq(1)*sk(1)*vcousq(1) ! sk(1,1,iw)=sks(1,1,iw)=H of Eq.(40).
c$$$     &            -vcousq(1)*vcousq(1)* sum( vcousq(2:ngb)*sk(2:ngb) * matmul(epstinv(2:ngb,2:ngb),sks(2:ngb)*vcousq(2:ngb)))
c$$$                write(*,"('mmmmzwp99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0), !(1d0/llw(iw,iq0)-1d0)*vcousq(1)**2,
c$$$c     &                  w_k(2:10:3)/llw(iw,iq0), w_k(63:70:3)/llw(iw,iq0)
c$$$     &                  w_ks(2:10:3)/llw(iw,iq0), w_ks(63:70:3)/llw(iw,iq0)
c$$$                write(*,"('mmmmzwp99x ')")
c$$$                endif
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

              do igb1=1+ix,ngb
                do igb2=1+ix,ngb
                  zw0(igb1,igb2) = vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                  if (igb1==igb2) zw0(igb1,igb2) = zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                enddo
              enddo
              zw(1:ngb,1:ngb) = zw0 ! copy because zw = zw(nblochpmx,nblochpmx)
#ifdef USE_MPI
              if (q_comm%ID == writer1) then
                if (q_comm%ID == writer1) write(ifrcw, rec= 1 ) zw !  WP = vsc-v
                if (q_comm%ID == writer1) write(*,'("[",I4,"] iq = ",I3," iw = ",I3," zw(1,1) = ",2F15.5)') rank, iq, iw, zw(1,1)
                call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
              endif
#else
              write(ifrcw,rec=iw-nw_i+1) zw !  WP = vsc-v
              call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
#endif
            endif

!           Takao Apr 2019 has localfieldcorrectionG always true ?
C           if (newaniso2.and.iq>nqibz) then
            if (iq>nqibz) then
              vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
              if (localfieldcorrectionG) then ! Full inversion to calculate eps with LFC
                ix=0
                do igb1=ix+1,ngb
                  do igb2=ix+1,ngb
                    if (igb1==1.and.igb2==1) then
                      epstilde(igb1,igb2) = 1d0 - vcou1*zxq(1,1,iw)
                      cycle
                    endif
                    epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                    if (igb1==igb2) then
                      epstilde(igb1,igb2) = 1d0 + epstilde(igb1,igb2)
                    endif
                  enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                if (ix == 0) then
                  call matcinv(ngb,epstinv)
                else
                  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                endif
                llw(iw,iq0)= 1d0/epstinv(1,1)
              else
                llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw)
              endif
              if (q_comm%ID == 0) write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8,0p,2x,2f12.6)') frr,iw,
     .          llw(iw,iq0),1d0-vcou1*zxq(1,1,iw),llw(iw,iq0)-1d0+vcou1*zxq(1,1,iw)
            endif

C              if (.not.newaniso2) then               ! Original mode
C                call rx('newaniso2 no longer supported')
Cc               call wcf( ngb, vcoul, zxq(1,1,iw), 1, zw0)  ! W= (1-v zxq)^{-1} v
C              endif

C              write(6,"('--- EigenValues for Im(W) ---')")
C              allocate(ebb(ngb))
C              call diagcvh2( (zw0-transpose(dconjg(zw0)))/2d0/img, ngb,ebb)
C              do ii=1,ngb
C                if ( abs(ebb(ii))>1d-8 .and. ebb(ii)>0) then
C                  write(6, "('hhhIWq : iw ii eb=',2i4,d13.5)") iw, ii,ebb(ii)
C                else
C                  write(6, "('hhhIWqxxx : iw ii eb=',2i4,d13.5)") iw,ii, ebb(ii)
C                endif
C              enddo
C              deallocate(ebb)
#ifdef USE_MPI
            if (q_comm%ID == writer1) ifrcw = iclose('WVR.'//charnum5(iq)//'.'//charnum5(recd))
#endif

!       ... Make dielectric function
            elseif (epsmode.and.(.not.chipm)) then !.and.ixc/=23) then ! No LFC (local field correction). It's better to use echo 4| hbasfp0.
C            if (newaniso2) then
!! there is difference of two vcmean below since we use (slightly) screened Coulomb (screenfac() in switch.F)
!! NOTE that we use vcoul with screening (screenfac() is used in hvccfp0.F
c             vcmean = fourpi/sum(q**2*tpioa**2) !aug2012
              vcmean = vcousq(1)**2
              epsi(iw,iqixc2)= 1d0/(1d0 - vcmean*zxq(1,1,iw))
              write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2d23.15,3x, 2d23.15,
     &        " vcmean x0mean =", 2d23.15,3x, 2d23.15)') iqixc2,iw,2*frr,
     &          1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean, zxq(1,1,iw) !x0mean(iw,1,1)
              write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)')
     &          q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              if (.not.nolfco) then
                ix = 0
                do igb1 = ix+1,ngb
                do igb2 = ix+1,ngb
! ! DMT: the branching is pointless for vcmean=vcousq(1)**2 (which is the one used) but i'm leaving it for the commented out case see the note 15 lines above
                  if (igb1==1.and.igb2==1) then
                    epstilde(igb1,igb2) = -vcmean*zxq(igb1,igb2,iw) !aug2012
                  else
                    epstilde(igb1,igb2) = -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                  endif
                  if (igb1==igb2) epstilde(igb1,igb2) = 1+epstilde(igb1,igb2)
                enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb) = epstilde(ix+1:ngb,ix+1:ngb)
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                epsi(iw,iqixc2) = epstinv(1,1)
                write(6,'( " iq iw omega eps epsi  wLFC="
     &         ,2i6,f8.3,2d23.15,3x, 2d23.15)')
     &          iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                write(6,*)
                write(ifepsdat,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)')
     &          q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              endif
C              else
C                 print *,'not support epsmode and newaniso=F mode now.'
Cc$$$                vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
Cc$$$                if (.not.nolfco) then
Cc$$$                 x0mean(iw,1,1) = sum( dconjg(zzr(:,1))* matmul(zxq(:,:,iw),zzr(:,1)))
Cc$$$                endif
Cc$$$                epsi(iw,iqixc2) = 1d0/(1- vcmean * x0mean(iw,1,1))
Cc$$$                write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2d23.15,3x, 2d23.15,
Cc$$$     &          " vcmean x0mean =", 2d23.15,3x, 2d23.15)') iqixc2,iw,2*frr,
Cc$$$     &          1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean,x0mean(iw,1,1)
Cc$$$                write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)')
Cc$$$     &          q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
Cc$$$                if (.not.nolfco) then ! With LFC
Cc$$$                  imode=2
Cc$$$                  call wcf( ngb, vcoul, zxq(1,1,iw), imode,
Cc$$$     &                zw0) !  write(6,"('ssschk1=',3d13.5)") sum(abs(zw0)) sum(abs(gbvec))
Cc$$$                  epsi(iw,iqixc2)= sum( dconjg(gbvec) * matmul(zw0,zzr(:,1)) )
Cc$$$                  write(6,'( " iq iw omega eps epsi  wLFC="
Cc$$$     &            ,2i6,f8.3,2d23.15,3x, 2d23.15)')
Cc$$$     &            iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
Cc$$$                  write(6,*)
Cc$$$                  write(ifepsdat,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)')
Cc$$$     &            q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
Cc$$$                endif
C              endif

!       ... Make transverse susceptibility
            elseif (epsmode.and.chipm) then
              allocate( x0meanx(nmbas,nmbas) )
              if (nolfco) then  ! ChiPM mode without LFC
                if (legas) then
                   call rx('hx0fp0.m.F need to implement homogeneous gas here if required')
! --- three lines below may work for test purpose for legas. But not sure.
c                  vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
c                  write(ifchipmn,'(3f12.8,2x,f8.5,2x,2e23.15)')
c     &            q, 2*schi*frr, 1d0-vcmean*2*x0mean(iw,1,1)  !4*pi*alat**2/sum(q**2)/4d0/pi**2*x0mean(iw)
                else
                  x0meanx = x0mean(iw,:,:)/2d0  !in Ry unit.
                endif
C         ... ChiPM mode with LFC... NoLFC part
              else
                zxq(1:ngb,1:ngb,iw) = zxq(1:ngb,1:ngb,iw)/2d0  ! in Ry.
                do imb1=1,nmbas
                  do imb2=1,nmbas
                    x0meanx(imb1,imb2)=
     &              sum( svec(1:nbloch,imb1)*
     &              matmul(zxq(1:nbloch,1:nbloch,iw),svec(1:nbloch,imb2)))  !/ mmnorm**2  I removed mmnorm may2007
                  enddo
                enddo
!       x0meanx= <m|chi^+-(\omega)|m>/<m|m>**2
              endif
              do imb1=1,nmbas
                do imb2=1,nmbas
                  x0meanx(imb1,imb2) = x0meanx(imb1,imb2)/mmnorm(imb1)/mmnorm(imb2)
                enddo
              enddo
              write(ifchipmn_mat,'(3f12.8,2x,f20.15,2x,5000d23.15)') q, 2*schi*frr, x0meanx(:,:)
              if (.not.nolfco) write(ifchipm_fmat) q, 2*schi*frr, zxq(1:nbloch,1:nbloch,iw)
              deallocate(x0meanx)
            endif
#ifdef USE_MPI
            if (q_comm%ID == writer1) ifrcw = iclose('WVR.'//charnum5(iq)//'.'//charnum5(recd))
#endif
 1015     continue !iw
#ifdef USE_MPI
          if (iq==1) then !MLsend notification that file for iq=1 has been written:
            i_send = 1
            if (first_q_group /= last_q_group) then
              call MPI_Isend(i_send, 1, MPI_INTEGER, last_q_group, 1, q_inter%communicator, i_request2, ierror)
            endif
!           write(0,'("[",I4,"] Finished writing file for iq=1:",2I10)') rank, q_comm%group_ID, q_inter%ID
          endif
#endif
          if (allocated(zzr)) deallocate(zzr)
          if (allocated(x0mean)) deallocate(x0mean)
          if (allocated(gbvec)) deallocate(gbvec)
        endif                   ! end of RealOmega block

C   ... Close files for epsmode
        if (epsmode) then !iepsmode/=0) then      ! only calculate iq>nqibz
          if (chipm) then
            ifchipmn_mat=iclose('ChiPM'//charnum4(iqixc2)//'.nlfc.mat')
            if (.not.nolfco) then
              ifchipm_fmat=iclose( 'ChiPM'//charnum4(iqixc2)//'.fmat')
            endif
          else
            filepsnolfc ='EPS'//charnum4(iqixc2)//'.nolfc.dat'
            ifepsdatnolfc = iclose( filepsnolfc)
            if (.not.nolfco) then
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
              ifepsdat  = iclose(fileps)
            endif
          endif
        endif

C   --- Imaginary Omega ---
        if (imagomega) then
          print *,' ... start imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
c     ... Antiferro acceleration.
          if (anfexist) then
            do  iw  = 1, niw
              call anfx0k(natom,nclass,nblocha,iclass,pos,nbloch,ngc,
     i          q, ngveccB,qbas, ! for q+G
     i          anfvec,ibasf,    ! these are antiferro informations.
     m          zxqi(1,1,iw))    ! i/o
            enddo
          endif
          if (ecorr_on>0)then !ixc==101.or.(sergeyv.and.imagonly)) then
            imode=0
          else
            imode=1
          endif

          do  1016 iw = 1, niw  ! loop on the imaginary axis
            if (iq<=nqibz) then ! Eqs.(37),(38) in PRB81 125102
              if (iq==1) then
                ix = 1
                zw0(:,1) = 0d0
                zw0(1,:) = 0d0
              else
                ix = 0
              endif
              do igb1 = ix+1, ngb
                do igb2 = ix+1, ngb
                  epstilde(igb1,igb2) = -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                  if (igb1==igb2) epstilde(igb1,igb2) = 1+epstilde(igb1,igb2)
                enddo
              enddo
              epstinv = epstilde
              if (ix == 0) then
                call matcinv(ngb,epstinv)
              else
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
              endif
              do igb1 = ix+1,ngb
                do igb2 = ix+1,ngb
                  zw0(igb1,igb2) = vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                  if (igb1==igb2) zw0(igb1,igb2) = zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                enddo
              enddo
              zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
              if (q_comm%ID == writer2) then
                write(ifrcwi, rec=iw)  zw !  WP = vsc-v
                call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
              endif

            else ! iq>nqibz
              vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2 = fourpi/sum(q**2*tpioa**2-eee)
              if (localfieldcorrectionG) then ! Full inversion to calculate eps with LFC.
                ix = 0
                do igb1 = ix+1,ngb
                  do igb2 = ix+1,ngb
                    if (igb1==1.and.igb2==1) then
                      epstilde(igb1,igb2) = 1d0 - vcou1*zxqi(1,1,iw)
                      cycle
                    endif
                    epstilde(igb1,igb2) = -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                    if (igb1==igb2) then
                      epstilde(igb1,igb2) = 1d0 + epstilde(igb1,igb2)
                    endif
                  enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb) = epstilde(ix+1:ngb,ix+1:ngb)
                if (ix == 0) then
                  call matcinv(ngb,epstinv)
                else
                  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                endif
                llwI(iw,iq0) = 1d0/epstinv(1,1)
              else  ! llwI without LFC. LFC contribution is added in
                llwI(iw,iq0) = 1d0 - vcou1*zxqi(1,1,iw) !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
              endif
              if (q_comm%ID == 0) write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8,0p,2x,2f12.6)')
     .                 frr,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw),llwI(iw,iq0)-1d0+vcou1*zxqi(1,1,iw)
            endif

C              if (.not.newaniso2) then ! original mode
C                call rx( 'not checked here')
Cc               call wcf( ngb, vcoul,zxqi(1,1,iw),imode,  zw0)
C              endif

C              if (newaniso2.and.iq>nqibz) then
Cc               zw(1:ngb,1:ngb) = 0d0 ! zw(nblochpmx,nblochpmx)
Cc               write(ifrcwi, rec=(iq-iqxini)*niw + iw)  zw    !  WP = vsc-v
C              else
C                zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
Cc               write(ifrcwi, rec=(iq-iqxini)*niw + iw)  zw    !  WP = vsc-v
C                write(ifrcwi, rec=iw)  zw    !  WP = vsc-v
C                call tr_chkwrite("freq_i iq iw  imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
C              endif

!! --- Miyake's total energy branch !Nov2004. not maintained now... need to fix this maybe(2012takao) -----------
C            if (ecorr_on>0 .and. (.not.newaniso2)) then !I did not modified this for newaniso2 2012takao
C              if (debug) write(*,*)'ip,ix=',iq,iw,'  niw=',niw
C              call getwk(iq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)
C              call ecorq2 (vcoul, zw0, ngb, iq,iw,ieceig,
C     o        erpaqw, trpvqw, trlogqw)
CC --- Integration along imaginary axis.
C! omit k and basis index for simplicity
C! wint = -(i/4pi) < [w'=-inf,inf] Q(w') >
C!
C! When w' ==> iw', w' is now real,
C!   wint =  (1/2pi) < [w'=0,inf] Q(iw') >
C!
C! transform: x = 1/(1+w')
C! this leads to a denser mesh in w' around 0 for equal mesh x
C! which is desirable since Q is peaked around w'=0
C!    wint =  (1/2pi) < [x=0,1] Q(iw') / x^2 >
C              faca   =  wk4ec* wiw(iw)
C              trpv(iecut)    = trpv(iecut)  + faca* trpvqw
C              trlog(iecut)   = trlog(iecut) + faca* trlogqw
C              totexc(iecut)  = totexc(iecut)+ faca* erpaqw !  = trpv+ trlog
Cc            ecqw(iq,iw) = erpaqw
C              if (iw==1) then
C                write(ieclog,*)
C              endif
C              if (iw==1.and.iq==iqxini) then
C                write(ieclog,
C     &         "('   iq   iw   omega/i(Ry)        IntWgt',
C     &         '    trpvqw(eV)      ecqw(eV)   ecqw*IntWgt',
C     &         ' :  ecut   ecuts')")
C              endif
C              write(ieclog,"( 2i5,3f14.6,3f14.6,2f8.3)")
C     &       iq,iw, 2d0*freq_i(iw), faca, trpvqw*hartree, erpaqw*hartree,
C     &       faca*erpaqw*hartree, ecut(iecut),ecuts(iecut)
C              close(ieclog)
C              open(ieclog,file="ecorr.chk",access='append')
Cccccccccccccccccccccccccccccccccccccccccccccccccc
Cc           allocate( ovlpc(ngb,ngb),evall(ngb),
Cc     &          evecc(ngb,ngb))
Cc           evall=0d0
Cc           ovlpc=0d0
Cc           do i=1,ngb
Cc             ovlpc(i,i)=1d0
Cc           enddo
Cc           nmx=ngb
Ccc1           call diagcv(ovlpc,zw0/2d0+transpose(dconjg(zw0))/2d0,evecc,ngb, evall,nmx,1d99, nev)
Cc           call diagcv(ovlpc,zw0,evecc,ngb, evall,nmx,1d99, nev)
Cc           write(6,"('ngb nev=',2i5)") ngb,nev
Cc           write(6,"('chk eigen of zw0 Max Min=',2d13.6)")maxval(evall),minval(evall)
Cc           do i=1,3
Cc             write(6,*) i, evall(i)
Cc           enddo
Cc           do i=ngb-3,ngb
Cc            write(6,*) i, evall(i)
Cc           enddo
Cc           deallocate( ovlpc,evall,evecc)
Cccccccccccccccccccccccccccccccccccccccccccccccc
C            endif
 1016     continue              ! loop over imaginary frequencies

C         if (newaniso2) then
c            if (iq==1)    ifepstinv = iclose('EPS0inv') !iq==1 close write mode.
c            if (iq>nqibz) ifepstinv = iclose('EPS0inv') !iq>nqibz --> readin mode
            deallocate(epstinv)
            if (allocated(epstilde)) deallocate(epstilde)
C         endif
        endif                   ! Imag omega

        if (allocated(vcoul)) deallocate(vcoul)
        if (allocated(zw0)) deallocate(zw0)
        iw = idalloc('zxq',iverbose+3,nmbas1*nmbas,(niw-nw_i+1+niw)*2)
        if (allocated(zxq)) deallocate(zxq)
        if (allocated(zxqi)) deallocate(zxqi)

        if (q_comm%ID == writer2 .and. (ixc==101.or.normalm)) then
          ifrcwi = iclose('WVI.'//charnum5(iq))
        endif
        if (normalm) then
#ifdef USE_MPI
          write(*,'(" [",I4,"](",I3,",",I3,") End of q-loop: iq = ",I3," Time = ",F15.5)')
     .      rank, q_comm%group_ID, q_comm%ID,iq,(dwtime()-walltimeq)/60
#else
          ifrcw = iclose('WVR.'//charnum5(iq))
          write(*,*) "End of q-loop: iq = ", iq
#endif
        endif
        if (iq == iqxend) last_q = .true.

 1001 continue                  ! end of big iq loop
C!! === Receive llw and llwI at node 0, where q=0(iq=1) is calculated. ===
C      if(MPI__size/=1) then
C        do iq=nqibz+1,iqxend
C          iq0 = iq - nqibz
C          if(MPI__ranktab(iq)/=0) then !jan2012
C            if(MPI__ranktab(iq) == MPI__rank) then
C              dest=0
C              call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
C              call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
C            elseif(MPI__root) then
C              src=MPI__ranktab(iq)
C              call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
C              call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
C            endif
C          endif
C        enddo  
C      endif

#ifdef USE_MPI
!     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      write(0,'(" [",I4,"](",I3,",",I3,") Still alive after q-loop:",L5)') rank, q_comm%group_ID, q_comm%ID, last_q
#endif

C --- W(0) divergent part and non-analytic constant part ---
C     Special treatment of the coulomb singularity at q=0 (iq=1)
C     The singularity goes as 1/q^2, but it is integrable resulting in finite G*W
C     This code is based on Eq. 33 in J. Phys. Soc. Jpn. 83, 094711 (2014)
C     Coefficients w_L are defined in Eq. 28 for F_L, Eq. 39.
C     w_L and the other contents of EPSwklm are generated by qg4gw (see mkqg.F)
      if (last_q .and. q_comm%ID == 0 .and. ixc == 11) then
        print *, '==== Divergent part of W(q=0) ===='

!!      Return w0 and w0i = diagonal element at Gamma point and llmat = L(omega=0) matrix.
!!      Outputs w0,w0i,llmat
        allocate(w0(nw_i:nw),w0i(niw))
        call w0w0i(q_comm%ID,ixc.ne.1011,llw,llwI,nw_i,nw,nq0i,niw,freq_r,freq_i,q0i,w0,w0i,llmat)

#ifdef USE_MPI
!ML     this MPI_Wait call should only finish when iq=1 has send the i_ready
        if (first_q_group /= last_q_group) then
          write(0,'("[",I4,"] ",I4,":  Waiting for message from ",I5)') rank, q_inter%ID, first_q_group
          call MPI_Wait(i_request,i_status,ierror)
        endif
#endif

!      w0,w0i are stored in zw for q=0; w_ks*wk are stored in zw for iq>nqibz
!      Assume iq=1 belongs to rank=0
        do  iq = 1, 1           !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
          q = qibze(:,iq)
          do  ircw = 1, 2       ! Real and imaginary frequencies
            if (ircw==1) then
              nini=nw_i; nend=nw
#ifdef USE_MPI
              ifrcwx = iopen('WVR.'//charnum5(iq)//'.'//charnum5(iw-nini+1),0,-1,mrecl)
              read(ifrcwx, rec = 1 ) zw !(1:ngb,1:ngb)
#else
              ifrcwx = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
#endif
            elseif (ircw==2) then
              nini=1; nend=niw
              ifrcwx = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
            endif
            do iw=nini,nend
c             if (iq<=nqibz) read(ifrcwx, rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
              read(ifrcwx, rec= iw-nini+1) zw !(1:ngb,1:ngb)
              if (ircw==1) zw(1,1) = w0(iw)
              if (ircw==2) zw(1,1) = w0i(iw)
              write(ifrcwx,rec=iw-nini+1) zw !(1:ngb,1:ngb)
            enddo
            if (ircw==1) then
              ifrcwx = iclose('WVR.'//charnum5(iq))
            elseif (ircw==2) then
              ifrcwx = iclose('WVI.'//charnum5(iq))
            endif
          enddo ! Real, imaginary branches
        enddo
      endif ! last_q
      endif ! Faleev's normal mode

!! --- legas mode is not working now. Need fixing... voltot ntot are not given.
C      if (epsmode.and.legas) then
C        call rx( ' LEGAS mode is not maintained well. Need some fixing.')
C        voltot=0d0
C        ntot=0d0
C        print *,' Find LEGAS. legas =',legas
C        iflegas = 2101
C        open (iflegas,file='LEGAS')
C        read(iflegas,*)rs
C        close(iflegas)
C        alpha  = (9*pi/4d0)**(1d0/3d0)
C        qfermi = alpha/rs
C        efx  = qfermi**2
C        valn = efx**1.5d0*voltot/3d0/pi**2
C        write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
C        write (6,*)'     Exact Fermi momentum  qf  =', qfermi
C        write (6,*)'     Exact Fermi energy    Ef  =', efx
C        do iq = iqxini,iqxend ! q=(0,0,0) is omitted!
C          if (iq<=nqibz) cycle
C          print *,' iq=',iq
C          iqixc2 = iq- (nqibz+nq0ix)
C          filele ='EPSEG'//charnum4(iqixc2)//'.dat'
C          ife = iopen ( filele,1,3,0)
C          write(ife,"(a)")
C     &          ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
C          q = qibze(:,iq)
C          qt= sqrt(sum(qibze(1:,iq)**2))*2d0*pi/alat
C          qs= qt/qfermi
C          write(6,"(' qs qfermi=',2d13.5)"    ) qs,qfermi
C          write(6,"(' q-q^2/2 q+q^2=',2d13.5)") qs-qs**2/2d0,qs+qs**2/2d0
C          do iw  = nw_i,nw
C            ww  = freq_r(iw)
C            muu = ww/qfermi**2
C            if (qs<2d0 .and. muu < qs-qs**2/2d0) then
C              x0mx= -img*qfermi/(4*pi*qs)*2*muu
C            elseif ( qs<2d0 .and. muu < qs+qs**2/2d0) then
C              x0mx= -img*qfermi/(4*pi*qs)*( 1d0-(muu/qs-.5d0*qs)**2 )
C            else
C              x0mx=0d0
C            endif
C            vcmmmm= 4*pi/qt**2
C            epsi(iw,iqixc2) = 1d0/(1- vcmmmm * x0mx)
Cc            epsi(iw,iqixc2) = 1d0/(1- vcmmm(iq) * x0meanx)
C            write(ife,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)')
C     &        q, 2*ww,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
C          enddo
C        enddo
C        print *,' ----------legas end--------'
C      endif

C--- Write TEECOR ecorr_on mode
      if (imagomega.and.ecorr_on>0) then
        hartree=2d0*rydberg()
        ifcor   = iopen('TEECORR2',1,-1,0) ! output files
        do iecut=1,necut
          write(6,"( ' RPA Ec =' 3f23.15,'   ecut ecuts (Ry)=',2d12.4)")
     &   totexc(iecut)*hartree,trpv(iecut)*hartree, trlog(iecut)*hartree
     &    ,ecut(iecut),ecuts(iecut)
          write(ifcor,*) '============================'
          write(ifcor,*) 'Correlation energy Erpa (eV)'
          write(ifcor,*) '============================'
          write(ifcor,*)' ### '
          write(ifcor,"(5e23.15)")
     &     totexc(iecut)*hartree,trpv(iecut)*hartree,trlog(iecut)*hartree
     &    ,ecut(iecut),ecuts(iecut)
        enddo
C... output ecqw !    write(ifcor,*)'### ecqw(q,w) ###'
        write(ifcor,*)' nqibz =',nqibz
        write(ifcor,*)' nq0i  =',nq0i
        write(ifcor,*)' niw   =',niw
        write(ifcor,*)' --- See details of Ec in ecor.chk ---'
c        nqitot = nqibz + nq0i
c        call wecqw(ifcor,
c     d           nqibz,nqbz,nq0i,nqitot,niw,
c     o           wibz,wqt,wx,freqx,ecqw)
C... Write electron gas correlation energy
        legas = .false.
        INQUIRE (FILE = 'LEGAS', EXIST = legas)
        if (legas) then !!! test for electron gas case.
          call rx( ' LEGAS mode is not maintained well. Need some fixing.')
          voltot=0d0
          ntot=0d0
          print *,' find LEGAS. legas =',legas
          iflegas = 2101
          open (iflegas,file='LEGAS')
          read(iflegas,*)rs
          close(iflegas)
          alpha = (9*pi/4d0)**(1d0/3d0)
          qfermi = alpha/rs
          efx  = qfermi**2
          valn = efx**1.5d0*voltot/3d0/pi**2
          write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
          write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
          write (6,*)' egas  Exact Fermi energy    Ef  =', efx
          if (tetra) call rx( 'legas You have to give ef of  tetrahedron')
          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
          qfermi= dsqrt(efz)
          alpha = (9*pi/4d0)**(1d0/3d0)
          rs    = alpha/qfermi
          write (ifcor,*)' --- electron gas ---'
          write (ifcor,*)' density parameter rs= ', rs
          write (ifcor,*)' kf= ',qfermi
          write (ifcor,*)' ### Barth-Hedin formula'
          ecelgas = eclda_bh(rs) * hartree * ntot
          write (ifcor,*)ecelgas
          write (ifcor,*)' ### Perdew-Zunger formula'
          ecelgas = eclda_pz(rs) * hartree * ntot
          write (ifcor,*)ecelgas
          write (ifcor,*)' ### Gell-Mann and Brueckner formula'
          ecelgas = (-0.0311d0 * dlog(rs) -0.048d0) * hartree * ntot
          write (ifcor,*)ecelgas
        endif
      endif
C...close files
c        close(ifdpin)
      is = iclose('hbe.d')

      is = idalloc(' ',1,1,1)
      write(*,"(/' Max dynamic allocation :',i7,' MB'/)") is
      call cputid(0)
C     call MPI__Finalize

      call mpi_end

      if (nthreads > 1) then
        print 131, (dwtime() - walltime0)/60, cpusec()/60, nthreads
      else
        print 131, (dwtime() - walltime0)/60, cpusec()/60
      endif
  131 format(' hx0fp0 :  Wall clock time',f9.1,'m','  Processor time',f9.1,'m':,i4,' threads')

C      if (ixc==101) stop ' OK! hx0fp0 mode=101 Ecor '
C      if (ixc==1)   stop ' OK! hx0fp0 mode=1 normal'
C      if (ixc==2)   stop ' OK! hx0fp0 mode=2  eps NoLFC'
c      if (ixc==4)   stop ' OK! hx0fp0 ixc=4 test'
c      if (ixc==5)   stop ' OK! hx0fp0 ixc=5 test'
c      if (ixc==6)   stop ' OK! hx0fp0 ixc=6 test'
C      if (ixc==22)  stop ' OK! hx0fp0 mode=22 chi+- NoLFC'
C      if (ixc==23)  stop ' OK! hx0fp0 mode=23 chi+- LFC  '
C      if (ixc==3)   stop ' OK! hx0fp0 mode=3 eps with LFC'
      if (ixc==101) call rx0( ' OK! hx0fp0 mode=101 Ecor ')
      if (ixc==202) call rx0( ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC')
      if (ixc==203) call rx0( ' OK! hx0fp0 mode=203 sergeyv eps LFC ')
      if (ixc==222) call rx0( ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv')
      if (ixc==223) call rx0( ' OK! hx0fp0 mode=223 chi+- LFC sergeyv')
      if (ixc==11)  call rx0( ' OK! hx0fp0 mode=11 normal sergeyv')
      if (ixc==12)  call rx0( ' OK! hx0fp0 mode=12 Ecor sergeyv mode')
      end

c--------------------------------------------------------------------
      real*8 function eclda_bh(rs)
      real(8) :: rs,cp,rp,z
      cp       = 0.0504d0*0.5d0 ! 0.5 changes unit from Ry to Hartree
      rp       = 30.d0
      z        = rs / rp
      eclda_bh = -cp * ( (1.d0+z**3)*dlog(1.d0+1.d0/z)
     .                 + 0.5d0*z - z**2 - 0.33333333d0 )
      end
c--------------------------------------------------------------------
      real*8 function eclda_pz(rs)
      real(8) :: rs
      if (rs >= 1.d0) then
        eclda_pz = -0.1423d0 / (1.d0 + 1.0529d0*dsqrt(rs) + 0.334d0*rs)
      else
        eclda_pz = -0.0480d0 + 0.0311d0*dlog(rs) - 0.0116d0 * rs
     .         + 0.0020d0*rs*dlog(rs)
      endif
      end
c--------------------------------------------------------------------
      subroutine wecqw(ifcor,
     d                 nqibz,nqbz,nq0i,nqitot,niw,
     o                 wibz,wqt,wx,freqx,ecqw)

      implicit double precision (a-h,o-z)
      dimension   wibz(nqibz),wqt(nq0i),wx(niw),
     .            freqx(niw),ecqw(nqitot,niw)
      real(8):: rydberg
      write(ifcor,*)'### ecqw(q,w) ###'
      write(ifcor,*)'nqibz =',nqibz
      write(ifcor,*)'nq0i  =',nq0i
      write(ifcor,*)'niw   =',niw
      do ip = 2,nqitot
        if (ip <= nqibz) then
          wk = wibz(ip)*0.5d0 ! 0.5 for the normalization of wibz
        else
c        wk = wqt(ip-nqibz)*wibz(1)*0.5d0 ! 0.5 for the normalization of wibz
          wk = wqt(ip-nqibz)* 1d0/dble(nqbz)
        endif
        write(ifcor,*)'### iq,wq = ',ip,wk
        sume=0d0
        do ix = 1,niw
          write(ifcor,*)freqx(ix),ecqw(ip,ix),wx(ix)
          sume=sume+  wx(ix)/(freqx(ix)*freqx(ix)) * ecqw(ip,ix)
        enddo
        write(ifcor,*) '  sum ecqw*wx=', wk*sume*2d0*rydberg()
! end of ip-loop
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine getsqovlp(q,ngc,ngb,sqovlp)
!! == Get sqrt of ppovl ==
      implicit none
      real(8)::q(3)
      integer:: ngc,ngb,nbloch,i,nmxx,ix,iy,nev
      complex(8):: sqovlp(ngb,ngb)
      complex(8),allocatable:: ooo(:,:),ppo(:,:),ppovl(:,:)
      complex(8),allocatable:: evec(:,:)
      real(8),allocatable:: eval(:)
      nbloch = ngb-ngc
      if (ngc==0) goto 888

      allocate(ppovl(1:ngc,1:ngc))
      call readppovl0(q,ngc,ppovl)
      allocate(ooo(ngc,ngc),ppo(ngc,ngc),evec(ngc,ngc),eval(ngc))
      ooo= 0d0
      do ix=1,ngc
        ooo(ix,ix)=1d0
      enddo
      ppo = ppovl
      deallocate(ppovl)
      nmxx = ngc
      evec = 0d0
      eval = 0d0
      call diagcv(ooo, ppo,
     &     evec, ngc, eval, nmxx, 1d99, nev)
      print *,' diagcv overlap ngc nev=',ngc,nev
      deallocate(ooo,ppo)
c
 888  continue
      sqovlp=0d0
      do i=1,nbloch
        sqovlp(i,i)=1d0
      enddo
      do i=1,ngc
        if (eval(i)<0d0) then
          call rx( 'getsqovlp:  eval(i) <0d0')
        endif
        do ix=1,ngc;  do iy=1,ngc
          sqovlp(ix+nbloch,iy+nbloch)=
     &     sqovlp(ix+nbloch,iy+nbloch)
     &      + evec(ix,i)* sqrt(eval(i))* dconjg(evec(iy,i))
        enddo ;      enddo
      enddo
      if (allocated(evec)) deallocate(evec)
      if (allocated(eval)) deallocate(eval)
      print *,' end of getsqovlp'
c         sqovlpi = sqovlp
c         call matcinv(ngb,sqovlp)     !  inverse
c         ovlpi=ovlp
c         deallocate(ppovl,ovlp)
      end

      subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
C- Check write for zw
      implicit none
      integer:: nblochpmx,nbloch,ngb,iw,i,iq
      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
      real(8):: freqq
      logical :: smbasis
      character*(*)::tagname
      trwv=0d0
      if (.not.smbasis()) then
        do i = 1,nbloch
          trwv = trwv + zw(i,i)
        enddo
      endif
      trwv2 = 0d0
      do i = 1,ngb
        trwv2 = trwv2 + zw(i,i)
      enddo
      write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8,0p,2x,2f12.6)') freqq,iw,trwv,trwv2,trwv2-trwv
c     write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
c      do i = 1,ngb
c        write(6,'("iii i=",i4,a,f10.4,2i5,4d22.14)')i,tagname,freqq,iq,iw,zw(i,i)
c      enddo
      end

c--------------------------------------------------------------------
C      subroutine test_xxx(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
C      implicit none
C      integer:: nblochpmx,nbloch,ngb,iw,i,iq
C      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
C      real(8):: freqq
C      logical :: smbasis
C      character*(*)::tagname
C      trwv2 = 0d0
C      forall( i = 1:ngb)
C        trwv2 = trwv2 + zw(i,i)
C      end forall
C      end
c--------------------------------------------------------------------

      function matcinvf(a) result(b)
!!== Test routine for Inversion ==
      implicit none
      integer :: info,n,n2(2)
      integer,allocatable :: ipiv(:)
      complex(8):: a(:,:), b(1)
      complex(8),allocatable:: work(:)
      n2= SHAPE(a)
      n=n2(1)
      call zcopy(n,b,1,a,1)
      call zgetrf(n,n,a,n,ipiv,info)
      if (info/=0) then
        print *,' matcinv: zegtrf info=',info
        call rx( ' matcinv: zegtrf ')
      endif
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if (info/=0) then
        print *,'matcinv: zegtri info=',info
        call rx( 'matcinv: zegtri ')
      endif
      end

c--------------------------------------------------------------------
      subroutine diagno00(nbloch,wpvc,eval)
!! == ontain eigenvalue only for input complex matrix wpvc(nbloch,nbloch)
      implicit none
      integer:: nbloch,nmx,nev,i
      complex(8),allocatable:: ovlpc(:,:),evecc(:,:),wpvcc(:,:)
      real(8)::eval(nbloch)
      complex(8):: wpvc(nbloch,nbloch)
      allocate( ovlpc(nbloch,nbloch),evecc(nbloch,nbloch),wpvcc(nbloch,nbloch))
      wpvcc= wpvc
      ovlpc= 0d0
      do i=1,nbloch
        ovlpc(i,i)=1d0
      enddo
      eval=0d0
      nev  = nbloch
      nmx  = nbloch
      call diagcv(ovlpc,wpvcc, evecc, nbloch, eval, nmx, 1d99, nev)
      deallocate(ovlpc,evecc,wpvcc)
      end
