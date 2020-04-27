      program hx0fp0_sc
!! Calculate W-V for QSGW
!! ------------------------------------------------
!! We calculate chi by the following three steps.
!!  tetwt5: tetrahedron weights
!!  x0kf_v4h: Accumulate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
!!  dpsion5: calculate real part by Hilbert transformation from the Im part
Cl Local variables
Cl  nmb     : total number of mixed product basis functions I,J in which <I|v|J> is calculated
Cl          : (formerly ngb)
Cl  nmbas   : number of mixed product basis functions used for susceptibility
Cl          : It may be nmb, but may be smaller in some contexts
Cl          : (formerly nmbas_in, nmbas, nmbas1, nmbas2)
Cl  nmbas1  : nmbas1=nmbas, except when offset-gamma points when local fields are omitted.
Cl          : Then nmbas1=1
Cl  This version requires PPOVL and PPOVL0 for matrix elements.
Cl  No other PPOVL files are needed.
!! ------------------------------------------------
      use m_readqg,only: readngmx,readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval
      use m_read_bzdata,only: read_bzdata,
     &  ngrp2=>ngrp,nqbz,nqibz,nqbzw,ntetf,n1,n2,n3,qbas,ginv,
     &  dq_,qbz,qibz,qbzw,
     &  idtetf,ib1bz
      use m_genallcf_v3,only: genallcf_v3,
     &     nclass,natom,nspin,nl,nn,nnc,ngrp,
     &     nlmto,nlnmx,nctot,niw,nw,
     &     alat,ef,dw,delta,iclass,
     &     il,in,im,nlnm,nlnmv,
     &     nlnmc,icore,ncore,
     &     pos,ecore,symgg ! symgg=w(igrp)

      use keyvalue,only: getkeyvalue
      use m_pbindex,only: PBindex
      use m_readqgcou,only: readqgcou

      use mpi_mod, q_comm => comm_1, q_inter => inter_1, rank => mpi_rank

!     Takao Apr 2019 new module
      use m_readq0p,only: readq0p,
     & wqt,q0i,nq0i ,nq0iadd,ixyz


      implicit none

      integer num_q_groups
      integer q_start, q_step
      logical last_q,modified_response

c     use m_spacegrouprot,only: Spacegrouprot
c     use m_rotMPB,only: rotMPB
Cr
Cr Be careful with the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
Cr      The indexing starts with core first and then valence on top of core
Cr      So n-index in "in" for valence electron is different from "inv".
Cr      Search in and inv.
      integer::iopen,maxocc2,iclose, ixc,iqxini,iqxend,
     &     nw_input, ifhbe,  nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &     i,nq0ix,ngcmx,nqbze,nqibze,ini,ix,ngrpx
     &     ,mdimx,nbloch,nblochpmx,is,nwp, !ifvcfpout,
     &     nxx,ifepsfile,ifepsfileI
     &     ,iw,nwhis,nw2,noccxv,noccx
     &     ,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx
     &     ,irot,iq,iqv,ngb
     &     ,kx,isf,job,nbnbx,nhwtot,noccxvx(2)=-9999,nwmax !,ifev1,ifev2
     &     ,jhwtot,ik,ibib,ib1,ib2
     &     ,ngpmx ,nwmin
      real(8)::  omg_c,omg2, q(3)
      real(8):: ua=1d0          ! this is a dummy.
      integer::  ndble=8
      procedure(integer):: idxk,nword_in_dble
C     Moved to module m_readq0p
C     integer:: nq0i
C     real(8),allocatable::   wqt(:),q0i(:,:)
      integer,allocatable :: ngveccB(:,:),nx(:,:),nblocha(:),lx(:)
      complex(8) :: Pbar
      complex(8),allocatable:: vcoul(:,:),
     &     zw(:,:),zw0(:,:), zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),
     &     symope(:,:),
     &     ppb(:,:), qbze(:,:),qibze(:,:) !,ecore(:,:)
      complex(8),allocatable :: rcxq(:,:,:,:)
      logical :: tetra=.true., crpa=.false.
      integer,allocatable :: n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:)

      logical:: w4pmode=.false.
      complex(8),allocatable:: wmu(:,:),wmuk(:,:)
      integer:: ifw4p,ngbq0,igb
      real(8):: qv(3,3)

      real(8):: ebmx
      integer:: nbmx,mtet(3)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:) !qbzw(:,: )
      real(8) ::anfvec(3)
      logical :: anfexist,debug
      integer,allocatable:: ibasf(:)
      real(8),allocatable ::  freq_r(:),freq_i(:),freqx(:),wx(:),expa(:),frhis(:) !ppovl(:,:)
      logical :: realomega, imagomega
      complex(8),allocatable:: zzr(:,:),ppovl(:,:),ppovlz(:,:)
      logical :: mpi__root=.true.
!     integer :: mpi__rank=0
      character(5) :: charnum5
      real(8) :: Emin, Emax      ,emax2,emin2
      real(8) :: ebasis = 0d0   ! Nonzero => make rcxq in E-basis.  ebasis = truncation
c     integer :: iSigma_en  !sf..21May02  !iSigma_en is integer
                                !parameter stored in GWIN_V2
                                !determines which approximation for mapping self-energy to QSGW.
                                !Self-energy should be made hermitian for energies to be real
                                !iSigma_en==0 SE_nn'(ef)+img integral:delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
                                !iSigma_en==1 SE_nn'(ef)+delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
                                !iSigma_en==2 [SE_nn'((e_n+e_n')/2)+h.c.]/2
                                !iSigma_en==3 [(SE_nn'(e_n)+SE_nn'(e_n'))/2+h.c.]/2
      real(8) :: omg2max
      real(8), allocatable :: freqr2(:)
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      real(8),allocatable :: whw(:)
      integer,allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      integer:: nwin,incwfin, verbose,nbcut,nbcut2
      integer::iverbose,idalloc
      real(8)::efin
      logical :: nolfco=.false.
      integer:: ngc,mrecg
      real(8)::  quu(3),deltaq(3)
      logical,allocatable :: iwgt(:,:,:,:)
      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg,smbasis
      logical   :: timereversal
      integer:: npm ,jpm,ncc,nw_i
      integer:: nw_ixx,nwxx  ! ixc=1011
      real(8)   :: frr !, sciss
      logical :: newaniso,newaniso2 !,z1offd
      integer :: ngb0,ifvcoud,idummy,igb1,igb2,nmbas1,nmbas,iq0
      complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:)
      real(8),allocatable:: vcousq(:)
      real(8):: fourpi,sqfourpi,tpioa,vcou1
!! Eq.(40) in PRB81 125102
c     complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:)
      complex(8),allocatable:: llw(:,:),llwI(:,:),llwmod(:,:),llwImod(:,:)
      integer:: ifrcwx,ircw,nini,nend
C     Outputs generated by w0w0i
      complex(8),allocatable :: w0(:),w0i(:)
      complex(8) :: llmat(3,3)

!     complex(8):: vc1vc2
      integer,allocatable:: nwgt(:,:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
      integer:: l2nl,iqxendx, GWversion
      logical:: eibz4x0,tiii,iprintx,chipm=.false.,iqinit,localfieldcorrectionllw,localfieldcorrectionG
      real(8)::qvv(3),ecut,ecuts,hartree,rydberg,pi
      character(128):: vcoudfile
      integer,parameter :: NULLI=-99999
      integer:: ifw0w0i
      logical :: eibzmode
      real(8):: schi=-9999 !dummy

      complex(8), parameter :: cmplx0 = 0, cmplx1 = 1

      integer nthreads,omppid
      real(8):: cpusec,walltime0,dwtime, walltimeQ

!! -------------------------------------------------------------------
!     For MPI
      logical MPIinitexists
      integer :: writer1, writer2, recd, i_group, first_q_group, last_q_group
      integer :: i_ready, i_request, i_request2, i_status, i_send

! --- First executable statement ---
      call mpi_start
!     write(0,'(" [",I4,"] Starting HX0FP0_SC:")') rank
      walltime0 = dwtime()
C$OMP parallel
      nthreads = omppid(0)
C$OMP end parallel

      allocate( zzr(1,1)) !dummy
      hartree=2d0*rydberg()
      pi = 4d0*datan(1d0)
      fourpi=4d0*pi
      sqfourpi=sqrt(fourpi)
      localfieldcorrectionG = mod(GWversion(),10) /= 0 .or. localfieldcorrectionllw()
      call getkeyvalue("GWinput","mod_x0" ,modified_response, default=.false. )

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

C      write(outfilename,'("out_",I1)') rank
C      open(6,file=trim(outfilename))
C
C      write(outfilename,'("out-",I1)') rank
C      open(0,file=trim(outfilename))

C --- computational mode select
      call getjobmode(ixc)
      if (ixc == NULLI .and. rank == 0) then
         write(6,*) '--- hx0fp0_sc Choose one of the following modes ---'
         write(6,*) '    11    make W'
         write(6,*) '    1011  Add W0W0I part for q=0'
         write(6,*) '    10011 constrained RPA'
         call readin5(ixc,iqxini,iqxend)
       endif
#ifdef USE_MPI
       call MPI_Bcast(ixc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
      if (ixc<0) call rx('Choose computational mode with --job=#')
      if (rank == 0) then
        call headver('hx0fp0_sc',ixc)
C        call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
C        write(6,"(' ScaledGapX0=',f5.3)") sciss

        if (ixc==11) then
          write(6,*) "ixc=11 s.f.'s normal mode "
        elseif(ixc==1011) then
          write(6,*) 'ixc=1011 Add W0W0I part at q=0'
        elseif(ixc==10011) then
          write(6,*) " OK ixc=10011 crpa mode "
          crpa=.true.
        else
          call rx('illegal mode')
        endif
      endif

!! newaniso2 is now fixed to be .true.
      newaniso2=.true.
C      if (newaniso()) then
C        newaniso2=.true.
C      else
C        call rx(' hx0fp0_sc: ony newaniso supported  now')
C      endif

C     Extract ecut, ecuts, ebasis from GWinput
      if (rank == 0) then
        call getkeyvalue("GWinput","ecut_p" ,ecut, default=1d10 )
        call getkeyvalue("GWinput","ecuts_p",ecuts,default=1d10 )
        call getkeyvalue("GWinput","ebasis",ebasis,default=0d0)
      endif
#ifdef USE_MPI
      call MPI_Bcast(ecut,1,MPI_Double_Precision,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(ecuts,1,MPI_Double_Precision,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(ebasis,1,MPI_Double_Precision,0,MPI_COMM_WORLD,ierror)
#endif

      if (ebasis <= -1d0) ebasis = -nint(-ebasis)
      if (rank == 0) then
        write(*,"(1x)")
        if (ebasis <= -1d0) then
          write(*,"( a,i4,a)") ' rcxq generated in E basis, remove',nint(-ebasis),' states'
        elseif (ebasis /= 0d0) then
          write(*,"( a,1p,E11.3)") ' rcxq generated in E basis, Ecut =', ebasis
        endif
        if (.not. timereversal()) write (*,"(' Time reversal suppressed')")
      endif

      call cputid (0)

C --- Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
      call read_BZDATA()

!! Use off-regular mesh (not including Gamma) for qbzreg()=F
!! This must be consistent with qg4gw.F-mkqg.F
      if (.not.qbzreg()) then
        deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
        do i=1,nqbz
          qbz(:,i) = qbz(:,i) - deltaq/2d0
          write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
         enddo
      endif
      if (rank == 0) then
         print *,' nqbz qbz =',nqbz
         print *,' nqibz ngrp=',nqibz,ngrp2
      endif

!! === Readin by genallcf (MPI stuff inside it) ===
      nwin   = 0                !Readin nw from NW file
      incwfin= 0                !use ForX0 for core in GWIN
      efin =  0d0               !readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if (ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      tpioa=2d0*pi/alat

!...  These are allocated and setted by genallcf_v3
!     integer::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
!     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
!     real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
!     character(120):: symgrp
!     character(6),allocatable :: clabl(:)
!     integer,allocatable:: iclass(:)
!     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
!     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
!     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
!     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
!     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
!     &    occv(:,:,:),unoccv(:,:,:)
!     &   ,occc(:,:,:),unoccc(:,:,:),
!     o    nocc(:,:,:),nunocc(:,:,:)
!     real(8), allocatable::
!     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)

      nw_input = nw ; print *, 'nw delta=',nw_input,delta
      debug=.false.; if (verbose()>=100) debug=.true.
      iverbose = (verbose()/30)*10
      if (debug) print *,' end of genallc'


!!!!  WE ASSUME iclass(iatom)= iatom,  nclass = natom.  !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom) call rx( ' nclass /= natom ')

C --- tetra or not
      if (delta <= 0d0) then
         tetra =  .true.
         delta = -delta
        if (rank == 0) print *,' hx0fp0: tetrahedron mode delta=',delta
      else
         tetra = .false.        ! switch for tetrahedron method for dielectric functions
      endif

!! --- read dimensions of h,hb
      if (rank == 0) then
        ifhbe = iopen('hbe.d',1,0,0)
        read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
        if (nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
        if (nqbz /=nqbzt) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')
        ifhbe      = iclose('hbe.d')
      endif
#ifdef USE_MPI
      call MPI_Bcast(nprecb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(mrecb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(mrece,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(nlmtot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(nqbzt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(nband,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(mrecg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

#endif

C --- Readin Offset Gamma ---
      call readq0p()
      write(6,"(' ### nqibz nq0i nq0iadd=', 3i5)") nqibz,nq0i,nq0iadd
c     --- readin by rdpp ; Radial integrals ppbrd and plane wave part
c     call getsrdpp( nclass,nl,
c     o               ngpmx,ngcmx,nxx )

      if (rank == 0) call getsrdpp2(nclass,nl,nxx)
#ifdef USE_MPI
      call MPI_Bcast(nxx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif

C     if (rank == 0) then
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      if (rank == 0) print *,' ngcmx ngpmx=',ngcmx,ngpmx
C     endif
C#ifdef USE_MPI
C      call MPI_Bcast(ngpmx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
C      call MPI_Bcast(ngcmx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
C#endif
      nqbze  = nqbz *(1 + nq0i+nq0iadd)
      nqibze = nqibz + nq0i+nq0iadd
      allocate( qbze(3, nqbze), qibze(3, nqibze))
      call dcopy(3*nqbz, qbz,  1, qbze,1)
      call dcopy(3*nqibz,qibz, 1, qibze,1)
      do i = 1,nq0i+nq0iadd
        qibze(:,nqibz+i)  = q0i(:,i)
        ini = nqbz*(1 + i -1)
        do ix=1,nqbz
          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
        enddo
      enddo
      ngrpx = 1                 ! dummy ngrpx=1
      l2nl=2*(nl-1)
!     Takao Apr 2019 reorganized allocations
      allocate(symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
!     write(0,'(" [",I4,"] nl = ",I5,", ngrpx = ",I5)') rank, nl, ngrpx
      allocate( nx(0:l2nl,nclass),  nblocha(nclass) ,lx(nclass),
     &     ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:l2nl,nxx, nspin*nclass),
     &     cgr(nl**2,nl**2,(2*nl-1)**2,ngrpx))


      call rdpp_v3(nxx, nl,ngrpx, nn, nclass, nspin,symope,qbas,
     o     nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
!--- Product basis information is now read

      if (.not.smbasis()) nblochpmx = nbloch + ngcmx
      allocate(ngveccB(3,ngcmx))
      iqxend = nqibz + nq0i
      if (rank == 0)  print *, ' nqibz nqibze=',nqibz,nqibze

!!... initialization of readEigen !readin m_hamindex
!     MPI considerations: Both calls populate a number of module data elements.
!     It seems best to call them on all nodes and handle
!     single node file access from within.
      call init_readeigen(ginv,nspin,nband,mrece) ! EVU EVD are readin in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)

c     --- ecore ---
c     allocate(ecore(nctot,nspin)) !core energies
c     do  is = 1,nspin
c     if (nctot > 0) then
c     call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c     print *,' ecore is=',is,ecore(:,is)
c     endif
c     enddo

c     --- set realomega, imagomega tetra nw niw nwp ifgb0vec ---------------------
!     nwp, freq_r,  frhis(1:nwhis+1)
c     if   ( ixc==1 ) then !old imagw = 2 case
c     realomega =.true.
c     imagomega =.true.
c     stop 'hsfp0sc: ixc==1 is not implemented'
ccccccccccccccccccccfaleev 21May02,  use only ixc=1,11 modes cccccccccc
c     elseif( ixc==2.or.ixc==3 ) then
c     realomega =.true.
c     imagomega =.false.
c     niw = 0
c     ifepscond = 2102
c     open (ifepscond,file='EPScond')
c     read (ifepscond,*) epsrng, dwry !epsrng dw in Ry
c     dw = dwry/2d0
c     close(ifepscond)
c     if (dw==0d0) then
c     nw = 1
c     else
c     nw = (epsrng/2d0 - 1d-10)/(dw/2d0) + 2 !epsrng/2d0 corresponds to in a.u.
c     endif
c     allocate(epsi(nw,neps))
c     if (paralellx0) then
c     ifgb0vec = iopen ( "Mix0vec."
c     &     //xxt(iqxini,iqxend),1,3,0)
c     c          "//charnum5(iqxini)//'to'//charnum5(iqxend),1,3,0)
c     else
c     ifgb0vec = iopen ( "Mix0vec",1,3,0)
c     endif
c     elseif(ixc==4.or.ixc==5.or.ixc==6) then
c     ! ... These are test modes.
c     ! ixc=4 tetrahedron weight test. tetwt5.vs.tetwt5. Write tethis.chk
c     ! ixc=5 Spectrum function (Img part) along the Real axis with tetwt4
c     ! ixc=6 Spectrum function (Img part) along the Real axis with tetwt5. Histgram method.
c     realomega = .true.
c     imagomega = .false.
c     tetra     = .true.
c     niw = 0
c     ! ---  For tetwt5 ---  the tetrahedron weight for spectrum function (imaginary part)
c     !   Histogram bins are specified by freq_r(1:nwp)
c     !     nwp=nw+1; frhis(1)=0
c     !     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
c     !     The last  bin  is     [frhis(nw), frhis(nwp)].
c
c     ! ... These parameters specifies a test histgram bins;Sergey's mesh just for test modes.
c     nw0 = 200    !100    800
c     dwh = 0.01d0 !0.02d0 0.0025d0 !in hartree
c     ! ...
c     call findemaxmin(ifev,nband,nqbz,nspin,emax,emin)
c     if (nctot > 0) Emin = minval(ecore)
c     omg2max = (Emax-Emin)*.5+.2d0    !(in Hartree) covers  all relevant omega, +.2 for margin
c     omg1max = dwh*(nw0-1)
c     nwp = int(sqrt(omg2max*(2*nw0-1d0)/dwh-(nw0**2-3*nw0+1d0)))+1 ! + 1 for margin
c     nw  = nwp-1
c     write(*,*) Emax,Emin,nw0,nw  ! nwp is new max number in frequency array
c     write(6,'(a32,2i7,2d15.3)')'hx0fp1: nw0,nw,omg1max,omg2max='
c     &              , nw0,nw,  omg1max,omg2max
c     if (nw <= nw0) stop 'hx0fp0:ixc==[456] nw2 <= nw'
c     allocate(freq_r(nwp))
c     do iw=1,nwp  !This is a test mesh by Sergey.Faleev
c     if (iw<=nw0) then;  freq_r(iw)=dwh*(iw-1)
c     else;  freq_r(iw)=dwh*(iw**2+nw0**2-3*nw0+1)/(2*nw0-1d0)
c     endif
c     enddo !freq_r(iw) is linear for iw<=nw and quadratic for nw<iw<=nw2
c     !freq_r(iw) chosen in such a way that it is continues with
c!!!  nw nwp=nw+1 freq_r(1:nwp) are used after here.
c     allocate(frhis(nwp))
c     frhis=freq_r(1:nwp)
c     nwhis=nw

!! We first accumulate Imaginary parts. Then it is K-K transformed to obtain real part.
      realomega = .true.
      imagomega = .true.
      tetra     = .true.
      if (rank == 0) call getkeyvalue("GWinput","omg_c",omg_c)
#ifdef USE_MPI
      call MPI_Bcast(omg_c,1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
#endif

      if (rank == 0) write(*,"(/' Parameters for real frequency mesh: dw, omg_c',2f12.6)") dw, omg_c
      call findemaxmin(nband,qbze,nqbze,nspin,emax,emin)
      if(.not.qbzreg()) then
         allocate(qbz2(3,nqbz))
         do iq=1,nqbz
            qbz2(:,iq)=qbz(:,iq)+dq_
         enddo
         call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
         emax=max(emax,emax2)
         emin=min(emin,emin2)
         deallocate(qbz2)
      endif
!     Takao Apr 2019 packaged into getwemax, getfreq
      if (nctot > 0) emin = minval(ecore(:,1:nspin))
      omg2max = (emax-emin)*.5d0+.2d0 !(in Hartree) covers  all relevant omega, +.2 for margin
      if (rank == 0) write(*,"(' emin emax',2f12.6,'  range',f12.6)") emin, emax, omg2max*2
      nw2 = int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
      allocate(freqr2(nw2))     !+1 b/c (iw-1)
      do  iw = 1, nw2           !linear + quadratic term
        freqr2(iw) = dw*(iw-1) + dw**2/2/omg_c*(iw-1)**2
      enddo
      if (nw2 < 2) call rx( "hx0fp0: nw2 < 2")
      if (dw*(nw-2) > freqr2(nw2-1)) call rx("hx0fp0: dw*(nw-2) > freqr2(nw2-1)")
      nwhis = nw2-1
      allocate(frhis(1:nwhis+1))
      frhis = freqr2(1:nwhis+1)
      if (rank == 0) write(*,"(' number of histogram bins =',i4/)") nwhis

      nw = nw2-1                ! nw+1 is how many points of real omega we use
                                ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
                                ! maximum nw=nw2-1 because nwhis=nw2-1
      do iw = 3, nw2-1
                                !nw is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
                                !Here frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
                                !nw was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
                                ! for correlation Coulomb Wc(omg),
                                ! and one more point omg=dw*(nw-1) needed for extrapolation.
                                ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
                                ! and one more point omg=frhis_m(nw) needed for extropolation
                                ! used in subroutine alagr3z in  sxcf.f.
        omg2 = (freqr2(iw-2)+freqr2(iw-1))/2d0
        if (omg2 > dw*(nw_input-3)) then
          nw=iw
          exit                  ! 'nw_input' is only used to get maximum frequency for
                                !  dressed coulomb line
        endif
      enddo
      allocate(freq_r(0:nw))
      freq_r(0)=0d0
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
      if (.not.timereversal()) then
        npm = 2
        nw_i = -nw
        if (.not.tetra) call rx('require tetra=T for timereversal=off')
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

C --- Imaginary frequencies
C     Points x for Gaussing quadrature between (0,1) and w=(1-x)/x
      if (imagomega) then
        allocate(freq_i(niw),freqx(niw),wx(niw),expa(niw) )
        call freq01 (niw,ua,freqx,freq_i,wx,expa)
        deallocate(freqx,wx,expa)
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
      noccxv = maxocc2 (nspin,ef, nband, qbze,nqbze)
      if (noccxv>nband) call rx( 'hx0fp0: all the bands filled! too large Ef')
      noccx = noccxv + nctot
      nprecx = ndble            !We use double precision arrays only.

C ... ppb
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass,nspin) )
C     ppb = 0d0

C ... This is just to get nblochpmx
      if (smbasis()) then
        call getngbpomat(nqibz+nq0i,nblochpmx,i)
      endif
      mrecl = nprecx*2*nblochpmx*nblochpmx/nword_in_dble()

C     Open WV.d and allocate zw
      if (rank == 0) ifwd = iopen('WV.d',1,-1,0)
      if (rank == 0) write (ifwd,"(1x,i3,i10,2i6,3i5)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz+nq0i-1,nw_i

      allocate(zw(nblochpmx,nblochpmx))
      nspinmx = nspin

!! --- READIN ANFcond. we don't use this path so much ---
!     Takao Apr 2019 discarded
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

C     This branch doesn't work now ... comment out until it works again
C      if (anfexist) then
C         allocate(ibasf(natom),transaf(1:3,natom))
C         ifianf = 211
C         open(ifianf,file='ANFcond')
C         read(ifianf,*) soflag
C         if (soflag=='SO') then  !SO case is not implemented
C            if (rank == 0) print *,' ANFcond found. but SO case is not implemented. '
C     &           ,'Not using ANF symm'
C            anfexist=.false.
C            goto 889
C         endif
C         read(ifianf,*)
C         read(ifianf,*) anfvec(1:3)
C         if (rank == 0) write(6,'(" ANFcond=",3f13.6)') anfvec
C         do ibas=1,natom
C            read(ifianf,*) ibas1, ibasf(ibas)
C            if (ibas1/=ibas) call rx( ' readin ANFcond error:ibas1/=ibas')
C            if (rank == 0) write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
Cc     &     ,transaf(1:3,ibas)
C         enddo
C         close(ifianf)
C         nspinmx =1
C      endif
C 889  continue
C      if (smbasis() .and. anfexist) then
C         print *,' smbasis=T & anfexist=T is not yet. '//
C     &        ' At leaset anfx0k shold be changed'
C         call rx( ' smbasis=T & anfexist=T is not yet. ')
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
      iqxini=1 !for newaniso
      eibzmode = eibz4x0()

C --- ppb ---
      do  is = 1,nspin          ! <Phi(R,L,n) Phi(R,L',n') B(R,i)> for each spin
        irot = 1
!       Parallelized with OMP and can run on all nodes: no I/O
!     Takao Apr 2019 allocates ppb for all rotated functions (array ppbir)
        call ppbafp_v2(irot,ngrpx,is,nspin,
     .    il,in,im, nlnm,
     .    nl,nn,nclass,nlnmx,
     .    mdimx,lx,nx,nxx,  !Bloch wave
     .    cgr,nl-1,         !rotated CG
     .    ppbrd,            !radial integrals
     o    ppb(1,is))
      enddo
      if (rank == 0) print *, ' end of ppbafp_v2'

!     Takao Apr 2019 Added some parts from tetwt5
      eibzmode = eibz4x0() !.and.(.not.nolfco)
!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
!!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
!     Takao Apr 2019 (eibz4x0()) -> (eibzmode)
      if (eibz4x0()) then
!! commentout block inversion Use iqxendx=iqxend because of full inversion
         if (verbose()>=30) call cputid2(' start eibzmode block',0)
         iqxendx=iqxend
         allocate(nwgt(nqbz,iqxini:iqxendx),igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
     .     eibzsym(ngrp,-1:1,iqxini:iqxendx))
         iprintx=.false.
         if (mpi__root) iprintx=.true.
         call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,timereversal(),ginv,iprintx,
     o        nwgt,igx,igxt,eibzsym,tiii)
         if (rank == 0) write(6,"(' Use TimeReversal with EIBZ = ',l1)") tiii

c$$$
c$$$         write(6,"('TimeRevesal switch = ',l1)") timereversal()
c$$$         call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,timereversal(),ginv,iprintx,
c$$$     o        nwgt,igx,igxt,eibzsym)
c$$$!! Check timereversal is required for symmetrization operation or not. If not tiii=timereversal=F is used.
c$$$!! this is because the symmetrization is a little time-consuming.
c$$$         tiii=timereversal()
c$$$         if (minval(igxt)==1) tiii=.false.
c$$$         iprintx=.true.
c$$$cccccccccccccccccccc
c$$$c     tiii=.true.
c$$$cccccccccccccccccccc
c$$$         write(6,"('=== goto eibzgen === used timereversal=',l1)")tiii
c$$$         call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,tiii,ginv,iprintx,
c$$$     o        nwgt,igx,igxt,eibzsym)

        call PBindex(natom,lx,l2nl,nx) ! stores indices in m_pbindex for product basis rotation
            ! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
        call readqgcou() !no input. Read QGcou and store date into variables.
!!      call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.
      else !dummy allocation to overlaid -check bound !sep2014
         iqxendx=iqxend
         allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &    ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      endif

!! == Calculate x0(q,iw) and W == main loop 1001 for iq.
!! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
!! Thus it is necessary to do iq=1 in advance to perform iq >nqibz, modify do 1001 loop.
      allocate(llw(nw_i:nw,nq0i),llwI(niw,nq0i))
      llw=1d99; llwI=1d99 ! for debugging
      allocate(llwmod(nw_i:nw,nq0i),llwImod(niw,nq0i))
      llwmod=1d99; llwImod=1d99 ! for debugging

      if (ixc==1011) then ! ixc==11 is a debug mode to test contribution at \Gamma point.
        last_q = .true.
        goto 1191
      endif

!! for w4phonon. all nodes have wmu array.
      w4pmode=.false.
C     if(sum(ixyz)/=0) w4pmode=.true.
      if(w4pmode) then
        allocate( wmuk(2:nblochpmx,3))
        wmuk=1d99
      endif


      if (rank == 0) then
        write(*,321) iqxend-iqxini+1, nq0i
  321   format(/' --- Start Q - loop ...',i3,' mesh points including',i2,' offset Gamma points')
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
      do 1001 iq = iqxini + q_start, iqxend, q_step ! q=(0,0,0) is omitted!  !==== Big q-loop ====!

        walltimeq = dwtime()

#ifdef USE_MPI
C        if (q_comm%ID == 0) write (*,*) "q point ",iq," running on",q_comm%group_ID,"(",rank,")"
        write(*,*) "Beginning of q-loop:", rank, q_comm%group_ID, q_comm%ID, " iq = ",iq
!       write (0,'(" [",I4,"](",I3,",",I3,") Running iq = ",I4)') rank, q_comm%group_ID, q_comm%ID, iq
!       if (q_comm%ID == writer1) ifrcw  = iopen('WVR.'//charnum5(iq),0,-1,-mrecl)
        if (q_comm%ID == writer2) ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
#else
        ifrcw = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
        if (q_comm%ID == writer2) ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
#endif
        call cputid (0)
        q = qibze(:,iq)
        if (q_comm%ID == 0) then
          write(6,"(' q vec =',3f14.6)") q
          call readqg('QGcou', q, ginv,  quu,ngc,ngveccB) ! q was qq
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

        if (newaniso2.and.iq==1) then ! *sanity check
          if (sum(q**2)>1d-10) then
            call rx( ' hx0fp0: sanity check. |q(1)| /= 0')
          endif
        endif

C   --- Read diagonalized Coulomb interaction ---
        ngb = nbloch + ngc
#ifdef USE_MPI
         write(0,"('[',I4,'](',I3,',',I3,')  do 1001: iq=',i3,' of',i3,'  q=',3f12.6)")
     &     rank,q_comm%group_ID, q_comm%ID, iq,iqxend-iqxini+1,q
#else
         write(6,"(' do 1001: iq=',i3,' of',i3,'  q=',3f12.6)") iq,iqxend-iqxini+1,q
#endif

!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
        if (newaniso2) then
C         Each q-point has its own file.  BUT, only the head of each q-group should do the reading!!
          iqv=iq
          if (iqv>nqibz) iqv = 1  ! Use vcoul at q=(0,0,0) for offset gamma points
          if (q_comm%ID == 0) then
            vcoudfile='Vcoud.'//charnum5(iqv) !this is closed at the end of do 1001.
            ifvcoud = iopen(trim(vcoudfile),0,-1,0)
C           rewind ifvcoud
            read(ifvcoud) ngb0
            if (ngb0/=ngb) then
              call rx( 'hx0fp0.m.f:ngb0/=ngb')
            endif
            read(ifvcoud) qvv
            if (sum(abs(qvv-q))>1d-10 .and. iq<=nqibz) then
              print *,'qvv =',qvv
              call rx( 'hx0fp0: qvv/=0 hvcc is not consistent')
            endif
          endif                 ! q_comm%ID == 0
#ifdef USE_MPI
          call MPI_Bcast(ngb0,1,MPI_INTEGER,0,q_comm%communicator,ierror)
          call MPI_Bcast(qvv,3,MPI_REAL8,0,q_comm%communicator,ierror)
#endif
C         Allocation needs to be done on all nodes
          if (allocated(zcousq)) deallocate( zcousq,vcousq )
          allocate(zcousq(ngb0,ngb0),vcousq(ngb0))
          if (q_comm%ID == 0) then
            read(ifvcoud) vcousq
            read(ifvcoud) zcousq
          endif                 ! qcomm%ID == 0
#ifdef USE_MPI
!         if (q_comm%ID == 0) write(0,'(" [",I4,"] ngb0 = ",I6,"l ngb = ",I6)') rank ,ngb0,ngb
          call MPI_Bcast(vcousq,ngb0,MPI_REAL8,0,q_comm%communicator,ierror)
          call MPI_Bcast(zcousq,ngb0*ngb0,MPI_Double_Complex,0,q_comm%communicator,ierror)
#endif
!         call zprm('zcousq',zcousq,ngb0,ngb0,ngb0)

          idummy = iclose(trim(vcoudfile))

          do  igb1 = 1, ngb0
            if (vcousq(igb1) <= 0) then
              if (q_comm%ID == 0) write(*,"(' warning!: vcousq negative for igb =',i6,' vcousq =',1pe11.3)")
     .        igb1,vcousq(igb1)
              vcousq(igb1) = 1d-6
            endif
            vcousq(igb1) = sqrt(vcousq(igb1))
          enddo
        endif                   ! newaniso2
        nmbas = ngb; nmbas1 = ngb; nolfco = .false.

!     Takao Apr 2019 Update new localfieldcorrectionllw())
        if (newaniso2.and. iq>nqibz.and.(.not.localfieldcorrectionG)) then ! For eps at q=0
          if (ngb0/=ngb) then
            write(0,'(" [",I4,"] iq = ",I5,", ngb0 = ",I5,"; ngb =",I5,": STOP!")') rank, iq, ngb0, ngb
#ifdef USE_MPI
            call MPI_ABORT(MPI_COMM_WORLD,-1,ierror)
#endif
            call rx( 'hx0fp0.m.f:ngb0/=ngb')
          endif
          nolfco =.true.  ! For offset points, turn off local fields
!         nmbas=ngb
          nmbas1 = 1  !nmbas1=1 is for P\mu\nu (\mu=1) components. (nmbas1=2 was tested for \nu=1, but need trick after x0kf_v4h)
        endif

        if (newaniso2) then
          if (allocated(ppovlz)) deallocate(ppovlz)
          if (allocated(ppovl)) deallocate(ppovl)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
           if (q_comm%ID == 0) call readppovl0(q,ngc,ppovl)
#ifdef USE_MPI
           call MPI_Bcast(ppovl,ngc*ngc,MPI_Double_Complex,0,q_comm%communicator,ierror)
#endif
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
!         ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
          if (ngb /= ngb0) stop 'ngb /= ngb0'
          call zgemm('n','n',ngc,ngb,ngc,cmplx1,ppovl,ngc,zcousq(nbloch+1,1),ngb0,cmplx0,ppovlz(nbloch+1,1),ngb)
c         if (allocated(zzr)) zzr = ppovlz
          if (q_comm%ID == 0) write(*,*) 'nbloch ngc ngb =',nbloch,ngc,ngb
        endif
!       call zprm('ppovlz',ppovlz,ngb,ngb,ngb)


c$$$        if (ixc==11) then
c$$$          print *," xxx2: memsize 8*ngb*ngb*nwhis=", 8*ngb*ngb*nwhis,' ngb nwhis=',ngb,nwhis
c$$$          allocate( rcxq(ngb,ngb,nwhis,npm) )
c$$$          rcxq=(0d0,0d0)
c$$$        else
c$$$          if (onceww(2)) print *," xxx2:allocate zxq zxqi memsize 16*ngb*ngb*(nwp+niw)=",
c$$$     &    16*ngb*ngb*(1+nwp+niw),' ngb nwp niw=',ngb,nwp,niw
c$$$          allocate(
c$$$     &    zxq (ngb,ngb,nw_i:nw),   !,nwp) feb2006
c$$$     &    zxqi(ngb,ngb,niw))
c$$$          zxq=0d0; zxqi=0d0
c$$$        endif

!     Takao Apr 2019 rcxq(nmbas1,nmbas2,nwhis,npm)
        is = idalloc('rcxq',iverbose+2,nmbas1*nmbas,nwhis*npm*2)
        allocate(rcxq(nmbas1,nmbas,nwhis,npm) )
        allocate(zw0(ngb,ngb))
!       allocate(zxq(ngb,ngb,nw_i:nw),zxqi(ngb,ngb,niw))
        rcxq = 0d0

C   --- For each spin, do ---
        do 1003 is = 1, nspinmx
C         if (q_comm%ID == 0) write(6,"(' ### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ### ')") iq, is, nqibz + nq0i,nspin
C         if (debug) print *,' niw nw=',niw,nw
          isf=is

!! --- tetrahedron method --- preparation
!!  nbnb    : total number of pairs (ib,jb) with non-zero weight.
!!  n1b  = band index for occ.   1\ge n1b \ge nband+nctot.
!!  Valence-core order(Core index follows valence index).
!!  n2b  = band index for unocc. 1\ge n2b \ge nband
          if (tetra) then
            if (q_comm%ID == 0)
     .        write(6,"(/' Frequency integration by tetrahedron method on real axis : niw nw nwhis=',3i5)") niw,nw,nwhis
!!          ekxx1  for rk, ekxx2 for q+rk See tetwt4
!     Takao Apr 2019 call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
            do  kx = 1, nqbz
              if (q_comm%ID == 0) then
                call readeval(qbz(:,kx),is,ekxx1(1:nband, kx) )
                call readeval(q+qbz(:,kx),is,ekxx2(1:nband, kx) )
              endif
            enddo
#ifdef USE_MPI
            call MPI_Bcast(ekxx1, nband*nqbz, MPI_REAL8, 0, q_comm%communicator, ierror)
            call MPI_Bcast(ekxx2, nband*nqbz, MPI_REAL8, 0, q_comm%communicator, ierror)
#endif

C          call gettetwt(q,iq,is,isf,nwgt(:,iq),frhis,nwhis,npm,
C     i     qbas,ginv,ef,nqibz,nband,ekxx1,ekxx2,nctot,ecore,
C     i     nqbz,qbz,nqbzw,qbzw,ntetf,idtetf,ib1bz,
C     i     nbmx,ebmx,mtet,eibzmode) !nov2016

C       ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
C           The pairs are not dependent on the energy omega in the denominator of the dielectric function.
C           if (q_comm%ID == 0) write(6,"(' ... Initial call to tetwt5 for array dimensions')")
            job = 0
            if (npm==1) then
              ncc=0
            else
              ncc=nctot
            endif

!           Takao Apr 2019 call tetwt5x_dtet4 removed here
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
     i          iq,is,is,nqibz,eibzmode,nwgt(1,iq),
     i          nbmx,ebmx,mtet)
!             write(0,'(" [",I4,"](",I3,",",I3,") tetwt5x_dtet4 initialization done.")') rank, q_comm%group_ID, q_comm%ID
              deallocate(ibjb,ihw,jhw,nhw,whw)

!             write(0,'(" [",I4,"](",I3,",",I3,") nband = ",I5)') rank, q_comm%group_ID, q_comm%ID, nband
!             write(0,'(" [",I4,"](",I3,",",I3,") nctot = ",I5)') rank, q_comm%group_ID, q_comm%ID, nctot
!             write(0,'(" [",I4,"](",I3,",",I3,") ncc   = ",I5)') rank, q_comm%group_ID, q_comm%ID, ncc
!             write(0,'(" [",I4,"](",I3,",",I3,") nqbz  = ",I5)') rank, q_comm%group_ID, q_comm%ID, nqbz
!             write(0,'(" [",I4,"](",I3,",",I3,") npm   = ",I5)') rank, q_comm%group_ID, q_comm%ID, npm

            endif               ! q_comm%ID == 0
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

            if (debug .and. q_comm%ID == 0) then
              do kx = 1, nqbz
                do jpm = 1, npm
                  write(6,"('jpm kx  minval n1b n2b=',4i5)") jpm,kx,
     &                   minval(n1b(1:nbnb(kx,jpm),kx,jpm)),
     &                   minval(n2b(1:nbnb(kx,jpm),kx,jpm))
                enddo
              enddo
            endif

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

c=======TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======
C --- METHOD (tetwt5) for the tetrahedron weight
!   Histogram secstions are specified by frhis(1:nwp)
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].
!     nwp=nw+1; frhis(1)=0
          if (frhis(1)/=0d0) call rx('hx0fp0: require frhis(1)=0d0')

C     ... make index sets
          allocate(ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm))
          jhwtot = 1
          do jpm =1,npm
            do ik   = 1,nqbz
              do ibib = 1,nbnb(ik,jpm)
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

C   ...   Populate ibjb
          whw = 0d0; ibjb = 0
          do jpm=1,npm
            do ik   = 1,nqbz
              do ibib = 1,nbnb(ik,jpm)
                ib1  = n1b(ibib,ik,jpm)
                ib2  = n2b(ibib,ik,jpm)
                ibjb(ib1,ib2,ik,jpm) = ibib
              enddo
            enddo
          enddo

C       ... Generate the histogram weights whw
          job = 1
C         if (q_comm%ID == 0) write(6,*) 'begin making histogram weights'
          allocate(iwgt(1,1,1,1),demin(1,1,1,1),demax(1,1,1,1)) ! Dummies to avoid pointer problems

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
     i            iq,is,is,nqibz,eibzmode,nwgt(1,iq),
     i          nbmx,ebmx,mtet)
          if (q_comm%ID == 0) write(*,"(' completed histogram weights'$)")
          call cputid(0)
          print *, ' '
          deallocate(iwgt,demin,demax)
!! ======TetrahedronWeight_5 block end =========

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

          call x0kr(npm,ihw,nhw,jhw,whw,nhwtot,n1b,n2b,nbnbx,nbnb,q,is,isf,qbas,qbz,
     i      nblocha,nlnmv,nlnmc,iclass,ppb(1,is),icore,ncore,nlmto,nqbz,noccx,noccxv,
     i      nl,nclass,natom,nnc,nlnmx,mdimx,nbloch,nwhis,iq,ngb,ngc,nband,
     i      nolfco,zzr,nmbas1,nmbas,ppovlz,ppovl,chipm,eibzmode,nwgt(1,iq),ebasis,
     o      rcxq)

C         Restore n1b,n2b to its original (valence,core) precedence
C         call alignn12(0,nbnbx,nqbz,0,npm,nband,nctot,ncc,nbnb,n1b,n2b)

          if (is==nspinmx .and. eibzmode .and. nmbas1 == ngb .and. nmbas == ngb) then
            if (ebasis /= 0) then
              call rcxq2Ebas(1,nmbas,nmbas,nwhis,npm,ngb,ppovlz,rcxq)
            endif
            call x0k_sym(npm,q,ginv,nbloch,nwhis,ngb,ngc,nolfco,zzr,nmbas,zcousq,chipm,ngrp,eibzsym(:,:,iq),rcxq)
          elseif (is==nspinmx .and. nmbas1 == ngb .and. nmbas == ngb) then
            call x0k_sym(npm,q,ginv,nbloch,nwhis,ngb,ngc,nolfco,zzr,nmbas,zcousq,chipm,1,[0d0],rcxq)
          endif

! ! dmt: this is probably a serial killer on xfi so i'd commenting it out.
!           write(*,'(" I am ",2I4,":  end of x0kr sum rcxq(",I4,",",I4,")=",2F15.5)')
!      &         q_comm%ID, rank, iq, is, sum(abs(rcxq))

          deallocate(ihw,nhw,jhw,whw,ibjb)
          if (tetra) deallocate(n1b,n2b)
 1003   continue ! end of spin loop
        if (q_comm%ID == 0) print *, 'end of spin-loop nwp=',nwp !end of spin-loop
C       call zprm('rcxq',rcxq,nmbas1,nmbas1,nmbas*nwhis)

C   --- Hilbert transform.  Generate Real part from Imaginary part ---
        if (allocated(zxq)) deallocate(zxq,zxqi)
        iw = idalloc('zxq',iverbose+2,nmbas1*nmbas,(niw-nw_i+1+niw)*2)
        allocate(zxq(nmbas1,nmbas,nw_i:nw), zxqi(nmbas1,nmbas,niw))
        if (q_comm%ID == 0) write(6,'(" ... begin dpsion5: generate Re chi from Im chi. nwhis nw_i niw  nmbas1 nmbas=",6i5)')
     .    nwhis,nw_i,nw,niw,nmbas1,nmbas
        call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega,imagomega,
     i    rcxq,npm,nw_i,nmbas1,nmbas, ! rcxq is overwritten --- used as work
     o    zxq,zxqi,
     i    chipm,schi,is,ecut,ecuts)
! ! dmt: see note on serial stuff above
!         if (q_comm%ID == 0) print *,'... end dpsion5 checksum W-v',sum(abs(zxq)),sum(abs(zxqi))
        if (allocated(rcxq)) then
          iw = idalloc('rcxq',iverbose+4,nmbas1*nmbas,nwhis*npm*2)
          deallocate(rcxq)
        endif
        is = idalloc(' ',11,1,1)

C   --- RealOmega ---
        if (realomega) then
          if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
          nwmax = nw
          nwmin = nw_i
          if (anfexist) then    !Antiferro acceleration.
            do iw  = nwmin, nwmax
              call anfx0k(natom,nclass,nblocha,iclass,pos,nbloch,ngc,
     i          q, ngveccB,qbas, ! for q+G
     i          anfvec,ibasf,   ! these are antiferro informations.
     .          zxq(1,1,iw))    ! i/o
            enddo
          endif

C     ... prepare for iq0.
          iq0 = iq - nqibz
          if (newaniso2) then
C            if (iq==1) then
C              ifepstinv = iopen('EPS0inv',0,-1,0)
C              write(ifepstinv) ngb
C            endif
            if (iqinit) then
              allocate(sk(ngb,nwmin:nwmax,nq0i),sks(ngb,nwmin:nwmax,nq0i) )
              allocate(skI(ngb,niw,nq0i),sksI(ngb,niw,nq0i))
              iqinit=.false.
            endif
            allocate(epstilde(ngb,ngb))
            allocate(epstinv(ngb,ngb))
          endif

          if (q_comm%ID == 0 .and. newaniso2.and.iq<=nqibz) then
            write(6,"(' === trace check for W-V, q =',3f12.6,' iq =',i3,' ===')") q,iq
            write(6,"(4x,' freq    iw',14x,'W-v(site)',25x,'W-v(tot)',24x,'diff')")
          endif
          if (q_comm%ID == 0 .and. newaniso2.and.iq>nqibz) then
            write(6,"(' === Contribution to eps(q=0), q =',3f12.6,' iq =',i3,' ===')") q,iq
            if(modified_response)then
            write(6,"(4x,' freq    iw',14x,'eps(LFC)',25x,'eps(MOD)',25x,'eps(noLFC)',23x,'diff')")
            else
            write(6,"(4x,' freq    iw',14x,'eps(LFC)',25x,'eps(noLFC)',23x,'diff')")
            endif
          endif

          if(iq>nqibz)ifepsfile=iopen('epsR.'//charnum5(iq),0,-1,mrecl)
          if(iq>nqibz)ifepsfileI=iopen('epsI.'//charnum5(iq),0,-1,mrecl)

C     --- iw loop for real axis ---
          do  1015  iw = nwmin, nwmax
#ifdef USE_MPI
            writer1 = mod(iw,q_comm%size)
            recd = iw-nw_i+1
            if (q_comm%ID == writer1) ifrcw  = iopen('WVR.'//charnum5(iq)//'.'//charnum5(recd),0,-1,mrecl)
#endif
            frr = dsign(freq_r(abs(iw)),dble(iw))
            if (newaniso2.and.iq<=nqibz) then !for mmmw
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

!             Takao Apr 2019 call tetw5x_dtet4 removed here w4p writing eps
                if(iw==0.and.w4pmode) then 
                !static epstinv is saved. For q=0 epstilde (mu=1 skipped). For q/=0 full matrix inversion.
                             !(ix=1 is set for q=0) 
!                ifw4p = ifile_handle()
                  open(ifw4p,file='W4PHONON.'//charnum5(iq),form='unformatted')
                  write(ifw4p) iq,q,ngb,ix !ix=0, or ix=1 for q=0 (iq=1)
                  write(ifw4p) epstinv(ix+1:ngb,ix+1:ngb) 
                  close(ifw4p)
                endif
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
C             if (iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
              if (q_comm%ID == writer1) then
                if (q_comm%ID == writer1) write(ifrcw, rec= 1 ) zw !  WP = vsc-v
                if (q_comm%ID == writer1) write(*,'("[",I4,"] iq = ",I3," iw = ",I3," zw(1,1) = ",2F15.5)') rank, iq, iw, zw(1,1)
                call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
              endif
#else
C             if (iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
!               print *, 'valgrind flags this: write(ifrcw,rec=iw-nw_i+1) zw'
              write(ifrcw,rec=iw-nw_i+1) zw !  WP = vsc-v
              call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
#endif
            endif

            if (newaniso2.and.iq>nqibz) then
              vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
              if (localfieldcorrectionG) then ! Full inversion to calculate eps with LFC
C###########################################
                if(modified_response)then
                do igb2=1,ngb
                  do igb1=1,ngb
                    epstilde(igb1,igb2)=-zxq(igb1,igb2,iw)*(vcousq(igb2))**2
                  if(igb2==1)epstilde(igb1,igb2)=0.0d0
                  if(igb1==igb2)epstilde(igb1,igb2)=1.0d0+epstilde(igb1,igb2)
                  enddo
                enddo
                call matcinv(ngb,epstilde)
                Pbar=0.0d0
                do igb1=1,ngb
                  Pbar=Pbar+epstilde(1,igb1)*zxq(igb1,1,iw)
                enddo
                llwmod(iw,iq0)=1.0d0-Vcou1*Pbar
                endif
C############################################
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
!               Takao Apr 2019 reworked somewhat
                write(ifepsfile,rec=iw-nw_i+1)epstilde
C                if(.not.modified_response)then
                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                if (ix == 0) then
                  call matcinv(ngb,epstinv)
                else
                  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                endif
                llw(iw,iq0)= 1d0/epstinv(1,1)
C                endif
              else
c commentout block inversion
c$$$                     sk  (1:ngb,iw,iq0)= zxq(1,1:ngb,iw)
c$$$                     sks (1:ngb,iw,iq0)= zxq(1:ngb,1,iw)
c$$$c                     sks (1:ngb,iw,iq0)= zxq(2,1:ngb,iw) !nmbas1=2 see z1stcol in x0kf_v4h.
c$$$                     vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
                llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw)
              endif
              if(modified_response)then
                if (q_comm%ID == 0) write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8,0p,2x,2f12.6)') frr,iw,
     .          llw(iw,iq0),llwmod(iw,iq0),llwmod(iw,iq0)-llw(iw,iq0)
                llw(iw,iq0)=llwmod(iw,iq0)
              else
                if (q_comm%ID == 0) write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8,0p,2x,2f12.6)') frr,iw,
     .          llw(iw,iq0),1d0-vcou1*zxq(1,1,iw),llw(iw,iq0)-1d0+vcou1*zxq(1,1,iw)
              endif
            endif

C            if (.ynot.newaniso2) then               ! Original mode
C              call rx('newaniso2 no longer supported')
Cc             call wcf( ngb, vcoul, zxq(1,1,iw), 1, zw0)  ! W= (1-v zxq)^{-1} v
C            endif

C              write(6,"('--- EigenValues for Im(W) ---')")
C              allocate(ebb(ngb))
C              call diagcvh2( (zw0-transpose(dconjg(zw0)))/2d0/img, ngb,ebb)
C              do ii=1,ngb
C                if (abs(ebb(ii))>1d-8 .and. ebb(ii)>0) then
C                  write(6, "('hhhIWq : iw ii eb=',2i4,d13.5)") iw, ii,ebb(ii)
C                else
C                  write(6, "('hhhIWqxxx : iw ii eb=',2i4,d13.5)") iw,ii, ebb(ii)
C                endif
C              enddo
C              deallocate(ebb)

c            if (newaniso2.and.iq>nqibz) then
cc                zw(1:ngb,1:ngb) = 0d0
cc                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
c            else
c                zw(1:ngb,1:ngb) = zw0
cc                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
c                write(ifrcw, rec= iw-nw_i+1) zw   !  WP = vsc-v
c                call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
c            endif

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
        endif                   ! end of RealOmega block

C   --- Imaginary Omega ---
        if (imagomega) then
          print *,' ... start imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
!     Takao Apr 2019 remove anfexist
c     ... Antiferro acceleration.
          if (anfexist) then
            do  iw  = 1, niw
              call anfx0k(natom,nclass,nblocha,iclass,pos,nbloch,ngc,
     i          q, ngveccB,qbas, ! for q+G
     i          anfvec,ibasf,    ! these are antiferro informations.
     m          zxqi(1,1,iw))    ! i/o
            enddo
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
#ifdef USE_MPI
c             if (iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
#else
C             if (iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
#endif
              if (q_comm%ID == writer2) then
                write(ifrcwi, rec=iw)  zw !  WP = vsc-v
                call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
              endif

            else ! iq>nqibz
              vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2 = fourpi/sum(q**2*tpioa**2-eee)
              if (localfieldcorrectionG) then ! Full inversion to calculate eps with LFC.
C########################################
                if(modified_response)then
                do igb2=1,ngb
                  do igb1=1,ngb
                    epstilde(igb1,igb2)=-zxqi(igb1,igb2,iw)*(vcousq(igb2))**2
                  if(igb2==1)epstilde(igb1,igb2)=0.0d0
                  if(igb1==igb2)epstilde(igb1,igb2)=1.0d0+epstilde(igb1,igb2)
                  enddo
                enddo
                call matcinv(ngb,epstilde)
                Pbar=0.0d0
                do igb1=1,ngb
                  Pbar=Pbar+epstilde(1,igb1)*zxqi(igb1,1,iw)
                enddo
                llwImod(iw,iq0)=1.0d0-Vcou1*Pbar
                endif
C###################################################
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
                write(ifepsfileI,rec=iw)epstilde
C                if(.not.modified_response)then
                epstinv(ix+1:ngb,ix+1:ngb) = epstilde(ix+1:ngb,ix+1:ngb)
                if (ix == 0) then
                  call matcinv(ngb,epstinv)
                else
                  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                endif
                llwI(iw,iq0) = 1d0/epstinv(1,1)
C                endif
              else  ! llwI without LFC. LFC contribution is added in
                llwI(iw,iq0) = 1d0 - vcou1*zxqi(1,1,iw) !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
              endif
              if(modified_response)then
              if (q_comm%ID == 0) write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8,0p,2x,2f12.6)')
     .                 frr,iw,llwI(iw,iq0),llwImod(iw,iq0),llwImod(iw,iq0)-llw(iw,iq0)
              llwI(iw,iq0)=llwImod(iw,iq0)
              else
              if (q_comm%ID == 0) write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8,0p,2x,2f12.6)')
     .                 frr,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw),llwI(iw,iq0)-1d0+vcou1*zxqi(1,1,iw)
              endif
c              write(ifisk) iw,iq0,q
c              write(ifisk) vcousq(2:ngb)*w_kI(2:ngb),vcousq(2:ngb)*w_ksI( 2:ngb)
            endif
 1016     continue              ! loop over imaginary frequencies
C         if (iq==1) ifepstinv = iclose('EPS0inv') !iq==1 close write mode.
          if (newaniso2) then
            deallocate(epstinv)
            if (allocated(epstilde)) deallocate(epstilde)
          endif
        endif                   ! Imag omega

        if (allocated(vcoul)) deallocate(vcoul)
        if (allocated(zw0)) deallocate(zw0)
        iw = idalloc('zxq',iverbose+3,nmbas1*nmbas,(niw-nw_i+1+niw)*2)
        if (allocated(zxq)) deallocate(zxq)
        if (allocated(zxqi)) deallocate(zxqi)

        if (q_comm%ID == writer2) ifrcwi = iclose('WVI.'//charnum5(iq))

#ifdef USE_MPI
!       if (q_comm%ID == writer1) ifrcw = iclose('WVR.'//charnum5(iq))
        write(*,'(" [",I4,"](",I3,",",I3,") End of q-loop: iq = ",I3," Time = ",F15.5)')
     .    rank, q_comm%group_ID, q_comm%ID,iq,(dwtime()-walltimeq)/60
#else
        ifrcw = iclose('WVR.'//charnum5(iq))
        write(*,*) "End of q-loop: iq = ", iq
#endif
        if(iq>nqibz)ifepsfile=iclose('epsR.'//charnum5(iq))
        if(iq>nqibz)ifepsfileI=iclose('epsI.'//charnum5(iq))
        if (iq == iqxend) last_q = .true.

 1001 continue                  ! end of big iq loop

C!! === Receive llw and llwI at node 0, where q=0(iq=1) is calculated. ===
CYY!$OMP parallel
CYY!$OMP master
C         print *,' eeeeeeeeeeeeeeeee111 mpi__rank=',mpi__rank
C         if (mpi__size/=1) then
C            do iq=nqibz+1,iqxend
C               iq0 = iq - nqibz
Cc     print *,' iq iq0 mpi_rank mpi_ranktab(iq)=',iq, iq0,mpi__rank,mpi__ranktab(iq),mpi__root,nw,nw_i,niw
C               if (mpi__ranktab(iq)/=0) then !jan2012
C               if (mpi__ranktab(iq) == mpi__rank) then
Cc     print *,' mpi_send iq from',iq,mpi__ranktab(iq)
Cc     print *,' send llw sum=',sum(abs(llw(:,iq0))),nw,nw_i
Cc     do i=nw_i,nw
Cc     write(*,*)'sendxxx',i,llw(i,iq0)
Cc     enddo
Cc     print *,' send llwI sum=',sum(abs(llwI(:,iq0))),niw
C                  dest=0
C                  call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
C                  call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
C               elseif (mpi__root) then
Cc     print *,' mpi_recv iq from',iq,mpi__ranktab(iq),nw,nw_i,niw
C                  src=mpi__ranktab(iq)
C                  call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
C                  call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
Cc     do i=nw_i,nw
Cc     write(*,*)'recivxxx',i,llw(i,iq0)
Cc     enddo
Cc     print *,' recv llw sum=',sum(abs(llw(:,iq0))),nw,nw_i
Cc     print *,' recv llwI sum=',sum(abs(llwI(:,iq0))),niw
C               endif
C               endif
C            enddo
C         endif
C         print *,' eeeeeeeeeeeeeeeee222 mpi__rank=',mpi__rank
CCYY!$OMP end master
CCYY!$OMP end parallel

#ifdef USE_MPI
!     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      write(0,'(" [",I4,"](",I3,",",I3,") Still alive after q-loop:",L5)') rank, q_comm%group_ID, q_comm%ID, last_q
#endif

c commentout block inversion
c$$$!! Add LFC (local field correction) to llw and llwI
c$$$         if (newaniso2 .and. mpi__rank == 0) then ! only on root node
c$$$            iq=1 !for q=0
c$$$            vcoudfile='Vcoud.'//charnum5(iq)
c$$$            ifvcoud = iopen(trim(vcoudfile),0,-1,0)
c$$$            read(ifvcoud) ngb0
c$$$            read(ifvcoud) qvv
c$$$            if (sum(abs(qvv))>1d-10) then
c$$$               print *,'qvv =',qvv
c$$$               stop 'hx0fp0: qvv/=0 hvcc is not consistent'
c$$$            endif
c$$$            if (allocated(zcousq0)) deallocate( zcousq0,vcousq0 )
c$$$            allocate( zcousq0(ngb0,ngb0),vcousq0(ngb0))
c$$$            read(ifvcoud) vcousq0
c$$$            read(ifvcoud) zcousq0
c$$$            idummy=iclose(trim(vcoudfile))
c$$$            vcousq=sqrt(vcousq)
c$$$            allocate(epstinv(ngb0,ngb0),w_k(ngb0),w_ks(ngb0),w_kI(ngb0),w_ksI(ngb0),eemat(ngb0,ngb0))
c$$$
c$$$            do iq0=1,nq0i
c$$$              iq = iq0 + nqibz
c$$$              q = qibze(:,iq)
c$$$
c$$$              vcoudfile='Vcoud.'//charnum5(iq)
c$$$              ifvcoud = iopen(trim(vcoudfile),0,-1,0)
c$$$              read(ifvcoud) ngb
c$$$              read(ifvcoud) qvv
c$$$              if (sum(abs(qvv-q))>1d-10) then
c$$$               print *,'qvv =',qvv
c$$$               stop 'hx0fp0: qvv/=0 hvcc is not consistent'
c$$$              endif
c$$$              if (allocated(zcousq)) deallocate(zcousq)
c$$$              if (allocated(vcousq)) deallocate(vcousq)
c$$$              allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
c$$$              read(ifvcoud) vcousq
c$$$              read(ifvcoud) zcousq
c$$$              idummy=iclose(trim(vcoudfile))
c$$$              vcousq=sqrt(vcousq)
c$$$
c$$$              ifepstinv = iopen('EPS0inv',0,0,0)
c$$$              read(ifepstinv) ngb
c$$$
c$$$               ngc=ngb-nbloch
c$$$               if (allocated(ppovlz)) deallocate(ppovlz)
c$$$               if (allocated(ppovl)) deallocate(ppovl)
c$$$               allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
c$$$               call readppovl0(q,ngc,ppovl) !q was qq
c$$$               ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
c$$$               ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
c$$$
c$$$!  eemat: Z\mu_i(\bfk=0)^* <i|j> Z\nu_j(\bfk)
c$$$               eemat =matmul(transpose(dconjg(zcousq0)),matmul(ppovlz,zcousq))
c$$$               vcou1  = fourpi/sum(q**2*tpioa**2) ! test-->vcousq(1)**2 !fourpi/sum(q**2*tpioa**2-eee)
c$$$               vcou1sq = vcou1**.5
c$$$               write(*,*)
c$$$
c$$$              do iw=nwmin,nwmax
c$$$                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$                epstinv(2:ngb,2:ngb) = matmul( transpose(dconjg(eemat(2:ngb,2:ngb))),
c$$$     &                                matmul(epstinv(2:ngb,2:ngb),eemat(2:ngb,2:ngb)) )
c$$$                if (iw/=iwx) then
c$$$                print *,'iw iwx=',iw,iwx
c$$$                stop 'hx0fp0_sc: iw/=iwx'
c$$$                endif
c$$$                w_k(2:ngb) = vcou1sq*matmul( epstinv(2:ngb,2:ngb), sk(2:ngb,iw,iq0)*vcousq(2:ngb))
c$$$                epslfc = -vcou1sq*sum( sks(2:ngb,iw,iq0) * w_k(2:ngb) *vcousq(2:ngb) )
c$$$                llw(iw,iq0) = llw(iw,iq0)  + epslfc
c$$$                write(6,"('eps(on real) iq iw',2i4,2f9.3,2x,2f9.3)") iq0,iw, llw(iw,iq0)-epslfc,llw(iw,iq0)
c$$$              enddo
c$$$              do iw=1,niw
c$$$                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$                if (iw/=iwx) then
c$$$                 print *,'iw iwx=',iw,iwx
c$$$                 stop 'hx0fp0_sc: iw/=iwx'
c$$$                endif
c$$$                w_kI(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), skI(2:ngb,iw,iq0)*vcousq(2:ngb))
c$$$                epslfc=- vcou1sq*sum( sksI(2:ngb,iw,iq0)* w_kI(2:ngb)*vcousq(2:ngb) )
c$$$                llwI(iw,iq0)= llwI(iw,iq0)+epslfc
c$$$                write(6,"('eps(on img ) iq iw',2i4,2f9.3,2x,2f9.3)")iq0,iw, llwI(iw,iq0)-epslfc,llwI(iw,iq0)
c$$$              enddo
c$$$              ifepstinv = iclose('EPS0inv')
c$$$           enddo
c$$$         endif

C --- W(0) divergent part and non-analytic constant part ---
C     see comments in m_w0w0i.F
 1191 continue
      if (last_q .and. q_comm%ID == 0) then

!!      ix=1011 is a special mode to overwrite llw and llwI for test purpose
!!      File W0W0I is generated by calling w0w0i, but not used except for this mode
        if (ixc==1011) then
          ifw0w0i = iopen('W0W0I',0,-1,0)
          read(ifw0w0i) nw_ixx,nwxx,niw,nq0ix
          write(6,*)'w0w0i: n=',nw_ixx,nwxx,niw,nq0ix
          if(nq0i/=nq0ix)  call rx('nq0i/=nq0ix')
          if(nw_i/=nw_ixx) call rx(nw_i/=nw_ixx)
          if(nw/=nwxx) call rx(nw/=nwxx)
          read(ifw0w0i) llw(nw_i:nw,1:nq0i)
          read(ifw0w0i) llwI(1:niw,1:nq0i)
          ifw0w0i = iclose('W0W0I')
        endif

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

!      is = iclose('hbe.d')

#ifdef USE_MPI
      call mpi_end
#endif

      is = idalloc(' ',1,1,1)
      write(*,"(/' Max dynamic allocation :',i7,' MB'/)") is
      call cputid(0)

! Check whether the next lines should be before mpi_end !!
      if (nthreads > 1) then
        print 131, (dwtime() - walltime0)/60, cpusec()/60, nthreads
      else
        print 131, (dwtime() - walltime0)/60, cpusec()/60
      endif
  131 format(' hx0fp0_sc :  Wall clock time',f9.1,'m','  Processor time',f9.1,'m':,i4,' threads')

!---------------------------------------------------------
      if (ixc==11) call rx0( ' OK! hx0fp0_sc ixc=11 Sergey F. mode')
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
