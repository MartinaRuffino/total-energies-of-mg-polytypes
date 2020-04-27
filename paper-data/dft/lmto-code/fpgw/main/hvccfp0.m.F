      program hvccfp0
c- Coulomb matrix. <f_i | v| f_j>_q. ------------------------
c input files
c    HVCCIN      : some inputs by hbg0.
c    PLN         : plane wave expansion data by nbg0.
c    BASFP//atom : product basis by hbasfp0. ic=iatom should be kept!
c output
c    VCCFP : the coulomb matrix vcoul(nblochpmx,nblochpmx) for all qibz.
c-------------------------------------------------------------
c int
c    strx: structure constant for e=0 (means 1/|r-r'| )
c
      use m_readqg
      use keyvalue
C      use m_mpi,only: MPI__hx0fp0_rankdivider2,mpi__task,MPI__Initialize,MPI__Finalize,mpi__root,
C     &   MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,mpi__rank,mpi__size,
C     &   mpi__ranktab,MPI__consoleout,mpi__iend,mpi__iini,mpi__getrange

      use mpi_mod, q_comm => comm_1, q_inter => inter_1, rank => mpi_rank

      use omp_lib

      implicit none
      integer(4) :: ifvcfpout,iopen,ifhvccfp,is, nqbz, nbas,lmxcg,
     &  nband, ifplane,ngpmx, ngcmx, nblochpmx, nbloch,
     &  ibas,ic,lxx,nxx,nrx,l,n,k,isx,kdummy, iclose,
     &  nkdmx,nkqmx,lmax,nkdest,nkrest,ngp,ngc,nlxx,i,lnjcg,lnxcg,
     &  nkd,nkq ,ibas1,ibas2,nlx1,nlx2, nqibz,iqibz
      real(8) :: alat, plat(3,3),qlat(3,3),q(3),p(3),voltot,
     &  pi,fpi,tripl,alat0,epsx,
     &  tol,as,tpiba,qb0(3,3),vol0,rdist0,qdist0,radd,qadd,
     &  a0,awald,alat1,tol1,r0,q0,awald0,qg(3),   absqg2,aaa,aaa12
      integer(4),allocatable :: jcg(:),indxcg(:),
     &  lx(:),kmx(:),nblocha(:),nr(:),ificrb(:),
     &  nx(:,:),ngvecp(:,:),ngvecc(:,:),ngvecci(:,:,:),iqibzx(:)
      real(8),allocatable :: qbz(:,:),qibz(:,:),bas(:,:),rmax(:),
     &  cg(:),rprodx(:,:,:,:),dlv(:,:),qlv(:,:),work(:),
     &  rojb(:,:,:), sgbb(:,:,:,:),aa(:),bb(:),rofit(:),phi(:),psi(:),
     &  wqt(:), q0i(:,:)
      complex(8) ,allocatable :: vcoul(:,:),geig(:,:),strx(:,:,:,:),
     &  sgpb(:,:,:,:),sgpp(:,:,:,:),
     &  fouvb(:,:,:,:),fouvp(:,:,:,:),vcoul0(:,:),
     &  s(:,:),sd(:,:),rojp(:,:,:) , vcoulnn(:,:)
      character*7,allocatable :: filename(:)
      character(20) :: xxt
      logical :: mpi__root = .true.
      integer mpi__iini,mpi__iend

c dimension work space
c      integer(4),parameter :: nsize= 1 000 000
c      real w(nsize)
c      common /w/ w

C     For managing negative eigenvalues of vcoul
      real(8), parameter :: toleb = 1d-10
      integer nposeb  ! number of positive eigenvalues
      real(8) :: maxposeb

      complex(8):: phasep,img=(0d0,1d0)
      integer(4)::ir,ig1,n1,n2
c check
      complex(8),allocatable :: hh(:,:),oox(:,:),ooxi(:,:),oo(:,:),zz(:,:),zzr(:)
      real(8),allocatable    :: eb(:)

      complex(8),allocatable :: matp(:),matp2(:)
      complex(8) :: xxx,trwv
      integer(4) :: ngb,nev,nmx,iqx,ipl1,ipl2,nq0i,igx1,igx2
      logical checkeig
      logical:: besseltest=.false. !test
      real(8) :: sss1,sss2,dnorm
c

      complex(8),allocatable:: gbvec(:), ppovl(:,:), b0mat(:)

      integer(4) ::igc,igc0,ifgb0vec,ifgb0vec1,ix, iy

      integer(4) :: iqxini, iqxend,imode
      logical :: allochk=.false. !paralellx0=.true.,

      integer :: idalloc,iverbose,iw ! Parameters for printout

      complex(8),allocatable:: hh1(:,:),oo1(:,:)
      integer(4):: nqnumc,ifiqgc
c      character(5)  :: charnum5
c      integer(4),allocatable:: iqok(:)
      real(8):: qqq(3),QpGcut_Cou       !,qq(3)

      integer(4),allocatable:: ngvecc0(:,:)
      integer(4):: ngc0

      real(8):: ginv(3,3),quu(3),det

c---
      real(8),allocatable :: rkpr(:,:,:),rkmr(:,:,:),rofi(:,:)
      real(8):: eee,eees, q_org(3),screenfac
      integer(4):: ifvcfporg,nqbz_in,nblochpmx_in
      complex(8),allocatable:: vcoul_org(:,:)

      logical :: smbasis,debug=.false.,smbb
      integer(4)   :: ifprodmt,nl_r,lx_,nxx_r,nxdim,ibl1,nn,no,ngbnew,jobBes,
     &  nmatch,ifpmatch,nmatch_q,ifpmatch_q,m,ifpomat,nbln,ibln,ngb_in,nnr,igc2
      character(3) :: charnum3
      character(5) :: charnum5
      character(11):: filenamep
      integer(4),allocatable:: nx_r(:), ibl(:,:,:,:)
     &  ,imatcho(:),imatchn(:),imatcho_q(:),imatchn_q(:)
      real(8),allocatable:: prodmt(:,:,:,:),rdmatch(:,:,:,:)
      complex(8),allocatable:: ppmt(:,:,:,:),pmat(:,:),pomat(:,:),oon(:,:)
      complex(8):: pval,pslo,phasex
      real(8)::absqq,qqx(3), epsmx,aaaa
      integer(4):: nnmx ,ngcnn,ngbo
cki      integer(4):: is_mix0vec ,ifgb0vec_a,ifgb0vec_b
      integer(4):: ifgb0vec_a,ifgb0vec_b , ifvcoud,idummy
      logical:: is_mix0vec,wvcc
      character(128):: vcoudfile
      integer(4):: verbose
      integer nthreads,omppid
      real(8):: walltime0


      integer num_q_groups
      integer iqxstart, iqxstep, myID
      integer, allocatable :: first(:), last(:), displace(:), recvcount(:)
      double precision a_coeff, b_coeff

      complex(8), allocatable :: pmat_local(:,:)

      logical :: MPIinitexists
      integer,parameter :: NULLI=-99999
      procedure(integer) :: GWversion
      procedure(real(8)) :: cpusec,dwtime

      walltime0 = dwtime()
! C$OMP parallel
!       nthreads = omppid(0)
! C$OMP end parallel
      nthreads = omp_get_max_threads()

      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      nposeb = 0; maxposeb = 0
      jobBes = mod(GWversion()/1000,10)

!     MPI initialization
      call mpi_start()
      mpi_N1 = mpi_size; mpi_N2 = 1
#ifdef USE_MPI
      if(rank == 0) then
        inquire(file='MPI_VCC.init',EXIST=MPIinitexists)
        if (MPIinitexists) then
          open(100,file="MPI_VCC.init")
          read(100,*) mpi_N1, mpi_N2
          close(100)
        endif
      endif
      call MPI_Bcast(mpi_N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(mpi_N2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      num_q_groups = mpi_N1
      call mpi_make_groups(comm_world,iqxend-iqxini,num_q_groups,q_comm,q_inter)
#ifdef USE_MPI
      write(0,'("MPI: rank = ",I4,", q_comm%ID  = ",I4,", q_comm%group_ID  = ",I4)')
     &          rank, q_comm%ID, q_comm%group_ID

      write(0,'("MPI: rank = ",I4,", q_inter%ID = ",I4,", q_inter%group_ID = ",I4)')
     &          rank, q_inter%ID, q_inter%group_ID
#endif

      call getjobmode(imode)
      if (imode == NULLI) then
#ifdef USE_MPI
        call rx('run command with --job=#')
#else
        write(6,"(a)") ' enter mode: (only 0 for now)'
        read(5,*) imode
#endif
      endif

      call headver('hvccfp0: start',imode)
      call cputid (0)
      iverbose = (verbose()/30)*10
      if(imode==202 ) then
C       write(6,*)' hvccfp0: imode=',imode
c      elseif(imode==101) then
c        write(6,*)' hvccfp0: imode=',imode
c        write(6,*)' remove_r0c is effective'
c        write(6,*)' Generate VCCFP = VCCFP.ORG - new_VCCFP'
c        ifvcfporg = iopen( "VCCFP.ORG",0,-1,0)
c      elseif(imode==102) then
c        write(6,*)' hvccfp0: imode=',imode
c        write(6,*)' remove_r0c is effective'
      elseif(imode==0) then
      elseif(imode==3) then
      else
        call rx( 'hvccfp0: now hvccfp0 support just normal mode=0 3 202 101')
      endif
c      if(iqxini< 2) paralellx0=.false.
c      if(paralellx0) then
c        write(6,"(' PARALELL.X0 mode: iqxini iqxend=',2i3)")
c     &  iqxini, iqxend
c      endif

C --- q, nqbz, alat, qlat, nbas, bas
      if (rank == 0) then
         ifhvccfp = iopen('HVCCIN',0,-1,0)
         read (ifhvccfp) alat, plat,qlat, nqbz, nbas, nband
      endif
#ifdef USE_MPI
      call MPI_Bcast(alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(plat, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(qlat, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nqbz,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nbas,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nband, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      if(allochk)
     & write(*,*) 'allocate(qbz(3,nqbz),bas(3,nbas),rmax(nbas))'

      allocate(qbz(3,nqbz),bas(3,nbas),rmax(nbas))

      if(rank == 0) then
         read(ifhvccfp) qbz, bas,rmax
         read(ifhvccfp) nqibz
      endif
#ifdef USE_MPI
      call MPI_Bcast(qbz,   3*nqbz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(bas,   3*nbas, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(rmax,  nbas,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nqibz,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      if(allochk)
     &  write(*,*)'allocate(qibz(3,nqibz),iqibzx(nqbz))'

      allocate(qibz(3,nqibz),iqibzx(nqbz))


      if(rank == 0) read(ifhvccfp) qibz(1:3,1:nqibz)
#ifdef USE_MPI
      call MPI_Bcast(qibz,   3*nqibz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
#endif

      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))

      if(rank == 0) then
        print *,' cell volume =',voltot
        write(6,"(i4,3f13.6)")(i,qibz(1:3,i),i=1,nqibz)
      endif

      call readngmx('QGcou',ngcmx)
      allocate(ngvecc(3,ngcmx))
      if (besseltest) allocate(ngvecci(3,ngcmx,nqibz))

c --- Readin BASFP//atom. The product basis functions.
      if(allochk) write(*,*)'allocte(lx,kmx,nblocha,nr,aa,bb,filename,ificrb'
      allocate(lx(nbas),kmx(nbas),nblocha(nbas),
     &         nr(nbas),aa(nbas),bb(nbas),filename(nbas),
     &         ificrb(nbas) )


      if(rank == 0) then
        do ibas = 1,nbas
          ic = ibas !
          filename(ibas)= 'BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10))
          ificrb(ibas)  = iopen ( filename(ibas),1,3,0)
          read(ificrb(ibas),"(4i6,2d24.16)") lx(ibas), kmx(ibas), nblocha(ibas), nr(ibas),aa(ibas),bb(ibas)
        enddo
      endif
#ifdef USE_MPI
      call MPI_Bcast(lx,      nbas,  MPI_INTEGER,  0,  MPI_COMM_WORLD, ierror)
      call MPI_Bcast(kmx,     nbas,  MPI_INTEGER,  0,  MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nblocha, nbas,  MPI_INTEGER,  0,  MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nr,      nbas,  MPI_INTEGER,  0,  MPI_COMM_WORLD, ierror)
      call MPI_Bcast(aa,   nbas, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(bb,   nbas, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
#endif

      lxx = maxval(lx)
      if(allochk) write(*,*) 'allocate( nx(0:lxx,nbas) )'
      allocate( nx(0:lxx,nbas) )

      if(rank == 0) then
        do ibas = 1,nbas
          read(ificrb(ibas),"(i5)") nx(0:lx(ibas),ibas)
        enddo
      endif
#ifdef USE_MPI
      call MPI_Bcast(nx, (lxx+1)*nbas, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      nxx = maxval(nx)
      nrx = maxval(nr)
      if(allochk) write(*,*) 'allocate( rprodx(nrx,nxx,0:lxx,nbas) )'
      allocate( rprodx(nrx,nxx,0:lxx,nbas) )

      if (rank == 0) then
      do ibas = 1,nbas
        do l = 0, lx(ibas)
          do n = 1, nx(l,ibas)
            read(ificrb(ibas),"(3i5)"   ) k, kdummy,kdummy
            read(ificrb(ibas),"(d23.15)") (rprodx(i,n,l,ibas),i=1,nr(ibas))
          enddo
        enddo
c      isx = iclose(filename(ibas))
      enddo
      endif
#ifdef USE_MPI
      call MPI_Bcast(rprodx, nrx*nxx*(lxx+1)*nbas, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
#endif

cccccccccccccccccccccccccccccccccccccccccccccc
cccc TEST cccccccccccccccc
c      open(117, file='xin')
c      do i=1,nr(1)
c        read(117,"(d24.16)") rprodx(i,1,0,1)
c      enddo
cccccccccccccccccccccccccccccccccccccccccccccccc

C       This test doesn't do anything now.  It should write to fort.196
C       You must set switches in vcoulq also
c TEST cccccccccccccccccccccccccccccccccccccc
C      if (besseltest .and. rank == 0) then
C        print *
C        print *
C        print *, ' *** TEST case ***  rprodx is given by Bessel.'
Cccc test G, corresponding <q+G|v|q+G> should be exact. e.g. ig1=1 and  ig1=35 for iqx=2
Cccc You can change these values for tests. ccccccccccccc
C        iqx  = 2
C        igx1 = 1
C        igx2 = 35
Cc
C
C        do  iqx = 1, nqibz
C          call readqg('QGcou',qibz(1,iqx),ginv,  quu,ngc, ngvecci(:,:,iqx))
C        enddo
C
C        write(6,"(' iqx=',i3,' ig1 ig2=',2i3)") iqx,igx1,igx2
C        write(6,"(a)")
C     & ' <q+G|v|q+G> for the corresponding iqx ig1 ig2 should be exact!'
C        write(6,"(a)") ' See fort.196'
C        write(6,"(a)")
C     & ' Errors will be from the radial function integrals !!!'
C        write(6,"(a)") ' You can slso so similar test from hbasfp0.'
C        write(6,"(a)") ' See test1 in basnfp0.'
Cc
C        deallocate(rprodx,nx)
C        tpiba=8.d0*datan(1.d0)/alat
C        lx = 4
C        nr = nr(1)
C        aa = aa(1)
C        bb = bb(1)
C        lxx = maxval(lx)
C        if(allochk) write(*,*)'allocate( nx(0:lxx,nbas) )'
C        allocate( nx(0:lxx,nbas) )
C        kmx= 1
C        nx = 2
C        nxx = maxval(nx)
C        nblocha= nxx *(lxx+1)**2
C        nrx = maxval(nr)
C        if(allochk) write(*,*)'allocate(rprodx,rofi ,phi,psi) '
C        allocate(rprodx(nrx,nxx,0:lxx,nbas),rofit(nrx),phi(0:lxx),psi(0:lxx))
C        rofit(1) = 0d0
C        do ir   = 1, nrx
C          rofit(ir) = bb(1)*( exp(aa(1)*(ir-1)) - 1d0)
C        enddo
C        do n = 1, nxx
C          if(n==1) ig1 = igx1
C          if(n==2) ig1 = igx2
C          qg(1:3) =
C     &    tpiba * (qibz(1:3,iqx)+ matmul(qlat, ngvecci(1:3,ig1,iqx)))
C          absqg2  = sum(qg(1:3)**2)
Cc
C          do ir =1,nrx
CC           call bessl(absqg2*rofit(ir)**2,lxx,phi,psi)
C            call besslr(absqg2*rofit(ir)**2,10*jobBes+0,0,lxx,phi,psi)
C            do ibas=1,nbas
C              do l = 0, lx(ibas)
C                rprodx(ir,n,l,ibas) = phi(l)* rofit(ir) **(l +1 )
C              enddo
C            enddo
C          enddo
C        enddo
Cc --- orthogonalized rprodx.
C        do ibas=1,nbas
C          do l = 0, lx(ibas)
C            rprodx(1:nr(ibas),1,l,ibas)=
C     &     rprodx(1:nr(ibas),1,l,ibas)
C     &   + rprodx(1:nr(ibas),2,l,ibas)
C            n = 1
C            call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas)
C     &    ,aa(ibas),bb(ibas),nr(ibas), aaa )
C            aaa = 1d0/sqrt(aaa)
C            rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
C            if(nxx==1) cycle
C            n1=1
C            n2=2
C            call gintxx(rprodx(1,n1,l,ibas),rprodx(1,n2,l,ibas)
C     &    ,aa(ibas),bb(ibas),nr(ibas), aaa12 )
C            rprodx(1:nr(ibas),n2,l,ibas) = rprodx(1:nr(ibas),n2,l,ibas)
C     &      - aaa12*rprodx(1:nr(ibas),n1,l,ibas)
C            n = 2
C            call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas)
C     &    ,aa(ibas),bb(ibas),nr(ibas), aaa )
C            aaa = 1d0/sqrt(aaa)
C            rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
C          enddo
C        enddo
C        call rx0('done bessel test')
C
C      endif ! besseltest

      nbloch    = sum(nblocha)
      nblochpmx = nbloch + ngcmx

c --- CG coefficients. <LM3|lm1 lm2>
c inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
c Injcg = indxcg(inxcg) to indxcg(inxcg)-1
c cg(inxcg)  : = <lm3|lm1 lm2>
c jcg(lnjcg) : = lm3
      lmxcg = lxx
      call scg_sizechk(lmxcg,lnjcg,lnxcg) !(lmax,c,cindx,js)
      allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
      call scg(lmxcg,cg,indxcg,jcg)
      call minv33(qlat,ginv)

c --- Get real-space vectors and reciprocal-space vectors for Ewald sum.
C defaults values for ewald sum
c      call lattc(awald0,tol,alat,alat,plat0,gx,gy,gz,gam,plat,qlat,
c     .   lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx,w(owork))
c- taken from lattc.f

c default values ok?
      awald0 = 2d0   !See p_lat_0
      tol    = 1d-9
      nkdmx  = 800
      nkqmx  = 800
      lmax   = 2*lxx  !lxx or lmax=6 ???

      vol0= abs(tripl(plat,plat(1,2),plat(1,3)))
      as   = awald0
      alat1= alat
      tpiba = 8.d0*datan(1.d0)/alat
      call cross_x(plat(1,2),plat(1,3),qb0)
      call cross_x(plat(1,3),plat(1,1),qb0(1,2))
      call cross_x(plat(1,1),plat(1,2),qb0(1,3))
      qb0(1:3,1:3) = qb0(1:3,1:3)/vol0

      rdist0=vol0**(1.d0/3.d0)
      qdist0=1.d0/rdist0
      radd=.7*rdist0
      qadd=.7*qdist0
      a0=as/rdist0
      awald=a0/alat
cccccccccccccccccccccccccccccccccccc
c takao
c      tol1= tol*alat**(lmax+1) *0.01
cccccccccccccccccccccccccccccccccccccc
      tol1= tol*alat**(lmax+1)
      if(allochk) write(*,*) 'allocate(dlv, qlv, work) '
      allocate(dlv(3,nkdmx), qlv(3,nkqmx), work(max0(nkdmx,nkqmx)) )
      call lctoff(a0,vol0,lmax,tol1,r0,q0)
      nkdest =4.18879*(r0+radd)**3/vol0+.5
      nkrest =4.18879*(q0+qadd)**3*vol0+.5

      if(rank == 0) write(6,340) as,tol,lmax,awald,vol0,alat1,nkdest,nkrest
  340 format(/' lattc:  as=',f6.3,'   tol=',1p,e9.2,'   lmax=',i2,
     .  '   awald=',0p,f7.4,'   v0=',f10.3/' alat1=',f9.5,
     .  '   estimates:   nkd',i6,'   nkr',i6)
      call lgen(plat,r0+radd,nkd,nkdmx,dlv,work)

      if(rank == 0) write(6,342) r0,r0*alat,radd,nkd
  342 format('  r0=',f9.4,'   rc=',f9.4,'   radd=',f9.4,'   nkd=', i7)
      call lgen(qb0,q0+qadd,nkq,nkqmx,qlv,work)
      if(rank == 0) write(6,341) q0,q0*tpiba,qadd,nkq
  341 format('  q0=',f9.4,'   qc=',f9.4,'   qadd=',f9.4,'   nkr=', i7)
      if(allochk) write(*,*) 'deallocate(work)'
      deallocate(work)

C... readin r0c
C     if(newaniso()) then
        eee=screenfac() !takao feb2012
c      elseif(imode==101.or.imode==102) then
c        eee = eees()
C      else
C        eee=0d0
C      endif

!! for eps_lmf and epsPP_lmf mode,
!! even the small eee=1d-4 can affect to dielectric function near q=0 when its values is large as one-hundred or more.
!! Thus we set eee=0d0 to avoid this.
      if(imode==202) then !
         eee=0d0
      endif

      if(rank == 0) write(6,"(' Coulomb is exp(sqrt(-eee)*r)/r. eee=',1pd13.6)")eee
      if(verbose()>=30) call cputid2('mkbj start',0)
C--- bessel and hankel for the expansion of exp(-r/r_0)/r.
c bessel and hankel is renomarized so that its behaves as r^l and r^{-l-1} near r=0.
c  rkpr means r^l*r for e=0 (r0c =infinity) case
      allocate(rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas),rofi(nrx,nbas))
      do ibas=1,nbas
        call genjh(eee,nr(ibas),aa(ibas),bb(ibas),lx(ibas), nrx,lxx,
     o     rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas))
      enddo

C--- onsite integrals <j(e=0)|B> and <B|v(onsite)|B>
cc      if(allochkw) write(*,*) ' allocate  rojb, sgbb '
c      allocate( rojb(nxx, 0:lxx, nbas), sgbb(nxx,  nxx,  0:lxx, nbas))
c      do ibas = 1,nbas
c        call mkjb( lxx, lx(ibas),nxx, nx(0:lxx,ibas),
c     i                 aa(ibas),bb(ibas), nr(ibas), nrx,
c     i                 rprodx(1,1,0,ibas),
c     o        rojb(1,0,ibas), sgbb(1,1,0,ibas))
c      enddo


      allocate( rojb(nxx, 0:lxx, nbas), sgbb(nxx,  nxx,  0:lxx, nbas))
      do ibas = 1,nbas
        call mkjb_4( lxx, lx(ibas),nxx, nx(0:lxx,ibas),
     i                 aa(ibas),bb(ibas), nr(ibas), nrx,
     i                 rprodx(1,1,0,ibas),
     i     rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas),
     o        rojb(1,0,ibas), sgbb(1,1,0,ibas))
      enddo

C#ifdef DEBUG
      print "(' sumcheck rojb',f16.10,'  sgbb',f16.10)", sum(rojb),sum(sgbb)
C#endif

      if(verbose()>=30) call cputid2('mkbj done',0)
c----------------
C--- coulomb matrix for each q = qibz
c----------------
      nlxx= (lxx+1)**2
c      ngb = nbloch + ngcn(1)
      if(rank == 0) then
         allocate(ngvecc0(3,ngcmx))
         call readqg('QGcou',(/0d0,0d0,0d0/),ginv,  quu,ngc0, ngvecc0)
         deallocate(ngvecc0)
      endif
#ifdef USE_MPI
      call MPI_Bcast(ngc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
      ngb = nbloch + ngc0
      iw = idalloc('vcoul',iverbose+2,2*nblochpmx,nblochpmx)
      allocate(vcoul(nblochpmx,nblochpmx) )
      vcoul  = 0d0

C... q near zero
      if(rank == 0) then
        write(6,*) '--- readin Q0P -------'
        open (101,file='Q0P')
        read (101,"(i5)") nq0i
      endif
#ifdef USE_MPI
      call MPI_Bcast(nq0i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
      if(allochk) write(*,*)'allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )'
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      if (rank == 0) then
        read (101,"(d24.16,3x, 3d24.16)") (wqt(i),q0i(1:3,i),i=1,nq0i)
        write (6,"(d13.5,3x, 3d13.5)") (wqt(i),q0i(1:3,i),i=1,nq0i)
        close(101)
      endif
#ifdef USE_MPI
      call MPI_Bcast(wqt,   nq0i, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(q0i, 3*nq0i, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
#endif

      if(rank == 0) write(6,*) ' nqibz nq0i=',nqibz,nq0i

C --- Check PARALELL.X0
c      INQUIRE (FILE = 'PARALELL.X0', EXIST = paralellx0)
C      wvcc=.true.
C      if(newaniso()) wvcc=.false.
      wvcc=.false.
      if(rank == 0) write(6,'(a)') " Mix0vec.XXX is not empty only when"
     &   //" the corresponding q is in Q0P with zero weight."
c      if(paralellx0) then
c        if(wvcc) ifvcfpout = iopen( "VCCFP." //xxt(iqxini,iqxend),0,-1,0)
c        ifgb0vec = iopen ( "Mix0vec."//xxt(iqxini,iqxend),1,3,0)
c        ifgb0vec1 = iopen ( "Mix0vec1."//xxt(iqxini,iqxend),1,3,0)
c      else
      iqxend = nqibz + nq0i

      if(rank == 0) then
        if(wvcc) ifvcfpout = iopen('VCCFP',0,-1,0)
      endif

!ML MPI remark: these files will be read by each q-group
      if(q_comm%ID == 0) then
        ifgb0vec = iopen ( "Mix0vec",1,3,0)
        ifgb0vec1 = iopen ( "Mix0vec1",1,3,0)
      endif

      if(imode==202) then
        iqxini= nqibz + 1
c      elseif(paralellx0) then
c     &  !skip
C      elseif(bzcase()==1) then
C        iqxini = 1 !oct2005
      else
        iqxini = 1
      endif
!!
C      if(newaniso().and.imode==0) then
      if(imode==0) then
        iqxini=1
c        iqxend=nqibz ! comment out at 18nov2012
      endif
      if(rank == 0) print *,'iqxini iqxend=',iqxini,iqxend

c qibz loop
c      epsx =  0.01d0
C      if(bzcase()==1) then
        if(abs(sum(qibz(:,1)**2))/=0d0) call rx('hvccfp0: sum(q**2)==0d0')
C      endif
      if(wvcc.and.(rank == 0)) write(ifvcfpout) nqbz, nblochpmx

c      if(imode==101) then
c        read(ifvcfporg) nqbz_in, nblochpmx_in
c        if(nqbz /= nqbz_in) stop 'nqbz /= nqbz_in VCCFP.ORG'
c        if(nblochpmx /= nblochpmx_in)
c     &     stop 'nblochpmx /= nblochpmx_in VCCFP.ORG'
c      endif



C... Readin PRODMT into prodmt. oct2005
      smbb = smbasis()
      if(rank == 0) write(6,*) ' smooth mixed basis=',smbb

      if(smbasis()) then
        allocate( prodmt(2,nxx,0:lxx,nbas))
        allocate( nx_r(0:lxx))

        do ibas =1,nbas
          if(rank == 0) then
            filenamep = 'PRODMT_'//charnum3(ibas)
            ifprodmt  = iopen(filenamep,0,-1,0)
            read(ifprodmt) nl_r
            read(ifprodmt) nxx_r
            write(6,"(' nxx =',100i3)")nxx_r
          endif
#ifdef USE_MPI
          call MPI_Bcast(nl_r,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
          call MPI_Bcast(nxx_r, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

          if( 2*(nl_r-1) /= lxx ) then
             write(6,*) 2*(nl_r-1),lxx
             call rx('2*nl_r-1 /= lxx ')
          endif

          if(nxx_r>nxx) call rx('nxx_r>nxx')

          if(rank == 0) then
             read(ifprodmt) nx_r(0:lxx)
             write(6,"(' nx_r=',100i3)") nx_r(0:lxx)
             lx_ = lx(ibas)
             if(sum(abs(nx(0:lx_,ibas)-nx_r(0:lx_))) /=0) then
                write(6,*)' debug: nx  =',nx(0:lx_,ibas)
                write(6,*)' debug: nx_r=',nx_r(0:lx_)
                call rx('nx /=nx_r')
             endif
             read(ifprodmt) prodmt(1:2, 1:nxx_r, 0:lxx, ibas)
             write(6,*)' sumcheck prodmt=',sum(abs(prodmt(:,:,:,ibas)))
             isx = iclose(filenamep)
          endif
#ifdef USE_MPI
          call MPI_Bcast(nx_r, lxx, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
          call MPI_Bcast(prodmt, 2*nxx*(lxx+1)*nbas, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
#endif
         enddo

ccccccccccccccccccccccccccccccccccccccc
C... Check write for radial part of the product basis
!ML
!ML MPI remark: this block has not yet been adapted to MPI
!ML
         if(.false.) then
           do ibas= 1,1         !1,nbas
             do l   =  0,lx(ibas)
               open(1011,file='ProdOld_ibas'//charnum3(ibas)//'_l'//charnum3(l))
c         open(2011,file='ProdNew_ibas'//charnum3(ibas)//'_l'//charnum3(l))
               nxdim = nx(l,ibas)
               do ix=1,nxdim
                 write(1011,"(' -- -- -- ',3i3,' --- ' )") ix,l,ibas
c     write(2011,"(' -- -- -- ',3i3,' --- ' )") ix,l,ibas
                 do ir =1,nr(ibas)
                   write(1011,"(d13.5,2x,2d18.8)") rofi(ir,ibas), rprodx(ir,ix,l,ibas),
     .               rprodx(ir,ix,l,ibas) /rofi(ir,ibas)
c     write(2011,"(d13.5,2x,2d18.8)")
c     &   rofi(ir,ibas), sum(rprodx(ir,1:nxdim,l,ibas)*rdmatch(1:nxdim,ix,l,ibas))
c     &   , sum(rprodx(ir,1:nxdim,l,ibas)*rdmatch(1:nxdim,ix,l,ibas))/rofi(ir,ibas)
                 enddo
               enddo
               close(1011)
c     close(2011)
             enddo
           enddo
c     stop 'text end'
         endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
         allocate( rdmatch(nxx,nxx,0:lxx,nbas) )
         do ibas= 1, nbas
           do l   = 0, lx(ibas)
             nxdim = nx(l,ibas)
             if(nxdim<=1)write(6,*) 'hvccfp0:smbasis case error nxdim <=1'
!       pval  = prodmt(1, 1:nxdim, l,ibas)
!       pslo  = prodmt(2, 1:nxdim, l,ibas)
!       prod(r, inew) = \sum_iold rrmat(inew,iold) * prod(r,iold)
             if(rank == 0) write(6,"('goto mkradmatch ibas lnxdim =',3i4)")ibas,l,nxdim
             call mkradmatch(prodmt(1:2, 1:nxdim, l,ibas), nxdim,
     o         rdmatch(1:nxdim,1:nxdim,l,ibas) )
           enddo
         enddo

! index (mx,nx,lx,ibas) ordering: taken from voul_4
         allocate(ibl(-lxx:lxx,nxx,0:lxx,nbas))
         ibl1 = 0
         ibl=999999
         do ibas= 1, nbas
           do l   = 0, lx(ibas)
             do n   = 1, nx(l,ibas)
               do m   = -l, l
                 ibl1  = ibl1 + 1
                 ibl(m,n,l,ibas) = ibl1
!       write(6,*)ibl1,n,l,m,lmbl(ibl1)
               enddo
             enddo
           enddo
         enddo
         if(ibl1/= nbloch) then
           write(6,*)' ibl1 nbloch',ibl1, nbloch
           call rx(' hvccfp0:smbasis mode  error ibl1/= nbloch')
         endif
! index (mx,nx,lx,ibas) ordering
ctttt
         nnr = 2                ! =2 new
        ! =0 equivalence with original mixed basis
         if(rank == 0) write(6,*)' sss:nbas lx=',nbas,lx(1:nbas)

         nbln=0
         do ibas= 1, nbas
           do l   = 0, lx(ibas)
             if(rank == 0) write(6,"('sss: nx=',3i4)") ibas,l,nx(l,ibas)
             if(nx(l,ibas)<=0) cycle
             if(nx(l,ibas)==1) call rx( 'nx(l,ibas) =1')
ccccccccccccccccd
ctttt
c         nnr = 2 ! =2 new
c               ! =0 equivalence with original mixed basis
c         if(l<=3) nnr=0
cccccccccccccccccc
             nbln = nbln + (2*l+1)*(nx(l,ibas)-nnr)
           enddo
         enddo
         allocate( pmat(nbloch+ngcmx, nbln+ngcmx) )
         pmat=0d0
         ibln = 0
         do ibas= 1, nbas
           do l   = 0, lx(ibas)
ccccccccccccccccccc
ctttt
c         nnr = 2 ! =2 new
c               ! =0 equivalence with original mixed basis
c         if(l<=3) nnr=0
cccccccccccccccccc
             do nn  = nnr+1, nx(l,ibas) !nn=1 and nn=2 corresponds to non-zero val sol
               do m   = -l, l
                 ibln = ibln +1
                 nxdim = nx(l,ibas)
                 pmat( ibl(m,1:nxdim,l,ibas), ibln) =  rdmatch(1:nxdim, nn, l,ibas)
ctttt
c          pmat( ibl(m,nn,l,ibas), ibln)
c     &    =  1d0
ccccccccccccccccccccccccccccccccccccccccccccccccccc
               enddo            ! m
             enddo              ! nn
           enddo                ! l
         enddo                  ! ibas

C... Store matting matrix (imatchn,imatcho,pmatch)
         ifpomat   = iopen('POmat',0,-1,0)
c       write(6,*)'ttt= sumchk pmat(b)=',sum(abs(pmat(1:nbloch, 1:nbln)))

       endif                    ! smbasis()


!! === open file Vcoud ===
!! This contains E(\nu,I), given in PRB81,125102

!! == main loop for iqx ==
      write(*,'(" [",I4,"] Ready for q-loop:",F15.5)') rank, cpusec()/60

      iqxstart = iqxini + q_comm%group_ID
      iqxstep  = num_q_groups

      do 1001 iqx = iqxstart, iqxend, iqxstep  ! q in IBZ. avoid q=0 case for iqx=1

        if(q_comm%ID == 0) then
        write (*,'(" do 1001 ip",i3,"  out of ",i3,"  k-points "$)') iqx,iqxend
        call cputid(0)
!       write(0,'("[",I4,"] opening file Vcoud.",I5.5)') rank,iqx
        vcoudfile='Vcoud.'//charnum5(iqx)  !this is closed at the end of do 1001
        ifvcoud = iopen(trim(vcoudfile),0,-1,0)
      endif

        if(iqx > nqibz) then  !       iq = 1
          q  = q0i(:,iqx-nqibz)
c          qq = 0d0
        else                  !       iq = iqx
          q  = qibz(:,iqx)
c          qq = q
        endif
cccccccccccccccccccccccc
c        if(imode==202) then !for iqx>nqibz
c           qq=q
c        endif
cccccccccccccccccccccccc

C        if(.not. newaniso() ) then !this is for fe_epsPP_lmfh_chipm feb2012
C          if(sum(q**2)<1d-12) q=(/1d-4,0d0,0d0/) !takao oct2006
C        endif
!! ==== q+G vector ====
        if(q_comm%ID == 0) call readqg('QGcou',q,ginv,  quu,ngc, ngvecc )
#ifdef USE_MPI
         call MPI_BCast(ngc, 1, MPI_INTEGER, 0, q_comm%communicator, ierror)
         call MPI_BCast(quu, 9, MPI_REAL8, 0, q_comm%communicator, ierror)
         call MPI_BCast(ngvecc, 3*ngc, MPI_INTEGER, 0, q_comm%communicator, ierror)
#endif

        ngb = nbloch + ngc  !it was ngcnn(iq)
        if(q_comm%ID == 0) write(6,'(" iqx q ngc =",i5,3f10.4,i5)') iqx,q,ngc

C        if(newaniso()) then
C          continue
C        elseif(bzcase()==1.and.iqx==1) then
C          goto 1101
C        endif

c        ngc = ngcn(iq)
c        ngvecc(1:3,1:ngc) = ngvecci(1:3,1:ngc,iq)
c        write(6,*)' iq ngc=',iq, ngc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  q test
c       q=(/ 0.09d0,0.09d0,0.09d0/)
c      q = q+(/ 0.01d0,0.01d0,0.01d0/)
c        q=q/4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if(verbose()>=30) call cputid2(' start strxq',0)

C--- strxq structure factor.
        if(q_comm%ID == 0) print "(/' ... make structure constants with screening energy',1pd12.4)", eee
        iw = idalloc('strx',iverbose+2,nlxx*nbas,nlxx*nbas)
        allocate(strx(nlxx,nbas,nlxx,nbas))
        do ibas1 =1,nbas
          do ibas2 =1,nbas
            p = bas(:,ibas2)-bas(:,ibas1)
            phasep =exp(img*2*pi*sum(q*p))
            nlx1 = (lx(ibas1)+1)**2
            nlx2 = (lx(ibas2)+1)**2
            if(allochk) write(*,*) 'allocate( s(nlx1,nlx2))'
            allocate( s(nlx1,nlx2),sd(nlx1,nlx2)) !kino add sd----but sd is dummy
            s = 0
            call strxq(1,eee,q,p,nlx1,nlx2,nlx1,alat,voltot,awald,nkd,nkq,dlv,qlv,cg,indxcg,jcg,
     o        s,sd)
C            print "(' sumcheck strx', 2f18.6)", sum(s)
C            call zprm('strx',s,nlx1,nlx1,nlx2)
            strx(1:nlx1,ibas1,1:nlx2,ibas2) = fpi*s      !!! *phasep
            if(allochk) write(*,*)'deallocate( s )'
            deallocate( s,sd )
          enddo                 ! ibas2
        enddo                   ! ibas1
C#ifdef DEBUG
        print "(' sumcheck strx/1e6', 2f18.6)", sum(strx)/1d6
C#endif

        if(verbose()>=30) call cputid2(' done strxq',0)
ccccccccccccccccccccccccccc
c        strx=0d0
cccccccccccccccccccccccccccc

C--- onsite integrals <j(e=0)|P^(q+G)_L> and <B|v(onsite)|B>
c$$$      if(.true.) then !==New version without sgpp and fouvp allocation June2004=====
        if(allochk) write(*,*)'allocate(rojp,sgpb,fouvb)'
        iw = idalloc('fouvb',iverbose+2,max(2*ngc*nxx,1),2*nlxx*nbas)
        allocate(  rojp(ngc,      nlxx, nbas),
     &            sgpb(ngc, nxx, nlxx, nbas),
     &           fouvb(ngc, nxx, nlxx, nbas))
c     &            sgpp(ngc, ngc, nlxx, nbas),
c     &           fouvp(ngc, ngc, nlxx, nbas) )

!$omp parallel private(ibas)
!$omp do
        do ibas = 1,nbas
          call mkjp_4(q,ngc,ngvecc,alat,qlat,
     i      lxx,lx(ibas),nxx,nx(0:lxx,ibas),
     i      bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas),
     i      nr(ibas),nrx,rprodx(1,1,0,ibas),
     i      eee,rofi(1,ibas),rkpr(1,0,ibas),rkmr(1,0,ibas),
     o      rojp(1,1,ibas),sgpb(1,1,1,ibas),
     o      fouvb(1,1,1,ibas))
        enddo                   ! ibas
!$omp end do
!$omp end parallel

       if(verbose()>=30) call cputid2(' complete onsite integrals',0)

C  ---  Coulomb matrix
       if((verbose()>=30).and.(q_comm%ID == 0)) call cputid2(' start vcoulq',0)
        call vcoulq(q,nbloch,ngc,nbas,lx,lxx,nx,nxx,alat,qlat,voltot,ngvecc,
     i     strx,rojp,rojb,sgbb,sgpb,fouvb,nblochpmx,bas,rmax,eee,aa,bb,nr,nrx,rkpr,rkmr,rofi,
     o     vcoul)

C        if((verbose()>=30).and.(q_comm%ID == 0)) call cputid2(' start vcoulq_4',0)
C        call vcoulq_4(q,nbloch,ngc,nbas,lx,lxx,nx,nxx,alat,qlat,voltot,ngvecc,
C     i     strx,rojp,rojb,sgbb,sgpb,fouvb,nblochpmx,bas,rmax,eee,aa,bb,nr,nrx,rkpr,rkmr,rofi,
C     o     vcoul)

        iw = idalloc('fouvb',iverbose+4,2*ngc*nxx,2*nlxx*nbas)
        iw = idalloc('strx',iverbose+4,nlxx*nbas,nlxx*nbas)
        deallocate(strx,rojp,sgpb,fouvb)

!        write(0,'("[",I4,"] After vcoulq_p4. iqx = ",I4)') rank, iqx

c$$$      else !===old version (allocation of sgpp and fouvp are required) ====
c$$$
c$$$       if(allochk) write(*,*) 'allocate(rojp,sgpb,sgpp,fouvb,fouvp)'
c$$$       allocate( rojp(ngc,      nlxx, nbas),
c$$$     &            sgpb(ngc, nxx, nlxx, nbas),
c$$$     &           fouvb(ngc, nxx, nlxx, nbas),
c$$$     &            sgpp(ngc, ngc, nlxx, nbas),
c$$$     &           fouvp(ngc, ngc, nlxx, nbas) )
c$$$       do ibas = 1,nbas
c$$$          print *,' xxx goto mkjp',ibas
c$$$          call mkjp2(q,ngc, ngvecc, alat, qlat,
c$$$     i      lxx, lx(ibas),nxx, nx(0:lxx,ibas),
c$$$     i      bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas),
c$$$     i      nr(ibas), nrx, rprodx(1,1,0,ibas),
c$$$     o      rojp(1,1,ibas),  sgpb(1,1,1,ibas),
c$$$     o      fouvb(1,1,1,ibas),
c$$$     o      sgpp(1,1,1,ibas),fouvp(1,1,1,ibas) )
c$$$       enddo
c$$$c--- the Coulomb matrix
c$$$       print *,' goto vcoulq'
c$$$       call vcoulq(q, nbloch, ngc,
c$$$     i                  nbas, lx,lxx, nx,nxx,
c$$$     i                  alat, qlat, voltot, ngvecc,
c$$$     i          strx, rojp,rojb, sgbb,sgpb,sgpp, fouvb,fouvp, nblochpmx,
c$$$     o          vcoul)
c$$$       if(allochk)
c$$$     &   write(*,*)'deallocate(strx, rojp,sgpb,sgpp, fouvb,fouvp)'
c$$$       deallocate( strx, rojp,sgpb,sgpp, fouvb,fouvp)
c$$$
c$$$      endif !=============================================================

c----check write
        trwv = 0d0
        do i = 1,nbloch
          trwv = trwv + vcoul(i,i)
        enddo

        if(q_comm%ID == 0) then
          write(6,'(" vcoul trwi=",i6,2d22.14)') iqx,trwv
          write(6,'(" sum vcoul(1:ngb,      1:ngb) ",1p,2d20.12,2x,d20.12)')
     &      sum(vcoul(1:ngb,1:ngb)), sum(abs(vcoul(1:ngb,1:ngb)))
          write(6,'(" sum vcoul(1:nbloch,1:nbloch) ",1p,2d20.12,2x,d20.12)')
     &      sum(vcoul(1:nbloch,1:nbloch)),sum(abs(vcoul(1:nbloch,1:nbloch)))
          write(6,*)
        endif

ccccccccccccccccccccccccccccccc
c      vcoul(:, nbloch+1:ngb)=0d0
c      vcoul(nbloch+1:ngb,:)=0d0
ccccccccccccccccccccccccccccccc

 1101   continue
        ngbo=ngb
C   ... Generate ppmt matrix oct2005
        if(smbasis()) then
          allocate( ppmt(2,(lxx+1)**2,nbas,ngc) )
          ppmt = 0d0
          call mkppmt(alat,plat,qlat, q,
     i    ngc, ngvecc,
     i    rmax, nbas,  bas, lx, lxx,
     o    ppmt) ! ppmt contains value and slove of e(i q+G r) at MT boundaries.
          ! ppmt(2,lmxaa,nbas)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c        write(6,*) 'lxx ppmtsum=',lxx, sum(abs(ppmt))
            if(q_comm%ID == 0) write(6,*) 'nbln ngc',nbln,ngc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

C... Matching matrix pmtch. ppmt and prodmt
          pmat(:, nbln+1:nbln+ngc)=0d0
c        write(6,*) 'sss nbln ngc',nbln,ngc

#ifdef USE_MPI

            allocate( first(1:q_comm%size) )
            allocate( last(1:q_comm%size) )
            allocate( displace(1:q_comm%size) )

            a_coeff = 1d0
            b_coeff = 1d0

!ML MPI remark: looks like one can parallelize the igc loop within the q_comm group

            call thdist1(0,a_coeff,b_coeff,q_comm%size,1,ngc,first,last)

            displace(:) =  (nbloch+ngcmx) * (first(:)-1)
            recvcount(:) = (nbloch+ngcmx) * (last(:) - first(:)+1)

            myID = q_comm%ID+1

! allocate local slice of pmat:

            allocate( pmat_local(nbloch+ngcmx, first(myID):last(myID) ) )

            do igc=first(myID),last(myID)

               pmat_local(nbloch+igc, igc) = 1d0

               do ibas= 1, nbas
                  do l  =  0, lx(ibas)
                     do m  = -l, l

                        pval= ppmt(1, l**2 + l+1 +m, ibas,igc)
                        pslo= ppmt(2, l**2 + l+1 +m, ibas,igc)

                        do n = 1,nx(l,ibas)
                           if(n==1.and.debug) write(6,"('ttt2: ')")

                           pmat_local(ibl(m,n,l,ibas), igc)
     &                          =  rdmatch(n,1,l,ibas) * pval
     &                          +  rdmatch(n,2,l,ibas) * pslo

                           if(debug.and.abs(pmat_local(ibl(m,n,l,ibas), igc))/=0d0)
     &                          write(6,"('ttt2: i1 i2 pmat=',2i5,2d13.5)")
     &                          ibl(m,n,l,ibas), nbln+igc, pmat_local(ibl(m,n,l,ibas), igc)
                        enddo   ! n
                     enddo      ! m
                  enddo         ! l
               enddo            ! ibas

            enddo               ! igc

            deallocate(ppmt)
            nn = nbln  +ngc     ! number for new smooth mixed basis.
            no = nbloch+ngc     ! number for original size of mixed basis.



! pmat(nbloch+ngcmx, nbln+ngcmx) )

            call MPI_AllGatherV(pmat_local,recvcount(myID), MPI_DOUBLE_COMPLEX,
     &           pmat(:,nbln+1),recvcount, displace,  MPI_DOUBLE_COMPLEX, q_comm%communicator, ierror)

            deallocate (first, last, displace, recvcount, pmat_local)

#else
          do igc=1,ngc
            pmat(nbloch+igc, nbln+igc) = 1d0
            do ibas= 1, nbas
              do l  =  0, lx(ibas)
                do m  = -l, l
                  pval= ppmt(1, l**2 + l+1 +m, ibas,igc)
                  pslo= ppmt(2, l**2 + l+1 +m, ibas,igc)
                  do n = 1,nx(l,ibas)
                    if(n==1.and.debug) write(6,"('ttt2: ')")
                    pmat(ibl(m,n,l,ibas), nbln+igc)
     &         =  rdmatch(n,1,l,ibas) * pval
     &         +  rdmatch(n,2,l,ibas) * pslo
                    if(debug.and.abs(pmat(ibl(m,n,l,ibas), nbln+igc))/=0d0)
     &        write(6,"('ttt2: i1 i2 pmat=',2i5,2d13.5)")
     &        ibl(m,n,l,ibas), nbln+igc, pmat(ibl(m,n,l,ibas), nbln+igc)
                  enddo         ! n
                enddo           ! m
              enddo             ! l
            enddo               ! ibas
          enddo                 ! igc

          deallocate(ppmt)
          nn = nbln  +ngc ! number for new smooth mixed basis.
          no = nbloch+ngc ! number for original size of mixed basis.

#endif
          if(debug .and. (q_comm%ID == 0) ) write(6,*) 'end of pmat'

C... oo(no,no). The original overlap matrix.
          allocate( pomat(nn,no) )
          allocate( ppovl(ngc,ngc),oo(no,no))
          call mkppovl2(alat,plat,qlat,
     i          ngc,  ngvecc,
     i          ngc,  ngvecc,
     i          nbas, rmax, bas,
     o          ppovl)
          oo = 0d0
          do ipl1 = 1,nbloch
            oo(ipl1,ipl1) = 1d0
          enddo
          do ix= 1,ngc
            do iy= 1,ngc
              oo(nbloch+ix, nbloch+iy) = ppovl(ix,iy)
            enddo
          enddo
          if(debug) write(6,*) 'end of oo'


C... oon(nn,nn) is the overlap matrix with new basis
          allocate(oon(nn,nn))
          oon = matmul( dconjg(transpose(pmat(1:no,1:nn)))
     &               ,matmul(oo,pmat(1:no,1:nn)) )


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Reduction of pmat by SVD ... not meaningful
c      nnmx = 1000000
c      epsmx= 1D20
cc      print *,' sumchk pmat=',sum(abs(pmat(1:no,1:nn)))
c      call  zgesvdnn2(nn,nn, nnmx,epsmx,
c     o   oon, ! pmat is reduced to pmat(1:no,1:nnn) by SVD.
c     o   ngcnn)
cc      call  zgesvdnn2(no,ngc, nnmx,epsmx,
cc     i   pmat(1:no,nbln+1:nbln+ngc), ! pmat is reduced to pmat(1:no,1:nnn) by SVD.
cc     o   ngcnn)
cc      nn= nbln+ngcnn
c      write(6,*)' svd ngc ngcnn=',ngc, ngcnn
c      stop 'test end xxxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccc
          if(.false.) then
          if(q_comm%ID == 0) then
            open(3011,file='oontest'//charnum5(iqx))
            do ix=nbln+1,nn
              igc=ix-nbln
              qqx(1:3) = (q(1:3)+ matmul(qlat, ngvecc(1:3,igc)))
              absqq  =  sqrt(sum(qqx(1:3)**2))
c           absqg2x(ix) =sum( (2*pi/alat *q0i(1:3,nq0i))**2)
              do iy=nbln+ 1,nn
                igc2=iy-nbln
                if(ix==iy) then
                  write(3011,"('on : ',2i8,3i3,2x,3i3,f13.5,3x,2f20.10)")
     &       ix,iy, ngvecc(1:3,igc),ngvecc(1:3,igc2),absqq, oon(ix,iy)
                else
                  write(3011,"('off:', 2i8,3i3, 2f20.10)")ix,iy,
     &       ngvecc(1:3,igc)-ngvecc(1:3,igc2),  oon(ix,iy)
                endif
              enddo
            enddo
            close(3011)
          endif
          endif

ccccccccccccccccccccccccccccccc


C... Generat pomat
!    zmelt_new(K, ij) = \sum_I pomat(K,I)* zmelt(I, ij)
!    means <psi_i psi_j | K> where |K> denote new mixed basis.
!  See sxcf_fal2 and x0kf.
!   Be careful its transpose procedure---it is a little confusing...
          call pmatorth(oo,oon, pmat(1:no,1:nn), no, nn,
     o    pomat)

ccccccccccccccccccccccccccccccccccccccccccccc
ctttt
c        pomat=0d0
c        do ix= 1,ngb
c          pomat(ix,ix)=1d0
c        enddo
c        do ix= 1,ngb
c        do iy= 1,ngb
c          if(pmat(ix,iy)/=0d0 )
c     &   write(6,"(' ttt: pmat=',2i3,2d13.6)")
c     &   ix,iy,pmat(ix,iy)
c        enddo
c        enddo
c        write(6,"(' ttt:sumchk=',2d13.6,2i4)")
c     &     sum(pomat(:,:)), no,nn
cccccccccccccccccccccccccccccccccccccccccccccc

          if( iqx <= nqibz ) deallocate(oon)
          deallocate(ppovl,oo)
C... Store matching matrix

          if(q_comm%ID == 0) then
            write(ifpomat) q,nn,no,iqx
            write(ifpomat) pomat
          endif
          deallocate(pomat)
        endif

C        if(newaniso()) then
C          continue
C        elseif(bzcase()==1.and.iqx==1)then
C          cycle
C        endif

!! == Write out VCCFP ==
        if(debug .and. (q_comm%ID == 0) ) write(6,*) 'write out vcoul'
        if(smbasis()) then
          ngb= nn
          allocate(vcoulnn(ngb,ngb))
          vcoulnn= matmul(transpose(dconjg(pmat(1:no,1:nn)))
     &               ,matmul(vcoul(1:no,1:no),pmat(1:no,1:nn)))
          vcoul(1:ngb,1:ngb)= vcoulnn
          deallocate(vcoulnn)
        endif
         if(wvcc .and. (q_comm%ID == 0) ) then
          write(ifvcfpout) ngb
          write(ifvcfpout) vcoul(1:ngb,1:ngb),q
        endif
         if(q_comm%ID == 0) write(6,"(' ngc ngb/ngbo=',6i6)") ngc,ngb,ngbo

c Mix0vec ---------------------------------
!! diagonalize the Coulomb matrix
        if(verbose()>=30) call cputid2(' start diagonalizing vcoul',0)

c       if( iqx > nqibz .or. iqx==1) then !feb2012 add iqx==1 for newansio()=T
        if(.true.) then
C       if(allochk) write(*,*) 'allocate( ppovl(ngc,ngc))'
        iw = idalloc('coulomb',iverbose+2,2*ngb,5*ngb) ! for oo, oox, vcoul0, hh, zz
        allocate( oo(ngb,ngb) )
        allocate( ppovl(ngc,ngc) )
        call mkppovl2(alat,plat,qlat,
     &        ngc,  ngvecc,
     &        ngc,  ngvecc,
     &        nbas, rmax, bas,
     o        ppovl)
C#ifdef DEBUG
        print "(' sumcheck ppovl', 2f18.6)", sum(ppovl)
C#endif
C       call zprm('ppovl',ppovl,ngc,ngc,ngc)

        if(smbasis()) then
          oo = oon
          deallocate(oon)
        else
          oo = 0d0
          do ipl1=1,nbloch
            oo(ipl1,ipl1) = 1d0
          enddo
          do ix=1,ngc
            do iy=1,ngc
              oo(nbloch+ix, nbloch+iy) = ppovl(ix,iy)
            enddo
          enddo
        endif

        allocate( oox(ngb,ngb) )
        oox = oo
        allocate(  vcoul0(ngb,ngb) )
        vcoul0 = vcoul(1:ngb,1:ngb)
        if(allochk)
     &  write(*,*) 'allocate(hh(ngb,ngb),oo(ngb,ngb),oox,zz,eb,zzr)'
        allocate(hh(ngb,ngb),zz(ngb,ngb),eb(ngb),zzr(ngb))
        hh  = - vcoul0
c        nmx = 15
        nmx = ngb
C       call zprm('ovlp',oo,ngb,ngb,ngb)
C       call zprm('vcoul',hh/1d5,ngb,ngb,ngb)
        call diagcv(oo,hh,zz,ngb, eb,nmx,1d99,nev)
C       call prmx('evals of zcoul',-eb,nev,nev,1)
C       call zprm('zcoul',zz,ngb,ngb,ngb)

C       Printout
        if(q_comm%ID == 0) then
          write(*,"(' ... eigenvalues of Coulomb matrix, q=',3f10.6,' dimension =',i5)") q,ngb
          do ipl1=1,nev
            if(ipl1==11) write(6,*)' ... '
            if(ipl1>10.and.ipl1<nev-5) cycle
            write(6,'(i5, 1pd20.12)') ipl1, -eb(ipl1)
          enddo
        endif

C       Zero out positive eigenvalues
        if(q_comm%ID == 0) then
          write(*,"(' ... eigenvalues of Coulomb matrix, q=',3f10.6,' dimension =',i5)") q,ngb
          do  ipl1 = 1, nev
            if (-eb(ipl1) > toleb) cycle
            maxposeb = max(maxposeb,eb(ipl1))
            nposeb = nposeb+1
            eb(ipl1) = -toleb
          enddo
        endif


c$$$!! Modify -eb
c$$$          if(iqx<=nqibz) then
c$$$             do iqbz=1,nqbz     !! check
c$$$                if(sum(abs(qbzwww(:,iqbz)-q))<1d-6) then
c$$$                   iqbzx=iqbz
c$$$                   goto 888
c$$$                endif
c$$$             enddo
c$$$             stop 'hvccfp0:sum(abs(qbzwww(:,iq)-qbz(:,iq)))>1d-6'
c$$$ 888         continue
c$$$             if(abs((-eb(1)+eb(2))/eb(2))<1d-2) then
c$$$!! Center. touching case. Respect smoothness when we change n1n2n3 division.
c$$$                eb(1)=eb(1)*wqfac(iqbzx)
c$$$                eb(2)=eb(2)*wqfac(iqbzx)
c$$$             else
c$$$                eb(1)=eb(1)*wqfac(iqbzx)
c$$$             endif
c$$$          endif

!! === save zz === apr2012takao
c          if( newaniso().and.iqx==1 ) then
c            if(sum(q**2)>1d-10) then
c               stop ' hvccfp0: sanity check. |q(iqx)| /= 0'
c            endif


!       write(0,'("[",I4,"] Before writing. iqx = ",I4)') rank, iqx
        if(q_comm%ID == 0) then
          write(*,'("[",I4,"] Writing to ",A16)') rank,vcoudfile
          write(ifvcoud) ngb
          write(ifvcoud) q
          write(ifvcoud) -eb
          write(ifvcoud) zz
        endif

c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$        print *,' dddddddddddddddd q=',q
c$$$        do ix=1,ngb
c$$$        do iy=1,ngb
c$$$          aaaa=  sum(  dconjg(zz(1:ngb,ix))*matmul( oox,zz(1:ngb,iy))  )
c$$$           if(ix==iy .and.  abs(aaaa-1d0) >1d-8 ) then
c$$$            write(*,*)' dddd zcousum check',ix,iy,aaaa
c$$$           endif
c$$$           if(ix/=iy .and.  abs(aaaa) >1d-8 ) then
c$$$            write(*,*)' dddd zcousum check',ix,iy,aaaa
c$$$           endif
c$$$        enddo
c$$$        enddo
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if(q_comm%ID == 0) then
          print *
          write(6,'(" eig0 must be equal to the largest =", 1p2d20.12)')
     &       sum(  dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1))  )
          write(6,'(" zz norm check=",1pd20.12)')
     &    sum( dconjg(zz(1:ngb,1))*matmul(oox,zz(1:ngb,1)) )
          print *
c          write(6,'(" --- vcoul(exact  no eee)=",1pd14.6," absq2=",1pd20.12)')
c     &    fpi*voltot/(sum(tpiba**2*q(1:3)**2))
c     &             , (sum(tpiba**2*q(1:3)**2))
          if (all(q /= 0)) write(6,'(" --- vcoul(exact)=",1pd14.6," absq2=",1pd20.12)')
     &    fpi*voltot/(sum(tpiba**2*q(1:3)**2)-eee)
     &             , (sum(tpiba**2*q(1:3)**2)-eee)
          write(6,'(" --- vcoul(cal ) =",1p2d14.6)')
     &    sum( dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1)) )*voltot
        endif

ccccccccccccccccccccccccccccccccccccccccc
c          do igc=1,ngb
c          qqx(1:3) = (q(1:3)+ matmul(qlat, ngvecc(1:3,igc)))
c          write(6,'(" --- vcoul(exact) xxx =",1pd14.6," absq2=",1pd20.12)')
c     &    fpi*voltot/(sum(tpiba**2*(qqx(1:3)**2)-eee))
c     &             , (sum(tpiba**2*(qqx(1:3)**2)-eee))
c          write(6,'(" --- vcoul(cal ) xxx =",1p2d14.6)')
c     &    sum( dconjg(zz(1:ngb,igc))*matmul( vcoul0,zz(1:ngb,igc)) )*voltot
c          enddo
cccccccccccccccccccccccccccccccccccccccccc
        deallocate(vcoul0)

        if( iqx-nqibz>=1 ) then
          if( wqt(iqx-nqibz)==0d0) then ! MIZUHO-IR

C --- To get the vector <Mixed basis| q=0> --------------
            if(.not.is_mix0vec()) then !used original befor oct2006
! See switch.F ---> this is not used now.
              ifgb0vec_a =ifgb0vec1
              ifgb0vec_b =ifgb0vec
            else
! ismix0vec=1 is to avoid problem at BZ boundary when is_mix0vec()=0.
              ifgb0vec_a =ifgb0vec
              ifgb0vec_b =ifgb0vec1
            endif
c1... Case1 to write ifgb0vec -------------------------------------
            if(q_comm%ID == 0) write(6,*) 'cell volume =',voltot

C      this would otherwise not be initialized.
            igc0 = 0


            if(ngc==0) then
              continue
            else
              do igc=1,ngc
                if( sum(abs( ngvecc(1:3,igc) ))==0 ) then
                  igc0=igc
                  exit
                endif
              enddo
              if(q_comm%ID == 0) print *,' igc0=',igc0,ngvecc(1:3,igc0)
              zzr(nbloch+1:nbloch+ngc) = ppovl(1:ngc,igc0)
            endif

            allocate( gbvec(ngb), b0mat(nbloch) )

C ... get a vector <Product Basis| q+0>
            call mkb0(q,lxx,lx,nxx,nx,aa,bb,nr,nrx,rprodx,alat,bas,nbas,nbloch,b0mat)
            zzr(1:nbloch) = b0mat(1:nbloch)
cccccccccccccccccccccccccccc
c           do igc=1,ngb
c             write(6,"('ssss: ',i5,2d14.6)") igc, zzr(igc)
c           enddo
cccccccccccccccccccccccccccc
            allocate(ooxi(ngb,ngb))
            ooxi=oox
            call matcinv(ngb,ooxi)
            gbvec = matmul(ooxi, zzr)

cccccccccccccc
c           do igc=1,ngb
c             write(6,"('ssss: ',i5,2d14.6)") igc, gbvec(igc)
c           enddo
ccccccccccccccc
            deallocate(ooxi)
            dnorm = sqrt( sum(dconjg(gbvec)*zzr) )
! remove /dnorm at 14June2008. See main/hx0fp0.
! dnorm corresponds to volume (or sum of MT volume if no IPW).
c            gbvec = gbvec /dnorm
c            zzr   = zzr   /dnorm
! Not dnorm=1 at 14June2008. See main/hx0fp0.
c            dnorm=1

            if(q_comm%ID == 0) then
              write(ifgb0vec_a,"(3d24.16,2i10,d24.16)") q, ngb,igc0,dnorm
              write(ifgb0vec_a,"(4d24.16)") (gbvec(i),zzr(i),i=1,ngb)
            endif
            deallocate(gbvec,b0mat)
c1----------------------------------------------------

c2... --- Case2 to write ifgb0vec c2 is problematic at BZ boundary...------
            dnorm  = 1d0
            zzr(:) = matmul (oox, zz(:,1))
            igc0 = 999999 !dummy now
c phasex ---just to clean. this is irrelevant
            phasex =1d0
            do i=1,ngb
              if(abs(zz(i,1)) > 1d-3) phasex = abs(zz(i,1))/zz(i,1)
            enddo
            do i=1,ngb
              zz(i,1)= phasex * zz(i,1)
              zzr(i) = phasex * zzr(i)
            enddo
            if(q_comm%ID == 0) then
              write (ifgb0vec_b,"(3d24.16,2i10,d24.16)") q, ngb,igc0,dnorm
              write (ifgb0vec_b,"(4d24.16)") (zz(i,1),zzr(i),i=1,ngb)
            endif
          endif                 ! wgt==0
        endif                   ! (iqx>nqibz)              ! ML


        if (verbose()>=40) call cputid2(' complete diagonalizing vcoulq',0)
        iw = idalloc('coulomb',iverbose+4,2*ngb,5*ngb)  ! for oo, oox, vcoul0, hh, zz
        deallocate(hh,oo,zz,eb,oox,zzr)
        deallocate(ppovl)
        endif
        idummy=iclose(trim(vcoudfile))

C        print *, '!! ... abort q loop for debugging'
C        goto 999

 1001 continue ! iqx;  end of q-loop !

C  999 continue


      call mpi_end()

      iw = idalloc('vcoul',iverbose+4,2*nblochpmx,nblochpmx)
      deallocate(vcoul)

      deallocate(ngvecc)
      call cputid(0)

      if (nposeb > 0) then
        print 132, nposeb, -maxposeb
  132   format(/' hvccfp0 (warning): ',i4,' negative eigenvalues found.  Min val = ',1pg10.3)
      endif

      if (verbose()>=40) then
        iw = idalloc(' ',11,1,1)
      else
        iw = idalloc(' ',1,1,1)
        write(*,"(/' Max dynamic allocation hvccfp0:',i7,' MB'/)") iw
      endif

      if (nthreads > 1) then
        print 131, (dwtime() - walltime0)/60, cpusec()/60, nthreads
      else
        print 131, (dwtime() - walltime0)/60, cpusec()/60
      endif
  131 format(' hvccfp0 :  Wall clock time',f9.1,'m','  Processor time',f9.1,'m':,i4,' threads')

      if(imode==202) call rx0( ' OK! hvccfp0 imode=202 only for Q0P')
c      if(imode==101) stop ' OK! hvccfp0 imode=101 overwirte VCCFP'
c      if(imode==102) stop ' OK! hvccfp0 imode=102 remove_r0c'
      if(imode==0) call rx0( ' OK! hvccfp0 imode=0')
      if(imode==3) call rx0( ' OK! hvccfp0 imode=3')

      end

      subroutine checkagree(a,b,char)
      real(8):: a(3),b(3)
      character*(*) :: char
      if(sum(abs(a-b))>1d-6) then
        call rxs('Error in checkagree:',char)
      endif
      end

      subroutine mkradmatch(p,nxdim,rdmatch)
C- make rdmatch
C----------------------------------------------------
Ci  p(1,i): phi     at mt for i-th basis
Ci  p(2,i): dphi/dr at mt for i-th basis
Co rdmatch(nxdim,nxdim)
C-------
Cr    phinew_j(r) =sum_i phi_i(r)* rdmatch (i,j)
Cr     phinew_1(rmt)    =1      phinew_2(rmt)   =0
Cr   d phinew_1(rmt)/dr =0    d phinew_2(rmt)/dr=1
Cr for k >=3
Cr     phinew_k(rmt)    =0
Cr   d phinew_k(rmt)/dr =0
C----------------------------------------------------
      implicit none
      integer(4):: nxdim,lbas,i,i1,i2,ix
      real(8):: p(1:2, 1:nxdim), rdmatch(1:nxdim,1:nxdim)
      real(8):: pd,p1,p1d,p2,p2d,s,t, eps=1d-3,delta
Cr                                       old     new
c      write(6,"('mkradmatch: nxdim=',i4)") nxdim
      if(nxdim <=0) return
      if(nxdim ==1) call rx( 'mkradmatch err nxdim==1')
      rdmatch=0d0
C... pivot--- get better set of phi for augmentation
      do
        i1= nxdim
        i2= nxdim-1
        p1 = p(1, i1)
        p2 = p(1, i2)
        p1d= p(2, i1)
        p2d= p(2, i2)
        write(6,"('mkradmatch: i1 p1 p1d=',i3,2d13.6)") i1,p1,p1d
        write(6,"('mkradmatch: i2 p2 p2d=',i3,2d13.6)") i2,p2,p2d
        delta = p1*p2d-p2*p1d
        if(abs(delta) <eps*p1*p2) then
          if(i2==1) then
            write(6,"(' i1 i2=',2i5,2d13.6)") i1,i2,p1d/p1,p2d/p2
            call rx( 'mkradmatch: err poor linear dep')
          endif
          i2=i2-1
        endif
        exit
      enddo
C...
      call phimatch(1d0,0d0,  p1,p1d,p2,p2d, s,t)
      rdmatch(i1, 1)=  s
      rdmatch(i2, 1)=  t
      write(6,"('mkradmatch: 1 0    st=',2d13.5)") s,t
      call phimatch(0d0,1d0,  p1,p1d,p2,p2d, s,t)
      rdmatch(i1, 2)=  s
      rdmatch(i2, 2)=  t
      write(6,"('mkradmatch: 0 1    st=',2d13.5)") s,t

      ix=2
      do i= 1,nxdim
        if(i==i1.or.i==i2) cycle
        ix=ix+1
c        write(6,"('mkradmatch: i p pd=',i3,2d13.5)") i,p(1,i),p(2,i)
        call phimatch(p(1,i),p(2,i),  p1,p1d,p2,p2d, s,t)
        rdmatch(i,  ix)=  1d0
        rdmatch(i1, ix)=  -s
        rdmatch(i2, ix)=  -t
        write(6,"('mkradmatch: ix st=',i3,2d13.5)") ix,s,t
      enddo
      end

      subroutine phimatch(p,pd, p1,p1d,p2,p2d, s,t)
C --- match for given p and pd
c   phi = s phi1 + t phi2 !slope and value are at MT
c     p  = s p1  + t p2
c     pd = s pd1 + t pd2
      implicit none
      real(8):: matinv(2,2),p,pd,p1,p1d,p2,p2d,s,t,delta,ddd1,ddd2
      delta = p1*p2d-p2*p1d
      matinv(1,1) = 1/delta *  p2d
      matinv(1,2) = 1/delta * (-p2)
      matinv(2,1) = 1/delta * (-p1d)
      matinv(2,2) = 1/delta *  p1
      s = matinv(1,1) *p  + matinv(1,2) *pd
      t = matinv(2,1) *p  + matinv(2,2) *pd
C... check
      ddd1 = abs(s*p1  + t*p2   -  p )
      if(  ddd1 >1d-8 ) call rx( 'phimatch: ddd1 err')
      ddd2 = abs(s*p1d + t*p2d  -  pd)
      if(  ddd2 >1d-8 ) call rx( 'phimatch: ddd2 err')
      end

      subroutine pmatorth(oo,oon,pmat,no,nn, pomat)
C get conversion matrix from old mixed basis(no) to augmented mixed basis(nn).
C pmatorth contains
c   oo^{-1}_IJ
      implicit none
      integer(4):: no,nn,io,in,i
      complex(8):: pmat(no,nn),pomat(nn,no),oo(no,no),oon(nn,nn)
      complex(8),allocatable:: ooninv(:,:)
      real(8),allocatable:: eb(:)
      allocate(ooninv(nn,nn))
      ooninv = oon
      call matcinv(nn,ooninv) !generate ooninv
c      pomat = matmul(ooninv, matmul(dconjg(transpose(pmat)),oo))
      pomat = transpose (matmul( oo, matmul(pmat,ooninv)))
      deallocate(ooninv)
      end
c      allocate(pp(nn,nn),ppin(nn,nn),eb(nn),zz(nn,nn),zze(nn,nn))
c      ppin = pp
c      call diagcvh(ppin,nn,eb,zz)
c      do i=1,nn
c        zze(:,i) =  zz(:,i)* sqrt(eb(i))
c      enddo
c      pomat = matmul(pmat, matmul(zze,dconjg(transpose(zz))))

      subroutine diagcvh(hh,ngb,eb,zz)
      implicit none
      integer(4):: nmx,nev,i,ngb
      complex(8):: hh(ngb,ngb),oo(ngb,ngb),zz(ngb,ngb)
      real(8):: eb(ngb)
      nmx=ngb
      oo = 0d0
      do i=1,ngb
        oo(i,i) = 1d0
      enddo
      call diagcv(oo,hh,zz,ngb, eb,nmx,1d99,nev)
      print *,' diagcvv: ngb,nev=',ngb,nev
      do i=1,nev
        write(6,'(i4,d23.16)')i, eb(i)
      enddo
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zgesvdnn2(no,nn, nnmx,epsmx,
     i   pmat,
     o   nnn)
c pmat(no,nn) ---> pmat(no,nnn)
Cio input          pmat(no,nn)
Cio output reduced pmat(no,nnn)
      implicit none
      integer(4):: lwork,info,nn,no,nnn,nnmx,i
      complex(8)::  pmat(no,nn),uu(no,no),vt(nn,nn)
      real(8):: ss(nn),epsmx
      real(8),allocatable:: rwork(:)
      complex(8),allocatable:: work(:),vtt(:,:),pmatx(:,:)
c      print *,' sumchk pmat=',sum(abs(pmat(1:no,1:nn)))
      lwork=4*no
      allocate(work(LWORK),rwork(5*no),pmatx(no,nn))
      pmatx =pmat
      call zgesvd('A','A',no,nn,pmat,no,SS,UU,no,VT,nn,work,lwork,rwork,info)
      nnn=-999
      do i=1,nn
        write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
!         vtt(i,:)=ss(i)*vt(i,:)
        if(nnn==-999.and.ss(i)<epsmx) nnn = i-1
      enddo
c      write(6,*) 'nnn=',nnn
      if(nnn==0) call rx( 'strange: nnn=0')
      if(nnn>nnmx) nnn=nnmx
      pmat=pmatx
c      pmat(:,1:nnn) = uu(:,1:nnn)
!      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
!     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
!      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk)))
!     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
!      deallocate(vtt)
      end


c---------------------------------------------------------------------
      subroutine mkb0( q, lxx,lx,nxx,nx, aa,bb, nrr,nrx,rprodx,
     i        alat,bas,nbas,nbloch,
     o        b0mat)
C--make the matrix elementes < B_q | exp(iq r)>
      implicit none
      integer(4) :: nlx,l,n,m,nr,ir,lm,ibl1,ibas,nrx,nbloch,jobBes

      integer(4) :: nbas,lxx, lx(nbas), nxx, nx(0:lxx,nbas),nrr(nbas)
      real(8)    :: rprodx(nrx,nxx,0:lxx,nbas),aa(nbas),bb(nbas),
     &   phi(0:lxx),psi(0:lxx), bas(3,nbas),
     &   alat,
     &   pi,fpi,tpiba,qg1(3),q(3),absqg,r2s,a,b
c
      complex(8) :: b0mat(nbloch),img=(0d0,1d0) ,phase
c
      integer(4),allocatable:: ibasbl(:), nbl(:), lbl(:), lmbl(:)
      real(8),allocatable :: ajr(:,:),rofi(:),rob0(:,:,:)
      real(8),allocatable::cy(:),yl(:)
      complex(8),allocatable :: pjyl(:,:)
      procedure(integer) :: GWversion
#ifdef COMMONLL
      integer(4) ll(51**2)
      common/llblock/ll
#else
      integer(4) ll
#endif

c-----
      write(6,*)'mkb0:'
      pi   = 4d0*datan(1d0)
      fpi  = 4*pi
      nlx  = (lxx+1)**2
      jobBes = mod(GWversion()/1000,10)
c
      tpiba = 2*pi/alat
      qg1(1:3) = tpiba * q(1:3)
      absqg    = sqrt(sum(qg1(1:3)**2))
c
      allocate(ajr(1:nrx,0:lxx), pjyl(nlx,nbas),rofi(nrx),
     &  ibasbl(nbloch), nbl(nbloch), lbl(nbloch), lmbl(nbloch),
     &  cy(nlx),yl(nlx),rob0(nxx,0:lxx,nbas))
c
      call sylmnc(cy,lxx)
!spherical factor Y( q+G )
      if (absqg /= 0) then
        call sylm( qg1/absqg,yl,lxx,r2s)
      else
        call sylm( qg1,yl,lxx,r2s)
      end if
c
      do ibas = 1,nbas
        a = aa(ibas)
        b = bb(ibas)
        nr= nrr(ibas)
        rofi(1)    = 0d0
        do ir      = 1, nr
          rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
C         call bessl(absqg**2*rofi(ir)**2,lx(ibas),phi,psi)
          call besslr(absqg**2*rofi(ir)**2,10*jobBes,0,lx(ibas),phi,psi)
          do l  = 0,lx(ibas)
c ... bessel function
            ajr(ir,l) = phi(l)* rofi(ir) **(l +1 )
            ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
          enddo
        enddo

c ... Coefficients for j_l yl  on MT  in the expantion of of exp(i q r).
        phase = exp( img*sum(qg1(1:3)*bas(1:3,ibas))*alat  )
        do lm = 1,(lx(ibas)+1)**2
          l = ll(lm)
          pjyl(lm,ibas) = fpi *img**l *cy(lm)*yl(lm) *phase  *absqg**l
        enddo
c ... rob0
        do l = 0,lx(ibas)
          do n = 1,nx(l,ibas)
            call gintxx( ajr(1,l), rprodx(1,n,l,ibas), a,b,nr,
     o                 rob0(n,l,ibas) )
          enddo
        enddo
      enddo

c ... index (mx,nx,lx,ibas) order.
      ibl1 = 0
      do ibas= 1, nbas
        do l   = 0, lx(ibas) ! write(6,'(" l ibas nx =",3i5)') l,nx(l,ibas),ibas
          do n   = 1, nx(l,ibas)
            do m   = -l, l
              ibl1  = ibl1 + 1
              ibasbl(ibl1) = ibas
              nbl   (ibl1) = n
              lbl   (ibl1) = l
              lmbl  (ibl1) = l**2 + l+1 +m ! write(6,*)ibl1,n,l,m,lmbl(ibl1)
            enddo
          enddo
        enddo
      enddo
c ... pjyl * rob0
      do ibl1= 1, nbloch
        ibas= ibasbl(ibl1)
        n   = nbl  (ibl1)
        l   = lbl  (ibl1)
        lm  = lmbl (ibl1)
        b0mat(ibl1) = pjyl(lm,ibas) * rob0(n,l,ibas)
      enddo
      deallocate(ajr, pjyl,rofi,
     &  ibasbl, nbl, lbl, lmbl,
     &  cy,yl,rob0)
      end
