      program h_uumatrix

c-------------------------------------------------
c  Calculate <u|u> matrix . u_kj(r) is the perodic part of eigencuntion.
c
c ixc=1: <u(k) | u(q_in)>
c ixc=2: <u(k) | u(k+b)>
c ixc=3: <u(k) | u(k+q0)>
c ixc=4: <u(k) | u(k+b)> for 1D (k=unshifted, shidted, or unshifted+q0)
c
c modified by Takahi Miyake, Aug, 2005 (ixc=4)
c modified by Takahi Miyake, Jul, 2005 (ixc=3)
c modified by Takahi Miyake, May, 2004
c originally written by Takao Kotani, April, 2004
c-------------------------------------------------
      use m_readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,readcphi,readgeig,
     & init_readeigen_phi_noeval,init_readeigen_mlw_noeval1d,init_readeigen_mlw_noeval
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue

      use mpi_mod

      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     & , ngvecpf1(:,:), ngvecpf2(:,:),
     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:),
     &phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:),
     &qbzs(:,:),qbz2(:,:)
      complex(8),parameter:: img=(0d0,1d0)
c,nw,incwf,natom,nclass,ipos,igrp,
c     & iinvg,nspin,nl,nn,nnv,nnc,
c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
c     o                   izdummy,
c     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
c     o   iecore,ikonf,iicore,incore,nctot,             !core
c     o   imagw_dummy,niw,idummy,
      integer(4)
     &   nw_input,
     &   ifhbe,
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
     &   ifepscond,nxx,ifvxcpout,ifgb0vec
     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax  !,ifev1,ifev2
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &   ,ngpmx

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini,
     & dwry,dwh,omg_c,omg2,xxx
      integer(4)::nwin, incwfin,  verbose
      real(8)::efin
      integer(4):: bzcase, mrecg,ifphi,
     & nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc,
     & ldim2,ixx,ngp1,ngp2,nq0it
      real(8):: qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
      real(8),parameter::  pi =     3.14159265358979323846d0
      real(8),parameter::  fpi =    4d0*pi

c      logical:: test_qbzreg
      logical:: qbzreg
!-------------------------------------------------------------------------
      integer(4),allocatable:: ncindx(:,:),
     &           lcindx(:,:),
     &           nrad(:),
     &           nindx_r(:,:),
     &           lindx_r(:,:),
     &           nc_max(:,:),
     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
      real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:)
     &                     ,phitotr(:,:,:,:,:),
     &        bb(:),zz(:),rmax(:),cy(:),yl(:)


      complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:)
     & ,uum(:,:,:), ppovl(:,:)
      complex(8):: ppj,phaseatom
      real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
      integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin
     & ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx
     &, iopen,ico,lxd,lx !,ll
      real(8):: ylk
#ifdef COMMONLL
      integer(4):: ll(51**2)
      common/llblock/ll
#else
      integer(4) ll
#endif
ctm
      integer(4) :: ixc,idummy,idummy2,i1,i2,i3,nbbloop,
     &              ifq0p,ifuu(2),
     &              ifbb,nbb,iko_ixs(2),iko_fxs(2),noxs(2),
     &              iqibz,iqbz,ibb,itmp,itmp2,iti,itf,
     &              nqibz2,nqbz2,iqb,ibb2,iqtmp,ibbtmp,ndg(3),
     &              nb1d,iq0i
      integer(4),allocatable:: ikidx(:),ikbidx(:,:),
     &              ibidx(:,:),ibidxs(:,:),ibidx0(:,:,:),
     &              ij1idx(:),ij2idx(:)
      real(8),allocatable :: bbv(:,:)
      real(8) :: uunorm,dqx(3),dqx0(3),dq0(3),dg(3),dqmin(3),adq0,
     &           q0wf(3),wgt
      integer(4) :: iqindx
c-------------------------


      stop 'huumat is out of date'
C
C
C      call mpi_start()
C
C      if(mpi_size > 1) then
C         call mpi_end()
C         stop 'MPI not yet implemented.'
C      endif
C
Cctm
Cc mode switch. --------------
C      write(6,*) ' --- Choose omodes below -------------------'
C      write(6,*) '  (1) (q,Qin),  (2) (q,q+b), (3) (q,q+q0)'
C      write(6,*) '  or (4) (q,q+b) for 1D'
C      write(6,*) ' --- Put number above ! --------------------'
C      call readin5(ixc,idummy,idummy2)
C      write(6,*) ' ixc=',ixc
C      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
C
C
Cc---  readin BZDATA. See gwsrc/rwbzdata.f
Cc--------readin data set when you call read_BZDATA ---------------
Cc       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
Cc     &   n_index_qbz
Cc       integer(4):: n1,n2,n3
Cc       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
Cc       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
Cc     &    ,wibz(:),qbzw(:,:)
Cc       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
Cc     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
Cc-----------------------------------------------------------------
C      call read_BZDATA()
C
Cc--- Use regular mesh even for bzcase==2
C      allocate(qbzs(3,nqbz))
C      if(bzcase()==2.and.qbzreg()) then
C        deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
C        do i=1,nqbz
C          qbzs(:,i) = qbz(:,i) -deltaq/2d0
C          write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
C        enddo
Cctm 050802
C        if (ixc <= 4) stop "huumat: bzcase=2 is not implemented."
C      endif
C      print *,' ======== nqbz qbz  =',nqbz
C      print *,' ======== nqibz ngrp=',nqibz,ngrp2
C      print *,  qbz
C      print *,'============================'
C      print *
C
CC--- readin GWIN and LMTO, then allocate and set datas.
C      nwin = -999    !readin condition. Not readin NW file
C      incwfin= 0     !readin condition. use ForX0 for core in GWIN
C      efin =  -999d0 !readin condition. Not readin EFERMI
C      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
C      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
Cc---  These are allocated and setted by genallcf_v3
Cc      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
Cc     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
Cc      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
Cc      character(120):: symgrp
Cc      character(6),allocatable :: clabl(:)
Cc      integer(4),allocatable:: iclass(:)
Cc     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
Cc     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
Cc     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
Cc     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
Cc     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
Cc     &    occv(:,:,:),unoccv(:,:,:)
Cc     &   ,occc(:,:,:),unoccc(:,:,:),
Cc     o    nocc(:,:,:),nunocc(:,:,:)
Cc      real(8), allocatable::
Cc     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)
C
C!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
C      if(nclass /= natom) stop ' nclass /= natom '
C
Cc --- read dimensions of h,hb
C      ifhbe      = iopen('hbe.d',1,0,0)
C      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
C      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
C      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'
C
Cc --- read by rdpp ; Radial integrals ppbrd and plane wave part
C      call getsrdpp2(nclass,nl,nxx)
C      call readngmx('QGpsi',ngpmx)
C      print *,' ngpmx=',ngpmx
C
Cc --- read radial functions PHIVC   (taken from hasfp0)
C      print *,' Go to readining phivc'
C      ifphi  = iopen('PHIVC', 0,-1,0)     ! PHIV+PHIC augmentation wave and core
C      read(ifphi) nbas, nradmx, ncoremx,nrx
C      if( nbas/=natom ) stop ' nbas(PHIVC) /= natom '
C      deallocate(ncore)
C      allocate(  ncindx(ncoremx,nbas),
C     &           lcindx(ncoremx,nbas),
C     &           nrad(nbas),
C     &           nindx_r(1:nradmx,1:nbas),
C     &           lindx_r(1:nradmx,1:nbas),
C     &        aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) ,
C     &        phitoto(nrx,0:nl-1,nn,nbas,nspin),
C     &        phitotr(nrx,0:nl-1,nn,nbas,nspin),
C     &        nc_max(0:nl-1,nbas),ncore(nbas),rmax(nbas) )
C      print *,' end of allocation'
C      read(ifphi) nrad(1:nbas)
C      read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
C      nc_max=0
C      do ibas=1,nbas
C        ic = ibas
C        write(6,*)' --- read PHIVC of ibas nrad=',ibas,nrad(ic)
C        read(ifphi) ncore(ic), ncoremx                            !core
C        read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
C        print *,' xxx0'
C        read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
C
C        print *, 'ic icx=',ic,icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
C        if(ic/=icx) then
C          stop ' h_uu: ic/=icx'
C        endif
C        print *,' xxx1 ncoremx ncore(ic)=',ncoremx,ncore(ic)
C        read(ifphi) rr(1:nrofi(ic),ic)
C        print *,' xxx2 ncoremx ncore(ic)=',ncoremx,ncore(ic)
C
C        print *,' xxx2 nspin=',nspin
C        rmax(ic) = rr(nrofi(ic),ic)
C        do isp = 1, nspin
C          write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
C          do ico = 1, ncore(ic) !core
C            l =  lcindx(ico,ic)
C            n =  ncindx(ico,ic)
C            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
C            phitotr(1:nrofi(ic),l,n, ic,isp)=              !core raw= core orthgonal
C     &      phitoto(1:nrofi(ic),l,n, ic,isp)               !
C            if(n>nc_max(l,ic)) nc_max(l,ic)=n
C            print *,' sss1c=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
C          enddo
C          do irad = 1, nrad(ic)   !valence
C            l = lindx_r (irad,ic)
C            n = nindx_r (irad,ic) + nc_max(l,ic)
C            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
C            read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
C            print *,' sss1=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
C            print *,' sss2=',sum(abs(phitotr(1:nrofi(ic),l,n, ic,isp)))
C          enddo
C        enddo
C      enddo
C
Cc--- cg coefficient.  y = cg y y ; y is the real spherical harmonics
C      ngrpx=1
C      allocate( cg(nl**2,nl**2,(2*nl-1)**2), symope(3,3) )
C      symope(1:3,1) = (/1d0,0d0,0d0/)
C      symope(1:3,2) = (/0d0,1d0,0d0/)
C      symope(1:3,3) = (/0d0,0d0,1d0/)
C      cg = 0d0 !for sanity check
C      call rotcg(nl-1,symope,ngrpx,cg)
C
Cc --- initiallization to get eigenfunctions
C      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
C      call init_readeigen2(mrecb,nlmto,mrecg)
C      call readngmx('QGpsi',ngpmx)
C      allocate( geig1(ngpmx,nband),geig2(ngpmx,nband))
C      print *, 'end of initialization'
C
Cc --- Readin nlam index
C      ifoc = iopen('@MNLA_CPHI',1,0,0)
C      ldim2 = nlmto
C      read(ifoc,*)
C      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
C      do ix =1,ldim2
C        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
C        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
C      enddo
C
Cc ---  q near zero
C      print *, 'reading QOP'
C      open (101,file='Q0P')
C      read (101,"(i5)") nq0i
C!      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
C      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
C      nq0it = nq0i
C      allocate( q0i(1:3,1:nq0i) ) !wqt(1:nq0i),
C!      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
C      nq0ix = nq0i
C      do i=1,nq0i
C        read (101,* ) xxx,q0i(1:3,i)
C        if(xxx==0d0 ) nq0ix = i-1
C      enddo
C      nq0i = nq0ix ! New nq0i July 2001
C      write(6,*) ' Used k number in Q0P =', nq0i
C      write(6,"(i3,f14.6,2x, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
C      close(101)
C
Cctm
Cc read BBVEC
C      ifbb = iopen('BBVEC',1,0,0)
C      if (ixc <= 3) then
C        read(ifbb,*)
C        read(ifbb,*)nbb,nqbz2
C        if (nqbz /= nqbz2) stop 'readbb: nqbz is wrong!'
C        allocate (bbv(3,nbb),ikbidx(nbb,nqbz))
C        call  readbb(ifbb,nqbz,nspin,nbb,
C     o             bbv,
C     o             ikbidx,
C     o             iko_ixs,iko_fxs,noxs)
C      elseif (ixc == 4) then
C        read(ifbb,*)
C        read(ifbb,*)nbb,nqbz2,nb1d
C        if (nbb /= 1) stop 'readbb: nbb is wrong!'
C        if (nqbz /= nqbz2) stop 'readbb: nqbz is wrong!'
C        allocate (bbv(3,nbb),qbz2(3,nqbz))
C        allocate(ibidx(nb1d,nqbz),ibidxs(nb1d,nqbz),
C     >           ibidx0(nb1d,nqbz,nq0i),ij1idx(nb1d),ij2idx(nb1d))
C        call  readbb1D(bzcase(),ifbb,nqbz,nspin,nbb,nb1d,nq0i,
C     i             qbz,qbzs,q0i,ginv,
C     o             bbv,
C     o             ibidx,ibidxs,ibidx0,
C     o             iko_ixs,iko_fxs,noxs)
C      endif
C
Cc open output files
Cc mode 1
C      if (ixc == 1) then
C
C        ifq0p = iopen('QPwf',1,0,0)
C        read(ifq0p,*)
C        read(ifq0p,"(d24.16,3x, 3d24.16)" ) wgt,q0wf(1:3)
C
C        ifuu(1) = iopen('UU0U',0,-1,0)
C        write(ifuu(1))'nqbz,nbb,iko_ixs(1),iko_fxs(1)'
C        write(ifuu(1))nqbz,nbb,iko_ixs(1),iko_fxs(1)
C        if (nspin == 2) then
C          ifuu(2) = iopen('UU0D',0,-1,0)
C          write(ifuu(2))'nqbz,nbb,iko_ixs(2),iko_fxs(2)'
C          write(ifuu(2))nqbz,nbb,iko_ixs(2),iko_fxs(2)
C        endif
C
Cc mode 2
C      elseif (ixc == 2) then
C
C        ifuu(1) = iopen('UUU',0,-1,0)
C        write(ifuu(1))'nqbz,nbb,iko_ixs(1),iko_fxs(1)'
C        write(ifuu(1))nqbz,nbb,iko_ixs(1),iko_fxs(1)
C        if (nspin == 2) then
C          ifuu(2) = iopen('UUD',0,-1,0)
C          write(ifuu(2))'nqbz,nbb,iko_ixs(2),iko_fxs(2)'
C          write(ifuu(2))nqbz,nbb,iko_ixs(2),iko_fxs(2)
C        endif
C
Cc mode 3
C      elseif (ixc == 3) then
C
C        ifuu(1) = iopen('UUq0U',0,-1,0)
C        write(ifuu(1))'nqbz,nq0i,iko_ixs(1),iko_fxs(1)'
C        write(ifuu(1))nqbz,nq0i,iko_ixs(1),iko_fxs(1)
C        if (nspin == 2) then
C          ifuu(2) = iopen('UUq0D',0,-1,0)
C          write(ifuu(2))'nqbz,nq0i,iko_ixs(2),iko_fxs(2)'
C          write(ifuu(2))nqbz,nq0i,iko_ixs(2),iko_fxs(2)
C        endif
C
Cc mode 4
C      elseif (ixc == 4) then
C        nbbloop = nq0i+1
C        if (bzcase()==2) nbbloop = nq0i+2
C        ifuu(1) = iopen('UU1dU',0,-1,0)
C        write(ifuu(1))nqbz,nq0i,nb1d,nbbloop
C        if (nspin == 2) then
C          ifuu(2) = iopen('UU1dD',0,-1,0)
C          write(ifuu(2))nqbz,nq0i,nb1d,nbbloop
C        endif
C
Cc end of ixc-if
C      endif
C
C
Cc======================================================================
Cc --- Set q1(j1range) q2(j2range)
Cc======================================================================
C! Note that the true q when we generate eigenfunctions are q1x and q2x.
C! q1-q1x should be a G vector.
C! So you may need to take into account the phase shift to <u|u> vectors.
C!
C! --- I inserted checkagree to make sure that q1=q1x and q2=q2x ...
C!
Cctm      q1 = qbz(:, 12)           ;  j1min=1;   j1max=8
Cctm      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
Cc      q1 = qbz(:, 5);  j1min=1;   j1max=8
Cc      q2 = qbz(:, 5);  j2min=1;   j2max=8
Cc      q1 = qbz(:, 12)+q0i(1:3,1);  j1min=1;   j1max=8
Cc      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
Cctm
C      j1min = iko_ixs(1)
C      j1max = iko_fxs(1)
C      if (nspin == 2) then
C        if (iko_ixs(2) < j1min) j1min = iko_ixs(2)
C        if (iko_fxs(2) > j1max) j1max = iko_fxs(2)
C      endif
C      j2min = j1min
C      j2max = j1max
Cc======================================================================
C
C      allocate( uum(j1min:j1max,j2min:j2max,nspin) )
C
Cctm
C      do 1070 iqbz = 1,nqbz
C        write(*,*)'iqbz =',iqbz, 'out of',nqbz
C
C        if (ixc == 1) then
C          nbbloop = 1
C        elseif (ixc == 2) then
C          nbbloop = nbb
C        elseif (ixc == 3) then
C          nbbloop = nq0i
C        elseif (ixc == 4) then
C          nbbloop = nq0i+1
C          if (bzcase()==2) nbbloop = nq0i+2
C        endif
C
C        do 1080 ibb = 1,nbbloop
C
C          if (ixc == 1) then
C            q1(:) = qbz(:,iqbz)
C            q2(:) = q0wf(:)
C
C          elseif (ixc == 2) then
C
C            write(*,*)'ibb  =',ibb, 'out of',nbb
C
C            iqb = ikbidx(ibb,iqbz)
C            q1(:) = qbz(:,iqbz)
Cc            q2 = qbz(:,iqb)
C            q2(:) = q1(:) + bbv(:,ibb)
C
C            if (iqb < iqbz) then
C              iqtmp = iqb
C              do ibb2 = 1,nbb
C                itmp = ikbidx(ibb2,iqtmp)
C                if (itmp == iqbz) then
C                  ibbtmp = ibb2
C                  goto 1200
C                endif
C              enddo
C              stop 'huumat: (iq,ib) error'
C 1200         continue
C
C              do ispin = 1,nspin
C                write(ifuu(ispin))-20
C                write(ifuu(ispin))iqbz,ibb,iqtmp,ibbtmp
C              enddo
C              goto 1090
C            endif
C
C          elseif (ixc == 3) then
C
C            write(*,*)'ibb  =',ibb, 'out of',nbb
C
C            q1(:) = qbz(:,iqbz)
Cc            q2 = qbz(:,iqb)
C            q2(:) = q1(:) + q0i(:,ibb)
C
C          elseif (ixc == 4) then
C
C            write(*,*)'iq0i  =',ibb, 'out of',nbbloop
C
C            if (bzcase()==2) then
C              if (ibb == 1) then
C                q1(:) = qbz(:,iqbz)
C                q2(:) = q1(:) + bbv(:,1)
C                ij1idx(:) = ibidx(:,iqbz)
C                iqtmp = iqindx(q2(:), ginv,qbz,nqbz)
C                ij2idx(:) = ibidx(:,iqtmp)
C              elseif (ibb == 2) then
C                q1(:) = qbzs(:,iqbz)
C                q2(:) = q1(:) + bbv(:,1)
C                ij1idx(:) = ibidxs(:,iqbz)
C                iqtmp = iqindx(q2(:), ginv,qbzs,nqbz)
C                ij2idx(:) = ibidxs(:,iqtmp)
C              else
C                iq0i = ibb-2
C                q1(:) = qbz(:,iqbz)+q0i(:,iq0i)
C                q2(:) = q1(:) + bbv(:,1)
C                ij1idx(:) = ibidx0(:,iqbz,iq0i)
C                do iq = 1,nqbz
C                  qbz2(:,iq) = qbz(:,iq) + q0i(:,iq0i)
C                enddo
C                iqtmp = iqindx(q2(:), ginv,qbz2,nqbz)
C                ij2idx(:) = ibidx0(:,iqtmp,iq0i)
C              endif
C            else
C              if (ibb == 1) then
C                q1(:) = qbz(:,iqbz)
C                q2(:) = q1(:) + q0i(:,ibb)
C                ij1idx(:) = ibidx(:,iqbz)
C                iqtmp = iqindx(q2(:), ginv,qbz,nqbz)
C                ij2idx(:) = ibidx(:,iqtmp)
C              else
C                iq0i = ibb-1
C                q1(:) = qbz(:,iqbz)+q0i(:,iq0i)
C                q2(:) = q1(:) + bbv(:,1)
C                ij1idx(:) = ibidx0(:,iqbz,iq0i)
C                do iq = 1,nqbz
C                  qbz2(:,iq) = qbz(:,iq) + q0i(:,iq0i)
C                enddo
C                iqtmp = iqindx(q2(:), ginv,qbz2,nqbz)
C                ij2idx(:) = ibidx0(:,iqtmp,iq0i)
C              endif
C            endif
Cc end of ixc-if
C          endif
C
Cc debug
Cc         q2 = q1
Cc         write(*,*)'q1 =',q1
Cc         write(*,*)'q2 =',q2
C
CC --- q1x and q2x
C          call readqg0('QGpsi',q1,ginv,  q1x, ngp1)
C          call readqg0('QGpsi',q2,ginv,  q2x, ngp2)
Cc      call checkagree(q1,q1x,' q1 ne q1x')
Cc      call checkagree(q2,q2x,' q2 ne q2x')
C          write(6,"('q1 q1x=',3f9.4,3x,3f9.4)") q1,q1x
C          write(6,"('q2 q2x=',3f9.4,3x,3f9.4)") q2,q2x
C
Cc ... dq
C          dqx = q1x-q2x
C          call q2q0(dqx,plat,dqx0)
C          call q02q0g0(dqx0,dq0,ndg)
C          ndg = -ndg
C          call q02q(dq0,qbas,dq)
C          dg = dq - dqx
C
Cc debug ***************
C          write(96,*)'plat'
C          write(96,*)plat
C          write(96,*)'qbas'
C          write(96,*)qbas
C
C          write(96,*)'iq=',iqbz
C          write(96,*)'dqx'
C          write(96,"(3f12.6)")dqx
C          write(96,*)'dqx0'
C          write(96,"(3f12.6)")dqx0
C          write(96,*)'dq0'
C          write(96,"(3f12.6)")dq0
C          write(96,*)'ndg'
C          write(96,*)ndg
C          write(96,*)'dq'
C          write(96,"(3f12.6)")dq
C          write(96,*)'dg'
C          write(96,"(3f12.6)")dg
Cc debug ***************
C
C
C          if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
C
C          absdq = sqrt(sum(dq**2))
C          absqg2 = (2*pi/alat)**2 *sum(dq**2)
C          absqg =sqrt(absqg2)
C
CC --- ppovl= <P_{q1+G1}|P_{q2+G2}>
C          allocate( ngvecpf1(3,ngp1), ngvecpf2(3,ngp2), ppovl(ngp1,ngp2) )
C          call readqg('QGpsi',q1,ginv, q1x, ngp1, ngvecpf1)
C          call readqg('QGpsi',q2,ginv, q2x, ngp2, ngvecpf2)
Cc      call checkagree(q1,q1x,' q1 ne q1x xxx2')
Cc      call checkagree(q2,q2x,' q2 ne q2x xxx2')
Cc      print *,' ngp1,ngp2 sum check=',ngp1,ngp2,sum(abs(ngvecpf1)),sum(abs(ngvecpf2))
Cc 19June2004 T.Miyake  ccccccccccccccccccccccccccccccc
C          do i = 1,3
C            ngvecpf2(i,1:ngp2) = ngvecpf2(i,1:ngp2) + ndg(i)
C          enddo
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C          call mkppovl2(alat,plat,qbas, !22April2004
C     &    ngp1, ngvecpf1,
C     &    ngp2, ngvecpf2,
C     &    nbas, rmax, pos,
C     o    ppovl)
Cc      print *,' end of mkppovl2'
C
Cc ... lxx and allocations
C          lxx=2*(nl-1)
C          allocate( ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nspin,nbas),
C     &   rprodx(nrx,0:lxx),
C     &   phij(0:lxx),psij(0:lxx),rphiphi(nrx))
C
Cc ... dq
Cc      dq = q1x-q2x
Cc      if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
Cc
Cc      absdq = sqrt(sum(dq**2))
Cc      absqg2 = (2*pi/alat)**2 *sum(dq**2)
Cc      absqg =sqrt(absqg2)
C
Cc ... YL(dq)
C          allocate(cy((lxx+1)**2),yl((lxx+1)**2))
C          call sylmnc(cy,lxx)
C          call sylm(dq/absdq,yl,lxx,r2s) !spherical factor Y(dq)
C
CC --- radial integral  ppbrd = <phi phi j_l>
C          ppbrd=0d0
C          do 900 ibas = 1,nbas
C            ic = ibas
Cc        write(6,"(' nindx=',10i3)") nindx(1:nl,ic)
Cc        write(6,*)' radial integral ibas=',ibas
C            do ir =1,nrofi(ic)
C              call bessl(absqg2*rr(ir,ibas)**2,lxx,phij,psij)
Cc  phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
C              do lx = 0, lxx
C                if(rr(ir,ibas)==0d0) then
C                  rprodx(ir,lx)=0d0
C                else
C                  rprodx(ir,lx) = rr(ir,ibas)* phij(lx)* (absqg*rr(ir,ibas))**lx
C                endif
C                ! = r \times j_l(|dq|r)  !bessel function
C              enddo
Cccccccccccccccccccccccc
Cc          write(1100,"(10d13.6)")rr(ir,ibas),rprodx(ir,0:lxx) !,phij(0:lxx)
Cccccccccccccccccccccccc
C            enddo
C            do 125 isp = 1,nspin
C              do 25 l1 = 0, nl-1
C              do 25 n1 = 1, nindx(l1+1,ic)
C              do 25 l2 = 0, nl-1
C              do 25 n2 = 1, nindx(l2+1,ic)
C                rphiphi(1)       = 0d0
C                rphiphi(2:nrofi(ic)) = phitoto(2:nrofi(ic),l1,n1,ic,isp)
C     &                          *phitoto(2:nrofi(ic),l2,n2,ic,isp)/rr(2:,ic) ! phi = u = r \phi
C              do 25 lx = 0, 2*(nl-1)
C                if(lx <abs(l1-l2) .or. l1+l2<lx) cycle
C                call gintxx( rprodx(1,lx), rphiphi,aa(ic),bb(ic),nrofi(ic),
C     &        ppbrd(l1, n1,l2, n2, lx, isp,ibas) )
Cc          if(l1==l2.and.n1==n2.and.lx==0)
Cc         print *, ' ppbrd=',l1,n1,ppbrd(l1, n1,l2, n2, lx, isp,ibas)
C   25         continue
C 125        continue
C 900        continue
C
CC --- Calcuate <u{q1x j1} | u_{q2x j2}>
Cc              = < exp(i(q1x-q2x)r) psi^*{q1x j1} psi_{q2x j2} >
Cc ... MT part
Ccr   ldim2 = nlmto
Ccr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
Ccr   l_indx   (1:ldim2) : l index
Ccr   ibas_indx(1:ldim2) : ibas index.
C            uum = 0d0
C            do 1050 ispin=1,nspin
C              allocate(cphi1 (nlmto,nband),cphi2(nlmto,nband) )
C              call readcphi(q1, nlmto, ispin, quu, cphi1)
Cc       call checkagree(q1,q1x,' q1 ne quu')
C              call readcphi(q2, nlmto, ispin, quu, cphi2)
Cc       call checkagree(q2,quu,' q2 ne quu')
C
C              do 1020 ia1 = 1,nlmto
C                ibas1= ibas_indx(ia1)
C                l1   = l_indx    (ia1)
C                m1   = m_indx    (ia1)
C                n1   = n_indx    (ia1) + nc_max(l1,ibas1)
C                lm1  = l1**2+l1+1  + m1
C                do 1010 ia2 = 1,nlmto
C                  ibas2 = ibas_indx(ia2)
C                  if(ibas2/=ibas1) cycle
C                  phaseatom = exp( img* 2d0*pi*sum(dq*pos(:,ibas1)) )
C                  m2   = m_indx    (ia2)
C                  l2   = l_indx    (ia2)
C                  n2   = n_indx    (ia2) + nc_max(l2,ibas2)
C                  lm2= l2**2 +l2+1 + m2
Ccccccccccccccccccccccccccccc
Cc Norm check test.
Cc          do j1= j1min,j1max
Cc          do j2= j2min,j2max
Cc            if(ia1==ia2) uum(j1,j2,ispin) = uum(j1,j2,ispin)
Cc     &        + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2)
Cc          enddo
Cc          enddo
Cccccccccccccccccccccccccccccc
C                  do lm3= (l1-l2)**2+1, (l1+l2+1)**2 ! l3 can take |l1-l2|,...l1+l2
C                    l3 = ll(lm3)
Cccccccccccccccccccccccccccccccc
Cc          ylk=0d0;   if(l3==0) ylk=1d0/sqrt(4*pi) !Y_00 only test
Cccccccccccccccccccccccccccccccc
C                    ylk= cy(lm3)*yl(lm3)
C                    ppj = ppbrd(l1,n1,l2,n2,l3,ispin,ibas1) *cg(lm1,lm2, lm3)
C     &          * fpi* img**l3* phaseatom * ylk
C
C! cg(lm1,lm2,lm3)= \int Y_lm3(\hat(r)) Y_lm2(\hat(r)) Y_lm1(\hat(r)) \frac{d \Omega}{4\pi}
C! This is based on inverse expansion. See Rose.Eq.3.8.
C                    if (ixc /= 4) then
Cc            do j1= j1min,j1max
Cc            do j2= j2min,j2max
C                      do j1= iko_ixs(ispin),iko_fxs(ispin)
C                        do j2= iko_ixs(ispin),iko_fxs(ispin)
C                          uum(j1,j2,ispin) =
C     &       uum(j1,j2,ispin) + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2)
C     &                        * ppj
C                        enddo
C                      enddo
C                    else
C                      do j1= iko_ixs(ispin),iko_fxs(ispin)
C                        do j2= iko_ixs(ispin),iko_fxs(ispin)
C                          uum(j1,j2,ispin) =  uum(j1,j2,ispin)
C     &  + dconjg(cphi1(ia1,ij1idx(j1)))*cphi2(ia2,ij2idx(j2)) * ppj
C                        enddo
C                      enddo
C                    endif
C                  enddo
C 1010           continue
C 1020         continue
C
Cc ... Interstitial Plane Wave part
C              call readgeig(q1, ngpmx, ispin, quu, geig1)
Cc       call checkagree(q1,quu,' q1 ne quu eig')
C              call readgeig(q2, ngpmx, ispin, quu, geig2)
Cc       call checkagree(q2,quu,' q1 ne quu eig')
C              if (ixc /= 4) then
Cc         do j1=j1min,j1max
Cc         do j2=j2min,j2max
C                do j1= iko_ixs(ispin),iko_fxs(ispin)
C                  do j2= iko_ixs(ispin),iko_fxs(ispin)
C                    uum(j1,j2,ispin)= uum(j1,j2,ispin) +
C     &      sum( dconjg(geig1(1:ngp1,j1))
C     &    * matmul(ppovl,geig2(1:ngp2,j2)) )
C                  enddo
C                enddo
C              else
C                do j1= iko_ixs(ispin),iko_fxs(ispin)
C                  do j2= iko_ixs(ispin),iko_fxs(ispin)
C                    uum(j1,j2,ispin)= uum(j1,j2,ispin) +
C     &      sum( dconjg(geig1(1:ngp1,ij1idx(j1)))
C     &    * matmul(ppovl,geig2(1:ngp2,ij2idx(j2))) )
C                  enddo
C                enddo
C              endif
C              deallocate(cphi1, cphi2)
C 1050       continue
C
C
Cctm
C            do ispin = 1,nspin
C              iti = iko_ixs(ispin)
C              itf = iko_fxs(ispin)
C              write(ifuu(ispin))-10
C              if (ixc == 1) then
C                write(ifuu(ispin)) iqbz
C              elseif (ixc == 2) then
C                write(ifuu(ispin))
C     &            iqbz,ibb,ikbidx(ibb,iqbz)
C              elseif (ixc == 3.or.ixc == 4) then
C                write(ifuu(ispin)) iqbz,ibb
C              else
C              endif
C              write(ifuu(ispin))
C     &         ((uum(j1,j2,ispin),j1=iti,itf),j2=iti,itf)
C            enddo
C
Ccheck write
C            if (mod(iqbz,10) == 1) then
C              do ispin = 1,nspin
C                iti = iko_ixs(ispin)
C                itf = iko_fxs(ispin)
C                do j1 = iti,itf
C                  uunorm = 0d0
C                  do j2 = iti,itf
C                    uunorm = uunorm + abs(uum(j1,j2,ispin))**2
C                  enddo
C                  write(*,980)ispin,iqbz,ibb,j1,uunorm
C                enddo
C              enddo
C            endif
C
Cc--- write resutlt
Cc      print *,' ============ result --- diagonal --- =============='
Cc      do ispin = 1,nspin
Cc      do j1=j1min,j1max
Cc      do j2=j2min,j2max
Cc       if(j1==j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)")
Cc     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
Cc      enddo
Cc      enddo
Cc      print *,'--- off diagonal ----------------------'
Cc      do j1=j1min,j1max
Cc      do j2=j2min,j2max
Cc       if(j1/=j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)")
Cc     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
Cc      enddo
Cc      enddo
Cc      enddo
C
C            deallocate(ngvecpf1, ngvecpf2, ppovl, ppbrd,
C     &   rprodx, phij, psij, rphiphi,
C     &   cy, yl)
C
C 1090       continue
C
Cc end of bb-loop
C 1080       continue
Cc end of iqbz-loop
C 1070       continue
C
C
C 980        format('uu: norm check',i4,i5,2i4,f12.6)
C 990        format(3i6)
C
C
C            call mpi_end()
C
C            print *, ' ====== end ========================================'
Cc      stop ' ====== end ========================================'
            end

            subroutine checkagree(a,b,char)
            real(8):: a(3),b(3)
            character*(*) :: char
            if(sum(abs(a-b))>1d-6) then
              print *,' Error in checkagree:',char
              stop ' Error in checkagree:'
            endif
            end
c-----------------------------------------------------------------------
            subroutine  readbb(ifbb,nqbz,nspin,nbb,
     o                   bbv,
     o                   ikbidx,
     o                   iko_ixs,iko_fxs,noxs)

            implicit real*8(a-h,o-z)
            parameter (eps = 1d-4)

            real (8) :: u(3),bbv(3,nbb)
            integer (4) :: iopen,
     i           iko_ixs(2),iko_fxs(2),noxs(2)
            integer(4):: ikbidx(nbb,nqbz)
c      integer(4),allocatable:: ikidx(:),ikbidx(:,:)

c      ifbb = iopen('BBVEC',1,0,0)
c      read(ifbb,*)
c      read(ifbb,*)nbb,nqibz2,nqbz2
c      if (nqibz /= nqibz2) stop 'readbb: nqibz is wrong!'
c      if (nqbz /= nqbz2) stop 'readbb: nqbz is wrong!'
c
c      allocate (ikidx(nqibz),ikbidx(nbb,nqibz))

            do i = 1,nbb
              read(ifbb,*)bbv(1,i),bbv(2,i),bbv(3,i),dummy4
            enddo
            do iq = 1,nqbz
              read(ifbb,*)itmp,u(1:3)
              do ib = 1,nbb
                read(ifbb,*)itmp,itmp2,ikbidx(ib,iq),u(1:3)
              enddo
            enddo
            read(ifbb,*)
            read(ifbb,*)nspin2
            if (nspin /= nspin2) stop 'nspin is wrong!'
            do is = 1,nspin
              read(ifbb,*)iko_ixs(is),iko_fxs(is),noxs(is)
            enddo


            return
            end
c-----------------------------------------------------------------------
            subroutine  readbb1D(ibzcase,ifbb,nqbz,nspin,nbb,nband,nq0i,
     i                   qbz,qbzs,q0i,ginv,
     o                   bbv,
     o                   ibidx,ibidxs,ibidx0,
     o                   iko_ixs,iko_fxs,noxs)

            implicit real*8(a-h,o-z)
            parameter (eps = 1d-4)

            real (8) :: u(3),bbv(3,nbb),ginv(3,3),qbz(3,nqbz),
     >            qbzs(3,nqbz),q0i(3,nq0i),q(3),qbz2(3,nqbz)
            integer (4) :: iopen,iqindx,itmp,i,
     i           iko_ixs(2),iko_fxs(2),noxs(2)
            integer(4):: ibidx(nband,nqbz),ibidxs(nband,nqbz)
     >            ,ibidx0(nband,nqbz,nq0i)
     >            ,ibidx2(nband)
c      integer(4),allocatable:: ikidx(:),ikbidx(:,:)

c read bb-vector
            do i = 1,nbb
              read(ifbb,*)bbv(1,i),bbv(2,i),bbv(3,i)
            enddo
c read q-vector

c unshifted grid
            do iq = 1,nqbz
              read(ifbb,*)q(1:3)
              read(ifbb,*)ibidx2(1:nband)
              itmp = iqindx(q(:), ginv,qbz,nqbz)
              ibidx(:,itmp) = ibidx2(:)
            enddo

c shifted grid
            if (ibzcase == 2) then
              do iq = 1,nqbz
                read(ifbb,*)q(1:3)
                read(ifbb,*)ibidx2(1:nband)
                itmp = iqindx(q(:), ginv,qbzs,nqbz)
                ibidxs(:,itmp) = ibidx2(:)
              enddo
            endif

c +q0 grid
            do iq0=1,nq0i
              do iq=1,nqbz
                qbz2(:,iq) = qbz(:,iq)+q0i(:,iq0)
              enddo
              do iq = 1,nqbz
                read(ifbb,*)q(1:3)
                read(ifbb,*)ibidx2(1:nband)
                itmp = iqindx(q(:), ginv,qbz2,nqbz)
                ibidx0(:,itmp,iq0) = ibidx2(:)
              enddo
            enddo

c iko_ixs and iko_fxs
            do is = 1,nspin
              iko_ixs(is) = 1
              iko_fxs(is) = nband
              noxs(is) = nband
            enddo

            return
            end
c-----------------------------------------------------------------------
            subroutine q2q0(q,plat,q0)

            implicit real*8(a-h,o-z)
            parameter (eps=1d-4)
            real(8) :: q(3),q0(3),plat(3,3)

            do ii = 1,3
              q0(ii) = sum(q(:)*plat(:,ii))
            enddo

            return
            end
c-----------------------------------------------------------------------
            subroutine q02q(q0,qbas,q)

            implicit real*8(a-h,o-z)
            real(8) :: q(3),q0(3),qbas(3,3)

            q(:) = 0d0
            do ii = 1,3
              q(:) = q(:) + qbas(:,ii)*q0(ii)
            enddo

            return
            end
c-----------------------------------------------------------------------
            subroutine q02q0g0(qin,qout,ng)
            implicit real*8(a-h,o-z)
c
c qin = qout + ng

            parameter (eps=1d-6)
            real(8) :: qin(3), qout(3)
            integer(4) :: ng(3)

            do i = 1,3
              q = qin(i) - eps
              ng(i) = nint(q)
c         g = anint(q)
              qout(i) = q - dble(ng(i)) + eps
            enddo

            return
            end
c-----------------------------------------------------------------------
