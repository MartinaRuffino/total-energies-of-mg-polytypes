      program h_uumatrix
c-------------------------------------------------
c  Calculate <u|u> matrix . u_kj(r) is the perodic part of eigencuntion.
c This is a test routine. A bit detailed comment---so this code is a kind of manual
c  to see the treatment of eigenfunctions.
c
c-------------------------------------------------
      use m_readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,readcphi,readgeig
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue

C      implicit none
C      real(8):: q(3),  qgbin(3),qx(3)
C      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
C     & , ngvecpf1(:,:), ngvecpf2(:,:),
C     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
C      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:),
C     &phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:)
C      complex(8),parameter:: img=(0d0,1d0)
Cc,nw,incwf,natom,nclass,ipos,igrp,
Cc     & iinvg,nspin,nl,nn,nnv,nnc,
Cc     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions
Cc     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
Cc     o                   izdummy,
Cc     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
Cc     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
Cc     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
Cc     o   iecore,ikonf,iicore,incore,nctot,             !core
Cc     o   imagw_dummy,niw,idummy,
C      integer(4)
C     &   nw_input,
C     &   ifhbe,jobBes,
C     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
C     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
C     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
C     &   ifepscond,nxx,ifvxcpout,ifgb0vec
C     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
C     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
C     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
C     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax  !,ifev1,ifev2
C     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
C     &   ,ngpmx
C
C      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini,
C     & dwry,dwh,omg_c,omg2,xxx
C      integer(4)::nwin, incwfin,  verbose
C      real(8)::efin
C      integer(4):: bzcase, mrecg,ifphi,
C     & nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc,
C     & ldim2,ixx,ngp1,ngp2,nq0it
C      real(8):: qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
C      real(8),parameter::  pi =     3.14159265358979323846d0
C      real(8),parameter::  fpi =    4d0*pi
C
C      logical:: qbzreg
C!-------------------------------------------------------------------------
C      integer(4),allocatable:: ncindx(:,:),
C     &           lcindx(:,:),
C     &           nrad(:),
C     &           nindx_r(:,:),
C     &           lindx_r(:,:),
C     &           nc_max(:,:),
C     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
C      real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:)
C     &                     ,phitotr(:,:,:,:,:),
C     &        bb(:),zz(:),rmax(:),cy(:),yl(:)
C
C
C      complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:)
C     & ,uum(:,:,:), ppovl(:,:)
C      complex(8):: ppj,phaseatom
C      real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
C      integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin
C     & ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx
C     &, iopen,ico,lxd,lx
C      real(8):: ylk
C      procedure(integer) :: GWversion
C#ifdef COMMONLL
C      integer(4):: ll(51**2)
C      common/llblock/ll
C#else
C      integer(4) ll
C#endif
Cc-------------------------
C      call headver('h_uumatrix',0)
C      jobBes = mod(GWversion()/1000,10)
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
C      if(bzcase()==2.and.qbzreg()) then
C        deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
C        do i=1,nqbz
C          qbz(:,i) = qbz(:,i) -deltaq/2d0
C          write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
C        enddo
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
CCstop2rx 2013.08.09 kino      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
C      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
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
CCstop2rx 2013.08.09 kino      if(nclass /= natom) stop ' nclass /= natom '
C      if(nclass /= natom) call rx( ' nclass /= natom ')
C
Cc --- read dimensions of h,hb
C      ifhbe      = iopen('hbe.d',1,0,0)
C      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
CCstop2rx 2013.08.09 kino      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
C      if(nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
CCstop2rx 2013.08.09 kino      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'
C      if(nqbz /=nqbzt ) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')
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
CCstop2rx 2013.08.09 kino      if( nbas/=natom ) stop ' nbas(PHIVC) /= natom '
C      if( nbas/=natom ) call rx( ' nbas(PHIVC) /= natom ')
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
CCstop2rx 2013.08.09 kino          stop ' h_uu: ic/=icx'
C          call rx( ' h_uu: ic/=icx')
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
CCstop2rx 2013.08.09 kino        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
C        if(ixx/=ix) call rx( 'failed to readin @MNLA_CPHI')
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
C      write(6,"(i3, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
C      close(101)
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
C      q1 = qbz(:, 1) ;  j1min=1;   j1max=8
C      q2 = qbz(:, 2) ;  j2min=1;   j2max=8
Ccc      q1 = qbz(:, 12)           ;  j1min=1;   j1max=8
Ccc      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
Cc      q1 = qbz(:, 5);  j1min=1;   j1max=8
Cc      q2 = qbz(:, 5);  j2min=1;   j2max=8
Cc      q1 = qbz(:, 12)+q0i(1:3,1);  j1min=1;   j1max=8
Cc      q2 = qbz(:, 12)+q0i(1:3,1);  j2min=1;   j2max=8
Cc======================================================================
C      allocate( uum(j1min:j1max,j2min:j2max,nspin) )
C
CC --- ppovl= <P_{q1+G1}|P_{q2+G2}>
C      call readqg0('QGpsi',q1,ginv,  q1x, ngp1)
C      call readqg0('QGpsi',q2,ginv,  q2x, ngp2)
C      call checkagree(q1,q1x,' q1 ne q1x') ! make sure q1=q1x
C      call checkagree(q2,q2x,' q2 ne q2x')
C      write(6,"('q1 =',3f9.4,3x,3f9.4)") q1
C      write(6,"('q2 =',3f9.4,3x,3f9.4)") q2
C      allocate( ngvecpf1(3,ngp1), ngvecpf2(3,ngp2), ppovl(ngp1,ngp2) )
C      call readqg('QGpsi',q1,ginv, q1x, ngp1, ngvecpf1)
C      call readqg('QGpsi',q2,ginv, q2x, ngp2, ngvecpf2)
C      call checkagree(q1,q1x,' q1 ne q1x xxx2') !make sure q1 =q1x
C      call checkagree(q2,q2x,' q2 ne q2x xxx2')
C      print *,' ngp1,ngp2 sum check=',ngp1,ngp2,sum(abs(ngvecpf1)),sum(abs(ngvecpf2))
C      call mkppovl2(alat,plat,qbas, !22April2004
C     &    ngp1, ngvecpf1,
C     &    ngp2, ngvecpf2,
C     &    nbas, rmax, pos,
C     o    ppovl)
C      print *,' end of mkppovl2'
C
Cc ... lxx and allocations
C      lxx=2*(nl-1)
C      allocate( ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nspin,nbas),
C     &   rprodx(nrx,0:lxx),
C     &   phij(0:lxx),psij(0:lxx),rphiphi(nrx))
C
Cc ... dq
C      dq = q1x-q2x
C      if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
C
C      absdq = sqrt(sum(dq**2))
C      absqg2 = (2*pi/alat)**2 *sum(dq**2)
C      absqg =sqrt(absqg2)
C
Cc ... YL(dq)
C      allocate(cy((lxx+1)**2),yl((lxx+1)**2))
C      call sylmnc(cy,lxx)
C      call sylm(dq/absdq,yl,lxx,r2s) !spherical factor Y(dq)
C
CC --- radial integral  ppbrd = <phi phi j_l>
C      ppbrd=0d0
C      do 900 ibas = 1,nbas
C        ic = ibas
C        write(6,"(' nindx=',10i3)") nindx(1:nl,ic)
C        write(6,*)' radial integral ibas=',ibas
C        do ir =1,nrofi(ic)
CC         call bessl(absqg2*rr(ir,ibas)**2,lxx,phij,psij)
C          call besslr(absqg2*rr(ir,ibas)**2,10*jobBes+0,0,lxx,phij,psij)
Cc  phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
C          do lx = 0, lxx
C            rprodx(ir,lx) = rr(ir,ibas)* phij(lx)* (absqg*rr(ir,ibas))**lx
C            ! = r \times j_l(|dq|r)  !bessel function
C          enddo
Cccccccccccccccccccccccc
Cc          write(1100,"(10d13.6)")rr(ir,ibas),rprodx(ir,0:lxx) !,phij(0:lxx)
Cccccccccccccccccccccccc
C        enddo
C        do 125 isp = 1,nspin
C          do 25 l1 = 0, nl-1
C          do 25 n1 = 1, nindx(l1+1,ic)
C          do 25 l2 = 0, nl-1
C          do 25 n2 = 1, nindx(l2+1,ic)
C            rphiphi(1)       = 0d0
C            rphiphi(2:nrofi(ic)) = phitoto(2:nrofi(ic),l1,n1,ic,isp)
C     &                          *phitoto(2:nrofi(ic),l2,n2,ic,isp)/rr(2:,ic) ! phi = u = r \phi
C          do 25 lx = 0, 2*(nl-1)
C            if(lx <abs(l1-l2) .or. l1+l2<lx) cycle
C            call gintxx( rprodx(1,lx), rphiphi,aa(ic),bb(ic),nrofi(ic),
C     &        ppbrd(l1, n1,l2, n2, lx, isp,ibas) )
Cc          if(l1==l2.and.n1==n2.and.lx==0)
Cc         print *, ' ppbrd=',l1,n1,ppbrd(l1, n1,l2, n2, lx, isp,ibas)
C   25     continue
C 125    continue
C 900    continue
C
CC --- Calcuate <u{q1x j1} | u_{q2x j2}>
Cc              = < exp(i(q1x-q2x)r) psi^*{q1x j1} psi_{q2x j2} >
Cc ... MT part
Ccr   ldim2 = nlmto
Ccr   n_indx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
Ccr   l_indx   (1:ldim2) : l index
Ccr   ibas_indx(1:ldim2) : ibas index.
C        uum = 0d0
C        do 1050 ispin=1,nspin
C          allocate(cphi1 (nlmto,nband),cphi2(nlmto,nband) )
C          call readcphi(q1, nlmto, ispin, quu, cphi1) !quu is true q, q1-q can be G vectors. But need to check it again!
C                                                      !Need to examine readcphi again.
C          call checkagree(q1,q1x,' q1 ne quu') !check for safe quu is true q
C          call readcphi(q2, nlmto, ispin, quu, cphi2)
C          call checkagree(q2,quu,' q2 ne quu')
C
C          do 1020 ia1 = 1,nlmto
C            ibas1= ibas_indx(ia1)
C            l1   = l_indx    (ia1)
C            m1   = m_indx    (ia1)
C            n1   = n_indx    (ia1) + nc_max(l1,ibas1)
C            lm1  = l1**2+l1+1  + m1
C            do 1010 ia2 = 1,nlmto
C              ibas2 = ibas_indx(ia2)
C              if(ibas2/=ibas1) cycle
C              phaseatom = exp( img* 2d0*pi*sum(dq*pos(:,ibas1)) )
C              m2   = m_indx    (ia2)
C              l2   = l_indx    (ia2)
C              n2   = n_indx    (ia2) + nc_max(l2,ibas2)
C              lm2= l2**2 +l2+1 + m2
Ccccccccccccccccccccccccccccc
Cc Norm check test.
Cc          do j1= j1min,j1max
Cc          do j2= j2min,j2max
Cc            if(ia1==ia2) uum(j1,j2,ispin) = uum(j1,j2,ispin)
Cc     &        + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2)
Cc          enddo
Cc          enddo
Cccccccccccccccccccccccccccccc
C              do lm3= (l1-l2)**2+1, (l1+l2+1)**2 ! l3 can take |l1-l2|,...l1+l2
C                l3 = ll(lm3)
Cccccccccccccccccccccccccccccccc
Cc          ylk=0d0;   if(l3==0) ylk=1d0/sqrt(4*pi) !Y_00 only test
Cccccccccccccccccccccccccccccccc
C                ylk= cy(lm3)*yl(lm3)
C                ppj = ppbrd(l1,n1,l2,n2,l3,ispin,ibas1) *cg(lm1,lm2, lm3)
C     &          * fpi* img**l3* phaseatom * ylk
C
C! cg(lm1,lm2,lm3)= \int Y_lm3(\hat(r)) Y_lm2(\hat(r)) Y_lm1(\hat(r)) \frac{d \Omega}{4\pi}
C! This is based on inverse expansion. See Rose.Eq.3.8.
C                do j1= j1min,j1max
C                  do j2= j2min,j2max
C                    uum(j1,j2,ispin) =
C     &      uum(j1,j2,ispin) + dconjg(cphi1(ia1,j1))*cphi2(ia2,j2) * ppj
C                  enddo
C                enddo
C              enddo
C 1010       continue
C 1020     continue
C
Cc ... Interstitial Plane Wave part
C          call readgeig(q1, ngpmx, ispin, quu, geig1)
C          call checkagree(q1,quu,' q1 ne quu eig')
C          call readgeig(q2, ngpmx, ispin, quu, geig2)
C          call checkagree(q2,quu,' q1 ne quu eig')
C          do j1=j1min,j1max
C            do j2=j2min,j2max
C              uum(j1,j2,ispin)= uum(j1,j2,ispin) +
C     &    sum( dconjg(geig1(1:ngp1,j1)) * matmul(ppovl,geig2(1:ngp2,j2)) )
C            enddo
C          enddo
C        deallocate(cphi1,cphi2)
C 1050   continue
C
Cc--- write resutlt
C        print *,' ============ result --- diagonal --- =============='
C        do ispin = 1,nspin
C          do j1=j1min,j1max
C            do j2=j2min,j2max
C              if(j1==j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)")
C     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
C            enddo
C          enddo
C          print *,'--- off diagonal ----------------------'
C          do j1=j1min,j1max
C            do j2=j2min,j2max
C              if(j1/=j2) write(6,"(' ispin=',i2,' j1j2=',2i4,' <u|u>=',2d13.5,' abs=',f13.5)")
C     & ispin,j1,j2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
C            enddo
C          enddo
C        enddo
C        print *, ' ====== end ========================================'
Cc      stop ' ====== end ========================================'
        end

        subroutine checkagree(a,b,char)
        real(8):: a(3),b(3)
        character*(*) :: char
        if(sum(abs(a-b))>1d-6) then
          print *,' Error in checkagree:',char
Cstop2rx 2013.08.09 kino          stop ' Error in checkagree:'
          call rx( ' Error in checkagree:')
        endif
        end
