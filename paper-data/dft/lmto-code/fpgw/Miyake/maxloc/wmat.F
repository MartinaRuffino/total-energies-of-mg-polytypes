      subroutine wmatq(kount,irot,ef,ef2,esmr,esmr2,tr,
ctm, 070501
c     i          iatomp,rsite,nsp,isp, !ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     i          iatomp,
c     i          rws,irws,nrws,
     i          rws1,rws2,nrws1,nrws2,nrws,
     i          nsp,isp, !ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     i                  ifrcw,ifrcwi,
     i                  qbas,ginv,
     i         qibz,qbz,wk,nstbz,wik,nstar,irk,  !koun,,iindxk
c
     i                  iclass,mdim,nlnmv,nlnmc,
     i                  icore,ncore,imdim,
     i                  ppb, !pdb,dpb,ddb,
     i                  freq_r,freqx,wx,expa,ua,dw, !deltaw,freq
     i                  ecore,
c
     d                  nlmto,nqibz,nqbz,nctot,
c     i          index_qbz, n_index_qbz,  !jan2004
     d                  nl,nnc,nclass,natom,
     d                  nlnmx,mdimx,nbloch,ngrp,nw_i,nw,nrw,niw,niwx,nq,
c
c     &     nblochpmx ,ngpn,ngcni,ngpmx,ngcmx,geigB,ngvecpB,ngveccBr,
     &     nblochpmx ,ngpmx,ngcmx, !ngveccBr,!Jan2004
     &     wgt0,wqt,nq0i,q0i,symope,alat, shtv,nband, ifvcfpout,
c     &     shtw,
     &     exchange, !tote,screen,cohtest, ifexsp,
ctetra
ctetra     &     wtet,wtetef,
ctetra    &     ntqx,ibzx,tetraex,
c
c     i omega,iwini,iwend,
     i     nbmx,ebmx, !takao 18June2003
     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !oct2005 for pomat
     i     nwf,
     o          rw_w,cw_w,rw_iw,cw_iw)
c
c 2006 May Takashi Miyake, updated for new fpgw
c 2004 Sep Takashi Miyake, off-site W
c 2004 Jul Takashi Miyake,
c 2004 Apr Takashi Miyake, from sxcf_fal2.f
c
c 2001 Sep. esec=omega(itp,iw). Genral iw mode for exchange =F
c
c 2000 takao kotani. This sxcf is modified from sec.f F.Aryasetiawan.
c
c  exchange=T : Calculate the exchange self-energy
c  exchange=F : Calculate correlated part of the self-energy
c
c---- correlation case documents by ferdi.Aryasetiawan.  -----------------
c 92.02.24
c 93.10.18 from sec.f modified to take into account equivalent atoms
c
c the screened coulomb potential
c Wc(r,r';w)  = W(r,r';w) - v(|r-r'|)
c             = < [r1,r2] v(|r-r1|) X(r1,r2;w) v(|r2-r'|) >
c W(r,r';w)   = < [r''] ei(r,r'';w) v(|r''-r'| >
c ei          = e^(-1), inverse dielectric matrix
c             = 1 + vX
c e           = 1 - vX0 in RPA

c expand Wc(r,r';w) in optimal product basis B
c Wc(r,r';w)  = S[k=FBZ] S[i,j=1,nbloch]
c               B(k,i,r) Wc(k,w)(i,j) B(k,j,r')^*
c Wc(k,w)(i,j) are  the matrix elements of Wc in B
c
c q       = q-vector in SEc(q,t)
c itq     = states t at q
c ntq     = no. states t
c eq      = eigenvalues at q
c ef      = fermi level in Rydberg
c tr      = translational vectors in rot*R = R' + T
c iatomp(R) = R'
c ifrw,ifcw,ifrwi,ifcwi
c   = direct access unit files for Re and Im coulomb matrix
c     along real and imaginary axis
c ifrb,ifcb,ifrhb,ifchb
c         = direct access unit files for Re and Im b,hb
c qbas    = base reciprocal lattice vectors
c ginv    = inverse of qbas s. indxrk.f
cxxxxx ippb,ipdb,idpb,iddb = pointers to work array w for
c  ppb     = <phi(RLn) phi(RL'n') B(R,i)>
c  pdb     = <phi(RLn) phidot(RL'n') B(R,i)>
c  dpb     = <phidot(RLn) phi(RL'n') B(R,i)>
c  ddb     = <phidot(RLn) phidot(RL'n') B(R,i)>
c freq    = frequencies along real axis
c freqx   = gaussian frequencies x between (0,1)
c freqw   = (1-freqx)/freqx
c wx      = weights at gaussian points x between (0,1)
c ua      = constant in exp(-ua^2 w'^2) s. wint.f
c expa    = exp(-ua^2 w'^2) s. wint.f
c dw      = frequency mesh along real axis
c deltaw  = energy mesh in SEc(qt,w) ---Not used now
c iclass  = given an atom, tells the class
c wk      = weight for each k-point in the FBZ
c indexk  = k-point index
c qbz     = k-points in the 1st BZ
c nstar   = no. stars for each k
c irk(k,R) = gives index in the FBZ with k{IBZ, R=rotation
c mdim    = dimension of B(R,i) for each atom R
c work arrays:
c rbq,cbq     = real and imaginary part of b(q)
c rhbq,chbq   = real and imaginary part of hb(q)
c rbkq,cbkq   = real and imaginary part of b(q-k)
c rhbkq,chbkq = real and imaginary part of hb(q-k)
c   b is the eigenvector of the LMTO-Hamiltonian
c ekq     = eigenvalues at q-k
c rmel,cmel = real and imaginary part of
c             <psi(q,t') | psi(q-k,t) B(k,R,i)>
c wr1 ... = work arrays
c dimensions:
c nqibz   = number of k-points in the irreducible BZ
c n1,n2,n3= divisions along base reciprocal lattice vectors
c natom   = number of atoms
c nctot   = no. allowed core states
c nbloch  = total number of Bloch basis functions
c nlnmx   = maximum number of l,n,m
c nlmto   = total number of LMTO basis functions
c ngrp    = no. group elements (rotation matrices)
c niw     = no. frequencies along the imaginary axis
c nw      = no. frequencies along the real axis
c niwx    = max(niw,nw)
c
c----------------------------------------------------------------------
c      implicit real*8(a-h,o-z)
      use m_readqg
      use keyvalue
      use m_readeigen,only: readeval,readcphi,readgeig,
     & init_readeigen_phi_noeval,init_readeigen_mlw_noeval1d,init_readeigen_mlw_noeval
      implicit none
      integer(4) :: ntq, natom,nqbz,nqibz,ngrp,nq,nw_i,nw,niw,
     &  nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,nstatex,
     &   irot,  iqisp,ikpisp,isp,nsp,  nlnmx, iq, ! idxk,
c     &  iwr1,iwr2,iwr3,iwr4,iwc1,iwc2,iwc3,iwc4,
     &  ip, it,itp, it2, itp2, !ifcphi,    ! ifrb,ifcb,ifrhb,ifchb,
     i  iiclass,mdim(*),
     i  ifrcw,ifrcwi, !iindxk,
     &  ifvcfpout,ndummy1,ndummy2,kx,kr,kr2,kr3,ngc,ngb,nbloch, !n1,n2,n3, k,
     &  kp,nt0,nocc, nt0p,nt0m,irkp,i,nt0org,nmax,nt,ntp0,
     &  nbmax,nclass,nl,nnc, nblochpmx,ix,nx,iw,iwp,ixs,ixsmx, !iwp sf 21May02
     &  mdimx, nwx,niwx,
     &  iatomp(natom),
     &  nstar(nqibz),irk(nqibz,ngrp),kount(nqibz),nwf,nrw,nrws !,iclose
c
      real(8) :: q(3),qbas(3*3),ginv(3*3),tr(3,natom),
     &  wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz),
     &  freqx(niw),wx(niw),expa(niw),
     &  eq(nband),
     &   ekq(nband), ekc(nctot+nband),
     &  tpi,ef,ef2,esmr,esmr2,efp,efm,wtx,wfac,wfacx,we,esmrx,ua,
     &  dw,wtt,wexx,www,exx,exxq,weight
c      complex(8) :: zsec(-1:1,ntq,nq)
c      real(8)    ::  shtw
c                       ! This shft is  to avoid some artificial resonance effects.
c                       ! shtw can be zero for esmr/=0 given by takao.
c
      integer(4):: ngpmx, ngcmx,  !ngcni(nqibz), !ngpn(nqbz),
     &   igc, !ngvecpB(3,ngpmx,nqbz),ngveccBr(3,ngcmx,nqibz),
     &  nadd(3)
      real(8) :: wgt0(nq0i,ngrp),wqt(nq0i),qk(3), !qfbz(3),
     &  qbasinv(3,3), det,qdiff(3),add(3),symope(3,3),
     &  qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot),
     &   ppb(1) !pdb(1),dpb(1),ddb(1)
      complex(8)::   alagr3z,wintz !geigB  (ngpmx,nband,nqbz),

c
      real(8),allocatable:: !rmel(:,:,:),cmel(:,:,:),
     &                      rmelt(:,:,:),cmelt(:,:,:),
     &                      rmelt2(:,:,:),cmelt2(:,:,:),
c     &                      rmelt3(:,:,:),cmelt3(:,:,:)
     &                      rmelt3(:,:,:,:),cmelt3(:,:,:,:)
      complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:),
     &     zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:),
     &     zmelt(:,:,:),zmelc(:,:,:,:)
c for exchange --------------------
      logical :: exchange,screen,cohtest,tote
      real(8),allocatable::
     &   w1p(:,:,:),w2p(:,:,:),w3p(:,:)
      complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)

c- debug write ---------------------
      logical :: debug=.false.

cccccccccccccc
c tetra
c      integer(4) :: ntqx
c      integer(4) :: ibzx(nqbz)
c      real(8)    :: wtet  (nband,nqibz,1:ntqx), wtetef(nband,nqibz)
c                ! where the last index is 3*itq+iw-1,itq=1,ntq,iw=-1,1
c     logical    :: tetraex
cccccczzcccccc

      complex(8) :: wintzav,wintzsg

      integer(4) :: ibl,iii,ivsumxxx,ifexsp ,iopen
      integer(4),save::ifzwz=-999

      integer(4) :: iwini, iwend, ia
c      real(8)    :: esec, omega(ntq, iwini:iwend)
      real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws,0:nrw),
     &           cw_w(nwf,nwf,nwf,nwf,nrws,0:nrw),
     &           rw_iw(nwf,nwf,nwf,nwf,nrws,niw),
     &           cw_iw(nwf,nwf,nwf,nwf,nrws,niw)
      complex(8),allocatable:: expikt(:)
      complex(8):: img=(0d0,1d0)
ctakao
      complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)
     .           , cphiqtmp(nlmto,nband)

cccccccccccccccccccccccccccccccccccccccccccccccc faleev 2002
      integer(4) :: nt_max, igb1,igb2,iigb, nw_w
      complex(8),allocatable:: zmel1(:)
      complex(8), allocatable :: zw_(:,:) !,zzmel(:,:)
      complex(8), allocatable :: zwz2(:,:),zw2(:,:,:,:) !0 variant
      complex(8) ::  zz2 ,zwz3(3)
      real(8) :: dd,omg_c,dw2
      real(8) :: freq_r(nw_i:nw)
      complex(8), allocatable :: zw3(:,:,:)


      real(8)::weavx,wfaccut=1d-10

      logical :: GaussSmear
      real(8) :: ebmx,ddw
      integer(4):: nbmx,nbmxe,nstatetot

c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

      integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

      integer(4)::verbose,nstbz(nqbz),iqini,iqend
      real(8):: wgtq0p

      integer(4):: iqindx,nrec,kxx
      real(8)::quu(3),qibz_k(3),qbz_kr(3)
      logical :: onlyQ0P, onlyimagaxis ,noq0p !,noq0p,test_omitq0p,

      logical ::zwz3mode
c      logical ::testimx=.false.

      real(8):: ua_,expa_(niw),ua2,freqw,freqw1,ratio,ua2_(niw)
      logical :: ua_auto
      integer(4):: icc=0
      real(8),allocatable:: uaa(:,:)

c      logical ::testimx=.false.
ccccc zvz test cccccccccccccccccccccccccc
      integer(4):: ngbx
c      complex(8):: vcoul(ngbx,ngbx)
      complex(8),allocatable:: vzz(:,:,:),aaa(:)
      complex(8):: zvz,zvz1
      integer(4):: ib1,ib2,ifix
ccccccccccccccccccccccccccccccccccc
      integer(4) ::nbcut,nbcutc
      logical ::iww2=.true., oncew


C...
      logical::smbasis
      integer(4):: nn,no,ifpomat,iclose,isx,iqx
      complex(8),allocatable:: pomat(:,:)
      real(8):: q_r(3)
      integer(4):: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo)
      complex(8):: pomatr(nnmx,nomx,nkpo)
      real(8):: qrr(3,nkpo)

      real(8):: elxx,ehxx,ekxx,efxx
      integer(4):: ixsmin,iwm,iir,nwxi
      real(8)   :: fffr(3)
      complex(8):: zwzz(3)

ctm
      integer(4) :: nqbz2,nwf2,iko_ix,iko_fx,iqtmp,ifmlw,nko,iqk
     .             ,ifi,in1,in2,imp,ilp,ii,jj,nrws1,nrws2
     .             ,ir1,ir2,ir3,ir
      real(8) :: norm2,qtmp(3),rws1(3,nrws1),rws2(3,nrws2),tmp
      complex(8) :: ztmp,expiqR1(nrws1),expiqR2
      complex(8),allocatable :: cnk(:,:,:),zmel2(:,:,:),zmel3(:,:,:)
      integer(4) :: itq(nwf)
      complex(8) :: weightc(nrws1)

c--------------------------------------------------------------------
      debug=.false.
      if(verbose()>=90) debug=.true.

c oct2005
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
      nbcutc=nctot+nbcut

      tpi         = 8d0*datan(1.d0)
c      iq         = idxk (q,qbz,nqbz) ! index for q
c      write(6,"(' iq q  =',i4,3f8.4)")iq,q
cccc      iq          = idxk (q,qbze,nqbze) ! index for q
      ekc(1:nctot)= ecore(1:nctot)  ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
C     call dinv33(qbas,0,qbasinv,det)
      call minv33(qbas,qbasinv)

c work arrays for psi2br.f
      if(debug) print *, ' sxcf: 1'
      allocate(expikt(natom))

C     if(bzcase()==1) then
        if(abs(sum(qibz(:,1)**2))/=0d0) stop ' sxcf assumes 1st qibz/=0 '
        if(abs(sum( qbz(:,1)**2))/=0d0) stop ' sxcf assumes 1st qbz /=0 '
C      endif

      do it = 1,nwf
        itq(it) = it
      enddo

ctm debug
c        write(*,*)'isp,',isp,rw_w(1,1,1,1,1,0),rw_w(2,2,2,2,1,0)

c-----
      if(exchange) then
        rewind  ifvcfpout
        read(ifvcfpout) ndummy1, ndummy2
      endif

c===============================
c loop over irreducible k-points
c===============================
ccccccccccccccccccccccccccccccccc
c      iii = ivsumxxx(irk,nqibz*ngrp)
c      print *,' sxcf:sum non-zero irk=',iii
c      stop "sss"
c
ccccccccccccccccccccccccccccccccc

C      if(bzcase()==1) then
        kx = 1  ! qibz(:,1)=0 contribution for kcount
        if(irk(kx,irot)/=0) kount(kx)= kount(kx) + 1
c      kount(kx)= kount(kx) + 1
C      endif

        stop 'this routine is out of date'

Cc --- main loop start
C      iqini=2
CC     if(bzcase()==2) iqini=1
C      iqend=nqibz+nq0i
C
Ccccccccccccccccccccccccccccccc
C      call getkeyvalue("GWinput","TestOnlyQ0P",onlyq0p,default=.false.)
C      call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
C      if(noq0p) iqend=nqibz
Cc      iqend=nqibz
Cc      if(test_omitq0p()) then
Cc        iqend=nqibz
Cc        print *,'iqend=',iqend
Cc      endif
Ccccccccccccccccccccccccccccc
C      do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
Cc debug:
Cc      do 1100 kx = iqini,iqini !kx=1 corresponds to q=0 is omitted.
Cc        write (6,"(i3,'  out of ',i3,$)") kx,iqend
C        if(debug)  print *, ' sxcf: goto loop kx=',kx
C
C        ! write(*,'("1  begin k-cycle",$)')
Cc         call cputid(0)
Cc          write(*,*)'kx, ip, irot=',kx, ip,irot
C
C        if( kx <= nqibz ) then
Cc          k  = kx
C          kr = irk(kx,irot) ! index for rotated k in the FBZ
C          qibz_k= qibz(:,kx)
Cc          qbz_kr= qbz (:,kr)
C          if(kr/=0) qbz_kr= qbz (:,kr) !feb2006
C        else
Cc          k = 1  ! corresponds to q=0
Cc          kr= 1  ! corresponds to q=0
Cc          k = iqindx((/0d0,0d0,0d0/), ginv, qibz,nqibz)
Cc          kr= iqindx((/0d0,0d0,0d0/), ginv, qbz,  nqbz)
C          kr=-99999 !for sanity check
C          qibz_k= 0d0
C          qbz_kr= 0d0
C        endif
Cc        ngc = ngcni(k)  ! k-points in IBZ
Cc        print *, ' k ngc=',k,ngc
Cc        ngb = nbloch + ngcni(k)
C
C        call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
Cc        ngc = ngcni(k)  ! k-points in IBZ
C        ngb = nbloch + ngc
C
Cc --- exchange
C        if(exchange) then
Cc          allocate(vcoul(ngb,ngb))
Cc          read(ifvcfpout) vcoul(1:ngb,1:ngb)
C          read(ifvcfpout) nn !oct2005
C          allocate(vcoul(nn,nn))
C          read(ifvcfpout) vcoul(1:nn,1:nn)
C        endif
C
Cc- weight check for cycle or not.
C        if( kx <= nqibz ) then
C          if (kr == 0)    then
Cc            stop 'wmat: kr=0'
C            if(exchange) deallocate(vcoul)
C            cycle
C          endif
C          kount(kx)= kount(kx) + 1     ! count the no. times k
C          ! appears in the 1st BZ
Ccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc         print *,' irot,ip, k, kount in  =',irot, ip, k, kount(k,ip)
Cc          deallocate(vcoul)
Cc          cycle
Cc          print *,' kount out =',kount(k)
Cccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc          if (kount(kx) > nstar(kx)) stop 'sexc: too many stars'
C          if (kount(kx) > nstar(kx)) stop 'wmat: kount > nkstar'
Cc          if (kount(kx) > 1) stop 'wmat: kount > 1'
C        else
C          if( wgt0(kx-nqibz,irot) == 0d0 ) then
C            if(exchange) deallocate(vcoul)
C            cycle
C          endif
C        endif
C
Cc---test
C        if(OnlyQ0P.and.kx<=nqibz) then
C          if(exchange) deallocate(vcoul)
C          cycle
C        endif
C
Cc phase factor for off-site W
C        do ir1=1,nrws1
C          expiqR1(ir1) = exp(-img*tpi* sum(qbz_kr(:)*rws1(:,ir1)))
C        enddo
Cc      tmp =  sum(qbz_kr*rsite)
Cc      write(*,*)'tmp',tmp
Cc      write(*,*)'qbz_kr',qbz_kr
Cc      write(*,*)'rsite',rsite
C
Cc loop over FBZ
C        allocate( rmelt3(ngb,nwf,nwf,nrws2),cmelt3(ngb,nwf,nwf,nrws2))
C        rmelt3 = 0d0
C        cmelt3 = 0d0
C
C        do iq = 1,nqbz
Cctm debug      do iq = 1,0
C
C          q(:) = qbz(:,iq)
C
Cc read eval and eigen
C          call readeval (qbz(:,iq), isp, eq)
C          call readcphi (qbz(:,iq), nlmto,isp, quu, cphiq)
Cctm debug
Cc      write(7400,*)iq,nlmto,nwf
Cc      write(7400,*)cphiq
C          if(debug) print *, ' sxcf: 2'
C
Cc qk = q - rk, rk is inside 1st BZ, not restricted to
Cc the irreducible BZ
Cc      call vminv   (q,qbz(1,kr),3,qk)
C          qk =  q - qbz_kr ! qbz(:,kr)
Cc        iqk = iqindx(qk, ginv,qbz,nqbz)
C
Cc bring q-rk into the 1st BZ ==> qfbz labelled by kp
Cc         call fbz2(qk, ginv,index_qbz,n_index_qbz, qbz,nqbz, !qbz=rk
Cc     o                qfbz,kp)
C
Cc        if(debug)  print *,'xxx2',kp,qk,n1,n2,n3,qbas,ginv
Cc read eigenvalues,b and hb at q-rk
Cc        call rwdd1   (ifev, kp,nband,ekq)
C          call readeval(qk, isp, ekq)
C
Cc        ikpisp=isp + nsp*(kp-1)
Cc        call rwdd1  (ifcphi, ikpisp,2*nlmtobnd, cphikq)
C          call  readcphi(qk, nlmto,isp, quu, cphikq)
C
C
Cc====================================================
Cc matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)>
Cc including the phase factor exp(ik.T)
Cc B(rot*k,r) = B(k,invrot*r)
Cc====================================================
C
Cc> phase factors exp(ik.T(R))
Cc        call cv      (tpi,qibz(1,k),3,qfbz)
Cc        call expkt   (qfbz,tr,natom,
Cc     o              w(icoskt),w(isinkt))
Cc equivalent with above phase
C          do ia = 1,natom
C            expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) !  write(6,'(" phase ",i3,2d12.4)')ia,expikt(ia)
C          end do
Cc
C          if(debug) print *, ' sxcf: tr=',tr
C          if(debug) print *, ' sxcf: goto psicb2'
Cc> matrix elements
Cc>> core
C          nbmax = nwf
C          nt   = nctot + nbmax ! = nstate for the case of correlation
C          ntp0 = nwf
C          allocate( zzmel(nbloch,nt,ntp0))
C          zzmel = 0d0
Cc>> valence
C          call psi2b_v2 (nbmax, ntp0, iclass,
C     i    dreal(expikt(1:natom)),dimag(expikt(1:natom)),
C     i              cphikq, !rbkq,cbkq,rhbkq,chbkq, !  q-rk nko
C     i              cphiq,  !rbq,cbq,rhbq,chbq,     !  q    nko
C     i              ppb,  !pdb,dpb,ddb,
C     i              nlnmv,nlnmc,mdim,nctot,
C     i              imdim,iatomp,
C     d              mdimx,nlmto,nbloch,nlnmx, nband, nt,ntp0,
C     d              natom,nclass,
C     o              zzmel) ! rmel,cmel)
C
C          if(debug) write(6,"('sum of zmel abszmel=',4d23.16)")
C     &                 sum(zzmel),sum(abs(zzmel) )
C
Cc-plane
Cc        allocate( rmelt(ngb, nctot+nbmax, ntp0), ! nstate= nctot+nband
Cc     &            cmelt(ngb, nctot+nbmax, ntp0))
C          allocate( zmelt(ngb, nctot+nbmax, ntp0) )
C          if(debug) print *, ' sxcf_fal2: goto drvmelp'
Cc        call drvmelp( iq,   ntp0, ! q in FBZ
Cc     i                kp, nbmax,  ! q-rk
Cc     i                k,          ! k in IBZ for e-product basis
C
C
Cc        call drvmelp2( q,   ntp0, ! q in FBZ
Cc     i                 q-qbz_kr, nbmax,  ! q-rk
Cc     i                 qibz_k,           ! k in IBZ for e-product basis
Cc     i                 isp,ginv,
Cc     i       ngc,ngcmx,ngpmx,nband,itq,
Cc     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
Cc     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
Cc     o       rmelt,cmelt)
C
C          call drvmelp3( q,   ntp0, ! q in FBZ
C     i                 q-qbz_kr, nbmax,  ! q-rk
C     i                 qibz_k,           ! k in IBZ for e-product basis
C     i                 isp,ginv,
C     i       ngc,ngcmx,ngpmx,nband,itq,
C     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
C     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
C     o       zmelt)
C
C          deallocate(zzmel) !rmel,cmel)
Cc
C          if(debug) print *, ' sxcf: goto wtt'
C          if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)") sum( rmelt),sum(cmelt)
C
C          do ir2=1,nrws2
C            expiqR2 = exp( img*tpi* sum(q(:)*rws2(:,ir2)))
Cc rmelt3,cmelt3
C            rmelt3(:,:,:,ir2) = rmelt3(:,:,:,ir2)
C     &                     + wk(iq) * dreal(zmelt(:,:,:)*expiqR2)
C            cmelt3(:,:,:,ir2) = cmelt3(:,:,:,ir2)
C     &                     + wk(iq) * dimag(zmelt(:,:,:)*expiqR2)
C          enddo  ! ir2
C          deallocate(zmelt)
C
Cc end of loop over FBZ
C        end do
C
Cctm debug
Cc      do jj = 1,nwf
Cc        norm2 = 0d0
Cc        do ii = 1,ngb
Cc           norm2 = norm2 + (rmelt3(ii,jj,jj)**2+cmelt3(ii,jj,jj)**2)
Cc        enddo
Cc        write(*,*)'norm check',jj,norm2
Cc      enddo
C
C        if(bzcase()==2)then
C          if(kx<=nqibz) then
C            wtt = wk(kr)
C            if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
Cctm nqibz            wtt = wik(kx)*0.5d0
Cctm nqibz            if(nstbz(kr)/=0) wtt = wtt*(1d0-wgtq0p()/nstbz(kr))
C          elseif(kx>nqibz) then  !   wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
C            wtt= wgt0(kx-nqibz,irot)
Cctm nqibz            wtt= wqt(kx-nqibz)*wgtq0p()/dble(nqbz)
C          endif
C        else
C          if(kx<= nqibz) then
C            wtt = wk(kr)                   !         wtx = 1d0
Cctm nqibz            wtt = wik(kx)*0.5d0             !
C          else
C            wtt = wk(1)*wgt0(kx-nqibz,irot) ! wtx = wgt0(kx-nqibz,irot)
Cctm nqibz            wtt = wk(1)*wqt(kx-nqibz)  !
C            if(abs(wk(1)-1d0/dble(nqbz))>1d-10) stop 'sxcf:wk(1) inconsistent'
C          endif
C        endif
C
C        weight = wtt
Cctm        weight = weight * wik(kx) * 0.5d0 / wk(kr)
C
C        if(debug) then
C          write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
C        endif
C
Cctm, phase factor arising from off-site W
Cc        weight = weight
Cc        weightc = weight
C        do ir1=1,nrws1
C          weightc(ir1) = weight*expiqR1(ir1)
C        enddo
C
Cc--------------------------------------------------------
Cc --- bare Coulomb section ---
Cc--------------------------------------------------------
Cc
Cc S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
Cc                        v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
Cc
Cc> z1p(j,t,t') = S[i=1,nbloch] <psi(q,t') | psi(q-rk,t) B(rk,i)> v(k)(i,j)
Cc
Cc
Cc      write(6,*)' vcoulsum=',sum(vcoul)
Cc      if(debug) write(6,*)'  sumz=',dcmplx(rmelt,cmelt),sum(vcoul)
Cc
C        if(exchange) then
C
C          if (debug) write(*,*) 'bare coulomb section begins'
Cc
C          allocate(z1p(ngb,nwf,nwf))
C          allocate(zmel1(ngb))
C          allocate(zmel(ngb, nwf, nwf),vcoult(1:ngb,1:ngb))
C          do ir2=1,nrws2
Cc> zmel(j,t,t') = <psi(q,t') | psi(q-rk,t) B(rk,j)>
C            zmel  = dcmplx (rmelt3(:,:,:,ir2),cmelt3(:,:,:,ir2))
C            vcoult= transpose(vcoul)
Cc
Cc          deallocate(vcoul)
C
Cc          allocate( z1p(ngb,nwf,nwf))
C            call matm( vcoult, zmel,
C     &           z1p, ngb,ngb,nwf*nwf )
Cc          deallocate(vcoult, zmel)
Cc          allocate(zmel1(ngb))
C
Cc         print *,' nctot nbmax=',nctot,nbmax
C
C
C            do ir3=1,nrws2
Cc
C              do itp2 = 1,nwf
C                do it2  = 1,nwf
C                  do it   = 1,nwf
C                    do itp  = 1,nwf
C
C                      zmel1(:)=dcmplx(rmelt3(:,it,itp,ir3),-cmelt3(:,it,itp,ir3))
C                      ztmp = sum ( z1p(:,it2,itp2)*zmel1 )
C                      do ir1=1,nrws1
C                        ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
C                        rw_w(itp2,it2,it,itp,ir,0) = rw_w(itp2,it2,it,itp,ir,0)
C     &                                 + real(ztmp*weightc(ir1))
C                        cw_w(itp2,it2,it,itp,ir,0) = cw_w(itp2,it2,it,itp,ir,0)
C     &                                 + imag(ztmp*weightc(ir1))
C                      enddo ! ir1
C                    enddo
C                  enddo
C                enddo
C              enddo
Cc
C            enddo ! ir3
C          enddo ! ir2
C          deallocate(vcoul)
C          deallocate(vcoult, zmel)
C          deallocate(z1p,rmelt3,cmelt3,zmel1)
C
C          if (debug) write(*,*) 'bare coulomb section finished'
C
Cc-- End of bare-Coulomb section --------------
C
C        else
Cc--------------------------------------------------------------------------
Cc--- screening effect section----------------------------------------------
Cc--------------------------------------------------------------------------
C
Cc--------------------------------------------------------------
Cc S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
Cc                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n') |psi(q,t')>
Cc--------------------------------------------------------------
C
C
Cc--- The matrix elements zmel.
C          allocate( zmelc (ngb, nwf, nwf,nrws2),
C     c            zw (nblochpmx,nblochpmx),
C     c            zw2(nwf,nwf,nwf,nwf) )
Cc> zmelc(j,t,t') = < B(rk,j) psi(q-rk,t) | psi(q,t') >
C          zmelc = dcmplx (rmelt3,-cmelt3)
C          deallocate(rmelt3,cmelt3)
C          if(debug) print *,' end of zmel'
C
C
Cc====================================================================
Cc Wc(qt,w) along the imaginary axis
Cc====================================================================
C
Cc------------------------------------------------
Cc loop over w' = (1-x)/x, frequencies in Wc(k,w')
Cc {x} are gaussian points between (0,1)
Cc------------------------------------------------
Cc
C          nx  = niw
Cc        if(niw <1) stop " sxcf:niw <1"
Cc
C          do ix = 1,nx     ! imaginary frequency w'-loop
C            nrec=(kx-2)*niw+ix
C            if(bzcase()==2) nrec= (kx-1)*niw+ix
C            read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis
Cc          read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! Readin W-v on imag axis
C
Cc zwz= S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
Cc                Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
Cc        do itp = 1,ntp0
Cc        do  it = 1,nstate
Cc          zwz(ix,it,itp) = sum(
Cc     &   dconjg(zmel(:,it,itp)),matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)) )
Cc        enddo
Cc        enddo
C
C            do ir3=1,nrws2
C              do ir2=1,nrws2
C                call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2),
C     i                                   zmelc(:,:,:,ir3),
C     i      nwf,nwf,ngb,
C     o      zw2)
C
C                do ir1=1,nrws1
C                  ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
C                  rw_iw(:,:,:,:,ir,ix)
C     &  = rw_iw(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
C                  cw_iw(:,:,:,:,ir,ix)
C     &  = cw_iw(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
C                enddo ! ir1
C              enddo ! ir2
C            enddo ! ir3
C
C          enddo
C
C
Cc====================================================================
Cc Wc(qt,w) along the real axis
Cc====================================================================
C
C          if(debug) print *,' go to poles'
C
Cc        do      ix = 0,nw                    ! real frequency w'-loop
C          do      ix = 0,nrw                    ! real frequency w'-loop
C            nrec=(kx-2)*(nw+1-nw_i)+ ix-nw_i+1
C            if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
C            read(ifrcw,rec=nrec) zw
C
Cc zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
C            do ir3=1,nrws2
C              do ir2=1,nrws2
C                call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2),
C     i                                   zmelc(:,:,:,ir3),
C     i      nwf,nwf,ngb,
C     o      zw2)
C
C                do ir1=1,nrws1
C                  ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
C                  rw_w(:,:,:,:,ir,ix)
C     &  = rw_w(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
C                  cw_w(:,:,:,:,ir,ix)
C     &  = cw_w(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
C                enddo ! ir1
C              enddo ! ir2
C            enddo ! ir3
C
Cc debug: >
Cc          write(*,*)'nw,nrw=',nw,nrw
Cc          write(*,*)'ix=',ix
Cc          write(*,*)'nrws=',nrws
Cc          do i=1,5
Cc          write(*,*)i,rw_w(i,i,i,i,1,ix)
Cc          enddo
Cc          stop 'debug: 1'
Cc debug: <
C
C          enddo
C
C          deallocate(zmelc,zw,zw2)
C
C          if(debug) print *,' end of screen-if'
C
Cc end of if (exchange)
C        endif
C
C
C 1100 continue  ! end of k-loop
C

      return
      end

c -------------------------------------------------------------------
      subroutine matzwz2(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,itp2,it2,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),
     c              zwz(ntp0,nstate,nstate,ntp0)
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,nstate,ntp0) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
      do itp2 = 1,ntp0
        do  it2 = 1,nstate
          do  it  = 1,nstate
            do itp  = 1,ntp0
              zwz(itp,it,it2,itp2)
     & = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
            enddo
          enddo
        enddo
      enddo
      deallocate(CC)
      end
c -------------------------------------------------------------------
      subroutine matzwz3(zw,zmel1,zmel2, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,itp2,it2,ngb
      complex(8) :: zw(ngb,ngb),zmel1(ngb,nstate,ntp0),
     c                          zmel2(ngb,nstate,ntp0),
     c              zwz(ntp0,nstate,nstate,ntp0)
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,nstate,ntp0) )
      call matm(zw,zmel2,cc, ngb, ngb, nstate*ntp0)
      do itp2 = 1,ntp0
        do  it2 = 1,nstate
          do  it  = 1,nstate
            do itp  = 1,ntp0
              zwz(itp,it,it2,itp2)
     & = sum( dconjg(zmel1(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
            enddo
          enddo
        enddo
      enddo
      deallocate(CC)
      end
c -------------------------------------------------------------------
      subroutine drvmelp3( q, ntp0, q_rk,nt0, qik, isp,ginv,
     i        ngc,ngcmx,ngpmx,nband,itp0,
     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,noccx,nctot,
     o       zmelt)
      use m_readqg,only: readqg
      use m_readeigen,only: readgeig
c ! q1=q    ---> iq 1:ntp0 q-point
c ! q2=q-rk ---> kp 1:nt0  occupied
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c  zmelpl
c-----------------------------------------------------
      implicit none
      real(8):: q(3),q_rk(3),qik(3),ginv(3,3)
      integer(4):: ngp1, ngp2, ngpmx,nqbz,nqibz,  !iq,kp,k,
     &      ngcmx ,nctot,nband,
     &      ntp0,nt0,nbloch,noccx,  itx, ngc,nnum,inum,ig1,ig2,igc,
     &             ngvecpB1(3,ngpmx),
     &             ngvecpB2(3,ngpmx),
     &             ngveccBr(3,ngcmx), itp0(ntp0),
     &     nadd(3),isp  !,ngpn(nqbz)
      complex(8),allocatable::  zmelpl(:,:,:),geigq(:,:)
      real(8):: qbas(3,3),shtv(3),qdiff(3),add(3)
     &         ,qibz(3,nqibz),qbz(3,nqbz),qbasinv(3,3),symope(3,3)
     &         ,rmel(nbloch,noccx,ntp0)
     &         ,cmel(nbloch,noccx,ntp0)
c     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
c     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     & ,pi=3.1415926535897932D0
!      complex(8) :: geigB  (ngpmx_in,nband,nqbz)
      complex(8):: geig1(ngpmx,nband),geig2(ngpmx,nband)
     &            ,zmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      integer(4)::verbose
c
c      print *,' drvmelp',ngpn(iq),ngpn(kp),ngc
      if(verbose()>50) print *,' drvmelp2:'
c      ngp1 = ngpn(iq)  ! q        in FBZ
c      ngp2 = ngpn(kp)  ! q-rk     in FBZ
      call readqg('QGpsi', q,   ginv, qt,   ngp1, ngvecpB1)
      call readqg('QGpsi', q_rk,ginv, q_rkt,ngp2, ngvecpB2)
      call readgeig(q,    ngpmx, isp, qu1, geig1)
      call readgeig(q_rk, ngpmx, isp, qu2, geig2)
c debug
c      geig1=0d0
c      geig2=0d0

      if(sum(abs(qt-qu1))>1d-10) stop 'drvmelp2;qu1/=qu1x'
      if(sum(abs(q_rkt-qu2))>1d-10) stop 'drvmelp2;qu2/=qu2x'
      if(verbose()>=100) print *,' end of read geig '

c     ngc  = ngcni(k)  ! k-points in IBZ
c      allocate( geigq(ngp1,ntp0) )
c      do itx = 1,ntp0
c        geigq(1:ngp1,itx) = geigB(1:ngp1,itp0(itx), iq)
c      enddo

c     print *,' drvmelp x1',symope,qibz(1:3,k),qbz(:,iq),qbz(:,kp)
c      qdiff = matmul(symope,qibz(1:3,k)) - qbz(:,iq) + qbz(:,kp)
      qdiff = matmul(symope,qik)  - qt + q_rkt
      ! rk    -q  +(q-rk) is not zero.
      add  = matmul(qbasinv,qdiff)
      nadd = idint( add + dsign(.5d0,add))
      if(sum(abs(add-nadd))>1d-7) stop "drvmlp: abs(add-nadd))>1d-10"
c     print *,' drvmelp qdiff',qdiff
c      write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)')
c     &  add,nadd,abs(add-nadd)

      allocate( zmelpl(ngc,nt0,ntp0) )

cccccccccccccccccccccccccccccccccccccccccccccccc
c      if(.false.) then !old version
c      call melpln( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           , ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           , ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c     i       shtv, matmul(symope,qibz(1:3,k)), qbas,
c     o       zmelpl)
c      else             !new version from Sep-10 2000
c       call melpln2( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           ,  ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           ,  ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c     i       shtv, matmul(symope,qibz(1:3,k)),qibz(1:3,k), symope, qbas,
c     o       zmelpl)
c       endif
      if (nt0 /= ntp0) stop 'drvmelp3: nt0 /= ntp0'
      call melpln2( ngp1, ngvecpB1
     &           ,   ngp2, ngvecpB2
c     &           ,   ngc,  ngveccBr, nadd,
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,  ! q1=q    ---> iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,itp0(1:nt0)), nt0,     ! q2=q-rk ---> kp ngp2 1:nt0  occupied
     i       shtv, matmul(symope,qik),qik, symope, qbas,
     o       zmelpl)
c      deallocate(geig1,geig2)
c      rmelt=0d0
c      cmelt=0d0
      zmelt=0d0
c      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
c     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
c      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
c     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
      zmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = dcmplx
     &(rmel (1:nbloch, 1:nctot+nt0, 1:ntp0),
     & cmel (1:nbloch, 1:nctot+nt0, 1:ntp0))
c extensiton to nbloch+ngc
c      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c     & = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
c      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
      zmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     &=zmelpl(1:ngc,               1:nt0,            1:ntp0)
      deallocate(zmelpl)
      end



c------------------------------------------------------
