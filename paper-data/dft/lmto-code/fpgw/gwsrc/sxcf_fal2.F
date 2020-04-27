C Remove: kount, nstar, alat, niwx, ifvcfpout
      subroutine sxcf_fal3z(kount,ixc,deltaw,shtw,qip,itq, ntq,ef,ef2,esmr,esmr2,
     i          tiat,miat, nsp,isp,
c     i          ifrcw,ifrcwi,
     i                  qbas,ginv,
     i                  qibz,qbz,wk,nstbz,wik,
     i                  nstar,irkip,
     i                  iclass,mdim,nlnmv,nlnmc,
     i                  icore,ncore,imdim,
     i                  ppb,
     i                  freq_r,freqx,wx,
     i                  dw,ecore,
     d                nlmto,nqibz,nqbz,nctot,
     d                  nl,nnc,nclass,natom,
     d                  nlnmx,mdimx,nbloch,ngrp, nw_i,nw ,niw,niwx,nq,
     &     nblochpmx ,ngpmx,ngcmx,
     &     wgt0,nq0i,q0i,symgg,alat, shtvg,nband, ifvcfpout,
     &     exchange,tote,screen,cohtest, ifexsp,
     i     omega,iwini,iwend,
     i     nbmx,ebmx,
     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo,
     i   wklm,lxklm, !mar2012 !wklm apr2012takao
     i invg,il,in,im,nlnm,nn_,lx,nx_,nxx_,cgr,ppbrd, !sep 2012
     o zsec,coh,exx)
      use m_readqg
      use m_readeigen,only: readeval,readcphi
      use keyvalue
      use mpi_mod, q_comm => comm_1, q_inter => inter_1,
     &             k_comm => comm_2, k_inter => inter_2,
     &             rank => mpi_rank
      implicit none
CMay2006 TimeReversal off. when nw_i is not zero.

C- Calcualte diagonal part simga_ii(e_i)= <i|Re[S](e)|i> -------
Co zsec: S_ij= <i|Re[S](e)|i> where e=e_i and e_i \pm deltaw
Co
Cr  exchange=T : Calculate the exchange self-energy
Cr          =F : Calculate correlated part of the self-energy
Cr
Cr
Cr---- 2001 Sep. esec=omega(itp,iw). Genral iw mode for exchange =F
Cr 2000 takao kotani. This sxcf is starting from sec.f F.Aryasetiawan.
C---------------------------------------------------------------


c---- original document for sce.f (correlation case) by ferdi.Aryasetiawan.
c 92.02.24
c 93.10.18 from sec.f modified to take into account equivalent atoms
c calculates the correlated part of the self-energy SE
c SEc(q,t,t') = <psi(q,t) |SEc| psi(q,t'>
c SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >

c the zeroth order Green function
c G(r,r';w)   = S[occ]   psi(kn,r) psi(kn,r')^* /(w-e(kn)-i*delta)
c             + S[unocc] psi(kn,r) psi(kn,r')^* /(w-e(kn)+i*delta)

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

c SEc(q,t,t') = S[k=FBZ] S[n=occ]   S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)-i*delta)>
c
c             + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)+i*delta)>

c the analytic structure of GWc for w <= ef
c                               |
c                               |   o = pole of G
c                               ^   x = pole of Wc
c                               |
c                               |   ef-w
c                               |----<-----
c                               |          |
c                 o  o  o  o  o |o  o  o   ^
c               x  x  x  x  x  x|          |
c  -----------------------------|---->------------------------------
c                               |x  x  x  x  x  x  x  x
c                               |              o  o  o  o  o
c                               |       <----->
c                               ^        gap in insulator
c                               |
c                               |

c the analytic structure of GWc for w > ef
c                               |
c                               |   o = pole of G
c                               |   x = pole of Wc
c                               |
c         gap in insulator      ^
c                <----->        |
c      o  o  o  o               |
c         x  x  x  x  x  x  x  x|
c  ------------------------>----|-----------------------------------
c                   |           |x  x  x  x  x  x  x  x
c                   ^   o  o  o  o  o  o  o
c                   |           |
c                    ------<----|
c                       w-ef    |
c                               ^
c                               |
c integration along the real axis from -inf to inf is equivalent to
c the integration along the path shown
c------------------------------------------------------------
c integration along the imaginary axis: wint (s. also wint.f) (takao ->wintz)
c   (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
c the i*delta becomes irrelevant
c------------------------------------------------------------
c
c omit k and basis index for simplicity and denote e(q-k,n) = e
c wint = (i/2pi) < [w'=-inf,inf] Wc(w')/(w+w'-e) >
c
c w' ==> iw', w' is now real
c wint = - (1/pi) < [w'=0,inf] Wc(iw') (w-e)/{(w-e)^2 + w'^2} >
c
c transform: x = 1/(1+w')
c this leads to a denser mesh in w' around 0 for equal mesh x
c which is desirable since Wc and the lorentzian are peaked around w'=0
c wint = - (1/pi) < [x=0,1] Wc(iw') (w-e)x^2/{(w-e)^2 + w'^2} >
c
c the integrand is peaked around w'=0 or x=1 when w=e
c to handel the problem, add and substract the singular part as follows:
c wint = - (1/pi) < [x=0,1] { Wc(iw') - Wc(0)exp(-a^2 w'^2) }
c                          * (w-e)/{(w-e)^2 +w'^2}x^2 >
c        - (1/2) Wc(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c
c the second term of the integral can be done analytically, which
c results in the last term
c a is some constant
c
c when w = e, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w') and
c the integral becomes -Wc(0)/2
c this together with the contribution from the pole of G (s.u.)
c gives the so called static screened exchange -Wc(0)

c--------------------------------------------
c contribution from the poles of G: SEc(pole)
c--------------------------------------------
c
c for w <= ef
c SEc(pole) = - S[k=FBZ] S[n=occ] S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c             Wc(k,e(q-k,n)-w)(i,j) theta(e(q-k,n)-w)
c
c for w > ef
c SEc(pole) = + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c             Wc(k,w-e(q-k,n))(i,j) theta(w-e(q-k,n))
c
c theta(x)  = 1   if x > 0
c           = 1/2 if x = 0
c           = 0   if x < 0

c FBZ = 1st BZ
c NOTE: the routine only calculates the diagonal elements of the SE
c       i.e. SEc(q,t)

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
c irk(k,R,nq) = gives index in the FBZ with k{IBZ, R=rotation
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
c secq(t) = <psi(q,t) |SEc| psi(q,t)>
c----------------------------------------------------------------------
c      implicit real*8(a-h,o-z)
      integer :: ntq, natom,nqbz,nqibz,ngrp,nq,nw,niw,
     &  nband,  nlmto, nq0i,nctot,nlmtobnd,nstate,nstatex,
     &   irot,  isp,nsp,  nlnmx,
     &  ip, it,itp,
     i  mdim(*),
     i  ifrcw,ifrcwi, !iindxk,
     &  ifvcfpout,kx,kr,ngc,ngb,nbloch, !n1,n2,n3, k,
     &  nt0,nocc, nt0p,nt0m,i,nt,ntp0,
     &  nbmax,nclass,nl,nnc, nblochpmx,ix,nx,iw,iwp,ixs,ixsmx, !iwp sf 21May02
     &  mdimx, nwx,niwx,
     &  itq(ntq),iatomp(natom),miat(natom,ngrp),
     &  nstar(nqibz),irkip(nqibz,ngrp,nq),kount(nqibz,nq)
c
      real(8) :: q(3),qbas(3*3),ginv(3*3),tiat(3,natom,ngrp),tr(3,natom),
     &  wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz),
     &  freqx(niw),wx(niw),!expa(niw),
     &  eq(nband,nq),
     &   ekq(nband), ekc(nctot+nband),
     &  tpi,ef,ef2,esmr,esmr2,efp,efm,wfac,wfacx,we,esmrx, !ua,
     &  dw,wtt,wexx,exx,exxq ,wfacx2,weavx2,xx
c      complex(8) :: zsec(-1:1,ntq,nq)
c      real(8)    ::  shtw
c                       ! This shft is  to avoid some artificial resonance effects.
c                       ! shtw can be zero for esmr/=0 given by takao.
c
      integer(4):: ngpmx, ngcmx,ntqxx
      real(8) :: wgt0(nq0i,ngrp),qk(3), !qfbz(3),
     &  det,symgg(3,3,ngrp),symope(3,3),
     &  qxx(3),q0i(1:3,1:nq0i),shtvg(3,ngrp),shftv(3),alat,ecore(nctot),
     &   ppb(1), !pdb(1),dpb(1),ddb(1), !*
     &  coh(ntq,nq)   !, pos(3,natom)
      complex(8)::   alagr3zz,zdotc !geigB  (ngpmx,nband,nqbz),

c
C      real(8),allocatable:: !rmel(:,:,:),cmel(:,:,:),
C     &                      rmelt(:,:,:),cmelt(:,:,:)
      complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:),zmelt(:,:,:),
     &     zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:)
c for exchange --------------------
      logical :: exchange,screen,cohtest,tote
      real(8),allocatable::
     &   w1p(:,:,:),w2p(:,:,:),w3p(:,:)
      complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)
      integer:: invg(ngrp),invrot,invr
      integer:: il(*),in(*),im(*),nlnm(*),nn_,lx(*),nx_(*),nxx_
      real(8)::  cgr(*),ppbrd(*)

c- debugwrite ---------------------
      logical :: debug=.false.
      logical :: initp=.true.

cccccccccccccc
c tetra
c      integer(4) :: ntqx
c      integer(4) :: ibzx(nqbz)
c      real(8)    :: wtet  (nband,nqibz,1:ntqx), wtetef(nband,nqibz)
c                ! where the last index is 3*itq+iw-1,itq=1,ntq,iw=-1,1
c     logical    :: tetraex
cccccczzcccccc

      complex(8) :: wintzav,wintzsg_npm

      integer(4) :: iii,ifexsp ,iopen
C     integer(4),save::ifzwz=-999

      integer(4) :: iwini, iwend, ia
      real(8)    ::  omega(ntq, iwini:iwend)
      complex(8) :: zsec(iwini:iwend,ntq,nq)
      complex(8),allocatable:: expikt(:)
      complex(8):: img=(0d0,1d0)
ctakao
      complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)

      integer(4) :: nt_max,igb2,  nw_i  !nw_i is at feb2006 TimeReversal off case
      complex(8),allocatable:: zmel1(:),zmel3(:)
      complex(8), allocatable :: zw_(:,:) !,zzmel(:,:)
      complex(8), allocatable :: zwz2(:,:),zmel2(:,:) !0 variant
      complex(8) ::  zz2 ,zwz3(3)
      real(8) :: omg
      real(8) :: freq_r(nw_i:nw)
      complex(8), allocatable :: zw3(:,:,:)


      real(8)::wfaccut=1d-10

      logical :: GaussSmear,gass
      real(8) :: ebmx,ddw
      integer(4):: nbmx,nbmxe,nstatetot

c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

      integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

      integer(4)::verbose,iverbose,idalloc
      integer(4)::nstbz(nqbz),bzcase,iqini,iqend
      real(8):: wgtq0p

      integer(4):: nrec
      real(8)::quu(3),qibz_k(3),qbz_kr(3)
      logical :: onlyQ0P,noq0p,ltmp

      logical ::zwz3mode


      real(8):: ua_,expa_(niw),freqw,freqw1,ratio,ua2_(niw)
      logical :: ua_auto
      real(8),allocatable:: uaa(:,:)

c      logical ::testimx=.false.
ccccc zvz test cccccccccccccccccccccccccc
c      complex(8):: vcoul(ngbx,ngbx)
      complex(8),allocatable::  zwzs(:)
ccccccccccccccccccccccccccccccccccc
      integer(4) ::nbcut,nbcutc
      logical :: oncew


C...
      logical::smbasis
      integer(4):: nn,no,iclose
      complex(8),allocatable:: pomat(:,:)
      real(8):: q_r(3)
      integer(4):: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo)
      complex(8):: pomatr(nnmx,nomx,nkpo)
      real(8):: qrr(3,nkpo)

      integer(4):: ixsmin,iir,nwxi, itini,itend, npm
      real(8)   :: fffr(3)
      complex(8):: zwzz(3)

      real(8),allocatable:: ebb(:)
      integer(4):: ii
      logical ::evaltest !, imgonly

      integer:: lxklm,ivc,ifvcoud,ngb0
      real(8):: wklm((lxklm+1)**2),pi,fpi,vc,qvv(3)
C     logical :: newaniso !,testomitq0
      complex(8),allocatable:: ppovl(:,:),ppovlz(:,:),zcousq(:,:)
      real(8),allocatable::vcoud(:)
      integer:: mrecl,nprecx,ifwd
      character(5):: charnum5

      integer:: ixc
      real(8):: qip(3,*),deltaw,shtw,eqx(nband)

      integer :: num_q_groups, num_k_groups, ipstart, ipend, ipstep, ikstart, ikend, ikstep

      print *,'sxcf_fal3z'
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      debug=.false.
      if (verbose()>=90) debug=.true.
!     debug=.true.
      if(.not.exchange) then !screened exchange
        if (rank == 0) then
          ifwd = iopen('WV.d',1,-1,0)
          read (ifwd,*) nprecx,mrecl
          ifwd = iclose('WV.d')
        endif
#ifdef USE_MPI
        call MPI_Bcast(nprecx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(mrecl, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
      endif ! .not.exchange
      if (rank == 0) call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
#ifdef USE_MPI
      call MPI_Bcast(nbcut,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif

      nbcutc = nctot+nbcut
      if(debug) then
        if(rank == 0) then
C         if(GaussSmear()) then
          gass=GaussSmear()
C         write(6,"(' sxcf_fal2z: q  =',3f8.4)")q
          write(6,*)' sxcf_fal2z: new version: gaussian smearing=',gass
C         endif
        endif
      endif
      tpi = 8d0*datan(1d0)
      if (nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
      nlmtobnd  = nlmto*nband
      nstatetot = nctot + nband
      allocate(expikt(natom))

!!== ip loop uses symmetry to accelerate integration over external q ==
      num_q_groups = mpi_N1
      call mpi_make_groups(comm_world,nq,num_q_groups,q_comm,q_inter)

#ifdef USE_MPI
      write(*,'("MPI: rank = ",I4,", q_comm%ID = ",I4,", q_comm%group_ID = ",I4)')
     &          rank, q_comm%ID, q_comm%group_ID

      write(*,'("MPI: rank = ",I4,", q_inter%ID = ",I4,", q_inter%group_ID = ",I4)')
     &          rank, q_inter%ID, q_inter%group_ID
#endif

      num_k_groups = mpi_N2
      call mpi_make_groups(q_comm, iqend-iqini, num_k_groups, k_comm, k_inter)
      ipstart = 1 + q_comm%group_ID
      ipend   = nq
      ipstep  = num_q_groups

      do  1001 ip = ipstart, ipend, ipstep
        if (sum(irkip(:,:,ip))==0) cycle
        q(1:3) = qip(1:3,ip)
        if (q_comm%ID == 0) then
          write(*,"(' do 1001: iq=',i3,' of',i3,'  q=',3f12.6)") ip,nq,q
          call cputid(0)
        endif

        if(ixc==2) then
c     iq  = idxk (q(1,ip),qbz,nqbz)
c     call rwdd1 (ifev(is),iq, nband,eq)
        call readeval(q,isp,eqx)
        do iw = iwini,iwend
          do i  = 1,ntq
            omega(i,iw) = eqx(itq(i)) + 2d0*(dble(iw)-shtw)*deltaw
          enddo
        enddo
        endif
c       read eval and eigen
        call readeval(q,isp,eq(1,ip))
        call readcphi(q,nlmto,isp,quu,cphikq) ! ecal Apr 19 no readcphi
        if (tote) exxq = 0.d0

c ecal Apr 19
C        if(ixc==4) then
Cc        dwplot=0.01
C        do iw = iwini,iwend
C          omega(1:ntq,iw) =  dwplot* iw + ef
C        enddo
C        endif

        do  it = 1, ntq ! augmentation part of eigenfunction for itp at q
          itp = itq(it)
          cphiq(1:nlmto,it) = cphikq(1:nlmto,itp)
        enddo

C       if (bzcase()==1) then
        if (abs(sum(qibz(:,1)**2))/=0d0) call rx(' sxcf assumes 1st qibz/=0 ')
        if (abs(sum( qbz(:,1)**2))/=0d0) call rx(' sxcf assumes 1st qbz /=0 ')
C       endif

C        Not used any more
C        if(exchange.and.(.not.newaniso())) then
C          rewind  ifvcfpout
C          read(ifvcfpout) ndummy1, ndummy2
C        endif

C   --- gauss_img : interpolation gaussion for W(i \omega).
        if (q_comm%ID == 0) call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
#ifdef USE_MPI
        call MPI_Bcast(ua_,1,MPI_DOUBLE_PRECISION,0,q_comm%communicator,ierror)
#endif

        if (ua_ <= 0d0) then
          ua_auto = .true.
          write(6,"(' ua_auto=T')")
        else
          ua_auto = .false.
          do  ix = 1, niw
            freqw     = (1d0 - freqx(ix))/ freqx(ix)
            expa_(ix) = exp(-(ua_*freqw)**2)
          enddo
        endif

!! == Big loop for kx ==
!! kx is for irreducible k points, kr=irk(kx,irot) runs all k points in the full BZ.
!       call getkeyvalue("GWinput","TestOnlyQ0P",onlyq0p,default=.false.)
        iqini=1
        iqend=nqibz             !no sum for offset-Gamma points.
        ikstart = iqini + k_comm%group_ID
        ikend   = iqend
        ikstep  = num_k_groups
        do 1100 kx = ikstart, ikend, ikstep !  kx=1 corresponds to q=0 is omitted.
          if (sum(irkip(kx,:,ip))==0) cycle

          iverbose = (verbose()/30)*10
          if (kx /= iqini .or. ip > 1) iverbose = 0
          if(k_comm%ID == 0) then
            write (*,'(" do 1100 kx",i3,"  out of ",i3,"  k-points "$)') kx,iqend
            call cputid(0)
          endif
          if (kx <= nqibz) then
            qibz_k = qibz(:,kx)
          else
            qibz_k = 0d0
          endif

          call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
          ngb = nbloch + ngc
          if (debug .and. k_comm%ID == 0) print *, ' sxcf: ngb=',ngb,nbloch

!! === Read diagonalized Coulomb interaction from file Vcoud.ID  where ID=kx ===
!! note sep102012takao
!! Vcoud is a sequential file containing the coulomb matrix for qibz_k.
!! It is is written in hvccfp0.m.F.
!! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)

!         Always newaniso now
C         if(newaniso()) then
          if(kx<=nqibz) qxx=qibz_k
          if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
          if (k_comm%ID == 0) ifvcoud = iopen('Vcoud.'//charnum5(kx),0,0,0)
          do
            if(k_comm%ID == 0) then ! only k-head reads
              read(ifvcoud) ngb0
              read(ifvcoud) qvv
            endif
#ifdef USE_MPI
            call MPI_Bcast(ngb0, 1, MPI_INTEGER, 0, k_comm%communicator, ierror)
            call MPI_Bcast(qvv,  3, MPI_DOUBLE_PRECISION, 0, k_comm%communicator, ierror)
#endif
            if(allocated(vcoud)) deallocate(vcoud)
            allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
            if(k_comm%ID == 0) then
              read(ifvcoud) vcoud
              read(ifvcoud) zcousq
            endif
#ifdef USE_MPI
            call MPI_Bcast(vcoud,  ngb0,    MPI_DOUBLE_PRECISION, 0, k_comm%communicator, ierror)
            call MPI_Bcast(zcousq, ngb0**2, MPI_DOUBLE_COMPLEX, 0, k_comm%communicator, ierror)
#endif

            if(sum(abs(qvv-qxx))<1d-6) goto 1133
          enddo                 ! done reading vcoud, zcousq

          if(sum(abs(qvv-qxx))>1d-6) then
            print *,'qvv =',qvv
            print *,'qxx=',qxx,kx
            call rx( 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent')
          endif
 1133     continue
          if (ngb0/=ngb ) then  !sanity check
            write(0,'("[",I4,"] qxx = ",I5,", ngb0 = ",I6", ngb = ",I6)') rank,qxx,ngb0,ngb
#ifdef USE_MPI
            call MPI_abort(MPI_COMM_WORLD,-1,ierror)
#endif
            call rx( 'hsfp0.m.f:ngb0/=ngb')
          endif

!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
!! ppovlz= O Zcousq
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
          if (allocated(ppovl)) deallocate(ppovl,ppovlz)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
          if(k_comm%ID == 0) call readppovl0(qibz_k,ngc,ppovl)
#ifdef USE_MPI
          call MPI_Bcast(ppovl,ngc*ngc,MPI_Double_Complex,0,k_comm%communicator,ierror)
#endif
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
!         ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
          call zgemm('N','N',ngc,ngb,ngc,(1d0,0d0),ppovl,ngc,zcousq(nbloch+1,1),ngb0,(0d0,0d0),ppovlz(nbloch+1,1),ngb)
!         call zprm('ppovzl',ppovlz,ngb,ngb,ngb)
          deallocate(zcousq)
!         endif

!! read Coulomb interaction for exchange for newaniso=F
C        if(exchange) then
Cc          if(testomitq0()) then
Cc            print *,'xxxxxxxxxxx takao test omit q0i aaaaa'
Cc            if(allocated(vcoul)) deallocate(vcoul)
Cc          endif
C          if(.not.newaniso()) then
C            read(ifvcfpout) nn !oct2005
C            allocate(vcoul(nn,nn))
C            read(ifvcfpout) vcoul(1:nn,1:nn)
C          endif
C        endif

!! === open WVR,WVI for correlation mode ===
          if (.not.exchange) then
            if(k_comm%ID == 0) then  ! MPI:  files are only opened on head node of the k_groups.
              ifrcw  = iopen('WVR.'//charnum5(kx),0,-1,mrecl)
              ifrcwi = iopen('WVI.'//charnum5(kx),0,-1,mrecl)
            endif
          endif

!!    --- loop over rotations ---
          do 1000 irot = 1,ngrp
            if (kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if(kr==0) cycle   ! next irot
              qbz_kr = qbz (:,kr)
            else
              kr = -99999         !for sanity check
              qbz_kr = 0d0
              if (wgt0(kx-nqibz,irot) == 0d0) cycle ! next irot
            endif
            if(debug) write(*,"('ip,kx irot=',3i5,' out of',2i4)") ip,kx,irot,iqend,ngrp

c        if(testomitq0())then
c          print *,'xxxxxxxxxxx takao test omit q0i aaaaa'
c          if(kx==1.or.kx>nqibz) cycle
c        endif

!!          Rotate atomic positions invrot*R = R' + T
            invr  = invrot (irot,invg,ngrp)
            tr    = tiat(:,:,invr)
            iatomp= miat(:,invr)
            symope= symgg(:,:,irot)
            shftv  = matmul(symope,shtvg(:,invr))
c     -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
            call ppbafp_v2 (irot,ngrp,isp,nsp,
     i           il,in,im,nlnm, !w(i_mnl),
     d           nl,nn_,nclass,nlnmx,
     i           mdimx,lx,nx_,nxx_, !Bloch wave
     i           cgr, nl-1,  !rotated CG
     i           ppbrd,      !radial integrals
     o           ppb)

c        if( kx <= nqibz ) then
c          kr = irk(kx,irot) ! index for rotated k in the FBZ
c          qibz_k= qibz(:,kx)
c          if(kr/=0) qbz_kr= qbz (:,kr) !feb2006
c        else !!   k = 1  ! corresponds to q=0
c          kr=-99999 !for sanity check
c          qibz_k= 0d0
c          qbz_kr= 0d0
c        endif

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! check zcousq^* * ppovl * zcousq gives delta function --->
c$$$        print *,' ddddd zcou qvv qibz_k=',sum(abs(qvv-qibz_k)),qvv
Ci          :NEW: I needed to replace mkppolv2 in rdata4gw_v2. (takao2012may).
c$$$        do ix=1,ngb
c$$$        do iy=1,ngb
c$$$           aaaa=sum(dconjg(zcousq(:,ix))*ppovlz(:,iy))
c$$$           if( ix/=iy.and.abs(aaaa) >1d-3 ) then
c$$$            write(*,*)' dddd zcousum check',ix,iy,aaaa
c$$$           endif
c$$$           if( ix==iy.and.abs(aaaa-1) >1d-3 ) then
c$$$            write(*,*)' dddd zcousum check',ix,iy,aaaa
c$$$           endif
c$$$        enddo
c$$$        enddo
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cc---test
c        if(OnlyQ0P.and.kx<=nqibz) then
c          if(exchange) deallocate(vcoul)
c          cycle
c        endif

!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
c          call vminv(q,qbz(1,kr),3,qk)
            qk = q - qbz_kr     ! qbz(:,kr)

C        if(debug) then
C        write(6,"(' q      =',3f13.5)") q
C        write(6,"(' qbz_kr =',3f13.5)") qbz_kr
C        write(6,"(' qk     =',3f13.5)") qk
C        endif

c bring q-rk into the 1st BZ ==> qfbz labelled by kp
c         call fbz2(qk, ginv,index_qbz,n_index_qbz, qbz,nqbz, !qbz=rk
c     o                qfbz,kp)

c        if(debug)  print *,'xxx2',kp,qk,n1,n2,n3,qbas,ginv
c read eigenvalues,b and hb at q-rk
c        call rwdd1   (ifev, kp,nband,ekq)
!!          Read eigenvalue and Cphi at q-k
            call readeval(qk, isp, ekq)
c        ikpisp=isp + nsp*(kp-1)
c        call rwdd1  (ifcphi, ikpisp,2*nlmtobnd, cphikq)
c        call  readcphi2(pos,natom,qk, nlmto,isp, quu, cphikq)
            call readcphi(qk, nlmto,isp, quu, cphikq)



c        if(debug) write(6,"('111111 end of readcphi: sum of cphikq=',4d23.16)")
c     &                 sum(cphikq),sum(abs(cphikq) )

c no. occupied (core+valence) and unoccupied states at q-rk
c      call dcopy   (nband,ekq,1,ekc(nctot+1),1)
            ekc(nctot+1:nctot+nband) = ekq(1:nband)
            nt0 = nocc(ekc,ef,.true.,nstatetot)
            ddw = .5d0
            if (GaussSmear()) ddw = 10d0
            efp = ef+ddw*esmr
            efm = ef-ddw*esmr
            nt0p = nocc(ekc,efp,.true.,nstatetot)
            nt0m = nocc(ekc,efm,.true.,nstatetot)

ctetra
c        if(tetraex) then
c         irkp = ibzx(kp) kp??? ! IBZ corresponding to qk
c         do i=1,nband
c          if( wtetef(i,irkp) == 0d0 ) exit
c         enddo
c        if(sum(abs(wtetef(i:nband,irkp)))/=0d0)
c     &   stop ' sxcf: wtetef sum err1'
c         nt0org = nocc (ekc,ef,.true.,nstatetot)
c         nt0 = i-1 + nctot
c        endif

            if(exchange) then
              nbmax = nt0p-nctot
            else
              nbmax = nband  ! intermediate states => all states initially
              nbmxe = nocc (ekc,ebmx,.true.,nstatetot)-nctot !? ebmx(2)
              nbmax  = min(nband,nbmx,nbmxe) !?nbmx(2)
              if(initp) then
                if(k_comm%ID == 0) print *,' nbmax=',nbmax
                initp=.false.
              endif
            endif
            nstate = nctot + nbmax ! = nstate for the case of correlation
            if(debug) print *,'ngb, nstate, ntp0=',ngb,nstate,ntp0

c====================================================
c matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)>
c including the phase factor exp(ik.T)
c B(rot*k,r) = B(k,invrot*r)
c====================================================

c> phase factors exp(ik.T(R))
c        call cv      (tpi,qibz(1,k),3,qfbz)
c        call expkt   (qfbz,tr,natom,
c     o              w(icoskt),w(isinkt))
c equivalent with above phase
        do ia = 1,natom
          expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) !  write(6,'(" phase ",i3,2d12.4)')ia,expikt(ia)
        end do
c
        if(debug) print *, ' sxcf: tr=',tr
        if(debug) print *, ' sxcf: goto psicb2'

        ntp0 = ntq

        nt   = nctot + nbmax    ! = nstate for the case of correlation .... note nt identical to nstate above
C       Allocate zzmel once only; fill in parts.  No drealzzmel,dimagzzmel,rmelt,cemlt
        ia = idalloc('zzmel',iverbose+2,ngb*nt,ntq*2)
        allocate(zzmel(ngb,ntq,nt))
        ltmp = verbose()>=40 .and. kx == iqini .or. verbose()>=41

        call mkmelt(1,isp,isp,qbas,q,qbz_kr,qibz_k,cphiq,cphikq,ppb,icore,ncore,
     .    symope,shftv,tr,nl,natom,nclass,iclass,iatomp,
     .    nlnmx,nlnmv,nlnmc,imdim,mdim,mdimx,nlmto,ngc,nbloch,
     .    1,ntq,itq,1,1,nbmax,1,ngb,1,nt,1,ntq,nnc,ltmp,zzmel)

        if (debug) write(6,"('mkmet sum zzmel, abszzmel=',4d23.12)") sum(zzmel),sum(abs(zzmel))
C       call zprm('zzmel',zzmel,ngb,ngb,nstate*ntq)

! need pomatr, nnr, nor
!oct2005
        if(smbasis()) then !
          nn= nnr(kx)
          no= nor(kx)
          allocate( pomat(nn,no) )
          pomat= pomatr(1:nn,1:no,kx)
          if( sum(abs(qibz_k-qrr(:,kx)))>1d-10 .and.kx <= nqibz ) then
            call rx( 'qibz/= qrr')
          endif
          if(no /= ngb.and.kx <= nqibz) then
!A bit sloppy check only for kx<nqibz because qibze is not supplied...
            write(6,"(' q  ngb  ',3d13.5,3i5)")  qibz_k,ngb
            write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
            call rx( 'x0kf_v2h: POmat err no/=ngb')
          endif
ccccccccccccccccccccccccccc
ctttt
c       write(6,*)
c         write(6,*)'kx=',kx
c       write(6,"('ngb nn=',3i4)") ngb,nn,no
c         write(6,"(' q  ngb  ',3d13.5,3i5)")  qibz_k,ngb
c         write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
c       if(ngb/=nn) stop ' sxcf_fal2 ngb/=nn'
c       goto 1113
ccccccccccccccccccccccccccc

          call rx('rewrite smbasis matrix elements')
C          ngb = nn       ! Renew ngb !!!
C          allocate ( zmel  (nn, nctot+nbmax, ntp0) )
C          call matm( pomat, dcmplx(rmelt,cmelt), zmel,
C     &                nn, no, (nctot+nbmax)*ntp0 )
C          deallocate(rmelt, cmelt)
C          allocate( rmelt(ngb, nctot+nbmax, ntp0), !ngb is reduced.
C     &           cmelt(ngb, nctot+nbmax, ntp0) )
C          rmelt = dreal(zmel)
C          cmelt = dimag(zmel)
C          deallocate(zmel,pomat)


c       print *,' smooth mixed basis : augmented zmel'
c 1113 continue
        else
          nn=ngb
          no=ngb
        endif

        if( oncew() ) then
          write(6,"('  ngb nn no=',3i6)") ngb,nn,no
        endif

ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
c      allocate( zmel (ngb, nctot+nbmax, ntp0), aaa(ngb))
c      zmel = dcmplx (rmelt,-cmelt)
c      ifix=9090
c      open (ifix,file='zvz_test')
c      print *,' ======= goto zvz test  ========='
c      do ib1=1,nctot+nbmax
c      do ib2=1,ntp0
c! zmelt*vcoul
c        zvz1= sum( dconjg(zmel(1:nbloch,ib1,ib2))
c     &          * matmul(vcoul(1:nbloch,1:nbloch),zmel(1:nbloch,ib1,ib2)) )
c        aaa= matmul(vcoul,zmel(:,ib1,ib2) )
c        zvz= sum( dconjg(zmel(:,ib1,ib2))*aaa(:) )
c        write(ifix,"(2i4,  2d13.5,'  ',2d13.5,'  ',2d13.5)") ib1,ib2,zvz ,zvz1, zvz-zvz1
c      enddo
c      enddo
c      close(ifix)
c      stop '--- zvz test end------------------'
ccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccc
c A matrix element check in Si. See the corresponding section of x0kf.f
c      write(6,"(3f13.5)") q         !ntp0
c      write(6,"(3f13.5)") qbz(:,kr) !product basis
c      print *,' sxcf: sumcheck rmelt cmelt='
c     &   ,sum(rmelt(1:ngb,4,1))
c     &   ,sum(cmelt(1:ngb,4,1))
c     &   ,sum(abs(rmelt(1:ngb,4,1)))
c     &   ,sum(abs(cmelt(1:ngb,4,1)))
c      stop ' a matrix ele check test end'
ccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
c testmelt0
c      rmelt(1: nbloch, :,:)=0d0
c      cmelt(1: nbloch, :,:)=0d0
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccc
c      do  it  = 1,nstate
c      itp=1
c      it =1
c     do  itp = 1,ntp0
c      print *,' ngb nbloch ngc =',ngb,nbloch,ngcni(k)
c      do  ibl  = 1,ngb
c        write(6,"(' rmelt cmelt =',3i4,4d14.6)")it,itp,ibl
c     &  ,rmelt(ibl,it,itp) ,cmelt(ibl,it,itp)
c     &  ,rmel (ibl,it,itp) ,cmel (ibl,it,itp)
c      enddo
cccccccccccccccccccccccccccccccccccccccccccccc
C        deallocate(zzmel) !rmel,cmel)
c
        if(debug) print *, ' sxcf: goto wtt'
C       if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)") sum(rmelt),sum(cmelt)

C        if(bzcase()==2)then
C          if(kx<=nqibz) then
C            wtt = wk(kr)
C            if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
C          elseif(kx>nqibz) then  !   wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
C            wtt= wgt0(kx-nqibz,irot)
C          endif
C        else
          if(kx<= nqibz) then
            wtt = wk(kr)        !         wtx = 1d0
          else
            wtt = wk(1)*wgt0(kx-nqibz,irot) !       wtx = wgt0(kx-nqibz,irot)
            if(abs(wk(1)-1d0/dble(nqbz))>1d-10)call rx( 'sxcf:wk(1)inconsistent')
          endif
C       endif

        if(debug) then
          write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
        endif


!!--------------------------------------------------------
!! === exchange section ===
!!--------------------------------------------------------
c
c S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c                        v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
c
c> z1p(j,t,t') = S[i=1,nbloch] <psi(q,t') | psi(q-rk,t) B(rk,i)> v(k)(i,j)
!! NOTE: zmel(igb, nctot+nbmax, ntp0) ---> <phi phi  |igb>
        if(exchange) then
c --- screened exchange case ----
c          if(screen) then
c            ix = 1
c            nrec=(kx-iqini)*nw+ix
c            if(bzcase()==2) nrec= (kx-1)*nw+ix
c            read(ifrcw,rec=nrec) zw  ! Readin W(0) - v
c            vcoul = vcoul + zw(1:ngb,1:ngb) !c  screen test
c          endif

C         stop 'sxcf_fal3z zmel -> zzmel+'
C         allocate( zmel(ngb, nctot+nbmax, ntp0), w3p( nctot+nbmax,ntp0))
          allocate(w3p(nctot+nbmax,ntp0))
C         zmel  = dcmplx (rmelt,cmelt)

cccccccccccccccccccccccccccccc
          if(debug) then
          do it=1,nctot+nbmax
             write(6,"('wwwww0 ',i5,2f10.4)") it,sum(abs(zzmel(:,1,it)))
          enddo
          write(6,*)'sumcheck ppovlz=',sum(abs(ppovlz(:,:)))
          write(6,*)'sumcheck vcoud=',sum(abs(vcoud(:)))
          endif
ccccccccccccccccccccccccccccc

!! ==== ngb=1 at q=0 is for exp(iqr) basis. ====
!!
c$$$      if(kx>=2 .and. kx <=nqibz) then !for testxxxxxx
c$$$          print *,'test skip kx=',kx
c$$$          do itp = 1,ntp0
c$$$          do it  = 1,nctot+nbmax
c$$$            w3p(it,itp) =0d0
c$$$          enddo
c$$$          enddo
c$$$          deallocate(vcoult, zmel)
c$$$      else

C         if(newaniso()) then !apr2012takao
C Original
C            do 992 itp = 1,ntp0
C            do 993 it  = 1,nctot+nbmax
C              w3p(it,itp) = 0d0
C              do 994 ivc=1,ngb
C                  if(ivc==1.and.kx==1) then
C                    vc= wklm(1)* fpi*sqrt(fpi) /wk(kx)
Cc                    print *,'wklm(1) vc=',wklm(1),vc
C                  else
C                    vc= vcoud(ivc)
C                  endif
C                  zmelt1 = sum( zmel(:,it,itp) *ppovlz(:,ivc) )
C                  w3p(it,itp) = w3p(it,itp)+ vc * abs(zmelt1)**2
Cccccccccccccccccccccc test cccccccccccccccccccccccc
Cc              if(ivc==1) w3p(it,itp) = 0d0
Cccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc              if(kx==1.and. ivc==1.and.abs(ekc(it)- eq(itq(itp),ip)) < 1d-6 ) then
Cc              voltot=270.01139400000000
Cc              print *,'zzzzzz abszmelt=',abs(zmelt1)**2*voltot
Cc              endif
Cc
Cc$$$!! may2012takao test mode. this is useful to see what zmel means
Cc$$$!! be careful for complex conjugate...  zmel(ib,it,itp)= <psi_itp psi_it|ib>
Cc$$$              if(.false.) then
Cc$$$              if(kx==1.and. ivc==1.and.abs(ekc(it)- eq(itq(itp),ip)) < 1d-6 ) then
Cc$$$                voltot = 302.838874d0   !test for gas_pw_gw_lmfh case
Cc$$$                zmelt0 = sum( zmel(1:nbloch,it,itp)*ppovlz(1:nbloch,ivc) )*sqrt(voltot)
Cc$$$c                zmelt1 = sum( zmel(nbloch+1:nbloch+ngc,it,itp)*
Cc$$$c     &                        matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,ivc)))
Cc$$$c     &                        *sqrt(voltot)
Cc$$$c                print *,'it itp zmelt0 1 =',it,itq(itp),ekc(it),abs(zmelt0),abs(zmelt0+zmelt1)
Cc$$$                print *,'it itp zmelttot =',it,itq(itp),ekc(it),abs(zmelt0)
Cc$$$              endif
Cc$$$              endif
C 994          continue
C 993        continue
C 992        continue
C            if(verbose()>=30) call cputid2(' complete w3p',0)
C            call dprm('inline w3p',w3p,(nctot+nbmax),(nctot+nbmax),ntp0)
C vectorized
            ia = idalloc('zmelt',iverbose+2,ngb*nstate,ntp0*2)
            allocate(zmelt(ngb,ntp0,nstate))
            call zgemm('T','N',ngb,nstate*ntp0,ngb,(1d0,0d0),
     .        ppovlz,ngb,zzmel,ngb,(0d0,0d0),zmelt,ngb)
            if(verbose()>=40) call cputid2(' complete zzmel * overlap .',0)
C           call zprm('zmel*ppovlz+, vectorized',zmelt,ngb,ngb,nstate*ntp0)
            vc = vcoud(1)
            if (kx==1) vcoud(1) = wklm(1)* fpi*sqrt(fpi) /wk(kx)
            do  itp = 1, ntp0
              do  it  = 1, nstate
                w3p(it,itp) = 0d0
                do  ivc = 1, ngb
                  w3p(it,itp) = w3p(it,itp) + vcoud(ivc) *
     .              (dble(zmelt(ivc,itp,it))**2+dimag(zmelt(ivc,itp,it))**2)
                enddo
              enddo
            enddo
            ia = idalloc('zmelt',iverbose+4,ngb*nstate,ntp0*2)
            deallocate(zmelt)
            vcoud(1) = vc
            if(verbose()>=30) call cputid2(' complete w3p',0)
C           call dprm('vectorized w3p',w3p,(nctot+nbmax),(nctot+nbmax),ntp0)

C         endif
          it = idalloc('zzmel',iverbose+4,ngb*nstate,ntqxx*2)
          deallocate(zzmel)

          if(debug) then
            do  it  = 1,nctot+nbmax
            do  itp = 1,ntp0
              write(6,"(' w3p =',2i4,2d14.6)") it,itp,w3p(it,itp)
            enddo
            enddo
          endif

c Write the Spectrum function for exchange May. 2001
          if(ifexsp/=0) then
            do it  = 1, nctot+nbmax
              do itp = 1,ntp0
                write(ifexsp,"(3i4, 3f12.4, ' ',d23.15,'  ',d23.15)")
     &      ip,itp,it, qbz_kr, ekc(it), -wtt*w3p(it,itp)
              enddo
            enddo
          endif

!! --- Correct weights wfac for valence by esmr
          do it = nctot+1, nctot+nbmax
cgaussian
            wfac = wfacx(-1d99, ef, ekc(it), esmr)
            if(debug) print *, ' it wfac=',it, wfac,ef, ekc(it), esmr
            w3p(it,1:ntp0) = wfac * w3p(it,1:ntp0)
          enddo

          if (.not.tote) then !total energy mode tote
            do itp = 1,ntp0  !S[j=1,nbloch]  z1p(j,t,t') <B(rk,j) psi(q-rk,n) |psi(q,t')>
              zsec(iwini,itp,ip) = zsec(iwini,itp,ip)
     &        - wtt * sum( w3p(:,itp) )
            enddo
          else
            do itp = 1,ntp0
cgaussian
              wfac = wfacx(-1d99, ef2, eq(itq(itp),ip), esmr2) !june-2002 takao
              w3p(1:nctot+nbmax,itp) = wfac * w3p(1:nctot+nbmax,itp)
              exxq = exxq - wtt * sum( w3p(:,itp) )
            enddo
          endif
C         if(.not.newaniso()) deallocate(vcoul)
C         deallocate( w3p,rmelt,cmelt)
          deallocate(w3p)
          cycle
        endif
c-- End of exchange section --------------


c--------------------------------------------------------------------------
c--- correlation section --------------------------------------------------
c--------------------------------------------------------------------------
c--- The matrix elements zmel.
        allocate( zmel (ngb, ntp0, nstate) )
C Original
C        zmel = dcmplx (rmelt,-cmelt)
Cc        print *,'ddddddddddddd original like mode ddddddddddddddddd'
C        if(newaniso()) then
C          do itp=1,ntp0
C            do it=1,nstate
Ccooo
C              zmel(:,it,itp) =  matmul(zmel(:,it,itp),dconjg(ppovlz(:,:)))
C            enddo
C          enddo
C          if(verbose()>=30) call cputid2(' complete zmel * overlap+',0)
C        endif
C        call zprm('zmel*ppovlz+, direct',zmel,ngb,ngb,nstate*ntp0)
C Vectorized
C        if(newaniso()) then
C         allocate(zmelt(ngb,nstate,ntp0))
C         zmelt = dcmplx(rmelt,-cmelt)
C         call zprm('zmelt',zmelt,ngb,ngb,nstate*ntp0)

!         Use CC of zzmel as calculated in sxcf_mvs
          zzmel = dconjg(zzmel)

          call zgemm('C','N',ngb,nstate*ntp0,ngb,(1d0,0d0),
     .      ppovlz,ngb,zzmel,ngb,(0d0,0d0),zmel,ngb)
C         deallocate(zmelt)
          if(verbose()>=30) call cputid2(' complete zmel * overlap+ .',0)
C         call zprm('zmel*ppovlz+, vectorized',zmel,ngb,ngb,nstate*ntp0)
C        else
CC         zmel = dcmplx(rmelt,-cmelt)
C          zmel = dconjg(zzmel)
C        endif
        ia = idalloc('zzmel',iverbose+4,ngb*nstate,ntqxx*2)
        deallocate(zzmel)
        if(debug) print *,' sum abs zmel=',sum(abs(zmel(:,:,:)))

        if(debug) print *,' end of zmel'
c================================================================
c The correlated part of the self-energy:
c S[n=all] S[i,j=1,nbloch]
c <psi(q,t) |psi(q-rk,n) B(rk,i)>
c  < [w'=0,inf] (1/pi) (w-e)/{(w-e)^2 + w'^2} Wc(k,iw')(i,j) >
c                                <B(rk,j) psi(q-rk,n) |psi(q,t)>
c e = e(q-rk,n), w' is real, Wc = W-v
c================================================================

        allocate( zw (nblochpmx,nblochpmx) )
c        if(imgonly) goto 3333 !jan2008


c====================================================================
c contribution to SEc(qt,w) from integration along the imaginary axis
c====================================================================

c------------------------------------------------
c loop over w' = (1-x)/x, frequencies in Wc(k,w')
c {x} are gaussian points between (0,1)
c------------------------------------------------
        allocate( zwz0(ntp0,nstate) )
        ix = 1  - nw_i !at omega=0

c        nrec=(kx-iqini)*(nw-nw_i+1) +ix ! 2---> iqini
c        if(bzcase()==2) nrec= (kx-1)*(nw-nw_i+1) +ix
        nrec=ix

        if(debug) print *,' wvr nrec kx nw nw_i ix=',nrec,kx,nw,nw_i,ix
        read(ifrcw,rec=nrec) zw  ! direct access read Wc(0) = W(0) - v
ccccccccccccccccccccccccccccccccccccccc
c        print *,'kkkkkk sumcheck zw=',kx,sum(abs(zw(1:ngb,1:ngb)))
ccccccccccccccccccccccccccccccccccccc
c1
c1        call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,
c1     o    zwz0) ! zwz0 = zmel*(W(0)-v)*zmel
c2
C Se

C       In-line
C        zwz0=0d0
Cc        allocate(zmel1(ngb))
C!! this loop looks complicated but just in order to get zwz0=zmel*zwz0*zmel
C!! Is this really efficient???
C!$OMP parallel do private(itp,it,igb2,zz2)
C        do itp=1,ntp0
C          do it=1,nstate
Cc          zmel1(:)=dconjg(zmel(:,it,itp))
C            do igb2=2,ngb
Cc            zz2 = sum( zmel1(1:igb2-1)*zw(1:igb2-1,igb2) )
C              zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )
C              zwz0(it,itp) = zwz0(it,itp)+zz2*zmel(igb2,it,itp)*2d0+
C     &       dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp)
C            enddo !igb2
C            zwz0(it,itp) = zwz0(it,itp)+
C     &    dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
C          enddo !it
C        enddo !itp
C!$OMP end parallel do
C        zwz0 = dreal(zwz0)
C        if(verbose()>=30) call cputid2(' complete zmel+ * (Wc-v) * zmel ',0)
C       call zprm('inline zwz0',zwz0,nstate,nstate,ntp0)
c2end

C       Vectorized zwz_t,t' = sum_i zz2_i_t,t'*zm_i,t,t'
C       where zz2_i,t,t' = sum_i' zm*_i' zw_i',i
C       call zprm('zmel',zmel,ngb,ngb,ntp0*nstate)
C       call zprm('zw',zw,nblochpmx,ngb,ngb)
C       ia = idalloc('zmelt',iverbose+2,ngb*nstate,ntp0*2)
        allocate(zmelt(ngb,ntp0,nstate))
        call zgemm('N','N',ngb,nstate*ntp0,ngb,(1d0,0d0),
     .    zw,nblochpmx,zmel,ngb,(0d0,0d0),zmelt,ngb)
C       call zprm('W * z',zmelt,ngb,ngb,nstate*ntp0)
C        do itp = 1, ntp0
C          do  it = 1, nstate
C            zwz0(it,itp) = zdotc(ngb,zmel(1,it,itp),1,zmelt(1,it,itp),1)
C          enddo
C        enddo
        do  it = 1, nstate
          do itp = 1, ntp0
            xx = 0
            do  igb2 = 1, ngb
              xx = xx + dble(zmel(igb2,itp,it))*dble(zmelt(igb2,itp,it))
     .                + dimag(zmel(igb2,itp,it))*dimag(zmelt(igb2,itp,it))
            enddo
            zwz0(itp,it) = xx
          enddo
        enddo
C       ia = idalloc('zmelt',iverbose+4,ngb*nstate,ntp0*2)
        deallocate(zmelt)

        if(verbose()>=30) call cputid2(' complete z * W(0) * z. ',0)
C       call zprm('vectorized zwz0',zwz0,nstate,nstate,ntp0)


cccccccccccccccccccccccc old variant ccccccccccccccccccccccccccccccccc
c        read(ifrcw,rec=((kx-iqini)*nw+ix)) zw  ! direct access read Wc(0) = W(0) - v
c        call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,
c     o    zwz0) ! zwz0 = zmel*(W(0)-v)*zmel
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c COH term test ----- The sum of the all states for zwz00 gives the delta function.
        if(cohtest) then
          do itp = 1,ntq
            coh(itp,ip)  = coh(itp,ip)
     &       + .5d0*wtt*sum(dreal(zwz0(itp,1:nstate)))
          enddo
          deallocate(zw,zwz0,zmel)
          cycle
        endif
c
        nx  = niw
        if(niw <1) call rx( " sxcf:niw <1")

c May2006 zwz
        npm = 1
        if(nw_i/=0) npm = 2  ! TimeReversal off case.

        if(allocated(zwz)) deallocate(zwz)
        if(allocated(zwzi)) deallocate(zwzi)
        allocate( zwz(niw*npm,ntp0, nstate),  zwzi(nstate,ntp0) )
        if(screen) allocate( zwz00(ntp0,nstate) )
c
        if(verbose()>50) write(*,'("6 before matzwz in ix cycle ",$)')
        if(verbose()>50) call cputid(0)

**************************** BEGIN 1 ************************************************
c         do ix=1,nx
c            nrec= (kx-iqini)*niw*npm + ix
c            read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis
c            nrec= (kx-iqini)*niw*npm + ix + nx
c            read(ifrcwi,rec=nrec) zw2  ! Readin W-v on imag axis
c           write(6,"(' diff =',2i4,d13.6)") nx,ix, sum(abs(zw-zw2))
c         enddo

        zwz=0d0            !...........................faleev, temp
        ! allocate(zmel1(ngb))
        ia = idalloc('zmelt',iverbose+2,ngb*nstate,ntp0*2)
        allocate(zmelt(ngb,ntp0,nstate))
        do ix = 1,nx  !*npm            ! imaginary frequency w'-loop

           nrec= ix
c          nrec= (kx-iqini)*niw  + ix
c          print *,' ix nrec=',ix,niw,nrec
c          if(bzcase()==2) nrec= (kx-1)*niw + ix

          if(debug) print *,' wvi nrec=',nrec
          read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis
ccccccccccccccccccccccccccccccccccccccc
c          print *,'kkkkkk sumcheck zwi=',kx,sum(abs(zw(1:ngb,1:ngb)))
ccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
c          write(6,"(' trans =',i4,2d13.6)") ix
c     &    ,sum(abs(zw-dconjg(transpose(zw)))),sum(abs(zw))
ccccccccccccccccccccccccccccccccc

C In-line
C          if(npm==1) then !then zwz is real so, we can use mode c2.
Cc$$$c1          call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,
Cc$$$c1     o      zwz(ix,1:nstate,1:ntp0))
C            do itp= 1,ntp0
C              do it = 1,nstate
C                ppp=0d0
C                do igb2 = 2,ngb
Cc               zz2 = sum( zmel1(1:igb2-1)*zw(1:igb2-1,igb2) )
C                  zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )
C! only take real part
C                  ppp = ppp + dreal(zz2*zmel(igb2,it,itp)) * 2
C     &           + dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp)
C                enddo !igb2
C                zwz(ix,it,itp) = ppp +
C     &       dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
C              enddo !it
C            enddo!itp
C!$OMP end parallel do
Cc2end
C          else !we need to use mode2 because zwz is not real now.
C            call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,
C     o          zwz(ix,1:nstate,1:ntp0))
C          endif
C          if(verbose()>=50) call cputid2(' complete zmel+ * (Wc-v) * zmel ',0)
C          call zprm('inline zwz',zwz(ix,1:nstate,1:ntp0),nstate,nstate,ntp0)


          call zgemm('N','N',ngb,nstate*ntp0,ngb,(1d0,0d0),
     .      zw,nblochpmx,zmel,ngb,(0d0,0d0),zmelt,ngb)
C         call zprm('W * z',zmelt,ngb,ngb,nstate*ntp0)
          do  it = 1, nstate
            do itp = 1, ntp0
                zwz(ix,itp,it) = zdotc(ngb,zmel(1,itp,it),1,zmelt(1,itp,it),1)
            end do
          enddo
          if(verbose()>=50) call cputid2(' complete zmel+ * (Wc(om)-v) * zmel ',0)
C         call zprm('vectorized zwz',zwz(ix,1:nstate,1:ntp0),nstate,nstate,ntp0)

        enddo !ix
        ia = idalloc('zmelt',iverbose+4,ngb*nstate,ntp0*2)
        deallocate(zmelt)
        if(verbose()>=30) call cputid2(' complete z * W(omega) * z. ',0)
        if(debug) print *,' sumzmel=',ngb, nstate, ntp0,sum(abs(zmel))
******************************* END 1 ***********************************************
!  zwz      symmetric part zwz(1:nx,        :,:)
!      anti-symmetric part zwz(nx+1:nx*npm, :,:)
c         do ix=1,nx
c           write(6,"(' diff =',i4,d13.6)") ix, sum(abs((zwz(ix+nx,:,:)-zwz(ix,:,:))))
c         enddo


c********************************* BEGIN 3 ***************************************
c       if(.false.) then
c        do ix = 1,nx     ! imaginary frequency w'-loop         !********* faleev
c         nrec=(kx-iqini)*niw+ix
c         if(bzcase()==2) nrec= (kx-1)*niw+ix
c         if(debug) print *,' wvi2 nrec=',nrec
c         read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis
c
cc zwz= S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
cc                Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
cc        do itp = 1,ntp0
cc        do  it = 1,nstate
cc          zwz(ix,it,itp) = sum(
cc     &   dconjg(zmel(:,it,itp)),matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)) )
cc        enddo
cc        enddo
c          call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb, !********* faleev
c     o      zwz(ix,1:nstate,1:ntp0))
c        enddo !ix
c        endif                        !********* faleev
cccccccccccccccccccccccccccccccccccccccc
cc zwz is the diagonal term and real.
cc      write(6,*) ' img sumzwz=',sum(abs(imag(zwz(ix,1:nstate,1:ntp0))))
cccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc
cc      write(6,*) ' sumzwz=',sum(zmel)
cc     & ,sum(zw(1:ngb,1:ngb)),sum(zwz(ix,1:nstate,1:ntp0))
cccccccccccccccccccccccccccccccccc
c************************************ END 3 *************************************
c
c        if(verbose()>50) write(*,'("7 after matzwz in ix cycle ",$)')
c        if(verbose()>50) call cputid(0)

cccccccccccccccccccccccccccccccccccccccccccccc
c      do  it  = 1,nstate
c      itp=1
c      ix =1
c     do  itp = 1,ntp0
c      do  ix = 1,nx
c        write(6,"(' zwz on img =',3i4,2d14.6)")ix,it,itp,zwz(ix,it,itp)
c      enddo
c      enddo
c     enddo
c      stop " *** TEST END *** "
cccccccccccccccccccccccccccccccccccccccccccccc


c--------------------------------------------------------------
c S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
c needed to take care of the singularity in the w' integration
c when w-e(q-rk,n) is small
c--------------------------------------------------------------
        if(screen) then
          zwz00 = zwz0
          zwz0  = 0d0
          do ix = 1,nx
            zwz(ix,:,:)=zwz(ix,:,:) - zwz00
          enddo
        endif

ccccccccccccccccccccccccccc
c       if(ifzwz==-999) ifzwz = iopen('zwz.hsfp0',1,-1,0)
c        do itp = 1,ntq
c        do it  = 1,nstate
c           do iw=1,niw
c             write(ifzwz,"(4i5,2d24.16)") iw,it,itp,kx,zwz(iw,it,itp)
c           enddo
c           iw = 0
c           write  (ifzwz,"(4i5,2d24.16)") iw,it,itp,kx,zwz0(it,itp)
c           write(ifzwz,*)
c        enddo
c        enddo
ccccccccccccccccccccccccccc

c------------------------------------------------
c loop over w in SEc(qt,w)
c------------------------------------------------
ccccccccccccccccccccccccccccccccccccc
ctest imx
c        testimx=.true.
c        if(testimx) zwzi=0d0
ccccccccccccccccccccccccccccccccccccc

        if(ua_auto) then
          allocate(uaa(nstate,ntq))
          do itp = 1,ntq
            do  it = 1,nstate
              ratio = abs(zwz(niw,itp,it)/zwz0(itp,it))
              call gen_uaa(ratio,freqx(niw),  uaa(it,itp))
              if(verbose()>45) then
                write(6,"(' it itp uaa=',2i4,12f8.4)") it,itp,uaa(it,itp)
              elseif(verbose()>40.and.mod(it,10)==1.and.mod(itp,10)==1) then
                write(6,"(' it itp uaa=', 2i4,12f8.4)") it,itp,uaa(it,itp)
              endif
            enddo
          enddo
        endif

        allocate(zwzs(npm*nx))
        do      iw = iwini,iwend
c frequency integration along the imaginary axis, s. wint.f
c for each e(q-rk,n) and w in SEc(qt,w)
          do 1385  itp = 1,ntq
            do 1387 it = 1,nstate
              we =.5d0*( omega(itp,iw) -ekc(it)) != .5d0*( eq(itq(itp),ip)+2d0*(dble(iw)-shtw)*deltaw-ekc(it))
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ctest imx
c         if(testimx) then
c           if(ef-omega(itp,iw)>0) then
c             if(  +0.001<-we .and. ekc(it)< ef) cycle
c           endif
c           if(ef-omega(itp,iw)<=0) then
c             if(  -0.001>-we .and. ekc(it)> ef) cycle
c           endif
c         endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccc
              if(verbose()>50) then
                do  ix = 1,niw
                  ratio  = abs(zwz(ix,itp,it)/zwz0(itp,it))
                  freqw1 = (1d0 - freqx(ix))/ freqx(ix)
                  ua2_(ix) = sqrt(- 1d0/freqw1*log(ratio))
c            write(6,"(' ix freqw1=',i4,f10.4)") ix,freqw1
                enddo
                write(6,"(' sxcf_fal2: ua=sqrt(1/w1*log(v0/v1))=',
     &         12f8.4)") ua2_(1:niw)
              endif

c          if(ua_auto) then
c            call gen_ua(abs(zwz(niw,it,itp)/zwz0(it,itp)), niw,freqx, expa_,ua_)
c            if(iw==ini) then
c            if(verbose()>45) then
c              write(6,"(' it itp ua_=',2i4,12f8.4)")it,itp,ua_
c            elseif(verbose()>40.and.mod(it,20)==1.and.mod(itp,20)==1) then
c              write(6,"(' it itp ua_=',3i4,12f8.4)")it,itp,ua_
c            elseif(irot==1.and.mod(it,10)==1.and.itp==it) then
c              write(6,"(' it itp ua_=',3i4,12f8.4)")it,itp,ua_
c            endif
c            endif
c          endif
              if(ua_auto) then
                ua_ = .5d0*uaa(it,itp)
                call gen_expa(niw,freqx,ua_,  expa_)
              endif

c----------------------
              esmrx = esmr
              if(it <= nctot) esmrx = 0d0
c
              do ix=1,nx
                zwzs(ix   ) = dreal( zwz(ix,itp,it)) ! w(iw) + w(-iw) symmetric part
                if(npm==2) then
                  zwzs(ix+nx) = dimag( zwz(ix,itp,it)) ! w(iw) - w(-iw)
                endif
              enddo

              if(GaussSmear()) then
                zwzi(it,itp) = wintzsg_npm(npm, zwzs,zwz0(itp,it),freqx,wx,ua_,expa_,we,nx, esmrx)
              else
                if(npm==2)
     &            call rx( ' ###Not implement wintzav for npm=2. Use Gausssmear.')
                zwzi(it,itp) = wintzav( zwzs,zwz0(itp,it),freqx,wx,ua_,expa_,we,nx, esmrx)
              endif
c    .    wintz (zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx)
ccccccccccccccccccccccccccccccc
c          if(verbose()>45) then
c          if(it==50.and.itp==1) then
c          write(6,"(' it itp abs(zwzi)=',2i4,12d13.5)")it,itp,abs( zwzi(it,itp))
c          icc=icc+1
c          if(icc==10) stop 'test end'
c          endif
c          endif
cccccccccccccccccccccccccccccc
 1387       continue
 1385     continue

c sum over both occupied and unoccupied states and multiply by weight
          do     itp = 1,ntq
            zsec(iw,itp,ip)  = zsec(iw,itp,ip) + wtt*sum(zwzi(:,itp))
          enddo

c end of SEc w-loop
        enddo
        deallocate(zwzs)
        if(debug) then
          print *,' ntq nstate sum(zwzi)=',ntq,nstate,sum(zwzi)
          print *,' ntq nstate sum(zwz )=',ntq,nstate,sum(zwz)
          do itp = 1,ntq
            write(6,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
          enddo
        endif
        deallocate(zwz,zwz0,zwzi)

        if(ua_auto) deallocate(uaa)
ccccccccccccccccccccccccccccccccccccc
c        if(onlyimagaxis()) then
c          print *,' sxcf_fal2: onlyimagaxis cycle'
c          deallocate(freq_r,zmel,zw,zmel1)
c         goto 1100
c        endif
c        if(testimx) then
c          print *,' sxcf_fal2: testimx cycle'
c          deallocate(freq_r,zmel,zw,zmel1)
c         goto 1100
c        endif
ccccccccccccccccccccccccccccccccccccc
c        do iw=nw_i,nw
c        write(6,*) iw,freq_r(iw)
c        enddo
c        stop 'cccccccccccccccccccccccccccccccccccccccccccc'



c 3333   continue
c===============================================================================
c contribution to SEc(qt,w) from the poles of G
c===============================================================================
ccccccccccccccccccccccc
c        debug=.true.
c        if(debug) print *,' --------- go to poles'
cccccccccccccccccccccccc

C... Check freq_r for nw_i case.
!    We assume freq_r(i) == -freq_r(-i) in this code. feb2006
        if(nw_i/=0) then
          if(nw/= -nw_i)        call rx( "sxcf_fal2z: nw/=-nw_i")
          if(freq_r(0)/=0d0)    call rx( "sxcf_fal2z: freq_r(0)/=0")
          if( sum(abs( freq_r(1:nw)+freq_r(-1:-nw:-1)))/=0)
     &       call rx( "sxcf_fal2z: freq_r /= -freq_r")
        endif

c---------------------------------------
c maximum ixs finder
c---------------------------------------
c      print *,' ekc at nt0p nt0m+1=', ekc(nt0p),ekc(nt0m+1)
c      print *,'  nt0p nt0m+1=', nt0p, nt0m+1
        ixsmx =0
        ixsmin=0
        do 3001 iw  = iwini,iwend
          do 3002 itp = 1,ntq
            omg = omega(itp,iw)
c          write(6,"('3001: iw itp omg-ef=',2i5,d13.6)")iw,itp,omg-ef
            if (omg < ef) then
              itini= 1
              itend= nt0p
            else
              itini= nt0m+1
              itend= nstate
            endif
            do 3011 it= itini,itend
              esmrx = esmr
              if(it<=nctot) esmrx = 0d0
              wfac = wfacx2(omg,ef, ekc(it),esmrx)
              if(GaussSmear()) then
                if(wfac<wfaccut) cycle
                we = .5d0*(omg-weavx2(omg,ef,ekc(it),esmr))
              else
                if(wfac==0d0) cycle
                if(omg>=ef) we = max( .5d0*(omg-ekc(it)), 0d0) ! positive
                if(omg< ef) we = min( .5d0*(omg-ekc(it)), 0d0) ! negative
              endif
              do iwp  = 1,nw      ! may2006
                ixs = iwp         ! ixs = iwp= iw+1
c                write (*,*) 'xxx freq we=',freq_r(iwp),abs(we)
                if(freq_r(iwp) > abs(we)) exit
              enddo
c This change is because G(omega-omg') W(omg') !may2006
c             if(ixs>ixsmx  .and. omg<=ef ) ixsmx  = ixs
c             if(ixs>ixsmin .and. omg> ef ) ixsmin = ixs
              if(ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
              if(ixs>ixsmin .and. omg< ef ) ixsmin = ixs
              wexx  = we
              if(ixs+1 > nw) then
                write (*,*) ' nw_i ixsmin',nw_i, ixsmin
                write (*,*) ' wexx, dw ',wexx,dw
                write (*,*) ' omg ekc(it) ef ', omg,ekc(it),ef
                call rx( ' sxcf 222: |w-e| out of range')
              endif
 3011       continue
 3002     continue  !end of SEc w and qt -loop
 3001   continue  !end of SEc w and qt -loop
        if(nw_i==0) then
          nwxi = 0
          nwx  = max(ixsmx+1,ixsmin+1)
        else
          nwxi = -ixsmin-1
          nwx  =  ixsmx+1
        endif
        if (nwx > nw   ) then
          call rx( ' sxcf nwx check : |w-e| > max(w)')
        endif
        if (nwxi < nw_i) then
          call rx( ' sxcf nwxi check: |w-e| > max(w)')
        endif
        if(debug) print *,' nwxi nwx nw=',nwxi,nwx,nw

C... Find nt_max ------------------------------------
        nt_max=nt0p !initial nt_max
        do 4001 iw  = iwini,iwend
          do 4002 itp = 1,ntq
            omg     = omega(itp,iw)
            if (omg > ef) then
              do  it = nt0m+1,nstate    ! nt0m corresponds to efm
                wfac = wfacx2 (ef,omg, ekc(it),esmr)
                if( (GaussSmear().and.wfac>wfaccut)
     &         .or.(.not.GaussSmear().and.wfac/=0d0)) then
                  if (it > nt_max) nt_max=it ! nt_max is  unocc. state
                endif                        ! that ekc(it>nt_max)-omega > 0
              enddo
            endif
4002      continue
4001    continue

ccccccccccccccccccccccc
c        debug=.false.
cccccccccccccccccccccccc

C... Set zw3 or zwz -----------------------------------
        zwz3mode=.true.
        if(iwend-iwini>2) then
          zwz3mode=.false.
        endif
ccccccccccccccccccccccccccccccccccc
c        write(6,*)'zwz3mode=',zwz3mode
cccccccccccccccccccccccccccccccccccccc

        if(zwz3mode) then
          allocate( zw3(ngb,ngb,nwxi:nwx))
          do ix = nwxi,nwx                   ! real frequency w'-loop
!           write(*,*)'r_om', ix,nwx,nw_w
c            nrec=(kx-iqini)*(nw+1-nw_i)+ ix-nw_i+1
c            if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
            nrec=ix-nw_i+1

            if(debug) print *,' wvr3 nrec=',nrec,nblochpmx,kx,ix,nw
            read(ifrcw,rec=nrec) zw
            zw3(1:ngb,1:ngb,ix) = zw(1:ngb,1:ngb)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            itp = 1
c            it = 1
c            zwz3x=0d0
c            do igb2=1,ngb
c              zwz3x = zwz3x
c     &                + sum(dconjg(zmel(1:ngb,it,itp))*zw3(1:ngb,igb2,ix))
c     &                  * zmel(igb2,it,itp)
c            enddo
c            write(6,"(' zwz3xtest: ix zwz3x=',i4,2d13.5)") ix, zwz3x
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(evaltest()) then
              write(6,"('iii --- EigenValues for zw --------')")
              allocate(ebb(ngb))
              call diagcvh2((zw(1:ngb,1:ngb)-transpose(dconjg(zw(1:ngb,1:ngb))))/2d0/img,
     &        ngb, ebb)
              do ii=1,ngb
                if(abs(ebb(ii))>1d-8.and.ebb(ii)>0) then
                  write(6,"('iii1xxx:  iw ii eb=',2i4,d13.5)") ix,ii,ebb(ii)
                else
                  write(6,"('iii1:  iw ii eb=',2i4,d13.5)") ix,ii,ebb(ii)
                endif
              enddo
              deallocate(ebb)
            endif
          enddo
          deallocate(zw)
        else
          nstatex= max(ntp0,nt_max)
          if(allocated(zwz)) deallocate(zwz)
          allocate( zwz(nwxi:nwx,ntp0,1:nstatex) )
          do      ix = nwxi,nwx
            nrec= ix-nw_i+1
c            nrec=(kx-iqini)*(nw+1-nw_i)+ ix-nw_i+1
c            if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1

            read(ifrcw,rec=nrec) zw  ! Readin (W-v)(k,w')(i,j) at k and w' on imag axis
c zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
            call matzwz(zw(1:ngb,1:ngb), zmel(1:ngb,1:ntp0,1:nstatex), ntp0,nstatex,ngb,
     o      zwz(ix,1:ntp0,1:nstatex))
! zmel (ngb, nstate, ntp0)
          enddo
          deallocate(zmel)
          deallocate(zw)
        endif

c---------------------------------------------
        if(screen) then
          if(zwz3mode) call rx( ' screen mode is not implemented')
          do ix = nw_i,nwx
            zwz(ix,:,:)=zwz(ix,:,:) - zwz00
          enddo
          deallocate(zwz00)
        endif

c-------------------------------
c loop over w and t in SEc(qt,w)
c-------------------------------
C        if(verbose()>50) write(*,'("10 wfacx  iw,itp,it cycles ",$)')
C        if(verbose()>50) call cputid(0)
        if(verbose()>50) call cputid2(' start loop over omega, band pairs',0)
        do 2001 iw  = iwini,iwend
        do 2001 itp = 1,ntq
c          if(debug) print *,'2011 0 zmel=',sum(abs(zmel(:,:,:)))
          omg = omega(itp,iw)
c          write(6,"('2001: iw itp omg=',2i5,d13.6)")iw,itp,omg-ef
          if (omg >= ef) then
            itini= nt0m+1
            itend= nt_max
            iii=  1
          else
            itini= 1
            itend= nt0p
            iii= -1
          endif

          do 2011 it= itini,itend
c           if(debug) print *,'2011 1 loop--- it=',iw,itp,it,sum(abs(zmel(:,:,:)))
            esmrx = esmr
            if(it<=nctot) esmrx = 0d0
            wfac = wfacx2(omg,ef, ekc(it),esmrx)
            if(GaussSmear()) then
              if(wfac<wfaccut) cycle
              we = .5d0*abs(omg-weavx2(omg,ef, ekc(it),esmr))
            else
              if(wfac==0d0) cycle
              if(omg>=ef) we = 0.5d0* abs(max(omg-ekc(it), 0d0)) ! positive
              if(omg< ef) we = 0.5d0* abs(min(omg-ekc(it), 0d0)) ! negative
            endif

            wfac= iii* wfac*wtt


cccccccccccccccccccccccccccccccc
c            if(debug)
c            if(iw==0)  write(6,"( ' xxx1',i3, 10d11.3)") iw,omg,ef, ekc(it),wfac
cccccccccccccccccccccccccccccccc



c three-point interpolation for Wc(we)
            do iwp = 1,nw
              ixs=iwp
              if(freq_r(iwp)>we) exit
            enddo
            if(nw_i==0) then
              if(ixs+1>nwx) then
                 print *,' ixs,nwx, we =',ixs,nwx,we
                 call rx( ' sxcf: ixs+1>nwx xxx2')
              endif
            else !   write(6,*)" ixs nwxi=",ixs,nwxi,freq_r(ixs-1),we,freq_r(ixs)
              if(omg >=ef .and. ixs+1> nwx ) then
                write(6,*)'ixs+1 nwx=',ixs+1,nwx
                call rx( ' sxcf: ixs+1>nwx yyy2a')
              endif
              if(omg < ef .and. abs(ixs+1)> abs(nwxi) ) then
                write(6,*)'ixs+1 nwxi=',ixs+1,nwxi
                call rx( ' sxcf: ixs-1<nwi yyy2b')
              endif
            endif

            iir=1
            if(omg < ef .and. nw_i/=0) iir = -1  !May2006 because of \int d omega' G(omega-omega') W(omega')

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(debug) then
              fffr = freq_r( iir*(ixs-1):iir*(ixs+1):iir )
              write(6,*)
              write(6,*) ' itp it iw=',itp,it,iw
              write(6,"(' uuu: itp it iw omega=',3i4,14d13.6)")
     &        itp,it,iw,omega(itp,iw),ekc(it)
              write(6,"(' uuu: weavx omg ef ekc esmr=',14d13.6)")
     &       weavx2(omg,ef, ekc(it),esmr), omg, ef, ekc(it), esmr
              write(6,"(' uuu: fffr=',d13.6,2x,3d13.6)") we,fffr
              print *,'uuu zwz3mode=',zwz3mode
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(zwz3mode) then
              zwz3=(0d0,0d0)
CCC!$OMP parallel do private( ix,igb2,zz2)
              if(debug) write(6,"('wwwwwww ixs=',10i4)") ixs,igb2,it,itp
              if(debug) print *,'2011 www zmel aaa=',sum(abs(zmel(:,:,:)))

              do ix = ixs, ixs+2
                do igb2=1,ngb
c                  print *,' ix igb2=',ixs,ix,igb2,iir
c                  print *,' it itp=',it,itp
c                  print *,' zw3=', sum(zw3(1:ngb,igb2,iir*(ix-1)) )
c                  print *,' zmel=',sum(dconjg(zmel(1:2,it,itp)))
                  zz2 = sum(dconjg(zmel(1:ngb,itp,it))*zw3(1:ngb,igb2,iir*(ix-1)) )
                  zwz3(ix-ixs+1) = zwz3(ix-ixs+1)+zz2 *zmel(igb2,itp,it)
                enddo !igb2
              enddo !ix
CCC!$OMP end parallel do
              if(debug) write(6,"('w xxxxxxxxxxxxx ixs loopend=',i4)") ixs
              if(debug) print *,zwz3(1:3)          !,freq_r(ixs-1),zwz3(1:3)
              if(debug) print *,'we frez zwz3=', we,ixs,freq_r(ixs-1:ixs+1)
              if(debug) print *,'2011 bbb www zmel=',sum(abs(zmel(:,:,:)))

              zsec(iw,itp,ip) = zsec(iw,itp,ip)
     &      + wfac *alagr3zz(we,freq_r(ixs-1),zwz3) !faleev

              if(debug) print *,'2011 ccc www zmel=',sum(abs(zmel(:,:,:)))
              if(debug) write(6,"('wwwwwww eo zsecsum')")

cccccccccccccccccccccccccccccc
c            if(iw==0) then
c              write(6,"('zzz1:')")
c              qqqq= dimag(alagr3zz(we,freq_r(ixs-1),zwz3))
c              if(omg>=ef) wex = 0.5d0* abs(max(omg-ekc(it), 0d0)) ! positive
c              if(omg< ef) wex = 0.5d0* abs(min(omg-ekc(it), 0d0)) ! negative
c               write(6,"('zzz1: ',2d13.5,2x,4d13.5)")
c     &          we,wex, freq_r(ixs-1:ixs+1)
c               write(6,"('zzz1: ',26x,2x,4d13.5)")
c     &           dimag(zwz3(1:3)), qqqq
c              write(6,*)
c            endif
cccccccccccccccccccccccccccccc



c           else
c            zsec(iw,itp,ip) = zsec(iw,itp,ip)
c     &      + wfac*alagr3(we,freq_r(ixs-1),zwz3) !faleev
c           endif
            else
              zwzz(1:3) = zwz(iir*(ixs-1):iir*(ixs+1):iir,itp, it)
c           if(npm==1) then
              zsec(iw,itp,ip) = zsec(iw,itp,ip)
     &        + wfac*alagr3zz(we,freq_r(ixs-1),zwzz)
c            else
cccccccccccccccccccccccccccccc
c            write(6,"('zzz: npm ',i3,6d13.5)") npm,zwzz(1:3)
c            write(6,"('zzz: ',6d13.5)") dimag(zwzz(1:3))
c            if(abs(iw)==30) then
c            write(6,"('zzz: omega=',2i4,6d13.5)") itp,iw,omega(itp,iw)
c            write(6,"('zzz: ',3i4,6d13.5)") iw,itp,ip,wfac*alagr3zz(we,freq_r(ixs-1),zwzz)
c            endif
cccccccccccccccccccccccccccccc
c            zsec(iw,itp,ip) = zsec(iw,itp,ip)
c     &      + wfac*alagr3(we,freq_r(ixs-1),zwzz)
c            endif
            endif
c------------
cccccccccccccccccccccccccccccccccccccccc
            if(debug) write(6,"('wwwwwww bottom of 2011 loop')")
cccccccccccccccccccccccccccccccccccccccc
 2011     continue
cccccccccccccccccccccccccccccccccccccccc
c         write(6,"('www1: ',4d13.5)") zsec(30,1,1)
cccccccccccccccccccccccccccccccccccccccc
 2001   continue  !end of SEc w and qt -loop
cccccccccccccccccccccccccccccccccccccccc
c         write(6,"('www222: ',4d13.5)") zsec(30,1,1),zsec(-30,1,1)
cccccccccccccccccccccccccccccccccccccccc
        if(debug) print *,' end of do 2001'
        if(verbose()>50) call cputid2(' end loop over omega, band pairs',0)
C        if(verbose()>50) then
C          write(*,'("11 after alagr3zz iw,itp,it cycles ",$)')
C          call cputid(0)
C        endif
        if(debug) then
          do itp = 1,ntq
            write(6,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
          enddo
        endif
        if (allocated(zmel1))  deallocate(zmel1) !,zwz)
        if(zwz3mode) then
          deallocate(zmel,zw3)
        else
          deallocate(zwz)
        endif
 1000 continue
!     if(newaniso()) ifvcoud =iclose('Vcoud.'//charnum5(kx))
      ifvcoud =iclose('Vcoud.'//charnum5(kx))
      if(.not.exchange) then
      ifrcw  = iclose('WVR.'//charnum5(kx))
      ifrcwi = iclose('WVI.'//charnum5(kx))
      endif
 1100 continue  ! end of k-loop
      if(verbose()>=40) call cputid2(' completed 1100 k-point loop',0)

c     if(smbasis()) isx = iclose("POmat")
c      if(irot==1) write(6,"('  sum(abs(zsec))=',d23.15)")sum(abs(zsec))

c      call rlse (iwr1)
ccccccccccccccccccccccccccccc
c      print *," end kount ip irot=", irot, ip, sum(kount(:,ip))
c       print *," sumkount 1=",irot,ip,sum(kount(:,:))
cccccccccccccccccccccccccccc
      if (tote) then
c         if (tetraex) then
c            exx = exx + nstar(ip) * exxq * 0.5d0
c         else
        exx = exx + wik(ip) * exxq * 0.25d0
c         endif
      endif
c      call rlse (iwr1)
ccccccccccccccccccccccccccccc
c      print *," end kount ip irot=", irot, ip, sum(kount(:,ip))
c       print *," sumkount 1=",irot,ip,sum(kount(:,:))
cccccccccccccccccccccccccccc
C      if (allocated(rmelt)) deallocate(rmelt)
C      if (allocated(cmelt)) deallocate(cmelt)
      if (allocated(zz)) deallocate(zz)
      if (allocated(zmel)) deallocate(zmel)
      if (allocated(zzmel))deallocate(zzmel)
      if (allocated(zw)) deallocate(zw)
      if (allocated(zwz)) deallocate(zwz)
      if (allocated(zwz0)) deallocate(zwz0)
      if (allocated(zwzi)) deallocate(zwzi)
      if (allocated(zwz00)) deallocate(zwz00)
      if (allocated(w1p)) deallocate(w1p)
      if (allocated(w2p)) deallocate(w2p)
      if (allocated(w3p)) deallocate(w3p)
      if (allocated(z1p)) deallocate(w1p)
      if (allocated(vcoul)) deallocate(vcoul)
      if (allocated(vcoult)) deallocate(vcoul)
      if (allocated(zmel1)) deallocate(zmel1)
      if (allocated(zmel3)) deallocate(zmel3)
      if (allocated(zw_)) deallocate(zw_)
      if (allocated(zwz2)) deallocate(zwz2)
c      if (allocated(zw2)) deallocate(zw2)
      if (allocated(zmel2)) deallocate(zmel2)
      if (allocated(zw3)) deallocate(zw3)
      if (allocated(uaa)) deallocate(uaa)
 1001 continue
      if (allocated(expikt)) deallocate(expikt)

      if (verbose()>=40) then
        itp = idalloc(' ',11,1,1)
      else
        itp = idalloc(' ',1,1,1)
        write(*,"(/' Max dynamic allocation sxcf_fal3z:',i7,' MB'/)") itp
      endif

      end

