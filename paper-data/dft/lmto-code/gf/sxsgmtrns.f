      subroutine sxsgmtrns(isp,nsp,sig,ldim,lihdim,pph,sk)
C- Generate ASA part of the  Overlap
C ---------------------------------------------------------------------
Ci Inputs
Ci   sig   :sigma in gamma representation
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lihdim:dimensions ccd and pph
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :pph(1..5,i,is): parms for ith RL and is(th) spin channel.
Ci         :pph(1) : enu
Ci         :pph(2) : calpha
Ci         :pph(3) : sqrdel
Ci         :pph(4) : palpha
Ci         :pph(5) : oalp
Ci   sk    :structure constants, s^beta
Ci   wk    :work array of length ldim
Co Outputs
Co  okl,okr: transformations left and right for sigma
Co  sig    :transformed sigma
Cr Remarks
Cr
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | psi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   The first form uses phidot = phidot^alpha; the second form uses
Cr     phidot^alpha = phidot^gamma + o phi; and calls 'dot'=phidot^gamma
Cr   Note that <phi|phi>=1 <phidot|phi>=o <phidot|phidot>=p
Cr   Considering the LL block only, the ASA part of the overlap is:
Cr    <psi|psi>_ASA = <phi|phi> + h<phidot|phi> + h.c.
Cr                    + h <phidot|phidot> h
Cr                  = 1 + ho + oh + hph
Cr
Cr   To work directly with  D = srdel S srdel, rather
Cr   that h = C-enu + D, the diagonal parts connecting C-enu
Cr   in the one-, two-, three-center terms are reshuffled.  Thus
Cr    <psi|psi>_ASA = 1 + ho + oh + hph
Cr                  = 1 + (C-e+D)o + o(C-e+D) + (C-e+D) p (C-e+D)
Cr                  = 1 + 2(C-e)o + (C-e)^2 p      (one center)
Cr                  + D(o+p(C-e)) + (o+p(C-e))D    (two center)
Cr                  + D p D                        (three center)
Cr
Cr   The hamiltonian <psi|H|psi>_ASA has a corresponding structure
Cr   with similar 1-, 2- and 3- center terms; but the diagonal parts
Cr   are calculated from <phi|H|phi>, <phidot|H|phi>, <phidot|H|phidot>
Cr   and are passed in array ccd.
Cr
Cr   Also, In the notation below, <k|k> = <kappa|kappa>
Cu   Updates
Cu Dec 04 04 made from the overlap and matrix Hamiltonian (T.Sandu)
Cu
Cu
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer lihdim,ldim,isp,nsp
      double precision pph(5,lihdim,nsp),wk(ldim)
      double precision sk(ldim,ldim*2),sig(ldim,ldim*2),sg(ldim,ldim*2)
      double precision okl(ldim,ldim*2),okr(ldim,ldim*2)
C Local parameters
      integer i,j,l2,ofpr,ofpl
C External calls
      external yygemm

      call tcn('sxsgmtrns')
      print*,'into sxsgmtrns'
      l2 = ldim**2

cc      print *,'isp nsp lihdim ldim',isp,nsp,lihdim,ldim
      ofpr = isp*lihdim
      ofpl = isp*lihdim
cc      print*,'ofpr,ofpl',ofpr,ofpl
C --- Make (calpha-enu)*dij+d*S*d with d = sqrt(delta) and
C ---  dij is Kronecker delta function ---
cc      call yprm('salpha',12,sk,ldim*ldim,ldim,ldim,ldim)
      call makdsd(1,ldim,ldim,ldim,ldim,0,0,pph(1,1,isp),sk,sk)
ccc      call yprm('halpha',12,sk,ldim*ldim,ldim,ldim,ldim)
ccc      call yprm('sigsxtrns',2,sig,ldim*ldim,ldim,ldim,ldim)


c      do  44  j = 1, ldim
c      do  44  i = 1, ldim
c
c        sk(i,j) = 0d0

c  44  continue



C --- O += o*D & D*o
      do  1  i = 1, ldim
    1 wk(i) = pph(5,i,isp)

      do  2  j = 1, ldim
      do  2  i = 1, ldim

        okr(i,j) = sk(i,j)*(wk(i) )
        okr(l2+i,j) = sk(l2+i,j)*(wk(i))
        okl(i,j) = sk(i,j)*(wk(j) )
        okl(l2+i,j) = sk(l2+i,j)*(wk(j))
    2 continue
C     call yprm('okr-O+= 2c',12,okr,ldim*ldim,ldim,ldim,ldim)
C     call yprm('O+= 2c',12,okl,ldim*ldim,ldim,ldim,ldim)
C --- O += 1+o*D & 1+D*o
      do  3  i = 1, ldim
        okr(i,i) = okr(i,i) + 1d0
        okl(i,i) = okl(i,i) + 1d0
    3 continue

c      call yprm('okr-O+= 2c',12,okr,ldim*ldim,ldim,ldim,ldim)
c      call yprm('Okl= 2c',12,okl,ldim*ldim,ldim,ldim,ldim)

c      do  4  j = 1, ldim
c      do  4  i = 1, ldim

c        sg(i,j) = 0d0
c        sig(i,j) = 1d0*sig(i,j)

c    4 continue

C ---transform sigma from gamma to alpha representation
C ----(1+D*o)*sigma*(1+o*D)
      call yygemm('N','N',ldim,ldim,ldim,1d0,sig(1,1),sig(l2+1,1),
     .   ldim,okr,okr(l2+1,1),ldim,0d0,sg(1,1),sg(l2+1,1),ldim)
      call yygemm('N','N',ldim,ldim,ldim,1d0,okl(1,1),okl(l2+1,1),
     .   ldim,sg,sg(l2+1,1),ldim,0d0,sig(1,1),sig(l2+1,1),ldim)

C.....make -sigma in alpha representation
c       print*,'before ------'
c       do 555 lj = 1,ldim
c       do 555 li = 1,ldim*2

c 555   sig(lj,li)= (-1d0)*sig(lj,li)

      call tcx('sxsgmtrns')
      end

      subroutine sqrdppc(s_ctrl,s_site,s_ham,s_pot,s_spec,vshft,zp,
     .  sqrdpp)
C- Calculate sqrt of p-dot for Green's function
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lrel nl nbasp nspin lham nbas nclass nspec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkpotf makpfz
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pnu spec class clabel v0
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkpotf makpfz
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  hord ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  mkpotf makpfz
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz ves
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pp pprel sop
Cio    Passed to:  mkpotf makpfz
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr z rmt lmxa hcr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkpotf makpfz
Ci Inputs
Ci   zp    :complex energy
Cio Inputs/ Outputs
Cio  spot  : pot. information is packed in structure spot; see Remarks
Cr Remarks
Cr   This routine generates energy-dependent potential
Cr   parameters needed to create the hamiltonian for a Green's function.
Cr
Cr   For 2nd-generation LMTO-ASA, these are potential functions
Cr    pot->{pf,dpf,ddpf,dddpf}.
Cr
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      double precision vshft(*),zp(2)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: sqdp(:)
      integer, allocatable :: ipc(:),ipdc(:)
      real(8), allocatable :: potp(:)
C     real(8), allocatable :: ovl(:)
C ... Local parameters
c      logical bittst
      integer hord,nbasp,nl,nsp,ldham(16),ldim,lidim,lihdim,
     .  iopt,lham,nspc,lrel
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (nspc,ldham(4))
      integer nbas,nclass,nlspc
      double precision xx
      double complex sqrdpp(1)

      lrel = mod(s_ctrl%lrel,10)

      nl = s_ctrl%nl
      nbasp = s_ctrl%nbasp
      nsp = s_ctrl%nspin
      lham = s_ctrl%lham
      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      hord = s_ham%hord
      ldham = s_ham%ldham
      allocate(sqdp(lihdim*nsp))

c      print *,'nclass',nclass,'nbasp',nbasp,'nclass',nclass

      if (lrel == 2) call rx('sqdpp:not ready for relativistic case')

      allocate(ipc(nbasp))
      call rx('sxsgmtrns: update sid')
      iopt = 0
      if (hord == 3) iopt = 1
      if (hord == 1) iopt = 2
      if (hord == 4) iopt = 3 + 3400000

C      print *, '!zp'
C      zp(1) = -.25d0
C      zp(2) = .05d0

C ... Make sure that potential parameters are in alpha reperesentation
      nlspc = nl*nsp*nclass
      allocate(potp(6*nlspc))
c      print*,'sqrtpp'
c      print*,'zp1 zp2',zp(1),zp(2)
c      print*,'pp'
c      call ppprnt1(nl,nsp,nclass,s_pot%pp)
      call dcopy(6*nlspc,s_pot%pp,1,potp,1)
      call getnewpp(nl,nsp,nclass,potp,s_pot%pp)
C      allocate(ovl(nlspc))
ccc      call pptrns(0,nl,ipc,nclass,nsp,w(oalph),nbas,s_pot%pp,ovl)
c      print*,'sqrtpp'
c      call pptrns(0,nl,ipc,nclass,nsp,w(oalph),nbas,s_pot%pp,ovl)
c      print*,'opotp'
c      call ppprnt1(nl,nsp,nclass,s_pot%pp)
C      deallocate(ovl)

C       Make P, P-dot, -1/2 P-dotdot/P-dot, (-1/2 P-dotdot/P-dot)-dot
C       The last is used for linear response
c        iopt = iopt + 10*(2**0+2**1+2**2+2**4)
C  ...with sqrt(dpf)
c        iopt = iopt + 10*(2**3)
C  ...with dpf
c        print*,'iopt before assign',iopt
      iopt = iopt + 10*(2**1)
C  trasf to gamma for debugging
c        iopt = iopt + 20000
c        print*,'iopt',iopt

c        print*,'into alpha and before mkpotf'
C        call mkpotf(iopt,sctrl,
C     .    s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,ipc,xx,vshft,zp,w(opf),
C     .    w(odpf),w(oddpf),xx,w(odddpf),xx,xx,xx)
C ... calculate sqrt(dpf)
C        call mkpotf(iopt,sctrl,
C     .    s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,
C     .    ipc,xx,vshft,zp,xx,
C     .    sqdp,xx,w,xx,xx,xx,xx)
C...calculate dpf
      allocate(ipdc(nbasp)); ipdc = 0
      call rx('update call to mkpotf')
C      call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,pfdim,ipc,ipdc,xx,vshft,
C     .  zp,xx,sqdp,xx,xx,xx,xx,xx,xx)
c        call pfprnt1(lihdim,1,sqdp)
      deallocate(ipdc)
      call dppvect(lidim,sqdp,sqrdpp)

cc       call zprm('dpf',2,sqdp,lihdim,lihdim,nsp)
cc       call zprm('dpf-trans',2,sqrdpp,lihdim,lihdim,nsp)

C      call zprm('palp',2,w(opalp),lihdim,lihdim,nsp)
C      call zprm('pfr',2,w(opfr),lihdim,lihdim,4)

      call dcopy(6*nlspc,potp,1,s_pot%pp,1)
      deallocate(sqdp,ipc,potp)
      end

      subroutine dppvect(ldg1,sqdpin,sqdpout)
C-

C ----------------------------------------------------------------------
      implicit none
      integer ldg1
      double complex sqdpin(ldg1,2)
      double complex sqdpout(ldg1,2)
      integer i,j

        do  10  j = 1, 2
        do  10  i = 1, ldg1
          sqdpout(i,j) = sqdpin(i,j)
   10   continue
      end

      subroutine sclsgm(isp,ldg1,sqdpf,gf)
C- Convert 2nd generation g_ij = (P-S)^-1 to G_ij by energy scaling

C ----------------------------------------------------------------------
      implicit none
      integer ldg1,isp
      double precision gf(ldg1,ldg1,2)
      double complex sqdpf(ldg1,2)
C ... Passed parameters
      double complex xxc
      integer i,j

C --- sigma_ij -> -(1/(sqrt(P^dot_i))* sigma_ij* 1/(sqrt(P^dot_j)) ---
        do  10  j = 1, ldg1
        do  10  i = 1, ldg1
C ...first sigma_ij --> -sigma_ij
          gf(i,j,1) = (-1d0)*gf(i,j,1)
          gf(i,j,2) = (-1d0)*gf(i,j,2)
C... make the scaling
          xxc =(1d0/sqdpf(i,isp))*(1d0/sqdpf(j,isp))
     .     *dcmplx(gf(i,j,1),gf(i,j,2))
          gf(i,j,1) = dble(xxc)
          gf(i,j,2) = dimag(xxc)
   10   continue
      end


      subroutine sigmadd(ldg1,sigma,gf1)
C- Add (-Sigma) to P-S

C ----------------------------------------------------------------------
      implicit none
      integer ldg1
      double precision gf1(ldg1,ldg1,2)
      double precision sigma(ldg1,ldg1,2)
C ... Passed parameters
C     double complex xxc
      integer i,j

C --- sigma_ij -> -(1/(sqrt(P^dot_i))* sigma_ij* 1/(sqrt(P^dot_j)) ---
        do  10  j = 1, ldg1
        do  10  i = 1, ldg1
C ...(P-S)_ij -sigma_ij

c          gf1(i,j,1) = gf1(i,j,1) - sigma(i,j,1)
c          gf1(i,j,2) = gf1(i,j,2) - sigma(i,j,2)
          gf1(i,j,1) = gf1(i,j,1) + sigma(i,j,1)
          gf1(i,j,2) = gf1(i,j,2) + sigma(i,j,2)
   10   continue
      end




      subroutine ppprnt1(nl,nsp,nclass,PP0)
C- Writes out the response matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nRLc  :dimension of P0
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of k-points for which P0 is available
Ci   P0    :response matrix = 1/i G delta v G
Co Outputs
Co    P0 written in a file
Cl Local variables
Cr Remarks
Cr Made Feb 20 2004
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nsp,nclass
      double precision PP0(6,nl,nsp,nclass)
C ... Local parameters
      integer i1,ik,ii,ij,fisi,fopna
      character*(10) fmt
      character a
      fmt = '(12f16.8)'
C      fisi = fopna('p0out',29,0)
C      print *, 'fisi=',fisi
C      rewind fisi
C      call prjrsp(nsp,1,2,nRLc,nkp,P0)
C --- Writes P0(q)
      fisi = fopna('PP0out',29,0)
      rewind fisi
      do  90  ik = 1, nclass

      do  90  i1 = 1, nsp
      do  88  ij = 1, nl
   88 write (fisi,fmt)(PP0(ii,ij,i1,ik), ii=1,6)
c      call fclose(fisi)
c      print *, 'kp= ',ik

   90 continue

      call fclose(fisi)
      read *, a

      end


c      subroutine sclfact(isp,ldg1,sqdpf,scfact)
C- Convert 2nd generation sigma_ij = to (P-S)^-1-like object

C ----------------------------------------------------------------------
C     implicit none
c      integer ldg1,isp
c      double precision scfact(ldg1,ldg1,2)
c      double complex sqdpf(ldg1,2)
C ... Passed parameters
c      double complex xxc
c      integer i,j

C --- sigma_ij <- 1/(sqrt(P^dot_i))* sigma_ij* 1/(sqrt(P^dot_j)
c        do  10  j = 1, ldg1
c        do  10  i = 1, ldg1
c          xxc = (1d0/sqdpf(i,isp))*(1d0/dpf(j,isp))
c          scfact(i,j,1) = dble(xxc)
c          scfact(i,j,2) = dimag(xxc)
c   10   continue
c      end


      subroutine sclsgm1(isp,ldg1,dpf,gf)
C- Convert 2nd generation sigma_ij to (P-S)^-1 like object

C ----------------------------------------------------------------------
      implicit none
      integer ldg1,isp
      double precision gf(ldg1,ldg1,2)
      double complex dpf(ldg1,2)
C ... Passed parameters
      double complex xxc
      integer i,j

C --- sigma_ij <- 1/(sqrt(P^dot_i))* sigma_ij* 1/(sqrt(P^dot_j)) ---
        do  10  j = 1, ldg1
        do  10  i = 1, ldg1
C ...first sigma_ij --> -sigma_ij
          gf(i,j,1) = (-1d0)*gf(i,j,1)
          gf(i,j,2) = (-1d0)*gf(i,j,2)
          xxc = (sqrt(dpf(i,isp)))*(sqrt(dpf(j,isp)))
     .         *dcmplx(gf(i,j,1),gf(i,j,2))
          gf(i,j,1) = dble(xxc)
          gf(i,j,2) = dimag(xxc)
   10   continue
      end

      subroutine sclsgm2(isp,ldg1,dpf,gf)
C- Convert 2nd generation sigma_ij to (P-S)^-1 like object

C ----------------------------------------------------------------------
      implicit none
      integer ldg1,isp
      double precision gf(ldg1,ldg1,2)
      double complex dpf(ldg1,2)
C ... Passed parameters
      double complex xxc
      integer i,j

C --- sigma_ij <- 1/(sqrt(P^dot_i))* sigma_ij* 1/(sqrt(P^dot_j)) ---
        do  10  j = 1, ldg1
        do  10  i = 1, ldg1
C ...first sigma_ij --> -sigma_ij
C....commented out since sigma has reversed sign
cccccc          gf(i,j,1) = (-1d0)*gf(i,j,1)
cccccc          gf(i,j,2) = (-1d0)*gf(i,j,2)
c          xxc = (sqrt(dpf(i,isp)))*(sqrt(dpf(j,isp)))
c     .         *dcmplx(gf(i,j,1),gf(i,j,2))
c          print*, 'sclsgm2 aici'
          xxc = (sqrt(dpf(i,isp)))*dcmplx(dble(sqrt(dpf(j,isp))),
     .          -dimag(sqrt(dpf(j,isp))))
     .         *dcmplx(gf(i,j,1),gf(i,j,2))

          gf(i,j,1) = dble(xxc)
          gf(i,j,2) = dimag(xxc)
   10   continue
      end

      subroutine sclsgm3(isp,ldg1,dpf,gf)
C- Convert 2nd generation sigma_ij to (P-S)^-1 like object

C ----------------------------------------------------------------------
      implicit none
      integer ldg1,isp
      double precision gf(ldg1,ldg1,2)
      double complex dpf(ldg1,2)
C ... Passed parameters
      double complex xxc
      integer i,j

C --- sigma_ij <- 1/(sqrt(P^dot_i))* sigma_ij* 1/(sqrt(P^dot_j)) ---
        do  10  j = 1, ldg1
        do  10  i = 1, ldg1
C ...first sigma_ij --> -sigma_ij
cccc          gf(i,j,1) = (-1d0)*gf(i,j,1)
cccc          gf(i,j,2) = (-1d0)*gf(i,j,2)
          xxc = (1d0/dpf(i,isp))*(1d0/dpf(j,isp))
     .         *dcmplx(gf(i,j,1),gf(i,j,2))
          gf(i,j,1) = dble(xxc)
          gf(i,j,2) = dimag(xxc)
   10   continue
      end

      subroutine mksqdela(nsp,lhdim,pph,sqdela)

      integer nsp,lhdim
      double precision pph(5,lhdim,nsp),sqdela(2,lhdim,nsp),xx
      integer i,j
      do 1 j = 1, nsp
      do 1 i = 1,lhdim

       xx = pph(3,i,j)
c       if (xx < 0) xx = - xx
       sqdela(1,i,j) = xx
       sqdela(2,i,j) = 0d0

 1    continue

      end


      subroutine pfprnt1(lihdim,isp,PF0)
C- Writes out the response matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nRLc  :dimension of P0
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of k-points for which P0 is available
Ci   P0    :response matrix = 1/i G delta v G
Co Outputs
Co    P0 written in a file
Cl Local variables
Cr Remarks
Cr Made Feb 20 2004
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lihdim,isp
      double precision PF0(2,lihdim,2)
C ... Local parameters
      integer ii,fisi,fopna
      character*(10) fmt
      character a
      fmt = '(12f16.8)'
C      fisi = fopna('p0out',29,0)
C      print *, 'fisi=',fisi
C      rewind fisi
C      call prjrsp(nsp,1,2,nRLc,nkp,P0)
C --- Writes P0(q)
      fisi = fopna('PFout',29,0)
      rewind fisi

      write (fisi,fmt)(PF0(1,ii,isp), ii=1,lihdim)
      write (fisi,fmt)(PF0(2,ii,isp), ii=1,lihdim)
c      call fclose(fisi)
c      print *, 'kp= ',ik

      call fclose(fisi)
      read *, a
      end

      subroutine getntabidx(nbas,ntabs,nttabs)
      integer nbas
      integer ntabs(nbas+1)
      integer nttabs

        nttabs = ntabs(nbas+1)
      end
