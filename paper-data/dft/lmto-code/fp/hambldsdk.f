      subroutine hambldsdk(mode,s_site,s_spec,s_lat,s_ham,nbas,isp,ispc,q,
     .  k1,k2,k3,qsig,nqsig,smpot,vconst,lcplxp,lso,alfa,ndimh,z,nev,dsdk)
C- Numerical k-differentation of dSigma/dk
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos sighh sighk sigkk pihh pihk pikk sighhx
Ci                 sighkx sigkkx pihhx pihkx pikkx
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hambl augmbl bstrux smhsbl hsibq hsubblock hsibq2
Cio                hsibq4
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt rsma pz name orbp ngcut
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hambl augmbl bstrux uspecb smhsbl hsibq tbhsi
Cio                hsubblock hsibq2 hsibq4
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:   alat  vol awald tol nkd nkq gmax nabc ng kv
Ci                 kv2
Co     Stored:     *
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:cg indxcg jcg cy qlv dlv gv igv kv igv2
Cio    Passed to:  hambl augmbl bstrux hxpbl ghibl hklbl gklbl hxpgbl
Cio                ghigbl hklgbl smhsbl hhibl phhibl hsmbl hsibq sugvec
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  oveps   nprs
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb iaxs hrs
Cio    Passed to:  hambl sopert3
Ci Inputs
Ci   mode  :Passed to hambls
Ci         :Routine does nothing if 1s digit does not contain 4
Ci   nbas  :size of basis
Ci   isp   :spin index
Ci   ispc  :index where to store h,s. (for noncollinear mode)
Ci   q     :Bloch vector (k-point)
Ci   k1,k2,k3 dimensions of smpot
Ci   qsig  :q-points at which sigma can be computed
Ci   nqsig :number of qsig
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   vconst:additional constant potential
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   osig,otau,oppi  augmentation matrices
Ci   alfa  :add alfa * overlap to hamiltonian
Ci         :This is for stability in evals.  Preferably alfa=0
Ci   ndimh :dimension of hamiltonian and overlap h,s
Co Outputs
Co   mode  :if mode returns <0, hambls finds (approximate)
Co         :evals and evecs internally; see also description of h,s
Co   h     :Hamiltonian matrix, unless 1s digit or 1000s digit mode is nonzero.
Co         :In that case mode is returned as -1, and
Co         :what is returned in h depends on 1000s digit mode, e.g.
Co         :if it is 1, h is diagonalized internally, and on exit
Co         :h is returned as a diagonal matrix, whose elements are eigenvalues;
Co         :s contains eigenvectors.  The method for doing this, and the
Co         :motivation for it, are described in Remarks.
Co   s     :overlap matrix, unless 1s digit or 1000s digit mode is nonzero.
Co         :See description for h.
Cl Local variables
Cr Remarks
Cu Updates
Cu   28 Mar 17
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,isp,ispc,k1,k2,k3,lcplxp,ndimh,nqsig,lso,nev
      double precision qsig(3,nqsig),alfa
      double precision q(3),vconst
      double complex smpot(k1,k2,k3),z(ndimh,ndimh)
      double complex dsdk(ndimh,ndimh,3,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Dynamically allocated arrays
      complex(8),allocatable:: h(:,:),sigq(:,:,:)
C ... Local parameters
      character outs*256
      logical ldsdk
      integer i,j,ii,lrsig,nevl,mode0,model,lwsig,ismidb
      real(8) :: dq(3),qs(3),xx(1)
      logical, parameter :: T=.true., F=.false.
      integer, parameter :: NULLI=-99999
      real(8), parameter :: tol=1d-7, delq=1d-3
      complex(8), parameter :: zer=(0d0,0d0),one=(1d0,0d0)
      procedure(integer) :: nglob
      integer hreal,nl,mxorb,nsp,nttabs,ndhrs,nlmto

C --- Setup ---
      mode0 = mod(mode,10); lrsig = mod(mode/10,10); lwsig = mod(mode/1000,10)
      ldsdk = mode0 >= 4
      if (.not. ldsdk .or. lrsig == 0 .or. lwsig /= 0) return
      if (lso /= 0) call rx('hamblsdsdk not ready for so')

!     call zprm('z',2,z,ndimh,ndimh,nev)

      call tcn('hamblsdsdk')

      hreal = 0
      i = 100000 + 4000 + 000 + 40*(1-hreal) + 10
      nl = nglob('nl'); mxorb = nglob('mxorb'); nsp  = nglob('nsp')
      nttabs = s_ham%nprs(nbas+1)
      nlmto = s_ham%nlmto       ! Number of LMTO's ... independent of k
      ndhrs = s_ham%ndhrs
      call shorps(1,s_lat%qlat,(/72,2,2/),q,qs)
      allocate(h(ndimh,ndimh))
      do  j  = 1, 3
        call bloch(i+j*1000000,qs,nl,s_lat%plat,mxorb,s_ham%iprmb,1,nttabs,s_ham%iaxs,s_ham%hrs,
     .    ndhrs,isp,nsp,nlmto,nlmto,0,ndimh,0,ndimh,0,dsdk(1,1,j,ispc),xx,xx)
C       print *, j; call zprm('raw dsdk, orb bas',2,dsdk(1,1,j,ispc),ndimh,nlmto,nlmto)
        call phmbls(1,ndimh,nev,xx,xx,h,dsdk(1,1,j,ispc),z,z,dsdk(1,1,j,ispc))
C       print *, j; call zprm('<z dS/dk z>',2,dsdk(1,1,j,ispc),ndimh,nev,nev)
      enddo
      deallocate(h)

C ... dS/dk By brute force, including renormalization of high lying sigma
C      allocate(h(ndimh,ndimh),sigq(ndimh,ndimh,2))
C      do  j = 1, 3
C        dq = 0
C        do  i = 1, 2
C          ii = 3-2*i
C          dq(j) = delq*ii
C          qs = q + dq
C          print *, j,ii,sngl(qs)
C          model = mode + 5000 - 4 ! Return sigma in orbital basis, do not calculate dS/dk
C          call hambls(model,s_site,s_spec,s_lat,s_ham,nbas,isp,ispc,qs,
C     .      k1,k2,k3,qsig,nqsig,smpot,vconst,lcplxp,lso,alfa,ndimh,h,sigq(1,1,i),dsdk,ismidb,nevl)
C!          print *, ii,j; call zprm('sigm(orb)',2,sigq(1,1,i),ndimh,ndimh,ndimh)
C        enddo
C        call daxpy(2*ndimh**2,-1d0,sigq(1,1,2),1,sigq(1,1,1),1)
C        call dscal(2*ndimh**2,1d0/(2*delq),sigq(1,1,1),1)
C!        print *, j; call zprm('dS/dk(orb)',2,sigq,ndimh,ndimh,ndimh)
C        call phmbls(1,ndimh,nev,xx,xx,h,sigq,z,z,sigq)
C        print *, j; call zprm('<z dS/dk z>',2,sigq,ndimh,nev,nev)
C        if (nev /= ndimh) call rx('hamblsdsdk: oops')
C        call dcopy(2*ndimh**2,sigq,1,dsdk(1,1,j,ispc),1)
C      enddo

C ... dH/dk by brute force, no renormalization of high lying sigma
C      allocate(h(ndimh,ndimh),sigq(ndimh,ndimh,2))
C      do  j = 1, 3
C        dq = 0
C        do  i = 1, 2
C          ii = 3-2*i
C          dq(j) = delq*ii
C          qs = q + dq
C          print *, j,ii,sngl(qs)
C          model = mode + 5000 - 4 ! Return sigma in orbital basis, do not calculate dS/dk
C          model = mod(model,100)
C          call hambls(model,s_site,s_spec,s_lat,s_ham,nbas,isp,ispc,qs,
C     .      k1,k2,k3,qsig,nqsig,smpot,vconst,lcplxp,lso,alfa,ndimh,sigq(1,1,i),h,dsdk,ismidb,nevl)
C!          print *, ii,j; call zprm('h(orb)',2,sigq(1,1,i),ndimh,ndimh,ndimh)
C        enddo
C        call daxpy(2*ndimh**2,-1d0,sigq(1,1,2),1,sigq(1,1,1),1)
C        call dscal(2*ndimh**2,1d0/(2*delq),sigq(1,1,1),1)
CC       print *, j; call zprm('dH/dk(orb)',2,sigq,ndimh,ndimh,ndimh)
C        call phmbls(1,ndimh,nev,xx,xx,h,sigq,z,z,sigq)
CC       print *, j; call zprm('<z dH/dk z>',2,sigq,ndimh,nev,nev)
C        if (nev /= ndimh) call rx('hamblsdsdk: oops')
C        call dcopy(2*ndimh**2,sigq,1,dsdk(1,1,j,ispc),1)
C      enddo


      call tcx('hamblsdsdk')
      end
