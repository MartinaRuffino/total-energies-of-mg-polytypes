      subroutine optinq(s_ctrl,s_optic,s_bz,eband,nbmax,nsp,nspc,efermi,
     .  ipq,alat,plat,nfilm,nempm,nptmx,npol,njdosw,optmt,jdosw,opt)
C- BZ integration of Im(eps) or joint DOS by tetrahedron method
C ----------------------------------------------------------------------
Cio Structures
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic lpart dw window ocrng unrng esciss iq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  optin2
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp lmet lshft ntet n range w
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  optin2
Ci Inputs:
Ci   eband :energy bands for irreducible part of BZ
Ci   nbmax :leading dimension of eband
Ci   nsp   :number of spins if nspc=1, 1 if nspc=2
Ci   nspc  :2 for noncoll or S/O; otherwise 1
Ci   efermi:Fermi level
Ci         :NB: efermi=-99 is a flag, indicating that the Fermi
Ci         :level is not set.
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nfilm,nempm: dimensions optmt.  Should be nfiup-nfilo+1,nemup-nemlo+1
Ci   nptmx :leading dimension of opt
Ci   npol  :number of polarizations
Ci         :Case njdosw=0 or loptic>0: npol should be 3
Ci         :Case njdosw>0 and loptic<0: npol should be 1+njdosw
Ci   optmt :<i|grad|j> connecting occ i with unocc j
Ci         :optmt is not used when making joint DOS (loptc=-1)
Ci         :only states i between nfilo and nfiup, and states j
Ci         :between nemlo and nemup are supplied
Co Outputs:
Co   popt  :opt, resolved in one of the following ways.
Co   opt   :Im(epsilon), unscaled by weight 16*pi**2/vol/e**2
Co         :opt is either :
Co         :(1) generated and the scaled epsilon written to file 'opt'
Co          (when information for all k-points is available; see ikp);
Co         :(2) partial contribution added to the unscaled epsilon
Co          (when information for 1 k-points is available; see ikp);
Co         :NB: output may be joint DOS instead of Im eps; see jdos
Co         :    below
Co         :NB: opt is NOT initialized by this routine if accumulating
Co         :    only 1 kp (ikp>0)
Co         :Caller should set opt=0 before first call to optinq.
Cl Local variables and switches
Cl  loptic  What optical property to calculate
Cl          Should be same as s_ctrl%loptc
Cl          1 generate Im(eps) (spin 1 only, if partial decompsn)
Cl          2 generate Im(eps) spin 2 (needed for partial decompsn)
Cl         -1 generate joint DOS
Cl         -2: generate joint DOS, spin 2 (use with LTET=3)
Cl         -3: generate up-down joint DOS (LTET=3)
Cl         -4: generate down-up joint DOS (LTET=3)
Cl         -5: generate spin-up single DOS (LTET=3)
Cl         -6: generate spin-dn single DOS (LTET=3)
Cl         Add 20 for Second Harmonic Generation
Cl   lpart :Switch controlling how to resolve JDOS or Im eps into parts
Cl         :1s digit
Cl         :0  No decomposition
Cl         :1  Resolve eps or dos into (occ,unocc) parts
Cl         :2  Resolve eps or dos by k
Cl         :3  Both 1 and 2
Cl         :10s digit
Cl         :0  save as ASCII file
Cl         :1  save as binary file
Cl   nfilo :first occupied band to include in Im eps or JDOS
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   nfiup :last  occupied band to include in Im eps or JDOS
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   nemlo :first unoccupied band
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   nemup :last  unoccupied band
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   liq   :True if to compute JDOS or Im eps with unocc shifted by q
Cb Bugs
Cb   Set of nfilo,nfiup,nemlo,nemup has duplicate nfilox,nfiupx,nemlox,nemupx
Cb   2nd set is used in tetwtq, but not clear why it should be distinct.
Cb   Possibly affects tetrahedron weights?
Cr Remarks
Cr   Uses tetwtq to makes joint density of states or Im(eps)
Cr   tetwtq calculates a (occ,unocc,k) resolved set of weights
Cr   for sequence of histograms 1,2,3,...
Cr     whw(i,ib,jb,k) =
Cr        \int_om_i^om_i+1 d\omega \times
Cr        \int d^3k f(e(k))(1-f(e(q+k)))\delta(omg - e(q+k) + e(k))
Cr
Cr   tetwtt uses these weights to assemble Im(eps) or JDOS.
Cr
Cu Updates
Cu   27 Aug 14 Better printout for memory management
Cu   21 Jun 14 s_optic%loptic=1 handles both spins unless partial decompsn
Cu   06 Jan 14 Implemented s_optic%loptic=2
Cu   12 Nov 12 migrated structures to f90 pointers
Cu   09 Aug 12 jdosw is spin polarized
Cu   10 Nov 11 Begin migration to f90 structures
Cu   11 Sep 09 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbmax,nsp,nspc,nfilm,nempm,nptmx,ipq(*),npol,njdosw
      double precision efermi,alat,plat(3,3)
      double precision opt(nptmx,npol,nsp),jdosw(max(njdosw,1),nfilm+nempm,*)
      double precision optmt(3,nfilm,nempm,nsp,*),
     .                 eband(nbmax,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_optic):: s_optic
      type(str_bz)::    s_bz
C ... Local parameters
      logical jdos,metal,ltop,liq,ltmp,memcheck
      integer i,ilo,ipr,iq1,iup,j,jpm,k,lgunit,isp,jsp,ipol,lmet,loptic,
     .  lpart,mxbnds,nband,nemlo,nemup,nfilo,nfiup,nkp,npts,ntet,nwhis,
     .  stdo,nfilmx,nempmx,nkpw,idum,nfilmw,nempw,nfilox,nfiupx,nemlox,
     .  nemupx,ibw,jbw,lpart0,lpart1,ksp,npass
      integer idalloc,allocvb
      integer ocrng(2),unrng(2)
      integer nkabc(3),n1,n2,n3,lshft(3),i1,i2,i3,nqbz,nqbzw,iqshft(3)
      integer npm,nctot,ncc,nwgtx,jhwtot,nhwtot,iq,jj1,jj2,jj3,ifac(3)
      double precision emin,emax,pi,esciss,vol,de(3),saveef,fac,
     .  optrng(4),q(3),qlat(3,3),tiny,qb(3,3),xv(10)
      real demaxx
      equivalence (emin,optrng(1)),(emax,optrng(2))
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      real(8),allocatable :: qbzw(:,:),ecore(:,:)
      real(8),allocatable,target :: efbz1(:,:,:,:),efbz2(:,:,:,:)
      real(8),pointer:: eocc(:,:,:,:),eunocc(:,:,:,:)
      real(8),allocatable :: emesh(:),frhis(:),popt(:,:,:,:,:)
      real(8),allocatable :: wtthis(:,:),wtphis(:,:,:,:,:),whw(:)
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      integer(4),allocatable::
     .  ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      integer,allocatable :: ib1bz(:),idtet(:,:),nwgt(:,:),nwgttt(:,:),
     .  nib(:,:,:),njb(:,:,:)
      integer :: noccxvx(2)=-9999
      logical,allocatable :: iwgt(:,:,:,:)
      character strn*32
      logical, parameter :: T=.true., F=.false.
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
C      double precision qk
C      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
C     .                    (jj2*ifac(2)-1)*qb(k,2) +
C     .                    (jj3*ifac(3)-1)*qb(k,3)
C ... External calls
      external cspline,csplint,dinv33,dosmsh,dpzero,dscal,getpr,
     .         hisrange,info2,info5,info8,optin2,optind,pair_index,
     .         pralloc,rxi,rxx,tcn,tcx,tetfbz,tetwtq,tetwtt

      call rx('OPTINQ: OPTICS package not installed')

      end
