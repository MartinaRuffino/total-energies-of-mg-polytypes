      subroutine optint(s_optic,s_bz,lteto,eband,nbmax,ksp,nsp,
     .  nspc,ikp,efermi,idtet,wgts,vol,nfilm,nempm,nptmx,npol,optmt,opt)
C- BZ integration of Im(eps) or joint DOS
C ----------------------------------------------------------------------
Ci Inputs:
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: loptc
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: optin2
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read: lpart dw window ocrng unrng esciss
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: optin2
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: nkabc nkp ntet lmet n range w
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: optin2
Ci   lteto :Integration mode
Ci         : 0 sampling
Ci         : 1 original tetrahedron integration
Ci         : 2 tetrahedron, calling bzjnos
Ci         : 3 tetrahedron, calling tetwtq
Ci         :-1 on-the-fly sampling
Ci   eband :energy bands for irreducible part of BZ
Ci   nbmax :leading dimension of eband
Ci   ksp   :starting spin index passed to mkjdos.  Used for
Ci         :on-the-fly sampling integration (partial contribution),
Ci         :when only information for second-spin is available.
Ci         :See Remarks
Ci   nsp   :number of spins if nspc=1, 1 if nspc=2
Ci   nspc  :2 for noncoll or S/O; otherwise 1
Ci   ikp   :current k-point (sampling integration only)
Ci         :If ikp>0, partial contribution from a single
Ci         :k-point (ikp) and spin (ksp) accumulated (see Remarks)
Ci         :ikp<0 signifies that information from all k-points is
Ci         :available to be accumulated.
Ci   efermi:Fermi level
Ci         :NB: efermi=-99 is a flag, indicating that the Fermi
Ci         :level is not set.
Ci   idtet :(0,i) no. of tetrahedra of the i'th kind
Ci         :(1-4,i) identifies the i'th tetrahedron in terms of
Ci         :four irreducible k-points
Ci   vol   :volume
Ci   nfilm,nempm: dimensions optmt.  Should be nfiup-nfilo+1,nemup-nemlo+1
Ci   nptmx :max number of bins; dimension for following 3 work arrays
Ci   npol  :number of polarizations; needed for dimensioning.
Ci         :npol should be 1 for jdos, 3 for eps.
Ci   optmt :<i|grad|j> connecting occ i with unocc j
Ci         :optmt is not used when making joint DOS (loptc<0)
Co Outputs:
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
Co         :Caller should set opt=0 before first call to optint.
Cr Remarks
Cr   Adapted from bzints to make joint density of states or Im(eps)
Cr   All energy differences between states below ef and states
Cr   above ef+emin are summed, and integrated over the BZ
Cr   Treatment near the critical points (ef and ef+emin) handled crudely
Cr   Optics package adapted from Sergey Rashkeev with Walter Lambrecht,
Cr   which was adapted from an earlier version by V. Antropov.
Cr
Cr   (Sampling only) The partial contribution to Im(eps) may be made from
Cr   a single k-point and spin.  This routine accumulates contributions from
Cr   all k-points uness ikp>0, in which case it accumulates a
Cr   partial contribution from (ikp,ksp) only.
Cl Local variables and switches
Cl   popt  :opt resolved by (unocc,occ) pairs.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   31 Dec 10 optmt from optdme is now calculated for all (nfilm,nempm)
Cu   11 Sep 09 Some reworking to accomodate future tetwtq
Cu   15 Jan 04 (Uk-Jin Roh) add partial sampling contribution to Im(eps)
Cu   20 Nov 02 (jek) extended to metals
Cu   02 Mar 01 Added scissors operator
Cu   20 Dec 00 (wrl) extended to noncollinear case
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lteto
      integer nbmax,nsp,nspc,nfilm,nempm,nptmx,ksp,ikp,npol
      integer idtet(0:4,*)
      double precision efermi,vol
      double precision opt(nptmx,npol,nsp)
      double precision optmt(3,nfilm,nempm,nsp,*),wgts(*),
     .                 eband(nbmax,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_optic):: s_optic
      type(str_bz)::    s_bz
C ... Local parameters
      logical jdos,metal,lonekp
      integer nfilo,nfiup,nemlo,nemup,nkp,ntet,ibf,ibm,ipts,k,ip,lmet
      integer fopn,lgunit,isw,loptic,lpart,lpart0,lpart1,ifil,ipr,stdo
      integer ib,jb,iq1,iq2,iq3,iq4,isp,itet,ocrng(2),unrng(2),npts
      integer nkabc(3),n1,n2,n3
      integer mpsord,iq,ifmin,ifmax,iwmin,iwmax,npopt
      real(8),allocatable :: yy(:),optpar(:,:),popt(:,:,:,:),emesh(:)
      double precision swidth,srnge,del
      double precision eci(4),ecj(4),volwgt,wtm,wt,optrng(4),emin,
     .  emax,fac,pi,ds2,ds3,de(3)
      double precision esciss
      logical cmdopt,samp,bonly
      character*120 strn
      equivalence (emin,optrng(1)),(emax,optrng(2))
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))

      call rx('OPTINT: OPTICS package not installed')


      end

      subroutine optin2(s_ctrl,s_optic,s_bz,
     .  lfac,nsp,nspc,vol,nptmx,
     .  npol,lbin,npopt,npts,emesh,popt,opt)
C- File output for Im(eps) or joint DOS, generated by optint or optinq
C ----------------------------------------------------------------------
Ci Inputs:
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: loptc
Ci   soptic:struct for parameters for optical ME; see routine uoptic
Ci     Elts read: window ocrng unrng esciss
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: nkabc nkp ntet lmet n range w
Ci   lfac  :scaling factor
Ci         :1s digit controls constant factor
Ci         :1 scale opt by 1
Ci         :2 scale opt by 2
Ci         :3 scale opt by 32d0*pi**2/vol  (Im eps)
Ci         :10s digit:
Ci         :1 scale opt by 1/E^2
Ci         :2 print data including zero and negative frequencies
Ci   nsp   :number of spins if nspc=1, 1 if nspc=2
Ci   nspc  :2 for noncoll or S/O; otherwise 1
Ci   vol   :volume
Ci   nptmx :leading dimension of opt,popt
Ci   npol  :number of polarizations; needed for dimensioning.
Ci         :npol must be 1 for jdos, 3 for eps.
Ci   lbin  :0 save popt as ASCII file
Cl         :1 save popt as binary file
Ci   npopt :number of (ib,jb) pairs eps or jdos is resolved into
Ci         :If 0, no partial decomposition is written
Ci   npts  :number of energy points
Ci   emesh :energy mesh
Ci   popt  :Im(eps) or joint DOS, resolved by (occ,unocc) or k, or both
Ci   opt   :Im(eps) or joint DOS
Co Outputs:
Co   opt:  :written to disk
Co   popt  :written to disk if npopt>0
Cr Remarks
Cr   Sanity checks:
Cr   mc -qr popt.ogan -split a 1,nr+1 1,2,nc+1 a11 a12 -csum -ccat opt.ogan -e2 x1 x2+x3+x4 --
Cr   mc popt.ogan -split a 1,nr+1 1,2,nc+1 a11 a12 -coll 1:84:3 -csum -ccat opt.ogan -e2 x1 x2 --
Cr   or:
Cr   mc popt.ogan -split a 1,nr+1 1,2,nc+1 a11 a12 -csum -ccat jdos.ogan --
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   11 Sep 09 Rearrangement to accomodate tetwtq.  Now writes popt.
Cu   15 Jan 04 (Uk-Jin Roh) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lbin
      integer lfac
      integer nsp,nspc,npts,nptmx,npol,npopt
      double precision vol,emesh(npts)
      double precision opt(nptmx,npol,nsp),popt(nptmx,npol,npopt)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_optic):: s_optic
      type(str_bz)::    s_bz
C ... Local parameters
      logical jdos
      integer nfilo,nfiup,nemlo,nemup,nkp,ntet,ipts,i,j,ip,lmet,ic
      integer ifio,ifiop,fopn,fopna,loptic,nchan
      integer ocrng(2),unrng(2)
      integer nkabc(3),n1,n2,n3
      integer mpsord
      double precision swidth,srnge
      double precision wt,ei,optrng(4),emin,emax,fac,pi
      double precision esciss
      real(8),allocatable :: pfile(:,:)
C     character*120 strn
      equivalence (emin,optrng(1)),(emax,optrng(2))
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))

      call rx('OPTIN2: OPTICS package not installed')

      end
