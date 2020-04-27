      subroutine amagnc(nbas,nl,ipc,rhos,nrhos,qnu,eula,neul,mode,
     .  amag,aamom,bxc)
C- Printout magnetic moments in unit cell
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   rhos  :spin density matrix (for mode=0)
Ci         :rhos should be hermitian in spin space, but may not be owing
Ci         :to energy integration errors in the complex plane.
Ci         :amagnc uses a symmetrized form to minimize errors.
Ci   qnu   :moments (for mode=1)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l- and m-independent,
Ci         :nl if Euler are l-dependent and m-independent
Ci         :nl**2 if Euler are l- and m-dependent
Ci   mode  :0 use spin density matrix to make moments along T
Ci         :1 use qnus to make moments along qnu
Ci         :Add 2 to generate average magnetization direction bxc
Co Outputs
Co   amag(1..3): net system magnetic moment
Co   aamom :local magnetic moments
Co   bxc   :average magnetization direction for each class
C          :(made when mode>=2)
Cr Remarks
Cr   Definition of rho in terms of M: (standard definition of sigma)
Cr      rho = M . sigma/2
Cr   Pauli matrices sigma:
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr   Given rho, M can be obtain from:
Cr     M_x =  2 Re(rho21) = Re (rho12+rho21)
Cr     M_y =  2 Im(rho21) = Im (rho21-rho12)
Cr     M_z =  (rho11)-(rho22)
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be properly hermitian, e.g. when rhos is generated
Cr   by a Green's function technique.
Cb Bugs
Cb   This routine does not symmetrize over class.
Cb   The first member of the class is assumed to be representative.
Cu Updates
Cu   09 Apr 11 can generate bxc; also small bug fixes
Cu   21 Apr 04 Revised to properly accomodate m-dependent Euler angles
Cu   17 Feb 03 Revised amagnc; cleaner and bug fixes.
C ----------------------------------------------------------------------
      implicit none
C Passed variables
      integer nbas,neul,nl,nrhos,ipc(nbas),mode
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,*),
     .  qnu(3,nl,2,*),amag(3),aamom(nbas),bxc(3,*)
C Local variables
      logical lrhol
      integer i,ib,ic,lp1,lgunit,ipr,k,stdo,ilm,l,m
      double precision alphan,betan,gamman,arhol,arhom,sarho,
     .  rotg(3,3),amom,amlm(3),amgm(3),saml(3),samg(3),amgl(3),aml(3)
      integer PRT1,PRT2
      parameter (PRT1=30,PRT2=40)

C     character*1 cmode(2)
C     data cmode /'o','i'/

      call rx('amagnc requires noncollinear package')

      end
      subroutine amagn2(nl,nclass,nbas,ipc,eula,neul,bxc,qnu,nrhos,rhos)
C- Renormalize magnetic part of qnu, corresponding to average direction
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   nrhos :number of channels which spin density-matrix is stored
Ci   bxc   :average direction for XC field (amagnc.f)
Ci   rhos  :spin density-matrix.  rhos++ and rhos-- are modified.
Co Outputs
Co   qnu   :energy-weighted moments of the sphere charges
Cl Local variables
Cl         :
Cr Remarks
Cr   If bxc is scaled by -1, the magnetic part of qnu is scaled by -1.
Cr   Only the product is relevant.  This routine scales bxc by -1 if
Cr   the moment from the given qnu is negative.
Cb Bugs
Cb   This routine does not symmetrize over class.
Cb   The first member of the class is assumed to be representative.
Cu Updates
Cu   09 Apr 11 Completely rewritten
Cu   09 Apr 04  First created.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nclass,nbas,neul,nrhos,ipc(nbas)
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,nclass)
      double precision qnu(3,nl,2,1),bxc(3,*)
C ... Local parameters
      integer ilm,l,j,m,k,i,ic,ib,iclbsj
      double precision mbar(3,0:2),amg(3,0:2),ammi(3),ddot,rotm(3,3),
     .  alpha,beta,gamma,ql,amli(0:2),qp,qm,eulat(neul,3),ph,th,
     .  hatbxc(3),pi,rotb(3,3),rotmb(3,3),rotmb3,rhom(3),qli(0:2),fac

      call rx('amagnc require noncollinear package')

      end
