      subroutine gfenint(s_ctrl,s_ham,s_pot,nbas,offH,nzp,zp,wz,nzpi,
     . lidim,lhh,nspc,nsp,isp,nkp,qp,plat,gint)
C- Contribution to density-matrix from file Green's function for one energy
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham lncol nl lrel
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 lncol ldham offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  gfg2g
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:papg gmar palp pf gma dpfr ddpfr dpf ddpf
Cio    Passed to:  gfg2g
Ci Inputs
Ci   nbas  :size of basis
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   nzp   :number of energy points
Ci   zp    :complex energy
Ci   wz    :weight of complex energy
Ci   nzpi  :current point of energy index corresponding to zp
Ci   lidim :dimension of the Hamiltonian
Ci   lhh   :leading dimension of higher orbitals block
Ci   nspc  : equals 1 collinear, =2 non-collinear spins
Ci   nsp   :=1 spin polarized, =2 spin polarized
Ci   isp   :current spin =1 or =2;
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   qp    :irreducible k-points
Ci   zp    :complex energy
Ci   plat  :primitive lattice vectors, in units of alat
Co Outputs
Co   gint  :integrated GF for all irreducible q-points.
Cr Remarks
Cr Made July 04 (T.Sandu).
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   20 Oct 04 updated to spin pol. case(T. Sandu)
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nkap0,n0H,lrel,ierr,nzp
      parameter (nkap0=4,n0H=5)
      integer nbas,nkp,ifi,lhh,nsp,nspc,nzpi,lidim,isp,
     .  offH(n0H,nkap0,1)
      double precision plat(3,3),qp(3,nkp)
      double complex zp,wz,gint(lidim,lidim,nkp,nsp)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      real(8), allocatable :: ggii(:)
C ... Local parameters
      integer iq,lham,lncol,nl,iprint,ld11,ld21,ldrgn
      integer ogi
C      integer ld2,lds
C     integer jj1,jj2,jj3,k
C     double precision qk
C ... For file I/O of gf
      integer clp(9*2),iogfrs,lio
C     Bits for matching conditions (10-100s digit of lio)
c      integer MCD,MZP
c      parameter (MCD=32,MZP=64)


      integer MNC,MNB,MNT,MLB,MPL,MZP,MCD
      parameter (MNC=1,MNB=2,MNT=4,MLB=8,MPL=16,MCD=32,MZP=64)

C     Bits for selecting which data to I/O, iogfrs
      integer ISH
      parameter (ISH=1)
ccccccccccccccccccccccccc
C     character*(10) fmt1
      integer fisi,fopnx
      integer kcplx,lrd
      double complex zp0
      double precision gfe(lidim,nspc,lidim,nspc,2),xx,ghh(lhh,lhh,2)
      character*8 fisnam
C
      if (nspc == 2) call rx('gfenint:not ready for non-colinear')
C     fmt1 = '(12f14.9)'
      kcplx = 0
      zp0 = dcmplx(0d0,0d0)

C --- Setup ---
C     call wkfast(.false.)
      call tcn('gfenint')
      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      nl = s_ctrl%nl

      lidim = offH(4,1,nbas+1)
C     ld2 = lidim*lidim*nspc**2
C      lds = lidim*nspc
      kcplx = 0
C     call defdc(oggii,(lidim*nspc)**2)
      allocate(ggii(2*(lidim*nspc)**2))
C     print *, 'into gf-energ-integ'

      call iinit(clp,9*2)
      clp(3) = lidim*nspc
      clp(4) = lidim*nspc
      lrel = mod(s_ctrl%lrel,10)

c      print *,'lrel2',lrel
ccccccccccccccccccccccccccccccccccccccc

C     ....rewind file read header
C      print *, 'before opening'
       ifi = fopnx('gfqp',100,16+8+4+0,-1)

C   ... rewind file, read header
       rewind ifi
       call pshpr(1)
       zp0 = dcmplx(0d0,0d0)
       lio = 2 + 10*(MNC*1+MNB*1+MNT*0+MLB*1+MPL*0+MCD*1+MZP*0)
       if (iogfrs(lio,0,0,' ',ifi,1,1,nbas,xx,zp,qp,plat,xx,clp,
     .   xx,xx,0,0,0,xx) /= 0) c
     .   all rx('gf-energ-i failed to read header')
       call poppr


C --- For each irreducible qp, make G in its star ---
      do  iq = 1, nkp
c      print *, 'iq',iq
C   ... Read from disk g for this iq.
        call pshpr(iprint()-10)
        if (iogfrs(10000*ISH+0*MZP+10*MCD+6,1,0,' ',ifi,1,1,0,0,zp,
     .    qp(1,iq),xx,xx,clp,xx,xx,0,0,kcplx,ggii) /= 0)
     .    call rxi('gf-energ-integ failed to read gf for qp no.',iq)
C       Make big GF from small GF
        call gfg2g(s_ham,s_pot,100,lrel,isp,1,nbas,1,nbas,
     .    lidim,lidim,ggii,ghh,ierr)
        call poppr

        call cplxdm(kcplx,lidim*nspc,lidim*nspc,ld11,ld21,ldrgn,ogi)
c        print *, 'iq zp',iq,zp
c       call yprm('ggii',2,ggii,ogi,lidim*nspc,lidim*nspc,lidim*nspc)
        call gfenint1(iq,nkp,wz,lidim,nspc,nsp,isp,ggii,gint)

      enddo

C   ... Save int of GF for all irreducible q-points on disk, spins 1+2
      if (nzpi == nzp) then
        fisnam = 'gfdm1'
        if (isp == 2) fisnam = 'gfdm2'
        fisi = fopnx(fisnam,100,16+8+4+0,-1)
C       Rewind file, write header
        rewind fisi
        call pshpr(1)
        lrd = iogfrs(3,0,0,fisnam,fisi,1,1,nbas,0,zp0,qp,xx,xx,clp,
     .    xx,xx,0,0,0,xx)
        call poppr

        do  iq = 1, nkp
C          call zprm('gf-int',2,gint(1,1,iq,isp),lidim*nspc,
C     .      lidim*nspc,lidim*nspc)

C         Change int-of-GF to kcplx=0
          call gfcnv000(iq,nkp,lidim,isp,nsp,gint,gfe)

C-....close the integration contour in the upper complex-energy plane by
C-....making the result hermitian, i.e. the path is symmetric with respect to
C-....real axis
          call p0hermrl(lidim*nspc,1,1,gfe)
C         call yprm('ggii',2,gfe,ogi,lidim*nspc,lidim*nspc,lidim*nspc)

C    ........write the result (dens. matrix)
          call pshpr(iprint()-10)
          lrd = iogfrs(10000*ISH+7,1,0,' ',fisi,1,1,0,0,zp0,
     .      qp(1,iq),xx,xx,clp,xx,xx,0,0,kcplx,gfe)

          call poppr
        enddo
        call fclose(fisi)
      endif

      deallocate(ggii)
      call tcx('gfenint')

      end

      subroutine gfenint1(iq,nkp,wz,lidim,nspc,nsp,isp,gf,gint)
C-
C ----------------------------------------------------------------------
Ci Inputs
Ci   iq
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   wz    :weights for complex energy integration
Ci   lidim :number of lower+intermediate orbitals
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   gf
Ci   gint
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   29 Sep 04
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lidim,nsp,nspc,nkp,iq,isp
      double precision gf(lidim,nspc,lidim,nspc,2)
      double complex gint(lidim,lidim,nkp,nsp),wz
C ... Local parameters
      double precision pi,rnsp
      double complex gfc,wt
      integer i,j

C ... wt is 2*1/(2*pi) * energy weight
      pi = 4*datan(1d0)
      rnsp = 1d0
      wt = (2d0/rnsp)*wz/pi
C      srm1 = (0d0,1d0)



c      do  10  is2 = 1, nspc
      do  10  j = 1, lidim
c      do  10  is1 = 1, nspc
      do  10  i = 1, lidim
c        gfc = dcmplx(gf(i,is1,j,is2,1),gf(i,is1,j,is2,2))
C         arra = (-1)*gf(i,is1,j,is2,2)
C         gfc = dcmplx(arra,gf(i,is1,j,is2,1))
         gfc = dcmplx(-gf(i,1,j,1,2),gf(i,1,j,1,1))

         gint(i,j,iq,isp) = gint(i,j,iq,isp) + wt*gfc
   10 continue
c      call zprm('gfc-int-1',2,gfc,lds,lds,lds)
c      call yprm('ggiint',2,gf,lidim*lidim*nspc**2,lidim*nspc,
c     .      lidim*nspc,lidim*nspc)

      end

      subroutine gfcnv000(iq,nkp,lidim,isp,nsp,gint,gf)
      implicit none
      integer lidim,nsp,isp
      integer i,j,iq,nkp
ccccc deleted the spin part from gfbz at the end
      double precision gf(lidim,lidim,2)
      double precision gint(2,lidim,lidim,nkp,nsp)

c      print *, 'replace gfcnv000 by standard call'


      do  100  j = 1, lidim
      do  100  i = 1, lidim



          gf(i,j,1) = gint(1,i,j,iq,isp)
          gf(i,j,2) = gint(2,i,j,iq,isp)


 100   continue
c      call yprm('gf-int-1',2,gf,lidim*lidim,lidim,lidim,lidim)

      end
