      subroutine mkorbm(s_site,s_spec,isp,nsp,nspc,nlmax,
     .  ndham,nphimx,nev,wtkp,iq,nbas,ppnl,aus,nl,nkp,orbtm)
C- Orbital moments inside augmentation spheres
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa rmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nlmax :leading dimension of aus
Ci   ndham :dimensions aus,wtkp
Ci   nphimx:dimensions aus:max number of partial waves of a given l channel at any site
Ci   nev   :number of eigenvectors to accumulate orbital moment
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   iq    :current k-point
Ci   nbas  :size of basis
Ci   ppnl  :NMTO potential parameters; see eg potpus.f
Ci   aus   :values of (phi,phidot) MT sphere boundary; see makusq
Ci   nl    :(global maximum l) + 1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Co Outputs
Co   orbtm :orbital moments accumulated for this qp
Cl Local variables
Cl   ispc  :the current spin index in the coupled spins case.
Cl         :Some quantities have no separate address space for each
Cl         :spin in the indepedent-spins case (evec,evl,ewgt) but do
Cl         :in the coupled-spins case.  A separate loop ispc=1..nspc
Cl         :must be added for the latter case
Cl         :ispc is the appropriate index for objects which distinguish
Cl         :spins in the spin-coupled case only
Cl   isp   :isp  is the appropriate index for objects which distinguish
Cl         :spins in the spin-uncoupled case only
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Cl         :ksp is appropriate spin index for quantities that have
Cl         :separate address space for each spin in every case
Cl         :(potential- and density-like objects).
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   16 Jul 09 move transformation to spherical harmonics to makusq
Cu   25 Apr 05 (A. Chantis) extended to local orbitals
Cu   24 Dec 04 Extended to spin-coupled case
Cu   30 Aug 04 (A. Chantis) first written, adapted from mkpdos
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nsp,nspc,nlmax,ndham,nphimx,nbas,nev
      integer n0,nppn,nl,nkp,nab,iq
      parameter (n0=10,nppn=12,nab=9)
      double precision wtkp(ndham*nspc,nsp/nspc,nkp)
      double precision ppnl(nppn,n0,nsp,*),orbtm(nl,nsp,*)
      double complex aus(nlmax,ndham*nspc,nphimx,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer lmxa,lmxax,ib,is,iv,ilm,l,m,
     . ll,lc,em,ispc,ksp
      double precision suml(11),s11,s22,s12,s33,s31,s32,s13,s23,
     .  rmt,sab(nab,n0,2)
      double complex au,as,az,iot

      lmxax = ll(nlmax)
      iot = dcmplx(0d0,1d0)
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        rmt = s_spec(is)%rmt
        lmxa = min(lmxa,lmxax)
        if (lmxa == -1) cycle

C   ... Convert potential parameters to hab,vab,sab
        call phvsfp(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)

C       In noncollinear case, isp=1 always => need internal ispc=1..2
C       ksp is the current spin index in both cases:
C       ksp = isp  in the collinear case
C           = ispc in the noncollinear case
C       ispc=1 for independent spins, and spin index when nspc=2
        do  ispc = 1, nspc
        ksp = max(ispc,isp)
        do  iv = 1, nev
          call dpzero(suml,1+lmxa)
          ilm = 0

C  ....  Rotate from real to spherical harmonics (order assumed: m,...,-m).
C        |Psi>_l = \Sum_{m}(A_l,m * u_l + B_l,m * s_l)*R_l,m --->
C        |Psi>_l = \Sum_{m}(C_l,m * u_l + D_l,m * s_l)*Y_l,m
C        R_l,m and Y_l,m are the real and spherical harmonics respectively.
C
C              | (-1)^m/sqrt(2)*A_l,-m + i*(-1)^m/sqrt(2)*A_l,m , m>0
C        C_l,m=|  A_l,m                                         , m=0
C              |  1/sqrt(2)*A_l,-m -  i*1/sqrt(2)*A_l,m         , m<0
C
C       Same relationships are valid between D and B.

          do  l = 0, lmxa
          lc = (l+1)**2 - l
          do  m = -l, l
            em = abs(m)
            ilm = ilm+1

CC    ...    m,...,-m order
CC            if (m < 0) then
CC            au = 1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) -
CC     .           iot/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
CC            as = 1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) -
CC     .           iot/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
CC            az = 1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) -
CC     .           iot/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
CC            else if (m > 0) then
CC            au = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) +
CC     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
CC            as = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) +
CC     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
CC            az = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) +
CC     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
CC            else
CC    ...   -m,...,m order
C            if (m < 0) then
C            au = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) +
C     .           1d0/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
C            as = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) +
C     .           1d0/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
C            az = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) +
C     .           1d0/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
C            else if (m > 0) then
C          au = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) +
C     .           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
C          as = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) +
C     .           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
C          az = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) +
C     .           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
C            else
C            au = aus(ilm,iv,1,ksp,ib)
C            as = aus(ilm,iv,2,ksp,ib)
C            az = aus(ilm,iv,3,ksp,ib)
C            endif


            au = aus(ilm,iv,1,ksp,ib)
            as = aus(ilm,iv,2,ksp,ib)
            az = 0
            if (nphimx > 2)
     .      az = aus(ilm,iv,3,ksp,ib)


            s11 = dconjg(au)*au*sab(1,l+1,ksp)
            s12 = 2*dconjg(au)*as*sab(2,l+1,ksp)
            s22 = dconjg(as)*as*sab(4,l+1,ksp)
            s33 = dconjg(az)*az*sab(5,l+1,ksp)
            s31 = dconjg(az)*au*sab(6,l+1,ksp)
            s32 = dconjg(az)*as*sab(7,l+1,ksp)
            s13 = dconjg(au)*az*sab(6,l+1,ksp)
            s23 = dconjg(as)*az*sab(7,l+1,ksp)

            orbtm(l+1,ksp,ib) = orbtm(l+1,ksp,ib) +
     .        m*(s11+s12+s22+s33+s32+s23+s31+s13)*wtkp(iv,isp,iq)

C           print *, m,l,iv,ib,sngl(orbtm(l+1,ksp,ib))

          enddo
          enddo

        enddo
        enddo
C        print*, isp, ib, 'ORB.MOMNT=',orbtm(isp,ib)
      enddo
      end
