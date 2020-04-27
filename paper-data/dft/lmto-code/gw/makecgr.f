      subroutine makecgr(s_lat,lmax,symgr,cgr)
C- Order and possibly rotate Clebsch-Gordan coefficients <lm3 | lm1 lm2>
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  indxcg jcg cg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   lmax :l-cutoff
Ci   symgr :rotation matrix.  If symgr(1,1) = NULLI, no rotation
Co Outputs
Co   cgr   :Clebsch Gordon coefficients in matrix form, possibly rotated by symgr
Cr Remarks
Cr
Cu Updates
Cu   20 Aug 18 Adapted from GW code rotcg
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lmax
      real(8):: cgr((lmax+1)**2,(lmax+1)**2,(2*lmax+1)**2)
      real(8) :: symgr(3,3)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer nlmax,ilma,la,ilmb,lh,ii,indx,icg,icg1,icg2,ilm,lm,l,m,lx,lm2,lm1,lcut,nlcut,ldr
      real(8) dlmm(-2*lmax:2*lmax,-2*lmax:2*lmax,0:2*lmax)
      real(8) cgloc((lmax+1)**2,(lmax+1)**2,(2*lmax+1)**2)
      real(8) rmat((2*lmax+1)**2,(2*lmax+1)**2)
C#ifdefC COMMONLL
C      integer(4)::ll(51**2)
C      common/llblock/ll
C#else
      procedure(integer) :: ll
C#endif

C     call dbgcg(s_lat,lmax)


      ldr = (2*lmax+1)**2
      nlmax = (lmax+1)**2

C --- Arrange CG in matrix order ---
      call dpzero(cgloc,size(cgloc))
      do ilma = 1, nlmax
        la = ll(ilma)
        do ilmb = 1, nlmax
          lh = ll(ilmb)
          ii = max0(ilma,ilmb)
          indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
          icg1 = s_lat%indxcg(indx)
          icg2 = s_lat%indxcg(indx+1)-1
          do  icg = icg1, icg2
            ilm = s_lat%jcg(icg)
            cgloc(ilma,ilmb,ilm) = s_lat%cg(icg)
          enddo
        enddo
      enddo
C      call prmx('cg(unrotated)',cgloc,nlmax**2,nlmax**2,ldr)

C     Special case symgr is a unit matrix: nothing to rotate
      if (sum(abs(symgr)) == 3 .and. symgr(1,1)==1 .and. symgr(2,2)==1 .and. symgr(3,3)==1) then
        call dcopy(size(cgloc),cgloc,1,cgr,1)
        return
      endif

CC --- Rotation matrix, matrix product by l ---
C      call rotdlmm(symgr,1,2*lmax+1,dlmm)
C
CC --- Rotate Clebsch Gordan coefficients ---
C      do  lm = 1, ldr
C        l = ll(lm)
C        m = lm - l**2 - l - 1
C        lx = l**2 + l + 1
C        do  lm2 = 1, nlmax
C          do  lm1 = 1, nlmax
C            cgr(lm1,lm2,lm) = sum(cgloc(lm1,lm2,lx-l:lx+l)*dlmm(-l:l,m,l))
C          enddo
C        enddo
C      enddo

C --- Rotate lm3.  Divide into blocks to reduce dgemm time ---
      call ylmrtg(ldr,symgr,rmat)  ! Rotates Y_L
C     call prmx('rmat',rmat,ldr,ldr,ldr)
      call dpzero(cgr,size(cgr))
      nlcut = min(25,ldr)
      call dgemm('N','N',nlmax**2,nlcut,nlcut,1d0,cgloc,nlmax**2,rmat,ldr,0d0,cgr,nlmax**2)
      if (nlcut < ldr) then
        call dgemm('N','N',nlmax**2,ldr-nlcut,ldr-nlcut,1d0,cgloc(1,1,nlcut+1),nlmax**2,
     .    rmat(nlcut+1,nlcut+1),ldr,0d0,cgr(1,1,nlcut+1),nlmax**2)
      endif

C      call prmx('cgr',cgr,nlmax**2,nlmax**2,ldr)

      end
      subroutine dbgcg(s_lat,lmax)
C- Debugging routines for CG coefficients
      use structures
      implicit none
C ... Passed parameters
      integer lmax
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      real(8), allocatable :: ylp(:,:),ylwp(:,:),ylylp(:),ylyll(:,:,:)
      real(8), allocatable :: gylp(:,:,:),gylyll(:,:,:,:)
      real(8), allocatable :: cgm(:,:,:)
C ... Local parameters
      logical :: lgrad = .true.
      integer, parameter :: nnn = 122
      integer np,nth,nlm,nlm2,nlmg,ip,ilm,jlm,icg,icg1,icg2,ii,ilma,ilmb,indx,la,lh
      double precision p(3,nnn),wp(nnn),p2(nnn,3),r2(nnn),fpi
      procedure(integer) :: ll

      nlm = (lmax+1)**2
      nlm2 = (2*lmax+1)**2
      nlmg = nlm2
      if (lgrad) nlmg = (2*lmax+2)**2      ! Only (2*lmax+1)**2 if grad not needed
      nth = -nnn
      fpi = 16d0*datan(1d0)
      call fpiint(nth,0,np,p,wp)
C ... Spherical harmonics (unit radius) ... to be overwritten by Ylp * wp
      call dmcpy(p,1,3,p2,np,1,np,3)
      allocate (ylp(np,nlmg),ylwp(np,nlmg))
      call ropyln(np,p2,p2(1,2),p2(1,3),ll(nlmg),np,ylp,r2)
      call dcopy(size(ylp),ylp,1,ylwp,1) ! Not needed?
      if (lgrad) then
        allocate (gylp(np,nlm2,3))
        call ropylg(1,ll(nlmg)-1,nlm2,np,np,p2,p2(1,2),p2(1,3),r2,ylp,gylp)
      endif

C ... Resolve product YLa YLb into YM by numerical integration on the sphere
      allocate (ylylp(np),ylyll(nlm2,nlm,nlm))
      if (lgrad) allocate (gylyll(nlm2,nlm,nlm,3))
      do  ilma = 1, nlm
      do  ilmb = 1, nlm
        forall (ip=1:np) ylylp(ip) = ylp(ip,ilma)*ylp(ip,ilmb)
        call fp2yl(1,nlm2,1,np,wp,ylylp,ylp,0d0,ylyll(1,ilma,ilmb))
        if (lgrad) then
          do  ii = 1, 3
            forall (ip=1:np) ylylp(ip) = ylp(ip,ilma)*gylp(ip,ilmb,ii)
            call fp2yl(1,nlm2,1,np,wp,ylylp,ylp,0d0,gylyll(1,ilma,ilmb,ii))
          enddo
        endif
      enddo
      enddo

C --- Arrange CG in matrix order, cg(K,L,M) ---
C     Expansion theorem: Y_L Y_M = sum_K cg(K,L,M) Y_K
      allocate(cgm(nlm2,nlm,nlm))
      call dpzero(cgm,size(cgm))
      do  ilma = 1, nlm
        la = ll(ilma)
        do ilmb = 1, nlm
          lh = ll(ilmb)
          ii = max0(ilma,ilmb)
          indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
          icg1 = s_lat%indxcg(indx)
          icg2 = s_lat%indxcg(indx+1)-1
          do  icg = icg1, icg2
            ilm = s_lat%jcg(icg)
            cgm(ilm,ilma,ilmb) = s_lat%cg(icg)
          enddo
        enddo
      enddo

C --- CG for expansion Y_L grad Y_M ---
C     grad f(r) Y_M = C_N fbar(r) Y_M has two kinds of contributions:
C     Coefficients C_N returned by ylmbyr in compact form as cgcof(:,:,1:3,M) for x,y,z
C     1. ll(N) = ll(M)+1: fbar(r) = (df/dr - l*f/r)  = g1(r) and C_N given by cgcof(:,1,1:3,M)
C     1. ll(N) = ll(M)-1: fbar(r) = (df/dr + (l+1)*f/r) = g2(r); C_N given by cgcof(:,2,1:3,M)
C     Divide grad f(r) Y_M into two terms:
C       grad f(r) Y_M = g1(r) sum_N C_N Y_N  +  g2(r) sum_N C_N Y_N
C       where index N described in ylmbyr.
C     Expansion theorem:
C       Y_L grad f(r) Y_M = Y_L grad f(r) Y_M = Y_L [g1(r) sum_N C_N Y_N + g2(r) sum_N C_N Y_N]
C     Each sum over N runs over 2 values, which can be determined from M and p=1,2,3 for x,y,z
C     Write N as N1(1:2,p,M) for ll(N) = ll(M)+1 and N2(1:2,p,M) for ll(N) = ll(M)-1
C       Y_L grad f(r) Y_M = g1(r) cgcof(1:2,1,p,M) Y_L Y_N1(1:2,p,M) + g2(r) cgcof(1:2,2,p,M) Y_L Y_N2(1:2,p,M)
C                         = sum_K cg(K,L,N1(1:2,p,M)) g1(r) cgcof(1:2,1,p,M) Y_K +
C                         = sum_K cg(K,L,N2(1:2,p,M)) g2(r) cgcof(1:2,2,p,M) Y_K
C     Write gcgm(K,L,M,1,p) = cg(K,L,N1(1:2,p,M)) g1(r) cgcof(1:2,1,p,M)
C     and   gcgm(K,L,M,2,p) = cg(K,L,N2(1:2,p,M)) g1(r) cgcof(1:2,2,p,M)

      do  ilma = 1, nlm
      do  ilmb = 1, nlm
        do  ilm = 1, min(nlm2,((15-lmax*2)+1)**2)  ! num int exact only up to l=15
          if (abs(cgm(ilm,ilma,ilmb)) + abs(ylyll(ilm,ilma,ilmb)) > 1d-10) then
            print 333, ilm, ilma, ilmb, cgm(ilm,ilma,ilmb), ylyll(ilm,ilma,ilmb), cgm(ilm,ilma,ilmb)-ylyll(ilm,ilma,ilmb)
  333       format(3i4,3f14.8)
          endif
        enddo
      enddo
      enddo

      stop


      end
