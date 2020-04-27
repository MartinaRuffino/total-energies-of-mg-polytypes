      subroutine fixplat(opt,plat,tol,ng,g)
C- Adjusts site positions to agree with given symmmetry to machine precision
C ----------------------------------------------------------------------
Ci Inputs:
Ci   opt   :1s digit
Ci         :nonzero => shift any component of plat
Ci         :whose value, or squre, is an integer multiple of 1/12,
Ci         :to the integer multiple
Ci         :10s digit
Ci         :nonzero => shift plat to conform to each of ng point group ops,
Ci                     followed by shift in g
Ci         :See Remarks
Ci   plat  :primitive lattice vectors, in units of alat
Ci   tol   :used when 1s digit is nonzero:
Ci          upper bound to allowed shift in any component of plat
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   istab: atom transformation table; see symtab
Co Inputs/Outputs:
Cio  plat  :primitive lattice vectors, in units of alat
Ci         :plat may be modified
Cr Remarks:
Cr   Lattice vectors may be supplied with limited precision.
Cr   This routine can modify plat in one of 2 ways:
Cr   1s digit opt>0: lattice vectors are near, but not exactly coincident with
Cr      numbers (or square root of numbers) that are integral multiples of 1/12,
Cr      e.g. -0.50000001 or 0.86602539 ~ sqrt(3/4)
Cr   this routine will shift any component within tol of such a multiple.
Cr   10s digit opt>0: attempt to make lattice vectors exactly comply with given
Cr       point group operations.  For each of ng point group operations g
Cr       fixplat will adjust plat to correspond to g to machine precision.
Cr       This routine runs through each g sequentially, finding a shift for each
Cr       g makes, by solving a linear equation (Sylvester equation).
Cr       Note: this does not guarantee that the final plat will conform to all point
Cr       group operations.  The conditions of two separate g's may be incompatible.
Cr       Example:  This lattice
Cr          plat= 1.0 0.0 0.0 -0.50000001 0.86602539 0.0 0.0 0.0 2.96724119
Cr       has 12 space group ops given by generators
Cr          r6z::(0.0000001,0,1/2) my::(0.0000001,0,1/2)
Cr       Both of these two point group operations are slightly out of compliance:
Cr         my
Cr         m(sqrt(3)/2,1/2,0)
Cr       The shift that makes one symop comply, undoes compliance of the other if the
Cr       ops themselves are slightly in error (generated with a slightly misshapen plat)
Cr       The remedy chosen here is to make g conform to plat as a second step.
Cr       A full remedy requires that the g's and plat be adjusted as coupled problem.
Cr       This is a big job, and it has not been attempted.
Cu Updates
Cu   03 Apr 16 First created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer opt,ng
      double precision plat(3,3),g(3,3,ng),tol
C Local parameters:
      integer ig,i,j
      double precision delp(3,3),plat0(3,3),xx,sdlength
      procedure(real(8)) :: dlength
      integer, parameter :: IPRT=60

      call info2(20,1,-1*0,' FIXPLAT: mode=%i, ng=%i',opt,ng)

      call dcopy(9,plat,1,plat0,1)

C --- For each component, seek values that are multiples of 1/12 ---
      if (mod(opt,10) > 0) then
        do  i = 1, 3
        do  j = 1, 3
          xx = (12*plat(i,j))**2
          if (abs(xx-nint(xx)) < tol) then
            xx = nint(xx)
            xx = dsign(dsqrt(xx)/12d0,plat(i,j))
            plat(i,j) = xx
          endif
        enddo
        enddo
      endif

C --- For each g, find and shift by delp ---
      if (mod(opt/10,10) > 0 .and. ng > 1) then
        call info0(IPRT,0,0,'   ig     delta p')
        do  ig = 1, ng
          call pfixplat(plat,g(1,1,ig),delp)
          xx = dlength(9,delp,1)
          call info2(IPRT,0,0,'  %,3i   %:-2,3;3g',ig,xx)
C         call prmx('delp',delp,3,3,3)
          plat = plat + delp
        enddo
      endif
      call info2(20,0,-1,' delta |plat|=%,3;3g',dlength(9,plat-plat0,1),0)

      sdlength = 0
      if (mod(opt/10,10) > 0 .and. ng > 1) then
        call info0(IPRT,1,0,'   ig     delta g')
        do  ig = 1, ng
          call pfixgrp(plat,g(1,1,ig),delp)
          xx = dlength(9,delp,1)
          sdlength = sdlength + xx
          call info2(IPRT,0,0,'  %,3i   %:-2,3;3g',ig,xx)
          g(:,:,ig) = g(:,:,ig) + delp
          do  i = 1, 3
            call info2(IPRT+20,0,0,'%3;14,9D',g(1,i,ig),0)
          enddo
        enddo
        call info2(20,0,-1,'  delta |g|=%,3;3g',sdlength,0)
      endif

      sdlength = 0
      do  ig = 1, ng
        call pfixgrp(plat,g(1,1,ig),delp)
        sdlength = dlength(9,delp,1)
      enddo
      call info2(20,0,0,'  error in updated g=%,3;3g',sdlength,0)

      call info0(IPRT+20,0,0,' Modified plat:')
      do  i = 1, 3
        call info2(IPRT+20,0,0,'%3;14,9D',plat(1,i),0)
      enddo

      end

      subroutine pfixplat(plat,g,delp)
C- Adjust lattice vectors to conform to given symop
C ----------------------------------------------------------------------
Ci Inputs:
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   g     :point group operator
Co Outputs:
Co   delp  :shift in plat
Cr Remarks:
Cr   This routine uses the following algorithm:
Cr   If plat conforms to point group g
Cr      (plat)^-1 g (plat) = pigp0 where pigp0 = some 3x3 integer array
Cr   There will be some deviation from integer.
Cr   Let pigp0 be the integer array closest to (plat)^-1 g (plat)
Cr   The elements of plat are adjusted by dp so that
Cr      (plat+dp)^-1  g (plat+dp) = pigp0    (Eq 1)
Cr   (Eq 1) is equivalent to the Sylvester equation
Cr      g dp - dp pigp0 = -C where C = g plat - plat pigp0
Cr   This routine determines dp from (Eq 1).
Cr   Sylvester equation for solution of AX + XB = -C may be written
Cr     L vec X = -vec C
Cr   where vec transforms a matrix into a vector and
Cr     L = I \otimes A + B^T \otimes I where \otimes means the Kroneker product
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      double precision plat(3,3),g(3,3),delp(3,3)
C ... Local parameters
      integer i
      integer, parameter :: nplat=3
      procedure(real(8)) :: dlength
      double precision platr(3,3),qlat(3,3),vol,abnrm
      double precision pigp(3,3),pigp0(3,3),Lr(9,9),Cr(9)
      double complex evl(9),z(9,9),invz(9,9),rhs(9),zdelp(3,3)

      call mkqlat(plat,qlat,vol)

C     call prmx('g',g,3,3,3)
C     call prmx('qlat',qlat,3,3,nplat)

C ... Integer part pigp0 of plat^-1 g plat
      call dgemm('N','N',3,3,3,1d0,g,3,plat,3,0d0,platr,3)
C     call prmx('rotated plat',platr,3,3,nplat)
      call dgemm('T','N',3,3,3,1d0,qlat,3,platr,3,0d0,pigp,3)
C     call prmx('p^-1 g p',pigp,3,3,nplat)
      pigp0 = nint(pigp)

C ... Make L
C     call prmx('integer part of p^-1 g p',pigp0,3,3,nplat)
      call makdkron(2,3,3,3,3,-transpose(pigp0),[vol],Lr)
C     call prmx('- transpose(int[p^-1 g p] otimes I)',Lr,9,9,9)
      call makdkron(11,3,3,3,3,[vol],g,Lr)
C     call prmx('I otimes g - int[p^-1 g p] otimes I',Lr,9,9,9)

C ... C = g plat - plat pigp0
      call dgemm('N','N',3,3,3,-1d0,plat,3,pigp0,3,0d0,Cr,3)
      call dgemm('N','N',3,3,3,1d0,g,3,plat,3,1d0,Cr,3)
C     call prmx('vec C',Cr,9,9,1)
      if (dlength(9,Cr,1) < 1d-12) then ! Matrix nearly singular
        delp = 0
        return
      endif

C ... Diagonalize L: get eigenvalues eigenvectors and inverse
      call zgeevs('P','N','V','N',9,dcmplx(Lr),9,evl,z,9,z,9,abnrm,i)
C     call zprm('eigenvalues',2,evl,9,9,1)
C     call zprm('eigenvectors z',2,z,9,9,9)
      call zcopy(9*9,z,1,invz,1); call zinv(9,invz,9)
C     call zprm('inverse of z',2,invz,9,9,9)

C ... Make vec delp = z e^-1 z^-1 (-vec C), substituting 0 for small e^-1
C     Eigenvalues of L should be of order 1, or 0 for the three components parallel to rotation axis
      call zgemm('N','N',9,1,9,(-1d0,0d0),invz,9,dcmplx(Cr),9,(0d0,0d0),rhs,9)
      do i = 1, 9
        if (cdabs(evl(i)) > .001d0) then
          rhs(i) = rhs(i) / evl(i)
        else
          rhs(i) = 0
        endif
      enddo
      call zgemm('N','N',9,1,9,(1d0,0d0),z,9,rhs,9,(0d0,0d0),zdelp,9)
      delp = zdelp ! Real part
      if (sum(cdabs(delp-zdelp)) > 1d-6) call rx('bug in pfixplat')
C     call prmx('delp',delp,3,3,3)

      end

      subroutine pfixgrp(plat,g,delg)
C- Adjust lattice vectors to conform to given symop
C ----------------------------------------------------------------------
Ci Inputs:
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   g     :point group operator
Co Outputs:
Co   delg  :shift in g
Cr Remarks:
Cr   This routine uses the following algorithm:
Cr   If plat conforms to point group g
Cr      (plat)^-1 g (plat) = pigp where pigp = some 3x3 integer array
Cr   There will be some deviation from integer.
Cr   Let pigp0 be the integer array closest to pipg.
Cr   The elements of g are adjusted by dg so that
Cr      (plat)^-1  (g+dg) plat = pigp0
Cr      dg = plat (pigp0 - pipg) (plat)^-1
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      double precision plat(3,3),g(3,3),delg(3,3)
      integer, parameter :: nplat=3
      double precision platr(3,3),qlat(3,3),vol
      double precision pigp(3,3),pigp0(3,3)
      procedure(real(8)) :: dlength

      call mkqlat(plat,qlat,vol)

C     call prmx('g',g,3,3,3)
C     call prmx('qlat',qlat,3,3,nplat)

C ... Integer part pigp0 of plat^-1 g plat
      call dgemm('N','N',3,3,3,1d0,g,3,plat,3,0d0,platr,3)
C     call prmx('rotated plat',platr,3,3,nplat)
      call dgemm('T','N',3,3,3,1d0,qlat,3,platr,3,0d0,pigp,3)
C     call prmx('p^-1 g p',pigp,3,3,nplat)
      pigp0 = nint(pigp)
C     print *, dlength(9,pigp-pigp0,1)
      call dgemm('N','N',3,3,3,1d0,plat,3,pigp0-pigp,3,0d0,platr,3)
      call dgemm('N','T',3,3,3,1d0,platr,3,qlat,3,0d0,delg,3)

C ... debug: check that p^-1 (g+dg) p is integer
      call dgemm('N','N',3,3,3,1d0,g+delg,3,plat,3,0d0,platr,3)
C     call prmx('rotated plat',platr,3,3,nplat)
      call dgemm('T','N',3,3,3,1d0,qlat,3,platr,3,0d0,pigp,3)
C     call prmx('p^-1 g p',pigp,3,3,nplat)
      pigp0 = nint(pigp)
C     print *, dlength(9,pigp-pigp0,1)

      end
