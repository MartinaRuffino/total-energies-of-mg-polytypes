      subroutine fpiint(nx,np,nxp,x,w)
C- Points and weights for integration on a sphere surface
C ----------------------------------------------------------------------
Ci Inputs
Ci   nx    :number of points in polar angle (Legendre integration)
Ci          Use nx<0 for special points; see Remarks.
Ci   np    :number of points in phi (uniform points)
Ci         :np=0 => np depends on nx making dphi approx constant.
Ci          for nx<0, np is not used.
Co Outputs
Co   nxp   :total number of number of points in quadrature
Co   x     :cartesian coordinates of points on unit sphere
Co   w     :weights corresponding to x
Cr Remarks
Cr   fpiint generates a mesh of points on a unit sphere for angular
Cr   integration, either using a set of special generated from the
Cr   Platonic solids, or by integrating the polar angle with Legendre
Cr   gaussian quadrature and the azimuthal angle with a set of evenly
Cr   spaced points on a circle.
Cr   For special points, invoke fpiint with one of the following:
C     nx= -4 integrates any ylm<=2 exactly (tetrahedron)
C     nx= -6 integrates any ylm<=3 exactly (faces of cube)
C     nx= -8 integrates any ylm<=3 exactly (cube)
C        -12 integrates any ylm<=5 exactly (icosahedron)
C        -20 integrates any ylm<=5 exactly (faces of icosahedron)
C        -30 integrates any ylm<=5 exactly (sides of icosahedron)
C        -60 integrates any ylm<=5 exactly (buckeyball)
C        -32 integrates any ylm<=9 exactly  (combination of 12,20)
C        -62 integrates any ylm<=11 exactly (combination of 12,20,30)
C        -92 integrates any ylm<=11 exactly (combination of 12,20,60)
C       -122 integrates any ylm<=15 exactly (combination of 12,20,30,60)
C ----------------------------------------------------------------------
      implicit none
      integer nx,np,nxp
      double precision x(3,*),w(*),fpi
      integer i,iprint

      if (nx >= 0) then
        call nintsp(nx,np,nxp,x,w)
      else
        fpi = 16*datan(1d0)
        nxp = -nx
        if (nx == -32) then
          call platsl(x,12)
          call platsl(x(1,13),20)
          do  10  i = 1, 12
   10     w(i) = 5d0*fpi/(14*12)
          do  12  i = 13, 32
   12     w(i) = 9d0*fpi/(14*20)
        elseif (nx == -62) then
          call platsl(x,12)
          call platsl(x(1,13),20)
          call platsl(x(1,33),30)
          do  110  i = 1, 12
  110     w(i) = 125d0*fpi/(14*12*33)
          do  112  i = 13, 32
  112     w(i) = 81d0*fpi/(14*20*33)
          do  114  i = 33, 62
  114     w(i) = 256d0*fpi/(14*30*33)
        elseif (nx == -92) then
          call platsl(x,12)
          call platsl(x(1,13),20)
          call platsl(x(1,33),60)
          do  210  i = 1, 12
  210     w(i) = 1/12.34817490904537d0*fpi/12
          do  212  i = 13, 32
  212     w(i) = 2.986997567806883d0/12.34817490904537d0*fpi/20
          do  214  i = 33, 92
  214     w(i) = 8.361177341238484d0/12.34817490904537d0*fpi/60
        elseif (nx == -122) then
          call platsl(x,12)
          call platsl(x(1,13),20)
          call platsl(x(1,33),30)
          call platsl(x(1,63),60)
          do  310  i = 1, 12
  310     w(i) = (0.0939463041645901d0)*fpi/12
          do  312  i = 13, 32
  312     w(i) = (0.2373458837681504d0)*fpi/20
          do  314  i = 33, 92
  314     w(i) = (0.0378880378880377d0)*fpi/30
          do  316  i = 63, 122
  316     w(i) = (0.6308197741792218d0)*fpi/60
        else
          call platsl(x,nxp)
          do  20  i = 1, nxp
   20     w(i) = fpi/nxp
        endif
      endif

C --- Printout ---
      if (iprint() < 80) return
      print '(/'' fpiint:'',i5,'' points generated:'')', nxp
      do  90  i = 1, nxp
   90 print 333, i, x(1,i), x(2,i), x(3,i), w(i)
  333 format(i3,4f20.15)

      end

      subroutine fp2yl(nr,nlm,nsp,np,wp,fp,yl,fac,fl)
C- Add to fl Yl-projection of function fp tabulated on an angular mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlm   :(lmax+1)**2 where lmax is l-cutoff for fl.  2nd dim for fl and yl
Ci   nsp   :number of functions
Ci   np    :fp is held on np angular mesh points
Ci   wp    :quadrature weights
Ci   fp    :representation of function on angular mesh
Ci   yl    :Spherical harmonics
Ci   fac   :Replace fl with (YL projection of fp) + fac*fl
Co Outputs
Co   fl    :project fp as sum_L fl(r) YL,
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   03 Aug 13  Folded into fpiint
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nlm,np,nsp
      double precision fac
      double precision fl(nr,nlm,nsp),fp(nr,np,nsp),yl(np,nlm),wp(np)
C ... Local parameters
      integer ip,ilm,i
C     integer ir
      double precision ylwp(np,nlm)
C     double precision xx

C     Scale yl by wp for fast multiplication
      do  ilm = 1, nlm
        do  ip = 1, np
          ylwp(ip,ilm) = yl(ip,ilm)*wp(ip)
        enddo
      enddo
      do  i = 1, nsp
        call dgemm('N','N',nr,nlm,np,1d0,fp(1,1,i),nr,ylwp,np,fac,
     .    fl(1,1,i),nr)
      enddo

C      do  i = 1, nsp
C      do  ip = 1, np
C      do  ilm = 1, nlm
C        xx = wp(ip)*yl(ip,ilm)
C        do  ir = 1, nr
C          fl(ir,ilm,i) = fl(ir,ilm,i) + fp(ir,ip,i)*xx
C        enddo
C      enddo
C      enddo
C      enddo
      end

      subroutine fyl2p(nr,nlm,nsp,np,yl,fl,fp)
C- Tabulate function represented in Yl expansion on an angular mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlm   :(lmax+1)**2 where lmax is l-cutoff for fl
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   np    :fp is held on np angular mesh points
Ci   yl    :Spherical harmonics
Ci   fl    :given function represented as sum_L fl(r) YL
Co Outputs
Co   fp    :tabulation of fl on angular mesh
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   03 Aug 13  Folded into fpiint
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nlm,np,nsp
      double precision fl(nr,nlm,nsp),fp(nr,np,nsp),yl(np,nlm)
C ... Local parameters
      integer i

      do  i = 1, nsp
        call dgemm('N','T',nr,np,nlm,1d0,fl(1,1,i),nr,yl,np,0d0,
     .    fp(1,1,i),nr)
      enddo
      end

C      subroutine fmain
CC- Tests fpiint
C      implicit none
C      logical cmdstr
C      character*(160) first,fnam
C      double precision pi,srfpi,sumd
C      integer nnn,lmax,nxl(0:7),nth,nph,np,nlm,ilm,ifi,fopng,i,rdm,
C     .  ip,jlm
C      parameter(nnn=122)
C      double precision p(3,nnn),wp(nnn),p2(nnn,3),r2(nnn)
C      real(8), allocatable :: yl(:,:),fl(:),fr(:,:),fbak(:)
C      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/
C
C      print *, 'check platsl for 4 points'
C      call platsl(p,4)
C      do  ip = 1, 4
C        print '(i4,3f15.10)', ip, p(:,ip)
C      enddo
C      print *, 'check platsl for 6 points'
C      call platsl(p,6)
C      do  ip = 1, 6
C        print '(i4,3f15.10)', ip, p(:,ip)
C      enddo
C      print *, 'check platsl for 8 points'
C      call platsl(p,8)
C      do  ip = 1, 8
C        print '(i4,3f15.10)', ip, p(:,ip)
C      enddo
C
C      call info0(0,1,1,' Decompose function tabulated on an angular'//
C     .  'mesh into YL expansion.')
C
CC      call setpr(81)
C
C      pi = 4d0*datan(1d0)
C      srfpi = dsqrt(4*pi)
C      lmax = 7
C      nlm = (lmax+1)**2
C      nth = nxl(min(lmax,7))
C      nth = nxl(7)
C      nph = 0
C      call fpiint(nth,nph,np,p,wp)
C      call dmcpy(p,1,3,p2,np,1,np,3)
CC     call prmx('p',p2,np,np,3)
C      allocate (yl(np,nlm),fr(np,4),fl(nlm),fbak(np))
C      call ropyln(np,p2,p2(1,2),p2(1,3),lmax,np,yl,r2)
C
CC ... Numerically integrate YL(jlm)
C      jlm = 7
C      call info2(0,0,0,' Test 1. Decompose YL,L=%i.'//
C     .  '%N   ilm%5ffl',jlm,0)
C      fl = 0
C      call fp2yl(1,nlm,1,np,wp,yl(1,jlm),yl,0d0,fl)
C      do  ilm = 1, nlm
C        if (abs(fl(ilm)) > 1d-15) then
C          call info2(0,0,0,'%,6i   %,15;15g',ilm,fl(ilm))
C        endif
C      enddo
C
CC ... Numerically integrate YL(jlm)
C      if (.not. cmdstr(1,first)) goto 99
C      if (first(1:4) /= '-fn=') goto 99
C      fnam = first(5:)
C      ifi = fopng(fnam,-1,0)
C      call info0(0,1,0,
C     .  '#Test 2 ... attempting to read (122,4) array from file '//
C     .  trim(fnam)//' ...')
C      i = rdm(ifi,0,np*4,first,fr,122,4)
C      if (i /= 1) call rx('failed to read data from file')
C      call info0(0,0,0,
C     .  '#%6f ... checking that cols 1..3 match mesh points')
C      do  ip = 1, np
C        if (abs(p2(ip,1)-fr(ip,1))+
C     .      abs(p2(ip,2)-fr(ip,2))+
C     .      abs(p2(ip,3)-fr(ip,3)) > 1d-10)
C     .    call rx('fr does not match mesh')
C      enddo
C      call info0(0,0,0,
C     .  '#%6f ... doing numerical integration'//
C     .  '%N#   ilm%5ffl')
C
CC     fr(:,4) = yl(:,11) ! replace f with yl(11) to show exact
C      call fp2yl(1,nlm,1,np,wp,fr(1,4),yl,0d0,fl)
C      do  ilm = 1, min(nlm,16)
C        if (abs(fl(ilm)) > 1d-15) then
C          call info2(0,0,0,'%,6i   %,10;10g',ilm,fl(ilm))
C        endif
C      enddo
C      fbak = 0
C      sumd = 0
C      call fyl2p(1,nlm,1,np,yl,fl,fbak)
C      print *
C      print *, 'Compare against original'
C      print *, ' pt     original      integrated      diff'
C      do  ip = 1, 122
C        write(6,101) ip, fr(ip,4),fbak(ip), fbak(ip)-fr(ip,4)
C  101   format(i4,3f15.10)
C        sumd = wp(ip)*(fbak(ip)-fr(ip,4))
C      enddo
C      call info2(0,0,0,' integral of diff = %,10;10g',sumd,0)
C
C      return
C
C   99 continue
CC     return
C      print *, 'Projects function tabulated on 122 special points',
C     .  ' in a YL expansion.'
C      print *, 'Usage:  fpiint -fn=filename'
C
C      end
