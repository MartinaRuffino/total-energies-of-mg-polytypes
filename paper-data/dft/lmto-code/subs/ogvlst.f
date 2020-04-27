      subroutine ogvlst(alat,plat,tau,n1,n2,n3,gmax,ngmx,ng,gv,kv)
C- Set up a list of reciprocal G vectors within cutoff gmax.
C ----------------------------------------------------------------------
Ci Inputs
Ci   job      1s digit
Ci              0 return ng only
Ci              1 return kv and igv
Ci              2 return kv and igv2
Ci              4 return kv and gv
Ci              8 return kv and gv, and sort list
Ci                any combination of 1,2,4 is allowed
Ci   alat,plat  Real-space lattice vectors
Ci   tau        offset subtracted from gv to measure length
Ci   n1,n2,n3   no. divisions along the three lattice vectors
Ci   gmax       Energy cutoff; generate it with gvctof
Ci   ngmx       Leading dimension of gv,kv; must be at least ng
Ci              as calculated by ogvlst; generate it with gvctof
Co Outputs
Co   ng         Number of lattice vectors
Co   gv         list of reciprocal lattice vectors G
Co   kv         indices for gather/scatter operations.
Co              kv(ig,i=1,2,3) for vector ig point to which entry
Co              (i1,i2,i3) vector ig belongs
Cr Remarks
Cr   Collects a list of Fourier entries from a 3D array into a list of
Cr   G-vectors that lie within a cutoff gmax.  List is sorted by length.
Cr   Criterion for cutoff is |(gv-tau)|<gmax
Cr   Use with the following:
Cr     call gvgetf(ng,1,kv,k1,k2,k3,c,c0)
Cr       to collect elements into list c0 from 3D array c
Cr     call gvputf(ng,1,kv,k1,k2,k3,c0,c)
Cr       to poke elements from list c0 into 3D array c
Cr     call gvaddf(ng,kv,k1,k2,k3,c0,c)
Cr       to add elements from list c0 into 3D array c
Cr
Cr   The following routines are designed to work either with a
Cr   specified uniform mesh of G vectors, specified by primitive
Cr   lattice vectors and numbers of divisions n1,n2,n3 (gvcutof,
Cr   gvlist) or without it (gvlst2).
Cr
Cr     gvctof takes as input a uniform mesh of points n1..n3, and
Cr            generates an appropriate energy cutoff gmax, and counts
Cr            the number of G vectors within the cutoff.  The list of
Cr            vectors is generated by looping over all points
Cr            0..n1,0..n2,0..n3, shortens each vector, and retaining
Cr            those for which G<Gmax.  Because only points on the
Cr            specified mesh are considered, there is no guarantee
Cr            that all vectors G<Gmax will be included in the list.
Cr
Cr     gvlist takes as input uniform mesh of G vectors specified by
Cr            n1..n3 and an energy cutoff Gmax, and creates a list of
Cr            G-vectors that lie within Gmax.  ogvlst operates in the
Cr            same way as gvctof, and generates the same list of
Cr            vectors.  The list of vectors is sorted by increasing
Cr            length.  NB: offset tau in ogvlst corresponds to -q in
Cr            gvlst2.
Cr
Cr     gvlst2 is designed to work without any constraint that the list
Cr            of vectors map onto a specified mesh n1..n3.  It takes
Cr            as input an energy cutoff Gmax, and returns a list of
Cr            all G vectors whose length G is G<Gmax.  There is no
Cr            shortening of vectors; the number of divisions n1..n3
Cr            needed to encompass all vectors cannot be specified, but
Cr            is output by gvlst2.  Thus, this routine is only
Cr            suitable in cases where there is no need for the vectors
Cr            to map onto a specified mesh.  NB:  offset q in gvlst2
Cr            corresponds to -tau in ogvlst.
Cr
Cb Bugs
Cb   This is an outdated version of gvlist.  Use gvlist instead.
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   07 Feb 01 changed gmax tolerance to be consistent with gvlst2
Cu   10 Jun 00 added extra argument to gvgetf and gvputf
Cu    2 Sep 98 Adapted from nfp gvlist.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,ng,ngmx,kv(ngmx,3)
      double precision alat,gmax,gv(ngmx,3),plat(3,3),tau(3)
C ... Dynamically allocated local arrays
      real(8), allocatable :: ggv(:)
      real(8), allocatable :: kk(:)
      real(8), allocatable :: iwk(:)
C ... Local parameters
      integer m,ig,j1,j2,j3,iprint
      double precision qlat(3,3),plat1(3,3),qlat1(3,3),g(3),gs(3),pi,
     .  tpiba,tol,gg,vol,gmax2
C ... External calls
      external dinv33,dvshel,pvogls,rx,shorbz

      pi = 4d0*datan(1d0)
      call dinv33(plat,1,qlat,vol)
      tpiba = 2*pi/alat
      tol = 1d-8

C ... Basis vectors for real-space mesh and recip-space supercell
      do  m = 1, 3
        plat1(m,1) = plat(m,1)/n1
        plat1(m,2) = plat(m,2)/n2
        plat1(m,3) = plat(m,3)/n3
        qlat1(m,1) = qlat(m,1)*n1
        qlat1(m,2) = qlat(m,2)*n2
        qlat1(m,3) = qlat(m,3)*n3
      enddo

C ... Loop through g vectors, shorten, keep if within gmax
      ig = 0
      gmax2 = (gmax-tol)**2
      do  j1 = 0, n1-1
      do  j2 = 0, n2-1
      do  j3 = 0, n3-1
        g(1) = j1*qlat(1,1)+j2*qlat(1,2)+j3*qlat(1,3) - tau(1)
        g(2) = j1*qlat(2,1)+j2*qlat(2,2)+j3*qlat(2,3) - tau(2)
        g(3) = j1*qlat(3,1)+j2*qlat(3,2)+j3*qlat(3,3) - tau(3)
        call shorbz(g,gs,qlat1,plat1)
C        print 333, j1,j2,j3,gs(1)+tau(1),gs(2)+tau(2),gs(3)+tau(3)
C  333   format(3i4,3f12.6)
        gg = (tpiba*tpiba)*(gs(1)**2+gs(2)**2+gs(3)**2)
        if (gg <= gmax2) then
C         if (abs(gs(1)+5) < 1d-4 .and. gs(3) < -10.8) then
C           print *, j1,j2,j3
C           print *, gs
C         endif
          ig = ig+1
          if (ig > ngmx) call rx('ogvlst: ng exceeds ngmx')
          gv(ig,1) = gs(1) + tau(1)
          gv(ig,2) = gs(2) + tau(2)
          gv(ig,3) = gs(3) + tau(3)

          kv(ig,1) = j1+1
          kv(ig,2) = j2+1
          kv(ig,3) = j3+1
        endif
      enddo
      enddo
      enddo

      ng = ig
      if (iprint() >= 60) print 861, gmax,ng,n1*n2*n3
  861 format(/' ogvlst: cutoff radius ',f7.3,' gives',i7,
     .  '   recips of max',i7)

C ... Sort the list of vectors
      allocate(ggv(ng),kk(ng),iwk(ng))
      call pvogls(ggv,kk,iwk,tau,ngmx,ng,gv,kv)
      deallocate(ggv,kk,iwk)

      end
      subroutine pvogls(gg,kk,iwk,p,ngmx,ng,gv,kv)
C- Kernel called by ogvlst to sort gv and kv
Ci   p          offset subtracted from g to measure length
      implicit none
      integer ngmx,ng,kk(ng),iwk(ng),kv(ngmx,3)
      double precision gv(ngmx,3),gg(ng),p(3)
C Local variables
      integer ig,m,jg

      do  ig = 1, ng
        gg(ig) =(gv(ig,1)-p(1))**2+(gv(ig,2)-p(2))**2+(gv(ig,3)-p(3))**2
      enddo

      call dvshel(1,ng,gg,iwk,1)
C     call dvheap(1,ng,gg,iwk,0d0,11)

C ... Rearrange gv,kv
      do  m = 1, 3
        do  ig = 1, ng
          jg = iwk(ig)+1
          gg(ig) = gv(jg,m)
          kk(ig) = kv(jg,m)
        enddo
        do  ig = 1, ng
          gv(ig,m) = gg(ig)
          kv(ig,m) = kk(ig)
        enddo
      enddo

C      print 333, p
C  333 format(' ogvlst: g vectors after sorting: p=',3f12.6)
C      do  30  ig = 1, min(ng,50)
C        gg(ig)=(gv(ig,1)-p(1))**2+(gv(ig,2)-p(2))**2+(gv(ig,3)-p(3))**2
CC       print 550, ig,gv(ig,1)-p(1),gv(ig,2)-p(2),gv(ig,3)-p(3),
C        print 550, ig,gv(ig,1),gv(ig,2),gv(ig,3),
C     .     kv(ig,1),kv(ig,2),kv(ig,3),sqrt(gg(ig))
C  550   format(i5,3f11.5,3i6,f11.5)
C   30 continue
C      pause

      end