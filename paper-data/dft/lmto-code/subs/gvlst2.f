      subroutine gvlst2(alat,plat,q,n1,n2,n3,gmin,gmax,mshlst,job,ngmx,ng,kv,gv,igv,igv2)
C- Set up a list of recip vectors within cutoff |q+G| < gmax
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat     Lattice constant
Ci   plat     Real-space primitive lattice vectors
Ci   q        vectors |q+G|<gmax are included in list
Ci   mshlst   if first entry is nonzero, a list of allowed values
Ci            n1,n2,n3 may take.  First entry is the size of
Ci            the list; then follows the list itself.
Ci   job      1s digit
Ci              0 return ng only
Ci              1 return kv and igv
Ci              2 return kv and igv2
Ci              4 return kv and gv
Ci              8 return kv and gv, and sort list
Ci                any combination of 1,2,4 is allowed
Ci            100s digit
Ci              1 to return internally generated values for n1,n2,n3,
Ci                as determined by call to gvlstn.
Ci                See description of n1,n2,n3 below.
Ci              2 to return internally generated values for n1,n2,n3,
Ci                as determined by maximum values found in list
Ci            1000s digit
Ci              1  purge list of duplicate (kv1,kv2,kv3) (if 1st digit job > 0)
Ci              2  not used ... OLD Do not change input gmax if it is nonzero.
Ci   ngmx     Leading dimension of kv,gv,igv
Ci   gmin     Lower bound for reciprocal lattice vectors, in a.u.
Cio Inputs/Outputs
Cio   gmax    On input, cutoff for reciprocal lattice vectors, in a.u.
Cio           Energy cutoff is gmax**2.
Cio           If input gmax is zero, gvlst2 will generate it from n1..n3
Cio           (It is an error for both gmax and n1..n3 to be zero.)
Cio           On output, gmax may be altered; see Remarks.
Cio  n1..3    On input, max # divisions along the three lattice vectors.
Cio           (It is an error for both gmax and n1..n3 to be zero.)
Cio           Otherwise, input n1..n3 additionally constrain which
Cio           vectors are added to list; see Remarks.
Cio           On output, any n1..n3 initially zero are found from gmax
Cio           n1..n3 are only returned if 100s digit of job is set
Co Outputs
Co   ng       Number of lattice vectors
Co   gv       list of reciprocal lattice vectors G
Co   igv      list of reciprocal lattice vectors G, represented as
Co            three integers (the multiples of qlat)
Co            gv and igv are related by:
Co              gv(1:3,1:ng) = 2*pi/alat * (qlat * igv(1:ng))
Co   igv2     same as igv except first and second columns are permuted
Co   kv       indices for gather/scatter operations.
Co            kv(ig,i=1,2,3) for vector ig point to which entry
Co            (i1,i2,i3) vector ig belongs
Cr Remarks
Cr   Collects a list of q + reciprocal lattice vectors (G+q) that lie
Cr   within a cutoff gmax.  List is optionally sorted by length.
Cr   Vectors G are integer multiples of primitive r.l.v. qlat.
Cr
Cr   Cutoff gmax may be input (preferred mode), or generated from
Cr   input values n1..n3.
Cr
Cr   Additional constraints may be imposed, namely that the number of
Cr   multiples in each of axes 1..3 not exceed input values n1..n3.
Cr
Cr   For copying the list of points to/from a regular mesh, use:
Cr     call gvgetf(ng,1,kv,k1,k2,k3,c,c0)
Cr       to collect elements into list c0 from 3D array c
Cr     call gvputf(ng,1,kv,k1,k2,k3,c0,c)
Cr       to poke elements from list c0 into 3D array c
Cr     call gvaddf(ng,kv,k1,k2,k3,c0,c)
Cr       to add elements from list c0 into 3D array c
Cr
Cr   The following routines are designed to work either with a
Cr   specified uniform mesh of G vectors, specified by primitive
Cr   lattice vectors and numbers of divisions n1,n2,n3 (gvctof,
Cr   gvlist) or without it (gvlst2).
Cr
Cr     gvctof takes as input a uniform mesh of points n1..n3. gvctof
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
Cr            G-vectors that lie within Gmax.  gvlist operates in the
Cr            same way as gvctof, and generates the same list of
Cr            vectors.  The list of vectors is sorted by increasing
Cr            length.
Cr
Cr     gvlst2 is designed to work without any constraint that the list
Cr            of vectors map onto a specified mesh n1..n3.  It takes
Cr            as input an energy cutoff Gmax, and returns a list of
Cr            all G vectors whose length G is G<Gmax.  There is no
Cr            shortening of vectors; the number of divisions n1..n3
Cr            needed to encompass all vectors cannot be specified, but
Cr            is output by gvlst2.  Thus, this routine is only
Cr            suitable in cases where there is no need for the vectors
Cr            to map onto a specified (n1..n3).  Offset q in gvlst2
Cr            corresponds to -tau in gvlist.
Cr
Cu Updates
Cu   03 Mar 13 Modified 100s digit job
Cu   11 Feb 13 Rearranged to better handle shortening of G vectors
Cu             (modified with sugvec.f, gvlst2.f)
Cu   11 Jul 08 New argument gmin
Cu   01 Jun 01 revised gvlst2 together with gvlist.  They form similar
Cu             operations but with different functions; see gvlist.f
Cu   26 Mar 01 Another bug fix for input n1..n3 and gmax nonzero
Cu   06 Mar 01 Bug fix for input n1..n3 and gmax nonzero
Cu   Routine was adapted from T. Kotani, routine getgv2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,job,ng,ngmx,kv(ngmx,3),igv(ngmx,3),igv2(3,*),
     .  mshlst(0:*)
      double precision alat,gmin,gmax,gv(ngmx,3),plat(3,3),q(3)
C ... Dynamically allocated arrays
C ... Local parameters
      integer ig,n1max,n1min,n2max,n2min,n3max,n3min,i1,i2,i3,nn,i
      integer n1l,n2l,n3l,n1x,n2x,n3x
      integer PRTG,PRTG2,iset(3),stdo,ipr,lgunit,job0,job1,job2,job4,
     .  job8,k1,k2,k3
      double precision qlat(3,3),vol,pi,tpiba,qpg(3),q2
      double precision gmin0,gmax0,gmin2,gmax2,h1,h2,h3,ddot,tol
      character*256 outs
      parameter (PRTG=30,PRTG2=70,tol=1d-8)

C     call prmx('gv',gv,ngmx,ng,3)
C ... Setup
      stdo = lgunit(1)
      call getpr(ipr)
      call dinv33(plat,1,qlat,vol)
      pi = 4d0*datan(1d0)
      tpiba = 2*pi/alat
C     Dimensionless gmax
      gmin0  = gmin/tpiba
      gmax0  = gmax/tpiba
      if (gmin < 0) call rx('gvlst2: input gmin <= 0')
      if (gmax <= 0) call rx('gvlst2: input gmax <= 0')
      job0 = mod(job,100)
      job1 = mod(job0,2)
      job2 = mod(job0/2,2)
      job4 = mod(job0/4,2)
      job8 = mod(job0/8,2)
      if (job8 /= 0) job4 = 1
C     Save n1..n3 for printout
      iset(1) = n1
      iset(2) = n2
      iset(3) = n3

C --- Compute range n1min,n1max for each direction vector n1,n2,n3 ---
      nn = n1
      call gvlstn(qlat(1,1),qlat(1,2),qlat(1,3),q,mshlst,gmax0,nn,n1min,n1max)
      n1l = nn
      if (mod(job,1000) >= 100) n1 = nn
      nn = n2
      call gvlstn(qlat(1,2),qlat(1,3),qlat(1,1),q,mshlst,gmax0,nn,n2min,n2max)
      n2l = nn
      if (mod(job,1000) >= 100) n2 = nn
      nn = n3
      call gvlstn(qlat(1,3),qlat(1,1),qlat(1,2),q,mshlst,gmax0,nn,n3min,n3max)
      n3l = nn
      if (mod(job,1000) >= 100) n3 = nn

C --- Count number, get list of vectors ---
      ig = 0
      n1x = 0; n2x = 0; n3x = 0;
      gmax2 = (gmax0-tol)**2
      gmin2 = gmin0**2
      do i1 = n1min, n1max
      if (job0 /= 0) then
        k1 = mod(i1+n1l,n1l)+1
      endif
      do i2 = n2min, n2max
      if (job0 /= 0) then
        k2 = mod(i2+n2l,n2l)+1
      endif
      do i3 = n3min, n3max
      if (job0 /= 0) then
        k3 = mod(i3+n3l,n3l)+1
      endif
        do  i = 1, 3
          qpg(i)= q(i) + qlat(i,1)*i1 + qlat(i,2)*i2 + qlat(i,3)*i3
        enddo

C        if (abs(qpg(1)+5) < 1d-4 .and. qpg(3) < -10.8) then
C          print *, i1,i2,i3,k1,k2,k3
C          print *, qpg
C        endif

        q2 = qpg(1)**2+qpg(2)**2+qpg(3)**2
        if (q2 >= gmin2 .and. q2 < gmax2) then
          ig = ig+1
C          print 333, ig,i1,i2,i3,qpg,qpg(1)**2+qpg(2)**2+qpg(3)**2-gmax2
C  333     format(i6,3i4,3f15.9,1pe10.2)
          n1x = max(n1x,iabs(i1))
          n2x = max(n2x,iabs(i2))
          n3x = max(n3x,iabs(i3))
          if (job0 /= 0) then
            kv(ig,1) = k1
            kv(ig,2) = k2
            kv(ig,3) = k3
            if (job1 /= 0) then
              igv(ig,1) = i1
              igv(ig,2) = i2
              igv(ig,3) = i3
            endif
            if (job2 /= 0) then
              igv2(1,ig) = i1
              igv2(2,ig) = i2
              igv2(3,ig) = i3
            endif
            if (job4 /= 0) then
              gv(ig,1) = qpg(1)
              gv(ig,2) = qpg(2)
              gv(ig,3) = qpg(3)
            endif
          endif
        endif
      enddo
      enddo
      enddo
      ng = ig

C --- Printout ---
      if (ipr >= PRTG .and. n1l*n2l*n3l == 0) then
        call awrit2('%x%N GVLST2: gmax = %,1;3d a.u. created %i'//
     .    ' recip. lattice vectors',outs,len(outs),stdo,gmax0*tpiba,ng)
      elseif (ipr >= PRTG) then
        h1 = alat*sqrt(ddot(3,plat(1,1),1,plat(1,1),1))/n1l
        h2 = alat*sqrt(ddot(3,plat(1,2),1,plat(1,2),1))/n2l
        h3 = alat*sqrt(ddot(3,plat(1,3),1,plat(1,3),1))/n3l
        call awrit4('%x%N GVLST2: gmax = %,1;3d'//
     .    ' a.u. created %i vectors of %i (%i%%)',
     .    outs,len(outs),0,gmax,ng,n1l*n2l*n3l,(ng*100)/(n1l*n2l*n3l))
        i = iset(1)*iset(2)*iset(3)
        call awrit7('%a%N%9f%?#n#(input) ##mesh has %i x %i x %i'//
     .    ' divisions; length %,3;3d, %,3;3d, %,3;3d',outs,
     .    len(outs),0,i,n1l,n2l,n3l,h1,h2,h3)
        if (i == 0 .and. iset(1)**2+iset(2)**2+iset(3)**2 /= 0) then
          call awrit3('%a%N%9fgenerated from input mesh with ('//
     .      '%?#n#%-1j%i#*#,%?#n#%-1j%i#*#,'//
     .      '%?#n#%-1j%i#*#) divisions',
     .      outs,len(outs),0,iset(1),iset(2),iset(3))
        endif
        call awrit0('%a',outs,len(outs),-stdo)
      endif

C --- Purge duplicate vectors ---
      if (job0 /= 0 .and. mod(job/1000,2) == 1) then
C        allocate(iprm(ng))
C        call ivdup(3,ng,kv,iprm,11)
C        deallocate(iprm)
        call rx('gvlst2: not ready for 1000s digit job')
      endif

C --- Sort the list of vectors --
      if (job8 /= 0) then
        call dpzero(qpg,3)  ! a work array here
        call gvlsts(job1,job2,ngmx,ng,qpg,gv,kv,igv,igv2)
C        if (job1 /= 0) then
C          call prmx('gv',gv,ng,ng,3)
C          call dgemm('N','T',ng,3,3,1d0,dble(igv(1:ng,1:3)),ng,qlat,3,
C     .      0d0,gv,ng)
C          call prmx('gvx',gv,ng,ng,3)
C        endif
      endif
      if (job1+job2+job4 /= 0 .and. ipr >= PRTG2) then
        call dpzero(qpg,3)  ! a work array here
        call gvlstp(job1,job2,alat,qlat,ngmx,ng,q,qpg,gv,kv,igv,igv2)
      endif

C      if (ipr >= PRTG2 .and. ng > 0 .and. job1+job2 /= 0) then
C        write(stdo,333)
C  333   format(' G vectors (multiples of reciprocal lattice vectors)'/
C     .     '   ig    G1   G2   G3     E')
C        do  ig = 1, ng
C          if (job1 /= 0) then
C            i1 = igv(ig,1)
C            i2 = igv(ig,2)
C            i3 = igv(ig,3)
C          endif
C          if (job2 /= 0) then
C            i1 = igv2(1,ig)
C            i2 = igv2(2,ig)
C            i3 = igv2(3,ig)
C          endif
C          do  i = 1, 3
C            qpg(i)= q(i) + qlat(i,1)*i1 + qlat(i,2)*i2 + qlat(i,3)*i3
C          enddo
C          q2 = (qpg(1)**2+qpg(2)**2+qpg(3)**2)*tpiba**2
C          write(stdo,334)  ig,i1,i2,i3,q2
C  334     format(i5,1x,3i5,2x,f8.4)
C        enddo
C      endif

C     call prmx('gv',gv,ngmx,ng,3)
C     print *, 'ng,ngmx=',ng,ngmx
      if (mod(job,1000) >= 200) then
        n1 = n1x
        n2 = n2x
        n3 = n3x
      endif

      end

      subroutine gvlstn(q0,q1,q2,qp,mshlst,gmax0,nn,nmin,nmax)
C- Multiples of r.l.v. that bound cutoff gmax0
C ----------------------------------------------------------------------
Ci Inputs
Ci   q0    :r.l.v. for which bounds nmin and nmax are to be computed
Ci   q1    :first  r.l.v. different from q0
Ci   q2    :second r.l.v. different from q0
Ci   qp    :k-point added to G vectors: sphere is centered G=qp.
Ci   mshlst:An ordered list of integers used to restrict the assignment
Ci         :of nn to one of an allowed list of points, should a value
Ci         :of nn be assigned (see description for nn below).
Ci         :If first entry is nonzero, mshlst(1..) = list of allowed
Ci         :values nn may take.  mshlst(0) is the size of mshlst.
Ci         :Certain fourier transform routines have restrictions
Ci         :on the allowed mesh sizes; this constraint is designed
Ci         :to handle that restriction.
Ci Inputs/Outputs
Ci   nn    :On input, if nonzero, maximum number of mesh points
Ci         :allowed.  Thus nmax <= (nn-1)/2  and  -nmin <= (nn-1)/2
Ci         :On output, if input is zero, nn = 2*max(|nmin|,nmax)+1
Cio  gmax0 :On input, cutoff G
Cio        :On output, gmax0 may be reduced because constraints
Cio        :on nmin,nmax cause lattice vectors to be neglected
Cio        :that are smaller than input gmax0 (if input nn nonzero).
Co Outputs
Co   nmin  :search for lattice vectors limited to (nmin..nmax)*q0
Co   nmax  :search for lattice vectors limited to (nmin..nmax)*q0
Cr Remarks
Cr   q0,q1,q2,qp and G are all dimensionless (units of 2*pi/a)
Cu Updates
Cu   15 Apr 05 Bug fix when nn > max mshlst
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nmin,nmax,nn,mshlst(0:*)
      double precision q0(3),q1(3),q2(3),qp(3),gmax0
C ... Local parameters
      integer indx
      double precision qperp(3),ddot,qqperp

C ... qperp = q1 x q2 / |q1 x q2| ; qqperp = q . qperp
      qperp(1)  = q1(2)*q2(3) - q1(3)*q2(2)
      qperp(2)  = q1(3)*q2(1) - q1(1)*q2(3)
      qperp(3)  = q1(1)*q2(2) - q1(2)*q2(1)
C     anorm = sqrt(ddot(3,q1,1,q1,1))
C     bnorm = sqrt(ddot(3,q2,1,q2,1))
C     call dscal(3,1/anorm/bnorm,qperp,1)
      call dscal(3,1/sqrt(ddot(3,qperp,1,qperp,1)),qperp,1)
      qqperp = ddot(3,q0,1,qperp,1)

C   10 continue
C     print *, gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp + 1
C     print *, -gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp - 1

      nmax =  gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp + 1
      nmin = -gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp - 1

C ... Assign nn, if input value is zero
      if (nn == 0) then
        nn = 2*max(iabs(nmin),nmax)+1
        if (mshlst(0) /= 0) then
          indx = 1
          call hunti(mshlst(1),mshlst,nn,0,indx)
          nn = mshlst(min(indx+1,mshlst(0)))
        endif
      endif

      end

      subroutine gvlsts(job1,job2,ngmx,ng,t,gv,kv,igv,igv2)
C- Kernel called by gvlist,gvlst2 to sort gv,kv,igv,igv2
C ----------------------------------------------------------------------
Ci Inputs
Ci   job1  :sort igv also
Ci   job2  :sort igv2 also
Ci   ngmx  :dimensions gv,kv,igv
Ci   ng    :number of G-vectors
Ci   t     :offset added to gv = G+q to measure length
Cio Inputs/Outputs
Cio        :the following arrays are sorted:
Cio  gv    :list of reciprocal lattice vectors G+q (gvlist.f)
Cio  kv    :indices for gather/scatter operations (gvlist.f)
Cio  igv   :(job1>0) G vectors in integer format
Cio  igv2  :(job2>0) like igv, but array dimensioned in inverse order
Cr Remarks
Cr   Merged with pvglst
Cu Updates
Cu   08 Feb 13 Split off printout
Cu   05 May 12 Made compatible with gvlist
Cu   20 Jul 02 Bug fix for nonzero t
Cu   01 Jun 01 Added igv,igv2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngmx,ng,kv(ngmx,3),igv(ngmx,3),igv2(3,ngmx),job1,job2
      double precision gv(ngmx,3),t(3)
C ... Local parameters
      integer ig,jg,m,ioff
      integer kk(ngmx),iprm(ngmx)
      double precision gg(ngmx)
C     parameter (PRTG=80)

      do  ig = 1, ng
        gg(ig) =(gv(ig,1)+t(1))**2+(gv(ig,2)+t(2))**2+(gv(ig,3)+t(3))**2
      enddo

C#ifdefC DVHEAP
CC     dvheap can preserve original order of vectors of equal length
C      call dvheap(1,ng,gg,iprm,1d-7,111); ioff = 0
C#else
C     but dvshel seems to work a little faster for vectors this size
      call dvshel(1,ng,gg,iprm,1); ioff = 1
C#endif

C ... For each coordinate do
      do  m = 1, 3

C   ... Rearrange gv,kv
        do  ig = 1, ng
          jg = iprm(ig)+ioff
          gg(ig) = gv(jg,m)
          kk(ig) = kv(jg,m)
        enddo
        do  ig = 1, ng
          gv(ig,m) = gg(ig)
          kv(ig,m) = kk(ig)
        enddo
C   ... Rearrange igv
        if (job1 /= 0) then
        do  ig = 1, ng
          jg = iprm(ig)+ioff
          kk(ig) = igv(jg,m)
        enddo
        do  ig = 1, ng
          igv(ig,m) = kk(ig)
        enddo
        endif
C   ... Rearrange igv2
        if (job2 /= 0) then
        do  ig = 1, ng
          jg = iprm(ig)+ioff
          kk(ig) = igv2(m,jg)
        enddo
        do  ig = 1, ng
          igv2(m,ig) = kk(ig)
        enddo
        endif
      enddo

C     call prmx('gv',gv,ngmx,ng,3)

C      if (iprint() >= PRTG) then
C        call gvlstp(job1,job2,alat,qlat,ngmx,ng,q,t,gv,kv,igv,igv2)
C      endif

      end

      subroutine gvlstp(job1,job2,alat,qlat,ngmx,ng,q,t,gv,kv,igv,igv2)
C- Kernel called by gvlist to print out G-vector data
C ----------------------------------------------------------------------
Ci Inputs
Ci   job1  :igv is available
Ci   job2  :igv2 is available
Ci   alat  :lattice constant (for printout)
Ci   ngmx  :dimensions gv,kv,igv
Ci   ng    :number of G-vectors
Ci   q     :constant in gv=G+q
Ci   t     :offset added to gv = G+q to measure length
Ci   kv    :integer work array of length ng
Ci   igv   :(job1>0) G vectors in integer format
Ci   igv2  :(job2>0) like igv, but array dimensioned in inverse order
Co Outputs
Co   G-vectors are printed to stdout
Cr Remarks
Cu Updates
Cu   08 Feb 13 Adapted from gvlsts
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngmx,ng,kv(ngmx,3),igv(ngmx,3),igv2(3,ngmx),job1,job2
      double precision gv(ngmx,3),q(3),t(3),alat,qlat(3,3)
C ... Local parameters
      integer ig,jg,stdo,igvl(3),ipr
      double precision gg,pi,tpiba,qpg(3)
      procedure(integer) :: nglob
      integer, parameter :: PRT2=80

      call getpr(ipr)
      pi = 4d0*datan(1d0)
      tpiba = 2*pi/alat

      stdo = nglob('stdo')
      write(stdo,1) q
    1 format(' GVLST: G+q vectors, q =',3F10.6/
     .  '   ig',17x,'q+G',19x,'kv',9x,'|q+G|',7x,'E',10x,'ig')
      do  ig = 1, ng
        jg = ig
        if (ig == 50 .and. ipr < PRT2) then
          write(stdo,"('  ...')")
          jg = ng
        endif

        if (job1 /= 0 .or. job2 /= 0) then
          if (job2 /= 0) igvl(:) = igv2(:,jg)
          if (job1 /= 0) igvl(:) = igv(jg,:)
          qpg = q + matmul(qlat,igvl)
          gg = (qpg(1)+t(1))**2+(qpg(2)+t(2))**2+(qpg(3)+t(3))**2
          write(stdo,2) jg,qpg,
     .      kv(jg,1),kv(jg,2),kv(jg,3),sqrt(gg),gg*tpiba**2,igvl
        else
          gg = (gv(ig,1)+t(1))**2+(gv(ig,2)+t(2))**2+(gv(ig,3)+t(3))**2
          write(stdo,2) jg,gv(jg,1),gv(jg,2),gv(jg,3),
     .      kv(jg,1),kv(jg,2),kv(jg,3),sqrt(gg),gg*tpiba**2
        endif
    2   format(i5,3f11.5,3i4,2f11.5:3i4)

        if (jg /= ig) exit
      enddo

      end