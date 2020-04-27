      subroutine chgmsh(iopt,n,plat,m1,m2,m3,l1,l2,l3,f0,
     .                        plax,n1,n2,n3,k1,k2,k3,f)
C- Retabulate a function on a new real-space mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt  :0 Use default (smaller of iopt=1,2)
Ci         :1 use Nyquist cutoff
Ci         :2 use cutoff for largest sphere in BZ
Ci         :3 use all the points
Ci         :Add 10 to restore f0 to original value
Ci         :Add 100 to reverse the sense of copy : f -> f0
Ci   n     :number of functions to change
Ci   plat  :primitive lattice vectors, in units of alat
Ci   m1..m3:number of divisions for the original mesh
Ci   l1..l3:dimensions of f0
Ci   f0    :function on m1,m2,m3 mesh
Ci   plax :primitive lattice vectors of new mesh density
Ci         :plax must consist of integral multiples of plat
Ci         :Not used if iopt == 3
Ci   k1..k3:dimensions of f
Cio Inputs/Outputs
Cio  n1..n3:Destination mesh
Cio        :n1,n2,n3 are input, unless iopt eq 3
Cio        :If iopt eq 3  n1..n3 are output as twice (m1..m3)
Co Outputs
Co    f0        FT of f0 is returned in f0
Co    f         function on n1,n2,n3 mesh
Cr Remarks
Cr   f0 and f may occupy the same address space (but not if iopt=13)
Cr   A check is made that lattice vectors plax is contained in plat
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu             New capability for mapping to supercells with plax
Cu   25 Jun 00 added argument n
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iopt,n,m1,m2,m3,n1,n2,n3,l1,l2,l3,k1,k2,k3
      double precision plat(3,3),plax(3,3)
      complex(8) :: f0(l1,l2,l3,n),f(k1,k2,k3,n)
C ... Dynamically allocated arrays
      integer, allocatable :: mapG(:,:,:),kv1(:,:),kv2(:,:),linc(:)
      complex(8), allocatable :: cv1(:)
      real(8), allocatable :: gv1(:,:),gv2(:,:)
      integer, allocatable :: kv1x(:,:),kv2x(:,:)
C      real(8), allocatable :: gv1(:,:),gv2(:,:),xkv1(:,:),xkv2(:,:)
C      complex(8), allocatable :: cv1(:) !,cv2(:)
C ... Local parameters
      logical lsamep,linside(3)
      integer n123(3),m123(3),iplam(3,3),i,j,k,j1,j2,j3,iq123(3),lohi(3,2),ng,iopt0
      double precision plam(3,3),qlat(3,3),qlax(3,3),gmax,qlat1(3,3),qlat2(3,3)
      integer ng1,ng2,ngmx,mm,ig
      double precision gmax1,gmax2
      complex(8) cx
      real(8), parameter :: tol = 1d-6, pi = 4*datan(1d0)
      integer, parameter :: PRT0=30, PRTG=55, PRTG2=70
      procedure(logical) :: latvec
      procedure(integer) :: iprint,isw
      procedure(real(8)) :: dlength

      iopt0 = mod(iopt,10)
      m123 = [m1, m2, m3]
      n123 = [n1, n2, n3]

      call mkqlat(plat,qlat,cx)
      call mkqlat(plax,qlax,cx)
C     Superlattice vectors as multiples of original lattice
      call dgemm('T','N',3,3,3,1d0,plax,3,qlat,3,0d0,plam,3)

C     mc qlax -i qlat -x plax -t  --

C ... Printout
      lsamep = dlength(9,plat-plax,1) < tol ! If plat == plax
      call info8(PRT0,1,0,' CHGMSH: remake f (%i * %i * %i)  to '//
     .    'f (%i * %i * %i)%?#n#  in supercell##,  opt=%i',m1,m2,m3,n1,n2,n3,isw(.not.lsamep),iopt)
      if (.not. lsamep) then
        call info0(PRT0,0,0,'         basis vectors of supercell:%N'//
     .    '%14fplax%27fplax/plat%26fqlax')
        do  k = 1, 3
          call info5(PRT0,0,0,'%3;10,5D   %3;10,5D   %3;10,5D',
     .      plax(1,k),plam(1,k),qlax(1,k),0,0)
        enddo
      endif
      if (.not. latvec(3,tol,qlat,plax)) then
        call logwarn(2,' *** Warning!  supercell lattice '//
     .    'not contained in the original')
        call rx('chgmsh requires commensurate cell')
      endif

      if (iopt0 >= 3) then

      forall (i=1:3, j=1:3) iplam(i,j) = nint(plam(i,j)) ! integer version of plam
      forall (i=1:3) lohi(i,:) = [-n123(i)/2, (n123(i)-1)/2]

C     Make the mapping from small to large cell
      allocate(mapG(3,2,m1*m2*m3))
      ng = 0
      do  j3 = -m3/2, (m3-1)/2
      do  j2 = -m2/2, (m2-1)/2
      do  j1 = -m1/2, (m1-1)/2
        forall (i=1:3) iq123(i) = j1*iplam(i,1) + j2*iplam(i,2) + j3*iplam(i,3)
        forall (i=1:3) linside(i) = iq123(i) >= lohi(i,1) .and. iq123(i) <= lohi(i,2)
C        print 344, j1,j2,j3, iq123, linside, linside(1).and.linside(2).and.linside(3)
C  344   format(3i4,2x,3i4,2x,3L3,L4)
        if (linside(1).and.linside(2).and.linside(3)) then
          ng = ng+1
          mapG(1:3,1,ng) = [j1, j2, j3]
          mapG(1:3,2,ng) = iq123(:)
        endif
      enddo
      enddo
      enddo
      call info2(PRTG,0,0,' chgmsh found %i points in common out of %i in original lattice',ng,m1*m2*m3)
      allocate(kv1(ng,3),kv2(ng,3))

C     Translate from (-m/2,m/2) to (1,..m)
      forall (ig=1:ng) kv1(ig,1:3) = mod(mapG(1:3,1,ig)+m123(:),m123(:))+1
      forall (ig=1:ng) kv2(ig,1:3) = mod(mapG(1:3,2,ig)+n123(:),n123(:))+1

C ... G vectors for printout or shortening
      do  i = 1, 3
        qlat1(i,:) = qlat(i,:)*m123(:)
        qlat2(i,:) = qlat(i,:)*n123(:)
      enddo
      allocate(gv1(3,ng),gv2(3,ng),linc(ng))
      do  ig = 1, ng
        forall (j=1:3) gv1(j,ig) = kv1(ig,1)*qlat(j,1) + kv1(ig,2)*qlat(j,2) + kv1(ig,3)*qlat(j,3)
        forall (j=1:3) gv2(j,ig) = kv2(ig,1)*qlax(j,1) + kv2(ig,2)*qlax(j,2) + kv2(ig,3)*qlax(j,3)
        linc(ig) = 1
      enddo
      call shorps(1,qlat1,[72,2,2],gv1,gv1)
      call shorps(1,qlat2,[72,2,2],gv2,gv2)

C ... Optionally exclude vectors for g>gmax
      iopt0 = mod(iopt,10)
      if (iopt0 < 3) then
        call pshpr(1)
        call gvctof(iopt0,1d0,pi*plat,[0d0,0d0,0d0],m1,m2,m3,gmax,ig)
        call poppr
        do  ig = 1, ng
          if (dlength(3,gv1(1,ig),1) > gmax) linc(ig) = 0
        enddo
      endif

      else  ! End of iopt0 = 3 branch

C     The preceding branch is supposed to cover these cases.
C     For some reason the sorting is different.
C     If the sorting could be made to match the preceding branch could be the sole branch

C ... Lists of vectors for original and target mesh
      call pshpr(max(iprint()-30,1))
      call gvctof(iopt0,1d0,plat,[0d0,0d0,0d0],m1,m2,m3,gmax1,ng1)
      call gvctof(iopt0,1d0,plax,[0d0,0d0,0d0],n1,n2,n3,gmax2,ng2)
C      call gvlst2(1d0,plat,[0d0,0d0,0d0],m1,m2,m3,0d0,gmax1,
C     .  [0],0,ngmx,ng1,kv1,gv1,cx,cx)
C      call gvlst2(1d0,plat,[0d0,0d0,0d0],n1,n2,n3,0d0,gmax2,
C     .  [0],0,ngmx,ng2,kv2,gv2,cx,cx)

      gmax = dmin1(gmax1,gmax2)

C ... Create the two lists and align them
C     Case plat same as plax (eventually not needed)
C     Both this branch and the following should be merged with opt=3.
C     Required : gv1 be sorted by length
      if (lsamep) then
        ngmx = min0(ng1,ng2)
        allocate(gv1(ngmx,3),gv2(ngmx,3),kv1(ngmx,3),kv2(ngmx,3))
        call gvlist(1d0,plat,cx,m1,m2,m3,gmax,8,ngmx,ng1,kv1,gv1,cx,cx)
        call gvlist(1d0,plax,cx,n1,n2,n3,gmax,8,ngmx,ng2,kv2,gv2,cx,cx)
C        call gvlst2(1d0,plat,[0d0,0d0,0d0],m1,m2,m3,0d0,gmax1,
C     .    [0],1008,ngmx,ng1,kv1,gv1,cx,cx)
C        call gvlst2(1d0,plat,[0d0,0d0,0d0],n1,n2,n3,0d0,gmax1,
C     .    [0],1008,ngmx,ng2,kv2,gv2,cx,cx)

        if (ng1 /= ng2) call rx('chgmsh: ng1 /= ng2')
        call gvmtch(ng1,gv1,kv1,ng2,gv2,kv2)
C     Case plat different from plax
      else
        ngmx = max0(ng1,ng2)
        allocate(gv1(ngmx,3),gv2(ngmx,3),kv1x(ngmx,3),kv2x(ngmx,3))
        call gvlist(1d0,plat,cx,m1,m2,m3,gmax,8,ngmx,ng1,kv1x,gv1,cx,cx)
        call gvlist(1d0,plax,cx,n1,n2,n3,gmax,8,ngmx,ng2,kv2x,gv2,cx,
     .    cx)
        if (ng1 > ng2) call rx('chgmsh: ng1 > ng2')
        call gvmtch2(ngmx,ng1,gv1,ng2,gv2,kv2x)
        allocate(kv1(ng1,3),kv2(ng1,3))
        do  mm = 1, 3
          do  ig = 1, ng1
            kv1(ig,mm) = kv1x(ig,mm)
            kv2(ig,mm) = kv2x(ig,mm)
          enddo
        enddo
        deallocate(kv1x,kv2x)
      endif
      deallocate(gv1,gv2)
      call poppr
      ng = ng1
      allocate(linc(ng))
      forall (ig=1:ng) linc(ig)=1
      endif

C ... Copy f0 to f or f to f0
      allocate(cv1(ng*n)); call dpzero(cv1,2*ng*n)
      if (iopt < 100) then
        call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,-1)
C       call zprm3('initial mesh (recip space)',0,f0,l1,l2,l3)
        call gvgetf(ng,n,kv1,l1,l2,l3,f0,cv1)
        forall (ig=1:ng) cv1(ig) = cv1(ig)*linc(ig)
        call gvputf(ng,n,kv2,k1,k2,k3,cv1,f)
C       call zprm3('final mesh (recip space)',0,f,k1,k2,k3)
        call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,1)
      else
        call rx('chgmsh: this branch never checked')
        call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,-1)
        call gvgetf(ng,n,kv2,k1,k2,k3,f,cv1)
        forall (ig=1:ng) cv1(ig) = cv1(ig)*linc(ig)
        call gvputf(ng,n,kv1,l1,l2,l3,cv1,f0)
        call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,1)
      endif

C ... Printout
      call info0(PRTG,0,0,' GVLST: G vectors in original and superlatttice%N'//
     .  '   ig%6fiGlat%14fGlat%18fiGlax%14fGlax%19ff(G)')

      do  ig = 1, ng
        if (ig > 50 .and. iprint() < PRTG2 .and. ig /= ng .or. .not. allocated(gv1)) cycle
        cx = f0(kv1(ig,1),kv1(ig,2),kv1(ig,3),1)
        call info8(PRTG,0,0,
     .    ' %,4i %3,4i %3;9,5D %3,4i %3;9,5D  %?#n#%2:-1;3,3g# ---#',
     .    ig,kv1(ig,:),gv1(:,ig),kv2(ig,:),gv2(:,ig),linc(ig),cx,8)
      enddo

      if (mod(iopt,100) >= 10) then
        if (iopt < 100) call fftz3(f0,m1,m2,m3,l1,l2,l3,n,0,1)
        if (iopt >= 100) call fftz3(f,n1,n2,n3,k1,k2,k3,n,0,-1)
      endif

      deallocate(kv1,kv2,cv1)

      end
      subroutine pchmsh(f0,m1,m2,m3,l1,l2,l3,k1,k2,k3,n,f)
C- Copies Fourier transform on one mesh to a doubled mesh
      implicit none
      integer m1,m2,m3,l1,l2,l3,k1,k2,k3,n
      double complex f0(l1,l2,l3,n),f(k1,k2,k3,n)
      integer i,i1,i2,i3,i1m,i2m,i3m,j1m,j2m,j3m

C     call zprm3('initial mesh',0,f0,l1,l2,l3)

      call dpzero(f,2*k1*k2*k3*n)
      do  10  i = 1, n
      do  10  i3 = 1, (m3+1)/2
      i3m = m3+1-i3
      j3m = 2*m3+1-i3

        do  20  i2 = 1, (m2+1)/2
        i2m = m2+1-i2
        j2m = 2*m2+1-i2

          do  30  i1 = 1, (m1+1)/2
          i1m = m1+1-i1
          j1m = 2*m1+1-i1

            f(i1,i2,i3,i)   = f0(i1,i2,i3,i)
            f(i1,i2,j3m,i)  = f0(i1,i2,i3m,i)
            f(i1,j2m,i3,i)  = f0(i1,i2m,i3,i)
            f(i1,j2m,j3m,i) = f0(i1,i2m,i3m,i)
            f(j1m,i2,i3,i)  = f0(i1m,i2,i3,i)
            f(j1m,i2,j3m,i) = f0(i1m,i2,i3m,i)
            f(j1m,j2m,i3,i) = f0(i1m,i2m,i3,i)
            f(j1m,j2m,j3m,i)= f0(i1m,i2m,i3m,i)


   30     continue
   20   continue
   10 continue

C     call zprm3('final mesh',0,f,2*m1,2*m2,2*m3)
      end
      subroutine pchms2(f0,m1,m2,l1,l2,k1,k2,n,f)
C- 2D analog of pchmsh
      implicit none
      integer m1,m2,l1,l2,k1,k2,n
      double complex f0(l1,l2,n),f(k1,k2,n)
      integer i,i1,i2,i1m,i2m,j1m,j2m

C     call zprm3('initial mesh',0,f0,m1,1,m2)

      call dpzero(f,2*k1*k2*n)
      do  10  i = 1, n
        do  20  i2 = 1, (m2+1)/2
        i2m = m2+1-i2
        j2m = 2*m2+1-i2

          do  30  i1 = 1, (m1+1)/2
          i1m = m1+1-i1
          j1m = 2*m1+1-i1

            f(i1,i2,i)   = f0(i1,i2,i)
            f(i1,i2,i)   = f0(i1,i2,i)
            f(i1,j2m,i)  = f0(i1,i2m,i)
            f(i1,j2m,i)  = f0(i1,i2m,i)
            f(j1m,i2,i)  = f0(i1m,i2,i)
            f(j1m,i2,i)  = f0(i1m,i2,i)
            f(j1m,j2m,i) = f0(i1m,i2m,i)
            f(j1m,j2m,i) = f0(i1m,i2m,i)


   30     continue
   20   continue
   10 continue

C     call zprm3('final mesh',0,f,2*m1,1,2*m2)
      end
