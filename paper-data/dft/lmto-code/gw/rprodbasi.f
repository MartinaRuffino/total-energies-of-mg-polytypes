      integer function rprodbasi(mode,nl,nsp,nrmx,nr,ri,rwt,ndphi,nphi,
     .  nocc,nunocc,phipb,tolbas,lcut,off,nblraw,nblr,rprodb,fiterr)
C- Construct product basis for one site
C ----------------------------------------------------------------------
Ci  mode   : If NULLI, return version number.  Otherwise:
Ci         : 1s digit should be anything but 4 for now.
Ci         :    0 normal
Ci         :    The following are old conventions from basnfp_v2
Ci         :    0 normal
Ci         :    3 coremode (same as 0)
Ci         :    4 ptest (not implemented)
Ci         :    5 Excore (same as 0)
Ci         :    6 core-valence (same as 0)
Ci         :    7 val-val Ex (same as 0)
Ci         :    8 normal + <rho_spin|B> (not implemented)
Ci         :10s digit
Ci         :    0 no gradients returned
Ci         :    1 Use dg(l)/r for g(l)   for gradient (only one kind is made)
Ci         :    2 Use dg/dr - (l+1)/r*g  for gradient (grad+)
Ci         :    3 Use dg/dr + l/r*g      for gradient (grad-)
Ci         :    Add 4 to return (grad phi2) * phi1 instead of phi2 * (grad phi1)
Ci         :100s digit
Ci         :    1 check completeness of raw product basis (commented out for now)
Ci         :    2 check completeness of truncated product basis
Ci         :    4 check completeness of truncated product basis, verbose output
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  ndphi  :dimensioning parameter : max number of valence partial waves of a particular l
Ci         :ndphi called nn in old GW code
Ci  nrmx   :leading dimension of gval,gcore
Ci  nr     :Number of radial mesh points
Ci  ri     :ri(1:nr) = radial mesh
Ci  rwt    :rwt(1:nr) = radial mesh weights that numerically integrate int dr f(r)
Ci  nphi   :number of (core states + valence partial waves) for a particular l
Ci         :nphi formerly named nindx in old GW code
Ci  nocc   :0 or 1 for l=0..nl-1 and i = 1..ndphi
Ci         :1 => include as left function when assembling
Ci         :     possible combinations for product basis
Ci         :0 => do not include
Ci  nunocc :0 or 1 for l=0..nl-1 and i = 1..ndphi
Ci         :1 => include as right function when assembling
Ci         :     possible combinations for product basis
Ci         :0 => do not include
Ci  phipb  :basis of (r * partial waves) that expand MTO's
Ci         :Formerly named phitotr in old GW code
Ci  tolbas :product basis tolerance... one number for every l.  See Remarks
Ci  lcut   :l-cutoff in product basis
Ci  off    : <0 do not return rprodb or any of its gradients
Ci         :>=0 return rprodb(:,off+1:off+nradpbi)
Co Outputs
Co  nblraw :number of radial functions for each l, before orthonormalization
Co  nblr   :number of radial functions for each l, after orthonormalization
Cl         :Formerly named nxx in old GW code
Co  probasi:total number of radial functions
Co  rprodb :product basis functions constructed from r * (phi1/r) * (phi2/r)
Co         :Note that rprodb = r * true B(r).
Co         :rprodb(:,off+1:off+nradpbi) is returned
Co  off    :If non-negative, rprodb is returned and off is incremented by rprodbasi
Cl Local variables
Cl  irad   :index to current radial product function
Cl  nprad  :raw number of radial product functions before reduction by orthogonalization
Cl  nprod  :raw total total number of radial product functions
Cl  rprod  :a vector (1:irad) of all (phi1,phi2) pairs, tabulated on a radial mesh.
Cl         :See npr supplies indexing that connect a particular irad to a specific (phi1,phi2) pair
Cl  npr    :irad = npr(l1,n1,l2,n2) = which (l1,n1) in phi1 and (l1,n2) in phi2 make rprod(:,irad)
Cl  phiav  :spin-average of partial waves phipb, used to make product basis
Cl  ibl    :see Remarks.  Formerly named ixx in old GW code
Cr Remarks
Cr   All pairs of participating partial waves (see nocc, nunocc) are initially
Cr   created.  Their overlap matrix is computed and diagonalized.
Cr   The product basis is orthonormalized and functions with eigenvalues
Cr   of the overlap < tolbas are removed from the basis.
Cr
Cr   In more detail:
Cr   A complete list of product functions is initially tabulated, numbering irad functions
Cr   Indexing arrays keep track of attributes of each function:
Cr     nblr(l)   records how many radial functions are present for a particular l
Cr     ibl(l,i)  points to which of the irad radial functions is the ith function of a particular l
Cr
Cr   Adapted from T. Kotani (basnfp.F), modified from F. A. (indxbas.f)
Cr
Cr   Remarks on basnfp_v2, from which this code was derived.
Cr   basnfp_v2 makes both the product basis and matrix elements of partial waves with it.
Cr   rprodbasi makes only the product basis; see prodbasme for the matrix elements
Cr   phipb is a basis for the expansion of actual partial waves phime used in the GW suite.
Cr
Cr   phime :Orthogonalized form of phipb (see gwcphi)
Cr         :These are partial waves used for matrix elements <phi1 phi2 | B>
Cr         :For valence states, phime consist of valence part of phib only
Cr         :One-center expansion of eigenfunction is (apart from L convergence) : cphi x phime.
Cr         :Valence waves should be orthogonalized, while while core states unchanged from phipb.
Cr   basnfp stores product basis and their matrix elements in files:
Cr   BASFP//atom     :  product basis. used in hvccfp0.
Cr   PPBRD_V2_//atom : <phi phi |B> radial integral. Read in hx0fp0 hsfp0 through rdpp_v2.
Cu Updates
Cu   19 Aug 18 First cut at gradients
Cu   12 Dec 17 Adapted from old GW basnfp_v2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,ndphi,nrmx,nr,nsp,lcut,mode,off
      integer nocc(0:nl-1,ndphi),nunocc(0:nl-1,ndphi),nphi(0:nl-1),nblraw(0:2*(nl-1)),nblr(0:2*(nl-1))
      double precision ri(nrmx),rwt(nrmx),tolbas(0:2*(nl-1)),fiterr(2)
      real(8),target :: phipb(nrmx,0:nl-1,ndphi,nsp)
      double precision rprodb(nrmx,*)
C ... Dynamically allocated arrays
      integer,allocatable ::  ibl(:,:),nnode(:)
      real(8),allocatable,target :: rprodl(:,:,:),rprod(:,:)! ,rgprod(:,:,:) !,rgprodl(:,:,:,:)
      real(8),allocatable :: phi21(:,:,:,:,:,:) ! ,gphi21a(:,:,:,:,:,:),gphi21b(:,:,:,:,:,:)
      real(8),allocatable :: eb(:),ebx(:),z(:,:),zt(:,:),ovb(:,:),ovbl(:,:),rprbme(:,:,:,:,:,:),rwk(:)
      real(8),pointer ::  phi1(:,:,:,:),phi2(:,:,:,:)
C ... Local parameters
      logical :: ldebug = .false.
      integer,parameter:: NULLI = -99999
      integer,parameter:: ivsn = 1
      integer,parameter:: nblmx=300
      integer,parameter:: ndif=8  ! number of points in num. diff of phi
      logical :: lgphi1,lgphi2
      integer i,ib1,ib2,ierr,l1,l2,lb,mode0,mode1,n1,n2,nb,nblrmx,noo,nuu,nprod,nradpb
      integer irad              ! index to current radial product function; finally, the number generated
      integer offsave           ! Hang on to passed offset; needed at the check phase at the end
      integer npr(0:nl-1,ndphi,0:nl-1,ndphi),nblr0(0:2*(nl-1))
      double precision fv1(nblmx),fv2(nblmx)
      double precision ovphi(ndphi,ndphi,0:nl-1,nsp)
      procedure(integer) :: fopng,iprint,nodnum,iinear
      procedure(real(8)) :: dot3

C --- Setup ---
      if (mode == NULLI) then
        rprodbasi = ivsn
        return
      endif
      mode0 = mod(mode,10); mode1 = mod(mod(mode/10,10),4); offsave = off
      if (mode0 == 4) call rx('prodbas:  mode 4 has been removed')
      lgphi1 = .false.; lgphi2 = .false.; phi1 => phipb; phi2 => phipb
      allocate(phi21(nr,0:nl-1,ndphi,0:nl-1,ndphi,nsp))

      if (mode0 > 0) then       ! Special purpose gradients
        call rx('no special product basis available yet')
      elseif (mode1 == 0) then  ! Regular product basis
        i = 13
!       call prodphi(13,nl,nsp,ndphi,nrmx,nr,ri,rwt,nphi,phi2,phi1,ovphi,phi21) ! ovphi returned but not used
!       call prodphi(12,nl,nsp,ndphi,nrmx,nr,ri,rwt,nphi,phi2,phi1,fv1,phi21)
      else                      ! gradient of regular product basis
C       call snot
        if (mod(mode/10,10)>4) then
          i = 1000
          if (mod(mode/100,10) > 1) then
            allocate(phi2(nrmx,0:nl-1,ndphi,nsp)) ! phi2 is first product in phi2 phi1
            call dcopy(size(phipb),phipb,1,phi2,1)
            lgphi2 = .true.
          endif
        else
          i = 100
          if (mod(mode/100,10) > 1) then
            allocate(phi1(nrmx,0:nl-1,ndphi,nsp)) ! phi1 is second product in phi2 phi1
            call dcopy(size(phipb),phipb,1,phi1,1)
            lgphi1 = .true.
          endif
        endif
        i = i*mode1 + 12        ! Set type of gradient; return phi21; spin average
        if (mod(mode/100,10) > 1) i = i + 4 ! return grad phipb in phi1
      endif
      call prodphi(i,nl,nsp,ndphi,nrmx,nr,ri,rwt,nphi,phi2,phi1,fv1,phi21)

C --- Product basis construction ---
      allocate(rprod(1:nr,nblmx),ibl(0:2*(nl-1),nblmx))
      call dpzero(rprod,size(rprod))
      npr   = 0                 ! marks which (phi1,phi2) pair generated
      nblr  = 0                 ! Tallies how many radial functions belong to a particular l
C ... Make radial products phiav * phiav and indexing arrays npr,nblr,ibl
      irad = 0                  ! running index for all product partial waves
      do  l1 = 0, nl-1
      do  n1 = 1, nphi(l1)
        noo = nocc(l1,n1)
        do  l2 = 0, nl-1
        do  n2 = 1, nphi(l2)
          nuu = nunocc(l2,n2)

C         Nothing to do unless both occ and unocc in list, and (phi1,phi2) not already in list
          if (noo == 0 .or. nuu == 0 .or. npr(l1,n1,l2,n2) > 0) cycle
          irad = irad+1
          npr(l1,n1,l2,n2) = irad; if (mode1 == 0) npr(l2,n2,l1,n1) = irad
C         l's participating in YL expansion of phiav(l1)*phiav(l2)
          do  lb = abs(l1-l2), l1+l2 ! Radial function can be spread over these l's
            if (mod(lb+l1+l2,2) == 1) cycle
            if (lb > lcut) cycle  ! l cutoff for product basis
            nblr(lb) = nblr(lb)+1
            if (nblr(lb) > nblmx) call rx('Enlarge nblmx in prodbas')
            ibl(lb,nblr(lb)) = irad
          enddo
CC         Debugging ... isolate lb=1 cases
CC           if (ldebug) then
C            i = iinear(nblr(1),irad,ibl(1,1:nblr(1)),1)
C            if (i > 0) then
C              if (irad-ibl(1,i) == 0) then
C                print *, 'hi',l1,n1,l2,n2,irad, sngl(gphi21a(nr,l1,n1,l2,n2,1))
C              endif
C            endif
CC           endif
          rprod(1,irad) = 0d0
          rprod(2:nr,irad) = phi21(2:,l1,n1,l2,n2,1) ! rprod =  r * (phi2/r) * (phi1/r)
        enddo ! loop over n2
        enddo ! loop over l2
      enddo
      enddo  ! loop over l1

C ... Add r^l and r^l+1 for l=0,1.
C     Query: why not r^l and r^(l+2)?  This is energy derivative of Bessel.
C     May be an adequate compromise to make good <M|v|M>.
      do  lb = 0, 1
        irad = irad+2           !Increment counter for 2 new radial functions
        nblr(lb) = nblr(lb)+2
        ibl(lb,nblr(lb)-1) = irad-1; ibl(lb,nblr(lb)) = irad
        if (nblr(lb) > nblmx) call rx('Enlarge nblmx in rprodbasi')
        rprod(1:nr,irad-1) = ri(1:nr)*ri(1:nr)**lb ! Remember rprod = r*[true B]
        rprod(1:nr,irad)   = ri(1:nr)*ri(1:nr)**(lb+1)
        if (mode1 > 0) then
          allocate(rwk(nr))
          call dcopy(nr,rprod(1,irad-1),1,rwk,1)
          call prodphig(mode1,lb,nr,ri,rwk,rprod(1,irad-1))
          call dcopy(nr,rprod(1,irad),1,rwk,1)
          if (mode1 == 2) then
            irad = irad-1       ! grad+ r^l is 0; discard this function
            nblr(lb) = nblr(lb)-1
!           ibl(lb,nblr(lb)-1) = irad-1
          endif
          call prodphig(mode1,lb,nr,ri,rwk,rprod(1,irad))
          deallocate(rwk)
        endif
      enddo

C     Count number of functions
      nprod = 0
      do  lb = 0, 2*(nl-1)
        nprod = nprod + (2*lb+1)*nblr(lb)
      enddo
      nblraw = nblr
      nblrmx = maxval(nblr)
      call info5(50,1,0,' %i radial, %i total product functions.  Decomposed by l: %n:1i',irad,nprod,nl,nblr,5)

C --- Orthonormalize basis functions ---
      allocate(rprodl(nr,nblrmx,0:2*(nl-1)),ovbl(nblrmx*nblrmx,0:2*(nl-1))) ! To handle either pb or gbp
C ... For each lb, do
      do  lb = 0, 2*(nl-1)
        nb = nblr(lb)
        nblr0(lb) = nblr(lb)
        if (nb == 0) cycle
        allocate(ovb(nb,nb),zt(nb,nb),z(nb,nb),eb(nb),ebx(nb),nnode(nb))
        do  ib1 = 1, nb
          do  ib2 = 1, nb
            ovb(ib1,ib2) = dot3(nr,rprod(1,ibl(lb,ib1)),rprod(1,ibl(lb,ib2)),rwt)
          enddo
        enddo
C   ... Preserve in case checking for completeness of raw product basis
        call dcopy(nb*nb,ovb,1,ovbl(1,lb),1)

        if (ldebug) call yprmi('overlap in B l=%i',lb,0,1,ovb,0,nb,nb,nb)
        call rs(nb,nb,ovb,eb,1,z,fv1,fv2,ierr)
        if (ierr/=0) call rx('prodbas failed to diagonalize overlap matrix')
        if (ldebug) call yprmi('evecs of overlap l=%i',lb,0,1,z,0,nb,nb,nb)

        ib2 = 0; call dpzero(ovb,nb); call ivset(nnode,1,nb,-1)

C       Debugging
C       By hand for 11 functions.  put pb in file pb, e in file e, z in file z, rprodl(nr,1:ib2,lb) in dat
C       for gradient, put gpb in file, rgprodl(nr,1:ib2,lb) in datg
C       mcx -r:nc=11 pb -a pb -r:nc=11 e -a e pb z -x  e -t -e1 '1/sqrt(x1)' -t -xe -coll:mirr -coll 1:5 dat --
C       mcx -r:nc=11 gpb -a pb -r:nc=11 e -a e pb z -x  e -t -e1 '1/sqrt(x1)' -t -xe -coll:mirr -coll 1:5 datg --
C        if (mode1 > 0) then
C        call yprmi('evecs of overlap l=%i',lb,0,1,z,0,nb,nb,nb)
C        print *, 'pb'
C        print *, rprod(nr,ibl(lb,1:nb))
C        print *
C        print *, 'gpb'
C        print *, rgprod(nr,ibl(lb,1:nb),1)
C        print *
C        print *, 'eb'
C        print *, eb(1:nb)
C        print *
C        endif

C   ... Scale z so that z^T ovlb z = 1, rotate braw->borth in reverse order and cull eb<tolbas
C       Debugging : put : rprod(nr) into braw   scaled z into z0s
C       Then do: mc braw -t z0s -x -coll:mirr -t borth --
C       Also :   mc borth -t -coll:mirr z0s -i -x -t braw --
C                mc borth -t z0s -coll:mirr -i -x -t braw --
C                mc borth -t z0s -i -rowl:mirr -x -t braw --
        do  ib1 = nb, 1, -1
          ebx(nb-ib1+1) = eb(ib1)
          if (eb(ib1) < tolbas(lb)) cycle
          ib2 = ib2+1
          zt(ib2,1:nb) = z(1:nb,ib1)*sqrt(eb(ib1))
          z(1:nb,ib1) = z(1:nb,ib1)/sqrt(eb(ib1))
          rprodl(1:nr,ib2,lb) = matmul(rprod(1:nr,ibl(lb,1:nb)),z(1:nb,ib1))
C          if (mode1 > 0) then
C            rgprodl(1:nr,ib2,lb,1) = matmul(rgprod(1:nr,ibl(lb,1:nb),1),z(1:nb,ib1))
C            if (mode1 == 2) then
C            rgprodl(1:nr,ib2,lb,2) = matmul(rgprod(1:nr,ibl(lb,1:nb),2),z(1:nb,ib1))
C            endif
C          endif
          ovb(ib2,1) = dot3(nr,rprodl(1,ib2,lb),rprodl(1,ib2,lb),rwt)
          nnode(ib2) = nodnum(rprodl(1,ib2,lb),nr)

        enddo  ! Loop over ib1

C       Copy overlap to file ovb and check the following:
C       mc -f9f12.6 ovlb -evc ovlb -evl -e1 'sqrt(1/x1)' -v2dia -x out.cr3si6 --
        if (ldebug) call yprmi('scaled evecs of overlap l=%i',lb,0,1,z,0,nb,nb,nb)
        if (ldebug) call yprmi('raw  b(nr), l=%i',lb,0,1,rprod(nr,ibl(lb,1:nb)),0,nb,nb,1)
        if (ldebug) call yprmi('orth b(nr), l=%i',lb,0,1,rprodl(nr,1:nb,lb),0,nb,nb,1)

        if (iprint() >= 50) then
        call info8(50,0,0,' l=%i : keep %i, discard %i waves for ecut < %;3g : %n:;2g ',
     .    lb,ib2,nb-ib2,tolbas(lb),nb-ib2,ebx(min(ib2+1,nb)),7,8)
C       Probably this should not bother to print out norm.
        call arrprt(' ib  eval norm nnod','%,3i %;5,5F %;5,5F%,2i','Iddi',ib2,0,4,
     .      0,'  | ',ovb,ebx,ovb,nnode,ovb,ovb,ovb,ovb)
        endif
        nb = ib2; nblr(lb) = ib2

        deallocate(ovb,z,zt,eb,ebx,nnode)
      enddo  ! loop over lb

C ... Return rprodb, if off is not negative
      i = sum(nblr)
      rprodbasi = i
      if (off >= 0) then
        do  lb = 0, 2*(nl-1)
          if (nblr(lb) == 0) cycle
          do  i = 1, nblr(lb)
C           print *, 'lb,i,off,nr,nrmx',lb,i,off,nr,nrmx
            call dvset(rprodb(1,1+off),nr+1,nrmx,0d0)
            call dcopy(nrmx,rprodl(1,i,lb),1,rprodb(1,1+off),1)
C            if (mode1 > 0) then
C              call dcopy(nrmx,rgprodl(1,i,lb,1),1,rgprodb(1,1,1+off),1)
C              if (mode1 == 2)
C     .        call dcopy(nrmx,rgprodl(1,i,lb,2),1,rgprodb(1,2,1+off),1)
C            endif
            off = off + 1
          enddo
        enddo
      endif


C --- Debugging: check raw product basis for completeness ---
C      if (mod(mode/100,10) > 0) then
C        nradpb = sum(nblr0)
C        deallocate(rprodl)
C        allocate(rprodl(nrmx,nradpb,1))
C        allocate(rprbme(0:nl-1,ndphi,0:nl-1,ndphi,0:2*(nl-1),nblrmx))
C        n2 = 0
C        do  lb = 0, 2*(nl-1)
C          nb = nblr0(lb)
C          n1 = n2+1
C          n2 = n2+nb
C          rprodl(1:nr,n1:n2,1) = rprod(1:nr,ibl(lb,1:nb))
C        enddo
C        call rprodbaschki(20+mode/200,nl,1,1,nblr0,ndphi,nrmx,nr,ri,rwt,nphi,-1,
C     .    nocc,nunocc,rprodl,phi1,phi2,rprbme,fiterr)
C        deallocate(rprbme)
C        print *, '----- end debug ----'
C        stop
C      endif

C --- Check product basis for completeness ---
      if (mod(mode/100,10) > 1) then
        nradpb = sum(nblr)
        allocate(rprbme(0:nl-1,ndphi,0:nl-1,ndphi,0:2*(nl-1),nblrmx))
        call rprodbaschki(20+mode/200,nl,1,1,nblr,ndphi,nrmx,nr,ri,rwt,nphi,-1,
     .    nocc,nunocc,rprodb(1,1+offsave),phi1,phi2,rprbme,fiterr)
        deallocate(rprbme)
      endif
      deallocate(rprod,rprodl,ibl,ovbl,phi21)
      if (lgphi2) deallocate(phi2)
      if (lgphi1) deallocate(phi1)

      end
      subroutine rprodbaschki(mode,nl,isp1,isp2,nblr,ndphi,nrmx,nr,ri,rwt,nphiv,nphic,nocc,nunocc,rprodb,
     .  phi2,phi1,rprbme,fiterr)
C- Check expansion of partial wave products in terms of product basis functions
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :1s digit controls what to write to stdout
Ci         :0  write nothing
Ci         :2  write information about each (l1,n1,l2,n2) pair
Ci         :10s digit
Ci         :1  Check completeness of core waves
Ci         :2  Check completeness of valence waves
Ci         :3  Combination of 1+2
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  isp1   :spin index for first partial wave
Ci  isp2   :spin index for second partial wave
Ci  nblr   :number of radial functions for each l
Ci  ndphi  :dimensioning parameter : max number of valence partial waves of a particular l
Ci         :ndphiv formerly named nn in old GW code
Ci  nrmx   :leading dimension of gtoto, rprodb
Ci         :nrmx formerly named nrofi in old GW code
Ci  nr     :Number of radial mesh points
Ci  ri     :ri(1:nr) = radial mesh
Ci  rwt    :rwt(1:nr) = radial mesh weights
Ci  nphiv  :number of valence partial waves for this l and site
Ci         :nphiv formerly named nphitv in old GW code
Ci         :If nphv(0) = -1, nphiv is not used.  In that case
Ci         :phi contains no valence states or nphic encompasss all the states
Ci  nphic  :number of core levels for this l and site
Ci         :If nphic(0) = -1, nphic is not used.  In that case
Ci         :phi contains no core states or nphiv encompasss all the states
Ci         :nphic formerly named nphitc in old GW code
Ci  nocc   :0 or 1 for l=0..nl-1 and i = 1..ndphi
Ci         :1 => include as left function when assembling
Ci         :     possible combinations for product basis
Ci         :0 => do not include
Ci         :Not used if nocc(0,1) < 0.  If it is used,
Ci         :the ordering nocc(l,:) should synchronize phi
Ci  nunocc :0 or 1 for l=0..nl-1 and i = 1..ndphi
Ci         :1 => include as right function when assembling
Ci         :     possible combinations for product basis
Ci         :0 => do not include
Ci         :Not used if nunocc(0,1) < 0.  If it is used,
Ci         :the ordering nunocc(l,:) should synchronize phi
Ci  rprodb :radial parts of product basis functions B for this site
Ci         :there is no requirement for orthonormality
Ci  phi2   :(r * partial waves) that expand MTO's and one of (phi2,phi1) pairs making product
Ci         :There is no requirement they be orthogonal.
Ci         :phi2 can consist of a concatenation of core + valence partial waves;
Ci         :the ordering should synchronize with nphic + nphiv.
Ci  phi1   :Second of (phi2,phi1) pairs forming product.
Ci         :phi1 can point to same address space as phi2, which see.
Co Outputs
Co  rprbme :radial matrix elements < gtoto gtoto B> for one site
Cl Local variables
Cr Remarks
Cr   General statement: If B is complete, sum_J |B_J(r')><B_J(r)| = delta(r-r')
Cr   Then
Cr      phi2(r) phi1(r) = int dr' phi2(r') phi1(r') sum_J |B_J(r')><B_J(r)|
Cr                      = sum_J B_J(r) <phi2(r') phi1(r') |B_J(r')>
Cr   A more limited form of completeness:
Cr   minimize difference <[phi2 phi1 - sum_J C_J |B_J]^2>
Cr   If difference can be made 0, B spans phi2 phi1.  Otherwise C_J gives best fit.
Cr   Let
Cr     f_J = sum_I C_I O_IJ  with f_J = <B_J | phi2 phi1>  and   O_IJ = <B_I | B_J>
Cr   This routine computes
Cr     C_J = sum_J f_J O^-1_IJ
Cr   and compares
Cr     sum C_J B_J(r) to r * (phi2(r)/r) * (phi1(r)/r)
Cu Updates
Cu   05 Jul 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,isp1,isp2,nrmx,nr,ndphi
      integer nocc(0:nl-1,ndphi),nunocc(0:nl-1,ndphi)
      integer nblr(0:2*(nl-1)),nphiv(0:nl-1),nphic(0:nl-1)
      double precision rprodb(nrmx,*),ri(nr),rwt(nr),fiterr(2)
      double precision rprbme(0:nl-1,ndphi,0:nl-1,ndphi,0:2*(nl-1),maxval(nblr))
      double precision phi2(nrmx,0:nl-1,ndphi,max(isp1,isp2))
      double precision phi1(nrmx,0:nl-1,ndphi,max(isp1,isp2))
C ... Dynamically allocated local arrays
      real(8),allocatable :: cj(:),wk(:,:),ovb(:,:)
C ... Local parameters
      logical :: ldebug = .false.
      integer irad,lb,nb,l1,n1,l2,n2,offrpb,mode0,mode1,i,noo,nuu,ib,nrms
      integer nphi(0:nl-1),nphi2(0:nl-1),nphi1(0:nl-1)
      double precision rphiphi(nrmx),rpb(nrmx,2),xv(10),errl(0:nl-1,ndphi,0:nl-1,ndphi),avgrms,maxrms
      procedure(real(8)) :: dot3

      if (sum(nblr(0:2*(nl-1))) == 0) return
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      call dpzero(errl,size(errl))
      fiterr = 0

      if (mode0 >= 2) call info0(30,0,0,
     .  '  l1  n1  l2  n2  lb  nb%4fphi*phi(2)%3fB(2)%4fphi*phi(nr)%3fB(nr)%6f<B>%6f<B-gg>  sqr<(B-gg)^2>')

C ... Assemble nphi=nphic+nphiv, nphi2=lower bound to n1,n2,  nphi1=upper bound to n1,n2
      nphi = 0
      if (nphic(0) >= 0) nphi = nphi + nphic
      if (nphiv(0) >= 0) nphi = nphi + nphiv
      nphi2 = 1
      if (nphic(0) >= 0 .and. mod(mode1,2) == 0) nphi2 = nphic+1 ! phi contains cores but do not evaluate them
      nphi1 = nphi
      if (nphiv(0) >= 0 .and. mod(mode1/2,2) == 0) nphi1 = nphic ! phi contains val but do not evaluate them

C     Make f_J = <phi2 phi1 rprodb>
      call prodbasmei(0,nl,1,nblr,ndphi,nrmx,nr,ri,rwt,nphi,rprodb,phi2,phi1,rprbme)

      noo = 1 ; nuu = 1; offrpb = 0; avgrms = 0; maxrms = 0; nrms = 0
      do  lb = 0, 2*(nl-1)
        nb = nblr(lb)
        if (lb>0) offrpb = offrpb + nblr(lb-1)
        if (nb == 0) cycle
        allocate(ovb(nb,nb),wk(nb,nb),cj(nb))

C       Make O^-1 for this lb
        do  n1 = 1, nb
          do  n2 = 1, nb
            ovb(n1,n2) = dot3(nr,rprodb(1,offrpb+n1),rprodb(1,offrpb+n2),rwt)
          enddo
        enddo
C       call yprmi('rprodbaschki: PB overlap l=%i',lb,0,1,ovb,0,nb,nb,nb)
        call dqinv('s',ovb,nb,0,nb,wk,nb,i)

        do  l1 = 0, nl-1
        do  n1 = nphi2(l1), nphi1(l1)
          if (nocc(0,1) >= 0) noo = nocc(l1,n1)

          do  l2 = 0, nl-1
          if (lb<abs(l1-l2) .or. l1+l2<lb) cycle
          if (mod(lb+l1+l2,2) == 1) cycle
          do  n2 = 1, nphi(l2)

            if (nunocc(0,1) >= 0) nuu = nunocc(l2,n2)

C           Nothing to do unless both occ and unocc in list, and (phi2,phi1) not already in list
            if (noo == 0 .or. nuu== 0) cycle

C           Make coefficients cj
            call dgemm('N','N',nb,1,nb,1d0,ovb,nb,rprbme(l1,n1,l2,n2,lb,1:nb),nb,0d0,cj,nb)
            if (ldebug) call yprmi('f_J l=%i',lb,0,1,rprbme(l1,n1,l2,n2,lb,1:nb),0,nb,nb,1)

            rphiphi(1)    = 0d0
            rphiphi(2:nr) = phi2(2:nr,l1,n1,isp2) * phi1(2:nr,l2,n2,isp1)/ri(2:) ! phi = u = r \phi
            call dpzero(rpb,nr)
            do  ib = 1, nb
              irad = ib + offrpb
              call daxpy(nr,cj(ib),rprodb(1,irad),1,rpb,1)
            enddo
            call dcopy(nr,rphiphi,1,rpb(1,2),1)
            call daxpy(nr,-1d0,rpb,1,rpb(1,2),1)

C           Debugging
C            print *, 'debugging l1,l2=',l1,l2
C            rpb(1,1) = 0
C            rpb(2:nr,1) = rpb(2:nr,1)/ri(2:nr)
C            rpb(2:nr,2) = rphiphi(2:nr)/ri(2:nr)
C            print *, l1+l2, rphiphi(2)/ri(2), rphiphi(2)/ri(2)/ri(2)**(l1+l2)
C            call prrmsh('rphiphi, rpb',ri,rpb,nr,nr,2)

            xv(1) = rphiphi(2)/ri(2)
            xv(2) = rpb(2,1)/ri(2)
            xv(3) = rphiphi(nr)
            xv(4) = rpb(nr,1)
            xv(5) = dot3(nr,ri,rpb(1,1),rwt)
            xv(6) = dot3(nr,ri,rpb(1,2),rwt)
            xv(7) = sqrt(dot3(nr,rpb(1,2),rpb(1,2),rwt))
            nrms = nrms + 1
            maxrms = max(maxrms,xv(7))
            avgrms = avgrms + xv(6)
            if (mode0 >= 2) then
            call info2(30,0,0,'%6:1,3i %7:4;7,7F%?#n#',[l1,n1,l2,n2,lb,nb],xv)
            endif

C            errl(l1,n2,l2,n2) = xv(4)

          enddo               ! loop over n2
          enddo               ! loop over l2
        enddo
        enddo                 ! loop over l1

        deallocate(ovb,wk,cj)
      enddo
      fiterr(1) = avgrms/nrms
      fiterr(2) = maxrms

      end
C      subroutine snot
C      end
