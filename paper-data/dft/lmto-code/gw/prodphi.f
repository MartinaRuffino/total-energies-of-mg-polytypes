      subroutine prodphi(mode,nl,nsp,ndrphi,nrmx,nr,ri,rwt,nrphi,phi2,phi1,ovphi,phi21)
C- products of partial waves phi2 * phi1 for one site, or some gradient of phi2*phi1
C ----------------------------------------------------------------------
Ci  mode   : 1s digit governs what is returned
Ci         :    0 return version number in mode
Ci         :    1 return overlap in ovphi, possibly decomposed (see 1000s digit)
Ci         :    2 return products of partial waves in phi21
Ci         :    4 replace phi1 with gradient if 100s digit set
Ci         :      replace phi2 with gradient if 1000s digit set
Ci         :    Any combination of 1,2,4 allowed
Ci         :10s digit modifies kind partial waves used to make product basis
Ci         :    0 Use partial waves as given
Ci         :    1 Take the spin average of partial waves
Ci         :    2 Return phi1(isp=1) * phi2(isp=2) and phi1(isp=2) * phi2(isp=1)
Ci         :100s digit
Ci         :    0 use phi1 for g1 in constructing product basis g2 * g1 / r
Ci         :    1 Replace g1(l) with g1'(l,r) = dg1(l,r)/r
Ci         :    2 Replace g1(l) with g1'(l,r) - (l+1)g1(l,r)/r
Ci         :      Note that angular momentum associated with g1 is l+1, not l
Ci         :    4 Replace g1(l) with g1'(l,r) + l*g1(l,r)/r
Ci         :      Note that angular momentum associated with g1 is l-1, not l
Ci         :      Switches 1,2,4 are mutually exclusive
Ci         :1000s digit
Ci         :    0 use phi2 for g2 in constructing product basis g2 * g1 / r
Ci         :    1 Replace g2(l) with g2'(l,r) = dg2(l,r)/r
Ci         :    2 Replace g2(l) with g2'(l,r) - (l+1)g2(l,r)/r
Ci         :      Note that angular momentum associated with g2 is l+1, not l
Ci         :    4 Replace g2(l) with g2'(l,r) + l*g2(l,r)/r
Ci         :      Note that angular momentum associated with g2 is l-1, not l
Ci         :      Switches 1,2,4 are mutually exclusive
Ci         :10000s digit affects what is returned in ovphi, if 1s digit mode>0
Ci         :    0 return overlap
Ci         :    1 return Cholesky-decomposed overlap
Ci         :    2 return Cholesky-decomposed inverse of overlap
Ci         :    3 return eigenvectors of ovphi / sqrt(eigenvalue)
Ci         :    4 return eigenvectors of scaled ovphi / sqrt(eigenvalue),
Ci         :      using Kotani convention for scaling: sqrt(1d0+0.1d0*n1)*sqrt(1d0+0.1d0*n2)
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  ndrphi :dimensioning parameter : max number of partial waves of a particular l
Ci         :ndrphi called nn in old GW code
Ci  nrmx   :leading dimension of phi1, phi2
Ci  nr     :Number of radial mesh points for this site
Ci  ri     :ri(1:nr) = radial mesh
Ci  rwt    :rwt(1:nr) = radial mesh weights to integrate int dr f(r)
Ci  nrphi  :number of (core states + valence partial waves) for a particular l
Ci         :nrphi formerly named nindx in old GW code
Ci  phi1   :partial waves enter into product, multiplied by r, right factor
Ci         :True partial wave is psi_lm = (phi(r)/r)_l*Y_lm
Ci  phi2   :partial waves enter into product, multiplied by r, left factor
Ci         :True partial wave is psi_lm = (phi(r)/r)_l*Y_lm
Co Outputs
Co  ovphi  :(mode=1) Overlap between partial waves
Co  phi21  :(mode=2) products of partial waves = r * (phi2/r) * (phi1/r)
Cl Local variables
Cr Remarks
Cr   This routine returns:
Cr   1. Overlap matrix ovphi between phi1 and phi2
Cr   2. Overlap matrix is Cholesky decomposed or diagonalized
Cr   3. products phi21(1:nr,;l1,l2,n1,n2,isp) of phi1 and phi2, or their gradient
Cr
Cr   Remarks concerning gradients.
Cr   A product basis constructed from r g2(l')/r grad g1(l)/r generates two kinds of radial functions:
Cr   because grad g1(l)/r generate f+(l+1)  proportional to dg/dr - (l+1)g/r
Cr                             and f-(l-1)  proportional to dg/dr + (l)g/r
Cr   See https://www.questaal.org/docs/numerics/spherical_harmonics/#gradients-of-spherical-harmonics
Cr   Thus if 100's or 1000's digit mode is 2, it is shifted to one higher l
Cr        if 100's or 1000's digit mode is 3, it is shifted to one lower  l
Cu Updates
Cu   19 Aug 18 First cut at gradients
Cu   04 Jul 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,ndrphi,nrmx,nr,nsp,mode,nrphi(0:nl-1)
      double precision ri(nrmx),rwt(nrmx)
      real(8) :: phi1(nrmx,0:nl-1,ndrphi,nsp),phi2(nrmx,0:nl-1,ndrphi,nsp)
      real(8),target :: ovphi(ndrphi,ndrphi,0:nl-1,nsp)
      real(8) :: phi21(nr,0:nl-1,ndrphi,0:nl-1,ndrphi,nsp)
C ... Dynamically allocated arrays
      real(8),allocatable,target :: g1(:,:,:,:),g2(:,:,:,:)
      real(8),allocatable :: gx(:,:,:,:)
C     real(8),allocatable :: eb(:),ebx(:),z(:,:),ovb(:,:)
C ... Local parameters
      integer,parameter:: NULLI = -99999
      integer,parameter:: ivsn = 1
      integer,parameter:: ndif=8  ! number of points in num. diff of phi
      real(8),parameter:: ovtol = 1d-10
      integer isp,l1,l2,n1,n2,mode0,mode1,mode2,mode3,mode4,nspp,ir
      double precision ovv,fac,evl(ndrphi),fv1(ndrphi),fv2(ndrphi),z(ndrphi,ndrphi)
      procedure(integer) :: iprint
      procedure(real(8)) :: dot3

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      mode4 = mod(mode/10000,10)
      if (mode0 == 0) then
        mode = ivsn
        return
      endif

C --- Make local copy of phi, possibly spin-averaged, or "magnetic" phi ---
      nspp = nsp  ! Default number of spins in phi-phi products
      allocate(g1(nr,0:nl-1,ndrphi,nsp),g2(nr,0:nl-1,ndrphi,nsp))
      do  l1 = 0, nl-1
        do  n1 = 1, nrphi(l1)
          do  isp = 1, nsp
            call dcopy(nr,phi1(1,l1,n1,isp),1,g1(1,l1,n1,isp),1) ! Default
            call dcopy(nr,phi2(1,l1,n1,isp),1,g2(1,l1,n1,isp),1) ! Default
            if (mode1 == 2) then
              if (nsp == 1) call rx('no magnetic product basis with nsp=1')
              call rx('prodphi magnetic product basis not implemented yet')
            elseif (mod(mode1,2) > 0 .and. nsp == 2) then
              fac = dsign(1d0,phi1(2,l1,n1,1)*phi1(2,l1,n1,2))
              g1(1:nrmx,l1,n1,1) = phi1(1:nrmx,l1,n1,1)/2 + fac*phi1(1:nrmx,l1,n1,2)/2
              fac = dsign(1d0,phi2(2,l1,n1,1)*phi2(2,l1,n1,2))
              g2(1:nrmx,l1,n1,1) = phi2(1:nrmx,l1,n1,1)/2 + fac*phi2(1:nrmx,l1,n1,2)/2
              nspp = 1
              exit
            endif
          enddo
        enddo
      enddo

C --- Replace g1 with a gradient ---
C mch -vl=$l -f1p9g18.8 g1 -coll 1,l+2 -diff g1 -e2 'x1*0' '(l+1)*x{l+2}/x1' -- out.coo2 -coll 1,l+3 --
C mch -vl=$l -f1p9g18.8 g1 -coll 1,l+2 -diff g1 -e2 'x1*0' '(l)*x{l+2}/x1' -+  out.coo2 -inc 'x1>0' -coll 1,l+1 --
      if (mode2 > 0) then
        allocate(gx(nr,0:nl-1,ndrphi,nspp)); call dpzero(gx,size(gx))
        do  isp = 1, nspp
        do  l1 = 0, nl-1
        do  n1 = 1, nrphi(l1)
          call prodphig(mode2,l1,nr,ri,g1(1,l1,n1,isp),gx(1,l1,n1,isp))
        enddo
        enddo
        enddo
        call dcopy(size(gx),gx,1,g1,1)
        deallocate(gx)
      endif

C     call prrmsh('g1',ri,g1(1,0,1,1),nr,nr,nl)
C     call prrmsh('g1-2',ri,g1(1,0,2,1),nr,nr,nl)

C --- Replace g2 with a gradient ---
      if (mode3 > 0) then
        allocate(gx(nr,0:nl-1,ndrphi,nspp)); call dpzero(gx,size(gx))
        do  isp = 1, nspp
        do  l1 = 0, nl-1
        do  n1 = 1, nrphi(l1)
          call prodphig(mode3,l1,nr,ri,g2(1,l1,n1,isp),gx(1,l1,n1,isp))
        enddo
        enddo
        enddo
        call dcopy(size(gx),gx,1,g2,1)
        deallocate(gx)
      endif
C     call prrmsh('g2',ri,g2(1,0,1,1),nr,nr,nl)
C     call prrmsh('g2-2',ri,g2(1,0,2,1),nr,nr,nl)

C      print *, 'mode, sumcheck', mode, sum(g2), sum(g1)
C      print *, 'check g1(nr)', g1(nr,:,:,1)
C      print *, 'check g2(nr)', g2(nr,:,:,1)

C --- Replace g1 with phi1 and/or g2 with phi2 ---
      if (mode0 >= 4) then
        do  l1 = 0, nl-1
          do  n1 = 1, nrphi(l1)
            do  isp = 1, nsp
              if (mode2 > 0) call dcopy(nr,g1(1,l1,n1,isp),1,phi1(1,l1,n1,isp),1)
              if (mode3 > 0) call dcopy(nr,g2(1,l1,n1,isp),1,phi2(1,l1,n1,isp),1)
            enddo
          enddo
        enddo
      endif

C --- Overlap ---
      if (mod(mode0,2) == 0) goto 99

      call dpzero(ovphi,size(ovphi))
      do  isp = 1, nspp
      do  l1 = 0, nl-1
        do  n1 = 1, nrphi(l1)
        do  n2 = n1, nrphi(l1)
          ovv = dot3(nr,g1(1,l1,n1,isp),g2(1,l1,n2,isp),rwt)
          if (abs(ovv) < ovtol) ovv = 0
          if (mode4 == 4) ovv = ovv*sqrt(1d0+0.1d0*n1)*sqrt(1d0+0.1d0*n2)
          ovphi(n1,n2,l1,isp) = ovv
          ovphi(n2,n1,l1,isp) = ovv
        enddo
        enddo
      enddo
      enddo

C ... Printout overlap
      do  isp = 1, nspp
      if (iprint() >= 60) then
        call info2(60,1,0,' Overlap matrix of partial waves%?#(n>1)# (spin %i)##%N'//
     .    '  n1   n2   l=0        l=1 ...',nspp,isp)
        ir = maxval(nrphi)
        do  n1 = 1, ir
        do  n2 = n1, ir
          call info5(60,0,0,' %,3i  %,3i %n:2;9,9F',n1,n2,nl,ovphi(n1,n2,:,isp),5)
        enddo
        enddo
      endif
      if (iprint() >= 55) then
        call info2(55,1,0,' Normalization of partial waves <phi_n | phi_n>%?#(n>1)# (spin %i)##%N'//
     .    '   l   n=0        n=1 ...',nspp,isp)
        do  l1 = 0, nl-1
          forall (ir=1:nrphi(l1)) fv1(ir) = ovphi(ir,ir,l1,isp)
          call info5(55,0,0,' %,3i %n:2;9,9F',l1,nrphi(l1),fv1,4,5)
        enddo
      endif
      enddo

      if (mode4 == 0) goto 99

      do  isp = 1, nsp
      do  l1 = 0, nl-1

C   ... Diagonalize overlap matrix
C       Eigenvectors z of diagonalize overlap:
C       \sum_j ovphi_ij*z_jk = z_ik e_k with ovphi_ij = <phi_i | phi_j>
C       ovphi~_ij = (z+ ovphi z)_ij = e_i delta_ij
C       Then phi~_i = \sum_k phi_k z_ki / sqrt(e_i) forms orthonormal basis
        if (mode4 == 3 .or. mode4 == 4) then
          call dpzero(z,size(z))
          n1 = nrphi(l1)
          call rs(ndrphi,n1,ovphi(1,1,l1,isp),evl,1,z,fv1,fv2,ir)
          call rxx(ir/=0,'prodphi failed to diagonalize overlap matrix')
          do  ir = 1, n1
            call dscal(n1,1/sqrt(evl(ir)),z(1,ir),1)
          enddo
C         call yprmi('z for l=%i',l1,0,1,z,0,ndrphi,n1,n1)
          call dcopy(size(z),z,1,ovphi(1,1,l1,isp),1)

C   ... Cholesky decomposition of overlap matrix, or its inverse
        elseif (mode4 == 1 .or. mode4 == 2) then
          n1 = nrphi(l1)
C         call prmx('overlap',ovphi(1,1,l1,isp),ndrphi,n1,n1)
          call dschd(ndrphi,n1,ovphi(1,1,l1,isp),z,mode4==2,ir)
          call rxx(ir/=0,'gworthphi: error in dschd')
          do  n2 = 1, n1
            forall (l2=1:n2-1) ovphi(l2,n2,l1,isp) = 0
          enddo
C         call prmx('cd, L^-1',ovphi(1,1,l1,isp),ndrphi,n1,n1)
        endif

      enddo
      enddo

   99 continue ! End of overlap branch

C --- Make phi21 ---
      if (mod(mode0/2,2) == 0) goto 199

      do  isp = 1, nspp
      do  l1 = 0, nl-1
      do  n1 = 1, nrphi(l1)
        do  l2 = 0, nl-1
        do  n2 = 1, nrphi(l2)
          phi21(1,l1,n1,l2,n2,isp) = 0
          phi21(2:nr,l1,n1,l2,n2,isp) = g2(2:,l2,n2,isp)/ri(2:)*g1(2:,l1,n1,isp) ! phi21 = r * (phi2/r) * (phi1/r)
        enddo
        enddo
      enddo
      enddo
      enddo

  199 continue ! end of phi21 branch
      deallocate(g1,g2)

      end

      subroutine prodphig(mode,l,nr,ri,g,gradg)
C- Return some form of the radial gradient of g
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : 0 do nothing
Ci         : 1 gradg = dg/dr
Ci         : 2 gradg = dg/dr - (l+1)*g/r
Ci         : 3 gradg = dg/dr + l*g/r
Ci   nr    :number of radial mesh points
Ci   ri    :radial mesh
Ci   g     :normalized wave function times r
Co Outputs
Ci   gradg : gradient of g, according to mode
Cu Updates
Cu   18 Aug 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nr,l
      double precision ri(nr),g(nr),gradg(nr)
      integer,parameter:: ndif=8  ! number of points in num. diff of phi
C ... Local parameters
      integer ir
      double precision gg(nr)

      call poldvm(ri,g,nr,ndif,.false.,1d-8,ir,gg)

      select case (mode)
C     dg/dr
      case (1); call dcopy(nr,gg,1,gradg,1)
C     dg/dr - (l+1)/r*g -> g(l+1), but store in g(l)
      case (2)
        forall (ir = 2:nr) gg(ir) = gg(ir) - (l+1)/ri(ir)*g(ir)
        call dcopy(nr,gg,1,gradg,1)
C     dg/dr + l/r*g -> g(l-1), but store in g(l)
      case (3)
        forall (ir = 2:nr) gg(ir) = gg(ir) + l/ri(ir)*g(ir)
        call dcopy(nr,gg,1,gradg,1)
      case default; call rx('prodphig : bad mode')
      end select
      gradg(1) = (ri(3)*gradg(2)-ri(2)*gradg(3))/(ri(3)-ri(2))

      end
