      subroutine mstokm(mode,nl,ncomp,norb,gms,gkm,idx)
C- Convert from lms to l,kappa,mu rep or back depending on mode
C-----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 transform from lms to kappa-mu
Ci         :1 transform from kappa-mu to lms
Ci         :2 don't transform, instead make the index array
Ci         :4 don't transform, instead return in gkm coefficients C
Ci            such that C^T g(lms) C = g(kmu)
Ci         :5 don't transform, instead return in gms coefficients C
Ci            such that C^T g(kmu) C = g(lms)
Ci   nl    :(maximum l) + 1
Ci   ncomp :number of functions to transform (needed for CPA)
Ci   norb  :number of orbitals, (lmax+1)^2
Cio Inputs/Outputs (not used if mode=2)
Cio  gms   :Matrix in lms rep (input for mode 0, output for mode 1)
Cio  gkm   :Matrix in kmu rep (input for mode 1, output for mode 0)
Co Outputs
Co   idx   :Index array idx(i,l,imu) = index(kappa_i), where i=1,2
Co         :idx(1,:) = 0 => 1 soln and (kap1,kap2) block reduces to (2,2) element
Co         :Output for mode 2 only
Cl Local  variables
Cl   mu2   :mu*2 (an integer)
Cu Updates
Cu  18 Jun 18 Synchronize with updated spherical harmonics
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,ncomp,norb,idx(2,0:nl-1,2*nl)
      double complex gms(norb*2,norb*2,ncomp),gkm(norb*2,norb*2,ncomp)
C ... Local parameters
      integer i,i0,icomp,ifac,im1,im2,k1,k2,l,l1,l2,lmax,m1,m2,n1,n2,n3,mode0
      integer lrec(norb),mrec(norb),ip(4,norb),kappa(norb*2),mu2(norb*2)
      double complex wk(norb*2,norb*2)
      double precision cl(norb*2,norb*2),u,c1,c2
      procedure(integer) :: ll

      mode0 = mod(mode,10)

      call sanrg(.true.,mode0,0,4,'mstokm','mode')
      if (mod(mode0,4) == 0) gkm = 0 ; if (mod(mode0,4) == 1) gms = 0

C     Order of orbitals in kappa-mu basis
C     #   kappa  mu    l
C     1    -1   -1/2   0
C     2    -1    1/2   0
C     3     1   -1/2   1
C     4     1    1/2   1
C     5    -2   -3/2   1
C     6    -2   -1/2   1
C     7    -2    1/2   1
C     8    -2    3/2   1
C     ...
C     Clebsch coupling between k1+k2=-1, mu1=mu2

C --- Build the transformation matrix from Clebsch coefficients
      lmax = ll(norb)
C ... Initialize index arrays
      kappa = 0 ; mu2 = 0 ; i0 = 0
      if (mode0 == 2) idx = 0
      do l = 0, lmax
        n1 = 2*l ; n2 = 2*(l+1)
        if (n1 /= 0) then
          do i = 1, n1
            i0 = i0 + 1
            kappa(i0) = l ; mu2(i0) = - 2*l - 1 + 2*i ! mu = -(l+1/2) when i=1
            if (mode0 == 2) idx(1,l,i+1) = i0
          enddo
        endif
        do i = 1, n2
          i0 = i0 + 1
          kappa(i0) = - l - 1 ; mu2(i0) = -2*l - 1 + 2*(i-1)
          if (mode0 == 2) idx(2,l,i) = i0
        enddo
      enddo

      if (mode0 == 2) return

C ... l and m for each orbital in order of appearance. Next section requires m ordered l:-l
      call lmorder(0,i,lrec,mrec)
      if (i == 2) i = 1; if (i == 0 .or. i == 3) i = 2
      call lmorder(3-i,ll(norb),lrec,mrec)

C ... Make a list of orbital pairs (same l, mu) or max mu for one l
      ip = 0 ; i = 0
      do n1 = 1, 2*norb - 1
        k1 = kappa(n1) ; im1 = mu2(n1)
        l1 = k1 ; if (k1 < 0) l1 = -k1 - 1
        ifac = 1 ; if (-2*k1 - 1 == abs(im1)) ifac = - 1 ! look for -mu for same l
        do n2 = n1+1, 2*norb
          k2 = kappa(n2) ; l2 = k2 ; if (k2 < 0) l2 = -k2 - 1
          im2 = mu2(n2)
          if (l2 == l1 .and. im1 == ifac*im2) then
            i = i + 1 ; ip(1,i) = n1 ; ip(2,i) = n2
c           Find two orbitals with same l, m+s = mu
            do n3 = 1, norb
              if (lrec(n3) /= l1) cycle
              if (mrec(n3) == (im1 + 1)/2) ip(3,i) = n3
              if (mrec(n3) == (im2 - 1)/2) ip(4,i) = n3 + norb
            enddo
            exit
          endif
        enddo
      enddo

C ... Assemble the trafo matrix CL (lms = CL * kmu) block by block
      cl = 0
      do i = 1, norb ! there are exactly norb orbital pairs
        n1 = ip(1,i) ; n2 = ip(2,i) ; m1 = ip(3,i) ; m2 = ip(4,i)
        im1 = mu2(n1) ; im2 = mu2(n2) ; l = lrec(m1)
        if (im1 + im2 == 0) then ! maximal abs(mu) pair
          if (mode0 == 0) then
            cl(m1,n1) = 1d0 ; cl(m2,n2) = 1d0
          elseif (mode0 == 1) then
            cl(n1,m1) = 1d0 ; cl(n2,m2) = 1d0
          endif
          cycle
        endif
        u = (im1/2d0)/(l+5d-1) ; c1 = sqrt((1+u)/2) ; c2 = sqrt((1-u)/2)
C       By construction of the pairs list, n1 is kappa = l; n2 is -l-1
        if (mode0 == 0) then
          cl(m1,n1) =   c1 ; cl(m1,n2) =   c2
          cl(m2,n1) = - c2 ; cl(m2,n2) =   c1
        elseif (mode0 == 1) then
          cl(n1,m1) =   c1 ; cl(n1,m2) = - c2
          cl(n2,m1) =   c2 ; cl(n2,m2) =   c1
        endif
      enddo

C     call prmx('kappa-mu Clebsh Gordan coeffs',cl,norb*2,norb*2,norb*2)

      if (mode0 == 4) then
        gkm(:,:,1) = cl(:,:)
        return
      endif
      if (mode0 == 5) then
        gms(:,:,1) = cl(:,:)
        return
      endif

C ... Perform the trafo
      do icomp = 1, ncomp
        if (mode0 == 0) wk(:,:) = gms(:,:,icomp)
        if (mode0 == 1) wk(:,:) = gkm(:,:,icomp)
C       print *, sum(wk)
        wk = matmul(transpose(cl),matmul(wk,cl))
C       print *, sum(wk)
        if (mode0 == 0) gkm(:,:,icomp) = wk(:,:)
        if (mode0 == 1) gms(:,:,icomp) = wk(:,:)
      enddo

      end
