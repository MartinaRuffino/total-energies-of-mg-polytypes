      subroutine desxdn(n,p0,dedv,vsxl)
C-  Make dE_sx/dn(R') from dE_sx/dv(R).
C ----------------------------------------------------------------------
Ci Inputs
Ci   p0,n   response function and dimension.  Left intact on exit
Ci   dedv   dE_sx/dv_R.  Left intact unless vsxl and dedv point to
Ci          same space
Co Outputs
Co   vsxl   dE_sx/dn_R ... may occupy same address space as dedv.
Cr Remarks
Cr   vsxl is obtained from inverting (Kotani, PRL 74, 2989 (95))
Cr   dE_sx/dv(R) = \sum_R' p0(R,R') dE_sx/dn(R').
Cr   Inverse of p0 by singular value decomposition, because p0 has
Cr   a 0 eigenvalue corresponding to a uniform shift in potential.
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n
      double precision p0(n,n),dedv(n),vsxl(n)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk(:),z(:),p0l(:),e(:)
C ... Local parameters

C ... Allocate some local arrays
      allocate(wk(n*11))
      allocate(z(n*n))
      allocate(p0l(n*n))
      allocate(e(n))
      call dcopy(n*n,p0,1,p0l,1)
      call pdesxd(n,p0l,wk,z,e,dedv,vsxl)
      deallocate(wk,z,p0l,e)

      end
      subroutine pdesxd(n,p0,wk,z,e,dedv,vsxl)
C- Kernel called by desxdn
      implicit none
      integer n
      double precision p0(n,n),dedv(n),vsxl(n),
     .  wk(n,11),z(n,n),e(n)
C Local variables
      integer k,ierr
C ... This is the smallness parameter for singular values
      double precision tol,emax
      parameter (tol=1d-10)

C      call prmx('p0',p0,n,n,n)
C      call prmx('dedv',dedv,n,n,1)

C ... Eigenvalues and eigenvectors of p0
*     call dsev1(n,p0,w,wk,0,.false.,.false.,.false.,n,9d9,nn,z,e)
      call dtridx(n,n,p0,wk,wk(1,4),wk(1,5),wk(1,2))
      call imtqlv(n,wk,wk(1,4),wk(1,5),e,wk(1,11),ierr,wk(1,6))
      call rxx(ierr /= 0,'imtqlv cannot find all evals')
      call tinvit(n,n,wk(1,1),wk(1,4),wk(1,5),n,e,wk(1,11),z,
     .  ierr,wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10))
      call rxx(ierr /= 0,'tinvit cannot find all evecs')
      call dtribx(n,n,p0,wk(1,2),n,z)

C ... Zero any eigenvalues less than tol*largest value
      emax = 0d0
      do  13  k = 1, n
   13 emax = max(emax,dabs(e(k)))
      do  14  k = 1, n
   14 if (dabs(e(k)) < emax*tol) e(k) = 0

C ... Solve p0 vsxl = dedv in restricted subspace
      call svbksb(n,n,n,1,e,z,z,dedv,vsxl,wk)

      end

