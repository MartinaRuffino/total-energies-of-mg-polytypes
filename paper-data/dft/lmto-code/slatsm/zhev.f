C#define LAPACK
      subroutine zhevx(n,lh,h,s,lov,lx,nmx,emx,nev,wk,linv,e,lz,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   n:    order of h and s
Ci   lh:   leading dimension of h and s
Ci   h:    hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix, (used only if lov is true)
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci         (not used if LAPACK zhegv is invoked)
Ci   wk:   work array of length at least 11n
Ci         NB: If LAPACK version is used, and eigenvectors are sought
Ci         wk should be dimensioned max(11,2*n)*nmx
Ci   lov:  0 no overlap matrix
Ci         1 overlap matrix, return evecs of nonorthogonal H
Ci   lx:   if T, calls routines to exploit unit stride lengths (risc)
Ci         Not used if LAPACK zhegv is invoked.
Ci   linv: if T, using inverse iteration
Ci         Not used if LAPACK zhegv is invoked.
Ci   lz:   leading dimension of z
Co Outputs:
Co   e:    eigenvalues
Co   nev:  number of eigenvectors found
Co   z:    eigenvectors (1..nev)  (declared as z(n,*)
Co   s:    has been decomposed into and LL+ decomposition.
Co         You can call zhev2 to scale a vector by L
Cr Remarks:
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cr   Aborts on exit
Cp Procedures used:
Cp   (lapack)  zhegv
Cp   (eispack) htribk, htridx, imtql2, tqlrat
Cu Updates
Cu   24 Feb 07 Bug fix when nmx=0
Cu   17 May 03 Adapted from zhev, intended to supersede zhev.
Cu   14 Aug 02 Added zheev when lov is F; new zhev2.
Cu   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      logical linv,lx
      integer lov,lh,lz
      integer n,nev,nmx,nel
!       complex(8) :: h(lh,*),s(lh,*),z(lz,*)
      complex(8) :: h(*),s(*),z(*)
      real(8) :: e(n),wk(11*n),emx


C#ifdef LAPACK
C Local parameters
      character :: jobz
      integer :: lwork, info
      real(8) :: abstol, vu, vl
      integer, allocatable :: iwork(:), ifail(:)
      real(8), allocatable :: rwork(:)
      complex(8), allocatable :: work(:)


      call tcn('zhev')


      abstol = -1
      vl = 0
      vu = 0
      info = 0
      nev = 0
      lwork = -1
      jobz = 'v'
      if (nmx <= 0) jobz = 'n'
!! unievp often uses different lda and ldz, it seems this is not well tested in mkl and while it gives correct results it does trigger some uninitialised variable deep in mkl and crashes if in debug or optimised debug mode.
!       call unievp1(lov>0, .true., 1, jobz, 'a', 'u', n, h, lh, s, lh, vl, vu, 1, n, nel, nev, e, z, lz)

C     Workspace query
      allocate(work(1),rwork(7*n),iwork(5*n),ifail(n))
      work = 1
      if (lov > 0) then
        call zhegvx(1,jobz,'a','u',n,h,lh,s,lh,vl,vu,1,n,abstol,nev,e,z,lz,work,lwork,rwork,iwork,ifail,info)
      else
        call zheevx(jobz,'a','u',n,h,lh,vl,vu,1,n,abstol,nev,e,z,lz,work,lwork,rwork,iwork,ifail,info)
      end if

C     Diagonalization
      lwork = max(1,nint(real(work(1))))
      deallocate(work)
      allocate(work(lwork))
      if (lov > 0) then
        call zhegvx(1,jobz,'a','u',n,h,lh,s,lh,vl,vu,1,n,abstol,nev,e,z,lz,work,lwork,rwork,iwork,ifail,info)
      else
        call zheevx(jobz,'a','u',n,h,lh,vl,vu,1,n,abstol,nev,e,z,lz,work,lwork,rwork,iwork,ifail,info)
      end if

      deallocate(work,rwork,iwork,ifail)

      if (info > 0) write(*,'("zhevx: ",i0," vector(s) failed to converge. Try setting abstol=2*dlamch(s).")') info
      if (info /= 0) call rxi('zhevx: exited with nonzero status', info)

C#elseC
CC Local parameters
C      integer i,j,k,m,mdim,k1,k2,k3,k4,k5,k6,k7,k8,k9,n2,ier
C
C      call tcn('zhev')
C      n2 = n**2
C
CC --- Take care of n=1 case ---
C      if (n == 1) then
C         e(1) = h(1)
C         if (lov > 0) e(1) = h(1)/s(1)
C         if (1 > nev .or. e(1) > emx) return
C         z(1) = 1
C         z(2) = 0
C         return
C      endif
C
CC --- Separate real and imaginary parts ---
C      mdim = 2*lh
C      call ztoy(h,lh,n,n,0)
C      if (lov > 0) call ztoy(s,lh,n,n,0)
C
CC --- Debugging: eigenvalues of overlap ---
CC      call htridx(mdim,n,s(1),s(lh+1),e,wk,wk(3*n+1),wk(n+1))
CC      do  11  j = 1, n
CC   11 wk(j) = wk(j)**2
CC      call tqlrat(n,e,wk,ier)
CC      write(6,600) e
CC  600 format(' evl='/(1x,1p,5e14.6))
CC      call rx('eigenvalues of overlap')
C
CC --- H <- S^-1/2  H  S^-1/2 ---
C      if (lov > 0) then
C        call yyhchd(mdim,n,s,s(lh+1),wk,lx,.true.,ier)
C        call rxx(ier /= 0,'ZHEV: error in yyhchd')
C        if (lx) then
C          if (lz /= lh) call rx('zhev not ready for lz ne lh')
C          call yyhrdx(mdim,n,h,h(lh+1),s,s(lh+1),z,z(lz+1))
C        else
C          call yyhred(mdim,n,h,h(lh+1),s,s(lh+1),.true.)
C        endif
C      endif
C
CC --- Transform to tridiagonal matrix ---
C      if (linv) then
C        k1 = 1
C        k2 = k1 + 3*n
C        k3 = k2 + n
C        k4 = k3 + n
C        k5 = k4 + n
C        k6 = k5 + n
C        k7 = k6 + n
C        k8 = k7 + n
C        k9 = k8 + n
C#ifndefC GENERIC
C        call htridx(mdim,n,h(1),h(lh+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C#elseC
C        call htridi(mdim,n,h(1),h(lh+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C#endifC
C      else
C#ifndefC GENERIC
C        call htridx(mdim,n,h(1),h(lh+1),e,wk,wk(3*n+1),wk(n+1))
C#elseC
C        call htridi(mdim,n,h(1),h(lh+1),e,wk,wk(3*n+1),wk(n+1))
C#endifC
C      endif
C
CC --- Eigenvalues only ---
C      if (nmx <= 0) then
C        do  12  j = 1, n
C   12   wk(j) = wk(j)**2
C        call tqlrat(n,e,wk,ier)
C        call rxx(ier /= 0,'ZHEV: tqlrat cannot find all evals')
C        nev = 0
C        goto 100
C
CC --- Eigenvalues and eigenvectors ---
C      else if (linv) then
C        call imtqlv(n,wk(k1),wk(k2),wk(k3),e,wk(k9),ier,wk(k4))
C        call rxx(ier /= 0,'ZHEV: imtqlv cannot find all evals')
CC   --- Determine number of eigenvectors to be calculated ---
C        nev = 1
C        do  14  j = 2, n
C          if (j <= nmx .and. e(j-1) <= emx) nev = j
C   14   continue
C        call tinvit(2*lz,n,wk(k1),wk(k2),wk(k3),nev,e,wk(k9),z,ier,
C     .    wk(k4),wk(k5),wk(k6),wk(k7),wk(k8))
C        call rxx(ier /= 0,'ZHEV: tinvit cannot find all evecs')
C      else
C        do  17  j = 1, n
C          k = (j-1)*2*lz
C          m = k+n
C          do  16  i = k+1, m
C   16     z(i) = 0d0
C          z(k+j) = 1d0
C   17   continue
C        call imtql2(2*lz,n,e,wk,z,ier)
CC       call prmx('eval',e,n,n,1)
C        call rxx(ier /= 0,'ZHEV: imtql2 cannot find all evecs')
C
CC   --- Determine number of eigenvectors to be calculated ---
C        nev = 1
C        do  15  j = 2, n
C          if (j <= nmx .and. e(j-1) <= emx) nev = j
C   15   continue
C      endif
CC     call prmx('eval',e,n,n,1)
C
C      if (nev > 0) then
C#ifndefC GENERIC
C        if (lz /= lh) call rx('zhev not ready for lz ne lh')
C        call htribx(mdim,n,h(1),h(lh+1),wk(n+1),nev,z(1),z(lz+1))
C#elseC
C        call htribk(mdim,n,h(1),h(lh+1),wk(n+1),nev,z(1),z(lz+1))
C#endifC
C
CC --- Get the eigenvectors of H - E O ---
C        if (lov > 0) then
C          if (lx) then
C            call ymcpy(z,2*lz,1,lz,h,mdim,1,lh,n,n)
C            call yympy(s,s(lh+1),mdim,1,h,h(lh+1),mdim,1,z,z(lz+1),2*lz,
C     .        1,n,nev,n)
C          else
C            call yyhbak(mdim,n,s,s(lh+1),nev,z,z(lz+1),.true.)
C          endif
C        endif
C
CC   --- Convert eigenvectors to double complex storage ---
C        call ztoy(z,lz,n,nev,1)
CC       call zprm('evecs',2,z,lz,n,nmx)
C      endif
C
C#endif

  100 call tcx('zhev')
      end

      subroutine zhev(n,h,s,lov,lx,nmx,emx,nev,wk,linv,ltime,e,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   n:    dimension of h
Ci   h,n:  hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix, (used only if lov is true)
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci   wk:   work array of length at least 11n
Ci   lov:  if T, non-orthogonal
Ci   lx:   if T, calls routines to exploit unit stride lengths (risc)
Ci         Not used if LAPACK zhegv is invoked.
Ci   linv: if T, using inverse iteration
Ci         Not used if LAPACK zhegv is invoked.
Ci   ltime:Not used now
Co Outputs:
Co   e:    eigenvalues
Co   nev:  number of eigenvectors found
Co   z:    eigenvectors (1..nev)  (declared as z(n,*)
Co   s:    has been decomposed into and LL+ decomposition.
Co         You can call zhev2 to scale a vector by L
Cb Bugs:
Cb    The gnu library seems to require 4*N for RWORK in the call to zhegv.
Cb    Workaround: allocate local array work; add extra memory to speed it up.
Cr Remarks:
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cr   Aborts on exit
Cp Procedures used:
Cp   (lapack)  zhegv
Cp   (eispack) htribk, htridx, imtql2, tqlrat
Cu Updates
Cu   14 Aug 02 Added zheev when lov is F; new zhev2.
Cu   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      logical linv,lx
      integer n,nev,nmx,ltime
      double precision h(1),s(1),e(n),wk(11*n),z(2),emx
      logical lov
      complex(8), allocatable :: work(:)

C#ifdef LAPACK
C Local parameters
      integer ier,lwork
      character jobz

      call tcn('zhev')
      allocate(work(max(n*(n+1),4)/2))
      lwork = size(work)
      if (nmx <= 0) then
        jobz = 'N'
C       lwork = 4*n
        if (lov) then
C        call zhegv(1,jobz,'U',n,h,n,s,n,e,wwk,lwork,wk(1),ier)
         call zhegv(1,jobz,'U',n,h,n,s,n,e,work,lwork,wk(1),ier)
        else
          call zheev(jobz,'U',n,h,n,e,work,lwork,wk(1),ier)
        endif
        nev = 0
      else
        jobz = 'V'
C       lwork = n*min(n,nmx)
        if (lov) then
          call zhegv(1,jobz,'U',n,h,n,s,n,e,z,lwork,wk(1),ier)
        else
          call zheev(jobz,'U',n,h,n,e,z,lwork,wk(1),ier)
        endif
        call zcopy(n*min(n,nmx),h,1,z,1)
        nev = min(n,nmx)
      endif
      deallocate(work)
      call rxx(ier /= 0,'zhev: zhegv cannot find all evals')
C#elseC
CC Local parameters
C      integer i,j,k,m,mdim,k1,k2,k3,k4,k5,k6,k7,k8,k9,n2,ier
C
CC     call prm('(5f20.15)',s,n,n)
C      call tcn('zhev')
C      n2 = n**2
C
CC --- Take care of n=1 case ---
C      if (n == 1) then
C         e(1) = h(1)
C         if (lov) e(1) = h(1)/s(1)
C         if (1 > nev .or. e(1) > emx) return
C         z(1) = 1
C         z(2) = 0
C         return
C      endif
C
CCC --- Use diagno if lx ---
CC      if (lx) then
CC        call cplx2r(n2,0,h,z)
CC        if (lov) call cplx2r(n2,0,s,z)
CC        i = 0
CC        if (lov) i = 1
CC        call diagno(n,h,s,wk,lx,i,linv,nmx,emx,nev,z,e)
CC        call cplx2r(n2,1,z,h)
CC        goto 40
CC      endif
C
CC --- Separate real and imaginary parts ---
C      mdim = 2*n
C      call ztoy(h,n,n,n,0)
C      if (lov) call ztoy(s,n,n,n,0)
CC      mdim = 2*n
CC      do  10  j = 1, n
CC        k = (j-1)*mdim + 1
CC        l = k + n
CC        call dcopy(n,h(k+1),2,wk,1)
CC        call dcopy(n,h(k),2,h(k),1)
CC        call dcopy(n,wk,1,h(l),1)
CC        if (lov) then
CC          call dcopy(n,s(k+1),2,wk,1)
CC          call dcopy(n,s(k),2,s(k),1)
CC          call dcopy(n,wk,1,s(l),1)
CC        endif
CC   10 continue
C
C
CC --- Debugging: eigenvalues of overlap ---
CC      call htridx(mdim,n,s(1),s(n+1),e,wk,wk(3*n+1),wk(n+1))
CC      do  11  j = 1, n
CC   11 wk(j) = wk(j)**2
CC      call tqlrat(n,e,wk,ier)
CC      write(6,600) e
CC  600 format(' evl='/(1x,1p,5e14.6))
CC      call rx('eigenvalues of overlap')
C
CC --- H <- S^-1/2  H  S^-1/2 ---
C      if (lov) then
C        call yyhchd(mdim,n,s,s(n+1),wk,lx,.true.,ier)
C        call rxx(ier /= 0,'ZHEV: error in yyhchd')
C        if (lx) then
C          call yyhrdx(mdim,n,h,h(n+1),s,s(n+1),z,z(n+1))
C        else
C          call yyhred(mdim,n,h,h(n+1),s,s(n+1),.true.)
C        endif
C      endif
C
CC --- Transform to tridiagonal matrix ---
C      if (linv) then
C        k1 = 1
C        k2 = k1 + 3*n
C        k3 = k2 + n
C        k4 = k3 + n
C        k5 = k4 + n
C        k6 = k5 + n
C        k7 = k6 + n
C        k8 = k7 + n
C        k9 = k8 + n
C#ifndefC GENERIC
C        call htridx(mdim,n,h(1),h(n+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C#elseC
C        call htridi(mdim,n,h(1),h(n+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C#endifC
C      else
C#ifndefC GENERIC
C        call htridx(mdim,n,h(1),h(n+1),e,wk,wk(3*n+1),wk(n+1))
C#elseC
C        call htridi(mdim,n,h(1),h(n+1),e,wk,wk(3*n+1),wk(n+1))
C#endifC
C      endif
C
CC --- Eigenvalues only ---
C      if (nmx <= 0) then
C        do  12  j = 1, n
C   12   wk(j) = wk(j)**2
C        call tqlrat(n,e,wk,ier)
C        call rxx(ier /= 0,'ZHEV: tqlrat cannot find all evals')
C        nev = 0
C        goto 100
C
CC --- Eigenvalues and eigenvectors ---
C      else if (linv) then
C        call imtqlv(n,wk(k1),wk(k2),wk(k3),e,wk(k9),ier,wk(k4))
C        call rxx(ier /= 0,'ZHEV: imtqlv cannot find all evals')
CC   --- Determine number of eigenvectors to be calculated ---
C        nev = 1
C        do  14  j = 2, n
C          if (j <= nmx .and. e(j-1) <= emx) nev = j
C   14   continue
C        call tinvit(mdim,n,wk(k1),wk(k2),wk(k3),nev,e,wk(k9),z,ier,
C     .    wk(k4),wk(k5),wk(k6),wk(k7),wk(k8))
C        call rxx(ier /= 0,'ZHEV: tinvit cannot find all evecs')
C      else
C        do  17  j = 1, n
C          k = (j-1)*mdim
C          m = k+n
C          do  16  i = k+1, m
C   16     z(i) = 0d0
C          z(k+j) = 1d0
C   17   continue
C        call imtql2(mdim,n,e,wk,z,ier)
C        call rxx(ier /= 0,'ZHEV: imtql2 cannot find all evecs')
C
CC   --- Determine number of eigenvectors to be calculated ---
C        nev = 1
C        do  15  j = 2, n
C          if (j <= nmx .and. e(j-1) <= emx) nev = j
C   15   continue
C      endif
C
C      if (nev > 0) then
C#ifndefC GENERIC
C        call htribx(mdim,n,h(1),h(n+1),wk(n+1),nev,z(1),z(n+1))
C#elseC
C        call htribk(mdim,n,h(1),h(n+1),wk(n+1),nev,z(1),z(n+1))
C#endifC
C
CC --- Get the eigenvectors of H - E O ---
C        if (lov) then
C          if (lx) then
C            call dcopy(n2*2,z,1,h,1)
C            call yympy(s,s(n+1),mdim,1,h,h(n+1),mdim,1,z,z(n+1),mdim,1,
C     .        n,nev,n)
C          else
C            call yyhbak(mdim,n,s,s(n+1),nev,z,z(n+1),.true.)
C          endif
C        endif
C
CC   --- Convert eigenvectors to double complex storage ---
C        call ztoy(z,n,n,nev,1)
CC       call zprm('evecs',2,z,lz,n,min(nev,nmx))
C      endif
C
CC  40 if (nev > 0 .and. ltime >= 0) print 337, nev,nmx,n,emx
CC  337 format(' nev, nevmx, ndim=',3i4,' emx=',f10.5)
C#endif

  100 call tcx('zhev')
      end

      subroutine zhev2(mode,ldhs,n,m,lx,h,s)
C- Scales a hermitian matrix h by Cholesky decomposition of overlap s
C ----------------------------------------------------------------
Ci Inputs:
Ci   mode :For all modes, see Remarks.
Ci        :0 do nothing to h
Ci        :2  replace h to orthogonal form L^-1 h L+^-1 or U+^-1 h U^-1
Ci        :   h is assumed to be hermitian
Ci        :4  replace h to orthogonal by L h L+ or U+ h U
Ci        :   Modes 2 and 4 perform inverse functions of each other
Ci        :10 replace h to U^-1 h
Ci        :20 replace h to U h
Ci        :   Modes 10 and 20 perform inverse functions of each other
Ci        :
Ci        :Adding 1 to any of the following modes causes zhev2 to
Ci        :Cholesky-decompose hermitian matrix s.
Ci        :h is not used in the decomposition of s.
Ci        :
Ci   n    :dimension of h and s
Ci   m    :second dimension of h (mode>=10)
Ci   lx   :if T, calls custom routines for Cholesky decomp'sn
Ci        :that exploit unit stride lengths
Ci        :(not used in LAPACK implementation)
Co Inputs/Outputs:
Cio  s    :Overlap matrix.  s should be in Cholesky-decomposed form as
Cio       :computed by zhev or zhev2; or else 1 should be added to mode
Cio       :to decompose s.  s is unchanged apart from possible decomps'n
Cio  h    :On input:
Cio       :  h is a hermitian matrix of dimension (n,n) for mode<10
Cio       :  h is a matrix of dimension (n,m) for mode>=10
Cio       :On output, h is scaled depending on mode (or mode-1 if s is C.D.)
Cio       :  mode           h is transformed into:
Cio       :   0             h is unchanged
Cio       :   2             L^-1 h L+^-1 or U+^-1 h U^-1
Cio       :   4             L h L+       or U+ h U
Cio       :  10             U^-1 h
Cr Remarks:
Cr   zhev2 performs one of several operations that deal with overlap
Cr   matrices in a generalized eigenvalue problem.  A hermitian overlap
Cr   matrix S can be decomposed into L L+ (LAPACK implementation: U+ U)
Cr   Then the eigenvalues e of (H - e S) are the eigenvalues of
Cr     H' = L^-1 H L+^-1 (H' = U+^-1 H U^-1)
Cr   and eigenvectors Z of (H - e S) are related to vectors Z' of H' as
Cr     Z = L+^-1 Z' (Z = U^-1 Z')
Cu Updates
Cu   13 Aug 02 first created
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer mode,ldhs,n,m
      double precision h(1),s(1)
      logical lx
C Local parameters
      integer ier,itype,lmode
      double complex one
      parameter (one=(1.0d+0,0.0d+0))

      call tcn('zhev2')

      lmode = mode

C#ifndefC LAPACK
C      call rx('zhev2:  only implemented in LAPACK mode, sorry')
C#endif

C --- Cholesky decompose s ---
      if (mod(lmode,2) /= 0) then
        call zpotrf('U',n,s,ldhs,ier)
        if (ier /= 0) then
          call rxi('zhev2: failed to c.d. s: ier = ',ier)
        endif
        lmode = mode-1
      endif

C --- do nothing to h ---
      if (lmode == 0) then
        continue

C --- h <- U+^-1 h U^-1 ---
      elseif (lmode == 2) then
        itype = 1
        call zhegst(itype,'U',n,h,ldhs,s,ldhs,ier)
        call rxx(ier /= 0,'zhev: zhegv cannot find all evals')
        call z2herm('U',ldhs,n,h)

C --- h <- U+ h U ---
      elseif (lmode == 4) then
C       Overwrite h with U+ h
        call ztrmm('L','U','C','N',n,n,one,s,ldhs,h,ldhs)
C       Overwrite U+ h with U+ h U
        call ztrmm('R','U','N','N',n,n,one,s,ldhs,h,ldhs)
C       print '(a)', '#U+ h U'
C       call zprm('(5f12.6)',h,n,n)

C --- h <- U^-1 h ---
      elseif (lmode == 10) then
        call ztrsm('L','U','N','N',n,m,one,s,ldhs,h,ldhs)

C --- h <- U h ---
      elseif (lmode == 20) then
        call ztrmm('L','U','N','N',n,m,one,s,ldhs,h,ldhs)

      else
        call rxi('zhev2: bad mode',mode)
      endif


      call tcx('zhev2')
      end

      subroutine zhev3(mode,n,h,s,z)
C- Scales a hermitian matrix h by evec decomposition of overlap s
C ----------------------------------------------------------------
Ci Inputs:
Ci   mode :For all modes, see Remarks.
Ci        :0 do nothing to h
Ci        :2  replace h to orthogonal form L^-1 h L+^-1 or U+^-1 h U^-1
Ci        :   h is assumed to be hermitian
Ci        :4  replace h to orthogonal by L h L+ or U+ h U
Ci        :   Modes 2 and 4 perform inverse functions of each other
Ci        :10 replace h to U^-1 h
Ci        :20 replace h to U h
Ci        :   Modes 10 and 20 perform inverse functions of each other
Ci        :
Ci        :Adding 1 to any of the following modes causes zhev3 to
Ci        :Cholesky-decompose hermitian matrix s.
Ci        :h is not used in the decomposition of s.
Ci        :
Ci   n    :dimension of h and s
Ci   m    :second dimension of h (mode>=10)
Ci   lx   :if T, calls custom routines for Cholesky decomp'sn
Ci        :that exploit unit stride lengths
Ci        :(not used in LAPACK implementation)
Co Inputs/Outputs:
Cio  s    :Overlap matrix.  s should be in Cholesky-decomposed form as
Cio       :computed by zhev or zhev3; or else 1 should be added to mode
Cio       :to decompose s.  s is unchanged apart from possible decomps'n
Cio  h    :On input:
Cio       :  h is a hermitian matrix of dimension (n,n) for mode<10
Cio       :  h is a matrix of dimension (n,m) for mode>=10
Cio       :On output, h is scaled depending on mode (or mode-1 if s is C.D.)
Cio       :  mode           h is transformed into:
Cio       :   0             h is unchanged
Cio       :   2             L^-1 h L+^-1 or U+^-1 h U^-1
Cio       :   4             L h L+       or U+ h U
Cio       :  10             U^-1 h
Cr Remarks:
Cr   zhev3 performs one of several operations that deal with overlap
Cr   matrices in a generalized eigenvalue problem.  A hermitian overlap
Cr   matrix S can be decomposed into L L+ (LAPACK implementation: U+ U)
Cr   Then the eigenvalues e of (H - e S) are the eigenvalues of
Cr     H' = L^-1 H L+^-1 (H' = U+^-1 H U^-1)
Cr   and eigenvectors Z of (H - e S) are related to vectors Z' of H' as
Cr     Z = L+^-1 Z' (Z = U^-1 Z')
Cu Updates
Cu   13 Aug 02 first created
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer mode,n
      double complex h(n,n),s(n,n),z(n,n)
C Local parameters
      real(8),allocatable:: e(:),wk(:)
      logical lx
      integer i,lmode
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      parameter(lx=.true.)

      call tcn('zhev3')

C#ifndefC LAPACK
C      call rx('zhev3:  only implemented in LAPACK mode, sorry')
C#endif

      lmode = mode

C --- Eigenvectors and eigenvalues of s ---
      if (mod(mode,2) /= 0) then
      allocate(e(n),wk(11*n))
      call zhev(n,s,s,.false.,lx,n,9d9,i,wk,.false.,0,e,z)
C      if (lmode == 2) then
C      do  i = 1, n
C        do  j = 1, n
C          z(i,j) = z(i,j) / sqrt(e(j))
C        enddo
C      enddo
C      endif
C      if (lmode == 4) then
C      do  i = 1, n
C        do  j = 1, n
C          z(i,j) = sqrt(e(i)) * z(i,j)
C        enddo
C      enddo
C      endif
      deallocate(e,wk)
      lmode = mode-1
      endif

C --- do nothing to h ---
      if (lmode == 0) then
        continue

C --- h <- z+ h z ---
      elseif (lmode == 2) then
        call zgemm('N','N',n,n,n,one,h,n,z,n,zer,s,n)
        call zgemm('C','N',n,n,n,one,z,n,s,n,zer,h,n)
C       call zprm('z+ h z',2,h,n,n,n)

C --- h <- z h z+ ---
      elseif (lmode == 4) then
        call zgemm('N','C',n,n,n,one,h,n,z,n,zer,s,n)
        call zgemm('N','N',n,n,n,one,z,n,s,n,zer,h,n)
C       call zprm('z h z+',2,h,n,n,n)

      else
        call rxi('zhev3: bad mode',mode)

      endif

      call tcx('zhev3')
      end

      subroutine zhevo(n,lh,h,s,nmx,emx,epsovl,ovncut,nevl,nev,e,eo,lz,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix in
C- reduced hilbert space, eliminating projections onto small overlap
C ----------------------------------------------------------------
Ci Inputs:
Ci   n    :order of h and s
Ci   lh   :leading dimension of h and s
Ci   nmx  :maximum number of eigenvectors to be found
Ci   emx  :No longer used
Ci  epsovl:lower bound to eigenvalue of overlap
Ci  ovncut:discard or modify hilbert space belonging to lowest evals of s.
Ci        :See Remarks for description of how hilbert space is modified.
Ci         ovncut<=0 and epsovl<=0 : no special treatment.
Ci         ovncut=0 and epsovl>0 : modify space for which evals of s<epsovl
Ci         ovncut>0: modify space spanning first ovncut evals of s
Ci         Note: If ovncut>0 but epsovl=0, an internal default of epsovl
Ci               will be used (10**-12 * largest eval of s)
Ci   lz:   leading dimension of z
Cio Inputs/Outputs:
Cio  nevl: On input, governs manner in which hilbert space is modified
Cio      :   nevl=0: Reduce hilbert space (see Remarks) and return evecs
Cio      :           for reduced space
Cio      :   nevl=1: Modify s by replacing evals of diagonalized s
Cio      :           with epsovl (if nonzero), or fac*(input eo(1)) if epsovl=0
Cio      :           The modified s is rotated back to the original space
Cio      :           and the original eigenvalues returned in eo.
Cio      :           nmx, nev, e, h, z are not addressed
Cio      :   nevl=2: Initial steps are same as for nevl=1.
Cio      :           Generalized eigenvalues of (h,modified s) are returned
Cio      :           and h,s are destroyed.
Cio      :   nevl<0: Reduce hilbert space, then artificially restore it
Cio      :           to dimension n by assigning eigenvalues of the
Cio      :           eliminated part of s to epsovl (see Remarks)
Cio      : On output:
Cio      :  If input nevl=1, unchanged on exit
Cio      :  Otherwise, reduced dimension of hilbert space
Cio  eo  : Not addressed unless s modified: ovncut>0 or epsovl>0.
Cio      : If input nevl=0 => not used.  Otherwise:
Cio      : ecap=input eo(1): used to artifically set evals outside reduced part
Cio      : Case input nevl<0 and i>nevl: e(i) = ecap
Cio      : Case input nevl>1 and i>nevl: e(i) = min(ecap,e(i))
Cio      : On output:
Cio      : unused if input nevl=1 or no special treatment (ovncut=0 and epsovl<=0)
Cio      : Otherwise eo = eigenvalues of overlap matrix
Cio  h:  : Case nevl=1: h is never addressed.
Cio      : On input, hermitian matrix to be diagonalized.
Cio      : On output: h is destroyed
Cio  s:  : On input, hermitian overlap matrix.
Cio      : On output:
Cio      : Case nevl=1: s is modified to be positive definite (see Remarks)
Cio      : Otherwise s is destroyed
Co Outputs:
Co       : the following are generated if nevl<>1
Co   nev : number of eigenvectors found
Co   e   : eigenvalues of reduced or modified (h,s)
Co   z   : eigenvectors 1..nev of reduced or modified (h,s)
Cr Remarks:
Cr   Case input nevl > 1:
Cr     S is replaced by a positive definite matrix S'.
Cr     Let S = U+ d U where d are the eigenvalues of S.
Cr     Then S' = U+ d' U where d are modified eigenvalues of S:
Cr     d'(i) = max(d(i),epsovl)
Cr   Case input nevl < 0:
Cr   zhevo works as follows. Let T = evecs of S, such that
Cr     T+ S T = es
Cr   (es = eigenvalues of S).  For those states j for which es es_j > epsovl
Cr   (there are nevl such values, and ni-1 excluded values).
Cr   Define a matrix
Cr     T'_ij = T_ij / es_(j+ni-1), j=1,nevl
Cr   so that
Cr     T'+ S T = 1
Cr   Make a reduced hamiltonian
Cr     H' = T'+ H T'
Cr   H' has nevl eigenvalues E' and eigenvectors Z'. zhevo returns nevl
Cr   eigenvalues of H'.  When eigenvectors are calculated, zhev returns
Cr     Z = T' Z'
Cr
Cr * Case input nevl < 0:
Cr   1. the eigenvalues of the excluded hilbert space are assigned a constant
Cr   eigenvalue (whatever is supplied in input eo(1)).
Cr
Cr   2. If eigenvectors are also calculated, vectors z are also returned
Cr   for the excluded parts of the hilbert space.  That is to say, T'
Cr   is once again made into square matrix by appending eigenvectors
Cr   to T' for j = nevl+1 .. n
Cr     T'_i,j = T_i,(j+ni-1) / es_(j+ni-1), j=1,nevl   as before
Cr     T'_i,j = T_i,(nevl+j  / ecap,     j=nevl+1, n
Cr   The eigenvectors are orthogonal in the original sense, i.e.
Cr     I' = Z+ S Z
Cr   where I' is a diagonal matrix.  I' the unit matrix only in the
Cr   unexcluded part of the hilbert space.
Cr
Cp Procedures used:
Cp   (lapack) zheev
Cp   (blas)   zgemm
Cb Bugs
Cb   Probably possible to improve memory usage.
Cb   Replace call to zheev with call to zhevx
Cu Updates
Cu   18 Apr 12 New options, nevl=1 or 2
Cu   14 Apr 09 New ovncut
Cu   19 Jan 09 Input nevl < 0 generates evecs, large eval
Cu             for excluded part of hilbert space.
Cu   8  Jul 08 (T. Kotani) first implementation, adapted from zhev.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer lh,lz
      integer n,nevl,nev,nmx,ovncut
      double complex h(lh,n),s(lh,n),z(lz,n)
      double precision e(n),eo(n),emx,epsovl
C Local parameters
      logical lfullh
      integer ier,lwork,i,ni
      character jobz
      double precision ecap,xx,epsovli
      complex(8),allocatable:: wk(:),hh(:,:)
      real(8),allocatable:: rwork(:)
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))

C ... No special treatment: just call standard diagonalizer
      if (ovncut == 0 .and. epsovl <= 0) then
        lwork = n*n
        allocate(wk(lwork))
        call zhevx(n,lh,h,s,1,.true.,n,9d19,nev,wk,.false.,e,lz,z)
        return
      endif

      call tcn('zhevo')

C      call zprm('hamiltonian',2,h,lh,n,n)
C      call zprm('overlap',2,s,lh,n,n)

C     Preserve input eo(1), since eo will be overwritten
      ecap = eo(1)

C ... Eigenvalues and eigenvectors of s
      jobz = 'V'
      lwork = n*n
      allocate(wk(lwork),rwork(max(1,3*n-2)))
C     case nevl>0: preserve s for future use
      if (nevl > 0) then
C       call zprm('s',2,s,lh,n,n)
        allocate(hh(n,n))
        call zmcpy('N',s,lh,1,hh,n,1,n,n)
C       call zprm('s',2,hh,n,n,n)
C       hh is overwritten with evecs of s; eo holds evals of s
        call zheev(jobz,'U',n,hh,n,eo,wk,lwork,rwork,ier)
C     Case nevl<=0: s is overwritten with evecs of s; eo holds evals of s
      else
        call zheev(jobz,'U',n,s,lh,eo,wk,lwork,rwork,ier)
      endif
      deallocate(wk,rwork)

C ... Internal value of epsovl
      epsovli = epsovl
      if (epsovl <= 0) epsovli = eo(n)*1d-12

C ... Find 1st eigenvalue ni of S exceeding threshold
      if (ovncut > 0) then
        ni = ovncut+1
      else
        do  i = 1, n
          if (eo(i) > epsovli) then
            ni = i
            exit
          endif
        enddo
      endif

C --- Case nevl>0: generate modified s ---
      if (nevl > 0) then
C       Modify s unless no eigenvalues were changed
        if (ni /= 1) then
          if (epsovli <= 0)
     .      call rx1('zhevo: epsovli must be >0 but is %g',epsovli)
          allocate(rwork(n))
          call dcopy(n,eo,1,rwork,1)
          forall (i = 1:ni-1) rwork(i) = epsovli
          call zhdmpy(n,n,hh,n,rwork,(0d0,0d0),s,lh)
          deallocate(rwork)
        endif
        deallocate(hh)
C       Case nevl=1 : return with modified s
C       call zprm('modified s',2,s,lh,n,n)
        if (nevl == 1) return
        allocate(wk(lwork))
C       call zprm('h',12,h,lh,n,n)
        call zhevx(n,lh,h,s,1,.true.,n,9d19,nev,wk,.false.,e,lz,z)
C       call zprm('z',2,z,lz,n,n)
        deallocate(wk)
        nevl = n-ni+1
        if (ecap > 0) then
          do  i = ni, n
            e(i) = min(e(i),ecap)
          enddo
        endif
        return
      endif

C --- Case nevl<=0 ---
      lfullh = nevl < 0
      if (ni == 1) lfullh = .false.
C     If user doesn't supply ecap, use epsovli for it
      if (ecap == 0) ecap = 1/epsovli

C ... Map s to reduced Hilbert space (dimension = nevl)
C     For now, use temporary zz.  Better: use s in place of zz.
      nevl = n-ni+1

C ... Construct projection matrix
      allocate(hh(n,n))
      do  i = ni, n
        hh(:,i-ni+1) = s(:,i)/sqrt(eo(i))
      enddo
C     Put rejected part of hilbert space at the end; use epsovli for eval of s
      xx = sqrt(epsovli)
      do  i = 1, ni-1
        hh(:,nevl+i) = s(:,i)/xx
      enddo
      do  i = 1, n
        s(:,i) = hh(:,i)
      enddo
      deallocate(hh)

C ... Hamiltonian  <zz | H | zz>
      allocate(hh(nevl,n))
C     call zprm('h',2,h,lh,n,n)
C     call zprm('zz',2,s,lh,n,nevl)
      call zgemm('C','N',nevl,n,n,one,s,lh,h,lh,zer,hh,nevl)
C     call zprm('zz+*h',2,hh,nevl,nevl,n)
      call zgemm('N','N',nevl,nevl,n,one,hh,nevl,s,lh,zer,h,lh)
C     call zprm('zz+ * h * zz',2,h,lh,nevl,nevl)

C ... Diagonalize
      if (nmx <= 0) then
        jobz = 'N'
        lwork = 4*nevl
        allocate(wk(lwork),rwork(max(1,3*nevl-2)))
        call zheev(jobz,'U',nevl,h,lh,e,wk,lwork,rwork,ier)
        nev = 0
        if (lfullh) then
          do  i = nevl+1, n
            e(i) = ecap
          enddo
        endif
      else
        jobz = 'V'
        nev = min(nmx,nevl)
        lwork = max(nevl*nev,2*nevl)
        allocate(wk(lwork),rwork(max(1,3*nevl-2)))
C       h contains eigenvectors in reduced hilbert space
        call zheev(jobz,'U',nevl,h,lh,e,wk,lwork,rwork,ier)
C       Artifically enlarge hilbert space of h
        if (lfullh) then
          do  i = nevl+1, n
            e(i) = ecap
            h(i,:) = 0d0
            h(:,i) = 0d0
            h(i,i) = 1d0
          enddo
          nev = n
        endif
C       Rotate back to nonorthogonal evecs including overlap matrix
        if (lfullh .and. nevl < n) then
          call zgemm('N','N',n,n,n,one,s,lh,h,lh,zer,z,lz)
        else
          call zgemm('N','N',n,nev,nevl,one,s,lh,h,lh,zer,z,lz)
        endif
C       call prmx('evl',e,n,n,1)
C       call zprm('z',2,z,lz,n,nev)
      endif
      deallocate(wk,rwork,hh)

      call tcx('zhevo')

      end
!
! Test code for some of the diagonalisations above. The whole thing should be replaced with unievp pending removal of the memory doubling there.
!           real(8) :: tol
!           integer, allocatable :: iwork(:), ifail(:)
!           real(8), allocatable :: rwork(:)
!
!           procedure(real(8)) :: dlamch
!           tol = 2*dlamch('S')
!           ier = 0
!           allocate(ifail(n))
!           allocate(iwork(5*n))
!           allocate(rwork(7*n))
!
!           deallocate(work)
!           allocate(work(1))
!           lwork = -1
!
!           call zhegvx(1,jobz,'a','u',n,h,lh,s,lh,0.0_8,0.0_8,0,0,tol,nev,e,z,lz,work,lwork,wk,iwork,ifail,ier)
!
!           lwork  = max(4,int(work(1)))
!
!           deallocate(work)
!           allocate(work(lwork))
!
!           call zhegvx(1,jobz,'a','u',n,h,lh,s,lh,0.0_8,0.0_8,0,0,tol,nev,e,z,lz,work,lwork,wk,iwork,ifail,ier)
!
!           deallocate(work)
!           deallocate(rwork)
!           deallocate(iwork)
!           deallocate(ifail)
