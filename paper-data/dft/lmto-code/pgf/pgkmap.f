      subroutine pgkmap(s_ham,s_pot,s_str,s_site,plat,ipl,pgplp,qp,isp,mode,zp,rvec,r)
C- Generate evals and evecs of bulk GF, and optionally bulk GF
C ----------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 ldham lncol lham neula offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  pgbevl plham
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:palp pfnc
Cio    Passed to:  pgbevl plham
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  pgbevl plham
Ci Inputs:
Ci   plat  :primitive lattice vectors, in units of alat
Ci   ipl   :index to current principal layer
Ci          -1 => strux for left bulk PL; npl => strux for right bulk
Ci   npl   :number of principal layers (pgfset.f)
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   qp    :k-point
Ci   isp   :current spin channel (1 or 2)
Ci   mode  : 1s digit
Ci             0, Use strx (s00,s0L,s0L+)
Ci             1, Use strx (s00,s0R+,s0R)
Ci           10s digit:
Ci             0, r,rvec are passed externally
Ci             1, internal arrays used for r,rvec
Ci   zp    :complex energies at which to evaluate r
Co Outputs:
Co   r,rvec:eigenvalues and eigenvectors of bulk GF, layer ipl
Co          provided 10s digit of mode is set.  Else, left untouched.
Co   sbulk :s00,s0L,s0R that generate evals and evecs
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,ipl,mode,pgplp(6,-1:*)
      double precision plat(3,3),qp(3),rvec(*),r(*),zp(2,1)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ir,ld0,l2,ifi,npl,nspc
      double precision tol,pi,xv(1)
      double complex logk
      procedure(integer) :: fopn,nglob

      npl = nglob('npl')
      ld0 = pgplp(4,ipl)
      l2  = 2*ld0
      pi = 4*datan(1d0)
      ifi = fopn('BNDS')
      tol = 1d-3
      nspc = 1
      call pgbevl(s_ham,s_pot,s_str,s_site,plat,isp,nspc,ipl,npl,pgplp,qp,mode,-1,rvec,r,xv,xv,ld0,0,xv,xv)

      do  ir = 1, 2*ld0
        if (abs(abs(dcmplx(r(ir),r(ir+l2)))-1) < tol) then
          logk = log(dcmplx(r(ir),r(ir+l2)))
          write(ifi,333) dimag(logk)/(2*pi), zp, isp
  333     format(3f12.6,i3)
          print 333, dimag(logk)/(2*pi), zp, isp
        endif
      enddo

      end
