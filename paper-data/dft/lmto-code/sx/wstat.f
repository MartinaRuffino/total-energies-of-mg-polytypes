      subroutine wstat(nn,p0,v0,ws)
C- Screened potential ws := v0*(1- p0*v0)^-1
C  Output: ws
Cr T. Sandu Oct. 14 2004
Cr.Adapted to spin polarized case:
Cr ws := v0*(1- p01*v0-p02*v0)^-1, where
Cr p01 and p02 are the polarization functions for spin up and down.
C  Debugging
C  mc -f12f9.4 v0 -a v0 p0 -a p0 -1:12 p0 v0 -x -- -i  v0 -tog -x
C  set sx = "-subx 2,3,2,3 -subx 3,4,3,4 -subx 4,5,4,5 -subx 5,6,5,6"  n = 4
C  set n = 4
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nn
      double precision p0(nn,nn,2,2),v0(nn,nn,2),ws(nn,nn,2)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wkr(:)
      real(8), allocatable :: wki(:)
      real(8), allocatable :: wk(:)
C ... Local parameters
      integer ipr,ierr,nn2

      call getpr(ipr)
      nn2 = nn**2*2

      if (ipr >= 100) then
C        call o1r(v0,nn,' v0r   ')
C        call o1r(v0(1,1,2),nn,' v0i   ')
C        call o1r(p0,nn,' p0r   ')
C        call o1r(p0(1,1,2),nn,' p0i   ')
        call yprm('p01',2,p0(1,1,1,1),nn*nn,nn,nn,nn)
        call yprm('p02',2,p0(1,1,1,2),nn*nn,nn,nn,nn)
        call yprm('v0',2,v0,nn*nn,nn,nn,nn)
      endif

C --- Make wkr,wki := -(1 - p01 v0 -p02 v0) ---
C    ...spin non-polarized case p02 is set to 0 by default
      allocate(wkr(nn*nn))
      allocate(wki(nn*nn))
      call yygemm('N','N',nn,nn,nn,1d0,p0(1,1,1,1),p0(1,1,2,1),nn,
     .  v0,v0(1,1,2),nn,0d0,wkr,wki,nn)
      call yygemm('N','N',nn,nn,nn,1d0,p0(1,1,1,2),p0(1,1,2,2),nn,
     .  v0,v0(1,1,2),nn,1d0,wkr,wki,nn)
      call daxpy(nn,1d0,-1d0,0,wkr,nn+1)

C --- Make ws := v0 (1 - p0 v0)^-1 ---
      allocate(wk(nn*(nn+1)))
      call dpcopy(v0,ws,1,nn2,-1d0)
      call yyqnvb('t',wkr,wki,nn,nn,nn,wk,nn,wk,
     .  ws,ws(1,1,2),nn,ierr)
      if (ierr /= 0) call rx('wstat failed to invert v0')

      if (ipr >= 100) then
        call yprm('wstat',2,ws,nn*nn,nn,nn,nn)
      endif

      deallocate(wkr,wki,wk)
      end
