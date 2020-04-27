      integer function mxbnds(mode,eb,nde,nband,nk,ef,kmax)
C- Find the maximum number of states below Ef among a set of k-points
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0, find the largest  number of bands under Ef for any k
Ci             There should be NO bands eb(nmax+1,1:nk)<ef
Ci             There should be at LEAST one state eb(nmax,1:nk)<ef
Ci         :1, find the smallest number of bands under Ef for any k
Ci             There should be NO bands eb(nmax,1:nk)>ef
Ci             There should be at LEAST one state eb(nmax+1,1:nk)>ef
Ci         :2, Same as mode 0, but add 1 for padding
Ci   eb    :eigenvalues k-points 1..nk
Ci   nde   :leading dimension of eb
Ci   nband :maximum number of states to search
Ci   nk    :number of k-points
Ci   ef    :Fermi level
Co Outputs
Co   mxbnds: highest number of occupied states
Co   kmax  : the index in eb(:,k) which determined mxbnds
Cl Local variables
Cr Remarks
Cr   Adapted from noccx1 in tetwt4.F
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nde,nband,nk,kmax
      double precision eb(nde,nk),ef
C ... Local parameters
      integer nmax,k,it
C     integer kmax

      nmax = 0
      if (mode == 1) nmax = nband+1
      do  k = 1, nk
        do  it = 1, nband
          if (eb(it,k) > ef) goto 1111
        end do
        it = nband+1
 1111   continue
C   ... At this point it points to first band above ef
        if (mode == 0) then
          it = it-1
          if (it > nmax) then
            nmax = it
            kmax = k
          endif
        elseif (mode == 1) then
          it = it-1
          if (it < nmax) then
            nmax = it
            kmax = k
          endif
        elseif (mode == 2) then
          if (it > nmax) then
            nmax = min(it,nband)
            kmax = k
          endif
        endif
C       print *, it,nmax,kmax
      end do
      mxbnds = nmax

      end

