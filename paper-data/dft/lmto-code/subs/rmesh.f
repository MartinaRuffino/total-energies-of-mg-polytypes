      subroutine rmesh(z,rmax,lrel,lgrad,nrmx,a,nr)
C- Generate parameters for shifted logarithmic radial mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   lrel  :0 for non-relativistic
Ci         :1 for scalar relativistic
Ci         :2 for Dirac equation
Ci   lgrad :0 for LDA, nonzero for gradient corrections
Ci   nrmx  :maximum allowed number of points
Cio Inputs/Outputs
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio        :a is not altered if input a>0; otherwise a is set here.
Cio        :When a is set, it is independent of rmax and nr
Cio  nr    :number of radial mesh points
Cio        :nr is not altered if input nr>0; otherwise nr is set here.
Cio        :The calculated value of nr depends on both a and z
Cr Remarks
Cr   Uses input values for a,nr if >0; otherwise rmesh sets them
Cu Updates
Cu   01 Aug 16 Relativistic case no longer uses NR=1501; chooses a=.015 as default
Cu   19 Apr 12 Remove restriction nr must be odd
Cu   18 Mar 03 Default parameters for fully relativistic case
Cu   11 Oct 02 No longer uses a smaller number of points for
Cu             the nonrelativistic case.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nrmx,nr,lgrad,lrel
      double precision z,rmax,a
C ... Local parameters
      double precision b
      integer nr0
      integer,parameter:: NULLI =-99999, nrmax=5001, nrrel=1501

C     Hang on to given value
      nr0 = nr

      b = 1d0/(2*z+1)
      if (lrel == 2 .and. a < 1d-6) a = .015d0
      if (lgrad /= 0) then
        if (a < 1d-6) a = 0.015d0
        if (nr <= 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
      else
        if (a < 1d-6) a = 0.03d0
        if (nr <= 0) nr = 2*(.5d0+dlog(1+rmax/b)/a)
      endif
      nr = max0(51,((nr-1)/2)*2+1)
      if (iabs(nr-nr0) == 1) nr = nr0
      if (nrmx > 0) nr = min0(nr,nrmx)

      call info5(50,0,0,' RMESH:  Z=%d  a=%1;4d  nr=%i  rmax=%1;6d',
     .  z,a,nr,rmax,0)

      end
