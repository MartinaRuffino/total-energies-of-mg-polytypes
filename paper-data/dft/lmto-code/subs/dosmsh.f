      subroutine dosmsh(mode,rnge,de,nptmx,npts,emesh)
C- Return a mesh of points, depending on mode
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :Determines kind of mesh
Ci         :sign of mode: if <0, only npts is returned, not mesh itself
Ci         :1  uniform mesh of points spaced de, in rnge(1..2)
Ci         :   emesh(i) = rnge(1) + de*(i-1)
Ci         :11 quadratic mesh:
Ci         :   emesh(i) = rnge(1) + de*(i-1) + de**2/omg/2*(i-1)**2
Ci         :2-9 same as mode 1, but multiple windows
Ci         :rnge(1..2); rnge(2..3); ...; rnge(mode:mode+1)
Ci         :100s digit:
Ci         :1 Add points until npts == nptmx,
Ci         :  instead of using rnge(2) as upper bound.
Ci         :  Not compatible with 1s digit mode 2-9
Ci   rnge  :energy boundaries:
Ci         :rnge(i)   = lower bound to ith window
Ci         :rnge(i+1) = upper bound to ith window
Ci   de    :mesh spacing.  First argument de, 2nd argument omg
Ci         : Uniform mesh:
Ci         : emesh(i) = rnge(1) + de*(i-1)
Ci         : Faleev's quadratic mesh:
Ci         : emesh(i) = rnge(1) + de*(i-1) + de**2/omg/2*(i-1)**2
Ci         : Reduces to linear spacing as omg->infty
Ci         : The point at which the spacing approximately doubles
Ci         : that of de is i = omg/de.
Ci   nptmx :leading dimension of emesh (not used if mode<0)
Co Outputs
Co    npts :number of mesh points
Co   emesh :energy mesh (not touched if mode<0)
Co         :It is the caller's
Cl Local variables
Cl         :
Cr Remarks
Cr   Note that Takao's freq_r in hx0fp0 is in hartree, and consists
Cr   of midpoints of frhis, defined here as emesh:
Cr     freq_r(iw) = (frhis(iw)+frhis(iw+1))/2d0
Cu Updates
Cu   03 Apr 17 Slight redesign of 100s digit for consistency
Cu    6 Dec 13 swapped ie += de for ie = de*(i-1). (dmt)
Cu             This routine seems to often abused with aliased nptmx and npts
Cu   14 Sep 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nptmx,npts
      double precision de(*),rnge(*),emesh(nptmx)
C ... Local parameters
      logical lmakee,linc
      integer modea,j
      double precision ei,fuzz
      parameter (fuzz = 1d-10)

      call rx('DOSMSH not installed; optics library required')

      end
