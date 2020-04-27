      subroutine setcg(s_lat,lmxcg,lmxcy)
C- Allocate space for, and make Clebsch-Gordan coeffs
C ----------------------------------------------------------------------
Cio Structures
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  cg jcg indxcg cy
Cio    Elts passed:cy cg indxcg jcg
Cio    Passed to:  *
Ci Inputs
Ci   lmxcg :Compute cg,jcg,indxcg for 0<l<lmxcg
Ci   lmxcy :Compute cy for 0<l<lmxcy
Co Outputs
Cr Remarks
Cu Updates
Cu   23 Dec 17 calls scg0 to determine dimensions of lnjcg and lnxcg for any lmxcg
Cu   01 Sep 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lmxcg,lmxcy
C ... For structures
!      include 'structures.h'
      type(str_lat):: s_lat
C ... Local parameters
      integer lnjcg,lnxcg,nlm,ixx(1)
      double precision xx(1)

C ... Allocate and occupy the arrays
      nlm = (lmxcy+1)**2
      call scg0(0,lmxcg,xx,ixx,ixx,lnjcg,lnxcg)
      call ptr_lat(s_lat,8+1,'cg',lnjcg,0,0,xx)
      call ptr_lat(s_lat,8+1,'jcg',lnjcg,0,0,xx)
      call ptr_lat(s_lat,8+1,'indxcg',lnxcg,0,0,xx)
      call ptr_lat(s_lat,8+1,'cy',nlm,0,0,xx)

      call sylmnc(s_lat%cy,lmxcy)
      call scg(lmxcg,s_lat%cg,s_lat%indxcg,s_lat%jcg)

C      call ptr_lat(s_lat,4+1,'cg',lnjcg,0,0,w(ocg))
C      call ptr_lat(s_lat,4+1,'jcg',lnjcg,0,0,w(ojcg))
C      call ptr_lat(s_lat,4+1,'indxcg',lnjcg,0,0,w(oidxcg))
C      call ptr_lat(s_lat,4+1,'cy',nlm,0,0,w(ocy))

      end
