      subroutine gfrscp(mode,clp,icl1,icl2,gsrc,gdst)
C- Copy one R.S. Green's function to another, possibly changing format
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit: source matrix
Ci          0 matrix has real followed by imaginary
Ci          1 matrix is in complex*16 format
Ci          2 matrix is in transpose form (see Remarks)
Ci          3 combination of 1 and 2
Ci         10s digit: destination matrix
Ci          0 matrix has real followed by imaginary
Ci          1 matrix is in complex*16 format
Ci          2 matrix is in transpose form (see Remarks)
Ci          3 combination of 1 and 2
Ci   clp   :index and dimensioning information for cluster (suordn.f)
Ci   icl1,icl2:range of clusters to copy
Ci   gsrc  :source matrix
Co Outputs
Co   gdst  :destination matrix
Cr Remarks
Cr   g is a collection of ncl two-dimensional arrays g_mk all strung
Cr   together.  It would be simpler to have a vector of ncl pointers,
Cr   but fortran does not allow them, so they are strung together
Cr   instead.  The offset to the i'th array is offg = clp(6,i).
Cr   Relative to that offset the i'th array g is dimensioned
Cr   g(ndimg,ndimb); see the code to see how they are extracted.
Cr
Cr   Alternatively, g is a a series of arrays {g_km}, in which the
Cr   offsets to the start of array k is as above, but the array g_km is
Cr   stored as the transpose of g_mk.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nclp,icl1,icl2,mode
      parameter (nclp=9)
      integer clp(nclp,icl2)
      double precision gsrc(*),gdst(*)
C Local variables:
      integer nsc,icl,offxg,offg,ndimb,ndimg
      integer nrs,nrd,ncs,ncd,offs,offd,modes,moded

      modes = mod(mode,10)
      moded = mod(mode/10,10)

C --- For each cluster, copy gsrc to gdst ---
      do  10  icl = icl1, icl2
        nsc   = clp(1,icl+1) - clp(1,icl)
        offxg = clp(2,icl)
        offg  = clp(6,icl)
        ndimb = clp(3,icl)
        ndimg = (clp(6,icl+1) - clp(6,icl)) / ndimb
        print 333, nsc, offxg,offg,ndimg,ndimb
  333   format(6i7)

        nrs  = 1
        ncs  = ndimg
        offs = ndimg*ndimb
        if (modes >= 2) then
          nrs = ndimb
          ncs = 1
        endif
        if (mod(modes,2) == 1) then
          nrs = 2*nrs
          ncs = 2*ncs
          offs = 1
        endif

        nrd  = 1
        ncd  = ndimg
        offd = ndimg*ndimb
        if (moded >= 2) then
          nrd = ndimb
          ncd = 1
        endif
        if (mod(moded,2) == 1) then
          nrd = 2*nrd
          ncd = 2*ncd
          offd = 1
        endif

        call ymcpy(gsrc(1+2*offg),ncs,nrs,offs,
     .             gdst(1+2*offg),ncd,nrd,offd,ndimg,ndimb)

C#ifdefC DEBUG
C        if (moded == 0)
C     .  call yprm('gdst',2,gdst(1+2*offg),ndimg*ndimb,ndimg,ndimg,ndimb)
C        if (moded == 1)
C     .  call yprm('gdst',3,gdst(1+2*offg),ndimg*ndimb,ndimg,ndimg,ndimb)
C        if (moded == 2)
C     .  call yprm('gdst',2,gdst(1+2*offg),ndimb*ndimg,ndimb,ndimb,ndimg)
C        if (moded == 3)
C     .  call yprm('gdst',3,gdst(1+2*offg),ndimb*ndimg,ndimb,ndimb,ndimg)
C#endif

   10 continue

      end
