      subroutine hoffcl(ic1,ic2,iaxg,offsH,offc1,offcH,ndimH)
C- For one cluster, table of Hamiltonian offsets and ham. dimension
C ----------------------------------------------------------------------
Ci Inputs
Ci  ic1,ic2:first, final members in iaxg for which to calculate offcH
Ci   iaxg  :neighbor table for a cluster (paircl.f)
Ci   offsH :put table in offcH(offsH+1,...)
Ci          usually 0 unless offcH is part of a larger array.
Ci   offc1 :starting value, i.e. offcH(offsH+1)
Ci          usually 0 unless offcH is part of a larger array.
Co Outputs
Co   offcH :Table of hamiltonian offsets for this cluster
Co          starting at offcH(offsH+1)
Co   ndimH :dimension for this cluster
Cr Remarks
Cr   Uses iaxg(9)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax
      parameter (niax=10)
      integer ic1,ic2,iaxg(niax,1),offsH,offc1,offcH(*),ndimH
C Local variables:
      integer ic,kc

      offcH(offsH+1) = 0
      kc = offsH
      do  10  ic = ic1, ic2
        kc = kc+1
        offcH(kc+1) = offcH(kc) + iaxg(9,ic)
   10 continue

C     call yprm('offch',0,offch,0,ntab0+1,ntab0+1,2)

      end
