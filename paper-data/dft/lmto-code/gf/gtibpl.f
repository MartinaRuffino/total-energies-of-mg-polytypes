      subroutine gtibpl(ibl,nbl,pgplp,ib1,ib2)
C- Find range of basis atoms (first and last) in a principal layer or subblock
C ----------------------------------------------------------------------
Ci Inputs
Ci   ibl   :index to current principal layer or ham. subblock
Ci          subblocks range from 0 ... nbl-1
Ci          In the layer context, -1 => left bulk PL; nbl => right bulk
Ci   nbl   :number of principal layers or subblocks (susite)
Ci   pgplp :index and lattice dimensioning information for crystal subblocks.
Ci          The meaning of pgplp depends on the context; see susite.f
Co Outputs
Co   ib1,ib2:range of sites in subblock ibl
Cr Remarks
Cr   For principal-layer context, gtibpl uses the following:
Cr     nbas  = pgplp(1,nbl-1)
Cr     npadl = pgplp(1,0)
Cr     npadr = pgplp(1,nbl-1) - pgplp(1,nbl-2)
Cu Updates
Cu  13 Oct 98  convention for pgplp extended to non-layer case
C ----------------------------------------------------------------------
      implicit none
      integer ibl,nbl,ib1,ib2,pgplp(6,-1:nbl),jbl

      if (pgplp(1,-1) < 0)  then
        if (ibl < 0 .or. ibl > nbl) call rxi('gtibpl: bad ibl',ibl)
C   ... Permit ibl=nbl so call to gtibl(nbl,nbl,...) returns ib2 = nbas
        jbl = min(ibl,nbl-1)
C   ... max reqd because pgplp(1,-1) can't be both a flag and an offset
        ib1 = max(pgplp(1,jbl-1),0)+1
        ib2 = pgplp(1,jbl)
        return
      endif

      if (ibl <= -1) then
        ib1 = pgplp(1,nbl-1) + 1
        ib2 = ib1-1 + pgplp(1,0)
      else if (ibl >= nbl) then
        ib1 = pgplp(1,nbl-1) + pgplp(1,0) + 1
        ib2 = ib1-1 + pgplp(1,nbl-1) - pgplp(1,nbl-2)
      else
        ib1 = pgplp(1,ibl-1)+1
        ib2 = pgplp(1,ibl)
      endif
      end
