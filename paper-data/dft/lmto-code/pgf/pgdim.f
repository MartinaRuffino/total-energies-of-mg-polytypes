      integer function pgdim(mode,ipl,npl,npl1,npl2,pgplp)
C- Generate dimensions for layer Green's function connecting neighbor PL
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return dimension of g connecting layers npl1..npl2
Ci         :1 return (dimension of ipl) * (dimension of g)
Ci         :2 suppress ipl in counting (for offsets)
Ci         :3 return largest PL dimension between npl1...npl2
Ci   ipl   :index to current principal layer; used for mode=1,2
Ci   npl   :total number of PL, needed for indexing pgplp at endpoints
Ci   npl1,npl2:accumulate dimensions for PL npl1 .. npl2; see mode
Co Outputs
Co   pgdim :dimension of off-diagonal GF, according to mode
C ----------------------------------------------------------------------
      implicit none
      integer mode,ipl,npl,npl1,npl2,pgplp(6,-1:npl)
      integer jpl,ndim

      ndim = 0
      do  10  jpl = npl1, npl2
C   ... Case do not add ipl to offset (100s digit mode = 2)
        if (mode == 2 .and. jpl == ipl) goto 10
        if (mode == 3) then
          ndim = max(ndim,pgplp(4,max(min(jpl,npl),-1)))
        else
          ndim = ndim + pgplp(4,max(min(jpl,npl),-1))
        endif
   10 continue
      if (mode == 1) ndim = ndim*pgplp(4,ipl)
c     print *, 'pgdim: ipl,npl1,npl2,ndim=',ipl,npl1,npl2,ndim
      pgdim = ndim

      end
