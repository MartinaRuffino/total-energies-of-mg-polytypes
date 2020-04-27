      subroutine pggmal(mode,ip1,ip2,pgplp,kp1,gii)
C- Allocates memory for layer GF gii
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          1 gf is noncollinear; dimensions of g double
Ci   ip1,ip2 range of PL for which to allocate memory
Ci   kp1   :not used, but intention was to put array pointers in kp1,kp1+1,...
Ci         :Thus gii(kp1) allocated for gf(ip1),
Ci         :     gii(kp1+1) allocated for gf(ip1+1),  etc
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Co Outputs
Co    gii  : table of arrays containing GF by layer
Cr Remarks
Cu Updates
Cu   02 Sep 15 complete migration to f90 pointers
C ----------------------------------------------------------------------
      use structures, only : s_lgf
      implicit none
C ... Passed parameters
      integer mode,ip1,ip2,pgplp(6,-1:ip2),kp1
C     Pointer arrays to different layer GF's
      type (s_lgf) :: gii(-2:*)
C ... Local parameters
      integer ipl,ld0,ndg,nspc

      nspc = 1
      if (mod(mode,10) == 1) nspc = 2
      do  ipl = ip1, ip2
        ld0 = pgplp(4,ipl)*nspc
        ndg = pgplp(3,ipl)*nspc
        allocate(gii(ipl+kp1-ip1)%gll(ld0,ndg))
      enddo
      end

      subroutine pggfre(mode,ip1,ip2,kp1,gii)
C- Free Dynamic memory allocation for layer GF
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode
Ci   ip1,ip2 range of PL for which to free memory
Ci   pgplp :index and dimensioning information for each PL (pgfset.f,subas.f)
Ci   kp1   :Array pointers stored in kp1,kp1+1,...
Ci         :Thus gii(kp1) allocated for gf(ip1),
Ci         :     gii(kp1+1) allocated for gf(ip1+1),  etc
Co Outputs
Co    gii  : table of arrays containing GF by layer
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      use structures, only : s_lgf
      implicit none
C ... Passed parameters
      integer mode,ip1,ip2,kp1
C     Pointer arrays to different layer GF's
      type (s_lgf) :: gii(-2:*)
C ... Local parameters
      integer ipl

      do  ipl = ip2, ip1, -1
        deallocate(gii(ipl+kp1-ip1)%gll)
      enddo
      end

