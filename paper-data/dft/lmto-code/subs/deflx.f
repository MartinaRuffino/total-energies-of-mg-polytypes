      subroutine deflx(mode,modeb,z,lmxb,lmxa)
C- Generate default quantities for species
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 generate default lmxb
Ci         :2 generate default lmxa
Ci         : Combination 1+2 allowed
Ci   modeb :mode for autoset lmxb (see fadflb)
Ci   z     :nuclear charge
Co Outputs
Co   lmxb  :basis l-cutoff
Co   lmxa  :augmentation l-cutoff
Cr Remarks
Cr
Cu Updates
Cu   19 Oct 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,modeb,lmxb(*),lmxa
      double precision z
C ... Local parameters
      logical d3,d4,d5,rareE,actinide
C     logical nm
      integer iqocc(5)
C      double precision

      d3 = z >= 21 .and. z <= 29
      d4 =  z >= 39 .and. z <= 47
      d5 =  z == 57 .or. z >= 72 .and. z <= 79
C     nm  = z == 29 .or. z == 47 .or. z == 79
      rareE = z >= 58 .and. z <= 71
      actinide = z >= 90 .and. z <= 103

      if (mod(mode,2) /= 0) then
        call fadflb(modeb,99,z,lmxb,iqocc)
      endif

      if (mod(mode/2,2) /= 0) then
        lmxa = 3
        if (z >= 55 .or. d3 .or. d4 .or. d5) lmxa = 4
        if (rareE .or. actinide) lmxa = 6
        if (mod(mode,2) /= 0) lmxa = max(lmxa,lmxb(1))
      endif

      end
