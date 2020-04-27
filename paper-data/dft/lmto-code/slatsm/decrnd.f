      double precision function decrnd(x,n)
C- Truncate after n significant figures
C ----------------------------------------------------------------
Ci Inputs
Ci   x,n
Co Outputs
Co   decrnd
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x
      integer n
C Local parameters
      double precision scale
C This construction avoids taking integer part of negative number
      if (x == 0) then
        decrnd = 0
        return
      endif
      scale = dlog10(dabs(x))
      if (scale > 0) then
        scale = int(scale)
      else
        scale = -int(-scale) - 1
      endif
      scale = 10d0**(n-1-scale)
      decrnd = nint(x*scale)/scale
      end
C#ifdefC TEST
C      double precision decrnd,x
C      x = 12345678
CC     x = -x
C      do  20 j = 1,7
C        do  10  i = -20,20
C   10   print *, decrnd(x*(10d0**i),j), decrnd(x*(10d0**i),j)/
C     .      (x*(10d0**i))
C        print *, '----'
C   20 continue
C      end
C#endif
