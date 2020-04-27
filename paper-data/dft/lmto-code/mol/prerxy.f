      subroutine prerxy(bin,bout,nlm,ixy)
c  pre-multiplies bin by differential rotation matrix
c  iab=1: x-rotation. iab=2: y-rotation
      implicit real*8 (a-h,p-z), integer(o)
      dimension bin(1),bout(1)
      lmax=ll(nlm)
      call dpzero(bout,nlm)
c -------- here for ixy=1 -----------
      if(ixy == 1) then
      do 20 l=0,lmax
        mm=l*l+l+1
        k=0
        do 20 m=-l,l-1
        n=m+1
        k=k-m
        aa=dsqrt(0.5d0*k)
        if(n == 0) aa=0d0
        if(n == 1) aa=dsqrt(2d0)*aa
        if(m < 0) aa=-aa
        bout(n+mm)=bout(n+mm)+aa*bin(m+mm)
  20    bout(m+mm)=bout(m+mm)-aa*bin(n+mm)
      endif
c -------- here for ixy=2 -----------
      if(ixy == 2) then
      do 30 l=0,lmax
        mm=l*l+l+1
        k=0
        do 30 m=-l,l-1
        n=m+1
        k=k-m
        aa=dsqrt(0.5d0*k)
        if(n == 0) aa=0d0
        if(n == 1) aa=dsqrt(2d0)*aa
        if(m < 0) aa=-aa
        modm=mod(m+40,4)
        m1=m
        if(modm == 1.or.modm == 3) m1=-m
        sgm=1d0
        if(modm == 2.or.modm == 3) sgm=-1d0
        modn=mod(n+40,4)
        n1=n
        if(modn == 1.or.modn == 3) n1=-n
        sgn=1d0
        if(modn == 2.or.modn == 3) sgn=-1d0
        bout(n1+mm)=bout(n1+mm)+aa*sgm*sgn*bin(m1+mm)
  30    bout(m1+mm)=bout(m1+mm)-aa*sgm*sgn*bin(n1+mm)
      endif
      end
