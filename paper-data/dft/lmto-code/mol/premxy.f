      subroutine premxy(bin,bout,nlm,nb,ixy)
c  pre-multiplies bin by differential rotation matrix, add into bout
c  iab=1: x-rotation. iab=2: y-rotation
      implicit real*8 (a-h,p-z), integer(o)
      dimension bin(nlm,nb),bout(nlm,nb)
      lmax=ll(nlm)
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
        do 20 i=1,nb
        bout(n+mm,i)=bout(n+mm,i)+aa*bin(m+mm,i)
  20    bout(m+mm,i)=bout(m+mm,i)-aa*bin(n+mm,i)
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
        do 30 i=1,nb
        bout(n1+mm,i)=bout(n1+mm,i)+aa*sgm*sgn*bin(m1+mm,i)
  30    bout(m1+mm,i)=bout(m1+mm,i)-aa*sgm*sgn*bin(n1+mm,i)
      endif
      end
      subroutine pretxy(bin,bout,nb,nlm,ixy)
c  multiplies bin by tranpose of differential rot matrix, add into bout
c  iab=1: x-rotation. iab=2: y-rotation
      implicit real*8 (a-h,p-z), integer(o)
      dimension bin(nb,nlm),bout(nb,nlm)
      lmax=ll(nlm)
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
        do 20 i=1,nb
        bout(i,n+mm)=bout(i,n+mm)+aa*bin(i,m+mm)
  20    bout(i,m+mm)=bout(i,m+mm)-aa*bin(i,n+mm)
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
        do 30 i=1,nb
        bout(i,n1+mm)=bout(i,n1+mm)+aa*sgm*sgn*bin(i,m1+mm)
  30    bout(i,m1+mm)=bout(i,m1+mm)-aa*sgm*sgn*bin(i,n1+mm)
      endif
      end
