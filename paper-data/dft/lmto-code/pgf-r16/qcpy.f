      subroutine qcpy(tr,ti,dr,di,t1,t2)
C- complex multiply (t1,t2) = (tr,ti) * (dr,di)
      real*16 tr,ti,dr,di,t1,t2
      real*16 tmp
      tmp = tr*dr - ti*di
      t2  = tr*di + ti*dr
      t1 = tmp
      end
