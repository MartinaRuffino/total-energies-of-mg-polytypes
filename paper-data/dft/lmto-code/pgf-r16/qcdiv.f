      subroutine qcdiv(tr,ti,dr,di,t1,t2)
C- complex divide (t1,t2) = (tr,ti) / (dr,di)
Cr Remarks
Cr   Adapted from eispack.
Cr   It is permissible for (t1,t2) to occupy the same address space
Cr   as either (tr,ti) or (dr,di)
      real*16 tr,ti,dr,di,t1,t2
      real*16 rr,d,tmp

      if (abs(di) > abs(dr)) then
        rr = dr / di
        d = dr * rr + di
        tmp = (tr * rr + ti) / d
        t2 = (ti * rr - tr) / d
        t1 = tmp
      else
        rr = di / dr
        d = dr + di * rr
        tmp = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        t1 = tmp
      endif
      end