      subroutine gtqval(ib1,ib2,znuc,qc,ipc,zval)
C- Get the nuclear-core charge for a block of sites
      implicit none
      integer ib1,ib2,ipc(ib2),ib,ic
      double precision znuc(*),qc(*),zval
C     if (zval /= 0) return
      zval = 0
      do  12  ib = ib1, ib2
      ic = ipc(ib)
   12 zval = zval + znuc(ic)-qc(ic)
      end
