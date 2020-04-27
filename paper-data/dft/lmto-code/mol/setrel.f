C#define LINUX_PGI | INTEL_IFC | INTEL_IFORT
      integer function lrel()
      implicit none
C#ifndefC LINUX_PGI | INTEL_IFC | INTEL_IFORT
C      integer lsave,lr,setrel
C      save
C#else
      integer lsave,lr
      common /setrel/ lsave
C#endif
      lrel=lsave
      return
C#ifndefC LINUX_PGI | INTEL_IFC | INTEL_IFORT
C      entry setrel(lr)
C      lsave=lr
C      return
C#endif
      end

      subroutine relpr(ifi)
      if(ifi >= 1) then
         if (lrel() == 0) write(ifi,101)
         if (lrel() /= 0) write(ifi,102)
  101    format(' non-relativistic calculation')
  102    format(' scalar-relativistic calculation')
      endif
      end
