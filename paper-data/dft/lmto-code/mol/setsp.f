C#define LINUX_PGI | INTEL_IFC  | INTEL_IFORT
      integer function lsp()
      implicit none
C#ifndefC LINUX_PGI | INTEL_IFC | INTEL_IFORT
C      integer lsave,setsp,lspin
C      save
C#else
      integer lsave
      common /setsp/ lsave
C#endif
      lsp=lsave
      return
C#ifndefC LINUX_PGI | INTEL_IFC | INTEL_IFORT
C      entry setsp(lspin)
C      lsave=lspin
C#endif
      end
