C#define LINUX_PGI | INTEL_IFC | INTEL_IFORT
      integer function lxcf()
C- returns index to xc functional
      implicit none
C#ifndefC LINUX_PGI | INTEL_IFC | INTEL_IFORT
C      integer lsave,lsavg,lxcg,lxcfun,lxcgrd,setxcf,setxcg
C      save lsave,lsavg
C      data lsave /1/ lsavg /0/
C#else
      integer lsave,lsavg,lxcg,lxcfun,lxcgrd
      common /setxc/ lsave,lsavg
C#endif
      lxcf=lsave
      return

      entry lxcg()
C- returns index to xc gradient functional
      lxcg=lsavg
      return

C#ifndefC LINUX_PGI | INTEL_IFC | INTEL_IFORT
C      entry setxcf(lxcfun)
C      lsave=lxcfun
C      return
C
C      entry setxcg(lxcgrd)
C      lsavg=lxcgrd
C#endif
      end
      subroutine prxcf(s1,s2,ifi)
      implicit none
      character*(*) s1,s2
      integer lxcf,lxcg,ifi,ip
      double precision uk
      s1 = ' '
      if (lxcf() == 1) s1 = 'Ceperly-Alder'
      if (lxcf() == 2) s1 = 'Barth-Hedin'
      if (lxcf() == 3) s1 = 'Perdew-Wang ''91'
      s2 = ' '
      if (lxcg() == 1) s2 = 'Langreth-Mehl'
      if (lxcg() == 2) s2 = 'Perdew-Wang ''91'
      if (lxcg() == 3) then
        call getuk(uk)
        call awrit1("Perdew-Burke-Ernzerhof kappa=%d",s2,80,0,uk)
      endif
      if (lxcg() == 4) then
        call getuk(uk)
        call awrit1("Perdew-Burke-Ernzerhof kappa=%d, Becke exchange",
     .              s2,80,0,uk)
      endif
      if (ifi >= 1) write(ifi,333) s1,s2
  333 format(' xc ',a,1x,' gc ',a)
      end
