C#define OKA
      integer function ioka()
      implicit none
C- returns 0, or 1000 if OKA is defined
C#ifdef OKA
      ioka = 1000
C#elseC
C      ioka = 0
C#endif
      end
