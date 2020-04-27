C#define SP
      subroutine remhyb(ldim,nl,nbas,h,o)
C- Remove hybridisation blocks of H and O
C-----------------------------------------------------------------------
Ci ldim,nl,nbas,h,o
Co h and o changed on exit
Cr precompiler option SP retains sp hybridisation
C-----------------------------------------------------------------------
      implicit none
      integer ldim,nl,nbas
      double precision h(0:ldim-1,0:ldim-1,2),o(0:ldim-1,0:ldim-1,2)

      integer nl2,ip,jp,nlm1,nlm2

      if (nl == 1) return

      nl2 = nl**2
      do  nlm1 = 0, nl2-1
        if (nlm1 == 0) then
          do  nlm2 = 0, nl2-1
C#ifdef SP
            if (nlm2 > 3) then
C#elseC
C            if (nlm2 > 0) then
C#endif
              do  ip = 0, (nbas-1)*nl2, nl2
              do  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
              enddo
              enddo
            endif
          enddo
        elseif (nlm1 < 4) then
          do  nlm2 = 0, nl2-1
C#ifdef SP
            if (nlm2 > 3) then
C#elseC
C            if (nlm2 == 0 .or. nlm2 > 3) then
C#endif
              do  ip = 0, (nbas-1)*nl2, nl2
              do  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
              enddo
              enddo
            endif
          enddo
        elseif (nlm1 < 9) then
          do  nlm2 = 0, nl2-1
            if (nlm2 < 4 .or. nlm2 > 8) then
              do  ip = 0, (nbas-1)*nl2, nl2
              do  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
              enddo
              enddo
            endif
          enddo
        elseif (nlm1 < 16) then
          do  nlm2 = 0, nl2-1
            if (nlm2 < 9) then
              do  ip = 0, (nbas-1)*nl2, nl2
              do  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
              enddo
              enddo
            endif
          enddo
        endif
      enddo
      end


