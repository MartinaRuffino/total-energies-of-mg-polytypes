      subroutine gfasdc(ib1,ib2,vbar,vshft,nl,nsp,ipc,qnu,edcc)
C- Compute double-counting correction for total energy
      implicit none
      integer ib1,ib2,ipc(ib2),nl,nsp
      double precision vbar,vshft(ib2),qnu(3,nl,nsp,*),edcc
      integer ib,ic
      double precision qat,vat,dsum

      call rx('this routine never checked')
      edcc = 0
      do  10  ib = ib1, ib2
        ic = ipc(ib)
        qat = dsum(nl*nsp,qnu(1,1,1,ic),3)
        vat = vshft(ib) + vbar
        edcc = edcc + qat * vat
   10 continue

      end
