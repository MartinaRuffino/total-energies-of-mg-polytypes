      subroutine mkiaxd(npr,lmx,ips,iax)
C- Sets iax(9,:) to number of orbitals
C ----------------------------------------------------------------------
Ci Inputs
Ci   npr   :set iax(9,1:npr)
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   ips   :species table: site ib belongs to species ips(ib)
Cio Inputs/Outputs
Cio  iax   :neighbor table containing pair information (pairc.f)
Cio        :iax(9,:) set to (lmx+1)**2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer npr,niax,lmx(*),ips(*)
C ... Local parameters
      parameter (niax=10)
      integer iax(niax,npr)
      integer nli,is,i,jbas
c     ibas0 = 0
      do  i = 1, npr
c       ibas = iax(1,i)
        jbas = iax(2,i)
        is = ips(jbas)
        nli = (lmx(is)+1)**2
        iax(9,i) = nli
C        if (ibas /= ibas0) then
C          ibas0 = ibas
CC         iax(10,i) = nli
C        else
CC         iax(10,i) = iax(10,i-1) + nli
C        endif
C       print *, isum(50,iax(9,1),niax)
C      print *, isum(npr,iax(9,1),niax)
      enddo
      end
