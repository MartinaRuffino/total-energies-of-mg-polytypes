      subroutine wsloca(nsp,indxsh,lihdim,nbas,nl,ipc,qnu,nn,wscrl)
C- Make local part of ASA SX potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   nn    :dimensions wscrl
Cio Inputs/Outputs
Cio   wscrl:On input, site-diagonal part of ASA potential
Cio        :On output, wsloc replaced by charge-weighted average on
Cio        :each site
Cr Remarks
Cr
Cu Updates
Cu   07 Oct 04 (T.Sandu) Adapted to handle spin polarized case
Cu   12 Jun 01
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nn,ipc(1),nsp,indxsh(*),lihdim,nbas,nl
      double precision qnu(3,nl,nsp,*),wscrl(nn,nn,2)
C ... Local parameters
      integer i0,ibas,ic,isp,l,lmr,lmr0,i
      double precision avg,sumq,sumw


      lmr = 0
      isp = 1
      i = 0
      do  ibas = 1, nbas
        ic = ipc(ibas)
        sumw = 0
        sumq = 0
        i0 = i
        lmr0 = lmr
C   ... Make average wscrl for this site
        do  l = 0, nl-1
          if (indxsh(lmr+1) > lihdim) goto 2
          i = i+1
C         if ((qnu(1,l+1,1,ic)+qnu(1,l+1,nsp,ic))/nsp < 5) then
          do isp = 1, nsp
           sumw = sumw + qnu(1,l+1,isp,ic)*wscrl(i,i,1)
           sumq = sumq + qnu(1,l+1,isp,ic)
          enddo
C         endif
    2     lmr = lmr + 2*l+1
        enddo
        avg = sumw/sumq
C   ... Distribute average wscrl over this site
        i = i0
        lmr = lmr0
        do  l = 0, nl-1
          if (indxsh(lmr+1) > lihdim) goto 3
          i = i+1
          wscrl(i,i,1) = avg
    3     lmr = lmr + 2*l+1
        enddo
      enddo

      end
