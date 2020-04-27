      subroutine hoffs(nbas,lphi,lmxa,n0,ips,nel,ioffb,ioffa,ioffl,nlbx)
C- Make offset tables ioffb for basis, ioffl for linked, ioffa for aug.
C  Note: lphi(nel+1) must be set to max(lphi) for linked case.
      implicit none
      integer n0,nbas,nel,lphi(n0,2),ioffb(nbas,nel),ioffa(nbas+1),
     .  ioffl(nbas+1),ips(1),ib,ie,is,nlmb,ibe,lmxa(nbas),nlb,nlbx,ioff

C --- Offsets table for hamiltonian ---
      ioffb(1,1) = 0
      ioff = 0
      ibe = 1
      do  10  ie = 1, nel
      do  10  ib = 1, nbas
      is = ips(ib)
      nlmb = (lphi(ie,is)+1)**2
      ibe = ibe+1
      ioffb(ibe,1) = ioff + nlmb
      ioff = ioffb(ibe,1)
   10 continue
C --- Offsets table for augmentation ---
      ioffa(1) = 0
      do  12  ib = 1, nbas
      is = ips(ib)
   12 ioffa(ib+1) = ioffa(ib) + (lmxa(is)+1)**2
      nlbx = 0
      do  14  ie = 1, nel
      nlb = ioffb(nbas+1,ie)-ioffb(1,ie)
   14 nlbx = max0(nlbx,nlb)
C --- Offsets table for linking energy ---
      ie = nel+1
      ioffl(1) = 0
      do  20  ib = 1, nbas
      is = ips(ib)
      nlmb = (lphi(ie,is)+1)**2
      ioffl(ib+1) = ioffl(ib) + nlmb
   20 continue

      end
