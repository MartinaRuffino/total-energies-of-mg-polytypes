      subroutine hpair3(nbas,ips,ntab,iax,rham,rtab,
     .  ntab3,iax3,rtab3)
C- Neighbor table for 3c terms in Hamiltonian matrix
C  Inputs:  ntab,iax,rham,rtab (from hpair)
C  Outputs: ntab3,iax3,rtab3, analogs of ntab,iax,rtab for 3C
      implicit none
      integer nbas,ips(nbas),ntab(nbas+1),iax(10,1),
     .  iax3(40),ntab3(nbas+1)
      double precision rham(nbas),rtab(3,1),rtab3(3,1)
      integer ib,is,jb,js,ipr,nt,i,nttab
      double precision r1,r2,rr
      real w(1)
      common /w/ w

C --- Printout for this ib ---
C     call awrit2('xx ntab %n:1i',' ',80,6,nbas+1,ntab)

      call getpr(ipr)
      if (ipr > 40) print 332
  332 format(/' hpair3     ib   jb',9x,
     .  '------- dr --------',11x,'d      r1+r2')

      nttab = 0
      ntab3(1) = 0
      do  10  ib = 1, nbas
        do  12  nt = ntab(ib)+1, ntab(ib+1)
        rr = dsqrt(rtab(1,nt)**2 + rtab(2,nt)**2 + rtab(3,nt)**2)
        r1 = rham(ips(iax(1,nt)))
        r2 = rham(ips(iax(2,nt)))
        if (rr < (r1+r2)/2) then
C snot for testing...
C         if (ib == 2 .and. nt == 5) then
C         else
          nttab = nttab+1
          iax3(nttab) = nt - ntab(ib)
C         iax3(1,nttab) = iax(1,nt)
C         iax3(2,nttab) = iax(2,nt)
C         iax3(3,nttab) = iax(3,nt)
C         iax3(4,nttab) = iax(4,nt)
C         iax3(5,nttab) = iax(5,nt)
C         rtab3(1,nttab) = rtab(1,nt)
C         rtab3(2,nttab) = rtab(2,nt)
C         rtab3(3,nttab) = rtab(3,nt)
          if (ipr > 40)
     .    print 333, iax(1,nt),iax(2,nt),(rtab(i,nt),i=1,3), rr,r1+r2
  333     format(i14,i5,3f11.6,2f9.4,2x,3i3)
        endif
   12 continue
      ntab3(ib+1) = nttab
   10 continue
      if (ipr >= 20 .and. ipr <= 40) print
     . '('' hpair3:'',i5,'' pairs taken from'',i5)', nttab, ntab(1+nbas)

C      print *, 'muck up iax3 for testing ...'
C      iax3(2) = 3
C      iax3(3) = 2

      end
      subroutine hoff3(nbas,iax3,iax,lphi,lmxa,n0,ips,nel,ioffb)
C- Make offset tables ioffb for 3C basis, ioffl for linked
      implicit none
      integer n0,nbas,nel,lphi(n0,2),ioffb(nbas,nel),
     .  iax3(nbas),iax(10,nbas),
C    .  ioffl(nbas+1),
     .  ips(3),ib,ie,is,nlmb,ibe,lmxa(nbas),nlb,nlbx,ioff

C     call awrit2('xx iax3 %n:1i',' ',80,6,nbas,iax3)

C --- Offsets table for hamiltonian ---
      ioffb(1,1) = 0
      ioff = 0
      ibe = 1
      do  10  ie = 1, nel
      do  10  ib = 1, nbas
      is = ips(iax(2,iax3(ib)))
C     call awrit3('xx iax3 %,3i  iax %,3i  ips %,3i',' ',80,6,
C    .  iax3(ib),iax(2,iax3(ib)),ips(iax(2,iax3(ib))))
      nlmb = (lphi(ie,is)+1)**2
      ibe = ibe+1
      ioffb(ibe,1) = ioff + nlmb
      ioff = ioffb(ibe,1)
   10 continue
C     call awrit2('xx ioffb %n:1i',' ',80,6,ibe,ioffb)
CC --- Offsets table for linking energy ---
C      ie = nel+1
C      ioffl(1) = 0
C      do  20  ib = 1, nbas
C      is = ips(iax3(ib))
C      nlmb = (lphi(ie,is)+1)**2
C      ioffl(ib+1) = ioffl(ib) + nlmb
C   20 continue
C
      end
