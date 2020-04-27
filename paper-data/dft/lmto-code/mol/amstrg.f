      subroutine amstrg(lphi,n0,ib,nlma,nhs,pos,nbas,ips,alat,el,nel,
     .  cg,jcg,indxcg,cy,gb,gbp)
C- Gradient of molecular strux constants and energy derivative
C  Analog of amstrp, but expand about one site only
      implicit none
      integer n0,nel,lphi(n0,1),jcg(1),indxcg(1),ips(1),ib,nbas,nhs,nlma
      double precision el(nel),alat,pos(3,1),dr(3),cg(*),cy(*),
     .  gb(nlma,nhs,3),gbp(nlma,nhs,3)
      integer ilb,kb,ks,m,nlmb,i,ie

      ilb = 1
C --- Loop over energies and sites which to expand ---
      do  10  ie = 1, nel
      do  10  kb = 1, nbas
      ks = ips(kb)
      nlmb = (lphi(ie,ks)+1)**2
      if (kb /= ib) then
        do  14  m = 1, 3
   14   dr(m) = alat*(pos(m,kb)-pos(m,ib))
        call mstrgg(el(ie),dr,nlma,nlmb,nlma,nhs,cg,indxcg,jcg,cy,
     .    gb(1,ilb,1),gbp(1,ilb,1))
      else
        do  30  i = 1, 3
        call dpzero(gb(1,ilb,i),nlma*nlmb)
   30   call dpzero(gbp(1,ilb,i),nlma*nlmb)
      endif
      ilb = ilb+nlmb
   10 continue

C      call prmx(gb(1,1,3),nlma,nlma,nhs)
C      stop
      end
