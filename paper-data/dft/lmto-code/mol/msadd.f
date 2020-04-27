      subroutine msadd(n0,nel,rsm,rmt,llink,elink,tlink,lphi,el,ceadd)
C- Coefficients for linked basis of Hankel functions, KE input
C  ceadd: coefficients of linking energy to be added for linked.
      implicit none
      integer n0,nel,lphi(nel)
      double precision elink,tlink(n0,nel),el(nel),p(0:10),pl(0:10),
     .  ceadd(25,5),add,rmt,rsm,wk(6),tol
      integer lmax,ie,l,m,ilm,ipr,idx(2),llink

      call getpr(ipr)

C --- Put elink to el(nel+1) ---
      el(nel+1) = elink
      lmax = -1
      if (llink /= 0) then
        do  19  ie = 1, nel
        if (dabs(elink-el(ie)) < 1d-5) call rx('msadd: elink eq el')
   19   lmax = max0(lmax,lphi(ie))
      endif
   18 lphi(nel+1) = lmax

      if (llink == 0) return

C --- Stop with error if elink equals any of the elmto ---
      do  40  ie = 1, nel
   40 if (dabs(elink-el(ie)) < 1d-5)
     .    call rx('elink too close to basis energy')

C --- Coefficients ceadd for linked basis ---
      tol = 1d-10
      call hansrg(rsm,10,1,lmax,elink,rmt**2,1,1,tol,idx,wk,0,pl)
      do  30  ie = 1, nel
        call hansrg(rsm,10,1,lphi(ie),el(ie),rmt**2,1,1,tol,idx,wk,
     .    0,p)
        lmax = lphi(ie)
        ilm = 0
        do  34  l = 0, lmax
        if (dabs(tlink(l+1,ie)-elink) < 1d-7)
     .      call rx('msadd: tlink too close to elink')
        add = -p(l)/pl(l)*(tlink(l+1,ie)-el(ie))/(tlink(l+1,ie)-elink)
        do  34  m = -l, l
        ilm = ilm+1
   34   ceadd(ilm,ie) = add
   30 continue

C --- Printout ---
      if (ipr < 50) return
  346 format('  l   ceadd ...')
      lmax = lphi(nel+1)
      do  44  l = 0, lmax
        ilm = (l+1)**2
        print 333, l, (ceadd(ilm,ie), ie = 1, nel)
  333   format(i3,5f12.7)
   44 continue

      end
