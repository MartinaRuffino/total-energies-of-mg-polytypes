      subroutine msetpl(n0,nel,spid,nspec,nbas,rmt,nphi,ips,llink,
     .  elink,tphi,lphi,el,nelx,nhl,ceadd)
C- Setup for linked basis, all species
C  Outputs:
C    nelx= nel+1 if not linked basis, else nel.
C    nhl: dimension of strx for linked basis,
C    ceadd: coefficients of linking energy to be added
      implicit none
      integer n0,nel,nelx,nspec,nbas,nhl,lphi(n0,nspec),ips(3),nphi(1)
      double precision elink,tphi(n0,5,nspec),el(nel),
     .  ceadd(25,5,1),rmt(nspec)
      integer ie,is,ib,l,ilm,ipr,llink
      character*8 spid(nspec)

      call getpr(ipr)
      nelx = nel
      if (llink == 1) nelx = nel+1

C --- Coefficients for linked basis ---
      call pshpr(1)
      do  10  is = 1, nspec
   10 call msadd(n0,nel,0d0,rmt(is),llink,elink,tphi(1,1,is),
     .  lphi(1,is),el,ceadd(1,1,is))
      call poppr

C --- Dimensions for strx of linked basis ---
      nhl = 0
      do  20  ib = 1, nbas
      is = ips(ib)
   20 nhl = nhl+(lphi(nel+1,is)+1)**2

C --- Printout ---
      if (ipr >= 20) print 333, elink
  333 format(/' msetpl: linked basis, elink=',f8.4)
      if (ipr < 30) return
      print 332
  332 format(' spec       el     l      tl       ceadd')
      do  40  is = 1, nspec
      do  40  ie = 1, nel
C      print 335, is, (tphi(ie,is), ie=1,nel)
C  335 format(1x,i4,'  tl',5f12.7)
      do  44  l = 0, lphi(ie,is)
        ilm = (l+1)**2
        if (ie == 1 .and. l == 0) print 334,
     .    spid(is), el(ie), l, tphi(l+1,ie,is), ceadd(ilm,ie,is)
        if (ie > 1 .and. l == 0) print 334,
     .    '    ',   el(ie), l, tphi(l+1,ie,is), ceadd(ilm,ie,is)
        if (l > 0) print 335,
     .                      l, tphi(l+1,ie,is), ceadd(ilm,ie,is)
  334   format(2x,a4,f10.4,i4,f10.4,f12.7)
  335   format(2x,4x,10x,  i4,f10.4,f12.7)
   44 continue
   40 continue
      end
