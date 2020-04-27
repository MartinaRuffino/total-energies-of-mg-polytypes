      subroutine astrxj(nphi,lphi,nel,el,nlma,pos,n0,nclus,jc,ips,alat,
     .  cg,jcg,indxcg,cy,llink,ioffc,iofcl,ixp,dxp,b,bd,bl,bld)
C- All strux around site jc for Hamiltonian
      implicit none
      integer i1,ic,ie,is,jc,m,nclus,nel,nlma,nlmb,nlk,ll,lmxaj,
     .  nlmbx,nlmp,npow,ocf,oip,oikl,ojkl(0:20),ojj,ohl,ohd,
     .  ih,ndimc
      integer n0,llink,nphi(1),lphi(n0,1),jcg(1),indxcg(1),ips(1),
     .  ioffc(nclus,nel),iofcl(nclus+1),ixp(1)
      double precision el(1),pos(3,nclus),dr(3),cg(1),cy(1),
     .  b(nlma,1),bl(nlma,1),alat,e,dxp(1),
     .  bd(nlma,1),bld(nlma,1)
      real w(1)
      common /w/ w

C --- ic loop over function centers, jc = site where augmented ---
      ndimc = ioffc(nclus+1,nel)
      nlk = 0
      if (llink > 0) nlk = 1
      lmxaj = ll(nlma)
      call strxsu(nlma,nphi,lphi,n0,0,nclus,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      do  10  ic = 1, nclus
        is = ips(ic)
        do  14  m = 1, 3
   14   dr(m) = alat*(pos(m,ic)-pos(m,jc))
C  ...  Strx for all energies, this pair
        call defrr(ohl,nlmp*(nphi(is)+nlk))
        call defrr(ohd,nlmp*(nphi(is)+nlk))
        call rstr0(nphi(is)+nlk,lphi(1,is),el,nlmp,1,dr(1),dr(2),dr(3),
     .    lmxaj,0,w(ohl),w(ohd))
        do  12  ie = 1, nel
          nlmb = ioffc(ic+1,ie)-ioffc(ic,ie)
          if (nlmb <= 0) goto 12
          i1 = ioffc(ic,ie)+1
          if (ixp(1) == 0 .and. ic == jc) then
            call dpzero(b(1,i1), nlma*nlmb)
            call dpzero(bd(1,i1),nlma*nlmb)
          else
            e = el(ie)
C           call mstrxp(e,dr,b(1,i1),bd(1,i1),nlma,nlmb,nlma,
C    .        cg,indxcg,jcg,cy)
            ojj = ojkl(lphi(ie,is))
            ih = nlmp*(ie-1)
            call hstrud(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),b(1,i1),bd(1,i1))
          endif
   12   continue

C --- Strux for linked basis ---
        if (llink > 0) then
          ie = nel+1
          nlmb = iofcl(ic+1)-iofcl(ic)
          if (nlmb <= 0) goto 16
          i1 = iofcl(ic)+1
          if (ixp(1) == 0 .and. ic == jc) then
            call dpzero(bl(1,i1), nlma*nlmb)
            call dpzero(bld(1,i1),nlma*nlmb)
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmp*(ie-1)
            call hstrud(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),bl(1,i1),bld(1,i1))
          endif
   16     continue
        endif

        call rlse(ohl)
   10 continue
      call rlse(ocf)

c      call prmx('b from astrxj',b,nlma,nlma,ndimc)
C      call prmx('bd from astrxj',bd,nlma,nlma,ndimc)

      end
