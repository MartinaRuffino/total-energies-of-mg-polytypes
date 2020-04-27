      subroutine astrgj(nphi,lphi,nel,el,nlma,pos,n0,nbas,jb,ips,alat,
     .  cg,jcg,indxcg,cy,nhs,nhl,ioffb,ioffl,ixp,dxp,b,bd,bl,bld,
     .  gb,gbd,gbl,gbld)
C- All strux and their gradient around site jb
      implicit none
      integer i1,ib,ie,is,jb,js,m,n0,nbas,nel,nhs,nlma,nlmb,
     .  nlmbx,nlmp,npow,ocf,oip,oikl,ojkl(0:20),ojj,ohl,ohd,
     .  nlmj,ih,nhl
      integer nphi(1),lphi(n0,1),lmxaj,jcg(1),indxcg(1),ips(1),
     .  ioffb(nbas,nel),ioffl(nbas+1),ixp(1),nlk,ll,nlmh
      double precision el(1),pos(3,1),dr(3),cg(1),cy(1),
     .  alat,e,dxp(1),b(nlma,1),bl(nlma,1),gb(nlma,1),gbl(nlma,1),
     .  bd(nlma,1),bld(nlma,1),gbd(nlma,1),gbld(nlma,1)
      real w(1)
      common /w/ w

C --- ib loop over function centers, jb = site where augmented ---
      nlk = 0
      if (nhl > 0) nlk = 1
      lmxaj = ll(nlma)
      call strxsu(nlma,nphi,lphi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      do  10  ib = 1, nbas
      is = ips(ib)
        do  14  m = 1, 3
   14   dr(m) = alat*(pos(m,ib)-pos(m,jb))
C  ...  Strx for all energies, this pair ...
        nlmh = (ll(nlmp)+2)**2
        call defrr(ohl,nlmh*(nphi(is)+nlk))
        call defrr(ohd,nlmh*(nphi(is)+nlk))
        call rstr0(nphi(is)+nlk,lphi(1,is),el,nlmh,1,dr(1),dr(2),dr(3),
     .    lmxaj+1,0,w(ohl),w(ohd))
        do  12  ie = 1, nel
          nlmb = ioffb(ib+1,ie)-ioffb(ib,ie)
          if (nlmb <= 0) goto 12
          i1 = ioffb(ib,ie)
          if (ixp(1) == 0 .and. ib == jb) then
            call dpzero(b(1,i1+1), nlma*nlmb)
            call dpzero(bd(1,i1+1),nlma*nlmb)
            call dpzero(gb(1,3*i1+1), nlma*nlmb*3)
            call dpzero(gbd(1,3*i1+1),nlma*nlmb*3)
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmh*(ie-1)
            call hstrpg(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),b(1,i1+1),bd(1,i1+1),
     .        gb(1,3*i1+1),gbd(1,3*i1+1))
          endif
   12   continue

C   ... Strux for linked basis ...
        if (nhl > 0) then
          ie = nel+1
          nlmb = ioffl(ib+1)-ioffl(ib)
          if (nlmb <= 0) goto 16
          i1 = ioffl(ib)
          if (ixp(1) == 0 .and. ib == jb) then
            call dpzero(bl(1,i1+1), nlma*nlmb)
            call dpzero(bld(1,i1+1),nlma*nlmb)
            call dpzero(gbl(1,3*i1+1), nlma*nlmb*3)
            call dpzero(gbld(1,3*i1+1),nlma*nlmb*3)
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmh*(ie-1)
            call hstrpg(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),bl(1,i1+1),bld(1,i1+1),
     .        gbl(1,3*i1+1),gbld(1,3*i1+1))
          endif
   16     continue
        endif

        call rlse(ohl)

   10 continue
      call rlse(ocf)

C      call prmx('b from astrxj',b,nlma,nlma,nhs)
C      call prmx('bd from astrxj',bd,nlma,nlma,nhs)

      end
