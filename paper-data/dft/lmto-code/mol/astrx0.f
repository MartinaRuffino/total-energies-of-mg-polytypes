      subroutine astrx0(nphi,lphi,nel,el,lmxa,pos,n0,nbas,ips,alat,
     .  cg,jcg,indxcg,cy,nla,nhs,nhl,ioffa,ioffb,ioffl,
     .  ixp,dxp,b,bd,bl,bld)
C- Structure constants and energy derivatives for LMTO Hamiltonian
      implicit none
      integer i1,ib,ie,is,j1,jb,js,m,n0,nbas,nel,nhs,nla,nlma,nlmb,
     .  nlmbx,nlmp,npow,ocf,oip,oikl,ojkl(0:20),ojj,ohl,ohd,os,osd,
     .  id,jd,ndims,ndim,nlmj,ih,nhl
      integer nphi(1),lphi(n0,1),lmxa(1),jcg(1),indxcg(1),ips(1),
     .  ioffb(nbas,nel),ioffa(nbas+1),ioffl(nbas+1),ixp(1),nlk
      double precision el(1),pos(3,1),dr(3),cg(1),cy(1),
     .  b(nla,nhs),bl(nla,nhl),alat,e,dxp(1),bd(1),bld(1)
      real w(1)
      common /w/ w

C --- ib loop over function centers, jb over sites where expanded ---
      nlk = 0
      if (nhl > 0) nlk = 1
      do  11  jb = 1, nbas
      js = ips(jb)
      nlma = (lmxa(js)+1)**2
      call strxsu(nlma,nphi,lphi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      do  10  ib = 1, nbas
      is = ips(ib)
        do  14  m = 1, 3
   14   dr(m) = alat*(pos(m,ib)-pos(m,jb))
C  ...  Strx for all energies, this pair
        call defrr(ohl,nlmp*(nphi(is)+nlk))
        call defrr(ohd,nlmp*(nphi(is)+nlk))
        call rstr0(nphi(is)+nlk,lphi(1,is),el,nlmp,1,dr(1),dr(2),dr(3),
     .    lmxa(js),0,w(ohl),w(ohd))
        call defrr(os, nlma*nlmbx)
        call defrr(osd,nlma*nlmbx)
        ndims = 1
        do  12  ie = 1, nel
          ndim = ioffb(nbas+1,ie)-ioffb(1,ie)
          nlmb = ioffb(ib+1,ie)-ioffb(ib,ie)
          if (nlmb <= 0) goto 12
          j1 = ioffa(jb)+1
          i1 = ioffb(ib,ie)+1
          id = ioffb(ib,ie)-ioffb(1,ie)
          jd = ioffb(jb,ie)-ioffb(1,ie)
          if (ixp(1) == 0 .and. ib == jb) then
            call dpzero(w(os), nlma*nlmb)
            call dpzero(w(osd),nlma*nlmb)
          else
            e = el(ie)
C           call mstrxp(e,dr,b(j1,i1),bp(j1,i1),nlma,nlmb,nla,
C    .        cg,indxcg,jcg,cy)
            ojj = ojkl(lphi(ie,is))
            ih = nlmp*(ie-1)
            call hstrud(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd))
          endif
          nlmj = ioffb(jb+1,ie)-ioffb(jb,ie)
          call axtrx0(nlma,nlmb,nlmj,nla,ndim,id,jd,w(os),w(osd),
     .      b(j1,i1),bd(ndims))
   12   ndims = ndims + ndim**2

C --- Strux for linked basis ---
        if (nhl > 0) then
          ie = nel+1
          ndim = ioffl(nbas+1)-ioffl(1)
          nlmb = ioffl(ib+1)-ioffl(ib)
          if (nlmb <= 0) goto 16
          j1 = ioffa(jb)+1
          i1 = ioffl(ib)+1
          id = ioffl(ib)-ioffl(1)
          jd = ioffl(jb)-ioffl(1)
          if (ixp(1) == 0 .and. ib == jb) then
            call dpzero(w(os), nlma*nlmb)
            call dpzero(w(osd),nlma*nlmb)
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmp*(ie-1)
            call hstrud(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd))
          endif
          nlmj = ioffl(jb+1)-ioffl(jb)
          call axtrx0(nlma,nlmb,nlmj,nla,ndim,id,jd,w(os),w(osd),
     .      bl(j1,i1),bld)
   16     continue
        endif

        call rlse(ohl)
   10 continue
      call rlse(ocf)
   11 continue

C      call prmx(b,nla,nla,nhs)
C      call prmx(bp,nla,nla,nhs)

      end
      subroutine axtrx0(nlma,nlmb,nlmj,nla,nlb,id,jd,s,sd,b,bd)
C- Copy strx into proper block
      implicit none
      integer nlma,nlmb,nla,i,j,nlb,id,jd,nlmj
      double precision s(nlma,1),sd(nlma,1),b(nla,nlb),bd(nlb,nlb)

      do  20  i = 1, nlmb
      do  20  j = 1, nlma
   20 b(j,i) = s(j,i)

      do  30  i = 1, nlmb
      do  30  j = 1, nlmj
   30 bd(j+jd,i+id) = sd(j,i)

      end
