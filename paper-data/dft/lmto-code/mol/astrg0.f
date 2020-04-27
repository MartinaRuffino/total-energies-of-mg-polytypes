      subroutine astrg0(nphi,lphi,nel,el,nlma,jb,pos,n0,nbas,ips,alat,
     .  cg,jcg,indxcg,cy,nhs,nstate,ioffb,ixp,dxp,z,bz,gbz,gdz)
C- Grad structure constants * evec for LMTO forces (vectorizes)
C  All strx are expanded about site jb;
C  grad strx decomposed into kb for explicit decomposition of force.
      implicit none
      integer n0,nbas,nel,nhs,nphi(1),lphi(n0,1),jcg(1),indxcg(1),
     .  ips(1),ioffb(nbas,nel),ixp(1),nlma,nstate
      double precision el(*),pos(3,1),dr(3),cg(*),cy(*),alat,e,
     .  gbz(nlma,nstate,nbas,3,nel),z(nhs,nstate),dxp(*),
     .  gdz(nlma,nstate,3,nel),bz(nlma,nstate,nel)
      integer kb,ie,ks,jb,m,nlmb,nlmbx,nlmp,npow,nlmj,ih,nlmh,ll,iofz,
     .  ocf,oip,oikl,ojkl(0:20),ojj,ohl,ohd,os,osd,ogs,ogd
      real w(1)
      common /w/ w

      call strxsu(nlma,nphi,lphi,n0,0,nbas,ips,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      call defrr(os,  nlma*nlmbx)
      call defrr(osd, nlma*nlmbx)
      call defrr(ogs, nlma*nlmbx*3)
      call defrr(ogd, nlma*nlmbx*3)
      call dpzero(bz, nlma*nstate*nel)
      call dpzero(gdz,nlma*nstate*nel*3)

C --- kb loop over function centers; expand around (passed) jb ---
      do  10  kb = 1, nbas
      ks = ips(kb)
        do  14  m = 1, 3
   14   dr(m) = alat*(pos(m,kb)-pos(m,jb))
C  ...  Strx for all energies, this pair
        nlmh = (ll(nlmp)+2)**2
        call defrr(ohl,nlmh*nphi(ks))
        call defrr(ohd,nlmh*nphi(ks))
        call rstr0(nphi(ks),lphi(1,ks),el,nlmh,1,dr(1),dr(2),dr(3),
     .    ll(nlma)+1,0,w(ohl),w(ohd))
        do  12  ie = 1, nel
          nlmb = ioffb(kb+1,ie)-ioffb(kb,ie)
          if (nlmb <= 0) goto 12
          if (ixp(1) /= 0 .or. kb /= jb) then
            e = el(ie)
            ojj = ojkl(lphi(ie,ks))
            ih = nlmh*(ie-1)
            call hstrpg(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd),w(ogs),w(ogd))
            nlmj = ioffb(jb+1,ie)-ioffb(jb,ie)
            iofz = ioffb(kb,ie)+1
            call axtrg0(nlma,nlmb,nlmj,nhs,nbas,nstate,z(iofz,1),w(os),
     .       w(ogs),w(ogd),bz(1,1,ie),gbz(1,1,kb,1,ie),gdz(1,1,1,ie))
          endif
   12   continue
        call rlse(ohl)
   10 continue
      call rlse(ocf)
      end
      subroutine axtrg0(nlma,nlmb,nlmj,nhs,
     .  nbas,nstate,z,s,gs,gd,bz,gbz,gdz)
C- Grad strux * evec for proper block
      implicit none
      integer nlma,nlmb,nhs,nlmj,ix,nbas,nstate
      double precision s(nlma,nlmb),gs(nlma,nlmb,3),z(nhs,nstate),
     .  bz(nlma,nstate),gbz(nlma,nstate,nbas,3),
     .  gd(nlma,nlmb,3),gdz(nlma,nstate,3)

C --- b * z ---
      call dgemm('N','N',nlma,nstate,nlmb,1d0,s,nlma,z,nhs,1d0,bz,nlma)

C --- Grad b * z,  Grad b-dot z ---
      do  10  ix= 1, 3
        call dgemm('N','N',nlma,nstate,nlmb,1d0,gs(1,1,ix),nlma,
     .    z,nhs,0d0,gbz(1,1,1,ix),nlma)
        call dgemm('N','N',nlmj,nstate,nlmb,1d0,gd(1,1,ix),nlma,
     .    z,nhs,1d0,gdz(1,1,ix),nlma)
   10 continue

      end
