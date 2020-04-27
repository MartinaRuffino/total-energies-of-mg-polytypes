      subroutine astrgc(nphi,lphi,nel,el,nlma,jc,pos,n0,nbas,nclus,
     .  nclus3,iprs,alat,cg,jcg,indxcg,cy,nhs,nhl,nstate,ioffb,iofc,
     .  iofcl,iax,iax3,ixp,dxp,z,ceadd,bz,gbz,gdz,bzl,gbzl,gdzl)
C- Grad strux * evec for LMTO forces (cluster, linked basis)
C  All strx are expanded about cluster jc;
C  grad strx decomposed into ic for explicit decomposition of force.
      implicit none
      integer n0,nbas,nclus,nclus3,nel,nhs,nhl,nlma,nstate,
     .  nphi(1),lphi(n0,1),jcg(1),indxcg(1),iax(10,1),iax3(1),
     .  iprs(1),ioffb(nbas,nel),iofc(nclus,nel),iofcl(nclus+1),ixp(1)
      double precision el(1),pos(3,1),dr(3),cg(1),cy(1),alat,e,
     .  gbz(nlma,nstate,nbas,3,nel),z(nhs,nstate),dxp(1),ceadd(25,5,1),
     .  gdz(nlma,nstate,3,nel),bz(nlma,nstate,nel),
     .  gbzl(nlma,nstate,nbas,3,nel),
     .  gdzl(nlma,nstate,3,nel),bzl(nlma,nstate,nel)
      integer ic,ie,je,is,jc,m,nlmb,nlmbx,nlmp,npow,nlmj,ih,nlmh,ll,
     .  iofz,ocf,oip,oikl,ojkl(0:20),ojj,ohl,ohd,os,osd,ogs,ogd,owk,
     .  nlk,lmxaj,ib,l3c,ic3
      real w(1)
      common /w/ w
      character *80 strn

C --- Setup ---
      nlk = 0
      if (nhl > 0) nlk = 1
      lmxaj = ll(nlma)
      call strxsu(nlma,nphi,lphi,n0,0,nclus,iprs,cg,jcg,indxcg,
     .   nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      call defrr(os,  nlma*nlmbx)
      call defrr(osd, nlma*nlmbx)
      call defrr(ogs, nlma*nlmbx*3)
      call defrr(ogd, nlma*nlmbx*3)
      call defrr(owk, nlma*nlmbx)
      call dpzero(bz, nlma*nstate*nel)
      call dpzero(gdz,nlma*nstate*nel*3)

C --- ic loop over function centers; expand around jc ---
      do  10  ic = 1, nclus
        ib = iax(2,ic)
        is = iprs(ic)
C   ... l3c is 1 if ic is part of 3C cluster too
        l3c = 0
        do  21  ic3 = 1, nclus3
          if (ic == iax3(ic3)) then
            l3c = 1
            goto 22
          endif
   21   continue
   22   continue
C   ... Uncomment following line to get make b*z for all clusters!
C       l3c=1
        do  14  m = 1, 3
   14   dr(m) = alat*(pos(m,ic)-pos(m,jc))
C  ...  Strx for all energies, this pair
        nlmh = (ll(nlmp)+2)**2
        call defrr(ohl,nlmh*(nphi(is)+nlk))
        call defrr(ohd,nlmh*(nphi(is)+nlk))
        call rstr0(nphi(is)+nlk,lphi(1,is),el,nlmh,1,dr(1),dr(2),dr(3),
     .    lmxaj+1,0,w(ohl),w(ohd))
        do  12  ie = 1, nel
          nlmb = iofc(ic+1,ie)-iofc(ic,ie)
          if (nlmb <= 0) goto 12
          if (ixp(1) == 0 .and. ic == jc) then
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmh*(ie-1)
            call hstrpg(e,nlma,nlmb,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd),w(ogs),w(ogd))
            nlmj = iofc(jc+1,ie)-iofc(jc,ie)
            iofz = ioffb(ib,ie)+1
            call xxtrg0(nlma,nlmb,nlmj,nhs,l3c,nbas,nstate,z(iofz,1),
     .        w(os),w(ogs),w(ogd),bz(1,1,ie),gbz(1,1,ic,1,ie),
     .        gdz(1,1,1,ie))
          endif
   12   continue

C ---   b*z, grad b*z, bd*z grad bd*z for linked basis ---
        if (nhl > 0) then
          ie = nel+1
          nlmbx = iofcl(ic+1)-iofcl(ic)
          if (nlmbx <= 0) goto 16
          if (ixp(1) == 0 .and. ic == jc) then
          else
            e = el(ie)
            ojj = ojkl(lphi(ie,is))
            ih = nlmh*(ie-1)
            call hstrpg(e,nlma,nlmbx,nlmp,npow,ih,w(oikl),w(ojj),w(oip),
     .        w(ocf),w(ohl),w(ohd),w(os),w(osd),w(ogs),w(ogd))
            do  15  ie = 1, nel
              nlmb = iofc(ic+1,ie)-iofc(ic,ie)
              if (nlmb <= 0) goto 15
              if (ixp(1) == 0 .and. ic == jc) then
              else
                nlmj = 0
                do  17  je = 1, nel
   17           nlmj = max(nlmj,iofc(jc+1,je)-iofc(jc,je))
                iofz = ioffb(ib,ie)+1
                call xxtrgl(nlma,nlmbx,nlmb,nlmj,nhs,l3c,nbas,nstate,
     .            z(iofz,1),w(os),w(ogs),w(ogd),ceadd(1,ie,is),w(owk),
     .            bzl(1,1,ie),gbzl(1,1,ic,1,ie),gdzl(1,1,1,ie))
              endif
   15       continue
          endif
   16     continue
        endif

        call rlse(ohl)
   10 continue


      call rlse(ocf)

      end
      subroutine xxtrg0(nlma,nlmb,nlmj,nhs,l3c,
     .  nclus,nstate,z,s,gs,gd,bz,gbz,gdz)
C- Grad strux * evec for one block
      implicit none
      integer nlma,nlmb,nhs,nlmj,ix,nclus,nstate,l3c
      double precision s(nlma,nlmb),gs(nlma,nlmb,3),z(nhs,nstate),
     .  bz(nlma,nstate),gbz(nlma,nstate,nclus,3),
     .  gd(nlma,nlmb,3),gdz(nlma,nstate,3)

C --- Accumulate b * z from this block ---
      if (l3c == 1)
     .call dgemm('N','N',nlma,nstate,nlmb,1d0,s,nlma,z,nhs,1d0,bz,nlma)

C --- Grad b * z,  Grad b-dot z ---
      do  10  ix= 1, 3
        call dgemm('N','N',nlma,nstate,nlmb,1d0,gs(1,1,ix),nlma,
     .    z,nhs,0d0,gbz(1,1,1,ix),nlma)
        call dgemm('N','N',nlmj,nstate,nlmb,1d0,gd(1,1,ix),nlma,
     .    z,nhs,1d0,gdz(1,1,ix),nlma)
   10 continue

      end
      subroutine xxtrgl(nlma,nlmbx,nlmb,nlmj,nhs,l3c,
     .  nclus,nstate,z,s,gs,gd,ceadd,wk,bz,gbz,gdz)
C- Grad strux * evec for one block, linking basis
      implicit none
      integer nlma,nlmbx,nlmb,nhs,nlmj,nclus,nstate,ix,i,j,l3c
      double precision s(nlma,nlmbx),gs(nlma,nlmbx,3),z(nhs,nstate),
     .  bz(nlma,nstate),gbz(nlma,nstate,nclus,3),ceadd(nlmb),
     .  gd(nlma,nlmbx,3),gdz(nlma,nstate,3),wk(nlma,nlmb)

C --- b * ceadd * z ---
      if (l3c == 0) goto 12
      do  10  i = 1, nlma
      do  10  j = 1, nlmb
   10 wk(i,j) = s(i,j)*ceadd(j)
      call dgemm('N','N',nlma,nstate,nlmb,1d0,wk,nlma,z,nhs,1d0,bz,nlma)
   12 continue

C --- Grad b * ceadd * z,  Grad b-dot * ceadd * z ---
      do  20  ix= 1, 3
        do  24  i = 1, nlma
        do  24  j = 1, nlmb
   24   wk(i,j) = gs(i,j,ix)*ceadd(j)
        call dgemm('N','N',nlma,nstate,nlmb,1d0,wk,nlma,
     .    z,nhs,0d0,gbz(1,1,1,ix),nlma)
        do  28  i = 1, nlma
        do  28  j = 1, nlmb
   28   wk(i,j) = gd(i,j,ix)*ceadd(j)
        call dgemm('N','N',nlmj,nstate,nlmb,1d0,wk,nlma,
     .    z,nhs,1d0,gdz(1,1,ix),nlma)
   20 continue

      end
