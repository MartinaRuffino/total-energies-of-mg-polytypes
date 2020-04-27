      subroutine hsmder(nbas,isp,nphi,lphi,ips,n0,nhs,nhl,nstate,
     .  nla,ioffa,ioffb,ioffl,
     .  el,nel,evl,z,hvecs,svecs,hi,si,ceadd,
C    .  b,bd,bl,bld,gb,gbd,gbl,gbld,
     .  gh,gj,puu,pus,pss,iofc,iofcl,iofc3,ioffp,
     .  iprs,ntab,iax,rtab,ntab3,iax3,rtab3,ixp,dxp,
     .  pos,alat,cg,jcg,indxcg,cy,wgt,ewt,f3)
C- Change in sumeV from virtual displacements (cluster h, linked basis)
C  This version truncates clusters following hsmatr.
      implicit none
      integer n0,nbas,isp,nel,nhs,nla,ips(1),nphi(1),lphi(n0,1),
     .  indxcg(1),jcg(1),iprs(1),ntab(1),iax(10,1),ntab3(1),iax3(7),
     .  nstate,ixp(1),iofc(1),iofc3(1),iofcl(nbas+1),
     .  ioffp(1),ioffb(nbas,nel),ioffl(nbas+1),ioffa(1),nhl
      double precision hvecs(nla,4,10),svecs(nla,4,6),hi(6),si(6),el(1),
     .  cg(1),cy(1),pos(3,1),gh(nla,2,nel),gj(nla,2,nel),z(nhs,nstate),
     .  puu(1),pus(1),pss(1),evl(nhs,2),dxp(1),alat,f3(3,nbas),
     .  ewt(nhs,isp),wgt,ceadd(25,5,1),rtab(3,1),rtab3(3,1)
C    .  ,b(1),bl(1),bd(1),bld(1),gb(1),gbl(1),gbd(1),gbld(1)

      double precision xxh,xxs,xx2,xx3,xxt
      integer ie,iv,ivs,je,ipr,ib,is,nlm1,nlm2,nlma,nx,k0,io,ic,ic3,
     .  nclus,nclus3,nlbx,owka,owkb,owkal,owkbl,j,
     .  obz,ogbz,ogdz,obzl,ogbzl,ogdzl
      integer nsp,lsp,ivl,jvl,nvl,ivls,jvls,nvls,il,ldiag,ndimi,ndimj,
     .  ndim3j,oiofa,ofe1,ndim3i,ofe2,
     .  opkk,opkj,opjk,opjj,opkk2,opkj2,opjk2,opjj2,
     .  opkk3,opkj3,opjk3,opjj3,opkkt,opkjt,opjkt,opjjt
      parameter (nx=49)
      logical llink
      real w(1)
      common /w/ w

C --- Setup ---
      call tcn('hsmder')
      call getpr(ipr)
      llink = nhl > 0
      nsp = lsp()+1
      call defrr(opkk, nx*nx)
      call defrr(opkj, nx*nx)
      call defrr(opjk, nx*nx)
      call defrr(opjj, nx*nx)
      if (nhl > 0) then
        call defrr(opkk2, nx*nx)
        call defrr(opkj2, nx*nx)
        call defrr(opjk2, nx*nx)
        call defrr(opjj2, nx*nx)
        call defrr(opkkt, nx*nx)
        call defrr(opkjt, nx*nx)
        call defrr(opjkt, nx*nx)
        call defrr(opjjt, nx*nx)
        call defrr(opkk3, nx*nx)
        call defrr(opkj3, nx*nx)
        call defrr(opjk3, nx*nx)
        call defrr(opjj3, nx*nx)
        il = nel+1
        nvl = ((nel+1)*(nel+2))/2
        nvls= nsp*(nvl-1)+isp
      endif

C --- For virtual displacement of site ib ---
      do  10  ib = 1, nbas
        is = ips(ib)
C ...   Dimensions of b,bd are locally b(nlma,*)
        nlma = ioffa(ib+1)-ioffa(ib)
        k0 = ioffa(ib)+1
C ...   iofc,iofl: table of offsets for cluster strux
        ic  = 1+ntab(ib)
        ic3 = 1+ntab3(ib)
        do  11  j = ntab(ib)+1, ntab(ib+1)
   11   iprs(j-ntab(ib)) = ips(iax(2,j))
        nclus = ntab(ib+1) - ntab(ib)
        nclus3 = ntab3(ib+1)-ntab3(ib)
        call defrr(oiofa, nclus+1)
        call hoffs(nclus,lphi,w(oiofa),n0,iprs,nel,iofc,
     .    w(oiofa),iofcl,nlbx)
        call hoff3(nclus3,iax3(ic3),iax(1,ic),lphi,w(oiofa),
     .    n0,ips,nel,iofc3)
        call rlse(oiofa)

C ...   Work arrays allocated
        if (nlma > nx) call rx('hsmder: nlma gt nx')
        call defrr(owka,     nlma*nstate)
        call defrr(owkb,     nlma*nstate)
        call defrr(ogbz,     nlma*nstate*nbas*3*nel)
        call defrr(ogdz,     nlma*nstate*3*nel)
        call defrr(obz,      nlma*nstate*nel)
        call dpzero(w(ogbz), nlma*nstate*nbas*3*nel)
        call dpzero(w(ogdz), nlma*nstate*3*nel)
        call dpzero(w(obz),  nlma*nstate*nel)
        if (nhl > 0) then
          call defrr(ogbzl,    nlma*nstate*nbas*3*nel)
          call defrr(ogdzl,    nlma*nstate*3*nel)
          call defrr(obzl,     nlma*nstate*nel)
          call dpzero(w(ogbzl),nlma*nstate*nbas*3*nel)
          call dpzero(w(ogdzl),nlma*nstate*3*nel)
          call dpzero(w(obzl), nlma*nstate*nel)
          call defrr(owkal,    nlma*nstate)
          call defrr(owkbl,    nlma*nstate)
        endif

C --- b * z, grad b * z, grad bdot * z, all sites with tails in ib ---
        call astrgc(nphi,lphi,nel,el,nlma,1,rtab(1,ic),n0,nbas,nclus,
     .    nclus3,iprs,alat,cg,jcg,indxcg,cy,nhs,nhl,nstate,ioffb,iofc,
     .    iofcl,iax(1,ic),iax3(ic3),0,dxp,z,ceadd,
     .    w(obz),w(ogbz),w(ogdz),w(obzl),w(ogbzl),w(ogdzl))

C --- b, grad b, bdot, grad bdot, sites in cluster with tails in ib ---
C        call astrgj(nphi,lphi,nel,el,nlma,rtab(1,ic),n0,nclus,1,iprs,
C     .    alat,cg,jcg,indxcg,cy,nhs,nhl,iofc,iofcl,0,dxp,b,bd,bl,bld,
C     .    gb,gbd,gbl,gbld)
C        call defrr(owk, nlma*25)
C        call xxmgbz(ib,nlma,nel,nhs,nbas,nclus,ioffb,iofc,iax(1,ic),
C     .    .false.,w,w,iprs,nstate,z,gb,gbd,w(ogbz),w(ogdz))
C        call xxmbz(ib,nlma,nel,nhs,nbas,nclus,nclus3,ioffb,iofc,
C     .    iax(1,ic),iax3(ic3),.false.,w,w,iprs,nstate,z,b,w(obz))
C        if (llink) then
C          call xxmgbz(ib,nlma,nel,nhs,nbas,nclus,ioffb,iofcl,
C     .      iax(1,ic),llink,ceadd,w(owk),iprs,nstate,z,gbl,gbld,
C     .      w(ogbzl),w(ogdzl))
C          call xxmbz(ib,nlma,nel,nhs,nbas,nclus,nclus3,ioffb,iofcl,
C     .      iax(1,ic),iax3(ic3),llink,ceadd,w(owk),iprs,nstate,z,bl,
C     .      w(obzl))
C        endif
C        call rlse(owk)

C ---   Loop over energy pairs ---
        iv = 0
        do  20  ie = 1, nel
C ...   Local cluster dimensions for this (ie,kb)
        ivl = (ie*(il+il-ie+1))/2
        ivls= nsp*(ivl-1)+isp
        ofe1 = iofc(1+(ie-1)*nclus)
        ndimi = iofc(1+ie*nclus)   - ofe1
        nlm1 = iofc(2+(ie-1)*nclus) - ofe1
        ndim3i= iofc3(1+ie*nclus3) - iofc3(1+(ie-1)*nclus3)

        do  30  je = ie, nel
C ...   Local cluster dimensions for this (je,kb)
        jvl = (je*(il+il-je+1))/2
        jvls= nsp*(jvl-1)+isp
        ofe2 = iofc(1+(je-1)*nclus)
        nlm2 = iofc(2+(je-1)*nclus)- ofe2
        ndimj = iofc(1+je*nclus)-ofe2
        ndim3j= iofc3(1+je*nclus3) - iofc3(1+(je-1)*nclus3)
        iv = iv+1
        ivs = nsp*(iv-1)+isp
        ldiag = 0
        if (ie == je) ldiag = 1

C ---   perturbation matrices for (ib,ie,je) ---
        io = nsp*ioffp(ib)+1
        if (isp == 2) io = io + nlma**2
        xxh = 0d0
        if (ie /= je) xxh = hi(iv)/(el(je)-el(ie))
        call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,ivs),nla,xxh,
     .    gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,je),gh(k0,2,je),
     .    gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,je),gj(k0,2,je),
     .    nlm1,nlm2,nlma,w(opkk),w(opkj),w(opjk),w(opjj))
        xxs = 0
        if (ie /= je) xxs = si(iv)/(el(je)-el(ie))

C ---   Forces for this (ie,je,ib) block, no linked basis ---
        if (nhl == 0) then
          if (ndimi > 0 .and. ndimj > 0)
     .      call xxrde3(nbas,nclus3,nhs,nla,nstate,ie,je,nel,ib,nlm1,
     .      nlm2,nlma,ioffb,iax(1,ic),iax3(ic3),w(opjj),
     .      svecs(k0,1,ivs),evl(1,isp),w(owka),w(owkb),
     .      w(obz),w(ogbz),w(ogdz),wgt,ewt(1,isp),f3)
          call xxrde2(nbas,nclus,nhs,nla,nstate,ie,je,nel,ib,nlm1,
     .      nlm2,nlma,ixp,ioffb,iax(1,ic),hi(iv),w(opkj),w(opjk),
     .      si(iv),svecs(k0,1,ivs),xxs,evl(1,isp),w(owka),w(owkb),
     .      z,w(ogbz),w(ogdz),wgt,ewt(1,isp),f3)
        else
C ...     pert matrices for (ie,linked) block
          xxh = hi(ivl)/(el(il)-el(ie))
          call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,ivls),nla,xxh,
     .    gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,il),gh(k0,2,il),
     .    gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,il),gj(k0,2,il),
     .    nlm1,nlm2,nlma,w(opkk2),w(opkj2),w(opjk2),w(opjj2))
C ...     pert matrices for (je,linked) block
          xxh = hi(jvl)/(el(il)-el(je))
          call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,jvls),nla,xxh,
     .    gh(k0,1,je),gh(k0,2,je),gh(k0,1,il),gh(k0,2,il),
     .    gj(k0,1,je),gj(k0,2,je),gj(k0,1,il),gj(k0,2,il),
     .    nlm2,nlm1,nlma,w(opkkt),w(opkjt),w(opjkt),w(opjjt))
C ...     pert matrices for (linked,linked) block
          xxh = 0d0
          call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,nvls),nla,xxh,
     .    gh(k0,1,il),gh(k0,2,il),gh(k0,1,il),gh(k0,2,il),
     .    gj(k0,1,il),gj(k0,2,il),gj(k0,1,il),gj(k0,2,il),
     .    nlm1,nlm2,nlma,w(opkk3),w(opkj3),w(opjk3),w(opjj3))
          xx2 = si(ivl)/(el(il)-el(ie))
          xxt = si(jvl)/(el(il)-el(je))
          xx3 = 0d0
          if (ndimi > 0 .and. ndimj > 0)
     .      call xxrdl3(nbas,nclus3,nhs,nla,nstate,ie,je,nel,ib,nlm1,
     .      nlm2,nlma,ixp,ioffb,iax(1,ic),iax3(ic3),
     .      svecs(k0,1,ivs),xxs,svecs(k0,1,ivls),xx2,svecs(k0,1,nvls),
     .      xx3,svecs(k0,1,jvls),xxt,w(opjj),w(opjj2),w(opjj3),
     .      w(opjjt),evl(1,isp),w(owka),w(owkb),w(owkal),w(owkbl),
     .      z,w(obz),w(ogbz),w(obzl),w(ogbzl),wgt,ewt(1,isp),f3)
          call xxrdl2(nbas,nclus,nhs,nla,nstate,ie,je,nel,ib,nlm1,nlm2,
     .      nlma,ixp,ioffb,iax(1,ic),hi(iv),hi(ivl),si(iv),si(nvl),
     .      svecs(k0,1,ivs),xxs,svecs(k0,1,ivls),xx2,svecs(k0,1,nvls),
     .      xx3,svecs(k0,1,jvls),xxt,ceadd(1,ie,is),ceadd(1,je,is),
     .      w(opkj),w(opjk),w(opkj2),w(opjk2),w(opkj3),w(opjk3),w(opkjt)
     .      ,w(opjkt),evl(1,isp),w(owka),w(owkb),w(owkal),w(owkbl),
     .      z,w(ogbz),w(ogdz),w(ogbzl),w(ogdzl),wgt,ewt(1,isp),f3)
        endif
   30   continue
        if (nhl /= 0) iv = iv+1
   20 continue
      call rlse(owka)

   10 continue
      call rlse(opkk)

C ---   Printout -------------------
      if(ipr >= 40) then
      write(6,289)
      do 55 ib=1,nbas
  55  write(6,288) ib,f3(1,ib),f3(2,ib),f3(3,ib)
  288 format(i6,3f13.6)
C 288 format(i6,3f22.14)
  289 format(/'    ib     eigval-force from asa hamiltonian')
      endif
      call tcx('hsmder')
      end
      subroutine xxrde2(nbas,nclus,nhs,nla,nstate,ie,je,nel,ib,
     .  nlm1,nlm2,nlma,ixp,ioffb,iax,hi,pkj,pjk,si,svec,xxx,evl,
     .  wk2a,wk2b,z,gbz,gdz,wgt,ewt,f3)
C- Change in eV sum from displacing ib, one pair of energies, 2C terms
      implicit none
      integer nbas,nclus,nhs,ib,nstate,nla,nel,nlm1,nlm2,nlma,ie,je,
     .  ioffb(nbas,1),iax(10,1),ixp(1)
      double precision pkj(nlm1,nlma),xxx,sum,sumi,wgt,ewt(1),evl(1),
     .  pjk(nlma,nlm2),hi,si,svec(nla,4),f3(3,nbas)
      double precision z(nhs,1),wk2a(nlma,nstate),wk2b(nlma,nstate),
     .  gbz(nlma,nstate,nbas,3,nel),gdz(nlma,nstate,3,nel),ddot
      integer ix,jb,ioffz,nlmk,is,jlm,jc

C --- Contract pkj,pjk blocks with z ---
      call dgemm('T','N',nlma,nstate,nlm1,1d0,
     .  pkj,nlm1,z(1+ioffb(ib,ie),1),nhs,0d0,wk2b,nlma)
      do  66  is = 1, nstate
      do 166  jlm = nlm1+1, nlma
  166 wk2b(jlm,is) = wk2b(jlm,is) * ewt(is)
      do  66  jlm = 1, nlm1
   66 wk2b(jlm,is) = (wk2b(jlm,is) - z(jlm+ioffb(ib,ie),is)*evl(is)*
     .    (svec(jlm,2)+xxx)) * ewt(is)
      call dgemm('N','N',nlma,nstate,nlm2,1d0,
     .  pjk,nlma,z(1+ioffb(ib,je),1),nhs,0d0,wk2a,nlma)
      do  68  is = 1, nstate
      do 168  jlm = nlm2+1, nlma
  168 wk2a(jlm,is) = wk2a(jlm,is) * ewt(is)
      do  68  jlm = 1, nlm2
   68 wk2a(jlm,is) = (wk2a(jlm,is) - z(jlm+ioffb(ib,je),is)*evl(is)*
     .    (svec(jlm,3)-xxx)) * ewt(is)

      do  10  ix = 1, 3

C ---   2C terms <ib,ie | grad jb,je> ---
        if (nlm1 > 0) then
          do  20  jc = 2, nclus
            jb = iax(2,jc)
            if (jb > ib .or. ie /= je) then
              sum = ddot(nstate*nlma,wk2b,1,gbz(1,1,jc,ix,je),1)
              sum = sum*2
              f3(ix,ib) = f3(ix,ib) + sum*wgt
              f3(ix,jb) = f3(ix,jb) - sum*wgt
C              call prf('2C ',ix,ie,ib,je,jb,wk2b,gbz(1,1,jc,ix,je),
C     .          sum,f3)
            endif
   20     continue
        endif

C ---   2C terms <ib,je | grad jb,ie> ---
        if (nlm2 > 0) then
          do  30  jc = 2, nclus
            jb = iax(2,jc)
            if (jb < ib .or. ie /= je) then
              sum = ddot(nstate*nlma,wk2a,1,gbz(1,1,jc,ix,ie),1)
              sum = sum*2
              f3(ix,ib) = f3(ix,ib) + sum*wgt
              f3(ix,jb) = f3(ix,jb) - sum*wgt
C              call prf('2Cx',ix,ie,ib,je,jb,wk2a,gbz(1,1,jc,ix,ie),
C     .          sum,f3)
            endif
   30     continue
        endif

C ---   grad bdot ---
        if (nlm1 > 0 .and. ie == je) then
          sum = 0
          nlmk = ioffb(ib+1,ie)-ioffb(ib,ie)
          ioffz = ioffb(ib,ie)
          do  74  is = 1, nstate
          sumi = 0
          do  75  jlm = 1, nlmk
   75     sumi = sumi + z(jlm+ioffz,is)*gdz(jlm,is,ix,ie)
   74     sum = sum + sumi*(hi-evl(is)*si)*ewt(is)
          sum = sum*2
          f3(ix,ib) = f3(ix,ib) + sum*wgt
C          call prf('dot',ix,ie,ib,je,0,z,gdz,sum,f3)
        endif

   10 continue

      end
      subroutine xxrde3(nbas,nclus3,nhs,nla,nstate,ie,je,nel,ib,nlm1,
     .  nlm2,nlma,ioffb,iax,iax3,pjj,svec,evl,wk3a,wk3b,
     .  bz,gbz,gdz,wgt,ewt,f3)
C- Change in eV sum from displacing ib, one pair of energies, 3C terms
      implicit none
      integer nbas,nclus3,nhs,ib,nstate,nla,nel,nlm1,nlm2,nlma,
     .  ie,je,ioffb(nbas,1),iax(10,1),iax3(1)
      double precision sum,wgt,ewt(1),evl(1),
     .  pjj(nlma,nlma),svec(nla,4),f3(3,nbas)
      double precision bz(nlma,nstate,nel),
     .  gbz(nlma,nstate,nbas,3,nel),
     .  gdz(nlma,nstate,3,nel),wk3a(nlma,nstate),wk3b(nlma,nstate),ddot
      integer ix,jc,jc3,jb,is,jlm

C --- Contract pjj blocks with bz ---
      call dgemm('N','N',nlma,nstate,nlma,1d0,pjj,nlma,
     .  bz(1,1,je),nlma,0d0,wk3a,nlma)
      do  60  is = 1, nstate
      do  60  jlm = 1, nlma
   60 wk3a(jlm,is) = wk3a(jlm,is) -
     .               bz(jlm,is,je)*svec(jlm,4)*evl(is)
      if (ie == je) then
        do  62  is = 1, nstate
        do  62  jlm = 1, nlma
   62   wk3b(jlm,is) = wk3a(jlm,is)
      else
        call dgemm('T','N',nlma,nstate,nlma,1d0,
     .    pjj,nlma,bz(1,1,ie),nlma,0d0,wk3b,nlma)
        do  64  is = 1, nstate
        do  64  jlm = 1, nlma
   64   wk3b(jlm,is) = (wk3b(jlm,is) - evl(is) * bz(jlm,is,ie)*
     .      svec(jlm,4)) * ewt(is)
      endif

      do  10  ix = 1, 3

C ---   z+ grad b+(ie) pjj b(je) z ---
        do  40  jc3 = 2, nclus3
          jc = iax3(jc3)
          jb = iax(2,jc)
          sum = ddot(nstate*nlma,wk3a,1,gbz(1,1,jc,ix,ie),1)
          sum = sum*2
          f3(ix,ib) = f3(ix,ib) + sum*wgt
          f3(ix,jb) = f3(ix,jb) - sum*wgt
C          call prf('3C ',ix,ie,ib,je,jb,wk3a,gbz(1,1,jc,ix,ie),
C     .      sum,f3)
   40   continue

C ---  z+ b+(ie) pjj grad b(je) z  ---
       if (ie /= je) then
         do  50  jc3 = 2, nclus3
           jc = iax3(jc3)
           jb = iax(2,jc)
           sum = ddot(nstate*nlma,wk3b,1,gbz(1,1,jc,ix,je),1)
           sum = sum*2
           f3(ix,ib) = f3(ix,ib) + sum*wgt
           f3(ix,jb) = f3(ix,jb) - sum*wgt
C           call prf('3Cx',ix,ie,ib,je,jb,wk3b,gbz(1,1,jc,ix,je),
C     .       sum,f3)
   50     continue
        endif

   10 continue

      end
      subroutine xxrdl2(nbas,nclus,nhs,nla,nstate,ie,je,nel,ib,
     .  nlm1,nlm2,nlma,ixp,ioffb,iax,hi,hil,si,sil,
     .  svec,xxs,svc2,xx2,svc3,xx3,svct,xxt,addi,addj,
     .  pkj,pjk,pkj2,pjk2,pkj3,pjk3,pkjt,pjkt,
     .  evl,wk2a,wk2b,wk2al,wk2bl,z,gbz,gdz,gbzl,gdzl,wgt,ewt,f3)
C- Change in sumeV from ie,je and displacing ib, linked basis
      implicit none
      integer nbas,nclus,nhs,ib,nstate,nla,nel,nlm1,nlm2,nlma,ie,je,
     .  ioffb(nbas,1),iax(10,1),ixp(1)
      double precision sum,sumi,wgt,ewt(1),evl(1),
     .  pkj (nlm1,nlma),pjk (nlma,nlm2),hi,hil,si,sil,
     .  pkj2(nlm1,nlma),pjk2(nlma,nlm2),
     .  pkjt(nlm2,nlma),pjkt(nlma,nlm2),svct(nla,4),xxt,
     .  pkj3(nlm1,nlma),pjk3(nlma,nlm2),
     .  svec(nla,4),xxs,svc2(nla,4),xx2,svc3(nla,4),xx3,
     .  f3(3,nbas),addi(nlm1),addj(nlm2)
      double precision
     .  gbz(nlma,nstate,nbas,3,nel),gbzl(nlma,nstate,nbas,3,nel),
     .  z(nhs,1),gdz(nlma,nstate,3,nel),gdzl(nlma,nstate,3,nel),
     .  wk2a(nlma,nstate),wk2b(nlma,nstate),
     .  wk2al(nlma,nstate),wk2bl(nlma,nstate),ddot
      integer ix,jc,jb,ioffz,nlmk,is,ilm,jlm,ila

C --- Contract pkj,pjk blocks with z ---
      do  21  ila = 1, nlma
      do  21  ilm = 1, nlm1
   21 pkj(ilm,ila) = pkj(ilm,ila) + addi(ilm)*pjkt(ila,ilm)
      call dgemm('T','N',nlma,nstate,nlm1,1d0,
     .  pkj,nlm1,z(1+ioffb(ib,ie),1),nhs,0d0,wk2b,nlma)
      do  66  is = 1, nstate
      do 166  ilm = nlm1+1, nlma
  166 wk2b(ilm,is) = wk2b(ilm,is) * ewt(is)
      do  66  ilm = 1, nlm1
   66 wk2b(ilm,is) = (wk2b(ilm,is) - z(ilm+ioffb(ib,ie),is)*evl(is)*
     .    ((svec(ilm,2)+xxs) + addi(ilm)*(svct(ilm,3)-xxt))) * ewt(is)
      do  22  ila = 1, nlma
      do  22  ilm = 1, nlm1
   22 pkj3(ilm,ila) = addi(ilm)*pkj3(ilm,ila) + pkj2(ilm,ila)
      call dgemm('T','N',nlma,nstate,nlm1,1d0,
     .  pkj3,nlm1,z(1+ioffb(ib,ie),1),nhs,0d0,wk2bl,nlma)
      do  67  is = 1, nstate
      do 167  ilm = nlm1+1, nlma
  167 wk2bl(ilm,is) = wk2bl(ilm,is) * ewt(is)
      do  67  ilm = 1, nlm1
   67 wk2bl(ilm,is) = (wk2bl(ilm,is) - z(ilm+ioffb(ib,ie),is)*evl(is)*
     .    (svc2(ilm,2)+xx2 + addi(ilm)*(svc3(ilm,2)+xx3))) * ewt(is)
      do  31  ilm = 1, nlm2
      do  31  ila = 1, nlma
   31 pjk(ila,ilm) = pjk(ila,ilm) + pjk2(ila,ilm)*addj(ilm)
      call dgemm('N','N',nlma,nstate,nlm2,1d0,
     .  pjk,nlma,z(1+ioffb(ib,je),1),nhs,0d0,wk2a,nlma)
      do  68  is = 1, nstate
      do 168  jlm = nlm2+1, nlma
  168 wk2a(jlm,is) = wk2a(jlm,is) * ewt(is)
      do  68  jlm = 1, nlm2
   68 wk2a(jlm,is) = (wk2a(jlm,is) - z(jlm+ioffb(ib,je),is)*evl(is)*
     .    (svec(jlm,3)-xxs + addj(jlm)*(svc2(jlm,3)-xx2))) * ewt(is)
      do  33  ilm = 1, nlm2
      do  33  ila = 1, nlma
   33 pjk3(ila,ilm) = pjk3(ila,ilm)*addj(ilm) + pkjt(ilm,ila)
      call dgemm('N','N',nlma,nstate,nlm2,1d0,
     .  pjk3,nlma,z(1+ioffb(ib,je),1),nhs,0d0,wk2al,nlma)
      do  69  is = 1, nstate
      do 169  jlm = nlm2+1, nlma
  169 wk2al(jlm,is) = wk2al(jlm,is) * ewt(is)
      do  69  jlm = 1, nlm2
   69 wk2al(jlm,is) = (wk2al(jlm,is) - z(jlm+ioffb(ib,je),is)*evl(is)*
     .    (svct(jlm,2)+xxt + addj(jlm)*(svc3(jlm,3)-xx3))) * ewt(is)

      do  10  ix = 1, 3

C ---   2C terms <ib,ie+le | grad jb,je> ---
        if (nlm1 > 0) then
          do  20  jc = 2, nclus
            jb = iax(2,jc)
            if (jb > ib .or. ie /= je) then
              sum = ddot(nstate*nlma,wk2b,1,gbz(1,1,jc,ix,je),1) +
     .              ddot(nstate*nlma,wk2bl,1,gbzl(1,1,jc,ix,je),1)
              sum = sum*2
              f3(ix,ib) = f3(ix,ib) + sum*wgt
              f3(ix,jb) = f3(ix,jb) - sum*wgt
C              call prf('2C ',ix,ie,ib,je,jb,wk2b,gbz(1,1,jc,ix,je),
C     .          sum,f3)
            endif
   20     continue
        endif

C ---   2C terms <ib,je+le | grad jb,ie> ---
        if (nlm2 > 0) then
          do  30  jc = 2, nclus
            jb = iax(2,jc)
            if (jb < ib .or. ie /= je) then
              sum = ddot(nstate*nlma,wk2a,1,gbz(1,1,jc,ix,ie),1) +
     .              ddot(nstate*nlma,wk2al,1,gbzl(1,1,jc,ix,ie),1)
              sum = sum*2
              f3(ix,ib) = f3(ix,ib) + sum*wgt
              f3(ix,jb) = f3(ix,jb) - sum*wgt
C              call prf('2Cx',ix,ie,ib,je,jb,wk2a,gbz(1,1,jc,ix,ie),
C     .          sum,f3)
            endif
   30     continue
        endif

C ---   grad bdot ---
        if (nlm1 <= 0) goto 10
        if (ie == je) then
          sum = 0
          do  74  is = 1, nstate
          nlmk = ioffb(ib+1,ie)-ioffb(ib,ie)
          ioffz = ioffb(ib,ie)
          sumi = 0
          do  75  jlm = 1, nlmk
   75     sumi = sumi + z(jlm+ioffz,is)*gdz(jlm,is,ix,ie)
   74     sum = sum + sumi*(hi-evl(is)*si)*ewt(is)
          sum = sum*2
          f3(ix,ib) = f3(ix,ib) + sum*wgt
C          call prf('dot',ix,ie,ib,je,0,z(1+ioffz,1),
C     .      gdz(1,1,ix,ie),sum,f3)
        endif

C ...   grad bdot contribution from linked basis
        sum = 0
        do  76  is = 1, nstate
          nlmk = ioffb(ib+1,ie)-ioffb(ib,ie)
          ioffz = ioffb(ib,ie)
          sumi = 0
          do  77  jlm = 1, nlmk
   77     sumi = sumi + z(jlm+ioffz,is)*addi(jlm)*gdzl(jlm,is,ix,je)
          sum = sum + sumi*(hil-evl(is)*sil)*ewt(is)
          if (ie /= je) then
            nlmk = ioffb(ib+1,je)-ioffb(ib,je)
            ioffz = ioffb(ib,je)
            sumi = 0
            do  78  jlm = 1, nlmk
   78       sumi = sumi + z(jlm+ioffz,is)*addj(jlm)*gdzl(jlm,is,ix,ie)
            sum = sum + sumi*(hil-evl(is)*sil)*ewt(is)
          endif
   76   continue
        sum = sum*2
        f3(ix,ib) = f3(ix,ib) + sum*wgt
C       print 346, ix,'dtl ib,ie,je,sum=',ib,ie,je,sum,f3(1,2)
C 346   format(1x,i1,1x,a,3i4,2g24.16)

   10 continue
      end
      subroutine xxrdl3(nbas,nclus3,nhs,nla,nstate,ie,je,nel,ib,nlm1,
     .  nlm2,nlma,ixp,ioffb,iax,iax3,svec,xxs,svc2,xx2,svc3,
     .  xx3,svct,xxt,pjj,pjj2,pjj3,pjjt,evl,wk3a,wk3b,wk3al,wk3bl,
     .  z,bz,gbz,bzl,gbzl,wgt,ewt,f3)
C- Change in 3C sumeV from ie,je displacing ib, linked basis
      implicit none
      integer nbas,nclus3,nhs,ib,nstate,nla,nel,nlm1,nlm2,nlma,
     .  ie,je,ioffb(nbas,1),iax(10,1),iax3(1),ixp(1)
      double precision sum,wgt,ewt(1),evl(1),
     .  pjj(nlma,nlma),pjj2(nlma,nlma),pjjt(nlma,nlma),pjj3(nlma,nlma),
     .  svec(nla,4),xxs,svc2(nla,4),xx2,svc3(nla,4),xx3,svct(nla,4),xxt,
     .  f3(3,nbas)
      double precision bz(nlma,nstate,nel),bzl(nlma,nstate,nel),
     .  gbz(nlma,nstate,nbas,3,nel),gbzl(nlma,nstate,nbas,3,nel),
     .  wk3a(nlma,nstate),wk3b(nlma,nstate),z(nhs,1),
     .  wk3al(nlma,nstate),wk3bl(nlma,nstate),ddot
      integer ix,jc,jc3,jb,is,jlm

C --- Contract pjj blocks with bz ---
      call dgemm('N','N',nlma,nstate,nlma,1d0,pjj,nlma,
     .  bz(1,1,je),nlma,0d0,wk3a,nlma)
      call dgemm('N','N',nlma,nstate,nlma,1d0,pjj2,nlma,
     .  bzl(1,1,je),nlma,1d0,wk3a,nlma)
      do  60  is = 1, nstate
      do  60  jlm = 1, nlma
   60 wk3a(jlm,is) = (wk3a(jlm,is) - evl(is) * (svec(jlm,4)*
     .    bz(jlm,is,je) + svc2(jlm,4)*bzl(jlm,is,je))) * ewt(is)
      call dgemm('T','N',nlma,nstate,nlma,1d0,pjjt,nlma,
     .  bz(1,1,je),nlma,0d0,wk3al,nlma)
      call dgemm('N','N',nlma,nstate,nlma,1d0,pjj3,nlma,
     .  bzl(1,1,je),nlma,1d0,wk3al,nlma)
      do  61  is = 1, nstate
      do  61  jlm = 1, nlma
   61 wk3al(jlm,is) = ( wk3al(jlm,is) - evl(is) * (svct(jlm,4)*
     .    bz(jlm,is,je) + svc3(jlm,4)*bzl(jlm,is,je)) ) * ewt(is)

      call dgemm('T','N',nlma,nstate,nlma,1d0,
     .  pjj,nlma,bz(1,1,ie),nlma,0d0,wk3b,nlma)
      call dgemm('N','N',nlma,nstate,nlma,1d0,
     .  pjjt,nlma,bzl(1,1,ie),nlma,1d0,wk3b,nlma)
      do  62  is = 1, nstate
      do  62  jlm = 1, nlma
   62 wk3b(jlm,is) = (wk3b(jlm,is) - evl(is) * (bz(jlm,is,ie)*
     .    svec(jlm,4) + bzl(jlm,is,ie)*svct(jlm,4))) * ewt(is)
      call dgemm('T','N',nlma,nstate,nlma,1d0,
     .  pjj2,nlma,bz(1,1,ie),nlma,0d0,wk3bl,nlma)
      call dgemm('T','N',nlma,nstate,nlma,1d0,
     .  pjj3,nlma,bzl(1,1,ie),nlma,1d0,wk3bl,nlma)
      do  63  is = 1, nstate
      do  63  jlm = 1, nlma
   63 wk3bl(jlm,is) = (wk3bl(jlm,is) - evl(is) * (bz(jlm,is,ie)*
     .    svec(jlm,4) + bzl(jlm,is,ie)*svct(jlm,4))) * ewt(is)

      do  10  ix = 1, 3

C ---   z+ grad b+(ie) (pjj b(je) + pjj2 b(le)) z ---
C and   z+ grad b+(le) (pjjt b(je) + pjj3 b(le)) z
        do  40  jc3 = 2, nclus3
          jc = iax3(jc3)
          jb = iax(2,jc)
          sum = ddot(nlma*nstate,wk3a, 1,gbz(1,1,jc,ix,ie), 1) +
     .          ddot(nlma*nstate,wk3al,1,gbzl(1,1,jc,ix,ie),1)
          sum = sum*2
          f3(ix,ib) = f3(ix,ib) + sum*wgt
          f3(ix,jb) = f3(ix,jb) - sum*wgt
C          call prf('3C ',ix,ie,ib,je,jb,wk3a,gbz(1,1,jc,ix,ie),
C     .      sum,f3)
   40   continue

C ---   z+ (b+(ie) pjj  +  b+(le) pjjt)  grad b(je) z  ---
C and   z+ (b+(ie) pjj2 +  b+(le) pjj3)  grad b(le) z
        if (ie /= je) then
         do  50  jc3 = 2, nclus3
           jc = iax3(jc3)
           jb = iax(2,jc)
           sum = ddot(nlma*nstate,wk3b, 1,gbz(1,1,jc,ix,je), 1) +
     .          ddot(nlma*nstate,wk3bl,1,gbzl(1,1,jc,ix,je),1)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C            call prf('3Cx',ix,ie,ib,je,jb,wk3b,gbz(1,1,jc,ix,je),
C     .        sum,f3)
   50     continue
        endif

   10 continue

      end
      subroutine xxmgbz(kb,nlma,nel,nhs,nbas,nclus,ioffb,iofc,iax,
     .  llink,ceadd,wk,iprs,nstate,z,gb,gd,gbz,gdz)
C- grad strux * evec for all strx centered about kb
C  llink t: use strux for ie=1; but scale by ceadd(ie)
      implicit none
      integer kb,nlma,nel,nhs,nbas,nclus,
     .  ioffb(nbas,nel),iofc(nclus,1),iax(10,1),iprs(1),nstate
      logical llink
      double precision gb(nlma,1),gd(nlma,1),z(nhs,nstate),
     .  gbz(nlma,nstate,nbas,3,1),gdz(nlma,nstate,3,1),
     .  ceadd(25,5,1),wk(nlma,1)
      integer ie,ic,ib,nlmb,nlmkb,iofz,i1,ii,ix
C     character *80 strn
      integer i,j,is,le

      call dpzero(gdz, nlma*nstate*3*nel)

      do  10  ie = 1, nel
      le = ie
      if (llink) le=1
      nlmkb = iofc(2,le)-iofc(1,le)
      do  10  ic = 1, nclus
        ib = iax(2,ic)
        is = iprs(ic)
        nlmb = iofc(ic+1,le)-iofc(ic,le)
        if (nlmb <= 0) goto 10
        if (ib /= kb) then
          iofz = ioffb(ib,ie)+1
          i1 = iofc(ic,le)
C ---     Grad b * z, grad bdot * z or grad bl * ceadd * z ---
          ii = 1
          do  12  ix= 1, 3
            if (llink) then
              do  16  j = 1, nlmb
              do  16  i = 1, nlma
   16         wk(i,j) = gd(i,j-1+3*i1+ii)*ceadd(j,ie,is)
              call dgemm('N','N',nlmkb,nstate,nlmb,1d0,wk,
     .          nlma,z(iofz,1),nhs,1d0,gdz(1,1,ix,ie),nlma)
              do  14  j = 1, nlmb
              do  14  i = 1, nlma
   14         wk(i,j) = gb(i,j-1+3*i1+ii)*ceadd(j,ie,is)
              call dgemm('N','N',nlma,nstate,nlmb,1d0,wk,
     .          nlma,z(iofz,1),nhs,0d0,gbz(1,1,ic,ix,ie),nlma)
            else
              call dgemm('N','N',nlma,nstate,nlmb,1d0,gb(1,3*i1+ii),
     .          nlma,z(iofz,1),nhs,0d0,gbz(1,1,ic,ix,ie),nlma)
              call dgemm('N','N',nlmkb,nstate,nlmb,1d0,gd(1,3*i1+ii),
     .          nlma,z(iofz,1),nhs,1d0,gdz(1,1,ix,ie),nlma)
            endif
            ii = ii+nlmb
   12     continue

C          if (.not. llink) goto 10
c          ii = 1
C          print *, 'xgb for kb,ie,ib,ic=', kb,ie,ib,ic
C          call prmx(' ',wk,nlma,nlma,nlmb)
C          print *, 'xgb for kb,ie,ib,ic=', kb,ie,ib,ic
C          call prmx(' ',gb(1,3*i1+ii),nlma,nlma,nlmb)
C          print *, 'z for ie,ib=', ie,ib
C          call prmx(' ',z(iofz,1),nlmb,nlmb,nstate)
c          print *, 'gbz for kb,ie,ib,ic=', kb,ie,ib,ic
c         call prmx(' ',gbz(1,1,ic,1,ie),nlma,nlma,nstate)
C          call prmx(' ',gbz(1,1,ie,1,1),nlma,nlma,nstate)
C          strn = ' '
C          call awrit3('diff out.dat xx/bz.%i%i%i',strn,80,0,kb,ie,ib)
c          call awrit3('diff out.dat xx/xgbz.%i%i%i',strn,80,0,kb,ie,ib)
c          call awrit3('mc -qr out.dat -qr xx/xgbz.%i%i%i -- -px',strn,
c     .      80,0,kb,ie,ib)
C          print *, 'gdz for kb,ie,ib,ic=', kb,ie,ib,ic
C          call prmx(' ',gdz(1,1,1,ie),nlma,nlmkb,nstate)
C          strn = ' '
C          call awrit3('diff out.dat xx/xgbz.%i%i%i',strn,80,0,kb,ie,ib)
C          call fsystm(strn,i)
C          call awrit1('shell command '//strn//'%a returned %i',
C     .      strn,len(strn),-6,i)
        endif
   10 continue

      end
      subroutine xxmbz(kb,nlma,nel,nhs,nbas,nclus,nclus3,ioffb,iofc,
     .  iax,iax3,llink,ceadd,wk,iprs,nstate,z,b,bz)
C- strx * evec for all strx centered about kb
C  llink t: use strux for ie=1; but scale by ceadd(ie)
      implicit none
      integer kb,nlma,nel,nhs,nbas,nclus,nclus3,ioffb(nbas,nel),
     .  iofc(nclus,1),iax(10,1),iprs(1),iax3(1),nstate
      logical llink
      double precision b(nlma,1),z(nhs,nstate),
     .  bz(nlma,nstate,1),ceadd(25,5,1),wk(nlma,1)
      integer is,le,ie,ic,ic3,ib,nlmb,nlmkb,iofz,i1
C     character *80 strn
      integer i,j
      call dpzero(bz, nlma*nstate*nel)

      do  10  ie = 1, nel
      le = ie
      if (llink) le=1
      nlmkb = iofc(2,le)-iofc(1,le)
      do  10  ic3 = 2, nclus3
        ic = iax3(ic3)
        ib = iax(2,ic)
        is = iprs(ic)
        nlmb = iofc(ic+1,le)-iofc(ic,le)
        if (nlmb <= 0) goto 10
        if (ib /= kb) then
          iofz = ioffb(ib,ie)+1
          i1 = iofc(ic,le)
C ---     b * z  or  b * ceadd * z ---
          if (llink) then
            do  18  j = 1, nlmb
            do  18  i = 1, nlma
   18       wk(i,j) = b(i,j+i1)*ceadd(j,ie,is)
            call dgemm('N','N',nlma,nstate,nlmb,1d0,wk,nlma,
     .        z(iofz,1),nhs,1d0,bz(1,1,ie),nlma)
          else
            call dgemm('N','N',nlma,nstate,nlmb,1d0,b(1,1+i1),nlma,
     .        z(iofz,1),nhs,1d0,bz(1,1,ie),nlma)
          endif
        endif
   10 continue

      end
C      subroutine xxhm0(ldot,lgrad,llink,nlma,nlbe,nll,
C     .  ob,obd,ogb,ogbd,obl,obld,ogbl,ogbld)
CC- Allocates memory for strx about one site, and its derivatives
CC  nlbe is dimensioned for all neighbors and energies connecting strx.
CC  nll is for linked basis.
CC  ldot (lgrad, llink) ne 0: allocate arrays for dot (gradient, linked)
C      implicit none
C      integer nlma,nlbe,nll,ob,obd,ogb,ogbd,obl,obld,ogbl,ogbld
C      logical ldot,lgrad,llink
C
C      call defrr(ob, nlma*nlbe)
C      if (ldot)  call defrr(obd, nlma*nlbe)
C      if (lgrad) call defrr(ogb, nlma*nlbe*3)
C      if (ldot .and. lgrad) call defrr(ogbd, nlma*nlbe*3)
C      if (llink) then
C        call defrr(obl, nlma*nll)
C        if (ldot)  call defrr(obld, nlma*nll)
C        if (lgrad) call defrr(ogbl, nlma*nll*3)
C        if (ldot .and. lgrad) call defrr(ogbld, nlma*nll*3)
C      endif
C      end
C
C      subroutine xxrdmt(ewt,nhs,nst,z,wk,dmat)
CC- Makes density matrix
C      implicit none
C      integer nhs,nst,i,j
C      double precision z(nhs,nst),wk(nst,nhs),ewt(nst),
C     .  dmat(nhs,nhs)
C
C      do  10  j = 1, nst
C      do  10  i = 1, nhs
C   10 wk(j,i) = ewt(j)*z(i,j)
C      call dgemm('N','N',nhs,nhs,nst,1d0,z,nhs,wk,nst,0d0,dmat,nhs)
C      end
