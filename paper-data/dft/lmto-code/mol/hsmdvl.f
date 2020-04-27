      subroutine hsmdvl(nbas,isp,nphi,lphi,ips,n0,nhs,nhl,nstate,
     .  nla,ioffa,ioffb,ioffl,
     .  el,nel,evl,z,hvecs,svecs,hi,si,ceadd,
     .  gh,gj,puu,pus,pss,ioffp,
     .  ixp,dxp,pos,alat,cg,jcg,indxcg,cy,wgt,ewt,f3)
C- Change in sumeV from virtual displacements (real h, linked basis)
C  This version calls astrgl, decomposes force explicitly by pairs
      implicit none
      integer n0,nbas,isp,nel,nhs,nla,ips(1),nphi(1),lphi(n0,1),
     .  indxcg(1),jcg(1),nstate,ixp(1),
     .  ioffp(1),ioffb(nbas,nel),ioffl(nbas+1),ioffa(1),nhl
      double precision hvecs(nla,4,10),svecs(nla,4,6),hi(6),si(6),el(1),
     .  gh(nla,2,nel),gj(nla,2,nel),puu(1),pus(1),pss(1),evl(nhs,2),
     .  dxp(1),
     .  cg(1),cy(1),pos(3,1),alat,z(nhs,nstate),f3(3,nbas),
     .  ewt(nhs,isp),wgt,ceadd(25,5,1)
      double precision xxh,xxs,xx2,xx3,xxt
      integer ie,iv,ivs,je,ipr,ib,ib1,ib2,is,nlm1,nlm2,nlma,nx,k0,io,
     .  owk3a,owk3b,owk2a,owk2b,owk3al,owk3bl,owk2al,owk2bl,
     .  obz,ogbz,ogdz,obzl,ogbzl,ogdzl
      integer nsp,lsp,ivl,jvl,nvl,ivls,jvls,nvls,il,ndimi,ndimj,
     .  opkk,opkj,opjk,opjj,opkk2,opkj2,opjk2,opjj2,
     .  opkk3,opkj3,opjk3,opjj3,opkkt,opkjt,opjkt,opjjt,oadd1,oadd2
      parameter (nx=49)
      double precision pkk(nx*nx),pkj(nx*nx),pjk(nx*nx),pjj(nx*nx)
      real w(1)
      common /w/ w

C For MPI ...
      integer, dimension(:), allocatable :: bproc
      integer mpipid,lgunit,procid,master,numprocs,length,ierr
      logical mlog,cmdopt,MPI
      character outs*80, shortname*10, datim*26
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      call mpiprn(shortname,length)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      call tcn('hsmdvl')
      call getpr(ipr)
      nsp = lsp()+1
      call defrr(opkk, nx*nx)
      call defrr(opkj, nx*nx)
      call defrr(opjk, nx*nx)
      call defrr(opjj, nx*nx)
      il = 0
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
        call defrr(oadd1,25*nbas)
        call defrr(oadd2,25*nbas)
        il = nel+1
        nvl = ((nel+1)*(nel+2))/2
        nvls= nsp*(nvl-1)+isp
      endif

C --- For virtual displacement of site ib ---
      if (MPI) then
        allocate (bproc(0:numprocs), stat=ierr)
        call dstrbp(nbas,numprocs,1,bproc(0))
        ib1 = bproc(procid)
        ib2 = bproc(procid+1)-1
      else
        ib1 = 1
        ib2 = nbas
      endif
      do  10  ib = ib1, ib2
        if (MPI .and. mlog) then
        if (ib == bproc(procid)) then
          call gettime(datim)
          call awrit4(' hsmdvl '//datim//' Process %i of %i on '
     .        //shortname(1:length)//
     .        ' starting atoms %i to %i',' ',256,lgunit(3),
     .         procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        endif
        is = ips(ib)
        nlma = ioffa(ib+1)-ioffa(ib)

C ...   Work arrays allocated
        if (nlma > nx) call rx('hsmdvl: nlma gt nx')
        call defrr(owk3a,    nlma*nstate)
        call defrr(owk3b,    nlma*nstate)
        call defrr(owk2a,    nlma*nstate)
        call defrr(owk2b,    nlma*nstate)
        call defrr(ogbz,     nlma*nstate*nbas*3*nel)
        call defrr(ogdz,     nlma*nstate*3*nel)
        call defrr(obz,      nlma*nstate*nel)
        call dpzero(w(ogbz), nlma*nstate*nbas*3*nel)
        call dpzero(w(ogdz), nlma*nstate*3*nel)
        call dpzero(w(obz),  nlma*nstate*nel)
        obzl = 1; ogbzl = 1; ogdzl = 1
        if (nhl > 0) then
          call defrr(ogbzl,    nlma*nstate*nbas*3*nel)
          call defrr(ogdzl,    nlma*nstate*3*nel)
          call defrr(obzl,     nlma*nstate*nel)
          call dpzero(w(ogbzl),nlma*nstate*nbas*3*nel)
          call dpzero(w(ogdzl),nlma*nstate*3*nel)
          call dpzero(w(obzl), nlma*nstate*nel)
          call defrr(owk3al,   nlma*nstate)
          call defrr(owk3bl,   nlma*nstate)
          call defrr(owk2al,   nlma*nstate)
          call defrr(owk2bl,   nlma*nstate)
        endif

C --- b * z, grad b * z, grad bdot * z, all sites with tails in ib ---
        call astrgl(nphi,lphi,nel,el,nlma,ib,pos,n0,nbas,ips,alat,cg,
     .    jcg,indxcg,cy,nhs,nhl,nstate,ioffb,ioffl,ixp,dxp,z,
     .    ceadd,w(obz),w(ogbz),w(ogdz),w(obzl),w(ogbzl),w(ogdzl))

C --- Loop over energy pairs ---
      iv = 0
      do  20  ie = 1, nel
C        if (nhl > 0) call hmcadd(ips,nbas,ioffb,ie,ceadd,w(oadd1))
        ivl = (ie*(il+il-ie+1))/2
        ivls= nsp*(ivl-1)+isp
        ndimi = ioffb(nbas+1,ie)-ioffb(1,ie)
        do  30  je = ie, nel
        jvl = (je*(il+il-je+1))/2
        jvls= nsp*(jvl-1)+isp
        ndimj = ioffb(nbas+1,je)-ioffb(1,je)
C       if (nhl > 0) call hmcadd(ips,nbas,ioffb,je,ceadd,w(oadd2))
        iv = iv+1
        ivs = nsp*(iv-1)+isp
        nlm1 = (lphi(ie,is)+1)**2
        nlm2 = (lphi(je,is)+1)**2

C ---   perturbation matrices for (ib,ie,je) ---
        io = nsp*ioffp(ib)+1
        if (isp == 2) io = io + nlma**2
        k0 = ioffa(ib)+1
        xxh = 0d0
        if (ie /= je) xxh = hi(iv)/(el(je)-el(ie))
        call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,ivs),nla,xxh,
     .    gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,je),gh(k0,2,je),
     .    gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,je),gj(k0,2,je),
     .    nlm1,nlm2,nlma,pkk,pkj,pjk,pjj)
C       print *, 'chk hvec',ib,ie,je,hvecs(1,1,ivs)

C ---   Forces for this (ie,je,ib) block ---
        xxs = 0
        if (ie /= je) xxs = si(iv)/(el(je)-el(ie))
C ...   Case no linked basis
        if (nhl == 0) then
          if (ndimi > 0 .and. ndimj > 0) then
            call xxmde2(nbas,nhs,nla,nstate,ie,je,nel,ib,nlm1,nlm2,nlma,
     .        ixp,ioffb,ioffa,hi(iv),pkj,pjk,pjj,si(iv),svecs(k0,1,ivs),
     .        xxs,evl(1,isp),w(owk3a),w(owk3b),w(owk2a),w(owk2b),
     .        z,w(obz),w(ogbz),w(ogdz),wgt,ewt(1,isp),f3)
          endif
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
C         if (ndimi > 0 .and. ndimj > 0)
          call xxmdvl(nbas,nhs,nla,nstate,ie,je,nel,ib,nlm1,nlm2,nlma,
     .      ixp,ioffb,ioffa,hi(iv),hi(ivl),si(iv),si(nvl),
     .      svecs(k0,1,ivs),xxs,svecs(k0,1,ivls),xx2,svecs(k0,1,nvls),
     .      xx3,svecs(k0,1,jvls),xxt,ceadd(1,ie,is),ceadd(1,je,is),
     .      pkj,pjk,pjj,w(opkj2),w(opjk2),w(opjj2),
     .      w(opkj3),w(opjk3),w(opjj3),
     .      w(opkjt),w(opjkt),w(opjjt),
     .      evl(1,isp),w(owk3a),w(owk3b),w(owk2a),w(owk2b),
     .      w(owk3al),w(owk3bl),w(owk2al),w(owk2bl),
     .      z,w(obz),w(ogbz),w(ogdz),w(obzl),w(ogbzl),w(ogdzl),
     .      wgt,ewt(1,isp),f3)
        endif
   30   continue
        if (nhl /= 0) iv = iv+1
   20 continue
      call rlse(owk3a)
   10 continue
      call rlse(opkk)
      if (MPI) then
        call mpibc2(f3,3*nbas,4,3,mlog,'hsmdvl','f3')
        deallocate (bproc, stat=ierr)
      endif
C ---   Printout -------------------
      if (ipr >= 30) then
        write(6,289)
        do  ib = 1, nbas
          write(6,288) ib,f3(1,ib),f3(2,ib),f3(3,ib)
        enddo
      endif
  288 format(i6,3f13.6)
  289 format(/'    ib     eigval-force from asa hamiltonian')

      call tcx('hsmdvl')
      end
      subroutine xxmde2(nbas,nhs,nla,nstate,ie,je,nel,ib,nlm1,nlm2,nlma,
     .  ixp,ioffb,ioffa,hi,pkj,pjk,pjj,si,svec,xxx,evl,
     .  wk3a,wk3b,wk2a,wk2b,z,bz,gbz,gdz,wgt,ewt,f3)
C- Change in eV sum from displacing ib, one pair of energies
      implicit none
      integer nbas,nhs,ib,nstate,nla,nel,nlm1,nlm2,nlma,ie,je,
     .  ioffb(nbas,1),ioffa(1),ixp(1)
      double precision pkj(nlm1,nlma),xxx,sum,sumi,wgt,ewt(1),evl(1),
     .  pjk(nlma,nlm2),pjj(nlma,nlma),hi,si,svec(nla,4),f3(3,nbas)
      double precision bz(nlma,nstate,nel),z(nhs,1),
     .  gbz(nlma,nstate,nbas,3,nel),wk2a(nlma,nstate),wk2b(nstate,nlma),
     .  gdz(nlma,nstate,3,nel),wk3a(nlma,nstate),wk3b(nstate,nlma)
      integer ix,jb,ioffz,nlmk,is,jlm,jbot,jtop

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
   62   wk3b(is,jlm) = wk3a(jlm,is)
      else
        call dgemm('T','N',nstate,nlma,nlma,1d0,
     .    bz(1,1,ie),nlma,pjj,nlma,0d0,wk3b,nstate)
        do  64  is = 1, nstate
        do  64  jlm = 1, nlma
   64   wk3b(is,jlm) = wk3b(is,jlm) -
     .                 bz(jlm,is,ie)*svec(jlm,4)*evl(is)
      endif

C --- Contract pkj,pjk blocks with z ---
      if (nlm1 > 0)
     .  call dgemm('T','N',nstate,nlma,nlm1,1d0,
     .  z(1+ioffb(ib,ie),1),nhs,pkj,nlm1,0d0,wk2b,nstate)
      do  66  is = 1, nstate
      do  66  jlm = 1, nlm1
   66 wk2b(is,jlm) = wk2b(is,jlm) -
     .    z(jlm+ioffb(ib,ie),is)*(svec(jlm,2)+xxx)*evl(is)
      if (nlm2 > 0)
     .  call dgemm('N','N',nlma,nstate,nlm2,1d0,
     .  pjk,nlma,z(1+ioffb(ib,je),1),nhs,0d0,wk2a,nlma)
      do  68  is = 1, nstate
      do  68  jlm = 1, nlm2
   68 wk2a(jlm,is) = wk2a(jlm,is) -
     .    z(jlm+ioffb(ib,je),is)*(svec(jlm,3)-xxx)*evl(is)

      do  10  ix = 1, 3

C        print *, 'WW skip to jk,kj'
C        goto  999

C ---   z+ grad b+(ie) pjj b(je) z ---
        do  40  jb = 1, nbas
          sum = 0
          do  44  is = 1, nstate
          sumi = 0
          do  45  jlm = 1, nlma
   45     sumi = sumi + gbz(jlm,is,jb,ix,ie)*wk3a(jlm,is)
   44     sum = sum + sumi*ewt(is)
          sum = sum*2
          f3(ix,ib) = f3(ix,ib) + sum*wgt
          f3(ix,jb) = f3(ix,jb) - sum*wgt
C         print 345, ix,'3C   ib,jb,sum=',ib,jb,sum,f3(1,2)
   40   continue
C  345   format(1x,i1,1x,a,2i4,2g24.16)

C ---   z+ b+(ie) pjj grad b(je) z  ---
       if (ie /= je) then
          do  50  jb = 1, nbas
            sum = 0
            do  54  is = 1, nstate
            sumi = 0
            do  55  jlm = 1, nlma
   55       sumi = sumi + gbz(jlm,is,jb,ix,je)*wk3b(is,jlm)
   54       sum = sum + sumi*ewt(is)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'3C-2 ib,jb,sum=',ib,jb,sum,f3(1,2)
   50     continue
        endif

C        print *, 'WW do jj only'
C        goto 10

C ---   pkj * grad b(je), (2C terms with head in ib) --
  999   continue
        if (nlm1 > 0) then
          jbot = 1
          if (ie == je) jbot = ib
          do  20  jb = jbot, nbas
            if (jb == ib .and. ixp(1) == 0) goto 20
            sum = 0
            do  24  is = 1, nstate
            sumi = 0
            do  25  jlm = 1, nlma
   25       sumi = sumi + wk2b(is,jlm)*gbz(jlm,is,jb,ix,je)
   24       sum = sum + sumi*ewt(is)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'2C   ib,jb,sum=',ib,jb,sum,f3(1,2)
   20     continue
        endif

C ---   grad b+(je) pjk (2C terms with tail in ib)--
        if (nlm2 > 0) then
          jtop = nbas
          if (ie == je) jtop = ib
          do  30  jb = 1, jtop
            if (ib == jb .and. ixp(1) == 0) goto 30
            sum = 0
            do  34  is = 1, nstate
            sumi = 0
            do  35  jlm = 1, nlma
   35       sumi = sumi + wk2a(jlm,is)*gbz(jlm,is,jb,ix,ie)
   34       sum = sum + sumi*ewt(is)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'2Cx  ib,jb,sum=',ib,jb,sum,f3(1,2)
   30     continue
        endif

C        print *, 'WW skip bdot'
C        goto 10

C ---   grad bdot ---
        if (nlm1 > 0 .and. ie == je) then
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
C         print 345, ix,'dot ib,allj,sum=',ib,0,sum,f3(1,2)
        endif

   10 continue

      end
      subroutine xxmdvl(nbas,nhs,nla,nstate,ie,je,nel,ib,nlm1,nlm2,nlma,
     .  ixp,ioffb,ioffa,hi,hil,si,sil,
     .  svec,xxs,svc2,xx2,svc3,xx3,svct,xxt,addi,addj,
     .  pkj,pjk,pjj,pkj2,pjk2,pjj2,pkj3,pjk3,pjj3,pkjt,pjkt,pjjt,
     .  evl,
     .  wk3a,wk3b,wk2a,wk2b,wk3al,wk3bl,wk2al,wk2bl,
     .  z,bz,gbz,gdz,bzl,gbzl,gdzl,wgt,ewt,f3)
C- Change in sumeV from ie,je and displacing ib, linked basis
      implicit none
      integer nbas,nhs,ib,nstate,nla,nel,nlm1,nlm2,nlma,ie,je,
     .  ioffb(nbas,1),ioffa(1),ixp(1)
      double precision sum,sumi,wgt,ewt(1),evl(1),
     .  pkj (nlm1,nlma),pjk (nlma,nlm2),pjj (nlma,nlma),hi,hil,si,sil,
     .  pkj2(nlm1,nlma),pjk2(nlma,nlm2),pjj2(nlma,nlma),
     .  pkjt(nlm2,nlma),pjkt(nlma,nlm2),pjjt(nlma,nlma),
     .  pkj3(nlm1,nlma),pjk3(nlma,nlm2),pjj3(nlma,nlma),
     .  svec(nla,4),xxs,svc2(nla,4),xx2,svc3(nla,4),xx3,svct(nla,4),xxt,
     .  f3(3,nbas),addi(nlm1),addj(nlm2)
      double precision bz(nlma,nstate,nel),bzl(nlma,nstate,nel),
     .  gbz(nlma,nstate,nbas,3,nel),gbzl(nlma,nstate,nbas,3,nel),
     .  z(nhs,1),gdz(nlma,nstate,3,nel),gdzl(nlma,nstate,3,nel),
     .  wk2a(nlma,nstate),wk2b(nlma,nstate),
     .  wk3a(nlma,nstate),wk3b(nlma,nstate),
     .  wk2al(nlma,nstate),wk2bl(nlma,nstate),
     .  wk3al(nlma,nstate),wk3bl(nlma,nstate),ddot
      integer ix,jb,ioffz,nlmk,is,ilm,jlm,jbot,jtop,ila

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

C --- Contract pkj,pjk blocks with z ---
      do  21  ila = 1, nlma
      do  21  ilm = 1, nlm1
   21 pkj(ilm,ila) = pkj(ilm,ila) + addi(ilm)*pjkt(ila,ilm)
      if (nlm1 > 0)
     . call dgemm('T','N',nlma,nstate,nlm1,1d0,
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
      if (nlm1 > 0)
     .  call dgemm('T','N',nlma,nstate,nlm1,1d0,
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
      if (nlm2 > 0)
     .  call dgemm('N','N',nlma,nstate,nlm2,1d0,
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
      if (nlm2 > 0)
     .  call dgemm('N','N',nlma,nstate,nlm2,1d0,
     .  pjk3,nlma,z(1+ioffb(ib,je),1),nhs,0d0,wk2al,nlma)
      do  69  is = 1, nstate
      do 169  jlm = nlm2+1, nlma
  169 wk2al(jlm,is) = wk2al(jlm,is) * ewt(is)
      do  69  jlm = 1, nlm2
   69 wk2al(jlm,is) = (wk2al(jlm,is) - z(jlm+ioffb(ib,je),is)*evl(is)*
     .    (svct(jlm,2)+xxt + addj(jlm)*(svc3(jlm,3)-xx3))) * ewt(is)

      do  10  ix = 1, 3

C      print *, 'WW skip to jk,kj'
C      goto  999

C ---   z+ grad b+(ie) (pjj b(je) + pjj2 b(le)) z ---
C and   z+ grad b+(le) (pjjt b(je) + pjj3 b(le)) z
        do  40  jb = 1, nbas
          sum = ddot(nlma*nstate,wk3a, 1,gbz(1,1,jb,ix,ie), 1) +
     .          ddot(nlma*nstate,wk3al,1,gbzl(1,1,jb,ix,ie),1)
          sum = sum*2
          f3(ix,ib) = f3(ix,ib) + sum*wgt
          f3(ix,jb) = f3(ix,jb) - sum*wgt
C         print 345, ix,'3C   ib,jb,sum=',ib,jb,sum,f3(1,2)
   40   continue
C 345   format(1x,i1,1x,a,2i4,2g24.16)

C ---   z+ (b+(ie) pjj  +  b+(le) pjjt)  grad b(je) z  ---
C and   z+ (b+(ie) pjj2 +  b+(le) pjj3)  grad b(le) z
        if (ie /= je) then
          do  50  jb = 1, nbas
          sum = ddot(nlma*nstate,wk3b, 1,gbz(1,1,jb,ix,je), 1) +
     .          ddot(nlma*nstate,wk3bl,1,gbzl(1,1,jb,ix,je),1)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'3C-2 ib,jb,sum=',ib,jb,sum,f3(1,2)
   50     continue
        endif

C        print *, 'WW do jj only'
C        goto 10
C  999   continue

C ---   2C terms <ib,ie+le | grad jb,je> ---
        if (nlm1 > 0) then
          jbot = 1
          if (ie == je) jbot = ib
C         call dgemv('T',nstate*nlma,nbas-jbot+1,1d0,
C    .      gbz(1,1,jbot,ix,je),nstate*nlma,wk2b,1,0d0,sumjb(jbot),1)
          do  20  jb = jbot, nbas
            if (jb == ib .and. ixp(1) == 0) goto 20
C           sum = sumjb(jb)*2
            sum = ddot(nstate*nlma,wk2b,1,gbz(1,1,jb,ix,je),1)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'2C   ib,jb,sum=',ib,jb,sum,f3(1,2)
   20     continue
C ...     Linked strux terms <ib,ie+le | grad jb,le>
          do  26  jb = jbot, nbas
            if (jb == ib .and. ixp(1) == 0) goto 26
            sum = ddot(nstate*nlma,wk2bl,1,gbzl(1,1,jb,ix,je),1)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'2Cl  ib,jb,sum=',ib,jb,sum,f3(1,2)
   26     continue
        endif

C ---   2C terms <ib,je+le | grad jb,ie> ---
        if (nlm2 > 0) then
          jtop = nbas
          if (ie == je) jtop = ib
          do  30  jb = 1, jtop
            if (ib == jb .and. ixp(1) == 0) goto 30
            sum = ddot(nstate*nlma,wk2a,1,gbz(1,1,jb,ix,ie),1)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'2Cx  ib,jb,sum=',ib,jb,sum,f3(1,2)
   30     continue
C ...     Linked strux terms <ib,je+le | grad jb,le>)
          do  36  jb = 1, jtop
            if (ib == jb .and. ixp(1) == 0) goto 36
            sum = ddot(nstate*nlma,wk2al,1,gbzl(1,1,jb,ix,ie),1)
            sum = sum*2
            f3(ix,ib) = f3(ix,ib) + sum*wgt
            f3(ix,jb) = f3(ix,jb) - sum*wgt
C           print 345, ix,'2Clx ib,jb,sum=',ib,jb,sum,f3(1,2)
   36     continue
        endif

C        print *, 'WW skip bdot'
C        goto 10

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
C         print 345, ix,'dot ib,ie,sum=',ib,ie,sum,f3(1,2)
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
