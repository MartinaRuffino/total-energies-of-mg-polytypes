      subroutine ropfxc(n1,n2,n3,nbas,ips,n0,nxi,lxi,
     .  exi,rsm,alat,plat,pos,jobg,ioff,rhoi,rho0,
     .  vxci,vxc0,lforce,dvxci,dvxc0,rhoq,exc,vxc,f)
C- Exchange-correlation energy and potential of smoothed density
C  obtained by FFT of rho(q).  jobg: passed to ropffq
C  24 Dec spin pol version.  (Not set up for rho holding gaussians)
      implicit none
      integer nbas,n1,n2,n3,n0,ips(1),lxi(n0,1),nxi(1),jobg,lforce,
     .  ioff(1)
      double precision exi(n0,1),rsm(1),pos(3,1),rhoi(1),f(3,nbas),
     .  rho0(1),alat,exc,vxc,totq,vxci(1),vxc0(1),dvxci(1),dvxc0(1),
     .  repnl(2),rmunl(2)
      double precision rhoq(n1+2,n2,n3,*)
*     complex*16 rhoq(n1,n2,n3,1)
      integer oqp,oexc,ovxc,orho,orho2,ovxc2,ogrho,oidx,oibp,obp,ovnl,
     .  oexcx,oexcc,ovxcx,ovxcc,ovxcx2,ovxcc2,
     .  i1,i2,k,isp,m,np,npf,n,id,id2,nq,nqt,nri,ipr,j1,j2,j3,ib,i,
     .  nbp,nn,nsp,lsp,lxcf,lxcg,i1mach
      double precision plat(3,3),vol,tripl,xx,vol0,ecut,exc2,vxc2,totq2
      logical lfrc
      real w(1)
      common /w/ w

C --- Setup ---
      call tcn('ropfxc')
C     call wkprnt(1)
      call getpr(ipr)
      nsp = lsp()+1
      np = n1*n2*n3
      n = n2*n3
      id = n1/2+1
      id2 = n1+2
      if (2*id /= id2) call rx('ropfxc: n1 must be even')
* reset id
*     id  = n1
*     id2 = n1
      npf = 2*id*n2*n3
      lfrc = lforce == 1
      nri = ioff(nbas+1)
      call dpzero(rhoq,npf*nsp)
      vol = alat**3*tripl(plat,plat(1,2),plat(1,3))
      vol0 = vol/np
      call xcecut(alat,plat,n1,n2,n3,ecut)
C      if (ipr >= 30) print 333, n1,n2,n3,np,vol0,ecut
      if (ipr >= 30)
     .  call awrit6(' n1..n3:%i %i %i, np=%i, vol/cell=%d, ecut=%d',' ',
     .              120,i1mach(2),n1,n2,n3,np,vol0,ecut)
      if (ipr >= 40) print 337, plat
  337 format('         plat:',9f7.3)
C 333 format(' ropfxc: n1..n3,np=',3i3,i6,'  vol/cell,ecut=',2f12.6)
C jobg=1 and nsp needs rho0(1+nvi) ...
      if (jobg == 1)
     .  call rx('ropfxc not set up for forces or spin pol and jobg=1')

C --- Make rho(q) for all points, one plane at a time ---
      call tcn('make rho(q)')
      nn = 2*n
      call defrr(oqp, nn*3)
      call defcc(orho,nn+1)
      call defrr(oidx,nn)
      call defrr(oibp,nn)
      call defrr(obp, nn*3)
      nqt = 0
      do  10  i1 = 1, id
        call gmeshx(alat,plat,0,ecut,n1,n2,n3,i1,nn,
     .    nq,nbp,w(oibp),w(obp),w(oidx),w(oqp))
        nqt = nqt+nq+nbp
        call ropffq(nbas,ips,n0,nxi,lxi,exi,nri,rsm,alat,pos,1/vol,
     .    nn,nq+nbp,w(oqp),rhoi,rho0,w(orho),lfrc,w,dvxci,1,jobg)
        call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .    w(oqp),i1,1,w(orho),0,w,rhoq)
        if (nsp == 2) then
          call ropffq(nbas,ips,n0,nxi,lxi,exi,nri,rsm,alat,pos,1/vol,
     .      nn,nq+nbp,w(oqp),rhoi(1+nri),rho0,w(orho),lfrc,w,dvxci,
     .      1,jobg)
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,1,w(orho),0,w,rhoq(1,1,1,2))
        endif
   10 continue
      call rlse(oqp)
      if (ipr >= 40) print 336, nqt, id*n, (100*nqt)/(id*n)
  336 format(' ropfxc: retained',i6,' points out of',i6,' (',i3,'%)')
      call tcx('make rho(q)')

C     call zprx('rhoq',rhoq,id,n2,n3)
*     call zprm3('rhoq',rhoq,n1,n2,n3)

C --- Make rho(r) by backward FT rho(q) ---
      call tcn('make rho(r)')
      totq =  vol*rhoq(1,1,1,1)
      totq2 = vol*rhoq(1,1,1,nsp)
C ... if gradient corrections, let vxcnlp do the FFT
      if (lxcg() == 0) then
        call rfft3(rhoq,n1,n2,n3,1,1)
C       call prm3x('rho(r) after rfft3',0,rhoq,n1+2,n2,n3,n1,n2,n3)
        if (nsp == 2) call rfft3(rhoq(1,1,1,2),n1,n2,n3,1,1)
*       call cfft3(rhoq,n1,n1,n2,n3,1,1)
*       if (nsp == 2) call cfft3(rhoq(1,1,1,2),n1,n1,n2,n3,1,1)
      else
        call defrr(ovnl,   npf*nsp)
        call vxcnlp(n1,n2,n3,nsp,alat,plat,vol0,ecut,rhoq,repnl,rmunl,
     .    w(ovnl))
      endif
      call tcx('make rho(r)')
C     call prm3('rho(r)',rhoq,n1+2,n1,n2,n3)
*     call zprm3('rho(r)',rhoq,n1,n2,n3)

C --- Make local exc, vxc, and overwrite rhoq with vxc ---
      call tcn('make local xc')
      exc = 0d0
      vxc = 0d0
      exc2 = 0d0
      vxc2 = 0d0
      call defrr(orho,n)
      call defrr(oexc,n)
      call defrr(ovxc,n)
      if (nsp == 2) then
        call defrr(orho2,n)
        call defrr(ovxc2,n)
      endif
      do  20  i1 = 1, n1
        call xrhfxc(id2,n,i1,0,w(orho),rhoq)
        if (nsp == 1) then
          call defrr(oexcx,n)
          call defrr(oexcc,n)
          call defrr(ovxcx,n)
          call defrr(ovxcc,n)
          call defrr(ovxcx2,n)
          call defrr(ovxcc2,n)
C         call ropevx(w(orho),w(oexc),w(ovxc),n)
          if (lxcf() > 2) then
            call evxcp(w(orho),w(orho),n,nsp,lxcf(),w(oexcx),w(oexcc),
     .        w(oexc),w(ovxcx),w(ovxcx2),w(ovxcc),w(ovxcc2),w(ovxc),
     .        w(ovxc))
          else
            call evxcv(w(orho),w(orho),n,nsp,lxcf(),
     .        w(oexc),w(oexcx),w(oexcc),w(ovxc),w(ovxcx),w(ovxcc))
          endif
          call rlse(oexcx)
        else
          call xrhfxc(n1+2,n,i1,0,w(orho2),rhoq(1,1,1,2))
          call defrr(oexcx,n)
          call defrr(oexcc,n)
          call defrr(ovxcx,n)
          call defrr(ovxcc,n)
          call defrr(ovxcx2,n)
          call defrr(ovxcc2,n)
          if (lxcf() > 2) then
            call evxcp(w(orho),w(orho2),n,nsp,lxcf(),w(oexcx),w(oexcc),
     .        w(oexc),w(ovxcx),w(ovxcx2),w(ovxcc),w(ovxcc2),w(ovxc),
     .        w(ovxc2))
          else
            call dpadd(w(orho2),w(orho),1,n,1d0)
            call evxcv(w(orho2),w(orho),n,2,lxcf(),
     .        w(oexc),w(oexcx),w(oexcc),
     .        w(ovxc),w(ovxcx),w(ovxcc))
            call dpadd(w(orho2),w(orho),1,n,-1d0)
            call dpadd(w(orho),w(orho2),1,n,1d0)
            call evxcv(w(orho),w(orho2),n,2,lxcf(),
     .        w(oexc),w(oexcx),w(oexcc),
     .        w(ovxc2),w(ovxcx),w(ovxcc))
            call dpadd(w(orho),w(orho2),1,n,-1d0)
          endif
          call rlse(oexcx)
        endif
        call dpdot(w(orho),w(oexc),n,xx)
        exc = exc + xx*vol0
        if (nsp == 2) then
          call dpdot(w(orho2),w(oexc),n,xx)
          exc2 = exc2 + xx*vol0
        endif
        call dpdot(w(orho),w(ovxc),n,xx)
        vxc = vxc + xx*vol0
        if (nsp == 2) then
          call dpdot(w(orho2),w(ovxc2),n,xx)
          vxc2 = vxc2 + xx*vol0
        endif
        call xrhfxc(id2,n,i1,1,w(ovxc),rhoq)
        if (nsp == 2) call xrhfxc(id2,n,i1,1,w(ovxc2),rhoq(1,1,1,2))
   20 continue
      call rlse(orho)
      call tcx('make local xc')

C --- Combine terms ---
      if (nsp == 1) then
        if (ipr >= 20) print 334, totq,exc,vxc
  334   format(' ropfxc: q, local exc,vxc=', 3f12.6)
        if (lxcg() /= 0) then
          call dpadd(rhoq,w(ovnl),1,npf*nsp,1d0)
          exc = exc+repnl(1)
          vxc = vxc+rmunl(1)
          if (ipr >= 20) print 332, repnl(1),rmunl(1),exc,vxc
  332     format('         nonlocal exc,vxc=', 12x,2f12.6/
     .           '   total smoothed exc,vxc=', 12x,2f12.6)
        endif
      elseif (nsp == 2) then
        if (ipr >= 20) then
          print 335, totq,exc,vxc
          print 335, totq2,exc2,vxc2
        endif
        totq = totq+totq2
        exc = exc+exc2
        vxc = vxc+vxc2
        print 335, totq,exc,vxc
        if (lxcg() /= 0) call rx('ropfxc: not ready for spol and grd')
      endif
  335 format(' ropfxc: q,exc,vxc=', 3f12.6)

C --- Forward FT V(r) to make V(q) ---
C     call prm3('V(r)',rhoq,n1+2,n1,n2,n3)
*     call zprm3('V(r)',rhoq,n1,n2,n3)
      call rfft3(rhoq,n1,n2,n3,1,-1)
*     call cfft3(rhoq,n1,n1,n2,n3,1,-1)
      if (nsp == 2) then
        call rfft3(rhoq(1,1,1,2),n1,n2,n3,1,-1)
*       call cfft3(rhoq(1,1,1,2),n1,n1,n2,n3,1,-1)
      endif
C     call zprx('Vq',rhoq,id,n2,n3)
*     call zprm3('Vq',rhoq,n1,n2,n3)

C --- Make zetas in reciprocal space ---
      call defrr(oqp, nn*3)
      call defcc(orho,nn+1)
      ogrho = 1
      if (lfrc) call defcc(ogrho,3*(nn+1))
      call defrr(oidx,nn)
      do  30  i1 = 1, id
        call gmeshx(alat,plat,0,ecut,n1,n2,n3,i1,nn,
     .    nq,nbp,w(oibp),w(obp),w(oidx),w(oqp))
        call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .    w(oqp),i1,0,w(orho),lforce,w(ogrho),rhoq)
        xx = 2
        if (i1 == 1 .or. i1 == id) xx = 1
*       xx = 1
        call ropffq(nbas,ips,n0,nxi,lxi,exi,nri,rsm,alat,pos,xx,
     .    nn,nq+nbp,w(oqp),vxci,vxc0,w(orho),lfrc,w(ogrho),dvxci,2,jobg)
        if (nsp == 2) then
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,0,w(orho),lforce,w(ogrho),rhoq(1,1,1,2))
          call ropffq(nbas,ips,n0,nxi,lxi,exi,nri,rsm,alat,pos,xx,
     .      nn,nq+nbp,w(oqp),vxci(1+nri),vxc0,w(orho),lfrc,w(ogrho),
     .      dvxci(1+nri*3),2,jobg)
        endif
   30 continue
      call rlse(oqp)

C      call pshpr(41)
C      call prxcp(nri,nxi,lxi,n0,nbas,ips,rhoi,vxci,dvxci)
C      call prxcp(nri,nxi,lxi,n0,nbas,ips,rhoi(1+nri),
C     .  vxci(1+nri),dvxci(1+3*nri))

C --- Forces from rho * grad(zeta) ---
      if (lfrc) then
        do  40  isp = 1, nsp
        j1 = 0
        j2 = nri
        j3 = 2*nri
        k = nri*(isp-1)
C       if (ipr >= 30) write(6,289)
        do  40  ib = 1, nbas
          i1 = ioff(ib) + 1
          i2 = ioff(ib+1)
          if (isp == 1) then
            f(1,ib) = 0d0
            f(2,ib) = 0d0
            f(3,ib) = 0d0
          endif
          do  41  i = i1, i2
          f(1,ib) = f(1,ib) - rhoi(i+k)*dvxci(i+j1+k*3)
          f(2,ib) = f(2,ib) - rhoi(i+k)*dvxci(i+j2+k*3)
          f(3,ib) = f(3,ib) - rhoi(i+k)*dvxci(i+j3+k*3)
   41     continue
          if (ipr >= 30 .and. isp == nsp)
     .      write(6,288) ib,(f(m,ib),m=1,3)
  288     format(i6,3f13.6)
C 289     format(/'    ib     xc-force from smooth density')
   40   continue
      endif

C     call ppot3(nri,nxi,lxi,0,n0,nbas,ips,vxci,vxci,
C    .  rhoi,rho0,rho0,rho0,0)
      if (ipr >= 50) call tc('Exit ropfxc ...')
      call tcx('ropfxc')
      end
      subroutine zrhfxc(id,n,idx,nn,nq,nbp,ibp,G,i1,isw,rho,lforce,
     .  grho,rhoq)
C- Copy rhoq from/to plane to/from 3D array; average boundary points
C  isw=0: copy rhoq(i1,i) to rho(idx(i))  1: copy rho to rhoq
C  lforce=0: do nothing extra
C  lforce=1  isw=0: make grho(q) = G . rhoq  isw=1: copy grho to rhoq
C  lforce=2  isw=0: make grho(q) = G*G rhoq  isw=1: copy grho to rhoq
      implicit none
      integer id,n,i1,isw,idx(1),nn,nbp,ibp(1),nq,lforce
      double precision G(nn,3)
      complex*16 rho(0:nn), grho(0:nn,3), rhoq(id,n,1), cxx
      integer i,ib,ix,nqbp

C --- Copy 3D potential rho into 2D array; extra points for boundary ---
      if (isw == 0) then
        do  10  i = 1, n
   10   rho(idx(i)) = rhoq(i1,i,1)
C ...   Use average pot for boundary G vectors
        do  14  ib = 1, nbp
          i = ibp(ib)
          cxx = rho(i)/2
          rho(i) = cxx
          rho(ib+nq) = cxx
   14   continue
C ...   For integrals of potential gradient
        if (lforce == 0) return
        nqbp = nq+nbp
        if (lforce == 1) then
          do  16  ix = 1, 3
          call dpzero(grho(0,ix),2*(nqbp+1))
          do  16  i = 1, nqbp
   16     grho(i,ix) = rho(i)*dcmplx(0d0,G(i,ix))
        elseif (lforce == 2) then
          call dpzero(grho(0,1),2*(nqbp+1))
          do  18  i = 1, nqbp
   18     grho(i,1) = -rho(i)*(G(i,1)**2 + G(i,2)**2 + G(i,3)**2)
        endif
C --- Copy 2D rho array into 3D; avg extra points from boundary ---
      else
C ...   Use average rho for boundary G vectors
        do  24  ib = 1, nbp
          i = ibp(ib)
          rho(i) = (rho(i) + rho(ib+nq))/2
C         uncomment to get second set of points only
C         rho(i) = rho(ib+nq)
C         uncomment to copy second set of G vectors into first
C         G(i,3) = G(ib+nq,3)
   24   continue
C ...   Do the copy
        if (lforce == 0) then
C          ix = 2
C          print *, 'moving iG to rho ..,ix=',ix
C          do  99  i = 1, n
C   99     rho(i) = dcmplx(0d0,G(i,ix))
          do  20  i = 1, n
   20     rhoq(i1,i,1) = rho(idx(i))
        elseif (lforce == 1) then
          do  26  ix = 1, 3
          do  26  i = 1, n
   26     rhoq(i1,i,ix) = grho(idx(i),ix)
        elseif (lforce == 2) then
          do  28  i = 1, n
   28     rhoq(i1,i,1) = grho(idx(i),1)
        endif
      endif
      end
      subroutine xrhfxc(n1p2,n,i1,isw,rho,rhoq)
C- Copy rhoq(i1,*) to/from rho(*)
      implicit none
      integer n1p2,n,i1,isw
      double precision rho(n), rhoq(n1p2,n)
*     double precision rho(n)
*     double complex rhoq(n1p2,n)
      integer i

      if (isw == 0) then
        do  10  i = 1, n
   10   rho(i) = rhoq(i1,i)
      else
        do  20  i = 1, n
   20   rhoq(i1,i) = rho(i)
      endif
      end
      subroutine prxcp(nri,nxi,lxi,n0,nbas,ips,rhoi,zeta,gzet)
C- Prints out exchange zeta and gradient
C  For printout, hankels renormalized by dividing by (2l+1)!!
      implicit real*8 (a-h,p-z), integer (o)
      dimension zeta(1),ips(1),nxi(1),lxi(n0,1),rhoi(1),
     .  gzet(nri,*)
      sum1=0d0
      sum2=0d0
      sum3=0d0
      sum4=0d0
      do 1 i=1,nri
      sum1=sum1+rhoi(i)*zeta(i)
      sum2=sum2+rhoi(i)*gzet(i,1)
      sum3=sum3+rhoi(i)*gzet(i,2)
    1 sum4=sum4+rhoi(i)*gzet(i,3)
      if(iprint() >= 30) write(6,331) sum1,sum2,sum3,sum4
  331 format(/' prxcp: rhoi*(zetxc,gzet):',4f11.6)
      if (iprint() < 41) return
      write(6,334)
  334 format(/'    n ib ie ilm    Rhoi      Zetxc',17x,'Gzetxc')
      i=0
      do 10 ib=1,nbas
      is=ips(ib)
      do 10 ie=1,nxi(is)
      df=1d0
      ilm=0
      do 10 l=0,lxi(ie,is)
      df=df*(2*l+1)
      f=1d0/df
C     f = 1
      do 10 m=1,2*l+1
      i=i+1
      ilm=ilm+1
      top=dmax1(dabs(rhoi(i)),dabs(zeta(i)),dabs(gzet(i,1)),
     .  dabs(gzet(i,2)),dabs(gzet(i,3)))
      if (top > 1d-5) print 333,i,ib,ie,ilm,
     .  rhoi(i)/f,f*zeta(i),f*gzet(i,1),f*gzet(i,2),f*gzet(i,3)
  333 format(i5,3i3,6f11.6)
   10 continue
      end
