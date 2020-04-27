      subroutine sprho2(nbas,pos,ng,g,ag,nclas,ipclas,ips,lmxa,lmxl,
     .   quu,qus,qss,ioffp,z,rmt,nr,a,pnu,pnuz,ipan,lsc,n0,ov0,
     .   orho,qsph,cg,indxcg,jcg,cy)
c- assemble density in sphere; two-panel
C 6 Jan 95 spin polarized (MvS)
      parameter( ncmx=20, nlmx=49, npwr=5 )
      implicit real*8 (a-h,p-z), integer (o)
      integer intopt,nglob
      dimension ipbas(ncmx),posc(3,ncmx),g(9,1),ag(3,1),ipclas(1),
     .   ips(1),lmxl(1),pos(3,1),quu(1),qus(1),qss(1),ioffp(1),
     .   orho(1),lmxa(1),z(1),nr(1),a(1),cg(1),jcg(1),indxcg(1),cy(1),
     .   pnu(n0,2,1),ov0(1),rmt(1),pnuz(n0,2,1),smec(2),smtc(2),qcor(2)
      real w(1)
      common /w/ w
      call tcn('sprho2')
      call getpr(ipr)
      if(ipr >= 20) write(6,220) nclas,ng
  220 format(/' sprho2:  number of classes=',i3,'     group ops=',i3)
      srfpi=dsqrt(16d0*datan(1d0))
      qsph=0d0
      nsp=lsp()+1

c ------ start loop over classes: make list of atoms in class -----
      do 50 ic=1,nclas
      nrc=0
      do 1 ib=1,nbas
      if(ipclas(ib) == ic) then
        nrc=nrc+1
        if(nrc > ncmx) call rx('sprho2 incr ncmx')
        ipbas(nrc)=ib
        do 2 m=1,3
    2   posc(m,nrc)=pos(m,ib)
      endif
    1 continue
c ------ representative atom is first in the class -----
      ia=1
      ib=ipbas(ia)
      is=ips(ib)
      nlml=(lmxl(is)+1)**2
      lmaxa=lmxa(is)
      nr1=nr(is)
      if(ipr >= 35) write(6,910) ic,is,nrc,ib
  910 format(/' symmetry class',i3,'    species=',i3,'    atoms:',i3,
     .   '    first atom: ib=',i3)
      call defrr(orofi,  nr1)
      call defrr(orwgt,  nr1)
      call radmsh(rmt(is),a(is),nr1,w(orofi))
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rmt(is),a(is),nr1,w(orwgt))
      call defrr(orhol,  nr1*nlml)
      call defrr(orhosp, nr1*nlml*nsp)
      call defrr(orhoc,  nr1*nsp)

c ------ make core density --------------
      if(ipan == 1) then
C --- call getcore from lmf, single panel only
        call defrr(opz0,  -10*nsp)
        qcor(1) = 0
        qcor(2) = 0
        call getcor(0,z(is),a(is),pnu(1,1,ib),w(opz0),nr1,lmaxa,
     .    w(orofi),w(ov0(ib)),0,-1,qcor,smec,smtc,w(orhoc),ncore,
     .    0d0,0d0)
        call rlse(opz0)
      endif

c ------ setup for assembling symmetrized density ---------
      call defrr(osym,   nlml*nlml*nrc)
      call symprj(ia,nlml,posc,nrc,g,ag,ng,cy,w(osym))

c ------ assemble uu,us,ss parts of density for each spin ----------
      do 30 isp=1,nsp
        nrlm = nr1*nlml
        jsp = isp-1
        call defrr(ou1,    nr1*(lmaxa+1))
        call defrr(ou2,    nr1*(lmaxa+1))
        call defrr(othet,  nlml*nrc*(lmaxa+1)**2)
        call defrr(oh1,    nr1)
        call dpzero(w(orhol),  nrlm)
        call defrr(ovx,    nr1)
        call dpscop(w(ov0(ib)),w(ovx),nr1,1+nr1*jsp,1,1d0)
C --- call makusp from lmf, single panel only, lrel in common /cc/
        call defrr(opz0,  -10*nsp)
        call defrr(ogz,   nr1*(lmaxa+1)*nsp)
        call defrr(oruu,  nr1*(lmaxa+1)*2*nsp)
        call defrr(orus,  nr1*(lmaxa+1)*2*nsp)
        call defrr(orss,  nr1*(lmaxa+1)*2*nsp)
        call makusp(n0,z(is),1,1,rmt(is),lmaxa,w(ovx),a(is),nr1,
     .    0d0,0d0,pnu(1,isp,ib),w(opz0),w,w,w(ou1),w(ou2),w(ogz),
     .    w(oruu),w(orus),w(orss))
        call rlse(opz0)
        k=jsp*ioffp(nbas+1)
        call symchd(lmaxa,nlml,nr1,w(ou1),w(ou1),w(orhol),quu(1+k),
     .    ipbas,ioffp,w(osym),nrc,cg,indxcg,jcg,w(othet),w(oh1))
        call symchd(lmaxa,nlml,nr1,w(ou1),w(ou2),w(orhol),qus(1+k),
     .    ipbas,ioffp,w(osym),nrc,cg,indxcg,jcg,w(othet),w(oh1))
        call symchd(lmaxa,nlml,nr1,w(ou2),w(ou2),w(orhol),qss(1+k),
     .    ipbas,ioffp,w(osym),nrc,cg,indxcg,jcg,w(othet),w(oh1))
c|      call rhinfo(nr1,nlml,nsp,w(orhol),w(orofi),w(orwgt),a(is))
        call rlse(ou1)
        call dpdot(w(orhol),w(orwgt),nr1,sum)
        qsph=qsph+nrc*srfpi*sum
        call dpscop(w(orhol),w(orhosp),nrlm,1,1+jsp*nrlm,1d0)
        if (ipan == 1)
     .  call dpsadd(w(orhosp),w(orhoc),nr1,1+nrlm*jsp,1+nr1*jsp,1/srfpi)
        call rlse(ou1)
   30 continue

c ------ add rhol with proper symmetry into w(orho) ------
        do 26 ja=1,nrc
        jb=ipbas(ja)
        or=orho(jb)
        if(ipan == 1) call dpzero(w(or),nr1*nlml*nsp)
        call addrhl(ja,nr1,nlml,nsp,nrc,w(osym),w(orhosp),w(or))
        call rhinfo(nr1,nlml,nsp,w(or),w(orofi),w(orwgt),a(is))
   26   continue

      call rlse(orofi)
  50  continue
      call tcx('sprho2')
      end
c --------- sub rhinfo: print out info about sphere density -----
      subroutine rhinfo(nr1,nlml,nsp,rhol,rofi,rwgt,a)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhol(nr1,nlml,nsp),rofi(nr1),rwgt(nr1),
     .  ff(0:4),qsph(2),roav(2)
      call getpr(ipr)
      b=rofi(nr1)/(dexp(a*nr1-a)-1.d0)
      srfpi=dsqrt(16d0*datan(1d0))
      if(ipr >= 40) write(6,301)
      do 12 isp=1,nsp
      if(ipr >= 40 .and. nsp == 2) print '('' spin'',i2,'':'')', isp
      do 10 ilm=1,nlml
      l=ll(ilm)
      top=rhol(1,ilm,isp)
      bot=top
      qmom=0d0
      do 11 ir=2,nr1
      top=dmax1(top,rhol(ir,ilm,isp))
      bot=dmin1(bot,rhol(ir,ilm,isp))
   11 qmom=qmom+rwgt(ir)*rhol(ir,ilm,isp)*rofi(ir)**l
      rhov=rhol(nr1,ilm,isp)/rofi(nr1)**2
      do 29 i=0,4
   29 ff(i)=rhol(i+nr1-4,ilm,isp)/rofi(i+nr1-4)**2
      rhos=(25*ff(4)-48*ff(3)+36*ff(2)-16*ff(1)+3*ff(0))/12d0
      rhos=rhos/(a*(rofi(nr1)+b))
      xx=dmax1(dabs(top),dabs(bot),dabs(qmom),dabs(rhov),dabs(rhos))
      if(xx > 1d-5.and.ipr >= 35)
     .   write(6,300) ilm,bot,top,qmom,rhov,rhos
  300 format(i5,5f13.6)
  301 format(/'   ilm',6x,'bot',10x,'top',10x,'qmom',
     .   7x,'val(rmt)',5x,'slo(rmt)')
      if(ilm == 1) qsph(isp)=srfpi*qmom
   10 if(ilm == 1) roav(isp)=rhov/srfpi
      if(ipr >= 20) write(6,750) qsph(isp),roav(isp)
  750 format(' charge in sphere=',f11.6,'   avg density at rmt=',f11.6)
   12 continue
      if (nsp == 2 .and. ipr > 1) then
        print 751, qsph(1)+qsph(2),qsph(1)-qsph(2),roav(1)+roav(2)
      endif
  751 format('     total charge=',f11.6,'  moment=',f11.6,
     .  '  density=',f11.6)
      end
c --------- sub addrhl -------------------
      subroutine addrhl(ja,nr1,nlml,nsp,nrc,sym,rhol,rhut)
      implicit real*8 (a-h,p-z), integer (o)
      dimension sym(nlml,nlml,nrc),rhol(nr1,nlml,1),rhut(nr1,nlml,1)

      do 10 isp=1,nsp
      do 10 ilm=1,nlml
      do 10 jlm=1,nlml
      sym1=nrc*sym(jlm,ilm,ja)
      if(dabs(sym1) > 1d-9) then
        do 11 ir=1,nr1
   11   rhut(ir,ilm,isp)=rhut(ir,ilm,isp)+sym1*rhol(ir,jlm,isp)
      endif
  10  continue
      end
