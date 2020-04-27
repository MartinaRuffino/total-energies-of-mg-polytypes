      subroutine roisym(nbas,pos,ng,g,ag,nclas,ipclas,ips,nxi,lxi,
     .   n0,ioff,cy,rhoi,f)
c- symmetrize interstitial density and forces.
C 6 Jan 95 spin polarized (MvS)
      parameter( ncmx=20 )
      implicit real*8 (a-h,p-z), integer (o)
      dimension ipbas(ncmx),posc(3,ncmx),g(9,1),ag(3,1),ipclas(1),
     .  ips(1),pos(3,1),cy(1),nxi(1),lxi(n0,1),rhoi(1),f(3,1),
     .  ioff(1),x(200)
      real w(1)
      common /w/ w
      call tcn('roisym')
      call getpr(ipr)
      nsp=lsp()+1
      nri=ioff(nbas+1)
      if(ipr >= 20) write(6,220) nclas,ng
  220 format(/' roisym:  number of classes=',i3,'     group ops=',i3)
      if(ng <= 1) then
        call tcx('roisym')
        return
      endif
c ------ start loop over classes: make list of atoms in class -----
      do 50 ic=1,nclas
      nrc=0
      do 1 ib=1,nbas
      if(ipclas(ib) == ic) then
        nrc=nrc+1
        if(nrc > ncmx) call rx('roisym incr ncmx')
        ipbas(nrc)=ib
        do 2 m=1,3
  2     posc(m,nrc)=pos(m,ib)
      endif
  1   continue
c ------ representative atom is first in the class -----
      ia=1
      ib=ipbas(ia)
      is=ips(ib)
      lmax=-1
      do 3 ie=1,nxi(is)
  3   lmax=max0(lmax,lxi(ie,is))
      nlmx=(lmax+1)**2
      if(ipr >= 40) write(6,910) ic,is,nrc,ib
  910 format(/' symmetry class',i3,'    species=',i3,'    atoms:',i3,
     .   '    first atom: ib=',i3)

c ------ make symmetry projectors --------------
      call defrr(osym,   nlmx*nlmx*nrc)
      call symprj(ia,nlmx,posc,nrc,g,ag,ng,cy,w(osym))
c ... symmetrize forces
      call xxrsyf(nlmx,nrc,w(osym),ipbas,f)

c ------ symmetrize rhoi ---------------
      do 40 isp=1,nsp
      joff=1+(isp-1)*nri
      do 40 ie=1,nxi(is)
      lx=lxi(ie,is)
      nlm=(lx+1)**2
c ... make symmetrized coeffs on site ib
      do 41 ilm=1,nlm
  41  x(ilm)=0d0
      do 42 ja=1,nrc
      jb=ipbas(ja)
      koff=ioff(jb)+joff
  42  call xxrsy1(nlmx,nlm,w(osym),ja,rhoi(koff),x)
c ... put into rhoi
      do 22 ilm=1,nlm
  22  x(ilm)=nrc*x(ilm)
      do 44 ja=1,nrc
      jb=ipbas(ja)
      koff=ioff(jb)+joff
  44  call xxrsy2(nlmx,nlm,w(osym),ja,x,rhoi(koff))
  40  joff=joff+nlm
      call rlse(osym)
  50  continue
      call tcx('roisym')
      end
c ------ help subs ---------
      subroutine xxrsy1(nlmx,nlm,sym,ja,rhoi,x)
      implicit real*8 (a-h,p-z), integer (o)
      dimension sym(nlmx,nlmx,1),x(nlm),rhoi(nlm)
      do 42 ilm=1,nlm
      do 42 jlm=1,nlm
  42  x(ilm)=x(ilm)+sym(ilm,jlm,ja)*rhoi(jlm)
      end
      subroutine xxrsy2(nlmx,nlm,sym,ja,x,rhoi)
      implicit real*8 (a-h,p-z), integer (o)
      dimension sym(nlmx,nlmx,1),x(nlm),rhoi(nlm)
      do 1 ilm=1,nlm
  1   rhoi(ilm)=0d0
      do 42 ilm=1,nlm
      do 42 jlm=1,nlm
  42  rhoi(jlm)=rhoi(jlm)+sym(ilm,jlm,ja)*x(ilm)
      end
c ------ sub to symmetrize forces ----------
      subroutine xxrsyf(nlmx,nrc,s,ipbas,f)
      implicit real*8 (a-h,p-z), integer (o)
      dimension x(3),f(3,1),s(nlmx,nlmx,1),ipbas(1)
      x(1)=0d0
      x(2)=0d0
      x(3)=0d0
      do 2 ja=1,nrc
      jb=ipbas(ja)
      x(1)=x(1)+s(4,4,ja)*f(1,jb)+s(4,2,ja)*f(2,jb)+s(4,3,ja)*f(3,jb)
      x(2)=x(2)+s(2,4,ja)*f(1,jb)+s(2,2,ja)*f(2,jb)+s(2,3,ja)*f(3,jb)
  2   x(3)=x(3)+s(3,4,ja)*f(1,jb)+s(3,2,ja)*f(2,jb)+s(3,3,ja)*f(3,jb)
      do 3 ja=1,nrc
      jb=ipbas(ja)
      f(1,jb)=(s(4,4,ja)*x(1)+s(2,4,ja)*x(2)+s(3,4,ja)*x(3))*nrc
      f(2,jb)=(s(4,2,ja)*x(1)+s(2,2,ja)*x(2)+s(3,2,ja)*x(3))*nrc
  3   f(3,jb)=(s(4,3,ja)*x(1)+s(2,3,ja)*x(2)+s(3,3,ja)*x(3))*nrc
      end
