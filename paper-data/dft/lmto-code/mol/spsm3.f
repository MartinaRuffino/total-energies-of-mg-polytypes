      subroutine spsm3(ib0,r1,nr1,a1,nlml,nxi,lxi,exi,n0,rmt,rsm,
     .  rint,pos,nspec,nbas,ips,rhoi,poti,pot0,pot00,ioff,rho0,
     .  ioffv0,cy,cg,jcg,indxcg,ioffp,lmxa,nel,el,nx,nph,nrad,pjj,vab0,
     .  rhol,vesl,vval,nri,nvi,zets,zet0s,sxi0,sg0,rep,rmu,rph,qsm)
c  Sphere program for smooth interstitial density and potential.
c  Tail density and potential are made on mesh inside the sphere.
c  Output:  rhol = l-decomposed tail density
c           vesl = l-decomposed tail electrostatic potential
c           vval = values of head electrostatic potential on sphere
c           zets = integrals of smooth_pot* xi_m over sphere
c           zet0s= integrals of smooth_pot* gaussian over sphere
c           sxi0 = integral of xi excluding spheres
c           sg0  = integral of gaussians excluding spheres
c           qsm  = smooth-density charge inside sphere
c           rep,rmu,rph: integrals of density times exc,vxc,ves
      parameter( nnn=100, nlmx= 49, npwr=5 )
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(n0,1),nxi(1),ips(1),lxi(n0,1),ioff(1),ioffv0(1),
     .  rhol(nr1,nlml),rhoi(1),poti(1),pot0(1),pos(3,1),cy(1),rint(1),
     .  p(3,nnn),wp(nnn),rad(nnn+1),rmt(1),rsm(1),vesl(nr1,nlml),
     .  zets(nri,nbas),zet0s(nvi,nbas),el(1),ioffp(1),rho0(1),
     .  spk0(0:npwr),sxi0(nri),sg0(nvi),
     .  wrad(nnn),pot00(1),vval(nlml),rval(nlmx),pjj(1),vab0(9,nbas,1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('spsm3')
      pi=4d0*datan(1d0)
      y0=1d0/dsqrt(4d0*pi)
      if(nlml > nlmx) call rx('spsm3: increase nlmx')
      lmaxl=ll(nlml)
      js=ips(ib0)
      jr=ioff(ib0)+1
      jv=ioffv0(ib0)+1
      rsmg=0.5d0*r1
      call defrr(orofi,   nr1)
      call radmsh(r1,a1,nr1,w(orofi))
      do 22 i=1,nri
   22 zets(i,ib0)=0d0
      do 23 i=1,nvi
   23 zet0s(i,ib0)=0d0
      is0=ips(ib0)
c ------ set radial and angular meshes ----------
      call fpiint(nx,nph,np,p,wp)
      if(nrad > nnn.or.np > nnn) call rx('spsm3: incr nnn')
      call mklegw(nrad,rad,wrad,0)
      do  1  i = 1, nrad
      rad(i) = r1/2*rad(i) + r1/2
    1 wrad(i) = wrad(i)*r1/2
      if(iprint() >= 30) write(6,200) nrad,nx,nph,np
  200 format(/'  nrad,nx,nph=',3i5,'   gives   np=',i4)
      call defrr(oyl,   nlml*np)
      call setylm(nlml,np,p,cy,w(oyl))

c ------ define some arrays -------------------
      nlmbx=49
      call defrr(orhal,       nrad*nlml)
      call defrr(ovhal,       nrad*nlml)
      call defrr(orhokl,      (npwr+1)*nlml)
      call defrr(oveskl,      (npwr+1)*nlml)
      call defrr(ockl,        (npwr+1)*nlml*nlmbx)

c ------ loop over other sites; accumulate tail density and ves ----
      call dpzero(w(orhokl),  (npwr+1)*nlml)
      call dpzero(w(oveskl),  (npwr+1)*nlml)
      do 10 ib=1,nbas
      if(ib == ib0) goto 10
      is=ips(ib)
      dd=dsqrt((pos(1,ib0)-pos(1,ib))**2+(pos(2,ib0)-pos(2,ib))**2
     .         +(pos(3,ib0)-pos(3,ib))**2)
      if(dd > rint(is)+r1) goto 10
      ir=ioff(ib)+1
c ---------- this part: exi /= 0 --------------------
      lmax=-1
      do 11 ie=1,nxi(is)
      lmax=max0(lmax,lxi(ie,is))
      nlmb=(lxi(ie,is)+1)**2
      if(nlmb > nlmbx) call rx('spsm3: increase nlmbx')
      call hsmxpn(pos(1,ib),pos(1,ib0),exi(ie,is),rsm(is),rsmg,
     .   nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call xxsp11(npwr,nlml,nlmb,rhoi(ir),w(ockl),w(orhokl))
      call xxsp11(npwr,nlml,nlmb,poti(ir),w(ockl),w(oveskl))
  11  ir=ir+nlmb
c ---------- this part: exi == 0 and gaussians ------
      nlmb=(lmax+1)**2
      iv=ioffv0(ib)+1
      call hsmxpn(pos(1,ib),pos(1,ib0),0d0,rsm(is),rsmg,
     .   nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call xxsp11(npwr,nlml,nlmb,pot0(iv),w(ockl),w(oveskl))
      call hsmxpn(pos(1,ib),pos(1,ib0),0d0,0d0,rsmg,
     .   nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call xxsp11(npwr,nlml,nlmb,pot00(iv),w(ockl),w(oveskl))
      call gauxpn(pos(1,ib),pos(1,ib0),rsm(is),rsmg,
     .  nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call xxsp11(npwr,nlml,nlmb,rho0(iv),w(ockl),w(orhokl))
  10  continue

      call xxsp12(nrad,rad,rsmg,npwr,nlml,w(orhokl),w(ockl),w(orhal))
      call xxsp12(nrad,rad,rsmg,npwr,nlml,w(oveskl),w(ockl),w(ovhal))

c ------ interpolate tail rho and ves to fine radial mesh ----
      nn=10
      call splpol(nn,nrad,nlml,r1,w(orhal),nr1,w(orofi),rhol,001)
      call splpol(nn,nrad,nlml,r1,w(ovhal),nr1,w(orofi),vesl,000)

c ------ add smooth head density -----------
      call defrr(oxi,    nrad*(lmaxl+1))
      call sphed3(nrad,rad,wrad,lmaxl,nxi(js),lxi(1,js),exi(1,js),
     .   rsm(js),w(oxi),rhoi(jr),rho0(jv),poti(jr),pot0(jv),pot00(jv),
     .   w(orhal),w(ovhal))
      call rlse(oxi)

c ------ make values of rho,ves on sphere surface --------
      call dpzero(rval,   nlml)
      call dpzero(vval,   nlml)
      call spval3(nlml,nxi(js),lxi(1,js),exi(1,js),r1,rsm(js),
     .   rhoi(jr),rho0(jv),poti(jr),pot0(jv),pot00(jv),rval,vval)
      if(ipr >= 20) write(6,330) y0*vval(1),y0*vesl(nr1,1),
     .   y0*(vesl(nr1,1)+vval(1))
  330 format(/' vmad...   head=',f11.6,'   tail=',f11.6,'   tot=',f11.6)

c ------ assemble lmaxl-projected smooth density and pot on mesh ----
      n=np*nrad
      call defrr(orho,  n)
      call defrr(oves,  n)
      call spmkrh(nrad,nlml,np,w(orhal),w(oyl),w(orho))
      call spmkrh(nrad,nlml,np,w(ovhal),w(oyl),w(oves))

c ------ put lmaxl-projection of mesh pot (times weights) into ves ---
      call defrr(oexc,     n)
      call defrr(ovxc,     n)
      call ropevx(w(orho),w(oexc),w(ovxc),n)
C|    PRINT *, 'SPSM3, ZERO OUT VXC  FOR NOW'
C|    CALL DPZERO(W(OVXC), N)
c the next two lines project out l < lmaxl from vxc
      call sptlpj(nrad,nlml,np,w(ovxc),w(oyl),wp,w(ovhal))
      call spmkrh(nrad,nlml,np,w(ovhal),w(oyl),w(ovxc))
      call spmkxc(nrad,rad,wrad,np,wp,w(orho),w(oexc),w(ovxc),
     .   w(oves),rep,rmu,rph,qsm)
      call rlse(oexc)

c ------ on-site integrals zets,zet0s ---------------------
      call defrr(oxi,    nrad*(lmaxl+1))
      ir=ioff(ib0)+1
      call zthead(nrad,rad,nlml,nxi(js),lxi(1,js),exi(1,js),r1,rsm(is0),
     .  w(oxi),np,w(oves),w(oyl),w(ovhal),zets(ir,ib0),sxi0(ir))
      ir=ioffv0(ib0)+1
c ... z0head uses vhal already assembled (zthead).
      call z0head(nrad,rad,nlml,lmax,r1,rsm(is0),w(oxi),w(ovhal),
     .  zet0s(ir,ib0),sg0(ir))
      call rlse(oxi)

c ------ perturbation matrices (j vsm j): vhal scaled by weights -----
      npjj=ioffp(nbas+1)
      call spvls2(nlml,lmxa,nrad,nel,el,rad,nbas,w(ovhal),
     .   indxcg,jcg,cg,pjj,npjj,ioffp(ib0))
      call defrr(oxi,    nrad*(lmxa+1))
      call defrr(oxi2,   nrad*(lmxa+1))
      call defrr(oh1,    nrad)
      call defrr(oh2,    nrad)
      call spvab0(lmxa,ib0,nbas,nel,el,nrad,rad,
     .  w(ovhal),w(oh1),w(oh2),w(oxi),w(oxi2),vab0)
      call rlse(oxi)

c ------ loop over other sites; make integrals zets --------
      call defrr(orhokl,      (npwr+1)*nlml)
      call xxsp13(nrad,rad,wrad,rsmg,npwr,nlml,w(ockl),w(ovhal),
     .  w(orhokl),spk0)
      do 20 ib=1,nbas
      if(ib == ib0) goto 20
      is=ips(ib)
      dd=dsqrt((pos(1,ib0)-pos(1,ib))**2+(pos(2,ib0)-pos(2,ib))**2
     .         +(pos(3,ib0)-pos(3,ib))**2)
      if(dd > rint(is)+r1) goto 20
      ir=ioff(ib)+1
      do 21 ie=1,nxi(is)
      nlmb=(lxi(ie,is)+1)**2
      call hsmxpn(pos(1,ib),pos(1,ib0),exi(ie,is),rsm(is),rsmg,
     .  nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call xxsp14(npwr,nlml,nlmb,w(ockl),w(orhokl),zets(ir,ib0),
     .  spk0,sxi0(ir))
  21  ir=ir+nlmb
      nlmb=ioffv0(ib0+1)-ioffv0(ib0)
      call gauxpn(pos(1,ib),pos(1,ib0),rsm(is),rsmg,
     .  nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      iv=ioffv0(ib)+1
      call xxsp14(npwr,nlml,nlmb,w(ockl),w(orhokl),zet0s(iv,ib0),
     .  spk0,sg0(iv))
   20 continue

      sam=0d0
      do 24 i=1,nri
   24 sam=sam+rhoi(i)*zets(i,ib0)
      sum=0d0
      do 25 i=1,nvi
   25 sum=sum+rho0(i)*zet0s(i,ib0)
      if(ipr >= 50) write(6,986) sam,sum
  986 format('zets: L-truncated integral(rhoi,rho0*vtot):',2f12.6)

      call tcx('spsm3')
      end
c --------- sub sptlpj: projects out to lmaxl ----------
      subroutine sptlpj(nrad,nlml,np,rho,yl,wp,rhal)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rho(nrad,np),yl(nlml,np),wp(np),rhal(nrad,nlml)
      call dpzero(rhal,nrad*nlml)
      do 10 ilm=1,nlml
      do 10 ip=1,np
      do 10 ir=1,nrad
  10  rhal(ir,ilm)=rhal(ir,ilm)+rho(ir,ip)*wp(ip)*yl(ilm,ip)
      end
c --------- spmkrh: make point-wise density from yl-expansion ----
      subroutine spmkrh(nrad,nlml,np,rhal,yl,rho)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rho(nrad,np),yl(nlml,np),rhal(nrad,nlml)
      call dpzero(rho,nrad*np)
      do 10 ip=1,np
      do 10 ilm=1,nlml
      do 10 ir=1,nrad
  10  rho(ir,ip)=rho(ir,ip)+rhal(ir,ilm)*yl(ilm,ip)
      end
c --------- spmkxc: add xc potential to ves; multiply in weights ---
      subroutine spmkxc(nrad,rad,wrad,np,wp,rho,exc,vxc,ves,
     .   rep,rmu,rph,qsm)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rad(1),wp(1),wrad(1),rho(nrad,np),
     .   exc(nrad,np),vxc(nrad,np),ves(nrad,np)
      qsm=0d0
      rep=0d0
      rmu=0d0
      rph=0d0
      do 1 ip=1,np
      do 1 ir=1,nrad
      weight=(rad(ir)**2*wrad(ir))*wp(ip)
      qsm=qsm+rho(ir,ip)*weight
      rep=rep+exc(ir,ip)*rho(ir,ip)*weight
      rmu=rmu+vxc(ir,ip)*rho(ir,ip)*weight
      rph=rph+ves(ir,ip)*rho(ir,ip)*weight
  1   ves(ir,ip)=(ves(ir,ip)+vxc(ir,ip))*weight
      write(6,725) qsm,rep,rmu,rph
  725 format(/' spmkxc:  q=',f10.6,'   rep=',f10.6,'   rmu=',f10.6,
     .   '   rph=',f11.6)
      end
c --------- xxsp11: add together rhokl ------------
      subroutine xxsp11(npwr,nlml,nlmb,rhoi,ckl,rhokl)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhoi(1),rhokl(0:npwr,1),ckl(0:npwr,nlml,nlmb)
      do 14 ilmb=1,nlmb
      do 14 ilml=1,nlml
      do 14 k=0,npwr
  14  rhokl(k,ilml)=rhokl(k,ilml)+rhoi(ilmb)*ckl(k,ilml,ilmb)
      end
c --------- xxsp12: add together rhal -------------
      subroutine xxsp12(nrad,rad,rsmg,npwr,nlml,rhokl,pkl,rhal)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rad(1),rhal(nrad,nlml),pkl(0:npwr,0:1),
     .   rhokl(0:npwr,1)
      call dpzero(rhal,   nrad*nlml)
      lmaxl=ll(nlml)
      do 20 ir=1,nrad
      call pauskl(rad(ir),rsmg,pkl,lmaxl,npwr,npwr)
      do 21 ilm=1,nlml
      l=ll(ilm)
      do 21 k=0,npwr
  21  rhal(ir,ilm)=rhal(ir,ilm)+rhokl(k,ilm)*pkl(k,l)*rad(ir)**l
  20  continue
      end
c --------- xxsp13: make integrals vhal*pkl in rhokl -----
      subroutine xxsp13(nrad,rad,wrad,rsmg,npwr,nlml,pkl,vhal,
     .  rhokl,spk0)
c  makes spk0 = integral p_k0  in sphere and rhokl= int with pot.
      implicit real*8 (a-h,p-z), integer (o)
      dimension rad(1),wrad(1),vhal(nrad,nlml),pkl(0:npwr,0:1),
     .   rhokl(0:npwr,1),spk0(0:npwr)
      srfpi=dsqrt(16d0*datan(1d0))
      call dpzero(rhokl,  (npwr+1)*nlml)
      do 1 k=0,npwr
    1 spk0(k)=0d0
      lmaxl=ll(nlml)
      do 20 ir=1,nrad
      call pauskl(rad(ir),rsmg,pkl,lmaxl,npwr,npwr)
      do 2 k=0,npwr
    2 spk0(k)=spk0(k)+srfpi*wrad(ir)*pkl(k,0)*rad(ir)**2
      do 21 ilm=1,nlml
      l=ll(ilm)
      do 21 k=0,npwr
  21  rhokl(k,ilm)=rhokl(k,ilm)+vhal(ir,ilm)*pkl(k,l)*rad(ir)**l
  20  continue
      end
c --------- xxsp14: add to zets -------------------
      subroutine xxsp14(npwr,nlml,nlmb,ckl,rhokl,zets,spk0,sxi0)
      implicit real*8 (a-h,p-z), integer (o)
      dimension zets(1),rhokl(0:npwr,1),ckl(0:npwr,nlml,nlmb),
     .   spk0(0:npwr),sxi0(1)
      do 14 ilmb=1,nlmb
      do 14 ilml=1,nlml
      do 14 k=0,npwr
   14 zets(ilmb)=zets(ilmb)+rhokl(k,ilml)*ckl(k,ilml,ilmb)
      do 15 ilmb=1,nlmb
      do 15 k=0,npwr
   15 sxi0(ilmb)=sxi0(ilmb)-spk0(k)*ckl(k,1,ilmb)
      end
