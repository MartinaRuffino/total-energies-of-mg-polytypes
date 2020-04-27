      subroutine spsmt(ib0,r1,nr1,a1,nlml,nxi,lxi,exi,n0,rsm,rint,
     .  pos,nbas,ips,rhoi,ioff,rho0,ioffv0,cy,cg,jcg,indxcg,
     .  nx,nph,nrad,rhol)
C- Sphere program to obtain tail charge density on a single atom
c  Output:  rhol = l-decomposed tail density
      implicit none
      integer nbas,nnn,nlmx,npwr,n0,nr1,nlml,ib0,nph,nrad,nx
      parameter( nnn=100, nlmx= 49, npwr=5 )
      integer nxi(1),ips(1),lxi(n0,1),ioff(1),indxcg(1),ioffv0(1),jcg(1)
      double precision exi(n0,1),rho0(1),cg(1),r1,
     .  rhol(nr1,nlml),rhoi(1),pos(3,1),cy(1),rint(1),
     .  p(3,nnn),wp(nnn),rad(nnn+1),rsm(1),wrad(nnn)
      integer i,ib,ie,ipr,iprint,ir,is,iv,ll,
     .  lmax,nlmb,nlmbx,nn,np,orofi,oyl,ockl,orhokl,orhal
      double precision a1,dd,rsmg

      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('spsmt')
      if(nlml > nlmx) call rx('spsm3: increase nlmx')
      rsmg=0.5d0*r1
      call defrr(orofi,   nr1)
      call radmsh(r1,a1,nr1,w(orofi))

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
      call defrr(orhokl,      (npwr+1)*nlml)
      call defrr(ockl,        (npwr+1)*nlml*nlmbx)

c ------ loop over other sites; accumulate tail density ----
      call dpzero(w(orhokl),  (npwr+1)*nlml)
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
  11  ir=ir+nlmb
c ---------- this part: exi == 0 and gaussians ------
      nlmb=(lmax+1)**2
      iv=ioffv0(ib)+1
      call hsmxpn(pos(1,ib),pos(1,ib0),0d0,rsm(is),rsmg,
     .   nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call hsmxpn(pos(1,ib),pos(1,ib0),0d0,0d0,rsmg,
     .   nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call gauxpn(pos(1,ib),pos(1,ib0),rsm(is),rsmg,
     .  nlmb,nlml,npwr,cg,indxcg,jcg,cy,w(ockl))
      call xxsp11(npwr,nlml,nlmb,rho0(iv),w(ockl),w(orhokl))
  10  continue

      call xxsp12(nrad,rad,rsmg,npwr,nlml,w(orhokl),w(ockl),w(orhal))

c ------ interpolate tail rho to fine radial mesh ----
      nn=10
      call splpol(nn,nrad,nlml,r1,w(orhal),nr1,w(orofi),rhol,001)

      call rlse(orofi)
      call tcx('spsmt')
      end
