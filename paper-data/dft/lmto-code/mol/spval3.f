      subroutine spval3(nlml,nxi,lxi,exi,r1,rsm,
     .   rhoi,rho0,poti,pot0,pot00,rval,vval)
C  adds head-part to rho and ves on sphere surface
      implicit real*8 (a-h,p-z), integer (o)
      dimension exi(1),lxi(1),rhoi(1),poti(1),pot0(1),pot00(1),
     .   rval(1),vval(1),xi(0:20),xi0(0:20),rho0(1)
      call getpr(ipr)
      lmaxl=ll(nlml)
      asm=1d0/rsm
      pi=4d0*datan(1d0)
C ------ contribution from smooth Hankels, e /= 0 ------
      i0=0
      do 10 ie=1,nxi
      lx=min0(lxi(ie),lmaxl)
      nlmx=(lx+1)**2
      call hansmr(r1,exi(ie),asm,xi,lx)
      do 11 ilm=1,nlmx
      l=ll(ilm)
      rval(ilm)=rval(ilm)+rhoi(ilm+i0)*(xi(l)*r1**l)
  11  vval(ilm)=vval(ilm)+poti(ilm+i0)*(xi(l)*r1**l)
  10  i0=i0+(lxi(ie)+1)**2
c ------- contribution from smooth hankels with e=0 --------
      call hansmr(r1,0d0,asm,xi,lmaxl)
      call ropgau(rsm,lmaxl,1,r1,xi0,1)
      do 8 ilm=1,nlml
      l=ll(ilm)
      rval(ilm)=rval(ilm)+rho0(ilm)*xi0(l)
  8   vval(ilm)=vval(ilm)+pot0(ilm)*(xi(l)*r1**l)
c ------- contribution from normal hankels with e=0 --------
      xi(0)=1d0/r1
      do 23 l=1,lmaxl
  23  xi(l)=(2*l-1d0)*xi(l-1)/r1
      do 24 ilm=1,nlml
      l=ll(ilm)
  24  vval(ilm)=vval(ilm)+pot00(ilm)*xi(l)

c --------- printout ------------------
      if(ipr >= 40) then
      write(6,201)
      do 20 j=1,nlml
      top=dmax1(dabs(rval(j)),dabs(vval(j)))
  20  if(top > 1d-5) write(6,200) j,rval(j),vval(j)
  200 format(i12,3f14.6)
  201 format(/' spval2:  ilm',7x,'rval',10x,'vval       (heads only)')
      endif

      end
