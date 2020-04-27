      subroutine ovldum(nspec,spid,r,rsm,lmxl,nxi,lxi,exi,nr,a,f,n0,
     .   nbas,ips,orho,rhox,nrx,rhoi,nri)
c  dummy input routine for atomic charge densities
      implicit real*8 (a-h,p-z), integer (o)
      integer intopt,nglob
      dimension rhox(nrx,nspec),f(n0,1),exi(n0,1),nr(1),a(1),r(1),
     .   nxi(1),ips(1),rhoi(nri),lxi(n0,1),orho(1),ex(20),lmxl(1),
     .   rsm(1),xi(0:10),x0(0:10),xdum(0:10)
      character*4 spid(1),spix
      real w(1)
      common /w/ w
      call getpr(ipr)
      fuzz=1.d-8
      pi=4d0*datan(1d0)
      srfpi=dsqrt(4d0*pi)
      y0=1d0/srfpi
c ------ read free-atom densities --------
      if(ipr >= 30) write(6,*) ' '
      jfi=78
      rewind jfi
      do 33 k=1,nspec
      read(jfi,401) spix,rxx,a(k),nr(k),nx
      read(jfi,402) (ex(m),m=1,nx)
      if (ipr >= 30) write(6,843) spix,rxx,nr(k)
  843 format(' read fa density, sp=',a4,'   rx=',f12.6,'   nr=',i5)
      lok=1
      if(spix /= spid(k).or.dabs(r(k)-rxx) > fuzz) lok=0
      if(nx /= nxi(k)) lok=0
      do 3 m=1,nx
  3   if(dabs(ex(m)-exi(m,k)) > fuzz) lok=0
      if(lok == 0) call rx('ovldum: bad data, spec=<'//spid(k)//'>')
      if(nr(k) > nrx) call rx('ovldum: nrx exceeded')
      read(jfi,402) (f(m,k),m=1,nx)
  33  read(jfi,410) (rhox(i,k),i=1,nr(k))
  401 format(1x,a4,2x,2f15.9,2i8)
  402 format(1p,4d18.10)
  410 format(1p,4d18.10)
c ------ define/occupy arrays for atomic densities -------
      do 10 ib=1,nbas
      is=ips(ib)
      nr1=nr(is)
      nlml=(lmxl(is)+1)**2
      call defrr(orho(ib),     nr1*nlml)
      call dpzero(w(orho(ib)), nr1*nlml)
      call dpcopy(rhox(1,is),w(orho(ib)),1,nr1,y0)
          call defrr(orwgt,    nr1)
          intopt = 10*nglob('lrquad')
          call radwgt(intopt,r(is),a(is),nr1,w(orwgt))
          call dpdot(w(orho(ib)),w(orwgt),nr1,qq)
          if(ipr >= 30) write(6,881) qq/y0
  881     format(' charge in sphere=',f12.6)
          call rlse(orwgt)
  10  continue

c ------ interstitial density -------
      if(ipr >= 40) write(6,333)
  333 format(/'  ib  spec     i         exi          rhoi')
      call dpzero(rhoi,nri)
      i=0
      do 30 ib=1,nbas
      is=ips(ib)
      do 40 j=1,nxi(is)
      rhoi(i+1)=f(j,is)
      if (ipr >= 40) write(6,334) ib,spid(is),i+1,exi(j,is),rhoi(i+1)
  334 format(i4,2x,a4,i6,f12.3,f15.5)
   40 i=i+(lxi(j,is)+1)**2
   30 continue

      end
