      subroutine oncmak(lx,ex,nx,lp1,ep1,lp2,ep2,r,rsm,ri,ifi)
c  on-site expansion of products of smoothed hankels
      implicit real*8 (a-h,p-z), integer(o)
      dimension ex(50),lx(50),err(50),jp(50)
      real w(1)
      common /w/ w
      call getigv(1,ipr)
      lmax=0
      do 2 j=1,nx
  2   lmax=max0(lmax,lx(j))
      lmax=min0(lmax,lp1+lp2)
      ncof=(lmax+1)*nx*(lp1+1)*(lp2+1)
      call defrr(ocof,        ncof)
      call dpzero(w(ocof),    ncof)
c ------ print out specifications for functions -----
      if(ipr >= 20) then
      write(6,991) lmax,r,rsm,ri
      write(6,990) (lx(i),i=1,nx),(200,i=nx+1,6),(ex(i),i=1,nx)
      write(6,992) lp1,ep1,lp2,ep2
  991 format(/' hyfmk0:  lmax=',i3,'   r,rsm,ri=',3f8.4)
  990 format(' xi:   l= ',6(1x,i1),'   e=',6f7.2)
  992 format(' phi1:  lmx=',i3,'    e=',f8.4
     .  /' phi2:  lmx=',i3,'    e=',f8.4)
      endif
c ------ setup for radial integration --------------
      rmax=30d0
      aa=0.03d0
      nmax=401
      call defrr(or,   nmax)
      call defrr(ow,   nmax)
      call radint(ri,aa,rmax,nmax,w(or),w(ow),nr)
c ------ write to log -------------------
CL      write(71,710) r,rsm,lmax
  710 format(' ------- oncmak:  r',f8.4,'  rsm',f8.4,'  lmax',i3)
CL      write(71,993) ri,rmax,aa,nr
CL      write(71,996) (lx(i),i=1,nx),(200,i=nx+1,6),(ex(i),i=1,nx)
CL      write(71,997) lp1,lp2,ep1,ep2
  993 format(' ri',f8.4,'    rmax',f8.2,'    a',f7.3,'   nr',i5)
  996 format(' xi:   l= ',6(1x,i1),'   e=',6f7.2)
  997 format(' phi   lmx=',2i2,'   e=',2f7.2)
c ------ make radial parts of the phi1,phi2 -------
      call defrr(ophi1,    nr*(lp1+1))
      call defrr(ophi2,    nr*(lp2+1))
      call xxfmk0(ep1,lp1,rsm,w(or),nr,w(ophi1))
      call xxfmk0(ep2,lp2,rsm,w(or),nr,w(ophi2))
c ------ start big loop over lxi: get dimensions ---------
      do 80 lxi=0,lmax
      ndim=0
      do 1 j=1,nx
      if(lx(j) >= lxi) ndim=ndim+1
  1   if(lx(j) >= lxi) jp(ndim)=j
      if(ndim == 0) goto 80
      call defrr(ob,     ndim*(lp1+1)*(lp2+1))
      call defrr(os,     ndim*ndim)
c ------ make the radial integrals for overlap and rhs -----
      call defrr(oxi,    nr*ndim)
      call oncint(ri,rsm,w(or),w(ow),nr,ndim,jp,lxi,ex,w(oxi),
     .   lp1,lp2,w(ophi1),w(ophi2),w(os),w(ob))
      call rlse(oxi)
c ------ solve the least-squares problem -------------
      call defrr(owk,    ndim)
      call chlr2f(w(os),w(owk),ndim,ndef)
      if(ndef /= ndim) call rx('hyfmk0:  matrix not pos definite''')
      call chlr2s(w(os),w(ob),ndim,(lp1+1)*(lp2+1))
c ------ copy to cof -------------------
      call hk0cop(ndim,jp,lxi,lp1,lp2,lmax,nx,w(ob),w(ocof))
      call rlse(ob)
  80  continue
c ------ calculate fit errors ---------------
      call oncerr(ep1,ep2,lp1,lp2,nx,ex,lmax,w(ocof),r,rsm,
     .   w(or),w(ow),nr,err)
c ------ output on file ifi -----------------
      i=0
      z=0.d0
      nalf=1
      lsym=3
      call hyfout(r ,r ,rsm ,rsm ,z ,z   ,z ,i    ,nalf,lmax,
     .   ncof,nr,i ,ri ,ri ,i  ,lp1,ep1,lp2,ep2,lsym,lx ,ex ,nx ,
     .   lx ,ex ,nx ,w(ocof),err,ifi)
      call rlse(ocof)

      return
      end
c ------ sub xxfmk0: make radial parts of the phi's -----
      subroutine xxfmk0(ephi,lphi,rsm,r,nr,phi)
      implicit real*8 (a-h,p-z), integer (o)
      dimension r(nr),f(21),phi(nr,1),ph(21),ps(21)
c   this part for smooth hankels
c|    asm=1.d0/rsm
c|    do 10 ir=1,nr
c|    call hansmr(r(ir),ephi,asm,f,lphi)
c|    rl=1.d0
c|    do 10 k=1,lphi+1
c|    phi(ir,k)=f(k)*rl
c|10  rl=rl*r(ir)
c   this part for normal hankels
      do 10 ir=1,nr
      call bessl(ephi*r(ir)**2,lphi,ph,ps)
      xx=1d0/r(ir)
      do 10 k=1,lphi+1
      phi(ir,k)=ps(k)*xx
  10  xx=xx*(1d0/r(ir))
      return
      end
c ------ sub oncint -----------------------------
      subroutine oncint(ri,rsm,r,w,nr,ndim,jp,lxi,ex,xi,
     .   lp1,lp2,phi1,phi2,s,b)
      implicit real*8 (a-h,p-z), integer (o)
      dimension f(0:20),jp(1),r(nr),w(nr),xi(nr,ndim),ex(1),
     .   phi1(nr,1),phi2(nr,2),s(1),b(ndim,lp1+1,lp2+1)
      asm=1.d0/rsm
c --------- set up radial parts of the xi -----------
      do 10 i=1,ndim
      ix=jp(i)
      do 10 ir=1,nr
      call hansmr(r(ir),ex(ix),asm,f,lxi)
  10  xi(ir,i)=f(lxi)*r(ir)**lxi
c --------- make overlap matrix ---------
      is=0
      do 20 i=1,ndim
      do 20 j=1,i
      is=is+1
      sum=0.d0
      do 21 ir=1,nr
  21  sum=sum+xi(ir,i)*xi(ir,j)*w(ir)*r(ir)**2
  20  s(is)=sum
c --------- make integrals for rhs ------------
      do 30 i=1,ndim
      ix=jp(i)
      do 30 k1=1,lp1+1
      do 30 k2=1,lp2+1
      sum=0.d0
      do 31 ir=1,nr
  31  sum=sum+phi1(ir,k1)*phi2(ir,k2)*xi(ir,i)*w(ir)*r(ir)**2
  30  b(i,k1,k2)=sum
      return
      end
c ------ sub hk0cop -----------------------------
      subroutine hk0cop(ndim,jp,lxi,lp1,lp2,lmax,nx,b,cof)
      implicit real*8 (a-h,p-z), integer (o)
      dimension jp(1),b(ndim,1),cof(lmax+1,nx,1)
      do 10 i=1,ndim
      ix=jp(i)
      do 10 kkk=1,(lp1+1)*(lp2+1)
  10  cof(lxi+1,ix,kkk)=b(i,kkk)
      return
      end
