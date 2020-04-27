      subroutine ftcsev(mmax,ndim,nsdmx,s,lpr)
C- look at evals of fit overlap matrix
      implicit real*8 (a-h,p-z), integer (o)
      dimension s(nsdmx,0:mmax),ndim(0:20),evl(800),rat(0:20)
      real w(1)
      common /w/ w
      call getpr(ipr)
      if(ipr < 20.and.lpr == 0) return
      if(ipr >= 20) write(6,500)
  500 format(/' ftcsev:  eigenvalues of overlap matrix:')

c --------- loop over m -------------------
      do 10 m=0,mmax
      nd=ndim(m)
      if(nd > 800) call rx('ftcsev: increase dim evl')
      call defrr(oa,    nd*nd)
      call defrr(os0,   nd*nd)
      call defrr(ot,    nd*nd)
      call xxftse(nd,s(1,m),w(oa),w(os0))
      call pshpr(1)
      call diagrr(nd,nd,w(os0),w(oa),w(ot),evl)
      call poppr
      ebot=evl(1)
      etop=evl(nd)
      if(ipr >= 40) write(6,300) m,nd,evl(1),evl(nd),evl(nd)/evl(1)
  300 format(9x,i2,'   dim',i4,'    min,max,rat',
     .   1p,3d11.2)
      rat(m)=evl(nd)/evl(1)
      if(ipr >= 60) write(6,301) (evl(i),i=1,nd)
  301 format(' evl',1p,5d15.5)
      call rlse(oa)
  10  continue
      if(ipr >= 20.and.ipr < 40) write(6,710) (rat(m),m=0,mmax)
CL      if(lpr == 1) write(71,710) (rat(m),m=0,mmax)
  710 format(' evn/ev1 ',1p,6d10.2)
      end

      subroutine xxftse(nd,s,a,s0)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s(1),a(nd,nd),s0(nd,nd)
      is=0
      do 11 j=1,nd
      do 11 i=1,j
      is=is+1
      a(i,j)=s(is)
  11  a(j,i)=s(is)
      do 20 i=1,nd
      do 21 j=1,nd
  21  s0(i,j)=0d0
  20  s0(i,i)=1d0
      end
