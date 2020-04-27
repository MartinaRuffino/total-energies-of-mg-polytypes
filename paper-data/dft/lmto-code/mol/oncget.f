      subroutine oncget(r,e1,nlm1,e2,nlm2,x,nf,nxi,lxi,exi,
     .   tspec,tdata,cg,jcg,indxcg)
c  Searches through tspec to find the correct table, then gets out
c  the tabulated 1-center expansion. Expansion coefficients for all
c  nlm1*nlm2 products are returned in x.
c  jxi connects table and calculation xi's: for i running over the
c  table xi energies i:1...nx, j=jxi(i) lies between 1 and nxi
c  and exi(j) equals ex(i).
      implicit real*8 (a-h,p-z), integer (o)
      dimension tspec(100,1),lxi(1),exi(1),cg(1),jcg(1),indxcg(1),
     .   jxi(10),tdata(1),lx(10),x(nf,nlm1,nlm2)
      real w(1)
      common /w/ w
      call getpr(ipr)
      l1=ll(nlm1)
      l2=ll(nlm2)
c ------- locate the correct table -----------------
      call oncloc(r,l1,e1,l2,e2,nxi,lxi,exi,tspec,itbl,jxi)
      if(itbl == 0) call rx('oncget:  no suitable table found''')
      it=iabs(itbl)
      ndt0=idnint(tspec(25,it))
      if(ipr >= 80) write(6,650) r,l1,l2,e1,e2,itbl,ndt0
  650 format(' oncget:   r=',f8.3,'   l=',2i3,'   e=',2f9.4,
     .   /' itbl=',i3,'    data offset=',i10)
c ------- get specifications of table ---------
      lmax=idnint(tspec( 8,it))
      ncof=idnint(tspec( 9,it))
      lsym=idnint(tspec(13,it))
      lp1 =idnint(tspec(15,it))
      lp2 =idnint(tspec(17,it))
      nx  =idnint(tspec(19,it))
      do 14 i=1,nx
  14  lx(i)=idnint(tspec(40+i,it))
c ------- call oncxpn ------------
      mlm1=(lp1+1)**2
      mlm2=(lp2+1)**2
      mf=0
      do 17 k=1,nx
  17  mf=mf+(lx(k)+1)**2
      call defrr(oa,    mf*mlm1*mlm2)
      call oncxpn(lx,lmax,nx,mf,lp1,lp2,mlm1,mlm2,tdata(ndt0),
     .   w(oa),cg,jcg,indxcg)
c ------- copy over into x1,x2 -----------------
      call xxcget(w(oa),mf,mlm1,mlm2,x,nf,nlm1,nlm2,
     .   nxi,lxi,nx,lx,jxi,itbl)
      call rlse(oa)

      return
      end
c ------- sub xxcget --------------------------
      subroutine xxcget(a,mf,mlm1,mlm2,x,nf,nlm1,nlm2,
     .   nxi,lxi,nx,lx,jxi,itbl)
c  copies the final 1-c-f over to array x.
      implicit real*8 (a-h,p-z), integer (o)
      dimension x(nf,nlm1,nlm2),a(mf,mlm1,mlm2),
     .    jxi(10),lxi(1),lx(1),j0(10)
      call dpzero(x,nf*nlm1*nlm2)
c --------- set up offsets for x ----------
      j0(1)=0
      do 1 j=2,nxi
  1   j0(j)=j0(j-1)+(lxi(j-1)+1)**2
c --------- this part if no switching needed ---------
      if(itbl > 0) then
      do 10 ilm1=1,nlm1
      do 10 ilm2=1,nlm2
      ii0=0
      do 15 i=1,nx
      j=jxi(i)
      jj0=j0(j)
      nlma=(lx(i)+1)**2
      nlmb=(lxi(j)+1)**2
      nlmi=min0(nlma,nlmb)
      do 16 ilm=1,nlmi
  16  x(jj0+ilm,ilm1,ilm2)=a(ii0+ilm,ilm1,ilm2)
  15  ii0=ii0+nlma
  10  continue
      endif
c --------- this part to switch phi1,phi2 ---------
      if(itbl < 0) then
      do 20 ilm1=1,nlm1
      do 20 ilm2=1,nlm2
      ii0=0
      do 25 i=1,nx
      j=jxi(i)
      jj0=j0(j)
      nlma=(lx(i)+1)**2
      nlmb=(lxi(j)+1)**2
      nlmi=min0(nlma,nlmb)
      do 26 ilm=1,nlmi
  26  x(jj0+ilm,ilm1,ilm2)=a(ii0+ilm,ilm2,ilm1)
  25  ii0=ii0+nlma
  20  continue
      endif
c --------- output -------------
c|99  ilm1=99
c|    write(6,*) 'enter ilm1,ilm2'
c|    read (5,*)  ilm1,ilm2
c|    if(ilm1 == 99) return
c|    do 80 i=1,nf
c|    xxx=x(i,ilm1,ilm2)
c|80  if(dabs(xxx) > 1.d-8) write(6,800) i,xxx
c|800 format(i6,f14.6)
c|    if(.true.) goto 99

      return
      end
