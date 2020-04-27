c --------- sub tlstat ---------------------
      subroutine tlstat(nla,ndim,nel,b,t,d)
c  adds together tails on all sites, by different lmto energies
      implicit real*8 (a-h,p-z), integer(o)
      dimension ndim(1),b(nla,1),t(1),d(nla,nel)
      ne=0
      do 20 ie=1,nel
      do 18 k=1,nla
  18  d(k,ie)=0d0
      do 19 i=1,ndim(ie)
      do 19 k=1,nla
  19  d(k,ie)=d(k,ie)+b(k,i+ne)*t(i+ne)
  20  ne=ne+ndim(ie)
      end
c --------- sub hhadqj ---------------------------------
      subroutine hhadqj(wgt,nlma,nla,nel,d,njj,qjjs)
      implicit real*8 (a-h,p-z)
      dimension d(nla,nel),qjjs(njj,1)
C|    do 30 ilm=1,nlma
C|30  write(6,300) ilm,(d(ilm,ie),ie=1,nel)
C|300 format(i6,'   d=',3f12.6)
      iv=0
      do 10 ie=1,nel
      do 10 je=ie,nel
      iv=iv+1
      www=wgt
      if(ie /= je) www=2d0*wgt
      iq=0
      do 12 jlm=1,nlma
      do 12 ilm=1,nlma
      iq=iq+1
  12  qjjs(iq,iv)=qjjs(iq,iv)-www*d(ilm,ie)*d(jlm,je)
  10  continue
      end
c --------- sub hhadj0 ---------------------------------
      subroutine hhadj0(wgt,lmaxa,nla,nel,d,n0,nv,qjj0)
      implicit real*8 (a-h,p-z)
      dimension d(nla,nel),qjj0(nv,0:1)
      iv=0
      do 10 ie=1,nel
      do 10 je=ie,nel
      iv=iv+1
      www=wgt
      if(ie /= je) www=wgt*2d0
      i=0
      qjj=0d0
      do 12 l=0,lmaxa
      do 12 m=-l,l
      i=i+1
  12  qjj0(iv,l)=qjj0(iv,l)-www*d(i,ie)*d(i,je)
  10  continue
      end
