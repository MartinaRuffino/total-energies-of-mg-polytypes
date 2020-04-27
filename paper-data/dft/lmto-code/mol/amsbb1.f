      subroutine amsbb1(nel,nbas,lphi,lmxa,ips,n0,ndim,nla,bb,bb1)
C- Extracts needed strux energy derivatives from full array.
      implicit real*8 (a-h,p-z), integer(o)
      dimension lmxa(1),ips(1),lphi(n0,1),ndim(1),bb(nla,1),bb1(1)
      call tcn('amsbb1')
      ne1=1
      je1=1
      do 20 ie=1,nel
      call xxsbb1(lmxa,nbas,ips,n0,lphi(ie,1),nla,ndim(ie),
     .   bb(1,ne1),bb1(je1))
      ne1=ne1+ndim(ie)
   20 je1=je1+ndim(ie)**2

      call tcx('amsbb1')
      end
c ---------- xxsbb1 -----------------------
      subroutine xxsbb1(lmxa,nbas,ips,n0,lphi,nla,nlb,bb,bb1)
      implicit real*8 (a-h,p-z), integer(o)
      dimension lmxa(1),ips(1),lphi(n0,1),bb1(nlb,nlb),bb(nla,nlb)
      if(nlb <= 0) return
      koff=0
      iof1=0
      do 36 ib=1,nbas
      is=ips(ib)
      nlmb=(lphi(1,is)+1)**2
      nlma=(lmxa(is)+1)**2
      do 37 ilm=1,nlmb
      i=iof1+ilm
      k=koff+ilm
      do 37 j=i,nlb
   37 bb1(i,j)=bb(k,j)
      iof1=iof1+nlmb
   36 koff=koff+nlma
      end
