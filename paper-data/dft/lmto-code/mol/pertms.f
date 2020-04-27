      subroutine pertms(ak,aj,lmxa,nlmb,nlmr,v,nlmf,rofi,a,nr,
     .   p1,p2,p3,cg,jcg,indxcg,nsp)
C- makes matrices p1,p2,p3 for pertubation treatment of nonspherical
c  potential in spheres. input= radial fcts ak,aj and pot v.
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
      implicit real*8 (a-h,p-z), integer(o)
      integer intopt,nglob
      dimension v(nr,nlmf,nsp),ak(nr,lmxa+1,nsp),aj(nr,lmxa+1,nsp),
     .  rofi(1),p1(nlmb,nlmb,nsp),p2(nlmb,nlmr,nsp),p3(nlmr,nlmr,nsp),
     .  cg(1),jcg(1),indxcg(1)

      do  10  i = 1, nsp
   10 call pertm1(ak(1,1,i),aj(1,1,i),nlmb,nlmr,v(1,1,i),nlmf,rofi,a,nr,
     .   p1(1,1,i),p2(1,1,i),p3(1,1,i),cg,jcg,indxcg)
      end
      subroutine pertm1(ak,aj,nlmb,nlmr,v,nlmf,rofi,a,nr,
     .   p1,p2,p3,cg,jcg,indxcg)
C- makes matrices p1,p2,p3 for pertubation treatment of nonspherical
c  potential in spheres. input= radial fcts ak,aj and pot v.
      implicit real*8 (a-h,p-z), integer(o)
      dimension rofi(1),v(nr,1),ak(nr,1),aj(nr,1),cg(1),jcg(1),
     .   p1(nlmb,nlmb),p2(nlmb,nlmr),p3(nlmr,nlmr),indxcg(1)
      real w(1)
      common /w/ w
      lmaxb = ll(nlmb)
      lmaxr = ll(nlmr)
      ndim = max0(lmaxb+1,lmaxr+1)
      call defrr(osum,   nlmf*ndim*ndim)
      call defrr(owt,    nr)
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rofi(nr),a,nr,w(owt))
      call xxper1(nlmb,nlmb,p1,ak,ak,1,w(osum),ndim,v,nlmf,
     .   rofi,w(owt),a,nr,cg,jcg,indxcg)
      call xxper1(nlmb,nlmr,p2,ak,aj,2,w(osum),ndim,v,nlmf,
     .   rofi,w(owt),a,nr,cg,jcg,indxcg)
      call xxper1(nlmr,nlmr,p3,aj,aj,1,w(osum),ndim,v,nlmf,
     .   rofi,w(owt),a,nr,cg,jcg,indxcg)
      call rlse(osum)
      end

c ----------- xxper1 ----------------------
      subroutine xxper1(nlm1,nlm2,p,a1,a2,job,sum,ndim,v,nlmf,
     .   rofi,h,a,nr,cg,jcg,indxcg)
C- makes a nlm1*nlm2 matrix of p pertubation integrals.
c  job=1 for symmetric matrix, job=2 for asymmetric matrix.
c  use array sum to save and re-use radial integrals.
      implicit real*8 (a-h,p-z), integer(o)
      dimension p(nlm1,nlm2),a1(nr,1),a2(nr,1),sum(nlmf,ndim,ndim),
     .   v(nr,1),h(nr),cg(1),jcg(1),indxcg(1),rofi(nr)

C      b=rofi(nr)/(dexp(a*nr-a)-1.d0)
C      do  2  ir = 1, nr
C      wgt=2.d0*(mod(ir+1,2)+1)
C      if (ir == 1.or.ir == nr) wgt=1.d0
C    2 h(ir)=a*(rofi(ir)+b)*wgt/3.d0
C     call prmx('h',h,nr,nr,1)

      lmax1 = ll(nlm1)
      lmax2 = ll(nlm2)
      do  3  l1p1 = 1, lmax1+1
      do  3  l2p1 = 1, lmax2+1
      do  3  mlm = 1, nlmf
    3 sum(mlm,l1p1,l2p1) = 1.1d10
      nint = 0
      do  10  ilm1 = 1, nlm1
      l1p1 = ll(ilm1)+1
      ilm2tp = nlm2
      if (job == 1) ilm2tp = ilm1
      do  10  ilm2 = 1, ilm2tp
      l2p1 = ll(ilm2)+1
      ix = max0(ilm1,ilm2)
      ix = (ix*(ix-1))/2+min0(ilm1,ilm2)
      psum = 0.d0
      do  11 icg = indxcg(ix),indxcg(ix+1)-1
      mlm = jcg(icg)
      if (mlm > 1.and.mlm <= nlmf) then
        sam = sum(mlm,l1p1,l2p1)
        if (sam > 1.d10) then
          nint = nint+1
          call ffvint(sam,a1(1,l1p1),a2(1,l2p1),v(1,mlm),h,nr)
          sum(mlm,l1p1,l2p1) = sam
        endif
        psum = psum+cg(icg)*sam
      endif
   11 continue
      p(ilm1,ilm2) = psum
      if (job == 1) p(ilm2,ilm1) = psum
   10 continue
c ----- printout ----------------
c|    write(6,842) nlm1,nlm2,job,nint
c|842 format(' pertm1   nlm1,nlm2=',2i3,'  job=',i2,'   nint=',i5)
c|    do  40  i = 1, nlm1
c|    do  40  j = 1, nlm2
c|40  if (dabs(p(i,j)) > 1.d-5) write(6,400) i,j,p(i,j)
c|400 format(2i5,f12.6)
      return
      end
c -------------------------------------------
      subroutine ffvint(sum,a1,a2,v,h,nr)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a1(1),a2(1),v(1),h(1)
      sam = 0.d0
      do  10  ir = 2, nr
   10 sam = sam + h(ir)*a1(ir)*a2(ir)*v(ir)
      sum = sam
      end
