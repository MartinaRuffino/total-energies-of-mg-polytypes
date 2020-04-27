      subroutine pertmw(nlm1,nlm2,f1,f2,v,nlmf,job,
     .  sum,ndim,nr,cg,jcg,indxcg,p)
C- Makes nlm1*nlm2 pertubation integrals
C  given tabulated functions f1 and f2.
C  v must be scaled by radial weights, ie sum_r v(r) = integral v
C  job=1 for symmetric matrix.  Work array sum saves radial integrals.
      implicit none
      integer jcg(1),indxcg(1),nr,job,ndim,nlm1,nlm2,nlmf
      double precision p(nlm1,nlm2),f1(nr,1),f2(nr,1),
     .  sum(nlmf,ndim,ndim),v(nr,1),cg(1),psum,sam
      integer i,icg,ilm1,ilm2,ilm2tp,ix,j,l1p1,l2p1,ll,lmax1,lmax2,
     .  mlm,nint,ir

C --- Large value of sum flags it has yet to be generated ---
      lmax1 = ll(nlm1)
      lmax2 = ll(nlm2)
      do  3  l1p1 = 1, lmax1+1
      do  3  l2p1 = 1, lmax2+1
      do  3  mlm = 1, nlmf
    3 sum(mlm,l1p1,l2p1) = 1.1d10

C --- Integrals by products of V and CG expansion of f1*f2 ---
      nint = 0
      do  10  ilm1 = 1, nlm1
      l1p1 = ll(ilm1)+1
      ilm2tp = nlm2
      if (job == 1) ilm2tp = ilm1
      do  10  ilm2 = 1, ilm2tp
      l2p1 = ll(ilm2)+1
      ix = max0(ilm1,ilm2)
      ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
      psum = 0d0
      do  11  icg = indxcg(ix), indxcg(ix+1)-1
      mlm = jcg(icg)
      if (mlm > 1 .and. mlm <= nlmf) then
         sam = sum(mlm,l1p1,l2p1)
         if (sam > 1d10) then
           nint = nint+1
           sam = 0d0
           do  20  ir = 1, nr
   20      sam = sam + f1(ir,l1p1)*f2(ir,l2p1)*v(ir,mlm)
           sum(mlm,l1p1,l2p1) = sam
         endif
         psum = psum + cg(icg)*sam
      endif
   11 continue
      p(ilm1,ilm2) = psum
      if (job == 1) p(ilm2,ilm1) = psum
   10 continue

C --- Printout ---
C      write(6,842) nlm1,nlm2,job,nint
C  842 format(' pertmw   nlm1,nlm2=',2i3,'  job=',i2,'   nint=',i5)
C      do  40  i = 1, nlm1
C      do  40  j = 1, nlm2
C   40 if (dabs(p(i,j)) > 1.d-5) write(6,400) i,j,p(i,j)
C  400 format(2i5,f12.6)

C      do  30  i = 1, nlm1
C   30 print 333, (p(i,j), j=1,nlm2)
C  333 format(9f8.5)


      end
