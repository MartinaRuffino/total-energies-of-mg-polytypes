      subroutine han2cg(tau1,tau2,e1,e2,r1,r2,nlm1,nlm2,b1,b2,g1,g2,
     .  cg,indxcg,jcg,cy,s,gs)
C- Gradient of integrals of products of hankel fcts on two sites,
C  heads subtracted out. Strux & gradient returned in b1,b2; g1,g2.
      implicit none
      integer indxcg(1),jcg(1),nlm1,nlm2
      double precision cy(1),cg(1),dr(3),tau1(3),tau2(3),
     .  b1(nlm1,nlm2),b2(nlm1,nlm2),g1(nlm1,nlm2,3),g2(nlm1,nlm2,3),
     .  s(nlm1,1),gs(nlm1,nlm2,3)
      integer i,i1,i2,j,l1,l2,ll,ll1,ll2,m
      double precision e,e1,e2,fac1,fac2,fuzz,r1,r2,xx,
     .  fkj1(10),fkj2(10),fkk(10),fjj(10),fjk(10)
      parameter (fuzz = 1d-8)
      do  1  m = 1, 3
    1 dr(m) = tau2(m)-tau1(m)
      l1 = ll(nlm1)
      l2 = ll(nlm2)
C --- For integrals outside spheres ---
      call wronkj(e1,e2,r1,l1,fkk,fkj1,fjk,fjj)
      call wronkj(e2,e1,r2,l2,fkk,fkj2,fjk,fjj)
C --- Case e1 == e2 ---
      if (dabs(e1-e2) < fuzz) then
        e = (e1+e2)/2
        call mstrpg(e,dr,nlm1,nlm2,nlm1,nlm2,cg,indxcg,jcg,cy,
     .     b1,b2,g1,g2)
        do  10  i1 = 1, nlm1
          ll1 = ll(i1)+1
          fac1 = fkj1(ll1)
          do  10  i2 = 1, nlm2
          ll2 = ll(i2)+1
          fac2 = fkj2(ll2)
          s(i1,i2) = b1(i1,i2)*(fac1+fac2) + b2(i1,i2)
          do  12  i = 1, 3
   12     gs(i1,i2,i) = g1(i1,i2,i)*(fac1+fac2) + g2(i1,i2,i)
   10   continue
        do  15  j = 1, nlm2
        do  15  i = 1, nlm1
   15   b2(i,j) = b1(i,j)
C --- Case e1 /= e2 ---
      else
        call mstrug(e1,dr,nlm1,nlm2,nlm1,cg,indxcg,jcg,cy,b1,g1)
        call mstrug(e2,dr,nlm1,nlm2,nlm1,cg,indxcg,jcg,cy,b2,g2)
        do  20  i1 = 1, nlm1
          ll1 = ll(i1)+1
          fac1 = fkj1(ll1) + 1d0/(e2-e1)
          do  20  i2 = 1, nlm2
          ll2 = ll(i2)+1
          fac2 = fkj2(ll2) + 1d0/(e1-e2)
          s(i1,i2) = b1(i1,i2)*fac2 + b2(i1,i2)*fac1
          do  22  i = 1, 3
   22     gs(i1,i2,i) = g1(i1,i2,i)*fac2 + g2(i1,i2,i)*fac1
   20   continue
        endif
      end
