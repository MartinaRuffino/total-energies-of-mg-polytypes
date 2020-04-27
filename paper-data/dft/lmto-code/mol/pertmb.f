      subroutine pertmb(puu,pus,pss,hv,nla,d0,h1,sh1,h2,sh2,
     .  b1,sb1,b2,sb2,nlm1,nlm2,nlma,pkk,pkj,pjk,pjj)
C-  Multiplies puu,pus,pss right and left with vals and slopes
c   and adds diagonal part to make one perturbation matrix.
      implicit none
      integer nlm1,nlm2,nlma,i1,i2,nla
      double precision h1(nlm1),sh1(nlm1),h2(nlm2),sh2(nlm2),d0,
     .  b1(nlm1),sb1(nlm1),b2(nlm2),sb2(nlm2),
     .  puu(nlma,nlma),pus(nlma,nlma),pss(nlma,nlma),hv(nla,4),
     .  pkk(nlm1,nlm2),pkj(nlm1,nlma),pjk(nlma,nlm2),pjj(nlma,nlma)

C --- pkk ---
      do 10 i2=1,nlm2
      do 10 i1=1,nlm1
      pkk(i1,i2)=h1(i1)*( puu(i1,i2)*h2(i2) + pus(i1,i2)*sh2(i2) )
     .         +sh1(i1)*( pus(i2,i1)*h2(i2) + pss(i1,i2)*sh2(i2) )
   10 continue
      do 11 i1=1,min0(nlm1,nlm2)
   11 pkk(i1,i1) = pkk(i1,i1) + hv(i1,1)

C --- pkj ---
      do 20 i2=1,nlma
      do 20 i1=1,nlm1
      pkj(i1,i2)=h1(i1)*( puu(i1,i2)*b2(i2) + pus(i1,i2)*sb2(i2) )
     .         +sh1(i1)*( pus(i2,i1)*b2(i2) + pss(i1,i2)*sb2(i2) )
   20 continue
      do 21 i1=1,nlm1
   21 pkj(i1,i1) = pkj(i1,i1) + hv(i1,2) + d0

C --- pjk ---
      do 30 i2=1,nlm2
      do 30 i1=1,nlma
      pjk(i1,i2)=b1(i1)*( puu(i1,i2)*h2(i2) + pus(i1,i2)*sh2(i2) )
     .         +sb1(i1)*( pus(i2,i1)*h2(i2) + pss(i1,i2)*sh2(i2) )
   30 continue
      do 31 i1=1,nlm2
   31 pjk(i1,i1) = pjk(i1,i1) + hv(i1,3) - d0

C --- pjj ---
      do 40 i2=1,nlma
      do 40 i1=1,nlma
      pjj(i1,i2)=b1(i1)*( puu(i1,i2)*b2(i2) + pus(i1,i2)*sb2(i2) )
     .         +sb1(i1)*( pus(i2,i1)*b2(i2) + pss(i1,i2)*sb2(i2) )
   40 continue
      do 41 i1=1,nlma
   41 pjj(i1,i1) = pjj(i1,i1) + hv(i1,4)

C       print *, 'pertmb keep jk,kj chuck kk,jj'
c       print *, 'pertmb keep kk'
c       call dpzero(pkk,nlm1*nlm2)
C       call dpzero(pkj,nlm1*nlma)
C       call dpzero(pjk,nlma*nlm2)
C       call dpzero(pjj,nlma*nlma)

      end
