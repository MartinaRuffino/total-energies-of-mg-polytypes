      subroutine qifix1(n,si,six,wk,res)
C-  Make (six-si)*inv(si)*(six-si).  res can be same as si
      implicit real*8 (a-h,p-z), integer (o)
      dimension si(n,n),six(n,n),wk(n,n),res(n,n)
      real w(1)
      common /w/ w
c ----- make wk=inv(si)*(six-si) -------
      do 10 j=1,n
      do 10 i=1,n
      six(i,j)=six(i,j)-si(i,j)
  10  wk(i,j)=six(i,j)
      call defrr(opiv,  n)
      call defrr(oww,   n)
      call leqr2f(si,n,n,w(opiv),w(oww))
      call leqr2s(si,n,n,wk,n,w(opiv))
      call rlse(opiv)
c -------- put (six-si)*wk into res ----------
      do 11 i=1,n
      do 11 j=1,n
      sum=0d0
      do 12 k=1,n
   12 sum=sum+six(i,k)*wk(k,j)
   11 res(i,j) = sum

      end
