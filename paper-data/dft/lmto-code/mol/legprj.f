      subroutine legprj(n,x,w,f,nn,nf,cof,err)
C- Projects coffs to Legendre polynomial approximating functions
C  n=number of points, nn=order of approximation, nf= number of f's.
C  Points, weights as defined by mklegw
C  Uses the orthogonality L_n L_n' = 2/(2n+1) delta(n,n')
      implicit real*8 (a-h,p-z), integer (o)
      dimension x(1),f(n,nf),w(1),cof(0:nn,nf),p(0:50)
      do  10  i = 1, nf
      do  10  j = 0, nn
   10 cof(j,i) = 0
      sw=0d0
      do 4 k=1,n
  4   sw=sw+w(k)
      if(nn > 50) call rx('legprj: incr dim p')

      do  20  k = 1, n
      call legpol(nn,x(k),p)
      do  20  j = 0, nn
      xx = p(j)*w(k)*dble(2*j+1)/sw
C|    xx = p(j)*w(k)*dble(2*j+1)/2
      do  20  i = 1, nf
   20 cof(j,i) = cof(j,i) + xx*f(k,i)

c ------- compare approximation with function ------
      err=0d0
      do 25 i=1,n
      call legpol(nn,x(i),p)
      do 25 k=1,nf
      sum=0d0
      do 21 j=0,nn
  21  sum=sum+cof(j,k)*p(j)
  25  err=dmax1(err,dabs(sum-f(i,k)))

      end
