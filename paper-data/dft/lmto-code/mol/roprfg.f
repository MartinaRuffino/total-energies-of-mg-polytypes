      subroutine roprfg(m,f,g,cm,sm,rho,n)
c  adds f*cm and g*sm to density, rho
      implicit real*8 (a-h,p-z), integer (o)
      dimension f(n),g(n),cm(n),sm(n),rho(n)
      do 2 i=1,n
  2   rho(i)=rho(i)+f(i)*cm(i)
      if(m == 0) return
      do 4 i=1,n
  4   rho(i)=rho(i)+g(i)*sm(i)
      end
