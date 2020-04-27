      subroutine fixcor(n,nstate,t,fix,evl)
c  undo effect of hamiltonian fix on eigenvalues
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(n,n),fix(n,n),evl(n)
      do 10 is=1,nstate
      cor=0d0
      do 11 i=1,n
      do 11 j=1,n
  11  cor=cor-t(i,is)*fix(i,j)*t(j,is)
      write(6,991) is,cor,evl(is),evl(is)+cor
  991 format(i5,'   corr,evl,evl+corr=',3f11.5)
      evl(is)=evl(is)+cor
  10  continue
      end
