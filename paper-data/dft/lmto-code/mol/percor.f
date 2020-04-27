      subroutine percor(n,nstate,t,fix,dqi,vint,vintf,evl)
C- two pertubation corrections to eigenvalues:
C  (1) undo hamiltonian fix  (2) change from vint to vintf
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(n,n),fix(n,n),dqi(n,n),evl(n)
      call getpr(ipr)
      demax=0d0
      deavg=0d0
      cormax=0d0
      coravg=0d0
      if(ipr >= 30) write(6,992)
      do 10 is=1,nstate
      cor=0d0
      dq=0d0
      do 11 i=1,n
      do 11 j=1,n
      dq=dq+t(i,is)*dqi(i,j)*t(j,is)
  11  cor=cor-t(i,is)*fix(i,j)*t(j,is)
      de=dq*(vint-vintf)
      if(ipr >= 30) write(6,991) is,dq,de,cor,evl(is),evl(is)+cor+de
  991 format(i5,5f11.6)
  992 format(' istate   dq',9x,'de1',8x,'de2',8x,'evl',7x,'evl+corr')
      evl(is)=evl(is)+de+cor
      demax=dmax1(demax,dabs(de))
      deavg=deavg+dabs(de)/nstate
      cormax=dmax1(cormax,dabs(cor))
  10  coravg=coravg+dabs(cor)/nstate
      if(ipr >= 10) print 200, demax,deavg,deavg*nstate,
     .  cormax,coravg,coravg*nstate
  200 format(/' percor: max dq*vi=',f10.6,'  avg dq*vi=',f10.6,
     .  '  sum dq*vi=',f10.6
     .       /'         max a-fix=',f10.6,'  avg a-fix=',f10.6,
     .  '  sum a-fix=',f10.6)
CL      write(71,710) nstate,demax,deavg,cormax,coravg
  710 format(' px  nst',i4,'   dq*vi',2f10.6,'   fix',2f10.6)

      end
