      subroutine perco2(n,nstate,t,dqi,vint,vintf,evl)
C- Change from vint to vintf by pertubation
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(n,n),dqi(n,n),evl(n)
      call tcn('perco2')
      call getpr(ipr)
      if(ipr >= 35) write(6,992)
      cormax=0d0
      coravg=0d0
      sde=0d0
      do 10 is=1,nstate
      cor=0d0
      do 11 j=1,n
      sam=0d0
      do 12 i=1,n
  12  sam=sam+t(i,is)*dqi(i,j)
  11  cor=cor+sam*t(j,is)
      de=cor*(vint-vintf)
      sde=sde+de
      if(ipr >= 35) write(6,991) is,cor,de,evl(is),evl(is)+de
  991 format(i5,5f13.6)
  992 format(/' istate     dq',11x,'de',11x,'evl',9x,'evl+de')
      evl(is)=evl(is)+de
      cormax=dmax1(cormax,dabs(cor))
  10  coravg=coravg+dabs(cor)/nstate
      fac=dabs(vint-vintf)
      if(ipr >= 10) write(6,200)
     .  cormax,coravg,sde,fac*cormax,fac*coravg,fac*sde
  200 format(
     .  ' max dq=   ',f12.6,'   avg dq=   ',f12.6,'   sum dq=   ',f12.6,
     . /' max dq*vi=',f12.6,'   avg dq*vi=',f12.6,'   sum dq*vi=',f12.6)
CL      write(71,710) nstate,fac*cormax,fac*coravg,cormax,coravg
  710 format(' mce  nst',i4,'   demx,avg',2f10.6,'  dqmx,avg',2f10.6)

      call tcx('perco2')
      end
