      subroutine perco1(n,nstate,nst1,t,fix,evl)
C- Undo hamiltonian fix by pertubation
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(n,n),fix(n,n),evl(n)
      call tcn('perco1')
      call getpr(ipr)
      if(ipr >= 45) write(6,992)
      demax=0d0
      deavg=0d0
      sde=0d0
      do 10 is=1,nstate
      cor=0d0
      do 11 j=1,n
      sam=0d0
      do 12 i=1,n
  12  sam=sam+t(i,is)*fix(i,j)
  11  cor=cor+sam*t(j,is)
      de=-cor
      if(ipr >= 45) write(6,991) is,de,evl(is),evl(is)+de
  991 format(i5,5f13.6)
  992 format(/' istate     de',11x,'evl',9x,'evl+de')
      evl(is)=evl(is)+de
      if(is <= nst1) demax=dmax1(demax,dabs(de))
      if(is <= nst1) sde=sde+de
  10  if(is <= nst1) deavg=deavg+dabs(de)/nst1
      if(iprint() > 10)
     .call awrit3(' max h-fix=%#8,8d avg h-fix=%#8,8d sum h-fix=%#8,8d',
     .            ' ',120,i1mach(2),demax,deavg,sde)
      call tcx('perco1')
      end
