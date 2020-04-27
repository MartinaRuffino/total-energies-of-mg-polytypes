      subroutine move(it,nit0,lrx,gamrx,ldyn,tau,fric,temp,nbas,ips,
     .   z,amass,f,pos,fsav,psav)
C- Moves atoms. lrx=1: relax along forces. ldyn=1: dynamics.
      implicit real*8 (a-h,p-z), integer (o)
      dimension ips(1),z(1),amass(1),f(3,nbas),pos(3,nbas),
     .   fsav(3,nbas),psav(3,nbas)
      if(lrx /= 0.and.ldyn /= 0) call rx('move: lrx and ldyn nonzero')
c ------- return if it <= nit0 ------
      if(it <= nit0) then
        write(71,711) it,nit0
  711   format(' px  it',i4,'   mv=0   nit0=',i3)
        return
      endif
c ------- relax along forces -------
      if(lrx == 1) then
      do 10 ib=1,nbas
      do 11 m=1,3
      psav(m,ib)=pos(m,ib)
      fsav(m,ib)=f(m,ib)
  11  pos(m,ib)=pos(m,ib)+gamrx*f(m,ib)
      write(6,101) ib,(psav(m,ib),m=1,3),(pos(m,ib),m=1,3)
  101 format(i5,'  old',3f9.3,'   new',3f9.3)
  10  continue
      endif
c ------- dynamics -----------------
      if(ldyn == 1) then
      write(6,955) tau
  955 format(' tau=',f12.6)
c ... first step: approx v=0
      if(it == 1) then
        do 20 ib=1,nbas
        is=ips(ib)
        do 20 m=1,3
        psav(m,ib)=pos(m,ib)
        fsav(m,ib)=f(m,ib)
  20    pos(m,ib)=pos(m,ib)+f(m,ib)*(tau*tau/amass(is))
      endif
c ... subsequent steps
      if(it >= 2) then
        do 21 ib=1,nbas
        is=ips(ib)
        do 21 m=1,3
        pnew=2d0*pos(m,ib)-psav(m,ib)+f(m,ib)*(tau*tau/amass(is))
      write(6,926) f(m,ib),psav(m,ib),pos(m,ib),pnew
  926 format(' f,psav,pos,pnew=',4f12.6)
        psav(m,ib)=pos(m,ib)
        pos(m,ib)=pnew
  21    fsav(m,ib)=f(m,ib)
      endif
      endif
c ------- printout and log ---------
      shmx=0d0
      ftop=0d0
      dbot=1d10
      do 30 ib=1,nbas
      call dpdist(pos(1,ib),psav(1,ib),3,shift)
      shmx=dmax1(shmx,shift)
      ff=dsqrt(f(1,ib)**2+f(2,ib)**2+f(3,ib)**2)
      ftop=dmax1(ftop,ff)
      do 30 jb=ib+1,nbas
      call dpdist(pos(1,ib),pos(1,jb),3,dd)
      dbot=dmin1(dbot,dd)
  30  continue
      if(lrx == 1) then
      write(6,220) gamrx,ftop,shmx,dbot
  220 format(/' move:  relax with gam=',f5.2,
     .   '   ftop,shmx,dbot=',2f8.4,f9.4)
      write(71,710) it,gamrx,ftop,shmx,dbot
  710 format(' px  it',i4,'   mv=rx  gam',f7.3,'  fx',f8.4,
     .   '  shx',f8.4,'  dmn',f9.4)
      endif
      end
