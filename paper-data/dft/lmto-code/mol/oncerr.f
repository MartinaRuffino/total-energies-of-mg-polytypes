      subroutine oncerr(ep1,ep2,lp1,lp2,nx,ex,lmax,cof,
     .   rmt,rsm,r,w,nr,errs)
c  calculates fit errors for the one-center case
      implicit real*8 (a-h,p-z), integer (o)
      dimension f(0:20),r(nr),w(nr),ex(1),fac1(0:10),fac2(0:10),
     .   cof(0:lmax,nx,0:lp1,0:lp2),errs(0:lp1,0:lp2),
     .   erx(0:6,0:6,3),ph(0:20),ps(0:20)
      character*1 cl(0:8)
      data cl /'s','p','d','f','g','5','6','7','8'/
      call getigv(1,ipr)
      asm=1.d0/rsm
      fpi=16.d0*datan(1.d0)
      erralf=1.5d0
      lsym=0
      if(lp1 == lp2.and.dabs(ep1-ep2) < 1d-10) lsym=1
c ------ get avg values on sphere to scale errors -------
      call hrmsav(rmt,lp1,ep1,fac1)
      call hrmsav(rmt,lp2,ep2,fac2)
      if(ipr >= 40) write(6,330) 1,(fac1(l),l=0,lp1)
      if(ipr >= 40) write(6,330) 2,(fac2(l),l=0,lp2)
  330 format(' rms avg of phi',i1,' on sphere: ',6f9.4)
c --------- start loop over pairs phi1,phi2 ------
      if(ipr >= 40) write(6,992)
  992   format(/'  case      l1      l2',
     .    '    max abs val    lm    rms-err  max-err (%)')
      do 30 l1=0,lp1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 30 l2=0,ltop
      erx(l1,l2,1)=0.d0
      erx(l1,l2,2)=0.d0
      lm1=iabs(l1-l2)
      lm2=l1+l2
      do 30 lm=lm1,lm2,2
      errmax=0.d0
      errrms=0.d0
      errvol=0.d0
      valmax=0.d0
      do 35 ir=1,nr
      rr=r(ir)
c  this part for smooth hankels
c|    call hansmr(rr,ep1,asm,f,l1)
c|    prd=f(l1)
c|    call hansmr(rr,ep2,asm,f,l2)
c|    prd=prd*f(l2)* rr**(l1+l2)
c  this part normal hankels
      call bessl(ep1*rr*rr,l1,ph,f)
      prd=f(l1)
      call bessl(ep2*rr*rr,l2,ph,f)
      prd=prd*f(l2)/rr**(l1+l2+2)
      fit=0.d0
      do 10 ix=1,nx
      call hansmr(rr,ex(ix),asm,f,lm)
  10  fit=fit+cof(lm,ix,l1,l2)*f(lm)* rr**lm
      errmax=dmax1(errmax,dabs(fit-prd))
      valmax=dmax1(valmax,dabs(prd))
      if(rr <= erralf*rmt) errvol=errvol+w(ir)*rr*rr
      if(rr <= erralf*rmt) errrms=errrms+w(ir)*rr*rr*(fit-prd)**2
  35  continue
      fac=1.d0/(fac1(l1)*fac2(l2)*fpi)
      errmax=errmax*fac
      valmax=valmax*fac
      errrms=dsqrt(errrms/errvol)*fac
      if(ipr >= 50.and.lm == lm1) write(6,991) cl(l1),cl(l2),
     .   l1,l2,valmax,lm,errrms*100,errmax*100
  991 format(1x,a1,' * ',a1,2x,i6,2x,i6,f13.4,i8,1x,2f9.2)
      if(ipr >= 50.and.lm > lm1)
     .   write(6,993) lm,errrms*100,errmax*100
  993 format(37x,i6,1x,2f9.2)
      erx(l1,l2,1)=dmax1(erx(l1,l2,1),errrms)
      erx(l1,l2,2)=dmax1(erx(l1,l2,2),errmax)
  30  erx(l1,l2,3)=valmax
c ------ put into errs, output small table ----------
      if(ipr >= 20) write(6,671)
CL      write(71,671)
      do 65 l1=0,lp1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 65 l2=0,ltop
      px=erx(l1,l2,3)
      if(ipr >= 20)
     . write(6,670) cl(l1),cl(l2),px,100*erx(l1,l2,1),100*erx(l1,l2,2)
CL      write(71,670) cl(l1),cl(l2),px,100*erx(l1,l2,1),100*erx(l1,l2,2)
      errs(l1,l2)=erx(l1,l2,1)*100
      if(lsym == 1) errs(l2,l1)=errs(l1,l2)
  65  continue
  670 format(3x,a1,' * ',a1,f13.4,1x,2f11.2)
  671 format(/' oncchk:'/'    case',
     .   '    max abs val    rms-err   max-err (%)')

      return
      end
