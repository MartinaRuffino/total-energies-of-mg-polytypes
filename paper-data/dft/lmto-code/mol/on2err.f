      subroutine on2err(ep1,ep2,lp1,lp2,nx,lmax,lmaxl,cof,
     .  rmt,xie,phi1,phi2,r,wr,nr,errs)
C- Calculates fit errors for the one-center case
      implicit real*8 (a-h,p-z), integer (o)
      integer lp1,lp2,nx,lmax,nr,lmaxl
      double precision
     .  r(nr),wr(nr),fac1(0:10),fac2(0:10),
     .  cof(0:lmax,nx,0:lp1,0:lp2),errs(0:lp1,0:lp2),
     .  erx(0:6,0:6,3),
     .  xie(nr,0:lmax,nx),phi1(nr,0:1),phi2(nr,0:1)
      double precision arrmax,arrrms,arrvol,ep1,ep2,
     .  erralf,errmax,errrms,errvol,fac,fit,fpi,
     .  prd,px,rmt,rr,valmax
      integer ipr,ir,ix,l,l1,l2,lm,lm1,lm2,lsym,ltop
      character*1 cl(0:8)
      data cl /'S','P','D','F','G','5','6','7','8'/
      call getpr(ipr)
      fpi=16.d0*datan(1.d0)
      erralf=1.5d0
      lsym=0
      if(lp1 == lp2.and.dabs(ep1-ep2) < 1d-10) lsym=1
c ------ get avg values on sphere to scale errors -------
      call hrmsav(rmt,lp1,ep1,fac1)
      call hrmsav(rmt,lp2,ep2,fac2)
      if(ipr >= 40) write(6,330) 1,(fac1(l),l=0,lp1)
      if(ipr >= 40) write(6,330) 2,(fac2(l),l=0,lp2)
  330 format(' on2err: rms avg of phi',i1,' on sphere: ',6f9.4)
c ------ start loop over pairs phi1,phi2 ------
      if(ipr >= 40) write(6,992)
      if(ipr >= 40) write(6,994)
  992 format(/'  case   l1   l2','    max val   lm',
     .  '     error inside rmt    outside rmt (%)')
  994 format(39x,'rms      max        rms      max')
      do 30 l1=0,lp1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 30 l2=0,ltop
      erx(l1,l2,1)=0.d0
      erx(l1,l2,2)=0.d0
      lm1=iabs(l1-l2)
      lm2=min(l1+l2,lmax)
      do 30 lm=lm1,lm2,2
      errmax=0.d0
      errrms=0.d0
      errvol=0.d0
      arrmax=0.d0
      valmax=0.d0
      arrrms=0.d0
      arrvol=0.d0
      do 35 ir=1,nr
      rr=r(ir)
c ------ make phi1*phi2 and fit -----------
      prd=phi1(ir,l1)*phi2(ir,l2)
      fit=0.d0
      do 10 ix=1,nx
   10 fit = fit + cof(lm,ix,l1,l2)*xie(ir,lm,ix)
c      if(lm == 4) write(49,721) rr,prd,fit
c  721 format(3f12.6)
c ------ put into error arrays --------------
      valmax=dmax1(valmax,dabs(prd))
      if(rr > rmt) errmax=dmax1(errmax,dabs(fit-prd))
      if(rr <= erralf*rmt.and.rr >= rmt) then
        errvol=errvol+wr(ir)*rr*rr
        errrms=errrms+wr(ir)*rr*rr*(fit-prd)**2
        endif
      if(rr < rmt) then
        arrvol=arrvol+wr(ir)*rr*rr
        arrrms=arrrms+wr(ir)*rr*rr*(fit-prd)**2
        arrmax=dmax1(arrmax,dabs(fit-prd))
        endif
  35  continue
c ------ output long table ----------------
      fac=1.d0/(fac1(l1)*fac2(l2)*fpi)
      errmax=errmax*fac
      valmax=valmax*fac
      errrms=dsqrt(errrms/errvol)*fac
      arrmax=arrmax*fac
      arrrms=dsqrt(arrrms/dmax1(arrvol,1d0))*fac
      if (lm <= lmaxl) arrrms=0
      if (lm <= lmaxl) arrmax=0
      if(ipr >= 40.and.lm == lm1) write(6,991) cl(l1),cl(l2),l1,l2,
     .   valmax,lm,arrrms*100,arrmax*100,errrms*100,errmax*100
      if(ipr >= 40.and.lm > lm1) write(6,993) lm,arrrms*100,
     .   arrmax*100,errrms*100,errmax*100
  991 format(1x,a1,' * ',a1,2i5,f10.3,i6,1x,2f9.2,2x,2f9.2)
  993 format(26x,i6,1x,2f9.2,2x,2f9.2)
      erx(l1,l2,1)=dmax1(erx(l1,l2,1),errrms)
      erx(l1,l2,2)=dmax1(erx(l1,l2,2),errmax)
  30  erx(l1,l2,3)=valmax
c ------ put into errs, output small table ----------
      if(ipr >= 30.and.ipr < 40) write(6,671)
CL      write(71,671)
      do 65 l1=0,lp1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 65 l2=0,ltop
      px=erx(l1,l2,3)
      if(ipr >= 30.and.ipr < 40)
     . write(6,670) cl(l1),cl(l2),px,100*erx(l1,l2,1),100*erx(l1,l2,2)
CL      write(71,670) cl(l1),cl(l2),px,100*erx(l1,l2,1),100*erx(l1,l2,2)
      errs(l1,l2)=erx(l1,l2,1)*100
      if(lsym == 1) errs(l2,l1)=errs(l1,l2)
  65  continue
  670 format(3x,a1,' * ',a1,f13.4,1x,2f11.2)
  671 format(/' on2err:'/'    case',
     .   '    max abs val    rms-err   max-err (%)')

      return
      end
