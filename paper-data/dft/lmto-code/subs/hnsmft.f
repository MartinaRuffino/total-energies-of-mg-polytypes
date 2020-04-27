      subroutine hnsmft(rofi,rho,nr,qout,a,nrmt,e,qcst,xi,f,nf,err)
C- Fit the charge density tail to nf functions
C-----------------------------------------------------------------------
Ci Inputs
Ci   rofi  :radial mesh points
Ci   rho   :spherical valence charge density times 4*pi*r*r
Ci   nr    :number of radial mesh points
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nrmt  :number of points between 0..rmt
Ci   e     :vector of fitting function energies (printout only)
Ci   qcst  :vector of integrals of fitting functions from rmt...infty
Ci   xi    :vector of functions
Ci   nf    :number of functions to fit
Co Outputs
Co   qout  :true charge in tail outside rmt
Co   f     :fit coefficients (see remarks)
Co   err   :rms error in fit
Cr Remarks
Cr   The input charge density is defined on the mesh of points
Cr   1..nr, though only points nrmt..nr are included in the fitting.
Cr
Cr   rho is fit to the form rho(r) = sum_i=1,nf f_i H_i
Cr   where H_i is input as a function tabulated on the rofi mesh
Cr
Cr   Fits to the constraint of total charge conservation.
Cu Updates
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
C-----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nr,nrmt,nf
      double precision a,qout,err,rofi(nr),rho(nr),e(nf),f(nf+1),
     .  xi(nr,nf),qcst(nf)
C Local variables
      integer i,ip,ipr,ir,ire,irs,j,nfmax,info,lpr,stdo,intopt,nglob
      parameter (nfmax=10)
      integer ipiv(nfmax)
      double precision r,rhot,rhotol,dif,fpi,fit,qtf,qtr,qof,qor
      double precision s(nfmax,nfmax),wk(2*nfmax),wgt,wt(nr),ddot,dot3,
     .  sqfpi,wtr2(nr)

C --- Setup ---
      if (nf+1 > nfmax) call rx('HNSMFT: nf gt nfmax')
      fpi = 16*datan(1d0)
      sqfpi = dsqrt(fpi)
      rhotol = 1d-12
      call getpr(ipr)
      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')

C --- Find final mesh points for integrations of the tail density ---
      irs = nrmt
      ire = min(irs+3,nr)
C     Backward compatibility with 7.7
C     do  i = irs+5, nr-2, 2
C       ire = i
C       if (dabs(rho(i)) < rhotol) exit
C     enddo
      do  i = irs+3, nr
        ire = i
        if (dabs(rho(i)) < rhotol/2) exit
      enddo
      if (dabs(rho(ire)) > rhotol)
     .  call info0(10,0,0,' hnsmft (warning): rho(nr) > rhotol')

C     print *, 'ire',ire

C --- Calculate charge outside the MT sphere ---
      call radmwtp(intopt,rofi,a,irs,ire,wt)
      call radmwtp(intopt+1,rofi,a,irs,ire,wtr2)
      qout = ddot(ire-irs+1,rho(irs),1,wt(irs),1)

C --- Calculate F-integrals (rho,H_i) and S-integrals (H_i,H_j) ---
      call dpzero(s,nfmax**2)
      do  i = 1, nf
C       call fint(xi(nrmt,i),nrmt,rho,a,b,rofi,irs,ire,f(i))
        f(i) = dot3(ire-nrmt+1,xi(nrmt,i),rho(nrmt),wt(nrmt))/sqfpi
C       call awrit2('hnsmft f(%i) = %g',' ',80,stdo,i,f(i))
        do  j = 1, i
C         call sint(xi(nrmt,i),xi(nrmt,j),nrmt,a,b,rofi,irs,ire,s(i,j))
          s(i,j) = dot3(ire-nrmt+1,xi(nrmt,i),xi(nrmt,j),wtr2(nrmt))
          s(j,i) = s(i,j)
C         call awrit3('hnsmft s(%i,%i) = %g',' ',80,stdo,i,j,s(i,j))
        enddo
      enddo

C --- Add constraint that charge outsite rmt matches rho ---
      call dcopy(nf,qcst,1,s(1,nf+1),1)
      call dcopy(nf,qcst,1,s(nf+1,1),nfmax)
      s(nf+1,nf+1) = 0
      f(nf+1) = qout

C --- Least-squares fit ---
C     call prmx('s',s,nfmax,nf+1,nf+1)
C     call prmx('f',f,nfmax,nf,1)
      call dsytrf('L',nf+1,s,nfmax,ipiv,wk,2*nfmax,info)
      if (info /= 0) call rx('hnsmft: normal matrix not pos def')
      call dsytrs('L',nf+1,1,s,nfmax,ipiv,f,nf+1,info)

C --- Printout ---
      if (ipr > 30) then
        write(stdo,100) ire-irs+1,rofi(irs),rofi(ire),qout
        call awrit2(' E:%n:4;8F',' ',80,stdo,nf,e)
        call awrit2(' C:%n:4;8F',' ',80,stdo,nf,f)
      endif
  100 format(' HNSMFT:',i4,' points in interval',f9.5,f10.5,
     .  ';  q=',f10.6)
C     write(stdo,101) 'E',(e(ip), ip=1,nf)
C     write(stdo,101) 'C',(f(ip), ip=1,nf)
C 101 format(5x,a1,':',10f12.5)

      err = 0
      qof  = 0
      qor  = 0
      qtf  = 0
      qtr  = 0
      call radwgt(intopt+1,rofi(ire),a,ire,wt)
      if (ipr > 30) write (stdo,1)
    1 format(7x,' r ',9x,'rho',9x,'fit',9x,'diff')
      do  ir = 2, ire
        r = rofi(ir)
        rhot = rho(ir)/fpi/r/r
        fit = 0
        do  ip = 1, nf
          fit = fit + f(ip)*xi(ir,ip)
        enddo
        fit = fit/dsqrt(fpi)
        dif = rhot - fit
C       wgt  = (mod(ir+1,2)+1)*fpi*r*r*(r+b)*2*a/3d0
C       print *, wgt- wt(ir)*fpi,wtr2(ir)
        wgt = wt(ir)*fpi
        qtf = qtf + wgt*fit
        qtr = qtr + wgt*rhot
        if (ir >= irs) then
          err = err + wtr2(ir)*fpi*dif**2
          qof = qof + wtr2(ir)*fpi*fit
          qor = qor + wtr2(ir)*fpi*rhot
        endif
C       rho(ir) = fpi*r*r*fit
        if ( .not. (ipr <= 30.or.ir < irs.and.ipr < 50)) then
        lpr = 0
        if (ipr >= 80 .or. ir == irs) lpr = 1
        if (ipr >= 40 .and. ir == 2) lpr = 1
        if (dabs(rhot) > 1d-6  .and. (mod(ir-irs,10) == 0)) lpr = 1
C          if (lpr >= 1) write (stdo,2) r,rhot,fit,dif
C    2     format(4F12.6)
          if (lpr >= 1) call info5(10,0,0,'%;12,6D%;12,6D%;12,6D%;12,6D',
     .      r,rhot,fit,dif,0)
        endif
      enddo

      err = err/(rofi(ire)-rofi(irs))
      err = dsqrt(err)
      if (ipr >= 30) then
        write (stdo,3) qof,err
    3   format('    q(fit):',f13.6,4x,'rms diff:',f11.6)
        write (stdo,4) qof,qtf-qof,qtf,qor,qtr-qor,qtr
    4   format(4x,'fit: r>rmt',f10.6,'   r<rmt',f10.6,'   qtot',f10.6,
     .        /4x,'rho: r>rmt',f10.6,'   r<rmt',f10.6,'   qtot',f10.6)
      endif

      end

C      subroutine sint(xi,xj,nrmt,a,b,rofi,irs,ire,sum)
C      implicit none
C      integer ire,irs,ir,nrmt,jr
C      double precision rofi(ire),r,a,b,sum,xi(*),xj(*)
C      sum = 0d0
C      do  ir = irs+1, ire-1
C        jr = ir+1-nrmt
C        r = rofi(ir)
C        sum = sum + (mod(ir+1,2)+1)*(r+b)*r**2*xi(jr)*xj(jr)
C      enddo
C      r = rofi(irs)
C      sum = sum + .5d0*(r+b)*r**2*xi(irs+1-nrmt)*xj(irs+1-nrmt)
C      r = rofi(ire)
C      sum = sum + .5d0*(r+b)*r**2*xi(ire+1-nrmt)*xj(ire+1-nrmt)
C      sum = 2d0*sum*a/3d0
C      end
C      subroutine fint(xi,nrmt,rho,a,b,rofi,irs,ire,sum)
C      implicit none
C      integer ire,irs,ir,nrmt
C      double precision rho(*),rofi(ire),r,a,b,sqfpi,sum,xi(*)
C      sum = 0d0
C      sqfpi = dsqrt(16*datan(1d0))
C      do  ir = irs+1, ire-1
C        r = rofi(ir)
C        sum = sum + (mod(ir+1,2)+1)*rho(ir)*(r+b)*xi(ir+1-nrmt)
C      enddo
C      r = rofi(irs)
C      sum = sum + .5d0*rho(irs)*(r+b)*xi(irs+1-nrmt)
C      r = rofi(ire)
C      sum = sum + .5d0*rho(ire)*(r+b)*xi(ire+1-nrmt)
C      sum = 2d0*sum*a/3d0/sqfpi
C      end
