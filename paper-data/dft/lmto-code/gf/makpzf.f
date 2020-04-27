      subroutine makpzf(job,avw,kbak,l,hcr,slo,val,rmax,ptf,beta)
C- Makes potential function and downfolding beta, complex energy
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job   :1s digit : number of derivatives in ptf to compute
Ci         :Must be one of 0, 2 or 4.
Ci         :10s digit :
Ci         :0 kbak is independent of energy
Ci         :1 kbak is energy-dependent (not implemented)
Ci         :100s digit :
Ci         :0 ptf is fit to a tangent form:
Ci         :      ptf = (z'-(C-znu))/(del + gam(z'-(C-znu)))
Ci         :  where
Ci                z' = z-znu
Ci         :  See outputs for storages of coefficients C-znu,del,gam
Ci         :1 ptf is fit to a tangent form as in job 0,
Ci         :  but z' = z-znu + pgam (z-znu)**3  (third order form)
Ci         :2 ptf is fit to a polynomial;
Ci            makpzf returns the coefficients to a Taylor series
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   kbak  :kinetic energy of the envelope function
Ci   l     :angular momentum
Ci   hcr   :hard core sphere radius
Ci   slo   :slope of g at rmax, and optionally higher order derivatives,
Ci         :if 1s digit of job0>0.
Ci   val   :value of g at rmax, and optionally higher order derivatives,
Ci         :if 1s digit of job0>0.
Ci   rmax  :MT radius Wigner-Seitz sphere radius, in atomic units
Co Outputs:
Co   ptf   :potential function P^0 = W{phi,K}/W{phi,J}
Co         :In 2nd gen LMTO, where kbak=0, using Andersen's old defs:
Co         :P^0 -> 2(2l+1) (avw/rmax)^(2l+1) (D+l+1) / (D-l)
Co         :and (P^alpha)^-1 = (P^0)^-1 - alpha
Co         :If higher order derivatives are computed (1s digit job>0)
Co         :data about the expansion of ptf is return in ptf(2..)
Co         :The parameterization of ptf takes one of two forms,
Co         :depending on 100s digit of job.
Co         :Tangent form:
Co         :  ptf(2) = C-znu (znu = energy at which val,slo are made)
Co         :  ptf(3) = del
Co         :  ptf(4) = gam
Co         :  ptf(5) = pgam (100s digit job=1), otherwise zero
Co         :Polynomial form:
Co         :  ptf(i) = ith energy derivative of ptf
Co   beta  :finds the specific alpha that matches to w.f. at hcr,
Co          optimal for downfolding.
Cr Remarks:
Cr   the potential function is 1/beta, where beta is obtainable
Cr   from the back extrapolated function phi^0, which has the form
Cr
Cr     phi^0 = const * (j^0(kr) - beta n^0(kr))
Cr
Cr   See Lecture Notes in Physics 535, after Eq. 34.
Cr   This was adapted from the Stuttgart third-generation LMTO package.
Cr
Cr   The standard relations of Wronskians to parameters C, gam, del
Cr   don't work for complex energy.  Here the parameters are obtained
Cr   by a least-squares procedure.
Cu Updates
Cu   18 Dec 01 Adapted from makptf.f
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer job,l
      double precision avw,hcr,rmax
      double complex kbak,beta,slo(*),val(*),ptf(*)
C Local variables:
      logical pass1
      integer nlmax,job0,i,j,k,ipiv(5),ierr,iprint,n
      parameter(nlmax=20)
*     double precision ei(5),dele
      double complex ei(5),dele,z,de1,del1,de2,del2
      double complex vali,sloi,ptfi(5)
      double complex er2,fi(0:nlmax+1),gi(0:nlmax+1),wn,wj,dphi
      double complex a(5,5),rhs(5),f(5),c,del,gam,pgam
C External calls:
      external besslz

      job0 = mod(job,10)
      if (job < 0 .or. job0/2 > 2) call rxi('makpzf: illegal job',job)
      if (job/100 > 1) call rxi('makpzf: not implemented, job',job)
      if (job/100 == 1 .and. job0 /= 4)
     .   call rxi('makpzf: illegal job',job)


C     Do the parameterizations by generating a set of finite differences
C     and then fitting parameterized form, or making polynomial fit.
      dele = .002d0

      ei(1) = 0
      ei(2) = 2
      ei(3) = -2
      ei(4) = 4
      ei(5) = -4
      if (job0 == 2) then
        ei(2) = 2
        ei(3) = -2
      endif

      do  i = 1, 1+job0
        ei(i) = ei(i)*dele

        z = ei(i)
        if (job0 == 0) then
          vali = val(1)
          sloi = slo(1)
        elseif (job0 == 2) then
          vali = (val(3)*z+val(2))*z+val(1)
          sloi = (slo(3)*z+slo(2))*z+slo(1)
        else
          vali = (((val(5)/24*z+val(4)/6)*z+val(3)/2)*z+val(2))*z+val(1)
          sloi = (((slo(5)/24*z+slo(4)/6)*z+slo(3)/2)*z+slo(2))*z+slo(1)
        endif

        er2 = kbak*rmax*rmax
        call besslz(er2,11,0,l+1,fi(0),gi(0))
        dphi = rmax*sloi /vali - 1.d0
C   ... wj = (D{phi}-D{J}), wn = (D{phi}-D{K})
        wn = (dphi-l) + gi(l+1)/gi(l)*(l+l+1)
        wj = (dphi-l) +er2*fi(l+1)/fi(l)/(l+l+1)
        ptfi(i) = (avw/rmax)**(l+l+1)*gi(l)/fi(l)*wn/wj
        if (i == 1) then
          if (hcr /= 0) then
            er2 = kbak*hcr*hcr
            call besslz(er2,11,0,l,fi(0),gi(0))
            beta = fi(l)/gi(l)*(hcr/avw)**(l+l+1)
          else
            beta = fi(l)/gi(l)
          endif
        endif
      enddo

      ptf(1) = ptfi(1)
      if (job0 /= 0) then
        de1  = (ei(2) - ei(3))/2
        del1 = (ei(2) + ei(3))/2
        de2  = (ei(4) - ei(5))/2
        del2 = (ei(4) + ei(5))/2

C        print *, '!!'
C        do  i = 1, 1+job0
C          write(66,357) ei(i),ptfi(i)
C  357     format(4f15.10)
C        enddo


C   --- Least-squares fit to tangent form ---
C       Fit is:
C         y = (a + b*z + d*z^3) (1 + c*z + d*gam*z^3)^-1
C       or
C         y = a + b*z - c*y*z + d(1-gam*y)*z*3
C       The third-order p can't be fit until gam is known ...  so
C       First pass: fit
C          y = a + b*z - c*y*z
C       => a = -C/(del-C*gam)  b = 1/(del-C*gam)  c=gam/(del-C*gam)
C       Second pass includes fit with 2nd order estimate for gam and
C       => d = pgam/(del-C*gam)
C       Then
C        gam=c/b  C=-a/b  del=C*gam + 1/b  pgam=d/b
        if (job/100 <= 1) then
          gam = 0
          pgam = 0
          pass1 = .true.
C         Re-entry point for second pass
   20     continue
          n = 4
          if (pass1) n = 3
          call dpzero(a,2*5*5)
          call dpzero(rhs,2*5)
          do  i = 1, 1+job0
            f(1) = 1
            f(2) = ei(i)
            f(3) = -ptfi(i)*ei(i)
            f(4) = (1-ptfi(i)*gam)*ei(i)**3
            do  j = 1, 4
              do  k = 1, 4
C               a(j,k) = a(j,k) + dconjg(f(j))*f(k)
                a(j,k) = a(j,k) + f(j)*f(k)
              enddo
            enddo
            do  j = 1, 4
              rhs(j) = rhs(j) + f(j)*ptfi(i)
            enddo
          enddo

C         call zprm('normal',2,a,5,n,n)
C         call zprm('rhs',2,rhs,5,n,1)
          ierr = 1
          call zgetrf(n,n,a,5,ipiv,ierr)
          call rxx(ierr /= 0,'bug in makpzf')
          call zgetrs('N',n,1,a,5,ipiv,rhs,5,ierr)
          call rxx(ierr /= 0,'bug in makpzf')

          gam = rhs(3)/rhs(2)
          c = -rhs(1)/rhs(2)
          del = c*gam + 1/rhs(2)
          if (.not. pass1) pgam = rhs(4)/rhs(2)
          ptf(2) = c
          ptf(3) = del
          ptf(4) = gam
          ptf(5) = pgam
          if (iprint() >= 100) then
            print 333, ' c=',c,'gam=',gam,'del=',del,'p=',pgam
C            print 333,'check P=',ptf(1),'-c/(del-c*gam)-P=',
C     .        -c/(del-c*gam)-ptf(1)
  333       format(4(a,2f12.7,2x))
          endif
          pass1 = .not. pass1
          if (job/100 /= 0 .and. .not. pass1) goto 20

C   --- Polynomial fit ---
        else
          call dfphiz(de1,del1,de2,del2,1,ptfi,ptfi(2),job0 == 4)
          call zcopy(job0,ptfi(2),1,ptf(2),1)
C          do  i = 1, 4
C            write(66,357) ei(i),ptfi(i+1)
C  357       format(4f15.10)
C          enddo
C          close(66)
        endif
      endif

      end

