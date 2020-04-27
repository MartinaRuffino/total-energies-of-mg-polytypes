      subroutine asxdev(iq,sigma,ldim,evec,ecorr)
C- SX correction to eigenvalues, 1st order perturbation theory
Cr Routine needs to be optimized by breaking into
Cr separate multiplications.
      implicit none
      integer iq,ldim
      double precision sigma(ldim,ldim,2),ecorr(ldim,2),
     .  evec(ldim,ldim,2,*)
      integer nu,j,i,iprint,lgunit
      double complex c1,c2,c3

C --- de_nu = sum z+_i(nu) sig_ij z_j(nu) ---
      do  10  nu = 1, ldim
        ecorr(nu,1) = 0d0
        ecorr(nu,2) = 0d0
        do  20  j = 1, ldim
          c2 = dcmplx(evec(j,nu,1,iq),evec(j,nu,2,iq))
          do  30  i = 1, ldim
            c1 = dcmplx(evec(i,nu,1,iq),-evec(i,nu,2,iq))
            c3 = dcmplx(sigma(i,j,1),sigma(i,j,2))*c1*c2
            ecorr(nu,1) = ecorr(nu,1) - dble(c3)
            ecorr(nu,2) = ecorr(nu,2) - dimag(c3)
   30     continue
   20   continue
   10 continue

      if (iprint() >= 30) then
        write(*,'(9f9.4)') (ecorr(i,1),i=1,ldim)
        j = lgunit(2)
        write(j,'(9f9.4)') (ecorr(i,1),i=1,ldim)
        write(*,*) ' '
      endif
      if (iprint() >= 80) then
        call yprm('ecorr',2,ecorr,ldim,ldim,ldim,1)
      endif
      end
