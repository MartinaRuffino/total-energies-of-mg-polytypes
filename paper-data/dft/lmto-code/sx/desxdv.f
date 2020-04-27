      subroutine desxdv(ldim,nsp,nbas,sigma,evec,h3,wtkp,efermi,zdel,
     .  eval,nbmx,iq,dedv)
C- Calculate linear response matrix dE_sx/dv
C ----------------------------------------------------------------------
Ci Inputs
Ci   h3 is a work array
Co Outputs
Co   dedv  is accumulated for this iq.
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
Cu   03 Oct 01 Added nbmx to dimension eval; new arg list
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iq,ldim,nbas,nsp,nbmx
      double precision efermi,wtkp(iq),eval(nbmx,iq),dedv(ldim)
      double precision sigma(ldim,ldim,2),evec(ldim,ldim,2,iq),
     .  h3(ldim,ldim,2)
      double complex zdel
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk(:)
C ... Local parameters
      integer nu1,nuocc,nuunoc

      if (nsp == 2) call rx('desxdv not spin pol')

C      call yprm('sigma',2,sigma,ldim*ldim,ldim,ldim,ldim)
C      call yprm('evec',2,evec(1,1,1,iq),ldim*ldim,ldim,ldim,ldim)

C ... nuocc is highest occ evec; nuunoc is lowest occ
      do  2  nu1 = 1, ldim
      nuunoc = nu1
    2 if(eval(nu1,iq) > efermi) goto 3
    3 continue
      nuocc = nuunoc-1

C ... work space h3 = sigma+ * z(1..nu1)
      call yygemm('C','N',ldim,nuocc,ldim,1d0,sigma,sigma(1,1,2),ldim,
     .  evec(1,1,1,iq),evec(1,1,2,iq),ldim,0d0,h3,h3(1,1,2),ldim)
C     call yprm('sigma+ z',2,h3,ldim*ldim,ldim,ldim,nuocc)

C ... Make dE/dv
      allocate(wk(ldim*2))
      call pesxdv(ldim,nbas,evec(1,1,1,iq),h3,abs(wtkp(iq)),zdel,
     .  eval(1,iq),dedv,wk,nuocc)
      deallocate(wk)
C     call yprm('desx/dv',1,dedv,0,ldim,ldim,1)

      end

      subroutine pesxdv(ldim,nbas,z,h3,wtkp,zdel,eval,dedv,wk,nuocc)
C- Kernel called by desxdv
      implicit none
      integer ldim,nbas,nuocc
      double precision wtkp,eval(ldim),dedv(ldim)
      double precision z(ldim,ldim,2),h3(ldim,ldim,2),wk(ldim,2)
      double complex zdel
C Local variables
      double precision ez,xxr,xxi
      integer i,nu1,nu2

C --- Loop over (occ,unocc) band pairs ---
      do  10  nu1 = 1, nuocc
C ... 12, 14 loops written independently to encourage unrolling
      do  12  nu2 = nuocc+1, ldim

C   --- Make sum_ij z+(j,nu1) sigma(j,i) z(i,nu2) ---
C   ... Use h3(*,nuocc+1) to permit unrolling of nu2 loop
        wk(nu2,1) = 0
        wk(nu2,2) = 0
        do  24  i = 1, ldim
          wk(nu2,1) = wk(nu2,1) +
     .      h3(i,nu1,1)*z(i,nu2,1) + h3(i,nu1,2)*z(i,nu2,2)
          wk(nu2,2) = wk(nu2,2) +
     .      h3(i,nu1,1)*z(i,nu2,2) - h3(i,nu1,2)*z(i,nu2,1)
   24   continue
   12 continue

C --- dedv = z+(i,nu2) z(i,nu1) sum_ij z+(nu1) sigma z(nu2) ---
      do  14  nu2 = nuocc+1, ldim
        ez = 2*wtkp/(eval(nu1)-eval(nu2) + zdel)
        do  30  i = 1, ldim
          xxr = z(i,nu1,1)*z(i,nu2,1) + z(i,nu1,2)*z(i,nu2,2)
          xxi = z(i,nu1,2)*z(i,nu2,1) - z(i,nu1,1)*z(i,nu2,2)
          dedv(i) = dedv(i) + (xxr*wk(nu2,1) -xxi*wk(nu2,2))*ez
   30   continue
   14 continue

   10 continue

      end
