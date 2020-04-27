      subroutine hml2nc(nbas,nl,indxsh,qspirl,eula,neul,
     .  pph,ccor,lss,lnc,ccd,wk,vmtz,elin,ldim,sk,sd,hk)
C- Generate two-center ASA noncollinear hamiltonian
C ---------------------------------------------------------------------
Ci Inputs
Ci   ccor,lss,lnc switches for comb. corr, spin spiral and noncollinear
Ci   qspirl(4): rotation angle of spin spiral
Ci   eula:  euler angles of noncollinear hamiltonian
Ci   pph:  vector of potential parameters (see makpph)
Ci   ccd, diagonal matrices constant, linear and
Ci   bilinear in S^alpha: they are the terms in parentheses in eq.3.87
Ci   in the Kanpur notes multiplied by w^2; vmtz, muffin-tin zero;
Ci   ldim: dimension of the hamiltonian
Ci   sk,sd: structure constants (sdot used only if ccor is true).
Ci   wk: work array of length ldim*2
Co Outputs
Co   hk
Cr Remarks
Cr   downfolding not implemented.
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,neul,nl,indxsh(*),ldim
      logical ccor,lss,lnc
      double precision ccd(ldim,0:2),eula(nbas,neul,3),pph(5,ldim,2),
     .  sk(ldim,ldim,2,2),sd(ldim,ldim,2,2),hk(ldim,2,ldim,2,2),
     .  vmtz,elin,qspirl(4),wk(ldim,2)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wkl(:)
C ... Local parameters
      integer i,j,ncsw,i1
      double precision xx

      call tcn('hml2nc')
      ncsw = 2000
      if (lss) ncsw = ncsw + 20000

C --- For ccor:  add (vmtz-elin)*(ccd(i,1)+ccd(j,1)) to s_ij ---
      if (ccor) then
        do j = 1, ldim
          do i = 1, ldim
          xx = 1 + (vmtz-elin) * (ccd(i,1)+ccd(j,1))
          sk(i,j,1,1) = sk(i,j,1,1)*xx
          sk(i,j,2,1) = sk(i,j,2,1)*xx
          if (lss) then
            sk(i,j,1,2) = sk(i,j,1,2)*xx
            sk(i,j,2,2) = sk(i,j,2,2)*xx
          endif
          enddo
        enddo
      endif

C --- hk <- C + d*S*d + (vmtz-elin)d*Sdot*d with d = sqrt(delta) ---
      call dpzero(hk,ldim*2*ldim*2*2)
      if (ccor) then
        call dscal(ldim**2*2,vmtz-elin,sd,1)
        if (lss) call dscal(ldim**2*2,vmtz-elin,sd(1,1,1,2),1)
      endif
      if (lnc .or. lss) then

        do i1 = 1, 2
          do i = 1, ldim
            wk(i,i1) = pph(3,i,i1)
          enddo
        enddo

        call rotspn(ncsw,1,nbas,nbas,nl,indxsh,eula,neul,qspirl(4),
     .    wk,wk,ldim,ldim,ldim,ldim,ldim,sk,hk)

        if (ccor) then
          allocate(wkl(ldim*2*ldim*2*2))
          call rotspn(ncsw,1,nbas,nbas,nl,indxsh,eula,neul,qspirl(4),
     .      wk,wk,ldim,ldim,ldim,ldim,ldim,sd,wkl)
          call dpadd(hk,wkl,1,ldim*2*ldim*2*2,1d0)
          deallocate(wkl)
        endif

C       call yprm('hk',02,hk,(ldim*2)**2,ldim*2,ldim*2,ldim*2)

      endif

C --- H += C + (vmtz-eln)*<k|k>_constant ---
      call daxpy(ldim,1d0,pph(2,1,1),5,hk,2*ldim+1)
      call daxpy(ldim,1d0,pph(2,1,2),5,hk(1,2,1,2,1),2*ldim+1)
      if (ccor) then
        do  i = 1, ldim
          hk(i,1,i,1,1) = hk(i,1,i,1,1) +
     .      (vmtz-elin)*ccd(i,0)*pph(3,i,1)**2
          hk(i,2,i,2,1) = hk(i,2,i,2,1) +
     .      (vmtz-elin)*ccd(i,0)*pph(3,i,2)**2
      enddo
      endif

C     call prmx('h in hml2nc',hk,ldim*2,ldim*2,ldim*2)

      call tcx('hml2nc')

      end
