      module gcpa
      contains
      subroutine gc00(mode,norb,nspc,lrel,lso,th,pfi,pfr,omg,g00)

      implicit none
      integer mode,norb,nspc,lrel
      logical lso
      real(8) th(nspc)
      complex(8), intent(in), optional :: pfi(norb,2,1+lrel/2),
     .  pfr(norb,2,norb,2)
      complex(8), intent(in) :: omg(norb,nspc,norb,2)
      complex(8), intent(out) :: g00(norb,2,norb,2)

      complex(8), allocatable :: wk(:,:,:,:)
      real(8) pi,eula(3)
      real(8), allocatable :: wrk(:)
      integer nl,ierr
      logical lrot

      if (lso .and. .not. present(pfr))
     .  call rx('gcon00: lso requires pfr passed as argument')
      if (.not. lso .and. .not. present(pfi))
     .  call rx('gcon00: not lso and pfi not passed as argument')
      pi = 4*datan(1d0)
      nl = sqrt(norb + 0.1)
      allocate(wk(norb,2,norb,2))

C ... wk <- (-omg)
      if (nspc == 2) then
        wk = - omg
      else
        wk = 0
        wk(:,1,:,1) = - omg(:,1,:,1)
        wk(:,2,:,2) = - omg(:,1,:,2)
      endif

C ... Rotate Omega to frame of local Bxc and add Pa in that frame
      lrot = th(1) /= 0d0
      if (lrel == 2 .or. lso) then
        if (lrot) then
          eula(1) = th(2)
          eula(2) = th(1) ; eula(3) = - eula(1)
          call rot_LandS(11,eula,nl,norb,1,wk)
        endif
        if (lso) then
          call daxpy(2*norb*norb*4,1d0,pfr,1,wk,1)
        else
          call pfr2block(1,nl,0,0,0,0,norb,pfi,norb,1,wk)
        endif
      else
        if (lrot) call rotm(norb,-th(1),wk)
        call adrotp(0,norb,pfi(:,1,1),pfi(:,2,1),0d0,wk)
      endif

C ... Invert wk in place to make (Pa-Omega)^-1
      allocate(wrk(66*2*norb))
      call zqinv('N',wk,2*norb,-66*2*norb,2*norb,wrk,2*norb,ierr)
      if (ierr /= 0) call rx1('mkgcpa: matrix inversion,ierr=%i',ierr)
      deallocate(wrk)

      g00 = wk

      deallocate(wk)
      end subroutine gc00

      end module gcpa
