      subroutine tbfrc3(nbas,nlmq1,nlmq,ipc,lmxl,vm,qmpol,f)
C- TB forces due to multipole charges (lmf way)
C-----------------------------------------------------------------------
Ci Inputs:
Ci    nbas  : number of atoms in unit cell
Ci    nlmq  : L-cutoff for multipoles, leading dimension of qmpol
Ci    nlmq1 : leading dimension of vm, nlmq1=(ll(nlmq)+1)**2
Ci    ipc   : ipc(ib) class to which atom ib belongs
Ci    lmxl  : max angular momentum of multipoles Q_L for each class
Ci    vm    : Madelung potential due to charge transfer
Ci    qmpol : multipole moments on each site (from tbmpol.f)
Co Outputs:
Co   f(3,nbas): es forces on each atom
Cl Local:
Cl   vmg    : gradient of vm
Cr Remarks
Cu Updates
Cu   27 Feb 10 (SL) alternative version of tbfrc2.f, uses scglp1 instead
Cu                  of going through Gaunt coefficients
C-----------------------------------------------------------------------
      implicit none

C Passed parameters
      integer, intent(in) :: nbas,nlmq1,nlmq,ipc(nbas),lmxl(*)
      double precision, intent(in)    :: qmpol(nlmq,nbas)
      double precision, intent(inout) :: vm(nlmq1,nbas)
      double precision, intent(out)   :: f(3,nbas)
C Local variables
      integer :: ib,i,ilm,ilm1,ilm2,il,ic,nlm,nlm1
      integer :: ll,ierr
      double precision cc
C Automatic arrays
c     double precision :: vmg(nlmq,nbas,3)
      integer          :: kz(nlmq),kx1(nlmq),kx2(nlmq),
     .                    ky1(nlmq),ky2(nlmq)
      double precision :: cz(nlmq),cx1(nlmq),cx2(nlmq),
     .                    cy1(nlmq),cy2(nlmq)
C Allocatable arrays
      double precision, allocatable :: vmg(:,:,:)

c     call tcn('tbfrc3')

      if (allocated(vmg))
     .  call rx0(' tbfrc3: array vmg is already allocated')
      allocate(vmg(nlmq,nbas,3),stat=ierr)
      if (ierr /= 0) call rx0(' tbfrc3: failed to allocate vmg')

c...  Reshuffle vm to Jackson (lm)-indecies
      do ib = 1, nbas
        ic = ipc(ib)
        nlm1 = (lmxl(ic)+2)**2
        call tbshfl(1,nlmq1,nlm1,1,vm(1,ib))
      enddo

c...  Prepare Gaunt coefficients needed for gradients
      do  ilm = 1, nlmq
        call scglp1(ilm,kz(ilm),cz(ilm),kx1(ilm),kx2(ilm),
     .    cx1(ilm),cx2(ilm),ky1(ilm),ky2(ilm),cy1(ilm),cy2(ilm))
      enddo
c Probably unnecessary
c       kmax = max(kx1,kx2,ky1,ky2,kz)
c       if (kmax > nlmq1)
c    .    call rxi(' tbfrc3: insufficient size of vmm. Need ',kmax)

c...  Make vmg = grad(vm) for all atoms
      vmg(1:nlmq,1:nbas,1:3) = 0d0
      do ib = 1, nbas
        ic = ipc(ib)
        nlm = (lmxl(ic)+1)**2
        do  ilm = 1, nlm
          vmg(ilm,ib,1) = vmg(ilm,ib,1) + cx1(ilm)*vm(kx1(ilm),ib)
     .                      +cx2(ilm)*vm(kx2(ilm),ib)
          vmg(ilm,ib,2) = vmg(ilm,ib,2) + cy1(ilm)*vm(ky1(ilm),ib)
     .                      +cy2(ilm)*vm(ky2(ilm),ib)
          vmg(ilm,ib,3) = vmg(ilm,ib,3) + cz(ilm)*vm(kz(ilm),ib)
        enddo
      enddo

c...  reshuffle vm and grad(vm) back to Stone's convention
      do ib = 1, nbas
        ic = ipc(ib)
        nlm = (lmxl(ic)+1)**2
        nlm1 = (lmxl(ic)+2)**2
        call tbshfl(0,nlmq1,nlm1,1,vm(1,ib))
        do i = 1, 3
          call tbshfl(0,nlmq,nlm,1,vmg(1,ib,i))
        enddo
      enddo

c...  rescale grad(vm) to Stone's convention. This absorbs the
c     rescaling of vm by sqrt(4pi/2l+3), Eq.(A2) in Tony's notes,
c     which therefore needn't be done explicitly
      do il = 0, ll(nlmq)
        ilm1 = il*il+1
        ilm2 = (il+1)**2
        cc = dsqrt(dfloat((2*il+1)*(2*il+3)))
        vmg(ilm1:ilm2,1:nbas,1:3) = vmg(ilm1:ilm2,1:nbas,1:3)*cc
      enddo

C...  Make forces as F(R) = - sum_L grad_R(Vm_L)*Q_RL
      f(1:3,1:nbas) = 0d0
      do ib = 1, nbas
        ic = ipc(ib)
        nlm = (lmxl(ic)+1)**2
        do ilm = 1, nlm
          f(1:3,ib) = f(1:3,ib) - vmg(ilm,ib,1:3)*qmpol(ilm,ib)
        enddo
      enddo

      deallocate(vmg)

c     call tcx('tbfrc3')
      end
