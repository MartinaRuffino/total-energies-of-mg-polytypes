      subroutine tbfrc2(nbas,nlmq1,nlmq,ipc,lmxl,vm,qmpol,
     .  indxcg,jcg,cg,f)
C- TB forces from point multipoles
C-----------------------------------------------------------------------
Ci Inputs:
Ci    nbas  : number of atoms in unit cell
Ci    nlmq  : L-cutoff for multipoles, leading dimension of qmpol
Ci    nlmq1 : leading dimension of vm, nlmq1=(ll(nlmq)+1)**2
Ci    ipc   : ipc(ib) class to which atom ib belongs
Ci    lmxl  : max angular momentum of multipoles Q_L for each class
Ci    vm    : Madelung potential due to charge transfer
Ci    qmpol : multipole moments on each site (from tbmpol.f)
Ci    gaunt1: Gaunt coefficients
Co Outputs:
Co   f(3,nbas): es forces on each atom
Cl Local:
Cl   dVn    : gradient of vm
Cr Remarks
Cu Updates
Cu   24 Feb 10 (SL) split off from tbesel.f, L summation limits modified
C-----------------------------------------------------------------------
      implicit none

C Passed parameters
      integer, intent(in) :: nbas,nlmq1,nlmq,ipc(nbas),lmxl(*)
      double precision, intent(in)  :: vm(nlmq1,nbas),qmpol(nlmq,nbas)
      double precision, intent(out) :: f(3,nbas)
C Clebsch-Gordan coefficients
      integer, intent(in) :: jcg(*),indxcg(*)
      double precision, intent(in) :: cg(*)
C Local variables
      integer :: ib,i,ilm,ilmp,ilmpp,il,ilc,nlm,nlm1,ic
      integer :: ll
C Automatic arrays
      double precision :: gaunt1(9,nlmq,nlmq1)
      double precision :: dVn(nlmq),pi

c     call tcn('tbfrc2')

      pi = 4d0*datan(1d0)
      f(1:3,1:nbas) = 0d0

C...  Make Gaunt coefficients
      call makcgn(9,nlmq,nlmq1,indxcg,jcg,cg,gaunt1)

      do  ib = 1, nbas
        ic = ipc(ib)
        nlm1 = (lmxl(ic)+2)**2
        nlm  = (lmxl(ic)+1)**2
c       if (nlmq1 < 4) call
c    .    rxi(' tbfrc2: inconcistency in the force branch. nlmq1 =',
c    .    nlmq1)
        do  i = 1, 3
          ilm = i+1
          f(i,ib) = f(i,ib) - vm(ilm,ib)*qmpol(1,ib)
c ... sum_m' Q_1m' * (d/dx_i)V_1m'
          if (nlm1 < 9) cycle
          dVn(1:nlm) = 0d0
          do  ilmp = 2, 4
C ... make (d/dx_i)V_1m'
            do  ilmpp = 5, 9
              dVn(ilmp) = dVn(ilmp) +
     .             gaunt1(ilm,ilmp,ilmpp)*vm(ilmpp,ib)
            enddo
c ... coefficient below is sqrt[4pi*(2l+1)*(2l+3)/3], where l = ll(ilmp)
            il = ll(ilmp)
            ilc = 4*(2*il+1)*(2*il+3)
            f(i,ib) = f(i,ib) -
     .             dsqrt(dfloat(ilc)*pi/3d0)*qmpol(ilmp,ib)*dVn(ilmp)
c    .                  dsqrt(20d0*pi)*qmpol(ilmp,ib)*dV1(ilmp)
          enddo
C ... sum_m' Q_2m' * (d/dx_i)V_2m'
          if (nlm1 < 16) cycle
          do  ilmp = 5, 9
C ... make (d/dx_i)V_2m'
            do  ilmpp = 10, 16
              dVn(ilmp) = dVn(ilmp) +
     .             gaunt1(ilm,ilmp,ilmpp)*vm(ilmpp,ib)
            enddo
            il = ll(ilmp)
            ilc = 4*(2*il+1)*(2*il+3)
            f(i,ib) = f(i,ib) -
     .            dsqrt(dfloat(ilc)*pi/3d0)*qmpol(ilmp,ib)*dVn(ilmp)
          enddo
C ... sum_m' Q_3m' * (d/dx_i)V_3m'
          if (nlm1 < 25) cycle
          do  ilmp = 10, 16
C ... make (d/dx_i)V_3m'
            do  ilmpp = 17, 25
              dVn(ilmp) = dVn(ilmp) +
     .             gaunt1(ilm,ilmp,ilmpp)*vm(ilmpp,ib)
            enddo
            il = ll(ilmp)
            ilc = 4*(2*il+1)*(2*il+3)
            f(i,ib) = f(i,ib) -
     .            dsqrt(dfloat(ilc)*pi/3d0)*qmpol(ilmp,ib)*dVn(ilmp)
          enddo
C ... sum_m' Q_4m' * (d/dx_i)V_4m'
          if (nlm1 < 36) cycle
          do  ilmp = 17, 25
C ... make (d/dx_i)V_3m'
            do  ilmpp = 26, 36
              dVn(ilmp) = dVn(ilmp) +
     .             gaunt1(ilm,ilmp,ilmpp)*vm(ilmpp,ib)
            enddo
            il = ll(ilmp)
            ilc = 4*(2*il+1)*(2*il+3)
            f(i,ib) = f(i,ib) -
     .            dsqrt(dfloat(ilc)*pi/3d0)*qmpol(ilmp,ib)*dVn(ilmp)
          enddo
        enddo                                  ! end of loop over i
      enddo                                    ! end of loop over atoms

c     call tcx('tbfrc2')
      end
