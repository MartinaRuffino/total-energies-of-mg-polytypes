      subroutine addoff(tbc,nl,nbas,lmx,ipc,indxsh,ldim,pot0,sk,hk,typesz)
!C- Add self consistent off-site terms to H_k
!C ----------------------------------------------------------------------
!Ci Inputs:
!Ci   nl,nbas,lmx,ipc,indxsh,ldim,
!Ci   indH: work length nbas
!Ci   pot0: monopole potentials, from tbesel
!Ci   pot0: monopole potentials, from tbesel
!Ci   sk,hk: Overlap and H in k-representation (real part only)
!Co Outputs:
!Co   hk: off site terms added (eq 7.81 Mike's book)
!Cr Remarks
!C ----------------------------------------------------------------------
      use tbprl, only : tbc_t

      implicit none
!C Passed Parameters
      type(tbc_t), intent(in) :: tbc
      integer nl,nbas,indxsh(nl**2*nbas),lmx(1),ipc(nbas),ldim,indH(nbas),typesz
      real(8) :: pot0(nbas),hk(typesz,ldim,*),sk(typesz,ldim,*)

      interface
         function islocal(coords,blcksz,procsz,proccr)
            implicit none
            integer, intent(inout) :: coords(2) ! global coordinates of the element
            integer, intent(in), dimension(2) :: blcksz, & ! block size
                                                procsz, & ! process array size
                                                proccr    ! process coordinates
            logical :: islocal
         end function islocal
      end interface

!C Local Variables
      integer nlm,i,j,ib,jb,ixi,ixj,nlmi,nlmj,lmi,lmj,ipa,pid
      integer, dimension(2) :: ij,blcksz,procsz,proccr
!C Intrinsic function
      nlm(i) = (1+lmx(i))**2

!C --- Form table of indices that mark offset of Rth block in H(k) ---
      indH(1) = 0
      do  ib = 2, nbas
        indH(ib) = indH(ib-1) + nlm(ipc(ib-1))
      enddo


      blcksz = tbc % blcksz
      procsz = tbc % c2d % szs(0:1)
      proccr = tbc % c2d % crd(0:1)

      pid = tbc % c2d % id

      if (.not. tbc%sl) then
         do ib = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)
            do ipa = tbc % neighmap(ib)+1, tbc % neighmap(ib+1)
               jb = tbc % neighidx(ipa)
               if (ib == jb) cycle
               ixi = indH(ib)
               ixj = indH(jb)
               nlmi = nlm(ipc(ib))
               nlmj = nlm(ipc(jb))
               do lmi = 1, nlmi
                  do  lmj = 1, nlmj
                     i = indxsh(ixi+lmi)
                     j = indxsh(ixj+lmj)
                     if (i <= ldim .and. j <= ldim) then
!                         hk(1:typesz,i,j)=hk(1:typesz,i,j) &
!                               & + (pot0(ib) + pot0(jb))*sk(1:typesz,i,j)*0.5d0
                        hk(1,i,j) = hk(1,i,j) + (pot0(ib) + pot0(jb))*sk(1,i,j)*0.5d0
                     end if
                  end do
               end do
            end do
         end do
      else
!This is dangerous. What happens when d2amap(pid)+1:d2amap(pid+1) does not hold all atoms
! who happen to have elements in the local hk? reverting to ib=1:nbas //DMT
!          do ib = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)
         do ib = 1, nbas
            do ipa = tbc % neighmap(ib)+1, tbc % neighmap(ib+1)
               jb = tbc % neighidx(ipa)
               if (ib == jb) cycle
               ixi = indH(ib)
               ixj = indH(jb)
               nlmi = nlm(ipc(ib))
               nlmj = nlm(ipc(jb))
               do lmi = 1, nlmi
                  do  lmj = 1, nlmj
                     ij = [indxsh(ixi+lmi),indxsh(ixj+lmj)] - 1
                     if (all(ij < ldim)) then
                        if (islocal(ij,blcksz,procsz,proccr)) then
                           ij = ij + 1
!                            hk(1:typesz,ij(1),ij(2)) = hk(1:typesz,ij(1),ij(2)) &
!                                  & + (pot0(ib) + pot0(jb))*sk(1:typesz,ij(1),ij(2))*0.5d0
                           hk(1,ij(1),ij(2)) = hk(1,ij(1),ij(2)) &
                                 & + (pot0(ib) + pot0(jb))*sk(1,ij(1),ij(2))*0.5d0
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end if








      end subroutine addoff
