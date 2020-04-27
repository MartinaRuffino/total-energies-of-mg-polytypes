      subroutine tbadh2(nl,nbas,nspu,nsites,npr,dh,h0,hrs)
C- Add mixed increments to H to real-space H
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nl,nbas,nspu,nsites,npr
Ci   dh, dho: increments to on-site and off-site hamiltonian
Ci   h0 : input hamiltonian, H_in
Ci   ors : real space overlap matrix
Co Outputs:
Co   hrs : real space hamiltonian: H_0 + H'
Cr Remarks
Cr   This routine for density mixing branch. The density matrix, rhoc
Cr   or rhon, is mixed after bndtb and before tbesel. The components
Cr   of H' returned from tbesel are added here to H_in
C ----------------------------------------------------------------------
      implicit none
c     integer niax
c     parameter (niax=10)
C Passed Parameters
      integer nl,nbas,nspu,nsites,npr(0:1,nbas)
c    .        ipc(*),iax(niax,nsites)
      double precision dh(nl**2,nl**2,nbas,nspu),
c    .                 pos(3,*),plat(3,3),
     .                 h0(nl**2,nl**2,nsites,nspu),
     .                 hrs(nl**2,nl**2,nsites,nspu)
c    .                 ors(nl**2,nl**2,nsites)

C Local Variables
      integer ispu,ib,j,ilm,ilmp,iprint

      if (iprint() > 40) print 100
      do  ib = 1, nbas
        do ispu = 1, nspu
          if (iprint() > 40 .and. nspu == 2) write (*,150) ispu
          j = npr(1,ib) + 1
          do  ilm = 1, nl**2
            do  ilmp = 1, nl**2

C --- verbose output ---
              if ((ilm == ilmp .and. iprint() > 40)
     .          .or. iprint() > 50) then
                write (*,200)
     .          ib,ilm,ilmp,hrs(ilm,ilmp,j,ispu),
     .          h0(ilm,ilmp,j,ispu)+dh(ilm,ilmp,ib,ispu)
              endif
C ----------------------

              hrs(ilm,ilmp,j,ispu) = h0(ilm,ilmp,j,ispu)
     .                             + dh(ilm,ilmp,ib,ispu)
            enddo
          enddo
        enddo
      enddo

  100 format(' TBADH2: site ilm ilm''   H_old   H_new')
  150 format ('        spin ',i1)
  200 format (7x,3i4,2x,2f10.6)


      end
