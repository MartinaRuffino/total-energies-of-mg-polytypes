      subroutine asagqr(n1,n2,n3,nsp,idim,jdim,wz,ldiag,tij,tji,gR)
C- Makes G(R) from G(q) for a block of orbitals by the FT technique
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1,n2,n3  QP mesh
Ci   wz        integration weight
Ci   idim,jdim dimensions of G(i,j) subblock
Ci   tij,tji   T-matrix connecting orbitals ij and ji
Ci   ldiag     ij is diagonal block, in which case tji is not
Ci             explicitly available, but must be computed from tij
Co  Outputs
Co   gR        g(R,R') for R=0,R'=0,nn
Cr  Remarks
C ----------------------------------------------------------------------
      implicit none
      logical ldiag
      integer idim,jdim
      integer n1,n2,n3,nsp
      double complex wz
      double complex tij(n1,n2,n3,idim,jdim,nsp),gR(idim,jdim,2)
      double complex tji(n1,n2,n3,jdim,idim,nsp)
C Local variables
      integer id,jd,k
      double complex wpi,wpi2

C ... Weight for integration 1/pi de ...
      wpi = wz/(4*datan(1d0))
      wpi2 = wpi/2

C --- Loop over orbital pairs ---
      do  10  id = 1, idim
      do  10  jd = 1, jdim

C        print *, id,jd
C        call zprm3(' tij %i',id+100+jd,tij(1,1,1,id,jd,1),n1,n2,n3)
        call fftz3(tij(1,1,1,id,jd,1),n1,n2,n3,n1,n2,n3,1,0,-1)
        if (.not. ldiag) then
C          call zprm3(' tji %i',id+100+jd,tji(1,1,1,id,jd,2),n1,n2,n3)
          call fftz3(tji(1,1,1,jd,id,1),n1,n2,n3,n1,n2,n3,1,0,-1)
        endif

C       call zprm3(' tij %i',id+100+jd,tij(1,1,1,id,jd,1),n1,n2,n3)
        do  20  k = 1, 2
   20   gR(id,jd,k) = tij(k,1,1,id,jd,1)

   10 continue
      call zprm('gR00',2,gR,idim,idim,jdim)
      call zprm('gR01',2,gR(1,1,2),idim,idim,jdim)

      end
