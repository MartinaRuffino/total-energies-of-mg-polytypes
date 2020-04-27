      subroutine p0herm(nRLc,nsp,nkp,P0)
C- Makes the response matrix hermitian
C ----------------------------------------------------------------------
Ci Inputs
Ci   nRLc  :dimension of P0
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of k-points for which P0 is available
Ci   P0    :response matrix = 1/i G delta v G
Co Outputs
Co    P0 is made hermitian
Cl Local variables
Cr Remarks
Cr Made Feb 20 2004 (T. Sandu)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nRLc,nsp,nkp
      double complex P0(nRLc,nRLc,nkp,nsp)
C ... Local parameters
      integer i1,ik,ii,ij
      double complex B(nRLc,nRLc,nkp,nsp)
c      print *, 'into p0herm'
C --- Make P0(q) hermitian
      do  88  ik = 1, nkp
      do  88  i1 = 1, nsp
      do  88  ii = 1, nRLc
      do  88  ij = 1, nRLc
   88 B(ii,ij,ik,i1) = DCONJG(P0(ij,ii,ik,i1))
      do  100  i1 = 1, nsp
      do  100  ik = 1, nkp
      do  100  ii = 1, nRLc
      do  100  ij = 1, nRLc
      P0(ii,ij,ik,i1) = 0.5d0*(P0(ii,ij,ik,i1)
     .+B(ii,ij,ik,i1))
  100 continue
      end

      subroutine p0hermrl(nRLc,nsp,nkp,P0)
C- Makes the response matrix hermitian complex format=0
C ----------------------------------------------------------------------
Ci Inputs
Ci   nRLc  :dimension of P0
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of k-points for which P0 is available
Ci   P0    :response matrix = 1/i G delta v G
Co Outputs
Co    P0 is made hermitian
Cl Local variables
Cr Remarks
Cr
Cr Made June 17 2004 (T. Sandu)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nRLc,nsp,nkp
      double precision P0(nRLc,nRLc,nkp,nsp,2)
C ... Local parameters
      integer i1,ik,ii,ij
      double precision B(nRLc,nRLc,nkp,nsp,2)
c      print *, 'into p0herm'
C --- Make P0(q) hermitian
      do  88  i1 = 1, nsp
      do  88  ik = 1, nkp
      do  88  ii = 1, nRLc
      do  88  ij = 1, nRLc
      B(ii,ij,ik,i1,1) = P0(ij,ii,ik,i1,1)
      B(ii,ij,ik,i1,2) = -P0(ij,ii,ik,i1,2)
 88   continue
      do  100  i1 = 1, nsp
      do  100  ik = 1, nkp
      do  100  ii = 1, nRLc
      do  100  ij = 1, nRLc
      P0(ii,ij,ik,i1,1) = 0.5d0*(P0(ii,ij,ik,i1,1)
     .+B(ii,ij,ik,i1,1))
      P0(ii,ij,ik,i1,2) = 0.5d0*(P0(ii,ij,ik,i1,2)
     .+B(ii,ij,ik,i1,2))
  100 continue
      end

      subroutine p0c(nRLc,nsp,nkp,P0)
C- Makes complex-conjugate of P0
C ----------------------------------------------------------------------
Ci Inputs
Ci   nRLc  :dimension of P0
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of k-points for which P0 is available
Ci   P0    :response matrix = 1/i G delta v G
Co Outputs
Co    P0 is made hermitian
Cl Local variables
Cr Remarks
Cr Made Sept 10 2004(T. Sandu)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nRLc,nsp,nkp
      double complex P0(nRLc,nRLc,nkp,nsp)
C ... Local parameters
      integer i1,ik,ii,ij
c      print *, 'into p0c'
C --- Make P0(q) complec-conjugate
      do  88  ik = 1, nkp
      do  88  i1 = 1, nsp
      do  88  ii = 1, nRLc
      do  88  ij = 1, nRLc
   88 P0(ii,ij,ik,i1) = DCONJG(P0(ii,ij,ik,i1))

      continue
      end
