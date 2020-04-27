      subroutine yhmlso(ldim,sk,so3c,hk)
C- 3-center terms of the L-S+ hamiltonian, noncollinear case.
C ----------------------------------------------------------------
Ci Inputs
Ci   sk
Ci   so3c: quasidiagonal matrix in the (1,2) block (see mksod)
Co Outputs
Co   hk:  accumulate sk * so3c * sk into hk
Cr Remarks
Cr   matrix multiplication done by blocks for efficiency
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldim
      double precision sk(ldim,2,ldim,2,2),hk(ldim,2,ldim,2,2),so3c(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk(:)
      real(8), allocatable :: wk2(:)
C ... Local parameters
      integer nrow
      parameter (nrow=48)

      allocate(wk(nrow*ldim*2))
      allocate(wk2(nrow*ldim*2*2))
      call xyhmls(nrow,ldim,sk,so3c,wk,wk2,hk)
      deallocate(wk,wk2)
      end
      subroutine xyhmls(nrow,ldim,sk,so3c,wk,wk2,hk)
      implicit none
      integer ldim,mrow,irow,nrow
      double precision so3c(*),wk(nrow,ldim,2),wk2(nrow,ldim,2,2),
     .  sk(ldim,2,ldim,2,2),hk(ldim,2,ldim,2,2)
      integer i,j

C --- Matrix multiplication by blocks of size nrow ---
      do  irow = 1, ldim*2, nrow
      mrow = min(2*ldim-irow+1,nrow)

C ...   sk * soc
        do  j = 1, ldim
          do  i = 1, mrow
          wk(i,j,1) = sk(i+irow-1,1,j,1,1)*so3c(j)
          wk(i,j,2) = sk(i+irow-1,1,j,1,2)*so3c(j)
          enddo
        enddo

C       call yprm('s.sod.4',2,wk,nrow*ldim,nrow,mrow,ldim)

C ...   sk * soc * sk
        call yygemm('N','N',mrow,2*ldim,ldim-1,1d0,wk,wk(1,1,2),nrow,
     .  sk(2,2,1,1,1),sk(2,2,1,1,2),ldim*2,0d0,wk2,wk2(1,1,1,2),nrow)

C       call yprm('s.sod.4.s',2,wk2,nrow*ldim*2,nrow,mrow,ldim*2)

C ...   Add into hk and Hermitian congugate into hk for L+.S- block
        do  j = 1, ldim*2
          do  i = 1, mrow
          hk(i+irow-1,1,j,1,1) = hk(i+irow-1,1,j,1,1) + wk2(i,j,1,1)
          hk(i+irow-1,1,j,1,2) = hk(i+irow-1,1,j,1,2) + wk2(i,j,1,2)
          hk(j,1,i+irow-1,1,1) = hk(j,1,i+irow-1,1,1) + wk2(i,j,1,1)
          hk(j,1,i+irow-1,1,2) = hk(j,1,i+irow-1,1,2) - wk2(i,j,1,2)
          enddo
        enddo

      enddo

      end
