      subroutine shffle(scat,ndim,indxsh,vector,wk)
C- Scatter/Gather : re-orders a vector in or out of downfolding order
C-----------------------------------------------------------------------
Ci Inputs
Ci   scat, ndim, indxsh, vector (see remarks); wk work array
Co Outputs
Co   vector
Cr Remarks
Cr   indxsh contains the downfolding indices. If scat is .true.
Cr   the vector is returned scattered into downfolding order.
Cr   If scat is .false. the vector is already so ordered and is
Cr   returned gathered.
Cr----------------------------------------------------------------------
      implicit none
      logical scat
      integer ndim,indxsh(ndim)
      double precision vector(ndim),wk(ndim)
      integer i
      external dcopy

      call dcopy(ndim,vector,1,wk,1)
      if (scat) then
C --- scatter ---
        do  i = 1, ndim
          vector(indxsh(i)) = wk(i)
        enddo
      else
C --- gather ---
        do  i = 1, ndim
          vector(i) = wk(indxsh(i))
        enddo
      endif
      end
