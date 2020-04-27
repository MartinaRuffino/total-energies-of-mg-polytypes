      subroutine gfradp1(kcplx,ldg,ldim,idim,pffr,gf)
C- Add fully relativistic Potential functions to -S
C ----------------------------------------------------------------------
Ci Inputs
Ci   kcplx :distinguishes how gf stores complex matrices
Ci          0: real, imaginary separated: gf = gf(ldg,ldg,1..2)
Ci          1: complex*16: gf = gf(ldg,ldg)
Ci          2: real, imaginary in columns : gf = gf(ldg,1..2,ldg)
Ci   ldg   :dimension of gf and pffr
Ci   pffr  :Fully relativistic potential functions
Cio Inputs/Outputs
Cio   gf   :On input, -S, where S = structure constants
Cio        :On output, P is added to -S
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   22 Feb 03 first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kcplx,ldg,ldim,idim
      double complex pffr(ldg,2,2),gf(ldg,2,ldg,2)
C ... Local parameters
      integer i,ldgx

C       call dpzero(gf,ldg*2*ldg*2*2)

      ldgx = ldg*2

C ... Internally we make g double complex
      if (kcplx /= 1) then
C       call yprm('starting S',2,gf,(ldgx)**2,ldgx,ldgx,ldgx)
        call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,kcplx,1)
C       call zprm('ending S',2,gf,ldgx,ldgx,ldgx)
      endif

C ... Add P to -S, diagonal parts
      do  i = 1, ldg
        gf(i,1,i,1) = gf(i,1,i,1) + pffr(i,1,1)
        gf(i,2,i,2) = gf(i,2,i,2) + pffr(i,2,2)
      enddo
C ... Add P to -S, off-diagonal parts
      do  i = 2, ldg
        gf(i-1,1,i,2) = gf(i-1,1,i,2) + pffr(i-1,1,2)
        gf(i,2,i-1,1) = gf(i,2,i-1,1) + pffr(i,2,1)
      enddo

      call zprm('ending P-S',2,gf,ldgx,ldgx,ldgx)

C ... restore g to expected format
      if (kcplx /= 1) then
        call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,1,kcplx)
      endif

      end
