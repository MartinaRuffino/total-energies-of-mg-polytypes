      subroutine gfradpl(kcplx,ldim,idim,ldg,pffr,ogi,gf)
C- Add fully relativistic Potential functions to lower subblock -S
C  Overwrite the intermediate blocks by 1/P
C ----------------------------------------------------------------------
Ci Inputs
Ci   kcplx :distinguishes how gf stores complex matrices
Ci          0: real, imaginary separated: gf = gf(ldg,ldg,1..2)
Ci          1: complex*16: gf = gf(ldg,ldg)
Ci          2: real, imaginary in columns : gf = gf(ldg,1..2,ldg)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   ldg   :dimension of gf and pffr
Ci   pffr  :Fully relativistic potential functions in downfolding order
Ci   ogi   :not used (offset to Im g, if kcplx is not 1)
Cio Inputs/Outputs
Cio   gf   :On input, -S, where S = structure constants
Cio        :On output, P is added to -S
Cl Local variables
Cl         :
Cr Remarks
Cr   Ignore spin-orbital coupling for intermediate waves
Cr   Potential function and Structure constant matrices must both be in
Cr   downfolding order at the entry.
Cu Updates
Cu   22 Feb 03 first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kcplx,ldim,idim,ldg,ogi
      double complex pffr(ldg,2,2),gf(ldg,2,ldg,2)
C ... Local parameters
      integer i,ldgx,m,n,j
C ... External calls
      external ztoyy

C     call dpzero(gf,ldg*2*ldg*2*2)

      ldgx = ldg*2

C ... Internally we make g double complex
      if (kcplx /= 1) then
C       call yprm('starting S',2,gf,(ldgx)**2,ldgx,ldgx,ldgx)
        call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,kcplx,1)
C       call zprm('ending S',2,gf,ldgx,ldgx,ldgx)
      endif

C ... Add P to -S, spin diagonal parts for the lower block
      do  i = 1, ldim
        gf(i,1,i,1) = gf(i,1,i,1) + pffr(i,1,1)
        gf(i,2,i,2) = gf(i,2,i,2) + pffr(i,2,2)
      enddo

C ... Add P to -S, spin off-diagonal parts for the lower block
      do  i = 2, ldim
        gf(i-1,1,i,2) = gf(i-1,1,i,2) + pffr(i-1,1,2)
        gf(i,2,i-1,1) = gf(i,2,i-1,1) + pffr(i,2,1)
      enddo

C ... Initialize the intermediate block, all spin components
      do i = 1, idim
         do j = 1, idim
            do m = 1, 2
               do n = 1, 2
                  gf(i+ldim,m,j+ldim,n) = 0.0d0
               enddo
            enddo
         enddo
      enddo

C ... Overwrite the intermediate blocks by P, diagonal parts
C     Off-diagonal parts are ignored
      do  i = 1, idim
        gf(i+ldim,1,i+ldim,1) = pffr(i+ldim,1,1)
        gf(i+ldim,2,i+ldim,2) = pffr(i+ldim,2,2)
      enddo

C     call zprm('ending P-S dnf',2,gf,ldgx,ldgx,ldgx)

C ... Restore g to expected format
      if (kcplx /= 1) then
        call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,1,kcplx)
      endif

      end









