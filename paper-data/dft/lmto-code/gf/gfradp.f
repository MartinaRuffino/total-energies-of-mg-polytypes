      subroutine gfradp(mode,ldg,pffr,gf)
C- Add fully relativistic Potential functions to -S
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :kcplx distinguishes how gf stores complex matrices
Ci          0: real, imaginary separated: gf = gf(ldg,ldg,1..2)
Ci          1: complex*16: gf = gf(ldg,ldg)
Ci          2: real, imaginary in columns : gf = gf(ldg,1..2,ldg)
Ci         :10s digit
Ci          0: m is ordered m=l,..-l
Ci          1: m is ordered m=-l,..l
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
Cu  11 May 18 Added 100s digit option for m ordering
Cu  22 Feb 03 first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ldg
      double complex pffr(ldg,2,2),gf(ldg,2,ldg,2)
C ... Local parameters
      integer i,ldgx,kcplx,mod1,im

C     call dpzero(gf,ldg*2*ldg*2*2)
      ldgx = ldg*2
      kcplx = mod(mode,10)
      mod1 = mod(mode/10,10)
!     im = 1 ; if (mod1 == 2) im = -1
      im = -1 ; if (mod1 == 1) im = 1

C ... Internally we make g double complex
C     call yprm('starting S',2,gf,(ldgx)**2,ldgx,ldgx,ldgx)
      call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,kcplx,1)
C     call zprm('ending S',2,gf,ldgx,ldgx,ldgx)

C ... Add P to -S, diagonal parts
      do  i = 1, ldg
        gf(i,1,i,1) = gf(i,1,i,1) + pffr(i,1,1)
        gf(i,2,i,2) = gf(i,2,i,2) + pffr(i,2,2)
      enddo

C     call zprm('ending P-S (diag)',2,gf,ldgx,ldgx,ldgx)

C ... Add P to -S, off-diagonal parts
      do  i = 1, ldg
        if (i+im < 1 .or. i+im > ldg) cycle
C        if (mod1 == 1) then
C          gf(i+im,1,i,2) = gf(i+im,1,i,2) + pffr(i+im,1,2)
C          gf(i,2,i+im,1) = gf(i,2,i+im,1) + pffr(i,2,1)
C        else
          gf(i,1,i+im,2) = gf(i,1,i+im,2) + pffr(i+im,1,2)
          gf(i+im,2,i,1) = gf(i+im,2,i,1) + pffr(i,2,1)
C        endif
      enddo

C     call zprm('ending P-S',2,gf,ldgx,ldgx,ldgx)
C     Alternatively write in binary format
C     i = fopna('pms',-1,4)
C     call ywrm(1,' ',3,i,'(5f15.9)',gf,0,ldgx,ldgx,ldgx)
C     call fclose(i)
C     call rx('done')

C ... restore g to expected format
      if (kcplx /= 1) then
        call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,1,kcplx)
      endif

      end

      subroutine gtpffr(pfall,nzp,lhdim,iz,indx,pffr,ldg)
C indx=1->P  indx=2->Pdot indx=3-> -1/2 Pdotdot/pdot indx=4->sqrt(Pdot)
      integer iz,nzp,lhdim,ldg,indx
      double complex pfall(nzp,lhdim,4,2,2)
      double complex pffr(ldg,2,2)
      integer i,i1,i2

      do i1 = 1,2
      do i2 = 1,2
        do i = 1, ldg
          pffr(i,i1,i2) = pfall(iz,i,indx,i1,i2)
        enddo
      enddo
      enddo

C      print *, 'iz=',iz
C      call zprm('pf-rel',2,pffr,ldg,ldg,4)

      end
