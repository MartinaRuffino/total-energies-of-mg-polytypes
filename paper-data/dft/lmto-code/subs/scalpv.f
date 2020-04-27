      subroutine scalpv(strx,nlf,mlf,nkap,alphv,SIDE)
C- Scales structure constant matrix by alphv
C ----------------------------------------------------------------
Ci Inputs
Ci   strx,nlf,nlf,nkap,alphv
Ci   strx: structure constants to scale
Ci   SIDE: if LEFT (ie 0) scale strx on lhs, otherwise scale on rhs
Co Outputs
Co   strx are scaled by alphv, either lhs or rhs.
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nkap,nlf,mlf,SIDE
      double precision alphv(nlf,1),strx(nlf,mlf)
C Local parameters
      double precision tmp(2,2)
      integer nskip,i,j
      integer LEFT,RIGHT
      parameter (LEFT=0,RIGHT=1)

      if (nkap == 1) then
        if (SIDE .EQ. LEFT) then
          do  i = 1, nlf
          do  j = 1, mlf
            strx(i,j) = alphv(i,1)*strx(i,j)
          enddo
          enddo
        else
          do  i = 1, nlf
          do  j = 1, mlf
            strx(i,j) = strx(i,j)*alphv(j,1)
          enddo
          enddo
        endif
      else
        nskip = nlf*mlf
        if (SIDE .EQ. LEFT) then
          do  i = 1, nlf
          do  j = 1, mlf
            call dmpy(alphv(i,1),nlf+nlf,nlf,
     .                strx(i,j),nskip+nskip,nskip,
     .                tmp,2,1,2,2,2)
            call dcopy(4,tmp,1,strx(i,j),nskip)
          enddo
          enddo
        else
          do  i = 1, nlf
          do  j = 1, mlf
            call dmpy(strx(i,j),nskip+nskip,nskip,
     .                alphv(j,1),nlf+nlf,nlf,
     .                tmp,2,1,2,2,2)
            call dcopy(4,tmp,1,strx(i,j),nskip)
          enddo
          enddo
        endif
      endif
      end
