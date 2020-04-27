      subroutine invalp(nlb,iopt,alpha)
C- replace alpha by inverse of alpha
C ----------------------------------------------------------------
Ci Inputs
Ci   nlb:  nl*nbas, or nl**2*nbas, depending on how alpha is stored
Ci         (see gena2)
Ci   iopt: 1, single basis
Ci   iopt: 2: double basis, format as in gena2
Co Outputs
Co   alpha is replace by inverse
Cr Remarks
Cr   Alpha is assumed to be diagonal in R and L
C ----------------------------------------------------------------
      implicit none
      integer nlb,iopt
      double precision alpha(nlb,2,2)
C Local parameters
      integer i
      double precision tmp,det

      if (iopt == 1) then
        do  i = 1, nlb
          alpha(i,1,1) = 1/alpha(i,1,1)
        enddo
      else if (iopt == 2) then
        do  i = 1, nlb
          det = alpha(i,1,1)*alpha(i,2,2) - alpha(i,1,2)*alpha(i,2,1)
          tmp          = +alpha(i,1,1)/det
          alpha(i,1,2) = -alpha(i,1,2)/det
          alpha(i,2,1) = -alpha(i,2,1)/det
          alpha(i,1,1) = +alpha(i,2,2)/det
          alpha(i,2,2) =  tmp
        enddo
      else
        call fexit(-1,119, 'INVALP: bad iopt',0)
      endif
      end
