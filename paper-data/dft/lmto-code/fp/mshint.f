      subroutine mshint(vol,n,n1,n2,n3,k1,k2,k3,c,sum1,sum2)
C- Integrate a function tabulated on a real-space mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   vol   :cell volume
Ci   n     :number functions
Ci n1,n2,n3:size of mesh
Ci k1,k2,k3:dimensions of array c
Ci    c    :mesh of points to integrate
Co Outputs
Co   sum1  :real part of integral
Co   sum2  :imaginary part of integral
Cr Remarks
Cu Updates
Cu   21 Apr 00 Adpated from nfp meshint
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,n1,n2,n3,k1,k2,k3
      double precision vol,sum1,sum2
      double complex c(k1,k2,k3,n)
C ... Local parameters
      integer stdo,lgunit,i,i1,i2,i3
      double precision bot1,bot2,top1,top2
      double complex csum

      stdo = lgunit(1)
      csum = 0d0
      top1 = -1d10
      top2 = -1d10
      bot1 = 1d10
      bot2 = 1d10
      do  i = 1, n
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              csum = csum + c(i1,i2,i3,i)
              top1 = dmax1(top1,dble(c(i1,i2,i3,i)))
              top2 = dmax1(top2,dimag(c(i1,i2,i3,i)))
              bot1 = dmin1(bot1,dble(c(i1,i2,i3,i)))
              bot2 = dmin1(bot2,dimag(c(i1,i2,i3,i)))
            enddo
          enddo
        enddo
      enddo

      csum = csum*vol/(n1*n2*n3)
      sum1 = dble(csum)
      sum2 = dimag(csum)

C      if (iprint() >= 50) write(stdo,1) bot1,top1,sum1,bot2,top2,sum2
C    1 format(/' mshint: real part: min,max',2f10.6,'    integral',f12.6
C     .  /'         imag part: min,max',2f10.6,'    integral',f12.6)


      end
