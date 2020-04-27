      integer function comfac(ifac,nfac,N,nn)
C- Reduce vector of integers n by common factors ifac
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifac  :vector of factors to check
Ci   nfac  :number of elements ifac
Ci   N     :Numbers to seek common factor
Ci   nn    :number of elements N
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   13 Jun 16
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nfac,nn,ifac(nfac),N(nn)
C ... Local parameters
      integer jfac,i,j,k

      k = 1
      do  j = 1, nfac
        jfac = ifac(j)
        do while (.true.)
          do  i = 1, nn
            if (jfac*(N(i)/jfac) /= N(i)) goto 10
          enddo
          k = k*jfac
          do  i = 1, nn
            N(i) = N(i)/jfac
          enddo
        enddo
   10   continue
      enddo

      comfac = k
      end
C     test comfac
C      implicit none
C      integer iv(4),comfac,i
C
C      iv(1:4) = [12*5,16*5,20,180]
C      print *, iv
C      i = comfac([2,3,5],3,iv,4)
C      print *,i
C      print *, iv
C      print *, i*iv
C      end
