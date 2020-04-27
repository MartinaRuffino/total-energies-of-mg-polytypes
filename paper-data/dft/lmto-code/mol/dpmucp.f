      subroutine dpmucp(a,i1,i2,j1,j2,b,ldb)
C- Copies subblock of packed matrix to a normal matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   a is the input matrix (packed)
Ci   i1,i2,j1,j2: subblock to be copied
Co Outputs:
Ci   b,ldb is the subblock matrix and the row dimension of b
Cr Remarks:
C ----------------------------------------------------------------
C Passed Parameters
      implicit none
      integer ldb,ldc,nr,nc,l,i1,i2,j1,j2
      double precision a(1), b(ldb,ldb)
C Local parameters
      double precision sum
      integer i,j,k,offa,offb,jtop,ia,ja

      do  20  ia = i1, i2
        offa = (ia*(ia-1))/2
        jtop = min(ia,j2)
        do  22  ja = j1, jtop
C         PRINT *, I1,I2,J1,J2, IA-I1+1, JA-J1+1, '*', a(ja+offa)
   22   b(ia-i1+1,ja-j1+1) = a(ja+offa)
        do  24  ja = max(ia+1,j1), j2
C          PRINT *, I1,I2,J1,J2, ia,ja, IA-I1+1, JA-J1+1, '**',
C     .      a(ia+(ja*(ja-1))/2)
   24   b(ia-i1+1,ja-j1+1) = a(ia+(ja*(ja-1))/2)
   20 continue

      end
