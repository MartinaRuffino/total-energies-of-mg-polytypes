      subroutine prmsk2(lmx1,lmx2,ipc,nm,nRL1,nRL2,s,fmt0,string,
     .  id1,id2,indxsh)
C- Print out structure constant and/or hamiltonian matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   lmx1,nRL1: lmx, nRL for rows;
Ci   lmx2,nRL2: lmx, nRL for columns.
Ci   s:      matrix to be printed
Ci   nm:     number of successive matrices to print
Ci           (2 for complex, 4 for real double-kappa, etc)
Ci   ipc: ipc(i) holds pointer to lmx for ith block.
Ci   string: header, printed first
Ci   id1,id2,indxsh are dummies for now
Co Outputs
Co   s is printed in blocks of lmx
Cr Remarks
Cr   adapted from prmsk; ipc reverted to old ipc
Cr   Needs to be adapted for downfolded, permuted matrices
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nm,nRL1,nRL2,indxsh(*)
      integer lmx1(*),lmx2(*),ipc(*)
      character*(*) string,fmt0
      double precision s(0:nRL1-1,0:nRL1-1)
C Local parameters
      double precision scal
      integer i1h,i1l,i2h,i2l,ib1,ib2,im,lm1,lm2,n2
      character*(20) fmt
      integer id1,id2
      data scal /1d0/

      fmt = '1x,9f8.4'
      if (fmt0 /= '?') fmt = fmt0
      print *
      print *, string
      n2 = 0
      do  im = 1, nm
        i1l = 0
        ib1 = 0
C --- Loop through the rows ---
    5   continue
        ib1 = ib1+1
          i1h = i1l - 1 + (lmx1(ipc(ib1))+1)**2
          i2l = 0
          ib2 = 0
C ---   Loop through the columns ---
   10   continue
        ib2 = ib2+1
        i2h = i2l - 1 + (lmx2(ipc(ib2))+1)**2
        print 332, ib1,ib2, i1l,i2l, i1h-i1l+1, i2h-i2l+1,
     .             n2/(nRL1*nRL2), nRL1, nRL2
  332   format('sites:', 2i3, '   block:', 2i4, '   size:', 2i3,
     .         '   offset:', i3,'   nRL: (',2i4, ')')
        do  lm1 = i1l, i1h
          write(*,'('//fmt//')') (scal*s(n2+lm1,lm2),lm2=i2l,i2h)
    2     format(1x,9F8.4)
        enddo
        i2l = i2h+1
        if (i2l < nRL2) goto 10
        i1l = i1h+1
        if (i1l < nRL1) goto 5
        n2 = n2 + nRL1*nRL2
      enddo
      end
