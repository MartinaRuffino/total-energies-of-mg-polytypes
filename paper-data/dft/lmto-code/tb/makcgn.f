      subroutine makcgn(nlm1,nlm2,nlm3,indxcg,jcg,cg,cgn)
C- Make nxmxk matrix of "Clebsch Gordan" (actually Gaunt) coefficients
C ----------------------------------------------------------------------
Ci Inputs:
Ci nlm1,nlm2,nlm3: dimensions of cgn. Also the actual limits to make cgn
Ci indxcg,jcg,cg : CG (Gaunt) coefficients stored in condensed form
Co Outputs:
Co cgn    : Gaunt coefficients in Stone's convention stored as a table
Cr Remarks
Cr   The program is only used by tbfrc2.f. An alternative to the latter,
Cr   tbfrc3.f, doesn't need makcgn.
Cr
Cr   From demonstration program written by Michael.
Cb Bugs
Cb   At present, re-shuffling only works for nlm1 = 9.
Cu Updates
Cu   17 Feb 10 (SL) remake of makcg9.f for arbitrary nlm2 and nlm3.
Cu                  nlm1 still has to be 9 because of xxxinm (see Bugs).
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer, intent(in) :: nlm1,nlm2,nlm3
      integer, intent(in) :: jcg(*),indxcg(*)
      double precision, intent(in)  :: cg(*)
      double precision, intent(out) :: cgn(nlm1,nlm2,nlm3)
C Local Variables
      integer ii,ilm1,ilm2,ilm3,indx,icg,icg1,icg2,iprint

      cgn(1:nlm1,1:nlm2,1:nlm3) = 0d0

c ... loop over the first two indices
      do  ilm1 = 1, nlm1
        do ilm2 = 1, nlm2

c ...     Find top and bottom limits icg1,icg2 in the lists.
c         The min and max functions are used because cases like
c         (1,3) and (3,1) are only tabulated once.
c         What you normally want is
c            indx=(ilm1*(ilm1-1))/2+ilm2
c         where it is known that ilm1 is larger or equal to ilm2.
c         The max, min stuff always gets the proper icg1 and icg2,
c         no matter which one of ilm1 and ilm2 is larger.

          ii = max0(ilm1,ilm2)
          indx = (ii*(ii-1))/2 + min0(ilm1,ilm2)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1

c ...     loop over the relevant part of list, get out
c         the coefficient and the third index.
c         Note that ilm3 will run to higher values then either
c         ilm1 and ilm2 (twice the lmax value).
c         If you want ilm3 to stay below nlm also, you need
c         extra if..else  statements.

          do  icg = icg1, icg2
            ilm3 = jcg(icg)
            if (ilm3 <= nlm3) cgn(ilm1,ilm2,ilm3) = cg(icg)
          enddo
        enddo
      enddo

      if (nlm1 == 9) then
        call xxxinm(nlm2,nlm3,cgn)
      else
        call rxi(' makcgn: first dimension of cgn should be 9. nlm1 =',
     .    nlm1)
      endif

      if (iprint() >= 120) then
        print *, ' makcgn: Gaunt coefficients ...'
        do  ilm3 = 1, nlm3
          print *,' ilm3 = ',ilm3
          write (*,10) ((cgn(ilm1,ilm2,ilm3),ilm2=1,nlm2),ilm1=1,nlm1)
          print *
        enddo
   10   format (25f10.6)
      endif

      end subroutine makcgn

      subroutine xxxinm(n,m,cgn)
C- Shuffle 9 x n x m matrix from MSM lm indices to TBE lm indices
C ----------------------------------------------------------------------
Ci Inputs:
Ci   cgn : 9 x n x m matrix of Clebsch Gordans
Co Outputs:
Co   cgn reordered
Cr Remarks
Cr   Michael's structure constants are ordered according to the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     y     z     x    xy    yz 3z^2-r^2 zx   x^2-y^2
Cr   while the TBE programs use the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     x     y     z    xy    yz    zx  x^2-y^2  3z^2-r^2
Cr   This widget rearranges the matrix made by makcg9 into
Cr   the TBE order.
Cr   The l>2 ordering is unchanged
Cu Updates
Cu   17 Feb 10 (SL) modified from xxxxin (makcg9.f) wrt dimensions
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer, intent(in) :: n, m
      double precision, intent(inout) :: cgn(9,n,m)
C Local Variables
      double precision wk(9,n,m)
      integer i, j, k, ind(9), nn, mm
      data ind /1,4,2,3,5,6,8,9,7/

      nn = min(n,9)
      mm = min(m,9)
c... upper left corner
        forall (i=1:9,j=1:nn,k=1:mm)
     .    wk(i,j,k) = cgn(ind(i),ind(j),ind(k))
        cgn(1:9,1:nn,1:mm) = wk(1:9,1:nn,1:mm)
c... remaining columns
        if (n > 9) then
          forall (i=1:9,j=10:n,k=1:mm)
     .      wk(i,j,k) = cgn(ind(i),j,ind(k))
          cgn(1:9,10:n,1:mm) = wk(1:9,10:n,1:mm)
        endif
c... remaining rows
        if (m > 9) then
          forall ( i=1:9, j=1:nn, k=10:m)
     .      wk(i,j,k) = cgn(ind(i),ind(j),k)
          cgn(1:9,1:nn,10:m) = wk(1:9,1:nn,10:m)
        endif
c... and finally
        if (m > 9 .and. n > 9) then
          forall ( i=1:9, j=10:n, k=10:m)
     .      wk(i,j,k) = cgn(ind(i),j,k)
          cgn(1:9,10:n,10:m) = wk(1:9,10:n,10:m)
        endif

      end subroutine xxxinm
