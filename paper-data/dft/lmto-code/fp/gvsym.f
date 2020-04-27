      subroutine gvsym(ng,gv,ips0,bgv,c,csym)
C- Symmetrize a function c, given in the form of a list of plane waves
C ----------------------------------------------------------------------
Ci Inputs
Ci   gv,ng   Lattice vectors, and number
Ci   ips0    pointer to first vector in star of this vector; see sgvsym
Ci   bgv     phase factor sum; see sgvsym
Ci   c       unsymmetrized function
Co Outputs:
Co   csym    symmetrized function
Cr Remarks:
Cu    7 Sep 98 Adapted from nfp gvsym.f
C ----------------------------------------------------------------------
      implicit none
      integer ng,ips0(ng)
      double precision gv(ng,3)
      double complex bgv(ng),c(ng),csym(ng)
      integer i,j,i0,kstar,ipr,iprint,nstar

      call dpzero(csym,2*ng)

C ... Sum up coefficients for first vector in each star
      do  i = 1, ng
        j = ips0(i)
        csym(j) = csym(j) + bgv(i)*c(i)
      enddo

C ... Normalize : divide by the number of elements in the star
      do  i0 = 1, ng
        if (ips0(i0) == i0) then
          kstar = 0
          do  i = i0, ng
            if (ips0(i) == i0) kstar = kstar+1
          enddo
          csym(i0) = csym(i0)/kstar
        endif
      enddo

C ... Make all the coefficients
      do  i = 1, ng
        j = ips0(i)
        csym(i) = csym(j)*dconjg(bgv(i))
      enddo

C ... Printout
      ipr = iprint()
      if (ipr < 60) return
      print 255
      nstar = 0
      do  i0 = 1, ng
        if (ips0(i0) == i0) then
          nstar = nstar+1
          if (ipr >= 61) print *, ' '
          do  i = i0, ng
            if (ips0(i) == i0) then
              if (i == i0) then
                print 251, nstar,i,gv(i,1),gv(i,2),gv(i,3),c(i),csym(i)
              else if (ipr >= 61)  then
                print 252, i,gv(i,1),gv(i,2),gv(i,3),c(i),csym(i)
              endif
  251         format(i4,i5,3f6.1,2f12.8,1x,2f12.8)
  252         format(4x,i5,3f6.1,2f12.8,1x,2f12.8)
  255         format(/' star  ig',8x,'recip',17x,'c_in',20x,'c_sym')
            endif
          enddo
        endif
      enddo

      end

      subroutine gvsymafm(ng,ips0,bgv,c)
C- Special AFM symmetrization of a function c, given in the form of a list of plane waves
C ----------------------------------------------------------------------
Ci Inputs
Ci   gv,ng   Lattice vectors, and number
Ci   ips0    pointer to first vector in star of this vector; see sgvsym
Ci   bgv     phase factor sum; see sgvsym
Ci   c       unsymmetrized function
Co Outputs:
Co   csym    symmetrized function
Cr Remarks:
Cu  7 Sep 98 This routine seems to work, but more checks are needed.
Cu           It assumes that there are either no G vectors in the star of G
Cu           or there are two that come in pairs, for any G
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,ips0(ng)
      double complex bgv(ng),c(ng,2)
C ... Local parameters
      integer i,j
      double complex csum(2,2),cdif(2,2)

C ... Sum up coefficients for first vector in each star
      do  i = 1, ng
        j = ips0(i)
        csum(1,1) = c(i,1)+c(i,2)
        csum(2,1) = c(j,1)+c(j,2)
        cdif(1,1) = c(i,1)-c(i,2)
        cdif(2,1) = c(j,1)-c(j,2)

!       print 333, i,j, c(i,1),c(i,2),c(j,1),c(j,2)

! Re + Im Consistent at least for FeSe and FePt cases
        if (j == i) then
!           print 333, i,j, c(i,1),c(i,2),c(i,2)-dconjg(c(j,1))*bgv(i),dble(bgv(i))/100,dble(bgv(j))/100

          csum(1,2) = (csum(1,1) + csum(2,1))/2
          csum(2,2) = (csum(2,1) + csum(1,1))/2
          cdif(1,2) = (cdif(1,1) - cdif(2,1)*bgv(i))/2
          cdif(2,2) = (cdif(2,1) - cdif(1,1)*bgv(i))/2
!           print 333, i,j, csum(:,1), cdif(:,1), dble(bgv(i))/100
!           print 333, i,j, csum(:,2), cdif(:,2), dble(bgv(i))/100
C            csym(i,1) = (c(i,2)+dconjg(c(j,1))*bgv(i))/2
C            csym(i,2) = (c(i,2)+dconjg(c(j,1))*bgv(i))/2
        else

!            print 333, i,j, c(i,1),c(i,2),c(j,1),c(j,2)

          csum(1,2) = (csum(1,1) + csum(2,1)*bgv(i))/2
          csum(2,2) = (csum(2,1) + csum(1,1)*bgv(i))/2
          cdif(1,2) = (cdif(1,1) - cdif(2,1)*bgv(i))/2
          cdif(2,2) = (cdif(2,1) - cdif(1,1)*bgv(i))/2

!            print 333, i,j, c(i,1),c(i,2),c(i,2)-c(j,1)*bgv(i),dble(bgv(i))/100
!            print 333, i,j, csum(:,1), cdif(:,1), dble(bgv(i))/100
!             print 333, i,j, csum(:,2), cdif(:,2), dble(bgv(i))/100

          endif
          c(i,1) = (csum(1,2) + cdif(1,2))/2
          c(i,2) = (csum(1,2) - cdif(1,2))/2
          c(j,1) = (csum(2,2) + cdif(2,2))/2
          c(j,2) = (csum(2,2) - cdif(2,2))/2
!          print 333, i,j, c(i,1),c(i,2),c(j,1),c(j,2)

  333     format(2i5,6(2x,2p,2f12.6))

C          if (i == 100) exit

        enddo

      end
