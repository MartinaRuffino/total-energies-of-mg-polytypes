      subroutine gvsym(ng,gv,ips0,bgv,c,csym)
C- Symmetrize a function c, given in the form of a list
C ----------------------------------------------------------------------
Ci Inputs
Ci   gv,ng   Lattice vectors, and number
Ci   ips0    pointer to first vector in star of this vector; see sgvsym
Ci           ips0<0 => star not complete (skip)
Ci   bgv     phase factor sum; see sgvsym
Ci   c       unsymmetrized function
Co Outputs:
Co   csym    symmetrized function
Cr Remarks:
Cu   24 Jun 09 Skip vectors i when ips0(i)<0
Cu    7 Sep 98 Adapted from nfp gvsym.f
C ----------------------------------------------------------------------
      implicit none
      integer ng,ips0(ng)
      double precision gv(ng,3)
      double complex bgv(ng),c(ng),csym(ng)
      integer i,j,i0,kstar,ipr,iprint,nstar

C ... Sum up coefficients for first vector in each star
      do  i = 1, ng
        csym(i) = (0d0,0d0)
      enddo
      do  i = 1, ng
        j = ips0(i)
        if (j > 0) csym(j) = csym(j) + bgv(i)*c(i)
      enddo

C ... Normalize
      do  i0 = 1, ng
        if (ips0(i0) == i0) then
          kstar = 0
          do  i = i0, ng
            if (ips0(i) == i0) kstar = kstar+1
          enddo
          if (j > 0) csym(i0) = csym(i0)/kstar
        endif
      enddo

C ... Make all the coefficients
      do  i = 1, ng
        j = ips0(i)
        csym(i) = csym(j)*dconjg(bgv(i))
      enddo

C ... Printout
      ipr = iprint()
      if (ipr < 55) return
      print 255
      nstar = 0
      do  i0 = 1, ng
        if (ips0(i0) == i0) then
          nstar = nstar+1
          if (ipr >= 60) print *, ' '
          do  i = i0, ng
            if (ips0(i) == i0) then
              if (i == i0) then
                print 251, nstar,i,gv(i,1),gv(i,2),gv(i,3),
     .             c(i),csym(i)
              else
                if (ipr >= 60)
     .             print 252, i,gv(i,1),gv(i,2),gv(i,3),
     .             c(i),csym(i)
              endif
  251         format(i4,i5,3f6.1,2f12.8,1x,2f12.8)
  252         format(4x,i5,3f6.1,2f12.8,1x,2f12.8)
  255         format(/' star  ig',8x,'recip',17x,'c_in',20x,'c_sym')
            endif
          enddo
        endif
      enddo

      end
