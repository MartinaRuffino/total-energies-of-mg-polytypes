      subroutine sliny(volwgt,ec,emin,emax,dosi,nr)
C- Adds to number-of-states for one tetrahedron
C ----------------------------------------------------------------
Ci Inputs
Ci   volwgt, weight on tetrahedron; ec energies at corners of tethdn.;
Ci   emin, emax, energy window; nr, number of bins + 1
Co Outputs
Co   dosi(k), integrated density in kth bin from tethdn.
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nr
      double precision ec(4),dosi(nr),volwgt,emin,emax
C Local parameters
      integer i,i01,i02,i03,i04,i1,i2,i3,i4,j
      double precision c1,c2,c3,cc,de,e,e1,e2,e3,e4,x

      do  3  i = 1, 3
        do  3  j = 1, 4-i
        if (ec(j) > ec(j+1)) then
          e = ec(j)
          ec(j) = ec(j+1)
          ec(j+1) = e
        endif
    3 continue
      e1 = ec(1)
      e2 = ec(2)
      e3 = ec(3)
      e4 = ec(4)
      de = (emax-emin)/(nr-1)
      i01 = (e1-emin)/de + 1.9999999d0
      i02 = (e2-emin)/de + 1.9999999d0
      i03 = (e3-emin)/de + 1.9999999d0
      i04 = (e4-emin)/de + 1.9999999d0
C --------------------------------
      i1 = max0(i01,1)
      i2 = min0(i02-1,nr)
      if (i1 <= i2) then
        cc = volwgt*3.d0/((e2-e1)*(e3-e1)*(e4-e1))
        do  20  i = i1, i2
          x = emin - e1 + (i-1)*de
          dosi(i) = dosi(i) + cc*x**2
   20   continue
      endif
      i2 = max0(i02,1)
      i3 = min0(i03-1,nr)
      if (i2 <= i3) then
        c3 = volwgt*3.d0*(e1+e2-e3-e4)/((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2))
        c2 = volwgt*6.d0/((e3-e1)*(e4-e1))
        c1 = c2*(e2-e1)*0.5d0
        do  21  i = i2, i3
          x = emin - e2 + (i-1)*de
          dosi(i) = dosi(i) + c1 + x*(c2 + x*c3)
   21   continue
      endif
      i3 = max0(i03,1)
      i4 = min0(i04-1,nr)
      if (i3 <= i4) then
        cc = volwgt*3.d0/((e3-e4)*(e2-e4)*(e1-e4))
        do  22  i = i3, i4
          x = emin - e4 + (i-1)*de
          dosi(i) = dosi(i) - cc*x**2
   22   continue
      endif
      end
