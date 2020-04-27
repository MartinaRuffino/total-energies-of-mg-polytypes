      subroutine grpfnd(tol,g,ag,ig,pos,nbas,qlat,ia,ja,dlat)
C- Find index to site ja site into which g,ag transforms site ia
C ----------------------------------------------------------------------
Ci Inputs
Ci   tol   :tolerance in site positions
Ci   g     :rotation part of space group
Ci   ag    :translation part of space group
Ci   ig    :which group operation
Ci   pos   :basis vectors
Ci   nbas  :size of basis
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   ia    :site for which to find equivalent by group op (g,ag)
Ci         :If ia<0, absolute value of ia is used and additionally
Ci         :point group operation -g is used.
Co Outputs
Co   ja    :site that ia is transformed into by (g,ag)
Co         :i.e. R(ja) = g(ig) R(ia) + ag(ig) + dlat
Co         :where dlat is some lattice translation vector
Co         :if zero, no equivalent site was found.
Co   dlat  :lattice translation vector R(ja)-R(ia)
Cu Updates
Cu   15 Mar 11 Returns dlat
Cu   26 Jan 01 Add ability to operate with -g (ia<0)
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer ia,ja,ig
      double precision g(3,3,ig),ag(3,ig),pos(3,1),qlat(3,3),tol,dlat(3)
C Local parameters
      double precision d(3)
      logical latvec
      integer ka,nbas,m,k

      ka = iabs(ia)
      do  m = 1, 3
      d(m) = ag(m,ig)
        do  k = 1, 3
          d(m) = d(m) + g(m,k,ig)*pos(k,ka)
        enddo
      enddo
      if (ia < 0) call dscal(3,-1d0,d,1)

      ja = 0
      do  ka = 1, nbas
        dlat(1) = d(1) - pos(1,ka)
        dlat(2) = d(2) - pos(2,ka)
        dlat(3) = d(3) - pos(3,ka)
        if (latvec(1,tol,qlat,dlat)) then
          ja = ka
          return
        endif
      enddo
      end

