      subroutine getmom(nsp,nl,qnu,nclass,nrclas,idxdn,zval)
C- Get total number of valence electrons from moments
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp,nl,nclass,nrclas
Ci   qnu(1,l,isp,ic) - zeroth moments for each L, spin, and class
Ci   idxdn(l,ic) - downfolding switches (if > 1, exclude orbitals)
Co Outputs
Co   zval: total number of valence electrons
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nsp,nl,nclass,nrclas(nclass),idxdn(nl,nclass)
      double precision zval,qnu(3,nl,nsp,1)
C Local variables
      integer l,ic,isp
      double precision d1mach

      zval = 0
      do  ic = 1, nclass
        do  isp = 1, nsp
          do  l = 1, nl
            if (idxdn(l,ic) <= 1) then
            zval = zval + nrclas(ic)*qnu(1,l,isp,ic)
            endif
          enddo
        enddo
      enddo

C      print 10, zval
C      10 format('GETMOM: zval = ',f10.6)

      if (zval < 10*d1mach(3)) zval = 0d0

      end
