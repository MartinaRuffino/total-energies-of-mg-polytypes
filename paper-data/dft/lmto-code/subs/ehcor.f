      subroutine ehcor(nclass,nrc,dclabl,vdiff,nl,nsp,qold,qnu,de)
C- Term (V(nin)-V)(nout-nin)/2 for shifts in sphere potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   nrc   :number of atoms per class
Ci   dclabl:class name, packed as a real number
Ci   vdiff :(V(nin)-V)
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   qold  :input moments
Ci   qnu   :output energy-weighted moments of the sphere charges
Co Outputs:
Ce   de    :(V(nin)-V)(nout-nin)/2
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nrc(*),nl,nsp,nclass
      double precision vdiff(*),qold(3,nl*nsp,1),qnu(3,nl*nsp,1),de
      double precision dclabl(nclass)
C Local Variables
      integer ic,i,iprint
      double precision xx,dq
      character*8 clabl

      if (iprint() >= 40) print 1
    1 format(/' ehcor: class     dq         dv         de')
      de = 0
      do  ic = 1, nclass
        dq = 0
        do  i = 1, nl*nsp
          dq = dq+qnu(1,i,ic)-qold(1,i,ic)
        enddo
        xx = dq*vdiff(ic)/2
        if (iprint() >= 40) then
          call r8tos8(dclabl(ic),clabl)
          print 2,clabl,dq,vdiff(ic),xx
    2     format(9x,a4,3F11.6)
        endif
        de = de + xx*nrc(ic)
      enddo

      end
