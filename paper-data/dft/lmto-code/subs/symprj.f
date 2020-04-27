      subroutine symprj(nrclas,nlml,ngrp,nbas,istab,g,ag,plat,qlat,
     .  pos,sym)
C- Set up symmetry projectors for one class
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrclas:number of atoms in this class
Ci   nlml  :L-cutoff to which symmetry projectors are calculated
Ci   ngrp  :size of space group
Ci   nbas  :size of basis
Ci   istab :not used
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   plat  :primitive lattice vectors, in units of alat
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   pos   :basis vectors, in units of alat
Co Outputs
Co   sym   :symmetry projectors for each site within this class
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlml,nrclas,ngrp,nbas,istab(nbas,ngrp)
      double precision sym(nlml,nlml,nrclas),plat(3,3),qlat(3,3),
     .   g(3,3,ngrp),ag(3,ngrp),pos(3,nrclas),d(3)
C ... Local parameters
      integer lmxl,ll,ig,ja,m,ia,iprint,ilm,l,jlm,jlm1,jlm2
      double precision wgt
      real(8), allocatable :: rmat(:)

      lmxl = ll(nlml)
      allocate(rmat(nlml*nlml))
      wgt = 1d0/ngrp
      call dpzero(sym, nlml*nlml*nrclas)

C --- For each group operation, do ---
      do  ig = 1, ngrp

C ...   Find site mapped into first site under this operation
        do  ja = 1, nrclas
          do  m = 1, 3
            d(m) = g(m,1,ig)*pos(1,ja) + g(m,2,ig)*pos(2,ja)
     .           + g(m,3,ig)*pos(3,ja) + ag(m,ig) - pos(m,1)
          enddo
          call shorbz(d,d,plat,qlat)
          ia = ja
          if (d(1)*d(1)+d(2)*d(2)+d(3)*d(3) < 1d-9) goto 5
        enddo
        call rxi('symprj: no site mapped into first under op',ig)
    5   continue

C       print *, ig,ia,ja,istab(ja,ig),istab(ia,ig)

C ...   Make and add transformation matrix
        call ylmrtg(nlml,g(1,1,ig),rmat)

        call dpadd(sym(1,1,ia),rmat,1,nlml*nlml,wgt)
      enddo

      deallocate(rmat)
C ... Printout
      if(iprint() < 60) return
      do  ja = 1, nrclas
        write (6,1) ja
    1   format(/' projection matrices for ja=',i3)
        ilm = 0
        do  l = 0, lmxl
          jlm1 = l*l+1
          jlm2 = (l+1)**2
          do  m = 1, 2*l+1
            ilm = ilm+1
            write (6,2) (sym(ilm,jlm,ja),jlm = jlm1,jlm2)
    2       format(1x,9F12.6)
          enddo
        enddo
      enddo
      end
