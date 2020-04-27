      subroutine symtab(nbas,iclass,bas,g,ag,ng,rb,qb,iwk,istab)
C- Make symmetry transformation table for basis atoms; check classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas:  number of atoms in basis
Ci   iclass:  jth atom belongs to class iclass(j)
Ci   bas:     bas(i,j) are Cartesian coordinates of jth atom in basis
Ci   g,ag,ng: defines space group (see sgroup.f)
Ci            NB: sign of ng used as a flag to suppress checking of
Ci            atom classes.
Ci   rb:  real space lattice vectors
Ci   qb:  reciprocal space lattice vectors
Ci   iwk: work array.  Not used if input ng<0
Co Outputs
Co   istab:  site ib is transformed into istab(ib,ig) by grp op ig
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nbas,ng
      integer iclass(*),iwk(*),istab(nbas,*)
      double precision bas(3,1),g(9,1),ag(3,1),rb(9),qb(9)
C Local variables
      double precision dlat(3)
      integer ib,ic,ig,jb,jc,stdo,nglob

      if (ng == 0) return
      stdo = nglob('stdo')

C --- Make atom transformation table ---
      do  ig = 1, iabs(ng)
        do  ib = 1, nbas
C         call gpfndx(g(1,ig),ag(1,ig),ib,jb,bas,nbas,rb,qb)
          call grpfnd(1d-5,g,ag,ig,bas,nbas,qb,ib,jb,dlat)
          if (jb == 0) then
            write(stdo,100) ib,ig
            call rx('SYMTAB: bad group operations')
          endif
          if (ng > 0) then
            ic = iclass(ib)
            jc = iclass(jb)
            if (ic /= jc) then
              write(stdo,110) ib,jb,ig
              call rx('SYMTAB: invalid atom classes')
            endif
          endif
          istab(ib,ig) = jb
        enddo
      enddo

  100 format(/' ***ERROR*** Map of atom ib=',i4,' not found for ig=',i3)
  110 format(/' ***ERROR*** Atom ib=',i4,
     .  ' not in same class as mapped atom jb=',i4,' for ig=',i3)

      if (ng < 0) return

C --- Check atom classes ---
      do  ib = 1, nbas
        ic = iclass(ib)
        call iinit(iwk,nbas)
        do  ig = 1, ng
          iwk(istab(ib,ig)) = 1
        enddo
        do  jb = 1, nbas
          if (iwk(jb) == 1) cycle
          jc = iclass(jb)
          if (ic == jc) then
            write(stdo,120) ib,jb
            call rx('SYMTAB: bad atom classes')
          endif
        enddo
      enddo

  120 format(/' ***ERROR*** Atom ib=',i4,
     .  ' in same class as symmetry-inequivalent atom jb=',i4)

      end
