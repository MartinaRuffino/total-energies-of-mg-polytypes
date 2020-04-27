      subroutine pairsx(nbas,nbasp,alat,plat,rmax,baspp,ipsp,nd,iltab,
     .  pltab,nttab,ntab,oiax,mxcsiz)
C- Allocate memory for and create neighbor table for a crystal
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis (input)
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   rmax  :maximum range for connecting vector, in a.u.
Ci          All connecting vectors with length < rmax(i)+rmax(j)
Ci          are retained.  rmax may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ipsp(1);
Ci          see description of ipsp
Ci   baspp :basis vectors, doubly padded for planar geometry
Ci   ipsp  :index to which species each site belongs, for padded basis;
Ci          identifies which rmax is associated with each site. NB:
Ci          ipsp(1) = -1 => rmax is a global scalar, independent of site
Ci          ipsp(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ipsp is not used.
Ci   nd    :number of dimensions for which periodic boundary conditions
Ci          are used
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci   mxcsiz:if nonzero, use in place of internal formula for mxnbr
Co Outputs
Co   nttab   :total number of pairs in neighbor table and iax
Co   ntab    :ntab array; see pairc below where it is generated
Co   w(oiax) :iax array; see pairc, below where it is generated
Co   mxcsiz  :size of the largest cluster encountered
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
      integer nbas,nbasp,nttab,oiax,nd,iltab
      integer ntab(nbasp+1),ipsp(*),pltab(*)
      double precision alat,plat(9),rmax(*),baspp(3,1)
C ... Dynamically allocated local arrays
      integer, allocatable :: iaxl(:)
      real(8),allocatable :: wk(:)
C Local variables
      integer modep(3),nbaspp,mxcsiz,mxnbr,i,niax,isw,paird
      parameter (niax=10)
C heap:
      integer w(1)
      common /w/ w

C ... Set up input for call to pairc
      if (mxcsiz > 0) then
        mxnbr = mxcsiz
      else
        mxnbr = paird(nbas,nbasp,alat,plat,rmax,mxcsiz)
      endif
      do  i = 1, 3
        modep(i) = 2
        if (i > nd) modep(i) = 0
      enddo
      nbaspp = 2*nbasp - nbas
      nttab = mxnbr
      isw = 0

C ... Allocate iax array and make the neighbor table
      allocate(iaxl(niax*mxnbr)); call iinit(iaxl,niax*mxnbr)
      allocate(wk(3*mxnbr))
      call pairc(1,nbasp,nbaspp,modep,isw,ipsp,alat,plat,baspp,baspp,
     .  rmax,iltab,pltab,nttab,ntab,iaxl,wk,mxcsiz)
      call defi(oiax, niax*nttab)
      call icopy(niax*nttab,iaxl,1,w(oiax),1)
      deallocate(wk,iaxl)
      end
