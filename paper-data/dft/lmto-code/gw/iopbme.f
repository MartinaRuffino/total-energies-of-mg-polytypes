      subroutine iopbme(mode,iat1,iat2,nl,ndrphi,npb,off,rprbme)
C- File I/O of radial product basis matrix elements
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0  file READ dimensioning parameters
Ci         :1  file READ rprbme
Ci         :2  file WRITE
Ci         :10s digit
Ci         :File format old GW style (one file per site)
Ci   iat1  :write matrix elements for sites in range (iat1,iat2)
Ci   iat2  :write matrix elements for sites in range (iat1,iat2)
Ci   nl    :(global maximum l) + 1
Ci   ndrphi:(dimensioning) max number of core + valence radial functions of a particular l
Ci         :Formerly named nn in old GW code
Ci   npb   :indexing for rprodb.
Ci         :Functions npb(l,iat):npb(l+1,iat) are the family of radial B functions
Ci         :stored in rprodb for quantum number l and site iat.
Ci   off   :reserved in future for offset to rprbme. Not used now
Cio Inputs/Outputs
Cio rprbme :radial matrix elements < gtoto gtoto | B> for sites iat1 ... iat2
Cr Remarks
Cr   It is the caller's responsibility to ensure that rprbme is large enough
Cr   Structure of old GW code files PPBRD_V2_nn, nn is site index
Cr   Records:
Cr     1. ntpbi (formerly nblocha or mdim), 2*(nl-1), nblr (formerly nxx)
Cr     2. rprbme for functions belonging to site nn
Cu Updates
Cu   05 Jul 18
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: mode,iat1,iat2,nl,ndrphi,off
      integer :: npb(0:2*(nl-1),iat2+1)
      double precision rprbme(0:nl-1,ndrphi,0:nl-1,ndrphi,0:2*(nl-1),*)
C ... Dynamically allocated local arrays
C ... Local parameters
      integer mode0,mode1,nblr(0:2*(nl-1)),iat,ntpbi,ifi,mxrbli,offpr
!     double precision x
      character*11 filenamep
      procedure(integer) :: fopng

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      if (mode0 < 2) call rx('iopbme not ready for file READ')

      call info0(30,1,0,'')

C --- File WRITE ---
      if (mode1 == 1) then ! old GW style, separate files named PPBRD_V2_nn
        offpr = 0
        do  iat = iat1, iat2
          call nrpb2nbl(03,nl,iat,npb,nblr,mxrbli,ntpbi)

          filenamep = 'PPBRD_V2_'//char(48+iat/10)//char(48+mod(iat,10))
          call info5(30,0,0,' iopbme writing file '//filenamep//' site%,3i  ntpbi=%,3i  nblr=%s,%ni  max=%i',
     .      iat,ntpbi,2*nl-2,nblr,mxrbli)
          ifi = fopng(filenamep,-1,4)
          write(ifi) ntpbi, 2*(nl-1), nblr
          write(ifi) rprbme(:,:,:,:,:,offpr+1:offpr+mxrbli)
          call fclr(' ',ifi)
          offpr = offpr + mxrbli
        enddo
      endif

      end
      subroutine nrpb2nbl(mode,nl,nat,nrpb,nblr,mxrbls,ntpb)
C- Construct nblr for one site from nrpb
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do nothing; just return
Ci         :1 return mxrbls,ntpb
Ci         :2 return nblr
Ci         :3 combination of 1+2
Ci         :10s digit
Ci         :0 nat refers to a single site.
Ci         :  ntpb and mxrbls are calculated for site nat only
Ci         :  nblr (if it is to be returned) is returned in nblr(:,1)
Ci         :1 nat refers to a collection of sites 1:nat
Ci         :  ntpb and mxrbls are summed over all sites
Ci         :  nblr (if it is to be returned) is returned in nblr(:,1:nat)
Ci   nl    :dimensions nrpb and nblr
Ci   nat   :extract parameters for this site
Ci   nrpb  :indexing for rprodb.
Ci         :Functions nrpb(l,nat):nrpb(l+1,nat) are the family of radial B functions
Ci         :stored in rprodb for quantum number l and site nat.
Ci         :In former GW code, nrpb(l+1,nat)-nrpb(l,nat) was stored in nxx
Co Outputs
Co   nblr  :number of radial product basis functions for each l
Co         :Formerly named nxx in old GW code
Co   mxrbls:sum of maxval(nblr) for all sites considered
Co   ntpb  :total number of local product basis functions within MT for sites considered
Co         :If for one site, formerly named nblocha in old GW code
Co         :If for all sites, formerly named nbloch in old GW code
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   07 Jul 18
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,nat,nrpb(0:2*(nl-1),nat+1),nblr(0:2*(nl-1),nat),mxrbls,ntpb
C ... Local parameters
      integer lx,iat,iat1,nrbli(0:2*(nl-1)),mxrbli,ntpbi,mode0

      iat1 = nat
      mode0 = mod(mode,10)
      if (mod(mode/10,10) > 0) iat1 = 1
      if (mod(mode0,2) > 0) then
        mxrbls = 0; ntpb = 0
      endif

      do  iat = iat1, nat

        forall (lx=0:2*nl-3) nrbli(lx) = nrpb(lx+1,iat) - nrpb(lx,iat)
        nrbli(2*(nl-1)) = nrpb(0,iat+1) - nrpb(2*(nl-1),iat) ! last element done separately so compiler does not fuss

        ntpbi = 0
        do  lx = 0, 2*(nl-1)
          ntpbi = ntpbi + nrbli(lx)*(2*lx+1)
        enddo
        mxrbli = maxval(nrbli)
        if (mod(mode0,2) > 0) then
          mxrbls = mxrbls + mxrbli ! To dimension global rprbme
          ntpb = ntpb + ntpbi
        endif

        if (mod(mode0/2,2) > 0 .and. mod(mode/10,10) > 0) then
          call icopy(2*nl-1,nrbli,1,nblr(1,iat),1)
        elseif (mod(mode0/2,2) > 0) then
          call icopy(2*nl-1,nrbli,1,nblr,1)
        endif

      enddo

      end
