      subroutine mksod(ldim,lihdim,ndim,indxsh,pph,soph,sod)
C- Diagonal arrays for one-, two- and three- center terms in the SO ham
C ---------------------------------------------------------------------
Ci Inputs
Ci   ldim  :lower dimension of hamiltonian matrix (makidx.f)
Ci   lihdim:dimensions pph
Ci   ndim  :dimensions of indxsh and soph; total no. of orbitals.
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :pph(1..5,i,is): parms for ith RL and is(th) spin channel.
Ci         :pph(1) : enu
Ci         :pph(2) : calpha
Ci         :pph(3) : sqrdel
Ci         :pph(4) : palpha
Ci         :pph(5) : oalp
Ci
Ci   soph  :Spin-orbit parameters, in downfolding order (mksoph.f)
Ci         :soph(*,*,*,i) are spin-orbit parameters for ith RL channel (mksoph.f)
Ci         :soph(1,is1,is2,i) <phi Hso phi> for spin indices s1,s2
Ci         :soph(2,is1,is2,i) <phi Hso dot> for spin indices s1,s2
Ci         :soph(3,is1,is2,i) <dot Hso dot> for spin indices s1,s2
Ci         :soph(4:6,:,:,i) are not used by this routine
Ci         :soph(7,1,1,i) l quantum number
Ci         :soph(7,2,1,i) m quantum number.  m is ordered -l, -l+1, ... l-1, l
Co Outputs
Co   sod(*,1:3,isp): arrays for 1-, 2-, 3- center terms, ++,-- blocks
Co   sod(*,1,3):    for 1- center terms, L- S+ block
Co   sod(*,2:3,3):  for 2- center terms, L- S+ block
Co   sod(*,4,3):    for 3- center terms, L- S+ block
Cr Remarks
Cr   Makes diagonal matrices belonging to the 1-, 2- and 3-center
Cr   terms in the matrix elements of the SO hamiltonian.  Its form is
Cr   SO_{RL,RL'} = <chi_{RL}| I L . S | chi_{RL'}>.
Cr   where L.S = Lz.Sz + .5(L+.S- + L-.S+) and I = (2/(c^2) dV/dr*(1/r).
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | chi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   The first form uses phidot = phidot^alpha; the second form uses
Cr     phidot^alpha = phidot^gamma + o phi; and calls 'dot'=phidot^gamma
Cr
Cr   This subroutine combines radial matrix elements I between phi's
Cr   and dot's (phidot^gamma) and the L.S operator into diagonal
Cr   matrices associated with the 1-, 2- and 3-center terms comprising
Cr   matrix elements (H_so)_(i,j).
Cr   For the (L-S+) block, the SO coupling uses:
Cr     L+ |lm> = sqrt((l-m)*(l+m+1)) |l,m+1>
Cr     L- |lm> = sqrt((l+m)*(l-m+1)) |l,m-1>
Cr     <phi1_l,m | L-.S+ | phi2_l,m+1> = F/2 < phi1_l | I | phi2_l>
Cr   where F = sqrt((l+m+1)*(l-m))
Cr   L-.S+ couples l,m only to l,m+1 (2l elements for each l).  The
Cr   following evaluates (2l+1) elements foreach l, exploiting the
Cr   fact that <lm   | L-S+ | lm+1> is zero for m=l because F=0, and
Cr   similarly <lm-1 | L-S+ | lm> is zero for m=-l.
Cr
Cr   Matrix elements have 1-, 2- and 3- center terms composed of
Cr   matrix elements of partial waves inside augmentation spheres
Cr     xi1 = <phi|SO|phi>  xi2 = <phi|SO|dot>  xi3 = <dot|SO|dot>
Cr   <chi | SO | chi> =
Cr     (1+ho)<phi_L| + h<dot_L| SO  | phi_L>(1+oh + |dot_L> h =
Cr     (1+ho) xi1 (1+oh) + (1+ho) xi2 h + cc + h xi3 h
Cr   Divide into 1-, 2-, 3- center terms
Cr     w1      where w1 = xi1
Cr     h w2    where w2 = (o xi1 + xi2)
Cr     h w3 h  where w3 = (o^2 xi1 + 2 o xi2 + xi3)
Cr
Cr   To work directly with  D = srdel S srdel, rather than
Cr   h = C-enu + D, the diagonal parts w1 connecting C-enu
Cr   in the one-, two-, three-center terms are reshuffled.  Thus
Cr     w1 + w2 (C-e+D) + h.c. + (C-e+D) w3 (C-e+D)  becomes
Cr     w1+2*w2*(C-e)+w3(C-e)**2  +  (w2+(C-e)w3) D + h.c.  +  D w3 D.
Cr   For the spin-diagonal parts, terms proportional to D^0, D^1 and D^3
Cr   are stored in sod(:,1:3,isp)
Cu Updates
Cu   16 Aug 13 Bug fix initializing some elements of sod
Cu   08 May 13 Initialize sod (some elements were never referenced)
Cu   05 Jan 04 Fix upper limit of loop in downfolding case
Cu   14 Feb 03 Dimensions of sod, soph altered to accomodate B-field
C ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldim,lihdim,ndim,indxsh(ndim)
      double precision pph(5,lihdim,2),soph(7,2,2,ndim),sod(ldim,4,3)
C ... Local parameters
      double precision sp,sz,flmp,flmm
      parameter (sp=.5d0)
      integer i,j,k,isp,m
      double precision oalp,oalp2
C     flmp = coff for L+(m) = sqrt[(l-m)(l+m+1)] = L-(-m).  Note also L+(m-1) = L-(m)
C     flmm = coff for L-(m) = sqrt[(l+m)(l-m+1)]
      flmp(i) = dsqrt((soph(7,1,1,i)+soph(7,2,1,i)+1)*(soph(7,1,1,i)-soph(7,2,1,i)))
      flmm(i) = dsqrt((soph(7,1,1,i)+soph(7,2,1,i))*(soph(7,1,1,i)-soph(7,2,1,i)+1))

      call dpzero(sod,ldim*4*3)

C --- Components for up-up, down-down blocks of S-O hamiltonian ------
      do isp = 1, 2
        sz = 1.5d0 - isp

C   ... For 3C part  ho <phi Hso phi> oh + 2 ho <phi Hso dot> h + h <dot Hso dot> h
        do k = 1, ndim
          i = indxsh(k)
          if (i <= ldim) then
            oalp = pph(5,i,isp)
            m = soph(7,2,1,i)
            sod(i,3,isp) = m*sz*soph(1,isp,isp,i)*oalp**2 +
     .                   2*m*sz*soph(2,isp,isp,i)*oalp +
     .                     m*sz*soph(3,isp,isp,i)
          endif
        enddo

C   ... For 2C part  <phi Hso phi> oh + <phi Hso dot> h + h.c.
        do k = 1, ndim
          i = indxsh(k)
          if (i <= ldim) then
            m = soph(7,2,1,i)
            oalp = pph(5,i,isp)
            sod(i,2,isp) = m*sz*(soph(2,isp,isp,i) + soph(1,isp,isp,i)*oalp)
          endif
        enddo

C   ... For 1C  <phi Hso phi>
        do k = 1, ndim
          i = indxsh(k)
          if (i <= ldim) then
            m = soph(7,2,1,i)
            sod(i,1,isp) = m*sz*soph(1,isp,isp,i)
          endif
        enddo

C   ... Modify 1C term: add diagonal (C-enu) subtracted from 2- and 3-center
C       (See "working directly with D" above)
        do k = 1, ndim
          i = indxsh(k)
          if (i <= ldim) then
            sod(i,1,isp) = sod(i,1,isp) +
     .      2*sod(i,2,isp)*(pph(2,i,isp)-pph(1,i,isp)) +
     .      sod(i,3,isp)*(pph(2,i,isp)-pph(1,i,isp))**2
          endif
        enddo
        do k = 1, ndim
          i = indxsh(k)
          if (i <= ldim) then
            sod(i,2,isp) = sod(i,2,isp) +
     .      sod(i,3,isp)*(pph(2,i,isp)-pph(1,i,isp))
          endif
        enddo
      enddo

C --- Components for up-down blocks of S-O hamiltonian ------

C ... For 3C contribution to L-S+ block
      do k = 1, ndim-1
        i = indxsh(k)
        if (i < ldim) then
          oalp = pph(5,i,1)
          oalp2 = pph(5,i+1,2) ! Nonsensical if m=l but ok since flmp(i) = 0
          sod(i,4,3) = sp*flmp(i) * (oalp*(soph(2,1,2,i) + oalp2*soph(1,1,2,i))
     .                            + soph(2,2,1,i)*oalp2 + soph(3,1,2,i))
        endif
      enddo
      sod(ndim,4,3) = 0

C ... For 2C contribution to L-S+ block
      do k = 1, ndim-1
      i = indxsh(k)
      if (i <= ldim) then
        oalp = pph(5,i,1)
        sod(i,2,3) = sp*flmp(i)*(soph(2,2,1,i) + oalp*soph(1,1,2,i))
      endif
      enddo
      do k = 2, ndim
      i = indxsh(k)
      if (i <= ldim) then
        oalp = pph(5,i,2)
        sod(i,3,3) = sp*flmm(i)*(soph(2,1,2,i) + oalp*soph(1,1,2,i))
      endif
      enddo
      sod(ldim,2,3) = 0
      sod(1,3,3) = 0

C ... For 1C contribution to L-S+ block
      do k = 1, ndim-1
        i = indxsh(k)
        if (i <= ldim) then
          sod(i,1,3) = soph(1,1,2,i)*flmp(i)*sp
        endif
      enddo
      sod(ldim,1,3) = 0

C ... Add diagonal (C-enu) subtracted from 2- and 3-center
      if (.true.) then
      do k = 1, ndim-1
        i = indxsh(k)
!       Nonsensical when m=l but ok since sod(i,1,3) should be zero there
        if (i < ldim) then
          sod(i,1,3) = sod(i,1,3) +
     .    (pph(2,i+1,1)-pph(1,i+1,1))*sod(i,2,3) +
     .    (pph(2,i+1,2)-pph(1,i+1,2))*sod(i+1,3,3) +
     .    (pph(2,i,1)-pph(1,i,1))*(pph(2,i+1,2)-pph(1,i+1,2))*sod(i,4,3)
        endif
      enddo

      do k = 1, ndim-1
        i = indxsh(k)
!       Nonsensical if m=l but ok since sod(i,4,3) should be zero there
        if (i < ldim) then
          sod(i,2,3) = sod(i,2,3) + (pph(2,i+1,2)-pph(1,i+1,2))*sod(i,4,3)
        endif
      enddo

      do k = 2, ndim
        i = indxsh(k)
        if (i <= ldim) then
          j = i-1
          sod(i,3,3) = sod(i,3,3) + (pph(2,j,1)-pph(1,j,1))*sod(j,4,3)
        endif
      enddo
      endif

C      call prmx('soph(3,:,:,:)',soph(3,:,:,:),4,4,ldim)
C      call prmx('sod',sod,ldim,ldim,12)
      end

