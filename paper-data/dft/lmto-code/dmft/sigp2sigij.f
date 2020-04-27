      subroutine sigp2sigij(s_dmft,mode,noms,nsp,iomg,cix,ldcix,sigmapr,sigmapc,sigmaij)
C- Converts sigma stored in compact form to/from matrix form
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_dmft
Ci     Elts read:  nzsigi iasig
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:icix l
Cio    Passed to:  *
Ci   mode  : 1s digit
Ci         : 0 Return uncompressed sigma for one block, m block only in sigmaij
Ci         :   row and column indices for m range from 1 ... 2l+1
Ci         :   sigmapr is accessed; sigmapc is not touched
Ci         : 1 Return uncompressed sigma for one block, lm ordering, in sigmaij
Ci         :   row and column indices for m range from l**2+1 ... (l+1)**2
Ci         :   sigmapc is accessed; sigmapr is not touched
Ci         : 2 Add 2 to reverse the sense of copy : sigmij copied to sigmapr or sigmapc
Ci         : 10s digit
Ci         : 0 Compressed sigma read/written to sigmapr (real format)
Ci         : 1 Compressed sigma read/written to sigmapc (complex format)
Ci         :   (not implemented yet)
Ci         : 100s digit
Ci         : 1 Skip cix if not the first inequivalent one (NO LONGER USED)
Ci         : 1000s digit (not implemented yet)
Ci         : 1 Rotate sigma by group operation
Ci   noms  : Leading dimension of sigmapr, sigmapc
Ci   iomg  : Frequency index to sigmapr,sigmapc
Ci   cix   : which cix block to copy
Ci   ldcix : dimensions sigmaij
Cio Inputs/Outputs
Cio sigmapr: Compressed sigma, real form
Cio sigmapc: Compressed sigma, complex form
Cio sigmaij: Uncompressed sigma (complex form)
Cs Command-line switches
Cl Local variables
Cr Remarks
Cr   Compressed form addressable for frequency iomg.
Cr   Uncompressed form in one frequency only: it should correspond to frequency iomg
Cu Updates
Cu   18 May 18 (MvS) Enable equivalent cix blocks with spin flip
Cu   20 Feb 16  First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_dmft) :: s_dmft
C ... Passed parameters
      integer mode,iomg,cix,ldcix,noms,nsp
      real(8) :: sigmapr(noms,*)
      complex(8) :: sigmapc(noms,*),sigmaij(ldcix,ldcix,nsp)
C ... Local parameters
      logical ltmp
      integer ichan,isigind,mode0,mode1,m1,m2,isp,ilm1,ilm2,icix,offl,cixl
C     double precision fac
      procedure(logical) :: isanrg

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      ltmp = isanrg(mode0,0,3,'sigi2sigij','mode0',.true.)
      ltmp = isanrg(mode1,0,1,'sigi2sigij','mode1',.true.)
      if (mode0 < 2) call dpzero(sigmaij,2*size(sigmaij))

!     do  cixl = cix, cix
      cixl = cix
      icix = iabs(s_dmft%icix(cixl))
C      fac = 1
C      if (s_dmft%icix(cixl) < 0) then
C        fac = -1
C        if (mod(mode/100,10) == 1) goto 99
C      endif
      offl = 0 ; if (mode0 == 1 .or. mode0 == 3) offl = (s_dmft%l(icix))**2
      do  ichan = s_dmft%nzsigi(cixl-1)+1, s_dmft%nzsigi(cixl)
        m1  = s_dmft%iasig(1,ichan)
        m2  = s_dmft%iasig(2,ichan)
        isp = s_dmft%iasig(3,ichan)
        cix  = s_dmft%iasig(5,ichan)
        isigind = s_dmft%iasig(4,ichan) ! iasig(4) knows about afm indexing
C       if (s_dmft%icix(cix) < 0 .and. nsp == 2) isp = nsp+1-isp
        ilm1 = m1 + offl
        ilm2 = m2 + offl

C       Return in sigmaij uncompressed sigma for this block from sigmapr
        if (mode0 <= 1 .and. mode1 == 0) then
          sigmaij(ilm1,ilm2,isp) = dcmplx(sigmapr(iomg,2*isigind-1),sigmapr(iomg,2*isigind))
C       Return in sigmaij uncompressed sigma for this block from sigmapc
        elseif (mode0 <= 1 .and. mode1 == 1) then
          sigmaij(ilm1,ilm2,isp) = sigmapc(iomg,isigind)
C       Compress sigma for this block to sigmapr
        elseif (mode1 == 0) then
          sigmapr(iomg,2*isigind-1) = dble(sigmaij(ilm1,ilm2,isp))
          sigmapr(iomg,2*isigind)   = dimag(sigmaij(ilm1,ilm2,isp))
C       Compress sigma for this block to sigmapc
        else
          sigmapc(iomg,isigind) = sigmaij(ilm1,ilm2,isp)
        endif
      enddo
!      enddo ! cix

      end
