      subroutine lumo(lmol,zval,nste,ldim,nsp,nspc,nbmax,eband,wtkb)
C- Try and make a HOMO-LUMO excitation by manipulating band weights
C ----------------------------------------------------------------------
Ci Inputs:
Ci   zval, nsp, wtkb
Ci   use HOMO-nste(1) and LUMO+nste(2)
Co Outputs:
Co   wtkb
Cr Remarks
Cr   Try to make a HOMO-LUMO excitation by exchanging the weights of the
Cr   HOMO and LUMO.
Cb Bugs
Cb   If the HOMO or LUMO are degenerate then this routine wil fail
Cb   with an error exit. In that case the weights should have to be
Cb   remade as in evlwsm in mol; this is not yet implemented in tbe.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nsp,nspc,ldim,nbmax,nste(2)
      double precision zval,wtkb(ldim,2),eband(nbmax,2)
      logical lmol
C Local Variables
!       integer obot,otop,obmap,owk
      integer, allocatable :: bmap(:)
      real(8), allocatable :: bot(:), top(:), wk(:)
      integer i,isp,ibh,isph,ibl,ispl,nbpw,nkp,nevxx,nbmxx
      integer i1mach,iprint
      double precision elo,ehi,ehomo,elumo,wt,qv,dummy
      logical lfill,efrng2,fndh,fndl
      logical, parameter :: T = .true., F = .false.

C --- Sanity checks ---
      call rxx(.not.lmol, 'Use LUMO for molecules only')
      call rxx(nsp /= 2,  'LUMO needs spin pol')
      call rxx(nspc == 2, 'LUMO needs npsc=1')

      if (iprint() > 40) then
        print *,'  Bands and weights ..'
        print *,'  Eigenvalue            weight'
        write (*,10) ((eband(i,isp),wtkb(i,isp),i=1,ldim),isp=1,nsp)
      endif

      nkp=1
      nbpw = int(log(real(i1mach(9),8)+1d0)/log(2d0))
      allocate(bmap(ldim*nsp*nkp/nbpw+1)); bmap = 0
      allocate(wk(ldim*nsp))
C --- gather bands ---
      call ebcpl(0,nbmax,ldim,nsp,nspc,nkp,nbpw,bmap,wk,eband)
      allocate(bot(ldim*nsp))
      allocate(top(ldim*nsp))
      nevxx = ldim*nsp
      nbmxx = nbmax*nsp
C --- find HOMO-n
      qv = (zval - 2d0*nste(1)) * 2d0
      lfill = efrng2(1,nkp,nbmxx,nevxx,qv,eband,bot,top,
     .               elo,ehi,ehomo,dummy)
C --- find LUMO+n
      qv = (zval + 2d0*nste(2)) * 2d0
      lfill = efrng2(1,nkp,nbmxx,nevxx,qv,eband,bot,top,
     .               elo,ehi,dummy,elumo)
      if (iprint() >= 40) then
        call awrit4(' LUMO efrng2 found eHOMO(-%i)=%d, eLUMO(+%i)%d',
     .              ' ',120,i1mach(2),nste(1),ehomo,nste(2),elumo)
      endif
C --- scatter bands ---
      call ebcpl(1,nbmax,ldim,nsp,nspc,nkp,nbpw,bmap,wk,eband)
!       call rlse(obmap)
      deallocate(top,bot,wk,bmap)

      isph = 0
      fndh = F
      do  i = 1, ldim
        do  isp = 1, nsp
          if (dabs(eband(i,isp)-ehomo) < 1d-12) then
            if (fndh) then
              if (isp == isph) then
                call rx('LUMO: HOMO is degenerate')
              endif
            endif
            ibh = i
            isph = isp
            fndh = T
          endif
        enddo
      enddo
      fndl = F
      do  i = 1, ldim
C ... flip the spin: Hund's rule
        do  isp = nsp, 1, -1
          if (dabs(eband(i,isp)-elumo) < 1d-12) then
            if (fndl) then
              if (isp == ispl) then
                call rx('LUMO: LUMO is degenerate')
              endif
            endif
            ibl = i
            ispl = isp
            fndl = T
          endif
        enddo
      enddo

      if (iprint() >= 30) then
        call awrit4(' LUMO found eHOMO=%d (weight=%d), eLUMO=%d'//
     .              ' (weight=%d)',' ',120,i1mach(2),
     .              eband(ibh,isph),wtkb(ibh,isph),
     .              eband(ibl,ispl),wtkb(ibl,ispl))
      endif

C --- exchange weights ---
      wt = wtkb(ibh,isph)
      wtkb(ibh,isph) = wtkb(ibl,ispl)
      wtkb(ibl,ispl) = wt

      if (iprint() > 31) then
        print *,'  Modified bands and weights ..'
        print *,'  Eigenvalue            weight'
        write (*,10) ((eband(i,isp),wtkb(i,isp),i=1,ldim),isp=1,nsp)
      endif
   10 format (f10.6,4x,f10.6)

      end
