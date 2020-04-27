      subroutine pokepf(mode,kcplx,pfdim,ldp,lhmn,lhmx,iprmb,offp,ib1,ib2,nspc,nsp,Pfv,Pfout)
C- Write vector pfun in hamiltonian form for a range of sites.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  not used yet
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 Copy to Pfout
Ci         :1 Add to Pfout
Ci         :10s digit
Ci         :0 m is ordered m=l,..-l
Ci         :1 m is ordered m=-l,..l
Ci   kcplx :Ordering of complex arithmetic
Ci   pfdim :dimensions Pfv
Ci   ldp   :dimensions Pfout
Ci   lhmn,x:orbitals with index outside (lhmn+1,lhmx) are omitted from the copy
Ci   lhmx  :orbitals with index outside (lhmn+1,lhmx) are omitted from the copy
Ci   iprmb :If iprmb(1) = 0 then table is not used.  Orbitals are ordered 1,2,3,...
Ci         :permutation indices ordering orbitals in downfolding order
Ci   offp  :offset to first upper left corner of Pfout to be filled in
Ci   ib1,2 :Return Pfout for orbitals connected to range of sites (ib1:ib2)
Ci         :Each site in Pfv has indices for mxorb orbitals (mxorb is a global variable)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cio Inputs/Outputs
Cio  Pfv   :Global pfun, vector form
Cio  Pfout :Global pfun, matrix form
Cu Updates
Cu   05 Jun 18 (MvS) Redesign for new lmorder
Cu   25 Mar 16 Some redesign
Cu   07 Mar 16 (MvS) adapted from pokeg0
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer,intent(in) :: mode,kcplx,pfdim,ldp,lhmn,lhmx,iprmb(pfdim),offp,ib1,ib2,nspc,nsp
      complex(8),intent(in) ::  Pfv(pfdim,nspc*nsp)
      complex(8),intent(out) :: Pfout(ldp,nspc,ldp,nsp)
C ... Local parameters
      logical ldnfld
      integer morder,i,iv,ilm,im,m1,m2,l,mxorb,lmr,ib,lmax
      procedure(integer) :: nglob,iprint,ll

C ... Setup
      mxorb = nglob('mxorb'); lmax = ll(mxorb)
      if (mod(mode,10) == 0) call dpzero(Pfout,2*size(Pfout))
      morder = mod(mode/10,10)
      im = -1 ; if (morder==1 .or. morder==2) im = 1

      ldnfld = iprmb(1) /= 0
      if (kcplx /= 1) call ztoyy(Pfout,ldp*nsp,ldp*nsp,ldp*nsp,ldp*nsp,kcplx,1)

C ... Copy P associated with a range of sites
      i = offp
      do  ib = ib1, ib2
        lmr = mxorb*(ib-1) ! Offset to Pfv
        do  l = 0, lmax
        m1 = l**2 + 1 ; m2 = (l+1)**2
        do  ilm = m1, m2
          i = i+1; if (i > ldp) call rx('pokepf: dimensioning error')

          iv = lmr+ilm; if (ldnfld) iv = iprmb(lmr+ilm)
          if (iv <= lhmn .or. iv > lhmx) cycle ! Outside range of this block

C         Spin diagonal part
          Pfout(i,1,i,1) = Pfout(i,1,i,1) + Pfv(iv,1)
C         The next line merely duplicates the preceding when nsp=1
          Pfout(i,nspc,i,nsp) = Pfout(i,nspc,i,nsp) + Pfv(iv,nspc*nsp)

C         Spin off-diagonal part
          if (nspc == 1 .or. ilm+im < m1 .or. ilm+im > m2) cycle
          Pfout(i+im,2,i,1) = Pfout(i+im,2,i,1) + Pfv(iv,2)
          Pfout(i,1,i+im,2) = Pfout(i,1,i+im,2) + Pfv(iv+im,3)

        enddo
        enddo
      enddo

      if (iprint() >= 100/10) then
C        call zprm('pfv',2,pfv,pfdim,pfdim,4)
C        call zprm('pokepf: Pfout',2,Pfout,pfdim*nspc,lhmx*nspc,lhmx*nspc)
      endif

      if (kcplx /= 1) call ztoyy(Pfout,ldp*nsp,ldp*nsp,ldp*nsp,ldp*nsp,1,kcplx)
      end
