      subroutine getq(nsp,nl,lmx,nc,z,pnu,qnu,ics,s_spec,qc,qt,dq)
C- Gets the charge and related parameters for a set of atoms
C ----------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci         :Not used unless ics(1)>0
Ci     Elts read:  coreh coreq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  gtpcor
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :(global maximum l) + 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   nc    :number of classes
Ci   z     :nuclear charge
Ci   pnu   :pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci         :Integer part used to get princ. quant. number for counting
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci         :Set to zero to suppress information for core hole
Co Outputs
Co   qc:    core electronic charge
Co   qt:    total charge (nuclear + electronic) within sphere
Co   dq:    difference between spin up and spin down charge
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   04 Jan 06 Redesign to enable core hole
Cu   17 Feb 03 Changed dq to correct convention (q+ - q-)
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nsp,nc,lmx(nc),ics(nc)
      double precision pnu(nl,nsp,nc),qnu(3,nl,nsp,nc),
     .                 qc(nc),qt(nc),dq(nc),z(nc)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ic,kcor,lcor,is
      integer, parameter :: n0=10
      double precision qcor(2),xx(1),pz(n0)

      call dpzero(pz,n0)
      do  ic = 1, nc
        if (ics(1) == 0) then
          kcor = 0
          lcor = 0
        else
          is = ics(ic)
          call gtpcor(s_spec,is,kcor,lcor,qcor)
        endif
        call getqvc(nsp,nl,lmx(ic),z(ic),pnu(1,1,ic),qnu(1,1,1,ic),pz,
     .    0,0,kcor,lcor,qcor,qc(ic),qt(ic),dq(ic),xx,xx)
      enddo

      end
      subroutine getqvc(nsp,nl,lmx,z,pnu,qnu,pz,ncmx,nvmx,
     .                  kcor,lcor,qcor,qc,qt,dq,ec,ev)
C- Gets the charge and related parameters within an atom
C ----------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :dimensions pnu,qnu
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   z     :nuclear charge
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   pz    :boundary conditions for semicore levels, if they exist
Ci   ncmx: :0 do not initialize ec (see remarks)
Ci         :else, maximum allowed number of core orbitals
Ci   nvmx: :0 do not initialize ev (see remarks)
Ci         :else, maximum allowed number of valence orbitals
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Co Outputs
Co   qc:    core electronic charge
Co   qt:    total charge (nuclear + electronic) within sphere
Co   dq:    difference between spin up and spin down charge
Co   ec:    starting guess for core eigenvalues (see remarks)
Co   ev:    starting guess for valence eigenvalues (see remarks)
Cl Local variables
Cl   npnu:  number of valence princ. quant. numbers for a given l
Cl          npnu(l) = 1, or 2 pz(l) is a semicore state
Cr Remarks
Cr   Generates charges and (optionally) initializes ec and ev
Cr   ec is initialized if ncmx > 0 and is otherwise unused
Cr   ev is initialized if nvmx > 0 and is otherwise unused
Cu Updates
Cu   28 May 14 Extended to include possible semicore states: new pz
Cu   26 Jan 06 bug fix for core hole, case lmx < lcor
Cr   (ATP) adapted from getq for one sphere, enables core hole
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nsp,lmx,ncmx,nvmx,kcor,lcor
      double precision pnu(nl,nsp),qnu(3,nl,nsp),qcor(2),
     .                 qc,qt,dq,z,ec(*),ev(*),pz(nl)
C ... Local parameters
      integer l,k,konf,konfig,konfg(0:10),isp,lmaxc,ncore,nval
      integer iprint,npnu(0:10)
      double precision ecore0,eval0,evalsc0,deg,pl
      parameter (ecore0=-5.d0, eval0=-.5d0, evalsc0=-2d0)

      call config(pnu,lmx,z,konfg,lmaxc)
      if (kcor > 0) then
        lmaxc = max(lmaxc,lcor)
        konfg(lcor) = max(konfg(lcor),kcor+1)
      endif
      qc = 0
      qt = 0
      dq = 0
      ncore = 0
      nval  = 0
C --- Core charge and trial values for core levels
      do  l = 0, lmaxc
        do  isp = 1, nsp
          deg = 2*(2*l+1)
          konfig = konfg(l)
          if (l < nl) then
            pl = pnu(l+1,isp); npnu(l) = 1
            if (pz(l+1) > 0 .and. int(pl-1) == mod(int(pz(l+1)),10)) then
              pl = mod(pz(l+1),10d0); npnu(l) = 2 ! increase number of valence p.q.n
            endif
            if (l <= lmx) konfig = pl
          endif
          if (konfig <= l) then
            if (iprint() == 0) call setpr(10)
            call fexit2(-1,1,' Exit -1 GETQVC: '//
     .      'bad pnu(l=%i) (%,1d)',l,pl)
          endif
          do  konf = l+1, konfig-1
            if (ncmx /= 0) then
              ncore = ncore+1
              ec(ncore) = ecore0
            endif
            if (konf == kcor .and. l == lcor) then
              deg = deg + qcor(1)
              if (nsp > 1) dq = dq + qcor(2)/nsp
            endif
            qc = qc + deg/nsp
            if (konf == kcor .and. l == lcor) deg=deg-qcor(1)
          enddo                 ! Loop over core pqn
        enddo                   ! Loop over spins
      enddo                     ! Loop over l

C --- Valence charge and trial values for valence energies
      do  k = 1, lmx+1
        if (nsp > 1) dq = dq + qnu(1,k,1) - qnu(1,k,2)
        do  isp = 1, nsp
          if (nvmx /= 0) then
            if (npnu(k-1) == 2) then  ! Semicore state
              nval = nval+1
              ev(nval) = evalsc0
            endif
            nval = nval+1
            ev(nval) = eval0
          endif
          if (npnu(k-1) == 2) qt = qt + (4*k-2)/nsp ! Semicore state
          qt = qt + qnu(1,k,isp)
        enddo
      enddo

      if (ncore > ncmx .or. nval > nvmx) call rx('GETQVC: too many orbitals')

      qt = qc + qt - z
      end
