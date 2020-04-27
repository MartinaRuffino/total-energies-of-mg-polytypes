C#define NSPH
      subroutine asamad(s_ctrl,s_pot,s_lat,s_spec,lq,p,q,vrl,ves,emad,
     .  trumad,vmtz,etrms)
C- Calculate ASA Madelung potential, Madelung energy, vmtz
C ----------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nl nspec nspin lpgf lmet nclass zbak nclasp
Ci                 npadl npadr ldlm nccomp ipcp nrc ips ipc rmax
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ics ncomp lasa nrc ipcp dclabl rmax idcc ipc ips
Cio    Passed to:  asavqm asaqmp
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vconst qpp vmtz0
Co     Stored:     vconst
Co     Allocated:  *
Cio    Elts passed:vrmax gibbs mad qc qpp pmpol
Cio    Passed to:  asavqm asaqmp
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat as tol vol nkdmx nkqmx platl platr qlat
Ci                 awald nkd nkq nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:vol pos cg jcg indxcg qlv dlv symgr ag
Cio    Passed to:  asavqm asaqmp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa z coreh coreq name lmxl lmxf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  getq gtpcor dlmq asavqm
Ci Inputs
Ci   lq   :  0 input p,q are pnu,qnu; make qt from pnu,qnu
Ci        :    NB: Dirac case: make qt from s_pot%qnur
Ci           1 input q is qt; p is not used
Ci           2 Take ves as input; evaluate vmtz only
Ci          10s digit nonzero substitutes nonspherical ves for 0,1
Ci             iff qpp available
Ci             NB: qpp should be passed as argument, and lq more
Ci                 sensibly arranged.
Ci          100s digit applies to the layer branch (lpgf>0)
Ci               0 No special treatment: estat potential of all layers
Ci                 computed from sphere charges
Ci               1 Estat potential of sites in embedding region
Ci                 computed from Madelung, but the potential of
Ci                 sites in L- and R- regions are computed as though
Ci                 in bulk with periodic boundary conditions.
Ci               2 Display the deviation in the potential in L- and R-
Ci                 PL as computed by the two approaches above.
Ci               3 Add to vconst(1) a best average of this deviation.
Ci   p     :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci          Used only to compute nuclear contribution to sphere charge
Ci   q     :energy-weighted moments of the sphere charges.
Ci          Used to compute electronic contribution to sphere charge
Cl   vrl   :Difference in estat potential in the left- and right- electrodes (PGF)
Ci         :vrl = bias across device: + [vconst(3)-vconst(2)]
Ci         :First term = ef(R)-efermi, where ef(R) is fermi energy of right lead
Ci         :Second term is needed because of how Poisson's eqn is solved (see lmasa-gf.f)
Ci  sarray: requires arrays oz
Co Outputs (see Remarks for qualifications)
Co   ves:  Hartree potential on sphere surfaces
Co   emad: Madelung energy using convention as discussed in remarks
Co   trumad:True interatomic Madelung energy; see madpot
Co   vmtz: muffin-tin zero (vmtz(2) = vmtz(spin-up)-vmtz(spin-dn))
Co   etrms:Madelung energy for each class is stored
Cl Local variables
Cl   These variables only apply to the layer case
Cl   vavo :electrostatic potential of the layers computed from
Cl        :multilayer
Cl        :(1) = avg ves in L
Cl        :(2) = avg ves in R
Cl   vavn :electrostatic potential computed from charges for L- and R-
Cl        :endpoints individually , under the condition each endpoint
Cl        :is subject to its own periodic boundary conditions as
Cl        :L- and R- layers were separate (periodic) materials.
Cl        :(1) = avg ves in L
Cl        :(2) = avg ves in R
Cr Remarks
Cr   asamad generates the ASA electrostatic potential by Ewald
Cr   summation (see madpot.f, below).  There is an additional
Cr   potential that may be added; this has no fundamental signifance
Cr   but has a practical significance in various contexts.  First,
Cr   double counting of the potential shifts must be eliminated in
Cr   the total energy.  Second, in some cases the electrostatic
Cr   potential is not sought, for models or to facilitate convergence
Cr   to self-consistency.  Last, in cases when the energy mesh or
Cr   Fermi level is explicit, as in Green's function calculations,
Cr   the mesh is affected by potential shifts.
Cr
Cr  *Case 1s digit lq eq 2:
Cr   The input ves is used.
Cr
Cr  *The USUAL band structure case (lq ne 2):
Cr   Electrostatic potential calculated directly by Ewald summation.
Cr
Cr  *The LAYER Green's function case (lq eq 0,1):
Cr   Case separate bulk calculations (lpgf eq 2):
Cr   The electrostatic potentials are calculated for the left-
Cr   and right- bulk PL independently.
Cr
Cr   Case interface calculation (lpgf eq 1):
Cr   The electrostatic potentials initially calculated as in the usual
Cr   band-structure case.  Were the left and right bulk-like PL
Cr   really bulk-like, this would be sufficient, apart from a possibly
Cr   constant shift for insulating L or R.  (The shift is zero for
Cr   a metallic end PL).  But there may arise unwanted deviations owing
Cr   to the nonlocal contribution from other layers.  The potentials
Cr   in these layers are obtained from the electrostatic potentials of
Cr   the bulk, plus a constant shift, calculated as:
Cr       the average electrostatic potential of the interface -
Cr       the average electrostatic potential of the bulk PL
Cr   The caller can force either L or R, or both shifts be zero.
Cr   Establish this condition by ctrl->lmet, bit 2 (L) and bit 3 (R).
Cu Updates
Cu   17 Oct 17 asamad calls madmat for leads in serial mode
Cu   07 Aug 16 Charges taken from qnur in relativistic case
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   29 Nov 10 (P. Larson) extensions for disordered local moments
Cu   21 Jul 07 (pgf) vne->vrl (for inequivalent left- and right- end layeers)
Cu   10 Feb 04 (S.Faleev) sbz struc added to argument list
Cu              to handle non-equilibrium mode
Cu   22 Mar 03 Revised electrostatics for layer case;
Cu             estat potential is be shifted by pot->vconst
Cu   19 Feb 02 Some changes adapt to updates in layer code
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lq
      double precision ves(*),emad,trumad,vmtz(2),p(*),q(*),vrl,etrms(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: lmx(:)
      real(8), allocatable :: dq(:),qc(:),qt(:)
      real(8), allocatable :: z(:),pos2(:)
      real(8), allocatable :: vval(:),vh(:),vesl(:)
      real(8), allocatable :: qmp(:),mad2(:),qtotr(:,:)
      real(8), allocatable :: dlv(:),qlv(:)
C ... Local parameters
      logical lves,metL,metR,bittst
      integer nbas,nl,nsp,nspec,nclasp,nclspp,npadl,npadr,nkdmx,nkqmx,
     .  lpgf,iprint,nbaspp,lmxst,nkd,nkq,ic,ic1,ic2,jc1,jc2,ib,
     .  nclass,lmet,lnsph,isw,nlmf,stdo
      double precision rhrmx(2),platl(3,3),platr(3,3),alat,awald0,
     .  tol,vol,plat(9),qlat(9),awald,xx,vavo(5),vavn(5),dum(3),
     .  vglob,vconst(3),rmsv(3)
C ... Parameters for DLM
      integer ldlm,nangl,nclspd,nclsppd
      integer,parameter:: NULLI = -99999, maxth = 100
C ... Parameters for non-equilibrium mode (S.F.)
      double precision plat0(3,3)
      procedure(integer) :: nglob

C      call pshpr(30)
C      print *, '!! asamad',lq

C ... Setup
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      lpgf = s_ctrl%lpgf(1)
      lmet = s_ctrl%lmet
      nclass = s_ctrl%nclass
      rhrmx = s_ctrl%zbak
      call dscal(2,1/s_lat%vol,rhrmx,1)
      nclasp = s_ctrl%nclasp
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      vconst = s_pot%vconst
      alat = s_lat%alat
      plat0 = s_lat%plat

      ldlm = s_ctrl%ldlm
      nangl = s_ctrl%nccomp
      stdo = nglob('stdo')

      nclspp = 2*nclasp-nclass
      nbaspp = 2*(npadl+npadr) + nbas
      nclspd = nclasp + nangl
      nclsppd = nclspp + nangl
      allocate(qt(nclsppd))
      if (mod(lq,100) == 0) then
        allocate(dq(nclsppd))
        allocate(qc(nclsppd))
        allocate(lmx(nclsppd))
        allocate(z(nclsppd))
        call spec2class(s_spec,nclsppd,s_ctrl%ics,'lmxa',1,lmx,xx)
        call spec2class(s_spec,nclsppd,s_ctrl%ics,'z',1,xx,z)
        call getq(nsp,nl,lmx,nclsppd,z,p,q,s_ctrl%ics,s_spec,qc,qt,dq)
        if (mod(s_ctrl%lrel,10) == 2 .and. s_pot%qnur(1,1) /= NULLI) then
          allocate(qtotr(2,nclsppd))
          call qtotrel(nl,nclsppd,s_pot%qnur,qtotr)
          do  ic = 1, nclsppd
            dq(ic) = qtotr(1,ic) - (z(ic)-qc(ic))
          enddo
          if (iprint() >= 30) then
            call info0(1,1,0,'')
            call arrprt(' Class dq(SR)  Dirac','%,4i%,4;9D%,4;8D',
     .          'Idd',nclsppd,0,3,0,'  | ',xx,qt,dq,xx,xx,xx,xx,xx)
          endif
          qt(:) = dq(:)
        endif
        deallocate(dq,qc,lmx,z)
      elseif (mod(lq,100) == 1) then
        call dcopy(nclsppd,q,1,qt,1)
      else
        call dpzero(qt,nclsppd)
      endif
      if (ldlm > 0) call dlmq(nclspp,s_spec,s_ctrl%ics,qt,
     .  s_pot%vrmax,s_ctrl%ncomp,s_pot%gibbs)
      lves = mod(lq,100) == 2

C ... Copy charges and potential to doubly padded layer
C      if (nclasp /= nclass) then
C        call dpscop(qt,qt,nclasp-nclass,nclass+1,nclasp+1,1d0)
C        call dpscop(ves,ves,nclasp-nclass,nclass+1,nclasp+1,1d0)
C        call dpscop(s_pot%vrmax,s_pot%vrmax,2*(nclasp-nclass),
C     .    2*nclass+1,2*nclasp+1,1d0)
CC        call prmx('vrmax',s_pot%vrmax,2,2,nclspp)
CC        call prmx('qt',qt,nclspp,nclspp,1)
C      endif

C      call setpr(50)
C      call dvset(qt,nclass+1,nclass+1,.2d0)
*      call dvset(qt,nclass+2,nclass+2,-.2d0)

C#ifdef NSPH
C --- Make potential from multipole moments ---
C     lnsph = isw(lgors('ctrl lasa,32',sctrl))
      lnsph = isw(IAND(s_ctrl%lasa,32) /= 0)
      if (lnsph /= 0 .and. .not. lves) then
      if (s_pot%qpp(1) >= 0d0) then
        nlmf = (2*nl-1)**2
        allocate(vval(nlmf*nbas)); call dpzero(vval,nlmf*nbas)
        allocate(qmp(nlmf*nbas))
        allocate(vh(nbas))
        call asavqm(0,s_ctrl,s_pot,s_lat,s_spec,nlmf,vh,qmp,vval)
        call dpzero(ves,nclasp)
        do ib = 1,nbas+(npadl+npadr)
          ic = s_ctrl%ipcp(ib)
          xx = 1/dble(s_ctrl%nrc(ic))
          ves(ic) = ves(ic) + vh(ib)/s_ctrl%nrc(ic)
        enddo
        lves = .true.
        deallocate(vval,qmp,vh)
      endif
      endif
C#endif

C --- Electrostatic potential for charges as given ---
      if (lpgf /= 2) then
        if (lpgf /= 0) call info(20,1,0,
     .  ' Electrostatics for embedded L-C-R system:%?#n# vrl = %d##',isw(abs(vrl) >= 1d-6),vrl)
        vmtz(1) = s_pot%vmtz0
        call madpot(nbaspp,1,nclspp,s_ctrl%nrc,s_ctrl%ipcp,
     .    s_ctrl%dclabl,qt,vconst,rhrmx,s_ctrl%rmax,s_pot%mad,
     .    s_lat%pos,plat0,vrl,lves,s_pot%vrmax,ves,emad,trumad,vmtz,
     .    nangl,s_ctrl%idcc,s_pot%gibbs,etrms)
        if (abs(vrl) >= 1d-6) then
          allocate(vh(nbaspp))
          call dplmom(nbaspp,s_lat%pos,alat,plat0,s_ctrl%ipcp,qt,0d0,vh,dum(1),dum(2))
          deallocate(vh)
        endif
C       call snot('ves[qt]',nbaspp,ves,s_ctrl%ipcp)
C   ... Nothing further if potentials are input
        if (lves) goto 99
      endif

C --- Layer GF: es pot for end layers and pot shift for ch. neutrality
      if (lpgf /= 0 .and. mod(lq/100,100) > 0) then
        awald0 = s_lat%as
        tol = s_lat%tol
        vol = s_lat%vol
        nkdmx = s_lat%nkdmx
        nkqmx = s_lat%nkqmx
        platl = s_lat%platl
        platr = s_lat%platr
C       nbasp = nbas + npadl + npadr
C   ... Save original ves(Mad) for printout
        allocate(vesl(nclspd))
        call dpcopy(ves,vesl,1,nclspd,1d0)

C   ... Average potentials for interface L,R, to calculate shift
        ic1 = nclasp
        ic2 = 1
        do  ib = nbas, nbas+npadl-1
          ic1 = min(ic1,s_ctrl%ipcp(ib+1))
          ic2 = max(ic2,s_ctrl%ipcp(ib+1))
        enddo
        jc1 = nclasp
        jc2 = 1
        do  ib = nbas+npadl, nbas+npadl+npadr-1
          jc1 = min(jc1,s_ctrl%ipcp(ib+1))
          jc2 = max(jc2,s_ctrl%ipcp(ib+1))
        enddo
        call pvmadp(ic1,ic2,s_ctrl%nrc,0,xx,s_ctrl%rmax,ves,vavo(1),
     .    rmsv)
        call pvmadp(jc1,jc2,s_ctrl%nrc,0,xx,s_ctrl%rmax,ves,vavo(2),
     .    rmsv(2))

C        print *, ic1,ic2
C        print *, jc1,jc2
C       stop

C        call dvset(qt,2+ic1,2+ic1,.1d0)
C        call dvset(qt,2+ic1+1,2+ic1+1,-.1d0)

C   ... Electrostatic potential for bulk L
        allocate(dlv(3*nkdmx))
        allocate(qlv(3*nkqmx))
        allocate(pos2(max(npadl,npadr)*3))
        allocate(mad2(max(npadl,npadr)**2))
        lmxst = 6
        call pshpr(1)
        call lattc(awald0,tol,0d0,alat,alat,platl,0d0,0d0,1d0,1d0,plat,
     .    qlat,lmxst,vol,awald,dlv,nkd,qlv,nkq,nkdmx,nkqmx)
        call dpscop(s_lat%pos,pos2,3*npadl,3*nbas+1,1,1d0)
c       call prmx('pos',pos2,3,3,npadl)
        call madmat(npadl,pos2,awald,alat,-abs(vol),dlv,nkd,qlv,nkq,mad2)
        call poppr
c       call prmx('mad',mad2,npadl,npadl,npadl)
        vavn(1) = 0
        call info(20,1,0,
     .    ' ES potential for infinitely repeating L princ layer:',0,0)
        call madpot(npadl,ic1,ic2,s_ctrl%nrc,s_ctrl%ipcp(nbas+1),
     .    s_ctrl%dclabl,qt,vconst(2),rhrmx,s_ctrl%rmax,mad2,xx,
     .    xx,0d0,lves,s_pot%vrmax,ves,emad,xx,vavn(1),0,xx,xx,etrms)
C       Keep vavn(1) for sanity check (see if unchanged after bulk R)
        call pvmadp(ic1,ic2,s_ctrl%nrc,0,xx,s_ctrl%rmax,ves,dum(1),rmsv)

C   ... Electrostatic potential for bulk R
        call pshpr(1)
        call lattc(awald0,tol,0d0,alat,alat,platr,0d0,0d0,1d0,1d0,plat,
     .    qlat,lmxst,vol,awald,dlv,nkd,qlv,nkq,nkdmx,nkqmx)
        call dpscop(s_lat%pos,pos2,3*npadr,3*(nbas+npadl)+1,1,1d0)
        call madmat(npadr,pos2,awald,alat,-abs(vol),dlv,nkd,qlv,nkq,mad2)
        call poppr
        vavn(2) = 0
        call info(20,1,0,
     .    ' ES potential for infinitely repeating R princ layer:',0,0)
        call madpot(npadr,jc1,jc2,s_ctrl%nrc,
     .    s_ctrl%ipcp(nbas+npadl+1),s_ctrl%dclabl,qt,vconst(3),
     .    rhrmx,s_ctrl%rmax,mad2,xx,xx,0d0,lves,s_pot%vrmax,ves,emad,
     .    trumad,vavn(2),0,xx,xx,etrms)
        deallocate(dlv,qlv,pos2,mad2)

        call pvmadp(ic1,ic2,s_ctrl%nrc,0,xx,s_ctrl%rmax,ves,vavn(1),rmsv)
        call pvmadp(jc1,jc2,s_ctrl%nrc,0,xx,s_ctrl%rmax,ves,vavn(2),rmsv(2))

C       Sanity check: did <ves(L)> change after calc ves(R)?
        if (abs(dum(1)-vavn(1)) > 1d-6 .and. iprint() >= 10) then
          call info2(10,1,0,'  asamad (warning): adjusting vconst(L):'//
     .      ' average ves(L) altered by %,6;6d%N%8f after ves(R) '//
     .      'calculated (possible misuse of equivalent classes)',
     .      vavn(1)-dum(1),0)
          vconst(2) = vconst(2)+vavn(1)-dum(1)
          s_pot%vconst = vconst
        endif

C   ... RMS deviation in the two ways to compute ves[padded]
        if (lpgf /= 2) then
          call dpadd(ves,vesl,1,nclasp,-1d0)
          call pvmadp(ic1,ic2,s_ctrl%nrc,0,xx,s_ctrl%rmax,ves,dum(1),rmsv)
          call pvmadp(jc1,jc2,s_ctrl%nrc,0,xx,s_ctrl%rmax,ves,dum(2),rmsv(2))
          call dpadd(ves,vesl,1,nclasp, 1d0)

          if (iprint() >= 30) then
            metL = bittst(lmet,8)
            metR = bittst(lmet,4)
            write(stdo,334)
  334       format(/' Deviations in end potentials:'/' region met',
     .        '  <ves>Bulk','   <ves>layer     Diff',6x,'RMS diff')
            write(stdo,335) 'L', metL,
     .        vavn(1), vavo(1), vavo(1)-vavn(1), dsqrt(rmsv(1))
            write(stdo,335) 'R', metR,
     .        vavn(2), vavo(2), vavo(2)-vavn(2), dsqrt(rmsv(2))
  335       format(2x,a,L7,5f12.6)
          endif
          rmsv(3) = dsqrt((rmsv(1)**2*npadl+rmsv(2)**2*npadr)/(npadl+npadr))
          call info5(20,1,0,' RMS pot difference in '//
     .      'L PL = %,6;6d  in R PL = %,6;6d  total = %,6;6d',
     .      rmsv,rmsv(2),rmsv(3),0,0)
          vglob = vconst(1) + (vavn(1) - vavo(1) - vavo(2) + vavn(2))/2

          call info0(40,1,0,' vconst should be fixed by charge '//
     .      'neutrality requirements.'//
     .      '%N%8fIt can found iteratively in a GF pass using mode 1.')

          if (iprint() < 30)
     .      call info2(20,0,0,' vconst is now %,6;6d  vconst '//
     .      'minimizing RMS diff to end layers = %,6;6d',vconst,vglob)
          if (iprint() >= 30) then
          call info2(20,0,0,' vconst that minimizes potential '//
     .      'mismatch to end layers =%;10,6D',vglob,0)
          call info2(20,0,0,' vconst is now (estimate to'//
     .      ' satisfy charge neutrality)  =%;10,6D',vconst,0)
          call info2(20,0,0,' difference%45f=%;10,6D',vglob-vconst(1),0)
          endif

          if (iprint() >= 40) call query('vconst=',4,vconst)

          s_pot%vconst = vconst
        else
          call info2(20,1,0,' shifts added to bulk ES potential: '//
     .      ' vconst(L)=%,6;6d   vconst(R)=%,6;6d',vconst(2),vconst(3))
        endif


C   OLD
CC   --- Global potential shift, layer case ---
CC       We are permitted to shift all layers by a single constant
CC       An end layer that is metallic should not have any potential
CC       shifts to satisfy charge neutrality.  However, this requirement
CC       cannot be exactly satisfied if both end layers are metallic,
CC       because there is only one free parameter.  In that case we pick
CC       the average of the two layers.  If only one is metallic, it
CC       determines the global potential shift.  If neither is metallic,
CC       choose the shift from the central layer shift.
C        if (lpgf /= 2 .and. mod(lq/100,10) >= 2) then
C
CC         fac=0 if show shifts, but don't add to vshft; fac=1 to add
C          fac = 0
C          if (mod(lq/100,10) >= 3) fac = 1
C
CC     ... Compute global potential shift according to comments above
CC         and assign vavo(3..5) to shifts we will actually add
CC         to existing values.
CC         L layer: vglob would be vavn(1) - vavo(1)
CC         R layer: vglob would be vavn(2) - vavo(2)
CC         C layer: vglob would be dval(s_pot%vshft,6) ??
C          metL = bittst(lmet,8)
C          metR = bittst(lmet,4)
C
C          if (metL .and. metR) then
C            vglob   = (vavn(1) - vavo(1) - vavo(2) + vavn(2))/2
C            vavo(3) = 0*vshfl
C            vavo(4) = 0*vshfr
C            vavo(5) = vglob
C          elseif (metL) then
C            vglob   = (vavn(1) - vavo(1))
C            vavo(3) = 0*vshfl
C            vavo(4) = vavo(4) - (vavo(1) - vavn(1))
C            vavo(5) = vglob
C          elseif (metR) then
C            vglob   = (vavn(2) - vavo(2))
C            vavo(3) = vavo(3) + (vavn(2) - vavo(2))
C            vavo(4) = 0*vshfr
C            vavo(5) = vglob
C          else
C            vglob = vavo(5)
C            stop 'not ready'
C            vavo(3) = vavo(3) + vglob
C            vavo(4) = vavo(3) + vglob
C            vavo(5) = 0
C          endif
C
C          if (iprint() >= 20) then
C            write(stdo,334) metL,metR,vglob
C  334       format(/' Global potential shift to best satisfy charge',
C     .        ' neutrality:'/' met(L) = ',L1,'  met(R) = ',L1,
C     .        '  global vshft =',f11.6//
C     .        ' region <ves>Bulk  <ves>layer     Diff',
C     .        4x,'<vshft_in>',2x,'<vshft_out>')
C            write(stdo,335) 'L', vavn(1), vavo(1), vavo(1)-vavn(1),
C     .        vshfl, vshfl+vavo(3)*fac
C            write(stdo,335) 'R', vavn(2), vavo(2), vavo(2)-vavn(2),
C     .        vshfr, vshfr+vavo(4)*fac
C            write(stdo,336) 'C', vshfc,vshfc*(1-fac)+vavo(5)*fac
C  335       format(2x,a,1x,5f12.6)
C  336       format(2x,a,37x,2f12.6)
C          endif
C          vshfl = vshfl + vavo(3)*fac
C          vshfr = vshfr + vavo(4)*fac
CC         vshfc = vshfc + vavo(5)*fac
C          vshfc = vshfc*(1-fac)+vavo(5)*fac
C          call dvset(s_pot%vshft,2,2,vshfl)
C          call dvset(s_pot%vshft,4,4,vshfr)
C          call dvset(s_pot%vshft,6,6,vshfc)
C          if (fac /= 0) call dvset(s_pot%vshft,8,8,1d0)
C
C          if (fac /= 0) then
C            write(stdo,500)
C            write(lgunit(2),500)
C            do  20  ic = 1, jc2
C              vadd = vshfc
C              if (ic > nclass) vadd = vshfl
C              if (ic > ic2) vadd = vshfr
C              call r8tos8(dval(s_ctrl%dclabl,ic),clabl)
Cc              ves(ic) = ves(ic) + vadd
C              if (iprint() >= 30) then
C                do  21  i = 1, 2
C   21           write(lgunit(i),501) clabl,dval(qt,ic),
C     .              dval(vesl,ic),ves(ic),ves(ic)+vadd,ves(ic)+vadd+
C     .              dval(s_pot%vrmax,2*ic-1)
C  501           format(1x,a,5f11.6)
C  500           format(/' Class',8x,'Qtot    ',
C     .            'Vh(Iface)   Vh(Bulk)   Vh+Shift    V(Rmax)')
C              endif
C   20       continue
C          endif
C        endif

      endif

   99 continue
      deallocate(qt)
      if (allocated(vesl)) deallocate(vesl)

C      if (iprint() >= 20) pause
C      call poppr

      end
