      subroutine loctsh(mode,lvsn11,spid,z,a,nr,nrmt,nsp,lmxb,rofi,v,pnu,pnz,
     .  rs3,eh3,vmtz,vsel,rsml,ehl)
C- Fit value and slope of local orbitals to smoothed Hankel
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit concerns fit to low-lying local orbitals
Ci         : 0 do not fit low-lying local orbitals
Ci         : 1 Generate val,slo, and K.E. from potential
Ci         : 2 Fit low-lying local orbitals attempting to match rs,eh
Ci         :   to val, slope and K.E.
Ci         : 3 combination of 1+2
Ci         : 4  -----
Ci         :10s digit concerns fit to high-lying local orbitals
Ci         :   NB: this branch is not implemented
Ci         : 0 do not fit high-lying local orbitals
Ci         : 1 Generate val,slo, and K.E. from potential
Ci         : 2 -----
Ci         : 4 Vary rsm to fit log deriv of high-lying local orbitals
Ci         :   to phi using specified eh3
Ci         : 5 combination of 1+4
Ci         :100s digit concerns constraint on rsm
Ci         : 0 no constraint on maximum rsm
Ci         : 1 constrain rsm to be <= rmt
Ci         :1000s digit deals with spin-polarized case:
Ci         : 0 Treat spins separately
Ci         : 1 Compute ehl,rsml for average potential
Ci         :10000s digit 1 plot wave functions (not implemented)
Ci   lvsn11:T => Use spin 1 pnz to construct LO.
Ci         :F Use spin average pnz
Ci         :Prior to Oct 18, lvsn11=T => also used spin-1 potential
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
Ci   nr    :number of radial mesh points in potential sphere
Ci   nrmt  :number of points between 0..rmt
Ci   lmxb  :l-cutoff for basis
Ci   pl    :boundary conditions for valence wavefunctions.
Ci   pnz   :boundary conditions for local orbital. pnz=0 -> no loc. orb.
Ci         :10s digit controls how local orbital included in hamiltonian
Ci         :10s digit nonzero -> smooth Hankel tail is attached.
Ci   rs3   :minimum allowed smoothing radius in attaching Hankel tails
Ci         :to local orbitals
Ci   eh3   :Hankel energy when attaching Hankel tails to high-lying
Ci         :local orbitals
Ci   vmtz  :parameter used in attaching Hankel tails to local orbitals
Ci         :It is used as a constant shift to Hankel energies for the
Ci         :fitting of local orbitals to Hankel tails. Thus vmtz
Ci         :is an estimate for the potential at the MT radius.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   spid  :species label
Cio Inputs/Outputs
Cio  vsel  :value, slope, K.E. of radial w.f.
Cio        :output if 1s digit (or 10s digit) has 1's bit set
Cio        :input  if 1s digit (or 10s digit) has 1's bit clear
Co Outputs
Co   ehl   :energies of smoothed Hankel tail for local orbital
Co   rsml  :smoothing radii for smoothed Hankel tail of local orbital
Cl Local variables
Cl   rmt   :muffin-tin radius, in a.u.
Cl   modei :1s digit mode for low-lying local orbitals
Cl         :10s digit mode for high-lying local orbitals
Cb Bugs
Cb   Te fully rel 4d state shows kink near RMT, from GW cdtetest.tar.gz
Cr Remarks
Cu Updates
Cu   03 Aug 15 Bug fix: previously, vmine used to make rsml,ehl
Cu             is calculated from the spin-averaged potential.
Cu             To retain the old convention (spin-1 potential) set lvsn11=.true.
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   11 Jan 10 Patch phidx call for deep semicore states
Cu   06 Jul 05 generation of (val,slo,k.e.) and fitting split into
Cu             independent operations
Cu   27 Jun 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lvsn11
      integer lmxb,nr,nrmt,nsp,n0,mode
      double precision a,z,rs3,eh3,vmtz
      parameter (n0=10)
      double precision rofi(*),v(nr,nsp),pnz(n0,nsp),pnu(n0,nsp)
      double precision rsml(n0,2),ehl(n0,2),vsel(4,n0)
C ... Dynamically allocated local arrays
      real(8), allocatable :: g(:)
      real(8), allocatable :: gp(:)
C ... Local parameters
      character spid*8,orbit(2)*4,flg(2)*1
      integer ipr,iprint,i,j,l,info,nglob,stdo,
     .  mode0,mode1,mode2,mode3,mode4,modei,loclo,nfit,IPRT,isw
      integer nrx,nrbig,ir,lp1
      double precision dphi,dphip,eh,eval,p,phi,phip,
     .  pnul,rsm,rsmin,rsmax,rmt,ekin,h,sum1,sum2,qrmt
C     emin and emax are the maximum allowed ranges in Hankel energies
C     for the fitting of local orbitals to Hankel tails.
      double precision emin,emax,vmine
      parameter (nrx=1501)
      parameter (emin=-10d0,emax=-.02d0,IPRT=20)
      double precision rbig,rofix(nrx),rwgtx(nrx),xi(0:20),pnbar(n0,2)
C     For plotting wave functions
C     integer itab(n0,2)
C     double precision pl(n0,nsp)
      data orbit/'low','high'/
      data flg/'*',' '/

C     call prrmsh('loctsh v',rofi,v,nr,nr,nsp)

C --- Setup ---
      ipr = iprint()
      stdo = nglob('stdo')
C     stdl = nglob('stdl')
      if (lmxb > 8) call rx('loctsh:  lmax too large')
      rmt = rofi(nrmt)
C     b = rmt/(dexp(a*nrmt-a)-1d0)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      mode4 = mod(mode/10000,10)
      nfit = 0

      allocate(g(2*nr),gp(2*nr*4))

      pnbar(:,1) = (pnz(:,1) + pnz(:,nsp))/2
      pnbar(:,2) = (pnz(:,1) + pnz(:,nsp))/2
      if (lvsn11) pnbar = pnz

      do  i = 1, nsp

C --- Loop over local orbitals ---
      do  l = 0, lmxb

C       itab(l+1,i) = 0
        rsml(l+1,i) = 0
        ehl(l+1,i) = 0

C       pl(l+1,i) = pnu(l+1,i)
C       konfig = mod(pnz(l+1,i),10d0)

C       Skip all but local orbitals with tails attached
        if (mod(pnz(l+1,i),100d0) < 10) cycle

C       Case local orbital deeper than valence
        if (int(pnu(l+1,i)-1) == int(mod(pnz(l+1,i),10d0))) then
          loclo = 1
          if (mode0 == 0) cycle
          modei = mode0

C       Case local orbital higher than the valence state
        elseif (int(pnu(l+1,i)+1) == int(mod(pnz(l+1,i),10d0))) then
          loclo = 0
          if (mode1 == 0) cycle
          modei = mode1

C       Local orbital neither low nor high: error
        else
          call fexit3(-1,111,' Exit -1 loctsh, l=%i:  sc PZ=%d '//
     .      'incompatible with valence P=%;3d',l,pnbar(l+1,i),pnu(l+1,i))
        endif

        pnul = mod(pnbar(l+1,i),10d0) + 100*(int(pnbar(l+1,i))/100) ! For printout

C   ... Initial printout
        nfit = nfit + 1
        if (nfit == 1 .and. modei > 1) then
          call info2(20,1,0,' Fit local orbitals to sm hankels, species '
     .      //spid//'%a, rmt=%;7g',rmt,0)
          if (ipr >= IPRT) write (stdo,261)
        endif

C   ... Make value, slope, kinetic energy
        if (mod(modei,2) == 1) then
C     ... Overwrite V+, V- with (V+ + V-)/2, (V+ - V-)/2,
          if (nsp == 2 .and. mode3 == 1) then
            call dsumdf(nr,0.5d0,v,0,1,v(1,2),0,1)
          endif

C     ... Wave function and potential parameters at MT sphere
          j = 0
          if (int(mod(pnbar(l+1,i),10d0)) < int(pnu(l+1,1))) j = 100
          call makrwf(1000+j,z,rmt,l,v(1,i),a,nrmt,rofi,pnbar(l+1,i),-0.5d0,4,
     .      g,gp,eval,phi,dphi,phip,dphip,p)
C         call prrmsh('gz',rofi,g,nr,nr,2)

          vmine = eval - (v(nrmt,1)-2*z/rmt)
          vsel(1,l+1) = phi
          vsel(2,l+1) = dphi
          vsel(3,l+1) = vmine
          vsel(4,l+1) = eval

C     ... Restore V+, V-
          if (nsp == 2 .and. mode3 == 1) then
            call dsumdf(nr,1d0,v,0,1,v(1,2),0,1)
          endif

C     ... v11 Bug ... retain for backwards compatibility
C          if (lvsn11) then
C            vmine = eval - (v(nrmt,i)-2*z/rmt)
C            vsel(3,l+1) = vmine
C          endif

        else

          phi   = vsel(1,l+1)
          dphi  = vsel(2,l+1)
          vmine = vsel(3,l+1)
          eval  = vsel(4,l+1)

        endif

C   ... Set conditions on envelope functions
        if (modei > 1) then
        rsmin = rs3
        rsmax = 5
        if (mode2 == 1) rsmax = rmt
        if (eval < emin)
     .    call rx1('increase emin in loctsh: eval=%;4g',eval)

C   ... Match val,slo, and K.E. of Hankel to phi,dphi,vmine
        if (modei == 2 .or. modei == 3) then
          rsm = 0
          eh = min(eval-vmtz,emax)
C         call pshpr(110-10)
          call mtchre(103,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,
     .      vmine,dphi,rsm,eh,ekin,info)
C         call poppr

C   ... Vary rsm to match sm Hankel to phi,dphi
        elseif (modei == 4 .or. modei == 5) then
          call rx('this branch not checked')
          rsm = rsmin
          call mtchre(100,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,phi,
     .      dphi,rsm,eh,ekin,info)

        else
          call rxi('loctsh: not implemented fitting mode=',modei)
        endif

C   ... Printout of fit functions
        if (ipr >= IPRT) then
          call radext(11,nr,nrx,2d0,a,rmt,nrbig,rbig,rofix,rwgtx)
C         make sphere charge for r>rmt for printout
          lp1 = l+1
          sum1 = 0
          do  ir = 1, nr
C           h = r*radial part of sm. Hankel
            call hansmr(rofix(ir),eh,1/rsm,xi,l)
            h = xi(l)*(rofix(ir)**lp1)
            sum1 = sum1 + rwgtx(ir)*h**2
          enddo
          sum2 = 0
          do  ir = nr, nrbig
C           h = r*radial part of sm. Hankel
            call hansmr(rofix(ir),eh,1/rsm,xi,l)
            h = xi(l)*(rofix(ir)**lp1)
            sum2 = sum2 + rwgtx(ir)*h**2
          enddo
          qrmt = sum2/(sum1+sum2)

          write (stdo,260) l,orbit(2-loclo),pnul,eval,vmine,rsm,eh,qrmt,ekin,
     .      flg(2-isw(dabs(ekin-vmine) > 1d-5))
  261   format('  l  type    Pnu      Eval        K.E.',
     .    7x,'Rsm       Eh      Q(r>rmt)    Fit K.E.')
  260   format(i3,2x,a,f8.3,2f12.6,3f10.5,f12.6,a1)
        endif
C       pause

C       itab(l+1,i) = 1
        rsml(l+1,i) = rsm
        ehl(l+1,i) = eh
        endif


      enddo
      if (mode3 /= 0) exit
      enddo

C --- Make plot file ---
      if (mode4 == 1) then
        call rx('loctsh: plotting not ready')
C        if (nsp == 2) call rx('loctsh is not spinpol yet')
C        allocate(psi,nr*nlxb)
C        call popta5(lmxb,rsml,ehl,itab,z,pl,rmax,rmt,nr,nrmt,
C     .     rofi,psi,v,g,a,b,spid)
C        deallocate(psi)
      endif

      deallocate(g,gp)

      end

