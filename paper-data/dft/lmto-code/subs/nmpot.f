      subroutine nmpot(job,s_ctrl,s_spec,s_site,s_ham,s_pot,
     .  lidim,lihdim,iprmb,dclabl,ppn)
C- Generates NMTO potential parameters
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: nbasp nclass nspec nspin
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: a nr z rmt lmxa hcr
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec class v0
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read: nmto kmto
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read: vmtz ves
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   job   :1s digit
Ci         :0 Make ppar for each of ham->nmto energies
Ci         :1 Make ppar for one energy, b.c.'s supplied by pnu
Ci         :  (not implemented)
Ci         :10s  digit
Ci         :0 Read potential from class file (requires dclabl)
Ci         :1 Use potential from site->ov0
Ci         :4 add 4 to 10s digit if pot->ves(ic) should be added
Ci         :  to spherical potential
Ci         :100s  digit
Ci         :1 Use hcr = rmt, rather than spec->hcr
Ci         :  NB: 2nd gen LMTO uses this switch
Ci   lidim :number of lower+intermediate orbitals
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   dclabl:class name, packed as a real number
Ci         :dclabl(0) => classes are not named
Co Outputs
Co   ppn   :3rd generation potential parameters, in downfolding order
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer job,lidim,lihdim,iprmb(*)
      integer nppn
      parameter (nppn=12)
      double precision ppn(nppn,lihdim,*),dclabl(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      integer, allocatable :: ibc(:)
      real(8), allocatable :: rofi(:)
C ... Local parameters
      logical sw,aiopot
      integer fopna,ib,ic,ifi,ii,ik,iks,is,isp,jb,jj,l,lmr,lmx,
     .  m,n0,nbasp,nclass,nkaph,nl,nmto,nr,nsp,nspec,offji,job0,
     .  job1,job2
      integer iprint,lgunit,jobx
      parameter(n0=10)
      double precision z,avw,a,b,hcr(n0),rmt,kmto(10),vmtz,dglob
      double precision ppnl(nppn,n0,n0,2),ves,pnu(n0,2)
      real(8), pointer    :: p_v0(:,:)
C     double complex kmtoz(10)
      character*8 clabl
      procedure(integer) :: nglob

C --- Unpack variables and other setup ---
      job0 = mod(job,10)
      job1 = mod(job/10,10)
      job2 = mod(job/100,10)
      nbasp = s_ctrl%nbasp
      nclass = s_ctrl%nclass
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin

C     Unpack other job-specific parameters
      if (job0 == 0) then
        nmto = s_ham%nmto
        kmto(1:6) = s_ham%kmto
        jobx = 2
      elseif (job0 == 1) then
        call rxi('nmpot: ill defined job',job)
C       pnu = s_site(ib)%pnu
        jobx = 0
        nmto = 1
      else
        call rxi('nmpot: ill defined job',job)
      endif
      vmtz = s_pot%vmtz
      ves = 0

      nkaph = nglob('nkaph')
      avw = dglob('avw',0d0,0)
      nl = nglob('nl')
      nbasp = nglob('nbasp')
      call setcc(mod(nglob('lrel'),10),1d0)
      if (lidim /= lihdim) call rx('nmpot2 not set up for downfolding')
      call dpzero(ppn, nppn*lihdim*nmto*nsp)
      call dpzero(ppnl,nppn*n0*n0*2)
C     call defi(oibc,-nbasp)
      allocate(ibc(nbasp)); call iinit(ibc,nbasp)
      if (iprint() >= 20) then
        if (job0 == 0) call awrit2('%N NMPOT: parameters for %i%-1j'//
     .    ' kinetic energies :%n:2,1d',' ',80,lgunit(1),nmto,kmto)
        if (job0 == 1) call awrit0('%N NMPOT: parameters for given'//
     .    ' boundary conditions ',' ',80,lgunit(1))
        if (iprint() >= 30) then
          write(lgunit(1),301)
  301     format(' l,E',6x,
     .      ' hcr     <phi phi>       a*D',
     .      '        a*Ddot      phi(a)      phip(a)')
        else
          write(lgunit(1),'(1x)')
        endif
      endif

C --- For each class, get potential parameters ---
      lmr = 0
      do  ib = 1, nbasp
        is = s_site(ib)%spec
        ic = s_site(ib)%class

C   ... Make pp's for new class ic
        if (ibc(ic) == 0) then
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          z = s_spec(is)%z
          rmt = s_spec(is)%rmt
          lmx = s_spec(is)%lmxa
          hcr = s_spec(is)%hcr
          if (job2 == 1) then
            call dvset(hcr,1,n0,rmt)
          endif
          allocate(rofi(nr*2))
          call radmsh(rmt,a,nr,rofi)
          b = rmt/(dexp(a*(nr-1)) - 1)

C         Get potential for parameters, depending on 1s digit job
          clabl = ' '
          if (mod(job1,4) == 0) then
            call r8tos8(dclabl(ic),clabl)
            ifi = fopna(clabl,-1,0)
            allocate(p_v0(nr,nsp))
            p_v0(1,1) = 0 ! this means rewind the pot file before loading
            sw = aiopot(nr,nsp,a,rmt,-99d0,p_v0,ifi)
            call fclose(ifi)
            if (.not. sw) call rxi('nmpot: no potential for class',ic)
          elseif (mod(job1,4) == 1) then
            p_v0 => s_site(ib)%v0
          else
            call rxi('nmpot: not implemented for job',job)
          endif

          if (job1 >= 4) then
            ves = s_pot%ves(ic)
C           call dpscop(s_pot%ves,ves,1,ic,1,1d0)
          endif

          call nmpot2(jobx,a,avw,b,ib,ic,is,clabl,kmto,pnu,lmx,nmto,nl,
     .      nr,nsp,rofi,hcr,p_v0,vmtz,ves,z,ppnl)
C          kmtoz = kmto
C          call nmpot3(jobx,a,avw,b,ib,ic,is,clabl,kmtoz,pnu,lmx,nmto,nl,
C     .      nr,nsp,rofi,hcr,w(ov0),vmtz,ves,z,ppnl)

          do  l = 0, nl-1
          if (iprmb(lmr+1) > lihdim) then
            lmr = lmr + 2*l+1
          else
            do   m = -l, l
              lmr = lmr+1
              ii = iprmb(lmr)
              do  ik = 1, nmto
                do  isp = 1, nsp
                  iks = ik + nmto*(isp-1)
                  call dcopy(nppn,ppnl(1,l+1,ik,isp),1,ppn(1,ii,iks),1)
                enddo
              enddo
            enddo
          endif
          enddo
          if (mod(job1,4) == 0) deallocate(p_v0)
          deallocate(rofi)
          ibc(ic) = ib

C   ... Copy from prior site jb with equivalent class
        else
          jb = ibc(ic)
          offji = nl*nl*nkaph*(jb-ib)
          do  l = 0, nl-1
          if (iprmb(lmr+1) > lihdim) then
            lmr = lmr + 2*l+1
          else
            do   m = -l, l
              lmr = lmr+1
              ii = iprmb(lmr)
              jj = iprmb(lmr+offji)
              do  ik = 1, nmto
                do  isp = 1, nsp
                  iks = ik + nmto*(isp-1)
                  call dcopy(nppn,ppn(1,jj,iks),1,ppn(1,ii,iks),1)
                enddo
              enddo
            enddo
          endif
          enddo

        endif

      enddo

      deallocate(ibc)

      if (iprint() >= 90)
     .  call prm('ppn',0,ppn,nppn*lihdim,nppn*lihdim,nmto*nsp)

      end
      subroutine nmpot2(job,a,avw,b,ib,ic,is,clabl,kmto,pnu,lmx,ne,nl,
     .  nr,nsp,rofi,hcr,v,vmtz,ves,z,ppnl)
C- Generates potential parameters for pot. functions and downfolding
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job   :passed to phidx.  Use
Ci         :0, boundary conditions specified val,slo
Ci         :2, boundary conditions specified by energy e (see Remarks)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   ib    :site index, for printout only
Ci   ic    :class index, for printout only
Ci   is    :species index, for printout only
Ci   clabl :class name, used for printout only
Ci   kmto  :set of kinetic energies for which to make ppars
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci         :pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci         :(not used now)
Ci   lmx   :l-cutoff
Ci   ne    :number of kmto kinetic energies
Ci   nl    :number of l's (for dimension of arrays)
Ci   nr    :number of mesh points
Ci   nsp   :=1 spin degenerate, =2 non-degenerate
Ci   rofi  :radial mesh points
Ci   hcr   :hard sphere radii
Ci   v     :spherical potential = v_true
Ci   vmtz  :muffin-tin zero
Ci   ves   :constant potential shift
Ci   z     :nuclear charge
Co Outputs:
Co   ppnl  :for each energy,l,isp:
Co         :(1) = inverse potential function
Co         :(2) = normalization of phi
Co         :(3) = a * log derivative of phi at a = hcr
Co         :(4) = a * log derivative of phidot at a = hcr
Co         :(5) = value of wave function at a = hcr
Co         :(6) = value of phidot at a = hcr
Co         :(7) = normalization of phidot (not computed now)
Co         :(8) = not used now
Cr Remarks:
Cr   This was adapted from the Stuttgart LMTO package.
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
C  ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,ib,ic,is,ne,nl,nsp,lmx,nr,nppn,n0
      parameter (nppn=12,n0=10)
      double precision a,avw,b,ppnl(nppn,n0,n0,nsp),kmto(ne),pnu(n0,2),
     .  rofi(nr),hcr(nl),v(nr,nsp),vmtz,ves,z
      character*8 clabl
C ... Dynamically allocated local arrays
      real(8), allocatable :: g(:)
      real(8), allocatable :: gp(:)
C ... Local parameters
      logical liophi
      integer ie,k,iprint,isp,il,ir,l,nn,lgunit,awrite
      double precision alpha,dla,dlap(4),dum,ptf,s,vali(5),sloi(5),rmax,
     .  phia,phiap(4),phio,ei,tol,phi,dphi,phip,dphip,p,hcrl
      parameter (tol=1d-12)
      character*80 outs
C ... External calls
      external awrit0,defdr,gintsr,makptf,norm2g,phidx,rxi

C     Add constant potential to satisfy boundary conditions
      do  isp = 1, nsp
      do  ir = 1, nr
        v(ir,isp) = v(ir,isp) + ves
      enddo
      enddo

      liophi = .false.
      rmax = rofi(nr)
      allocate(g(2*nr),gp(8*nr))

C     outs = description of this class, depending on info available
      if (iprint() >= 30) then
        outs = ' '
        il = awrite(' %?#n#%-1jsite %i, ##'//
     .              '%?#n#%-1jclass %i, ##'//
     .              '%?#n#%-1jspecies %i, ##'//
     .              '%b%b :',outs,len(outs),0,ib,ic,is,0,0,0,0,0)
        if (clabl /= ' ') outs(il+2:) = clabl
        call awrit0(outs,' ',-len(outs),lgunit(1))
      endif
      do  isp = 1, nsp
      do  ie = 1, ne
      if (job == 2) then
        ei = kmto(ie)+vmtz
      elseif (job == 2) then
        call rx('not ready')
C       nn = int(pnu(l+1,isp)) - l - 1
      else
        call rxi('nmpot2: nonsensical job',job)
      endif
      do  l = 0, lmx
        il = l+1
        hcrl = hcr(il)

        call phidx(job,z,l,v(1,isp),hcrl,vmtz,rofi,nr,4,tol,ei,vali,
     .    sloi,nn,g,gp,phi,dphi,phip,dphip,p,phia,phiap,dla,dlap)
C       dphip = (sloi(2)-phip)/rmax
        call makptf(avw,kmto(ie),l,hcrl,sloi,vali,rmax,ptf,alpha)
C       Scale wave function to make phi^0(a)=1
        call norm2g(a,b,kmto(ie),g,l,liophi,nr,hcrl,sloi,
     .    vali,rmax,rofi,dum,phio)
C       s is resulting overlap
        call gintsr(g,g,a,nr,z,ei,l,v(1,isp),rofi,s)
        ppnl(1,il,ie,isp) = 1d0/ptf
        ppnl(2,il,ie,isp) = s
        ppnl(3,il,ie,isp) = hcrl*dla
        ppnl(4,il,ie,isp) = hcrl*dlap(1)
        ppnl(5,il,ie,isp) = phia
        ppnl(6,il,ie,isp) = phiap(1)

C       if (iprint() >= 30) write(lgunit(1),301)
C    .    l,ei,hcrl,(ppnl(k,il,ie,isp),k=1,4)
        if (iprint() >= 30) write(lgunit(1),301)
     .    l,ie,hcrl,(ppnl(k,il,ie,isp),k=2,6)
  301   format(2i2,2f12.6,f13.6,10f12.6)
C       call awrit5('%,2i%,2i    %,6#12d%4,6;7#12g',' ',80,
C    .    lgunit(1),l,ie,hcrl,ppnl(2,il,ie,isp),0)

C       call prrmsh('g',rofi,g,nr,nr,1)
        enddo
      enddo
      enddo
      deallocate(g,gp)

C     Undo constant potential shift
      do  isp = 1, nsp
      do  ir = 1, nr
        v(ir,isp) = v(ir,isp) - ves
      enddo
      enddo

      end
C      subroutine nmpot3(job,a,avw,b,ib,ic,is,clabl,kmto,pnu,lmx,ne,nl,
C     .  nr,nsp,rofi,hcr,v,vmtz,ves,z,ppnl)
CC- Potential parameters for pot. functions and downfolding, complex e
CC ----------------------------------------------------------------------
CCi Inputs:
CCi   job   :passed to phidz.  Use
CCi         :0, boundary conditions specified val,slo
CCi         :2, boundary conditions specified by energy e (see Remarks)
CCi   b     :                 -//-
CCi   ib    :site index, for printout only
CCi   ic    :class index, for printout only
CCi   is    :species index, for printout only
CCi   clabl :atom name, for printout only
CCi   lmx   :l-cutoff
CCi   nl    :number of l's (for dimension of arrays)
CCi   nr    :number of mesh points
CCi   nsp   :=1 spin degenerate, =2 non-degenerate
CCi   rofi  :radial mesh points
CCi   hcr   :hard sphere radii
CCi   v     :spherical potential = v_true
CCi   vmtz  :muffin-tin zero
CCi   ves   :constant potential shift
CCi   z     :nuclear charge
CCo Outputs:
CCo   ppnl  :for each energy,l,isp:
CCo         :(1) = inverse potential function
CCo         :(2) = normalization of phi
CCo         :(3) = a * log derivative of phi at a = hcr
CCo         :(4) = a * log derivative of phidot at a = hcr
CCo         :(5) = value of wave function at a = hcr
CCo         :(6) = value of phidot at a = hcr
CCo         :(7) = normalization of phidot (not computed now)
CCo         :(8) = not used now
CCr Remarks:
CCr   This is an analog of nmpot2 for complex energies.
CC  ---------------------------------------------------------------------
C      implicit none
CC Passed variables:
C      integer job,ib,ic,is,ne,nl,nsp,lmx,nr,nppn,n0
C      parameter (nppn=12,n0=10)
C      double precision a,avw,b,pnu(n0,2),
C     .  rofi(nr),hcr(nl),v(nr,nsp),vmtz,ves,z
C      double complex kmto(ne),ppnl(nppn,n0,n0,nsp)
C      character*8 clabl
CC Local variables:
C      logical liophi
C      integer ie,k,iprint,isp,il,ir,l,nn,lgunit,og,ogp,awrite
C      double precision dum,rmax,tol,hcrl
C      double complex ei,sloi,vali,phia,phiap(4),phio,phi,dphi,phip,
C     .  dphip,p,dla,dlap(4),ptf,alpha,s
C      parameter (tol=1d-12)
C      character*80 outs
CC External calls:
C      external awrit0,defcc,gintz,makpzf,nrm2gz,phidz
C
CC     Add constant potential to satisfy boundary conditions
C      do  isp = 1, nsp
C      do  ir = 1, nr
C        v(ir,isp) = v(ir,isp) + ves
C      enddo
C      enddo
C
C      liophi = .false.
C      rmax = rofi(nr)
C      call defcc(og,2*nr)
C      call defcc(ogp,8*nr)
C
CC     outs = description of this class, depending on info available
C      if (iprint() >= 30) then
C        outs = ' '
C        il = awrite(' %?#n#%-1jsite %i, ##'//
C     .              '%?#n#%-1jclass %i, ##'//
C     .              '%?#n#%-1jspecies %i, ##'//
C     .              '%b%b :',outs,len(outs),0,ib,ic,is,0,0,0,0,0)
C        if (clabl /= ' ') outs(il+2:) = clabl
C        call awrit0(outs,' ',-len(outs),lgunit(1))
C      endif
C      do  isp = 1, nsp
C      do  ie = 1, ne
C      if (job == 2) then
C        ei = kmto(ie)+vmtz
C      elseif (job == 2) then
C        call rx('not ready')
CC       nn = int(pnu(l+1,isp)) - l - 1
C      else
C        call rxi('nmpot2: nonsensical job',job)
C      endif
C      do  l = 0, lmx
C        il = l+1
C        hcrl = hcr(il)
C        call phidz(job,z,l,v(1,isp),hcrl,vmtz,rofi,nr,4,tol,ei,vali,
C     .    sloi,nn,w(og),w(ogp),phi,dphi,phip,dphip,p,phia,phiap,dla,
C     .    dlap)
C        call makpzf(avw,kmto(ie),l,hcrl,sloi,vali,rmax,ptf,alpha)
CC       Scale wave function to make phi^0(a)=1
C        call nrm2gz(a,b,kmto(ie),w(og),l,liophi,nr,hcrl,sloi,
C     .    vali,rmax,rofi,dum,phio)
CC       s is resulting overlap
C        call gintz(w(og),w(og),a,b,nr,z,ei,l,v(1,isp),rofi,s)
C        ppnl(1,il,ie,isp) = 1d0/ptf
C        ppnl(2,il,ie,isp) = s
C        ppnl(3,il,ie,isp) = hcrl*dla
C        ppnl(4,il,ie,isp) = hcrl*dlap(1)
C        ppnl(5,il,ie,isp) = phia
C        ppnl(6,il,ie,isp) = phiap(1)
C
CC       if (iprint() >= 30) write(lgunit(1),301)
CC    .    l,ei,hcrl,(ppnl(k,il,ie,isp),k=1,4)
C        if (iprint() >= 30) write(lgunit(1),301)
C     .    l,ie,hcrl,(ppnl(k,il,ie,isp),k=2,6)
C  301   format(2i2,2f12.6,f13.6,10f12.6)
CC       call awrit5('%,2i%,2i    %,6#12d%4,6;7#12g',' ',80,
CC    .    lgunit(1),l,ie,hcrl,ppnl(2,il,ie,isp),0)
C
CC       call prrmsh('g',rofi,w(og),nr,nr,1)
C        enddo
C      enddo
C      enddo
C      call rlse(og)
C
CC     Undo constant potential shift
C      do  isp = 1, nsp
C      do  ir = 1, nr
C        v(ir,isp) = v(ir,isp) - ves
C      enddo
C      enddo
C
C      end
