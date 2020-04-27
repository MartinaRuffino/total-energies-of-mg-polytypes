      subroutine makpfz(job,s_ctrl,s_spec,s_site,s_ham,s_pot,
     .  lidim,lihdim,iprmb,vshft,nz,zp,pfz)
C- Generates potential functions for ASA GF
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbasp nclass nspec nspin
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr z rmt lmxa hcr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pnu spec class clabel v0
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz ves
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   job   :1s digit
Ci         :0 Make ppar for each of nz energies
Ci         :1 Make ppar for one energy, b.c.'s supplied by pnu
Ci         :  (not implemented)
Ci         :10s  digit
Ci         :0 Read potential from class file (requires dclabl)
Ci         :1 Use potential from s_site%v0
Ci         :4 add 4 to 10s digit if pot->ves(ic) should be added
Ci         :  to spherical potential
Ci         :100s  digit
Ci         :1 Use hcr = rmt, rather than spec->hcr
Ci         :  NB: 2nd gen LMTO uses this switch
Ci         :2 Envelope function has fixed k.e.=0 (as in 2nd gen LMTO)
Ci         :3 combination of 1+2
Ci   lidim :number of lower+intermediate orbitals
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   dclabl:class name, packed as a real number
Ci         :dclabl(0) => classes are not named
Co Outputs
Co   pfz   :3rd generation potential parameters, in downfolding order
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer job,lidim,lihdim,nz,iprmb(*)
      double precision vshft(*)
      integer npfz
      parameter (npfz=8)
      double complex zp(nz),pfz(npfz,lihdim,*)
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
      real(8), pointer  :: p_v0(:,:)
C ... Local parameters
      logical sw,aiopot
      integer fopna,ib,ic,ifi,ii,iz,izs,is,isp,jb,jj,l,lmr,lmx,m,n0,
     .  nbasp,nclass,nglob,nl,nr,nsp,nspec,offji,job0,job1,job2
      integer iprint,lgunit,jobx
      parameter(n0=10)
      double precision z,avw,a,b,hcr(n0),rmt,vmtz,dglob
      double precision ves,pnu(n0,2)
      double complex pfzl(8,n0,nz,2)
      character*8 clabl

C      zp(1) = (-.25d0, .05d0)
C      zp(1) = zp(1) + 1d-4 * 0
C      zp(1) = (-.25d0, .00d0)

C     print *, 'z:',zp

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
        jobx = 2
      elseif (job0 == 1) then
        call rxi('makpfz: ill defined job',job)
        pnu = s_site(ib)%pnu
        jobx = 0
      else
        call rxi('makpfz: ill defined job',job)
      endif
      vmtz = s_pot%vmtz
      ves = 0

      avw = dglob('avw',0d0,0)
      nl = nglob('nl')
      call setcc(mod(nglob('lrel'),10),1d0)
      if (lidim /= lihdim) call rx('makpfz3 not set up for downfolding')
      call dpzero(pfz,npfz*lihdim*nz*2*nsp)
      call dpzero(pfzl,npfz*n0*nz*2*2)
      allocate(ibc(nbasp)); call iinit(ibc,nbasp)
      if (iprint() >= 20) then
        if (job0 == 0) call awrit3('%N MAKPFZ: parameters for '//
     .    '%?#n==1#z=#%-1j%i energies:#%n:2,1d',' ',80,lgunit(1),nz,
     .    2*nz,zp)
C        if (job0 == 1) call awrit0('%N MAKPFZ: parameters for given'//
C     .    ' boundary conditions ',' ',80,lgunit(1))
        if (iprint() >= 30) then
          write(lgunit(1),301)
  301     format(' l,E',6x,' hcr     <phi phi>',13x,'a*D',20x,'a*Ddot',
     .      19x,'phi(a)',18x,'phip(a)')
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
          if (mod(job2,2) == 1) then
            call dvset(hcr,1,n0,rmt)
          endif
          allocate(rofi(nr*2))
          call radmsh(rmt,a,nr,rofi)
          b = rmt/(dexp(a*(nr-1)) - 1)

          call rx('check pointers, makpfz')

C         Get potential for parameters, depending on 1s digit job
          clabl = ' '
          if (mod(job1,4) == 0) then
            clabl = s_site(ib)%clabel
            ifi = fopna(clabl,-1,0)
            allocate(p_v0(nr,nsp))
            sw = aiopot(nr,nsp,a,rmt,-99d0,p_v0,ifi)
            call fclose(ifi)
            if (.not. sw) call rxi('makpfz: no potential for class',ic)
          elseif (mod(job1,4) == 1) then
            p_v0 => s_site(ib)%v0
          else
            call rxi('makpfz: not implemented for job',job)
          endif

          if (job1 >= 4) then
            ves = s_pot%ves(ic)
C           call dpscop(s_pot%ves,ves,1,ic,1,1d0)
          endif
          ves = ves + vshft(ib)
          if (job2 >= 2) jobx = jobx+10

          call mkpfz2(jobx,a,avw,b,ib,ic,is,clabl,zp,lmx,nz,nl,
     .      nr,nsp,rofi,hcr,p_v0,vmtz,ves,z,pfzl)

          do  l = 0, nl-1
          if (iprmb(lmr+1) > lihdim) then
            lmr = lmr + 2*l+1
          else
            do   m = -l, l
              lmr = lmr+1
              ii = iprmb(lmr)
              do  iz = 1, nz
                do  isp = 1, nsp
                  izs = iz + nz*(isp-1)
                  call zcopy(npfz,pfzl(1,l+1,iz,isp),1,pfz(1,ii,izs),1)
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
          offji = nl*nl*(jb-ib)
          do  l = 0, nl-1
          if (iprmb(lmr+1) > lihdim) then
            lmr = lmr + 2*l+1
          else
            do   m = -l, l
              lmr = lmr+1
              ii = iprmb(lmr)
              jj = iprmb(lmr+offji)
              do  iz = 1, nz
                do  isp = 1, nsp
                  izs = iz + nz*(isp-1)
                  call zcopy(npfz,pfz(1,jj,izs),1,pfz(1,ii,izs),1)
                enddo
              enddo
            enddo
          endif
          enddo

        endif

      enddo

      deallocate(ibc)

C      if (iprint() >= 90)
C     .  call prm('pfz',0,pfz,npfz*lihdim,npfz*lihdim,nz*nsp)

      end
      subroutine mkpfz2(job,a,avw,b,ib,ic,is,clabl,zp,lmx,ne,nl,
     .  nr,nsp,rofi,hcr,v,vmtz,ves,z,pfzl)
C- Potential parameters for pot. functions and downfolding, complex e
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job   :1s digit passed to phidz.  Use
Ci         :0, boundary conditions specified val,slo
Ci         :2, boundary conditions specified by energy e (see Remarks)
Ci         :10s digit:
Ci         :1, envelope function is energy-independent, with k.e.=0
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   ib    :site index, for printout only
Ci   ic    :class index, for printout only
Ci   is    :species index, for printout only
Ci   clabl :class name, used for printout only
Ci   zp    :set of energies for which to make potential functions
Ci   lmx   :l-cutoff
Ci   ne    :number of zp energies
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
Co   pfzl  :for each energy,l,isp:
Co         :(1) = 'C' parameter for parameterization of P; see Remarks
Co         :(2) = normalization of phi
Co         :(3) = a * log derivative of phi at a = hcr
Co         :(4) = a * log derivative of phidot at a = hcr
Co         :(5) = value of wave function at a = hcr
Co         :(6) = value of phidot at a = hcr
Co         :(7) = 'del' parameter for parameterization of P; see Remarks
Co         :(8) = 'gam' parameter for parameterization of P; see Remarks
Cr Remarks:
Cr   This is an analog of nmpot2 for complex energies.
Cr   However, the potential function is returned as a set of
Cr   parameters (c,gam,del) so that derivatives of P may be taken.
Cr   P^0 is parameterized in the usual tangent form:
Cr       P^0 = (z-C)/(del + gam(z-C))
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C  ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,ib,ic,is,ne,nl,nsp,lmx,nr,npfz,n0
      parameter (npfz=8,n0=10)
      double precision a,avw,b,
     .  rofi(nr),hcr(nl),v(nr,nsp),vmtz,ves,z
      double complex zp(ne),pfzl(npfz,n0,ne,nsp)
      character*8 clabl
C ... Dynamically allocated local arrays
      complex(8), allocatable :: g(:)
      complex(8), allocatable :: gp(:)
C ... Local parameters
      logical liophi
      integer ie,k,iprint,isp,il,ir,l,nn,lgunit,awrite,job0,job1
      double precision dum,rmax,tol,hcrl
      double complex ei,sloi(5),vali(5),phia,phiap(4),phio,phi,dphi,
     .  phip,dphip,p,dla,dlap(4),ptf(5),alpha,s,kmto
      parameter (tol=1d-12)
      character*80 outs
C External calls:
      external awrit0,defcc,gintz,makpzf,nrm2gz,phidz

      job0 = mod(job,10)
      job1 = mod(job/10,10)

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
      if (job0 == 2) then
        ei = zp(ie)
        kmto = ei-vmtz
      else
        call rxi('mkpfz2: nonsensical job',job)
      endif
      if (job1 == 1) kmto = 0
      do  l = 0, lmx
        il = l+1
        hcrl = hcr(il)
        nn = 0
        nn = -1
        call phidz(job0,z,l,v(1,isp),hcrl,vmtz,rofi,nr,4,tol,ei,vali,
     .    sloi,nn,g,gp,phi,dphi,phip,dphip,p,phia,phiap,dla,dlap)
        call makpzf(104,avw,kmto,l,hcrl,sloi,vali,rmax,ptf,alpha)
C       Scale wave function to make phi^0(a)=1
        call nrm2gz(a,b,kmto,g,l,liophi,nr,hcrl,sloi,
     .    vali,rmax,rofi,dum,phio)
C       s is resulting overlap
        call gintz(g,g,a,b,nr,z,ei,l,v(1,isp),rofi,s)
        pfzl(1,il,ie,isp) = ptf(2) + ei
        pfzl(2,il,ie,isp) = s
        pfzl(3,il,ie,isp) = hcrl*dla
        pfzl(4,il,ie,isp) = hcrl*dlap(1)
        pfzl(5,il,ie,isp) = phia
        pfzl(6,il,ie,isp) = phiap(1)
        pfzl(7,il,ie,isp) = ptf(3)
        pfzl(8,il,ie,isp) = ptf(4)

C       if (iprint() >= 30) write(lgunit(1),301)
C    .    l,ei,hcrl,(pfzl(k,il,ie,isp),k=1,4)
        if (iprint() >= 30) write(lgunit(1),301)
     .    l,ie,hcrl,dble(s),(pfzl(k,il,ie,isp),k=3,6)
  301   format(2i2,2f12.6,f13.6,20f12.6)
C       call awrit5('%,2i%,2i    %,6#12d%4,6;7#12g',' ',80,
C    .    lgunit(1),l,ie,hcrl,pfzl(2,il,ie,isp),0)

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
