      subroutine elocp(nbas,nsp,s_site,s_spec,job)
C- Make envelope parameters for extended local orbitals
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:v0
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name a nr rmt z lmxa lmxb rs3 eh3 vmtz pz orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   job   :1s  digit
Ci         : 0 do nothing ; just return
Ci         : 1 Make rsml, ehl for low lying local orbitals (stored in orbp)
Ci         : 2 Same as 1, but retain pre 7-14 bug in spin polarized case :
Ci         :   (rsml,ehl) calculated using vmt(spin 1) only
Co Outputs
Cr Remarks
Cu Updates
Cu   25 Oct 11 Started migration to f90 structures
Cu   06 Jul 05 first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,job,n0
      parameter (n0=10)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      character spid*8
      logical eloc,lvsn11
      integer ib,ibs,iclbsj,ipr,iprint,is,k,l,lmxa,lmxb,nglob,nkap0,
     .  nr,nrmx,nrspec,nspec,stdo
      parameter (nrmx=5001, nkap0=4)
      integer lh(nkap0),nkape,nkaph,idamax
      double precision z,a,rmt,rs3,eh3,vmtz,xx
      double precision rofi(nrmx),vseli(4,n0),vsel(4,n0,nbas),
     .  pnu(n0,2),pnz(n0,2),eh(n0,nkap0),rsmh(n0,nkap0),
     .  pnui(n0,2),pnzi(n0,2),ehl(n0,nkap0),rsml(n0,nkap0)
      integer,allocatable :: ips(:)
C     real(8), pointer    :: p_v0(:)


C --- Setup ---
      if (mod(job,10) == 0) return
      lvsn11 = job == 2
      call tcn('elocp')
      stdo = nglob('stdo')
      nspec = nglob('nspec')
      ipr = iprint()
      nkaph = nglob('nkaph')
      eloc = .false.
      call dpzero(vsel,4*n0*nbas)

C --- Find val, slo, K.E. for all sites ---
      do  ib = 1, nbas
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
C       p_v0 = s_site(ib)%p0
        pnz = s_site(ib)%pz
        spid = s_spec(is)%name
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        z = s_spec(is)%z
        lmxa = s_spec(is)%lmxa
        lmxb = s_spec(is)%lmxb

        if (lmxa == -1) cycle
        if (pnz(idamax(lmxb+1,mod(pnz,100d0),1),1) < 10) cycle
        eloc = .true.

        call radmsh(rmt,a,nr,rofi)
        call loctsh(1101,lvsn11,spid,z,a,nr,nr,nsp,lmxa,rofi,s_site(ib)%v0,
     .    pnu,pnz,xx,xx,vmtz,vsel(1,1,ib),rsml,ehl)
      enddo

      if (.not. eloc) goto 999
      if (ipr >= 30) write(stdo,199)
  199 format(/' elocp:')

C --- Determine shape of smooth Hankel tails for local orbitals ---
C     call defi(oips,nbas)
      allocate(ips(nbas))
C     call spackv(10,'site spec',ssite,1,nbas,ips)
C     ips(:) = s_site(1:nbas)%spec
      call sitepack(s_site,1,nbas,'spec',1,ips,xx)
C ... Loop over species containing extended local orbitals
      do  is = 1, nspec
        spid = s_spec(is)%name
        z = s_spec(is)%z
        lmxa = s_spec(is)%lmxa
        lmxb = s_spec(is)%lmxb
        if (lmxa == -1) cycle
        nrspec = iabs(iclbsj(is,ips,-nbas,nbas))
        if (nrspec == 0) cycle
        ib = iclbsj(is,ips,nbas,1)
        pnui = s_site(ib)%pnu
        pnzi = s_site(ib)%pz
        if (pnzi(idamax(lmxb+1,mod(pnzi,100d0),1),1) < 10) cycle

C   ... Average over sites within this species
        call dpzero(vseli,4*n0)
        do  ibs = 1, nrspec
          ib = iclbsj(is,ips,nbas,ibs)
          pnz = s_site(ib)%pz
          if (pnzi(idamax(lmxb+1,mod(pnzi,100d0),1),1) < 10) cycle
          call dpadd(vseli,vsel(1,1,ib),1,4*n0,1/dble(nrspec))
        enddo

C   ... Printout of input for parameters
        if (ipr >= 90/1) then
        write(stdo,261) spid
  261   format(/'  l  site    Eval        Val         Slo         K.E.',5x,'species ',a)
        do  l = 0, lmxb
          if (mod(pnz(l+1,1),100d0) < 10) cycle
          do  ibs = 1, nrspec
            ib = iclbsj(is,ips,nbas,ibs)
            write (stdo,260) l,ib,vsel(4,l+1,ib),(vsel(k,l+1,ib),k=2,4)
  260       format(i3,i4,4f12.6:a)
  262       format(i3,' avg',4f12.6)
          enddo
          write (stdo,262) l,vseli(4,l+1),(vseli(k,l+1),k=2,4)
        enddo
        endif

!      write(*,*)" no extLO optim" !JJ
!      return
C   ... Make parameters for this species
        call dpzero(rsmh,n0*nkap0)
        call dpzero(eh,n0*nkap0)
        call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkape)
        rs3 = s_spec(is)%rs3
        eh3 = s_spec(is)%eh3
        vmtz = s_spec(is)%vmtz
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        call dpzero(rsml,n0*2)
        call dpzero(ehl,n0*2)
        call radmsh(rmt,a,nr,rofi)
        call loctsh(1102,lvsn11,spid,xx,a,nr,nr,nsp,lmxa,rofi,xx,
     .    pnui,pnzi,rs3,eh3,vmtz,vseli,rsml,ehl)
        call dcopy(n0,rsml,1,rsmh(1,nkaph),1)
        call dcopy(n0,ehl, 1,eh(1,nkaph),1)
        call uspecb(1,1,s_spec,is,is,lh,rsmh,eh,nkaph)

      enddo
      deallocate(ips)

  999 continue
      call tcx('elocp')

      end
