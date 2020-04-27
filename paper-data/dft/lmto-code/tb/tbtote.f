      subroutine tbtote(s_ctrl,s_lat,s_site,nl,nsp,idxdn,nclass,
     .  nrclas,dclabl,pv,force,nbas,iax,npr,iam,npm,lscale,V0,ppmode,
     .  poly,cutppm,cutpp,qnu,sumev,alpha,entrpy,emad,ecorr,r,stni,jh,
     .  mmom,emag,e,ppdip,thrpv,tpvq,f,fe,fmax,fnou,fnom,eatm,erep,
     .  etot,efree,emg,amgm,vol)
C- Calculates tight-binding total energy and pair potential terms
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  ltb lstonr nbas ipc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  relax
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nl,nsp,idxdn,nclass,nrclas,dclabl
Ci   iax, neighbor lists; npr, see tbham
Ci   iam,npm:  see remarks
Ci   V0, parameters for pair potential, see rdtbh
Ci   qnu(1),qnu(2): initial charge and site energies (1st and 2nd moms)
Ci   sumev, sum of occupied eigenvalues
Ci   alpha, see Phys Rev B, 53, 15381 (1996)
Ci   entrpy, entropy term (actually kTS)
Ci   emad, electrostatic energy (if U=T)
Ci   ecorr, correction to UL=T total energy
Ci   r,   s,p, and d Mulliken charges
Ci   mmom,emag, local moments and magnetic energies (Stoner model)
Ci   e,   band energies for each atom
Ci   eri, work array for atom decomposition of pair energy
Ci   ppdip: atomic contributions to dipole moment
Ci   thrpv,tpvq: 3pV from bands and multipoles respectively
Ci   f,   forces from bands; fe, forces from electrostatix
Ci   fnou: force from overlap dependence of monopoles (Hubbard part)
Ci   fnom: force from overlap dependence of monopoles (Madelung part)
Ci   ppdip, accumulated atomic dipole moments (from tbesel)
Ci          RELAX=0 in ctrl is misused
Ci          to permit the dipole to be accumulated over a subset
Ci          of the sites.
Ci   vol, volume of unit cell
Co Outputs
Co   f,    force on each atom, has pair and e'static contributions added
Co   fmax, maximum force
Co   eatm, ref energy for each atom, unchanged if non-zero on entry
Co   erep, total energy contribution from pair potential
Co   etot, total energy
Co   efree, free energy: \Phi, eq (5), Phys Rev B, 53, 15381 (1996)
Co   emg,  magnetic energy (Stoner model)
Co   amgm, magnetic moment
Co   thrpv, 3pV: pV is the "internal virial" ie, (1/3)\sum r.f
Co         on output thrpv has pair and electrostatic terms added
Co   sumev is overwritten with Tr[rho H_0] if TRH is set
Cs Command-line switches
Cs   --fmax= :
Cr Remarks
Cr   iam(1,kk) and iam(2,kk) are the two classes for the kkth ME pair
Cr   k = iam(3,kk) is a pointer to tabme(*,k) which holds the
Cr   matrix elements <iam(1,kk) | H | iam(2,kk)>
Cr   npm(0,i) is number matrix elts tabulated for ith class
Cr   npm(1,i) is number matrix elts in all sites preceding ith class
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   04 Jun 08 (ATP) multipole contribution to 3PV
C ----------------------------------------------------------------------
      use mpi, only : mpi_comm_world, mpi_comm_rank
      use structures
      implicit none
C ... Passed parameters
      integer nl,nsp,nclass,nbas,niax
      parameter (niax=10)
      integer idxdn(nl,nclass),nrclas(nclass),iax(niax,2),npr(0:1,2),
     .  iam(3,2),npm(0:1,1)
      integer ppmode(1),cutppm(1),poly(1)
      double precision alat,sumev,alpha,entrpy,emad,thrpv,tpvq,erep,
     .                 etot,efree,amgm,emg,ecorr,fmax,vol
      double precision plat(3,3),V0(9,2),qnu(3,nl,nsp,nclass),
     .  r(nl,nsp,nbas),e(nbas,nsp),eri(nbas),cutpp(2,1),
     .  f(3,nbas),fe(3,nbas),fnou(3,nbas),fnom(3,nbas),
     .  stni(nclass),jh(4,nbas),dclabl(nclass),
     .  eatm(nclass,nsp),mmom(nbas),emag(nbas),ppdip(3)
      double precision getavJ
      character*8 clabli,clablj
      logical pv,force,lscale
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
C ... Local parameters
      integer i,ip,k,n,iat,ic,l,ifc,ipr,iatmx,isp,lgunit,ltb,lstnr(3),
     .  nsgrp,i1mach,ifrlx(3),parg
      double precision wk(0:3),erepi,derepi,eref,dv(3),dt(3),trh0(2),
     .  q(2),tpvpr,dp,fabs2,fmax2,fmaxx,d1mach,dsum,qtot(2),pos(3),lerep
     .  ,pdip(3),ptot(3),dipole,dip1,dip2,q0,Dq,J,tpv,ftot(3),dderep,pb
      logical cmdopt,bittst,lov,trh,TBU,UL,charge,stoner,MOL,lftot
      character*120 outs
      integer, pointer :: ipc(:)
      integer :: ierr, iproc
      logical :: root

      call tcn('tbtote')

      call mpi_comm_rank(mpi_comm_world, iproc, ierr)
      root = iproc == 0

      ltb = s_ctrl%ltb
      lstnr = s_ctrl%lstonr
      nbas = s_ctrl%nbas
      alat = s_lat%alat
      plat = s_lat%plat
      nsgrp = s_lat%nsgrp
      ipc => s_ctrl%ipc

      lov    = bittst(ltb,1)
      trh    = bittst(ltb,2**10)
      charge = bittst(ltb,2**11)
      UL     = bittst(ltb,2**15)
      TBU    = bittst(ltb,2**13)
C     Uav = (.not. bittst(ltb,2**14))
      MOL    = bittst(ltb,2**18)
      stoner = lstnr(1) /= 0
      if (UL .or. TBU) then
        emad = ecorr
        trh = .true.
      endif
      amgm = 0d0
      call getpr(ipr)

C --- Get reference energy ---
      eref = 0d0
      do  ic = 1, nclass
        if (abs(eatm(ic,1)) < 10*d1mach(3) .and.
     .      abs(eatm(ic,nsp)) < 10*d1mach(3)) then
          eatm(ic,1) = 0d0
          if (nsp == 2) eatm(ic,2) = 0d0
          do  isp = 1, nsp
            do  l = 1, nl
              if (idxdn(l,ic) <= 1) eatm(ic,isp) = eatm(ic,isp)
     &                                + qnu(1,l,isp,ic)*qnu(2,l,isp,ic)
            enddo
          enddo
        endif
        eref = eref + eatm(ic,1)*nrclas(ic)
        if (nsp == 2) eref = eref + eatm(ic,2)*nrclas(ic)
      enddo

C --- Loop over all atoms ---
      fmax2 = 0d0
      erep = 0d0
      tpvpr = 0d0
      call dpzero(eri,nbas)
      if (ipr > 20 .and. ipr <= 30 .and. force) write (*,102)
      do  iat = 1, nbas
C --- Get pair potential contributions ---
        ic = ipc(iat-1+1)
        n = npr(1,iat)
        if (force) dv = 0d0
        lerep = 0d0
        do  i = 2, npr(0,iat)
          call dlmn(nbas,plat,s_lat%pos,iax(1,i+n),wk)
          call meptr(ipc(iax(1,i+n)),ipc(iax(2,i+n)),iam,npm,k)
          if (k == 0) cycle
          call makvpp(alat,wk(0)*alat,v0(1,k),.false.,lscale,
     .           ppmode(k),poly(k),cutppm(k),cutpp(1,k),
     .           erepi,derepi,dderep,dp)
          if (force) dv = dv + derepi*wk(1:3)
          lerep = lerep + erepi
          tpvpr = tpvpr + dp
        enddo
        eri(iat) = 0.5d0*lerep
        erep = erep + lerep
C --- Add to and print forces ---
        if (force) then
          dt = f(1:3,iat) + fe(1:3,iat) + fnou(1:3,iat)
     .                                            + fnom(1:3,iat) + dv
          fabs2 = sum(dt*dt)
          if (fabs2 > fmax2) then
            fmax2 = fabs2
            iatmx = iat
          endif

          if (ipr > 20) then
            ic = ipc(iat)
            call r8tos8(dclabl(ic),clabli)
            call dpscop(s_lat%pos,pos,3,3*iat-2,1,1d0)
            if (ipr > 20 .and. ipr <= 30)
     .                     print 101, iat,clabli,(pos(ifc),ifc=1,3),dt
            if (ipr > 30) then
               if (UL .or. TBU) then
                  if (lov) then
                     if (root) print 103,iat,clabli,(pos(ifc),ifc=1,3),
     .                  (f(ifc,iat),ifc=1,3),(fe(ifc,iat),ifc=1,3),
     .                  (fnou(ifc,iat),ifc=1,3),(fnom(ifc,iat),ifc=1,3),
     .                  dv,dt
                  else
                     if (root) print 104,iat,clabli,(pos(ifc),ifc=1,3),
     .                  (f(ifc,iat),ifc=1,3),(fe(ifc,iat),ifc=1,3),
     .                  dv,dt
                  endif
               else
                  if (root) print 100, iat,clabli,(pos(ifc),ifc=1,3),
     .                   (f(ifc,iat),ifc=1,3),dv,dt
               endif
            endif
          end if
          f(1:3,iat) = dt
        endif
      enddo
      tpvpr = -tpvpr
      fmax = sqrt(fmax2)

      if (force .and. fmax > 1d-6 .and. ipr > 20 .and. root) then
        call r8tos8(dclabl(ipc(iatmx)),clablj)
        call awrit2('%x= %;6,6d on atom %i ('//
     .    clablj//'%a)',outs,80,0,fmax,iatmx)
        call awrit0('%N Maximum force'//outs,' ',-80,lgunit(1))
      endif
      if (root) print *, ' '

  100 format(/' Forces on atom ',i4,4x,'Species: ',a4/
     .  '  Coordinates:',3f14.8/
     .  '  From bands :',3f14.8/'  From pairs :',3f14.8/
     .  '  Total      :',3f14.8)
  103 format(/' Forces on atom ',i4,4x,'Species: ',a4/
     .  '  Coordinates        :',3f14.8/
     .  '  From bands         :',3f14.8/
     .  '  From e''stx         :',3f14.8/
     .  '  From overlap (U)   :',3f14.8/
     .  '  From overlap (Mad) :',3f14.8/
     .  '  From pairs         :',3f14.8/
     .  '  Total              :',3f14.8)
  104 format(/' Forces on atom ',i4,4x,'Species: ',a4/
     .  '  Coordinates        :',3f14.8/
     .  '  From bands         :',3f14.8/
     .  '  From e''stx         :',3f14.8/
     .  '  From pairs         :',3f14.8/
     .  '  Total              :',3f14.8)
  102 format(/'   Site',16x,'pos',30x,'force')
  101 format(i4,' ',a,3f10.5,3x,3f10.5)

c Check the total force
      if (force) then
c       ftot = 0d0
c       do iat = 1, nbas
c         ftot = ftot + f(1:3,iat)
c       enddo
c       print '(a,3g12.4)',' tbtote: The sum of all forces =',ftot
        ftot = sum(f,dim=2)
        lftot = maxval(abs(ftot)) > 1d-8
        if (ipr >= 20 .or. lftot) then
          if(root)print '(a,4g12.4)',' tbtote: sum of all forces =',ftot
          if(lftot) call rx0(' tbtote: forces do not sum up to zero')
        endif
      endif

C --- Print atom- or class-specific charges and energies ---
      emg = 0d0
      if (trh .or. charge .or.  UL .or. TBU .or. stoner) then
        trh0 = 0d0
        q = 0d0
        pdip = 0d0
        do  iat = 1, nbas
          call dpscop(s_lat%pos,pos,3,3*iat-2,1,1d0)
          ic = ipc(iat)
          if (trh .or. charge .or. UL .or. TBU) then
            qtot(1) = dsum(nl,r(1,1,iat),1)
            q(1) = q(1) + qtot(1)
            qtot(2) = 0d0
            if (nsp == 2) then
              qtot(2) = dsum(nl,r(1,2,iat),1)
              q(2) = q(2) + qtot(2)
            endif
C  --- Dipole moment ---
            ifrlx = s_site(iat)%relax
            q0 = dsum(nl,qnu(1,1,1,ic),3)
            if (nsp == 2) q0 = q0 + dsum(nl,qnu(1,1,2,ic),3)
            Dq = (qtot(1)+qtot(2)) - q0
            if (ifrlx(1) == 1) pdip = pdip + Dq*alat*pos
          endif
          if (ipr > 20) then
          call r8tos8(dclabl(ic),clabli)
          if ((nsp == 2 .and. ipr > 20 .and. root)
     .      .or. ipr > 30) write (*,310) iat,clabli
C  --- Mulliken charges ---
          if ( (trh .or. charge .or. UL .or. TBU)
     .      .and. ((nsp == 2 .and. ipr > 20)
     .            .or. ipr > 30) .and. root ) then
            if (nsp == 2) then
              if (nl == 3)
     .        write (*,400) nl,((r(l,1,iat)+r(l,2,iat)), l = 1, nl),
     .           dsum(nl,r(1,1,iat),1)+dsum(nl,r(1,2,iat),1)
              if (nl == 2)
     .        write (*,401) nl,((r(l,1,iat)+r(l,2,iat)), l = 1, nl),
     .           dsum(nl,r(1,1,iat),1)+dsum(nl,r(1,2,iat),1)
              if (nl == 1)
     .        write (*,402) nl,((r(l,1,iat)+r(l,2,iat)), l = 1, nl),
     .           dsum(nl,r(1,1,iat),1)+dsum(nl,r(1,2,iat),1)
              write (*,410)    ((r(l,1,iat)-r(l,2,iat)), l = 1, nl)
              if (nl > 2 .and. (UL .or. TBU)) then
                call awrit2('  d-band magnetic energy: I=%d, E_X=%d',
     .            ' ',120,i1mach(2),jh(3,iat),
     .            -0.25d0*(r(3,1,iat)-r(3,2,iat))**2*stni(ic))
              endif
            else
              if (nl == 3)
     .        write (*,450) nl,(r(l,1,iat), l = 1, nl),
     .                      dsum(nl,r(1,1,iat),1)
              if (nl == 2)
     .        write (*,451) nl,(r(l,1,iat), l = 1, nl),
     .                      dsum(nl,r(1,1,iat),1)
              if (nl == 1)
     .        write (*,452) nl,(r(l,1,iat), l = 1, nl),
     .                      dsum(nl,r(1,1,iat),1)
            endif
          endif
          end if
C --- Band energies, pair energies ---
          trh0(1) = trh0(1) + e(iat,1)
          e(iat,1) = e(iat,1) - eatm(ic,1)
          if (nsp == 2) then
            trh0(2) = trh0(2) + e(iat,2)
            e(iat,2) = e(iat,2) - eatm(ic,2)
          endif
          if (ipr > 20 ) then
            if (root) write (*,560) pos, Dq
          endif
C --- Stoner model ---
          if (stoner .and. ipr >= 20) then
            amgm = amgm + mmom(iat)
            emg  = emg  + emag(iat)
            if (root) write (*,575) mmom(iat),emag(iat)
          endif
C  --- Magnetic energy and moment ---
          if (nsp == 2 .and. .not. stoner) then
            amgm = q(1) - q(2)
            if (UL .or. TBU) then
              J = getavJ(nl,jh(1,iat),idxdn,ic)
              emg = emg - 0.25d0 * J * mmom(iat)*mmom(iat)
            endif
          endif
        enddo
      endif

C --- dipole moment ---
      if ((charge .or.  UL .or. TBU) .and. MOL) then
        pdip  = -2.541748d0 * pdip
        ppdip = -2.541748d0 * ppdip

        ptot = pdip + ppdip

        dipole = sqrt(sum(ptot*ptot))
        ptot = ptot/dipole

        dip1 = sqrt(sum(pdip*pdip))
        if (dip1 > d1mach(3)) pdip = pdip/dip1

        dip2 = sqrt(sum(ppdip*ppdip))
        if (dip2 > d1mach(3)) ppdip = ppdip/dip2

        if (root) then
        print *, ' '
        print *, ' Molecular dipole moment in Debye ... (unit vector)'
        call awrit2('  From point charges:  %d (%3:1d)',' ',128,
     .              i1mach(2),dip1,pdip)
        if (dip2 > d1mach(3)) then
          call awrit2('  From atomic dipoles: %d (%3:1d)',' ',128,
     .                i1mach(2),dip2,ppdip)
          call awrit2('  Total moment:        %d (%3:1d)',' ',128,
     .                i1mach(2),dipole,ptot)
        endif
        end if
      endif

  310 format(/' Atom ',i4,'   Species ',a4)
  400 format('  Charges: NL=',i1,': ',3(1x,f10.6),' (Total: ',f10.6,')')
  401 format('  Charges: NL=',i1,': ',2(1x,f10.6),' (Total: ',f10.6,')')
  402 format('  Charges: NL=',i1,': ',1(1x,f10.6),' (Total: ',f10.6,')')
  410 format(9x,'Moment: ',4(1x,f10.6))
  450 format('  Charges: NL=',i1,': ',3(1x,f10.6),' (Total: ',f10.6,')')
  451 format('  Charges: NL=',i1,': ',2(1x,f10.6),' (Total: ',f10.6,')')
  452 format('  Charges: NL=',i1,': ',1(1x,f10.6),' (Total: ',f10.6,')')
  560 format('  POS=',3f10.6,'  Dq/e=',f10.6)
  575 format('  Stoner Model:  MMOM=',f14.8,'   EMAG=',f14.8)

C --- Print total charges ---
      if (root) then
      if (ipr > 10) print *
      if (ipr >= 20) then
        if (nsp == 2) then
          if (trh) then
            write (*,600) q(1)+q(2),amgm,trh0(1),trh0(2),trh0(1)+trh0(2)
          elseif (charge) then
            write (*,600) q(1)+q(2),amgm
          endif
        else
          if (trh) then
            if (stoner) then
              write (*,675) q(1),amgm,trh0(1)
            else
              write (*,650) q(1),trh0(1)
            endif
          elseif (charge .or. UL .or. TBU) then
            if (stoner) then
              write (*,675) q(1),amgm
            else
              write (*,650) q(1)
            endif
          endif
        endif
      endif
      end if
C --- Total energy ---
      erep = erep/2
      tpvpr = tpvpr/2
      if (UL .or. TBU) then
        etot = trh0(1) + emad - eref + erep
        if (nsp /= 1) etot = etot + trh0(2)
      else
        etot = sumev + emad - eref + erep
      endif

      if (stoner) etot = etot + emg
      efree = etot - (1d0 - alpha)*entrpy

      if (ipr >= 10 .and. root) then
        if (UL .or. TBU) then
          if (nsp == 2) then
            write(*,700) sumev,emad,erep,eref,emg,amgm
          else
            write(*,700) sumev,emad,erep,eref
          endif
        else
          if (stoner) then
            write(*,750) sumev,erep,eref,emg
          else
            write(*,750) sumev,erep,eref
          endif
        endif
        write(*,800) etot, efree
      endif
      if (trh) sumev = trh0(1) + trh0(2)


C --- Pressure ---
      tpv = tpvpr
      if (pv) tpv = tpv + thrpv + tpvq
C --- pressure in bar ---
      pb = 147.116d0 * tpv / (3d0 * vol)

      if (ipr >= 10 .and. root) then
        if (pv) then
          write(*,900) tpvpr,thrpv,tpvq,tpv,pb
        else
          write(*,900) tpvpr
        endif
      print *
      endif
      thrpv = tpv

  600 format( ' Tr[rho]         total   :  ',f16.8,
     .       /'                moment   :  ',f16.8,:,
     .       /' Tr[rho][H_0]       up   :  ',f16.8,
     .       /'                  down   :  ',f16.8,
     .       /'                 total   :  ',f16.8)
  650 format( ' Tr[rho]                 :  ',f16.8,:,
     .       /' Tr[rho][H_0]            :  ',f16.8)
  675 format( ' Tr[rho]                 :  ',f16.8,:,
     .       /' Stoner magnetic moment  :  ',f16.8,:,
     .       /' Tr[rho][H_0]            :  ',f16.8)
  700 format( ' band structure energy   :  ',f16.8,
     .       /' E_2                     :  ',f16.8,:,
     .       /' pair potential energy   :  ',f16.8,:,
     .       /' reference energy        :  ',f16.8,:,
     .       /' Stoner magnetic energy  :  ',f16.8,:,
     .       /' Magnetic moment         :  ',f16.8)
  750 format( ' band structure energy   :  ',f16.8,
     .       /' pair potential energy   :  ',f16.8,:,
     .       /' reference energy        :  ',f16.8,:,
     .       /' Stoner magnetic energy  :  ',f16.8)
  800 format( ' total energy            :  ',f16.8,:,
     .       /' free energy             :  ',f16.8)
  900 format( ' 3PV              pair   :  ',f16.8,:,
     .       /'                 bands   :  ',f16.8,
     .       /'                 charges :  ',f16.8,
     .       /'                 total   :  ',f16.8, ' (',f12.8 ,' bar)' )

C --- protect against relaxation blowing up ---
      if (cmdopt('--fmax=',7,0,outs)) then
        ip = 7
        call skipbl(outs,len(outs),ip)
        k = parg(' ',4,outs,ip,len(outs),' ',1,1,n,fmaxx)
        call rxx(k /= 1,' TBTOTE: error parsing --fmax=')
        if (fmax > fmaxx) then
          if (ipr > 10) then
            if (root) call awrit2(' TBTOTE fmax=%d > fmaxx=%d',' ',128,
     .        i1mach(2),fmax,fmaxx)
          endif
          call rx(' ')
        endif
      endif

      call tcx('tbtote')
      end
