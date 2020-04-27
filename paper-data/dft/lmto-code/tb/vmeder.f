      subroutine vmeder(memode,ider,ioff,ilme,ntermx,tabme,decay,dist,
     .                  il1,il2,vme,vmed,vmedd)
C- Matrix elements and their spatial derivatives (optional) for a single
C  pair of species at given distance
C ----------------------------------------------------------------------
Ci Inputs
Ci   memode: scaling law for matrix elements (see Remarks)
Ci   ider  : switch as to how many spatial derivatives of ME requested
Ci           0 no derivatives, only ME value
Ci           1 ME value + first derivative
Ci           2 ME value + first & second derivatives
Ci   ioff  : offset to tables with species parameters tabme and decay
Ci           (including spins)
Ci   ilme  : = nlme(nl), number of TB MEs
Ci   ntermx: leading dimension of tabme
Ci   tabme : table containing the ME parameters for each pair of species
Ci           the meaning of parameters depends on memode (see Remarks)
Ci   decay : decay of MEs with distance for memode = 2 or 3
Ci   dist  : distance at which MEs are sought
Ci   il1,il2: MEs evaluated between il1 and il2, eg, vme(il1:il2)
Cl Local variables
Co Outputs
Co   vme   : MEs at distance dist, V_i(dist)
Co   vmed  : V'_i(dist) if ider >= 1, otherwise not referenced
Co   vmedd : V"_i(dist) if ider = 2, otherwise not referenced
Cr Remarks
Cr   * The following types of ME scaling law are currently implemented:
Cr     memode = 0,  fixed MEs
Cr     memode = 1,  Harrison universal MEs
Cr     memode = 2,  exponential decay
Cr     memode = 3,  power decay
Cr     memode = 4,  ME = \sum_i=1,3 a_i d^b_i exp(-c_i d), the ordering is:
Cr                  a_1 b_1 c_1 a_2 ... c_3 for ss-sigma, then sp-sigma, etc
Cr     memode = 5,  Goodwin-Skinner-Pettifor,
Cr                  v_ij (r0/d)^n exp[n (-{d/rc}^nc + {r0/rc}^nc)]
Cr                  ordering: v_ij n nc r0 rc for ss-sigma, etc
Cr     memode = 6, (nl=1) Slater-Koster + Srini's extra term
Cr                  NB: NOT implemented for spin pol.
Cr     memode = 7,  a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc)
Cr                  ordering: a b c d0 for ss-sigma, etc
Cr    memode >= 10, use canonical TB Hamiltonian (under development)
Cu Updates
Cu   19 Apr 11 (SL)  first created from tbham.f
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer, intent(in) :: memode,ider,ntermx,ioff,ilme,il1,il2
      double precision, intent(in) :: decay(*),tabme(ntermx,*),dist
      double precision, intent(out) :: vme(ilme),vmed(ilme),vmedd(ilme)
C Local parameters
      integer il,im
      double precision wk1(ilme),wk2(ilme)
      double precision wk3(3,3,ilme),wk4(ntermx,ilme)
      double precision n,nc,r0,rc,A
      double precision d1,e,de,dle,frac
      double precision d1mach

      if (il1 < 1 .or. il2 > ilme .or. il1 > il2) then
        call info5(10,1,0,' vmeder: il1 = %i'//
     .      ' il2 = %i ilme = %i',il1,il2,ilme,0,0)
        call rx(' vmeder: check bounds')
      endif

      if (ider /= 0) d1 = 1d0/dist

C --- Fixed MEs ---
      if (memode == 0 .or. memode == 6) then
        vme(il1:il2) = tabme(1,ioff+il1-1:ioff+il2-1)
        if (ider >= 1) then
          vmed(il1:il2) = 0d0
          if (ider == 2) vmedd(il1:il2) = 0d0
        endif

C --- Harrison universal MEs ---
      elseif (memode == 1) then
        call rxx(ilme > 4,
     .    ' vmeder: No universal Hamiltonian for l > 1')
        call ume(ilme,dist,wk1)
        vme(il1:il2) = wk1(il1:il2)
        if (ider >= 1) then
          vmed(il1:il2) = -2d0*vme(il1:il2)*d1
          if (ider == 2) vmedd(il1:il2) = -3d0*vmed(il1:il2)*d1
        endif


C --- Exponential or power decay ---
      elseif (memode == 2 .or. memode == 3) then
c       call dcopy(ilme,tabme(1,ioff),1,wk1,1)
        wk1(il1:il2) = tabme(1,ioff+il1-1:ioff+il2-1)
        call dcopy(ilme,decay(ioff),1,wk2,1)
        do il = il1, il2
          if (memode == 2) then
            vme(il) = wk1(il) * dexp(-dist*wk2(il))
            if (ider >= 1) then
              vmed(il) = -wk2(il)*vme(il)
              if (ider == 2) vmedd(il) = -wk2(il)*vmed(il)
            endif
          else
            vme(il) = wk1(il) * dist**(-wk2(il))
            if (ider >= 1) then
              vmed(il) = -wk2(il)*vme(il)*d1
              if (ider == 2) vmedd(il) = -(wk2(il)+1d0)*vmed(il)*d1
            endif
          endif
        enddo

C --- ME = \sum_i=1,3 a_i d^b_i exp(-c_i d) ---
      elseif (memode == 4) then
        call dcopy(9*ilme,tabme(1,ioff),1,wk3,1)
        vme(il1:il2) = 0d0
        if (ider >= 1) then
          vmed(il1:il2) = 0d0
          if (ider == 2) vmedd(il1:il2) = 0d0
        endif
        do il = il1, il2
          do im = 1, 3
            if (dabs(wk3(1,im,il)) >= 1d-12) then
              e = wk3(1,im,il)*dist**wk3(2,im,il) *
     .            dexp(-dist*wk3(3,im,il))
              vme(il) = vme(il) + e
              if (ider >= 1) then
                dle = wk3(2,im,il)*d1 - wk3(3,im,il)   ! dle = log(e)'
                de = dle * e
                vmed(il) = vmed(il) + de
              if (ider == 2) vmedd(il) = vmedd(il) - ! de = e*dle =>
     .          e*wk3(2,im,il)*d1**2 + de*dle          ! de' = e*dle' + e'*dle
            endif
            endif
          enddo
        enddo

C --- Goodwin-Skinner-Pettifor: V (r0/d)^n exp[n(-{d/rc}^nc+{r0/rc}^nc)]
      elseif (memode == 5) then
        call dcopy(ntermx*ilme,tabme(1,ioff),1,wk4,1)
c       call dcopy(5*ilme,tabme(1,ioff),1,wk4,1)
        do il = il1, il2
          n  = wk4(2,il)
          nc = wk4(3,il)
          r0 = wk4(4,il)
          rc = wk4(5,il)
          if (nc > -d1mach(1)) then
            if (dabs(wk4(1,il)) > 1d-12) then
              A  = wk4(1,il)*r0**n*exp(n*(r0/rc)**nc)
              frac = (dist/rc)**nc
              vme(il) = A*exp(-n*frac)/dist**n
              if (ider >= 1) then
                dle = -n*d1*(1d0 + nc*frac)            ! dle = log(vme)'
                vmed(il) = vme(il) * dle
                if (ider == 2) vmedd(il) = vme(il) *
     .            (dle**2 + n*d1**2*(1d0 - nc*(nc-1d0)*frac))
              endif
            else
              vme(il) = 0d0
              if (ider >= 1) then
                vmed(il) = 0d0
                if (ider == 2) vmedd(il) = 0d0
              endif
            endif
          else
            call info0(10,1,0,' vmeder: for Znam et al. augmented'//
     .      ' power law use memode 3 or 4')
            call rxi(' vmeder : negative nc in GSP mode, nc = ',nc)
          endif
        enddo

C --- ME = a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc) ---
      elseif (memode == 7) then
        call dcopy(ntermx*ilme,tabme(1,ioff),1,wk4,1)
c       call dcopy(4*ilme,tabme(1,ioff),1,wk4,1)
c       call rxx(ider /= 0,
c    .    ' vmeder: ME derivatives are not implemented in memode 7')
        do il = il1, il2
          A = dexp(wk4(3,il)*(dist - wk4(4,il)))
          vme(il) = wk4(1,il) * dist**(-wk4(2,il))/(1d0 + A)
          if (ider >= 1) then
            dle = -wk4(2,il)*d1 - wk4(3,il)*A/(1 + A)                  ! dle = [log(vme)]'
            vmed(il) = vme(il) * dle
            if (ider == 2) vmedd(il) = vmed(il)*dle +                ! de = e*dle =>
     .          vme(il)*(wk4(2,il)*d1**2 - wk4(3,il)**2*A/(1 + A)**2)  ! de' = e*dle' + e'*dle
          endif
        enddo

      else
        call rxi(' vmeder: bad memode = ',memode)
      endif

      end
