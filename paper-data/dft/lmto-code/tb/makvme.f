      subroutine makvme(memode,ider,ioff,ilme,ntermx,tabme,decay,
     .   alat,lscale,dist,poly,cutmod,cutme,vmec,vmecd)
C- Makes matrix elements for a single pair of species and applies
C- a smooth cutoff
C ----------------------------------------------------------------------
Ci Inputs
Ci   memode: scaling law for matrix elements
Ci   ider  : switch as to how many spatial derivatives of ME requested
Ci           0 no derivatives, only ME value
Ci           1 ME value + first derivative
Ci   ioff  : offset to tables with species parameters tabme and decay
Ci           (including spins)
Ci   ilme  : = nlme(nl), number of TB MEs
Ci   ntermx: leading dimension of tabme
Ci   tabme : table containing the ME parameters for each pair of species
Ci           the meaning of parameters depends on memode (see rdtbh.f)
Ci   decay : decay of MEs with distance for memode = 2 or 3
Ci   alat  : lattice constant (a.u.)
Ci   lscale: scale (lscale=.true.)/do not scale (lscale=.false.)
Ci           cutoff distances with alat
Ci   dist  : distance at which MEs are sought (in a.u.)
Ci   poly  : order of the cutoff polynomial (n=4 or n=5, see poly45.f)
Ci   cutmod: cutoff mode
Ci           0  no cutoff
Ci           1  augmentative cutoff: vme is augmented with the cutoff polynomial
Ci           2  multiplicative cutoff: vme is multiplied by the cutoff polynomial
Ci           (see rdtbh.f)
Ci   cutme : set of cutoff distances r1 and rc (in units of alat or
Ci           in a.u. depending on lscale)
Co Outputs
Co   vmec  : MEs at distance dist, V_i(dist)
Co   vmecd : V'_i(dist) if ider \= 0, otherwise not referenced
Cr Remarks
Cr   makvme calls vmeder (which actually evaluates ME) and then applies
Cr   the specified cutoff
Cr
Cr   memode = 0 and memode = 6 do not imply any cutoff since these are
Cr   constant ME modes. Program stops if the respective cutmod is not 0.
Cu Updates
Cu   19 Apr 11 (SL) first created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical, intent(in) ::  lscale
      integer, intent(in) :: memode,ider,ntermx,ioff,ilme
      integer, intent(in) :: poly,cutmod
      double precision, intent(in) :: decay(ilme),tabme(ntermx,ilme)
      double precision, intent(in) :: cutme(2,ilme)
      double precision, intent(in)  :: alat,dist
      double precision, intent(out) :: vmec(ilme),vmecd(ilme)
C Local paramters
      integer mode0,il
      double precision r1,r2,e,de,dde,dummy(ilme)
      double precision val(ilme),slo(ilme),curv(ilme)
      double precision cut(2,ilme)

!       call tcn('makvme')

      mode0 = mod(memode,10)
      if (cutmod == 0) then
C ... case no cutoff
        call vmeder(mode0,ider,ioff,ilme,ntermx,tabme,decay,dist,
     .                  1,ilme,vmec,vmecd,dummy)

      else
c ...   stop if memode is a constant ME mode
        if ((mode0 == 0) .or. (mode0 == 6) .or. (mode0 == 7))
     .    call rxi(' makvme: set CUTMOD = 0 for memode = ',memode)

C ... Set cuttof distances
        if (lscale) then
          cut = 0d0
          call daxpy(2*ilme,alat,cutme,1,cut,1) ! cut = alat*cutme
        else
          cut = cutme
        endif

C ...   begin cycle over ilme
        do il = 1, ilme
          r1 = cut(1,il)
          r2 = cut(2,il)

          if (dist <= r1) then
C ...     case cutoff \= 0 but dist < r1  => no cutoff needed
            call vmeder(mode0,ider,ioff,ilme,ntermx,tabme,decay,dist,
     .                  il,il,vmec,vmecd,dummy)

          elseif (dist >= r2) then
C ...     case cutoff \= 0 and dist > r2 => set output to zero
            vmec(il) = 0d0
            if (ider /= 0) vmecd(il) = 0d0
          else
C ...     case cutoff \= 0 and r1 < dist < r2  => select the cutoff mode
            if (cutmod == 1) then
C ...       augmentative cutoff
              call vmeder(mode0,2,ioff,ilme,ntermx,tabme,decay,r1,
     .                    il,il,val,slo,curv)
              call pcut45(poly,dist,r1,r2,val(il),slo(il),curv(il),
     .                    e,de,dde)
              vmec(il) = e
              if (ider /= 0) vmecd(il) = de
            elseif (cutmod == 2) then
C ...       multiplicative cutoff
              call vmeder(mode0,ider,ioff,ilme,ntermx,tabme,decay,dist,
     .                    il,il,val,slo,dummy)
              call pcut45(poly,dist,r1,r2,1d0,0d0,0d0,e,de,dde)
              vmec(il) = val(il) * e
              if (ider /= 0) vmecd(il) = slo(il) * e + val(il) * de
            else
              call
     .        rxi(' makvme: cutoff mode cutmod = %i not implemented',
     .            cutmod)
            endif
          endif
c ...   end of cycle over ilme
        enddo
      endif

!       call tcx('makvme')
      end
