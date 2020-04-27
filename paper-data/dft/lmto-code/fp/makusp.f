      subroutine makusp(n0,z,nsp,nspc,rmax,lmxa,v,a,nr,rs3,vmtz,pnu,pnz,
     .  rsml,ehl,ul,sl,gz,ruu,rus,rss)
C- Augmentation fcts of pure val,slo (times r) from spherical V and b.c.
C ----------------------------------------------------------------------
Ci Inputs
Ci   n0    :leading dimension of pnu and pnz
Ci   z     :nuclear charge
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rmax  :augmentation radius, in a.u.
Ci   lmxa  :augmentation L-cutoff
Ci   v     :spherical potential (atomsr.f)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rs3   :minimum smoothing radius for extrapolation of MT potential
Ci         :Not used now
Ci   vmtz  :muffin-tin zero: subtracted from V in the fitting procedure.
Ci         :The asymptotic form of V-vmtz is taken to be zero.
Ci         :Not used now
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnz   :boundary conditions for (optional) second p.q.n.
Co Outputs
Co   ul    :r * linear combination of wave functions; see Remarks
Co   sl    :r * linear combination of wave functions; see Remarks
Co   gz    :r * state with 2nd p.q.n; see Remarks
Co         :If no pnz is nonzero, gz is not touched
Co   ruu   :diagonal product ul*ul, including small component
Co         :ruu = ruu(ir,l+1,1,isp)
Co         :If pnz for a given l channel is nonzero, ruu(ir,l+1,2,isp)
Co         :is assigned to   gz*ul
Co   rus   :diagonal product ul*sl, including small component
Co         :rus = rus(ir,l+1,1,isp)
Co         :If pnz for a given l channel is nonzero, rus(ir,l+1,2,isp)
Co         :is assigned to   gz*sl
Co   rss   :diagonal product sl*sl, including small component
Co         :rss = rss(ir,l+1,1,isp)
Co         :If pnz for a given l channel is nonzero, rss(ir,l+1,2,isp)
Co         :is assigned to  gz*gz
Cl Local variables
Cl   lpzi  :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz
Cl         :1 value and slope of gz constructed to be zero at rmax
Cl         :  by admixture of phi,phidot
Cl         :2 a smooth Hankel tail is attached (extended local orbital)
Cl         :  and it is included explicitly in the basis.
Cl         :3 a smooth Hankel tail is attached (extended local orbital)
Cl         :  and is coupled to the valence states in an extended atom
Cl         :  approximation.
Cb Bugs
Cb   lpzi=3 is not implemented
Cr Remarks
Cr   This routine makes linear combinations (u,s) of phi,phidot
Cr   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
Cr   ul and sl are r * u and r * s, respectively.
Cr
Cr   Construction of local orbital when val,slo=0 at rmax (lpzi=1):
Cr   Let phi_z be the w.f. corresponding to loc. orbital spec'd by pnz.
Cr   gz is made for any l for pnz is nonzero where:
Cr      gz = r* ( phi_z - phi_z(rmax) u - (phi_z)'(rmax) s )
Cr   By construction, gz/r has value, slope = 0 at rmax.
Cr
Cu Updates
Cu   19 Dec 13 Makes spin 12 part of ruu,rus,rss.  Altered argument list
Cu   11 Jan 10 Patch phidx call for deep semicore states
Cu   12 Aug 04 First implementation of extended local orbitals
Cu             Altered argument list.
Cu   21 May 04 (ATP) pass n0 as argument rather than local dimension
Cu   06 Mar 02 Added code to scale gz (see SCALEGZ)
Cu   21 Aug 01 Extended to local orbitals.  Altered argument list.
Cu   16 May 00 Adapted from nfp makusp, makusr and potpsr.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: lmxa,nr,nsp,nspc,n0
      real(8), intent(in) :: a,rmax,z,rs3,vmtz,rsml(0:lmxa),ehl(0:lmxa),
     .  v(nr,nsp),pnu(n0,nsp),pnz(10,nsp)
      real(8), intent(inout) :: ul(nr,lmxa+1,nsp),sl(nr,lmxa+1,nsp),gz(2*nr,lmxa+1,nsp),
     .  ruu(nr,lmxa+1,2,*),rus(nr,lmxa+1,2,*),rss(nr,lmxa+1,2,*)
C ... Local parameters
      integer :: i,j,l,k,lpz,lpzi,nrbig,idx,lwronsk
      integer, parameter :: nrx=1501
      real(8) :: dphi,dphip,enu,ez,p,phi,phip,phz,dphz
      real(8) :: rofi(nrx), rwgt(nrx), vbig(nrx*2), g(nrx*2), gp(nrx*8), gzl(nrx*2)
      real(8) :: ux(nr,2),sx(nr,2),rsu(nr),rzu(nr),rzs(nr)
      real(8) :: xi(0:n0),wk(2),fac1,e1,ez1,mu
      procedure(integer) :: nglob
C ... External calls
      external daxpy,dcopy,dpzero,dscal,hansr,isanrg,makus2,
     .         makus3,rxi,vxtrap

C --- Make rofi,rwgt, and possibly extended mesh ---
      call vxtrap(1,z,rmax,lmxa,v,a,nr,nsp,pnz,rs3,vmtz,nrx,
     .  lpz,nrbig,rofi,rwgt,vbig)

      call sanrg(.true.,nrbig,nr,nrx,'makusp:','nrbig')

      lwronsk = nglob('wronsk')
      mu = -0.5d0
C --- Loop over spins and l ---
      do  l = 0, lmxa
      k = l+1
      do  i = 1, nsp
        lpzi = 0
        if (mod(pnz(k,i),100d0) > 0)   lpzi = 1
        if (mod(pnz(k,i),100d0) >= 10) lpzi = 2
        if (mod(pnz(k,i),100d0) >= 20) lpzi = 3
        if (lpzi == 3) call rxi('makusp: not implemented lpzi=',lpzi)
        if (lpzi /= 0) then
          j = 10*lwronsk

C         Case local orbital deeper than valence
          if (int(mod(pnz(k,i),10d0)) < int(pnu(k,1))) j = j + 100

          j = j + 1000
          call makrwf(j,z,rmax,l,v(1,i),a,nr,rofi,pnz(k,i),mu,2,gzl,gp,ez,phz,dphz,phip,dphip,p)
C         Scale extended local orbital
          if (lpzi > 1) then
            call hansr(rsml(l),0,l,1,l,ehl(l),rmax**2,1,1,idx,wk,11,xi)
            fac1 = gzl(nr)/rmax/xi(l)
            call dscal(2*nr,1/fac1,gzl,1)
          endif
        endif

        call makrwf(10*lwronsk,z,rmax,l,v(1,i),a,nr,rofi,pnu(k,i),mu,2,g,gp,enu,phi,dphi,phip,dphip,p)

C        if (l == 1) then
C          call makrwf(10*lwronsk,z,rmax,l,v(1,i),a,nr,rofi,pnu(1,i),-0.5d0,2,g,gp,
C     .      enu,phi,dphi,phip,dphip,p)
C        endif


C   ... Scale gz so that <|gz-P(g,gp)|^2> = 1
C#ifdefC SCALEGZ
C        if (lpzi) then
C          call ortrwf(0,z,l,v(1,i),nr,nr,nr,rofi,rwgt,enu,enu,ez,
C     .      g,gp,gzl,D)
C          call dscal(nr*2,1/D,gzl,1)
C          phz = phz/D
C          dphz = dphz/D
C        endif
C#endif

        call makus2(lpzi,nr,rofi,g,gp,gzl,phi,dphi,phip,dphip,phz,dphz,
     .    l,enu,ez,z,v(1,i),ul(1,k,i),sl(1,k,i),ux(1,i),sx(1,i),
     .    ruu(1,k,1,i),rus(1,k,1,i),rss(1,k,1,i),
     .    ruu(1,k,2,i),rus(1,k,2,i),rss(1,k,2,i))

        if (pnz(k,i) > 0) call dcopy(2*nr,gzl,1,gz(1,k,i),1)

        if (nspc == 2 .and. i == 1) then
          e1 = enu
          ez1 = ez
        elseif (nspc == 2) then

          call makus3(lpzi,nr,rofi,l,e1,ez1,enu,ez,z,v,
     .      ul(1,k,1),sl(1,k,1),ux(1,1),sx(1,1),gz(1,k,1),
     .      ul(1,k,2),sl(1,k,2),ux(1,2),sx(1,2),gz(1,k,2),
     .      ruu(1,k,1,3),rus(1,k,1,3),rss(1,k,1,3),
     .      ruu(1,k,2,3),rus(1,k,2,3),rss(1,k,2,3),
     .      rsu,rzu,rzs)

C         For now, take average of (us,su), (uz,zu), (sz,zs)
          call dscal(nr,.5d0,rus(1,k,1,3),1)
          call daxpy(nr,.5d0,rsu,1,rus(1,k,1,3),1)
C          call prmx('ruu',ruu(1:nr,k,1,1:3),nr,nr,3)
C          call prmx('rus',rus(1:nr,k,1,1:3),nr,nr,3)
C          call prmx('rss',rss(1:nr,k,1,1:3),nr,nr,3)
          if (lpzi /= 0) then
            call dscal(nr,.5d0,ruu(1,k,2,3),1)
            call daxpy(nr,.5d0,rzu,1,ruu(1,k,2,3),1)
            call dscal(nr,.5d0,rus(1,k,2,3),1)
            call daxpy(nr,.5d0,rzs,1,rus(1,k,2,3),1)
          endif
        endif

      enddo
      enddo

C --- If at least one semicore state, zero out missing ones ---
      if (lpz /= 0) then
        do  i = 1, nsp
        do  l = 0, lmxa
          k = l+1
          if (pnz(k,i) == 0) then
            call dpzero(gz(1,k,i),nr)
          endif
        enddo
        enddo
      endif

C     call prrmsh('ul',rofi,ul,nr,nr,1+lmxa)
C      call wrhomt(1,'ul','ul',0,ul,ul,rofi,nr,1+lmxa,nsp)
C      call wrhomt(1,'sl','sl',0,sl,sl,rofi,nr,1+lmxa,nsp)
C      if (lpz /= 0) then
C      call wrhomt(1,'gz','gz',0,gz,gz,rofi,nr,1+lmxa,nsp)
C      endif

      end
