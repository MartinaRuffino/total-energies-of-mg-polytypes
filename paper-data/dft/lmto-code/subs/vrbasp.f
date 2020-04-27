      subroutine vrbasp(auto,z,a,nr,b,nrmt,lmxa,pl,ql,pz,
     .  nsp,v,rofi,itab,rtab,etab,pnuopt,pzopt,tolgv,gtop,eref1glob,eref2glob,eh)
C- Choose a smooth-Hankel basis from the free-atom potential
C ----------------------------------------------------------------------
Ci Inputs
Ci  ein : (negative) parameter for choosing RSMH1
Ci  eout : (negative) parameter for choosing RSMH2
Ci  eh : Hankel energy, must be set using HAM_AUTOBAS_EH
Cio Inputs/Outputs
Co Outputs
Co   rtab : two sets of smooth Hankel smoothing radii for all l upto lmxa
Co   etab : uniform eh
Co   gtop : gmax needed for describing basis (tol tolgv), occupied l
Co   pzopt : is defined elsewhere, printed here
Cl Local variables
Cr Remarks
Ci   simplified (rmt independent) scheme for setting up basis functions
Ci   if ein,eout are specified and negative (mode 1), these are used:
Ci     v(rsm1)=ein, v(rsm2)=eout
Ci   else use V at the muffin-tin (mode 2)
Ci     v(rsm1)=v(rmt)*2, v(rsm2)=v(rmt)/2
Cu Updates
Cu   17 May 2019, rewritten with second simple scheme using V(rmt)
Cu   02 Mar 2018, function copied from freeat.f
Cu   modified version of fabasp(), file freeat.f, see history there
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, parameter :: n0 = 10
      integer lmxa,nr,nrmt,nsp,auto,itab(n0,2)
      double precision rofi(nr),v(nr,nsp)
      double precision pl(n0,nsp),ql(3,n0,nsp),pz(n0),pnuopt(n0),pzopt(n0)
      double precision rtab(n0,2),etab(n0,2)
      double precision a,b,z,tolgv,gtop,eref1glob,eref2glob,eh

C ... Dynamically allocated arrays
      real(8), allocatable :: g(:)

C ... Local parameters
      integer autob,autop,autopz,irep,konfig,l,nn
      double precision eval,pnu,qvl,gam,gmax,veff(nr)
      double precision rlastzero,rlastpeak,rctp1,rctp2,rsm1,rsm2
      double precision vrmt,rmin,eref1,eref2

      if (lmxa .ge. 7) call rx('fabaspj: lmax too large') !this should be checked elsewhere
      allocate(g(2*nr))

      rmin = 0.6d0 !set from v7input?

C ... flags for determining output
      autob  = mod(auto,10)
      autopz = mod(auto/10,10)
      autop  = mod(auto/100,10)

      if (autob .ne. 2  .and. autob .ne. 4) then 
        call rx(" simple basis needs 2-kappa setup: MTO=2 or 4")
      endif

C ... store any existing/specified basis paramaeters and PZ values
      !call dcopy(n0*2,rtab,1,rtab0,1) !restore not implemented
      !call dcopy(n0*2,etab,1,etab0,1)
      call dcopy(n0,pz,1,pzopt,1)

C ... setup basis setup records
      call iinit(itab,n0*2)
      call dpzero(rtab,n0*2)
      call dpzero(etab,n0*2)

C ... effective potential which generates atomic wfns (spin averaged)
      if (nsp .eq. 1) then
        call dpscop(v,veff,nr,1,1,1d0)
      else
        call dpscop(v,veff,nr,1,1,0.5d0)
        call dpsadd(veff,v(1,2),nr,1,1,0.5d0)
      endif

C ... running maximum GMAX
      gtop = 0

      if ( eref1glob < 0.d0 .and. eref2glob < 0.d0 )then
C   ... mode 1: rsmh(l=0) given by ctp at specified energies
        eref1=eref1glob
        eref2=eref2glob
        call info(20,1,0,' Simple basis setup mode 1:'
     .    //' RSMH via c.t.p at specified energies: %;3d,%;3d',eref1,eref2)
      else 
C   ... mode 2: rsmh(l=0) given by v(rmt)
         vrmt = veff(nrmt)- 2*z/rofi(nrmt)
         eref1 = vrmt*2
         eref2 = vrmt/2
         call info(20,1,0,' Simple basis setup mode 2:'
     .     //' v(rmt)=%;3d',vrmt,vrmt)
         call info(20,0,0,' Simple basis setup mode 2:'
     .     //' RSMH via c.t.p at energies: %;3d,%;3d',eref1,eref2)
      endif

      call info(20,0,1,' Setup smoothed-Hankel basis to lmax=%;3i,'
     .  //' EH=%;3d',lmxa+1,eh)

C ... loop over all l, generate basis, pnus, pzs
      do  l = 0, lmxa

        konfig = int((pl(l+1,1)+pl(l+1,nsp))/2)
        nn = konfig-l-1
        qvl = (ql(1,l+1,1)+ql(1,l+1,nsp))/2

C   ... get exact fa wavefunction, eigval, pnu at rmt
        call popta3(0,l,z,nn,nr,nrmt,rofi,veff,a,b,eval,pnu,g)

C   ... get the c.t.p corresponding to eref1 and eref2
        call ppratf(eref1,z,nr,nr,rofi,a,b,veff,g,rlastzero,rlastpeak,rctp1)
        call ppratf(eref2,z,nr,nr,rofi,a,b,veff,g,rlastzero,rlastpeak,rctp2)

C   ... setup valance basis
        rsm1=rctp1/sqrt(dble(l+1))
        rsm2=rctp2/sqrt(dble(l+1))

C   ... assert that the rsm is not too small
        if (rsm1 < rmin) then
           rsm1 = rmin
           call info(20,0,0,' (Warning) RSMH1, l=%;i fixed to rmin=%;3d',l,rmin)
        endif
        if (rsm2 < rmin) then
           rsm2 = rmin
           call info(20,0,0,' (Warning) RSMH2, l=%;i fixed to rmin=%;3d',l,rmin)
        end if

C   ... estimate maximum FT G using HAM_TOL
        gam = min(rsm1,rsm2)**2/4.d0
        gmax = 1.0d0
        do  irep = 1, 10
          gmax = sqrt(-log(tolgv/gmax**l)/gam)
        enddo

C   ... correct charge (in case of possible pzs) and update gmax if q>0
        if (qvl.gt.0) gtop = max(gtop,gmax) 

C   ... save basis params
        itab(l+1,1) = 1
        itab(l+1,2) = 1
        rtab(l+1,1) = rsm1
        rtab(l+1,2) = rsm2
        etab(l+1,1) = eh 
        etab(l+1,2) = eh

C   ... setup initial PNU values
        if (autop .eq. 1) then !setup pnu
          pnuopt(l+1) = pnu
C        else if (autop .eq. 2) then !override manually
C          pnuopt(l+1) = pdef(l+1,1)  ! pdef never assigned
        end if

      enddo !loop over l

      call info2(20,0,0,' Autogenerated RSM1:%n:1,3;3d',lmxa+1,rtab(1,1))
      call info2(20,0,1,' Autogenerated RSM2:%n:1,3;3d',lmxa+1,rtab(1,2))
      call info2(20,0,0,' Autogenerated Pnu: %n:1,3;3d',lmxa+1,pnuopt)
      call info2(20,0,1,' Autogenerated PZ:  %n:1,3;3d',lmxa+1,pzopt)
      deallocate(g)
      end
