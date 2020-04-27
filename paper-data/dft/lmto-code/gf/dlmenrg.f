      subroutine dlmenrg(ic1,ic2,eterms,efermi)
C- Make total single-site Landau potentials (Etot-Q*Ef) for DLM classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   ic1,ic2:DLM class range to be processed
Ci   etrms  :array with energy terms for all classes (from madpot and asvsph)
Ci   efermi :Fermi energy
Co Outputs
Co   etrms :Landau potential packed into etrms(22,*)
Cr Remarks
Cr   Total energy is sum of K.E., Hartree energy, XC energy:
Cr      etot = ekin + utot + rhoexc
Cr   The kinetic energy is computed via double-counting terms
Cr     ekin = sev + sumec - rhov
Cb Bugs
Cb   The double-counting correction from vbardc is neglected
Cu Updates
Cu   27 Dec 11 (Belashchenko) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ic1,ic2
      double precision eterms(22,*),efermi
C ... Local parameters
      integer stdo,ipr,lgunit
      double precision sumev,sumec,rhov,ekin,utot,rhoexc,edcvxt,eref
      double precision muq,etot
      integer ic
      integer mpipid,procid,master

      stdo = lgunit(1)
      procid = mpipid(1)
      master = 0
      call getpr(ipr)

C ... Remnants from asetot
c     if (mode == 1) edcvxt = eterms(17)
c     if (mode == 2) edcvxt = eterms(18)
      edcvxt = 0d0
      eref = 0d0

      do  ic = ic1, ic2
        sumev = eterms(16,ic)
        sumec = eterms(8,ic)
        rhov = eterms(10,ic) + eterms(11,ic) + edcvxt
        ekin = sumev + sumec - rhov

C ...   Madelung term from slot 20 (see madpot) added to utot
        utot = eterms(3,ic) + eterms(20,ic)
        rhoexc = eterms(6,ic)
        muq = efermi * eterms(13,ic)
        etot = ekin + utot + rhoexc

        eterms(2,ic) = etot - muq

C ...   Printout
        if (ipr >= 45 .and. procid == master) then
          write(stdo,310) ic, sumev
          write(stdo,311) rhov,ekin,eref,rhoexc,utot,etot
          write(stdo,312) muq,etot-muq
        endif
      enddo

  310 format(/1x,' Omega potential for DLM class',i4,': sumev=',f10.6)
  311 format(' rhov=',  f17.6,'   ekin=',f17.6,'   eref=',f16.6
     .        /' rhoep= ',f15.6,'   utot=',f17.6,3x,'etot=',     f16.6)
  312 format(' mu*q=',  f17.6,'  Omega=',f17.6)

      end
