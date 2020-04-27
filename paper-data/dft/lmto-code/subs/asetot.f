      subroutine asetot(mode,s_ham,sev,etot)
C- Make ASA Harris energy or Kohn-Sham energy
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read: eterms seref
Co     Stored:    eterms
Ci         sham%eterms: various integrals for the total energy are used:
Ci         :(3)  utot   = total electrostatic energy
Ci         :(6)  rhoexc = rho * exc
Ci         :(8)  sumec  = sum-of-core eigenvalues
Ci         :(10) xcore  = rhoc * total potential
Ci         :(11) valvef = rhov * total potential
Ci         :(17) rinvxt = double counting terms, input rho * appl field
Ci         :(18) rouvxt = double counting terms, output rho * appl field
Ci         :(20) Bext.M = d.c. from external field.
Ci         :(21) Beff.M = d.c. from effective field
Ci         :              (Fixed spin moment method)
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   mode  :1 Harris Foulkes energy
Ci         :2 Kohn-Sham energy
Ci   sev   :sum of eigenvalues
Ci   eref  :reference energy to be subtracted from total
Co Outputs
Co   etot  :Harris energy (mode = 1)
Co         :Hohenberg-Kohn-Sham energy (mode = 2)
Co   sham->eterms various integrals for the total energy are stored:
Co         :(1 )  eh    = etot = Harris energy (mode=1)
Co         :(2 )  eks   = etot = Hohenberg-Kohn-Sham energy (mode=2)
Co         :(16)  sumev = sum of eigenvalues
Cl Local variables
Cl         :
Cr Remarks
Cr
Cr   Total energy is sum of K.E., Hartree energy, XC energy:
Cr      etot = ekin + utot + rhoexc
Cr   The kinetic energy is computed via double-counting terms
Cr     ekin = sev + sumec - rhov
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   07 Jan 11 Double counting terms from eterms(17,18,21) consistent
Cu   08 Mar 03 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode
      double precision sev,etot
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
C ... Local parameters
      integer stdo,ipr,lgunit,i,j
      double precision eterms(22),sumec,rhov,ekin,utot,rhoexc,eref,
     .  edcvxt
      character nam1*9,nam2*5

      call sanrg(.true.,mode,1,2,' asetot','mode')
      stdo = lgunit(1)
      call getpr(ipr)
C     ipr = 51

      eterms = s_ham%eterms
      eref = s_ham%seref

      eterms(16) = sev

      if (mode == 1) edcvxt = eterms(17)
      if (mode == 2) edcvxt = eterms(18)
      sumec = eterms(8)
C     sumtc = eterms(9)
C     rhov = vat*(rhoc+rhov) + Vapplied*rho
      rhov = eterms(10) + eterms(11) + edcvxt
      ekin = sev + sumec - rhov

      utot = eterms(3)
      rhoexc = eterms(6)
      etot = ekin + utot + rhoexc
C ... Double counting from bext * local moment
      if (eterms(21) /= -99) then
        etot = etot - eterms(21)
      endif
      etot = etot - eref

      if (mode == 1) then
        eterms(1) = etot
        nam1 = 'Harris'
        nam2 = 'ehar='
      else
        eterms(2) = etot
        nam1 = 'Kohn-Sham'
        nam2 = 'ehks='
      endif

      s_ham%eterms = eterms

C --- Printout ---
      if (ipr > 40) then
        call word(nam1,1,i,j)
        write(stdo,310) nam1(i:j),sev
        write(stdo,311) rhov,ekin,eref,rhoexc,utot,nam2,etot
  310   format(/1x,a,' energy:':'  sumev=',f10.6)
  311   format(' rhov=',  f17.6,'   ekin=',f17.6,'   eref=',f16.6
     .        /' rhoep= ',f15.6,'   utot=',f17.6,3x,a5,     f16.6)
      endif

      end

