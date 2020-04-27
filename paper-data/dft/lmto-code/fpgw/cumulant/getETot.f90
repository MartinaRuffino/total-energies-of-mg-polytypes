PROGRAM getETot
  IMPLICIT NONE
  REAL(8) ETot, E1, E2, E0, mu, mu0, tmp, EFermi, kFermi, rs, kT, norm
  REAL(8), ALLOCATABLE :: omega(:), omega_dos(:), dos(:), dos0(:), Ak(:), nocck(:), Ek(:), Ek1(:), Ek2(:), xk(:),epsk(:)
  REAL(8), PARAMETER :: fa = 1.919158292677512811d0
  INTEGER iomega, nomega, nomega_dos, ik, nk
  LOGICAL, PARAMETER :: NITest = .FALSE.
  REAL(8), EXTERNAL :: ChemicalPotential, fermi

  PRINT*, 'Enter Rs, temperature in hartree:'
  READ*, rs, kT
  kFermi = fa/rs
  EFermi = kFermi**2/2.d0
  PRINT*, 'Enter number of kpoints, number of frequency points'
  READ*, nk, nomega

  ! Allocate arrays.
  PRINT*, 'alllocating arrays.', nk, nomega
  ALLOCATE(omega(nomega), Ak(nomega), nocck(nk), Ek(nk), Ek1(nk), Ek2(nk), xk(nk), epsk(nk))
  PRINT*, 'done allocating arrays', nk, nomega
  IF(NITest) GOTO 20
  ! Open file
  OPEN(UNIT=16, FILE="dos.dat", STATUS="OLD")
  ! Get total number of points (nomega)
  nomega_dos = 0
  DO
     READ(16,*,END=10)
     nomega_dos = nomega_dos + 1
     PRINT*, nomega_dos
  END DO
10 CONTINUE
  REWIND(16)
  PRINT*, 'Allocating DOS array.', nomega_dos
  ALLOCATE(omega_dos(nomega_dos), dos(nomega_dos), dos0(nomega_dos))
  PRINT*, 'done allocating DOS array.'
  ! Read from dos.dat
  dos0 = 0.d0
  DO iomega = 1, nomega_dos
     READ(16,*) omega_dos(iomega), tmp, dos(iomega)
     IF(omega_dos(iomega).GE.0.d0) dos0(iomega) = SQRT(2.d0*omega_dos(iomega))*3.d0/kFermi**3
     PRINT*, dos0(iomega)
  END DO
  PRINT*, "Finished reading DOS."
20 CONTINUE
  IF(NITest) THEN
     nomega_dos = 10000
     ALLOCATE(omega(nomega),dos(nomega))
     DO iomega = 1, nomega
        omega_dos(iomega) = iomega*(1.84d0/nomega)
        dos(iomega) = SQRT(2.d0*omega(iomega))*3.d0/kFermi**3
     END DO
  ELSE
     dos = dos*3.d0/kFermi**3
  END IF
  ! Find mu
  mu = ChemicalPotential(kT,EFermi,omega_dos,dos,nomega_dos)
  mu0 = ChemicalPotential(kT,EFermi,omega_dos,dos0,nomega_dos)
  PRINT*, "mu = ", mu, mu0, EFermi
  ! Find total energy
  ETot = 0.d0
  norm = 0.d0

  OPEN(33,FILE='spfcn.dat',STATUS='OLD')
  OPEN(34,FILE='ek.dat',STATUS='OLD')
  OPEN(35,FILE='ETotk.dat',STATUS='REPLACE')
  DO ik = 1, nk
     ! Read epsk.
     READ(34,*) epsk(ik)
     xk(ik) = SQRT(2.d0*epsk(ik))

     ! Read spectral function.
     PRINT*, 'Reading Ak.', ik, nomega
     READ(33,*)
     DO iomega = 1, nomega
        READ(33, *) omega(iomega), tmp, Ak(iomega)
     END DO
     READ(33,*)
     ! Integrate to get nk and E2.
     nocck(ik) = 0.d0
     Ek2(ik) = 0.d0
     DO iomega = 1, nomega-1
        nocck(ik) = nocck(ik) + 0.5d0*(Ak(iomega+1)*fermi(omega(iomega+1),mu,kT) + &
             & Ak(iomega)*fermi(omega(iomega),mu,kT))*(omega(iomega+1) - omega(iomega))
        !nocck(ik) = nocck(ik) + 0.5d0*(Ak(iomega+1) + &
        !     & Ak(iomega))*(omega(iomega+1) - omega(iomega))
        Ek2(ik) = Ek2(ik) + 0.5d0*(Ak(iomega+1)*fermi(omega(iomega+1),mu,kT)*omega(iomega+1) + &
             & Ak(iomega)*fermi(omega(iomega),mu,kT)*omega(iomega))*(omega(iomega+1) - omega(iomega))
     END DO
     Ek1(ik) = epsk(ik)*nocck(ik)
     WRITE(35,'(5f20.10)') xk(ik), Ek1(ik), Ek2(ik), Ek1(ik) + Ek2(ik), nocck(ik)
  END DO

  ! Now integrate to get total energy.
  ETot = 0.d0
  E0   = 0.d0
  E1   = 0.d0
  E2   = 0.d0
  DO ik = 1, nk - 1
     E0 = E0 + 0.5d0*(epsk(ik+1)*fermi(epsk(ik+1),mu0,kT)*xk(ik+1)**2 + epsk(ik)*fermi(epsk(ik),mu0,kT)*xk(ik)**2)*(xk(ik+1)-xk(ik))
     E1 = E1 + 0.5d0*(Ek1(ik+1)*xk(ik+1)**2 + Ek1(ik)*xk(ik)**2)*(xk(ik+1)-xk(ik))
     E2 = E2 + 0.5d0*(Ek2(ik+1)*xk(ik+1)**2 + Ek2(ik)*xk(ik)**2)*(xk(ik+1)-xk(ik))
  END DO
  E0 = E0*3.d0/kFermi**3
  E1 = E1*3.d0/kFermi**3
  E2 = E2*3.d0/kFermi**3

  OPEN(36,FILE='etot.dat',STATUS='REPLACE')
  WRITE(36,'(5f20.10)') E1, E2, E1+E2, E0
END PROGRAM getETot

REAL(8) FUNCTION fermi(omega,mu,kT)
  REAL(8) omega, mu, kT

  IF(omega.LT.mu) THEN
     fermi = 1.d0/(EXP((omega-mu)/kT) + 1.d0)
  ELSE
     fermi = EXP((mu-omega)/kT)
     fermi = fermi/(1.d0 + EXP((mu-omega)/kT))
  END IF
  RETURN
END FUNCTION fermi

REAL(8) FUNCTION NElectron(omega,dos,nomega,mu,kT)
  REAL(8) omega(nomega), dos(nomega), mu, kT
  INTEGER iomega
  REAL(8), EXTERNAL :: fermi

  ! Integrate to get number of electrons.
  NElectron = 0.d0
  DO iomega = 1, nomega-1
     NElectron = NElectron + 0.5d0*(dos(iomega)*fermi(omega(iomega),mu,kT) + dos(iomega+1)*fermi(omega(iomega+1),mu,kT))*(omega(iomega+1) - omega(iomega))
  END DO
END FUNCTION NElectron

REAL(8) FUNCTION ChemicalPotential(kT,EFermi,omega,dos,nomega)
  REAL(8) kT, EFermi, mu, mu0, mu1, omega(nomega), dos(nomega)
  INTEGER iter
  REAL(8), PARAMETER :: tol = 1.d-5
  INTEGER, PARAMETER :: maxiter = 1000
  REAL(8), EXTERNAL :: NElectron

  ! mu = EFermi is upper bound.
  mu0 = EVermi
  mu1 = EFermi - kT
  DO iter = 1, maxiter
     mu = mu0*(NElectron(omega,dos,nomega,mu1,kT) - 1.d0) - mu1*(NElectron(omega,dos,nomega,mu0,kT) - 1.d0)
     mu = mu/(NElectron(omega,dos,nomega,mu1,kT) - NElectron(omega,dos,nomega,mu0,kT))
     mu0 = mu1
     mu1 = mu
     IF(ABS(mu1-mu0).LT.tol) EXIT
  END DO

  ChemicalPotential = mu
END FUNCTION ChemicalPotential

