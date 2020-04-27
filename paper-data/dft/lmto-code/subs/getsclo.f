      subroutine getsclo(autopz,pqrmx,elocmx,z,a,nr,b,nrmt,lmxa,pl,ql,pz,nsp,v,rofi)
C- Search for semicore states using atomic eigenvalues/charge beyond rmt
C ----------------------------------------------------------------------
Ci Inputs
Ci   autopz:flag controlling how to generate semicore LO
Ci         :1s digit controls precedence
Ci         :0 do not modify PZ
Ci         :1 or 2 Set PZ that satisfy QLOC or ELOC criterion
Ci         :1 If input PZ is nonzero, use it; otherwise choose LO if criterion met
Ci         :2 PZ satisfying criteria are always added as LO
Ci         :10s digit whether occupied states should consider LO anyway
Ci         :0 any l that for whih hi LO is considered
Ci         :1 restrict search to l with unoccupied valence orbitals
Ci            or are l=0 states of an alkali metal
Ci         :2 restrict search to l with unoccupied valence orbitals
Ci  pqrmx : inclusion parameter for semicore LO: minimum charge beyond rmt
Ci  elocmx : inclusion parameter for semicore LO: lower eigenvalue limit
Cio Inputs/Outputs
Co Outputs
Co   pz : (extended) semicore LO are added to pz table
Cl Local variables
Cr Remarks
Cr   separated from basis scheme to allow pz determination earlier in
Cr   free atom solution (to allow atm.ext file to have Qc corresponding
Cr   to inclusion of semicore LO)
Cr   - if pz on entry contains HLLO, these are kept
Cr     and no futher search is done for that l
Cr   - bug: if pz already contains SCLO, these are overwritten
Cu Updates
Cu   03 Apr 2018 separated from basis finder (fabaspj.f)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, parameter :: n0=10
      integer autopz,lmxa,nr,nrmt,nsp
      double precision a,b,z,pqrmx,elocmx
      double precision rofi(nr),v(nr,nsp),pl(n0,nsp),ql(3,n0,nsp),pz(n0)
C ... Dynamically allocated arrays
      real(8), allocatable :: g(:)
C ... Local parameters
      logical ltmp
      integer konfig,l,nn,nnlo
      double precision eval,pnu,pold,pnulo,qvl,q,qrmt,veff(nr)

      if (mod(autopz,10) == 0) return
      allocate(g(2*nr))

C ... effective potential which generates atomic wfns (spin averaged)
      if (nsp == 1) then
        call dpscop(v,veff,nr,1,1,1d0)
      else
        call dpscop(v,veff,nr,1,1,0.5d0)
        call dpsadd(veff,v(1,2),nr,1,1,0.5d0)
      endif

      call info0(21,1,0,' Search for possible local orbitals')
      call info2(41,0,0,' criteria QLOC: q(r>rmt) > %;3g  ELOC: eig < %,5,3d',
     .  pqrmx,elocmx)

C ... loop over all l, generate pz's
      do  l = 0, lmxa

        konfig = int((pl(l+1,1)+pl(l+1,nsp))/2)
        nn = konfig-l-1
        qvl = (ql(1,l+1,1)+ql(1,l+1,nsp))/2

C   ... get exact fa wavefunction, eigval, pnu at rmt
        call popta3(0,l,z,nn,nr,nrmt,rofi,veff,a,b,eval,pnu,g)

C  ...  if HLLO has already been specified
        if (mod(int(pz(l+1)),10) > int(pl(l+1,1))) then
          call info2(21,0,0,' high-lying LO l=%i: specified PZ=%,6;3d',l,pz(l+1))
          cycle ! that's enough for this l
        endif

C   ... setup semicore local orbitals
        if (nn == 0) cycle ! candidates must have a state below this one

C   ... Exclude from search states w/ charge --- cores are usually deep
        if (mod(autopz/10,10) > 0)  then
          ltmp = qvl > 0 .and. pz(l+1) == 0
          ! Exception: s states in column I metals
          if (mod(autopz/10,10) == 1)  then
          ltmp = ltmp .and.
     .      (l>0 .or. (z/=11 .and. z/=19 .and. z/=37 .and. z/=55 .and. z/=87))
          endif
          if (ltmp) cycle
        endif

        pnulo = konfig-1 + 0.5d0
        nnlo = nn-1
        call popta3(0,l,z,nnlo,nr,nrmt,rofi,veff,a,b,eval,pnulo,g)
        call gintsl(g,g,a,nr,rofi,q)
        call gintsl(g,g,a,nrmt,rofi,qrmt)

        pold = pz(l+1)
        pnulo = 10+pnulo
        if ((mod(autopz,10) == 2 .or. pz(l+1) <= 0) .and.
     .      (eval>elocmx .or. q-qrmt>pqrmx)) pz(l+1) = pnulo
        if (q-qrmt>pqrmx) then ! q(r>rmt) check
          call info5(21,0,0,' LO l=%i: qistl=%;3g exceeds QLOC,'//
     .      ' autogen PZ=%;5,3d use PZ=%;5,3d',l,q-qrmt,pnulo,pz(l+1),4)
        endif
        if (eval>elocmx) then ! e_i(SC) check
          call info5(21,0,0,' LO l=%i: e=%;3g exceeds ELOC,'//
     .      ' autogen PZ=%;5,3d use PZ=%;5,3d',l,eval,pnulo,pz(l+1),4)
        endif
        if (eval<=elocmx .and. q-qrmt<=pqrmx) then
          if (pz(l+1) > 0) then
          call info2(21,0,0,' LO l=%i: specified by user : PZ=%;5,3d',
     .        l,pz(l+1))
          else
          call info5(50,0,0,' LO l=%i: e=%;3g>ELOC and qistl=%;3g<QLOC,'//
     .      ' no LO',l,eval,q-qrmt,4,5)
          endif
        endif
      enddo ! loop over l

      end
