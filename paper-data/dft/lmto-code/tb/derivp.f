      subroutine derivp(nkp,nbas,nclass,nsite,nl,nsp,ldim,nvar,nfit,ip,
     .  ikp,nev,ifit,ipc,ivar,lmx,indxsh,iax,npr,qp,plat,dh,z,eband,
     .  lso,lov,cryf,ocryf,wk,iev,dpar)
C- Make derivatives of eigenvalues wrt TB parameters for one k-point
C-----------------------------------------------------------------------
Ci Inputs
Ci   nkp,nbas,nclass,nl,nsp,ipc,lmx,indxsh,plat
Ci   nsite: total number of neighbors in all clusters
Ci   ldim: dimension of "lower" block of hamiltonian matrix
Ci   nvar: number of parameters to vary out of all parameters
Ci   nfit: number of eigenvalues to fit out of ldim*nkp total bands
Ci   if ip >= 1 print timings
Ci   nev: number of eigenvectors found from diagno
Ci   ifit(1,i),ifit(2,i): range of bands to fit for ith k-point
Ci   ivar(1,i): points to the position in the full list of ith variable
Ci   ivar(2,i): parameter type of ith variable, with types 1 to 6:
Ci     1: Hamiltonian parameter
Ci     2: Hamiltonian crystal field parameter
Ci     3: overlap parameter
Ci     4: overlap crystal field parameter
Ci     5: diagonal Hamiltonian parameter
Ci     6: spin-orbit parameter
Ci   ivar2(1,i): points to position in full list of ith sticky parameter
Ci   ivar2(2,i): points to position in full list of parameter to stick to
Ci   iax: neighbor lists
Ci   npr: see tbham
Ci   ikp,qp: k-point index and k-point
Ci   dh: real space derivatives wrt parameters (Hamil., overlap, etc.)
Ci   z: eigenvectors for this k
Ci   eband: bands for this k
Ci   lso: include spin-orbit interactions
Ci   lov: include overlap matrix (non-orthogonal TB)
Ci   cryf: true if crystal field terms in Hamiltonian
Ci   ocryf: true if overlap crystal field terms
Ci   wk: complex work array of dimension ldim*ldim
Co Outputs
Co   iev: cumulative counter for total number of eigenvalues
Co   dpar: derivatives of eigenvalues wrt TB parameters for this k
C-----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nkp,nbas,nclass,nsite,nl,nsp,ldim,nvar,nfit,ip,ikp,nev,iev
      integer ifit(2,nkp),ipc(nbas),ivar(2,nvar),lmx(nclass),indxsh(1),
     .  iax(10,nsite),npr(2,nbas)
      double precision qp(3),plat(3,3),dh(nl**4*nsite*nsp**2,nvar),
     .  z(ldim,ldim,2),eband(ldim),wk(ldim,ldim,2),dpar(nvar,nfit)
      logical lso,lov,cryf,ocryf
C Local parameters
      integer n,n2,iprint,i,iev0,ib,nspu,isp
!       integer oiwk,odhk,odsk
      double precision dpib,ddot
      integer, allocatable :: iwk(:)
      complex(8), allocatable :: dhk(:), dsk(:)
C     call pshpr(120)

      if (ip >= 1 .and. ikp == 1) call tc('...   Entering derivp')
      n = ldim
      n2 = n*n
      if (nev < ifit(2,ikp)) call
     .  rxi('DERIVP: insufficient number of bands to fit, nev = ',nev)
!       call defi(oiwk, nbas)
      allocate(iwk(nbas))

      nspu = 1
      isp = 1

C --- Loop over variables: get derivatives ---
      do  20  i = 1, nvar
        iev0 = iev

C --- Hamiltonian parameter ---
        if (ivar(2,i) == 1 .or. ivar(2,i) == 5
     .                      .or. (ivar(2,i) == 6 .and. lso)) then
!           call defdc(odhk, n2)
          allocate(dhk(n2))
          call tbloch(.false.,qp,nl,nsp,nspu,isp,nbas,plat,lmx,ipc,
     .      indxsh,nsite,iax,npr,dh(1,i),0,0,.false.,ldim,dhk,
     .      iwk)
          call zmpy(dhk,n,1,n2,z,n,1,n2,wk,n,1,n2,n,nev,n)

C --- Printout ---
          if ((ikp == 1 .and. iprint() > 70) .or. (iprint()
     . >= 110)) then
            write(*,500) ivar(1,i)
  500       format(' Derivatives for Hamiltonian parameter: ',i4)
            call yprm('dH',2,dhk,ldim*ldim,ldim,ldim,ldim)
          endif
!           call rlse(odhk)
          deallocate(dhk)

C --- Overlap parameter ---
        elseif (ivar(2,i) == 3 .and. lov) then
!           call defdc(odsk, n2)
          allocate(dsk(n2))
          call tbloch(.false.,qp,nl,nsp,nspu,isp,nbas,plat,lmx,ipc,
     .      indxsh,nsite,iax,npr,dh(1,i),0,0,.false.,ldim,dsk,
     .      iwk)
          call zmpy(dsk,n,1,n2,z,n,1,n2,wk,n,1,n2,n,nev,n)

C --- Printout ---
          if ((ikp == 1 .and. iprint() > 70) .or. (iprint()
     . >= 110)) then
            write(*,510) ivar(1,i)
  510       format(' Derivatives for overlap parameter: ',i4)
            call yprm('dS',2,dsk,ldim,ldim,ldim,ldim)
          endif
!           call rlse(odsk)
          deallocate(dsk)
        elseif (ivar(2,i) /= 2 .and. ivar(2,i) /= 4) then
          call rx('DERIVP: bad ivar(2)')
        endif

C --- Sum over bands ---
        do  10  ib = ifit(1,ikp), ifit(2,ikp)
          iev0 = iev0 + 1
          call rxx(iev0 > nfit,'DERIVP: iev gt nfit')

C --- Crystal field parameter ---
          if (ivar(2,i) == 2 .and. cryf) then
            dpib = 0d0
            call cffor(nbas,nl,nsp,nsite,ldim,ib,lmx,ipc,iax,npr,indxsh,
     .        dh(1,i),z,iwk,.true.,dpib)
            if (iprint() >= 110) then
              write(*,520) ivar(1,i),ib,dpib
  520         format(' Deriv for crystal field parameter: ',i4,
     .          3x,'Band=',i4,4x,'Deriv=',f10.6)
            endif

C --- Overlap crystal field parameter ---
          elseif (ivar(2,i) == 4 .and. ocryf) then
            dpib = 0d0
            call cffor(nbas,nl,nsp,nsite,ldim,ib,lmx,ipc,iax,npr,indxsh,
     .        dh(1,i),z,iwk,.true.,dpib)
            dpib = -eband(ib)*dpib
            if (iprint() >= 110) then
              write(*,530) ivar(1,i),ib,dpib
  530         format(' Deriv for overlap crystal field parameter: ',i4,
     .          3x,'Band=',i4,4x,'Deriv=',f10.6)
            endif

C --- Hamiltonian or overlap parameter ---
          else
            if (ivar(2,i) == 3 .and. lov) then
              call daxpy(n,-1d0*eband(ib),wk(1,ib,1),1,wk(1,ib,1),1)
              call daxpy(n,-1d0*eband(ib),wk(1,ib,2),1,wk(1,ib,2),1)
            endif
            dpib = ddot(n,z(1,ib,1),1,wk(1,ib,1),1)
     .           + ddot(n,z(1,ib,2),1,wk(1,ib,2),1)
          endif

          dpar(i,iev0) = dpar(i,iev0) + dpib
   10   continue
   20 continue

      iev = iev0
!       call rlse(oiwk)
      deallocate(iwk)
      if (ip >= 1 .and. ikp == 1) call tc('Done making derivs')

      end
