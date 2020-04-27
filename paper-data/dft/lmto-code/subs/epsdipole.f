C#define CHECK
      subroutine epsdipole(mode,nbas,s_lat,eps00,zstar,mion,q)
C- Dipole-dipole contribution to the dielectric matrix
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  pos awald
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:nkd nkq awald alat vol qlv dlv
Cio    Passed to:  strxion2
Ci Inputs
Ci   mode  :1s digit
Ci         :0 make Cbar
Ci         :1 make Chat
Ci         :2 make eigenfunctions, eigenvalues of Chat
Ci         :3 make D^mat
Ci         :4 make eigenfunctions, eigenvalues of D^mat and convert to cm^-1
Ci   nbas  :size of basis
Ci   eps   :q=0 epsilon, 3x3 matrix
Ci   q     :wave number
Ci   zstar: Born charges
Ci   mion  : 1/sqrt(ion masses)
Cu Updates
Cu   20 Apr 16 (WRL) Add simple force constant matrix for MgO
Cu   12 Apr 16 (MvS) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas
      double precision eps00(3,3),q(3),zstar(3,3,nbas)
      double precision mion(3,3,nbas)
C ... For structures
!       include 'structures.h'
      type(str_lat):: s_lat
C ... Dynamically allocated arrays
      real(8), allocatable :: wk(:),eval(:), evalsr(:)
      complex(8), allocatable :: cij(:,:,:,:),cij0(:,:,:,:),evc(:,:)
      complex(8), allocatable :: csr(:,:,:,:),csr0(:,:,:,:),dij(:,:,:,:)
C ... Local parameters
      integer ib,jb,mod0,nev,i,j
      double precision pi, flst, knn
      complex(8) csr00(3,3)

      pi=acos(-1d0)
      mod0 = mod(mode,10)
      print *, 'epsdipole:',' pi=', pi, 'mod0=',mod0
      allocate(wk(11*3*nbas),eval(3*nbas),evc(3*nbas,3*nbas))
      allocate(cij(3,3,nbas,nbas),dij(3,nbas,3,nbas))

C     For Rocksalt MgO only:  nearest neighbor force constant in Ry/a0^2
      knn=0.148d0/2d0
      allocate(csr(3,3,nbas,nbas))
      do ib=1,nbas
         do jb=1,nbas
            do i=1,3
               do j=1,3
                  csr(i,j,ib,jb)=(0d0,0d0)
               enddo
            enddo
         enddo
       enddo
c      call dpzero(csr,size(csr))
      do ib=1,nbas
        do i=1,3
           csr(i,i,ib,ib)=2d0*knn
        enddo
      enddo

      do i=1,3
        csr(i,i,1,2)=-2d0*knn*cos(q(i)/2d0)
        csr(i,i,2,1)=-2d0*knn*cos(q(i)/2d0)
      enddo

      do ib=1,nbas
        do jb=1,nbas
          print *, ' ib,jb=',ib,jb, ' short-range force constant matrix'
          do  i = 1, 3
            call info2(2,0,0,'%3;11,6D',dble(csr(:,i,ib,jb)),0)
          enddo
        enddo
      enddo
      do ib=1,nbas
        do jb=1,nbas
          csr00 = matmul(matmul(mion(:,:,ib),csr(:,:,ib,jb)),mion(:,:,jb))
          csr(:,:,ib,jb) = csr00(:,:)
        enddo
      enddo

      do ib=1,nbas
        do jb=1,nbas
          print *, ' ib,jb=',ib,jb, ' SR dyn matrix'
          do  i = 1, 3
            call info2(2,0,0,'%3;11,6D',dble(csr(:,i,ib,jb)),0)
          enddo
          do i=1,3
            call info2(2,0,0,'%3;11,6D',dimag(csr(:,i,ib,jb)),0)
          enddo
        enddo
      enddo

      allocate(csr0(3,nbas,3,nbas))
      do ib = 1, nbas
        do jb = 1, nbas
          csr0(:,ib,:,jb) = csr(:,:,ib,jb)
        enddo
      enddo
c     call zprm('csr0',12,csr0,3*nbas,3*nbas,3*nbas)
      call zhev(3*nbas,csr0,csr0,.false.,.true.,3*nbas,1d99,nev,wk,.false.,-1,eval,evc)
      allocate(evalsr(3*nbas))
c       do i=1,3
c          eval(i)=(dble(csr(i,i,1,1))+dble(csr(i,i,2,2)))/2d0
c     .           -dsqrt(dabs((dble(csr(i,i,1,1))-dble(csr(i,i,2,2)))**2/4d0+dble(csr(i,i,1,2))))
c          eval(i+3)=(dble(csr(i,i,1,1))+dble(csr(i,i,2,2)))/2d0
c     .           +dsqrt(dabs((dble(csr(i,i,1,1))+dble(csr(i,i,2,2)))**2/4d0+dble(csr(i,i,1,2))))
c       enddo

      do i=1,3*nbas
        eval(i)=sign(sqrt(abs(eval(i)))*sqrt(2d0)*1E6*13.6057/6.58211889/2/pi/2.99792458,eval(i))
      enddo
      call info2(30,0,0,' evl(SR part)%16p%n;11,6D',3*nbas,eval)
      do i=1,3*nbas
        evalsr(i)=eval(i)
      enddo

      call strxion(min(mod0,1),nbas,s_lat,eps00,zstar,q,cij)
      call fc2dynmat(10,nbas,mion,cij,cij,dij)
      call zprm('cij',2,dij,3*nbas,3*nbas,3*nbas)
      if (mod0 == 0) goto 999

C     Subtract 2nd term, Eq. 71
      allocate(cij0(3,3,nbas,nbas))
      call strxion(1,nbas,s_lat,eps00,zstar,(/0d0,0d0,0d0/),cij0)
C     call zprm('cij0',12,cij0,3*nbas,3*nbas,3*nbas)
      do  ib = 1, nbas
        do  jb = 1, nbas
          cij(:,:,ib,ib) = cij(:,:,ib,ib) - cij0(:,:,ib,jb)
        enddo
      enddo
      deallocate(cij0)
      call fc2dynmat(10,nbas,mion,cij,cij,dij)
      call zprm('C~DD, Eq. 71',2,dij,3*nbas,3*nbas,3*nbas)

      call zhev(3*nbas,dij,dij,.false.,.true.,3*nbas,1d99,nev,wk,.false.,-1,eval,evc)
      call info2(30,0,0,' evals of C~DD%16p%n;11,6D',3*nbas,eval)
      if (mod(mode,10) <= 1) goto 999

C ... Convert C~DD into a dynamical matrix
      call fc2dynmat(1,nbas,mion,cij,cij,cij)
      call fc2dynmat(10,nbas,mion,cij,cij,dij)
      call zhev(3*nbas,dij,dij,.false.,.true.,3*nbas,1d99,nev,wk,.false.,-1,eval,evc)
      call info2(30,0,0,'1000*evals of D~DD%16p%n;11,6D',3*nbas,1000*eval)
      call zprm('dipole-dipole part of D',12,dij,3*nbas,3*nbas,3*nbas)

C     Add dipole force constants (divided by epsinf for MgO) to short range ones.
      do ib=1,nbas
         do jb=1,nbas
            do i=1,3
               do j=1,3
                  cij(i,j,ib,jb)=cij(i,j,ib,jb)/3.0189 + csr(i,j,ib,jb)
               enddo
            enddo
         enddo
      enddo

C --- Secular matrix for D^hat with SR part --
      call fc2dynmat(10,nbas,mion,cij,cij,dij)
      call zhev(3*nbas,dij,dij,.false.,.true.,3*nbas,1d99,nev,wk,.false.,-1,eval,evc)
      call info2(30,0,0,' evl(w/ dipoles)%16p%n;11,6D',3*nbas,eval)
C convert eigenvalues to cm^-1
      do i=1,3*nbas
        eval(i)=sign(sqrt(abs(eval(i)))*sqrt(2d0)*1E6*13.6057/6.58211889/2/pi/2.99792458,eval(i))
      enddo
      call info2(30,0,0,' evl (cm^-1)%16p%n;11,6D',3*nbas,eval)
C calculate LST factor from Optical frequencies with and without dipoles.
      flst=1d0
      do i=4,3*nbas
        flst=flst*(eval(i)/evalsr(i))**2
      enddo
      print *, ' flst=',flst, ' LST2=', (eval(nbas*3)/evalsr(nbas*3-2))**2
      deallocate(wk)

  999 continue
      deallocate(cij)
      call rx0('done')
      end

      subroutine strxion(mode,nbas,s_lat,eps00,zstar,q,cij)
C- Dipole-dipole contribution to force constant matrix for all sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  pos awald
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:nkd nkq awald alat vol qlv dlv
Cio    Passed to:  strxion2
Ci Inputs
Ci   mode  :0 return C^bar, Eq 73 in PRB 55, 10355.  Does not require Z*
Ci         :1 return C^hat, Eq 72 in PRB 55, 10355.  Requires Z*
Ci   nbas  :size of basis
Ci   eps   :q=0 epsilon, 3x3 matrix
Ci   q     :wave number
Ci   zstar :Born effective charges (mode=1 only)
Co Outputs
Co   cij   :C^bar_i,a;j,b(q) or C^hat
Co         :for sites i,j and Cartesian components a,b; see Remarks
Cr Remarks
Cr  This routine calculates the q-dependent dipole-dipole contribution
Cr  C^hat or C^bar to the force constant matrix as described by
Cr  Eqn. 72 (for C^hat) and 73 (for C^bar) in this paper:
Cr    Xavier Gonze and Changyol Lee, Phys. Rev. B55, 10355 (1997)
Cr  Any matrix C_0T has a Fourier transform C~(q), defined in Eq. 65,
Cr  where T is a lattice translation vector.
Cr  Eq. 70 is the dipole-dipole contribution to C in real space.
Cr  As the manuscript notes, it is a generalization of the standard
Cr  dipole-dipole term (Eq. 67) designed to reproduce the proper
Cr  nonanalytic behavior C^NA~(q) as q->0, Eq. 60.
Cr  C~(q), Eq. 71, is the Fourier transform of C.
Cr
Cr  This routine does not calculate C~(q), which depends on Z*,
Cr  but instead C^bar(q), Eq. 73, which does not, or C^hat, Eq. 73,
Cr  C^hat is obtained from C^bar and C~ from C^hat  (Eq. 71) :
Cr  1. C~_i,a;j,b(q) = C^hat_i,a;j,b(q) - delta_ij sum_k C^hat_i,a;k,b(q=0)
Cr     The second term is subtracted to fulfill the acoustic sum rule.
Cr  2. C^hat_i,a;j,b(q) = sum_a'b' Z*_i,a;i_a' Z*_i,b;i_b' [C^bar_i,a;j,b(q)
Cr  3. D^hat_i,a;j,b(q) = C^hat(i,a;j,b)/sqrt(mi,mj)
Cr  Here i,j,k are sites in the cell;  a, b are Cartesian components 1..3
Cu Updates
Cu   11 Apr 16  First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas
      double precision eps00(3,3),q(3),zstar(3,3,nbas)
      complex(8) cij(3,3,nbas,nbas)
C ... For structures
!       include 'structures.h'
      type(str_lat):: s_lat
C ... Local parameters
      real(8), allocatable :: wk(:),yl(:)
      integer ib,jb,i,nrx
      double precision tau(3),cy1,e
      integer, parameter :: IPRT=90/9
      complex(8) :: cij0(3,3),cijh(3,3)
      complex(8) :: hl(9),hlp(9),ghl(3,4)
      procedure(real(8)) :: dlength
      procedure(integer) :: iprint,isw

      call info2(IPRT,1,0,' STRXION: '//
     .  '%?;(n==0);Cbar;;%-1j'//
     .  '%?;(n==1);Force constant (Chat);;%-1j'//
     .  '%j matrices at q=%3:2,5;5d',mode,q)

C     For hsmqe0
C#ifdef CHECK
      nrx = max(s_lat%nkd,s_lat%nkq)
      allocate(wk(nrx*(2*2+10)),yl(nrx*9))
      cy1 = dsqrt(3/(16*datan(1d0)))
C#endif

C --- Generate Force constant matrix C~ ---
      do  ib = 1, nbas
        do  jb = 1, nbas
          do  i = 1, 3
            tau(i) = s_lat%pos(i,jb)-s_lat%pos(i,ib)
          enddo
C         call shortn(tau,tau,s_lat%dlv,s_lat%nkd)
          call strxion2(tau,s_lat%awald,s_lat%alat,s_lat%vol,s_lat%qlv,s_lat%nkq,s_lat%dlv,s_lat%nkd,eps00,q,cij(1,1,ib,jb))

C#ifdef CHECK
          call hsmqe0(2,0d0,10,q,tau,nrx,9,wk,yl,s_lat%awald,s_lat%alat,s_lat%qlv,s_lat%nkq,s_lat%dlv,s_lat%nkd,s_lat%vol,hl)
          e = 0
          if (dlength(3,q,1) == 0d0) then
            e = -1d-8
            call hsmq(1,0,2,e,0d0,10,q,tau,nrx,9,wk,yl,
     .      s_lat%awald,s_lat%alat,s_lat%qlv,s_lat%nkq,s_lat%dlv,s_lat%nkd,s_lat%vol,hl,hlp)
          endif
C         Alternate using rcnsl
C         call rcnsl(e,q,tau,s_lat%awald,2,s_lat%alat,s_lat%qlv,s_lat%nkq,s_lat%dlv,s_lat%nkd,s_lat%vol,cy,hl,hlp)
          call h2grad(0,2,e,1,hl,hl,ghl)
C         Park Y_1 in y,z,x order to cijh in x,y,z order; scale Y_1 -> (x,y,z)/r
          do  i = 3, 1, -1
            cijh(1:3,i) = ghl(1:3,i)/cy1               ! yl(2:3) are y,z
            if (i == 1) cijh(1:3,i) = ghl(1:3,4)/cy1 ! yl(4) is x
          enddo

C         Difference in two Ewald methods
          cijh(:,:) = cijh(:,:) - cij(:,:,ib,jb)
C#endif

          if (mod(mode,10) >= 1) then
            cij0 = matmul(matmul(zstar(:,:,ib),cij(:,:,ib,jb)),zstar(:,:,jb))
            cij(:,:,ib,jb) = cij0(:,:)
          endif

C     ... Printout
          if (iprint() < IPRT) cycle
C#ifdef CHECK
          call info5(IPRT,0,0,' pair  i=%i  j=%i  difference in two Ewald sum methods %;3,3g',ib,jb,dlength(9,cijh,1),0,0)
C#elseC
C          call info2(IPRT,0,0,' pair  i=%i  j=%i',ib,jb)
C#endif
          do  i = 1, 3
            call info2(2,0,0,'%3;11,6D',dble(cij(:,i,ib,jb)),0)
          enddo

        enddo
      enddo

C#ifdef CHECK
      deallocate(wk,yl)
C#endif
      end

      subroutine strxion2(tau,awald,alat,vol,glat,nkg,dlat,nkd,eps,q,cij)
C- Dipole-dipole contribution to force constant matrix for one pair of sites
C ----------------------------------------------------------------
Ci Inputs
Ci   tau   :vector connecting pairs in the unit cell
Ci   awald :Ewald parameter, called Lambda in Phys. Rev. B55, 10355
Ci         :awald is in (atomic units)^-1.
Ci         :awald = 0 => pure direct lattice sum, without Ewald
Ci         :awald < 0 => pure recip lattice sum, without Ewald
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   vol   :cell volume
Ci   glat  :reciprocal lattice vectors
Ci   nkg   :number of glat
Ci   dlat  :direct space lattice vectors
Ci   nkd   :number of dlat
Ci   eps   :q=0 epsilon, 3x3 matrix
Ci   q     :wave number
Co Outputs
Co   cij:  A term required for ionic contribution to force constant matrix,
Co         Eq. 73 in Phys. Rev. B55, 10355
Co         cij has units of the inverse cell volume.
Cr Remarks
Cr  First lattice vector must be the zero vector.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkg,nkd
      double precision tau(3),glat(3,nkg),dlat(3,nkd),q(3),eps(3,3)
      double precision awald,alat,vol
      complex(8) cij(3,3)
C ... Local parameters
      integer ir,ir1,i,j
      double precision pi,tpi,gamma,tpiba2,KepsK,tauK,gfac,expy,erfy,tsrpi,fac1,fac2,dete
      double precision epsinv(3,3),del(3),depsd,K(3),d(3)
      real(8), parameter :: tol=1d-6
      complex(8) gsum(3,3),rsum(3,3),eiphi
      procedure(real(8)) :: ddot,derfc,dlength

      call dinv33(eps,0,epsinv,dete) ! Inverse of epsilon
      pi = 4*datan(1d0)
      tpi = 2*pi
      tsrpi = 2/dsqrt(pi)
      tpiba2 = (tpi/alat)**2

C --- Reciprocal space sum reverse order to accumulate small numbers first ---
      gsum = 0
      if (awald /= 0) then  ! No g vectors if Ewald parameter is 0
      gamma = 0.25d0/awald**2 ! awald is Lambda in PRB55, 10355
      if (dlength(3,tpi/alat*q(1:3),1) > tol) then
        ir1 = 1
      else
        ir1 = 2
C       call rx('need to identify G=0 term at q=0')
      endif
      gfac = 1  ! gfact when there is no Ewald parameter
      do  ir = ir1, nkg, 1
C     do  ir = nkg, ir1, -1
        K(1:3) = tpi/alat * (glat(1:3,ir) + q(1:3))
        KepsK = 0
        do  i = 1, 3
        do  j = 1, 3
          KepsK = KepsK + K(i)*eps(i,j)*K(j)
        enddo
        enddo
        tauK = alat*(K(1)*tau(1) + K(2)*tau(2) + K(3)*tau(3))
        eiphi = cdexp(dcmplx(0d0,tauK))
        if (awald > 0) gfac = dexp(-gamma*KepsK)
        do  i = 1, 3
        do  j = 1, 3
          gsum(i,j) = gsum(i,j) + K(i)*K(j)/KepsK * eiphi * gfac
        enddo
        enddo
      enddo
      gsum = 4*pi/vol*gsum
      endif

C --- Real space sum reverse order to accumulate small numbers first  ---
      rsum = 0
      if (awald >= 0) then    ! No T vectors if Ewald parameter is < 0
      if (dlength(3,alat*tau,1) > tol) then
        ir1 = 1
      else
        ir1 = 2
        rsum = -4d0/(3*dsqrt(pi))*awald**3*epsinv
      endif
      do  ir = nkd, ir1, -1
        d(1:3) = alat*(tau(1:3)-dlat(1:3,ir)) ! Has dimension of length
        del = matmul(epsinv,d)
        depsd = dsqrt(ddot(3,del,1,d,1))
        expy = dexp(-(awald*depsd)**2)  ! awald*depsd is dimensionless
        erfy = derfc((awald*depsd))
C       First factor in Eq 74 is del(i)*del(j)/depsd**2 * fac1
        fac1 = 3*erfy/depsd**3 + 2/dsqrt(pi)*expy*(3*awald/depsd**2 + 2*awald**3) !
C       Second factor in Eq 74 is epsinv(i,j) * fac2
        fac2 = erfy/depsd**3 + 2/dsqrt(pi)*expy*awald/depsd**2            ! Second factor in Eq 74
        eiphi = cdexp(dcmplx(0d0,tpi) * (q(1)*dlat(1,ir)+q(2)*dlat(2,ir)+q(3)*dlat(3,ir)))
        do  i = 1, 3
        do  j = 1, 3
          rsum(i,j) = rsum(i,j) - eiphi*(del(i)*del(j)/depsd**2*fac1 - epsinv(i,j)*fac2)
        enddo
        enddo
      enddo
      endif

      cij = gsum + rsum

      end

      subroutine fc2dynmat(mode,nbas,mion,cij,dij1,dij2)
C- Convert force constant matrix into dynamical matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return the force constant matrix cij in dij, with indices permuted
Ci         :1 Scale cij by (1/sqrt(mi*mj)) to make the dynamical matrix
Ci         :10s digit
Ci         :0 store result into dij1(3,3,nbas,nbas)
Ci         :1 store result into dij2(3,nbas,3,nbas)
Ci   nbas  :size of basis
Ci   mion  :1/sqrt(ion masses)
Ci   cij   :Force constant matrix
Co Outputs
Co   dij   :either cij or dynamical matrix, depending on mode
Co         :If 10s digit mode=0, cij and dij1 can occupy the same address space
Co         :If 10s digit mode=1, dij2 cannot point to the same space as cij
Cu Updates
Cu   27 Apr 16  Split off from strxion
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas
      double precision mion(3,3,nbas)
      complex(8) cij(3,3,nbas,nbas),dij1(3,3,nbas,nbas),dij2(3,nbas,3,nbas)
C ... Local parameters
      integer ib,jb,i
      complex(8) cij0(3,3)

      do  ib = 1, nbas
        do  jb = 1, nbas
          if (mod(mode,10) == 0) then
            cij0(:,:) = cij(:,:,ib,jb)
          else
            cij0 = matmul(matmul(mion(:,:,ib),cij(:,:,ib,jb)),mion(:,:,jb))
          endif
          if (mod(mode/10,10) == 0) then
            dij1(:,:,ib,jb) = cij0(:,:)
          else
            do  i = 1, 3
            dij2(:,ib,i,jb) = cij0(:,i)
            enddo
          endif
        enddo
      enddo

      end
