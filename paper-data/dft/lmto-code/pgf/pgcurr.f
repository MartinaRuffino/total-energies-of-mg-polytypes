      subroutine pgcurr(zp,qp,s_ctrl,s_ham,s_pot,s_str,s_site,plat,lpgf,isp,nspc,nspcl,pgplp,ikp,
     .  npl,ld1n,ndg1n,gLR,gRL,gsL,gLL,ldl,ndgl,ofl,gsR,gRR,ldr,ndgr,ofr,jzk)
C- Generate the transmission or reflectance for a given energy and k point
C  from the Green's functions connecting the two layers
C ----------------------------------------------------------------------
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 ldham lncol lham neula offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  plham plhamso
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfnc palp
Cio    Passed to:  plham
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  plham
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  plham plhamso
Ci Inputs
Ci   --- info for particular E, k|| point ---
Ci    zp : Energy
Ci    qp : k point
Ci   ikp : Index for current momentum state
Ci   npl : index for rightmost layer
Ci   --- info for potential and calculating structure matrix ---
Ci s_ham : strux containing info on hamiltonian
Ci s_pot : strux containing info on potential
Ci  plat : crystal lattice vectors
Ci   isp : index for current spin state
Ci  nspc :2 if spin-up and spin-down channels are coupled; else 1.
Ci  nspcl:2 if spin-up and spin-down channels are coupled in the leads; else 1.
Ci  lncol:local copy of s_ctrl%lncol
Ci   --- info for Green functions ---
Ci   lpgf: 5 or 6  calculate transmittance
Ci       : 7       calculate reflectance
Ci   gsR : retarded surface Green function on right side
Ci   gRR : retarded interface Green function on right side
Ci         Used only to calculate reflectance
Ci   ldr : leading dimension of gsR
Ci  ndgr : second dimension of gsR
Ci   ofr : offset for right surface Green function
Ci   gsL : retarded surface Green function on left side
Ci   gLL : retarded interface Green function on left side
Ci         Used only to calculate reflectance
Ci   ldl : leading dimension of gsL
Ci  ndgl : second dimension of gsL
Ci   ofl : offset for left surface Green function
Ci   gLR : advanced GF connecting -1,npl (see Remarks and PRB71, 195422 Eq. 43)
Ci       : Not used now; retarded Green's function gLR(E+) is needed.
Ci       : Routine gets it from the relation Gr_LR = Ga_RL^dagger
Ci   gRL : advanced gf connecting npl,-1 gRL(E-) (see Remarks and PRB71, 195422 Eq. 43)
Ci       : NB: advanced Ga_RL = Gr_LR^dagger
Ci         Used only to calculate transmittance
Ci  ld1n : leading dim of GLR and 2nd dim of GRL
Ci       : It is maximum dimension diagonal GF, layers -1 and npl only
Ci ndg1n : leading dim of GRL and 2nd dim of GLR
Ci       : It is maximum dimension of all diagonal GF in layers -1..npl
Ci pgplp : index and dimensioning info for each PL (pgfset.f)
Ci
Co Outputs
Co   jzk : Current density for this energy zp and k||  (lpgf=5)
Co       : Transmittance (mode 5)
Co       : jzk(1:2,:,:,:,:,1) = Re, Im parts of j
Co       : jzk(:,i,j1,k,l1,1) = b^L_ij1 * g^RL_kj1 * b^R_kl1 * g^RL_l1i
Co       : If leads are collinear, j1=l1=1, and jzk = jzk(:,i,1,k,1,1)
Co       : Reflectance (mode 7)
Co       : jzk(1:2,:,:,:)  = Re, Im parts of j
Co       : jzk(:,isp1,isp2,:) = spin-resolved j (noncollinear case)
Co       : jzk(:,:,:,1) = transmittance (lpgf=5)
Cr Remarks
Cr  D. Stewart 11/12/2001
Cr
Cr   This treatment follows Kudrnovsky et al. Surf Sci 454-456
Cr                                           p. 918-924 (2000)
Cr
Cr   Transmission through the interface region in this case is given by
Cr
Cr   T(E,k||) = trace(B_-1(E,k||)*Gr_-1N(k||)*B_N(E,k||)*Ga_N1(k||)
Cr
Cr   where B_-1 = I*S_-1-2(k||)*[gsL_r - gsL_a]*S_-2-1(k||)
Cr
Cr   We assume that gsL(-2) = gsL(-1)  and  gsR(N+1) = gsR(N)
Cr
Cr        B_N = I*S_N,N+1(k||)*[gsR_r-gsR-a]*S_N+1,N(k||)
Cr
Cr   Reflectance from left side is given by
Cr
Cr    R(E,k||) = trace(B_-1(E,k||)*Gr_-1-1(k||)*B_-1(E,k||)*Ga_-1,-1(k||)
Cr
Cr   Reflectance from left side is given by
Cr
Cr     R(E,k||) = trace(B_N(E,k||)*Gr_NN(k||)*B_N(E,k||)*Ga_NN(k||)
Cr
Cr   Nomenclature
Cr
Cr   Gr_-1N - retarded Green function connecting layer -1 to layer N
Cr   Ga_N-1 - advanced Green function connecting layer N to layer -1
Cr   Retarded and advanced Green functions are related by
Cr     Ga_ij = (Gr_ji)*.  i and j refer to combined spin and orbital indices.
Cr
Cr   gsR_(r,a) - (retarded/advanced) surface Green function for Right
Cr                                   semi-infinite region
Cr   gsL_(r,a) - (retarded/advanced) surface Green function for Left
Cr                                   semi-infinite region
Cr   S_ij  - structure constant connecting layers i and j
Cr
Cr   The current density for a given E and k||, J(E,k||) is
Cr
Cr   J(E,k||) = (e^2/h)*T(E,k||)*(fe-fc)
Cr
Cr   where e is the electronic charge, h is Planck's constant, and
Cr   fe and fc are the Fermi factors for the emitter and collector
Cr   respectively.
Cr
Cr Notes on dimensioning gLR and gRL:
Cr     array       dimensioning              subblock used
Cr      gLR   (ld1n*nspc,ndg1n*nspc,1:2) (ld01*nspc,ld0n*nspc,1:2)
Cr      gRL   (ndg1n*nspc,ld1n*nspc,1:2) (ld0n*nspc,ld01*nspc,1:2)
Cu Updates
Cu   23 Nov 15 (MvS) Option to calculate reflectance.  Modified argument list
Cu   29 May 15 (K Belashchenko) resolve spin flip part of transmission
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Jun 04 (A Chantis) corrections for noncollinear case
Cu   14 Nov 03 (MvS) Extended to noncollinear case
Cu   01 Nov 02 (D. Stewart) bug fix in k-point weights for current
Cu   17 Dec 01 (S Faleev) bug fix in dimensioning
C  ---------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nspc,nspcl,lpgf
      integer ld1n,ndg1n,ldr,ndgr,ldl,ndgl
      integer ikp,npl,ofl,ofr
      integer pgplp(6,-1:10)
      double precision plat(3,3),qp(3)
      double precision zp(2,*)
      double precision jzk(2,2,2,2,2,2)
      double precision gLR(ld1n,nspc,ndg1n,nspc,2)
      double precision gRL(ndg1n,nspc,ld1n,nspc,2)
      double precision gsR(ldr,nspc,ndgr,nspc,2),gRR(ldr,nspc,ndgr,nspc,2)
      double precision gsL(ldl,nspc,ndgl,nspc,2),gLL(ldl,nspc,ndgl,nspc,2)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: s11(:)
      real(8), allocatable :: b1(:,:,:,:,:),bn(:,:,:,:,:)
      complex(8), allocatable :: snn(:)
      complex(8), pointer :: pfa(:,:)
C ... Local parameters
C     logical :: debug = .false.
C     logical :: debug = .true.
      logical lso
      integer lds1,ldsn,ld01,ld0n,i1,i2,j1,j2,stdo,ifi,ifi2,i11,n11,nbasp
      integer ld1nx,ndg1nx,ldlx,ldrx,ndglx,ndgrx,ld01x,lds1x,ld0nx,ldsnx,lhdimp,lidimp
      integer, parameter :: PRTG=40
      procedure(logical) bittst
      procedure(integer) :: fopna,iprint,nglob,rotspd
      character strn*(160)
      character(len=1), parameter :: ch(2) = ['l','r']
C     transmission matrix:
C     In transmssion mode, t11(:,:,:,1) = transmission
C     In reflectance mode, t11(:,:,:,1) = left reflectance
C                          t11(:,:,:,2) = right reflectance
      double precision tr11(nspc,nspcl,nspc,nspcl,2,2)
      double precision t11(nspc,nspc,2,2)
C     double precision t22(ldl,nspc,ldl,nspc,2)
C     double precision t33(ldl,nspc,ldl,nspc,2), jtot
C     Difference between emitter and collector Fermi factors
      double precision fdiff
C     Required for f90 compatibility
      interface
      subroutine pgtrns(idx,ld1n,ndg1n,ld01,ld0n,nspc,nspcl,GLR,GRL,bL,bR,t)
      implicit none
      integer idx,ld1n,ndg1n,ld01,ld0n,nspc,nspcl
      double precision GLR(ld1n,nspc,ndg1n,nspc,2)
      real(8), target :: GRL(ndg1n,nspc,ld1n,nspc,2)
      double precision bL(ld01,nspc,ld01,nspc,2)
      double precision bR(ld0n,nspc,ld0n,nspc,2)
      double precision t(nspc,nspc,2)
      end
      end interface

      ld1nx  = ld1n*nspc
      ndg1nx = ndg1n*nspc
      ldlx   = ldl*nspc
      ldrx   = ldr*nspc
      ndglx  = ndgl*nspc
      ndgrx  = ndgr*nspc
      if (nspcl == 2 .and. nspc /= 2) call rx('pgcurr: nspcl=2 but spins are independent')

      stdo = nglob('stdo')
      if ((lpgf == 5.or.lpgf == 6) .and. isp == 1) then
        ifi = fopna('jzk',-1,0)
      elseif (lpgf == 5.or.lpgf == 6) then
        ifi = fopna('jzk2',-1,0)
      elseif (lpgf == 7 .and. isp == 1) then
        if (nspcl == 2) call rx('pgcurr: nspcl=2 not implemented for reflection, sorry')
        ifi  = fopna('rzkl',-1,0)
        ifi2 = fopna('rzkr',-1,0)
      else
        ifi  = fopna('rzkl2',-1,0)
        ifi2 = fopna('rzkr2',-1,0)
      endif

      if (ikp == 1) then
C        call info2(PRTG,0,0,'        z         iq  %?#(n==7)#  Re R#Re T#%-1j'//
C     .    '%46fIm %?#(n==7)#R#T#%?#(n==2)#  (spin 2)##',lpgf,isp)
        call info5(PRTG,0,0,'        z         iq  %?#(n==7)#  Re R#Re T#'//
     .    '%8f%?#(n==2)#%38f##%-2jIm %?#(n==7)#R#T#%j%?#(n==2)#  (spin 2)##',lpgf,nspc,isp,4,5)
        call awrit2('#     Z%13fiq  Re %?#(n==7)#R#T#%-1j ... '//
     .  '%49fIm %?#(n==7)#R#T#',strn,len(strn),ifi,lpgf,isp)
      endif

C     for now just set fdiff to one (Difference in Fermi functions)
      fdiff = 1d0
C     By doing this, we are effectively calculating the generalized
C     transmission for each energy

C ... Make dimensions for left and right layers
C     ld01 dimensions strx s and self-energy b1 for left layer
C     ld0n dimensions strx s and self-energy bn for right layer
C     Use fact that layers -2,-1,0 all have same dimension
      ld01 = pgplp(4,-1)
      lds1 = 3*ld01
      ld0n = pgplp(4,npl)
      ldsn = 3*ld0n
      ld01x = ld01*nspc
      lds1x = lds1*nspc
      ld0nx = ld0n*nspc
      ldsnx = ldsn*nspc
      lso = bittst(s_ctrl%lncol,4)
      nbasp = s_ctrl%nbasp

C      if (debug) then
C        call yprm('gLR',2,gLR,ld1nx*ndg1nx,ld1nx,ld1nx,ld01x)
C        call yprm('gRL',2,gRL,ld1nx*ndg1nx,ndg1nx,ld01x,ld1nx)
C        call yprm('gsr',2,gsR,ldrx*ndgrx,ldrx,ldrx,ndgrx)
C        call yprm('gsl',2,gsL,ldlx*ndglx,ldlx,ldlx,ndglx)
CCC        if (lpgf == 7) then
CCC          call yprm('grr',2,gRR,ldrx*ndgrx,ldrx,ldrx,ndgrx)
CCC          call yprm('gll',2,gLL,ldlx*ndglx,ldlx,ldlx,ndglx)
CCC        endif
C      endif

C     Assemble s_pot%pfnc.  Similar code in pgsif, pgflu, pgemb
C     Note: potential functions are not used anyway since only inter-PL S are used
      if (lso) then             ! plhamso does rotation internally; not needed here
      elseif (nspc == 2) then
        lidimp = s_ham%ldham(2)
        lhdimp = s_ham%ldham(3)
        call ptr_pot(s_pot,8+1,'pfnc',lhdimp,nspc**2,[0d0])
        pfa => s_pot%palp
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1,1,1d0)
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1+2*lhdimp,1+2*3*lhdimp,1d0)
        i1 = rotspd(0)          ! local -> global z
        call rotspn(230010+100*i1,1,nbasp,nbasp,s_ctrl%nl,s_ham%iprmb,s_ham%eula,s_ham%neula,
     .    [0d0],[0d0],[0d0],lhdimp,lidimp,lidimp,lhdimp,1,[0d0],s_pot%pfnc)
C       call zprm('pfnc',2,s_pot%pfnc,lhdimp,lhdimp,4)
      endif

C --- Calculate the self energy for the left exterior layer ---
C     B1 = I*S_{-1,-2}*(GsL_r(-2)-GsL_a(-2))*S_{-2,-1}
      allocate(s11(ld01x*lds1x))
      if (nspc == 2) then
        allocate(b1(2,ld01x,lds1x,1,1)) ! b1 used as a temporary work array
        call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,-1,npl,ld01,lds1,pgplp,qp,b1,s11)
        deallocate(b1)
      else
        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,-1,pgplp,qp,lds1,ld01,s11)
      endif
C      if (debug) then
C      call yprm('s11',2,s11,lds1x*ld01x,ld01x,ld01x,lds1x)
C      endif

C ... Use fact that layer dimension for PL=-2 is equal to
C     layer dimension of layer -1 since this is the bulk region
      allocate(b1(ld01,nspc,ld01,nspc,2))
C     call pgbfac0(gsL,ldl,ndgl,ofl,nspc,s11,ld01,1,b1)
      call pgbfac(gsL,ldl,ndgl,ofl,nspc,s11,ld01,1,b1)
      deallocate(s11)
C      if (debug) then
C      call yprm('left self energy',2,b1,ld01x**2,ld01x,ld01x,ld01x)
C      endif

C --- Calculate the self energy for the right exterior layer ---
C     BN = I*S_{N-1,N}*(GsR_r(N)-GsR_a(N))*S_{N,N-1}
      allocate(snn(ld0nx*ldsnx))
      if (nspc == 2) then
        allocate(bn(2,ld0nx,ldsnx,1,1)) ! bn is a temporary work array
        call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,npl,npl,ld0n,ldsn,pgplp,qp,bn,snn)
        deallocate(bn)
      else
        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,npl,pgplp,qp,ldsn,ld0n,snn)
      endif
C      if (debug) then
C      call yprm('snn',2,snn,ldsnx*ld0nx,ld0nx,ld0nx,ldsnx)
C      endif

      allocate(bn(ld0n,nspc,ld0n,nspc,2))
C     call pgbfac0(gsR,ldr,ndgr,ofr,nspc,snn,ld0n,2,bn)
      call pgbfac(gsR,ldr,ndgr,ofr,nspc,snn,ld0n,2,bn)
      deallocate(snn)
C      if (debug) then
C        call yprm('right self energy',2,bn,ld0nx**2,ld0nx,ld0nx,ld0nx)
C      endif

      call dpzero(tr11,size(tr11))
      if (lpgf == 5 .or. lpgf == 6) then
        call pgtrns(1,ld1n,ndg1n,ld01,ld0n,nspc,nspcl,GLR,GRL,b1,bn,tr11)
        n11 = 1
      elseif (lpgf == 7) then
        call pgrefl(ldl,ndgl,ld01,nspc,gLL,b1,t11)
        call pgrefl(ldr,ndgr,ld0n,nspc,gRR,bn,t11(1,1,1,2))
        n11 = 2
      else
        call rxi('pgcurr: no mode for lpgf =',lpgf)
      endif
      deallocate(b1,bn)

C ... Take the trace of the final product and put it in jzk
      do  i11 = 1, n11
        do  i1 = 1, nspc
          do  i2 = 1, nspc
            if (lpgf == 7) then
              jzk(1,i1,1,i2,1,i11) = - t11(i1,i2,1,i11)
              jzk(2,i1,1,i2,1,i11) = - t11(i1,i2,2,i11)
            else
              do  j1 = 1, nspcl
                do  j2 = 1, nspcl
                  jzk(1,i1,j1,i2,j2,i11) = - tr11(i1,j1,i2,j2,1,i11)
                  jzk(2,i1,j1,i2,j2,i11) = - tr11(i1,j1,i2,j2,2,i11)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

C ... Subtract the trace of the tt+ without spin-coupled terms
C      do  j = 1, ldl
C          jzk(1,1,2) = jzk(1,1,2) - t11(j,1,j,1,1) + t22(j,1,j,1,1)
C          jzk(2,1,2) = jzk(2,1,2) - t11(j,1,j,1,2) + t22(j,1,j,1,2)
C          jzk(1,2,1) = jzk(1,2,1) - t11(j,2,j,2,1) + t22(j,2,j,2,1)
C          jzk(2,2,1) = jzk(2,2,1) - t11(j,2,j,2,2) + t22(j,2,j,2,2)
C        enddo
C ... Subtract the trace of the whole tt+ without spin-coupled terms
C      do  j = 1, ldl
C          jzk(1,1,2) = jzk(1,1,2) - ( t11(j,1,j,1,1) + t11(j,2,j,2,1)
C     .                 - t22(j,1,j,1,1) - t22(j,2,j,2,1))/2d0
C          jzk(2,1,2) = jzk(2,1,2) - t11(j,1,j,1,2) + t22(j,1,j,1,2)
C          jzk(1,2,1) = jzk(1,2,1) - t11(j,2,j,2,1) + t22(j,2,j,2,2)
C          jzk(2,2,1) = jzk(2,2,1) - t11(j,2,j,2,2) + t22(j,2,j,2,2)
C        enddo

C    ... The transmission computed another way
c        jtot=0d0
c        call pgcur3(1,ld1nx,ndg1nx,ld01x,ld0nx,GLR,GRL,b1,bn,
c     . ldlx,t33)
c        do i=1, ldl
c        do j=1, ldl
c        jtot = jtot +
c     .    t33(i,1,j,1,1)*t33(i,1,j,1,1)
c     .   + t33(i,1,j,1,2)*t33(i,1,j,1,2)
c     .   t33(i,1,j,2,1)*t33(i,1,j,2,1)
c     .    + t33(i,1,j,2,2)*t33(i,1,j,2,2)
c        enddo
c        enddo

C ... Take the trace of the final product and put it in jzk
C      call dpzero(jzk,2*nspc**2)
C      do  i1 = 1, nspc
C      do  i2 = 1, nspc
C        do  j = 1, ldl
C          if(i1 == i2) print*, t11(j,i1,j,i2,1), t11(j,i1,j,i2,2)
C          jzk(1,i1,i2) = jzk(1,i1,i2) - t11(j,i1,j,i2,1)
C          jzk(2,i1,i2) = jzk(2,i1,i2) - t11(j,i1,j,i2,2)
C        enddo
C      enddo
C      enddo

C ... Printout
      if (iprint() >= PRTG .and. lpgf /= 7) then
        write(stdo,129) zp(1,1),ikp,
     .    ((((jzk(1,i1,j1,i2,j2,1), i1=1,nspc), j1=1,nspcl), i2=1,nspc), j2=1,nspcl),
     .    ((((jzk(2,i1,j1,i2,j2,1), i1=1,nspc), j1=1,nspcl), i2=1,nspc), j2=1,nspcl)
      elseif (iprint() >= PRTG) then
        do  i11 = 1, 2
          write(stdo,130) zp(1,1),ikp,ch(i11),
     .      ((jzk(1,i1,1,i2,1,i11), i1=1,nspc), i2=1,nspc),
     .      ((jzk(2,i1,1,i2,1,i11), i1=1,nspc), i2=1,nspc)
        enddo
C       call info2(PRTG,0,0,' %;11,4D,%,8i %:-2,4;4e',...)
      endif
C 129 format (1X,F11.4,I8,1P,(4E12.4:2X))
  129 format (1X,F11.4,I8,1P,(4E12.4:2X,4E12.4))
  130 format (1X,F11.4,I8,1x,a1,1P,4E12.4:/22X,4E12.4)
      if (nspcl == 2) then
        write(ifi,28,advance='no') zp(1,1),ikp
        do  j2 = 1, nspcl
          do  i2 = 1, nspcl
            if (i2+j2 > 2) write(ifi,'(22x)',advance='no')
            write(ifi,30)
     .        ((jzk(1,i1,j1,i2,j2,1), i1=1,nspc), j1=1,nspcl),
     .        ((jzk(2,i1,j1,i2,j2,1), i1=1,nspc), j1=1,nspcl)
          enddo
        enddo
      else
      write(ifi,29) zp(1,1),ikp,
     .      ((jzk(1,i1,1,i2,1,1), i1=1,nspc), i2=1,nspc),
     .      ((jzk(2,i1,1,i2,1,1), i1=1,nspc), i2=1,nspc)
      endif
      if (lpgf == 7) then
      write(ifi2,29) zp(1,1),ikp,
     .    ((jzk(1,i1,1,i2,1,2), i1=1,nspc), i2=1,nspc),
     .    ((jzk(2,i1,1,i2,1,2), i1=1,nspc), i2=1,nspc)
      endif
   28 format (1X,F13.8,I8)
   30 format (1P,4E14.6:2X,4E14.6)
   29 format (1X,F13.8,I8,1P,4E14.6:2X,4E14.6)

      end
      subroutine pgtrns(idx,ld1n,ndg1n,ld01,ld0n,nspc,nspcl,GLR,GRL,bL,bR,t)
C- Kernel called by pgcurr to make transmission matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   idx   :1 calculate t = bL * GRL+ * bR * GRL
Ci         :2 calculate t = bR * GLR+ * bL * GLR
Ci   ld1n  : leading dim of GLR and 2nd dim of GRL
Ci         : It is maximum dimension diagonal GF, layers -1 and npl only
Ci   ndg1n : leading dim of GRL and 2nd dim of GLR
Ci         : It is maximum dimension of all diagonal GF in layers -1..npl
Ci   ld01  : Dimension of 0th PL; also dimensions bL
Ci   ld0n  : Dimension of nth PL; also dimensions bR
Ci   gLR   : advancd GF connecting L,R layers
Ci         : gRL = GLR(ld1n,ndg1n)  used when idx=2
Ci   gRL   : off-diagonal retarded GF connecting R,L layers
Ci         : advanced gf connecting npl,-1 gRL(E-) (see Remarks and PRB71, 195422 Eq. 43)
Ci   bL    : self-energy for left layer; see remarks in pgcurr
Ci         : bL = bL(ld01,ld01)
Ci   bR    : self-energy for left layer; see remarks in pgcurr
Ci         : bR = bR(ld0n,ld0n)
Co Outputs
Ci       : NB: advanced Ga_RL = Gr_LR^dagger
Co   t     :(idx=1) bL * GrLR+ * bR * GaRL = bL * (gRL)+ * bR * gRL
Co         :(idx=2)                          bR * (gLR)+ * bL * GLR (not calculated now)
Co         :where bL=I*S_0L*(GsL_r-GsL_a)*S_L0
Co         :      bR=I*S_nR*(GsR_r-GsR_a)*S_Rn
Cr Remarks
Cr
Cr   The transmission is
Cr
Cr   t = Tr{ b^L g^LR(z+)  b^R g^RL(z-) }
Cr       Tr{ b^L g^RL(z-)+ b^R g^RL(z-) }
Cr
Cr   g^LR(z+) and g^RL(z-) are the retarded (advanced) Green's function connecting
Cr   left and right leads.
Cr   g is the advanced Green's function with (ga^LR)+ = (gr^RL). + => hermition conjugate.
Cr
Cr   b are the surface self-energies,
Cr     b^L = I*S_0L*(gs^L_r - gs^L_a)*S_L0
Cr     b^R = I*S_nR*(gs^R_r - gs^R_a)*S_Rn
Cr
Cr   Write out spin and orbital indices (sum over orbitals is implied):
Cr
Cr   t = sum_i,j,k,l b^L_ij(1:ld01,1:ld01) *  g^LR_jk(1:ld0n,1:ld01)[retarded] *
Cr                   b^R_kl(1:ld0n,1:ld0n) *  g^RL_li(1:ld0n,1:ld01)
Cr       sum_i,j,k,l b^L_ij(1:ld01,1:ld01) *  g^RL_kj(1:ld01,1:ld0n) *
Cr                   b^R_kl(1:ld0n,1:ld0n) *  g^RL_li(1:ld0n,1:ld01)
Cr
Cr   Resolve t into spin components.
Cr
Cr      a. Collinear leads (nspcl=1)
Cr         b^L, b^R are spin diagonal; t can be resolved into 2x2 t_ik
Cr         Returned as t(i,1,k,1)
Cr
Cr      b. Noncollinear leads:
Cr         t can be resolved into a 2x2x2x2 matrix, t_ij,kl
Cr         Returned as t(i,j,k,l)
Cr
Cb Bugs
Cb   Dimensioning for gLR and gRL is wrong if L and R leads are differently dimensioned1
Cu Updates
Cu   24 Nov 15 (MvS) Routine renamed from pgcur2
Cu   17 Dec 01 (S Faleev) bug fix in dimensioning
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer idx,ld1n,ndg1n,ld01,ld0n,nspc,nspcl
      complex(8) :: gLR(ld1n*nspc,ndg1n*nspc)
      complex(8) :: gRL(ndg1n*nspc,ld1n*nspc)
      complex(8) :: bL(ld01*nspc,ld01*nspc)
      complex(8) :: bR(ld0n*nspc,ld0n*nspc)
      double precision t(nspc,nspcl,nspc,nspcl,2)
C ... Local parameters
C     logical :: debug = .true.
C     logical :: debug = .false.
      complex(8), allocatable :: wk1(:,:),wk2(:,:),wt(:,:)
C     double precision wk1(ld0n,ld01,2),wk2(ld01,ld01,2),wt(ld01,ld01,2)
C      real(8), target  :: GRL1(ndg1n,ld1n,2)
C      real(8), pointer :: GRL2(:,:,:)
C      double precision bL1(ld01,ld01,2),bR1(ld0n,ld0n,2)
C     double precision GRLp(ndg1n,nspc,ld1n,nspc,2)
C      real(8), pointer :: GRLp(:,:,:,:,:)
      integer sp1,sp2,s1,s2,sq1,sq2,i,j,k,l,offi,offj,offk,offl,ilm
      integer ld01x,ld0nx,ndg1nx,kcplx,j1,l1
      complex(8), parameter :: zone = 1d0, znull = 0d0
      character strn*80

      if (ld01 > ld1n) call rx('pgtrns: dimension mismatch')

      ld01x = ld01*nspc
      ld0nx = ld0n*nspc
      ndg1nx = ndg1n*nspc
      kcplx = 0

C     Convert to standard complex form
      call ztoyy(gRL,ndg1nx,ld01x,ndg1nx,ld01x,kcplx,1)
      call ztoyy(gLR,ld01x,ndg1nx,ld01x,ndg1nx,kcplx,1)
      call ztoyy(BL,ld01x,ld01x,ld01x,ld01x,kcplx,1)
      call ztoyy(BR,ld0nx,ld0nx,ld0nx,ld0nx,kcplx,1)

C      if (debug) then
C        call zprm('gLR',2,gLR,ld01x,ld0nx,ld01x)
C        call zprm('gRL',2,gRL,ndg1nx,ld01x,ld0nx)
C        call zprm('bL',2,BL,ld01x,ld01x,ld01x)
C        call zprm('bR',2,BR,ld0nx,ld0nx,ld0nx)
C      endif

      allocate(wk1(ld0n,ld01),wk2(ld01,ld01),wt(ld01,ld01))

C     Tr bL(i,j) * gLR(j,k) * bR(k,l) * gRL(l,i)
      do  i = 1, nspc
        do  j = 1, nspc
          if  (i /= j .and. nspcl == 1) cycle
          offi = (i-1)*ld01
          offj = (j-1)*ld01

          do  k = 1, nspc
            do  l = 1, nspc
              if (k /= l .and. nspcl == 1) cycle

              offk = (k-1)*ld0n
              offl = (l-1)*ld0n

C              if (debug) then
C                call info5(1,0,0,'set i=%i j=%i k=%i l=%i',i,j,k,l,5)
C              endif

C             wk1(1:ld0n,1:ld01) <- bR(k,l) * gRL(l,i)
              call zgemm('N','N',ld0n,ld01,ld0n,zone,
     .          bR(1+offk,1+offl),ld0nx,gRL(1+offl,1+offi),ndg1nx,
     .          znull,wk1,ld0n)

C              if (debug) then
C                call yprmi('BR%s,(kl=%2i) * gRL(li=%2i)',[k,l],[l,i],3,wk1,0,ld0n,ld0n,ld01)
C              endif

C             wk2(1:ld01,1:ld01) <- g_rLR(j,k) * bR(k,l) * gRL(l,i) = g+RL(k,j) * bR(k,l) * gRL(l,i)
              call zgemm('C','N',ld01,ld01,ld0n,zone,gRL(1+offk,1+offj),ndg1nx,wk1,ld0n,znull,wk2,ld01)

C              if (debug) then
C                call yprmi('gLR%s,(jk=%2i) * BR(kl=%2i) * gRL',[j,k],[k,l],3,wk2,ld01*ld01,ld01,ld01,ld01)
C              endif

C             wt(1:ld01,1:ld01) <- bL(i,j) * gLR(j,k) * bR(k,l) * gRL(l,i)
              call zgemm('N','N',ld01,ld01,ld01,zone,bL(1+offi,1+offj),ld01x,wk2,ld01,znull,wt,ld01)

C              if (debug) then
C                call zprm('BL(ij) * gLR(jk) * BR(kl) * gRL(li)',2,wt,ld01,ld01,ld01)
C              endif

              j1 = j ; if (nspcl == 1) j1 = 1; l1 = l ; if (nspcl == 1) l1 = 1
              do  ilm = 1, ld01
                t(i,j1,k,l1,1) = t(i,j1,k,l1,1) + dble(wt(ilm,ilm))
                t(i,j1,k,l1,2) = t(i,j1,k,l1,2) + dimag(wt(ilm,ilm))
              enddo

C              if (debug) then
C              call info5(1,0,0,' t(i=%i j=%i k=%i l=%i) =%2;12,7D',i,j,k,l,t(i,j1,k,l1,1:2))
C              endif

            enddo
          enddo
        enddo
      enddo

C     call yprm('t',2,t,nspc*nspcl*nspc*nspcl,nspc*nspcl,nspc*nspcl,nspc*nspcl)

      end
      subroutine pgtrnsOLD(idx,ld1n,ndg1n,ld01,ld0n,nspc,nspcl,GLR,GRL,bL,bR,t)
C- Kernel called by pgcurr to make transmission matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   idx   :1 calculate t = bL * GRL+ * bR * GRL
Ci         :2 calculate t = bR * GLR+ * bL * GLR
Ci   ld1n  : leading dim of GLR and 2nd dim of GRL
Ci         : It is maximum dimension diagonal GF, layers -1 and npl only
Ci   ndg1n : leading dim of GRL and 2nd dim of GLR
Ci         : It is maximum dimension of all diagonal GF in layers -1..npl
Ci   ld01  : Dimension of 0th PL; also dimensions bL
Ci   ld0n  : Dimension of nth PL; also dimensions bR
Ci   GLR   : off-diagonal retarded GF connecting L,R layers
Ci         : GLR = GLR(ld1n,ndg1n)  used when idx=2
Ci   GRL   : off-diagonal retarded GF connecting R,L layers
Ci         : GRL = GRL(ndg1n,ld1n)  used when idx=1
Ci   bL    : self-energy for left layer; see remarks in pgcurr
Ci         : bL = bL(ld01,ld01)
Ci   bR    : self-energy for left layer; see remarks in pgcurr
Ci         : bR = bR(ld0n,ld0n)
Co Outputs
Co   t     :(idx=1) bL * GRL+ * bR * GRL
Co         :(idx=2) bR * GLR+ * bL * GLR
Co         :where bL=I*S_0L*(GsL_r-GsL_a)*S_L0
Co         :      bR=I*S_nR*(GsR_r-GsR_a)*S_Rn
Cr Remarks
Cr
Cr   The transmission is
Cr
Cr   t = Tr{ b^L g^LR b^R g^RL }
Cr
Cr   g is the advanced Green's function with [g^LR]+ = g^RL.
Cr   b are the surface self-energies,
Cr     b^L = I*S_0L*(gs^L_r - gs^L_a)*S_L0
Cr     b^R = I*S_nR*(gs^R_r - gs^R_a)*S_Rn
Cr
Cr   Write out spin and orbital indices (sum over orbitals is implied):
Cr
Cr   t = sum_h,i,j,k b^L_hi(1:l1,1:l1) * [g^RL_ij(1:ln,1:l1)]+ *
Cr                   b^R_jk(1:ln,1:ln) *  g^RL_kh(1:ln,1:l1)
Cr
Cr   Resolve t into spin components.
Cr
Cr      a. Collinear leads (nspcl=1)
Cr         b^L, b^R are spin diagonal; t can be resolved into 2x2 t_hk
Cr
Cr      b. Noncollinear leads:
Cr         t can be resolved into a 2x2x2x2 matrix, t_hi,jk
Cr
Cu Updates
Cu   24 Nov 15 (MvS) Routine renamed from pgcur2
Cu   17 Dec 01 (S Faleev) bug fix in dimensioning
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer idx,ld1n,ndg1n,ld01,ld0n,nspc,nspcl
      double precision GLR(ld1n,nspc,ndg1n,nspc,2)
      real(8), target :: GRL(ndg1n,nspc,ld1n,nspc,2)
      double precision bL(ld01,nspc,ld01,nspc,2)
      double precision bR(ld0n,nspc,ld0n,nspc,2)
      double precision t(nspc,nspcl,nspc,nspcl,2)

C ... Local parameters
      logical :: debug = .true.
      double precision wk1(ld0n,ld01,2),wk2(ld01,ld01,2),wt(ld01,ld01,2)
      real(8), target  :: GRL1(ndg1n,ld1n,2)
      real(8), pointer :: GRL2(:,:,:)
      double precision bL1(ld01,ld01,2),bR1(ld0n,ld0n,2)
C     double precision GRLp(ndg1n,nspc,ld1n,nspc,2)
      real(8), pointer :: GRLp(:,:,:,:,:)
      integer l1,l2,sp1,sp2,s1,s2,sq1,sq2

      if (ld01 > ld1n) call rx('pgtrns: wk1 dimensioned badly')

      if (debug) then
        call yprm('GRL',2,GRL,ndg1n*ld1n*nspc**2,ndg1n*nspc,ndg1n*nspc,ld1n*nspc)
        call yprm('BL',2,BL,ld01*nspc*ld01*nspc,ld01*nspc,ld01*nspc,ld01*nspc)
        call yprm('BR',2,BR,ld0n*nspc*ld0n*nspc,ld0n*nspc,ld0n*nspc,ld0n*nspc)
      endif

C ... TEST: (no spin2 input channel: zero the spin 2,2 block for bR)
C      do i = ld01/2+1,ld01
C         do j = ld01/2+1,ld01
C            do k = 1, 2
C            bR(i,j,k) = 0d0
C            enddo
C         enddo
C      enddo

C      do sp2=1, 2
C      do sp1=1, 2
C      do l1 = 1, ld1n
C      do l2 = 1, ld1n
C      GRLp(l1,sp1,l2,sp2,1) =0d0
C      GRLp(l1,sp1,l2,sp2,1) =0d0
C      enddo
C      enddo
C      enddo
C      enddo

      if (idx == 1) then

        gRLp => GRL
C   ... Copy to local array in case L, R are differently dimensioned
C       (MvS) what is the point of this?  GRLp and GRL have the same dimensions
        if (nspc == 2) then
          allocate(GRLp(ndg1n,nspc,ld1n,nspc,2))
          do l1 = 1, ld1n
            do l2 = 1, ld1n
              GRLp(l1,1,l2,1,1) = GRL(l1,1,l2,1,1)       ! (1,1) block
              GRLp(l1,1,l2,2,1) = GRL(l1,1,l2,2,1)       ! (1,2) block
              GRLp(l1,2,l2,1,1) = GRL(l1+ld1n,1,l2,1,1)  ! (2,1) block
              GRLp(l1,2,l2,2,1) = GRL(l1+ld1n,1,l2,2,1)  ! (2,2) block

              GRLp(l1,1,l2,1,2) = GRL(l1,1,l2,1,2)       ! Ditto, Im G
              GRLp(l1,1,l2,2,2) = GRL(l1,1,l2,2,2)
              GRLp(l1,2,l2,1,2) = GRL(l1+ld1n,1,l2,1,2)
              GRLp(l1,2,l2,2,2) = GRL(l1+ld1n,1,l2,2,2)
            enddo
          enddo
        endif

C   --- Make bL * gN1_a * bR * gN1_r ---
C       wk1(1:ld0n,1:ld01) <- bR(1:ld0n,1:ld0n) * GN1_r(1:ld0n,1:ld01)
        call dpzero(t,2*nspc**2*nspcl**2)
C   ... Loop over spin blocks to decompose t into (sp1,sp2) blocks
        do sp1 = 1, nspc
          do sp2 = 1, nspc
          do s1 = 1, nspcl
          do s2 = 1, nspcl
            if (nspcl == 1) then
              sq1 = sp1
              sq2 = sp2
            else
              sq1 = s1
              sq2 = s2
            endif
            forall ( l1=1:ld1n, l2=1:ndg1n) GRL1(l2,l1,1) = GRLp(l2,sp2,l1,sp1,1)
            forall ( l1=1:ld1n, l2=1:ndg1n) GRL1(l2,l1,2) = GRLp(l2,sp2,l1,sp1,2)
            if ((sq1 == sp1 .and. sq2 == sp2) .or. nspcl == 1) then
              GRL2 => GRL1
            else
              allocate(GRL2(ndg1n,ld1n,2))
              forall ( l1=1:ld1n, l2=1:ndg1n) GRL2(l2,l1,1) = GRLp(l2,sq2,l1,sq1,1)
              forall ( l1=1:ld1n, l2=1:ndg1n) GRL2(l2,l1,2) = GRLp(l2,sq2,l1,sq1,2)
            endif
            forall ( l1=1:ld01, l2=1:ld01)  bL1(l2,l1,1) = bL(l2,sp1,l1,sq1,1)
            forall ( l1=1:ld01, l2=1:ld01)  bL1(l2,l1,2) = bL(l2,sp1,l1,sq1,2)
            forall ( l1=1:ld0n, l2=1:ld0n)  bR1(l2,l1,1) = bR(l2,sp2,l1,sq2,1)
            forall ( l1=1:ld0n, l2=1:ld0n)  bR1(l2,l1,2) = bR(l2,sp2,l1,sq2,2)

C            do l1 = 1, ld1n
C              do l2 = 1, ndg1n
C                GRL1(l2,l1,1) = GRLp(l2,sp2,l1,sp1,1)
C                GRL1(l2,l1,2) = GRLp(l2,sp2,l1,sp1,2)
C              enddo
C            enddo
C            do l1 = 1,ld01
C              do l2 = 1,ld01
C                bL1(l2,l1,1) = bL(l2,sp1,l1,sp1,1)
C                bL1(l2,l1,2) = bL(l2,sp1,l1,sp1,2)
C              enddo
C            enddo
C            do l1 = 1,ld0n
C              do l2 = 1,ld0n
C                bR1(l2,l1,1) = bR(l2,sp2,l1,sp2,1)
C                bR1(l2,l1,2) = bR(l2,sp2,l1,sp2,2)
C              enddo
C            enddo

C           wk1(1:ld0n,1:ld01) <- BR(1:ld0n,1:ld0n) * GN1_r(1:ld0n,1:ld01)
            call yygemm('N','N',ld0n,ld01,ld0n,1d0,bR1,bR1(1,1,2),ld0n,
     .        GRL1,GRL1(1,1,2),ndg1n,0d0,wk1,wk1(1,1,2),ld0n)

            if (debug) then
              call yprmi('BR * GN1_r, isp=%s,(%2i)',[sp1,sp2],0,2,wk1,ld0n*ld0n,ld0n,ld0n,ld01)
            endif

C           wk2(1:ld01,1:ld01) <- G1N_a * BR * GN1_r = GN1+_r (1:ld01,1:ld0n)* (BR * GN1_r)
            call yygemm('C','N',ld01,ld01,ld0n,1d0,GRL2,GRL2(1,1,2),ndg1n,
     .        wk1,wk1(1,1,2),ld0n,0d0,wk2,wk2(1,1,2),ld01)

C            if (debug) then
C              call yprmi('G1N_a * BR * GN1_r, isp=%s,(%2i)',[sp1,sp2],0,2,wk2,ld01*ld01,ld01,ld01,ld01)
C            endif

C          wt(1:ld01,1:ld01) <- BL(1:ld01,1:ld01) * (G1N_a * BR * GN1_r)(1:ld01,1:ld01)
            call yygemm('N','N',ld01,ld01,ld01,1d0,bL1,bL1(1,1,2),ld01,
     .        wk2,wk2(1,1,2),ld01,0d0,wt,wt(1,1,2),ld01)

C            if (debug) then
C              call yprmi('BL * G1N_a * BR * GN1_r, isp=%s,(%2i)',[sp1,sp2],0,2,wt,ld01*ld01,ld01,ld01,ld01)
C            endif

            do  l1 = 1, ld01
              t(sp1,s1,sp2,s2,1) = t(sp1,s1,sp2,s2,1)+wt(l1,l1,1)
              t(sp1,s1,sp2,s2,2) = t(sp1,s1,sp2,s2,2)+wt(l1,l1,2)
            enddo
            if ((sq1 == sp1 .and. sq2 == sp2) .or. nspcl == 1) then
            else
              deallocate(GRL2)
            endif
          enddo
            enddo
          enddo
        enddo
      else
        call rx('pgtrns: idx not implemented')
      endif

      if (nspc == 2) deallocate(GRLp)

      end
      subroutine pgrefl(ldg,ndg,ld0,nspc,g,b,r)
C- Kernel called by pgcurr to make reflection matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldg   : Leading dimension of g
Ci   ndg   : 2nd dim of g
Ci   ld0   : Dimension of surface PL; also dimensions b
Ci   nspc  : 2 if spin-up and spin-down channels are coupled; else 1.
Ci   g     : Retarded Green's function at end layer: g = g(ld0,nspc,ndg,nspc)
Ci   b     : self-energy for end layer: b = b(ld0,ld0)
Co Outputs
Co   r     : Tr[b * ga * b * g_r]
Co         :where b = I*S_0L*(GsL_r-GsL_a)*S_L0  (left lead) or
Co         :      b = I*S_nR*(GsR_r-GsR_a)*S_Rn  (right lead)
Co
Co         :Spin structure in the noncollinear case
Co         :  r_ij = b_i g_ij^a b_j g^r_ji
Co         :       = b_i g+_ji  b_j g_ji    where g = g^r = (g+)^a
Cr Remarks
Cr
Cu Updates
Cu    5 Dec 15 (MvS) adapted from pgtrns
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldg,ndg,ld0,nspc
      double precision g(ldg,nspc,ndg,nspc,2)
      double precision b(ld0,nspc,ld0,nspc,2)
      double precision r(nspc,nspc,2)
C ... Local parameters
      double precision wk1(ld0,ld0,2),wk2(ld0,ld0,2),wt(ld0,ld0,2)
      integer l1,i,j,ld0x,ldgx,ndgx

      ld0x = ld0*nspc
      ldgx = ldg*nspc
      ndgx = ndg*nspc

C      call yprm('g_r',2,g,ndgx*ldgx,ndgx,ndgx,ldgx)
C      call yprm('b',2,b,ld0x*ld0x,ld0x,ld0x,ld0x)

C ... TEST: (no spin2 input channel: zero the spin 2,2 block for bR)
C      do i = ld0/2+1,ld0
C         do j = ld0/2+1,ld0
C            do k = 1, 2
C            b(i,j,k) = 0d0
C            enddo
C         enddo
C      enddo

C --- Make Tr[b * g_a * b * g_r] ---
C     Trace is over orbital indices only.  See Remarks for spin structure.
C     Loop over spin blocks to decompose r into (i,j) blocks
      call dpzero(r,2*nspc**2)
      do  i = 1, nspc
        do  j = 1, nspc

C         wk1 <- b_j g_ji
          call yygemm('N','N',ld0,ld0,ld0,1d0,b(1,j,1,j,1),b(1,j,1,j,2),ld0x,
     .      g(1,j,1,i,1),g(1,j,1,i,2),ldgx,0d0,wk1,wk1(1,1,2),ld0)

C          print *, i,j
C          call yprm('B * g_r',2,wk1,ld0*ld0,ld0,ld0,ld0)

C         wk2 <- g+_ji  b_j g_ji
          call yygemm('C','N',ld0,ld0,ld0,1d0,g(1,j,1,i,1),g(1,j,1,i,2),ldgx,
     .      wk1,wk1(1,1,2),ld0,0d0,wk2,wk2(1,1,2),ld0)

C          print *, i,j
C          call yprm('g_a * b * g_r',2,wk2,ld0*ld0,ld0,ld0,ld0)

C         wt <- b_i g+_ji  b_j g_ji
          call yygemm('N','N',ld0,ld0,ld0,1d0,b(1,i,1,i,1),b(1,i,1,i,2),ld0x,
     .      wk2,wk2(1,1,2),ld0,0d0,wt,wt(1,1,2),ld0)

C         mch -f9f15.10 b $b22 g $g22 b11 g11 -cc -t b11 g11 -x -x -x out.femgo -- -abs -max:g -px
C         print *, i,j
C         call yprm('b * g_a * b * g_r',2,wt,ld0*ld0,ld0,ld0,ld0)

C         Trace b * g_a * b * g_r
          do  l1 = 1, ld0
            r(i,j,1) = r(i,j,1) + wt(l1,l1,1)
            r(i,j,2) = r(i,j,2) + wt(l1,l1,2)
          enddo

        enddo                 !j
      enddo                   !i

      end
      subroutine pgbfac0(gs,ld,ndg,of,nspc,sii,ld0,idx,bout)
C- Kernel to calculate the dressed surface spectral function
C  bout = I*Sij*(Gs_r-Gs_a)*Sji for the current following Kudrnovsky
C ----------------------------------------------------------------------
Ci Inputs:
Ci    gs: surface Green function Gs_r to dress.  Note Gs_a = Gs_r(+)
Ci    ld: row dimension of gs
Ci   ndg: col dimension of gs
Ci    of: offset for surface Green function
Ci   sii: the structure matrix for the adjacent layer
Ci   example : gs(-1) corresponding sii -> sii(0)
Ci             gs(N)  corresponding sii -> sii(N)
Ci   ld0: row dimension of sii and bout
Ci   idx: index marker 1 - left surface Green func
Ci                     2 - right surface Green func
Cr Remarks
Cr  This routine is written for collinear strux and should
Cr  be superseded by pgbfac
Co Outputs:
Co   bout: dressed surface spectral function
C ----------------------------------------------------------------------
      implicit none
C  Input variables
      integer idx,ld,ndg,of,ld0,nspc
      double precision gs(ld,nspc,ndg,nspc,2)
      double precision sii(ld0,ld0,3,2)
C     double precision sii(ld0,lds,2)
      double precision bout(ld0,nspc,ld0,nspc,2)
C  Local variables
C     double precision dfgs(ld,ndg,2)
C     double precision wk(ld,ld0,2)
C     double precision smnr(ld0,lmnr,2)
      integer i,j,k1,k2,i1,i2

C ... Set up the difference: Gs_r - Gs_a; store in bout as work array
      do  i1 = 1, nspc
      do  i2 = 1, nspc
      do  i = 1, ld0
      do  j = 1, ld0
        bout(i,i1,j,i2,1) = gs(i,i1,j+of,i2,1) - gs(j,i2,i+of,i1,1)
        bout(i,i1,j,i2,2) = gs(i,i1,j+of,i2,2) + gs(j,i2,i+of,i1,2)
      enddo
      enddo
      enddo
      enddo

C     call yprm('gs(r)-gs(a)',2,bout,(ld0*nspc)**2,ld0*nspc,ld0*nspc,ld0*nspc)

      do  i1 = 1, nspc
      do  i2 = 1, nspc

C ... Calculate sii(1) <- (Gs_r(j)-Gs_a(j))*sji
C     left surface case: (Gs_r(-2)-Gs_a(-2))*S(-2,-1)
      k1 = 3
      k2 = 2
C     right surface case: (Gs_r(-2)-Gs_a(-2))*S(-2,-1)
      if (idx /= 1) then
        k1 = 2
        k2 = 3
      endif

C     (Gs_r-Gs_a)*S
      call yygemm('N','N',ld0,ld0,ld0,1d0,
     .  bout(1,i1,1,i2,1),bout(1,i1,1,i2,2),ld0*nspc,
     .  sii(1,1,k1,1),sii(1,1,k1,2),ld0,0d0,
     .  sii,sii(1,1,1,2),ld0)

C ... Calculate bout <- I*sij*(Gs_r(j)-Gs_a(j))*sji
      call yygemm('N','N',ld0,ld0,ld0,1d0,
     .  sii(1,1,k2,1),sii(1,1,k2,2),ld0,
     .  sii,sii(1,1,1,2),ld0,0d0,
     .  bout(1,i1,1,i2,1),bout(1,i1,1,i2,2),ld0*nspc)

      enddo
      enddo


      end
      subroutine pgbfac(gs,ld,ndg,of,nspc,sij,ld0,idx,bout)
C- Kernel to calculate the dressed surface spectral function
C  bout = I*Sij*(Gs_r-Gs_a)*Sji for the current following Kudrnovsky
C ----------------------------------------------------------------------
Ci Inputs:
Ci    gs   : surface Green function Gs_r to dress.  Note Gs_a = Gs_r(+)
Ci    ld   : row dimension of gs
Ci   ndg   : col dimension of gs
Ci    of   : offset for surface Green function
Ci   sij   : 2D Bloch summed strux, containing Si,i+1, Si+1,i  (plhamnc)
Ci   ld0   : row dimension of sij and bout
Ci   idx   : index marker 1 - left surface Green func
Ci                        2 - right surface Green func
Co Outputs:
Co   bout: dressed surface spectral function I*Sij*(Gs_r-Gs_a)*Sji
Co       : where I = sqrt(-1)
Cr Remarks
Cr   Example : gs(-1) corresponding sij -> sij(0)
Cr             gs(N)  corresponding sij -> sij(N)
Cr
Cr   Onsite S is OVERWRITTEN by pgbfac
Cr
Cr   Debugging: store arrays gsll and psmsl (see calling routine), and bL.  Compare:
Cr   mc -f9f12.6 gsll -p -cc -t -- -a dg pmsl -split pms 1,nr+1 1:nc+1:nr -pop \
Cr   pms12 dg pms13 -x -x bL -- -px
Cu Updates
Cu   28 Aug 16 In the noncollinear case, pgfac requires sij to be noncollinear
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer idx,ld,ndg,of,ld0,nspc
      double precision gs(ld,nspc,ndg,nspc,2)
      double precision sij(ld0*nspc,3*ld0*nspc,2)
C     double precision sij(ld0,lds,2)
      double precision bout(ld0,nspc,ld0,nspc,2)
C ... Local parameters
      logical :: debug = .false.
      integer i,j,k1,k2,i1,i2,ldr,ldl,ldlx,ldrx,ld0x

C ... Set up the difference: Gs_r - Gs_a; store in bout as work array
      do  i1 = 1, nspc
      do  i2 = 1, nspc
      do  i = 1, ld0
      do  j = 1, ld0
        bout(i,i1,j,i2,1) = gs(i,i1,j+of,i2,1) - gs(j,i2,i+of,i1,1)
        bout(i,i1,j,i2,2) = gs(i,i1,j+of,i2,2) + gs(j,i2,i+of,i1,2)
      enddo
      enddo
      enddo
      enddo

      if (debug) then
      call yprm('gs(r)-gs(a)',2,bout,(ld0*nspc)**2,ld0*nspc,ld0*nspc,ld0*nspc)
      endif
      ldr  = ld0; ldl  = ld0    ! L- and R- right adjacent to end layers same dim
      ld0x = ld0*nspc; ldlx = ldl*nspc; ldrx = ldr*nspc

C ... Calculate sij(1) <- (Gs_r(j)-Gs_a(j))*sji
C     Left surface case: (Gs_r(-2)-Gs_a(-2))*S(-2,-1)
      k1 = 3
      k2 = 2
      k2 = ldlx+1       !  index to S(-1,-2)
      k1 = ldlx+ld0x+1  !  index to S(-2,-1)
C     Right surface case: (Gs_r(-2)-Gs_a(-2))*S(-2,-1)
      if (idx /= 1) then
        k1 = 2
        k2 = 3
        k1 = ldlx+1             !  index to S(-1,-2)
        k2 = ldlx+ld0x+1        !  index to S(-2,-1)
      endif

C ... S(onsite) <- (Gs_r-Gs_a)*S
      call yygemm('N','N',ld0x,ld0x,ld0x,1d0,
     .  bout,bout(1,1,1,1,2),ld0x,
     .  sij(1,k1,1),sij(1,k1,2),ld0x,0d0,
     .  sij,sij(1,1,2),ld0x)

      if (debug) then
      call yprm('(gs_r-gs_a)*S',2,sij,3*(ld0*nspc)**2,ld0*nspc,ld0*nspc,ld0*nspc)
      endif

C ... bout <- I*sij*(Gs_r(j)-Gs_a(j))*sji
      call yygemm('N','N',ld0x,ld0x,ld0x,1d0,
     .  sij(1,k2,1),sij(1,k2,2),ld0x,
     .  sij,sij(1,1,2),ld0x,0d0,
     .  bout,bout(1,1,1,1,2),ld0x)

      if (debug) then
      call yprm('b',2,bout,(ld0*nspc)**2,ld0*nspc,ld0*nspc,ld0*nspc)
      endif

      end

C      subroutine pgcur22(idx,ld1n,ndg1n,ld01,ld0n,GLR,GRL,bL,bR,ldt,t1)
CC- Kernel called by pgcurr to make transmission matrix
CC ----------------------------------------------------------------------
CCi Inputs
CCi   idx   :1 calculate t = bL * GRL+ * bR * GRL
CCi         :2 calculate t = bR * GLR+ * bL * GLR
CCi   ld1n :leading dimension of gLR, second dim. of gRL
CCi   ndg1n:leading dimension of gRL, second dim. of gLR
CCi   ld01 :Dimension of 0th PL; also dimensions bL
CCi   ld0n :Dimension of nth PL; also dimensions bR
CCi   GLR   :off-diagonal retarded GF connecting L,R layers
CCi         :GLR = GLR(ld1n,ndg1n)
CCi   GRL   :off-diagonal retarded GF connecting R,L layers
CCi         :GRL = GRL(ndg1n,ld1n)
CCi   bL    :self-energy for left layer; see remarks in pgcurr
CCi         :bL = bL(ld01,ld01)
CCi   bR    :self-energy for left layer; see remarks in pgcurr
CCi         :bR = bR(ld0n,ld0n)
CCi   ldt   :leading dimension of t
CCo Outputs
CCo   t     :(idx=1) bL * GRL+ * bR * GRL
CCo         :(idx=2) bR * GLR+ * bL * GLR
CCo         :where bL=I*S_0L*(GsL_r-GsL_a)*S_L0
CCo         :      bR=I*S_nR*(GsR_r-GsR_a)*S_Rn
CCr Remarks
CCr
CCu Updates
CCu   17 Dec 01 (S Faleev) bug fix in dimensioning
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer idx,ld1n,ndg1n,ld01,ld0n,ldt,i,j,k
C      double precision GLR(ld1n,ndg1n,2),bL(ld01,ld01,2),t1(ldt,ldt,2)
C      double precision GRL(ndg1n,ld1n,2),bR(ld0n,ld0n,2)
CC ... Local parameters
C      double precision wk1(ld0n,ld01,2),wk2(ld01,ld01,2)
C
C
C      if (ld01 > ld1n) call rx('pgcurr: wk1 dimensioned badly')
C
CC ... zero the spin 1,2 block for GLR
C      do i = 1, ld1n/2
C         do j =(ndg1n/2+1),ndg1n
C            do k = 1, 2
C            GLR(i,j,k) = 0d0
C            enddo
C         enddo
C      enddo
CC ... zero the spin 2,1 block for GLR
C      do i = (ld1n/2+1), ld1n
C         do j = 1, ndg1n/2
C            do k = 1, 2
C            GLR(i,j,k) = 0d0
C            enddo
C         enddo
C      enddo
CC ... zero the spin 1,2 block for GRL
C      do i = 1, ndg1n/2
C         do j = (ld1n/2+1), ld1n
C            do k = 1, 2
C            GRL(i,j,k) = 0d0
C            enddo
C         enddo
C      enddo
CC ... zero the spin 2,1 block for GRL
C      do i = (ndg1n/2+1), ndg1n
C         do j = 1, ld1n/2
C            do k = 1, 2
C            GRL(i,j,k) = 0d0
C            enddo
C         enddo
C      enddo
C
CC      call yprm('GRL',2,GRL,ldt*ldt,ldt,ld01,ld01)
CC      call yprm('GLR',2,GLR,ldt*ldt,ldt,ld01,ld01)
CC      call yprm('bR',2,bR,ldt*ldt,ldt,ld01,ld01)
CC      call yprm('bL',2,bL,ldt*ldt,ldt,ld01,ld01)
C
C      if (idx == 1) then
C
CC   --- Make bL * gN1_a * BR * gN1_r ---
CC       wk1(1:ld0n,1:ld01) <- BR(1:ld0n,1:ld0n) * GN1_r(1:ld0n,1:ld01)
C        call yygemm('N','N',ld0n,ld01,ld0n,1d0,bR,bR(1,1,2),ld0n,
C     .    GRL,GRL(1,1,2),ndg1n,0d0,wk1,wk1(1,1,2),ld0n)
C
CC   ... wk2(1:ld01,1:ld01) <- G1N_a * BR * GN1_r = GN1+_r (1:ld01,1:ld0n)* (BR * GN1_r)
C        call yygemm('C','N',ld01,ld01,ld0n,1d0,GRL,GRL(1,1,2),
C     .    ndg1n,wk1,wk1(1,1,2),ld0n,0d0,wk2,wk2(1,1,2),ld01)
C
CC   ... t1(1:ld01,1:ld01) <- BL(1:ld01,1:ld01) * (G1N_a * BR * GN1_r)(1:ld01,1:ld01)
C        call yygemm('N','N',ld01,ld01,ld01,1d0,bL,bL(1,1,2),ld01,
C     .    wk2,wk2(1,1,2),ld01,0d0,t1,t1(1,1,2),ldt)
C
CC       call yprm('BL * G1N_a * BR * GN1_r',2,t1,ldt*ldt,ldt,ld01,ld01)
C      else
C        call rx('pgcur22: idx not implemented')
C      endif
C
C      end
C      subroutine pgcur3(idx,ld1n,ndg1n,ld01,ld0n,GLR,GRL,bL,bR,ldt,t2)
CC- Kernel called by pgcurr to make transmission matrix
CC ----------------------------------------------------------------------
CCi Inputs
CCi   idx   :1 calculate t = bL * GRL+ * bR * GRL
CCi         :2 calculate t = bR * GLR+ * bL * GLR
CCi   ld1n :leading dimension of gLR, second dim. of gRL
CCi   ndg1n:leading dimension of gRL, second dim. of gLR
CCi   ld01 :Dimension of 0th PL; also dimensions bL
CCi   ld0n :Dimension of nth PL; also dimensions bR
CCi   GLR   :off-diagonal retarded GF connecting L,R layers
CCi         :GLR = GLR(ld1n,ndg1n)
CCi   GRL   :off-diagonal retarded GF connecting R,L layers
CCi         :GRL = GRL(ndg1n,ld1n)
CCi   bL    :self-energy for left layer; see remarks in pgcurr
CCi         :bL = bL(ld01,ld01)
CCi   bR    :self-energy for left layer; see remarks in pgcurr
CCi         :bR = bR(ld0n,ld0n)
CCi   ldt   :leading dimension of t
CCo Outputs
CCo   t     :(idx=1) bL * GRL+ * bR * GRL
CCo         :(idx=2) bR * GLR+ * bL * GLR
CCo         :where bL=I*S_0L*(GsL_r-GsL_a)*S_L0
CCo         :      bR=I*S_nR*(GsR_r-GsR_a)*S_Rn
CCr Remarks
CCr
CCu Updates
CCu   17 Dec 01 (S Faleev) bug fix in dimensioning
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer idx,ld1n,ndg1n,ld01,ld0n,ldt
C      double precision GLR(ld1n,ndg1n,2),bL(ld01,ld01,2),t2(ldt,ldt,2)
C      double precision GRL(ndg1n,ld1n,2),bR(ld0n,ld0n,2)
CC ... Local parameters
CC     double precision wk1(ld0n,ld01,2),wk2(ld01,ld01,2)
C
C      if (ld01 > ld1n) call rx('pgcurr: wk1 dimensioned badly')
C
CC      call yprm('GRL',2,GRL,ldt*ldt,ldt,ld01,ld01)
CC      call yprm('GLR',2,GLR,ldt*ldt,ldt,ld01,ld01)
CC      call yprm('bR',2,bR,ldt*ldt,ldt,ld01,ld01)
CC      call yprm('bL',2,bL,ldt*ldt,ldt,ld01,ld01)
C
C      if (idx == 1) then
C
CC   --- Make bL * gN1_a * BR * gN1_r ---
CC       wk1(1:ld0n,1:ld01) <- BR(1:ld0n,1:ld0n) * GN1_r(1:ld0n,1:ld01)
C        call yygemm('N','N',ld0n,ld01,ld0n,1d0,bR,bR(1,1,2),ld0n,
C     .    GRL,GRL(1,1,2),ndg1n,0d0,t2,t2(1,1,2),ld0n)
C
CC       call yprm('BL * G1N_a * BR * GN1_r',2,t1,ldt*ldt,ldt,ld01,ld01)
C      else
C        call rx('pgcur22: idx not implemented')
C      endif
C
C      end
