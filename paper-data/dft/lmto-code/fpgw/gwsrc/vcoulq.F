      subroutine vcoulq(q,nbloch,ngc,nbas,lx,lxx,nx,nxx,alat,qlat,vol,ngvecc,strx,
     .  rojp,rojb,sgbb,sgpb,fouvb,nblochpmx,bas,rmax,eee,aa,bb,nr,nrx,rkpr,rkmr,rofi,
     o  vcoul)
C- Coulomb matrix (Yukawa potential) for one q.
C ----------------------------------------------------------------------
Ci Inputs
Ci   q     :calculate vcoulq(q)
Ci   nbloch:number of site-based product functions
Ci   ngc   :number of IPW's
Ci   nbas  :size of basis
Ci   lx    :lx(ib) = l-cutoff for site ib
Ci   lxx   :dimensioning parameter.  Must be at least as large as max(lx).
Ci   nx    :nx(l,ibas) = number of radial function product functions
Ci         :Note that the definition is a bit different from nx in basnfp.
Ci   nxx   :max number of nx among all l and ibas.
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   vol   :cell volume
Ci   ngvecc:q+G vector
Ci   strx  :Structure factors
Ci   rojp  :rho-type onsite integral
Ci         :rojp = <j_aL(r) | P(q+G)_aL > where
Ci         :       |P(q+G)_aL> : the projection of exp(i (q+G) r) to aL channnel.
Ci         :       |j_aL>      : \def r^l/(2l+1)!! Y_L.  The spherical bessel functions near r=0.  Energy-dependence is omitted.
Ci   rojb  :rho-type onsite integral
Ci   sgbb  :sigma-type onsite integral
Ci   sgpb  :sigma-type onsite integral
Ci   fouvb :
Ci   nblochpmx
Ci   bas   :basis vectors, in units of alat
Ci   rmax  :augmentation radius, in a.u.
Ci   eee   :Convergence parameter (bare coulomb -> Yukawa with decay eee)
Ci   aa    :radial mesh coefficients : rofi(ir) = b*(exp(a*(ir-1))-1)
Ci   bb    :radial mesh coefficients :  rofi(ir) = b*(exp(a*(ir-1))-1)
Ci   nr    :nr(ib) number of radial mesh points for site ib
Ci   nrx   :dimensioning parameter.  Must be at least as large as maxval(nr).
Ci   rkpr  :
Ci   rkmr  :
Ci   rofi  :radial mesh points
Co Outputs
Co   vcoul :coulomb interaction, in Hartrees
Co         :To convert to Ry units, e^2=2 => vcoul (Ry) = 2*vcoul (Hartree)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 May 15
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: nbloch,nblochpmx,nbas,lxx,nrx,nxx,lx(nbas),nx(0:lxx,nbas),nr(nbas)
      integer :: ngc,ngvecc(3,ngc)
      real(8) :: alat,qlat(3,3),vol,q(3),eee
      real(8) :: aa(nbas),bb(nbas),rofi(nrx,nbas),bas(3,nbas),rmax(nbas)
      real(8) :: rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas)
      complex(8) :: strx((lxx+1)**2,nbas,(lxx+1)**2,nbas)
C     rho-type onsite integral
      real(8)    :: rojb(nxx,0:lxx,nbas)
      complex(8) :: rojp(ngc,(lxx+1)**2,nbas)
C     sigma-type onsite integral
      real(8)    :: sgbb(nxx,nxx,0:lxx,     nbas)
      complex(8) :: sgpb(ngc,nxx,(lxx+1)**2,nbas)
C     Fourier
      complex(8) :: fouvb(ngc,nxx,(lxx+1)**2,nbas)
C     Matrix elements
      complex(8) :: vcoul(nblochpmx,nblochpmx)
C ... Dynamically allocated arrays
      real(8), allocatable :: yl(:),ajr(:,:,:),a1(:,:,:),wjj(:,:,:),rijj(:,:,:)
      complex(8), allocatable :: rojpv(:,:),rojpvs(:,:)
      complex(8), allocatable :: pjyl(:,:),pjylp(:,:),phase(:,:)
C     Help arrays for making fast Bessel functions
      real(8), allocatable :: xi(:,:),y(:),h(:)
C     real(8), allocatable :: cy(:)
C ... Local parameters
!ML MPI related variables
#ifdef USE_MPI
      integer myID
#endif
!       integer, allocatable :: first(:), last(:) ! , displace(:), recvcount(:)
!     double precision a_coeff, b_coeff
!     complex(8), allocatable :: vcoul_local(:,:)
!MvS OMP related variables
!       integer omppid,idalloc,nthreads,tid
      procedure(integer) :: idalloc

      logical, parameter :: debug=.false.
      integer iverbose,verbose,ib1,lm1,ib2,lm2,ig,l,lm,m,n,i1,i2,nlmbas,ig1,ig2,ir
      integer ibl1,ibl2,l1,l2,m1,m2,n1,n2,ibloch,jobBes
      real(8):: pi,fpi,fpivol,tpiba,xx,wv12,wv21
      real(8) :: qg2(ngc),qg(3),qhat(3)
      real(8):: int1x(nrx),int2x(nrx),wt(nrx),sigx(0:lxx),radsig(0:lxx)
C     real(8):: phi(0:lxx),psi(0:lxx)
      complex(8) :: zfacl,p12,p21,v12
C     complex(8) :: fouvp12,fouvp21,sgpp12
      complex(8), parameter :: img=(0d0,1d0)
C     complex(8) :: rojpstrx((lxx+1)**2,nbas)
C     Table of indices for inverse mapping (ibas,l,n,m) for product basis ibloch
      integer :: ibasbl(nbloch),nbl(nbloch),lbl(nbloch),mbl(nbloch),lmbl(nbloch)
      integer, allocatable :: aidx(:), alm(:)
      complex(8), allocatable :: rojpc(:,:), rojbc(:,:), strxc(:,:), tmp(:,:)
      complex(8), parameter :: zone = 1.0_8, znull = 0.0_8

#ifdef ptest
      integer nblochngc
      complex(8),allocatable :: matp(:),matp2(:)
#endif
      procedure(integer) :: GWversion

C     double precision delwc,swalltime

C     if (verbose()>=40) call cputid2(' start vcoulq',0)
      iverbose = (verbose()/30)*10

C     swalltime = delwc()
      pi     = 4d0*datan(1d0)
      fpi    = 4*pi
      fpivol = 4*pi*vol
      tpiba = 2*pi/alat
      jobBes = mod(GWversion()/1000,10)
      call dpzero(vcoul,2*size(vcoul))

      if(verbose()>=30) call cputid2(' vcoul packing start',0)

      allocate(aidx(nbas+1), alm(nbas))

      aidx(1) = 0
      do  ib1 = 1, nbas
        alm(ib1)  = (lx(ib1)+1)**2
        aidx(ib1+1) = aidx(ib1) + alm(ib1)
      end do

      nlmbas = aidx(nbas+1)

      allocate(rojbc(nbloch,nlmbas), rojpc(ngc,nlmbas), strxc(nlmbas,nlmbas))

      do  ib1 = 1, nbas
        call zcopy(ngc*alm(ib1), rojp(1,1,ib1), 1, rojpc(1,aidx(ib1)+1), 1)
      end do

      do  ib2 = 1, nbas
        do  ib1 = 1, nbas
          strxc(aidx(ib1)+1:aidx(ib1+1), aidx(ib2)+1:aidx(ib2+1))
     &         = strx(1:alm(ib1), ib1, 1:alm(ib2), ib2)
        end do
      end do


C ... Indices for inverse mapping bloch function ibloch to (ibas,l,n,m)
      ibloch = 0
      do  ib2 = 1, nbas
        do  l  = 0, lx(ib2)
          do  n  = 1, nx(l,ib2)
            do  m  = -l, l
              ibloch  = ibloch + 1
              ibasbl(ibloch) = ib2
              nbl   (ibloch) = n
              lbl   (ibloch) = l
              mbl   (ibloch) = m
              lmbl  (ibloch) = l**2 + l+1 + m
            enddo
          enddo
        enddo
      enddo

      rojbc = znull
      ibloch = 0
      do  ib1 = 1, nbas
        do  l  = 0, lx(ib1)
          do  n  = 1, nx(l,ib1)
            do  m  = 1, 2*l+1
              ibloch  = ibloch + 1
              rojbc(ibloch,aidx(ib1)+l*l+m) = rojb(n, l, ib1)
            enddo
          enddo
        enddo
      enddo

      if (ibloch/= nbloch) call rx('vcoulq: product basis mismatch')

      if (verbose()>=40) call cputid2('vcoulq complete packing',0)

C --- <B|v|B> ---
!ML  Parallelise over ibl2
!ML      do ibl2 = 1,nbloch
!ML         do ibl1 = ibl2,nbloch
!ML  With this order, we can split Vcoul with respect to ibl2, which leaves contiguous
!ML  slices of memory.
!ML  vcoul = 0 at thus stage. Each processor can just add to its bit
!ML  and we can collect all with an AllReduce
C#ifdef USE_MPI
C
C      a_coeff = 1.d0
C      b_coeff = 0.d0
C
C      allocate( first(1:q_comm%size) )
C      allocate( last(1:q_comm%size) )
C      allocate( displace(1:q_comm%size) )
C      allocate( recvcount(1:q_comm%size) )
C
C      call thdist1(1,a_coeff,b_coeff,q_comm%size,1,nbloch,first,last)
C
C      displace(:) = nblochpmx*(first(:)-1)
C      recvcount(:) = (last(:) - first(:) + 1) * nblochpmx
C
C
C      myID = q_comm%ID+1
C
C
C      write(0,'("[",I4,"] I am ",I4," and run the slice ",I4,":",I4," nblochpmx = ",I5)')
C     &     rank,myID,first(myID),last(myID), nblochpmx
C
C!ML allocate local slice of vcoul:
C      allocate( vcoul_local(nblochpmx, first(myID):last(myID)))
C      vcoul_local = 0.d0
C
C      do ibl2= first(myID), last(myID)  ! 1, nbloch
C
C        ib2= ibasbl(ibl2)
C        n2   = nbl (ibl2)
C        l2   = lbl (ibl2)
C        m2   = mbl (ibl2)
C        lm2  = lmbl(ibl2)
C
C        do ibl1= ibl2, nbloch
C
C          ib1= ibasbl(ibl1)
C          n1   = nbl (ibl1)
C          l1   = lbl (ibl1)
C          m1   = mbl (ibl1)
C          lm1  = lmbl(ibl1)
C
C          vcoul_local(ibl1,ibl2) =
C     &     rojb(n1, l1, ib1) *strx(lm1,ib1,lm2,ib2)
C     &    *rojb(n2, l2, ib2)
C
C          if(ib1==ib2 .and. lm1==lm2) then
C            vcoul_local(ibl1,ibl2) = vcoul_local(ibl1,ibl2) + sgbb(n1,n2,l1, ib1)
C            ! sigma-type contribution. onsite coulomb
C          endif
C
C        enddo ! ibl1
C
C      enddo ! ibl2
C
C
C!      write(0,'("[",I4,"] trying to send   ",I7," bytes; ",I7)') rank, nblochpmx*(last(myID)-first(myID) +1),
C!     &     size(vcoul_local)
C!      write(0,'("[",I4,"] trying to receive",I7," bytes at ",I7)') rank, recvcount(myID),displace(myID)
C!
C      call MPI_AllGatherV(vcoul_local, recvcount(myID), MPI_DOUBLE_COMPLEX,
C     &     vcoul,recvcount, displace,  MPI_DOUBLE_COMPLEX, q_comm%communicator, ierror)
C
C      deallocate (first, last, displace, recvcount, vcoul_local)
C
C#else

      allocate(tmp(nbloch,nlmbas))
      call zgemm('n','n',nbloch,nlmbas,nlmbas,zone,rojbc,nbloch,strxc,nlmbas,znull,tmp,nbloch)
      call zgemm('n','t',nbloch,nbloch,nlmbas,zone,tmp,nbloch,rojbc,nbloch,znull,vcoul,nblochpmx)
      deallocate(tmp)

      ibl1 = 1
      i1 = 1
      do ib1 = 1, nbas
        do l = 0, lx(ib1)
          do n1 = 1, nx(l,ib1)
            do m1 = 1, 2*l+1
              ibl2 = i1
              do n2 = 1, nx(l,ib1)
                do m2 = 1, 2*l+1
                  if (m1 == m2) vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l,ib1)
                  ibl2 = ibl2 + 1
                end do
              end do
              ibl1 = ibl1 + 1
            end do
          end do
          i1 = ibl2
        end do
      end do

C      print "(' sumcheck vcoulq(bb)', 2f18.6)", sum(vcoul(1:nbloch,1:nbloch))
C      call zprm('vcoul(bb)',vcoul,nblochpmx,nbloch,nbloch)

! C$OMP parallel
!         nthreads = omppid(0)
! C$OMP end parallel
!        allocate(first(0:nthreads-1),last(0:nthreads-1))
!        call thdist1(0,0d0,0d0,nthreads,1,nbloch,first,last)
! !
! C     Comment out: May produce wrong answer, and saves little time
! CC$OMP  parallel private(tid,ib2,n2,l2,m2,lm2,ib1,n1,l1,m1,lm1) shared(vcoul,rojb,strx,sgbb)
!       tid = omppid(1)
! C     print *, 'tid',tid,nthreads,first(tid),last(tid)
!       do  ibl2 = first(tid), last(tid) ! 1, nbloch
!         ib2 = ibasbl(ibl2)
!         n2  = nbl (ibl2)
!         l2  = lbl (ibl2)
!         m2  = mbl (ibl2)
!         lm2 = lmbl(ibl2)
!
!         do  ibl1 = 1, nbloch
!           ib1 = ibasbl(ibl1)
!           n1  = nbl (ibl1)
!           l1  = lbl (ibl1)
!           m1  = mbl (ibl1)
!           lm1 = lmbl(ibl1)
!
! !           vcoul(ibl1,ibl2) = rojb(n1, l1, ib1) * strx(lm1,ib1,lm2,ib2) * rojb(n2, l2, ib2)
!           if(ib1==ib2 .and. lm1==lm2) then
!             write(334,'("ib1, l, n1, m1, ibl1, n2, m2, ibl2:",8(x,5i))') ib1, l1, n1, m1, ibl1, n2, m2, ibl2
!             vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ib1)
!             ! sigma-type contribution. onsite coulomb
!           endif
!         enddo                   ! ibl1
!       enddo                     ! ibl2
!
!
!
!              stop
!
! C$OMP parallel
!       nthreads = omppid(0)
! C$OMP end parallel
!       nthreads = max(nthreads,1)
!
!       allocate(first(0:nthreads-1),last(0:nthreads-1))
!       call thdist1(0,0d0,0d0,nthreads,1,nbloch,first,last)
!
! C     Comment out: May produce wrong answer, and saves little time
! CC$OMP  parallel private(tid,ib2,n2,l2,m2,lm2,ib1,n1,l1,m1,lm1) shared(vcoul,rojb,strx,sgbb)
!       tid = omppid(1)
! C     print *, 'tid',tid,nthreads,first(tid),last(tid)
!       do  ibl2 = first(tid), last(tid) ! 1, nbloch
!         ib2 = ibasbl(ibl2)
!         n2  = nbl (ibl2)
!         l2  = lbl (ibl2)
!         lm2 = lmbl(ibl2)
!
!         do  ibl1 = 1, nbloch
!           ib1 = ibasbl(ibl1)
!           n1  = nbl (ibl1)
!           l1  = lbl (ibl1)
!           lm1 = lmbl(ibl1)
!
!           if(ib1==ib2 .and. lm1==lm2) then
!             vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ib1)
!             ! sigma-type contribution. onsite coulomb
!           endif
!         enddo                   ! ibl1
!       enddo                     ! ibl2
!
CC$OMP   end parallel
      if (verbose()>=40) call cputid2(' complete vcoul(BB)',0)

C#endif

      if (ngc == 0) goto 99

C --- <P_G|v|B> ---
C     Comment out: May produce wrong answer, and saves little time
CC$OMP  parallel private(tid,ib2,n2,l2,m2,lm2,ig1,ib1,lm1) shared(vcoul,rojb,strx,sgpb)


      allocate(tmp(nlmbas,ngc))
      call zgemm('t','c',nlmbas,ngc,nlmbas,zone,strxc,nlmbas,rojpc,ngc,znull,tmp,nlmbas)
      call zgemm('t','t',ngc,nbloch,nlmbas,-zone,tmp,nlmbas,rojbc,nbloch,znull,vcoul(nbloch+1,1),nblochpmx)

      ibl1 = 1
      do ib1 = 1, nbas
        do l = 0, lx(ib1)
          do n = 1, nx(l,ib1)
            do m = 1, 2*l+1
              lm = l*l+m
              call zaxpy(ngc,zone,fouvb(1,n,lm,ib1),1,vcoul(nbloch+1,ibl1),1)
              call zaxpy(ngc,-zone,sgpb(1,n,lm,ib1),1,vcoul(nbloch+1,ibl1),1)
              ibl1 = ibl1 + 1
            end do
          end do
        end do
      end do

C     print "(' sumcheck vcoulq(pb)', 2f18.6)", sum(vcoul(nbloch+1:nbloch+ngc,1:nbloch))
C     call zprm('vcoul(pb)',vcoul(nbloch+1,1),nblochpmx,ngc,nbloch)

!       tid = omppid(1)
!       do ibl2 = first(tid), last(tid) ! 1, nbloch
!         ib2= ibasbl(ibl2)
!         n2   = nbl (ibl2)
!         l2   = lbl (ibl2)
!         m2   = mbl (ibl2)
!         lm2  = lmbl(ibl2)
!         do  ig1 = 1, ngc
!           i1 = nbloch + ig1
! !           vcoul(i1,ibl2) = fouvb(ig1,  n2, lm2, ib2)
!
!           do ib1 = 1, nbas
!             do lm1  = 1, (lx(ib1)+1)**2
!               vcoul(i1,ibl2) = vcoul(i1,ibl2) - dconjg(rojp(ig1,lm1,ib1)) * strx(lm1,ib1,lm2,ib2) * rojb(n2,l2,ib2)
! !               if(ib1==ib2 .and. lm1==lm2) then
! !                 vcoul(i1,ibl2) = vcoul(i1,ibl2) - sgpb(ig1,n2,lm2,ib2)
! !               endif
!             enddo               ! lm1
!           enddo                 ! ib1
!         enddo                   ! ig1
!       enddo                     ! ibl2
CC$OMP   end parallel
      if (verbose()>=40) call cputid2(' complete vcoul(PB)',0)

C --- <P_G|v|P_G> ---
C ... Make the following help arrays
C     qg2(ig)            = |q+G|^2
C     phase(ig,ibas)     = exp( i (q+G) R )
C     pjyl(lm,ig)        = <j_l Y_L | exp i q+G r> = projection of exp(i q+G r) to jl yl on MT
      allocate(pjyl((lxx+1)**2,ngc),phase(ngc,nbas),yl((lxx+1)**2))
!       allocate(pjyl((lxx+1)**2,ngc),pjylp(nlmbas,ngc),phase(ngc,nbas),yl((lxx+1)**2))

C     call sylmnc(cy,lxx)
      do  ig = 1, ngc
        qg(1:3) = tpiba * (q(1:3) + matmul(qlat,ngvecc(1:3,ig)))
        qg2(ig)  = sum(qg(1:3)**2)

C   ... pjyl
        if (qg2(ig) > 0) then
          qhat = qg/sqrt(qg2(ig))
        else
          qhat = 0d0
        endif
C       call sylm(qhat,yl,lxx,xx) ! cy*yl = Y_L(\hat{q+G})
        call ropyln(1,qhat(1),qhat(2),qhat(3),lxx,1,yl,xx) ! This one makes Y_L directly
        lm = 0
        do  l = 0, lxx
          zfacl = fpi*img**l * sqrt(qg2(ig))**l
          do  m = -l, l
            lm = lm+1
            pjyl(lm,ig) = zfacl * yl(lm) ! * cy(lm)
          enddo
        enddo

C   ... phase
        do  ib2 = 1, nbas
          phase(ig,ib2) = exp( img*sum(qg(1:3)*bas(1:3,ib2))*alat  )
        enddo

      enddo                     ! ig

C ... Make rojpv, rhojpvs = rhojp^+ * strx.  Contiguous elements for zgemm.
C     rojpstrx = 0d0

      if(verbose()>=30) call cputid2(' vcoul start GG',0)
      call zgemm('t','t',ngc,ngc,nlmbas,zone,tmp,nlmbas,rojpc,ngc,znull,vcoul(nbloch+1,nbloch+1),nblochpmx)
      if(verbose()>=30) call cputid2(' vcoul end GG end',0)

      deallocate(tmp,rojbc,rojpc,strxc,aidx,alm)

C     Diagonal term in interstitial vcoulq

!$omp parallel do
      do  ig = 1, ngc
        vcoul(nbloch+ig,nbloch+ig) = vcoul(nbloch+ig,nbloch+ig) + fpivol/(qg2(ig) - eee) !eee must be negative
      enddo
!$omp end parallel do

!       deallocate(rojpv,rojpvs)

! C ... Make pjylp = pjylp*phase(ib).  Contiguous elements for zgemm.
!       i2 = 0
!       do  ib2 = 1, nbas
!         do  lm1 = 1, (lx(ib2)+1)**2
!           i2 = i2+1
!           do  ig = 1, ngc
!             pjylp(i2,ig) = pjyl(lm1,ig)*phase(ig,ib2)
!           enddo
!         enddo
!       enddo
      deallocate(yl)
!       deallocate(pjyl,phase,yl)

C ... On-site contribution to <P_G|v|P_G>
C     ajr(ir,l,ibas,igc) = j_l(sqrt(e) r) * r / (sqrt(e))**l
      allocate(ajr(nrx,ngc,0:lxx)) ! Spherical Bessel function
      allocate(a1(nrx,ngc,0:lxx))  ! integral of (spherical bessel)*(rkp rkm)
      ig = idalloc('ajr+a1',iverbose+2,2*nrx*(lxx+1),ngc)

      allocate(rijj(ngc,ngc,0:lxx)) ! Integrals of bessel * v * bessel
C     ig = idalloc('rijj',iverbose+2,ngc*(lxx+1),ngc)
C     rijj = -999

      allocate(wjj(ngc,ngc,0:lxx)) ! Wronskians between W(j_l(ig1),j_l(ig2)) at rmt
C     ig = idalloc('wjj',iverbose+2,ngc*(lxx+1),ngc)

      i2 = 0
      do  ib2 = 1, nbas

        call radmwt(100,rmax(ib2),aa(ib2),nr(ib2),rofi(1,ib2),wt) ! Integration weights
        if(verbose()>=30) call cputid2(' vcoul radmwt ',0)
C   ... ajr = Spherical Bessel function
C   ... a1 for integrals like \int_(r)^(rmax) rk(p,m)r(r') j_l(q2,r') dr'
C       Scale a1 by radial integration weights to prepare for radial integration I(a1(ig1),ajr(ig2))
        if (eee /= 0d0) then

        allocate(xi(nr(ib2),0:lx(ib2)),y(nr(ib2)),h(nr(ib2)))
        do  ig1 = 1, ngc

          call ropbes(rofi(1,ib2),qg2(ig1),lx(ib2),y,h,xi,nr(ib2),10*jobBes+1)
          do  ir = 1, nr(ib2)
C           call bessl(qg2(ig1)*rofi(ir,ib2)**2,lxx,phi,psi)
            xx = rofi(ir,ib2)
            do  l = 0, lx(ib2)
C             ajr(ir,ig1,l) = phi(l) * rofi(ir,ib2)**(l+1)
              ajr(ir,ig1,l) = xi(ir,l) * xx
              xx = xx*rofi(ir,ib2)
            enddo ! l
          enddo ! ir

          do  l = 0, lx(ib2)
            call intn_smpxxx(rkpr(1,l,ib2),ajr(1,ig1,l),int1x,aa(ib2),bb(ib2),rofi(1,ib2),nr(ib2),0)
            call intn_smpxxx(rkmr(1,l,ib2),ajr(1,ig1,l),int2x,aa(ib2),bb(ib2),rofi(1,ib2),nr(ib2),0)
            a1(1,ig1,l) = 0d0
            do  ir = 1, nr(ib2)
              a1(ir,ig1,l) = wt(ir)*(rkmr(ir,l,ib2)*(int1x(1)-int1x(ir))+rkpr(ir,l,ib2)*int2x(ir))
            enddo ! ir
          enddo ! l

        enddo ! ig1

        deallocate(xi,y,h)
        if(verbose()>=30) call cputid2(' vcoul dgemm start',0)
C       Matrix elements of the screened coulomb interaction with (j(iq1),j(iq2)) (?)
        do  l = 0, lx(ib2)
C         call dqsmpy(2,'T','N',ngc,nr(ib2),a1(1,1,l),nrx,ajr(1,1,l),nrx,0d0,rijj(1,1,l),ngc)
          call dgemm('T','N',ngc,ngc,nr(ib2),1d0,a1(1,1,l),nrx,ajr(1,1,l),nrx,0d0,rijj(1,1,l),ngc)
        enddo

        endif ! eee /= 0
        if(verbose()>=30) call cputid2(' vcoul wronjje start',0)
C   ... Wronskian wjj = W{qg2(ig1),qg2(ig2)} at rmax(ib2) for l=0:lx(ib2)
        call wronjje(jobBes,qg2,qg2,rmax(ib2),ngc,ngc,lx(ib2),wjj)

        if(verbose()>=30) call cputid2(' vcoul PP1 start ',0)
C   ... Loop ig1 = 1, ngc to make vcoul(nbloch+ig1,nbloch+ig2)
        do  ig1 = 1, ngc
!$omp parallel  private(ig2,sigx,radsig,l,m,lm,wv12,wv21,p12,p21,v12)
!$omp do

CC    ... Wronskian kept for a single ig1 to save space (at slower execution speed)
C         call wronjje(jobBes,qg2(ig1),qg2,rmax(ib2),1,ngc,lx(ib2),wjj)

C     ... Loop ig2 = 1, ig1 to make vcoul(nbloch+ig1,nbloch+ig2).  Lower triangle since vcoul is hermitian
          do  ig2 = 1, ig1

C       ... On-site sigma integral
            if (eee == 0d0) then
C             sig(l) = \int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} * j_l(qg2 r)/qg2**l
              call sigintpp(sqrt(qg2(ig1)),sqrt(qg2(ig2)),lx(ib2),rmax(ib2),sigx)
            else
              do  l = 0, lx(ib2)
                sigx(l) = rijj(ig2,ig1,l)
C                sigx(l) = 0
C                do  ir = 1, nr(ib2)
C                  sigx(l) = sigx(l) + a1(ir,ig1,l)*ajr(ir,ig2,l)
C                enddo
              enddo
            endif
            do  l = 0, lx(ib2)
C              if (abs(sigx(l)-rijj(ig2,ig1,l)) > 1d-10) then
C                print *, ig1,ig2,sigx(l),sigx(l)-rijj(ig2,ig1,l)
C                stop 'oops'
C              endif
              radsig(l) = fpi/(2*l+1) * sigx(l)
            enddo

            lm = 0; v12 = 0
            do  l = 0, lxx
C             Scale wjj by fpi/(qg2(ig1)-eee)
              wv12 = fpi/(qg2(ig1)-eee)*wjj(ig1,ig2,l)
              wv21 = fpi/(qg2(ig2)-eee)*wjj(ig2,ig1,l) ! Could exploit wjj(ig1,ig2,l) = wjj(ig2,ig1,l)
              p12 = 0 ; p21 = 0
              do  m = -l, l
                lm = lm+1

C                fouvp =  pjylp(ig1,lm) | W{j(qg2(ig1),qg2(ig2),rmt(ibas),l} | pjylp(ig2,lm)
C                fouvp12 = -wv12 * (dconjg(pjylp(i2+lm,ig1)) * pjylp(i2+lm,ig2))
C                fouvp21 = -wv21 * (dconjg(pjylp(i2+lm,ig2)) * pjylp(i2+lm,ig1))
C                sgpp12  =  radsig(l) * (dconjg(pjylp(i2+lm,ig1)) * pjylp(i2+lm,ig2))
CC               print 355, lm,ig2,ig1,ib2,sgpp12
C  355           format(4i4,2f22.12)
C                vcoul(nbloch+ig1,nbloch+ig2) = vcoul(nbloch+ig1,nbloch+ig2) + sgpp12 - dconjg(fouvp21) - fouvp12
                p12 = p12 + dconjg(pjyl(lm,ig1)*phase(ig1,ib2)) * pjyl(lm,ig2)*phase(ig2,ib2)
!                 p12 = p12 + dconjg(pjylp(i2+lm,ig1)) * pjylp(i2+lm,ig2)
              enddo             ! m
              v12 = v12 + (radsig(l) + wv12 + wv21) * p12
            enddo               ! l
            vcoul(nbloch+ig1,nbloch+ig2) = vcoul(nbloch+ig1,nbloch+ig2) + v12
          enddo                 ! ig2

!$omp end do
!$omp end parallel
        enddo                   ! ig1
        if(verbose()>=30) call cputid2(' vcoul PP1 end ',0)
        i2 = i2 + (lxx+1)**2
      enddo                     ! ib2

C      print "(' sumcheck vcoulq(pp)/1000', 2f18.6)", sum(vcoul(nbloch+1:nblochpmx,nbloch+1:nblochpmx))/1000
C      call zprm('vcoul (pp)',vcoul(nbloch+1,nbloch+1),nblochpmx,ngc,ngc)
      if (verbose()>=40) call cputid2(' complete vcoul(PP)',0)
C     ig = idalloc('rijj',iverbose+4,ngc*(lxx+1),ngc)
C     ig = idalloc('wjj',iverbose+4,ngc*(lxx+1),ngc)
      ig = idalloc('ajr+a1',iverbose+4,(2*nrx+ngc)*(lxx+1),ngc)
      deallocate(ajr,a1,wjj)
!       deallocate(pjylp,ajr,a1,wjj)

C     Re-entry point if ngc=0
   99 continue

C --- Right-Upper triangle of vcoul ---
      if(verbose()>=30) call cputid2(' vcoul transpose strt',0)
!$omp parallel do private(i1,i2)
      do  i1 = 1, nbloch+ngc
        do  i2 = 1, i1-1
          vcoul(i2,i1) = dconjg(vcoul(i1,i2))
        enddo
      enddo
!$omp end parallel do
      if(verbose()>=30) call cputid2(' vcoul transpose end',0)

      do  i1 = 1, nbloch+ngc
        if (i1==11) write(6,"(' ... ')")
        if(i1>10 .and. i1<nbloch+ngc .and. mod(i1,max(10,(nbloch+ngc)/20)) /= 0) cycle
        write(6,"(' Diagonal Vcoul =',i5,1p,2d18.10)") i1,vcoul(i1,i1)
      enddo

C      print "(' sumcheck vcoul/1000', 2f18.6)", sum(vcoul(1:nbloch+ngc,1:nbloch+ngc))/1000
C     call zprm('vcoul',vcoul,nblochpmx,nbloch+ngc,nbloch+ngc)

C --- Check Coulomb by plane wave expansion ---
#ifdef ptest
      write(6,*) ' --- plane wave Coulomb matrix check --- '
      write(197,*) ' --- off diagonal ---- '
      nblochngc = nbloch+ngc
      allocate(matp(nblochngc),matp2(nblochngc))
      do  ig1 = 1, ngc
        matp = 0d0
        do ibl2= 1, nbloch
          ib2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp(ibl2) = fouvb(ig1, n2, lm2, ib2)*qg2(ig1)/fpi
        enddo
        matp(nbloch+ig1) = 1d0
        ig2=ig1
c      do ig2 = 1,ngc !off diagonal
        matp2 = 0d0
        do  ibl2= 1, nbloch
          ib2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp2(ibl2) = fouvb(ig2, n2, lm2, ib2)*qg2(ig2)/fpi
        enddo
        matp2(nbloch+ig2) = 1d0
        v12= sum(
     &    matmul(matp(1:nblochngc),vcoul(1:nblochngc,1:nblochngc))
     &             *dconjg(matp2(1:nblochngc))  )
        if(ig1/=ig2) then  !off diagonal
          if(abs(v12)>1d-1 ) then
            write(197,'(2i5, 2d13.6)') ig1,ig2, v12
            write(197,'("    matpp ", 2d13.6)')
     &      vcoul(nbloch+ig1,nbloch+ig2)
            write(197,*)
          endif
        else
          write(196,'(2i5," exact=",3d14.6," q ngsum=",3f8.4,i5)')
     &      ig1,ig2,fpi*vol/qg2(ig1)
     &    , fpi*vol/qg2(ig2),qg2(ig1), q(1:3)
     &    , sum(ngvecc(1:3,ig1)**2)
          write(196,'("           cal  =", 2d14.6)') v12
          write(196,'("           vcoud=", 2d14.6)')
     &       vcoul(nbloch+ig1,nbloch+ig2)
          write(196,*)
        endif
c      enddo !off diagonal
      enddo

      deallocate(matp,matp2)
      stop ' *** ptest end *** See fort.196 and 197'
#endif

      end
