      subroutine vcoulq_4(q,nbloch,ngc,nbas,lx,lxx,nx,nxx,alat,qlat,vol,ngvecc,
     &  strx,rojp,rojb,sgbb,sgpb,fouvb,
     i  nblochpmx,bas,rmax,
     i  eee,aa,bb,nr,nrx,rkpr,rkmr,rofi,
     o  vcoul)
C- Coulmb matrix for each q.
C  Superseded by vcoulq, which see.
Ci strx:  Structure factors
Ci nlx corresponds to (lx+1)**2 . lx corresponds to 2*lmxax.
Ci rho-type integral
Ci  ngvecc     : q+G vector
Ci  rojp rojb  : rho-type integral
ci  sigma-type onsite integral
ci  Fourier
Ci  nx(l,ibas) : max number of radial function index for each l and ibas.
Ci               Note that the definition is a bit different from nx in basnfp.
ci  nxx        : max number of nx among all l and ibas.
ci  lx(nbas)   : max number of l for each ibas.
ci  lxx        :
ci
ci  vol : cell vol
c
Co Vcoul
cr vcoul is in a.u. You have to multiply e~2=2 if you want to it in Ry,
cr    vcoul = 2d0*vcoul !  in Ry unit.
c---------------------------------------------------------------------------
c  rojp = <j_aL(r) | P(q+G)_aL > where
c         |P(q+G)_aL> : the projection of exp(i (q+G) r) to aL channnel.
c         |j_aL>      : \def r^l/(2l+1)!! Y_L.  The spherical bessel functions near r=0.  Energy-dependence is omitted.
c

      use mpi_mod, q_comm => comm_1, q_inter => inter_1, rank => mpi_rank


      implicit none

      integer :: nbloch, nblochpmx, nbas,
     &              lxx,lx(nbas), nxx, nx(0:lxx,nbas)
      real(8)    :: vol,q(3),fpi

ci structure con
      complex(8) :: strx((lxx+1)**2, nbas, (lxx+1)**2,nbas)
ci |q+G|**2
      integer :: ngc, ngvecc(3,ngc)
      real(8)    :: qlat(3,3),alat,absqg2(ngc),qg(3),qgn(3)

ci rho-type onsite integral
      complex(8) ::   rojp(ngc, (lxx+1)**2, nbas)
      real(8)    ::   rojb(nxx, 0:lxx, nbas)
ci sigma-type onsite integral
      real(8)    :: sgbb(nxx,  nxx,  0:lxx,      nbas)
      complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2, nbas)
c     &              ,sgpp(ngc,  ngc,  (lxx+1)**2, nbas)
ci Fourier
      complex(8) ::
     &             fouvb(ngc,  nxx, (lxx+1)**2, nbas)
Co
     &             ,vcoul(nblochpmx, nblochpmx)
c     &             ,fouvp(ngc,  ngc, (lxx+1)**2, nbas)

cinternals
      integer :: ibl1, ibl2,ig1,ig2,ibas,ibas1,ibas2,
     &              l,m,n, n1,l1,m1,lm1,n2,l2,m2,lm2,ipl1,ipl2,jobw
      integer :: ibasbl(nbloch), nbl(nbloch), lbl(nbloch),
     &              mbl(nbloch), lmbl(nbloch)
      real(8) :: pi, fpivol,tpiba
      complex(8) :: rojpstrx((lxx+1)**2,nbas)

c check
C      complex(8),allocatable :: hh(:,:),oo(:,:),zz(:,:)
C      real(8),allocatable    :: eb(:)

      complex(8),allocatable :: matp(:),matp2(:)
      complex(8) :: xxx
      integer :: nblochngc,ix
      logical :: ptest=.false. ! See ptest in basnfp.f

c-------------------
      real(8),   allocatable :: cy(:),yl(:)
      complex(8),allocatable :: pjyl_(:,:),phase(:,:)
      complex(8) :: img=(0d0,1d0)
      real(8):: bas(3,nbas),r2s,rmax(nbas)
      integer:: lm
#ifdef COMMONLL
      integer::ll(51**2)
      common/llblock/ll
#else
      integer :: ll
      external ll
#endif
      real(8)::  fkk(0:lxx),fkj(0:lxx),fjk(0:lxx),fjj(0:lxx),sigx(0:lxx),radsig(0:lxx)
      complex(8):: fouvp12, fouvp21, sgpp12

      integer:: nrx,nr(nbas),ir,ig,jobBes
      real(8):: eee , int1x(nrx),int2x(nrx),phi(0:lxx),psi(0:lxx)
     &  ,aa(nbas),bb(nbas),rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas)
     &  ,rofi(nrx,nbas)
      real(8), allocatable:: ajr(:,:,:,:), a1(:,:,:)
      logical :: debug=.false.
      integer:: verbose
      procedure(integer) :: GWversion

!ML MPI related variables
      integer myID
      integer, allocatable :: first(:), last(:), displace(:), recvcount(:)
      double precision a_coeff, b_coeff
      complex(8), allocatable :: vcoul_local(:,:)

      write(6,'(" vcoulq_4: nblochpmx  nbloch ngc=",3i6)') nblochpmx,nbloch,ngc
c     print *, ' sum fouvp=',sum(fouvp(:,:,:,1))
c     print *, ' sum fouvb=',sum(fouvb(:,:,:,1))
      pi    = 4d0*datan(1d0)
      fpi    = 4*pi
      fpivol = 4*pi*vol

c---for sgpp fouvp
      allocate(  !ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc),rkpr(nr,0:lx),rkmr(nr,0:lx),
     &  pjyl_((lxx+1)**2,ngc),phase(ngc,nbas) )
      allocate(cy((lxx+1)**2),yl((lxx+1)**2))
      call sylmnc(cy,lxx)

c=======================================================
      vcoul = 0d0
      phase = 0d0
      pjyl_ = 0d0
      jobBes = mod(GWversion()/1000,10)

c-gvec
      tpiba = 2*pi/alat

!      write(0,'("[",I4,"] vcoulq_p4. Check-point 1: nr = ",I5)') rank,nr
!      do ig1 = 1+q_comm%ID,ngc, q_comm%size
      do ig1 = 1,ngc
        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
        absqg2(ig1)  = sum(qg(1:3)**2)

c---for spgg fourvp ----------
        do ibas=1,nbas
          phase(ig1,ibas) = exp( img*sum(qg(1:3)*bas(1:3,ibas))*alat  )
        enddo

        if(sqrt(absqg2(ig1)) > 0) then
           qgn = qg/sqrt(absqg2(ig1))
        else
           qgn = 0d0
        endif
        call sylm(qgn,yl,lxx,r2s) !spherical factor Y( q+G )

        do lm =1,(lxx+1)**2
          l = ll(lm)
          pjyl_(lm,ig1) = fpi*img**l *cy(lm)*yl(lm)  * sqrt(absqg2(ig1))**l  !*phase
          ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
        enddo
      enddo ! ig1
c

#ifdef USE_MPI

!      call MPI_AllReduce(MPI_IN_PLACE, pjyl_, ngc*(lxx+1)**2,
!     &     MPI_DOUBLE_COMPLEX, MPI_SUM,
!     &     q_comm%communicator, ierror)

!      call MPI_AllReduce(MPI_IN_PLACE, phase, ngc*nbas,
!     &     MPI_DOUBLE_COMPLEX, MPI_SUM,
!     &     q_comm%communicator, ierror)
#endif

!      write(0,'("[",I4,"] vcoulq_p4. Check-point 2")') rank
!ML ----------------------------------------------------------------------------
!ML MPI remark:    This loop NEEDS to be executed in sequential order,
!ML                and the index arrays need to be known on all processors

c-- index (mx,nx,lx,ibas) order.
      ibl1 = 0
      do ibas= 1, nbas
        do l   = 0, lx(ibas)
c        write(6,'(" l ibas nx =",3i5)') l,nx(l,ibas),ibas
          do n   = 1, nx(l,ibas)
            do m   = -l, l
              ibl1  = ibl1 + 1
              ibasbl(ibl1) = ibas
              nbl   (ibl1) = n
              lbl   (ibl1) = l
              mbl   (ibl1) = m
              lmbl  (ibl1) = l**2 + l+1 +m
c        write(6,*)ibl1,n,l,m,lmbl(ibl1)
            enddo
          enddo
        enddo
      enddo
      if(ibl1/= nbloch) then
        write(6,*)' ibl1 nbloch',ibl1, nbloch
        call rx( ' vcoulq: error ibl1/= nbloch')
      endif


c-- <B|v|B> block
c      write(6,*)' vcoulq: bvb block xxx rojbsum='
c      write(6,*) sum(rojb(:,:,1))
c      write(6,*) sum(rojb(:,:,2))
c      write(6,*) sum(rojb(:,:,3))
c      write(6,*) sum(rojb(:,:,4))
c      write(6,*)' vcoulq: bvb block xxx sgbbbsum='
c      write(6,*) sum(sgbb(:,:,:,1))
c      write(6,*) sum(sgbb(:,:,:,2))
c      write(6,*) sum(sgbb(:,:,:,3))
c      write(6,*) sum(sgbb(:,:,:,4))




!ML This loop can be large !!
!ML
!ML  can be parallelised over ibl1
!ML
!ML  vcoul = 0 at thus stage. Each processor can just add to its bit
!ML  and we can collect all with an AllReduce


!ML  Turn around the loops:
!ML
!ML  Original form:
!ML
!ML      do ibl1= 1, nbloch
!ML         do ibl2= 1, ibl1
!ML
!ML  this means that we always have ibl2 <= ibl1.
!ML
!ML  new form:
!ML
!ML      do ibl2 = 1,nbloch
!ML         do ibl1 = ibl2,nbloch
!ML
!ML  also here ibl1 runs from 1 to nbloch, and we have ibl2 >= ibl1
!ML
!ML  WIth this order, we can split Vcoul with respect to ibl2, which leaves contiguous
!ML  slices of memory.

!      write(0,'("[",I4,"] vcoulq_p4. Check-point 4")') rank

#ifdef USE_MPI

      a_coeff = 1.d0
      b_coeff = 0.d0

      allocate( first(1:q_comm%size) )
      allocate( last(1:q_comm%size) )
      allocate( displace(1:q_comm%size) )
      allocate( recvcount(1:q_comm%size) )

      call thdist1(1,a_coeff,b_coeff,q_comm%size,1,nbloch,first,last)

      displace(:) = nblochpmx*(first(:)-1)
      recvcount(:) = (last(:) - first(:) + 1) * nblochpmx


      myID = q_comm%ID+1


      write(0,'("[",I4,"] I am ",I4," and run the slice ",I4,":",I4," nblochpmx = ",I5)')
     &     rank,myID,first(myID),last(myID), nblochpmx

!ML allocate local slice of vcoul:

      allocate( vcoul_local(nblochpmx, first(myID):last(myID)))

      vcoul_local = 0.d0

!ML original:      do ibl1= 1, nbloch

!      do ibl2= 1, nbloch
      do ibl2= first(myID), last(myID)

        ibas2= ibasbl(ibl2)
        n2   = nbl (ibl2)
        l2   = lbl (ibl2)
        m2   = mbl (ibl2)
        lm2  = lmbl(ibl2)

        do ibl1= ibl2, nbloch

          ibas1= ibasbl(ibl1)
          n1   = nbl (ibl1)
          l1   = lbl (ibl1)
          m1   = mbl (ibl1)
          lm1  = lmbl(ibl1)

          vcoul_local(ibl1,ibl2) =
     &     rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2)
     &    *rojb(n2, l2, ibas2)

          if(ibas1==ibas2 .and. lm1==lm2) then
            vcoul_local(ibl1,ibl2) = vcoul_local(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1)
            ! sigma-type contribution. onsite coulomb
          endif

        enddo ! ibl1

      enddo ! ibl2


!      write(0,'("[",I4,"] trying to send   ",I7," bytes; ",I7)') rank, nblochpmx*(last(myID)-first(myID) +1),
!     &     size(vcoul_local)
!      write(0,'("[",I4,"] trying to receive",I7," bytes at ",I7)') rank, recvcount(myID),displace(myID)
!
      call MPI_AllGatherV(vcoul_local, recvcount(myID), MPI_DOUBLE_COMPLEX,
     &     vcoul,recvcount, displace,  MPI_DOUBLE_COMPLEX, q_comm%communicator, ierror)


      deallocate (first, last, displace, recvcount, vcoul_local)



#else

      do ibl1= 1, nbloch
        ibas1= ibasbl(ibl1)
        n1   = nbl (ibl1)
        l1   = lbl (ibl1)
        m1   = mbl (ibl1)
        lm1  = lmbl(ibl1)

        do ibl2= 1, ibl1
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)

          vcoul(ibl1,ibl2) =
     &     rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2)
     &    *rojb(n2, l2, ibas2)
          if(ibas1==ibas2 .and. lm1==lm2) then
            vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1)
            ! sigma-type contribution. onsite coulomb
          endif
        enddo ! ibl2
      enddo ! ibl1
#endif

!     write(0,'("[",I4,"] vcoulq_p4. Check-point 5")') rank


ccccccccccccccccccccccccccc
c      goto 1112
ccccccccccccccccccccccccccc

c <P_G|v|B>
      if(debug) write(6,*)' vcoulq_4: pgvb block 1111'
      do ibl2= 1, nbloch
        ibas2= ibasbl(ibl2)
        n2   = nbl (ibl2)
        l2   = lbl (ibl2)
        m2   = mbl (ibl2)
        lm2  = lmbl(ibl2)
        do ig1 = 1,ngc
          ipl1 = nbloch + ig1
          vcoul(ipl1,ibl2) = fouvb(ig1,  n2, lm2, ibas2)

          do ibas1= 1, nbas
            do lm1  = 1, (lx(ibas1)+1)**2
              vcoul(ipl1,ibl2) = vcoul(ipl1,ibl2) -
     &      dconjg(rojp(ig1, lm1, ibas1)) *strx(lm1,ibas1,lm2,ibas2)
     &      *rojb(n2, l2, ibas2)
              if(ibas1==ibas2 .and.lm1==lm2) then
                vcoul(ipl1,ibl2) = vcoul(ipl1,ibl2) -
     &        sgpb(ig1, n2, lm2, ibas2)
              endif
            enddo               !  lm1
          enddo                 !  ibas1
        enddo                   ! ig1
      enddo                     ! ibl2

!     write(0,'("[",I4,"] vcoulq_p4. Check-point 6")') rank


C ... Make function ajr and a1.
C     ajr:spherical bessel, a1: integral of (sperical bseel)*(rkp rkm)
      if(debug) write(6,*)' vcoulq_4: ajr allocate'
      allocate( ajr(nrx,0:lxx,nbas,ngc), a1(nrx,0:lxx,nbas) )
      if(debug) write(6,*)' vcoulq_4: end ajr allocate'
      do ig1 = 1,ngc
        do ibas= 1,nbas
          if(debug) write(6,"('ccc: ',10i15)")ig1,ibas
          do ir = 1,nr(ibas)
C           call bessl(absqg2(ig1)*rofi(ir,ibas)**2,lxx,phi,psi)
            call besslr(absqg2(ig1)*rofi(ir,ibas)**2,10*jobBes,0,lxx,phi,psi)
            do l  = 0,lx(ibas)

              if(debug.and.ig==162.and.ibas==8) then
                write(6,"('ccc: ',10i15)")ig1,ibas,ir,l
                write(6,*)"ccc:", phi(l)
                write(6,*)"ccc:", rofi(ir,ibas)
              endif

              ajr(ir,l,ibas,ig1) = phi(l)* rofi(ir,ibas) **(l +1 )
              ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
            enddo
          enddo
        enddo
      enddo

!      write(0,'("[",I4,"] vcoulq_p4. Check-point 7")') rank

c --- <P_G|v|P_G> ---
      if(debug) write(6,*)' vcoulq_4: pgvpg block'
      if(verbose()>=40) call cputid2(' start pgvpg',0)

      do  ig1 = 1, ngc

!         write(0,'("[",I4,"] vcoulq_p4. Check-point 7a: ig1,ngc = ",2I4)') rank, ig1,ngc

        ipl1 = nbloch + ig1

C   ... Make rojpstrx
        rojpstrx = 0d0
        do  ibas1 = 1, nbas
          do  lm1  = 1, (lx(ibas1)+1)**2
            do  ibas2 = 1, nbas
              do  lm2  = 1, (lx(ibas2)+1)**2
                rojpstrx(lm2, ibas2) = rojpstrx(lm2, ibas2) +
     &      dconjg(rojp(ig1, lm1, ibas1)) *strx(lm1,ibas1,lm2,ibas2)
              enddo
            enddo
          enddo
        enddo

!         write(0,'("[",I4,"] vcoulq_p4. Check-point 7b")') rank

C   ... Make int1x,int2x,a1

!         write(0,'("[",I4,"] vcoulq_p4. aa(ibas) = ",10E10.2)') rank, (aa(ibas), ibas=1,min(5,nbas))
!         write(0,'("[",I4,"] vcoulq_p4. bb(ibas) = ",10E10.2)') rank, (bb(ibas), ibas=1,min(5,nbas))
!         write(0,'("[",I4,"] vcoulq_p4. nr(ibas) = ",10I8)') rank, (nr(ibas), ibas=1,min(5,nbas))

        do  ibas = 1, nbas

!            write(0,'("[",I4,"] vcoulq_p4. Check-point 7b: ibas,nbas = ",2I4)') rank, ibas, nbas
!            write(0,'("[",I4,"] vcoulq_p4. Check-point 7b: ibas,nr = ",2I4)') rank, nr(ibas)
!            write(0,'("[",I4,"] vcoulq_p4. Check-point 7b: ibas,lx = ",2I4)') rank, lx(ibas)
!
!            write(0,'("[",I4,"] vcoulq_p4. rofi(i",I4,") = ",10E10.2)') rank, ibas, (rofi(i,ibas), i=1,min(5,nr(ibas)))

          do  l = 0, lx(ibas)


!               write(0,'("[",I4,"] vcoulq_p4: rkpr(i,",I4,",",I4,") = ",10E10.2)') rank, l, ibas, (rkpr(i,l,ibas),i=1,5)
!               write(0,'("[",I4,"] vcoulq_p4: ajr(i,",I4,",",I4,",",I4,") = ",10E10.2)') rank, l, ibas, ig1,
!     &              (ajr(i,l,ibas,ig1),i=1,5)

            call intn_smpxxx( rkpr(1,l,ibas), ajr(1,l,ibas,ig1),int1x,aa(ibas),bb(ibas),rofi(1,ibas),nr(ibas),0)
            call intn_smpxxx( rkmr(1,l,ibas), ajr(1,l,ibas,ig1),int2x,aa(ibas),bb(ibas),rofi(1,ibas),nr(ibas),0)
            a1(1,         l,ibas) = 0d0
            a1(2:nr(ibas),l,ibas) =
     &             rkmr(2:nr(ibas),l,ibas)*(int1x(1)-int1x(2:nr(ibas)))
     &           + rkpr(2:nr(ibas),l,ibas)* int2x(2:nr(ibas))

          enddo                 ! l
        enddo                   ! ibas

!      write(0,'("[",I4,"] vcoulq_p4. Check-point 7c")') rank

C Something is ill defined in this loop ... OMP doesn't work.   Wasn't able to identify the problem
!!C$OMP parallel do default(shared) private(ig2,ipl2,ibas2,l,lm2,phi,psi,fkk,fkj,fjk,fjj,sigx,radsig,fouvp12,fouvp21,sgpp12)
        do  ig2 = 1, ig1
          ipl2 = nbloch + ig2
          if (ig1==ig2) vcoul(ipl1,ipl2) = fpivol/(absqg2(ig1) - eee) !eee is negative
          do  ibas2 = 1, nbas
C       ... for fouvp and sgpp
            call wronkj(jobBes,absqg2(ig1),absqg2(ig2),rmax(ibas2),lx(ibas2),
     o                 fkk,fkj,fjk,fjj)

            if(eee==0d0) then
              call sigintpp(absqg2(ig1)**.5,absqg2(ig2)**.5,lx(ibas2),rmax(ibas2),
     o             sigx)
            else
              do  l = 0, lx(ibas2)
                call gintxx(a1(1,l,ibas2),ajr(1,l,ibas2,ig2),
     .            aa(ibas2),bb(ibas2),NR(ibas2), sigx(l))
              enddo
            endif
            do  l = 0, lx(ibas2)
              radsig(l) = fpi/(2*l+1) * sigx(l)
            enddo

            do  lm2  = 1, (lx(ibas2)+1)**2

              l= ll(lm2)
c          ...fouvp sgpp
              fouvp12 = fpi/(absqg2(ig1)-eee) *
     .          dconjg(pjyl_(lm2,ig1)*phase(ig1,ibas2)) *
     .          (-fjj(l)) * pjyl_(lm2,ig2)*phase(ig2,ibas2)
              fouvp21 = fpi/(absqg2(ig2)-eee) *
     .          dconjg(pjyl_(lm2,ig2)*phase(ig2,ibas2)) *
     .          (-fjj(l)) * pjyl_(lm2,ig1)*phase(ig1,ibas2)
              sgpp12  = dconjg(pjyl_(lm2,ig1)*phase(ig1,ibas2))*
     .          radsig(l)*pjyl_(lm2,ig2)*phase(ig2,ibas2)

              vcoul(ipl1,ipl2) = vcoul(ipl1,ipl2)
     .       +  rojpstrx(lm2,ibas2)*rojp(ig2, lm2, ibas2)
c     &      -  dconjg( fouvp(ig2,  ig1, lm2, ibas2)) !BugFix Mar5-01 It was dcmplx.
c     &      -          fouvp(ig1,  ig2, lm2, ibas2)
c     &      +  sgpp(ig1, ig2, lm2, ibas2)
     .       -  dconjg( fouvp21 )
     .       -          fouvp12
     .       +  sgpp12

            enddo               ! lm2
          enddo                 ! ibas2
        enddo                   ! ig2
!!C$OMP end parallel do
      enddo                     ! ig1



C      call zprm('vcoul',vcoul,nblochpmx,nblochpmx,nblochpmx)
C      call rx('done')

ccccccccccccccccccccccccccccccc
c 1112 continue
ccccccccccccccccccccccccccccccc

!      write(0,'("[",I4,"] vcoulq_p4. Check-point 8")') rank


c-- Right-Upper part of vcoul.
      if(debug) write(6,*)' vcoulq_4: right-upper'
      do ipl1=1, nbloch+ngc
        do ipl2=1, ipl1-1
          vcoul(ipl2,ipl1) = dconjg(vcoul(ipl1,ipl2))
        enddo
      enddo

cccccccccccccccccccccccccccccc
c test.xxxxxxxxxx
c$$$      do ibl2= 1, nbloch
c$$$        ibas2= ibasbl(ibl2)
c$$$        n2   = nbl (ibl2)
c$$$        l2   = lbl (ibl2)
c$$$        m2   = mbl (ibl2)
c$$$        lm2  = lmbl(ibl2)
c$$$        if(l2==1.and.ibas2>2) then
c$$$          vcoul(nbloch+1:nbloch+ngc, ibl2) = 0d0
c$$$          vcoul(ibl2, nbloch+1:nbloch+ngc) = 0d0
c$$$        endif
c$$$      enddo
ccccccccccccccccccccccccccc

c vcoul is in a.u. You have to multiply e~2=2 if you want to it in Ry,
c     vcoul = 2d0*vcoul !  in Ry unit.
c

c check write
      do ix = 1,nbloch+ngc,20
        write(6,"(' Diagonal Vcoul =',i5,2d18.10)") ix,vcoul(ix,ix)
      enddo
      if( allocated(yl)   ) deallocate(yl)
      if( allocated(cy)   ) deallocate(cy)
      if( allocated(phase)) deallocate(phase)
      if( allocated(pjyl_)) deallocate(pjyl_)
      if(.not.ptest) return



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c! Below ia a plane-wave test.
c--- check! Coulomb by plane wave expansion.
      write(6,*) ' --- plane wave Coulomb matrix check 1---- '
      write(197,*) ' --- off diagonal ---- '
      nblochngc = nbloch+ngc
      allocate(matp(nblochngc),matp2(nblochngc))
      do ig1 = 1,ngc
        matp = 0d0
        do ibl2= 1, nbloch
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp(ibl2) = fouvb(ig1, n2, lm2, ibas2)*absqg2(ig1)/fpi
        enddo
        matp(nbloch+ig1) = 1d0
        ig2=ig1
c      do ig2 = 1,ngc !off diagnal
        matp2 = 0d0
        do ibl2= 1, nbloch
          ibas2= ibasbl(ibl2)
          n2   = nbl (ibl2)
          l2   = lbl (ibl2)
          m2   = mbl (ibl2)
          lm2  = lmbl(ibl2)
          matp2(ibl2) = fouvb(ig2, n2, lm2, ibas2)*absqg2(ig2)/fpi
        enddo
        matp2(nbloch+ig2) = 1d0
        xxx= sum(
     &    matmul(matp(1:nblochngc),vcoul(1:nblochngc,1:nblochngc))
     &             *dconjg(matp2(1:nblochngc))  )
        if(ig1/=ig2) then  !off diagnal
          if(abs(xxx)>1d-1 ) then
            write(197,'(2i5, 2d13.6)') ig1,ig2, xxx
            write(197,'("    matpp ", 2d13.6)')
     &      vcoul(nbloch+ig1,nbloch+ig2)
            write(197,*)
          endif
        else
          write(196,'(2i5," exact=",3d13.6,"q ngsum=",3f8.4,i5)')
     &       ig1,ig2,fpi*vol/absqg2(ig1)
     &    , fpi*vol/absqg2(ig2),absqg2(ig1), q(1:3)
     &    , sum(ngvecc(1:3,ig1)**2)
          write(196,'("           cal  =", 2d13.6)') xxx
          write(196,'("           vcoud=", 2d13.6)')
     &       vcoul(nbloch+ig1,nbloch+ig2)
          write(196,*)
        endif
c      enddo !off diagnal
      enddo
c
!      write(0,'("[",I4,"] vcoulq_p4. Check-point 9")') rank

      deallocate(matp,matp2)
c      stop ' *** ptest end *** See fort.196 and 197'
      end

      subroutine mkjp_4(q,ngc,ngvecc, alat, qlat, lxx,lx,nxx,nx,
     i                  bas, a,b,rmax,nr,nrx,rprodx,eee,rofi,rkpr,rkmr,
     o                  rojp,sgpb,fouvb)
C- Make integrals in each MT and the Fourier matrix fouvb
C ----------------------------------------------------------------------
Ci Inputs
Ci   q     :wave number
Ci   ngc   :number of G vectors
Ci   ngvecc:q+G vectors
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   lxx   :dimensons nx,rprodx,rojp,sgpb
Ci   lx    :l-cutoff for rprodx,rkpr,rojp
Ci   nxx   :dimensions fouvb,sgpb
Ci   nx    :l-dependent number of radial functions in product basis
Ci   bas   :basis vectors, in units of alat
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rmax  :augmentation radius, in a.u.
Ci   nr    :number of radial mesh points
Ci   nrx   :Leading dimension of rprodx
Ci   rprodx:Product basis functions * r: r*BRlμ(r)
Ci   eee   :Yukawa screening energy
Ci   rofi  :radial mesh points
Ci   rkpr  :r*Jl(E,r)
Ci   rkmr  :r*Hl(E,r)
Co Outputs
Co   rojp  :the analog of mkjb_4 in mkjb_4.  I think (unsure about prefactors)
Co         :rojp(G,l) = 1/(2l + 1)!! int dr r*Jl(E,r) r*Jl(G^^2,r)
Co   sgpb  :the analog of sgbb in mkjb_4.
Co   fouvb :4*pi/(G^2-E) * integral(J*B)
Cl Local variables
Cl         :
Cr Remarks
Cr the integrals rojp, fouvb, fouvp
Cr are for  J_L(r)= j_l(sqrt(e) r)/sqrt(e)**l Y_L,
Cr which behaves as r^l/(2l+1)!! near r=0.
Cr
Cr Onsite integral is based on
Cr   1/|r-r'| = \sum 4 pi /(2k+1) \frac{r_<^k }{ r_>^{k+1} } Y_L(r) Y_L(r')
Cr See PRB34 5512(1986) Eq. A12 for sigma type integral
Cu Updates
C ----------------------------------------------------------------------
      implicit none

      integer :: ngc,ngvecc(3,ngc), lxx, lx, nxx,nx(0:lxx),nr,nrx
      real(8)    :: q(3),bas(3), rprodx(nrx,nxx,0:lxx),a,b,rmax,alat,
     &              qlat(3,3)
      integer jobBes
      procedure(integer) :: GWversion

ci rho-type onsite integral
      complex(8) :: rojp(ngc, (lxx+1)**2)
ci sigma-type onsite integral
      complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2)
c     &              sgpp(ngc,  ngc,  (lxx+1)**2)
      real(8),allocatable::cy(:),yl(:)
ci Fourier
      complex(8) ::
     &              fouvb(ngc,  nxx, (lxx+1)**2)
c     &              fouvp(ngc,  ngc, (lxx+1)**2)
c internal
      integer :: nlx,ig1,l,n,ir,lm !, ibas
#ifdef COMMONLL
      integer::ll(51**2)
      common/llblock/ll
#else
      integer :: ll
      external ll
#endif
      real(8)    :: pi,fpi,tpiba, qg1(3),
     &  fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx),absqg1,
     &  radint,phi(0:lx),psi(0:lx)
     &   ,r2s,sig
      complex(8), parameter :: img =(0d0,1d0)
      complex(8) :: phase
      complex(8),allocatable :: pjyl(:,:)
      real(8),allocatable ::ajr(:,:,:),a1(:,:,:), !rkpr(:,:),rkmr(:,:),
     &  qg(:,:),absqg(:)
      real(8):: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx),qg1n(3),eee
      logical, parameter :: debug=.false.
C     Help arrays to make fast Bessel functions
      real(8), allocatable :: xi(:,:),y(:),h(:)

      if(debug) print *,' mkjp_4:'
      nlx = (lx+1)**2
      allocate(ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc),
     &  qg(3,ngc),absqg(ngc),
     &  pjyl((lx+1)**2,ngc) )

      pi    = 4d0*datan(1d0)
      fpi   = 4*pi
      tpiba = 2*pi/alat
      allocate(cy((lx+1)**2),yl((lx+1)**2))
      jobBes = mod(GWversion()/1000,10)
      call sylmnc(cy,lx)
!      print *,' mkjp_4: end of sylmnc'
C... q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
      do ig1 = 1,ngc
        qg(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
        qg1(1:3) = qg(1:3,ig1)
        absqg(ig1)  = sqrt(sum(qg1(1:3)**2))
        absqg1   = absqg(ig1)
        if (absqg1 > 0) then
           qg1n = qg1/absqg1
        else
           qg1n = 0.d0
        endif
        phase = exp( img*sum(qg1(1:3)*bas(1:3))*alat  )
        call sylm(qg1n,yl,lx,r2s) !spherical factor Y( q+G )
C        print *, qg1
C        print *, absqg1,r2s
C        do lm =1,nlx
C          l = ll(lm)
C          yl(lm) = cy(lm)*yl(lm)
C        enddo
C        print *, yl(1:9)
C        call ropyln(1,qg1(1)/absqg1,qg1(2)/absqg1,qg1(3)/absqg1,lx,1,yl,r2s)
C        print *, yl(1:9)
C       pjyl are scaling factors; see Eq 49 in ecalj manual (June 19, 2015)
        do lm =1,nlx
          l = ll(lm)
          pjyl(lm,ig1) = fpi*img**l *cy(lm)*yl(lm) *phase  *absqg1**l
          ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
        enddo
      enddo
C#ifdefC DEBUG
C      print "(' sumcheck pjyl', 2f18.6)", sum(pjyl)
C#endif
C     call zprm('pjyl',pjyl,(lx+1)**2,(lx+1)**2,ngc)

cc rofi and aj = r**l / (2l+1)!! \times r. Sperical Bessel at e=0.
c      rofi(1) = 0d0
c      do ir   = 1, nr
c        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
c      enddo
c      do l = 0,lx
c        rkpr(1:nr,l) = rofi(1:nr)**(l      +1 )
c        rkmr(2:nr,l) = rofi(2:nr)**(-l-1   +1 )
c        rkmr(1,l)    = rkmr(2,l)
c      enddo

c ... rojp(G,l) = pjyl * W{J(E), J(G^2)}
      if(debug) print *,' mkjp_4: rojp'
      do ig1 = 1,ngc
        call wronkj(jobBes, absqg(ig1)**2, eee, rmax,lx,
     o                  fkk,fkj,fjk,fjj)
        do lm = 1,nlx
          l = ll(lm)
          rojp(ig1,lm) = (-fjj(l))* pjyl(lm,ig1)
        enddo
      enddo
C#ifdefC DEBUG
C      print "(' sumcheck rojp', 2f18.6)", sum(rojp)
C#endif

C     ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
      allocate(xi(nr,0:lx),y(nr),h(nr))
      do ig1 = 1,ngc
        call ropbes(rofi,absqg(ig1)**2,lx,y,h,xi,nr,10*jobBes+1)
        do ir  = 1,nr
CC         call bessl(absqg(ig1)**2*rofi(ir)**2,lx,phi,psi)
C         call besslr(absqg(ig1)**2*rofi(ir)**2,10*jobBes,0,lx,phi,psi)
          r2s = rofi(ir)
          do l   = 0,lx
C           ajr(ir,l,ig1) = phi(l) * rofi(ir)**(l +1 )
            ajr(ir,l,ig1) = xi(ir,l) * r2s
            r2s = r2s*rofi(ir)
          enddo
        enddo
      enddo
C#ifdefC DEBUG
C      print "(' sumcheck jg', f18.6)", sum(ajr(nr,:,:))
C#endif
      deallocate(xi,y,h)

c-------------------------
      if(eee==0d0) then
c        print *,' mkjp_4: use sigintAn1 eee=0(r0c=infty) mode'
        do ig1 = 1,ngc
          call sigintAn1( absqg(ig1), lx, rofi, nr
     o                  ,a1(1:nr, 0:lx,ig1) )
        enddo
c      else
c We need to implement a version of sigintAn1 to treat eee/=0 case...
      endif

c-------------------------
c sgpb
      sgpb = 0
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)                      ! r jl        , r B(r)
            if(eee==0d0) then
              call gintxx(a1(1,l,ig1),rprodx(1,n,l),A,B,NR, sig )
ccccccccccccccccccc
c        write(6,"( ' sgpb= ',3i5,2d14.6)") ig1,n,lm, sgpb(ig1,n,lm)
ccccccccccccccccccc
            else !for a while, we use this version of sgpb
              call sigint_4(rkpr(1,l),rkmr(1,l), lx,a,b,nr, ajr(1,l,ig1),rprodx(1,n,l)
     &            , rofi, sig)
            endif
            sgpb(ig1,n,lm) = dconjg(pjyl(lm,ig1))* sig/(2*l+1)*fpi
ccccccccccccccccccc
c        write(6,"( ' sgpb= ',3i5,2d14.6)") ig1,n,lm, sgpb(ig1,n,lm)
c        write(6,*)
ccccccccccccccccccc
          enddo
        enddo
      enddo
C#ifdefC DEBUG
C      print "(' sumcheck sgpb', 2f18.6)", sum(sgpb)
C#endif

c---------------------------------------
c sgpp block------->removed
c---------------------------------------

c fouvb = 4*pi/(G^2-E) * integral(J*B)
      if(debug) print *,' mkjp_4: Four'
      fouvb=0d0
      do ig1 = 1,ngc
        do lm  = 1,nlx
          l = ll(lm)
          do n =1,nx(l)
            call gintxx(ajr(1,l,ig1), rprodx(1,n,l), a,b,nr,
     o                 radint )
            fouvb(ig1, n, lm) =
     &    fpi/(absqg(ig1)**2-eee) *dconjg(pjyl(lm,ig1))*radint !eee is supposed to be negative

          enddo
        enddo
      enddo
C#ifdefC DEBUG
C      print "(' sumcheck fouvb', 2f18.6)", sum(fouvb)
C#endif

c-----------------------------
c fouvp block --->removed
c-----------------------------

      deallocate(ajr,a1,   qg,absqg,   pjyl)
      if (allocated( cy )) deallocate(cy)
      if (allocated( yl )) deallocate(yl)
      end

      real(8) function fac2m(i)
cC A table of (2l-1)!!
c     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
      logical,save::  init=.true.
      real(8),save:: fac2mm(0:100)
      if(init) then
        fac2mm(0)=1d0
        do l=1,100
          fac2mm(l)=fac2mm(l-1)*(2*l-1)
        enddo
      endif
      fac2m=fac2mm(i)
      end
c=====================================================================
      subroutine genjh(eee,nr,a,b,lx, nrx,lxx,
     o      rofi,rkpr,rkmr)
C-- Generate radial mesh rofi, spherical bessel, and hankel functions
Cr  rkpr, rkmr are real fucntions --
ci   ee :  E <0 dumping wave
cr      rkpr = (2l+1)!! * j_l(i sqrt(abs(E)) r) * r / (i sqrt(abs(E)))**l
cr      rkmr = (2l-1)!! * h_l(i sqrt(abs(E)) r) * r * i*(i sqrt(abs(E)))**(l+1)
cr rkpr reduced to be r**l*r      at E \to 0
cr rkmr reduced to be r**(-l-1)*r at E \to 0
C-----------------------------------------------------------
      implicit none
      integer:: nr,lx, nrx,lxx,ir,l
      real(8):: a,b,eee,psi(0:lx),phi(0:lx)
      real(8):: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx),fac2m
      integer jobBes
      procedure(integer) :: GWversion

      jobBes = mod(GWversion()/1000,10)
      rofi(1)    = 0d0
      do ir      = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
      enddo
      if(eee==0d0) then
        do l = 0,lx
          rkpr(1:nr,l) = rofi(1:nr)**(l +1)
          rkmr(2:nr,l) = rofi(2:nr)**(-l-1 +1)
          rkmr(1,l)    = rkmr(2,l)
        enddo
      else
        do ir  = 1, nr
C         call bessl(eee*rofi(ir)**2,lx,phi(0:lx),psi(0:lx))
          call besslr(eee*rofi(ir)**2,10*jobBes,0,lx,phi(0:lx),psi(0:lx))
          do l = 0,lx    !fac2m(l)= (2l-1)!!
c            print *,' phi=',l,phi(l),phi(l)*fac2m(l+1)
c            print *,' psi=',l,psi(l),psi(l)/fac2m(l)
            rkpr(ir,l) = phi(l)* rofi(ir)**(l +1) *fac2m(l+1)
            if(ir/=1) rkmr(ir,l) = psi(l)* rofi(ir) **(-l ) /fac2m(l)
          enddo
        enddo
        rkmr(1,0:lx) = rkmr(2,0:lx)
      endif
      end

      subroutine mkjb_4(lxx,lx,nxx,nx,a,b,nr,nrx,rprodx,rofi,rkpr,rkmr,
     o  rojb,sgbb)
C- Make onsite integrals rojb and sgbb in one MT
C ----------------------------------------------------------------------
Ci Inputs
Ci   lxx   :dimensions nx,rprodx,rojb,sgbb
Ci   lx    :l-cutoff for rprodx,rkpr,rojbo
Ci   nxx   :dimensions rojb,sgbb
Ci   nx    :l-dependent number of radial functions in product basis
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nrx   :Leading dimension of rprodx
Ci   rprodx:Product basis functions * r: r*BRlμ(r)
Ci   rofi  :radial mesh points
Ci   rkpr  :r*Jl(E,r)
Ci   rkmr  :r*Hl(E,r)
Co Outputs
Co   rojb  :Onsite integral rhojb(BRlμ) = 1/(2l + 1)!! int dr r*Jl(E,r) r*BRlμ(r)
Co         :Needed to evaluate integrals <B^k_RLμ(r)|vk(r,r')|B^k_R'L'μ'(r')>
Co         :when (R,T) and (R',T') are different since 1-center exp. of v is prop Jl
Co   sgbb  :Onsite integrals <B|v(onsite)|B> corresponding to rojb when (R,T) = (R',T')
Co         :sgbb = sigma integral, given by Eq 43 in ecalj manual (June 19, 2015)
Co         :Follows from the 1-center expansion of the Yukawa potential
Co         :sgbb(B_Rlμ, B_Rlν) =
Co         :  4*pi/(2l+1) * int dr dr' (r<)Jl(E,r<) (r>)Hl(E,r>) (r B_Rlμ(r)) (r' B_Rlν(r'))
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr  Here the Jl and Hl are "bar Bessel and Hankel functions"
Cr  with the property Jl -> r^l  and Hl -> r^-l-1 for E=0.  They are
Cr  proportional to the radial part of customary sphereical Hankels and Bessels.
Cr
Cr  See Section 18.4 of Takao's ecalj manual (June 2015)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer :: lxx,lx,nxx,nx(0:lxx),nr,nrx
      real(8) :: rprodx(nrx,nxx,0:lxx),a,b
      real(8) :: rojb(nxx,0:lxx) ! rho-type onsite integral
      real(8) :: sgbb(nxx,nxx,0:lxx) ! sigma-type onsite integral
c internal
      integer :: l,n,n1,n2
      real(8) :: fac,fpi,pi,sig
      real(8) :: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx)
      pi  = 4d0*datan(1d0)
      fpi = 4*pi

c      real(8),allocatable :: rkpr(:,:),rkmr(:,:)
c
c      allocate(rkpr(nr,0:lx),rkmr(nr,0:lx))
c-------------------------------------------------
c rofi and aj = r**l / (2l+1)!! \times r. Sperical Bessel at e=0.
ccccccccccccccccccccccccccccccccc
c      do l = 0,lx
c      do n = 1,nx(l)
c      do n1 = 1,nx(l)
c        call gintxx(rprodx(1:nr,n,l), rprodx(1:nr,n1,l), a,b,nr,
c     o                 xxx )
c      write(6,*)' check rprodx =',l,n,n-n1,xxx
c      enddo
c      enddo
c      enddo
c      stop 'xxx'
cccccccccccccccccccccccccccccccc

c      rofi(1)    = 0d0
c      do ir      = 1, nr
c        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
c      enddo
c      do l = 0,lx
c        rkpr(1:nr,l) = rofi(1:nr)**(l +1)
c        rkmr(2:nr,l) = rofi(2:nr)**(-l-1) *rofi(2:nr)
c        rkmr(1,l)    = rkmr(2,l)
c      enddo

C --- initialize ---
      rojb = 0d0
      sgbb = 0d0
C ... rojb
      fac = 1d0
      do l = 0,lx
        fac = fac/(2*l+1)
        do  n = 1, nx(l)
          call gintxx(rkpr(1,l),rprodx(1,n,l),a,b,nr,rojb(n,l))
        enddo
        rojb(1:nx(l),l) = fac*rojb(1:nx(l),l)
      enddo
C ... sgbb
      do l  = 0,lx
        do n1 = 1,nx(l)
          do n2 = 1,nx(l)
            call sigint_4(rkpr(1,l),rkmr(1,l),lx,a,b,nr,rprodx(1,n1,l),rprodx(1,n2,l)
     &            , rofi,sig )
            sgbb(n1, n2, l)=sig/(2*l+1)*fpi
          enddo
        enddo
      enddo

c      write(6,*) ' rojbsum=', sum(rojb(:,:)),   sum(abs(rojb(:,:)))
c      write(6,*) ' sgbbsum=', sum(sgbb(:,:,:)), sum(abs(sgbb(:,:,:)))
ccccccccccccccccccccccccccccccccccccc
c      write(6,*)' sigint 1 1 0=',sgbb(1, 1, 0) !/(16d0*datan(1d0))
c      sgbb(1, 1, 0) =0d0
cccccccccccccccccccccccccccccccccccccccc
c      deallocate(rkpr,rkmr)
      end

      subroutine sigint_4(rkp,rkm,kmx,a,b,nr,phi1,phi2,rofi, sig)
C-
C ----------------------------------------------------------------------
Ci Inputs
Ci   rkp
Ci   rkm
Ci   kmx
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   phi1
Ci   phi2
Ci   rofi  :radial mesh points
Ci   sig
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   11 Dec 17
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
C ... Local parameters
      integer :: nr,kmx
      real(8):: a,b, a1(nr),a2(nr),b1(nr),rkp(nr),rkm(nr),
     &    int1x(nr),int2x(nr), phi1(nr), phi2(nr),rofi(nr),sig
      real(8),parameter:: fpi = 4d0*3.14159265358979323846d0
c
      a1(1) = 0d0;  a1(2:nr) = rkp(2:nr)
      a2(1) = 0d0;  a2(2:nr) = rkm(2:nr)
      b1(1:nr) = phi1(1:nr)
      call intn_smpxxx(a1,b1,int1x,a,b,rofi,nr,0)
      call intn_smpxxx(a2,b1,int2x,a,b,rofi,nr,0)
c
      a1(1) = 0d0; a1(2:nr) =
     &  rkm(2:nr) *( int1x(1)-int1x(2:nr) )+ rkp(2:nr) * int2x(2:nr)
      b1(1:nr) = phi2(1:nr)
      call gintxx(a1,b1,A,B,NR, sig )
      end

      subroutine intn_smpxxx(g1,g2,int,a,b,rofi,nr,lr0)
C- Integral of two wave functions, used in ppdf
C ----------------------------------------------------------------------
Ci Inputs
Ci   g1    :left function entering integral
Ci   g2    :right function entering integral
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   lr0   :not used now
Co Outputs
Co   int(r) = \int_(r)^(rmax) g1(r') g2(r') dr'
Cr Remarks
Cr  Simpson's rule, with higher rule for odd-numbered mesh points
Cu Updates
Cu   11 Dec 17
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, INTENT(IN) :: nr
      integer, INTENT(IN) :: lr0
      double precision, INTENT(IN)  :: g1(nr),g2(nr),a,b,rofi(nr)
      double precision, INTENT(OUT) :: int(nr)
C ... Local parameters
      integer ir
      real(8), parameter :: w1 =  5.0_8/12.0_8,
     &                      w2 =  8.0_8/12.0_8,
     &                      w3 = -1.0_8/12.0_8,
     &                      ooth = 1.0_8/3.0_8,
     &                      foth = 4.0_8/3.0_8

      if(mod(nr,2) == 0) call rx( ' INTN: nr should be odd for simpson integration rule')
c
      int(1)=0.0d0
      DO  10  IR = 3,NR,2
        int(ir)=int(ir-2)
     &             + ooth*G1(IR-2)*G2(IR-2)*( a*(b+rofi(ir-2)) )
     &             + foth*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )
     &             + ooth*G1(IR)*G2(IR)*( a*(b+rofi(ir)) )
   10 CONTINUE

c At the value for odd points, use the same interpolation above
      do 20 ir = 2,nr-1,2
        int(ir)=int(ir-1)
     &             + w1*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) )
     &             + w2*G1(IR)  *G2(IR)*  ( a*(b+rofi(ir)  ) )
     &             + w3*G1(IR+1)*G2(IR+1)*( a*(b+rofi(ir+1)) )
   20 continue
      do ir=1,nr
        int(ir)=int(nr)-int(ir)
      enddo
      END

c-------------------------------------------------------------------
C      subroutine sigintAn1( absqg, lx, rofi, nr,
C     o             a1int)
Cc a1int(r') = r' * \int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} *
Cc                j_l(absqg r)/absqg**l
C      implicit none
C      integer :: nr,l,ir,lx
C      real(8):: a1int(nr,0:lx), rofi(nr),absqg
C      real(8)::
C     &   ak(0:lx) ,aj(0:lx), dk(0:lx), dj(0:lx),
C     &   aknr(0:lx),ajnr(0:lx),dknr(0:lx),djnr(0:lx),
C     &   phi(0:lx),psi(0:lx)
Cc---
Cc      print *,' sigintAn1: absqg=',absqg
C      if(absqg<1d-10) then
Cc      if(absqg<1d-6) then !23jan2004 1d-10 ok?
C        stop "sigintAn1: absqg=0 is not supported yet. Improve here."
Cc This part for absqg=0 has not been checked yet!
Ccc       call bessl(0d0,lx,phi,psi)
Cc       call besslr(0d0,10*jobBes,0,lx,phi,psi)
Cc        do ir = 1,nr
Cc        do l  = 0,lx
Cc          a1int(ir,l) = .5d0* rofi(nr)**2     * rofi(ir)**l     * phi(l)
Cc     &                +(1d0/(2d0*l+3d0)-.5d0) * rofi(ir)**(l+2) * phi(l)
Cc       enddo
Cc        enddo
C      else
C        call  radkj(absqg**2, rofi(nr),lx,aknr,ajnr,dknr,djnr,100*job+0)
C        a1int(1,:) = 0d0
C        do ir = 2,nr
C          call radkj(absqg**2, rofi(ir),lx,ak,aj,dk,dj,100*job+0)
C          do l  = 0,lx
C            a1int(ir,l) = ( (2*l+1)* aj(l)
C     &     -((l+1)* ajnr(l)+ rofi(nr)*djnr(l) )*(rofi(ir)/rofi(nr))**l)
C     &    /absqg**2
C     &   *rofi(ir)
C          enddo
C        enddo
C      endif
Cc      print *,' sigintAn1: end'
C      end

      subroutine sigintAn1(absqg,lx,rofi,nr,a1int)
C- a1int(r') = r' * \int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} * j_l(absqg r)/absqg**l
C ----------------------------------------------------------------------
Ci Inputs
Ci   absqg :wave number
Ci   lx    :generate a1int for l=0:lx
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Co Outputs
Co   a1int :integral
Cr Remarks
Cr   Let:
Cr     aj(r,l) = j_l at r
Cr     dj(r,l) = dj_l/dr at r
Cr     R = rmax = rofi(nr)
Cr   Integral is:
Cr   I = r/q^2 [(2l+1)*aj(r,l) - (l+1)*aj(R,l)(r/R)^l - R*dj(R,l)(r/R)^l]
Cr   For small q (see radkj)
Cr     aj(r,l) -> r^l * [1/(2l+1)!! - ((q*r)^2/2)/(2*l+3)!!]
Cr     dj(r,l) -> r^(l-1)[l/(2l+1)!! - (l+2)*((q*r)^2/2)/(2*l+3)!!]
Cr   Expand brackets [...] in I as power series in q:
Cr     r^l/(2l+1)!! [(2l+1) - (l+1) - l] -
Cr     q^2/2 r^l/(2l+3)!! * [(2l+1)r^2/2 - (l+1)R^2/2 - (l+2)*R^2/2]
Cr   Term proportional to q^0 vanishes and integral becomes
Cr     I = r^(l+1)/(2l+1)!! [R^2/2 - (2l+1)r^2/2/(2l+3)]
Cu Updates
Cu   27 Apr 13 Correct treatment q=0 case.
C ----------------------------------------------------------------------
      implicit none
      integer :: nr,l,ir,lx
      real(8):: a1int(nr,0:lx), rofi(nr),absqg
      real(8)::
     &   ak(0:lx) ,aj(0:lx), dk(0:lx), dj(0:lx),
     &   aknr(0:lx),ajnr(0:lx),dknr(0:lx),djnr(0:lx)
      double precision dfac,aq0,rlp1,rmax2
      integer jobBes
      procedure(integer) :: GWversion

      jobBes = mod(GWversion()/1000,10)

C ... Small q case: if q<1d-5 the q->0 limit is correct to order absqg^2
C     Switch to q=0 case at point where loss in machine precision from
C     general formula exceeds higher order corrections to power series.
      if (absqg<5d-5) then

        do  ir = 1, nr
          dfac = 1
          rlp1 = 1
          rmax2 = rofi(nr)**2
          do  l = 0, lx
            rlp1 = rlp1*rofi(ir)
            dfac = dfac*(2*l+1)
            aq0 = rmax2/2  - (2*l+1)*((rofi(ir))**2/2/(2*l+3))
            a1int(ir,l) = aq0*rlp1/dfac
C           print *, ir,l, a1int(ir,l), aq0*rlp1/dfac
          enddo
        enddo

C ... General q case: accurate to 8 decimal places for q>5d-5
      else
        call radkj(absqg**2, rofi(nr),lx,aknr,ajnr,dknr,djnr,100*jobBes+0)
        a1int(1,:) = 0d0
        rmax2 = rofi(nr)**2
        do  ir = 1, nr
          if (rofi(ir) == 0d0) cycle
          call radkj(absqg**2, rofi(ir),lx,ak,aj,dk,dj,100*jobBes+0)
C         dfac = 1
C         rlp1 = 1
          do  l = 0, lx
C            rlp1 = rlp1*rofi(ir)
C            dfac = dfac*(2*l+1)
            a1int(ir,l) = ( (2*l+1)*aj(l) -
     .        ((l+1)*ajnr(l) + rofi(nr)*djnr(l))*(rofi(ir)/rofi(nr))**l)
     .        *rofi(ir)/absqg**2

C            aq0 = rmax2/2  - (2*l+1)*((rofi(ir))**2/2/(2*l+3))
C            print *, ir,l, a1int(ir,l), a1int(ir,l)- aq0*rlp1/dfac
C            print *, ir,l, a1int(ir,l), aq0*rlp1/dfac

          enddo
        enddo
      endif
      end

C      subroutine fmain
C
C      integer nr,lx
C      parameter (nr=5,lx=3)
C      integer ir,l
C      real(8):: a1int(nr,0:lx),rofi(nr),absqg
C
C
C      rofi(1) = 1d-3
C      do  ir = 2, nr
C        rofi(ir) = dble(ir)/2 - .1d0
C      enddo
C
C      absqg = .5d-6
C      read(*,*) absqg
C
C      call sigintAn1(absqg,lx,rofi,nr,a1int)
C      print "('#q=',1pe12.4)",absqg
C      do  ir = 1, nr
C        print 333, (a1int(ir,l), l=0,lx)
C  333   format(10f15.10)
C      enddo
C      end

c-------------------------------------------------
      subroutine sigintpp( absqg1, absqg2, lx, rmax,
     o             sig)
c sig(l)   =  \int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} *
c               j_l(absqg1 r)/absqg1**l
c               j_l(absqg2 r)/absqg2**l
c e1\ne0 e2\ne0
      implicit none
      integer :: l,lx,jobBes
      procedure(integer) :: GWversion
      real(8)::  rmax,sig(0:lx), absqg1,absqg2, e1,e2,
     &   ak1(0:lx) ,aj1(0:lx), dk1(0:lx), dj1(0:lx),
     &   ak2(0:lx) ,aj2(0:lx), dk2(0:lx), dj2(0:lx),
     &  fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx)
c---
      e1 = absqg1**2
      e2 = absqg2**2
      jobBes = mod(GWversion()/1000,10)
c
c      print *," sigintpp",e1,e2
c
      call wronkj(jobBes,e1,e2,rmax,lx,fkk,fkj,fjk,fjj)
      call  radkj( e1,    rmax,lx,   ak1,aj1,dk1,dj1,100*jobBes+0)
      call  radkj( e2,    rmax,lx,   ak2,aj2,dk2,dj2,100*jobBes+0)
c
      do l = 0,lx
        sig(l)= ( -l*(l+1)*rmax*aj1(l)*aj2(l)
     &            + rmax**3 * dj1(l)*dj2(l)
     &            + 0.5d0*rmax**2* (aj1(l)*dj2(l)+aj2(l)*dj1(l))
     &            - fjj(l)*(2*l+1)*(e1+e2)/2d0
     &           ) /(e1*e2)
      enddo
      end
