      subroutine makechi0loc(mode,s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,
     .     nspx,nkfbz,ldcix,nsp,ncix,sigdc,ef,nomgv,nomgw,iqtab_in,nqtab_in,chi0f,chi0_locf)
Ci    nomgw : nb of frequencies for w
C     i    nomgv : nb of frequencies for v
!     mode 0 : computes chi0 for the full bz and stores chi0 and chi0loc in h5,
!     chi0f, chi0 , iqtab not touched
!     1 : computes chi0 for point in iqtab and stores chi0 and chi0loc in h5
!     chi0f, chi0  not touched
!     2 : computes chi0 for point in iqtab and stores it in  chi0f and chi0locf
!     iqtab   : integer array of size 3 x nqchi0 which contains reduced coordonate of q where chi0 is computed
!     nqtab   : size of iqtab
      use structures
C #ifdef H5
      use h5

      use mpi
      use meshes, only : khmesh
      use vertex, only : read_vertex

      implicit none
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_dmft)::  s_dmft
      type(h5space) :: s_h5

!     input
      integer, intent(in) :: nspx, nkfbz, ldcix, nsp, ncix, nomgv, nomgw, mode
      real(8),intent(in) :: ef
      
!     output
!     complex(8),intent(out) ::chi0(nkfbz,nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix),chi0_loc(nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix)
      complex(8),allocatable :: chi0(:,:,:,:,:,:,:),chi0_loc(:,:,:,:,:,:)
      complex(8) :: chi0f(2*nomgv,nomgw,ldcix,ldcix,ncix*ncix,nsp,nqtab_in),chi0_locf(nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix)
!     local variables
      integer :: iomg, nomg, nlo, nhi, iqs, i, j, isp
      integer ::iomgw, iomgv
      integer :: iqfbz, ifiproj, iq, ig, nevn, ifac(3), iv(3), nk1, nk2, nk3, nqtab,nqtab_in, iqtab_in(3,nqtab_in)      
      integer, allocatable :: i123(:,:), iqtab(:,:) 
      integer cix1, cix2
      integer :: il1, il2

      integer :: icix, ispc,ksp, nev, nl1, nspc

      integer,allocatable :: params
      real(8) :: qp(3), qpr(3), qb(9)
      real(8),allocatable ::  evlk(:,:)
      complex(8),allocatable :: sigbark(:,:,:), dmftu(:,:,:,:), gij(:,:,:), tmp(:,:)
      complex(8), allocatable :: gkloc(:,:,:,:,:,:,:), tmp1(:,:), gloc(:,:,:,:,:,:)
      complex(8) :: sigdc(s_dmft%ndsig)
      complex(8) :: omfac       ! 1 or sqrt(-1) depending on s_dmft%lsigim

      complex(8) :: old,u,ud
      complex(8),allocatable :: gtmp(:,:,:,:,:,:,:),gtmp_loc(:,:,:,:,:,:)
      complex(8),allocatable ::chi0_tmpz(:,:,:,:,:,:,:),chi0_tmp_locz(:,:,:,:,:,:)
      real(8),allocatable ::chi0_tmp(:,:,:,:,:,:),chi0_tmp_loc(:,:,:,:,:,:)
      real(8),allocatable:: gkloc_tmp(:,:,:,:,:,:,:)
      real(8),allocatable:: gloc_tmp(:,:,:,:,:,:)
      character(len=256) :: strn
!     qlist
      integer :: nq, nkread
      integer, allocatable ::qlist(:,:)
C     For MPI
      integer, allocatable :: kpproc(:) ! q loop blocking indices for multithreaded case
      integer, allocatable :: displs(:), recvcounts(:)
      integer :: rank
      integer procid,master,nproc,err
      real(8), parameter :: ry2eV = 13.60569193d0

      procedure(integer) :: nglob,fopna,fopng,fopnx,iprint,mpipid
      procedure(logical) :: cmdopt

C     setup mpi
      nproc = mpipid(0)
      procid = mpipid(1)



      nk1=s_bz%nkabc(1)
      nk2=s_bz%nkabc(2)
      nk3=s_bz%nkabc(3)

      nlo=s_dmft%nlohi(1)
      nhi=s_dmft%nlohi(2)
      nevn=nhi-nlo+1            ! number of bands


      nomg=s_dmft%nomg

      allocate(sigbark(nlo:nhi,nlo:nhi,nspx))
      allocate(gij(nlo:nhi,nlo:nhi,nspx))
      allocate(dmftu(nlo:nhi,ldcix,nsp,ncix))
      allocate(evlk(s_ham%ndham,nsp))
      allocate(tmp(nlo:nhi,nlo:nhi))
      allocate(gkloc(ldcix,ldcix,nsp,ncix,ncix,nomg,nkfbz))
      allocate(gloc(ldcix,ldcix,nsp,ncix,ncix,nomg))






      omfac = 1 ; if (s_dmft%lsigim) omfac = (0d0,1d0)

      call dpzero(sigbark,2*size(sigbark))
      call bzmsh00(s_lat%plat,s_bz%lshft,0,s_bz%nkabc,ifac,qb) ! Needed for iqstar
      ifiproj  = fopna('proj',-1,4); rewind ifiproj
      nspc = 1; if (nsp > nspx) nspc = 2
      allocate(tmp1(ldcix,nevn))
      call dpzero(gij,2*size(gij))
      call dpzero(tmp1,2*size(tmp1))
      call dpzero(gkloc,2*size(gkloc))
      call dpzero(gloc,2*size(gloc))


      allocate(i123(3,nkfbz))

!...  (UgU^+)_(RL)(R'L')
      do  iqfbz = 1, nkfbz
         

         call iodmftu(s_dmft,.false.,s_ham%ndham,nlo,nhi,ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifiproj)
         

         qpr = qp
         if (ig == 1) iqs = 0   ! First point of new star
!     lqstar is assumepd to be true
         call iqstar(2,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,iv)
         
         i123(:,iqfbz) = iv(:)
         

         call dpzero(sigbark,2*size(sigbark))
         
         do iomg = 1,  nomgv+nomgw
!     make g_ij
            call embed_sigma(s_dmft,nlo,nhi,ldcix,nsp,nspx,dmftu,sigdc,s_dmft%sig(iomg,:),sigbark(:,:,:))
            do  isp = 1, nspx

!     --- Full Green's function (eigenfunction basis) and eigenvalues in energy subspace ---
!     One-body part of g is diagonal
!     [g_ij]^-1 = delta_ij [omega + ef - evl] - sigbark_ij
               do  j = nlo, nhi
                  do  i = nlo, nhi
                     tmp(i,j) = -sigbark(i,j,isp)
!     add omega+mu-eval (and broadening) on the diagonal
                     if (i==j) then
                        tmp(i,j) = tmp(i,j) + s_dmft%omg(iomg)*omfac + ef - evlk(i,isp)
                        if (.not. s_dmft%lsigim) call rx('not ready for real axis')
                     endif
                  enddo
               enddo
!     perform inversion on (nlo:nhi) subblock (It's not efficient at the 1st iteration)

               call zinv(nevn,tmp,nevn)
               do  j = nlo, nhi
                  do  i = nlo, nhi
                     gij(i,j,isp)=tmp(i,j)
                  enddo         !i
                  if (.not. s_dmft%lsigim) gij(j,j,isp) = gij(j,j,isp) + (0d0,1d0)*s_dmft%gammac
               enddo            !j


!     .....   Projection
               do  cix1 = 1, ncix
                  do  cix2 = 1, ncix

!     In noncollinear case, isp=1 always => need internal ispc=1..2
!     ksp is the current spin index in both cases:
!     ksp = isp  in the collinear case
!     = ispc in the noncollinear case
!     ispc=1 for independent spins, and spin index when nspc=2
                     do  ispc = 1, nspc
                        ksp = min(max(ispc,isp),nsp)

                        icix = iabs(s_dmft%icix(cix1))
                        nl1 = 2*s_dmft%l(icix)+1
                        if (nl1 > ldcix) call rx('bug in sudmft')

                        call zgemm('C','N',nl1,nevn,nevn,(1d0,0d0),
     .                       dmftu(nlo,1,ksp,cix1),nevn,
     .                       gij(nlo,nlo,isp),nevn,(0d0,0d0),
     .                       tmp1,ldcix)
                        call zgemm('N','N',nl1,nl1,nevn,(1d0,0d0),
     .                       tmp1,ldcix,
     .                       dmftu(nlo,1,ksp,cix2),nevn,(0d0,0d0),
     .                       gkloc(1,1,ksp,cix1,cix2,iomg,iqfbz),ldcix)
!     call daxpy(2*ldcix*ldcix,1d0,gkloc(1,1,ksp,cix1,cix2,iomg,iqfbz),1,gloc(1,1,ksp,cix1,cix2,iomg),1)
                        gloc(:,:,ksp,cix1,cix2,iomg)=gloc(:,:,ksp,cix1,cix2,iomg)+gkloc(:,:,ksp,cix1,cix2,iomg,iqfbz)/nkfbz
                     enddo      ! ispc
                  enddo         ! cix1
               enddo            ! cix2
            enddo               !isp
         enddo                  !iomg
      enddo                     !iqfbz
      call fclose(ifiproj)



!     Q index has to be the fast index for the convolution
      allocate(gtmp(nkfbz,nomg,nsp,ldcix,ldcix,ncix,ncix))
      call dpzero(gtmp,2*size(gtmp))
      allocate(gtmp_loc(nomg,nsp,ldcix,ldcix,ncix,ncix))
      call dpzero(gtmp_loc,2*size(gtmp_loc))
      do iq =1,nkfbz
         do iomg=1,nomg
            do isp =1, nsp
               do il1 =1,ldcix
                  do il2=1,ldcix
                     do cix1=1,ncix
                        do cix2=1,ncix
                           gtmp(iq,iomg,isp,il1,il2,cix1,cix2)=gkloc(il1,il2,isp,cix1,cix2,iomg,iq)
                           if(iq==1) then
                              gtmp_loc(iomg,isp,il1,il2,cix1,cix2)=gloc(il1,il2,isp,cix1,cix2,iomg)
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo                     !iqfbz
      
      if (mode == 0) then         
         nq = nkfbz
         allocate(iqtab(3,nq))
         iqtab(:,:) = i123(:,:)
      else
         nq = nqtab_in
         allocate(iqtab(3,nq))
         iqtab(:,:) = iqtab_in(:,:)
      endif
      
!     chi0(iomv) , v=(2*(iomv-1-nomv)+1)*pi/beta
!     if(iomv<=nomgv)  G(iomv)=conj(gtmp(nomv-iomv))
!     if(iomv>nomgv)  G(iomv)=gtmp(iomv+nomv)

      allocate(chi0(nq,nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix))

      allocate(chi0_loc(nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix))

      call dpzero(chi0,2*size(chi0))
      call dpzero(chi0_loc,2*size(chi0_loc))
      do cix1=1,ncix
         do cix2=1,ncix
            write(*,*)'chi0 for cix :',cix1,cix2,ncix*(cix2-1)+cix1
            call makechi0k(gtmp(1,1,1,1,1,cix1,cix2),gtmp(1,1,1,1,1,cix1,cix2),nkfbz,nk1,nk2,nk3,
     .           nomg,nsp,ldcix,nomgw,nomgv,i123,iqtab,nq,chi0(1,1,1,1,1,1,ncix*(cix2-1)+cix1))

            call  makechi0_loc( gtmp_loc(1,1,1,1,cix1,cix2),gtmp_loc(1,1,1,1,cix1,cix2),
     .           nomg,nsp,ldcix,nomgw,nomgv,chi0_loc(1,1,1,1,1,ncix*(cix2-1)+cix1))
            write(*,*)'chi0 for cix :',cix1,cix2
         enddo
      enddo


      if( mode < 2) then
         allocate(chi0_tmp_loc(nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix))      
         open(file='chi0.h5',unit=123)
         close(123,status='delete')
         write(*,*)'print chi0'

         allocate(chi0_tmp(nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix))
         call s_h5 % init( shape(chi0_tmp))
         do  iqfbz = 1, nq
            call awrit1('%x_%i',strn,len(strn),0,iqfbz)

            forall (iomgw=1:nomgw,iomgv=1:2*nomgv,isp=1:nsp,il1=1:ldcix,il2=1:ldcix,cix1=1:ncix*ncix)
     .           chi0_tmp(iomgw,iomgv,isp,il1,il2,cix1)=real(chi0(iqfbz,iomgw,iomgv,isp,il1,il2,cix1))/ry2eV/ry2eV

            call h5_write('chi0.h5:/chi0-r'//strn, chi0_tmp, find_h5dtype(chi0_tmp(1,1,1,1,1,1)))


            forall (iomgw=1:nomgw,iomgv=1:2*nomgv,isp=1:nsp,il1=1:ldcix,il2=1:ldcix,cix1=1:ncix*ncix)
     .           chi0_tmp(iomgw,iomgv,isp,il1,il2,cix1)=aimag(chi0(iqfbz,iomgw,iomgv,isp,il1,il2,cix1))/ry2eV/ry2eV
            call h5_write('chi0.h5:/chi0-i'//strn, chi0_tmp, find_h5dtype(chi0_tmp(1,1,1,1,1,1)))

         end do


         call h5_write('chi0.h5:/iqtab', iqtab,  find_h5dtype(iqtab(1,1)))



         call h5_write('chi0.h5:/beta', [s_dmft%beta], find_h5dtype(s_dmft%beta))



         write(*,*)'print chi0_loc'
         chi0_tmp_loc(:,:,:,:,:,:)=real(chi0_loc(:,:,:,:,:,:))/ry2eV/ry2eV

         call h5_write('chi0.h5:/chi0loc-r', chi0_tmp_loc, find_h5dtype(chi0_tmp_loc(1,1,1,1,1,1)))
         chi0_tmp_loc(:,:,:,:,:,:)=aimag(chi0_loc(:,:,:,:,:,:))/ry2eV/ry2eV
         call h5_write('chi0.h5:/chi0loc-i', chi0_tmp_loc, find_h5dtype(chi0_tmp_loc(1,1,1,1,1,1)))
      else

!         allocate(chi0_tmpz(nomgw,2*nomgv,nsp,ldcix,ldcix,ncix*ncix,nq))         

         forall (iomgw=1:nomgw,iomgv=1:2*nomgv,isp=1:nsp,il1=1:ldcix,il2=1:ldcix,cix1=1:ncix*ncix,iq=1:nq)
!     .        chi0_tmpz( iomgw, iomgv, isp, il1, il2, cix1,iq) = chi0(iqfbz, iomgw, iomgv, isp, il1, il2, cix1)
     .        chi0f(  iomgv,iomgw,  il1, il2, cix1,isp,iq) = chi0(iqfbz, iomgw, iomgv, isp, il1, il2, cix1)         
!         chi0f(:,:,:,:,:,:,:) = chi0_tmpz(:,:,:,:,:,:,:)

      endif



      deallocate(chi0)

      end subroutine






      subroutine makechi0k(gij1,gij2,nkfbz,nk1,nk2,nk3,
     .     nomg,nspx,dimij,nomgw,nomgv,i123,qlist,nq,chi0)
!     chi0(q,w,v)=\sum_k G(k+q,v) G(k+q,v+w)
!     i    nkfbz : nb of k-points in the full BZ
!     i    nk1,nk2,nk3 : division of the full BZ
!     i    nomg : nb of frequencis of gij1 and gij2
!     i     gij(nkfbz,nomg,nspx,i,j)
!     i      i123  : (mode>1) indices i1,i2,i3, to microcell in full BZ (bzmesh)
!     i      qlist  :  selection of i123 made by readind qlist
!     o      chi0(nkfbz,nomgw,2*nomgv,nspx,dimij,dimji)
!     o      Fermionic grid v=2*pi*(n-nomgv-1)/beta n=1..2*nomv

      use structures
      implicit none
      integer ,intent(in) :: nkfbz,nk1,nk2,nk3,nomg,nspx,dimij,nomgw,nomgv,i123(3,nkfbz),nq,qlist(3,nq)
      complex(8),intent(in) :: gij1(nkfbz,nomg,nspx,dimij,dimij)
      complex(8),intent(in) :: gij2(nkfbz,nomg,nspx,dimij,dimij)
      complex(8),intent(out) :: chi0(nq,nomgw,2*nomgv,nspx,dimij,dimij)
C     local variables
      integer :: k1,k2,k3,iomgw,iomgv,ind1,ind2
      integer ::i,j,i1,i2,i3,iqfbz,isp,iq, ik1,ik2,ik3
      complex(8) :: tmp1(nk1,nk2,nk3),tmp2(nk1,nk2,nk3)

      call fftz30(nk1,nk2,nk3,k1,k2,k3)
      if(nomg <nomgw +nomgv-1) then
         call rx('nomg <nomgw +nomgv-1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      endif

      do  j = 1, dimij
         do  i = 1, dimij
            do iomgw = 0, nomgw-1 !Bosonic Grid
               do iomgv = 1, 2*nomgv
                  ind1=iomgv-nomgv
                  ind2=iomgv-nomgv+iomgw

                  do  isp = 1, nspx
                     do  iqfbz = 1, nkfbz
                        ik1 = i123(1,iqfbz)
                        ik2 = i123(2,iqfbz)
                        ik3 = i123(3,iqfbz)
                        
                        if(ind1 >=1) then
                           tmp1(ik1,ik2,ik3) = gij1(iqfbz,ind1,isp,i,j)
                        else
                           tmp1(ik1,ik2,ik3) = conjg(gij1(iqfbz,-ind1+1,isp,i,j))
                        endif

                        if(ind2 >=1) then
                           tmp2(ik1,ik2,ik3) = gij2(iqfbz,ind2,isp,j,i)
                        else
                           tmp2(ik1,ik2,ik3) = conjg(gij2(iqfbz,-ind2+1,isp,j,i))
                        endif

                     enddo

                     call fftz3c(tmp1,tmp2,nk1,nk2,nk3,k1,k2,k3,000,0-1)

                     do  iq = 1, nq
                        chi0(iq,iomgw+1,iomgv,isp,i,j)=tmp1(qlist(1,iq),qlist(2,iq),qlist(3,iq))
                     enddo
                  enddo         !isp
               enddo
            enddo               !iomg
         enddo                  !i
      enddo                     !j
      end


      subroutine makechi0_loc(gij1,gij2, nomg,nspx,dimij,nomgw,nomgv,chi0)
C     i     g(nomg,nspx,i,j)
C     o      chi0(nomgw,nomgv,nspx,dimij,dimji)
      use structures
      implicit none
      integer ,intent(in) :: nomg,nspx,dimij,nomgw,nomgv
      complex(8),intent(in) :: gij1(nomg,nspx,dimij,dimij)
      complex(8),intent(in) :: gij2(nomg,nspx,dimij,dimij)
      complex(8),intent(out) :: chi0(nomgw,2*nomgv,nspx,dimij,dimij)
C     local variables
      integer :: iomgw,iomgv,ind1,ind2
      integer ::i,j,i1,i2,i3,isp
      complex(8) :: tmp1,tmp2
      if(nomg <nomgw +nomgv-1) then
         call rx('nomg <nomgw +nomgv-1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      endif

      do  j = 1, dimij
         do  i = 1, dimij
            do iomgw = 0, nomgw-1 !Bosonic Grid
               do iomgv = 1, 2*nomgv
                  ind1=iomgv-nomgv
                  ind2=iomgv-nomgv+iomgw
                  do  isp = 1, nspx
                     if(ind1 >=1) then
                        tmp1= gij1(ind1,isp,i,j)
                     else
                        tmp1= conjg(gij1(-ind1+1,isp,i,j))
                     endif

                     if(ind2 >=1) then
                        tmp2= gij2(ind2,isp,j,i)
                     else
                        tmp2= conjg(gij2(-ind2+1,isp,j,i))
                     endif
                     chi0(iomgw+1,iomgv,isp,i,j)=tmp1*tmp2
                  enddo         !isp
               enddo            !iomgv
            enddo               !iomg
         enddo                  !i
      enddo                     !j
      end







      subroutine makechikh(s_bz,s_ham,s_lat,s_pot,s_spec,s_site,s_dmft,
     .     gkloc,ldcix,nsp,ncix,nkfbz,nomg,i123,nomgw,nomgv)
C     - compute chi0 from gkloc as in suscept.py
      use structures
C     #ifdef H5
      use h5
C     #endif
      implicit none
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_dmft)::  s_dmft
      integer,intent(in) :: ldcix,nsp,ncix,nkfbz,nomg,i123(nkfbz,3)
      integer,intent(in) :: nomgw,nomgv
      complex(8),intent(in) :: gkloc(ldcix,ldcix,nsp,ncix,nkfbz,nomg)
C     local variables
      complex(8) :: tmp(nkfbz,nomg,nsp,ldcix,ldcix,ncix)
      complex(8),allocatable ::chi0(:,:,:,:,:,:,:)
      integer :: iq,iomg,il1,il2,isp,icix,nk1,nk2,nk3
      nk1=s_bz%nkabc(1)
      nk2=s_bz%nkabc(2)
      nk3=s_bz%nkabc(3)
      do iq =1,nkfbz
         do iomg=1,nomg
            do isp =1, nsp
               do il1 =1,ldcix
                  do il2=1,ldcix
                     do icix=1,ncix
                        tmp(iq,iomg,isp,il1,il2,icix)=gkloc(il1,il2,isp,icix,iq,iomg)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo                     !iqfbz
      allocate(chi0(nkfbz,nomgw,nomgv,nsp,ldcix,ldcix,ncix))
      call dpzero(chi0,2*size(chi0))
      do icix=1,ncix
         call makechi0k(tmp(1,1,1,1,1,icix),tmp(1,1,1,1,1,icix),nkfbz,nk1,nk2,nk3,
     .        nomg,nsp,ldcix,nomgw,nomgv,i123,chi0(1,1,1,1,1,1,icix))
         write(*,*)'chi0 for cix :',icix
      enddo
      write(*,*)'print chi0 in h5 !'
C     #ifdef H5

      call h5_write('chilockh.h5:/chi0', chi0, find_h5dtype(chi0(1,1,1,1,1,1,1)))
      call h5_write('chilockh.h5:/iqtab', i123,  find_h5dtype(i123(1,1)))

      end subroutine



      subroutine readqlist(mode,nq,nk,qlist)
!     read file qlist with the following format :
!     size of the nk mesh in the full BZ
!     q1(1) q1(2) q1(3)
!     q2(1) q2(2) q2(3)
!     ....
!     where q(1) =1..     size of the nk mesh in the full BZ
!     mode : 0 : read nk and the nb of q point in the List
!     1 : read the qlist(nq,3)
      implicit none
      integer,intent(in) ::mode
      integer,intent(inout) :: nk,nq
      integer,intent(out) ::qlist(3,nq)
!     local variables
      integer ::x,y,z,i

      open(file='qlist',unit=123)
      rewind 123
      read(123,*) x
      if (x /= nk) then
         call rx('size of the k-mesh in the file not compatible')
      endif

      if( mode==0) then
         nq=0
         do
            read(123,*,END=999)x,y,z
            nq=nq+1
         enddo
 999     continue
      else if (mode == 1) then
         do i=1, nq
            read(123,*)x,y,z
            qlist(1,i) = x + 1 
            qlist(2,i) = y + 1
            qlist(3,i) = z + 1
         enddo

         if (minval(qlist) <1) then
            call rx('qlist : min value of coordinates is  0, but should be 1 in fortran')
         else if (maxval(qlist) >nk) then
            call rx('qlist : max value of coordinates higher than the size of the mesh...')
         endif

      else
         call rx('readqlist: mode unknow')
      endif
      

      end subroutine
