C#define H5
      subroutine printgprt(s_dmft,s_bz,s_ham,s_lat,mode,ldcix,nkp,nspx,nsp,ncix,nkfbz,nomg,gkloc,gloc)
C     - This subroutine  prints gkloc and gloc  for gpmode=3,5
Ci      mode : 1 print diagonal gkloc and gloc
Ci             2 print gloc and gkloc as a matrix flattened
Ci             3 print g as a matrix flattened,skip orbitals without sig_dmft
C     gkloc and gloc also save in h5file
      use structures
C#ifdef H5
      use h5
C#endif
      implicit none
      type(str_dmft)      ::  s_dmft
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat

      integer,intent(in) :: mode,ldcix,nsp,ncix,nkfbz,nomg,nkp,nspx
      complex(8),intent(in) :: gkloc(ldcix,ldcix,nsp,ncix,nkfbz,nomg)
      complex(8),intent(in) :: gloc(ldcix,ldcix,nsp,ncix,nomg)


C     local variables
      integer ::cix_ind,il1,iqfbz,iomg,ifig,ifigk,dim_cix,dim_red
      integer ::nicixi(ncix),icix_ind
      procedure(integer) :: fopna
      procedure(logical) :: cmdopt
      character(len=256) :: strn
      real(8), parameter :: ry2eV = 13.60569193d0

C     info k_mesh
      integer ::nk1,nk2,nk3,ifac(3),ifiproj,iq,ig,iqs,iv(3)
      real(8),allocatable ::  evlk(:,:)
      complex(8),allocatable ::dmftu(:,:,:,:)
      real(8) :: qp(3),qpr(3),qb(9)
      allocate(dmftu(s_dmft%nlohi(1):s_dmft%nlohi(2),ldcix,nsp,s_dmft%ncix))
      allocate(evlk(s_ham%ndham,nsp))


      nk1=s_bz%nkabc(1)
      nk2=s_bz%nkabc(2)
      nk3=s_bz%nkabc(3)
      call bzmsh00(s_lat%plat,s_bz%lshft,0,s_bz%nkabc,ifac,qb) ! Needed for iqstar
      call ineqcix(s_dmft,ncix,nicixi)
      ifiproj = fopna('proj',-1,4); rewind ifiproj
      do cix_ind =1,ncix
       dim_cix = 2*s_dmft%l(iabs(s_dmft%icix(cix_ind)))+1 ! if the cix blocks have different sizes




         !if there is just 1 ind cix, the name of the file does not have cix number
         !This has to be removed because could be confusing. I put it for backward compatilibility
         if(s_dmft%nicix >1) then
            call awrit1('%x_%i',strn,len(strn),0,cix_ind)
            ifigk=fopna('gkloc'//strn,-1,2)
            ifig=fopna('gloc'//strn,-1,2)
         else
            ifigk=fopna('gkloc',-1,2)
            ifig=fopna('gloc',-1,2)
         endif
         rewind ifigk
         rewind ifig

C.....print gkloc
         if(mode==2) then
            write(ifigk,'("#",5(x,i0),a,/,a)') nkfbz, 1, nomg, dim_cix, 1
     &           ,"  # nkpt, nsymop, nom, cixdms norbitals"
     &           ,"#  0.00000  0.00000  0.00000    # actual position of correlated atoms in the unit cell" ! not really actual...
         else if (mode==3) then
            dim_red=0
            do il1=1,dim_cix
               icix_ind=abs(nicixi(cix_ind))
               if(s_dmft%sigind(il1,il1,icix_ind,1) /= 0)  dim_red=dim_red+1
            enddo
            write(*,*)dim_red
            write(ifigk,'("#",5(x,i0),a,/,a)') nkfbz, 1, nomg, dim_red, 1
     &           ,"  # nkpt, nsymop, nom, cixdms norbitals"
     &           ,"#  0.00000  0.00000  0.00000    # actual position of correlated atoms in the unit cell" ! not really actual...

         endif
         rewind ifiproj
         do iqfbz=1,nkfbz
            call iodmftu(s_dmft,.true.,s_ham%ndham,s_dmft%nlohi(1),s_dmft%nlohi(2),ldcix,
     .       nsp,s_dmft%ncix,iq,ig,qp,evlk,dmftu,ifiproj)
            qpr = qp
            if (ig == 1) iqs = 0 ! First point of new star
            call iqstar(2,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,iv)
            if(mode==1) then
               call awrit6('%% rows %i cols %i  iq %,5i  irr %,5i  i123 %3,4i  qp %3:2,6;6d',
     .              ' ',128,ifigk,nomg,1+2*nsp*(2*s_dmft%l(1)+1),iqfbz,iq,iv,qpr)

            endif
            do iomg=1,nomg
               write(ifigk,'(f14.8)',advance='no') s_dmft%omg(iomg)*ry2eV
               dim_cix = 2*s_dmft%l(iabs(s_dmft%icix(cix_ind)))+1 ! if the cix blocks have different sizes
               call printg(mode,s_dmft,ldcix,nkp,nspx,nsp,dim_cix,cix_ind,ifigk,gkloc(1,1,1,cix_ind,iqfbz,iomg))
               write(ifigk,*)
            enddo               ! iomg
         enddo                  ! iqfbz

         do iomg=1,nomg
            write(ifig,'(f14.8)',advance='no') s_dmft%omg(iomg)*ry2eV
            call printg(mode,s_dmft,ldcix,nkp,nspx,nsp,dim_cix,cix_ind,ifig,gloc(1,1,1,cix_ind,iomg))
            write(ifig,*)
         enddo                  ! iomg
         call fclose(ifig)
         call fclose(ifigk)
      enddo                     !ncix


C#ifdef H5
      open(file='gloc.h5',unit=123)
      close(123,status='delete')
      call h5_write('gloc.h5:/gloc', gloc, find_h5dtype(gloc))

      open(file='gkloc.h5',unit=123)
      close(123,status='delete')
      call h5_write('gkloc.h5:/gkloc', gkloc, find_h5dtype(gkloc))
C#endif

      end subroutine

      subroutine printg(mode,s_dmft,ldcix,nkp,nspx,nsp,dim_cix,cix_ind,ifi,g)
C     - print in 1 line g(l1,l2,nsp) depending of the mode of gprt
Ci      mode : 1 print diagonal of g
Ci             2 print g as a matrix flattened
Ci             3 print g as a matrix flattened,skip orbitals without sig_dmft
      use structures
      implicit none
      type(str_dmft)      ::  s_dmft
      complex(8),intent(in) :: g(ldcix,ldcix,nsp)
      integer,intent(in) ::mode,ldcix,nkp,nspx,nsp,dim_cix,cix_ind,ifi
      integer ::nicixi(s_dmft%ncix)
C local variables
      real(8), parameter :: ry2eV = 13.60569193d0
      integer :: isp,ispc,i,j,ksp,nspc,icix_ind
      call ineqcix(s_dmft,s_dmft%ncix,nicixi)
      icix_ind=abs(nicixi(cix_ind))
      nspc = 1; if (nsp > nspx) nspc = 2
      do  isp = 1, nspx
         do  ispc = 1, nspc
            ksp = min(max(ispc,isp),nsp)


C............mode 1
            if(mode==1) then
               do i = 1,ldcix
                  write(ifi,'(2(x,f14.8))',advance='no') g(i,i,ksp)/ry2eV
               enddo
C............mode 2
            else if (mode ==2) then
               if(nspx>1 ) then
                  call rx('printg : mode 2 not ready for non collinear case')
               endif
               do i = 1,ldcix
                  if (ksp == 2) then
                     do j = 1, dim_cix
                        write(ifi,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                     end do
                  end if
                  do j = 1, dim_cix
                     write(ifi,'(2(x,f14.8))',advance='no') g(j,i,ksp)/ry2eV
                  end do
                  if (ksp < nsp) then
                     do j = 1, dim_cix
                        write(ifi,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                     end do
                  end if
               enddo            !i
C............mode  3
            else if (mode == 3) then
               if(nspx>1 ) then
                  call rx('printg : mode 2 not ready for non collinear case')
               endif
               do i = 1,ldcix
                  if(s_dmft%sigind(i,i,icix_ind,1)==0)  cycle
                  if (ksp == 2) then
                     do j = 1, dim_cix
                        if(s_dmft%sigind(j,j,icix_ind,ksp)==0) cycle
                        write(ifi,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                     end do
                  end if
                  do j = 1, dim_cix
                     if(s_dmft%sigind(j,j,icix_ind,ksp)==0) cycle
                     write(ifi,'(2(x,f14.8))',advance='no') g(j,i,ksp)/ry2eV
                  end do
                  if (ksp < nsp) then
                     do j = 1, dim_cix
                        if(s_dmft%sigind(j,j,icix_ind,ksp)==0) cycle
                        write(ifi,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                     end do
                  end if
               enddo            !i
            else
               call rx('printg : unknown mode')
            endif
         enddo                  !ispc
      enddo                     !isp

      end subroutine
