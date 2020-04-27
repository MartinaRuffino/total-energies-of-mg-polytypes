
        module vertex
      use h5

      implicit none

      contains

      subroutine read_vertex(mode,beta,nvfl,nomv,nOm,norb,infile,outfile)
!     Translation in fortran of Suscept.py vertex reader
!     mode : 0  read the header
!            1   read all the file and put it in chis,chic (not implemented yet)
!            2   read all the file and print it in hdf5 format
!            3   1+2

      real(8),intent(inout) :: beta
      integer,intent(in) ::  mode
      character(len=20),intent(in):: infile, outfile
      integer,intent(inout) :: norb, nvfl, nomv, nOm
!     local variables
      type(h5space) :: s_h5
      integer:: i,j,ib,iom,dom,sm0,em0,sm1,em1,im0,ti0,ti1,im1,nomm
      integer :: iv1,iv2
      integer,allocatable :: bfl_index(:),orb_ud(:,:)
      real(8) :: Om,om1
      real(8) :: vhr,vhi,vfr,vfi
      complex(8),allocatable :: vertexf(:,:,:,:,:),vertexh(:,:,:,:,:)
!      complex(8),intent(out) :: chis(2*nomv,2*nomv,norb,norb,nOm)
!      complex(8),intent(out) :: chic(2*nomv,2*nomv,norb,norb,nOm)
      complex(8),allocatable :: chis(:,:,:,:,:),chic(:,:,:,:,:)
      real(8),allocatable :: chi_tmp(:,:,:,:,:)
      complex(8),allocatable :: vf_uu(:,:),vh_uu(:,:)
      complex(8),allocatable :: vh_ud(:,:)
      complex(8),allocatable :: v_uu(:,:),v_ud(:,:),v_dd(:,:)


      open(file=infile,unit=123)
      rewind 123
      read(123,*) !first line
      read(123,*) beta, nvfl, nomv, nOm, nomm


      write(*,*)'beta=',beta
      write(*,*)'nvfl=',nvfl
      write(*,*)'nomv=',nomv
      write(*,*)'nOm= ',nOm

      allocate(bfl_index(nvfl))

      do ib =1,nvfl
         read(123,*)i,j
         bfl_index(ib)=j
      enddo
      norb = nvfl / 2


      allocate(orb_ud(norb,2))


      allocate(vf_uu(2*nomv,2*nomv))
      call dpzero(vf_uu,2*size(vf_uu))

      allocate(vh_uu(2*nomv,2*nomv))
      call dpzero(vh_uu,2*size(vh_uu))

      allocate(vh_ud(2*nomv,2*nomv))
      call dpzero(vh_ud,2*size(vh_ud))

      allocate(v_uu(2*nomv,2*nomv))
      call dpzero(v_uu,2*size(v_uu))

      allocate(v_dd(2*nomv,2*nomv))
      call dpzero(v_dd,2*size(v_dd))

      allocate(v_ud(2*nomv,2*nomv))
      call dpzero(v_ud,2*size(v_ud))

      do ib =1,norb
         orb_ud(ib,1)=ib
         orb_ud(ib,2)=ib+norb
      enddo

      ! do ib =1,nvfl
      !    if (ib <= norb) then
      !       orb_ud(bfl_index(ib)+1,1)=ib
      !    else
      !       orb_ud(bfl_index(ib)+1,2)=ib
      !    endif
      ! enddo

      write(*,*) 'the vertex is read'
      read(123,*)               ! # comment # b0 b1 Om om1
      allocate(vertexh(nvfl,nvfl,2*nom-1,2*nomv,2*nomv))
      allocate(vertexf(nvfl,nvfl,2*nom-1,2*nomv,2*nomv)) ! needs to be set to zero
      call dpzero(vertexh,2*size(vertexh))
      call dpzero(vertexf,2*size(vertexf))
      do i =1,Nvfl
         do j =1,Nvfl
            write(*,*)i,j
            do iom=0,2*nom-2
               dom=nom-iom-1
               sm0=max(0,-dom)
               em0=min(2*nomv,2*nomv-dOm)

               do im0=sm0+1,em0
                  read(123,*)ti0,ti1,Om,om1
                  sm1 = max(0,-dOm)
                  em1 = min(2*nomv,2*nomv-dOm)
                  do im1=sm1+1,em1
                     read(123,*)Om,vhr,vhi,vfr,vfi
                     vertexh(i,j,iom+1,im0,im1) =cmplx(vhr,vhi)
                     vertexf(i,j,iom+1,im0,im1) =cmplx(vfr,vfi)
                  enddo
               enddo
            enddo
         enddo
      enddo


      allocate(chi_tmp(nvfl,nvfl,2*nom-1,2*nomv,2*nomv))
      write(*,*) 'the vertex is stored in chiS_lm.h5 file'

      open(file='chiS_lm.h5',unit=123)
      close(123,status='delete')

      chi_tmp(:,:,:,:,:)=real(vertexh(:,:,:,:,:))
      call h5_write('chiS_lm.h5:/vh-r', chi_tmp, h5t_native_real8)
      chi_tmp(:,:,:,:,:)=aimag(vertexh(:,:,:,:,:))
      call h5_write('chiS_lm.h5:/vh-i', chi_tmp, h5t_native_real8)

      chi_tmp(:,:,:,:,:)=real(vertexf(:,:,:,:,:))
      call h5_write('chiS_lm.h5:/vf-r', chi_tmp, h5t_native_real8)
      chi_tmp(:,:,:,:,:)=aimag(vertexf(:,:,:,:,:))
      call h5_write('chiS_lm.h5:/vf-i', chi_tmp, h5t_native_real8)
      deallocate(chi_tmp)

      allocate(chis(nom,norb,norb,2*nomv,2*nomv))
      allocate(chic(nom,norb,norb,2*nomv,2*nomv))
      call dpzero(chis,2*size(chis))
      call dpzero(chic,2*size(chic))

      do iom=1,nom
         do iv1 = 1,2*nomv
            do iv2 = 1,2*nomv

               do i =1,norb
                  vf_uu(iv1,iv2) = 0.5*(vertexf(orb_ud(i,1),orb_ud(i,1),iom,iv1,iv2) &
                      + vertexf(orb_ud(i,2),orb_ud(i,2),iom,iv1,iv2))
                  do j =1,norb
                     vh_uu(iv1,iv2) = 0.5*(vertexh(orb_ud(i,1),orb_ud(j,1),iom,iv1,iv2) &
                         +vertexh(orb_ud(i,2),orb_ud(j,2),iom,iv1,iv2))
                     vh_ud(iv1,iv2) = 0.5*(vertexh(orb_ud(i,1),orb_ud(j,2),iom,iv1,iv2) &
                         +vertexh(orb_ud(i,2),orb_ud(j,1),iom,iv1,iv2))
                     v_ud(iv1,iv2) = vh_ud(iv1,iv2)
                     if ( i== j ) then
                        v_uu(iv1,iv2) = vh_uu(iv1,iv2) - vf_uu(iv1,iv2)

                     else
                        v_uu(iv1,iv2) = vh_uu(iv1,iv2)
                     endif

                     chis(nom+1-iom,i,j,iv1,iv2) = v_uu(iv1,iv2)-v_ud(iv1,iv2)
                     chic(nom+1-iom,i,j,iv1,iv2) = v_uu(iv1,iv2)+v_ud(iv1,iv2)

                  enddo
               enddo
            enddo
         enddo
      enddo

      allocate(chi_tmp(nom,norb,norb,2*nomv,2*nomv))
      close(123,status='delete')

      chi_tmp(:,:,:,:,:)=real(chis(:,:,:,:,:))
      call h5_write(outfile//':/chis-r', chi_tmp, h5t_native_real8)
      chi_tmp(:,:,:,:,:)=aimag(chis(:,:,:,:,:))
      call h5_write(outfile//':/chis-i', chi_tmp, h5t_native_real8)

      chi_tmp(:,:,:,:,:)=real(chic(:,:,:,:,:))
      call h5_write(outfile//':/chic-r', chi_tmp, h5t_native_real8)
      chi_tmp(:,:,:,:,:)=aimag(chic(:,:,:,:,:))
      call h5_write(outfile//':/chic-i', chi_tmp, h5t_native_real8)




      call h5_write(outfile//':/norb', [norb], find_h5dtype(norb))
      call h5_write(outfile//':/nOm', [nom], find_h5dtype(nom))
      call h5_write(outfile//':/nomv', [ nomv ], find_h5dtype(nomv))
      call h5_write(outfile//':/beta', [ beta ], find_h5dtype(beta))

      end subroutine read_vertex

      end module vertex
