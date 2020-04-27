      subroutine gradv1(nlml, z, nr, nsp, nspc, rofi, v1, grv1)
C- Calculates the gradient of the potential v1
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlml  :L-cutoff for density expansion
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rofi  :radial mesh points for density and potential
Ci   v1    :true potential at this site, excluding nuclear contribution
Ci         :but including background contribution
Co Outputs
Co   grv1  :
Cl Local variables
Cl   nn    :number of points used to differentiate radial f
Cl   np    :number of points for angular integration
Cl   ll    :function that takes nlml as input and returns  lmxl
Cl          lmxl is the l-cutoff for density, potential on the radial
Cl   r2    : created by ropyln, where defined as rsq(i) square of length
Cof point i
Cmesh
Cr Remarks
Cu   24 Aug 15 (Scott McKecknie) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr, nsp, nspc, nlml
      double precision rofi(nr), v1(nr,nlml,nsp*nspc), z, grv1(nr,3)
C ... Dynamically allocated local arrays
      real(8),allocatable :: yl(:,:),p(:,:),p2(:,:),r2(:),wp(:),
     .  gyl(:,:),angrv1(:,:,:),vavg(:,:),testggp(:,:)
C ... Local parameters
      integer nn, np, i, lmin, ll, ilm, stdo, nglob, lx, lmax1
      double precision fpi, srfpi, pi

C ... Setup
      stdo = nglob('stdo')
      fpi = 16*datan(1d0)
      pi = 4*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      lmin=1
      nn = 6
      np = 122
      lx=10     ! 10 means extrapolate for rofi(1) since zero
                ! the zero part means don't scale by r**2
      lmax1 = ll(nlml)+1
      allocate(yl(np,(lmax1+1)**2),wp(np))
      allocate(p(3,np),p2(np,3),r2(np))
      i = -np
      call fpiint(i,0,np,p,wp)  ! np
      call sanrg(.true.,-i,np,np,'locpot','np') ! Check int value
      call dmcpy(p,1,3,p2,np,1,np,3) ! matrix copy
      call ropyln(np,p2,p2(1,2),p2(1,3),lmax1,np,yl,r2)

C ... Gradient of spherical harmonics
      allocate (gyl(np*nlml,3))
      call ropylg(1,ll(nlml),nlml,np,np,p2,p2(1,2),p2(1,3),r2,yl,gyl)

C ... Average over spin and add in nuclear term for l=0 (nlml=1)
      allocate(vavg(nr,nlml))
      vavg(1:nr,1:nlml) = (v1(1:nr,1:nlml,1) + v1(1:nr,1:nlml,nsp))/2
      vavg(2:nr,1) = vavg(2:nr,1) - 2*z*srfpi/rofi(2:nr)

      write(stdo,'(a)') 'gradv1: Calculating gradient of v1...'
      do i=1,5
        !write(stdo,'(f20.8)') 2*z/rofi(i)
        do ilm=1,5
          write(stdo,'(x,f20.8)',advance='no') vavg(i,ilm)
        enddo
        write(stdo,'()',advance='yes')
      enddo

C ... Test gradfl (see notes on gradfl for this, should be similar

C ... Get gradient of the average potential with angular dependence
      allocate(angrv1(nr,np,3),testggp(nr,np))
      call gradfl(ll(nlml),nlml,nr,np,1,nr,0,lx,nn,rofi,yl,gyl,
     .    vavg,angrv1,testggp)
      deallocate(testggp)

!      write(stdo,'(a)') 'angrv1'
!      do i=1,10
!        do ilm=1,5
!          write(stdo,'(x,f20.8)',advance='no') angrv1(i,ilm,1)
!        enddo
!        write(stdo,'()',advance='yes')
!      enddo

C ... Integrate over angle
      do i=1,nr
!        write(stdo,223) vavg(i,1), vavg(i,2), vavg(i,3), vavg(i,4)
        grv1(i,1) = sum(angrv1(i,:,1)*wp(:))
        grv1(i,2) = sum(angrv1(i,:,2)*wp(:))
        grv1(i,3) = sum(angrv1(i,:,3)*wp(:))
!      write(stdo,223) grv1(i,1), grv1(i,2), grv1(i,3)
      enddo

!      write(stdo,'(x,a,x,f20.8)') 'ROFI1= ', rofi(1)

      do i=1,5
        write(stdo,'(x,3(x,f20.8))') grv1(i,1), grv1(i,2),
     .  grv1(i,3)
      enddo

!  222 format(a)
!  223 format(nlml(1x,f14.8))

      end subroutine gradv1
