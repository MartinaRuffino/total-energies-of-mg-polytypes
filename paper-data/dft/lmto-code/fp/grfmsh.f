      subroutine grfmsh(isw,alat,ng,gv,kv,k1,k2,k3,n1,n2,n3,fn,fgrd,flap)
C- Gradient and Laplacian of a function tabulated on a uniform mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :1s digit
Ci         : 0 input function fn is real
Ci         : 1 input function fn is complex
Ci         :10s digit
Ci         : 0 Input fn is in real space
Ci         : 1 Input fn is in reciprocal space
Ci         :100s digit handles output
Ci         : 0 No output
Ci         : 1 Output fgrd, in reciprocal space, no flap
Ci         : 2 Output fgrd, in real space, no flap
Ci         : 3 Output flap, in reciprocal space, no fgrd
Ci         : 4 Output flap, in real space, no fgrd
Ci         : 5 Output fgrd and flap, in reciprocal space
Ci         : 6 Output fgrd and flap, in real space
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   vol   :cell volume
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   k1..k3:dimensions of fn,fgrd,flap
Ci   n1..n3:size of uniform mesh for which fn,fgrd,flap are tabulated
Ci   fn    :function on uniform mesh, either in real or recip space,
Ci         :depending on isw
Co Outputs
Co   fgrd  :gradient of fn either in real or recip space, depending on isw
Co   flap  :laplacian of fn either in real or recip space, depending on isw
Cl Local variables
Cl   gi    :i * G (for gradient)
Cl   g2    :(i * G)^2 (for Laplacian)
Cl   fg    :gradient of function, G space
Cl   fl    :Laplacian of function, G space
Cr Remarks
Cr Sample test:
Cr   cmfft -f9f15.10 -248 rho | mc -f9f15.10  . igx -xe | cmfft -f9f15.10 -i -248 .
Cu Updates
Cu   08 Apr 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,isw,k1,k2,k3,n1,n2,n3
      integer kv(ng,3)
      double precision alat,gv(ng,3)
      double complex fn(k1,k2,k3),fgrd(k1,k2,k3,3),flap(k1,k2,k3)
C ... Local parameters
      integer i,isw0,isw1,isw2
      double precision pi,tpiba,g2
      double complex tpibai,gi(3)
      complex(8),allocatable:: fg(:),fgg(:,:),fg2(:)

      call tcn('grfmsh')
      pi   = 4d0*datan(1d0)
      tpiba=2*pi/alat
      isw0 = mod(isw,10)
      isw1 = mod(isw/10,10)
      isw2 = mod(isw/100,10)
      if (isw0 /= 1)
     .  call rx('grfmsh: gradient of real function not implemented')

C ... FT of smooth function to reciprocal space
      if (isw1 == 0) then
C       call zprm3('fn(r)',0,fn,n1,n2,n3)
        call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,-1)
C       call zprm3('fn(G)',0,fn,n1,n2,n3)
      endif

C ... Gather function G coefficients
      allocate(fg(ng))
      call gvgetf(ng,1,kv,k1,k2,k3,fn,fg)

C ... Restore given function to real space
      if (isw1 == 0) then
        call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,1)
      endif

C ... Make iG * f(G)  and (iG)^2 * f(G)
      fg(1) = 0
      tpibai = dcmplx(0d0,1d0)*tpiba
      allocate(fgg(ng,3),fg2(ng))
      do  i = 1, ng

        gi(1) = tpibai*(gv(i,1))
        gi(2) = tpibai*(gv(i,2))
        gi(3) = tpibai*(gv(i,3))
        g2 = -tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)

        fgg(i,1) = gi(1) * fg(i)
        fgg(i,2) = gi(2) * fg(i)
        fgg(i,3) = gi(3) * fg(i)
C       fgg(i,1) = gi(1)
C       fgg(i,2) = gi(2)
C       fgg(i,3) = gi(3)
        fg2(i)   = g2    * fg(i)

      enddo

C ... Scatter gradient, laplacian into respective arrays
      if (isw2 == 1 .or. isw2 == 2 .or. isw2 == 5 .or. isw2 == 6) then
        call gvputf(ng,1,kv,k1,k2,k3,fgg(1,1),fgrd(1,1,1,1))
        call gvputf(ng,1,kv,k1,k2,k3,fgg(1,2),fgrd(1,1,1,2))
        call gvputf(ng,1,kv,k1,k2,k3,fgg(1,3),fgrd(1,1,1,3))
C       call zprm3('i Gx fn(G)',0,fgrd(1,1,1,1),n1,n2,n3)
C       call zprm3('i Gy fn(G)',0,fgrd(1,1,1,2),n1,n2,n3)
C       call zprm3('i Gz fn(G)',0,fgrd(1,1,1,3),n1,n2,n3)
      endif
      if (isw2 == 3 .or. isw2 == 4 .or. isw2 == 5 .or. isw2 == 6) then
        call gvputf(ng,1,kv,k1,k2,k3,fg2,     flap)
      endif

C ... Gradient, laplacian in real space
      if (isw2 == 2 .or. isw2 == 6) then
        call fftz3(fgrd(1,1,1,1),n1,n2,n3,k1,k2,k3,1,0,1)
        call fftz3(fgrd(1,1,1,2),n1,n2,n3,k1,k2,k3,1,0,1)
        call fftz3(fgrd(1,1,1,3),n1,n2,n3,k1,k2,k3,1,0,1)
C       call zprm3('gradx fn(r)',0,fgrd(1,1,1,1),n1,n2,n3)
C       call zprm3('grady fn(r)',0,fgrd(1,1,1,2),n1,n2,n3)
C       call zprm3('gradz fn(r)',0,fgrd(1,1,1,3),n1,n2,n3)
      endif
      if (isw2 == 4 .or. isw2 == 6) then
        call fftz3(flap         ,n1,n2,n3,k1,k2,k3,1,0,1)
C       call zprm3('lap fn(r)',0,flap,n1,n2,n3)
      endif

      deallocate(fg,fgg,fg2)

      call tcx('grfmsh')
      end

C      subroutine chk_grfmsh(s_lat,n1,n2,n3,smpot,smrho)
CC- Debugging check of grfmsh
CC  Input should be smrho, smpot = -8*pi lap^-1 rho
CC  Call this routine from smves.f, which see
C      implicit none
CC ... Passed parameters
C      integer n1,n2,n3
C      double complex smrho(n1,n2,n3),smpot(n1,n2,n3)
CC ... For structures
C      include 'structures.h'
C      type(str_lat)::   s_lat
CC ... Local parameters
C      integer i
C      double precision pi8,diff12,diffr,diffl
C      complex(8),allocatable:: zgrho1(:,:,:,:),zgrho2(:,:,:,:,:),zggrho(:,:,:),zggrhol(:,:,:),zggrho2(:,:,:),zdiff(:,:,:)
C      procedure(real(8)) :: dlength
C
C      pi8  = 8*4d0*datan(1d0)
C      allocate(zgrho1(n1,n2,n3,3),zgrho2(n1,n2,n3,3,3),zggrho(n1,n2,n3),zggrhol(n1,n2,n3),zggrho2(n1,n2,n3),zdiff(n1,n2,n3))
C
CC     call fftgrad(n1,n2,n3,alat,plat,smrho(1,1,1,i),zgrho1(1,1,i))
CC      call fftlap(n1,n2,n3,s_lat%alat,s_lat%plat,smpot,zggrho)
CC      call zprm3('lap smpot from fftlap',0,zggrho,n1,n2,n3)
C
C      call fftlap(n1,n2,n3,s_lat%alat,s_lat%plat,smpot,zggrhol)
C      call zprm3('lap smpot from fftlap',0,zggrhol,n1,n2,n3)
C
CC     Make Laplacian of smpot directly
C      call grfmsh(601,s_lat%alat,s_lat%ng,s_lat%gv,s_lat%kv,n1,n2,n3,n1,n2,n3,
C     .  smpot,zgrho1,zggrho)
C      call zprm3('lap smpot from grfmsh',0,zggrho,n1,n2,n3)
C
CC     Make Laplacian from grad.grad(smpot)
C      do  i = 1, 3
C        call grfmsh(601,s_lat%alat,s_lat%ng,s_lat%gv,s_lat%kv,n1,n2,n3,n1,n2,n3,
C     .    zgrho1(1,1,1,i),zgrho2(1,1,1,1,i),zggrho2)
C      enddo
C      call dpzero(zggrho2,2*size(zggrho2))
C      zggrho2(:,:,:) = zgrho2(:,:,:,1,1) + zgrho2(:,:,:,2,2) + zgrho2(:,:,:,3,3)
C      call zprm3('lap smpot from grad.grad(smpot)',0,zggrho2,n1,n2,n3)
C
C      zdiff = zggrho2 - zggrho
C      diff12 = dlength(2*size(zdiff),zdiff,1)
C
CC     Make -8*pi*(rho-avg rho)
C      call dcopy(2*size(zdiff),smrho,1,zdiff,1)
C      call dscal(n1*n2*n3*2,-pi8,zdiff,1)
C      call fftz3(zdiff,n1,n2,n3,n1,n2,n3,1,0,-1)
C      zdiff(1,1,1) = 0
C      call fftz3(zdiff,n1,n2,n3,n1,n2,n3,1,0,1)
C      call zprm3('-8*pi*smrho',0,zdiff,n1,n2,n3)
C
C      zggrhol = zdiff - zggrhol
C      diffl = dlength(2*size(zggrhol),zggrhol,1)
C
C      zggrho = zdiff - zggrho
C      diffr = dlength(2*size(zggrho),zggrho,1)
C
C      call info5(1,0,0,
C     .  ' difference [grad.grad(smpot) - lap(smpot)] = %,4g%N'//
C     .  ' difference [lap(smpot) - -8*pi*smrho] = %,4g%N'//
C     .  ' difference [lap(smpot,fftlap) - -8*pi*smrho] = %,4g',
C     .  diff12,diffr,diffl,0,0)
C      call rx0('done')
C
C      end

