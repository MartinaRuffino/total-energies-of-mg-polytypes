      subroutine rotspa(mode,ldr,ldz,nr,rhor,rhoc,z)
C- Rotate a vector of spin-density-like objects to/from a quantization axis
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :Specifies how spin components of the density are stored,
Ci         :and also whether density generates rotation matrix, or
Ci         :if the rotation matrix is applied to rotate the density.
Ci         :Density may be real or complex, depending on mode.
Ci         :Cases rho is REAL    - rhor is used; rhoc is not used.
Ci         :Cases rho is COMPLEX - rhoc is used; rhor is not used.
Ci
Ci         :rhor or rhoc takes one of the following formats (see Remarks)
Ci         : (a) standard format as a linear combination of Pauli matrices
Ci         :     (density must be the complex rhoc in this case)
Ci         : (b) compressed format with diagonal parts only (vector of length 2)
Ci         : (c) compressed format x and y components in v(:,3) and v(:,4)
Ci
Ci         :1s digit for rho = rhor or rhoc
Ci         :1 input rho (complex) in standard 2x2 spinor format (a)
Ci         :2 input rho (real) in collinear format (2 elements) (b)
Ci         :3 input rho (complex) in collinear format (2 elements) (b)
Ci         :4 input rho (real) in compressed noncollinear format (4 elements)
Ci         :      Note: if 10s digit mode=1 OUTPUT rho has 4 elements
Ci         :5 input rho (complex) in compressed noncollinear format (4 elements)
Ci         :      Note: if 10s digit mode=1 OUTPUT rho has 4 elements
Ci
Ci         :10s digit for rotation z (complex)
Ci         :0 diagonalize rhor or rhoc and store eigenvectors in z
Ci         :  z diagonalizes rho, i.e. z^-1 rho z is diagonal
Ci         :  z is OUTPUT; rho is INPUT
Ci         :1 Rotate rhor or rhoc by inverse of eigenvector z, i.e.
Ci         :  rho -> z rho z^-1 and store the result in rho
Ci         :  Note: the spin part of rhor or rhoc must contain
Ci         :  four elements in this case, even if mode=2 or 3
Ci         :  z is INPUT; rho is ROTATED
Ci         :2 Same as 1, but reverse the sense of rotation z
Ci
Ci         :100s digit not used
Ci         :1000s digit
Ci         :0 eigenvalues ordered to maximize |z(1,1)|; see zdia22.f
Ci         :1 eigenvalues ordered in descending order to M>0 along spin axis
Ci   ldr  :leading dimension of rhor or rhoc
Ci   ldz :leading dimension of routc
Ci   nr    :number of elements to rotate
Cio Inputs/Outputs
Cio  rhor  :On input, a vector of spinors in real compressed format
Cio        :Output if 10s digit mode=1: rhor is rotated by z
Cio  rhoc  :On input, a vector of spinors in some complex format
Cio        :Output if 10s digit mode=1: rhor is rotated by z
Cio  z     :eigenvectors of rotation
Cio        :output if 10s digit mode=0; otherwise input.
Cr Remarks
Cr   The three components of the spin densityr rho are stored as a
Cr   linear combination of Pauli matrices sigma:
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr
Cr   rho = sigma0 rho0(r) + sigmax rhox(r) + sigmay rhoy(r) + sigmaz rhoz(r)
Cr         rho0 is the charge density
Cr         rhox, rhoy, rhoz are the components of the spin density
Cr
Cr  *Standard spinor format stores the three components of the density as
Cr          density                             magnetization
Cr       (rho++  rho+-)             M_x =  2 Re(rho+-) = Re (rho+-+rho-+)
Cr       (            )             M_y = -2 Im(rho+-) = Im (rho-+-rho+-)
Cr       (rho-+  rho+-)             M_z =  (rho++)-(rho--)
Cr                                  charge rho =  (rho++)+(rho--)
Cr
Cr  *Compressed spinor format is designed to handle collinear
Cr   and noncollinear cases in a single form.
Cr   It is a vector of length 2 (input v collinear) or 4 (noncollinear)
Cr    v(1) = rho++  v(2) = rho--  v(3) = Re rho+-  v(4) = Im rho+-
Cr
Cr   The density transforms in the same way as the Pauli spin matrices.
Cr       rhoT = U rho U+ ... see rotspu.f for a detailed example.
Cr
Cr   This routine either diagonalizes each spinor returning z, or
Cr   uses z to rotate from the local quantization axis to the original one.
Cu Updates
Cu   05 Nov 13  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ldr,ldz,nr
      double precision rhor(ldr,4)
      double complex rhoc(ldr,2,2),z(ldz,2,2)
C ... Local parameters
      logical ltmp,isanrg
      integer i,mode0,mode1,mode3
      double precision er(2)
      double complex ec(2),r12,zwk(2,2),rwk(2,2),rwk2(2,2)
C ... External calls
      external rxi,zdia22,zinv22,zmpy22

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
C     mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      ltmp = isanrg(mode3,0,1,'rotspa:','1000s digit mode',.true.)

C ... Diagonalize rho; store in z
      select case (mode0)
C     rho (complex) in standard 2x2 spinor format
      case (1)
      do  i = 1, nr
        rwk(1,1) = rhoc(i,1,1)
        rwk(1,2) = rhoc(i,1,2)
        rwk(2,1) = rhoc(i,2,1)
        rwk(2,2) = rhoc(i,2,2)
        if (mode1 == 0) then
C         mode3=1: maj spin in column 1; mode3=0: axis closest to z in col 1
          call zdia22(4-2*mode3,rwk,er,ec,z(i,1:2,1:2))
          rhoc(i,1,2) = 0; rhoc(i,2,1) = 0
          rhoc(i,1,1) = er(1)
          rhoc(i,2,2) = er(2)
        else
          zwk(:,:) = z(i,:,:)
          if (mode1 == 2) call zinv22(zwk,zwk)
          call zmpy22(zwk,rwk,rwk2)
          call zinv22(zwk,zwk)
          call zmpy22(rwk2,zwk,rwk)
          rhoc(i,1,1) = rwk(1,1)
          rhoc(i,1,2) = rwk(1,2)
          rhoc(i,2,1) = rwk(2,1)
          rhoc(i,2,2) = rwk(2,2)
        endif
      enddo

C     rho (real) in compressed format
      case (2,4)
      do  i = 1, nr
        rwk(1,1) = rhor(i,1)
        rwk(2,2) = rhor(i,2)
        if (mode0 == 2) then  ! Use if input rhor || z
          rwk(1,2) = 0
          rwk(2,1) = 0
        else
          r12 = dcmplx(rhor(i,3),rhor(i,4))
          rwk(1,2) = r12
          rwk(2,1) = dconjg(r12)
        endif
        if (mode1 == 0) then  ! Not useful but include for completeness
C         mode3=1: maj spin in column 1; mode3=0: axis closest to z in col 1
          call zdia22(4-2*mode3,rwk,er,ec,z(i,1:2,1:2))
          rhor(i,1) = er(1)
          rhor(i,2) = er(2)
          rhor(i,3) = 0; rhor(i,4) = 0
        else
          zwk(:,:) = z(i,:,:)
          call zmpy22(zwk,rwk,rwk2)
          call zinv22(zwk,zwk)
          call zmpy22(rwk2,zwk,rwk)
          rhor(i,1) = rwk(1,1)
          rhor(i,3) = dble(rwk(1,2))
          rhor(i,4) = dimag(rwk(1,2))
          rhor(i,2) = rwk(2,2)
        endif
      enddo

C     rho (complex) in compressed format
      case (3,5)
      do  i = 1, nr
        rwk(1,1) = dble(rhoc(i,1,1))
        rwk(2,2) = dble(rhoc(i,2,1))
        if (mode0 == 3) then  ! Use if input rhoc || z
          rwk(1,2) = 0
          rwk(2,1) = 0
        else
          r12 = dcmplx(dble(rhoc(i,1,2)),dble(rhoc(i,2,2)))
          rwk(1,2) = r12
          rwk(2,1) = dconjg(r12)
        endif
        if (mode1 == 0) then  ! Not useful but include for completeness
C         mode3=1: maj spin in column 1; mode3=0: axis closest to z in col 1
          call zdia22(4-2*mode3,rwk,er,ec,z(i,1:2,1:2))
          rhoc(i,1,1) = er(1)
          rhoc(i,2,1) = er(2)
          rhoc(i,1,2) = 0; rhoc(i,2,2) = 0
        else
          zwk(:,:) = z(i,:,:)
          call zmpy22(zwk,rwk,rwk2)
          call zinv22(zwk,zwk)
          call zmpy22(rwk2,zwk,rwk)
          rhoc(i,1,1) = dble(rwk(1,1))
          rhoc(i,1,2) = dble(rwk(1,2))
          rhoc(i,2,2) = dimag(rwk(1,2))
          rhoc(i,2,1) = dble(rwk(2,2))
        endif
      enddo

      case default
        call rxi('rotspv: illegal value of 1s digit mode:',mode0)

      end select

      end

C     Test rotspa
C      subroutine fmain
C      implicit none
C      double precision rinr(2,2),routr(2,2),pi
C      double complex u(2,2),rinc(2,2),rout(2,2),z(2,2)
C
C      pi = 4d0*datan(1d0)
C
CC     call rotspu(0,1,1,(/0d0,0d0,0d0/),1,u)
CC     call zprm('u',2,u,2,2,2)
C
C      rinr = 0; rinr(1,1) = 5; rinr(2,1) = 3
CC     call prmx('rinr',rinr,2,2,2)
C
CC     Rotation for magnetization along x
C      call rotspu(0,1,1,(/0d0,pi/2d0,0d0/),1,u)
CC     Rotation for magnetization along y
C      call rotspu(0,1,1,(/pi/2,pi/2,0d0/),1,u)
C
CC     Put magnetization along y, compressed real format (use later)
C      routr = 0
C      call rotspv(42,u,1,1,1,rinr,rinc,routr,rout)
CC     call prmx('rroty, (r-comp)',routr,4,4,1)
C
CC     Put magnetization along y, standard spinor format
C      call rotspv(12,u,1,1,1,rinr,rinc,routr,rout)
CC     call zprm('rroty (r-comp -> std)',2,rout,2,2,2)
C
C      print *, 'diagonalize rout(y) -> rout(z), mode 1'
C      call rotspa(1,1,1,1,routr,rout,z)
CC     call zprm('z',2,z,2,2,2)
C
C      print *, 'rotate rout(z) to (y) by using z'
C      call rotspa(11,1,1,1,routr,rout,z)
C      call zprm('rrotz -> rroty',2,rout,2,2,2)
C
C      print *, 'rout(y) -> rout(y), mode 4'
C      call rotspa(4,1,1,1,routr,rout,z)
CC     call zprm('z',2,z,2,2,2)
C
C      print *, 'rotate rout(z) to (y) by using z, mode 2'
C      call rotspa(12,1,1,1,routr,rout,z)
C      call prmx('rrotz -> rroty, (r-comp)',routr,4,4,1)
C
C      print *, 'rout(y) -> rout(y), mode 5'
C      rinc = routr
C      call rotspa(5,1,1,1,routr,rinc,z)
CC     call zprm('z',2,z,2,2,2)
C
C      print *, 'rotate rout(z) to (y) by using z, mode 3'
C      call rotspa(13,1,1,1,routr,rinc,z)
C      call zprm('rrotz -> rroty, (c-comp)',2,rinc,2,4,1)
C
C      end
