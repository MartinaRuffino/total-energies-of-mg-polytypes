      subroutine rotspv(mode,u,ldin,ldout,nr,rinr,rinc,routr,routc)
C- Rotate a vector of spin-density-like objects by a fixed rotation
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :How rin and routc store spin components of the density.
Ci         :Cases rin is REAL     - data taken from rinr (rinc is not used)
Ci         :Cases rin is COMMPLEX - data taken from rinc (rinr is not used)
Ci         :Styles are one of the following (see Remarks)
Ci         : *standard format as linear combination of Pauli matrices
Ci         : *compressed format with diagonal parts only (vector of length 2)
Ci         : *compressed format x and y components in v(:,3) and v(:,4)
Ci         :1s digit for input rin = rinr or rinc
Ci         :1 rin (complex) in standard 2x2 spinor format
Ci         :2 rin (real) in collinear format (2 elements)
Ci         :3 rin (complex) in collinear format (2 elements)
Ci         :4 rin (real) in compressed noncollinear format (4 elements)
Ci         :5 rin (complex) in compressed noncollinear format (4 elements)
Ci         :10s digit for output routc (complex)
Ci         :1 routc (complex) stored in standard spinor format
Ci         :4 routr (real) stored compressed noncollinear format (4 elements)
Ci         :5 routc (complex) stored in compressed noncollinear format (4 elements)
Ci         :100 digit  Defines u.
Cii        :0 u is a spinor rotation matrix
Ci         :1 Reverse the sense of rotation u
Ci         :2 Use a 2x2 unit matrix for u (useful for transposing formats).
Ci         :  u is not used in this case
Ci         :4 u is a quantization axis (a REAL vector of length 3)
Ci         :  Rotate FROM local quantization axis TO global axis
Ci         :5 Same as 4, but rotate TO local axis FROM global axis
Ci         :  (reverse sense of rotation)
Ci         :6 u is a triplet of Euler angles (a REAL vector of length 3)
Ci         :7 Same as 6, but reverse sense of rotation
Ci         :1000s digit used in cases 4,5
Ci         :0 eigenvalues ordered to maximize |z(1,1)|; see zdia22.f
Ci         :1 eigenvalues ordered in descending order to M>0 along spin axis
Ci         :10000s digit
Ci         :1 suppress rotation on the left so rhoT = rho U+ ... see Remarks
Ci         :2 suppress rotation on the right so rhoT = U rho ... see Remarks
Ci   u     :usually, spin rotation matrix. See rotspu.f.
Ci         :But see 100s digit mode.
Ci   ldin  :leading dimension of rinr or rinc
Ci   ldout :leading dimension of routc
Ci   nr    :number of elements to rotate
Ci   rinr  :a vector of spinors to be rotated
Ci         :real in compressed format with spin parallel to z axis
Ci   rinc  :a vector of spinors to be rotated
Co Outputs
Co   rinc   :10s digit mode 4: rin rotated by u
Co   routc  :10s digit mode 1,5: rin rotated by u
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
Cr   Note:
Cr     It is assumed that the density matrix is hermitian,
Cr     or a hermitian matrix scaled by phase.
Cr     For now, the phase must be zero if 1s or 10s digit job is 4 or 5
Cu Updates
Cu   29 Oct 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ldin,ldout,nr
      double precision rinr(ldin,4),routr(ldout,2,2)
      double complex u(2,2),rinc(ldin,2,2)
      complex(8),target :: routc(ldout,2,2)
C ... Dynamically allocated arrays
      complex(8),allocatable :: rwk(:,:,:)
C ... Local parameters
      logical ltmp,isanrg
      integer i,mode0,mode1,mode2,mode3,mode4
      double precision r11,r22,saxis(3),eula(3)
      double complex cx(2,2),r12
      complex(8),target :: ul0(2,2),unity(2,2)
      complex(8), pointer :: ul(:,:)

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      mode4 = mod(mode/10000,10)
      ltmp = isanrg(mode3,0,1,'rotspv:','1000s digit mode',.true.)
      ltmp = isanrg(mode4,0,2,'rotspv:','10000s digit mode',.true.)
      allocate(rwk(ldout,2,2))
      unity(1,2) = 0 ; unity(2,1) = 0 ; unity(1,1) = 1; unity(2,2) = 1

C ... ul = u or u+
      select case (mode2)
      case (0); ul0 = u
      case (1); call zinv22(u,ul0)
      case (2); ul0 = unity
      case (4,5);
        call dcopy(3,u,1,saxis,1)
        cx(1,1) = saxis(3)/2
        cx(2,2) = -saxis(3)/2
        cx(1,2) = dcmplx(saxis(1),-saxis(2))/2
        cx(2,1) = dconjg(cx(1,2))
C       Rotation matrix that diagonalizes saxis
C       mode3=1: maj spin in column 1; mode3=0: axis closest to z in col 1
        call zdia22(4-2*mode3,cx,saxis,saxis,ul0) ! evals -> saxis
        if (mode2 == 5) call zinv22(ul0,ul0)
C     u passed to this routine is a set of 3 Euler angles
      case (6,7);
        call dcopy(3,u,1,eula,1)
        call rotspu(0,1,1,eula,1,ul0)
        if (mode2 == 7) call zinv22(ul0,ul0)
      case default
        call rxi('rotspv: illegal value of 10s digit mode:',mode2)
      end select
C     call zprm('ul',2,ul0,2,2,2)

      if (mode4 /= 1) then
        ul => ul0
      else
        ul => unity
      endif

C ... Make u v
      select case (mode0)
C     rin (complex) in standard 2x2 spinor format
      case (1)
        do  i = 1, nr
          rwk(i,1,1) = ul(1,1)*rinc(i,1,1) + ul(1,2)*rinc(i,2,1)
          rwk(i,1,2) = ul(1,1)*rinc(i,1,2) + ul(1,2)*rinc(i,2,2)
          rwk(i,2,1) = ul(2,1)*rinc(i,1,1) + ul(2,2)*rinc(i,2,1)
          rwk(i,2,2) = ul(2,1)*rinc(i,1,2) + ul(2,2)*rinc(i,2,2)
        enddo

C     rin (real) in compressed format (collinear, 2 elements)
      case (2)
        do  i = 1, nr
          rwk(i,1,1) = ul(1,1)*rinr(i,1)
          rwk(i,1,2) = ul(1,2)*rinr(i,2)
          rwk(i,2,1) = ul(2,1)*rinr(i,1)
          rwk(i,2,2) = ul(2,2)*rinr(i,2)
        enddo

C     rin (complex) in compressed format (collinear, 2 elements)
      case (3)
        do  i = 1, nr
          rwk(i,1,1) = ul(1,1)*rinc(i,1,1)  ! u11 r11
          rwk(i,1,2) = ul(1,2)*rinc(i,2,1)  ! u12 r22
          rwk(i,2,1) = ul(2,1)*rinc(i,1,1)  ! u21 r11
          rwk(i,2,2) = ul(2,2)*rinc(i,2,1)  ! u22 r22
        enddo

C     rin (real) in compressed format (noncollinear, 4 elements)
      case (4)
        do  i = 1, nr
          r11 = rinr(i,1)
          r22 = rinr(i,2)
          r12 = dcmplx(rinr(i,3),rinr(i,4))
          rwk(i,1,1) = ul(1,1)*r11 + ul(1,2)*dconjg(r12)
          rwk(i,1,2) = ul(1,1)*r12 + ul(1,2)*r22
          rwk(i,2,1) = ul(2,1)*r11 + ul(2,2)*dconjg(r12)
          rwk(i,2,2) = ul(2,1)*r12 + ul(2,2)*r22
        enddo

C     rin (complex) in compressed format (noncollinear, 4 elements)
      case (5)
        do  i = 1, nr
          r11 = dble(rinc(i,1,1))
          r22 = dble(rinc(i,2,1))
          r12 = dcmplx(dble(rinc(i,1,2)),dble(rinc(i,2,2)))
          rwk(i,1,1) = ul(1,1)*r11 + ul(1,2)*dconjg(r12)
          rwk(i,1,2) = ul(1,1)*r12 + ul(1,2)*r22
          rwk(i,2,1) = ul(2,1)*r11 + ul(2,2)*dconjg(r12)
          rwk(i,2,2) = ul(2,1)*r12 + ul(2,2)*r22
        enddo

      case default
        call rxi('rotspv: illegal value of 1s digit mode:',mode0)

      end select

C     return

      if (mode4 /= 2) then
        ul => ul0
      else
        ul => unity
      endif


C ... Make u v u+
      call zinv22(ul,ul)

      select case (mode1)

      case (1)

C       routc is stored in standard spinor format
        do  i = 1, nr
          routc(i,1,1) = rwk(i,1,1)*ul(1,1) + rwk(i,1,2)*ul(2,1)
          routc(i,1,2) = rwk(i,1,1)*ul(1,2) + rwk(i,1,2)*ul(2,2)
          routc(i,2,1) = rwk(i,2,1)*ul(1,1) + rwk(i,2,2)*ul(2,1)
          routc(i,2,2) = rwk(i,2,1)*ul(1,2) + rwk(i,2,2)*ul(2,2)
        enddo

C     routr is stored in compressed spinor format (4 elements)
      case (4)
        do  i = 1, nr
          cx(1,1) = rwk(i,1,1)*ul(1,1) + rwk(i,1,2)*ul(2,1)
          cx(1,2) = rwk(i,1,1)*ul(1,2) + rwk(i,1,2)*ul(2,2)
          cx(2,1) = rwk(i,2,1)*ul(1,1) + rwk(i,2,2)*ul(2,1)
          cx(2,2) = rwk(i,2,1)*ul(1,2) + rwk(i,2,2)*ul(2,2)

          routr(i,1,1) = cx(1,1) ! routr(i,1)
          routr(i,1,2) = dble(cx(1,2)+cx(2,1))/2 ! routr(i,3)
          routr(i,2,2) = dimag(cx(1,2)-cx(2,1))/2 ! routr(i,4)
          routr(i,2,1) = cx(2,2) ! routr(i,2)
        enddo

C     routc is stored in compressed spinor format (4 elements)
      case (5)
        do  i = 1, nr
          cx(1,1) = rwk(i,1,1)*ul(1,1) + rwk(i,1,2)*ul(2,1)
          cx(1,2) = rwk(i,1,1)*ul(1,2) + rwk(i,1,2)*ul(2,2)
          cx(2,1) = rwk(i,2,1)*ul(1,1) + rwk(i,2,2)*ul(2,1)
          cx(2,2) = rwk(i,2,1)*ul(1,2) + rwk(i,2,2)*ul(2,2)

          routc(i,1,1) = cx(1,1) ! routc(i,1)
          routc(i,1,2) = dble(cx(1,2)+cx(2,1))/2 ! routc(i,3)
          routc(i,2,2) = dimag(cx(1,2)-cx(2,1))/2 ! routc(i,4)
          routc(i,2,1) = cx(2,2) ! routc(i,2)
        enddo

      case default
        call rxi('rotspv: illegal value of 10s digit mode:',mode1)

      end select

      deallocate(rwk)

      end

C     Test rotspv
C      subroutine fmain
C      implicit none
C      double precision rinr(2,2),routr(2,2),pi
C      double complex u(2,2),rinc(2,2),rout(2,2)
C
C      pi = 4d0*datan(1d0)
C
CC     call rotspu(0,1,1,(/0d0,0d0,0d0/),1,u)
CC     call zprm('u',2,u,2,2,2)
C
C      rinr = 0; rinr(1,1) = 5; rinr(2,1) = 3
CC     call prmx('rinr',rinr,2,2,2)
C
C      call rotspu(0,1,1,(/0d0,pi/2d0,0d0/),1,u)
CC     Should put magnetization along y, standard spinor format
C      call rotspu(0,1,1,(/pi/2,pi/2,0d0/),1,u)
C      call rotspv(12,u,1,1,1,rinr,rinc,routr,rout)
C      call zprm('rroty (r-comp -> std)',2,rout,2,2,2)
C
C
CC     Should restore magnetization to
C      rinr = -1; rinc = -1
C      call rotspv(141,u,1,1,1,routr,rout,rinr,rinc)
C      call prmx('r-orig, (std -> r-comp)',rinr,4,4,1)
C
CC     Should put magnetization along y, compressed spinor format
C      call rotspv(52,u,1,1,1,rinr,rinc,routr,rout)
C      call zprm('rroty  (r-comp -> comp)',2,rout,4,4,1)
C
C
C
C
CC     Should restore magnetization to z, standard spinor format
C      call rotspv(115,u,1,1,1,rinr,rout,routr,rinc)
C      call zprm('r-orig (c-comp -> std)',2,rinc,2,2,2)
CC     Should restore magnetization to z, compressed spinor format
C      call rotspv(155,u,1,1,1,rinr,rout,routr,rinc)
C      call zprm('r-orig, (c-comp -> comp)',2,rinc,4,4,1)
C
CC     Should put magnetization along xy, compressed format
C      call rotspu(0,1,1,(/pi/4,pi/2,0d0/),1,u)
C      call rotspv(13,u,1,1,1,rinr,rinc,routr,rout)
C      call zprm('rrotxy (c-comp -> std)',2,rout,2,2,2)
C      call rotspv(53,u,1,1,1,rinr,rinc,routr,rout)
C      call zprm('rrotxy, (c-comp -> c-comp)',2,rout,4,4,1)
C      call rotspv(43,u,1,1,1,rinr,rinc,routr,rout)
C      call prmx('rrotxy, (c-comp -> r-comp)',routr,4,4,1)
C
CC     Should restore magnetization to z, standard spinor format
C      call rotspv(114,u,1,1,1,routr,rout,rinr,rinc)
C      call zprm('r-orig (c-comp -> std)',2,rinc,2,2,2)
CC     Should restore magnetization to z, compress real spinor
C      rinr = 0
C      call rotspv(144,u,1,1,1,routr,rout,rinr,rinc)
C      call prmx('r-orig, (c-comp -> r-comp)',rinr,4,4,1)
C
C      end
