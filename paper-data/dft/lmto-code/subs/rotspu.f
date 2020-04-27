      subroutine rotspu(mode,ib1,ib2,nbas,nl,eula,neul,u)
C- Sets up spinor rotation matrices for sites ib1..ib2
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci           1 if to reverse the sense of rotation
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1  if Euler angles are l-independent
Ci          nl if Euler angles are l-dependent
Ci          nl*nl if Euler angles are lm-dependent
Co Outputs
Co   u     :spinor rotation matrices
Cr Remarks
Cr   Rotation convention follows that of Landau and Lifshitz.
Cr   Let A=exp(i alpha/2)  B=exp(i beta/2)  C=exp(i gamma/2)
Cr   Spinor rotation matrix U is
Cr         (A*C* ReB   -A*C ImB)              ( AC ReB    A*C ImB)   (u22 -u12)
Cr     U = (                   )  and  U^-1 = (                  ) = (        )
Cr         ( AC* ImB     AC ReB)              (-AC* ImB  A*C* ReB)   (-u21 u11)
Cr
Cr   Example 1: a 180 degree rotation around x or y:
Cr          ( 0  -1)                     ( 1   0)
Cr      U = (      )    Then U*sigz*U+ = (      ) = -sigz   (flips spin)
Cr          ( 1   0)                     ( 0  -1)
Cr
Cr   Example 2: a 90 degree (active) rotation counterclockwise around z mapping
CR   (x->y and y->-x) has Euler angles (alpha=pi/2, beta=gamma=0)
Cr   A=(1+i)/sqrt(2) B=C=1 and U is
Cr                    ( 1-i  0  )
Cr      U = 1/sqrt(2) (         )                                (A)
Cr                    ( 0   1+i )
Cr
Cr   U transforms Pauli matrices sigx -> sigy and sigy -> -sigx.
Cr   Explicitly
Cr
Cr            (0  1)           (0 -i)         (1  0)
Cr     sigx = (    )    sigy = (    )  sigz = (    )             (B)
Cr            (1  0)           (i  0)         (0 -1)
Cr
Cr   It is easy to verify by direct multiplication that
Cr      U sigx U+ = sigy    and U sigy U+ = - sigx
Cr
Cr   mcx commands demonstrating this.  
Cr   mcx uses passive convention so the sign of rotation must be reversed:
Cr   mcx '-array[2,2]' 0,1,1,0 -s1,0 -a sigx '-array[2,2]' 0,1,-1,0 -s0,1 -a sigy \
Cr       -rot=-z:pi/2 -ylm~l=0~spin -a u u sigx u -i -x -x sigy -- 
Cr   mcx '-array[2,2]' 0,1,1,0 -s1,0 -a sigx '-array[2,2]' 0,1,-1,0 -s0,1 -a sigy \
Cr       -rot=z:-pi/2 -ylm~l=0~spin -a u u sigy u -i -x -x sigx -s-1 -- 
Cr
Cr   Spinors also rotate by 90 degrees.
Cr   Let chi_x, chi_y, chi_-x, chi_-y be spinors parallel to x,y,-x,-y.
Cr   Explicit representations:
Cr       chi_x           chi_y        chi_-x           chi_-y
Cr             (1)        (1-i)              (-1)        (-1-i)
Cr   1/sqrt(2) ( )    1/2 (   )    1/sqrt(2) (  )    1/2 (    )  (C)
Cr             (1)        (1+i)              ( 1)        (-1+i)
Cr   since  sigx chi_x = chi_x  and -sigx chi_-x = chi_-x ; similarly for sigy.
Cr   The following can be verified by explicit multiplication:
Cr     U chi_x  = chi_y    U chi_y =   chi_-x
Cr     U chi_-x = chi_-y   U chi_-y = -chi_x   (note 180 degree phase shift!)
Cr
Cr   Rotation of sigz by (alpha,beta,gamma) is independent of gamma:
Cr                   ( cos(beta)                sin(beta) e^(-i alpha) )
Cr     U sigz U^-1 = (                                                 )
Cr                   ( sin(beta) e^(i alpha)   -cos(beta)              )
Cr
Cr   Inverse rotation is independent of alpha:
Cr                   ( cos(beta)               -sin(beta) e^(i gamma) )
Cr     U^-1 sigz U = (                                                )
Cr                   (-sin(beta) e^(-i gamma)  -cos(beta)             )
Cr
Cr   Note: rotsp1, which performs rotations taking u from this routine,
Cr   performs the inverse  operation  U^-1 sigz U.  Thus forward rotations by
Cr   rotsp1 are passive, making them consistent with orbital rotations in roth.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ib1,ib2,nbas,neul,nl
      double complex u(2,2,nl*nl,nbas)
      double precision eula(nbas,neul,3)
C ... Local parameters
      double precision xxc,xxs,alpha,gamma
      integer ib,il,im,ilm

C --- Spinor rotation matrices for all sites ---
      do  ib = ib1, ib2
C ... Assume initially Euler angles are not l-dependent
      xxc = cos(eula(ib,1,2)/2)
      xxs = sin(eula(ib,1,2)/2)
      alpha = eula(ib,1,1)
      gamma = eula(ib,1,3)
      ilm = 0
      do  il = 1, nl
C   ... If euler angles are l dependent
        if (neul == nl) then
          xxc = dcos(eula(ib,il,2)/2)
          xxs = dsin(eula(ib,il,2)/2)
          alpha = eula(ib,il,1)
          gamma = eula(ib,il,3)
        endif
        do  im = -il+1, il-1
        ilm = ilm+1
C   ... If euler angles are lm dependent
        if (neul == nl*nl) then
          xxc = dcos(eula(ib,ilm,2)/2)
          xxs = dsin(eula(ib,ilm,2)/2)
          alpha = eula(ib,ilm,1)
          gamma = eula(ib,ilm,3)
        endif
        if (mode == 1) then
          u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
          u(1,2,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
          u(2,1,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
          u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
        else
          u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
          u(2,1,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
          u(1,2,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
          u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
        endif
C        print *, ilm,ib
C        if (ilm == 1)
C     .    call zprm('u in rotspu',2,u(1,1,ilm,ib),2,2,2)
        enddo
      enddo
      enddo
      end
C      subroutine fmain
C      double complex u(2,2),sigz(2,2),uz(2,2),uzu(2,2)
C      double precision  eula(3),pi
C      double complex A,B,C
C
C      pi = 4*datan(1d0)
C
C      eula(1) = pi/2; eula(2)=.0d0; eula(3)=.0d0
C      eula(1) = .9d0; eula(2)=.5d0; eula(3)=.2d0
C      eula(1) = .9d0; eula(2)=.5d0; eula(3)=.7d0
C
C      call rotspu(0,1,1,1,1,eula,1,u)
C      call zprm('u from rotspu',2,u,2,2,2)
C      A = cdexp(dcmplx(0d0,eula(1)/2))
C      B = cdexp(dcmplx(0d0,eula(2)/2))
C      C = cdexp(dcmplx(0d0,eula(3)/2))
C
C      u(1,1) = dble(B)*dconjg(A)*dconjg(C)
C      u(1,2) = -dimag(B)*dconjg(A)*C
C      u(2,1) = dimag(B)*A*dconjg(C)
C      u(2,2) = dble(B)*A*C
C      call zprm('u again, confirming rotpu',2,u,2,2,2)
C
CC      print *, dble(B)
CC      print *, '11', dble(B)**2 - dimag(B)**2
CC      print *, '12', -u(1,1)*u(1,2)*2
CC      print *, '12', dble(B)*dimag(B)*dconjg(A)*dconjg(C)*dconjg(A)*C*2
CC      print *, '12', dsin(eula(2))*cdexp(dcmplx(0d0,-eula(1)))
C
C      sigz(1,1) = 1; sigz(1,2) = 0; sigz(2,1) = 0; sigz(2,2) = -1
C      call zmpy22(u,sigz,uz)
C      call zinv22(u,u)
C      call zmpy22(uz,u,uzu)
C      call zprm('u sigmz u^-1',2,uzu,2,2,2)
C
C      uzu(1,1) = dcos(eula(2))
C      uzu(2,2) = -dcos(eula(2))
C      uzu(1,2) = dsin(eula(2))*cdexp(dcmplx(0d0,-eula(1)))
C      uzu(2,1) = dsin(eula(2))*cdexp(dcmplx(0d0,eula(1)))
C      call zprm(' u sigmz u^-1 analytically',2,uzu,2,2,2)
C
C      call zmpy22(u,sigz,uz)
C      call zinv22(u,u)
C      call zmpy22(uz,u,uzu)
C      call zprm('u^-1 sigmz u',2,uzu,2,2,2)
C
C      uzu(1,1) = dcos(eula(2))
C      uzu(2,2) = -dcos(eula(2))
C      uzu(1,2) = -dsin(eula(2))*cdexp(dcmplx(0d0,eula(3)))
C      uzu(2,1) = -dsin(eula(2))*cdexp(dcmplx(0d0,-eula(3)))
C      call zprm(' u^-1 sigmz u analytically',2,uzu,2,2,2)
C
C
C      u(1,1) = dble(B)*A*C
C      u(1,2) = dimag(B)*dconjg(A)*C
C      u(2,1) = -dimag(B)*A*dconjg(C)
C      u(2,2) = dble(B)*dconjg(A)*dconjg(C)
C      call zprm('u inverse',2,u,2,2,2)
C
C
C      end
