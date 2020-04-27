      subroutine cb2sph(lhaveu,u,nu,cub,nlma,nlmb,wk,nd,sph)
C- Rotate one block of strux from cubic to spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci   cub,nlma,nlmb: strux and dimensions in cubic harmonics
Ci   lhaveu: T, rotation matrix u already made (F, cb2sph makes u)
Ci   wk:     work array of length nlma*nlmb*2
Ci   nd:     dimensioning parameter for sph
Cio Inputs/Outputs
Ci   u,nu:   Rotation matrix and dimension to do u S_real u+ = S_sph
Ci           If lhaveu=T, u is input
Ci           If lhaveu=F, u is generated and output
Co Outputs
Co   sph:    strx in spherical harmonics
Co   lhaveu: set to T
Cr Remarks
Cr   This routine makes and uses a rotation matrix u corresponding to
Cr   r^-1 in rcub2sph.f, which see.  See also Standard definition in
Cr   https://www.questaal.org/docs/numerics/spherical_harmonics
Cr
Cr   The test case shows that u*Rlm = Ylm
Cr   Other remarks:
Cr   lhaveu may be T on subsequent calls if u is preserved.
Cr   If nlma or nlmb=0, only u is made; cub,wk,nd,sph are unused.
Cu Updates
Cu   28 Apr 98 Changed u->u*; definition compatible with e.g. Jackson
C ----------------------------------------------------------------------
      implicit none
C Passed variables
      logical lhaveu
      integer nu,nlma,nlmb,nd
      double precision u(nu,nu,2),cub(nlma,nlmb),wk(nlma,nlmb,2),
     .  sph(nd,nd,2)
C Local variables
      integer na,nb,nl,l,m,la,la0,nla,lb,lb0,nlb,ll,lc
      double precision a,sm

C --- Make rotation matrix u ---
      nl = ll(nu)+1
      if (.not. lhaveu) then
        a = 1/dsqrt(2d0)
        call dpzero(u,nu*nu*2)
        do  l = 0, nl-1
          lc = (l+1)**2-l
          u(lc,lc,1) = 1
          sm = -1
C         This is r^-1, see Eq.(2) in rcub2sph.f
          do  m = 1, l
            u(lc-m,lc-m,2) = a     ! (-m,-m) i/sr2
            u(lc-m,lc+m,1) = a     ! (-m,+m) 1/sr2
            u(lc+m,lc-m,2) = -sm*a ! (+m,-m) -i(-1)^m/sr2
            u(lc+m,lc+m,1) = sm*a  ! (+m,+m) (-1)^m/sr2
            sm = -sm
          enddo
        enddo
        lhaveu = .true.
      endif

C     call yprm('u in cb2sph',2,u,nu*nu,nu,nu,nu)
      if (nlma <= 0 .or. nlmb <= 0) return
      na = ll(nlma)+1
      nb = ll(nlmb)+1

C --- make cub u+ by l-blocks ---
      do  lb = 0, nb-1
        lb0 = lb**2+1
        nlb = 2*lb+1
        call dgemm('N','T',nlma,nlb,nlb,1d0,cub(1,lb0),nlma,
     .    u(lb0,lb0,1),nu,0d0,wk(1,lb0,1),nlma)
        call dgemm('N','T',nlma,nlb,nlb,-1d0,cub(1,lb0),nlma,
     .    u(lb0,lb0,2),nu,0d0,wk(1,lb0,2),nlma)
      enddo

C --- make u (cub u+) by la-blocks ---
      do  la = 0, na-1
        la0 = la**2+1
        nla = 2*la+1
        call yygemm('N','N',nla,nlmb,nla,1d0,u(la0,la0,1),u(la0,la0,2),
     .              nu,wk(la0,1,1),wk(la0,1,2),nlma,0d0,sph(la0,1,1),
     .              sph(la0,1,2),nd)
      enddo

C     call yprm('sph in cb2sph',2,sph,nd*nd,nd,nlma,nlmb)

      end

C     This test calls cb2sph to check rotation from real to spherical harmonics
C      subroutine fmain
C      implicit none
C      integer i,l,m,nl,nl2,lc,am
C      double precision u(16,16,2),theta,phi,ylm(16,2),ylms(16,2),sr2,pi
C      double precision scub(16,16), swk(16,2,16), ssph(16,2,16),x,y,z,rr
C      logical lf
C
CC ... Use cp2sph to generate ssph
C      nl = 4
C      nl2 = nl*nl
C      call dpzero(scub,nl2*nl2)
C      lf = .false.
C      call cb2sph(lf,u,nl2,scub,0,0,swk,1,ssph)
CC     call yprm('u',2,u,16*16,16,16,16)
CC     call yprm('u',2,u,16*16,16,4,4)
C
CC ... Test cp2sph
C      sr2 = 1/dsqrt(2d0)
C      pi = 4*datan(1d0)
C      theta = 1.57079633d0; phi = 1.57079633d0
C      call setpr(10)
C      call info2(1,0,0,'theta,phi (radians) (default=%;8d,%;8d)? ',theta,phi)
C      read(*,*) theta,phi
C      x = sin(theta)*cos(phi)
C      y = sin(theta)*sin(phi)
C      z = cos(theta)
C
C      print 332, x,y,z,cos(phi),sin(phi),cos(theta),sin(theta),
C     .  -exp((0d0,1d0)*phi)*sin(theta)*sqrt(3d0/8d0/pi),
C     .  cos(theta)*sqrt(3d0/4d0/pi),
C     .  exp((0d0,-1d0)*phi)*sin(theta)*sqrt(3d0/8d0/pi)
C      print 333,
C     .  exp((0d0,2d0)*phi)*sin(theta)**2*sqrt(15d0/32d0/pi),
C     . -exp((0d0,1d0)*phi)*sin(theta)*cos(theta)*2*sqrt(15d0/32d0/pi),
C     .  (3*(cos(theta))**2-1)*sqrt(5d0/16d0/pi)
C  332 format('  x,y,z=',3f10.6/
C     .       '  cos(phi)=',f10.6,'  sin(phi)=',f10.6,
C     .       '  cos(th)=',f10.6,'  sin(th)=',f10.6/
C     .       '  Y11=',2f12.6,'     Y10=',f12.6,12x,'    Y1-1=',2f12.6)
C  333 format('  Y22=',2f12.6,'     Y21=',2f12.6,'     Y20=',f12.6)
C
C      call dpzero(ylm,2)
C      call ropyln(1,x,y,z,3,1,ylm,rr)
C
C      call ygemm('N','N',
C     .  nl2,1,nl2,1d0,u,nl2*nl2,nl2,ylm,nl2,nl2,0d0,ylms,nl2,nl2)
C
C
CC      print *, ' OLD definitions (compile cb2sph with OLD)'
CC      print *, ' l  m  L       Y_real        Re Y_L      Im Y_L'
CCC     echo .5 .2 | a.out  | tee out
CCC     extract-lines --quiet _real END 1 out | mc .  -e2 x5-x7 x6-x8
CC      i = 0
CC      do  l = 0, 3
CC      do  m = -l, l
CC        i = i+1
CC        if (m < 0) then
CC          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
CC     .      sr2*(ylm(i-2*m,1)),sr2*(ylm(i,1))
CC  334     format(3i3,2x,f12.6,2x,2f12.6,2x,2f12.6)
CC        elseif (m > 0) then
CC          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
CC     .      sr2*(-1)**m*ylm(i,1), sr2*(-1)**m*(-1)*ylm(i-2*m,1)
CC        else
CC          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
CC     .      ylms(i,1), ylms(i,2)
CC        endif
CC
CC      enddo
CC      enddo
C
C      print *, ' STANDARD definitions'
C      print *, ' l  m  L        R_L          Re Y_L      Im Y_L         +/- R_+ + i R_-'
C
CC     for this table check that these cols are zero:
CC     echo .5 .2 | a.out  | tee out
CC     extract-lines --quiet R_L END 1 out | mc .  -e2 x5-x7 x6-x8
C      i = 0
C      do  l = 0, 3
C      do  m = -l, l
C        i = i+1
C        lc = (l+1)**2-l
C        am = iabs(m)
C        if (m < 0) then
C          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
C     .      sr2*(dcmplx(0d0,1d0)*ylm(lc-am,1)+ylm(lc+am,1))
C  334     format(3i3,2x,f12.6,2x,2f12.6,2x,2f12.6)
C        elseif (m > 0) then
C          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
C     .      sr2*(-1)**m*(dcmplx(0d0,-1d0)*ylm(lc-am,1)+ylm(lc+am,1))
C        else
C          print 334, l,m,i, ylm(i,1), ylms(i,1), ylms(i,2),
C     .      ylms(i,1), ylms(i,2)
C        endif
C
C      enddo
C      enddo
C
C
C      end
