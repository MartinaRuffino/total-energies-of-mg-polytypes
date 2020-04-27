      subroutine s2sph(opt,nla,nlb,s,lsa,lsb,lspha,lsphb,sph)
C- Rotate a structure matrix from real to spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit governs storage of the imaginary part of sph.
Ci             See slatsm/ztoyy.f for a more complete description.
Ci          0: double complex with imaginary following real
Ci          1: complex*16
Ci          2: double complex with imaginary following real in columns
Ci         10s digit
Ci          0: s is real
Ci          1: s is complex.
Ci             The 1s digit of opt controls storage format for both s and sph
Ci         100s digit
Ci          0: Standard spherical harmonics (in e.g. Jackson), m ordered l,l-1,...,-l
Ci          1: Standard spherical harmonics (in e.g. Jackson), m ordered -l,-l+1,...l
Ci          2: Like 1, but alternate definition of rotation, which takes conjugate u* of standard.
Ci             This returns Y*_lm, provided the input s is REAL.
Ci             This was used as the standard convention between Oct 2013 and April 2018.
Ci             Note: since Y_l-m = (-1)^m Y^*_lm, it is similar to reversing the
Ci             m ordering of the Y_lm, but not identically so because of the (-1)^m.
Ci             Some gf routines partially handle internal inconsistencies, but not all
Ci          3: Like 0, but alternate definition of rotation
Ci         1000s digit
Ci          0  convert s (real harmonics) to s(spherical harmonics)
Ci          1  convert s (spherical harmonics) to s(real harmonics)
Ci   nla   :number of l's to rotate in the 1st (augmentation) dimension
Ci   nlb   :number of l's to rotate in the 2st (basis) dimension
Ci         :if nlb=0, u*s is returned, for lsb columns (for rotation of eigenvectors)
Ci   s     :Structure matrix, in Slater Koster form with with real harmonics
Ci         :orbitals ordered as defined e.g. in mstrx2.f
Ci   lsa   :leading dimension of s
Ci   lsb   :second dimension of s; not used unless: either
Ci         :nlb=0 or:
Ci         :s is complex and stored in format 0 (imaginary following real)
Ci   lspha :leading dimension of sph
Ci   lsphb :second dimension of sph; not used unless complex
ci         :format is 0 (imaginary following real)
Co Outputs
Co   sph   :u s u+, where u rotates Ylm from real to spherical harmonics
Co         :Special treatment when nlb is 0: returns u s instead of u s u+
Co         :See below for the construction of u.
Co         :If s is dimensioned large enough, sph can use the same address space as s
Cr Remarks
Cr   This routine rotates s to sph through a rotation u, s.t.
Cr      sph = u s u+, where u is the rotation defined as follows.
Cr
Cr   Let Y_lm = spherical harmonic.  From Jackson, 3.53:
Cr   Ylm(th,ph) = sqrt[(2l+1)/(4*pi)(l-m)!!/(l+m)!!]
Cr              * P_lm(cos(th)) exp(i m phi)
Cr   Note that Y_l,-m = (-1)^m Y^*_lm
Cr
Cr   Let R_lm = real harmonic, Y_lm = spherical harmonic.
Cr   STANDARD definitions : see Questaal web page,
Cr   https://www.questaal.org/docs/numerics/spherical_harmonics/
Cr     Y_l,-m = 1/sqrt(2)      (-i R_l,-m + R_l,m)   m>0
Cr     Y_l,m  = (-1)^m/sqrt(2) (+i R_l,-m + R_l,m)   m>0
Cr     Y_l,m  = 1/sqrt(2)      (-i R_l,m  + R_l,-m)  m<0
Cr     Y_l,0  = R_l,0                                m=0
Cr
Cr   The inverse relation is
Cr     R_l,-m =  1/sqrt(2) (i Y_l,-m - (-1)^m i Y_l,m)  m>0
Cr     R_l,m  =  1/sqrt(2) (+ Y_l,-m + (-1)^m   Y_l,m)  m>0
Cr     R_l,0  =  Y_l,0                                  m=0
Cr
Cr   Particular instances (see the Questaal web page)
Cr    l   m   R_l,m                     Y_l,m
Cr    1  -1   sqrt(3/4/pi)*y            sqrt(3/8/pi)(x - i*y)
Cr    1   0   sqrt(3/4/pi)*z            sqrt(3/4/pi)z
Cr    1   1   sqrt(3/4/pi)*x           -sqrt(3/8/pi)(x + i*y)
Cr
Cr   Matrix form (m>0)
Cr     (Y_l,-m) =           (-i         1     ) (R_l,-m)
Cr     (      ) = 1/sqrt(2) (                 ) (      ) = u^T^-1 R_l = u^* R_l
Cr     (Y_l,m ) =           (+i(-1)^m  (-1)^m ) (R_l,m )
Cr
Cr                          ( i          1     )
Cr     with   u = 1/sqrt(2) (                  ) as defined on the Questaal web page
Cr                          (-i(-1)^m   (-1)^m )
Cr
Cr   Consider a single function in the two representations
Cr       f(r) = sum_L (r_L R_L) = sum_L' (y_L' Y_L')
Cr   Relation between expansion coefficients r_L and y_L
Cr       y_L' = sum_L u_L',L r_L   r_L = sum_L' (u^-1) L,L' y_L
Cr
Cr   ALTERNATIVE definition (used between Oct 2013 and April 2018)
Cr     (Y_lm) (alternative) = (Y_lm)^* (Jackson)
Cr
Cr   Rotation of the strux S.  S are one-center expansion coefficients
Cr   such as the expansion of function H around another site:
Cr      H_L(r-dr) = (h R_L) = sum_L' S_LL' J_L' = S^R_LL'(dr) (j R_L')
Cr   If H and J are to be rotated from real harmonics R_L
Cr   to spherical harmonics Y_L, we have
Cr     (h R_L) =  sum_L' S^R_LL' (j R_L')
Cr     (h Y_L) =  sum_L' S^Y_LL' (j Y_L')
Cr   Then
Cr     (h u R_L) =  S^Y_LL' (j u R_L')
Cr   So
Cr     (h R_L) =  u^-1 S^Y_LL' u (j R_L') -> u^-1 S^Y_LL' u = S^R_LL'
Cr   Therefore
Cr     S^Y_LL' = u S^R_LL' u^-1 ->  S^Y = u S^R u+
Cr
Cr   Compare original to standard definitions for spherical harmonic repsn:
Cr         S^standard_LL' = (-1)^(l+l')(S^alt_LL')^*
Cr                        = (-1)^(l+l')(S^alt_LL')^T
Cr
Cr   Other remarks.
Cr   s2sph performs the same function as cb2sph, but is much faster.
Cb Bugs
Cb   Between Oct 2013 and April 2018, the alternative definition was thought
Cb   to be the standard one.
Cb   Not checked for the case nla ne nlb
Cu Updates
Cu   28 May 18 Rechecked definitions, and enabled l..-l or -l..l ordering
Cu   10 Oct 03 Extended to s being complex (s in kcplx mode 0)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nla,nlb,lsa,lsb,lspha,lsphb
      double precision s(lsa,lsb),sph(lspha,lsphb,2)
C ... Local parameters
C     logical :: debug=.false.
      logical lscplx,linv
      integer l,lma,lmb,nl2,lc,m,lalt,m1m
      double precision sr2
      complex(8), allocatable, target :: us(:,:)
      complex(8), pointer :: usu(:,:)
      double complex u(2,2),uc(2,2),srm1

C      if (opt < 100) then
C       call rx('s2sph test this branch')
C      endif
      nl2 = nlb*nlb
      if (nlb == 0) nl2 = lsb
      allocate(us(nla*nla,nl2),usu(nla*nla,nl2))
      if (lsphb < lsb) call rx('s2sph: illegal value for lsphb')

      srm1 = dcmplx(0d0,1d0); sr2 = dsqrt(0.5d0)
      lscplx = mod(opt/10,10) /= 0
      linv = mod(opt/1000,10) /= 0
C     if (mod(opt/10,10)/=0) call rx('s2sph complex not implemented')
      lalt = mod(opt/100,10)

      if (.not. lscplx) then
C       call yprm('s before rotation',1,s,0,lsa,nla*nla,nl2)
C --- Convert storage format and copy s to local usu (work array) ---
      else
        call ymscpk(0,mod(opt,10),1,nla*nla,nl2,s,lsa,lsb,usu,nla*nla,nl2)
C       if (opt > 1000) call zprm('s in cplx fmt',2,usu,nla*nla,nla*nla,nl2)
      endif

C     This is 2x2 rotation  defined in Remarks.  Factor (-1)^m will be incorporated later
C                           ( -i          1    )
C         U^ALT =   1/sqrt2 (                  )
C                           ( i(-1)^m   (-1)^m )
      if (lalt > 1) then
        u(1,1) = -srm1*sr2
        u(1,2) =  sr2
        u(2,1) =  srm1*sr2
        u(2,2) =  sr2
C                               (  i          1     )
C         STANDARD: u = 1/sqrt2 (                   ) = uALT^*
C                               ( -i(-1)^m   (-1)^m )
      else
        u(1,1) =  srm1*sr2
        u(1,2) =  sr2
        u(2,1) = -srm1*sr2
        u(2,2) =  sr2
      endif
      if (linv) then
        call zinv22(u,u)
      endif

C ... Make us = u*s, where u is rotation from real to sph. harm, above
      do  l = 0, nla-1
        lc = (l+1)**2-l
        do  lmb = 1, nl2
          us(lc,lmb) = s(lc,lmb)
          if (lscplx) us(lc,lmb) = usu(lc,lmb)
          m1m = -1  ! -1^m ... First m is always 1
          if (.not. lscplx) then  ! Use s directly as given
            do  m = 1, l
              us(lc-m,lmb) =     u(1,1)*s(lc-m,lmb) ! u-- s-
     .                     +     u(1,2)*s(lc+m,lmb) ! u-+ s+
              us(lc+m,lmb) = m1m*u(2,1)*s(lc-m,lmb) ! u+- s-
     .                     + m1m*u(2,2)*s(lc+m,lmb) ! u++ s+
              m1m = -m1m
            enddo
          else if (.not. linv) then  ! Use s directly as given
            do  m = 1, l
              us(lc-m,lmb) =     u(1,1)*usu(lc-m,lmb) ! u-- s-
     .                     +     u(1,2)*usu(lc+m,lmb) ! u-+ s+
              us(lc+m,lmb) = m1m*u(2,1)*usu(lc-m,lmb) ! u+- s-
     .                     + m1m*u(2,2)*usu(lc+m,lmb) ! u++ s+
              m1m = -m1m
            enddo
          else ! Inverse rotation
            do  m = 1, l
              us(lc-m,lmb) = u(1,1)    *usu(lc-m,lmb) ! u-- s-
     .                     + u(1,2)*m1m*usu(lc+m,lmb) ! u-+ s+
              us(lc+m,lmb) = u(2,1)    *usu(lc-m,lmb) ! u+- s-
     .                     + u(2,2)*m1m*usu(lc+m,lmb) ! u++ s+
              m1m = -m1m
            enddo
          endif
        enddo
      enddo

C     Debugging
C      if (lscplx .and. opt>1000) then
C        call yprmi('u s opt=%i',opt,0,3,us,0,nla**2,nla**2,nl2)
C      endif

C ... If multiply u*s only, skip to end
      if (nlb == 0) then
        if (lalt==0 .or. lalt==3) call pmorderx(1,1,1,1,nla,[0],0,nla**2,nla**2,nl2,0,nl2,0,us)
        deallocate(usu)
        usu => us
        goto 100
      endif

C ... usu = us*u+ = u*s*u+
      uc = dconjg(transpose(u))
      nl2 = nla*nla
      do  l = 0, nlb-1
        do  lma = 1, nl2
          lc = (l+1)**2-l
          usu(lma,lc) = us(lma,lc)
          m1m = -1  ! -1^m ... First m is always 1
          if (.not. linv) then
            do  m = 1, l
            usu(lma,lc-m) = us(lma,lc-m)*uc(1,1) + us(lma,lc+m)*uc(2,1)
            usu(lma,lc+m) =(us(lma,lc-m)*uc(1,2) + us(lma,lc+m)*uc(2,2))
     .                    *m1m
            m1m = -m1m
            enddo
          else ! Inverse rotation
            do  m = 1, l
            usu(lma,lc-m) = us(lma,lc-m)*uc(1,1) +
     .                      us(lma,lc+m)*m1m*uc(2,1)
            usu(lma,lc+m) = us(lma,lc-m)*uc(1,2) +
     .                      us(lma,lc+m)*m1m*uc(2,2)
            m1m = -m1m
            enddo
          endif
        enddo
      enddo

      if (lalt==0 .or. lalt==3) then
C       call pmorder(11,0,0,1,nla,[0],0,nla**2,nla**2,nla**2,nl2,nl2,usu)
        call pmorderx(11,1,1,1,nla,[0],0,nla**2,nla**2,nla**2,nl2,nl2,0,usu)
C       if (debug) call yprm('u s u+',3,usu,0,nla**2,nla**2,nl2)
      endif

  100 continue

C      if (lscplx) then
C      call yprm('u s u+',3,usu,0,nla**2,nla**2,nl2)
C      endif

C --- Convert storage format and copy to sph ---
      call ymscpk(0,1,mod(opt,10),nla*nla,nl2,usu,nla*nla,nl2,sph,lspha,lsphb)

C      if (mod(opt,10) == 0) then
C        call ztoyy(usu,nla*nla,nl2,nla*nla,nl2,1,0)
C        call ymscop(0,nla*nla,nl2,nla*nla,lspha,0,0,0,0,
C     .    usu,(nla*nlb)**2,sph,lspha*lsphb)
C      elseif (mod(opt,10) == 1) then
C        call zmscop(0,nla*nla,nl2,nla*nla,lspha,0,0,0,0,usu,sph)
C      else
C        call ztoyy(usu,nla*nla,nl2,nla*nla,nl2,1,2)
CC       call yprm('again',4,usu,0,nla*nla,nla*nla,nl2)
C        call ymscop(0,nla*nla,nl2,nla*nla*2,lspha*2,0,0,0,0,
C     .    usu,nla**2,sph,lspha)
C      endif

C     Debugging printout
C      m = 2+mod(opt,10)
C      call yprmi('output u s u+ opt=%i',opt,0,m,sph,lspha*lsphb,lspha,nla*nla,nl2)


      deallocate(us)
      if (nlb /= 0) deallocate(usu)

      end

C#ifdefC TESTV
CC     Test rotation of evec
C      subroutine fmain
C
C      integer i,m,kcplx
C      integer, parameter :: nl2=9
C      real(8) :: sin(nl2),sr2
C      complex(8) :: zsin(nl2),ans(nl2,2),srm1
CC     complex(8) :: u(nl2,nl2)
C      procedure(real(8)) :: dlength
C
C      print *, 'rotate vector to spherical harmonics'
C      sin = 0
C      forall (i=1:9) sin(i) = i
C
CC      lhaveu = .false.
CC      call cb2sph(lhaveu,u,nl2,[0d0],0,0,[0d0],nl2,[0d0])
C
C      srm1 = dcmplx(0d0,1d0); sr2 = dsqrt(0.5d0)
C      ans(1,2) = sin(1)
C      ans(2,2) = (srm1*sin(2) + sin(4))*sr2
C      ans(3,2) = sin(3)
C      ans(4,2) = (srm1*sin(2) - sin(4))*sr2
C      ans(5,2) = (srm1*sin(5) + sin(9))*sr2
C      ans(6,2) = (srm1*sin(6) + sin(8))*sr2
C      ans(7,2) = sin(7)
C      ans(8,2) = (srm1*sin(6) - sin(8))*sr2
C      ans(9,2) = (-srm1*sin(5) + sin(9))*1*sr2
C
CC      print *, 's2sph on vector, real s, with -l...l order'
CC      ans(:,1) = 0
CC      kcplx = 1; i = 100+kcplx
CC      print *, 'call s2sph with nlb=0 and opt=',i
CC      call s2sph(i,3,0,sin,9,1,9,1,ans)
CC      m = 2+kcplx
CC      call yprm('u s from s2sph',m,ans,9,9,9,2)
CC      print *, 'difference relative to inline',dlength(2*9,ans(:,1)-ans(:,2),1)
C
C      print *
C      print *, 's2sph on vector, complex s, with -l...l order'
C      ans(:,1) = 0
C      kcplx = 1; i = 110+kcplx
C      print *, 'call s2sph with nlb=0 and opt=',i
C      zsin = sin
C      call s2sph(i,3,0,zsin,9,1,9,1,ans)
C      m = 2+kcplx
C      call yprm('u s from s2sph',m,ans,9,9,9,2)
C      print *, 'difference relative to inline',dlength(2*9,ans(:,1)-ans(:,2),1)
C
C      print *
C      print *, 's2sph on vector, complex s, with l...-l order'
C      i = i - 100
C      call s2sph(i,3,0,zsin,9,1,9,1,ans)
C      m = 2+kcplx
C      call pmorder(1,0,0,1,ll(nl2)+1,[0],0,9,9,1,0,1,ans(1,2))
C      call yprm('u s from s2sph permuted order',m,ans,9,9,9,2)
C      print *, 'difference relative to inline',dlength(2*9,ans(:,1)-ans(:,2),1)
C
C      end
C#endif
C#ifdefC TEST
CC     Test rotation for REAL structure matrix sin
C      subroutine fmain
C      implicit none
C      integer i,rdm,nr,nc,ifi,fopng,ll,ns,m,kcplx,nl2,nd
C      parameter (ns=81)
C      real(8) :: sin(ns,ns),sout(ns+1,ns+2,2),u(9,9,2),swk(9,9,2),sout2(10,10,2)
C      real(8), allocatable :: srd(:,:,:)
C      logical lhaveu
C      character*(80) strn
C      procedure(real(8)) :: dlength
C
C      print *, 'rotate REAL structure matrix to YL, m=-l,...l order, read from file sin'
C
C      nr = 0; nc = 0; sout2 = 0
C      ifi = fopng('sin',-1,1)
C      i = rdm(ifi,0,0,' ',srd,nr,nc)
C      if ((ll(nr)+1)**2 /= nr) stop 'oops'
C      if ((ll(nc)+1)**2 /= nr) stop 'oops'
C      if (nr*nc > 81*81) stop 'oops'
C      allocate(srd(nr,nc,2))
C      rewind ifi
C      i = rdm(ifi,0,nr*nc,' ',srd,nr,nc)
C      sin(1:nr,1:nc) = srd(:,:,1)
C
C      call yprm('s before rotation',1,sin,ns*ns,ns,nr,nc)
C
C      kcplx = 0
C      i = 1*100+10*0+kcplx
C      i = i - 100
C      if (i/100 == 1) strn = ' (standard def, m=-l,...l order)'
C      if (i/100 == 0) strn = ' (standard def, m=l,...-l order)'
C      if (i/100 == 2) strn = ' (alternate def)'
C      print *, 'call s2sph with opt=',i,trim(strn)
C      call s2sph(i,ll(nr)+1,ll(nc)+1,sin,ns,ns,ns+1,ns+2,sout)
C
C      m = 2+kcplx
C      call yprm('u s u+ from s2sph',m,sout,(ns+1)*(ns+2),ns+1,nr,nc)
C
C      lhaveu = .false.; nl2 = 9; nd = 10
C      call cb2sph(lhaveu,u,nl2,srd,nr,nc,swk,nd,sout2)
C      call yprm('sph from cb2sph',2,sout2,nd*nd,nd,nl2,nl2)
C      if (i/100 == 0) then
C        call ztoyy(sout2,nd,nd,nd,nd,0,1)
C        call pmorder(11,0,0,1,ll(nr)+1,[0],0,nl2**2,nd,nd,nl2,nl2,sout2)
C        call ztoyy(sout2,nd,nd,nd,nd,1,0)
C        call yprm('sph from cb2sph after perm',2,sout2,nd*nd,nd,nl2,nl2)
C      endif
C
C      sout2(1:nl2,1:nl2,1:2) = sout2(1:nl2,1:nl2,1:2) - sout(1:nl2,1:nl2,1:2)
C      call yprmi('difference (mean=%;3g)',dlength(size(sout2),sout2,1),0,2,sout2,nd*nd,nd,nl2,nl2)
C
C      end
C#endif
C#ifdefC TEST2
C     Test rotation for COMPLEX structure matrix sin
C      subroutine fmain
C      implicit none
C      integer i,rdm,nr,nc,ifi,fopng,ll,ns,m,isw,kcplx
C      parameter (ns=81)
C      real(8) :: sin(ns,ns,2),sout(ns+1,ns+2,2)
C      real(8), allocatable :: srd(:,:,:)
C      logical lscplx
C
C      nr = 0; nc = 0
C      ifi = fopng('sin',-1,1)
C      i = rdm(ifi,0,0,' ',srd,nr,nc)
C      if ((ll(nr)+1)**2 /= nr) stop 'oops'
C      if ((ll(nc)+1)**2 /= nr) stop 'oops'
C      if (nr*nc > 81*81) stop 'oops'
C      allocate(srd(nr,nc,2)); srd = 0
C      rewind ifi
C      i = rdm(ifi,30,nr*nc*2,' ',srd,nr,nc)
C      sin(1:nr,1:nc,1:2) = srd(:,:,1:2)
C
C      lscplx = .true.
C      m = 2 ; if (.not. lscplx) m = 1
CC     call yprm('s before rotation',m,sin,ns*ns,ns,nr,nc)
C
C      kcplx = 0
C      call ztoyy(sin,ns,ns,ns,ns,0,kcplx)
CC      call yprm('s in kcplx format',2+kcplx,sin,ns*ns,ns,nr,nc)
C
C      i = 100+10*isw(lscplx)+kcplx
C      print *, 'call s2sph with opt=',i
C      call s2sph(i,ll(nr)+1,ll(nc)+1,sin,ns,ns,ns+1,ns+2,sout)
C
CC     Debugging printout
C      m = 2+kcplx
C      call yprm('u s u+',m,sout,(ns+1)*(ns+2),ns+1,nr,nc)
C
C      end
C#endif
