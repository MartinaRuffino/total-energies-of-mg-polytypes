      subroutine scglp1(nlm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
C- Makes Clebsch-Gordan coefficients c_m' to grad (H_lm) = sum_m' [c_m' H_l+1m'], H_lm = Hankel fn of 0 energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlm   :mlm = |nlm|:  compound (l,m) index to Ylm, real harmonic functions
Ci         :If nlm<0:  set (kx1,kx2) and also (ky1,ky2) are ordered by increasing value
Co Outputs
Co   kz    :ilm of z component for ll(mlm)+1; see Remarks
Co   cz    :coefficient for z component for ll(mlm)+1
Co   kx1   :ilm of 1st x component for ll(mlm)+1; see Remarks
Co   kx2   :ilm of 2nd x component for ll(mlm)+1; see Remarks
Co   cx1   :coefficient for 1st x component for ll(mlm)+1
Co   cx2   :coefficient for 2nd x component for ll(mlm)+1
Co   ky1   :ilm of 1st y component for ll(mlm)+1; see Remarks
Co   ky2   :ilm of 2nd y component for ll(mlm)+1; see Remarks
Co   cy1   :coefficient for 1st y component for ll(mlm)+1
Co   cy2   :coefficient for 2nd y component for ll(mlm)+1
Cr Remarks
Cr   Coefficients for Y1(grad) H_L, omitting factor sqrt(3/4pi) from Y_1m = sqrt(3/4pi)(y,z,x)
Cr   In the special case of spherical hankel function H_L(E=0,r) = r^-l-1/(2l-1)!! Y_L,
Cr   the gradient operator maps  grad H_L  into a linear combination of H_K.
Cr   Coefficients and corresponding klm quantum numbers are returned for nonzero components.
Cr
Cr   This routine is a special case of the more general grad f(r) Y_L.  See ylmbyr below.
Cr   grad f(r) Y_L consists of a linear combination of four terms:
Cr   grad f(r) Y_(l,m) = (df/dr - l*f/r)     [c11 Y(l+1,m1) + c21 Y(l+1,m2)]
Cr                     + (df/dr + (l+1)*f/r) [c12 Y(l-1,m1) + c22 Y(l-1,m2)]
Cr   m1,m2 depend on which component (x,y,z). See ylmbyr below for fuller discussion
Cr   The top row vanishes when f(r) = r^l
Cr   The bottom row vanishes when f(r) = r^-l-1 (i.e. the case applicable to scglp1)
Cr
Cr   c12 and c22 can be obtained from c11 and c21, as done in ylmbyr.
Cu Updates
Cu   20 May 19 New mkcga,mkcgar
Cu   1  May 00 Adapted from nfp scglp1.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlm,kz,kx1,kx2,ky1,ky2
      double precision cz,cx1,cx2,cy1,cy2
C ... Local parameters
      integer mlm,l,lav,m,mm,isg,kav,ma,mb
      double precision bot,top,tap,cofa,cofb
      procedure(integer) :: ll

      mlm = iabs(nlm)
      l = ll(mlm)
      lav = l*l+l+1
      m = mlm-lav
      mm = iabs(m)
      isg = 1
      if (m < 0) isg = -1
      kav = (l+1)*(l+1)+(l+1)+1
      bot = (2*l+1)*(2*l+3)
C     z coefficient
      kz = kav+m
      cz = dsqrt((l+mm+1)*(l-mm+1)/bot)
C     x,y coefficients
      top = (l+1+mm)*(l+2+mm)/2d0
      tap = (l+1-mm)*(l+2-mm)/2d0
      if (mm /= 0) top = 0.5d0*top
      if (mm /= 1) tap = 0.5d0*tap
      cofa = dsqrt(top/bot)
      cofb = dsqrt(tap/bot)
      ma = isg*(mm+1)  ! m+1 if m>0; m-1 if m<0
      mb = isg*(mm-1)  ! m-1 if m>0; m+1 if m<0
      kx1 = kav+ma
      cx1 = cofa                ! CG(x) for : l+1,m+1 if m>0;  l+1,m-1 if m<0
      kx2 = kx1
      cx2 = 0d0
C     call info5(1,0,0,' scglp1 l,m%,2i%,3i  cx1 %1;12,6D',l,m,cx1*(l+2),4,5)
      if (m/=-1 .and. m/=0) then
        kx2 = kav+mb            ! kx2 < kx1
        cx2 = -cofb             ! CG(x) for : l+1,m-1 if m>0;  l+1,m+1 if m<0
      endif
      ky1 = kav-ma
      cy1 = isg*cofa
      ky2 = ky1
      cy2 = 0d0
      if (m/=0 .and. m/=1) then
        ky2 = kav-mb
        cy2 = isg*cofb
      endif

      if (nlm>0) return

C     Reorder cx1,cx2 so that cx1 corresponds to m-1, cx2 to m+1
      if (kx1>kx2 .or. kx1==kx2 .and. kx1-kav-m>0 .and. cx2==0) then
        mlm=kx1; kx1=kx2; kx2=mlm
        cofa=cx1; cx1=cx2; cx2=cofa
      endif
C     Reorder cy1,cy2 so that cy1 corresponds to m-1, cy2 to m+1
      if (ky1>ky2 .or. ky1==ky2 .and. ky1-kav+m>0 .and. cy2==0) then
        mlm=ky1; ky1=ky2; ky2=mlm
        cofa=cy1; cy1=cy2; cy2=cofa
      endif

      end

      subroutine ylmbyr(lmax,cgcof,mg)
C- Return coefficients C_K to decomposition of grad (Y_L) or \vec{r}*Y_L  into  sum_K (C_K Y_K)
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmax  : coefficients are computed for l=0..lmax
Co Outputs
Co   cgcof : cgcof(1:2,1:2,(x,y,z),(l,m)) = the family of coefficients derived from Y_(l,m)
Co         : For each ilm, there can be up to four coefficients coupling
Co         : l to k=l+/-1  and |mk| = |ml| +/-1 (for grad_x,y) or mk=ml (for grad_z)
Co         : cgcof(:,:,:,ilm) coffs C_L mapping grad Y_L to Y_K, apart from radial factors
Cr         : cgcof(:,:,1,:) = coffs entering into grad_x f(r) Y_L
Cr         : cgcof(:,:,2,:) = coffs entering into grad_y f(r) Y_L
Cr         : cgcof(:,:,3,:) = coffs entering into grad_z f(r) Y_L
Cr         : cgcof(:,1,:,:) = coffs to k=l+1 Y_k,mk components of grad f(r) Y_L; see Remarks
Cr         : cgcof(:,2,:,:) = coffs to k=l-1 Y_k,mk components of grad f(r) Y_L; see Remarks
Cr         : cgcof(1,:,:,:) = coffs mk=ml-1 (grad_x), mk=-ml-1 (grad_y) and mk=ml (grad_z) component
Cr         : cgcof(2,:,:,:) = coffs mk=ml+1 (grad_x), mk=-ml+1 (grad_y) component
Cr         :
Cr         : Write out each term explicitly. Let Y'(l,m,p) = grad_p f Y(l,m)
Cr         : Y'(l,m,x) = f1 [cgcof(1,1,1,(l,m)) Y(l+1,m-1)  + cgcof(2,1,1,(l,m)) Y(l+1,m+1)]
Cr         :           + f2 [cgcof(1,2,1,(l,m)) Y(l-1,m-1)  + cgcof(2,2,1,(l,m)) Y(l-1,m+1)]
Cr         : Y'(l,m,y) = f1 [cgcof(1,1,2,(l,m)) Y(l+1,-m-1) + cgcof(2,1,2,(l,m)) Y(l+1,-m+1)]
Cr         :           + f2 [cgcof(1,2,2,(l,m)) Y(l-1,-m-1) + cgcof(2,2,2,(l,m)) Y(l-1,-m+1)]
Cr         : Y'(l,m,z) = f1 [cgcof(1,1,3,(l,m)) Y(l+1,m)    + cgcof(2,1,3,(l,m)) Y(l+1,m)]
Cr         : where f1 = (df/dr - l*f/r) and f2 = (df/dr + (l+1)*f/r)
Cr         :
Cr    mg   : mg(1:2,p,L) = mk quantum number cgcof(1:2,:,p,L).  That is
Cr         : grad_p Y_L [L = (l,m)] is mapped into a linear combination of Y_K, [K = (lk,mk)]
Cr         : lk : lk = l+1 for cgcof(:,1,p,L); lk = l-1 for cgcof(:,2,p,L)
Cr         : mk : For either lk, mk corresponding to cgcof(1:2,:,p,L) is mg(1:2,p,L) where
Cr         :        x                  y                     z
Cr         : mg(1,1,L) = m-1    mg(1,2,L) = -m-1      mg(1,3,L) = m
Cr         : mg(2,1,L) = m+1    mg(2,2,L) = -m+1      mg(2,3,L) = m
Cr         :
Cr         : Then
Cr         : grad_p f(r) Y_L = f1 sum_i=1,2 [ cgcof(i,1,p,l;m) Y_l+1;mg(i,p,m) ]
Cr         :                 + f2 sum_i=1,2 [ cgcof(i,2,p,l;m) Y_l-1;mg(i,p,m) ]
Cr         : where f1 = (df/dr - l*f/r) and f2 = (df/dr + (l+1)*f/r)
Cr
Cr Remarks
Cr   ylmbyr returns coefficients to the gradient grad[f(r)*Y_L], apart from a radial contribution.
Cr   The Y_L here are real harmonics (real harmonic polynomials/r^l)
Cr   These coefficients also apply to the position operator \vec{r}*[f(r)*Y_L] (see below).
Cr
Cr   cgcof omits a factor from the radial dependence of f.   Let g(r) = r*f(r).
Cr   Including the radial dependence, the coefficients mapping grad f Y_L to Y_L+1 are
Cr     grad_p (f(r) Y_L) = f1 sum_L c+_pL Y_L+1  + f2 sum_L c-_pL Y_L-1
Cr     where f1(r) = (df/dr - l*f/r) = (dg/dr - (l+1)*g/r)/r
Cr           f2(r) = (df/dr + (l+1)*f/r) = (dg/dr + l*g/r)/r
Cr   The radial-independent coefficients are c+_pL = cgcof(:,1,p,L)  and c-_pL = cgcof(:,2,P,L)
Cr   There are 2 c+_pL and 2 c+_pL; which L they correspond to is given by mg (see above).
Cr   Special cases:
Cr     f(r) = r^l.     (df/dr - l*f/r) = 0  and (df/dr + (l+1)*f/r) = (2l+1)/r
Cr     f(r) = r^-l-1.  (df/dr - l*f/r) = -(2l+1)/r  and (df/dr + (l+1)*f/r) = 0
Cr     f(r) = 1.       (df/dr - l*f/r) = -l/r  and (df/dr + (l+1)*f/r) = (l+1)/r
Cr
Cr   Examples:
Cr    *coefficient  -l * c(1,1,1,2) : gradx Y2
Cr     gradx Y2  = gradx Y_1,-1 = sqrt(3/4pi) gradx (y/r)
Cr               = sqrt(3/4pi) (-xy/r^3)|r=1 = sqrt(3/4pi) * (sqrt(4pi/15) Y_2,-2
Cr               = - sqrt(1/10) Y_2,-2 = - 0.447214 Y_2,-2 =  -l * c(1,1,1,2)
Cr     Same result for  grady Y4 = -l * c(1,1,2,4) = - 0.447214 Y_2,-2
Cr
Cr    *coefficient  -l * c(2,1,1,3) : gradx Y3
Cr     gradx Y3  = gradx Y_1,0 = sqrt(3/4pi) gradx (z/r)
Cr               = sqrt(3/4pi) (-xz/r^3)|r=1 = sqrt(3/4pi) * (sqrt(4pi/15) Y_2,1
Cr               = - sqrt(1/10) Y_2,1 = - 0.447214 Y_2,1 =  -l * c(2,1,1,3)
Cr     Similar result for grady Y3 =              = - 0.447214 Y_2,1 =  -l * c(2,1,1,3)
Cr     grady Y3  = - sqrt(1/10) Y_2,-1 (i.e. ~yz) = - 0.447214 Y_2,-1 = -l * c(1,1,2,3)
Cr     Same result for  gradz Y2 = -l * c(1,1,3,2) = - 0.447214 Y_2,-1
Cr     Similar result for  gradz Y4 = -l * c(1,1,3,4) = - 0.447214 Y_2,1
Cr
Cr    *coefficient  (l+1) * c(1,2,1,4) = coff to gradx Y_1,1 = coff to grady Y_1,-1
Cr     gradx Y4  = gradx Y_1,1 = sqrt(3/4pi) gradx (x/r)
Cr               = sqrt(3/4pi)/r gradx (x) + sqrt(3/4pi) x gradx (1/r) | r=1
Cr     The second term yields coeffs to Y2,0 and Y2,2, which ignore.  1st term:
Cr               = sqrt(3/4pi) = sqrt(3) Y_0,0 = (?)  (l+1) * c(1,2,1,4)
Cr               Should be sqrt(4/3) Y_0,0 ... didn't investigate why.
Cr     Similar result for grady Y2 : has a component 1.154701 Y_0,0 =  (l+1) * c(1,2,2,2)
Cr     Similar result for gradz Y3 : has a component 1.154701 Y_0,0 =  (l+1) * c(1,2,3,3)
Cr     Uncomment printout in this routine to display or the test case below to check C_K
Cr
Cr   Relation position operator \vec{r}*f(r)*Y_L to gradient operator grad(f*Y_L)
Cr   Use this relation (note that r_p means p component of \vec{r})
Cr   grad (r^2 f Y_L) = grad (r^2) f Y_L + r^2 grad (f Y_L) = 2*r (f Y_L) + r^2 grad (f Y_L)
Cr   Therefore r_p*(f Y_L) = (grad (r^2 (f Y_L)) - r^2 grad ((f Y_L)))/2
Cr                         = 1/2 [grad, r^2] (f Y_L)
Cr   Work out coefficients C+ and C- for grad (r^2 (f Y_L)) and for (r^2 grad ((f Y_L)))
Cr     (A)  grad_p (r^2 f Y_L) = (d(r^2*f)/dr - l*r^2*f/r) sum_L c+_pL Y_L+1 +
Cr                             + (d(r^2*f)/dr + (l+1)*r^2*f/r) sum_L c-_pL Y_L-1
Cr     (B)  r^2 grad_p (f Y_L) = r^2*(df/dr - l*f/r) sum_L c+_pL Y_L+1 +
Cr                             + r^2*(df/dr + (l+1)*f/r) sum_L c-_pL Y_L-1
Cr     (grad_p (r^2 (f Y_L)) - r^2 grad_p ((f Y_L)))/2 = ((A) - (B))/2
Cr   The multiplicative terms (proportional to f/r) cancel in the commutator (A) - (B),
Cr   leaving only
Cr     (grad_p (r^2 (f Y_L)) = (d(r^2*f)/dr - r^2*(df/dr)) sum_L c+_pL Y_L+1 / 2
Cr                           + (d(r^2*f)/dr - r^2*(df/dr)) sum_L c-_pL Y_L-1 / 2
Cr                           = r*f * sum_L c+_pL Y_L+1 + r*f * sum_L c-_pL Y_L-1
Cr   Therefore the gradient and position operators are related as
Cr     grad_p (f(r) Y_L) = (df/dr - l*f/r) sum_L c+_pL Y_L+1  + (df/dr + (l+1)*f/r) sum_L c-_pL Y_L-1
Cr     r_p (f(r) Y_L)    = ( r*f         ) sum_L c+_pL Y_L+1  + ( r*f             ) sum_L c-_pL Y_L-1
Cr
Cr   The (mg,cgcof) pair contains the indexing and coefficients
Cr   needed to construct C+ and C- in the Questaal manuscript
Cu Updates
Cu   11 Apr 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax
      integer mg(2,3,(lmax+1)**2)
      double precision cgcof(2,2,3,(lmax+1)**2)
C ... Local parameters
      integer ilm,kx1,kx2,ky1,ky2,kz,l,m,kav,am
      double precision cx1,cx2,cy1,cy2,cz !,f,f1,f2
      integer nlm1,im
C      integer lav
C      procedure(integer) :: ll

      call dpzero(cgcof,size(cgcof))
      call dpzero(cgcof,size(cgcof))
      nlm1 = lmax*lmax
      ilm = 0
      do  l = 0, lmax
        do  m = -l, l
          ilm = ilm+1
          call scglp1(-ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)

          kav = (l+1)*(l+1)+(l+1)+1
          am = iabs(m)

C         Y_l+1m' components of grad Y_l,m
          cgcof(1,1,1,ilm) = cgcof(1,1,1,ilm) + cx1
          cgcof(2,1,1,ilm) = cgcof(2,1,1,ilm) + cx2
          cgcof(1,1,2,ilm) = cgcof(1,1,2,ilm) + cy1
          cgcof(2,1,2,ilm) = cgcof(2,1,2,ilm) + cy2
          cgcof(1,1,3,ilm) = cgcof(1,1,3,ilm) + cz

          mg(1,1,ilm) = m-1
          mg(2,1,ilm) = m+1
          mg(1,2,ilm) = -m-1
          mg(2,2,ilm) = -m+1
          mg(1,3,ilm) = m
          mg(2,3,ilm) = m

C         Y_lm' components of grad Y_l+1,m.  ky1,ky2=index to grad Y_l+1,m l+2 = ll(kx1)+1
          if (ilm > nlm1) cycle

          im = 1; if (kx1-kav<m) im = 2
          cgcof(im,2,1,kx1) = cgcof(im,2,1,kx1) + cx1
          im = 1; if (kx2-kav<m) im = 2
          cgcof(im,2,1,kx2) = cgcof(im,2,1,kx2) + cx2

C         Let mk = m[ky1]=ky1-kav.
C         -mk-1 couples to m (im=1) or -mk+1 couples to m (im=2)
          im = 1; if (kav-ky1+1 == m) im = 2
          cgcof(im,2,2,ky1) = cgcof(im,2,2,ky1) + cy1
          im = 1; if (kav-ky2+1 == m) im = 2
          cgcof(im,2,2,ky2) = cgcof(im,2,2,ky2) + cy2

          cgcof(1,2,3,kz)  = cgcof(1,2,3,kz)  + cz

        enddo
      enddo

C ... Debugging printout
C      call info0(1,1,0,'  C_k for gradx YL%21fk=l+1%32fk=l-1')
C      call info0(1,0,0,'  L = l  m ilm     k,mk1,mk2   l+1,m-1     l+1,m+1     k,mk1,mk2     l-1,m-1     l-1,m+1')
C      do  ilm = 1, nlm1
C        l = ll(ilm); lav = l*l+l+1; m = ilm-lav
C        call info5(1,0,0,'    %3,3i     %3,3i%2;12,6D    %3,3i%2;12,6D',[l,m,ilm],[l+1,m-1,m+1],
C     .    -l*cgcof(:,1,1,ilm),[l-1,m-1,m+1],(l+1)*cgcof(:,2,1,ilm))
C      enddo
C      call info0(1,0,0,'  C_k for grady YL')
C      call info0(1,0,0,'  L = l  m ilm     k,mk1,mk2   l+1,-m-1    l+1,-m+1    k,mk1,mk2    l-1,-m-1    l-1,-m+1')
C      do  ilm = 1, nlm1
C        l = ll(ilm); lav = l*l+l+1; m = ilm-lav
C        call info5(1,0,0,'    %3,3i     %3,3i%2;12,6D    %3,3i%2;12,6D',[l,m,ilm],[l+1,-m-1,-m+1],
C     .    -l*cgcof(:,1,2,ilm),[l-1,-m-1,-m+1],(l+1)*cgcof(:,2,2,ilm))
C      enddo
C      call info0(1,0,0,'  C_k for gradz YL')
C      call info0(1,0,0,'  L = l  m ilm      k,mk%8fl+1,m%20fk,mk%8fl-1,m')
C      do  ilm = 1, nlm1
C        l = ll(ilm); lav = l*l+l+1; m = ilm-lav
C        call info5(1,0,0,'    %3,3i     %2,3i   %1;12,6D%12f     %2,3i  %1;12,6D',[l,m,ilm],[l+1,m],
C     .    -l*cgcof(1,1,3,ilm),[l-1,m],(l+1)*cgcof(1,2,3,ilm))
C      enddo
C      stop

      end
      subroutine mkcga(nlmm,nlmk,nlml,indxcg,jcg,cg,cga)
C- Copy Clebsch-Gordan (Gaunt) coefficients C_KLM to uncompressed array format
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlmk  :Maximum K; also dimensions cga
Ci   nlml  :Maximum L; also dimensions cga
Ci   nlmm  :Leading dimension of cga.
Ci         :Must be dimensioned at least as large as maximum M
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Co Outputs
Co   cga   :cga(M,K,L) is Gaunt coefficient C_KLM (MSM definition)
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlmk,nlml,nlmm,indxcg(*),jcg(*)
      double precision cg(*),cga(nlmm,nlmk,nlml)
C ... Local parameters
      integer K,L,M,ii,indx,ic1,ic2,icg

      do  K = 1, nlmk
        do  L = 1, nlml
          ii = max0(K,L)
          indx = (ii*(ii-1))/2 + min0(K,L)
          ic1 = indxcg(indx); ic2 = indxcg(indx+1)-1
          do  icg = ic1, ic2
            M = jcg(icg)
            cga(M,K,L) = cg(icg)
          enddo
        enddo
      enddo
      end

      subroutine mkcgar(mode,nlmm,nlmk,nlml,indxcg,jcg,cg,cgar)
C- Expansion coefficients to YK YL (x,y,z) as coefficients to spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do not distinguish l+1 and l-1 contributions to Y1p YK in cgar
Ci         :1 do     distinguish l+1 and l-1 contributions to Y1p YK in cgar
Ci   nlmk  :Maximum K; also dimensions cga
Ci   nlml  :Maximum L; also dimensions cga
Ci   nlmm  :Leading dimension of cga.
Ci         :Must be dimensioned at least as large as maximum M,
Ci         :which is [(ll(nlml)+ll(nlmk)+1)+1]^2.
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Co Outputs
Co   cgar  :cgar(M,K,L,p,pmi) is coefficient to YM in YM expansion of (x,y,z) YK YL
Co         :p = 1,2,3 for x,y,z
Co         :pmi = 1 if mode=0; pmi loops between 1 and 2 if mode = 1
Cl Local variables
Cl   cg1cof and mj:
Cl         : cg1cof(1:2,1:2,p,klm) = the family of Gaunt or CG coefficients that expand
Cl         : Y_1,p Y_k,mk as linear combinations of Y_k+1m' and Y_k-1,m' for p = one of (x,y,z).
Cl         : Two coefficients represent the Y_k+1,m' terms and two the Y_k-1,m' terms.
Cl         : mj(1:2,p,L) = corresponding m indices m'.  That is
Cl         : Y_1,p Y_k,mk   sum_i=1,2 [ cg1cof(i,1,p,klm) Y_k+1;mj(i,p,klm) ]
Cl         :              + sum_i=1,2 [ cg1cof(i,2,p,klm) Y_k-1;mj(i,p,klm) ]
Cl         : klm = kav+mk  where kav = k*k+k+1
Cr Remarks
Cu Updates
Cu   19 May 19 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nlmk,nlml,nlmm,indxcg(*),jcg(*)
      double precision cg(*),cgar(nlmm,nlmk,nlml,3,2)
C ... Local parameters
      integer jav,klm,llm,mlm,lmaxk,lmaxl,p,pm,pmi,mp,j,jlm
      integer lmax,mmax,ml,mm
      integer, allocatable :: mj(:,:,:)
      real(8), allocatable :: cga(:,:,:),cg1cof(:,:,:,:)
      procedure(integer) :: ll

C     mode0 = mod(mode,10)
      lmaxk = ll(nlmk); lmaxl = ll(nlml)
      lmax = lmaxk+lmaxl
      mmax = min(ll(nlmm),lmax+1)**2
      allocate(cga(nlmm,nlmk,nlml),cg1cof(2,2,3,mmax),mj(2,3,mmax))
      call ylmbyr(lmax,cg1cof,mj)
      call mkcga(nlmm,nlmk,nlml,indxcg,jcg,cg,cga)

      call dpzero(cgar,size(cgar))

      do  p = 1, 3

        mlm = 0
        do  ml = 0, lmax
        do  mm = -ml, ml
          mlm = mlm+1

C         Loop over four possible Y_j,mj that expand Y_1p Y_k;mk
          pmi = 1               ! Default case mode=0 : do not distinguish l+1 and l-1 contr
          do  pm = 1, 2         ! pm=1,2 corresponds to j = ml+1 and j = ml-1
            j = ml + 3 - 2*pm   ! j = ml+1 or  ml-1
            if (j < 0) cycle    ! No contribution when j<0
            jav = j*j+j+1
            do  mp = 1, 2        ! pair of indices m' == mj(mp,p,mlm)
              if (p == 3 .and. mp == 2) cycle   ! z component has only a single m'
              if (cg1cof(mp,pm,p,mlm) == 0) cycle ! Skip if coefficient vanishes
              if (mode > 0) pmi = pm  ! when distinguishing raising, lowering operators
              jlm = jav+mj(mp,p,mlm)
              if (jlm > nlmm) cycle
C             if (jlm == 4 .and. p == 1) call snot(1)
C             if (mlm == 9 .and. p == 1) call snot(1)
              forall (klm = 1:nlmk, llm = 1:nlml)
     .          cgar(jlm,klm,llm,p,pmi) = cgar(jlm,klm,llm,p,pmi) + cga(mlm,klm,llm)*cg1cof(mp,pm,p,mlm)
C             call snot(-1)

            enddo ! mp
          enddo ! pm

          enddo ! mm
          enddo ! ml

        enddo  ! p

      deallocate(cga,cg1cof,mj)

      end

C#ifdefC TEST
C
CC a.out [--gradcg] [--mode=#]
C      subroutine fmain
C      implicit none
C
C      character*(120) strn,strn2,cxyz*3,dc*1
C      integer, parameter :: nnn=300
C      integer mode,ilmx,nlm,nlmf,lmax,j,ilm,l,lav,m,mm,lp,pm,kavp,kavm,nf1,nf2
C      integer k,llm,jlm,klm,mk,ml,nlmm
C      integer np,nph,nth
C      double precision pi,srfpi,err(3),xx1
C      double precision p(3,nnn),wp(nnn)
C      integer, allocatable :: mg(:,:,:)
C      real(8) flnum(16,16,3),flanl(16,16,3)
C      real(8), allocatable :: xp(:),yp(:),zp(:),rp(:),r2(:)
C      real(8), allocatable :: yl(:,:),gyl(:,:,:),frl(:)
C      real(8), allocatable :: grp(:,:),ggrp(:,:),coff(:,:,:,:)
C      real(8), allocatable :: rgrad(:,:,:),gradm(:,:,:,:,:)
C
CC     For Clebsch Gordan coefficients
C      integer lmxcg,lmxcy,lnxcg,lnjcg
C      integer, allocatable :: jcg(:),indxcg(:)
C      real(8), allocatable :: cg(:),cy(:),cga(:,:,:),cgar(:,:,:,:,:)
C      procedure(logical) :: cmdopt
C      procedure(integer) :: ll,cmdoptsw,a2vec
C      integer, parameter :: NULLI=-99999
C      data cxyz / 'xyz'/
C
CC     mode=4  => decompose r^-l grad*(r^l Ylm) into linear combination of Ylm, Ylm = sph. harm. (no r^l)
CC     mode=14 => decompose      grad*(Ylm) into linear combination of Ylm
CC     mode=24 => decompose r^(l+1) grad*(r^(-l-1) Ylm) into linear combination of Ylm
C
C      mode =  4
C      mode = 14
CC      mode = 24
C      if (cmdopt('--mode=',NULLI,0,strn)) then
C        dc = '='
C        j = 7
C        if (a2vec(strn,len_trim(strn),j,2,', '//dc,3,2,1,ilm,mode) < 1) call rx('failed to parse '//trim(strn))
C      endif
C      if (mode/=4 .and. mode/=14 .and. mode/=24) call rx('bad mode')
C
C      ilmx = 9
CC     ilmx = 2
C      lmax = 4
C      nlmf = (lmax+2)**2
C      nlm  = (lmax+1)**2
C      pi = 4*datan(1d0)
C      srfpi = dsqrt(4*pi)
C
CC     Setup CG coefficients
C      lmxcg = 9; lmxcy = 12
C      call scg0(0,lmxcg,err,err,err,lnjcg,lnxcg)
C      allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg),cy((lmxcy+1)**2))
C      call sylmnc(cy,lmxcy)
C      call scg(lmxcg,cg,indxcg,jcg)
C      allocate(cga((2*lmax+1)**2,nlm,nlm))
C      call mkcga((2*lmax+1)**2,nlm,nlm,indxcg,jcg,cg,cga)
C      nth = -122
C      call fpiint(nth,nph,np,p,wp)
C      allocate(xp(np),yp(np),zp(np),r2(np),rp(np),grp(np,3),ggrp(np,3))
C      allocate(frl(nlm),yl(np,nlmf),gyl(np,nlm,3))
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
C      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
C
C      nlmm = (2*lmax+2)**2
C      allocate(cgar(nlmm,nlm,nlm,3,2))
C      call mkcgar(1,nlmm,nlm,nlm,indxcg,jcg,cg,cgar)
C
C      if (cmdopt('--gradcg',NULLI,0,strn)) then
C
C      call info2(1,1,0,' Check structure matrix Y1 YK YL against num int (%i angular points)',np,mode)
C
C      err = 0
C      do  j = 1, 3
C        jlm = 1+1 + mod(j+1,3)  ! x for j=1, y for j=2,  z for j=3
C        print *
C        print *, "  k mk    l ml klm llm   L       "//cxyz(j:j)//" YK YL ..."
C        klm = 0
C        do  k = 0, 2
C        do  mk = -k, k
C        klm = klm+1
C
C          llm = 0
C          do  l = 0, 2
C          do  ml = -l, l
C            llm = llm+1
C
CC           YK YL Y1 at each point in np
C            if (abs(yl(1,jlm)*sqrt(4*pi/3)-p(j,1)) .gt. 1d-10) stop 'oops'
C            gyl(:,1,1) = yl(:,klm)*yl(:,llm) * yl(:,jlm)*sqrt(4*pi/3)
C            call xyl(nlm,np,wp,gyl,r2,yl,frl)
C
C            xx1 = sum(abs(cgar(1:nlm,llm,klm,j,1)-frl(1:nlm)))
C            xx1 = sum(abs(cgar(1:nlm,llm,klm,j,1)+cgar(1:nlm,llm,klm,j,2)-frl(1:nlm)))
C            err(j) = err(j) + xx1
C            call awrit6(' Yl(%i%,2i) Yl(%i%,2i)%,3i%,3i: ',strn,len(strn),0,l,ml,k,mk,klm,llm)
C            strn2 = '%120pnum-CG = '
C            call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C
C            call prylx(0,strn,llm,0,frl,nlm,1d0,strn2)
C
C
C          enddo
C          enddo
C
C        enddo
C        enddo
C
C      enddo
C
C      call info2(1,0,0,' cumulative diff between (analytic, num CG * (x,y,z)) = %3:1;3,3g',err,2)
C      deallocate(xp,yp,zp,r2,rp,grp,ggrp,frl,yl,gyl)
C      call rx0('done')
C
C      endif
C
CC --- CG coefficients of real harmonics by numerical integration ---
C      call info2(1,0,0,' %i angular points.  mode=%i'//
C     . '%-1j%?#n==4# (gradients of r^l Y_L, Y_L = real harmonics)##'//
C     . '%-1j%?#n==14# (gradients of real harmonics Y_L)##'//
C     . '%-1j%?#n==24# (gradients of r^-l-1 Y_L, Y_L = real harmonics)##',
C     .  np,mode)
C
C      allocate(coff(2,2,3,(lmax+1)**2),mg(2,3,(lmax+1)**2))
C
C      lp=0; if (mod(mode/10,10) == 1) lp=1; if (mod(mode/10,10) == 2) lp=-1
C      call ropylg(lp,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)  ! gradient and laplacian of spherical harmonic polynomials
C      call info0(1,1,0,
C     .  " ... CG coefficients C*Y_l'm' by tabulating grad*Ylm on angular and integrating (grad*Ylm) Yl'm'"//
C     .  "%N%18fdl d+-m  coff ...")
C
C      do  j = 1, 3
C      do  ilm = 1, 16
C        l = ll(ilm)
C        lav = l*l+l+1
C        m = ilm-lav
C        if (m == -l) call info0(1,1,0,'')
C        call xyl(nlm,np,wp,gyl(1,ilm,j),r2,yl,frl)
C        call awrit2(' '//cxyz(j:j)//' grad Yl(%i%,2i): ',strn,len(strn),0,l,m)
C        pm = 1; if (j == 2) pm = -1
C        call pryl(0,strn,ilm,pm,frl,nlm,1d0,'')
C        call dcopy(16,frl,1,flnum(1,ilm,j),1)
C      enddo
C      enddo
C
C      call info0(1,1,0,"... analytic CG coefficients C*Y_l'm' to grad*Ylm by calling ylmbyr%N%18fdl d+-m  coff ...")
C
C      call ylmbyr(lmax,coff,mg)
C
CC     f(r) = 1 and r = 1 => (df/dr - l*f/r) = -l  and (df/dr + (l+1)*f/r) = l+1
C      do  ilm = 1, lmax**2
C        l = ll(ilm)
C        if (mode == 4) then  ! grad r^l Y_L
C          coff(:,1,:,ilm) = 0*coff(:,1,:,ilm)
C          coff(:,2,:,ilm) = (2*l+1)*coff(:,2,:,ilm)
C        elseif (mode == 14) then
C          coff(:,1,:,ilm) = -l*coff(:,1,:,ilm)
C          coff(:,2,:,ilm) = (l+1)*coff(:,2,:,ilm)
C        elseif (mode == 24) then ! grad r^-l-1 Y_L
C          coff(:,1,:,ilm) = -(2*l+1)*coff(:,1,:,ilm)
C          coff(:,2,:,ilm) = 0*coff(:,2,:,ilm)
C        endif
C      enddo
C
C      err = 0
C      do  j = 1, 3
C      do  ilm = 1, lmax**2
C        l = ll(ilm)
C        lav = l*l+l+1
C        kavp = (l+1)*(l+1)+(l+1)+1
C        kavm = (l-1)*(l-1)+(l-1)+1 ! k- restricted to kmin=l^2+1 to kmax=(l+1)**2
C        m = ilm-lav
C        mm = iabs(m)
C        if (m == -l) call info0(1,1,0,'')
C        call dpzero(frl,size(frl))
C        pm = 1
C        if (j == 1) then
C          frl(kavp+m-1) = coff(1,1,1,ilm)
C          frl(kavp+m+1) = coff(2,1,1,ilm)
CC         print *, m-l>=0, 1-m<=l, m-l - ((kavm+m-1) - l**2) !,  '   ', (l+1)**2 - (kavm+m+1),  (- m)+3*l-1
CC         print *, 3*l-m>=1, 1+m<=l, (l+1)**2 - (kavm+m+1) - (3*l-m-1)
C          if (m-l - ((kavm+m-1) - l**2) /= 0) stop 'oops'
C          if ((l+1)**2 - (kavm+m+1) - (3*l-m-1) /= 0) stop 'oops'
C          if (l>0) then
C          if (+m-1 >= 1-l .and. +m-1 <= l-1) frl(kavm+m-1) = coff(1,2,1,ilm)
C          if (+m+1 >= 1-l .and. +m+1 <= l-1) frl(kavm+m+1) = coff(2,2,1,ilm)
C          endif
CC         Debugging
CC          xx1=0; if (kavm+m-1>0) xx1=coff(1,2,1,ilm)
CC          xx2=0; if (kavm+m+1>0) xx2=coff(2,2,1,ilm)
CC          call info8(1,0,0," l,m,ilm =%3,3i   l',m1',ilm1' =%3,3i%;11,6D   l',m2',ilm2' =%3,3i%;11,6D",
CC     .      [l,m,ilm],[l-1,m-1,kavm+m-1],xx1,[l-1,m+1,kavm+m+1],xx2,6,7,8)
C        endif
C        if (j == 2) then
CC         coff(1,1,2,ilm) = term C*Y(l+1,-mk-1)
CC         coff(2,1,2,ilm) = term C*Y(l+1,-mk+1)
CC         coff(1,2,2,ilm) = term C*Y(l-1,-mk-1)
CC         coff(2,2,2,ilm) = term C*Y(l-1,-mk+1)
C          frl(kavp-m-1) = coff(1,1,2,ilm)
C          frl(kavp-m+1) = coff(2,1,2,ilm)
C          if (l>0) then
C          if (-m-1 >= 1-l .and. -m-1 <= l-1) frl(kavm-m-1) = coff(1,2,2,ilm)
C          if (-m+1 >= 1-l .and. -m+1 <= l-1) frl(kavm-m+1) = coff(2,2,2,ilm)
C          endif
C          pm = -1
C        endif
C        if (j == 3) then
C          frl(kavp+m) = coff(1,1,3,ilm)
C          if (kavm+m>0) frl(kavm+m) = coff(1,2,3,ilm)
C        endif
C        call awrit2(' '//cxyz(j:j)//' grad Yl(%i%,2i): ',strn,len(strn),0,l,m)
C        pm = 1; if (j == 2) pm = -1
C        call pryl(0,strn,ilm,pm,frl,nlm,1d0,'')
C        call dcopy(lmax**2,frl,1,flanl(1,ilm,j),1)
C        err(j) = err(j) + sum(abs(flnum(1:lmax**2,ilm,j)-flanl(1:lmax**2,ilm,j)))
C      enddo
C      enddo
C
C      call info2(1,1,0,
C     .  ' cumulative diff between (analytic, num diff) ME (x,y,z) = %3:1;3,3g',err,2)
CC     call prmx('difference (num-analytical) projection of grad Yl to out.dat',flnum-flanl,lmax**2,lmax**2,lmax**2*3)
C
CC ... Same test but call rgrme
C      call info0(1,1,0,
C     .  " ... repeat test, calling grmes for matrix elements.  Compare diff (grme,inline)"//"%N%19fl  m    coff ...")
CC     Same result with call to grmes
C      nf1 = 1; nf2 = 1
C      allocate(gradm(nlm,nlm,3,nf1,nf2))
C
C      allocate(rgrad(2,2,0:lmax)); call dpzero(rgrad,size(rgrad))
C      do  l = 0, lmax-1
C        if (mode ==  4) rgrad(2,2,l) =  2*l+1   ! for (l-1, l)
C        if (mode == 14) rgrad(1,1,l) =  -l      ! for (l+1, l)
C        if (mode == 14) rgrad(2,2,l) =  l+1     ! for (l-1, l)
C        if (mode == 24) rgrad(1,1,l) = -2*l-1   ! for (l+1, l)
C      enddo
C
C      call grmes(20,lmax,rgrad,1,1,[0],[0],nlm,gradm,gradm)
C
C      err = 0
C      do  j = 1, 3
C      do  ilm = 1, 16
C        l = ll(ilm)
C        lav = l*l+l+1
C        m = ilm-lav
C        call awrit2(' '//cxyz(j:j)//' grad Yl(%i%,2i): ',strn,len(strn),0,l,m)
C        call pryl(0,trim(strn),ilm,0,gradm(1:lmax**2,ilm,j,1,1)-flanl(1:lmax**2,ilm,j),lmax**2,1d0,'')
C        err(j) = err(j) + sum(abs(gradm(1:lmax**2,ilm,j,1,1)-flanl(1:lmax**2,ilm,j)))
C      enddo
C      print *
C      enddo
C      call info2(1,0,1,' cumulative diff between (grme, inline) ME (x,y,z) = %3:1;3,3g',err,2)
C
C
C      call cexit(1,1)
C      end
C
CC      subroutine snot(i)
CC      print *, i
CC      end
C
C      subroutine xyl(nlm,np,wp,fp,r2,yl,fl)
CC- Yl-projection of function tabulated on an angular mesh
C      implicit none
C      integer nlm,np,ip,ilm
C      double precision fl(nlm),r2(np),fp(np),yl(np,nlm),wp(np)
C      double precision rl
C      procedure(integer) :: ll
C
C      call dpzero(fl,nlm)
C      do  ip = 1, np
C      do  ilm = 1, nlm
C        rl = dsqrt(r2(ip))**ll(ilm)
C        fl(ilm) = fl(ilm) + fp(ip)*wp(ip)*yl(ip,ilm)/rl
C      enddo
C      enddo
C      end
C      subroutine pryl(mode,strn,mlm,pm,fl,nlm,fac,endstr)
CC- Prints out nonzero (l,m) components of fl(ilm)*fac
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 print fac*fl
CCi         :1 print fac*sign(fl)*fl**2
CCi         :1 print fac*sign(fl)*fl**2*(2l+1)*l
CCi   strn
CCi   mlm   :ilm from which this function was derived
CCi   pm    :1   delta m = m(ilm) - m(fl)
CCi         :-1  delta m = m(ilm) - m(fl)
CCi         :0   print l and m not delta l and delta m
CCi   fl    :function as linear combination of Ylm
CCi   nlm   :print coffs to fl up to nlm
CCi   fac   :scale fl by fac
CCo Outputs
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   06 Apr 18
CC ----------------------------------------------------------------------
C      implicit none
C      character*(*) strn,endstr
C      integer mode,mlm,nlm,pm
C      double precision fl(nlm),fac
C      character*120 outs
C      integer ilm,l2,m2,lmax,l,m,lav,dl,dm
C      double precision f
C      real(8),parameter:: tol=1d-6
C      procedure(integer) :: ll
C
C      l = ll(mlm)
C      lav = l*l+l+1
C      m = mlm-lav
C      call awrit2(strn,outs,len(outs),0,l,m)
C
C      lmax = ll(nlm)
C      ilm = 0
C      do  l2 = 0, lmax
CC        ilm0 = l2**2 + l2 + 1   ! the m=0 element for this l
CC        ilmm = ilm0 - 2*l2      ! the m=0 element for this l-1
CC        ilmp = ilm0 + 2*(l2+1)  ! the m=0 element for this l+1
C        do  m2 = -l2, l2
C          ilm = ilm+1
C          if (abs(fl(ilm)) > tol) then
C            if (mode == 0) f = fl(ilm)*fac
C            if (mode == 1) f = fl(ilm)**2*fac*dsign(1d0,fl(ilm))
C            dl = l2-l; if (pm==0) dl = l2
C            dm = m2-m; if (pm<0) dm=m2+m; if (pm==0) dm = m2
C            call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,f)
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l-1))
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l+3))
C
CC           if (l2-l /= -1) stop 'oops'  ! True for grad r^l YL, YL = SH (not polynomials)
CC           if (l2-l /= 1) stop 'oops'  ! True for grad r^-l-1 YL, YL = SH (not polynomials)
C
C          endif
C        enddo
C      enddo
C
C      if (endstr /= ' ') outs = trim(outs) // ' ' // trim(endstr)
C
C      call info0(1,0,0,trim(outs))
C
C      end
C      subroutine prylx(mode,strn,mlm,pm,fl,nlm,fac,endstr)
CC- Prints out nonzero (l,m) components of fl(ilm)*fac
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 print fac*fl
CCi         :1 print fac*sign(fl)*fl**2
CCi         :1 print fac*sign(fl)*fl**2*(2l+1)*l
CCi   strn
CCi   mlm   :ilm from which this function was derived
CCi   pm    :1   delta m = m(ilm) - m(fl)
CCi         :-1  delta m = m(ilm) - m(fl)
CCi         :0   print l and m not delta l and delta m
CCi   fl    :function as linear combination of Ylm
CCi   nlm   :print coffs to fl up to nlm
CCi   fac   :scale fl by fac
CCo Outputs
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   06 Apr 18
CC ----------------------------------------------------------------------
C      implicit none
C      character*(*) strn,endstr
C      integer mode,mlm,nlm,pm
C      double precision fl(nlm),fac
C      character*160 outs
C      integer ilm,l2,m2,lmax,l,m,lav,dl,dm
C      double precision f
C      real(8),parameter:: tol=1d-6
C      procedure(integer) :: ll
C
C      l = ll(mlm)
C      lav = l*l+l+1
C      m = mlm-lav
C      call awrit2(strn,outs,len(outs),0,l,m)
C
C      lmax = ll(nlm)
C      ilm = 0
C      do  l2 = 0, lmax
CC        ilm0 = l2**2 + l2 + 1   ! the m=0 element for this l
CC        ilmm = ilm0 - 2*l2      ! the m=0 element for this l-1
CC        ilmp = ilm0 + 2*(l2+1)  ! the m=0 element for this l+1
C        do  m2 = -l2, l2
C          ilm = ilm+1
C          if (abs(fl(ilm)) > tol) then
C            if (mode == 0) f = fl(ilm)*fac
C            if (mode == 1) f = fl(ilm)**2*fac*dsign(1d0,fl(ilm))
C            dl = l2-l; if (pm==0) dl = l2
C            dm = m2-m; if (pm<0) dm=m2+m; if (pm==0) dm = m2
C            call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,f)
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l-1))
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l+3))
C
CC           if (l2-l /= -1) stop 'oops'  ! True for grad r^l YL, YL = SH (not polynomials)
CC           if (l2-l /= 1) stop 'oops'  ! True for grad r^-l-1 YL, YL = SH (not polynomials)
C
C          endif
C        enddo
C      enddo
C
C      if (endstr /= ' ') outs = trim(outs) // ' ' // trim(endstr)
C
C      call info0(1,0,0,trim(outs))
C
C      end
C      subroutine prmx(strn,s,ns,nr,nc)
CC- writes matrix into out file (for debugging)
C      implicit none
C      integer nr,nc,ns,ifi
C      double precision s(ns,nc,2)
C      character*(14) fmt, fmt0, strn*(*), outs*100
C      integer i,j,fopna,i1mach
C      save fmt
C      data fmt /'(9f20.15)'/
C      fmt = '(1p9e20.10)'
C      ifi = fopna('out',29,0)
C      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
C#elseC
C      call awrit2('%% rows %i cols %i real',' ',100,ifi,nr,nc)
C#endifC
C      do  i = 1, nr
C        write (ifi,fmt) (s(i,j,1),j=1,nc)
C      enddo
C      write(ifi,*)
CC      do  12  i = 1, nr
CC   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
C      call fclose(ifi)
C
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#elseC
C      outs = ' prm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-100,-i1mach(2))
C#endifC
C      read(*,'(a100)') outs
C
C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in prmx'
C#elseC
C      if (outs == 'q') call rx0('quit in prmx')
C#endifC
C      return
C
C      entry prmx0(fmt0)
C      fmt = fmt0
C      end
C      subroutine zprm(strn,icast,s,ns,nr,nc)
CC- Print complex matrix to file out
C      implicit none
C      integer icast,nr,nc,ns,ifi
C      double precision s(2,ns,nc)
C      character*(20) fmt, outs*80, strn*(*)
C      integer i,j,fopna,i1mach
C      character fmt0*(*)
C      save fmt
C      data fmt /'(9f18.11)'/
C
C      outs = ' '
C      if (icast == 1)  outs = ' real'   ! print real part only
C      if (icast == 11) outs = ' symm'
C      if (icast == 2)  outs = ' complex'
C      if (icast == 12) outs = ' herm'
C      ifi = fopna('out',29,0)
C      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc,outs(1:10)
C#elseC
C      call awrit2('%% rows %i cols %i'//trim(outs),' ',80,ifi,nr,nc)
C#endifC
C      do  i = 1, nr
C        write (ifi,fmt) (s(1,i,j),j=1,nc)
C      enddo
C      if (mod(icast,10) > 1) then
C      write(ifi,*)
C        do  i = 1, nr
C          write (ifi,fmt) (s(2,i,j),j=1,nc)
C        enddo
C      endif
C      call fclose(ifi)
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#elseC
C      outs = ' zprm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endifC
C      read(*,'(a80)') outs
C
C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in zprm'
C#elseC
C      if (outs == 'q') call rx0('quit in zprm')
C#endifC
C
C      return
C
C      entry zprm0(fmt0)
C      fmt = fmt0
C      end
C
C#endif
