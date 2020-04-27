      subroutine rgrme(job,nf,nf1,nf2,nsp,nspc,lmxs,lmx,nr,ri,wi,phi1,phi2,rgrad)
C- Matrix elements of radial wave functions, radial part of gradient operator
C ----------------------------------------------------------------
Ci Inputs
Ci   job   :0,  matrix elements of grad operator
Ci         :1   matrix elements of pos operator
Ci   nf    :Leading dimensions of rgrad
Ci   nf1   :Number of radial fns g1 for which to calc. <g2 | grad g1>
Ci   nf2   :Number of radial fns g2 for which to calc. <g2 | grad g1>
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmxs  :Used to dimension rgrad
Ci   lmx   :make integrals for l=0:lmx
Ci   nr    :number of radial mesh points
Ci   ri    :radial mesh points
Ci   wi    :radial mesh weights for int dr f(r) (r^2 should to be folded into f)
Ci   phi1  :family of g1 for l=0:lmx, isp=1,nsp, if1=1:nf1
Ci         :g1 = r*(large component of radial wf)
Ci   phi2  :family of g2 for l=0:lmx, isp=1,nsp, if2=1:nf2
Ci         :g2 = r*(large component of radial wf)
Co Outputs
Co   job=0
Co   rgrad(if1,if2,ll,l,nsp):  <g2_l' grad_r g1_l> with l' = l +/- 1
Co         :if1 = index to ket partial wave g1
Co         :if2 = index to bra partial wave g2
Co         :Example: if1 = 1,2 for phi,phidot:
Co         :    if1   if2  g1      g2
Co         :     1     1   phi     phi
Co         :     2     1   phidot  phi
Co         :     1     2   phi     phidot
Co         :     2     2   phidot  phidot
Co         :ll= 1 => l' = l+1  <g2(l') | grad_r g1(l)> - (l+1) < g2(l') | 1/r g1(l) >
Co         :    2 => l' = l-1 <g2(l') | grad_r g1(l)> +    l  < g2(l') | 1/r g1(l) > 
Co
Co         :rgrad(:,:,1,:,:):    lp = l+1, f+(r) = (g2(lp) grad_r g1(l)) - (l+1) (g2(lp) 1/r g1(l))
Co         :rgrad(:,:,2,:,:):    lp = l-1, f-(r) = (g2(lp) grad_r g1(l)) +    l  (g2(lp) 1/r g1(l))
Co         :rgrad(:,:,:,l,:):    l index to g1(l)
Co         :rgrad(:,:,:,:,isp):  spin index (compound index in noncollinear case)
Co         :r * grad [g(r)/r Y_lm] generates two kinds of terms:
Co         : one term raises  l by one and is proportional to dg/dr - (l+1)g/r
Co         : the other lowers l by one and is proportional to dg/dr + (l)g/r
Co         :See https://www.questaal.org/docs/numerics/spherical_harmonics/#gradients-of-spherical-harmonics
Co  Noncollinear case: last index becomes effectively a 2x2 spinor (ispc,isp)
Co
Co   job=1
Co   rgrad(if1,if2,ll,l,nsp):  <g2_l' r g1_l> with l' = l +/- 1
Cl Local variables
Cl   ksp   : isp  in the collinear case, isp = 1 ... nsp
Cl         : ispc in the noncollinear case
Cl   kspc  : isp  in the collinear case, isp = 1 ... nsp
Cl         :      compound spin index (isp,ispc) in the noncollinear case
Cr Remarks
Cr   This routine calculates the radial part of <g2 grad g1> or to to <g2 r g1>
Cr   Angular momentum contribution obtained via Wigner Eckhart theorem;
Cr   see gradme.f.
Cr   Adapted from an optics package by V. Antropov
Cr
Cr   Debugging for position operator.  Uncomment write to file to phi1.  For <phi1 | phi2>
Cr   mcx phi1 -a phi phi -coll 1 phi -coll 2:4 phi -coll 3:5 -xe -ccat  -int 0  2.661766
Cr   For <phi1 | r phi2>
Cr   mcx phi1 -a phi phi -coll 1 phi -e3 'x1*x2' 'x1*x3' 'x1*x4' phi -coll 3:5 -xe -ccat  -int 0 2.661766
Cr Bugs
Cr   rgrad is dimensioned unecessarily large by a factor of 2.  Future:  merge indices 3 and 4.
Cu Updates
Cu   04 May 19 Add ability to calculate rme (altered argument list)
Cu   22 Jun 14 Extended to the noncollinear case
Cu   31 Dec 13 Redesign for fp, with
Cu   05 Jan 04 dimensions of rgrad distinct from lmx
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,nr,nsp,nspc,nf,nf1,nf2,lmxs,lmx
      double precision ri(nr),wi(nr)
      double precision phi1(nr,0:lmx,nsp,nf1),phi2(nr,0:lmx,nsp,nf2)
C     Outputs
      double precision rgrad(nf,nf,2,0:lmxs,nsp*nspc)
C ... Local parameters
      logical lgrad,lpos
      integer lun,isp,loc,ir,ll,if1,if2,ispc,ksp,kspc
      double precision wbr(nr),wir(nr),xv(2)
C     double precision rmex
      procedure(real(8)) :: dot3

C     call prrmsh('phi1',ri,phi1,nr,nr,(lmx+1)*nsp*nf1)
C      call prrmsh('phi',ri,phi1(:,0,1,1),nr,nr,lmx+1)
C      call prrmsh('phid',ri,phi1(:,0,1,2),nr,nr,lmx+1)

C ... Extra weights for radial mesh
      wbr(1) = 0
      forall (ir=2:nr) wbr(ir) = wi(ir)/ri(ir)
      wir(1) = 0
      forall (ir=2:nr) wir(ir) = wi(ir)*ri(ir)
      call dpzero(rgrad,nf*nf*2*(lmxs+1)*nsp*nspc)
      lgrad = job == 0
      lpos  = job /= 0

C --- Loop over spins ---
      do  isp = 1, nsp
      do  ispc = 1, nspc

        ksp = isp ; if (nspc == 2) ksp = ispc
        kspc = isp + 2*(ispc-1)

C   ... Radial part of matrix elements <g2 grad g1>
        do  lun = 0, lmx

C   ... Save radial integrals g1 grad g2 for each l
        do  loc = 0, lmx
          if (abs(loc-lun) /= 1) cycle
          ll = 1
          if (loc-lun == 1) ll = 2

          do  if2 = 1, nf2
          do  if1 = 1, nf1
            if (lgrad) then
              call prgrme(phi1(1,loc,isp,if1),phi2(1,lun,ksp,if2),ri,wi,wbr,nr,loc,xv(1),xv(2))
              rgrad(if1,if2,ll,loc,kspc) = xv(ll)
            endif

            if (lpos) then
C           debugging: matrix elements <phi1 | phi2>
C           rme = dot3(nr,phi1(1,loc,isp,if1),phi2(1,lun,ksp,if2),wi); print *, loc,lun,rme
C           debugging:  matrix elements <phi1 | r phi2>
C           rmex = dot3(nr,phi1(1,loc,isp,if1),phi2(1,lun,ksp,if2),wir); print *, loc,lun,rme
            rgrad(if1,if2,ll,loc,kspc) = dot3(nr,phi2(1,lun,ksp,if2),phi1(1,loc,isp,if1),wir)
C           print *, loc,lun,rme(if1,if2,ll,loc,kspc)
            endif

          enddo
          enddo

        enddo !loc
        enddo !lun
      enddo ! spins
      enddo ! noncollinear spins

C      call prmx('rgrad',rgrad,nf**2,nf**2,(lmxs+1))
C      call prmx('rgrad',rgrad,nf**2*4,nf**2*2,(lmxs+1))
C      call prmx('rme',rme,nf**2,nf**2,2*(lmxs+1))
      end

      subroutine prgrme(g1,g2,ri,wi,wbr,nr,l,r1,r2)
C- Radial integrals of partial waves with the radial gradient operator
C ----------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   g1    :r*partial wave of ket, quantum number l; see Outputs
Ci   g2    :r*partial wave of bra, normally quantum number l+1 or l-1; see Outputs
Ci   ri    :radial mesh points
Ci   wi    :radial mesh weights
Ci   wbr   :wi/r
Ci   nr    :number of radial mesh points
Ci   l     :l quantum number of g1
Co Outputs
Co    r1  = <g2 | grad g1>   - (l+1) < g2 | 1/r g1 >  (matrix element coupling l+1 to l)
Co    r2  = <g2 | grad g1>   +     l < g2 | 1/r g1 >  (matrix element coupling l-1 to l)
Cr Remarks
Cr   Use r1,r2 when g2 and g1 refer to r*partial wave.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,l
      double precision g1(nr),g2(nr),ri(nr),wi(nr),wbr(nr)
      double precision r1,r2
C ... Local parameters
      integer lerr
      double precision sum1,sum2
      double precision dg(nr)
      procedure(real(8)) :: dot3

C     debugging
C     sum1 = dot3(nr,g1,g2,wi); print *, '1',sum1

C     call radgra(a,b,nr,ri,g1,dg)
      call poldvm(ri,g1,nr,6,.false.,1d-12,lerr,dg)
C     call prrmsh('dg',ri,dg,nr,nr,1)
      sum1 = dot3(nr,g2,dg,wi)  ! <g2 | d/dr | g1>
      sum2 = dot3(nr,g2,g1,wbr) ! <g2 | 1/r  | g1>
      r1 = sum1 - sum2*(l+1)
      r2 = sum1 + sum2*l

C     Use this one for g(r) = r*phi(r) = r^l+1
C     print *, 'rgrme for phi(r)=r^l'; r1=0; r2=2*l+1
C     Use this one for g(r) = r*phi(r) = r-l
C     print *, 'rgrme for phi(r)=r^-l-1'; r1=-(2*l+1); r2=0
C     Use this one for g(r) = r*phi(r) = r-l
C     print *, 'rgrme for phi(r)=1'; r1=-l; r2=l+1

      end
      subroutine grmes(mode,lmax,rgrad,nf1,nf2,np1,np2,nd,gradm,zgradm)
C- Matrix elements of partial waves and the gradient operator in one sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :  0 return gradients grad(x,y,z) in real harmonic basis
Ci         :  1 return gradients grad(x,y,z) in spherical harmonic basis, m=-l...l order
Ci         :  2 return gradients grad(x,y,z) in spherical harmonic basis, m=l...-l order
Ci         :  3 return gradients grad(x,y,z) in spherical harmonic basis, old conventions
Ci         :    DOES NOT WORK.
Ci         :10s digit controls form of gradm returned by grmes
Ci         :  0 Return real part in gradm, discarding the imaginary part
Ci         :    Often (e.g. for mode=0 or mode=201 or 202) the result should be real anyway.
Ci         :    In this case gradm = gradm(k,n2,n1,ilm2,ilm1), as discussed in Remarks
Ci         :  1 Return gradients as COMPLEX matrix zgradm, with indices same as mode/10=0.
Ci         :  2 Same as mode/10=0, but indices are permuted; see pgrmes below
Ci         :  3 Same as mode/10=1, but indices are permuted; see pgrmes below
Ci         :100s digit type of gradient
Ci         :  0 return gradients grad(x,y,z)
Ci         :  1 return gradients grad(+,-,z)
Ci         :  Add 2 to swap components 1 and 2 (x <-> y or + <-> -)
Ci         :Use 101 for SH basis, (+,-,z) gradient and 301 for SH basis, (-,+,z) gradient
Ci   lmax  :l-cutoff for partial waves
Ci   rgrad :rgrad(if1,if2,j,ll,l,nsp) = <g2_lp(if2) grad_r g1_l(if1)> with lp = l +/- 1
Ci         :Here g_l = r*phi_l and phi_l is a partial wave, e.g. (phi(r), phidot(r), phiz(r))
Ci         :rgrad(if1,if2,:,:,:,:) integrals between partial waves g2(if2),g1(if1)
Ci         :rgrad(:,:,1,:,:,:):    (g2(lp) grad_r g1(l)) - (l+1) (g2(lp) 1/r g1(l))
Ci         :rgrad(:,:,2,:,:,:):    (g2(lp) grad_r g1(l)) +    l  (g2(lp) 1/r g1(l))
Ci         :rgrad(:,:,:,1,:,:):    lp = l+1
Ci         :rgrad(:,:,:,2,:,:):    lp = l-1
Ci         :rgrad(:,:,:,:,l,:):    l index to g1(l)
Ci         :rgrad(:,:,:,:,:,isp):  spin index (compound index in noncollinear case)
Ci         :This routine requires only rgrad(:,:,1,1,:,:) and rgrad(:,:,2,2,:,:)
Ci         :If 4's bit of 100's digit is set, treat rgrad as compact dimension:
Ci         :rgrad(:,:,1,:,:):      lp = l+1
Ci         :rgrad(:,:,2,:,:):      lp = l-1
Ci         :rgrad(:,:,:,l,:):      l index to g1(l)
Ci         :To make rgrad, see subroutine rgrme
Ci   nf2   :leading dimension of rgrad, and number of g2 if np2 is zero
Ci   nf1   :second  dimension of rgrad, and number of g1 if np1 is zero
Ci   np1   :np1(0:lmax) = number of partial waves g1 of a particular l
Ci         :If np1(0) = 0, assume np1 = nf1 for all l
Ci   np2   :np2(0:lmax) = number of partial waves g2 of a particular l
Ci         :If np2(0) = 0, assume np2 = nf2 for all l
Ci   nd    :dimensions gradm
Co Outputs
Co   gradm :gradients gradm(k,n2,n1,ilm2,ilm1) = <g2 grad g1> for
Co         : g1 = g1_l1(n1,r) Y_ilm1(rhat) and g2 = g2(n2,l2,r) Y_ilm2(rhat)
Co         : The meaning of Y_lm depends on 1s digit mode mode
Co         : gradm(:,:,:,ilm2,ilm1) gradients for matrix elements (ilm2 | grad | ilm1)
Co         : gradm(:,n2,n1,:,:) gradients for matrix elements (g2(n2) | grad | g1(n1))
Co         : gradm(k,:,:,:,:) component k of gradm, k=1,2,3 for x,y,z, or +,-,z if 2s bit mode is set.
Co         : k=1,2 are swapped if 4s bit mode is set
Co         :
Co   zgradm: Result may be returned in COMPLEX zgradm, depending on 10's digit mode
Co         :
Co         : Note: For either gradm or zgradm, if 10s digit mode > 1,
Co         : indices (k,n2,n1,ilm2,ilm1) are reordered to (ilm2,ilm1,k,n2,n1)
Cl Local variables
Cl   cgcof :coefficients for coupling angular part of gradient operator acting on
Cl         :real harmonics. Subroutine ylmbyr makes cgcof.
Cl   zrcg  :full matrix element = radial factor * angular factor.
Cl         :Matrix is real, but is internally complex to facilitate blas calls.
Cr Remarks
Cr   See www.questaal.org/docs/numerics/spherical_harmonics/
Cr Debugging: uncomment file dumps for rcgx, rcgy, ycgx, ycgy, ycg-, ycg+
Cr The following rotations to cg in spherical harmonics, but xy components, should yield zero:
Cr   mc u rcgx u -i -x -x ycgx --  -px
Cr   mc u rcgy u -i -x -x ycgy --  -px
Cr Check rotations SH(y,x) => SH(-,+)
Cr   set u1 = 'u -rowl 2,4 -coll 2,4 -a u1'
Cr   set defs = 'ycg- -a:nc=1 ycm ycg+ -a:nc=1 ycp ycgy -a:nc=1 ycy ycgx -a:nc=1 ycx'
Cr   mc $u1 $defs ycy ycx -ccat u1 -t -x -coll 1 ycm -- -px
Cr   mc $u1 $defs ycy ycx -ccat u1 -t -x -coll 2 ycp -- -px
Cu Updates
Cu   25 Apr 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,lmax,nf1,nf2,np1(0:lmax),np2(0:lmax),nd
      real(8) :: rgrad(nf2,nf1,2,0:lmax),gradm(3,nf1,nf2,nd,nd)
      complex(8) :: zgradm(3,nf1,nf2,nd,nd)
C ... Local parameters
      integer ilm,j,kavp,kavm,l,m,nlm1,nlm2,n1,n2,mode0,mode2,morder,mg(2,3,(lmax+1)**2)
      double precision rp,rm,cgcof(2,2,3,(lmax+1)**2),saig(3)
      complex(8) :: ut(2,2),zrcg((lmax+1)**2,(lmax+1)**2,3),zycg((lmax+1)**2,(lmax+1)**2,3)
      complex(8), parameter :: sr2 = dcmplx(1d0/dsqrt(2d0),0d0), iot = dcmplx(0d0,1d0)/dsqrt(2d0)

      mode0 = mod(mode,10)
      mode2 = mod(mode/100,10)

      if (np1(0) > 0 .or. np2(0) > 0) call rx('grmes: not ready for np/=0')

      nlm1 = (lmax+1)**2
      nlm2 = (lmax+1)**2        ! For now require a square matrix
      saig = 0
      if (nlm1<=0 .or. nlm2<=0) return

      call ylmbyr(lmax,cgcof,mg)

C --- For each pair partial waves, do ---
      do  n2 = 1, nf2
      do  n1 = 1, nf1

        if (np1(0) > 0) then
          if (n1 > np1(l)) cycle
        endif
        if (np2(0) > 0) then
          if (n2 > np2(l)) cycle
        endif

C       ilm = 0; call dpzero(rcg,size(rcg))
        ilm = 0; call dpzero(zrcg,2*size(zrcg))
        do  l = 0, lmax

          kavp = (l+1)*(l+1)+(l+1)+1
          kavm = (l-1)*(l-1)+(l-1)+1
          rp = rgrad(n1,n2,1,l)  ! (phi_l+1 | gradr | phi_l)
          rm = rgrad(n1,n2,2,l)  ! (phi_l-1 | gradr | phi_l)

          do  m = -l, l
            ilm = ilm+1

C           Note: grad taking Yl to Yl+1 restricted to l<lmax only because zrcg is dimensioned to lmax.
C           Actual grad has in fact components Ylmax+1, which can be used if zrcg is enlarged.
            do  j = 1, 3
              if (abs(mg(1,j,ilm))<=l+1 .and. l<lmax+0) zrcg(kavp+mg(1,j,ilm),ilm,j) = rp*cgcof(1,1,j,ilm)
              if (abs(mg(1,j,ilm))<=l-1) zrcg(kavm+mg(1,j,ilm),ilm,j) = rm*cgcof(1,2,j,ilm)
              if (j<3) then ! there is a second m
                if (abs(mg(2,j,ilm))<=l+1 .and. l<lmax+0) zrcg(kavp+mg(2,j,ilm),ilm,j) = rp*cgcof(2,1,j,ilm)
                if (abs(mg(2,j,ilm))<=l-1) zrcg(kavm+mg(2,j,ilm),ilm,j) = rm*cgcof(2,2,j,ilm)
              endif
            enddo

C            Same as above, explicit m-indexing, without using mg
C            if (l < lmax) then
C              zrcg(kavp+m-1,ilm,1) = rp*cgcof(1,1,1,ilm)
C              zrcg(kavp+m+1,ilm,1) = rp*cgcof(2,1,1,ilm)
C              zrcg(kavp-m-1,ilm,2) = rp*cgcof(1,1,2,ilm)
C              zrcg(kavp-m+1,ilm,2) = rp*cgcof(2,1,2,ilm)
C              zrcg(kavp+m,ilm,3)   = rp*cgcof(1,1,3,ilm)
C            endif
C            if (l > 0) then
C              if (+m-1 >= 1-l .and. +m-1 <= l-1) zrcg(kavm+m-1,ilm,1) = rm*cgcof(1,2,1,ilm)
C              if (+m+1 >= 1-l .and. +m+1 <= l-1) zrcg(kavm+m+1,ilm,1) = rm*cgcof(2,2,1,ilm)
C              if (-m-1 >= 1-l .and. -m-1 <= l-1) zrcg(kavm-m-1,ilm,2) = rm*cgcof(1,2,2,ilm)
C              if (-m+1 >= 1-l .and. -m+1 <= l-1) zrcg(kavm-m+1,ilm,2) = rm*cgcof(2,2,2,ilm)
C              if (kavm+m>0) zrcg(kavm+m,ilm,3) = rm*cgcof(1,2,3,ilm)
C            endif
          enddo
        enddo

C   ... Rotate zrcg from real to spherical harmonics
        if (mode0 > 0) then

          nlm1 = (lmax+1)**2
          nlm2 = (lmax+1)**2    ! For now require a square matrix

          select case (mode0)
            case (1); morder = 100
            case (2); morder = 000
!           case (3); morder = 200 ! Does not work!
            case default; call rx('grmes: bad mode')
          end select

          call s2sph(morder+11,lmax+1,lmax+1,zrcg(1,1,1),nlm1,nlm1,nlm1,nlm1,zycg(1,1,1))
C          call zprm('rcg(x)',2,zrcg(1,1,1),nlm1,nlm1,nlm1)
C          call zprm('ycg(x)',2,zycg(1,1,1),nlm1,nlm1,nlm1)

          call s2sph(morder+11,lmax+1,lmax+1,zrcg(1,1,2),nlm1,nlm1,nlm1,nlm1,zycg(1,1,2))
C          call zprm('rcg(y)',2,zrcg(1,1,2),nlm1,nlm1,nlm1)
C          call zprm('ycg(y)',2,zycg(1,1,2),nlm1,nlm1,nlm1)

          call zcopy(nlm1*nlm2*2,zycg,1,zrcg,1)
        endif

C   ... Rotate  (gradx  grady)  =>  (grad+  grad-)
C                                        ( i   i)
C       (grad-  grad+) =  (grady  gradx) (      )/sr2  or grad(-,+) = grad(y,x) u(odd)^T
C                                        ( 1  -1)
C
C       Here (y,x) order is reversed. Swap rows and columns
C                                        (-1   1)
C       (grad+  grad-) =  (gradx  grady) (      )/sr2  or grad(+,-) = grad(x,y) u~(odd)^T
C                                        ( i   i)
        if (mod(mode2,2) > 0) then

          ut(1,1:2) = [-sr2,  sr2] ! Transpose of u for odd m, swapping rows and columns
          ut(2,1:2) = [ iot,  iot]
          call zcopy(nlm1*nlm2*2,zrcg,1,zycg,1)
          call zgemm('N','N',nlm1*nlm2,2,2,(1d0,0d0),zycg,nlm1*nlm2,ut,2,(0d0,0d0),zrcg,nlm1*nlm1)

C          call zprm('cg(+)',2,zrcg(1,1,1),nlm1,nlm1,nlm1)
C          call zprm('cg(-)',2,zrcg(1,1,2),nlm1,nlm1,nlm1)

        endif

C   ... Swap components 1 and 2
        if (mode2 >= 2) then
          call zcopy(nlm1*nlm2,zrcg(1,1,1),1,zycg(1,1,2),1)
          call zcopy(nlm1*nlm2,zrcg(1,1,2),1,zycg(1,1,1),1)
          call zcopy(nlm1*nlm2*2,zycg(1,1,1),1,zrcg(1,1,1),1)
        endif

C       Debugging
C        call zprm('cg(1)',2,zrcg(1,1,1),nlm1,nlm1,nlm1)
C        call zprm('cg(2)',2,zrcg(1,1,2),nlm1,nlm1,nlm1)
C        call zprm('cg(3)',2,zrcg(1,1,3),nlm1,nlm1,nlm1)

C   ... Copy to gradm in a particular format
        call pgrmes(mode/10,lmax,nlm1,nlm2,n1,n2,nf1,nf2,nd,zrcg,gradm,zgradm,gradm,zgradm)

      enddo
      enddo

      end

      subroutine pgrmes(mode,lmax,nlm1,nlm2,n1,n2,nf1,nf2,nd,zrcg,gradm,zgradm,gradm3,zgradm3)
C- Kernel called by grmes to copy zrcg to gradm in suitable format
Ci  mode   :1s digit controls form of gradm returned by grmes
Ci         :0 zrcg is copied to REAL gradm(k,n2,n1,ilm2,ilm1)
Ci         :  The imaginary part is not copied.
Ci         :1 zrcg is copied to COMPLEX zgradm(k,n2,n1,ilm2,ilm1)
Ci         :2 zrcg is copied to REAL gradm3(ilm2,ilm1,k,n2,n1)
Ci         :3 zrcg is copied to COMPLEX zgradm3(ilm2,ilm1,k,n2,n1)
      implicit none
C ... Passed parameters
      integer mode,lmax,nlm1,nlm2,n1,n2,nf1,nf2,nd
      real(8) :: gradm(3,nf1,nf2,nd,nd),gradm3(nd,nd,3,nf1,nf2)
      complex(8) :: zrcg((lmax+1)**2,(lmax+1)**2,3),zgradm(3,nf1,nf2,nd,nd),zgradm3(nd,nd,3,nf1,nf2)
C ... Local parameters
      integer ilm,jlm,j,mode0
      complex(8), parameter :: sr2 = dcmplx(1d0/dsqrt(2d0),0d0), iot = dcmplx(0d0,1d0)/dsqrt(2d0)

      mode0 = mod(mode,10)

      if (mode0 == 0) then
        forall (jlm=1:nlm1,ilm=1:nlm2,j=1:3) gradm(j,n2,n1,jlm,ilm) = zrcg(jlm,ilm,j)
      elseif (mode0 == 1) then
        forall (jlm=1:nlm1,ilm=1:nlm2,j=1:3) zgradm(j,n2,n1,jlm,ilm) = zrcg(jlm,ilm,j)
      elseif (mode0 == 2) then
        forall (jlm=1:nlm1,ilm=1:nlm2,j=1:3) gradm3(jlm,ilm,j,n2,n1) = zrcg(jlm,ilm,j)
      elseif (mode0 == 3) then
        forall (jlm=1:nlm1,ilm=1:nlm2,j=1:3) zgradm3(jlm,ilm,j,n2,n1) = zrcg(jlm,ilm,j)
      else
        call rx('pgrmes: bad mode')
      endif

      end

C#ifdefC TEST
CC Tests ylmbyr and grmes
CC Run one of:
CC  a.out | grep diff
CC  a.out [--job=cg] [--mode=#]
C      subroutine fmain
C      implicit none
C
C      character*(120) strn,strn2,cxyz*3,dc*1
C      integer, parameter :: nnn=300
C      integer mode,ilmx,nlm,nlmf,lmax,j,ilm,jlm,l,lav,m,mm,lp,pm,kavp,kavm,nf1,nf2
C      integer np,nph,nth
C      double precision pi,srfpi,err(3),xx1 !,xx2
C      double precision p(3,nnn),wp(nnn)
C      real(8) flnum(16,16,3),flanl(16,16,3)
C      integer, allocatable :: mg(:,:,:)
C      real(8), allocatable :: xp(:),yp(:),zp(:),rp(:),r2(:)
C      real(8), allocatable :: yl(:,:),gyl(:,:,:),frl(:)
C      real(8), allocatable :: grp(:,:),ggrp(:,:),coff(:,:,:,:)
C      real(8), allocatable :: rgrad(:,:),gradm(:,:,:,:,:)
C      procedure(logical) :: cmdopt
C      procedure(integer) :: cmdoptsw,a2vec
C
CC     For Clebsch Gordan coefficients
C      integer lmxcg,lmxcy,lnxcg,lnjcg
C      integer, allocatable :: jcg(:),indxcg(:)
C      real(8), allocatable :: cg(:),cy(:),cga(:,:,:)
C      procedure(integer) :: ll
C      integer, parameter :: NULLI=-99999
C      data cxyz / 'xyz'/
C
CC     mode=4  => decompose r^-l grad*(r^l Ylm) into linear combination of Ylm, Ylm = sph. harm. (no r^l)
CC     mode=14 => decompose      grad*(Ylm) into linear combination of Ylm
CC     mode=24 => decompose r^(l+1) grad*(r^(-l-1) Ylm) into linear combination of Ylm
C
C      mode =  4
CC      mode = 14
CC      mode = 24
C
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
C      lmxcg = 8; lmxcy = 12
C      call scg0(0,lmxcg,err,err,err,lnjcg,lnxcg)
C      allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnjcg),cy((lmxcy+1)**2),cga((2*lmax+1)**2,nlm,nlm))
C      call sylmnc(cy,lmxcy)
C      call scg(lmxcg,cg,indxcg,jcg)
C      call mkcga((2*lmax+1)**2,nlm,nlm,indxcg,jcg,cg,cga)
C      nth = -122
C      call fpiint(nth,nph,np,p,wp)
C
CC --- CG coefficients of real harmonics by numerical integration ---
CC ... Make analytic Clebsch Gordan coefficients, calling mkcga
C      if (cmdoptsw('--job',',','cg','') > 0) then
C
C      call info2(1,1,0,' Check Gaunt coffs (array cga) against num int (%i angular points)',np,mode)
C
C      allocate(xp(np),yp(np),zp(np),r2(np),rp(np),grp(np,3),ggrp(np,3))
C      allocate(frl(nlm),yl(np,nlmf),gyl(np,nlm,3))
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
C      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
C
C      jlm = 0; err(1) = 0
C      print *, "  l m     l'm' ilm jlm   L       CG ..."
C      do  lp = 0, 1
C      do  pm = -lp, lp
C        jlm = jlm+1
C        ilm = 0
C        do  l = 0, 2
C        do  m = -l, l
C          ilm = ilm+1
C          gyl(:,1,1) = yl(:,jlm)*yl(:,ilm)
C          call xyl(nlm,np,wp,gyl,r2,yl,frl)
C!         if (ilm == 3) frl(2) = frl(2) + .01d0 ! Debugging
C          xx1 = sum(abs(cga(1:nlm,ilm,jlm)-frl(1:nlm)))
C
C          call awrit6(' Yl(%i%,2i) Yl(%i%,2i)%,3i%,3i: ',strn,len(strn),0,l,m,lp,pm,ilm,jlm)
C          strn2 = '%100pnum-CG = '
C          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C          call pryl(0,strn,ilm,0,frl,nlm,1d0,strn2)
C
C          call awrit4(' Clebsch Gordan :    ',strn,len(strn),0,l,m,lp,pm)
C          if (xx1 > 1d-6) call pryl(0,strn,ilm,0,cga(1,ilm,jlm),nlm,1d0,'')
C          err(1) = err(1) + xx1
C        enddo
C        enddo
C      enddo
C      enddo
C      call info2(1,0,0,' cumulative diff between (analytic, num CG) = %1;3,3g',err,2)
C      deallocate(xp,yp,zp,r2,rp,grp,ggrp,frl,yl,gyl)
C      call rx0('done')
C      endif
C
C      nth = -122
C      call fpiint(nth,nph,np,p,wp)
C      call info2(1,1,0,' %i angular points.  mode=%i'//
C     . '%-1j%?#n==4# (gradients of r^l Y_L, Y_L = real harmonics)##'//
C     . '%-1j%?#n==14# (gradients of real harmonics Y_L)##'//
C     . '%-1j%?#n==24# (gradients of r^-l-1 Y_L, Y_L = real harmonics)##',np,mode)
C
C      allocate(xp(np),yp(np),zp(np),r2(np),rp(np),grp(np,3),ggrp(np,3))
C      allocate(frl(nlm),yl(np,nlmf),gyl(np,nlm,3))
C      allocate(coff(2,2,3,(lmax+1)**2),mg(2,3,(lmax+1)**2))
C
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
C
C      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
C
C      lp=0; if (mod(mode/10,10) == 1) lp=1; if (mod(mode/10,10) == 2) lp=-1
C      call ropylg(lp,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)  ! gradient and laplacian of spherical harmonic polynomials
C      call info0(1,1,0,
C     .  " ... Numerical expansion of grad Y_L as combination of Y_L"//
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
C      call info0(1,1,0,"... analytic expansion (ylmbyr) of grad Y_L as combination of Y_L"//
C     .  "%N%18fdl d+-m  coff ...")
C
C      call ylmbyr(lmax,coff,mg)
CC     f(r) = 1 and r = 1 => (df/dr - l*f/r) = -l  and (df/dr + (l+1)*f/r) = l+1
C      do  ilm = 1, lmax**2
C        l = ll(ilm)
C        if (mode == 4) then
C          coff(:,1,:,ilm) = 0*coff(:,1,:,ilm)
C          coff(:,2,:,ilm) = (2*l+1)*coff(:,2,:,ilm)
C        elseif (mode == 14) then
C          coff(:,1,:,ilm) = -l*coff(:,1,:,ilm)
C          coff(:,2,:,ilm) = (l+1)*coff(:,2,:,ilm)
C        elseif (mode == 24) then
C          coff(:,1,:,ilm) = -(2*l+1)*coff(:,1,:,ilm)
C          coff(:,2,:,ilm) = 0*coff(:,2,:,ilm)
C        endif
C      enddo
C
CC ... ORIGINAL, without mg
CC      err = 0
CC      do  j = 1, 3
CC      do  ilm = 1, lmax**2
CC        l = ll(ilm)
CC        lav = l*l+l+1
CC        kavp = (l+1)*(l+1)+(l+1)+1
CC        kavm = (l-1)*(l-1)+(l-1)+1 ! k- restricted to kmin=l^2+1 to kmax=(l+1)**2
CC        m = ilm-lav
CC        mm = iabs(m)
CC        if (m == -l) call info0(1,1,0,'')
CC        call dpzero(frl,size(frl))
CC        pm = 1
CC        if (j == 1) then
CC          frl(kavp+m-1) = coff(1,1,1,ilm)
CC          frl(kavp+m+1) = coff(2,1,1,ilm)
CCC         print *, m-1,m+1,mg(1,1,ilm),mg(2,1,ilm)
CCC         print *, m-l>=0, 1-m<=l, m-l - ((kavm+m-1) - l**2) !,  '   ', (l+1)**2 - (kavm+m+1),  (- m)+3*l-1
CCC         print *, 3*l-m>=1, 1+m<=l, (l+1)**2 - (kavm+m+1) - (3*l-m-1)
CC          if (m-l - ((kavm+m-1) - l**2) /= 0) stop 'oops'
CC          if ((l+1)**2 - (kavm+m+1) - (3*l-m-1) /= 0) stop 'oops'
CC          if (l>0) then
CC          if (+m-1 >= 1-l .and. +m-1 <= l-1) frl(kavm+m-1) = coff(1,2,1,ilm)
CC          if (+m+1 >= 1-l .and. +m+1 <= l-1) frl(kavm+m+1) = coff(2,2,1,ilm)
CC          endif
CCC         Debugging
CCC          xx1=0; if (kavm+m-1>0) xx1=coff(1,2,1,ilm)
CCC          xx2=0; if (kavm+m+1>0) xx2=coff(2,2,1,ilm)
CCC          call info8(1,0,0," l,m,ilm =%3,3i   l',m1',ilm1' =%3,3i%;11,6D   l',m2',ilm2' =%3,3i%;11,6D",
CCC     .      [l,m,ilm],[l-1,m-1,kavm+m-1],xx1,[l-1,m+1,kavm+m+1],xx2,6,7,8)
CC        endif
CC        if (j == 2) then
CCC         coff(1,1,2,ilm) = term C*Y(l+1,-mk-1)
CCC         coff(2,1,2,ilm) = term C*Y(l+1,-mk+1)
CCC         coff(1,2,2,ilm) = term C*Y(l-1,-mk-1)
CCC         coff(2,2,2,ilm) = term C*Y(l-1,-mk+1)
CC          frl(kavp-m-1) = coff(1,1,2,ilm)
CC          frl(kavp-m+1) = coff(2,1,2,ilm)
CC          if (l>0) then
CC          if (-m-1 >= 1-l .and. -m-1 <= l-1) frl(kavm-m-1) = coff(1,2,2,ilm)
CC          if (-m+1 >= 1-l .and. -m+1 <= l-1) frl(kavm-m+1) = coff(2,2,2,ilm)
CC          endif
CC          pm = -1
CC        endif
CC        if (j == 3) then
CC          frl(kavp+m) = coff(1,1,3,ilm)
CC          if (kavm+m>0) frl(kavm+m) = coff(1,2,3,ilm)
CC        endif
CC        call awrit2(' '//cxyz(j:j)//' grad Yl(%i%,2i): ',strn,len(strn),0,l,m)
CC        pm = 1; if (j == 2) pm = -1
CC        call pryl(0,strn,ilm,pm,frl,nlm,'')
CC        call dcopy(lmax**2,frl,1,flanl(1,ilm,j),1)
CC        err(j) = err(j) + sum(abs(flnum(1:lmax**2,ilm,j)-flanl(1:lmax**2,ilm,j)))
CC      enddo
CC      enddo
C
CC     Vectorized using mg
C      err = 0
C      do  j = 1, 3
C      ilm = 0
C      do  l = 0, lmax-1
C      call info0(1,1,0,'')
C      do  m = -l, l
C        ilm = ilm+1
CC       l = ll(ilm)
C        lav = l*l+l+1
C        kavp = (l+1)*(l+1)+(l+1)+1
C        kavm = (l-1)*(l-1)+(l-1)+1 ! k- restricted to kmin=l^2+1 to kmax=(l+1)**2
CC       m = ilm-lav
C        mm = iabs(m)
C        call dpzero(frl,size(frl))
C
C        if (abs(mg(1,j,ilm))<=l+1) frl(kavp+mg(1,j,ilm)) = coff(1,1,j,ilm)
C        if (abs(mg(1,j,ilm))<=l-1) frl(kavm+mg(1,j,ilm)) = coff(1,2,j,ilm)
C        if (j<3) then
C          if (abs(mg(2,j,ilm))<=l+1) frl(kavp+mg(2,j,ilm)) = coff(2,1,j,ilm)
C          if (abs(mg(2,j,ilm))<=l-1) frl(kavm+mg(2,j,ilm)) = coff(2,2,j,ilm)
C        endif
C
CC       Debugging printout
CC        xx1=0; if (abs(mg(1,j,ilm))<=l-1) xx1=frl(kavm+mg(1,j,ilm))
CC        xx2=0; if (abs(mg(2,j,ilm))<=l-1) xx2=frl(kavm+mg(2,j,ilm)); if (j==3) xx2=0
CC        call info8(1,0,0,' '//cxyz(j:j)//" l,m,ilm =%3,3i   l',m1',ilm1' =%3,3i%;11,6D   l',m2',ilm2' =%3,3i%;11,6D",
CC     .    [l,m,ilm],[l-1,mg(1,j,ilm),kavm+mg(1,j,ilm)],xx1,[l-1,mg(2,j,ilm),kavm+mg(2,j,ilm)],xx2,6,7,8)
C
C        call awrit2(' '//cxyz(j:j)//' grad Yl(%i%,2i): ',strn,len(strn),0,l,m)
C        pm = 1; if (j == 2) pm = -1
C        call pryl(0,strn,ilm,pm,frl,nlm,1d0,'')
C        call dcopy(lmax**2,frl,1,flanl(1,ilm,j),1)
C        err(j) = err(j) + sum(abs(flnum(1:lmax**2,ilm,j)-flanl(1:lmax**2,ilm,j)))
C      enddo
C      enddo
C      enddo
C
C      call info2(1,1,0,
C     .  ' cumulative diff between (analytic, num diff) ME (x,y,z) = %3:1;3,3g',err,2)
CC     call prmx('difference (num-analytical) projection of grad Yl to out.dat',flnum-flanl,lmax**2,lmax**2,lmax**2*3)
C
CC ... Same test but call grmes
C      call info0(1,1,0,
C     .  " ... repeat test, calling grmes for matrix elements.  Compare diff (grme,inline)"//"%N%19fl  m    coff ...")
CC     Same result with call to grmes
C      nf1 = 1; nf2 = 1
C      allocate(gradm(nlm,nlm,3,nf1,nf2))
C
C      allocate(rgrad(2,0:lmax)); call dpzero(rgrad,size(rgrad))
C      do  l = 0, lmax-1
C        if (mode ==  4) rgrad(2,l) =  2*l+1   ! for (l-1, l)
C        if (mode == 14) rgrad(1,l) =  -l      ! for (l+1, l)
C        if (mode == 14) rgrad(2,l) =  l+1     ! for (l-1, l)
C        if (mode == 24) rgrad(1,l) = -2*l-1   ! for (l+1, l)
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
