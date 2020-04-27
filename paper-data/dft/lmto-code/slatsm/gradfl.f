      subroutine gradfl(lmax,nd,nr,np,ir0,ir1,lgg,lx,nn,ri,yl,gyl,fl,gp,ggp)
C- Gradient, Laplacian of function point-wise through sphere from YL expansion
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmax  :function is expanded to L=1:(lmax+1)**2
Ci   nd    :dimensions gyl
Ci   nr    :number of radial mesh points
Ci   np    :number of angular mesh points
Ci   ir0   :gp and gpp are made for radial points between ir0,ir1
Ci   ir1   :gp and gpp are made for radial points between ir0,ir1
Ci   lgg   :if zero, make gradient gp only; gpp not addressed.
Ci   lx    :1s   digit: if 1, fl is scaled by r**2
Ci         :10s  digit: extrapolate 1st point (ir0=1) from others
Ci         :100s digit: rational function interpolation for radial deriv
Ci   nn    :nn: number of points used to differentiate radial f
Ci   ri    :vector of radial mesh points
Ci   yl    :Spherical harmonics for L=0:(lmax+1)^2 at each angular mesh point
Ci         :Generate with a call to ropyln, points on unit radius
Ci   gyl   :Gradient of YL.  Generate with a call to ropylg with lp=1
Ci   fl    :true radial part (no scaling by power of r) of function in Y_lm
Ci         :representation. Function to be differentiated is sum_lm( fl_lm * Y_lm)
Co Outputs
Co   gp    :gradient of fl, on the combined radial and angular mesh,
Co         :x,y,z components
Co   ggp   :Laplacian of fl, on the combined radial and angular mesh
Cl Local variables
Cl   gf    :Work array (used for radial derivative of fl)
Cl   ggf   :Work array (used for 2nd radial derivative of fl)
Cr Remarks
Cr
Cu Updates
Cu   02 Apr 09 Made gf,ggf local; fixed bug for 10s digit lx case
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer np,nr,nd,lx,nn,lmax,ir0,ir1
      double precision fl(nr,(lmax+1)**2),gp(ir0:ir1,np,3),ggp(ir0:ir1,np),
     .  ri(nr),yl(np,*),gyl(np,nd,3)
C ... Local parameters
      integer i0,ilm,ip,ir,j,l,m,l2,lerr,lgg,iprint,jx,nx
      double precision xx,cy1,tol,egf0,gf(nr),ggf(nr)
      logical lrat
      parameter (tol=1d-12)

      if (ir0 < 1) call rx('gradfl: illegal value of ir0')
      cy1 = dsqrt(3/(16*datan(1d0)))
      l2 = mod(lx,100)
      lrat = lx/100 /= 0
      i0 = 1
      nx = nr
      if (l2/10 /= 0) then
        i0 = 2
        nx = nr-1
      endif

C --- Radial derivative : (d fl/dr) yl ---
      call dpzero(gp,    (ir1-ir0+1)*np*3)
      if (lgg /= 0) call dpzero(ggp,   (ir1-ir0+1)*np)
      ilm = 0
      do  20  l = 0, lmax
      do  20  m = -l, l
      ilm = ilm+1
      if (mod(l2,10) == 0) then
        call poldvm(ri(i0),fl(i0,ilm),nx,nn,lrat,tol,lerr,gf(i0))
        if (lerr /= 0) goto 99
      else
        forall (ir = i0:nr) ggf(ir) = fl(ir,ilm)/ri(ir)**2
        call poldvm(ri(i0),ggf(i0),nx,nn,lrat,tol,lerr,gf(i0))
        if (lerr /= 0) goto 99
      endif
C     Extrapolate gf to first point
      if (l2/10 /= 0) then
        jx = 1
        call polint(ri(2),gf(2),nx,nn,ri,0d0,0,jx,gf,egf0)
        lerr = 1
        if (iprint() >= 40 .and. dabs(egf0) > 1d-3*max(dabs(gf(1)),dabs(gf(2)))) then
          call info5(40,0,0,' gradfl (warning): uncertainty in grad'//
     .    ' f(r=0,L=%i):  f=%;3g  est err= %;3g',ilm,gf(1),egf0,0,0)
        endif
      endif
      do  ip = 1, np
      do  ir = ir0, ir1
        gp(ir,ip,1) = gp(ir,ip,1) + gf(ir)*yl(ip,ilm)
      enddo
      enddo

C --- Laplacian: (nabla fl) Yl + fl (nabla Yl) ---
      if (lgg /= 0) then
        call poldvm(ri(i0),gf(i0),nx,nn,lrat,tol,lerr,ggf(i0))
        if (lerr /= 0) goto 99
        if (mod(l2,10) == 0) then
          do  ir = i0, nr
            xx = 1/ri(ir)
            ggf(ir) = ggf(ir) + 2d0*gf(ir)*xx - l*(l+1)*fl(ir,ilm)*xx*xx
          enddo
        else
          do  ir = i0, nr
            xx = 1/ri(ir)
            ggf(ir) = ggf(ir) + 2d0*gf(ir)*xx - l*(l+1)*fl(ir,ilm)*xx**4
          enddo
        endif
        if (i0 == 2) ggf(1)= (ri(3)*ggf(2)-ri(2)*ggf(3))/(ri(3)-ri(2))
        do  ip = 1, np
          do  ir = ir0, ir1
            ggp(ir,ip) = ggp(ir,ip) + ggf(ir)*yl(ip,ilm)
          enddo
        enddo
      endif
   20 continue

C ... Split unit vector r into x,y,z- components
      do  j = 3, 1, -1
        do  ip = 1, np
          xx = yl(ip,j)/cy1
          if (j == 1) xx = yl(ip,4)/cy1
          do  ir = ir0, ir1
            gp(ir,ip,j) = xx*gp(ir,ip,1)
          enddo
        enddo
      enddo

C --- Add fl(r) grad yl + (rhat grad fl) yl.  gf = work array ---
      ilm = 0
      do  l = 0, lmax
      do  m = -l, l
        ilm = ilm+1
C   ... Factor 1/r from grad Yl
        if (mod(l2,10) == 0) then
          do  ir = max(i0,ir0), nr
            gf(ir) = fl(ir,ilm)/ri(ir)
          enddo
        else
          do  ir = max(i0,ir0), nr
            gf(ir) = fl(ir,ilm)/ri(ir)**3
          enddo
        endif
        if (i0 > ir0) gf(1) = (ri(3)*gf(2)-ri(2)*gf(3))/(ri(3)-ri(2))
        do  j = 1, 3
          do  ip = 1, np
            xx = gyl(ip,ilm,j)
            do  ir = ir0, ir1
              gp(ir,ip,j) = gp(ir,ip,j) + gf(ir)*xx
            enddo
          enddo
        enddo
      enddo
      enddo
      return

C --- Error handling ---
   99 print *, 'gradfl: stopping at ilm=',ilm,'  point', lerr
      call rx('gradfl: can''t diff radial function')

      end

C#ifdefC TEST
CC     Optional switch to check radial commutator
CC     Run as:
CC     a.out
CC     a.out --comm~z=0
CC     a.out --commrl~z=#[~enu=#,#,#][~checkse]
CC     a.out --norm
CC     a.out --gradphir | grep cumulative
CC     a.out --rphir | grep cumulative
CC     a.out --phigradphi | grep cumulative
CC     a.out --phirphi | grep cumulative
CC     a.out --rhs     | grep cumulative
CC     a.out --phigradphi | grep cumulative
CC     a.out --phigradphi --grmes | grep cumulative
CC     a.out --phirphi --rmes | grep cumulative
CC     a.out --hrh     | grep cumulative  ... this one not completed
CC     a.out --hsgradhs
C      subroutine fmain
C      implicit none
C      integer ifi,ns,nr,nl
C      real(8),allocatable:: phir(:,:),rofi(:),wi(:,:)
C      double precision a
C
C      logical ltmp
C      integer j,ilm,jlm,nlm,nlmy,nlmg,lmax,l,ir,jr,ip,lerr,lav,kavp,kavm,m,nf1,nf2,iv(10)
C      integer lnjcg,lnxcg,lmxcg,lmxcy,nlmm,klm,k,mk,llm,ml,kav
C      double precision pi,srfpi,xx1,xx2,xx3,xx4,err(3),enu(0:10),csigl(0:10)
C      character*(120) strn,strn2,cxyz*3,dc*1
C      integer, parameter :: nnn=300, npoly=6
C      double precision p(3,nnn),wp(nnn),z
CC     For smooth hankels
C      double precision e1,rsm1,e2,rsm2,rmax,b,hs(0:10),hsp(0:10),dum(1)
C      real(8), parameter :: tol=1d-12
C      integer np,nph,nth
C      procedure(logical) :: cmdopt
C      procedure(integer) :: ll,wordsw,a2vec,fopng,isw
C      procedure(real(8)) :: dot3,dlength
C      integer, parameter :: NULLI=-99999
C      data cxyz / 'xyz'/
C
C      integer pm
CC      real(8) flnum(16,16,3)
C
C      real(8), allocatable :: xp(:),yp(:),zp(:),rp(:),r2(:),res(:,:),resa(:,:)
C      real(8), allocatable :: yl(:,:),gyl(:,:,:),fl(:),phil(:,:),gphir(:,:),philm(:,:),frl(:,:)
C      real(8), allocatable :: wrp(:,:),phip(:,:,:),gphip(:,:,:),rphip(:,:,:,:),r2phip(:,:,:),lphip(:,:),rgrad(:,:),rme(:,:)
C      real(8), allocatable :: phi1(:,:),phip1(:,:,:),rhs1(:,:),rhs2(:,:)
C      real(8), allocatable :: coff(:,:,:,:),flj(:,:),fljn(:,:), gradm(:,:,:,:,:)
C      real(8), allocatable :: ggphir(:),v(:),rix(:),hg(:,:),hrg(:,:),commhg(:,:),vnl(:,:),sigllp(:,:)
C!     real(8), allocatable :: vphi(:,:)
C      integer, allocatable :: mg(:,:,:)
C      real(8), allocatable :: cg(:),cy(:),cga(:,:,:),cgar(:,:,:,:,:)
C      integer, allocatable :: jcg(:),indxcg(:)
C
C      print *, 'Test gradients of g=r*phi(l) and matrix elements (g | grad | g)'
C      print *, 'To run this test, you must supply phir as a disk file'
C
C      ifi = fopng('phir',-1,0) ! contains phi * r
C      nr = 0; nl = 0; a = 0
C      rewind ifi
C      call rdrmsh(ifi,[a],[a],nr,nr,nl,a)
C      if (nr < 2) call rx('file phir has < 2 rows')
C      if (nl < 1) call rx('file phir must have at one l')
C      ns = nr+1
C      allocate(phir(ns,0:nl),rofi(ns),wi(ns,3))
C      rewind ifi
C      call rdrmsh(ifi,rofi,phir,ns,nr,nl,a) ! read phir_l
CC     call prrmsh('test read of phir',rofi,phir,ns,nr,nl)
C      call radwgt(20,rofi(nr),a,nr,wi)      ! weights for int dr f(r)
C      call radwgt(21,rofi(nr),a,nr,wi(1,2)) ! weights for int dr r^2 f(r)
C      forall (ir = 2:nr) wi(ir,3) = wi(ir,2)/rofi(ir) ! weights for int dr r f(r)
C
CC      print *, '!! phir = r/r^(l+1)'
CC      forall (ir = 2:nr, l=0:lmax) phir(ir,l) = 1/rofi(ir)**l
C
C      lmax = nl-1
C      nlm  = (lmax+1)**2
C      nlmg = (lmax+2)**2
C      nlmy = (lmax+3)**2
C      pi = 4*datan(1d0)
C      srfpi = dsqrt(4*pi)
C
CC     Clebsch Gordan coefficients
C      lmxcg = 9; lmxcy = 12
C      call scg0(0,lmxcg,dum,dum,dum,lnjcg,lnxcg)
C      allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg),cy((lmxcy+1)**2))
C      call sylmnc(cy,lmxcy)
C      call scg(lmxcg,cg,indxcg,jcg)
C      allocate(cga((2*lmax+1)**2,nlm,nlm))
C      call mkcga((2*lmax+1)**2,nlm,nlm,indxcg,jcg,cg,cga)
C      nlmm = (2*lmax+2)**2
C      allocate(cgar(nlmm,nlm,nlm,3,2))
C      call mkcgar(1,nlmm,nlm,nlm,indxcg,jcg,cg,cgar)
C
CC     Angular mesh
C      nth = -122
C      call fpiint(nth,nph,np,p,wp)
C      allocate(coff(2,2,3,(lmax+1)**2),mg(2,3,(lmax+1)**2))
C      call ylmbyr(lmax,coff,mg)
C      allocate(xp(np),yp(np),zp(np),r2(np),rp(np))
C      allocate(fl(nlmg),yl(np,nlmy),gyl(np,nlmg,3))
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
C
CC     Spherical harmonics on angular mesh
C      call ropyln(np,xp,yp,zp,lmax+2,np,yl,r2)
C      call ropylg(1,lmax+1,nlmg,np,np,xp,yp,zp,r2,yl,gyl)  ! gradient and laplacian of spherical harmonic polynomials
C
CC --- Debugging check : numerical integration of grad Y_l
CC      call info0(1,1,0,
CC     .  " ... CG coefficients C*Y_l'm' by tabulating grad*Ylm on angular and integrating (grad*Ylm) Yl'm'")
CC      do  j = 1, 3
CC      do  ilm = 1, nlm
CC        l = ll(ilm)
CC        lav = l*l+l+1
CC        m = ilm-lav
CC        if (m == -l) call info0(1,1,0,'')
CC        call xyl(nlmg,np,wp,gyl(1,ilm,j),r2,yl,fl)
CC        if (j == 1) call awrit2(' gradx Yl(%i%,2i): dl,dm,cof = ',strn,len(strn),0,l,m)
CC        if (j == 2) call awrit2(' grady Yl(%i%,2i): dl,d-m,cof= ',strn,len(strn),0,l,m)
CC        if (j == 3) call awrit2(' gradz Yl(%i%,2i): dl,dm,cof = ',strn,len(strn),0,l,m)
CC        pm = 1; if (j == 2) pm = -1
CC        call pryl(0,strn,ilm,pm,fl,nlmg,1d0,'')
CCC       call dcopy(nlm,frl,1,flnum(1,ilm,j),1)
CC      enddo
CC      enddo
CC      stop
C
CC --- Tabulate the (true) phi_l point-wise through sphere for numerical integration
C      allocate(phip(nr,np,nlm),rphip(nr,np,nlm,3),r2phip(nr,np,nlm))
C      do  ilm = 1, nlm
C        do  ip = 1, np
C          l = ll(ilm)
C          xx1 = phir(2,l)/rofi(2) * Yl(ip,ilm)
C          xx2 = phir(3,l)/rofi(3) * Yl(ip,ilm)
C          phip(1,ip,ilm) = (rofi(3)*xx1-rofi(2)*xx2)/(rofi(3)-rofi(2))
C          rphip(1,ip,ilm,:) = 0
C          do  ir = 2, nr
C            phip(ir,ip,ilm) = phir(ir,l)/rofi(ir) * Yl(ip,ilm)
CC           r2phip(ir,ip,ilm) = rofi(ir)**2 * phir(ir,l)/rofi(ir) * Yl(ip,ilm)
C            do  j = 1, 3
C              rphip(ir,ip,ilm,j) = p(j,ip)*phir(ir,l) * Yl(ip,ilm)
C            enddo
C          enddo
C        enddo
C      enddo
CC     Weights for numerical integration (include r^2 in weights)
C      allocate(wrp(nr,np))
C      forall (ir = 1:nr, ip=1:np) wrp(ir,ip) = wi(ir,2) * wp(ip)
C
CC     Confirm the phil = phir/r are normalized (check both wi and wi(1,2))
CC     and make radial derivatives of r*phi
C      if (cmdopt('--norm',NULLI,0,strn)) then
C      call info0(1,1,0," ... Integrals int r (phi_l phi_l) from numerical radial integration")
C      endif
C      allocate(phil(nr,0:lmax),gphir(nr,0:lmax)); call dpzero(phil,size(phil))
C      do  l = 0, lmax
C        forall (ir = 2:nr) phil(ir,l) = phir(ir,l)/rofi(ir) ! phil = true partial wave
C        phil(1,l) = (rofi(3)*phil(2,l)-rofi(2)*phil(3,l))/(rofi(3)-rofi(2))
C        xx1 = dot3(nr,phir(1,l),phir(1,l),wi)        ! integrate (r*phi) * (r*phi) dr
C        xx2 = dot3(nr,phil(1,l),phil(1,l),wi(1,2))   ! integrate phi * phi * (r*r dr)
C        if (cmdopt('--norm',NULLI,0,strn)) then
C        call info5(1,0,0,' l = %i   <phi_l|phi_l>= %;12,6D%;12,6D',l,xx1,xx2,4,5)
C        endif
C        call poldvm(rofi(2),phir(2,l),nr-1,npoly,.false.,tol,lerr,gphir(2,l))
C        gphir(1,l) = (rofi(3)*gphir(2,l)-rofi(2)*gphir(3,l))/(rofi(3)-rofi(2))
C      enddo
CC     call prrmsh('d(phir)/dr',rofi,gphir,ns,nr,nl)
CC     call prrmsh('phip(ilm=8)',rofi,phip(1,1,8),nr,nr,np)
C
CC     Debugging: confirm the phip are orthonormal
C      allocate(res(nlm,nlm))
C      if (cmdopt('--norm',NULLI,0,strn)) then
C      call info0(1,1,0,
C     .  " ... Integrals int d^3r (phi_ilm phi_jlm) from pointwise integration in sphere%N%9fjl jm ...")
C      do  ilm = 1, nlm
C        do  jlm = 1, nlm
C          res(jlm,ilm) = dot3(nr*np,phip(1,1,ilm),phip(1,1,jlm),wrp)
C        enddo
C        call awrit2(' ilm=%,2i',strn,len(strn),0,ilm,2)
C        call pryl(0,strn,ilm,0,res(1,ilm),nlm,1d0,'')
C      enddo
CC     call prmx('res',res,nlm,nlm,nlm)
C      endif
C
CC     Radial integrals <phil | grad phil> and compare against in-line integrals
CC     mcx -vl=0 -vrmax=2.624229 phir -a g g -coll 1,2+l -diff g -coll 1,3+l -tog -ccat -e2 x1 'x2*x4' -int 0 rmax
C      call info0(1,1,0," ... Radial integrals (phi_l'|d/dr|phi_l), compare to rgrme")
C      allocate(rgrad(2,0:lmax),rme(2,0:lmax))
C      call rgrme(0,1,1,1,1,1,lmax,lmax,nr,rofi,wi,phir,phir,rgrad)
C      call rgrme(1,1,1,1,1,1,lmax,lmax,nr,rofi,wi,phir,phir,rme)
C      do  l = 0, lmax
C        xx1 = 0; xx2 = 0; xx3 = 0; xx4 = 0
C        if (l<lmax) then
CC          xx1 = dot3(nr,phil(1,l+1),gphir(1,l),wi(1,3)) !  < r*phi | grad(r*phi) >
CC          xx2 = dot3(nr,phil(1,l+1),phil(1,l),wi(1,3))  !  < r*phi | r*phi/r >
C          xx1 = dot3(nr,phir(1,l+1),gphir(1,l),wi) ! repeat w/ (r*phi) fns only < r*phi | grad(r*phi) >
C          xx2 = dot3(nr,phir(1,l+1),phil(1,l),wi)  ! <r*phi | r*phi/r >
C        endif
C        if (l>0) then
C          xx3 = dot3(nr,phil(1,l-1),gphir(1,l),wi(1,3)) ! int dr r phi_l-1(r) d(r*phi_l(r))/dr
C          xx4 = dot3(nr,phil(1,l-1),phil(1,l),wi(1,3)) ! int dr (r phi_l-1(r)) 1/r (r phi_l(r))
C        endif
C        call info5(1,0,0,' l = %i (l+1|d/dr|l) (l-1|d/dr|l)%;11,6D%;11,6D'//
C     .    ' diff rgrme%;11,6D  %;11,6D',l,xx1,xx3,
C     .    xx1-xx2*(l+1)-rgrad(1,l),xx3+xx4*l-rgrad(2,l))
C      enddo
C
CC --- Compare H g to eps g
CC     mcx -vl=2 -vZ=0 -veps=0.72623440837834485 -f6f20.10 phir -coll 1,2+l -a g g -diff -diff -e2 x1 -x2 v -inc 'x1>0' -e2 x1 'x2-2*Z/x1' g -ccat -e1 'x2*x4' g -e1 'l*(l+1)/x1/x1*x2'  -ccat -ccat -e2 x1 x2+x3+x4 g -e2 x1 'x2*eps' --
CC ... Check commutator.  Try r H g - H r g :
CC     mcx -vl=2 -vZ=0 -veps=0.72623440837834485 -f6f20.10 phir -coll 1,2+l -e2 x1 'x1*x2' -a g g -diff -diff -e2 x1 -x2 \
CC     v -inc 'x1>0' -e2 x1 'x2-2*Z/x1' g -ccat -e1 'x2*x4' g -e1 'l*(l+1)/x1/x1*x2'  -ccat -ccat -e2 x1 x2+x3+x4 -w hrg \
CC     mcx -vl=2 -vZ=0 -veps=0.72623440837834485 -f6f20.10 phir -coll 1,2+l -a g g -diff -diff -e2 x1 -x2 \
CC     v -inc 'x1>0' -e2 x1 'x2-2*Z/x1' g -ccat -e1 'x2*x4' g -e1 'l*(l+1)/x1/x1*x2'  -ccat -ccat -e2 x1 x2+x3+x4 -w hg
CC     mcx hrg hg -e2 x1 'x1*x2' -- -w commhg
CC     mcx hrg out.dat -inc 'x1>0' -coll 1,4 --
CC     Switch --commrl~z=#[~enu=#,#,#][~checkse]
C      j = 8
C      if (cmdopt('--commrl',NULLI,0,strn)) then
C      dc = strn(j+1:j+1)
C
C      stop 'oops'
C      call info0(1,1,0," ... Check commutator (-1/2) [h,r]g against grad_r g")
C      call info0(1,0,0,"     User must supply the local potential in file v")
C
CC     Read v
C      ifi = fopng('v',-1,0)     ! contains potential
C      ilm = nr; z = a           ! hang onto nr, a to confirm that v has same values as phir
C      nr = 0; jlm = 0; a = 0
C      rewind ifi
C      call rdrmsh(ifi,wi,wi,nr,nr,jlm,a)
C      if (nr < 2) call rx('file v has < 2 rows')
C      if (jlm /= 1) call rx('file v must have two columns')
C      if (nr /= ilm) call rx('radial mesh for v and phir are different')
C      ns = nr+1
C      allocate(v(ns),rix(ns))
C      rewind ifi
C      call rdrmsh(ifi,rix,v,ns,nr,jlm,a) ! read v
C      if (dlength(nr,rix-rofi,1) /= 0) call rx('radial mesh for v and phir are different')
C      deallocate(rix)
C
CC     Read Z to add 2*Z/r
C      j = wordsw(strn,dc,'z=','',l) + 2
C      if (j <= 2) call rx('invoke as -commrl~z=#')
C      if (a2vec(strn,len_trim(strn),j,4,', '//dc,3,2,1,ir,z) < 1) call rx('failed to parse '//trim(strn))
CC     Read enu
C      enu = 0
C      j = wordsw(strn,dc,'e=','',l) + 2
C      if (j > 2) then
C        if (a2vec(strn,len_trim(strn),j,4,', '//dc,3,2,1+lmax,iv,enu) < 1+lmax)
C     .    call rx('failed to parse '//trim(strn))
C      endif
C      call info5(1,0,0,"     using Z=%d%?;(n>2);  enu=%n:1d;;",z,j,1+lmax,enu,5)
C
CC     Make H g where g = r*phitrue
C      allocate(ggphir(nr),hg(nr,0:lmax),hrg(nr,0:lmax),commhg(nr,0:lmax))
C      do  l = 0, lmax
C        call poldvm(rofi(2),gphir(2,l),nr-1,npoly,.false.,tol,lerr,ggphir(2))
C        forall (ir = 2:nr) hg(ir,l) = -ggphir(ir) + (v(ir) - 2*z/rofi(ir) + l*(l+1)/rofi(ir)**2)*phir(ir,l)
C        hg(1,l) = (rofi(3)*hg(2,l)-rofi(2)*hg(3,l))/(rofi(3)-rofi(2))
C
CC       Debug: subtract enu*g
C        if (wordsw(strn,dc,'checkse','',j) > 0) then
C          if (enu(0) == 0) call rx('checkse requires enu=...')
C          forall (ir = 2:nr) hg(ir,l) = hg(ir,l) - enu(l)*phir(ir,l)
C        endif
C
C        forall (ir = 2:nr) ggphir(ir) = rofi(ir)*phir(ir,l) ! r*g
C        call poldvm(rofi(2),ggphir(2),nr-1,npoly,.false.,tol,lerr,hrg(2,l)) ! grad (r*g)
C        call poldvm(rofi(2),hrg(2,l),nr-1,npoly,.false.,tol,lerr,ggphir(2)) ! grad grad (r*g)
C        forall (ir = 2:nr) hrg(ir,l) = -ggphir(ir) + (v(ir) - 2*z/rofi(ir) + l*(l+1)/rofi(ir)**2)*rofi(ir)*phir(ir,l)
C        hrg(1,l) = (rofi(3)*hrg(2,l)-rofi(2)*hrg(3,l))/(rofi(3)-rofi(2))
C        forall (ir = 2:nr) commhg(ir,l) = hrg(ir,l) - rofi(ir)*hg(ir,l)
C
C      enddo
C
C      call prrmsh('H g',rofi,hg,nr,nr,nl)
C      call prrmsh('H rg',rofi,hrg,nr,nr,nl)
CC     Write this file to out.dat
CC     mcx -vl=2 -vZ=0 -veps=0.72623440837834485 -f6f20.10 phir  -a g g -diff out.dat -s-1/2 --
C      print *, 'copy this file to commhg'
C      call prrmsh('[H,r]g',rofi,commhg,nr,nr,nl)
C      print *, 'do either of the following:'
C      print *, "mcx -f6f20.10 commhg -s-1/2 out.dat -inc 'x1>0' --"
C      print *, "mcx -f6f20.10 phir -a g g -diff commhg -inc 'x1>0' -s-1/2 --"
CC     mcx -vl=2 -vZ=0 -veps=0.72623440837834485 -f6f20.10 phir  -a g g -diff out.dat -s-1/2 --
C      call prrmsh('grad g',rofi,gphir,nr,nr,nl)
C
C      call rx0('done')
C      endif
C
CC     Switch --commrnl~z=#~csigl=#,#,#]
CC     Nonlocal potential given by sig(r,r') = phi(r) csig(l) sig_ll' csig(l') phi(r')
C      j = 9
C      if (cmdopt('--commrnl',j,0,strn)) then
C      dc = strn(j+1:j+1)
C
C      call info0(1,1,0," ... Check commutator (-1/2) [h,r]g against grad_r g + < g(r') (r'-r) vnl(r,r') g(r)>")
C      call info0(1,0,0,"     User must supply the local potential in file v")
C
CC     Read v
C      ifi = fopng('v',-1,0)     ! contains potential
C      ilm = nr; z = a           ! hang onto nr, a to confirm that v has same values as phir
C      nr = 0; jlm = 0; a = 0
C      rewind ifi
C      call rdrmsh(ifi,wi,wi,nr,nr,jlm,a)
C      if (nr < 2) call rx('file v has < 2 rows')
C      if (jlm /= 1) call rx('file v must have two columns')
C      if (nr /= ilm) call rx('radial mesh for v and phir are different')
C      ns = nr+1
C      allocate(v(ns),rix(ns))
C      rewind ifi
C      call rdrmsh(ifi,rix,v,ns,nr,jlm,a) ! read local v
C      if (dlength(nr,rix-rofi,1) /= 0) call rx('radial mesh for v and phir are different')
C      deallocate(rix)
C
CC     Read Z to add 2*Z/r
C      j = wordsw(strn,dc,'z=','',l) + 2
C      if (j <= 2) call rx('invoke as -commrnl~z=#~csigl=...')
C      if (a2vec(strn,len_trim(strn),j,4,', '//dc,3,2,1,ir,z) < 1) call rx('failed to parse '//trim(strn))
CC     Read csigl : note sig(r,r') = phi(r) csig(l) sig_ll' csig(l') phi(r')
C      csigl = 0
C      j = wordsw(strn,dc,'csigl=','',l) + 6
C      if (j > 6) then
C        if (a2vec(strn,len_trim(strn),j,4,', '//dc,3,2,1+lmax,iv,csigl) < 1+lmax)
C     .    call rx('failed to parse '//trim(strn))
C      else
C        call rx('invoke as -commrnl~z=#~csigl=...')
C      endif
C      call info5(1,0,0,"     using Z=%d%?;(n>6);  csigl=%n:1d;;",z,j,1+lmax,csigl,5)
C
CC     Make sig(r,r')
C      allocate(vnl(nr,nr)); call dpzero(vnl,size(vnl))
C      do  ilm = 0, nl
C        do  jlm = 0, nl
C          do  ir = 1, nr
C            do  jr = 1, nr
C              vnl(ir,jr) = vnl(ir,jr) + phil(ir,ilm)*csigl(ilm)*csigl(jlm)*phil(jr,jlm)
C            enddo
C          enddo
C        enddo
C      enddo
CC     Debugging: confirm hermitian
CC      do  ir = 1, nr
CC        do  jr = 1, nr
CC          if (abs(vnl(ir,jr)-vnl(jr,ir)) > 1d-15) call rx('vnl not hermitian')
CC        enddo
CC      enddo
C
CC     mcx -f10e18.10 phir -e6 x1 x2/x1 x3/x1 x4/x1 x5/x1 x6/x1 > phil
CC     show phil is normalized in weight wi
CC     mch -f10f18.10 phil -coll 1 phil phil wi -xe -xe -show -coll 2:nc -ccat  -rsum
CC     call prrmsh('wi',rofi,wi(1,2),nr,nr,1)
C
CC     mch -f10e18.10 '-array[5,1]' "1 2 1.5 0.1 0.05" -a csigl phil -coll 2:6 csigl -v2dia -x -csum -p -t -x -w vnl -qr out.dat -sub 2,nr,2,nr --  -px
C      call prmx('vnl',vnl,nr,nr,nr)
C
CCC     Debugging ... make  <sig(r,r') phi(r')>
CC      allocate(vphi(nr,0:lmax)); call dpzero(vphi,size(vphi))
CC      do  jlm = 0, lmax
CC        do  ir = 1, nr
CC          do  jr = 1, nr
CC            vphi(ir,jlm) = vphi(ir,jlm) + vnl(ir,jr)*wi(jr,2)*phil(jr,jlm)
CC          enddo
CC        enddo
CC      enddo
CCC     mch -f10e18.10 phil -coll 1 vnl phil wi -xe -coll 2:nc -x -ccat -w vnlphi out.dat -inc 'x1>0' -- -px
CC      call prrmsh('<vnlphi>',rofi,vphi,nr,nr,nl)
C
CC     Equivalently make <phi(r) sig(r,r') phi(r')>
C      allocate(sigllp(0:lmax,0:lmax)); call dpzero(sigllp,size(sigllp))
C      do  ilm = 0, lmax
C        do  jlm = 0, lmax
C          do  ir = 1, nr
C            do  jr = 1, nr
C              sigllp(ilm,jlm) = sigllp(ilm,jlm) + phil(ir,ilm)*wi(ir,2)*vnl(ir,jr)*wi(jr,2)*phil(jr,jlm)
C            enddo
C          enddo
C        enddo
C      enddo
CC     mcx -f10e18.10 phil wi -xe -coll 2:nc -t vnlphi -coll 2:nc  -x -w phivnlphi out.dat -- -px
C      call prmx('<phivnlphi>',sigllp,lmax+1,lmax+1,lmax+1)
C
C      stop "here ... next make <phi (r-r') vnl phi>"
C
CC     Make Hloc g where g = r*phitrue
C      allocate(ggphir(nr),hg(nr,0:lmax),hrg(nr,0:lmax),commhg(nr,0:lmax))
C      do  l = 0, lmax
C        call poldvm(rofi(2),gphir(2,l),nr-1,npoly,.false.,tol,lerr,ggphir(2))
C        forall (ir = 2:nr) hg(ir,l) = -ggphir(ir) + (v(ir) - 2*z/rofi(ir) + l*(l+1)/rofi(ir)**2)*phir(ir,l)
C        hg(1,l) = (rofi(3)*hg(2,l)-rofi(2)*hg(3,l))/(rofi(3)-rofi(2))
C
C
C        forall (ir = 2:nr) ggphir(ir) = rofi(ir)*phir(ir,l) ! r*g
C        call poldvm(rofi(2),ggphir(2),nr-1,npoly,.false.,tol,lerr,hrg(2,l)) ! grad (r*g)
C        call poldvm(rofi(2),hrg(2,l),nr-1,npoly,.false.,tol,lerr,ggphir(2)) ! grad grad (r*g)
C        forall (ir = 2:nr) hrg(ir,l) = -ggphir(ir) + (v(ir) - 2*z/rofi(ir) + l*(l+1)/rofi(ir)**2)*rofi(ir)*phir(ir,l)
C        hrg(1,l) = (rofi(3)*hrg(2,l)-rofi(2)*hrg(3,l))/(rofi(3)-rofi(2))
C        forall (ir = 2:nr) commhg(ir,l) = hrg(ir,l) - rofi(ir)*hg(ir,l)
C
C      enddo
C
C      call prrmsh('Hloc g',rofi,hg,nr,nr,nl)
C      call prrmsh('Hloc rg',rofi,hrg,nr,nr,nl)
CC     Write this file to out.dat
CC     mcx -vl=2 -vZ=0 -veps=0.72623440837834485 -f6f20.10 phir  -a g g -diff out.dat -s-1/2 --
C      print *, 'copy this file to commhg'
C      call prrmsh('[Hloc,r]g',rofi,commhg,nr,nr,nl)
C      print *, 'do either of the following:'
C      print *, "mcx -f6f20.10 commhg -s-1/2 out.dat -inc 'x1>0' --"
C      print *, "mcx -f6f20.10 phir -a g g -diff commhg -inc 'x1>0' -s-1/2 --"
CC     mcx -vl=2 -vZ=0 -veps=0.72623440837834485 -f6f20.10 phir  -a g g -diff out.dat -s-1/2 --
C      call prrmsh('grad g',rofi,gphir,nr,nr,nl)
C
C      call rx0('done')
C      endif
C
CC --- Compare grad [phi_l Y_lm] numerical to analytical angular derivative at point nr
C      if (cmdopt('--gradphir',NULLI,0,strn)) then
C      call info2(1,1,0," ... compare gradfl to analytical grad [phi_l Y_lm] at r(%i)=%d"//
C     .    "%N%7film dl dm ...",nr,rofi(nr))
C      endif
C
C      allocate(gphip(nr,np,3),lphip(nr,np),frl(nr,nlm))  ! lphip is not calculated
C      allocate(philm(nr,nlm))
C      allocate(flj(nlmg,3),fljn(nlmg,3))
C      err = 0
C      do  ilm = 1, nlm
C        call dpzero(philm,size(philm))
C        l = ll(ilm)
C        lav = l*l+l+1; m = ilm-lav
C        forall (ir = 2:nr) philm(ir,ilm) = phir(ir,l)/rofi(ir)  ! True phi_l(r) * Y_lm
C        philm(1,ilm) = (rofi(3)*philm(2,ilm)-rofi(2)*philm(3,ilm))/(rofi(3)-rofi(2))
C        call gradfl(lmax,nlmg,nr,np,1,nr,0,10,npoly,rofi,yl,gyl,philm,gphip,lphip) ! And gradient
C
C        do  j = 1, 3
C          call fp2yl(nr,nlm,1,np,wp,gphip(1,1,j),yl,0d0,frl)
C          fljn(1:nlm,j) = frl(nr,1:nlm)
CC         call prrmsh('gradient of phir from gradfl, Y_l repsn',rofi,frl,ns,nr,nlm)
C        enddo
C
C        kavp = (l+1)*(l+1)+(l+1)+1
C        kavm = (l-1)*(l-1)+(l-1)+1
C        xx1 = (gphir(nr,l)-phil(nr,l)*(l+1))/rofi(nr)
C        xx2 = (gphir(nr,l)+phil(nr,l)*l)/rofi(nr)
C        flj = 0
C
C        do  j = 1, 3
C        if (abs(mg(1,j,ilm))<=l+1) flj(kavp+mg(1,j,ilm),j) = xx1*coff(1,1,j,ilm)
C        if (abs(mg(1,j,ilm))<=l-1) flj(kavm+mg(1,j,ilm),j) = xx2*coff(1,2,j,ilm)
C        if (j<3) then
C          if (abs(mg(2,j,ilm))<=l+1) flj(kavp+mg(2,j,ilm),j) = xx1*coff(2,1,j,ilm)
C          if (abs(mg(2,j,ilm))<=l-1) flj(kavm+mg(2,j,ilm),j) = xx2*coff(2,2,j,ilm)
C        endif
C        enddo
C
CC        if (l < lmax) then
CC          flj(kavp+m-1,1) = xx1*coff(1,1,1,ilm)
CC          flj(kavp+m+1,1) = xx1*coff(2,1,1,ilm)
CC          flj(kavp-m-1,2) = xx1*coff(1,1,2,ilm)
CC          flj(kavp-m+1,2) = xx1*coff(2,1,2,ilm)
CC          flj(kavp+m,3)   = xx1*coff(1,1,3,ilm)
CC        endif
CC        if (l > 0) then
CC          if (kavm+m-1 > 0) flj(kavm+m-1,1) = xx2*coff(1,2,1,ilm)
CC          flj(kavm+m+1,1) = xx2*coff(2,2,1,ilm)
CC          if (kavm-m-1>0) flj(kavm-m-1,2) = xx2*coff(1,2,2,ilm)
CC          flj(kavm-m+1,2) = xx2*coff(2,2,2,ilm)
CC          flj(kavm+m,3) = xx2*coff(1,2,3,ilm)
CC        endif
CC
C
C        if (cmdopt('--gradphir',NULLI,0,strn)) then
C        do  j = 1, 3
C          call awrit2(' grad'//cxyz(j:j)//' %,2i',strn,len(strn),0,ilm,2)
C          xx1 = sum(abs(flj(1:nlm,j)-fljn(1:nlm,j)))
C          strn2 = '%100pnum-an ='
C          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,strn,ilm,pm,fljn(1,j),nlm,1d0,strn2)
C          err(j) = err(j) + xx1
C          if (xx1 >= 1d-6)
C     .    call pryl(0,strn(1:6)//' diff analytic,numeric',ilm,0,flj(1:nlm,j)-fljn(1:nlm,j),nlm,1d0,'')
C        enddo
C        endif
C
C      enddo
C      if (cmdopt('--gradphir',NULLI,0,strn)) then
C        call info2(1,0,0,
C     .  ' cumulative diff between (analytic, num diff) at nr (x,y,z) = %3:1;3,3g',err,2)
C      endif
C
CC --- Compare r [phi_l Y_lm] numerical to analytical expansion
C      if (cmdopt('--rphir',NULLI,0,strn)) then
C      call info2(1,1,0," ... compare r*phi to analytical 1/2 [grad,r^2] phi_l Y_lm] at r(%i)=%d"//
C     .    "%N%7film dl dm ...",nr,rofi(nr))
C      err = 0
C      do  ilm = 1, nlm
C        l = ll(ilm)
C        lav = l*l+l+1; m = ilm-lav
C        do  j = 1, 3
C          call fp2yl(nr,nlm,1,np,wp,rphip(1,1,ilm,j),yl,0d0,frl)
C          fljn(1:nlm,j) = frl(nr,1:nlm)
CC         call prrmsh('gradient of phir from gradfl, Y_l repsn',rofi,frl,ns,nr,nlm)
C        enddo
C
C        kavp = (l+1)*(l+1)+(l+1)+1
C        kavm = (l-1)*(l-1)+(l-1)+1
C        xx1  = rofi(nr)*phil(nr,l)
C        xx2 =  xx1
C        flj = 0
C
C        do  j = 1, 3
C        if (abs(mg(1,j,ilm))<=l+1) flj(kavp+mg(1,j,ilm),j) = xx1*coff(1,1,j,ilm)
C        if (abs(mg(1,j,ilm))<=l-1) flj(kavm+mg(1,j,ilm),j) = xx2*coff(1,2,j,ilm)
C        if (j<3) then
C          if (abs(mg(2,j,ilm))<=l+1) flj(kavp+mg(2,j,ilm),j) = xx1*coff(2,1,j,ilm)
C          if (abs(mg(2,j,ilm))<=l-1) flj(kavm+mg(2,j,ilm),j) = xx2*coff(2,2,j,ilm)
C        endif
C        enddo
C
C        do  j = 1, 3
C          call awrit2(' phi*'//cxyz(j:j)//' %,2i',strn,len(strn),0,ilm,2)
C          xx1 = sum(abs(flj(1:nlm,j)-fljn(1:nlm,j)))
C          strn2 = '%100pnum-an ='
C          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,strn,ilm,pm,fljn(1,j),nlm,1d0,strn2)
C          err(j) = err(j) + xx1
C          if (xx1 >= 1d-6)
C     .    call pryl(0,strn(1:6)//' diff analytic,numeric',ilm,pm,flj(1:nlm,j)-fljn(1:nlm,j),nlm,1d0,'')
C        enddo
C
C      enddo
C      call info2(1,0,0,
C     .  ' cumulative diff between (analytic, num diff) at nr (x,y,z) = %3:1;3,3g',err,2)
C
C      call rx0('done')
C      endif
C
CC --- Integrals phi_l Y_lm grad [phi_l' Y_l'm']
CC     Note: this branch has been bundled into routine grmes and the calculation is repeated below
C      if (cmdopt('--phigradphi',NULLI,0,strn)) then
C      call info0(1,1,0," ... Compare integration int d^3r (phi_l Y_lm grad [phi_l' Y_l'm']"//
C     .  " in sphere to analytic angular integration")
C
CC     Debugging ... replace phi with sm hankels
C      e1 = -2.1d0; rsm1 = 2.0d0; e2 = -1.1d0; rsm2 = 1.6d0
C      do  ir = 1, nr
C        if (ir == nr) call snot
CC       Functions r * hs(1)
C        call hansmd(10,rofi(ir),e1,rsm1,lmax+1,hs,dum,dum,hsp,dum,dum)
C        do  l = 0, lmax
C          phir(ir,l) = rofi(ir)*hs(l) * 100
C        enddo
C      enddo
C      call rgrme(0,1,1,1,1,1,lmax,lmax,nr,rofi,wi,phir,phir,rgrad)
C      do  ilm = 1, nlm
C        do  ip = 1, np
C          l = ll(ilm)
C          xx1 = phir(2,l)/rofi(2) * Yl(ip,ilm)
C          xx2 = phir(3,l)/rofi(3) * Yl(ip,ilm)
C          phip(1,ip,ilm) = (rofi(3)*xx1-rofi(2)*xx2)/(rofi(3)-rofi(2))
C          rphip(1,ip,ilm,:) = 0
C          do  ir = 2, nr
C            phip(ir,ip,ilm) = phir(ir,l)/rofi(ir) * Yl(ip,ilm)
C            do  j = 1, 3
C              rphip(ir,ip,ilm,j) = p(j,ip)*phir(ir,l) * Yl(ip,ilm)
C            enddo
C          enddo
C        enddo
C      enddo
C
C      err = 0
C      do  j = 1, 3
C
C        print *
C        print *, "  l m   ilm  dl dm   phi_K grad"//cxyz(j:j)//" phi_L ..."
C
C        do  ilm = 1, nlm
C          l = ll(ilm)
C          m = ilm-l*l
C
CC         Make grad phi(ilm) on mesh of nr*np points
C          call dpzero(philm,size(philm))
C          forall (ir = 2:nr) philm(ir,ilm) = phir(ir,l)/rofi(ir)
C          philm(1,ilm) = (rofi(3)*philm(2,ilm)-rofi(2)*philm(3,ilm))/(rofi(3)-rofi(2))
CC         call prrmsh('phi, Y_l repsn',rofi,philm,ns,nr,nlm)
C          call gradfl(lmax,nlmg,nr,np,1,nr,0,10,npoly,rofi,yl,gyl,philm,gphip,lphip)
CC         Debugging
CC         call fp2yl(nr,nlm,1,np,wp,gphip,yl,0d0,frl)
CC         call prrmsh('gradient of phir, Y_l repsn',rofi,frl,ns,nr,nlm)
C
C          do  jlm = 1, nlm
CC           res(ilm,jlm) = dot3(nr*np,phip(1,1,jlm),phip(1,1,ilm),wrp)
C            res(jlm,ilm) = dot3(nr*np,phip(1,1,jlm),gphip(1,1,j),wrp)
C          enddo
C
C          call awrit3(' Y(%i%,2i) %,3i: ',strn,len(strn),0,l,m,ilm)
C
CC         pm = 1; if (j == 2) pm = -1
CC         call pryl(0,' num'//strn,ilm,pm,res(1,ilm),nlm,1d0,'')
C
CC         Redo with analytic angular derivative and radial derivative from rgrme
CC         grad1 = (g2(lp) grad_r g1(l)) - (l+1) (g2(lp) 1/r g1(l))
CC         grad2 = (g2(lp) grad_r g1(l)) +    l  (g2(lp) 1/r g1(l))
C          kavp = (l+1)*(l+1)+(l+1)+1
C          kavm = (l-1)*(l-1)+(l-1)+1
C          xx1 = rgrad(1,l)   ! (g2(lp) grad1 g1(l)
C          xx2 = rgrad(2,l)   ! (g2(lm) grad2 g1(l)
C          lav = l*l+l+1; m = ilm-lav
C          flj = 0
C
C          if (abs(mg(1,j,ilm))<=l+1) flj(kavp+mg(1,j,ilm),j) = xx1*coff(1,1,j,ilm)
C          if (abs(mg(1,j,ilm))<=l-1) flj(kavm+mg(1,j,ilm),j) = xx2*coff(1,2,j,ilm)
C          if (j<3) then
C            if (abs(mg(2,j,ilm))<=l+1) flj(kavp+mg(2,j,ilm),j) = xx1*coff(2,1,j,ilm)
C            if (abs(mg(2,j,ilm))<=l-1) flj(kavm+mg(2,j,ilm),j) = xx2*coff(2,2,j,ilm)
C          endif
C
CC          if (l < lmax) then
CC            flj(kavp+m-1,1) = xx1*coff(1,1,1,ilm)
CC            flj(kavp+m+1,1) = xx1*coff(2,1,1,ilm)
CC            flj(kavp-m-1,2) = xx1*coff(1,1,2,ilm)
CC            flj(kavp-m+1,2) = xx1*coff(2,1,2,ilm)
CC            flj(kavp+m,3)   = xx1*coff(1,1,3,ilm)
CC          endif
CC          if (l > 0) then
CC            if (kavm+m-1 > 0)
CC     .      flj(kavm+m-1,1) = xx2*coff(1,2,1,ilm)
CC            flj(kavm+m+1,1) = xx2*coff(2,2,1,ilm)
CC            if (kavm-m-1>0)
CC     .      flj(kavm-m-1,2) = xx2*coff(1,2,2,ilm)
CC            flj(kavm-m+1,2) = xx2*coff(2,2,2,ilm)
CC            flj(kavm+m,3)   = xx2*coff(1,2,3,ilm)
CC          endif
C
C          xx1 = sum(abs(flj(1:nlm,j)-res(1:nlm,ilm)))
C          err(j) = err(j) + xx1
C          strn2 = '%100pnum-an ='
C          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,strn,ilm,pm,flj(1,j),nlm,1d0,strn2)
C          if (xx1 >= 1d-6)
C     .      call pryl(0,'         num',ilm,pm,res(1:nlm,ilm),nlm,1d0,'')
C
C        enddo
C        print *
C      enddo
C      call info2(1,0,0,
C     .  ' cumulative diff between (analytic, num) matrix elements = %3:1;3,3g',err,2)
C      if (.not. cmdopt('--grmes',NULLI,0,strn)) call rx0('done')
C      endif
C
CC --- Integrals phi_l Y_lm r [phi_l' Y_l'm']
CC     Note: this branch has been bundled into routine grmes and the calculation is repeated below
C      if (cmdopt('--phirphi',NULLI,0,strn)) then
C      call info0(1,1,0," ... Compare integration int d^3r (phi_l Y_lm r [phi_l' Y_l'm']"//
C     .  " in sphere to analytic angular integration")
C      err = 0
C
C      call rgrme(1,1,1,1,1,1,lmax,lmax,nr,rofi,wi,phir,phir,rgrad,rme)
C
C      do  j = 1, 3
C
C        print *
C        print *, "  l m   ilm  dl dm   phi_K "//cxyz(j:j)//" phi_L ..."
C
C        do  ilm = 1, nlm
C
C          l = ll(ilm)
C          m = ilm-l*l
C
C          do  jlm = 1, nlm
C            res(jlm,ilm) = dot3(nr*np,phip(1,1,jlm),rphip(1,1,ilm,j),wrp)
C          enddo
C
C          call awrit3(' Y(%i%,2i) %,3i: ',strn,len(strn),0,l,m,ilm)
CC         call pryl(0,' num'//strn,ilm,0,res(1,ilm),nlm,1d0,'')
C
CC         Redo with analytic angular derivative and radial derivative from rgrme
C          kavp = (l+1)*(l+1)+(l+1)+1
C          kavm = (l-1)*(l-1)+(l-1)+1
C          xx1 = rme(1,l)        ! <phi_l+1 | r | phi_l>
C          xx2 = rme(2,l)        ! <phi_l-1 | r | phi_l>
C          lav = l*l+l+1; m = ilm-lav
C          flj = 0
C
C          if (abs(mg(1,j,ilm))<=l+1) flj(kavp+mg(1,j,ilm),j) = xx1*coff(1,1,j,ilm)
C          if (abs(mg(1,j,ilm))<=l-1) flj(kavm+mg(1,j,ilm),j) = xx2*coff(1,2,j,ilm)
C          if (j<3) then
C            if (abs(mg(2,j,ilm))<=l+1) flj(kavp+mg(2,j,ilm),j) = xx1*coff(2,1,j,ilm)
C            if (abs(mg(2,j,ilm))<=l-1) flj(kavm+mg(2,j,ilm),j) = xx2*coff(2,2,j,ilm)
C          endif
C
C          xx1 = sum(abs(flj(1:nlm,j)-res(1:nlm,ilm)))
C          err(j) = err(j) + xx1
C          strn2 = '%100pnum-an ='
C          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,strn,ilm,pm,flj(1,j),nlm,1d0,strn2)
C          if (xx1 >= 1d-6)
C     .      call pryl(0,' num'//strn,ilm,pm,res(1:nlm,ilm),nlm,1d0,'')
C        enddo
C        print *
C      enddo
C      call info2(1,0,0,
C     .  ' cumulative diff between (analytic, num) matrix elements = %3:1;3,3g',err,2)
C      endif
C
CC --- repeat test for grad but call grmes
C      if (cmdopt('--grmes',NULLI,0,strn)) then
C      call info0(1,1,0,' ... repeat test --phigradphi, calling grmes for matrix elements%N%11film  dl dm ...')
C
C      nf1 = 1; nf2 = 1
C      allocate(gradm(3,nf1,nf2,nlm,nlm))
C      call grmes(0,lmax,rgrad,1,1,[0],[0],nlm,gradm,gradm)
C
C      err = 0
C      do  j = 1, 3
C
C        print *
C        print *, "  l m   ilm  dl dm   phi_K "//cxyz(j:j)//" phi_L ..."
C
C        do  ilm = 1, nlm
C
C          l = ll(ilm)
C
C          call dpzero(philm,size(philm))
C          forall (ir = 2:nr) philm(ir,ilm) = phir(ir,l)/rofi(ir)
C          philm(1,ilm) = (rofi(3)*philm(2,ilm)-rofi(2)*philm(3,ilm))/(rofi(3)-rofi(2))
CC         call prrmsh('phi, Y_l repsn',rofi,philm,ns,nr,nlm)
C
C          call gradfl(lmax,nlmg,nr,np,1,nr,0,10,npoly,rofi,yl,gyl,philm,gphip,lphip)
CC         Debugging
CC         call fp2yl(nr,nlm,1,np,wp,gphip,yl,0d0,frl)
CC         call prrmsh('gradient of phir, Y_l repsn',rofi,frl,ns,nr,nlm)
C
C          do  jlm = 1, nlm
CC           res(ilm,jlm) = dot3(nr*np,phip(1,1,jlm),phip(1,1,ilm),wrp)
C            res(jlm,ilm) = dot3(nr*np,phip(1,1,jlm),gphip(1,1,j),wrp)
C          enddo
C        enddo
C
C        do  ilm = 1, nlm
C          l = ll(ilm)
C          m = ilm-l*l
C          call awrit3(' Y(%i%,2i) %,3i: ',strn,len(strn),0,l,m,ilm)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,trim(strn),ilm,pm,gradm(j,1,1,1:nlm,ilm)-res(1:nlm,ilm),nlm,1d0,'')
C          err(j) = err(j) + sum(abs(gradm(j,1,1,1:nlm,ilm)-res(1:nlm,ilm)))
C        enddo
C        print *
C      enddo
C      call info2(1,0,1,
C     .  ' cumulative diff between (analytic, num) matrix elements = %3:1;3,3g',err,2)
C      call rx0('done')
C      endif
C
CC --- Repeat test for integral phi r phi  but call grmes
C      if (cmdopt('--rmes',NULLI,0,strn)) then
C      call info0(1,1,0,' ... repeat test --phirphi, calling grmes for matrix elements')
C
C      nf1 = 1; nf2 = 1
C      allocate(gradm(3,nf1,nf2,nlm,nlm))
C      call grmes(0,lmax,rme,1,1,[0],[0],nlm,gradm,gradm)
C
C      err = 0
C      do  j = 1, 3
C
C        print *
C        print *, "  l m   ilm  dl dm   phi_K "//cxyz(j:j)//" phi_L ..."
C
C        do  ilm = 1, nlm
C!         l = ll(ilm)
C          do  jlm = 1, nlm
C            res(jlm,ilm) = dot3(nr*np,phip(1,1,jlm),rphip(1,1,ilm,j),wrp)
C          enddo
C        enddo
C
C        do  ilm = 1, nlm
C          l = ll(ilm)
C          m = ilm-l*l
C          call awrit3(' Y(%i%,2i) %,3i: ',strn,len(strn),0,l,m,ilm)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,trim(strn),ilm,pm,gradm(j,1,1,1:nlm,ilm)-res(1:nlm,ilm),nlm,1d0,'')
C          err(j) = err(j) + sum(abs(gradm(j,1,1,1:nlm,ilm)-res(1:nlm,ilm)))
C        enddo
C        print *
C      enddo
C      call info2(1,0,1,
C     .  ' cumulative diff between (analytic, num) matrix elements = %3:1;3,3g',err,2)
C      call rx0('done')
C      endif
C
CC --- Compare numerical integration Y_lm r Y_l'm' numerical to gradcg
C      if (cmdopt('--gcklm',NULLI,0,strn)) then
C      endif
C
CC --- Setup for --rhs and --hsrhs and --hsgradhs
C      if (cmdopt('--rhs',NULLI,0,strn) .or. cmdopt('--hsrhs',NULLI,0,strn) .or. cmdopt('--hsgradhs',NULLI,0,strn)) then
C        e1 = -2.1d0; rsm1 = 2.0d0; e2 = -1.1d0; rsm2 = 1.6d0
C
C        call info5(1,1,0,' Setup for tests with smooth Hankels use rsm1=%d e1=%d  rsm2=%d e2=%d',rsm1,e1,rsm2,e2,5)
C
C        if (cmdopt('--hsrhs',NULLI,0,strn) .or. cmdopt('--hsgradhs',NULLI,0,strn)) then
C!         print *, '!! debug'; rmax=rofi(nr)
C          rmax = -dlog(1d-16); a = .01d0; b = .1d0; ns = (1+dlog(1+rmax/b)/a); nr=ns
C!          print *, '!! debug'; e2 = -2.1d0; rsm2 = 2.0d0
C
C          deallocate(rofi); allocate(rofi(nr))
C          call radmsh(rmax,a,nr,rofi)
C          call info2(1,0,1,' Now using nr=%i and rmt=%d',nr,rofi(nr))
C        endif
C
C        allocate(rhs1(ns,0:nl),rhs2(ns,0:nl))
CC       allocate(fl(nlmg),yl(np,nlmy))
CC       call ropyln(np,xp,yp,zp,lmax+2,np,yl,r2)
CC       Remake weights for numerical integration
C        deallocate(phil,phir,wi,wrp,phip,gphip,rphip,philm)
C        allocate(phi1(nr,0:lmax),phil(nr,0:lmax),wi(ns,3),wrp(nr,np),philm(nr,nlm))
C        allocate(phip(nr,np,nlm),rphip(nr,np,nlm,3),gphip(nr,np,3),phip1(nr,np,nlm))
C        call radwgt(20,rofi(nr),a,nr,wi) ! weights for int dr f(r)
C        call radwgt(21,rofi(nr),a,nr,wi(1,2)) ! weights for int dr r^2 f(r)
C        forall (ir = 2:nr) wi(ir,3) = wi(ir,2)/rofi(ir) ! weights for int dr r f(r)
C        forall (ir = 1:nr, ip=1:np) wrp(ir,ip) = wi(ir,2) * wp(ip)
C
C        do  ir = 1, nr
C          if (ir == nr) call snot
CC         Functions r * hs(1)
C          call hansmd(10,rofi(ir),e1,rsm1,lmax+1,hs,dum,dum,hsp,dum,dum)
C          if (ir == nr) print *, '!! scale hs -> 10 hs'; hs = hs*10; hsp = 10*hsp
C          do  l = 0, lmax
C            phi1(ir,l) = hs(l)          ! phi1 = true partial wave, sm hankel hs1
C            rhs1(ir,l) = 2*hsp(l+1)     ! Analytic r*phi1
C          enddo
CC         Functions r * hs(2)
C          call hansmd(10,rofi(ir),e2,rsm2,lmax+1,hs,dum,dum,hsp,dum,dum)
C          if (ir == nr) print *, '!! scale hs -> 10 hs'; hs = hs*10; hsp = 10*hsp
C          do  l = 0, lmax
C            phil(ir,l) = hs(l)          ! phil = true partial wave, sm hankel hs2
C            rhs2(ir,l) = 2*hsp(l+1)     ! Analytic r*phil(l)
C          enddo
C        enddo
C
CC   ... Tabulate phil and rphil on radial + angular mesh
C        do  ilm = 1, nlm
C          do  ip = 1, np
C            l = ll(ilm)
C            phip1(1,ip,ilm) = 0
C            phip(1,ip,ilm) = 0
C            rphip(1,ip,ilm,:) = 0
C            do  ir = 2, nr
C              phip1(ir,ip,ilm) = phi1(ir,l) * Yl(ip,ilm)
C              phip(ir,ip,ilm) = phil(ir,l) * Yl(ip,ilm)
C              do  j = 1, 3
C                rphip(ir,ip,ilm,j) = p(j,ip)*rofi(ir)*phil(ir,l) * Yl(ip,ilm)
C              enddo
C            enddo
C          enddo
C        enddo
C
C      endif
C
CC --- Compare r [hs2 Y_lm] numerical to analytical expansion
C      if (cmdopt('--rhs',NULLI,0,strn)) then
C
C      err = 0
C      do  ilm = 1, nlm
C        l = ll(ilm)
C        lav = l*l+l+1; m = ilm-lav
C        do  j = 1, 3
C          call fp2yl(nr,nlm,1,np,wp,rphip(1,1,ilm,j),yl,0d0,frl)
C          fljn(1:nlm,j) = frl(nr,1:nlm)
C        enddo
C
C        kavp = (l+1)*(l+1)+(l+1)+1
C        kavm = (l-1)*(l-1)+(l-1)+1
C        xx1  = rhs2(nr,l)
C        xx2 =  xx1
C        flj = 0
C
C        do  j = 1, 3
C        if (abs(mg(1,j,ilm))<=l+1) flj(kavp+mg(1,j,ilm),j) = xx1*coff(1,1,j,ilm)
C        if (abs(mg(1,j,ilm))<=l-1) flj(kavm+mg(1,j,ilm),j) = xx2*coff(1,2,j,ilm)
C        if (j<3) then
C          if (abs(mg(2,j,ilm))<=l+1) flj(kavp+mg(2,j,ilm),j) = xx1*coff(2,1,j,ilm)
C          if (abs(mg(2,j,ilm))<=l-1) flj(kavm+mg(2,j,ilm),j) = xx2*coff(2,2,j,ilm)
C        endif
C        enddo
C
C        do  j = 1, 3
C          call awrit2(' hs*'//cxyz(j:j)//' %,2i',strn,len(strn),0,ilm,2)
C          xx1 = sum(abs(flj(1:nlm,j)-fljn(1:nlm,j)))
C          strn2 = '%100pnum-an ='
C          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,strn,ilm,pm,fljn(1,j),nlm,1d0,strn2)
C          err(j) = err(j) + xx1
C          if (xx1 >= 1d-6)
C     .    call pryl(0,strn(1:6)//' diff analytic,numeric',ilm,pm,flj(1:nlm,j)-fljn(1:nlm,j),nlm,1d0,'')
C        enddo
C
C      enddo
C      call info2(1,0,0,
C     .  ' cumulative diff between (analytic, num diff) at nr (x,y,z) = %3:1;3,3g',err,2)
C      call rx0('done')
C      endif
C
CC --- Compare numerical integration [HK ] grad [HL] or [HK ] r [HL] numerical to analytical expression
C      if (cmdopt('--hsrhs',NULLI,0,strn) .or. cmdopt('--hsgradhs',NULLI,0,strn)) then
C
C      ltmp = cmdopt('--hsrhs',NULLI,0,strn)
C
C      call info2(1,1,0," Compare integration int d^3r (H_K %?;(n);r;grad; H_L)"//
C     .  " in sphere to analytic angular integration",isw(ltmp),2)
C
C      do  ir = 2, nr
C        do  l = 0, lmax
C          phi1(ir,l) = rofi(ir)*phi1(ir,l)
C          phil(ir,l) = rofi(ir)*phil(ir,l)
C        enddo
C      enddo
C      call rgrme(0,1,1,1,1,1,lmax,lmax,nr,rofi,wi,phi1,phil,rgrad)
C      if (ltmp) then
C      call rgrme(1,1,1,1,1,1,lmax,lmax,nr,rofi,wi,phi1,phil,rgrad)
C      endif
C      do  ir = 2, nr
C        do  l = 0, lmax
C          phi1(ir,l) = phi1(ir,l)/rofi(ir)
C          phil(ir,l) = phil(ir,l)/rofi(ir)
C        enddo
C      enddo
C
C      call info0(1,1,0," sanity check: int d^3r (H_K H_L)")
C
CC     Compare against analytic overlap integrals of sm hankels  int d^3r (hs1 Y_lm [hs2 Y_l'm']"//
C      allocate(resa(nlm,nlm))
CC      call hhiml([0d0,0d0,0d0],[0d0,0d0,0d0],rsm1,rsm2,e1,e2,nlm,nlm,0,
CC     .  nlm,nlm,cg,indxcg,jcg,cy,resa)
C      call hhiml(0,[0d0,0d0,0d0],[0d0,0d0,0d0],[rsm1,99d0,99d0,99d0,99d0,99d0,99d0],
C     .  [rsm2,99d0,99d0,99d0,99d0,99d0,99d0],
C     .  [e1,99d0,99d0,99d0,99d0,99d0,99d0],
C     .  [e2,99d0,99d0,99d0,99d0,99d0,99d0],nlm,nlm,0,
C     .  nlm,nlm,cg,indxcg,jcg,cy,resa)
CC     Radial integrals done numerically
C      do  l = 0, lmax
C        hs(l) = dot3(nr,phi1(1,l),phil(1,l),wi(1,2))
C      enddo
C      write(*,'(" hhiml",25f10.6)') (resa(ilm,ilm)*100, ilm=1,16)
C      write(*,'(" num  ",25f10.6)') (hs(ll(ilm)), ilm=1,16)
C      write(*,'(" diff ",1p25e10.1)') (resa(ilm,ilm)*100-hs(ll(ilm)), ilm=1,16)
C
CC     Radial integrals <h1 r h2> done numerically
C      call info0(1,1,0," Sanity check compare radial integral < h1 | r | h2 > against analytic ...")
C      do  l = 0, lmax
C        hs(l) = dot3(nr,phi1(1,l),rhs2(1,l),wi(1,2))
C      enddo
CC     Radial integrals done numerically
C      write(*,'(" rgrme",25f10.6)')  (rgrad(1,l), l=0,lmax)
C      write(*,'("  num1",25f10.6)')  (hs(l), l=0,lmax)
C
C
C      stop
C
C
C
C      nf1 = 1; nf2 = 1
C      allocate(gradm(3,nf1,nf2,nlm,nlm))
C      call grmes(0,lmax,rgrad,1,1,[0],[0],nlm,gradm,gradm)
C
C      err = 0
C      do  j = 1, 3
C
C        print *
C        print *, "  l m   ilm  dl dm   H_K "//cxyz(j:j)//" H_L ..."
C
C        do  ilm = 1, nlm
C          l = ll(ilm)
C
C          call dpzero(philm,size(philm))
C!         forall (ir = 2:nr) philm(ir,ilm) = phir(ir,l)/rofi(ir)
C          forall (ir = 2:nr) philm(ir,ilm) = phil(ir,l)
C          philm(1,ilm) = (rofi(3)*philm(2,ilm)-rofi(2)*philm(3,ilm))/(rofi(3)-rofi(2))
CC         call prrmsh('phi, Y_l repsn',rofi,philm,ns,nr,nlm)
C
C          call gradfl(lmax,nlmg,nr,np,1,nr,0,10,npoly,rofi,yl,gyl,philm,gphip,lphip)
CC         Debugging
CC         call fp2yl(nr,nlm,1,np,wp,gphip,yl,0d0,frl)
CC         call prrmsh('gradient of phir, Y_l repsn',rofi,frl,ns,nr,nlm)
C
C          do  jlm = 1, nlm
C            res(jlm,ilm) = dot3(nr*np,phip1(1,1,jlm),gphip(1,1,j),wrp)
C            if (ltmp)
C     .      res(ilm,jlm) = dot3(nr*np,phip1(1,1,jlm),phip(1,1,ilm),wrp)
C          enddo
C        enddo
C
C        do  ilm = 1, nlm
C          l = ll(ilm)
C          m = ilm-l*l
C          call awrit3(' H(%i%,2i) %,3i: ',strn,len(strn),0,l,m,ilm)
C
C
C          xx1 = sum(abs(gradm(j,1,1,1:nlm,ilm)-res(1:nlm,ilm)))
C          err(j) = err(j) + xx1
C          strn2 = '%100pnum-an ='
C          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
C          pm = 1; if (j == 2) pm = -1
C          call pryl(0,strn,ilm,pm,res(1,ilm),nlm,1d0,strn2)
C        enddo
C        print *
C      enddo
C      call info2(1,0,1,
C     .  ' cumulative diff between (analytic, num) matrix elements = %3:1;3,3g',err,2)
C      call rx0('done')
C      endif
C
C
CC          do  jlm = 1, nlm
CC            res(jlm,ilm) = dot3(nr*np,phip(1,1,jlm),rphip(1,1,ilm,j),wrp)
CC          enddo
CC        enddo
CC
CC        do  ilm = 1, nlm
CC          l = ll(ilm)
CC          m = ilm-l*l
CC          call awrit3(' Y(%i%,2i) %,3i: ',strn,len(strn),0,l,m,ilm)
CC          pm = 1; if (j == 2) pm = -1
CC          call pryl(0,trim(strn),ilm,pm,gradm(j,1,1,1:nlm,ilm)-res(1:nlm,ilm),nlm,1d0,'')
CC          err(j) = err(j) + sum(abs(gradm(j,1,1,1:nlm,ilm)-res(1:nlm,ilm)))
CC        enddo
CC        print *
CC      enddo
CC      call info2(1,0,1,
CC     .  ' cumulative diff between (analytic, num) matrix elements = %3:1;3,3g',err,2)
CC      call rx0('done')
C
C
CC --- redundant
CC      call grmes(0,lmax,rme,1,1,[0],[0],nlm,gradm,gradm)
CC
CCC     Compare against analytic overlap integrals of sm hankels  int d^3r (hs1 Y_lm [hs2 Y_l'm']"//
CC      allocate(resa(nlm,nlm))
CC
C
C!!
CCC      call hhiml([0d0,0d0,0d0],[0d0,0d0,0d0],rsm1,rsm2,e1,e2,nlm,nlm,0,
CCC     .  nlm,nlm,cg,indxcg,jcg,cy,resa)
CC      call hhiml(0,[0d0,0d0,0d0],[0d0,0d0,0d0],[rsm1,99d0,99d0,99d0,99d0,99d0,99d0],
CC     .  [rsm2,99d0,99d0,99d0,99d0,99d0,99d0],
CC     .  [e1,99d0,99d0,99d0,99d0,99d0,99d0],
CC     .  [e2,99d0,99d0,99d0,99d0,99d0,99d0],nlm,nlm,0,
CC     .  nlm,nlm,cg,indxcg,jcg,cy,resa)
CC
CC      call info0(1,1,0," Sanity check compare radial integral < H1 | H2 > against analytic hhiml")
CCC     Radial integrals done numerically
CC      do  l = 0, lmax
CC        hs(l) = dot3(nr,phi1(1,l),phil(1,l),wi(1,2))
CC      enddo
CC      write(*,'(" hhiml",25f10.6)') (resa(ilm,ilm), ilm=1,16)
CC      write(*,'(" num  ",25f10.6)') (hs(ll(ilm)), ilm=1,16)
CC      write(*,'(" diff ",1p25e10.1)') (resa(ilm,ilm)-hs(ll(ilm)), ilm=1,16)
C!!
C
CC     Do again on 3D mesh ... don't learn anything new
CC      j  = 1
CC        do  ilm = 1, nlm
CC
CC          l = ll(ilm)
CC
CC          do  jlm = 1, nlm
CC            res(jlm,ilm) = dot3(nr*np,phip1(1,1,jlm),phip(1,1,ilm),wrp)
CC          enddo
CC
CC          call awrit1(' '//' '//'     %,3i',strn,len(strn),0,ilm)
CC
CC          flj(1:nlm,j) = resa(1:nlm,ilm)
CC
CC          xx1 = sum(abs(flj(1:nlm,j)-res(1:nlm,ilm)))
CC          err(j) = err(j) + xx1
CC          strn2 = '%100pnum-an ='
CC          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
CC          pm = 1; if (j == 2) pm = -1
CC          call pryl(0,' ana'//strn,ilm,pm,flj(1,j),nlm,1d0,strn2)
CC          if (xx1 >= 1d-6)
CC     .      call pryl(0,' num'//strn,ilm,pm,res(1:nlm,ilm),nlm,1d0,'')
CC
CC        enddo
C
CC      do  l = 0, lmax
CC        hs(l) = dot3(nr,phi1(1,l),rhs2(1,l),wi(1,2))
CC      enddo
CC
CC      call info0(1,1,0," Integrals of position operator")
CCC     \[{h_1}{Y_K}\vec{r}{h_2}{Y_L}\] = int (h_1 r h_2) int dOmega cgar(nlmm,K,L,p,2))
CCC      print *, cgar(1,:,:,1,1)
CCC      stop
CC
CC
CC      err = 0
CC      do  j = 1, 3
CC        jlm = 1+1 + mod(j+1,3)  ! x for j=1, y for j=2,  z for j=3
CC        print *
CC        print *, "  k mk    l ml klm llm   L       "//cxyz(j:j)//" YK YL ..."
CC
CC        klm = 0
CC        do  k = 0, ll(nlm)
CC        do  mk = -k, k
CC        klm = klm+1
CC
CC          llm = 0
CC          do  l = 0, ll(nlm)
CCC         Radial integral hsk r hsl
CC          hs(l) = dot3(nr,phi1(1,k),rhs2(1,l+1),wi(1,2))
CC          do  ml = -l, l
CC            llm = llm+1
CC
CCC           integral hsk YK r hsl YL
CC            res(llm,klm) = dot3(nr*np,phip1(1,1,klm),rphip(1,1,llm,j),wrp)
CC
CC          enddo
CC          enddo
CC
CC          do  klm = 1, 9
CC          print 345, (res(llm,klm), llm=1,9)
CC  345     format(16f12.6)
CC          enddo
CC          stop
CC
CC          call awrit6(' H(%i%,2i) H(%i%,2i)%,3i%,3i: ',strn,len(strn),0,l,ml,k,mk,klm,llm)
CCC         call pryl(0,' num'//strn,klm,0,res(1,klm),nlm,1d0,'')
CC
CCC         Redo with analytic angular derivative and radial derivative from rgrme
CC          kavp = (k+1)*(k+1)+(k+1)+1
CC          kavm = (k-1)*(k-1)+(k-1)+1
CC          xx1 = rme(1,k)        ! <phi_k+1 | r | phi_k>
CC          xx2 = rme(2,k)        ! <phi_k-1 | r | phi_k>
CC          kav = k*k+k+1; m = klm-kav
CC          flj = 0
CC
CC          if (abs(mg(1,j,klm))<=k+1) flj(kavp+mg(1,j,klm),j) = xx1*coff(1,1,j,klm)
CC          if (abs(mg(1,j,klm))<=k-1) flj(kavm+mg(1,j,klm),j) = xx2*coff(1,2,j,klm)
CC          if (j<3) then
CC            if (abs(mg(2,j,klm))<=k+1) flj(kavp+mg(2,j,klm),j) = xx1*coff(2,1,j,klm)
CC            if (abs(mg(2,j,klm))<=k-1) flj(kavm+mg(2,j,klm),j) = xx2*coff(2,2,j,klm)
CC          endif
CC
CC          xx1 = sum(abs(flj(1:nlm,j)-res(1:nlm,klm)))
CC          err(j) = err(j) + xx1
CC          strn2 = '%100pnum-an ='
CC          call awrit1('%a %;2g',strn2,len(strn2),0,xx1)
CC          pm = 1; if (j == 2) pm = -1
CC          call pryl(0,strn,klm,pm,flj(1,j),nlm,1d0,strn2)
CCC          if (xx1 >= 1d-6)
CCC     .      call pryl(0,' num'//strn,klm,pm,res(1:nlm,klm),nlm,1d0,'')
CC
CC        enddo
CC        enddo
CC        print *
CC      enddo
CC      call info2(1,0,0,
CC     .  ' cumulative diff between (analytic, num) matrix elements = %3:1;3,3g',err,2)
CC      endif
C      call rx0('done')
C
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
CCi   pm    :for printout
CCi         :1   print delta l, delta m = m(ilm) - m(fl)
CCi         :-1  print delta l, delta m = m(ilm) + m(fl)
CCi         :0   print l and m not delta l and delta m
CCi   fl    :function as linear combation of Ylm
CCi   nlm   :print coffs to fl up to nlm
CCi   fac
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
C      logical ldiff
C      procedure(integer) :: ll
C
C      l = ll(mlm)
C      lav = l*l+l+1
C      m = mlm-lav
C      call awrit2(strn,outs,len(outs),0,l,m)
C
C      lmax = ll(nlm)
C      ilm = 0
C      ldiff = .false.
C      do  l2 = 0, lmax
CC        ilm0 = l2**2 + l2 + 1   ! the m=0 element for this l
CC        ilmm = ilm0 - 2*l2      ! the m=0 element for this l-1
CC        ilmp = ilm0 + 2*(l2+1)  ! the m=0 element for this l+1
C        do  m2 = -l2, l2
C          ilm = ilm+1
C          if (abs(fl(ilm)) > tol) then
C            ldiff = .true.
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
C      if (trim(outs) == trim(strn))
C     .  call awrit1('%a*',outs,len(outs),0,tol)
CC    .  call awrit1('%a : no values exceeding tol = %g',outs,len(outs),0,tol)
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
C      subroutine snot
C      end
C#endif


C#ifdefC TEST2
CC Test program to check gradfl against analytic results for grad(Hankel) and grad(Bessel)
C      subroutine fmain
C      implicit none
C      integer nr,nlmx,lmax,ilm,nsize,ll
C      integer nlmf,nnn,np,nph,nth,ilm2
C
C      real(8), allocatable :: xp(:)
C      real(8), allocatable :: yp(:)
C      real(8), allocatable :: zp(:)
C      real(8), allocatable :: r2(:)
C      real(8), allocatable :: yl(:)
C      real(8), allocatable :: gyl(:)
C      real(8), allocatable :: rp(:)
C      real(8), allocatable :: gp(:)
C      real(8), allocatable :: ggp(:)
C      real(8), allocatable :: ggpb(:)
C      real(8), allocatable :: gfl(:)
C      real(8), allocatable :: ggfl(:)
C      real(8), allocatable :: rhol(:)
C      real(8), allocatable :: frl(:)
C
C
C      parameter (nr=250,nlmx=49,nsize=500000,nnn=144)
C      double precision p(3,nnn),wp(nnn),rofi(nr),a,b,rmax,scl,scl2
C      integer nxl(0:7)
C      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/
C      real w(nsize)
C      common /w/ w
C      call wkinit(nsize)
C
C      print *, 'test evaluates grad(Hankel,Bessel) analytically, then numerically'
C
C      ilm = 2
C      scl = 1
C      ilm2 = 5
C      scl2 = 0
C      print *, 'Function is scl1*j*yl(ilm1) + scl2*h*yl(ilm2).'
C      print *, 'Note! as test is now compiled, Function is pure Hankel, with ilm=ilm1'
C      print *, '      scl,ilm2,scl2 are not used.'
C      call info5(0,0,0,' Enter ilm1, scl1, ilm2, scl2 '//
C     .  '(default = %i %d %i %d): ',ilm, scl, ilm2, scl2,0)
C      read(*,*) ilm, scl, ilm2, scl2
C      lmax = ll(max(ilm,ilm2))+2
C
CC ... Angular mesh
C      nth = lmax+1
C      nph = 2*nth
C      nlmf = (lmax+1)**2
C      if (lmax > 6) then
C        nth = 2*lmax+2
C        nph = 0
C      else
C        nth = nxl(lmax)
C        nph = 0
C      endif
C      call fpiint(nth,nph,np,p,wp)
C      print *, nr, 'radial points;', np, ' angular points'
CC     call defrr(oxp,     np)
CC     call defrr(oyp,     np)
CC     call defrr(ozp,     np)
CC     call defrr(or2,     np)
C      allocate(xp(np))
C      allocate(yp(np))
C      allocate(zp(np))
C      allocate(r2(np))
CC ... 3 necessary if two derivatives taken ?
CC     call defrr(oyl,     (lmax+3)**2*np)
C      allocate(yl((lmax+3)**2*np))
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
CC     call defrr(ogyl,     nlmf*np*3)
C      allocate(gyl(nlmf*np*3))
C
CC ... Radial mesh
C      rmax = .1d0
C      a = 1d-6
C      a = .001
C      b = rmax/(dexp(a*(nr-1))-1d0)
C      call radmsh(rmax,a,nr,rofi)
CC ... Other setup
CC     call defrr(orp,     nr*np)
CC     call defrr(ogp,     nr*np*3)
CC     call defrr(oggp,    nr*np)
CC     call defrr(oggpb,   nr*np*3*3)
CC     call defrr(ogfl,    nr)
CC     call defrr(oggfl,   nr)
CC     call defrr(orhol,   nr*nlmf)
CC     call defrr(ofrl,    nr*nlmf)
C      allocate(rp(nr*np))
C      allocate(gp(nr*np*3))
C      allocate(ggp(nr*np))
C      allocate(ggpb(nr*np*3*3))
C      allocate(gfl(nr))
C      allocate(ggfl(nr))
C      allocate(rhol(nr*nlmf))
C      allocate(frl(nr*nlmf))
C
C      call testg(rofi,nr,np,ilm,scl,ilm2,scl2,nlmf,rp,
C     .  gp,ggp,ggpb,frl,yl,gyl,wp,rhol,
C     .  gfl,ggfl,xp,yp,zp,r2)
C
C      end
C      subroutine radmsh(r,a,nr,rofi)
C      implicit real*8 (a-h,p-z), integer (o)
C      dimension rofi(nr)
C      b=r/(dexp(a*nr-a)-1.d0)
C      do 1 ir=1,nr
C    1 rofi(ir)=b*(dexp(a*ir-a)-1d0)
C      end
C      subroutine prmr(strn,nr,nrx,rofi,pow,f,nl)
C      implicit none
C      integer nr,nrx,nl,ir,j,fopna,ifi
C      double precision rofi(nrx),f(nrx,nl),pow
C      character*(10) strn*(*)
C      ifi = fopna('out',19,0)
C      write(ifi,"('% rows',i4,' cols',i4,' scaled by r^',f8.4)") nr-1,nl+1,pow
C      do  10  ir = 2, nr
C        write(ifi,333) rofi(ir),
C     .    (f(ir,j)*(rofi(ir)+1d-12)**pow, j=1, nl)
C  333   format(f12.5/(1p10e14.6))
CC  333   format(f12.5,(7g18.10:/12x))
CC  333   format(f12.5,(17f14.6:/12x))
C   10 continue
C      call fclose(ifi)
CC     call snot
C      print *, trim(strn) // '( file out.dat)'
C      pause
C      end
C      subroutine testg(rofi,nr,np,ilm1,scl1,ilm2,scl2,nlm,rp,gp,ggp,
C     .  ggpb,frl,yl,gyl,wp,rl,gfl,ggfl,xp,yp,zp,r2)
CCi ilm1 : function is f(r) * YL(ilm1)
C      implicit none
C      integer nr,np,ilm1,ilm2,nlm
C      double precision rofi(1),yl(np,nlm),gyl(np,nlm,3),
C     .  wp(np),rl(nr,nlm),frl(nr,nlm),gfl(nr),ggfl(nr),
C     .  rp(nr,np),gp(nr,np,3),ggp(nr,np),ggpb(nr,np,3,3),x2,
C     .  xp(1),yp(1),zp(1),r2(1),xx,scl1,scl2,phi(0:20),psi(0:20),e,pow
C      integer ip,ir,ll,lmax,nlmx,lx1,lx2,lr2,ir0,ir1,ilm
C      logical lrat
C
C      ir0 = 247
C      ir1 = 249
C      lrat = .false.
C      lr2 = 10
C      lmax = ll(nlm)
C      lx1 = ll(ilm1)
C      lx2 = ll(ilm2)
C      nlmx = min(nlm,16)
C      if (ll(nlm) < lx1+2) call rx('testg: need bigger nlm')
C      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
C      call ropylg(1,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
C
CC ... Show by brute force for rl = 1, laplacian is -l(l+1)/r**2
CC      print *,
CC     .'Show by brute force for rl = 1, laplacian is -l(l+1)/r**2'
CC      call dpzero(rl, nr*nlm)
CC      call dcopy(nr,scl1,0,rl(1,ilm1),1)
CC      do  ip = 1, np
CC      do  ilm = 1, nlm
CC        do  ir = 1, nr
CC          rp(ir,ip) = rp(ir,ip) + rl(ir,ilm)*yl(ip,ilm)
CC        enddo
CC      enddo
CC      enddo
CC      gp = 0
CC      call blap(rofi,nr,np,nlm,rp,gp,ggp,ggpb,1,nr,
CC     .  lrat,lr2,frl,yl,gyl,wp,gfl,ggfl,0d0,lmax,nlmx)
C
CC...  Gradient, laplacian of Hankel or Bessel tabulated on a mesh
C      e = -.7d0
C      print *, 'nabla of hankel or bessel, brute force, e=',e
C      do  110  ir = 2, nr
C      call bessl(e*rofi(ir)**2,max(lx1,lx2),phi,psi)
CC ... here for Bessel
CC      lrat = .false.
CC      pow = -lx1
CC      xx = rofi(ir)**lx1
CC      x2 = rofi(ir)**lx2
CC      do  110  ip = 1, np
CC  110 rp(ir,ip) = scl1*phi(lx1)*xx*yl(ip,ilm1) +
CC     .            scl2*phi(lx2)*x2*yl(ip,ilm2)
C
CC ... here for Hankel (note for near origin, need rational f interp.
C      lrat = .true.
C      pow = lx1+1
C      xx = rofi(ir)**(-lx1-1)
C      do  110  ip = 1, np
C  110 rp(ir,ip) = psi(lx1)*xx*yl(ip,ilm1)
C
C      call makghl(rofi,nr,np,ilm1,nlm,gp,frl,yl,xp,yp,zp,wp,e,pow)
C
C      if (mod(lr2,10) /= 0) then
C        do  60  ip = 1, np
C        do  60  ir = 1, nr
C   60   rp(ir,ip) = rp(ir,ip)*rofi(ir)**2
C      endif
C      print *, 'numerical derivatives of Hankels ...'
C      call blap(rofi,nr,np,nlm,rp,gp,ggp,ggpb,
C     .  ir0,ir1,lrat,lr2,frl,yl,gyl,wp,gfl,ggfl,pow,lmax,nlmx)
C
C
C      end
C      subroutine blap(ri,nr,np,nlm,rp,gp,ggp,ggpb,
C     .  ir0,ir1,lrat,lr2,fl,yl,gyl,wp,gf,ggf,pow,lmax,nlmx)
CC- Gradient and Laplacian by brute force of function tabulated on a mesh
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ri
CCi   nr    :number of radial mesh points
CCi   np
CCi   nlm
CCi   rp    :function tabulated on a combined (radial, angular mesh)
CCi   ir0   :first point to differentiate for grad(grad) only
CCi   ir1   :last point to differentiate for grad(grad) only
CCi   lrat
CCi   lr2
CCi   fl
CCi   yl
CCi   gyl
CCi   wp
CCi   gf
CCi   ggf
CCi   pow
CCi   lmax  :maximum l for a given site
CCi   nlmx
CCo Outputs
CCo   gp    :gradient of rp, calculated by gradfl
CCo   ggp   :laplacian of rp, calculated by gradfl
CCo   ggpb
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   20 Aug 15  Patched for 2015 version of codes
CC ----------------------------------------------------------------------
C      implicit none
C      integer nr,np,nlm,lr2,ir0,ir1
C      logical lrat
C      double precision ri(nr),yl(np,nlm),gyl(np,nlm,3),
C     .  wp(np),fl(nr,nlm),gf(nr),ggf(nr),
C     .  rp(nr,np),gp(nr,np,3),ggp(nr,np),ggpb(ir0:ir1,np,3,3),xx,pow
C      integer ip,ir,lmax,j,nlmx,lx,itwo(2)
C
C      lx = 0
C      if (lrat) lx = 100
C      call fp2yl(nr,nlm,1,np,wp,rp,yl,0d0,fl)
C      call csmaln(fl,nr*nlm,1d-8,-1,itwo,itwo)
C      call prmr('Yl expansion of function rp ...',nr,nr,ri,pow,fl,min(nlmx,9))
CC ... Gradient and Laplacian
C      call gradfl(lmax,nlm,nr,np,1,nr,1,lx+lr2,8,ri,yl,gyl,fl,gp,ggp)
C      call fp2yl(nr,nlm,1,np,wp,ggp,yl,0d0,fl)
C      call csmaln(fl,nr*nlm,1d-8,-1,itwo,itwo)
C      call prmr('Yl expansion of Laplacian of rp from gradfl ...',nr,nr,ri,pow,fl,min(nlmx,9))
CC ... grad (gradient) ... points ir0:ir1 only
C      do  j = 1, 3
C        print *, 'Make Yl expansion of gradient, component', j
C        call fp2yl(nr,nlm,1,np,wp,gp(1,1,j),yl,0d0,fl)
C        call csmaln(fl,nr*nlm,1d-8,-1,itwo,itwo)
C        call prmr('Yl expansion of gradient of rp from gradfl',nr,nr,ri,pow+1,fl,nlmx)
CC        Uncomment the following lines if to calculate laplacian from grad(grad)
CC        call gradfl(lmax,nlm,nr,np,ir0,ir1,0,lx+10,8,ri,yl,gyl,fl,
CC     .      ggpb(ir0,1,1,j),ggp)
CC        call fp2yl(ir1-ir0+1,nlm,1,np,wp,ggpb(ir0,1,j,j),yl,0d0,fl)
CC        call prmr('ggpb',ir1-ir0+1,ir1-ir0+1,ri(ir0),pow+2,fl,nlmx)
C      enddo
C      return
C
CC     These lines make Laplacian by grad(grad)
C      do  16  ip = 1, np
C      do  16  ir = ir0, ir1
C      xx = ggpb(ir,ip,1,1) + ggpb(ir,ip,2,2) + ggpb(ir,ip,3,3)
C   16 ggpb(ir,ip,1,1) = xx
C      call fp2yl(ir1-ir0+1,nlm,1,np,wp,ggpb,yl,0d0,fl)
C      call prmr('Laplacian by grad(grad)',ir1-ir0+1,ir1-ir0+1,ri(ir0),
C     .  pow,fl,min(nlm,16))
C
C      end
C      subroutine makghl(ri,nr,np,ilm1,nlm,gp,fl,yl,xp,yp,zp,wp,e,pow)
CC- Gradient of solid Hankel functions
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ri    :radial mesh
CCi   nr    :number of radial mesh points
CCi   np    :number of points on angular mesh
CCi   ilm1  :Return hankel and gradient for ilm1 component
CCi   nlm
CCi   yl
CCi   xp,yp,zp : x-, y-, and z- coordinates of points on angular mesh
CCi   wp    :angular mesh weights
CCi   e     :Energy of solid Hankel
CCi   pow   :not used
CCo Outputs
CCo   gp    :gradient of solid hankel, on combined (radial, angular) mesh
CCo   fl    :representation of gp in spherical harmonics
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   20 Aug 15
CC ----------------------------------------------------------------------
C      implicit none
C      integer nr,np,nlm,ilm1
C      double precision ri(nr),gp(nr,np,3),
C     .  fl(nr,nlm),yl(np,nlm),xp(np),yp(np),zp(np),wp(np),e,pow
C      integer ndim,ir,ip,j,lmax,ll,nlmx,itwo(2)
C      parameter (ndim=200)
C      double precision dr(3),cy(16**2),hl(ndim),ghl(ndim,3),
C     .   hd(ndim),ghd(ndim,3)
C      common /static/ cy
C
C      nlmx = min(nlm,9)
C      call sylmnc(cy,15)
C      lmax = ll(ilm1)
C      do  ir = 2, nr
C      do  ip = 1, np
C        dr(1) = xp(ip)*ri(ir)
C        dr(2) = yp(ip)*ri(ir)
C        dr(3) = zp(ip)*ri(ir)
C        call solhpg(e,dr,lmax,ndim,hl,ghl,hd,ghd,cy)
C        do  j = 1, 3
C          gp(ir,ip,j) = ghl(ilm1,j)
C        enddo
C      enddo
C      enddo
C      print *, 'Yl expansion of analytic Hl, 3 components written in succession (file out.dat)'
C      do  j = 1, 3
C        print *, 'component',j
C        call fp2yl(nr,nlm,1,np,wp,gp(1,1,j),yl,0d0,fl)
C        call csmaln(fl,nr*nlm,1d-8,-1,itwo,itwo)
C        call prmr('Yl expansion of analytic grad hl',nr,nr,ri,pow+1,fl,nlmx)
C      enddo
C
C      end
C      subroutine solhpg(e,dr,lmax,ndim,hl,ghl,hd,ghd,cy)
CC- Solid Hankel functions with energy derivatives and gradients
C      implicit real*8 (a-h,p-z), integer (o)
C      dimension cy(1),dr(3),hl(ndim),ghl(ndim,3),phi(0:30),psi(0:30),
C     .   hd(ndim),ghd(ndim,3)
C      nlm=(lmax+1)**2
C      if((lmax+2)**2>ndim) call rx('solhgp: ndim too small')
C
CC --- Make solid Hankel functions HL ---
C      call sylm(dr,hl,lmax+1,r2)
C      call bessl(e*r2,lmax+2,phi,psi)
C      ilm=0
C      fac=dsqrt(r2)
C      do 10 l=0,lmax+1
C        fac=fac/r2
C        psidot=((l+l+1)*psi(l)-psi(l+1))/(e+e)
C        do 10 m=-l, l
C        ilm=ilm+1
C        hd(ilm)=fac*psidot*cy(ilm)*hl(ilm)
C  10    hl(ilm)=fac*psi(l)*cy(ilm)*hl(ilm)
C
CC ------ make gradients ----------
C      do 20 m=1,3
C      do 20 ilm=1,nlm
C      ghd(ilm,m)=0d0
C  20  ghl(ilm,m)=0d0
C
C      nlm1=lmax*lmax
C      do 22 ilm=1,nlm
C      call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
C      ghl(ilm,1)=ghl(ilm,1)-cx1*hl(kx1)-cx2*hl(kx2)
C      ghl(ilm,2)=ghl(ilm,2)-cy1*hl(ky1)-cy2*hl(ky2)
C      ghl(ilm,3)=ghl(ilm,3)-cz*hl(kz)
C      ghd(ilm,1)=ghd(ilm,1)-cx1*hd(kx1)-cx2*hd(kx2)
C      ghd(ilm,2)=ghd(ilm,2)-cy1*hd(ky1)-cy2*hd(ky2)
C      ghd(ilm,3)=ghd(ilm,3)-cz*hd(kz)
C      if(ilm<=nlm1) then
C        xx=e*hl(ilm)
C        ghl(kx1,1)=ghl(kx1,1)+cx1*xx
C        ghl(kx2,1)=ghl(kx2,1)+cx2*xx
C        ghl(ky1,2)=ghl(ky1,2)+cy1*xx
C        ghl(ky2,2)=ghl(ky2,2)+cy2*xx
C        ghl(kz,3)=ghl(kz,3)+cz*xx
C        xx=hl(ilm)+e*hd(ilm)
C        ghd(kx1,1)=ghd(kx1,1)+cx1*xx
C        ghd(kx2,1)=ghd(kx2,1)+cx2*xx
C        ghd(ky1,2)=ghd(ky1,2)+cy1*xx
C        ghd(ky2,2)=ghd(ky2,2)+cy2*xx
C        ghd(kz,3)=ghd(kz,3)+cz*xx
C        endif
C  22  continue
C      end
C
C      subroutine snot
C      return
C      end
C
C
C#endif

C#ifdefC TEST3
CC Decomposes a function of angle into SH
C      subroutine fmain
C      implicit none
C      integer ifi,ns,nr,nl
C      real(8),allocatable:: phir(:,:),rofi(:),wi(:,:)
C      double precision a
C      procedure(integer) :: fopng
C
C      integer j,ilm,jlm,nlm,nlmy,nlmg,lmax,l,ir,ip,lerr,lav,kavp,kavm,m,nf1,nf2,i
C      double precision pi,srfpi,xx1,xx2,xx3,xx4,err(3)
C      character*(120) strn,cxyz*3
C      integer, parameter :: nnn=300, npoly=6
C      double precision p(3,nnn),wp(nnn)
C      real(8), parameter :: tol=1d-12
C      integer np,nph,nth
C      procedure(integer) :: ll
C      procedure(real(8)) :: dot3
C      procedure(logical) :: a2bin
C      data cxyz / 'xyz'/
C
C
CC     integer j,ilm,lav,m,pm
CC     real(8) flnum(16,16,3),flanl(16,16,3)
C
C      real(8), allocatable :: xp(:),yp(:),zp(:),rp(:),r2(:),res(:,:)
C      real(8), allocatable :: yl(:,:),gyl(:,:,:),fl(:),phil(:,:),gphir(:,:),philm(:,:),frl(:,:)
C      real(8), allocatable :: wrp(:,:),phip(:)
C      real(8), allocatable :: coff(:,:,:,:),flj(:,:),fljn(:,:), gradm(:,:,:,:,:)
C
C      print *, 'Decompose function of angle into real harmonics'
C      write(*,3,advance='no')
C    3 format(' enter function of x,y,z : ')
C      read(*,'(a120)') strn
C
C      nl = 6
C      lmax = nl-1
C      nlm  = (lmax+1)**2
C      nlmg = (lmax+2)**2
C      nlmy = (lmax+3)**2
C      pi = 4*datan(1d0)
C      srfpi = dsqrt(4*pi)
C      nth = -122
C      call fpiint(nth,nph,np,p,wp)
C
C      allocate(xp(np),yp(np),zp(np),r2(np),yl(np,nlmy),fl(nlmg))
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
C
C      call ropyln(np,xp,yp,zp,lmax+2,np,yl,r2)
C
CC ... Tabulate the function point-wise through sphere for numerical integration
C      allocate(phip(np))
C      do  ip = 1, np
C        call lodsyv('x',1,xp(ip),i)
C        call lodsyv('y',1,yp(ip),i)
C        call lodsyv('z',1,zp(ip),i)
C        i = 0
C        if (.not. a2bin(strn,phip(ip),4,0,' ',i,-1)) call rx('failed to parse expr')
C      enddo
C      call xyl(nlmg,np,wp,phip,r2,yl,fl)
C      call pryl(0,strn,1,0,fl,nlmg,1d0,'')
C
C
C      end
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
C
C      subroutine pryl(mode,strn,mlm,pm,fl,nlm,fac)
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
CCi   fl    :function as linear combation of Ylm
CCi   nlm   :print coffs to fl up to nlm
CCi   fac
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
C      character*(*) strn
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
C      call info0(1,0,0,trim(outs))
C
C      end
