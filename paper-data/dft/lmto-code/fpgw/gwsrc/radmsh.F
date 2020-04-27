      subroutine radmsh(rmax,a,nr,rofi)
C- Makes shifted logarithmic mesh
      implicit none
      integer nr,ir
      double precision rmax,a,rofi(nr),b,rpb,ea
C     if (mod(nr,2) == 0) call rx('radmsh: nr is even')
      b = rmax / (dexp(a*nr-a)-1d0)

      ea = dexp(a)
      rpb = b
      do  ir = 1, nr
        rofi(ir) = rpb - b
        rpb = rpb*ea
      enddo
      end
      subroutine radwgt(opt,rmax,a,nr,wt)
C- Weights for numerical integration on shifted log mesh
C ----------------------------------------------------------------------
Ci   opt   :1 + 10's digit same as radmwt
Ci         :0 for uniform weight, i.e. int dr f(r)
Ci         :1 for r^2 weight,     i.e. int dr r^2 f(r)
Ci         :10s digit
Ci         :0 for 3-point quadrature (Simpson's rule)
Ci         :1 for 5-point quadrature
Ci         :2 for 7-point quadrature
Ci         :See radmwt for special treatment of endpoints
Ci   rmax  :augmentation radius, in a.u.
Ci   a     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Co Outputs
Co   wt    :weights for numerical integration on shifted logarithmic radial mesh
Cl Local variables
Cl         :
Cr Remarks
Cr   Routine has identical function to calling radmwt with 100's digit opt set
Cr
Cu Updates
Cu   04 Apr 12 Extended for 7-point quadrature, and for nr even
C ----------------------------------------------------------------------
      implicit none
      integer opt,nr
      double precision rmax,a,wt(nr)
      double precision xx

C     if (mod(nr,2) == 0) call rx('radwgt: nr is even')
      call radmwt(100+mod(opt,100),rmax,a,nr,xx,wt)

      end

      subroutine mshnrl(opt,rmax,a,nr,nrmx,nrl)
C- Estimate mesh size for transposing to shifted to standard log mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :determines criteria for nrl
Ci         :0  Return nrl = nr
Ci         :1  Find what nrl makes rofil(1) = rofi(2)
Ci   rmax  :augmentation radius, in a.u.
Ci   a     :radial mesh points are given by
Ci         :shifted log mesh:  rofi(i) = b [e^(a(i-1)) -1]
Ci         :            where  b = rmax / (dexp(a*nr-a)-1d0)
Ci         :standard log mesh: rofil(i) = [rmax/e^(a(nr-1))] e^(a(i-1))
Ci   nr    :number of radial mesh points
Ci   nrmx  :upper limit to number of mesh points
Co Outputs
Co   nrl   :Estimate number of mesh points for standard mesh
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   25 Apr 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nr,nrmx,nrl
      double precision rmax,a
C ... Local parameters
      double precision b,rofi2

      if (opt == 0) then
        nrl = nr
        return
      endif

      b = rmax / (dexp(a*nr-a)-1d0)
      rofi2 = b*(exp(a)-1)
      nrl = dlog(rmax/rofi2)/a
      nrl = min(max(((nrl-1)/2)*2+1,nr),nrmx)

C      call info8(31,0,0,
C     .  ' mshnrl: rmax=%d  a=%;6d  nr=%i -> log mesh, nr=%i'//
C     .  '%?;n==0;*;;  rofil(1)=%;3g',rmax,a,nr,nrl,nrl-nrmx,rofi2,0,0)

      end
