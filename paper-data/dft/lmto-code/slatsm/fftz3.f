C#define FFTW3
      subroutine fftz30(n1,n2,n3,k1,k2,k3)
C- Three-dimensional double complex-to-complex FFT
C ----------------------------------------------------------------------
Ci Inputs and parameters for fft3, fft30, fft3s
Ci    c(k1,k2,k3)   array to transform; overwritten on output (fft3)
Ci    n1,n2,n3      number of elements of the FFT
Ci    k1,k2,k3      dimensions of array c
Ci    isig         -1 for forward transform (real to recip)
Ci                  1 for inverse (assemble an real-space mesh)
Ci                 10s digit add 1 to suppress  / (n1*n2*n3) scaling
Ci    iset          a switch specifying information about setup;
Ci                  see Remarks.
Co Outputs
Cr   See Remarks.
Co   Forward transform: f(q) = sum_x f(x) exp(-i q x) / (n1*n2*n3)
Co   Reverse transform: f(q) = sum_x f(x) exp( i q x)
Cr Remarks
Cr  This code is a generic interface for several FFT routines.
Cr    fftz30(n1,n2,n3, k1,k2,k3)
Cr      Select the minimum dimensions k1,k2,k3 needed for fftz3.
Cr      This was included as an insurance in case some future
Cr      implementation requires k1,k2,k3 different from n1,n2,n3.
Cr      None of the existing ones here do, though some require
Cr      k1,k2,k3 are equal to n1,n2,n3.
Cr
Cr    fftz3s(n1,n2,n3,k1,k2,k3, iset)
Cr      This is an optional, implementation-specific setup routine.
Cr      Some implementations an initialization of internal coefficients
Cr      which need not repeated for subsequent calls.  For repeated
Cr      FFTs, the setup need not be repeated.
Cr      Calling fftz3s with iset=0 allocates memory for, and
Cr      initializes these coeficients.  To free this memory,
Cr      call fftz3s with the value of iset returned by the initial call.
Cr
Cr    fftz3(c,n1,n2,n3,k1,k2,k3,nfft,iset,isig)
Cr      This is the FFT routine.
Cr      Call with isig=-1 for forward, isig=1 for reverse transform:
Cr      Call with isig=-11 for forward, isig=11 for reverse transform,
Cr      swapping order of normalization (n1*n2*n3)
Cr      Forward transform: f(q) = sum_x f(x) exp(-i q x) / (n1*n2*n3)  isig=-1
Cr      Reverse transform: f(q) = sum_x f(x) exp( i q x)               isig=1
Cr      Forward transform: f(q) = sum_x f(x) exp(-i q x)               isig=-11
Cr      Reverse transform: f(q) = sum_x f(x) exp( i q x) / (n1*n2*n3)  isig=11
Cr      Call with iset=0 if fftz3s has not been invoked.
Cr
Cr    fftz3c(c,f,n1,n2,n3,k1,k2,k3,isw,isig)
Cr      Convolution of function c with f using the FFT routine.
C ----------------------------------------------------------------------
      implicit none
      integer n1,n2,n3,k1,k2,k3
C#ifdefC GOEDECKER
C#ifdefC SGI | SGI8
C      k1 = n1+1
C      k2 = n2+1
C      k3 = n3*2
C#elseifC DECA
C      k1 = n1
C      k2 = n2
C      k3 = n3*2
C#elseC
C      k1 = n1
C      k2 = n2
C      k3 = n3*2
C#endifC
C#elseifC ESSL
Cc ... For info, here are the strides suggested by essl
CC      call stride(n2, n1, inc2x, 'c', 0)
CC      call stride(n3, n2*inc2x, inc3x, 'c', 0)
CC      write(6,300) n1,n2,n3,inc2x,inc3x
CC  300 format('n1,n2,n3=',3i5,'    stride  suggests inc2,inc3=',2i5)
Cc ... but we just take k1=n1+1
C      k1 = n1+1
C      k2 = n2
C      k3 = n3
C#else
      k1 = n1
      k2 = n2
      k3 = n3
      call save_fftk(k1,k2,k3)
C#endif
      end

      subroutine fftz3(c,n1,n2,n3,k1,k2,k3,nfft,iset,isig)
      implicit none
      integer n1,n2,n3,k1,k2,k3,nfft,iset,isig
      double complex c(k1,k2,k3,nfft)
C Local variables
      integer i1,i2,i3,ifft
      double precision scale
C     save ow1,ow2,ow3,oiwk

C#ifdefC ESSL
C      integer w(1)
C      common /w/ w
C      integer naux,isign,inc2,inc3
C      call tcn('fftz3')
C
C      call rx('FFTZ3: this version does not seem to work')
Cc ... Work array
C      naux = 60000
C      if (max(n2,n3) >= 252) call rx('fftz3(essl): bad naux')
C      call defrr(ow1, naux)
C
Cc ... Do the FFT
C      isign = -isig
C      inc2 = k1
C      inc3 = k1*k2
C      scale = 1d0
C      if (isig < 0) scale = 1d0/dble(n1*n2*n3)
C      if (isig < -10) scale = 1
C      do  10  ifft = 1, nfft
C      call dcft3(c(1,1,1,ifft),inc2,inc3,c,inc2,inc3,n1,n2,n3,isign,
C     .           scale,w(ow1),naux)
C 10   continue
C      call rlse(ow1)
C
CC ... Renormalize forward transform ... already done by scale
CC      if (isig == 1 .or. isig == -11) goto 99
CC      scale = 1/dble(n1*n2*n3)
CC      do  20  ifft = 1, nfft
CC      do  20  i3 = 1, n3
CC      do  20  i2 = 1, n2
CC      do  20  i1 = 1, n1
CC   20 c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)
C
C#elseif FFTW | FFTW3
C
C --- A public-domain fft package.  See http://www.fftw.org/ ---
C ... Start of include file fftw_f77.i that comes with the fftw package
c     This file contains PARAMETER statements for various constants
c     that can be passed to FFTW routines.  You should include
c     this file in any FORTRAN program that calls the fftw_f77
c     routines (either directly or with an #include statement
c     if you use the C preprocessor).
      integer FFTW_FORWARD,FFTW_BACKWARD
      parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

C#ifdef FFTW3
      integer FFTW_ESTIMATE,FFTW_MEASURE
      parameter (FFTW_ESTIMATE=64,FFTW_MEASURE=0)
      INTEGER FFTW_PRESERVE_INPUT
      parameter (FFTW_PRESERVE_INPUT=16)
C#endif

      integer FFTW_THREADSAFE
      parameter (FFTW_THREADSAFE=128)
C ... End of include file fftw_f77.i that comes with the fftw package
      integer plan(2),jsig

      call tcn('fftz3')
      if (n1 == 1 .and. n2 == 1 .and. n3 == 1) goto 99

C#ifdef FFTW3
      jsig = 0
      if (isig == -1 .or. isig == -11) jsig = FFTW_FORWARD
      if (isig == 1  .or. isig == 11) jsig = FFTW_BACKWARD
      do  10  ifft = 1, nfft
      if (n2 == 1 .and. n3 == 1) then
        call dfftw_plan_dft_1d(plan,n1,c(1,1,1,ifft),c(1,1,1,ifft),
     .    jsig,FFTW_ESTIMATE)
      elseif (n3 == 1) then
        call dfftw_plan_dft_2d(plan,n1,n2,c(1,1,1,ifft),c(1,1,1,ifft),
     .    jsig,FFTW_ESTIMATE)
      else
        call dfftw_plan_dft_3d(plan,n1,n2,n3,c(1,1,1,ifft),c(1,1,1,ifft),
     .    jsig,FFTW_ESTIMATE)
      endif
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)
   10 continue
C#endif

C ... Renormalize forward transform
      if (isig == 1 .or. isig == -11) goto 99
      scale = 1/dble(n1*n2*n3)
      if (n1 == k1 .and. n2 == k2 .and. n3 == k3) then
        call zdscal(n1*n2*n3*nfft, scale, c, 1)
      else
        do  20  ifft = 1, nfft
        do  20  i3 = 1, n3
        do  20  i2 = 1, n2
        do  20  i1 = 1, n1
   20   c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)
      endif

C#endif

C     call zprm3('c',0,c,n1,n2,n3)
   99 call tcx('fftz3')
      return

      entry fftz3s(n1,n2,n3,k1,k2,k3,iset)

C#ifdef FFTW | FFTW3
C     Nothing needed using FFTW
C#elseC
C      if (iset == 0) then
C        id = k1
C        iopt = 0
C        if (n3 < 32) iopt = 1
C        iord = 1
C        call defrr(ow1, 4*n1+14)
C        call defrr(ow2, 4*n2*((1-iopt)+iopt*(id+1))+14)
C        call defrr(ow3, 4*n3+14)
C        call defrr(oiwk,max(n1,n2,n3))
C        call c3fft(c,id,n1,n2,n3,w(ow1),w(ow2),w(ow3),iopt,0,
C     .    iord,w(oiwk),ierr)
C        if (ierr == 2) call rx('c3fft needs prime factors 2,3,5')
C        if (ierr /= 0) call rxi('c3fft called improperly, err=',ierr)
C        iset = 1
C      else
C        call rlse(ow1)
C        iset = 0
C      endif
C#endif

      end

      subroutine gtpfac(npfac,pfac)
C- Returns allowed prime factors in integers for uniform mesh
      implicit none
      integer npfac,pfac(5)
C#ifdef SGI | SGI8 | FFTW | FFTW3
      npfac = 5
      pfac(1) = 2
      pfac(2) = 3
      pfac(3) = 5
      pfac(4) = 7
      pfac(5) = 11
C#elseC
C      npfac = 3
C      pfac(1) = 2
C      pfac(2) = 3
C      pfac(3) = 5
C#endif
      end

      subroutine fftz3c(c,f,n1,n2,n3,k1,k2,k3,isw,isig)
C- Convolution of double complex function with a filter by FFT
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   1s digit FT of input f,c
Ci         1  do not FT input f
Ci         2  do not FT input c
Ci         3  combination of 1+2
Ci        10s digit FT of output f,c
Ci         1  on output, do not inverse FT f, but return f~
Ci         2  on output, do not inverse FT c, but return c~
Ci         3  combination of 1+2
Ci        100s digit FFT type; see Remarks
Ci         0  type 0
Ci         1  type 1
Ci isig    sense of 'forward' and 'inverse' transforms
Cr         Use isig=-1 for usual definition of 'forward'
Cr            f~ = sum_x f(x) exp(-i q x) / (n1*n2*n3)
Cr         Use isig=1 for reverse definition:
Cr            f~ = sum_x f(x) exp( i q x)
Co Outputs
Co   c is overwritten by (c~(q) * f~(-q))~, where ~ is FT
Cr Remarks
Cr   Definition of Convolution in 1D and discrete analog with n points.
Cr   fftz3 handles two kinds: Basic definition (type 1) and Alternative definition (type 0)
Cr   Let F, G, H be the FT of f, g, h.
Cr
Cr                                   --- Basic definition ---
Cr                          Convolution                                    FT
Cr    Analytic  h(x)              = int[ du f(u) g(x-u) ]              H(k) G(k)
Cr    Discrete  h(1+cyclic(-i,n)) = sum_j f(1+j) g(1+cyclic(i+j,n))    H(1+i) G(1+i)
Cr              i and j loop between 0 and n-1.  Offset of 1 required by fortran convention
Cr
Cr                               --- Alternative definition ---
Cr                      Convolution                                        FT
Cr    Analytic  h(x)      = int[ du f(u) g(u-x) ]                      H(k) G(-k)
Cr    Discrete  h(i)      = sum_j f(1+j) g(1+cyclic(-i+j,n))           H(1+i) G(1+cyclic(i))
Cr
Cr  +Definition of cyclic = cyclic modulus of length n:
Cr   If i>=0 returns mod(i,n). If i<0 and not a multiple of n, returns n-mod(-i,n)
Cr   Cycles i in a continuous loop between 0 and n-1 for all i.
Cr      i cyclic(i,n)
Cr     -n   0
Cr      ...
Cr     -1   n-1
Cr      0   0
Cr      1   1
Cr      ...
Cr      n   0
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,k1,k2,k3,isw,isig
      double complex c(k1,k2,k3),f(k1,k2,k3)
C ... Local parameters
      integer isw0,isw1,isw2,i1,i1m,i2,i2m,i3,i3m
C     integer i,j,cyclic
C     cyclic(i,j) = i-floor(dble(i)/j)*j  ! Cyclic modulus: value betw/ 0 and j

      isw0 = mod(isw,10)
      isw1 = mod(isw/10,10)
      isw2 = mod(isw/100,10)

      if (isw0/2 == 0)      call fftz3(c,n1,n2,n3,k1,k2,k3,1,0,isig)
      if (mod(isw0,2) == 0) call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,isig)

C ... Alternative definition
      if (isw2 == 0) then
        do  i3 = 1, n3
C         i3m = 1+cyclic(1-i3,n3)  ! Construction below is is equivalent to cyclic modulus
          i3m = n3+2-i3; if (i3 == 1) i3m = 1
          do  i2 = 1, n2
C           i2m = 1+cyclic(1-i2,n2)
            i2m = n2+2-i2; if (i2 == 1) i2m = 1
            do  i1 = 1, n1
C             i1m = 1+cyclic(1-i1,n1)
              i1m = n1+2-i1; if (i1 == 1) i1m = 1
              c(i1,i2,i3) = c(i1,i2,i3)*f(i1m,i2m,i3m)
            enddo
          enddo
        enddo
C ... Basic definition
      elseif (isw2 == 1) then
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              c(i1,i2,i3) = c(i1,i2,i3)*f(i1,i2,i3)
            enddo
          enddo
        enddo

      else
        call rx('fftz3 : unknown FT type')
      endif

      if (isw1/2 == 0)      call fftz3(c,n1,n2,n3,k1,k2,k3,1,0,-isig)
      if (mod(isw1,2) == 0) call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-isig)
      end

      subroutine suconvolve1D(opt,npts,io,lp,npad,n1,src,dest)
C- Cyclic copy, possibly with pad
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt    :1s digit
Ci          : 0 src is real(8)
Ci          : 1 src is complex(8)
Ci          :10s digit
Ci          : 0 Forward map: copy src -> dest
Ci          : 1 Reverse map: copy dest -> src
Ci          :100s digit
Ci          : 0 order points as given
Ci          : 1 reverse order of elements in dest:
Ci          :   dest(i) <-> dest(npts+npad+1-i)
Ci          :   For forward map this is done on exit
Ci          :   For reverse map this is done on entry
Ci   npts   :number of points
Ci   io     :Shift of origin; see Remarks
Ci   lp     :Number of zero or positive elements on abscissa; see Remarks.
Ci          :npts-lp  is the number of negative elements
Ci          :Used only to insert padding after point lp
Ci   npad   :number of padding elements
Ci   n1     :Not used now.  Will become leading dimension of src,dest.
Ci   fsrc   :source function, to copy, shift and pad
Co Outputs
Co   dest  :destination function, copy and shift
Cr Remarks
Cr   This routine prepares arrays for convolution h(x)=∫du f(u) g(x−u) by FFT
Cr   Note that if F, G, H are the FT of f, g, h, then H(k) = F(k) G(k).
Cr   To manage convolution by FFT, there are several considerations:
Cr     1. First pont may not be origin.  The origin  must be shifted to first point
Cr        io = offset to element containing origin (i.e. origin at point io+1)
Cr     2. Array may not be periodic.  To mitigate this, pad region far from origin (npad)
Cr        lp is number of positive points
Cr        If npad>=npts, there is no error from replica?
Cr     3. May require convolution h(x)=∫du f(u) g(x+u)
Cr        Then g(x+u) should be mapped to g(-x-u).  Convolution will yield h(-x).
Cr
Cr   To accomplish this, src must be copied to dest with possible shift, pad, and swap x->0x
Cr   dest is copied from src with cyclic shift, and padding with npad pts after point lp
Cr
Cr                              Cyclic shift without padding
Cr      a        b          c              npts    b          c            npts
Cr       ____________________________________       ________________________
Cr      |        |          |                |     |          |             |
Cr      |<--io-->|<-- lp -->|<--npts-lp-io-->| =>  |<-- lp -->|<--npts-lp-->|
Cr      |        |          |                |     |          |             |
Cr       ------------------------------------       ------------------------
Cr         neg        pos          neg                 pos          neg
Cr
Cr                              Cyclic shift with padding
Cr      a        b          c              npts    b          c                      npts+npad
Cr       ____________________________________       _____________________________________
Cr      |        |          |                |     |          |            |             |
Cr      |<--io-->|<-- lp -->|<--npts-lp-io-->| =>  |<-- lp -->|<-- npad -->|<--npts-lp-->|
Cr      |        |          |                |     |          |            |             |
Cr       ------------------------------------       -------------------------------------
Cr         neg        pos          neg                 pos          neg
Cu Updates
Cu   29 May 17  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,npts,io,npad,lp,n1
      complex(8) :: src(npts),dest(npts+npad)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: wk(:)
C ... Local parameters
      integer i,is,id,opt0,opt1,opt2

C ... Setup
      if (npad < 0) call rx('suconvolve1D given bad padding')
      if (n1 /= 1) call rx('suconvolve1D not ready for n1')
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      opt2 = mod(opt/100,10)

C ... Reverse order of dest
      if (opt2 == 1 .and. opt1 == 1) then
        allocate(wk(npts+npad)); call zcopy(npts+npad,dest,1,wk,1)
        forall (id = 1:npts+npad) dest(id) = wk(npts+npad+1-id)
        deallocate(wk)
      endif

C ... Copy src to work array
      allocate(wk(npts)); call dpzero(wk,2*npts)
      if (opt0 == 0 .and. opt1 == 0) then
        call dcopy(npts,src,1,wk,2)
      elseif (opt1 == 0) then
        call dcopy(2*npts,src,1,wk,1)
      endif

C ... Copy with shift + pad
      id = 0
      do  i = 1, npts
        is = mod(i-1+io,npts)+1
        id = id+1
C        print 333, i,is,id
C  333   format(3i4)
        if (opt1 == 0) then
          dest(id) = wk(is)
        else
          wk(is) = dest(id)
        endif
        if (i == lp) then  ! Padding region
          if (opt1 == 0) then
          forall (id = lp+1:lp+npad) dest(id) = 0
          endif
          id = lp+npad
        endif
      enddo

C ... Reverse order of dest
      if (opt2 == 1 .and. opt1 == 0) then
        deallocate(wk)
        allocate(wk(npts+npad)); call zcopy(npts+npad,dest,1,wk,1)
        forall (id = 1:npts+npad) dest(id) = wk(npts+npad+1-id)
        deallocate(wk)
      endif

C ... Copy work array to src
      if (opt0 == 0 .and. opt1 == 1) then
        call dcopy(npts,wk,2,src,1)
      elseif (opt1 == 1) then
        call dcopy(2*npts,wk,1,src,1)
      endif
      if (allocated(wk)) deallocate(wk)

      end
