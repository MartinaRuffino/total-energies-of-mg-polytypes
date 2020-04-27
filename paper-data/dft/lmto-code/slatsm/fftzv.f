C#define FFTW3
      subroutine fftzv(c,n1,n2,n3,n4,nfft,iset,isig)
C-  Multidimensional FFT based on the FFTW package
C ----------------------------------------------------------------------
Ci Inputs
Ci    n1,..n4:number of elements of the FFT
Ci           :n4 = 1 => 3D FFT
Ci           :n3 = 1 => 2D FFT
Ci           :n2 = 1 => 1D FFT
Ci    iset   :Switch specifying information about setup
Ci           :1s digit:
Ci           :0 reset FFTW "wisdom" to original state
Ci           :1 make a plan for given parameters; execute no FFT
Ci           :  plan is NOT DESTROYED on exit
Ci           :2 execute FFT using already existing plan
Ci           :  plan is NOT DESTROYED on exit
Ci           :3 make a plan for given parameters; execute FFT
Ci           :  plan is DESTROYED on exit
Ci           :4 destroy plan
Ci           :10s digit: tells planner what  FFTW "wisdom" to build
Ci           :0 use FFTW_ESTIMATE (fast setup, less efficient execution)
Ci           :1 use FFTW_MEASURE (slower setup, more efficient execution)
Ci           :2 use FFTW_PATIENT (slow setup, most efficient execution)
Ci    isig   :-1 for forward transform (real to recip)
Ci           :1 for inverse (assemble an real-space mesh)
Ci           :10s digit add 1 to suppress  / (n1*n2*n3) scaling
Cio Inputs/Outputs
Ci    c      :c(n1..n4,nfft) = array to transform; nfft FTs are carried out
Cl Local variables
Cl         :
Cr Remarks
Cr   Multiple plans:
Cr   fftzv can be called with multiple plans, without reinitializing plans
Cr   You can extract or set a plan with fftzvp.
Cr
Cr   Calling fftv3 with multiple FFT's:
Cr   You can create a plan to generate many fft's at one blow,
Cr   or create a plan for one FFT, than use the plan with repeated
Cr   calls to fftzv.  There doesn't seem to be much efficiency gain
Cr   in generating all FFTs in one shot.
Cr
Cr   Compatibility when using FFT without FFTW:
Cr   If FFTW is not used, fftzv will merely call fftz3 with iset=0
Cr   if 1s digit of iset is 2 or 3; otherwise it will do nothing
Cr   It is assumed that fftz3 args k1,k2,k3 are n1,n2,n3.
Cr
Cu Updates
Cu   06 Jul 10
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,n4,nfft,iset,isig
      double complex c(n1,n2,n3,n4,nfft)
C ... Local parameters
      integer i1,i2,i3,i4,ifft,N,plan(2),iset0,iset1
      integer isw,plan0(2),k1,k2,k3,k1s,k2s,k3s
      double precision scale
      integer flags
      integer jsig,rank,n1234(4)
      save plan,k1s,k2s,k3s
C ... This is fftw3.f from FFTW package
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
      INTEGER FFTW_WISDOM_ONLY
      PARAMETER (FFTW_WISDOM_ONLY=2097152)
C ... end of fftw3.f

C#ifndefC FFTW3
C      iset0 = mod(iset,10)
C      if (iset0 /= 2 .and. iset0 /= 3) return
C      call fftz3(c,n1,n2,n3,k1s,k2s,k3s,nfft,0,isig)
C      return
C#else

      iset0 = mod(iset,10)
      iset1 = mod(iset/10,10)
      if (iset0 == 0) then
        call dfftw_forget_wisdom()
        return
      elseif (iset0 == 4) then
        call dfftw_destroy_plan(plan)
        return
      endif
      if (n1 == 1) return

      call tcn('fftzv')

      rank = 4
      if (n4 == 1) rank = 3
      if (n3 == 1) rank = 2
      if (n2 == 1) rank = 1
      n1234(1) = n1
      n1234(2) = n2
      n1234(3) = n3
      n1234(4) = n4
      N = n1*n2*n3*n4

C --- Make a plan ---
      if (mod(iset0,2) /= 0) then
        jsig = 0
        if (isig == -1) jsig = FFTW_FORWARD
        if (isig == 1) jsig = FFTW_BACKWARD

        flags = FFTW_ESTIMATE
        if (iset1 == 1) flags = FFTW_MEASURE
        if (iset1 == 2) flags = FFTW_PATIENT

        call dfftw_plan_many_dft(plan,rank,n1234,nfft,
     .                           c,n1234,1,N,
     .                           c,n1234,1,N,
     .                           jsig, flags)
      endif
      if (iset0 == 1) goto 99

C --- Execute plan ---
      if (iset0 == 2) then
        call dfftw_execute_dft(plan,c,c)
      elseif (iset0 == 3) then
        call dfftw_execute(plan)
        call dfftw_destroy_plan(plan)
      endif

C --- Renormalize forward transform --
      if (isig > 0 .or. isig < -10) goto 99
      scale = 1/dble(n1*n2*n3*n4)
      do  20  ifft = 1, nfft
      do  20  i4 = 1, n4
      do  20  i3 = 1, n3
      do  20  i2 = 1, n2
      do  20  i1 = 1, n1
   20 c(i1,i2,i3,i4,ifft) = scale*c(i1,i2,i3,i4,ifft)

   99 call tcx('fftzv')
      return
C#endif

      entry fftzvp(isw,plan0)
C isw = 0 : retrieve internal plan, return in plan0
C isw > 0 : set internal plan from passed plan0
      if (isw > 0) then
        plan(1) = plan0(1)
        plan(2) = plan0(2)
      else
        plan0(1) = plan(1)
        plan0(2) = plan(2)
      endif
      return

C ... This entry serves no function for fftw
C     For non-fftw implementations, k1,k2,k3 are preserved
C     to enable compatibility with fftz3.
      entry save_fftk(k1,k2,k3)
      k1s = k1
      k2s = k2
      k3s = k3

      end
