C#define EXIT
      subroutine mrqstp(nfit,ncof,nvar,ivar,dat,y,dy,sig,wk,
     .  cof,alsc,alp,alam,cov,chi,iter)
C- One Levenberg-Marquardt iteration for least-squares fit of function to data
C ----------------------------------------------------------------------
Ci Inputs
Ci   nfit  :number of points to fit
Ci   ncof  :total number of coefficients in fitting function;
Ci         :leading dimension of alp and cov
Ci   nvar  :number of coefficients to vary out of ncof parameters
Ci         :defining function.  The rest are frozen.
Ci         :ivar below points to list of 1:nvar parameters to be varied.
Ci         :Otherwise, mrqstp will check for consistency in array ivar.
Ci   ivar  :ivar(i) = index to ith parm to be varied in the cof vector
Ci   dat   :data to be fit
Ci   y     :current values of fitting function
Ci   dy    :derivatives of y wrt each of 1..nvar fitting coefficients
Ci   sig   :standard deviation sigma for each data point
Ci         :sig(i)=0 => point gets no weight (same as sigma = infinity)
Cio Inputs/Outputs
Cio  wk    :work array of dimension (ncof,5).  wk contains:
Cio        :wk(:,1) = beta vector; NR Eqns 15.5.8 and 15.5.6
Cio        :wk(:,2) = tryc: trial vector of cof
Cio        :wk(:,3) = dcof: change in cof
Cio        :wk(:,4) = pivot array in matrix inversion
Cio        :wk(:,5) = work array in matrix inversion
Cio        :wk should be preserved between successive calls to mrqstp
Cio  cof   :cof(1:ncof) contain coefficients to the fitting function
Cio        :Only a subset (nvar<ncof) will be varied; the elements
Cio        :that are varied is set by ivar above.
Cio        :On input, trial values for cofficients
Cio        :On output, updated trial values for coefficients
Cio        :cof must be preserved between successive calls to mrqstp
Cio  alsc  :alam scaling factor (Numerical Recipes uses alsc=10)
Cio        :If initially zero, alsc is set to 10
Cio  alp   :curvature (alpha) matrix, NR 15.5.8 and 15.5.11
Cio        :alp must be preserved between successive calls to mrqstp
Cio  alam  :L-M lambda parameter (see Remarks)
Cio        :On the first call: alam = initial value for lamdba
Cio        :If initially zero, alam starts with initial value of 0.001
Cio        :alam should be preserved between successive calls to mrqstp
Cio  iter  :number of iterations.
Cio        :On the first call, iter must be zero, which generates
Cio        :an initialization step.
Cio        :iter is incremented after each successive call.
Cio        :iter returned < 0 means an error has occured.
Cio        :iter = -1 => something wrong with ivar array
Cio        :iter = -2 => mrqstp could not invert alpha' matrix
Co Outputs
Co   cov   :covariance matrix; see Remarks
Co   chi   :chi(1) chi-squared value for fit with current parameters
Co         :chi(2) prior chi-squared value
Cl Local variables
Cr Remarks
Cr   This routine uses the Levenberg-Marquardt algorithm for
Cr   least-squares fitting of a nonlinear function of coefficients cof
Cr   to given data.  It closely follows Numerical Recipes 2nd edition,
Cr   Chapter 15.5.
Cr   Relation between variables in NR and those here:
Cr      NR      NR fortran       here
Cr     lamdba    alamda          alam
Cr     C         covar           cov (intermediate iter: holds alpha')
Cr     a         a               cof
Cr     alpha     alpha           alp
Cr     beta      beta            wk(:,1)
Cr               atry            wk(:,2)
Cr               da              wk(:,3)
Cr     N         ndata           nfit
Cr               nca             ncof
Cr     M         mfit            nvar
Cr               10              alsc
Cr               chisq           chi(1)
Cr               ochisq          chi(2)
Cr   Algorithm proceeds as follows.  For each iteration:
Cr     step      function
Cr      First iteration :
Cr      1        y,dy -> initial alpha, beta, chi
Cr      2        (no step 2)
Cr      3        alpha -> alpha' = alpha(1+lambda) (NR 15.5.13)
Cr               NB: alpha' stored in cov
Cr      4        alpha', beta -> trial cof (see mrqcof)
Cr               exit seeking y,dy of trial cof
Cr      Subsequent iterations have new y,dy, preserved :
Cr      1        y,dy -> trial alpha, beta, chi
Cr      2        If trial chi < old chi:
Cr                 replace cof,alpha,beta,chi with trial values
Cr                 decrease lambda
Cr               Otherwise:
Cr                 increase lambda
Cr               End of step 2
Cr      3        alpha -> alpha' = alpha(1+lambda) (NR 15.5.13)
Cr      4        alpha', beta -> trial cof (see routine mrqcof)
Cr   Calling this routine:
Cr     Iterate calls to this routine until convergence criteria is met.
Cr     (The caller must specify what that criterion is; see below)
Cr     A new set of trial cofs will be returned after each call.
Cr     Every iteration requires nfit,ncof,nvar,ivar,dat,y,dy,sig.
Cr     Set iter = 0 for the first iteration.
Cr     On each successive call, update y,dy; keep unchanged all other quantities
Cr
Cr     Example criteria:
Cr       1.  alam increases after e.g. 5 consecutive calls
Cr       2.  chi falls below a certain tolerance
Cu Updates
Cu   30 Jul 09  First created: adapted from NR mrqmin
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nfit,ncof,nvar,ivar(ncof),iter
      double precision alam,alsc,chi(2),alp(ncof,ncof),cov(ncof,ncof),
     .  cof(ncof),wk(ncof,5),sig(nfit),dat(nfit),y(nfit),dy(nvar,nfit)
C ... Local parameters
      integer j,k,kk,ihit,iprint,ierr,i1mach,stdo
      double precision rmsdel,dif
C ... External calls (calls to the latter three may be commented out)
C     external covsrt,dcopy,dgemm,dsidi,dsifa,dswap,mrqcof
C     external i1mach,info5,rx

C --- First iteration ---
      if (iter == 0) then

C   ... Setup
        if (alsc <= 0) alsc = 10
        if (alam <= 0) alam = 1d0
C   ... Ensure ivar contains a proper permutation of the parameters
        kk = nvar + 1
        do  j = 1, ncof
          ihit = 0
          do  k = 1, nvar
            if (ivar(k) == j) ihit = ihit + 1
          enddo
          if (ihit == 0) then
            ivar(kk) = j
            kk = kk + 1
          elseif (ihit > 1) then
C#ifdef EXIT
            call rx('MRQSTP: duplicate parameter in ivar')
C#endif
            iter = -1
            return
          endif
        enddo
        if (kk /= ncof+1) then
C#ifdef EXIT
          call rx('MRQSTP: Improper permutation in ivar')
C#endif
          iter = -1
        endif

C   ... 1. Make initial chi, alp, beta
C       call mrqcof(nfit,nvar,ncof,sig,dat,y,dy,chi,alp,beta)
        call mrqcof(nfit,nvar,ncof,sig,dat,y,dy,chi(1),alp,wk(1,1))
        chi(2) = chi(1)
        rmsdel = dsqrt(chi(1) / nfit)
C       call dcopy(ncof,cof,1,tryc,1)
        call dcopy(ncof,cof,1,wk(1,2),1)

C   ... Printout
C#ifdefC PRINT
C        call info5(30,1,0,'%N MRQSTP: Initial values: '
C     .    //'chi^2=%;3g  rms=%;4g  lambda=%;4g',chi,rmsdel,alam,0,0)
C#endif

C --- Subsequent iterations ---
      else
C       Re-entry: restore tryc, cof
C       call dswap(ncof,tryc,1,cof,1)
        call dswap(ncof,wk(1,2),1,cof,1)

C   ... 1. Make chi, cov, dcof for trial cof
        call mrqcof(nfit,nvar,ncof,sig,dat,y,dy,chi(1),cov,wk(1,3))

C   ... 2 Test whether trial case improved chi
        if (chi(1) < chi(2)*0.99999d0) then

C     ... Success, accept new solution and decrease alam
          alam = alam / alsc
C#ifdefC PRINT
C          call info5(30,1,0,' MRQSTP iter %i(-): chi^2=%,4;4g'//
C     .      '  dchi^2/chi^2=%,4;4g  lambda=%g',iter,chi,chi(1)/chi(2)-1,
C     .      alam,0)
C#endif
          chi(2) = chi(1)
          do  j = 1, nvar
            do  k = 1, nvar
              alp(j,k) = cov(j,k)
            enddo
C           beta(j) = dcof(j)
            wk(j,1) = wk(j,3)
C           cof(ivar(j)) = tryc(ivar(j))
            cof(ivar(j)) = wk(ivar(j),2)
          enddo
        else if (alam /= 0) then

C     ... Failure, increase alam
          alam = alsc*alam
C#ifdefC PRINT
C          call info5(30,1,0,' MRQSTP iter %i(+): chi^2=%,4;4g'//
C     .      '  dchi^2/chi^2=%,4;4g  lambda=%g',iter,chi,chi(1)/chi(2)-1,
C     .      alam,0)
C#endif
          chi(1) = chi(2)
        endif

      endif

C ... Printout
C#ifdefC PRINT
C      stdo = i1mach(2)
C      if (iprint() >= 50) then
C        write(stdo,500)
C        do  j = 1, nfit
C          dif = y(j) - dat(j)
C          write(stdo,510) j,dat(j),y(j),dif
C        enddo
C  500   format(/4x,'i',6x,'data',12x,'fit',13x,'diff')
C  510   format(i5,2f16.10,g20.10)
C      endif
C#endif

C --- 3. Scale linearized fitting matrix : C = alpha (1 + lambda) ---
      do  j = 1, nvar
        do  k = 1, nvar
          cov(j,k) = alp(j,k)
        enddo
        cov(j,j) = alp(j,j) * (1+alam)
C       Only needed if calling gaussj below
C       dcof(j) = beta(j)
      enddo

C --- 4. Solve linearized eqns for trial cof ---
C     Solver given in Numerical Recipes:
C     call gaussj(nvar,ncof,1,1,wk(1,4),wk(1,5),iwk3,cov,dcof)
C     call prmx('mrqstp after inversion: cov',cov,ncof,nvar,nvar)
C     Linpack analog:
      call dsifa(cov,ncof,nvar,wk(1,4),ierr)
      if (ierr /= 0) then
C#ifdef EXIT
        call rx('MRQSTP: failed to invert covariance matrix')
C#endif
        iter = -2
        return
      endif
      call dsidi(cov,ncof,nvar,wk(1,4),dif,j,wk(1,5),1)
      do  k = 1, nvar
      do  j = 1, k
        cov(k,j) = cov(j,k)
      enddo
      enddo
C     call dgemm('N','N',nvar,1,nvar,1d0,cov,ncof,beta,ncof,0d0,
C    .  dcof,ncof)
      call dgemm('N','N',nvar,1,nvar,1d0,cov,ncof,wk(1,1),ncof,0d0,
     .  wk(1,3),ncof)
C      call prmx('mrqstp after inversion: cov',cov,ncof,nvar,nvar)
C      call prmx('mrqstp: beta',wk(1,1),ncof,nvar,1)
C      call prmx('mrqstp: dcof',wk(1,3),ncof,nvar,1)

C ... If converged, evaluate covariance matrix with alam=0
      if (alam == 0) then
        call covsrt(nvar,ncof,ivar,cov)
        return
      endif

C ... Trial parameters from linearized solution
      do  j = 1, nvar
C       tryc(ivar(j)) = cof(ivar(j)) + dcof(j)
        wk(ivar(j),2) = cof(ivar(j)) + wk(j,3)
      enddo

C --- Prepare for next iteration, exit ---
      iter = iter+1
C     Swap new values into cof for caller to make new function call
C     call dswap(ncof,tryc,1,cof,1)
      call dswap(ncof,wk(1,2),1,cof,1)

      end

      subroutine mrqcof(nfit,nvar,ncof,sig,dat,y,dy,chi,alp,beta)
C- Linearized fitting matrix alpha and beta, Levenberg-Marquardt algorithm
C ----------------------------------------------------------------------
Ci Inputs
Ci   nfit  :number of data points to fit
Ci   nvar  :number of parameters to vary
Ci   ncof  :leading dimension of alp
Ci         :ncof=0 => just calculate chi; skip alp,beta
Ci   dat   :values to be fit by function
Ci   sig   :standard deviation for each data point
Ci   y     :current values of fitting function
Ci   dy    :derivatives of y wrt each of 1..nvar fitting parameters
Co Outputs
Co   chi   :chi-squared value for fit with current parameters, NR 15.5.5
Co   alp   :alpha matrix, NR 15.5.8 and 15.5.11
Co   beta  :beta vector, NR 15.5.8 and 15.5.6
Cr Remarks
Cr   Kernel called by mrqmin.
Cr   See Section 15.5, 'Nonlinear Models' in Numerical Recipes 2nd edition
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nfit,nvar,ncof
      double precision sig(nfit),dat(nfit),y(nfit),dy(nvar,nfit),
     .  chi,alp(ncof,ncof),beta(ncof)
C Local parameters
      integer i,j,k
      double precision del,sig2i,wt

C --- Initialize beta and (symmetric) alpha ---
      if (ncof > 0) then
        do  20  j = 1, nvar
          do  10  k = 1, j
            alp(j,k) = 0d0
   10     continue
          beta(j) = 0d0
   20   continue
      endif

C --- Loop over data and find chi-squared ---
      chi = 0d0
      do  50  i = 1, nfit
        del = dat(i) - y(i)
        if (sig(i) /= 0) then
          sig2i = 1d0 / (sig(i)*sig(i))
        else
          sig2i = 0
        endif
        if (ncof > 0) then
          do  40  j = 1, nvar
            wt = dy(j,i)*sig2i
            do  30  k = 1, j
              alp(j,k) = alp(j,k) + wt*dy(k,i)
   30       continue
            beta(j) = beta(j) + del*wt
   40     continue
        endif
        chi = chi + del*del*sig2i
   50 continue

C --- Fill in symmetric side of alpha ---
      if (ncof > 0) then
        do  70  j = 2, nvar
        do  60  k = 1, j - 1
          alp(k,j) = alp(j,k)
   60   continue
   70   continue
      endif

      end

      subroutine covsrt(nvar,ncof,ivar,cov)
C- Repack the covariance matrix to the true order
C ----------------------------------------------------------------------
Ci Inputs
Ci   nvar: number of parameters to vary out of ncof total parameters
Ci   ncof: total number of parameters (including fixed parameters)
Ci   ivar: ivar(i) = index to ith parm to be varied in the cof vector
Co Outputs
Co   cov: covariance matrix, repacked to true order on output
Cr Remarks
Cr   Elements of cov associated with fixed parameters set to zero.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nvar,ncof
      integer ivar(nvar)
      double precision cov(ncof,ncof)
C Local parameters
      integer i,j
      double precision swap

C --- Zero all elements below diagonal ---
      do  20  j = 1, ncof-1
        do  10  i = j+1, ncof
          cov(i,j) = 0d0
   10   continue
   20 continue

C --- Repack off-diagonal elements of fit into correct locations ---
      do  40  i = 1, nvar-1
        do  30  j = i+1, nvar
          if (ivar(j) > ivar(i)) then
            cov(ivar(j),ivar(i)) = cov(i,j)
          else
            cov(ivar(i),ivar(j)) = cov(i,j)
          endif
   30   continue
   40 continue

C --- Temporarily store original diagonal elements in top row ---
      swap = cov(1,1)
      do  50  j = 1, ncof
        cov(1,j) = cov(j,j)
        cov(j,j) = 0d0
   50 continue
      cov(ivar(1),ivar(1)) = swap

C --- Sort elements into proper order on diagonal ---
      do  60  j = 2, nvar
        cov(ivar(j),ivar(j)) = cov(1,j)
   60 continue

C --- Fill in above diagonal by symmetry ---
      do  80  j = 2, ncof
        do  70  i = 1, j-1
          cov(i,j) = cov(j,i)
   70   continue
   80 continue

      end


C#ifdefC TEST
CC     Test routine mrqstp
C      subroutine fmain
C      implicit none
C
CC     test1 declarations: read from file 77
CC      integer nfit,npar
CC      integer,allocatable:: ivar(:)
CC      real(8),allocatable :: efit(:),par(:),epar(:),dpar(:,:)
CC      real(8),allocatable :: wk(:,:),alp(:,:),cov(:,:),sigx(:)
CC      double precision alam,chi(2)
C
C
C      integer nvar
C      double precision alsc
C
CC     for test supplied by NR
C      INTEGER NPT,MA
C      double precision SPREAD
C      PARAMETER(NPT=100,MA=6,SPREAD=0.001d0)
C      INTEGER i,ia(MA),idum,ipass,itst,j,k,mfit,iter
C      double precision alamda,chisq(2),fgauss,gasdev,ochisq,x(NPT),
C     .  y(NPT),fy(NPT),sig(NPT),a(MA),covar(MA,MA),alpha(MA,MA),
C     .  gues(MA),dfyda(MA,NPT)
C      double precision wk2(MA,5)
C      EXTERNAL fgauss
C      DATA a /5d0,2d0,3d0,2d0,5d0,3d0/
C      DATA gues/4.5d0,2.2d0,2.8d0,2.5d0,4.9d0,2.8d0/
C
CC ... Test1: data is read from fort.77
CC      read(77) nfit,npar,nvar
CC      allocate(ivar(npar),efit(nfit),par(npar))
CC      read(77) ivar
CC      read(77) efit
CC      read(77) par
CC
CC      allocate(wk(npar,5),alp(npar,npar),cov(npar,npar),sigx(nfit))
CC      allocate(epar(nfit),dpar(nvar,nfit))
CC      alsc = 5; alam = 1; sigx = 1; iter = 0
CC
CC      do  i = 1, 5
CC      read(77) epar,dpar
CC      call mrqstp(nfit,npar,nvar,ivar,efit,epar,dpar,sigx,wk,par,
CC     .  alsc,alp,alam,cov,chi,iter)
CC      enddo
CC
C
CC     call setpr(60)
C
C      idum = -911
CC     first try a sum of two Gaussians
C      do  12  i = 1, NPT
C        x(i) = 0.1d0*i
C        y(i) = 0d0
C        do  11  j = 1, MA,3
C          y(i) = y(i)+a(j)*dexp(-((x(i)-a(j+1))/a(j+2))**2)
C11      continue
C        y(i) = y(i)*(1d0 + SPREAD*gasdev(idum))
C        sig(i) = SPREAD*y(i)
C12    continue
C      mfit = MA
C      do  13  i = 1, mfit
C        ia(i) = i
C13    continue
C      do  14  i = 1, MA
C        a(i) = gues(i)
C14    continue
C      alsc = 10; nvar = MA
C
CC     1st pass: no constraints; 2nd pass: constraints
C      do  16  ipass = 1, 2
CC       alamda = -1
C        iter = 0; alamda = 0.001d0
C        call getdyda(NPT,MA,nvar,x,a,ia,fy,dfyda,fgauss)
C        call mrqstp(NPT,MA,nvar,ia,y,fy,dfyda,sig,wk2,a,
C     .    alsc,alpha,alamda,covar,chisq,iter)
C
C        k = 1
C        itst = 0
CC       re-entry for iteration
C    1   write(*,'(/1x,a,i2,t18,a,f10.4,t43,a,e9.2)') 'Iteration #',k,
C     .       'Chi-squared:',chisq(1),'ALAMDA:',alamda
C        write(*,'(1x,t5,a,t13,a,t21,a,t29,a,t37,a,t45,a)') 'A(1)',
C     .       'A(2)','A(3)','A(4)','A(5)','A(6)'
C        write(*,'(1x,6f8.4)') (a(i),i=1,6)
C        k = k+1
C        ochisq = chisq(1)
C
C        call getdyda(NPT,MA,nvar,x,a,ia,fy,dfyda,fgauss)
C        call mrqstp(NPT,MA,nvar,ia,y,fy,dfyda,sig,wk2,a,
C     .    alsc,alpha,alamda,covar,chisq,iter)
C
C        if (chisq(1) > ochisq) then
C          itst = 0
C        else if (dabs(ochisq-chisq(1)) < 0.1d0) then
C          itst = itst+1
C        endif
C        if (itst < 4) then
C          goto 1
C        endif
C
CC       Final call with almda = 0 for covariance matrix
C
C        alamda = 0d0
C        call getdyda(NPT,MA,nvar,x,a,ia,fy,dfyda,fgauss)
C        call mrqstp(NPT,MA,nvar,ia,y,fy,dfyda,sig,wk2,a,
C     .    alsc,alpha,alamda,covar,chisq,iter)
C
C        write(*,*) 'Uncertainties:'
C        write(*,'(1x,6f8.4/)') (dsqrt(covar(i,i)),i=1,6)
C        write(*,'(1x,a)') 'Expected results:'
C        write(*,'(1x,f7.2,5f8.2/)') 5d0,2d0,3d0,2d0,5d0,3d0
C
C        if (ipass == 1) then
CC          write(*,*) 'press return to continue with constraint'
CC          read(*,*)
C          write(*,*) 'New test, holding a(2) and a(5) constant'
C          do  15  j = 1, MA
C            a(j) = a(j)+.1d0
C   15     continue
C          a(2) = 2d0
C          a(5) = 5d0
C          ia(1) = 1
C          ia(2) = 3
C          ia(3) = 4
C          ia(4) = 6
CC          ia(5) = 2
CC          ia(6) = 5
C          nvar = 4
C        endif
C   16 continue
C      END
C      subroutine getdyda(nfit,na,nvar,x,a,ia,y,dyda,funcs)
C      implicit none
C      integer nfit,na,nvar,ia(na)
C      double precision a(na),x(nfit),y(nfit),dyda(nvar,nfit)
C      double precision dydai(na)
C      integer i,j
C      EXTERNAL funcs
C
C      do  i=1,nfit
C        call funcs(x(i),a,y(i),dydai,na)
C        do j = 1, nvar
C          dyda(j,i) = dydai(ia(j))
C        enddo
C      enddo
C      end
C
C      SUBROUTINE fgauss(x,a,y,dyda,na)
C      INTEGER na
C      double precision x,y,a(na),dyda(na)
C      INTEGER i
C      double precision arg,ex,fac
C      y=0d0
C      do 11 i=1,na-1,3
C        arg=(x-a(i+1))/a(i+2)
C        ex=dexp(-arg**2)
C        fac=a(i)*ex*2d0*arg
C        y=y+a(i)*ex
C        dyda(i)=ex
C        dyda(i+1)=fac/a(i+2)
C        dyda(i+2)=fac*arg/a(i+2)
C11    continue
C      return
C      END
C
C      FUNCTION gasdev(idum)
C      INTEGER idum
C      double precision gasdev
CCU    USES ran1x
C      INTEGER iset
C      double precision fac,gset,rsq,v1,v2
C      REAL ran1x
C      SAVE iset,gset
C      DATA iset/0/
C      if (idum < 0) iset=0
C      if (iset == 0) then
C1       v1=2.*ran1x(idum)-1.
C        v2=2.*ran1x(idum)-1.
C        rsq=v1**2+v2**2
C        if(rsq >= 1d0. or.rsq == 0d0)goto 1
C        fac=dsqrt(-2d0*log(rsq)/rsq)
C        gset=v1*fac
C        gasdev=v2*fac
C        iset=1
C      else
C        gasdev=gset
C        iset=0
C      endif
C      return
C      END
C
C      FUNCTION ran1x(idum)
C      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
C      REAL ran1x,AM,EPS,RNMX
C      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
C     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
C      INTEGER j,k,iv(NTAB),iy
C      SAVE iv,iy
C      DATA iv /NTAB*0/, iy /0/
C      if (idum <= 0.or.iy == 0) then
C        idum=max(-idum,1)
C        do 11 j=NTAB+8,1,-1
C          k=idum/IQ
C          idum=IA*(idum-k*IQ)-IR*k
C          if (idum < 0) idum=idum+IM
C          if (j <= NTAB) iv(j)=idum
C11      continue
C        iy=iv(1)
C      endif
C      k=idum/IQ
C      idum=IA*(idum-k*IQ)-IR*k
C      if (idum < 0) idum=idum+IM
C      j=1+iy/NDIV
C      iy=iv(j)
C      iv(j)=idum
C      ran1x=min(AM*iy,RNMX)
C      return
C      END
C
C#endif
