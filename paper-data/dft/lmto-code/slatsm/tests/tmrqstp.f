C#ifdef TEST
C     Test routine mrqstp
      subroutine fmain
      implicit none

C     test1 declarations: read from file 77
C      integer nfit,npar
C      integer,allocatable:: ivar(:)
C      real(8),allocatable :: efit(:),par(:),epar(:),dpar(:,:)
C      real(8),allocatable :: wk(:,:),alp(:,:),cov(:,:),sigx(:)
C      double precision alam,chi(2)


      integer nvar
      double precision alsc

C     for test supplied by NR
      INTEGER NPT,MA
      double precision SPREAD
      PARAMETER(NPT=100,MA=6,SPREAD=0.001d0)
      INTEGER i,ia(MA),idum,ipass,itst,j,k,mfit,iter
      double precision alamda,chisq(2),fgauss,gasdev,ochisq,x(NPT),
     .  y(NPT),fy(NPT),sig(NPT),a(MA),covar(MA,MA),alpha(MA,MA),
     .  gues(MA),dfyda(MA,NPT)
      double precision wk2(MA,5)
      EXTERNAL fgauss
      DATA a /5d0,2d0,3d0,2d0,5d0,3d0/
      DATA gues/4.5d0,2.2d0,2.8d0,2.5d0,4.9d0,2.8d0/

C ... Test1: data is read from fort.77
C      read(77) nfit,npar,nvar
C      allocate(ivar(npar),efit(nfit),par(npar))
C      read(77) ivar
C      read(77) efit
C      read(77) par
C
C      allocate(wk(npar,5),alp(npar,npar),cov(npar,npar),sigx(nfit))
C      allocate(epar(nfit),dpar(nvar,nfit))
C      alsc = 5; alam = 1; sigx = 1; iter = 0
C
C      do  i = 1, 5
C      read(77) epar,dpar
C      call mrqstp(nfit,npar,nvar,ivar,efit,epar,dpar,sigx,wk,par,
C     .  alsc,alp,alam,cov,chi,iter)
C      enddo
C

C     call setpr(60)

      idum = -911
C     first try a sum of two Gaussians
      do  12  i = 1, NPT
        x(i) = 0.1d0*i
        y(i) = 0d0
        do  11  j = 1, MA,3
          y(i) = y(i)+a(j)*dexp(-((x(i)-a(j+1))/a(j+2))**2)
11      continue
        y(i) = y(i)*(1d0 + SPREAD*gasdev(idum))
        sig(i) = SPREAD*y(i)
12    continue
      mfit = MA
      do  13  i = 1, mfit
        ia(i) = i
13    continue
      do  14  i = 1, MA
        a(i) = gues(i)
14    continue
      alsc = 10; nvar = MA

C     1st pass: no constraints; 2nd pass: constraints
      do  16  ipass = 1, 2
C       alamda = -1
        iter = 0; alamda = 0.001d0
        call getdyda(NPT,MA,nvar,x,a,ia,fy,dfyda,fgauss)
        call mrqstp(NPT,MA,nvar,ia,y,fy,dfyda,sig,wk2,a,
     .    alsc,alpha,alamda,covar,chisq,iter)

        k = 1
        itst = 0
C       re-entry for iteration
    1   write(*,'(/1x,a,i2,t18,a,f10.4,t43,a,e9.2)') 'Iteration #',k,
     .       'Chi-squared:',chisq(1),'ALAMDA:',alamda
        write(*,'(1x,t5,a,t13,a,t21,a,t29,a,t37,a,t45,a)') 'A(1)',
     .       'A(2)','A(3)','A(4)','A(5)','A(6)'
        write(*,'(1x,6f8.4)') (a(i),i=1,6)
        k = k+1
        ochisq = chisq(1)

        call getdyda(NPT,MA,nvar,x,a,ia,fy,dfyda,fgauss)
        call mrqstp(NPT,MA,nvar,ia,y,fy,dfyda,sig,wk2,a,
     .    alsc,alpha,alamda,covar,chisq,iter)

        if (chisq(1) > ochisq) then
          itst = 0
        else if (dabs(ochisq-chisq(1)) < 0.1d0) then
          itst = itst+1
        endif
        if (itst < 4) then
          goto 1
        endif

C       Final call with almda = 0 for covariance matrix

        alamda = 0d0
        call getdyda(NPT,MA,nvar,x,a,ia,fy,dfyda,fgauss)
        call mrqstp(NPT,MA,nvar,ia,y,fy,dfyda,sig,wk2,a,
     .    alsc,alpha,alamda,covar,chisq,iter)

        write(*,*) 'Uncertainties:'
        write(*,'(1x,6f8.4/)') (dsqrt(covar(i,i)),i=1,6)
        write(*,'(1x,a)') 'Expected results:'
        write(*,'(1x,f7.2,5f8.2/)') 5d0,2d0,3d0,2d0,5d0,3d0

        if (ipass == 1) then
C          write(*,*) 'press return to continue with constraint'
C          read(*,*)
          write(*,*) 'New test, holding a(2) and a(5) constant'
          do  15  j = 1, MA
            a(j) = a(j)+.1d0
   15     continue
          a(2) = 2d0
          a(5) = 5d0
          ia(1) = 1
          ia(2) = 3
          ia(3) = 4
          ia(4) = 6
C          ia(5) = 2
C          ia(6) = 5
          nvar = 4
        endif
   16 continue
      END
      subroutine getdyda(nfit,na,nvar,x,a,ia,y,dyda,funcs)
      implicit none
      integer nfit,na,nvar,ia(na)
      double precision a(na),x(nfit),y(nfit),dyda(nvar,nfit)
      double precision dydai(na)
      integer i,j
      EXTERNAL funcs

      do  i=1,nfit
        call funcs(x(i),a,y(i),dydai,na)
        do j = 1, nvar
          dyda(j,i) = dydai(ia(j))
        enddo
      enddo
      end

      SUBROUTINE fgauss(x,a,y,dyda,na)
      INTEGER na
      double precision x,y,a(na),dyda(na)
      INTEGER i
      double precision arg,ex,fac
      y=0d0
      do 11 i=1,na-1,3
        arg=(x-a(i+1))/a(i+2)
        ex=dexp(-arg**2)
        fac=a(i)*ex*2d0*arg
        y=y+a(i)*ex
        dyda(i)=ex
        dyda(i+1)=fac/a(i+2)
        dyda(i+2)=fac*arg/a(i+2)
11    continue
      return
      END

      FUNCTION gasdev(idum)
      INTEGER idum
      double precision gasdev
CU    USES ran1x
      INTEGER iset
      double precision fac,gset,rsq,v1,v2
      REAL ran1x
      SAVE iset,gset
      DATA iset/0/
      if (idum < 0) iset=0
      if (iset == 0) then
1       v1=2.*ran1x(idum)-1.
        v2=2.*ran1x(idum)-1.
        rsq=v1**2+v2**2
        if(rsq >= 1d0. or.rsq == 0d0)goto 1
        fac=dsqrt(-2d0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

      FUNCTION ran1x(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1x,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum <= 0.or.iy == 0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum < 0) idum=idum+IM
          if (j <= NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum < 0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1x=min(AM*iy,RNMX)
      return
      END

C#endif
