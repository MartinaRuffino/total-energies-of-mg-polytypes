      subroutine mixpqb(nclass,nl,lmx,nsp,pnu,qnu,pold,qold,fnam,alph,
     .  a,namx,rms2,rmsdel,jac,beta,gamma,iter,broy,wc)
C- Mix to make a new vector of P and Q with modified Broyden mixing
C ------------------------------------------------------------------
Ci Inputs:
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   nl    :number of l's
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   nsp   :=1 spin degenerate, =2 non-degenerate
Cio  pnu   :principal quantum number
Cio  qnu   :moments
Ci   pold  :old principal quantum number
Ci   qold  :old moments
Ci   alph  :mixing parameter for linear mixing (only first iteration)
Ci   a     :workspace containing p and q
Ci   namx  :upper limit to number of elements in the vector
Ci   rmsdel:convergence criterion
Ci   broy  :T for modified Broyden mixing, F for Anderson mixing
Ci   wc    :mixing parameter (see routine Broyden).
Ci          wc<0: set wc_true = abs(wc) * RMS DQ, where RMS DQ is the
Cr          inverse of the RMS difference of the input and output vector
Co Outputs:
Cio  pnu,qnu: principal quantum number and moments
Co   rms2  :rms change in Q0
Co   jac   :Jacobian matrix
Co   beta  :define jac
Co   gamma :define jac
Co   iter  :total number of iterations
Cr Remarks:
Cr   The method requires a first estimation for the Jacobian matrix
Cr   if no information from previous iterations is available. Therefore,
Cr   for the first iteration, linear mixing is applied to calculate a
Cr   first jacobian matrix J^(1) as well as beta^(1) and gamma^(1).
Cr   Notations (according to PRB 30, 6118 (1984), App.):
Cr     x^(m): input vector for iteration m
Cr     F^(m): difference between output and input vector in iteration m
Cr     For the linear mixing, also the terms "in" and "out" must be
Cr     used.
Cr   Original version written by C. T. Chan, Ames Laboratory.
C   (D. Vanderbilt, S. G. Louie in Phys. Rev B 30, 6118 (1984))
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nclass,nl,nsp,namx,iter,lmx(nclass)
      character*(*) fnam
      double precision pnu(0:nl-1,nsp,nclass),pold(0:nl-1,nsp,nclass),
     .               qnu(3,0:nl-1,nsp,nclass),qold(3,0:nl-1,nsp,nclass),
     .               a(0:namx-1,0:1,2),jac(0:namx-1,0:namx-1),
     .               beta(0:namx-1,0:namx-1),gamma(0:namx-1,0:namx-1),
     .               alph,rmsdel(2),rms2,rmswc,wc
      logical broy
C ... Dynamically allocated local arrays
      integer, allocatable :: kpvt(:)
      real(8), allocatable :: bt(:)
      real(8), allocatable :: btx(:)
C ... Local parameters
      integer na,ic,isp,l,i,j,n1,nelts,stdo,
     .        fopn,fopnn,iprint,idum,info,iter0,lgunit,awrite
      double precision wctrue,diff,ddot
      character*80 outs
      logical lmix
      procedure(integer) :: nglob

C --- Convergence test (preserve Broyden file if convergence met) ---
      stdo = nglob('stdo')
      na = 0
      rms2 = 0
      do  ic = 1, nclass
        do  isp = 1, nsp
          do  l = 0, lmx(ic)
            na = na+1
            rms2 = rms2 + (qold(1,l,isp,ic) - qnu(1,l,isp,ic))**2
          enddo
        enddo
      enddo
      rms2 = dsqrt(rms2/na)
      if (4*na > namx) call rx('MIXPQB: too many elements')

C --- Determine wc_true if wc < 0 ---
      if (wc < 0) then
        rmswc = 0
        do  ic = 1, nclass
          do  isp = 1, nsp
            do  l = 0, lmx(ic)
              rmswc = rmswc + (pold(l,isp,ic) - pnu(l,isp,ic))**2
              do  i = 1, 3
                rmswc = rmswc + (qold(i,l,isp,ic) - qnu(i,l,isp,ic))**2
              enddo
            enddo
          enddo
        enddo
        wctrue = -wc/dsqrt(rmswc)
      else
        wctrue = wc
      endif

      iter = 1
      n1 = fopn(fnam)
      rewind n1
      read(n1,err=30,end=30) lmix, iter, idum, nelts
      if (nelts /= namx) then
        print *, 'MIXPQB: file mismatch ... do linear mix'
        goto 30
      endif

C---- Open file as new if incompatible ---
      if (.not. (lmix .and. broy)) then
        call fclose(n1)
        n1 = fopnn(fnam)
        rewind n1
        iter = 1
      else
        if (nelts /= namx) call rx('MIXPQB:  file mismatch')
      endif
   30 continue

C --- Printout ---
      if (iprint() >= 20) then
        print *, ' '
        j = awrite(' MIXPQB:  iteration %i read from file',
     .    outs,len(outs),0,iter-1,0,0,0,0,0,0,0)
        outs(j+2:len(outs)) = fnam
        call awrit1('%a.  RMS delta=%1,3;3e',outs,80,0,rms2)
        if (rmsdel(1) /= 0)
     .  call awrit1('%a  last it=%1,3;3e',outs,80,0,rmsdel)
        do  j = 1, 2
          call awrit0('%a',outs,-len(outs),-lgunit(j))
        enddo
        if (iter == 1) call query('beta',4,alph)
      endif
      rmsdel(1) = rms2
      rmsdel(2) = alph

C --- Modified Broyden mixing ---
      if (iter > 1) then
C ...   a(na,1,1): F^(m)  a(na,1,1): x^(m)
C ...   Save iter, set to 1 to revert to linear mix in case error
        iter0 = iter
        iter = 1
        read (n1,err=15,end=15) a
        read (n1,err=15,end=15) beta
        read (n1,err=15,end=15) gamma

C       F^(m+1) -> a(na,0,1)  x^(m+1) -> a(na,0,2)
C       dF^(m)  -> a(na,1,1)  dx^(m)  -> a(na,1,2)
        na = 0
        do  ic = 1, nclass
          do  isp = 1, nsp
            do  l = 0, lmx(ic)
              a(na,0,1) = pnu(l,isp,ic) - pold(l,isp,ic)
              a(na,0,2) = pold(l,isp,ic)
              a(na,1,1) = a(na,0,1) - a(na,1,1)
              a(na,1,2) = a(na,0,2) - a(na,1,2)
              do  i = 1, 3
                a(na+i,0,1) = qnu(i,l,isp,ic) - qold(i,l,isp,ic)
                a(na+i,0,2) = qold(i,l,isp,ic)
                a(na+i,1,1) = a(na+i,0,1) - a(na+i,1,1)
                a(na+i,1,2) = a(na+i,0,2) - a(na+i,1,2)
              enddo
              na = na+4
            enddo
          enddo
        enddo
        allocate(bt(namx)); call dpzero(bt,namx)
        allocate(btx(namx)); call dpzero(btx,namx)
C  ...  Broyden updates J, making J^(m+1)
        call broydn(jac,a(0,1,1),a(0,1,2),beta,gamma,bt,btx,
     .               namx,wctrue)
        deallocate(bt,btx)
C  ...  mv F^(m+1),x^(m+1) to old F^(m),x^(m) posns [a(*,1,1),a(*,1,2)]
        call dcopy (namx,a(0,0,1),1,a(0,1,1),1)
        call dcopy (namx,a(0,0,2),1,a(0,1,2),1)
C   ... Solve the linear equations (A6), storing y in a(*,0,1)
C       J^(m+1) * y = F^(m+1)  with  y = x^(m+2) - x^(m+1)
C       call defi(okpvt,namx)
        allocate(kpvt(namx))
        call dgefa(jac,namx,namx,kpvt,info)
        if (info /= 0) then
          print *, 'MIXPQB: Jacobian matrix singular ... do linear mix'
          goto 20
        endif
        call dgesl(jac,namx,namx,kpvt,a(0,0,1),0)
        deallocate(kpvt)

C   ... Make x^(m+2) = y + x^(m+1), store in a(na,0,1)
        na = 0
        do  ic = 1, nclass
          do  isp = 1, nsp
            do  l = 0, lmx(ic)
              pnu(l,isp,ic) = a(na,0,1) + a(na,1,2)
              a(na,0,1) = pnu(l,isp,ic)
              a(na,0,2) = a(na,1,1) + a(na,1,2)
              do  i = 1, 3
                qnu(i,l,isp,ic) = a(na+i,0,1) + a(na+i,1,2)
                a(na+i,0,1) = qnu(i,l,isp,ic)
                a(na+i,0,2) = a(na+i,1,1) + a(na+i,1,2)
              enddo
              na = na+4
            enddo
          enddo
        enddo
        iter = iter0
        alph = 1
C ---   Handle error; continue with linear mixing
   15   continue
        if (iter == 1)
     .    print *, 'MIXPQB: file read error ... do linear mix'
      endif
   20 continue

C --- Linear mixing with mixing alpha ---
C     x_out^(1)->a(na,0,2)  x_in^(1)->a(na,1,2)
C     x-in^(2) ->a(na,0,1)  F^(1)   ->a(na,1,1)
      if (iter == 1) then
        na = 0
        do  ic = 1, nclass
          do  isp = 1, nsp
            do  l = 0, lmx(ic)
              a(na,0,2) = pnu(l,isp,ic)
              a(na,1,1) = pnu(l,isp,ic) - pold(l,isp,ic)
              a(na,1,2) = pold(l,isp,ic)
              a(na,0,1) = alph*a(na,0,2) + (1-alph)*a(na,1,2)
              pnu(l,isp,ic) = a(na,0,1)
              do  i = 1, 3
                a(na+i,0,2) = qnu(i,l,isp,ic)
                a(na+i,1,1) = qnu(i,l,isp,ic) - qold(i,l,isp,ic)
                a(na+i,1,2) = qold(i,l,isp,ic)
                a(na+i,0,1) = alph*a(na+i,0,2) + (1-alph)*a(na+i,1,2)
                qnu(i,l,isp,ic) = a(na+i,0,1)
              enddo
              na = na+4
            enddo
          enddo
        enddo


C ---   Calculation of the first Jacobian matrix ---
C ...   See Eq. A6; for beta, eq. A14; for gamma, eq. A15.
C ...   x^(2) = x^(1) + [J^(1)]**(-1) * F^(m), so
C       J^(1) = F^(1) / [x^(2) - x^(1)]
C       Also, beta^(1) = 1 and gamma^(1) = J^(1)
        do  i = 0, namx-1
          diff = dabs(a(i,0,1) - a(i,1,2))
          if (diff > 1.d-12 .and. a(i,1,1) /= 0) then
            gamma(i,i) = a(i,1,1) / (a(i,0,1) - a(i,1,2))
          else
            gamma(i,i) = 1
            a(i,1,1) = 0
          endif
          beta(i,i) = 1
        enddo

      endif

C --- Printout ---
      if (iprint() > 30) then
        if (iter == 1) then
          call awrit3(' MIXPQB:  linear mixing  beta=%1;6d  nelts=%i'//
     .      '  rmsdel=%1;3e',' ',80,stdo,alph,na,rmsdel)
        else
          call awrit4(' MIXPQB:  Broyden mixing  iter=%i  nelts=%i'//
     .      '  rmsdel=%1;3e  wc=%1;3g',' ',80,stdo,iter,nelts,
     .      rmsdel,wctrue)
        endif
        if (iprint() > 40) then
          write (*,1)
    1     format(14x,'Old',11x,'New',10x,'Mixed',10x,'Diff'/)
          do  i = 0, na-1
            write (*,2) i,a(i,1,2),a(i,0,2),a(i,0,1),a(i,0,1)-a(i,1,2)
    2       format(i5,4F14.6)
            if (mod(i,4) == 3) write(stdo,*)
          enddo
        endif
      endif

      rms2 =  dsqrt(dabs(ddot(namx,a,1,a,1)
     .  -2*ddot(namx,a,1,a(0,1,2),1)
     .  + ddot(namx,a(0,1,2),1,a(0,1,2),1))/namx)

C --- Store newly mixed data in mixing file ---
      if (iprint() >= 40) call query('continue',-1,0)
      rewind n1
      iter = iter+1
      write(n1) broy, iter, idum, namx
      write(n1) a
      write(n1) beta
      write(n1) gamma
      call fclose (n1)

      end
