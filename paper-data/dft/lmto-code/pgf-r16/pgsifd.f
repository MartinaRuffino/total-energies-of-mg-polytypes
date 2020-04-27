      subroutine qpgsifd(ld,sl,spl,wkl,wk2l,lg00,maxit,ndg,g00l)
C- Makes semi-infinite Green's function by decimation
C ----------------------------------------------------------------
Ci Inputs
Ci   ld    :leading dimension of s,sp,g00
Ci   s     :-1 * layer hamiltonian
Ci   sp    :complex work array of same dimension as s
Ci   wk    :complex work array of dimension ld*ld
Ci   wk2   :real*16 work array of dimension 3*ld
Ci   lg00   1 Generate left  semi-infinite GF from s00,s0L,s0R
Ci          2 Generate right semi-infinite GF from s00,s0L,s0R
Ci   maxit :maximum number of iterations to attempt
Ci   ndg   :second dimension of g00
Co Outputs
Ci   s     :is OVERWRITTEN on output
Co   g00   :semi-infinite GF
Cr Remarks
Cr   Decimation works iteratively by solving for all odd g_2n+1,0
Cr   in terms of the even g_2n,0 and g_2n+2,0.  Thus the equation
Cr     (P - S_00) g_00 - S_01 g_10 = 1
Cr   becomes
Cr    (P - S'_00) g_00  - S'_01 g_20  = 1
Cr    (P - S''_00) g_00 - S''_01 g_40 = 1
Cr   ...
Cr   and is repeated until the second term becomes negligible.
Cr   The scheme starts with the equations
Cr     (P - S_00) g_00 - S_01 g_10 = 1
Cr     -S_10 g_n-1,0 + (P - S_11) g_n0 - S_01 g_n+10 = 0 for all n
Cr   At the start, S_11 = S_00.
Cr   Every odd equation is used to knock out all the g_2n+1:
Cr     g+2n+1,0 = (P - S_11)^-1 (S_10 g_2n,0 + S_01 g_2n+2,0)
Cr   The equations involving only the g_2n look just as those with g_n
Cr   provided that the S_00, S_11, S_01 and S_10 are redefined:
Cr     S_00' =  S_00 - S_01 (P - S_11) S_10
Cr     S_11' =  S_11 + S_10 (P - S_11) S_01 + S_01 (P - S_11) S_10
Cr     S_10' =  S_10 (P - S_11) S_10
Cr     S_01' =  S_01 (P - S_11) S_01
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ld,maxit,ndg,lg00
      real(8) :: wk2l(ld,3),sl(ld,ld,3,2),g00l(ld,ndg,2),spl(ld,ld,3,2),wkl(ld,ld,2)
C ... Local parameters
C     logical :: debug = .false.
      integer i,j,iter,info,lgunit,iprint,ld2,L,R
      real*16 det(4),tol,zero,one
      real(8) :: ddet(4)
      real*16, allocatable :: wk2(:,:),s(:,:,:,:),g00(:,:,:),sp(:,:,:,:),wk(:,:,:)
      parameter (tol=1d-16)

C     print *, '!!' ; call setpr(70)

C     Copy to local arays with higher precision
      allocate(wk2(ld,3),s(ld,ld,3,2),g00(ld,ndg,2),sp(ld,ld,3,2),wk(ld,ld,2))
      s = sl
      g00 = g00l
      sp = spl

      call tcn('pgsifd')
      ld2 = ld*ld
      zero = 0; one = 1

C --- Put (P-S) in g_00 and S_11 ---
      call qpcopy(s,g00,1,ld2,-one)
      call qpcopy(s(1,1,1,2),g00(1,1,2),1,ld2,-one)
      call qpcopy(g00,s,1,ld2,one)
      call qpcopy(g00(1,1,2),s(1,1,1,2),1,ld2,one)
C      if (debug) then
C      call yprm('(P-s00)',2,s,ld*ld*3,ld,ld,ld)
C      endif

      L = 2
      R = 3

C --- For right GF, transpose S(2) and S(3) ---
      if (lg00 == 2) then
        do  j = 1, ld
          call qswap(ld,s(1,j,2,1),1,s(1,j,3,1),1)
          call qswap(ld,s(1,j,2,2),1,s(1,j,3,2),1)
        enddo
      endif

C --- Start of decimation loop ---
C     Start with P-S_00 in s(1) S_01 in s(2), S_10 in s(3)
      do  iter = 1, maxit
C   ... sp(1,*) <- (P-S11) (temporary storage)
        do  i = 1, ld
        do  j = 1, ld
          sp(i,j,1,1) = s(i,j,1,1)
          sp(i,j,1,2) = s(i,j,1,2)
        enddo
        enddo
C       call yprm('(P-s11)',2,sp,ld*ld*3,ld,ld,ld)

C   ... sp(1,*) <- (P-s11)^-1 (temporary storage)
        call qyygefa(sp,sp(1,1,1,2),ld,ld,wk2,info)
        if (info /= 0) call fexit(-1,111,' Exit -1 PGSIFD: '//
     .    'matrix singular, iteration %i',iter)
        call qyygedi(sp,sp(1,1,1,2),ld,ld,wk2,det,wk2(1,2),wk2(1,3),1)
C       call yprm('(P-s11)^-1',2,sp,ld*ld*3,ld,ld,ld)

C   ... sp(3,*) <- s01 (P-s11)^-1 (temporary storage)
        call qyygemm('N','N',ld,ld,ld,one,s(1,1,L,1),s(1,1,L,2),ld,
     .    sp,sp(1,1,1,2),ld,zero,sp(1,1,3,1),sp(1,1,3,2),ld)
C       call yprm('s01 (P-s11)^-1',2,sp(1,1,3,1),ld*ld*3,ld,ld,ld)

C   ... wk <- s01 (P-s11)^-1 s10 (renormalization of S_00)
        call qyygemm('N','N',ld,ld,ld,one,sp(1,1,3,1),sp(1,1,3,2),ld,
     .    s(1,1,R,1),s(1,1,R,2),ld,zero,wk,wk(1,1,2),ld)
C       call yprm('s01 (P-s11)^-1 s10',2,wk,ld*ld,ld,ld,ld)

C   ... sp(2,*) <- s01 (P-s11)^-1 s01 (updated s01)
        call qyygemm('N','N',ld,ld,ld,one,sp(1,1,3,1),sp(1,1,3,2),ld,
     .    s(1,1,L,1),s(1,1,L,2),ld,zero,sp(1,1,2,1),sp(1,1,2,2),ld)
C       call yprm('s01 (P-s11)^-1 s01',2,sp(1,1,2,1),ld*ld*3,ld,ld,ld)

C   ... sp(3,*) <- s10 (P-s11)^-1 (temporary storage)
        call qyygemm('N','N',ld,ld,ld,one,s(1,1,R,1),s(1,1,R,2),ld,
     .    sp,sp(1,1,1,2),ld,zero,sp(1,1,3,1),sp(1,1,3,2),ld)
C       call yprm('s10 (P-s11)^-1',2,sp(1,1,3,1),ld*ld*3,ld,ld,ld)

C   ... sp(1,*) <- s10 (P-s11)^-1 s01 (temporary storage)
        call qyygemm('N','N',ld,ld,ld,one,sp(1,1,3,1),sp(1,1,3,2),ld,
     .    s(1,1,L,1),s(1,1,L,2),ld,zero,sp,sp(1,1,1,2),ld)
C       call yprm('s01 (P-s11)^-1 s10',2,sp,ld*ld*3,ld,ld,ld)

C   ... RMS value of s01 (P-s11)^-1 s10
C   ... sp(1,*) <- s11 - s01 (P-s11)^-1 s10 - s10 (P-s11)^-1 s01
C   ... g00 <- s00 - s01 (P-s11)^-1 s10 (updated s11, and s00)
C   ... wk  <- sp(3,*) = s10 (P-s11)^-1 (temporary storage)
        call qyydotc(ld2,wk,wk(1,1,2),1,wk,wk(1,1,2),1,det,det(2))
        do  j = 1, ld
        do  i = 1, ld
          sp(i,j,1,1) = s(i,j,1,1) - sp(i,j,1,1) - wk(i,j,1)
          sp(i,j,1,2) = s(i,j,1,2) - sp(i,j,1,2) - wk(i,j,2)
          g00(i,j,1) = g00(i,j,1) - wk(i,j,1)
          g00(i,j,2) = g00(i,j,2) - wk(i,j,2)
          wk(i,j,1) = sp(i,j,3,1)
          wk(i,j,2) = sp(i,j,3,2)
        enddo
        enddo
C       call yprm('updated s11',2,sp,ld*ld*3,ld,ld,ld)

C   ... sp(3,*) <- s10 (P-s11)^-1 s10 (updated s10)
        call qyygemm('N','N',ld,ld,ld,one,wk,wk(1,1,2),ld,
     .    s(1,1,R,1),s(1,1,R,2),ld,zero,sp(1,1,3,1),sp(1,1,3,2),ld)
C       call yprm('s10 (P-s11)^-1 s10',2,sp(1,1,3,1),ld*ld*3,ld,ld,ld)

C   ... RMS s10, s11
        det(2) = 0
        det(3) = 0
        do  j = 1, ld
        do  i = 1, ld
          det(2) = det(2) + sp(i,j,2,1)**2 + sp(i,j,2,2)**2
          det(3) = det(3) + sp(i,j,3,1)**2 + sp(i,j,3,2)**2
        enddo
        enddo
        do  i = 1, 3
          ddet(i) = sqrt(det(i)/ld2)
          det(i) = sqrt(det(i)/ld2)
        enddo

        if (iprint() >= 70 .or. iter == maxit) then
          ddet = det
          call awrit3(' pgfsid iter %i:  rms delta S00 = %1;4g, '//
     .    'rms S01,S10 =%2:1;4g',' ',80,lgunit(1),iter,ddet(1),ddet(2))
        endif
C   ... Copy updated s10,s01,s11 to new positions (g00 already updated)
        call qcopy(ld2*3*2,sp,1,s,1)
        if (det(1)+det(2)+det(3) < tol) exit
        if (iter == maxit) call info0(2,0,0,' PGSIFD (warning): GF not converged to tolerance ...')
      enddo

C --- Semi-infinite Green's function ---
      call qyygefa(g00,g00(1,1,2),ld,ld,wk2,info)
      if (info /= 0) call fexit(-1,111,' Exit -1 PGSIFD: '//
     .    'Green''s function matrix singular',0)
      call qyygedi(g00,g00(1,1,2),ld,ld,wk2,det,wk2(1,2),wk2(1,3),1)

C      if (debug) then
C        call yprm('surface GF',2,g00,ld*ndg,ld,ld,ld)
C      endif

      g00l = g00
      call tcx('pgsifd')
      deallocate(wk2,s,g00,sp,wk)

      end
