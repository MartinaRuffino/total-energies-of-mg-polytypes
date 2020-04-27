C test yqinv
C     To check many cases, eg
C     foreach n ( `mlist 150:188` )
C     a.out << END | grep -E 'overwritten|factor'
C
C     $n
C
C     END
C     end

      subroutine fmain
C     implicit none
      INTEGER   lda, ldb, ldc, ldw, n, i, j, i0, j0,
     .  opw, opw2, wksize, nlev, nmin, i1mach
      logical a2bin,cmdopt
      double precision mflops
      PARAMETER (lda=901,ldb=901,ldc=901,wksize=6000000)
      character*1 strn*10, cs, outs*80
      double precision
     .  sa( lda, lda, 2 ), pa( lda, lda, 2), pw( lda, lda, 2)
      double precision s_time, temp, s_mflops, q_time, q_mflops
      double precision cpusec,xx,errmx
      integer w(wksize)
      common /w/ w
      common /static/ sa, pa, pw

C --- Setup ---
      call finits(2,0,0,i)
      call wkinit(lda*lda*4)
      call tcinit(0,4)
      call nada
      n = 500
      nlev = 2
      cs = 's'
      print *, 'cs= (n or h)?'
      read(*,'(a1)') cs
      call initqu(.true.)
      print *, 'test yqinv: cs=',cs
      call query('n=?',2,n)
      if (n > lda) call rx('n gt lda')
      mflops = 2.0d-6 * n * n * n * 4

      j = 6
      if (cmdopt('-nmin=',j,0,outs)) then
        if (.not. a2bin(outs,nmin,2,0,' ',j,72))
     .    call rxs('cannot parse',outs)
        call yqinv0(nmin,1)
      endif

C#ifdefC gfortran
C      nmin = 32
C#else
      call yqinv0(nmin,0)
C#endif

C --- Matrix inversion, yygefa,yygedi ---
      CALL init ( sa,lda,n, cs )
C     call ywrm(0,' ',2,6,'(8f16.10)',sa,lda*lda,lda,n,n)
      temp = cpusec()
      if (cs == 'h') then
C        call yyhifa(sa,sa(1,1,2),lda,n,pw,i)
C        if (i /= 0) call rx(' matrix is singular')
C        call yyhidi(sa,sa(1,1,2),lda,n,pw,xx,i,pw(1,2,1),pw(1,3,1),1)
C        do  16  i = 1, n
C        do  16  j = 1, i
C        sa(i,j,1) =  sa(j,i,1)
C   16   sa(i,j,2) = -sa(j,i,2)

        call ztoyy(sa,lda,lda,n,n,0,1)
        temp = cpusec()
        call tcn('zhetr*')
        call zhetrf('L',n,sa,lda,pw,pa,lda*lda,i)
        if (i /= 0) call rx(' matrix is singular')
        call zhetri('L',n,sa,lda,pw,pa,i)
        call tcx('zhetr*')
        xx = cpusec()
        call ztoyy(sa,lda,lda,n,n,1,0)
        temp = temp + cpusec() - xx
        do  16  i = 1, n
        do  16  j = 1, i
        sa(j,i,1) =  sa(i,j,1)
   16   sa(j,i,2) = -sa(i,j,2)
        strn = 'zhetrf,i'
      else
        call ztoyy(sa,lda,lda,n,n,0,1)
        temp = cpusec()
        call tcn('zgetr*')
        call zgetrf(n,n,sa,lda,pw,i)
        if (i /= 0) call rx(' matrix is singular')
        call zgetri(n,sa,lda,pw,pa,lda*lda,i)
        call tcx('zgetr*')
        xx = cpusec()
        call ztoyy(sa,lda,lda,n,n,1,0)
        temp = temp + cpusec() - xx
        strn = 'zgetrf,i'
C        call yygefa(sa,sa(1,1,2),lda,n,pw,i)
C        if (i /= 0) call rx(' matrix is singular')
C        call yygedi(sa,sa(1,1,2),lda,n,pw,xx,pw(1,2,1),pw(1,3,1),1)
      endif
      s_time = cpusec() - temp
      if (s_time == 0d0) s_time = 1
      s_mflops = mflops / s_time
*     call yprm(.false.,'a^-1',2,6,'(5f16.10)',sa,lda,n,lda,n)

C --- Inversion by yqinv ---
      CALL init ( pa,lda,n, cs )
C ... (debugging) check usage of w
      ldw = min(n+5,lda)
      write(*,'('' using lda,ldw,n,nlev,nmin='',4i4)') lda,ldw,nlev,nmin
      call defrr(opw,  ldw*(n+2))
      call defrr(opw2, n*n)
      call dvset(w(opw),1,ldw*(n+2),-99d0)
      call dvset(w(opw2),1,n*n,-99d0)

      temp = cpusec()
      call yyqinv(cs,pa,pa(1,1,2),lda,nlev,n,w(opw),ldw,i)
C     call yqinv('0'//cs,pa,lda*lda,lda,nlev,n,w(opw),ldw,i)

      q_time = cpusec() - temp
      if (q_time == 0d0) q_time = 1
      q_mflops = mflops / q_time

C  ...Check explicitly hermitian
      if (cs == 'h') then
        errmx = 0d0
        i0 = 0
        j0 = 0
        do 10 i = 1, n
        do 10 j = 1, n
          if (max(dabs(pa(i,j,1)-pa(j,i,1)),
     .            dabs(pa(i,j,2)+pa(j,i,2))) > errmx) then
            i0 = i
            j0 = j
          endif
   10   errmx = max(errmx,dabs(pa(i,j,1)-pa(j,i,1)),
     .                    dabs(pa(i,j,2)+pa(j,i,2)))
        print 100, 'checking hermitian errmax ',errmx,i0,j0
  100   format(1X,a,':  errmx=',g10.3,2x,2i4)
      endif

C --- Check quality of inversion ---
      CALL diff('compare '//strn//'to yqinv', sa, pa, lda, n )
      CALL init ( sa,lda,n, cs )
      temp = cpusec()
      call yygemm('N','N',n,n,n,1d0,sa,sa(1,1,2),lda,pa,pa(1,1,2),lda,
     .  0d0,pw,pw(1,1,2),lda)
      temp = cpusec() - temp
      if (temp == 0d0) temp = 1
      do  40  i = 1, n
   40 pw(i,i,1) = pw(i,i,1) - 1
c     call yprm(.false.,'a^1 * a - 1',2,6,'(16f8.2)',pw,ldw,n,ldw,n)
      errmx = 0
      do 42 i = 1, n
      do 42 j = 1, n
   42 errmx = max(errmx,dabs(pw(i,j,1)),dabs(pw(i,j,2)))
      print 100, 'inv(a) * a - 1     errmax ',errmx,nint(mflops/temp)

      print 110, s_time, s_mflops, q_time, q_mflops, s_time/q_time,n

  110 format(/1X,'zgetrf time: ',F7.3,'  zgetrf MFlops:',F7.0,
     .       /1X,' yqinv time: ',F7.3,'   yqinv MFlops:',F7.0
     .       /1X,'     factor: ',F7.3,'  for n =',i4)


C ... Make sure pw not overwritten beyond ldw*(n+1)
      call wkchk('check if pw2 still intact')
      call dcopy(n*n,w(opw2),1,pw,1)
      j = 1
      do  50  i = 1, n*n
        if (pw(i,j,1) /= -99d0) call fexit2(-1,111,
     .    ' STOP: pw overwritten, i,j=%i %i',i,j)
   50 continue
C     This will fail if yqinv calls LAPACK
      if (ldw > n) then
        call dmcpy(w(opw),ldw,1,pw,lda,1,ldw,n+1)
        do  52  i = 1, ldw
        do  52  j = 1, n+1
            if (i <= n .and. j <= n+1) goto 52
            if (pw(i,j,1) /= -99d0) call fexit2(-1,111,
     .        ' STOP: outer block of pw overwritten, i,j=%i %i',i,j)
   52   continue
      endif
      call tcprt(i1mach(2))

      end

      SUBROUTINE init ( a, lda, n, cs)
C- Initialize arrays
      integer lda, n
      double precision a( lda, lda, 2)
      character *1 cs
      integer i,j

      real ran1
      call ran1in(1)

      do 10 i = 1, n
      do 10 j = 1, n
      a(i,j,1) = ran1()
   10 a(i,j,2) = ran1()

      if (cs == 'h') then
        do 20 i = 1, n
        do 20 j = 1, i
        a(i,j,1) =  a(j,i,1)
   20   a(i,j,2) = -a(j,i,2)
        do  22 i = 1, n
   22   a(i,i,2) = 0
      endif

      end
      SUBROUTINE diff (strn, sc, pc, ldc, n)
C- Compare the two arrays for differences
C     implicit none
      character*(*) strn
      integer ldc, n, i, j
      double precision sc(ldc,ldc,2), pc(ldc,ldc,2), errmx
      errmx = 0d0
      do 10 i = 1, n
      do 10 j = 1, n
   10 errmx = max(errmx,dabs(sc(i,j,1)-pc(i,j,1)),
     .                  dabs(sc(i,j,2)-pc(i,j,2)))
      print 100, strn,errmx
  100 format(1X,a,':  errmx=',g10.3)
      end
