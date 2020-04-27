      subroutine fmain
C-- tpzhev.f test pzhev.f
C-----------------------------------------------------------------------
C   Main program to test MPI parallel implementation of the BLACS and
C   SCALAPCK.
C
C   tpzhev.o needs to be linked against pzhev.o and the MPI slatsm.a
C   as well as LIBLOC, BLACS and SCALAPACK.
C
C   If the MPI lines are commented out using ccomp, then this routine
C   can be used to compare times for the serial routines in
C   LAPACK and EISPACK. Link against slatsm.a and LIBLOC
C
C   Options are,
C   -n=#      dimension of matrix
C   -nb=#     blocking factor for the BLACS process configuration
C             I've always found 16 significantly the optimum, and
C             this is the default
C   -nprow=# -npcol=# over ride the process configuration
C             not recommended (let MPI do it for you)
C   --inv     real matrix inversion instead of diagonalisation
C   --real    real matrix, otherwise complex
C   --g       generalised (ie, non orthogonal) eigenproblem
C   --test    check the parallel code gives the right answer
C   --mlog    each process # to write a log file to mlog.dat_#
C             gives details of process configuration and
C             memory alllocation. Process 0 (master) writes to
C             mlog.dat; look at that file to see wall times for
C             compute and distribution.
C   --serial  repeat the diagonalisation with serial LAPACK and
C             EISPACK (unit stride and generic) to compare the time
C   The following flags apply to the --serial option
C   --tinvit  Use inverse iteration in EISPACK
Cb  --test    does not apply to serial code. Only the LAPACK info is
Cb            printed out as a measure of success
C-------------------------------------------------------------------
      use mpi
      implicit none
      integer nbl,nrc,i,j,n,nevmx,nev,i1mach,ndmx,lwork,info,iov,iinv
      integer procid,master,ILAENV,lgunit,fopn,fext,fextg
      double precision eevmx,d1mach,cpusec,t1s,t2s,t1p,t2p,t1l,t2l,
     .                 t1x,t2x,work(3),abstol,DLAMCH
      logical test,twice,serial,l,lov,lc,invit,lx,linv,
     .        cmdopt,a2bin,T,F
      data T / .true. / F / .false. /
      character*6 Lname
      character*72 outs
      character*256 strn

C#ifdefC MPE
C      include "mpef.h"
C      integer MPE_LOG_GET_EVENT_NUMBER,MPE_DESCRIBE_STATE,
C     .        MPE_LOG_EVENT,MPE_INIT_LOG,MPE_FINISH_LOG
C#endif


      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS, dims(2)
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      logical mlog


C Heap allocation
      integer wksize
C     2GB memory request (4Byte per integer) ..
      parameter(wksize= 500 000 000)
      integer w(wksize)
C     Next two lines guarantee w is aligned along a d.p. boundary
      double precision ws
      equivalence (ws,w(1))
      common /w/ w

C Heap pointers
      integer oe,oh,os,ot,oz,ochk1,ochk2,owk,orwk,ozwk,oiwk,oifail

C MPI process configuration
      integer nb,nprow,npcol

C Initialise command line and heap
      call finits(2,0,0,i)
      call pshpr(0)
      call wkinit(wksize)
      call poppr
      call pshpr(30)

      procid = 0
      master = 0

      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs > 1) then
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      if (procid /= master) then
        call pshpr(0)
        call pshpr(0)
        do  i = 1, 4
          call sprt(i,0)
        enddo
      else
        call pshpr(30)
      endif
C#ifdefC MPE
C      ierr = MPE_INIT_LOG()
C#endifC
      mlog = cmdopt('--mlog',6,0,strn)
      i = fextg(ext)
      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
      if (procid == master) then
        call gettime(datim)
        if (mlog) i = fopn('MLOG')
        if (mlog) then
          call awrit2(' tpzhev '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' is master',' ',256,
     .      lgunit(3),procid,numprocs)
        endif
      else
        call strcat(ext,20,' ','_',1,' ',i)
        call bin2a(' ',0,0,procid,2,0,20,ext,i)
        ierr = fext(ext(1:i+1))
        if (mlog) ierr = fopn('MLOG')
        ierr = fextg(ext)
        call gettime(datim)
        if (mlog) then
          call awrit2(' tpzhev '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' file extension is '//ext(2:i+1),' ',
     .      256,lgunit(3),procid,numprocs)
        endif
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )

C MPI process configuration
      nb = 16
      nprow = -1
      npcol = -1
      endif

C Get dimension
      j = 3
      if (cmdopt('-n=',j,0,outs)) then
        if (a2bin(outs,n,2,0,' ',j,72)) goto 2
      endif
      if (procid == master) print *,
     .  'Usage: tpzhev --inv --real --g --test --twice --serial '//
     .  '--tinvit -n=# -nb=# -nprow=# -npcol=# [ext]'
      goto 3
    2 continue
      linv = F
      if (cmdopt('--inv',5,0,outs))    linv = T
      lc = T
      if (cmdopt('--real',6,0,outs))   lc = F
      test = F
      if (cmdopt('--test',6,0,outs))   test = T
      twice = F
      if (cmdopt('--twice',7,0,outs))  twice = T
      if (numprocs > 1) then
      serial = F
      if (cmdopt('--serial',8,0,outs)) serial = T
      else
      serial = T
      endif
      invit = F
      if (cmdopt('--tinvit',8,0,outs))  invit = T
      lov = F
      if (cmdopt('--g',3,0,outs))      lov = T
      if (numprocs > 1) then
      j = 4
      if (cmdopt('-nb=',j,0,outs))     l = a2bin(outs,nb,2,0,' ',j,72)
      j = 7
C MPI process configuration
      if (cmdopt('-nprow=',j,0,outs)) l = a2bin(outs,nprow,2,0,' ',j,72)
      j = 7
      if (cmdopt('-npcol=',j,0,outs)) l = a2bin(outs,npcol,2,0,' ',j,72)
      if (procid == master) then
        call awrit4('%N Starting .., n=%i nb=%i nprow=%i npcol=%i',' ',
     .              128,lgunit(1),n,nb,nprow,npcol)
      endif
      endif
      nevmx = n
      eevmx = 1d12
    4 continue
      if (linv) then
        lc = F
        lov = F
      endif
      if (lc) then
        nrc = 2
      else
        nrc = 1
      endif
      call defrr(oe, n)
      call defrr(oh, nrc*n*n)
      if (lov) then
        call defrr(os, nrc*n*n)
      else
        os = 5
      endif
      if (numprocs > 1) then
      if (linv) then
        call mkinv(n,w(oh))
      else
        if (lc) then
          call mkhsz(lov,n,w(oh),w(os))
        else
          call mkhsd(lov,n,w(oh),w(os))
        endif
      endif
      if (procid == master) t1p = MPI_WTIME()
      call pzhev(linv,lc,lov,n,oh,os,nb,nprow,npcol,eevmx,nevmx,nev,
     .           w(oe),ot)
      if (procid == master) t2p = MPI_WTIME()
      if (procid == master .and. .not. linv) then
        call awrit1(' pzhev found %i eigenvectors',' ',120,i1mach(2),
     .              nev)
      endif

C Test the result
      if (test .and. procid == master) then
        if (linv) then
          call mkinv(n,w(oh))
          call defrr(ochk1, n*n)
          call dmpy(w(oh),n,1,w(ot),n,1,w(ochk1),n,1,n,n,n)
          call checki(n,w(ochk1))
        else
          if (lc) then
            call mkhsz(lov,n,w(oh),w(os))
            call defcc(ochk1, n*n)
            call zmpy(w(oh),2*n,2,1,w(ot),2*n,2,1,w(ochk1),2*n,2,1,n,
     .                n,n)
            if (lov) then
              call defcc(ochk2, n*n)
              call zmpy(w(os),2*n,2,1,w(ot),2*n,2,1,w(ochk2),2*n,2,1,n,
     .                  n,n)
            else
              ochk2 = ot
            endif
            call checkz(n,w(ochk1),w(ochk2),w(oe))
          else
            call mkhsd(lov,n,w(oh),w(os))
            call defrr(ochk1, n*n)
            call dmpy(w(oh),n,1,w(ot),n,1,w(ochk1),n,1,n,n,n)
            if (lov) then
              call defrr(ochk2, n*n)
              call dmpy(w(os),n,1,w(ot),n,1,w(ochk2),n,1,n,n,n)
            else
              ochk2 = ot
            endif
            call checkd(n,w(ochk1),w(ochk2),w(oe))
          endif
        endif
      endif
      call rlse(ot)
      endif
C --- do serial calculations ---
      if (serial .and. procid == master) then
        if (linv) then
          call mkinv(n,w(oh))
          t1s = cpusec()
          call DPOTRF('U',n,w(oh),n,info)
          call awrit1(' DPOTRF returned info=%i',' ',96,i1mach(2),info)
          call DPOTRI('U',n,w(oh),n,info)
          call awrit1(' DPOTRI returned info=%i',' ',96,i1mach(2),info)
          call symm(n,w(oh))
          t2s = cpusec()
        else
          abstol = 2*DLAMCH('S')
          call defi(oiwk, 5*n)
          call defi(oifail, n)
          if (lov) then
            iov = 1
          else
            iov = 0
          endif
          if (invit) then
            iinv = 1
          else
            iinv = 0
          endif
          call defrr(owk, n*11)
          if (lc) then
            call defcc(oz, n*n)
            call mkhsz(lov,n,w(oh),w(os))
            call ztoyy(w(oh),n,n,n,n,1,0)
            if (lov) then
              call ztoyy(w(os),n,n,n,n,1,0)
            endif
            lx = F
            t1s = cpusec()
            call diagno(n,w(oh),w(os),w(owk),lx,iov,iinv,nevmx,eevmx,
     .                  nev,w(oz),w(oe))
            t2s = cpusec()
            call mkhsz(lov,n,w(oh),w(os))
            call ztoyy(w(oh),n,n,n,n,1,0)
            if (lov) then
              call ztoyy(w(os),n,n,n,n,1,0)
            endif
            lx = T
            t1x = cpusec()
            call diagno(n,w(oh),w(os),w(owk),lx,iov,iinv,nevmx,eevmx,
     .                  nev,w(oz),w(oe))
            t2x = cpusec()
            call mkhsz(lov,n,w(oh),w(os))
            t1l = cpusec()
            if (lov) then
              nbl = ILAENV(1,'ZHETRD','U',n,-1,-1,-1)
              lwork = n*(nbl+1)
              call defcc(ozwk, lwork)
              call defrr(orwk, 7*n)
              Lname = 'ZHEGVX'
              call ZHEGVX(1,'V','I','U',n,w(oh),n,w(os),n,0d0,0d0,1,
     .                    nevmx,abstol,nev,w(oe),w(oz),n,w(ozwk),
     .                    lwork,w(orwk),w(oiwk),w(oifail),info)
            else
              nbl = ILAENV(1,'ZHETRD','U',n,-1,-1,-1)
              nbl = max( nbl, ILAENV(1,'ZUNMTR','U',n,-1,-1,-1) )
              lwork = n*(nbl+1)
              call defcc(ozwk, lwork)
              call defrr(orwk, 7*n)
              Lname = 'ZHEEVX'
              call ZHEEVX('V','I','U',n,w(oh),n,0d0,0d0,1,nevmx,
     .                    abstol,nev,w(oe),w(oz),n,w(ozwk),lwork,
     .                    w(orwk),w(oiwk),w(oifail),info)
            endif
            t2l = cpusec()
          else
            call defrr(oz, n*n)
            call mkhsd(lov,n,w(oh),w(os))
            lx = F
            t1s = cpusec()
            call dsev1(n,w(oh),w(os),w(owk),0,lx,lov,iinv,nevmx,
     .                 eevmx,nev,w(oz),w(oe))
            t2s = cpusec()
            call mkhsd(lov,n,w(oh),w(os))
            lx = T
            t1x = cpusec()
            call dsev1(n,w(oh),w(os),w(owk),0,lx,lov,iinv,nevmx,
     .                 eevmx,nev,w(oz),w(oe))
            t2x = cpusec()
            call mkhsd(lov,n,w(oh),w(os))
            t1l = cpusec()
            if (lov) then
              nbl = ILAENV(1,'DSYTRD','U',n,-1,-1,-1)
              lwork = max( n*(nbl+3) , 8*n )
              call defrr(orwk, lwork)
              Lname = 'DSYGVX'
              call DSYGVX(1,'V','I','U',n,w(oh),n,w(os),n,0d0,0d0,
     .                    1,nevmx,abstol,nev,w(oe),w(oz),n,w(orwk),
     .                    lwork,w(oiwk),w(oifail),info)
            else
              nbl = ILAENV(1,'DSYTRD','U',n,-1,-1,-1)
              nbl = max( nbl, ILAENV(1,'DORMTR','U',n,-1,-1,-1) )
              lwork = max( n*(nbl+3) , 3*n )
              call defrr(orwk, lwork)
              Lname = 'DSYEVX'
              call DSYEVX('V','I','U',n,w(oh),n,0d0,0d0,1,nevmx,
     .                    abstol,nev,w(oe),w(oz),n,w(orwk),lwork,
     .                    w(oiwk),w(oifail),info)
            endif
            t2l = cpusec()
          endif
          call awrit2(' '//LNAME//': lwork = %i, info = %i',' ',
     .                128,i1mach(2),lwork,info)
        endif
      endif
      call rlse(oe)
      if (twice) then
        twice = F
        goto 4
      endif
C --- write timing results ---
      if (numprocs > 1) then
      if (procid == master) then
        if (.not. linv .and. serial) then
          call awrit4(
     .      ' Parallel Wall time: %;3ds.'//
     .      ' eispack CPU times: generic %;3ds, RISC %;3ds;'//
     .      ' lapack CPU time: %;3ds.%N',
     .      ' ',256,i1mach(2),t2p-t1p,t2s-t1s,t2x-t1x,t2l-t1l)
        elseif (linv .and. serial) then
          call awrit2(
     .      ' Parallel Wall time: %;3ds.'//
     .      ' lapack CPU time: %;3ds.%N',
     .      ' ',256,i1mach(2),t2p-t1p,t2s-t1s)
        else
          call awrit1(' Parallel Wall time: %;3ds.%N',
     .      ' ',256,i1mach(2),t2p-t1p)
        endif
      endif
      else
      if (linv) then
        call awrit1(' lapack CPU time: %;3ds.%N',' ',256,i1mach(2),t2s-t1s)
      else
        call awrit3(
     .      ' eispack CPU times: generic %;3ds, RISC %;3ds;'//
     .      ' lapack CPU time: %;3ds.%N',
     .      ' ',256,i1mach(2),t2s-t1s,t2x-t1x,t2l-t1l)
      endif
      endif
    3 continue
C#ifdefC MPE
C      ierr = MPE_FINISH_LOG('tpzhev')
C#endif
      if (numprocs > 1) then
C#ifdefC SCALI
C      call MPI_FINALIZE(ierr)
C      if ( procid == master ) then
C        call fexit(0,111,
C     .   'EXIT tpzhev on '//shortname(procid)(1:namelen(procid)),0)
C      endif
C#elseC
C      if ( procid == master ) then
C        call fexit(0,1,
C     .   ' EXIT tpzhev on '//shortname(procid)(1:namelen(procid)),0)
C      else
C        call fexit(0,0,' ',0)
C      endif
C#endifC
        continue
        else
        call fexit(0,0,' ',0)
        endif
      end

      subroutine mkhsz(lov,n,h,s)
      implicit none
      logical lov
      integer n,i,j
      double complex h(n,n),s(n,n)
      do i = 1, n
        do j = 1, n
          if (lov) s(i,j) = dcmplx(n-abs(i-j))
          if (i == j) then
            h(i,j) = dcmplx(1d0/(dble(i+j)-1d0))+dcmplx(1d0)
          else
            h(i,j) = dcmplx(1d0/(dble(i+j)-1d0),dble(j-i))
          endif
        enddo
      enddo
      end

      subroutine mkhsd(lov,n,h,s)
      implicit none
      logical lov
      integer n,i,j
      double precision h(n,n),s(n,n)
      do i = 1, n
        do j = 1, n
          if (lov) s(i,j) = n-abs(i-j)
          if (i == j) then
            h(i,j) = 1d0/(dble(i+j)-1d0)
          else
            h(i,j) = 1d0/(dble(i+j)-1d0)
          endif
        enddo
      enddo
      end

      subroutine mkinv(n,h)
      implicit none
      integer n,i,j
      double precision h(n,n)
      real r,ran1
      call ran1in(1)
      do  i = 1, n
        do j = i, n
          r = 0.01 * ran1()
          if (i. eq. j) then
            h(i,j) = 1d0 + r
          else
            h(i,j) = 0.5d0 + r
            h(j,i) = h(i,j)
          endif
        enddo
      enddo
      end

      subroutine checkz(n,c,zz,e)
      implicit none
      integer n,i,j,i1mach,iprint
      complex*16 c(n,n),zz(n,n)
      double precision dsqrt,DLAMCH,e(1),tol,diff,sum,rms

      tol = 4*n*DLAMCH('P')
      sum = 0d0
      do  i = 1, n
        do  j = 1, n
          diff = c(i,j) - e(j)*zz(i,j)
          if (abs(diff) > tol) then
            if (iprint() > 60) then
              call awrit3('H Z  -  E O Z != 0: i=%i, j=%i,'//
     .        ' c(i,j) - e(j)*zz(i,j)=%d',' ',256,i1mach(2),i,j,
     .        c(i,j) - e(j)*zz(i,j))
            endif
          endif
          sum = sum + diff**2
        enddo
      enddo
      rms = dsqrt(sum/n**2)
      if (rms > tol) then
        call awrit2(' Test failed. rms error: %;1g, tol=%;1g',' ',
     .              96,i1mach(2),rms,tol)
      else
        call awrit2(' Test passed. rms error: %;1g, tol=%;1g',' ',
     .              96,i1mach(2),rms,tol)
      endif
      end

      subroutine checkd(n,c,zz,e)
      implicit none
      integer n,i,j,i1mach,iprint
      double precision dsqrt,DLAMCH,c(n,n),zz(n,n)
      double precision e(1),tol,diff,sum,rms

      tol = n*DLAMCH('P')
      sum = 0d0
      do  i = 1, n
        do  j = 1, n
          diff = c(i,j) - e(j)*zz(i,j)
          if (abs(diff) > tol) then
            if (iprint() > 60) then
              call awrit3('H Z  -  E O Z != 0: i=%i, j=%i,'//
     .        ' c(i,j) - e(j)*zz(i,j)=%d',' ',256,i1mach(2),i,j,
     .        c(i,j) - e(j)*zz(i,j))
            endif
          endif
          sum = sum + diff**2
        enddo
      enddo
      rms = dsqrt(sum/n**2)
      if (rms > tol) then
        call awrit2(' Test failed. rms error: %;1g, tol=%;1g',' ',
     .              96,i1mach(2),rms,tol)
      else
        call awrit2(' Test passed. rms error: %;1g, tol=%;1g',' ',
     .              96,i1mach(2),rms,tol)
      endif
      end

      subroutine checki(n,a)
      implicit none
      integer n,i,j,k,i1mach
      double precision a(n,n),tol,diff,sum,rms,DLAMCH,dsqrt

      tol = n*DLAMCH('P')
      sum = 0d0
      do  i = 1, n
        do  j = 1, n
          if (i == j) then
            diff = a(i,j) - 1d0
          else
            diff = a(i,j)
          endif
          sum = sum + diff**2
        enddo
      enddo
      rms = dsqrt(sum/n**2)
      if (rms > tol) then
        call awrit2(' Test failed. rms error: %;1g, tol=%;1g',' ',
     .              96,i1mach(2),rms,tol)
      else
        call awrit2(' Test passed. rms error: %;1g, tol=%;1g',' ',
     .              96,i1mach(2),rms,tol)
      endif
      end

      subroutine symm(n,a)
C Copy upper triangle into lower triangle
      implicit none
      integer n
      double precision a(n,n)
      integer i, j
      do  i = 2, n
        do j = 1, i-1
          a(i,j) = a(j,i)
        enddo
      enddo
      end
