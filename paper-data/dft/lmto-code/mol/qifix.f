      subroutine qifix(alfa,n,si,six,ux)
C- Returns hamiltonian fix alfa*(qi^-qi)*inv(qi)*(qi^-qi) in ux.
C  At end, six contains qi^-qi and si is destroyed.
      implicit real*8 (a-h,p-z), integer (o)
      dimension si(n,n),six(n,n),ux(n,n)
      real w(1)
      common /w/ w

C For MPI ...
      integer, dimension(:), allocatable :: bproc
      double precision, dimension(:,:), allocatable :: sb
      integer mpipid,procid,master,numprocs,length,ierr
      integer i1mach,iprint
      logical mlog,cmdopt,MPI
      character outs*80, shortname*10, datim*26
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      call mpiprn(shortname,length)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      call tcn('qifix')
      if (iprint() >= 30) write(6,220) alfa
  220 format(' qifix:  alfa=',f10.4)
c ----- overwrite six with qi^-qi, make inv(qi)*(qi^-qi) in si ---
      call dpcopy(si,ux,1,n*n,1d0)
      call dpadd(six,si,1,n*n,-1d0)
      call dpcopy(six,si,1,n*n,1d0)
      call defrr(opiv,  n)
      call dsifa(ux,n,n,w(opiv),info)
      if (info /= 0) call rx('qifix:  cannot factor qi')
      if (MPI) then
        allocate (sb(n,n), stat=ierr)
        call dcopy(n*n,0d0,0,sb,1)
        allocate (bproc(0:numprocs), stat=ierr)
        call dstrbp(n,numprocs,1,bproc(0))
        n1 = bproc(procid)
        n2 = bproc(procid+1)-1
      else
        n1 = 1
        n2 = n
      endif
      do   j = n1, n2
        if (MPI .and. mlog) then
        if (j == bproc(procid)) then
          call gettime(datim)
          call awrit4(' qifix '//datim//' Process %i of %i on '
     .      //shortname(1:length)//
     .      ' starting columns %i to %i',' ',256,lgunit(3),
     .      procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        endif
        call dsisl(ux,n,n,w(opiv),si(1,j))
        if (MPI) then
          call dcopy(n,si(1,j),1,sb(1,j),1)
        endif
      enddo
      call rlse(opiv)
      ierr = mpipid(2)
      if (MPI) then
        call mpibc2(sb,n*n,4,3,mlog,'qifix','sb')
        call dcopy(n*n,sb,1,si,1)
        deallocate(bproc)
        deallocate(sb)
      endif
      call dgemm('N','N',n,n,n,1d0,six,n,si,n,0d0,ux,n)
      call dscal(n*n,alfa,ux,1)
      call tcx('qifix')
      end
