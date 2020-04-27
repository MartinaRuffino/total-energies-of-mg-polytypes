      subroutine dstrbp(nloop,np,nblk,index)
C- Distribute loop counter over processes
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nloop: number of times the loop is repeated
Ci   np   : number of processors
Ci   nblk : if  1, each processor gets contiguous loop elements
Ci        : if  0, a blocking size is chosen
Ci        : if -1, then do double loop balancing (see remarks)
Co Outputs:
Co   nblk : vector blocking size
Co   index: list of pointers, whose meaning depends on nblk.  See Remarks
Cr Remarks
Cr   A loop is of the form
Cr       do  i = 1, nloop, nblk
Cr   Each process will do a contiguous number of passes. Index(n)
Cr   is the counter that process n starts with. That is process n does
Cr       do i = index(n), index(n+1)-1, nblk
Cr   nblk is chosen (in the range 6-16) to give the best distribution.
Cr   The work is divided as evenly as possible over the processes.
Cr   Double loop balancing is for the case where the serial code is,
Cr       do i = 1, n
Cr       do j = i, n
Cr   In that case index is a 2D array and the loop is transformed into
Cr       do iloop = 1, index(procid,0)
Cr       i = index(procid,iloop)
Cr       do j = i, n
Cr   This distributes the work over i so that process 0 does 1 and n
Cr   process 1 does 2 and n-1 process 2 does 3 and n-2 and so on.
Cr   (Thanks to Andrew Canning)
Cu Updates
Cu   Written by Karen Johnston
C ----------------------------------------------------------------------
      implicit none
C Passed
      integer :: nloop,np,nblk
      integer, dimension(0:np) :: index
C Local
      integer :: i,step
      integer, dimension(0:np-1) :: xnode,inode
      integer iprint, lgunit

      if (nblk == -1) then
        call pdstlb(nloop,np,index)
        return
      endif
      step = nblk
! Initialising arrays
      do i=0,np-1
        inode(i)=0
        xnode(i)=0
      end do
      if (step==1) then
        call single(nloop,np,xnode,inode)
      else if (step == 0) then
! Find optimum step nblk 6-16
        call optimise(nloop,np,step,xnode,inode,nblk)
      else
        call multiple(nloop,np,step,xnode,inode)
      end if
      do i=0,np-1
        index(i)=inode(i)+1
      end do
      index(np)=nloop+1
      if (iprint() < 41) return
! Output indices
      write (*,'(" DSTRBP: i   inode   xnode   block=",i5)') nblk
      do i=0,np-1
        write (*,'(3x,3i7)') i,inode(i),xnode(i)
      end do

      end subroutine dstrbp

      subroutine single(nloop,np,xnode,inode)

      implicit none
      integer :: i=0,j=0,nloop,np
      integer :: min=0,rem=0
      integer :: times=0,rem2=0,total=0
      integer, dimension(0:np-1) :: inode,xnode
      integer iprint,lgunit

C      if (iprint() > 40) write(lgunit(1),*) 'DSTRBP SINGLE:'
! Split nblock evenly with nodes
      times = nloop/np
      rem2 = mod(nloop,np)
      call info5(41,0,0,
     .  ' DSTRBP, single:  nloop = %i  np = %i  times = %i  rem2 = %i',
     .  nloop,np,times,rem2,0)
C      if (iprint() > 40)
C     .  write(lgunit(1),*) 'nloop=',nloop,'np=',np,'times=',times
C     .  ,'rem2=',rem2
! Even no. of kpts per node without remainder
      inode(0)=0
      do i=1,np-1
        inode(i)=inode(i-1)+times
      end do
! Spread over the remainder
      if (rem2 /= 0) then
        do i=1,rem2
          do j=0,np-1
            if (j>i) then
              inode(j) = inode(j)+1
            end if
          end do
        end do
      end if
! Number of blocks in each node
      total = 0
      do i=0,np-2
        xnode(i) = inode(i+1)-inode(i)
        total = total + xnode(i)
      end do
      xnode(np-1) = nloop - total
      do i=1,np-1
        inode(i) = inode(i-1) + xnode(i-1)
      end do

      end subroutine single



      subroutine optimise(nloop,np,step,xnode,inode,size)

      implicit none
      integer, parameter :: stepmin=6,stepmax=16
      integer :: i=0,j=0,nloop,np,step
      integer :: size,rem=0
      integer :: nblock=0,rem2=0,total=0
      integer :: times=0,rem3=0
      integer, dimension(0:np-1) :: inode,xnode
      integer, dimension(0:nloop-1) :: iblock,xblock
      integer iprint,lgunit

      if (iprint() > 40) write(lgunit(1),*) 'DSTRBP OPTIMISE:'
! Optimises the block size
      size = stepmin
      do i=stepmin+1,stepmax
        if ( abs(mod(nloop,i)) <= abs(mod(nloop,size))) then
          size = i
        end if
      end do
      rem = mod(nloop,size)
! Split nloop into blocks
      nblock = nloop/size
      rem2 = mod(nloop,size)
      if (iprint() > 40) then
        write(lgunit(1),*) 'size=',size,'nblock=',nblock,'rem2=',rem2
        if (nblock < np) then
          write(lgunit(1),*) ""
          write(lgunit(1),*) "*** WARNING: No. blocks < no. nodes ***"
          write(lgunit(1),*) ""
        end if
      endif
! Even no. of kpts per block without remainder
      iblock(0)=0
      do i=1,nblock-1
        iblock(i)=iblock(i-1)+size
      end do
! Spread over the remainder
      if (rem2 /= 0) then
        do i=0,rem2-1
          do j=0,nblock-1
            if (j>i) then
              iblock(j) = iblock(j)+1
            end if
          end do
        end do
      end if
! Number of nloop in each block
      do i=0,nblock-2
        xblock(i) = iblock(i+1)-iblock(i)
        total = total + xblock(i)
      end do
      xblock(nblock-1) = nloop - total
      do i=1,nblock-1
        iblock(i) = iblock(i-1) + xblock(i-1)
      end do
! Split blocks over nodes
      times = nblock/np
      rem3 = mod(nblock,np)
      if (iprint() > 40)
     .  write(lgunit(1),*) 'nblock=',nblock,'np=',np,'times=',times
     .  ,'rem3=',rem3

! Even no. of blocks per node without remainder
      inode(0)=iblock(0)
      do i=1,np-1
        inode(i)=iblock(i*times)
      end do
! Spread over the remainder
      if (rem3 /= 0) then
        do i=0,rem3-2
          do j=np-1-i,np-1
            inode(j) = inode(j)+size
          end do
        end do
      end if
      if (iprint() > 40) then
        do i=0,np-1
          write(lgunit(1),*) 'inode(',i,')=',inode(i)
        end do
      endif
      do i=0,np-2
        xnode(i)=inode(i+1)-inode(i)
      end do
      xnode(np-1) = nloop-inode(np-1)

      end subroutine optimise


      subroutine multiple(nloop,np,step,xnode,inode)

      implicit none
      integer :: i=0,j=0,nloop,np,step
      integer :: size=0,rem=0
      integer :: nblock=0,rem2=0,total=0
      integer :: times=0,rem3=0
      integer, dimension(0:np-1) :: inode,xnode
      integer, dimension(0:nloop-1) :: iblock,xblock
      integer iprint,lgunit

      if (iprint() > 40) write(lgunit(1),*) 'DSTRBP MULTIPLE:'
      size = step
      rem = mod(nloop,size)
! Split nloop into blocks
      nblock = nloop/size
      rem2 = mod(nloop,size)
      if (iprint() > 40) then
        write(lgunit(1),*) 'size=',size,'nblock=',nblock,'rem2=',rem2
        if (nblock < np) then
          write(lgunit(1),*) ""
          write(lgunit(1),*) "*** WARNING: No. blocks < no. nodes ***"
          write(lgunit(1),*) ""
        end if
      endif
! Even no. of kpts per block without remainder
      iblock(0)=0
      do i=1,nblock-1
        iblock(i)=iblock(i-1)+size
      end do
! Spread over the remainder
      if (rem2 /= 0) then
        do i=0,rem2-1
          do j=0,nblock-1
            if (j>i) then
              iblock(j) = iblock(j)+1
            end if
          end do
        end do
      end if
! Number of nloop in each block
      do i=0,nblock-2
        xblock(i) = iblock(i+1)-iblock(i)
        total = total + xblock(i)
      end do
      xblock(nblock-1) = nloop - total
      do i=1,nblock-1
        iblock(i) = iblock(i-1) + xblock(i-1)
      end do
! Split blocks over nodes
      times = nblock/np
      rem3 = mod(nblock,np)
      if (iprint() > 40)
     .  write(lgunit(1),*) 'nblock=',nblock,'np=',np,'times=',times ,'rem3=',rem3

! Even no. of blocks per node without remainder
      inode(0)=iblock(0)
      do i=1,np-1
        inode(i)=iblock(i*times)
      end do
! Spread over the remainder
      if (rem3 /= 0) then
        do i=0,rem3-2
          do j=np-1-i,np-1
            inode(j) = inode(j)+size
          end do
        end do
      end if
      do i=0,np-2
        xnode(i)=inode(i+1)-inode(i)
      end do
      xnode(np-1) = nloop-inode(np-1)

      end subroutine multiple

      subroutine pdstlb(nloop,np,index)
      implicit none
      integer nloop, np, index(0:np-1,0:*)
      integer n, p, m, i, j, sum1, sum2, iprint, lgunit
      n = nloop
      do i = 0, np-1
        index(i,0) = 0
      enddo
      p = 0
      do  i = 0, n/2-1
        j = index(p,0) + 1
        index(p,j)   = 1 + i
        index(p,j+1) = n - i
        index(p,0) = index(p,0) + 2
        p = mod(p+1,np)
      enddo
      if (mod(n,2) /= 0) then
        j = index(p,0) + 1
        index(p,j) = n/2+1
        index(p,0) = index(p,0) + 1
      endif
      sum1 = 0
      sum2 = 0
      do  i = 1, n
        sum1 = sum1 + i
      enddo
      do i = 0, np-1
        do j= 1, index(i,0)
          sum2 = sum2 + index(i,j)
        enddo
      enddo
      if (iprint() > 40) then
        write (*,1)
        do i = 0, np-1
          write (lgunit(1),2) i,(index(i,j),j=0,index(i,0))
        enddo
      endif
      if (sum1 /= sum2) then
        if (iprint() > 0)
     .    call awrit2('Bug in pdstlb: sum1=%i sum2=%i',' ',128,
     .                lgunit(1),sum1,sum2)
        call rx0(' ')
      endif
    1 format (/' PDSTLB:'/' proc   total    loop number')
    2 format (i5,3x,i5,3x,256i5)
      end subroutine pdstlb

C      subroutine fmain
C      implicit none
C      integer procid,iloop,nloop,np,index0(1000000),i,j
C      integer,allocatable :: index(:,:)
C      write(*, 10)
C   10 format ('Enter nloop,np: '$)
C      read (*,*) nloop, np
C
CCC     write (*,20) nloop, np
CC   20 format ('nloop=',i5,' np=',i5)
CC     call pdstlb(nloop,np,index0)
C
C      allocate (index(0:np-1,0:nloop-1))
C      call dstrbp(nloop,np,-1,index(0,0))
C      do procid = 0, np-1
C        print "('procid=',i3'  nloop',i4)",
C     .    procid,index(procid,0)
C      enddo
C
C      do procid = 0, np-1
C        print "('pid',i3,2x,100(i3))", procid,
C     .    (index(procid,iloop), iloop = 1, index(procid,0))
C      enddo
C
C      end

