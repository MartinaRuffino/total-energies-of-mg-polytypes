      subroutine info(jpr,l1,l2,string,a1,a2)
C- Printout when ipr>=jpr, with some blank lines
C ----------------------------------------------------------------------
Ci Inputs
Ci   jpr   :print with verbosity is at least abs(jpr)
Ci         :jpr<0 -> print only if master node (MPI)
Ci   l1    :number of blank lines to print out before string
Ci   l2    :number of blank lines to print out after string
Ci         :l2=-1 => write without carriage return
Ci         :l2=-2 => print nothing except l1 blank lines (if any)
Ci   string:string to print
Co Outputs
Co   string is printed on stdo
Cr Remarks
Cu Updates
Cu   12 Apr 19 call awriteb to turn off buffer mode before writing to stdo
Cu   30 Jan 13 New l2=-2
Cu   16 Apr 12 Uses %$ construct in awrite to suppress newline (l2<0)
Cu   24 Nov 02 New info2, info5, info8
Cu   26 Jun 02 Add arguments a1,a2 so info can write data
Cu   02 Nov 01 Use awrite instead of write statement
Cu             Special handling of l2=-1
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer jpr,l1,l2,a1,a2,a3,a4,a5,a6,a7,a8
      character*(*) string
C ... Local parameters
      integer i,stdo,recl,lens,nargs,ibuf
      procedure(integer) :: iprint,lgunit,mpipid,awriteb
      parameter (recl=1024)  ! should be same as lens in awrite.f
      character*(recl) lstr

      nargs = 2
      goto 990

      entry info0(jpr,l1,l2,string)
      nargs = 0
      goto 990

      entry info2(jpr,l1,l2,string,a1,a2)
      nargs = 2
      goto 990

      entry info5(jpr,l1,l2,string,a1,a2,a3,a4,a5)
      nargs=5
      goto 990

      entry info8(jpr,l1,l2,string,a1,a2,a3,a4,a5,a6,a7,a8)
      nargs = 8

  990 continue ! All entry points meet here, with nargs defined
      if (iprint() < iabs(jpr)) return

      if (jpr < 0) then
        i = mpipid(1)
        if (i /= 0) return
      endif
      stdo = lgunit(1)
      do  i = 1, l1
        write(stdo,100)
      enddo
      if (l2 == -2) return

      ibuf = awriteb(0)
C     ibuf = 0  ! debugging to catch points where info should not be used
      if (ibuf /= 0) then
        i = awriteb(1)
C       print *, 'line 67',i
      endif

      lens = len(string)
      if (l2 >= 0) then

         if(nargs==0) call awrit0(string,' ', recl, stdo)
         if(nargs==2) call awrit2(string,' ', recl, stdo,a1,a2)
         if(nargs==5) call awrit5(string,' ', recl, stdo,a1,a2,a3,a4,a5)
         if(nargs==8) call awrit8(string,' ', recl, stdo,a1,a2,a3,a4,a5,a6,a7,a8)

CML        call awrit8(string,' ',recl,stdo,a1,a2,a3,a4,a5,a6,a7,a8)

C       call awrit2(string,' ',recl,stdo,a1,a2)
      else
        lstr = '%$'//string

         if(nargs==0) call awrit0(lstr(1:2+lens),' ', recl, stdo)
         if(nargs==2) call awrit2(lstr(1:2+lens),' ', recl, stdo,a1,a2)
         if(nargs==5) call awrit5(lstr(1:2+lens),' ', recl, stdo,a1,a2,a3,a4,a5)
         if(nargs==8) call awrit8(lstr(1:2+lens),' ', recl, stdo,a1,a2,a3,a4,a5,a6,a7,a8)


C        call awrit8(lstr(1:2+lens),' ',
C     .    recl,stdo,a1,a2,a3,a4,a5,a6,a7,a8)
C       i = awrite(string,lstr,recl,0,a1,a2,i,i,i,i,i,i)
C       i = awrite(string,lstr,recl,0,a1,a2,a3,a4,a5,a6,a7,a8)
C        call cwrite(lstr,0,i-1,0)
      endif
      do  i = 1, l2
        write(stdo,100)
      enddo
  100 format(1x)

      if (ibuf /= 0) then
        i = awriteb(-1)
C       print *, 'line 103',i
      endif

      return

C      entry info8(jpr,l1,l2,string,a1,a2,a3,a4,a5,a6,a7,a8)
CC- Same as info, but will handle up to 8 arguments.
C      stdo = lgunit(1)
C
C      do  110  i = 1, l1
C  110 write(stdo,100)
C      if (l2 >= 0) then
C        call awrit8(string,' ',recl,stdo,a1,a2,a3,a4,a5,a6,a7,a8)
C      else
C        i = awrite(string,lstr,recl,0,a1,a2,a3,a4,a5,a6,a7,a8)
C        call cwrite(lstr,0,i-1,0)
C      endif
C      do  120  i = 1, l2
C  120 write(stdo,100)

      end
C      subroutine fmain
C
C      call info2(20,0,0,'ab%9ftest indent and integer %i to %i',-1,1111)
C      call info2(20,0,-1,'ab%9ftest%N%11ftest no newline ..',-1,1111)
C      call info2(20,0,0,'. done',-1,1111)
C
C      end
