      subroutine zwrmdmft(opt,filel,fnam,fmt,s,ns,nr,nc)
C- Write to disk or stdout a complex matrix, Haule's format with one line per row
C ----------------------------------------------------------------
Ci Inputs:
Ci   opt   :1s digit:
Ci           0 writes in ascii mode
Ci           1 writes in binary mode
Ci         10s digit: extension
Ci           0 Append extension
Ci           1 fnam is full file name
Ci        100s digit: file pointer.  File opening through fopnx
Ci           0 Do nothing after opening through fopnx
Ci             If file is freshly opened, position is rewound
Ci             If file is already open, position is left unchanged.
Ci           1 rewind file after opening through fopnx
Ci             File pointer is always rewound
Ci           2 Put file pointer to end-of-file if opened through fopnx
Ci  filel  : Label to be put on the first line of the file
Ci   fnam  : file name
Ci           Pass a blank string => write to stdout
Ci    fmt  : Fortran format to write ascii string, e.g. '(5f12.6)'
Ci           You can use awrite format, e.g. '(%5,6g)'
Ci           Pass a blank string => program uses internal default
Ci     s   : matrix to be printed out
Ci    ns   : leading dimension of s
Ci    nr   : Number of rows to write
Ci    nc   : Number of columns to write
Cr Remarks
Cr   Binary write first record is
Cr     nr  nc  cast  optional-string-length
Cr   If a optional-string-length > 0, second record contains string
Cr   Next record is entire array, written as nr*nc elements
Cu Updates
Cu   19 Oct 14 Adapted from ywrm
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,ns,nr,nc
      complex(8), intent(in) :: s(ns,nc)
      character(len=*), intent(in) :: fmt, filel, fnam
C ... Local parameters
      integer :: ifi,i,j,fopnx,stdo,nglob,ista,mode
      character(len=32) :: lfmt

      lfmt = 'f12.8'; if (fmt /= ' ') lfmt = fmt
      stdo = nglob('stdo')

      if (fnam /= ' ') then
        ista = 0; if (mod(opt,10) == 1) ista = ista+4
        if (mod(opt/100,10) == 2) ista = ista+32
        mode = 0; if (mod(opt/10,10) == 1) mode = 2
        ifi = fopnx(trim(fnam),mode,ista,-1)
      else
        ifi = stdo
        if (mod(opt,10) == 1) call rx('cannot write binary file to stdout')
      endif

      if (mod(opt,10) == 1) then
        call ywrm(1,filel,3,ifi,fmt,s,1,ns,nr,nc)
      else
        if (filel /= ' ') write(ifi,'(a)') trim(filel)
        do  i = 1, nr
          do j = 1, nc
            write(ifi,'(2x,2(1x,'//trim(lfmt)//'))', advance='no') s(i,j)
          end do
          write(ifi,'("")')
        end do
      endif

      if (ifi /= stdo) call fclose(ifi)

      end

