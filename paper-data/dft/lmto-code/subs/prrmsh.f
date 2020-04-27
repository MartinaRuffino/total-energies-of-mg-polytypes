      subroutine prrmsh(strn,rofi,s,ns,nr,nc)
C- Print YL expansion to file
C ----------------------------------------------------------------------
Ci Inputs
Ci   strn  :printout string
Ci   rofi  :radial mesh points
Ci   s     :YL expansion to print
Ci   ns    :leading dimension of s
Ci   nr    :number of radial mesh points
Ci   nc    :number of columns
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cr   The following integrates the L=0 functions
Cr      mc out -coll 1,2 -av:nr,1 rmax -int 0 rmax
Cr   The following integrates r^2* L=1 function
Cr      mc out -e2 x1 'x2*x1*x1' -av:nr,1 rmax -int 0 rmax
Cu Updates
Cu   11 Nov 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nc,ns,ifi
      double precision rofi(nr),s(ns,nc)
C ... Local parameters
      character*(20) fmt, strn*(*), outs*120
      integer i,j,fopna,i1mach,i1
C     fmt = '(f10.6,10f15.10)'
C     Printout in 1 line, but only good for nc<=16
      fmt = '(f11.8,1p,16g18.10)'
C     Otherwise printout r in 1 line, followed by s
      if (nc>16) fmt = '(1p,(10e18.10))'
      ifi = fopna('out',-1,0)
      i1 = 1
      if (s(1,1) == 0) i1 = 2
      call awrit2('%% rows %i cols %i',' ',80,ifi,nr-i1+1,nc+1)
      do  i = i1, nr
        if (nc > 16) then
          write(ifi,'(f11.8)') rofi(i)
          write(ifi,fmt) (s(i,j), j=1,nc)
        else
          write(ifi,fmt) rofi(i), (s(i,j), j=1,nc)
        endif
      enddo
      call fclose(ifi)
      outs = 'prrmsh: done writing to file out data '//strn
      call awrit0('%a.  Continue?',outs,-len(outs),-i1mach(2))
      read(*,'(a80)') outs
      if (outs == 'q') call rx0('quit in prrmsh')
      end
      subroutine rdrmsh(ifi,rofi,s,ns,nr,nlm,a)
C- Read YL expansion from file
C ----------------------------------------------------------------------
Ci Inputs
Ci   rofi  :radial mesh points
Ci   ns    :leading dimension of s
Cio Inputs/Outputs
Cio  nlm   :An input zero value acts as f flag for rdrmsh to find nr, nlm and exit
Cio        :On output, nlm = number of columns in s
Co Outputs
Co   nr    :number of radial mesh points
Co   s     :s on radial mesh, if input nlm is nonzero
Co   a     :a mesh parameter, if input nlm is nonzero
Cr Remarks
Cu Updates
Cu   04 Apr 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ns,nr,nlm,ifi
      double precision rofi(nr),s(ns,nlm),a
C ... Dynamically allocated arrays
      real(8), allocatable :: dat(:,:)
C ... Local parameters
      integer i,j,k,nc,rdops
      real(8) xx,aold,gsign,tola
      procedure(integer) :: fopng,rdm
      procedure(real(8)) :: dmach

      rdops = 0; nr = 0; nc = 0; xx = 0
      rewind ifi
      i = rdm(ifi,rdops,0,' ',xx,nr,nc)
      if (nlm == 0) then
        nlm = nc-1
        return
      endif

      if (nlm < nc-1) call rx('rdrmsh: given array dimensioned too small to read')
      nlm = nc-1
      allocate(dat(nr,nc))
      rewind ifi
      i = rdm(ifi,rdops,nr*nc,' ',dat,nr,nc)
      k = 0
      if (dat(1,i) > 0) then
        nr = nr+1; k = 1; rofi(1) = 0; s(1,:) = 0
      endif
      if (ns<nr) call rx('rdrmsh: given array dimensioned too small to read')
      do  i = 1+k, nr
        rofi(i) = dat(i-k,1)
        forall (j = 2:nc) s(i,j-1) = dat(i-k,j)
      enddo

      a = log(rofi(nr)/rofi(nr-1))
      tola = 8*dmach(1)
      gsign = 1
      do   i = 1, 100
        aold = a
        a = log(rofi(nr)/rofi(nr-1)/(1-exp(-(nr-1)*a))*(1-exp(-(nr-2)*a)))
C       if (i > 95) write(stdo,'(i4,1p,2e15.5)') i,a,a-aold
        if (abs(a-aold) <= tola) exit
        if (i < 100) cycle
        call rx('rdrmsh failed to determine ''a'' parameter')
      enddo

      end
