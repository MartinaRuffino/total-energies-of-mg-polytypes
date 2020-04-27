      logical function aiocor(nr,nsp,a,rmax,rhoc,sumec,sumtc,ifi)
C- File I/O for core charge density.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci   nr,nsp,a,rmax
Ci   rhoc, if file write
Co Outputs
Co   rhoc, if file read
Cr Remarks
Cr    Format for core in atomic file begins with
Cr    category 'CORE:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, aiocor returns the value true and rhoc if:
Cr       the category is found, and
Cr       the file's value of a and nr match input and rmax is close to
Cr         file's value, and
Cr       the density is read without error.
Cu Updates
Cu   26 Apr 03 Added MPI calls.  Does not broadcast file read
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,rhoc(nr,nsp),sumec,sumtc
C Local parameters
      integer i,isp,nr2,nsp2
      integer mpipid,procid
      double precision a2,rmax2
      logical scat

      procid = mpipid(1)

      if (ifi > 0) then
        aiocor = .false.
        if (procid == 0) then
        if (scat(ifi,'CORE:',':',.true.)) then
          read(ifi,102,err=15) nr2,nsp2,a2,rmax2,sumec,sumtc
          if (nr == 0) nr = nr2
          if (nsp == 0) nsp = nsp2
          if (a == 0) a = a2
          if (rmax == 0) rmax = rmax2
          if (a2 == a .and. nr == nr2 .and.
     .      dabs(rmax2-rmax) <= 1D-5) then
          do isp = 1, min0(nsp2, nsp)
            read(ifi,101) (rhoc(i,isp), i=1,nr)
          enddo
          if (nsp < nsp2) call dscal(nr,2d0,rhoc,1)
          if (nsp > nsp2) then
            call dscal(nr,.5d0,rhoc,1)
            call dcopy(nr,rhoc,1,rhoc(1,2),1)
          endif
          aiocor = .true.
          endif
        endif
        endif
      else
      if (procid == 0) then
        write(-ifi,'(''CORE:'')')
        write(-ifi,102) nr,nsp,a,rmax,sumec,sumtc
          do isp = 1, nsp
          write(-ifi,101) (rhoc(i,isp),i = 1,nr)
          enddo
      endif
      aiocor = .true.
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f12.5,2f18.8)
   15 continue
      end
