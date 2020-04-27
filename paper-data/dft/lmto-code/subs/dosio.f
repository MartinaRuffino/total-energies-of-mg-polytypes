      subroutine dosio(dos,nemx,ndmx,ndos,nld,emin,emax,nspin,
     .  eferm,del,ifmt,ifile)
C- I/O for dos, MSM's format
C ----------------------------------------------------------------------
Ci Inputs
Ci   nemx  :leading dimension of dos array (must be >ndos)
Ci   ndmx  :max number of dos allowed for input
Ci   ndos  :number of energy mesh points
Ci   ifile :>0 => read, <0 => write
Cio Inputs/Outputs
Cio  The following header is read from or written to disk
Cio  emin  :dos ranges from (emin,emax)
Cio  emax  :dos ranges from (emin,emax)
Cio  nld   :number dos per spin channel in interval (emin,emax)
Cio  nspin :number of spin channels
Cio  eferm :Fermi level
Cio  del   :Sampling delta
Cio  ifmt  :DOS format style
Cio  After the header, this is read from or written to disk
Cio  dos   :density of states
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   18 Mar 10
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifile,ifmt,ndmx,ndos,nemx,nld,nspin
      double precision del,eferm,emax,emin
      double precision dos(nemx,1)
C ... Local parameters
      integer ie,ild,iedum

C --- Write branch ---
      if (ifile < 0) then
        write (-ifile,2) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ifmt == 0) then
          do  ie = 1, ndos
            write (-ifile,1) ie,(dos(ie,ild),ild=1,nld*nspin)
          enddo
        else
          do  ild = 1, nld*nspin
            write (-ifile,3) (dos(ie,ild),ie=1,ndos)
          enddo
        endif
      endif
C --- Read branch ---
      if (ifile > 0) then
        read (ifile,2) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ndos > nemx) call rx('dosio: ndos > nemx')
        if (nld*nspin > ndmx) call rx('dosio: nld > ndmx')
        if (ifmt == 0) then
          do  ie = 1, ndos
            read (ifile,1) iedum,(dos(ie,ild),ild=1,nld*nspin)
          enddo
        elseif (ifmt == 1) then
          do  ild = 1, nld*nspin
            read (ifile,3) (dos(ie,ild),ie=1,ndos)
          enddo
        else
          call rx('dosio: bad fmt')
        endif
      endif
    1 format(i5,6F12.5/(5x,6F12.5))
    2 format(2F10.5,3I5,2F10.5,i5)
    3 format(5F14.6)
      end
