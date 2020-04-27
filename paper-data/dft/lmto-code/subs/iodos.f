      subroutine iodos(mode,ifi,dos,nemx,ndmx,ndos,nchan,emin,emax,nspin,eferm,ifmt)
C- I/O for density-of-states or related quantity
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do nothing
Ci         :1 read/write header only
Ci         :2 read/write dos only
Ci         :3 read/write both header and DOS
Ci         :10s digit, read only
Ci         :1 => requires file ndos,nchan,nspin to match passed values
Ci         :100s digit
Ci         :0 DOS I/O with 6 decimal places
Ci         :1 DOS I/O with 8 decimal places
Ci   ifi   :>0 => read, <0 => write
Ci   nemx  :leading dimension of dos array (must be >ndos)
Ci         :Not referenced unless mode>=2
Ci   ndmx  :max number of dos channels in dos array:
Ci         :(must be > nchan*nspin)
Ci         :Not referenced unless mode>=2
Ci   ndos  :number of energy mesh points
Cio Inputs/Outputs
Cio  The following header is read from or written to disk
Cio  emin  :dos ranges from (emin,emax)
Cio  emax  :dos ranges from (emin,emax)
Cio  nchan :number dos channels per spin
Cio  nspin :number of spin channels
Cio  eferm :Fermi level
Cio  del   :Sampling delta
Cio  ifmt  :DOS format style
Cio        :0 dos read/written in rdm-compatible format
Cio        :1 dos read/written in Methfessel format
Cio  After the header, this is read from or written to disk
Cio  dos   :density of states ... not referenced unless mode>=2
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   31 Jul 17 Header in MSM format now free format
Cu   24 Feb 15 New long option, 100s digit mode
Cu   18 Mar 10 Adapted from dosio
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ifi,ifmt,ndmx,ndos,nemx,nchan,nspin
      double precision del,eferm,emax,emin
      double precision dos(nemx,*)
C ... Local parameters
      integer ie,ild,iedum,mode0,mode1,mode2,ndosl,nchanl,nspinl,ifmtl
      double precision de,omg
C ... External calls
      external awrit3,isanrg,rxi

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)

      if (mode0 == 0) return

C --- Write branch ---
      if (ifi < 0) then
        if (mod(mode0,2) == 1) then
          if (ifmt == 0) then
          elseif (ifmt == 1) then
C           write (-ifi,4) emin,emax,ndos,nchan,nspin,eferm,0d0,ifmt
            call awrit8('%2;10,5D  %i  %i  %i %;10,5D  %;1,1d  %i',' ',120,-ifi,[emin,emax],ndos,nchan,nspin,eferm,0d0,ifmt,8)
          endif
        endif
        if (mode0 < 2) return
        if (ifmt == 0) then
          call awrit3('%% rows %i  cols %i  eferm=%d',' ',80,-ifi,ndos,nchan*nspin+1,eferm)
C          write (-ifi,1) ndos,nchan*nspin+1,eferm
C    1     format('% rows',i6,'  cols',i5,'  eferm=',f11.6)
          de = (emax-emin)/(ndos-1)
          do  ie = 1, ndos
            omg = emin + de*(ie-1)
            if (mode2 == 0) write (-ifi,2)  omg,(dos(ie,ild),ild=1,nchan*nspin)
            if (mode2 == 1) write (-ifi,12) omg,(dos(ie,ild),ild=1,nchan*nspin)
    2       format(f12.5,6F13.6:/(12x,6F13.6))
   12       format(f14.8,6F18.8:/(14x,6F18.8))
          enddo
        elseif (ifmt == 1) then
          do  ild = 1, nchan*nspin
            if (mode2 == 0) write (-ifi,5) (dos(ie,ild),ie=1,ndos)
            if (mode2 == 1) write (-ifi,15) (dos(ie,ild),ie=1,ndos)
          enddo
        else
          call rxi('iodos does not recognize ifmt =',ifmt)
        endif
      endif

C --- Read branch ---
      if (ifi > 0) then
        if (mod(mode0,2) == 1) then
        if (ifmt == 0) call rx('iodos: not ready for format 0')
        if (mode1 == 0) then
            read (ifi,*) emin,emax,ndos,nchan,nspin,eferm,del,ifmtl
        else
          read (ifi,*) emin,emax,ndosl,nchanl,nspinl,eferm,del,ifmtl
          call sanrg(.true.,ndosl,ndos,ndos,'IODOS','file ndos')
          call sanrg(.true.,nchanl,nchan,nchan,'IODOS','file nchan')
          call sanrg(.true.,nspinl,nspin,nspin,'IODOS','file nspin')
        endif
        endif
        if (mode0 < 2) return
C       Sanity checks
        if (ndos > nemx) call rx('iodos: ndos > nemx')
        if (nchan*nspin > ndmx) call rx('iodos: nchan > ndmx')
        if (ifmt == 0) then
          do  ie = 1, ndos
            read (ifi,3) iedum,(dos(ie,ild),ild=1,nchan*nspin)
    3       format(i5,6F12.5/(5x,6F12.5))
          enddo
        elseif (ifmt == 1) then
          do  ild = 1, nchan*nspin
            if (mode2 == 0) read(ifi,5) (dos(ie,ild),ie=1,ndos)
            if (mode2 == 1) read(ifi,15) (dos(ie,ild),ie=1,ndos)
          enddo
        else
          call rx('iodos: bad format')
        endif
      endif
    4 format(2F10.5,3I5,2F10.5,i5)
    5 format(5F14.6)
   15 format(5F18.8)
      end
