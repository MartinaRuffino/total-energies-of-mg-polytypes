C Convolves a column of points
      subroutine fmain
      implicit none
      character*80 first,fmt
      character*1 f2(20)
      integer iarg,n,nr,nc,ifi,i,j,k,rdops,mode
      real(8), allocatable :: sin(:,:),sout(:,:)
      double precision gc,e0,ef,ep,sigma,xx
      procedure(logical) :: cmdstr,lsequ,a2bin
      procedure(integer) :: fopng,nargf,garg,rdm,i1mach
      equivalence (f2,first)

      nr = 0; nc = 0; mode = 2
      rdops = 0
      gc = .01d0
      ef = 0
      e0 = 0
      ep = 0.5d0*0
      sigma = 0

      fmt = '(5f12.6)'
      iarg = 1
   15 continue
      if (.not. cmdstr(iarg,first)) goto 20
      if (.not. lsequ(first,'-',1,' ',n)) goto 30

      if (lsequ(first,'-f',2,' ',n)) then
        fmt = '('
        call strcat(fmt,1,' ',f2(3),99,' ',n)
        call strcat(fmt,99,' ',')',1,' ',n)
      elseif (garg('-pr',iarg,2,' ',1,1,i,j,k) /= 0) then
        call pshpr(k)
      else if (lsequ(first,'-gc',3,' ',n)) then
        i = 4
        if (.not. a2bin(first,gc,4,0,' ',i,-1)) goto 20
      else if (lsequ(first,'-ef',3,' ',n)) then
        i = 4
        if (.not. a2bin(first,ef,4,0,' ',i,-1)) goto 20
      else if (lsequ(first,'-e0',3,' ',n)) then
        i = 4
        if (.not. a2bin(first,e0,4,0,' ',i,-1)) goto 20
      else if (lsequ(first,'-ep',3,' ',n)) then
        i = 4
        if (.not. a2bin(first,ep,4,0,' ',i,-1)) goto 20
      else if (lsequ(first,'-sig=',5,' ',n)) then
        i = 5
        if (.not. a2bin(first,sigma,4,0,' ',i,-1)) goto 20
      else if (lsequ(first,'-nc=',3,' ',n)) then
        i = 4
        if (.not. a2bin(first,nc,2,0,' ',i,-1)) goto 20
      elseif (lsequ(first,'- ',2,' ',n)) then
        iarg = iarg+1
        if (.not. cmdstr(iarg,first)) goto 20
        goto 30
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15
   30 continue

      if (first == '.' .or. nargf() == iarg) then
        ifi = i1mach(1)
      else
        ifi = fopng(first,-1,1)
      endif
      i = rdm(ifi,rdops,0,' ',[xx],nr,nc)

      if (mode == 1) then
        call rx('not ready')
      else
        call info5(10,0,0,' convolve : file '//trim(first)//' has %i rows, %i cols '//
     .    ' ... gaussian broadening with sigma=%d',nr,nc,sigma,0,0)
      endif
      allocate(sin(nr,nc),sout(nr,nc))
      rewind ifi
      i = rdm(ifi,rdops,nr*nc,' ',sin,nr,nc)
      call fclose(ifi)
      call dcopy(nr,sin,1,sout,1)
      do  i = 2, nc
        call convolve(2,gc,ef,e0,ep,sigma,nr,sin,sin(1,i),sout(1,i))
      enddo
      call snot

      iarg = iarg+1
      if (iarg >= nargf()) then
        ifi = i1mach(1)
      else
        if (.not. cmdstr(iarg,first)) goto 20
        ifi = fopng(first,-1,0)
      endif
      call prm(ifi,fmt,sout,nr,nr-1,nc)
      return

   20 print 331
  331 format('Usage: convolve [-switches] InFileName [OutFileName]'
     .   /7x,'-gc=Gamma_c :  (Lorenzian broadening) Core width'
     .   /7x,'-ef=E_F     :  (Lorenzian broadening) Fermi energy'
     .   /7x,'-e0=E_0     :  (Lorenzian broadening) Bottom of CB'
     .   /7x,'-ep=E_p     :  (Lorenzian broadening) Plasma energy'
     .   /7x,'-sig=sigma  :  (Gaussian broadening)  Instrumental Gaussian width')
      call fexit(0,0,' ',0)

      end
      subroutine prm(ifi,fmt,s,ns,nr,nc)
      integer nr,nc
      double precision s(ns,nc)
      character*(*) fmt
      integer i,j
      call awrit2('%% rows %i cols %i',' ',80,ifi,nr,nc)
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(i,j), j=1,nc)
      end
      subroutine snot
      end
