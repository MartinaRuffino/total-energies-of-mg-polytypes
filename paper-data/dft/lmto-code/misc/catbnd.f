C  Combines sets of bands into a single file
      subroutine fmain
      implicit none
      integer nb,nbmx,iq1,nq1,nline,ifile
      double precision de,ef,alat,q1,q2,q3,eferm,ebot,etop
      parameter( nbmx=500 )
      double precision eb(nbmx)

      integer iarg,ifmt,iargc,ifi,i,n
      character*20 fname
      character*1 f2(20)
      equivalence (f2,fname)
      logical cmdstr,lsequ,a2bin
      common /zz/ eb

      ifile = 0
      ef = -9999
      iarg = 0
      goto 12
   10 stop 'usage: catbnd [-ef=#] file ...'

   12 iarg = iarg+1
      if (.not. cmdstr(iarg,fname)) goto 10
      if (lsequ(fname,'-ef=',4,' ',n)) then
        i = 0
        if (.not. a2bin(f2(5),ef,4,0,' ',i,-1)) goto 10
        goto 12
      endif

      open(9, file='bnds.dat', err=10)

C ---- read until files exhausted -----
      print 345
  345 format('   File',16x,
     .  '  nq    emin      emax     efermi     shift')
      nline = 0
      iarg = iarg-1
   15 iarg = iarg+1
        if (.not. cmdstr(iarg,fname)) goto 20
        if (fname == 'bnds.dat')
     .    call rx('input file cannot be named "bnds.dat"')
        ifi = 10
        open(ifi, file=fname, status='OLD', err=10)
        ifile = ifile+1
        read(ifi,100) nb,eferm,alat
  100   format(i5,2f10.5)
        if (nb == 0) then
          close(ifi)
          goto 15
        endif
        if (ef == -9999) ef = eferm
        de = ef - eferm
        if (nb > nbmx) call rx('nb > nbmx')
        if (ifile == 1) write(9,100) nb,ef,alat

C ---   For each line this file ---
   91   read(ifi,*) nq1
        write(9,400) nq1
        ebot = 9999
        etop =-9999
        do  3  iq1 = 1, nq1
          read(ifi,500) q1,q2,q3,(eb(i),i=1,nb)
  500     format(3f10.5/(10f8.4))
          call daxpy(nb,1d0,de,0,eb,1)
          ebot = dmin1(ebot,eb(1))
          etop = dmax1(etop,eb(nb))
          write(9,500) q1,q2,q3,(eb(i),i=1,nb)
    3   continue
        print 334, fname, nq1, ebot, etop, ef, de
  334   format(2x,a20,i5,4f10.5)
        close(ifi)
      goto 15

C --- Clean up ---
   20 continue
      write(9,400) 0
  400 format(i5)
      close(9)
      end
      subroutine dositp(dos,ndos,emin,emax,ndos2,emin2,emax2,d2)
C- Interpolate density (or number) of states to another mesh
C ----------------------------------------------------------------
Ci Inputs
Ci   d2,ndos2,emin2,emax2: starting number of states,
Ci     number of mesh points and energy range of dos
Ci   ndos,emin,emax: output number of states,
Ci     number of mesh points and energy range of dos
Co Outputs
Co   dos: interpolated dos
Cr Remarks
C ----------------------------------------------------------------
C Passed parameters
      integer ndos,ndos2
      double precision dos(1),d2(1),emin,emax,emin2,emax2
C Local parameters
      integer i,iold
      double precision xiold,enew,oldslo,sloold,slonew

c      if (emin < emin2) stop 'dositp: cannot extrapolate emin'
c      if (emax > emax2) stop 'dositp: cannot extrapolate emax'

      sloold = (emax2-emin2)/(ndos2-1)
      slonew = (emax -emin )/(ndos -1)
      do  10  i = 1, ndos
        enew = slonew * (i-1) + emin
        xiold = (enew-emin2+.0000001)/sloold + 1
        iold  = xiold
        if (iold <= 0) iold=1
        if (iold >= ndos2) iold=ndos2-1
        xiold = xiold-iold
        dos(i) = d2(iold) + xiold*(d2(iold+1)-d2(iold))
   10 continue
      end
      SUBROUTINE DOSIO(DOS,NEMX,NDMX,NDOS,NLD,EMIN,EMAX,NSPIN,
     .   EFERM,DEL,IFMT,IFILE)
C- I/o for DOS, MSM's format
C ----------------------------------------------------------------
Ci Inputs
Ci   ifile
Co Outputs
Cp External procedures
Cp   none
Cr Remarks
Cr   dos is assumed to be declared as dos(nemx,nemx), but
Cr   dos(1..ndos,1..nld) are input/output
Cr
Cr   ifmt=0 for old format, ifmt=1 for new format
C ----------------------------------------------------------------
C      implicit none
C Passed parameters
      integer ifile,ifmt,ndmx,ndos,nemx,nld,nspin
      double precision del,eferm,emax,emin
      double precision dos(nemx,1)
C Local parameters
      integer ie,ild,iedum

C ----- WRITE BRANCH -----------------
      IF(IFILE.LT.0) THEN
      WRITE(-IFILE,760) EMIN,EMAX,NDOS,NLD,NSPIN,EFERM,DEL,IFMT
      IF(IFMT.EQ.0) THEN
        DO 10 IE=1,NDOS
  10    WRITE(-IFILE,761) IE,(DOS(IE,ILD),ILD=1,NLD)
      ELSE
        DO 11 ILD=1,NLD
  11    WRITE(-IFILE,762) (DOS(IE,ILD),IE=1,NDOS)
      ENDIF
  761 FORMAT(I5,6F12.5/(5X,6F12.5))
  760 FORMAT(2F10.5,3I5,2F10.5,I5)
  762 FORMAT(5F14.6)
      ENDIF
C ----- READ BRANCH -----------------
      IF(IFILE.GT.0) THEN
      READ(IFILE,760) EMIN,EMAX,NDOS,NLD,NSPIN,EFERM,DEL,IFMT
      IF(NDOS.GT.NEMX) STOP '*** DOSIO-READ: NDOS.GT.NEMX'
      IF(NLD .GT.NDMX) STOP '*** DOSIO-READ: NLD.GT.NDMX'
      IF(IFMT.EQ.0) THEN
        DO 20 IE=1,NDOS
  20    READ(IFILE,761) IEDUM,(DOS(IE,ILD),ILD=1,NLD)
      ELSE
        DO 21 ILD=1,NLD
  21    READ(IFILE,762) (DOS(IE,ILD),IE=1,NDOS)
      ENDIF
      ENDIF
      RETURN
      END
