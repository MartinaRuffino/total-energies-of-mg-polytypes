      subroutine zdia22(opt,a,er,ec,z)
C- Diagonalizes a 2x2 complex matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0  => hermitian matrix, eigenvalues returned in ascending order
Ci         :1  => complex general matrix
Ci         :2  => hermitian matrix, eigenvalues returned in descending order
Ci         :4  => hermitian matrix, eigenvalues returned in whatever order
Ci                maximizes |z(1,1)|
Ci         :10 => a is a real symmetric matrix; evecs are returned as real
Ci   a     :matrix to invert
Co Outputs
Co   er    :Eigenvalues if opt=0 or opt=10
Co   ec    :Eigenvalues if opt not zero
Co    z    :Eigenvectors
Cr Remarks
Cu Updates
Cu   04 Nov 13 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt
      double precision er(2)
      double complex a(2,2),ec(2),z(2,2)
C ... Local parameters
      integer i,opt0
      double complex zx(2),epm,ratio,zl(2,2),al(2,2)
      double precision tol,xx,areal(2,2)
      parameter (tol=1d-16)

C ... Copy a to al so as not to overwrite a in case it is real
      opt0 = mod(opt,10)
      if (opt == 10) then
        call dcopy(4,a,1,areal,1)
        al(:,:) = areal(:,:)
      else
        al(:,:) = a(:,:)
      endif

C ... Case one nondiagonal matrix element vanishes
      if (cdabs(al(1,2)*al(2,1)) == 0) then
        zl(1,1) = 1; zl(2,2) = 1; zl(1,2) = 0; zl(2,1) = 0
        zx(1) = al(1,1)
        zx(2) = al(2,2)

C ... Eigenvalues in general case
      else
        epm = (al(1,1) - al(2,2))**2/4 + al(1,2)*al(2,1)
        if (mod(opt0,2) == 0) then
          epm = dsqrt(dble(epm))
        else
          epm = cdsqrt(epm)
        endif
        zx(1) = (al(1,1) + al(2,2))/2 - epm
        zx(2) = (al(1,1) + al(2,2))/2 + epm
      endif

C     Copy to er (if hermitian) else ec
      if (mod(opt0,2) == 0) then
        er = zx
      else
        ec = zx
      endif

C     Diagonal matrix
      if (cdabs(al(1,2))+cdabs(al(2,1)) == 0 .or.
     .    cdabs(zx(1)-zx(2)) == 0) goto 90

C     Eigenvectors in zl, so as not to overwrite z in case z is real
      do i = 1, 2
        if (cdabs(al(1,1)-zx(i)) > cdabs(al(2,2)-zx(i))) then
          ratio = al(1,2)/(al(1,1)-zx(i)) ! z1/z2
          xx = dsqrt(1+dble(ratio*dconjg(ratio))) !normalization
          zl(1,i) = -ratio/xx
          zl(2,i) = 1/xx
        else
          ratio = al(2,1)/(al(2,2)-zx(i)) ! z2/z1
          xx = dsqrt(1+dble(ratio*dconjg(ratio))) !normalization
          zl(1,i) = 1/xx
          zl(2,i) = -ratio/xx
        endif
      enddo

   90 continue

C     Reorder eigenvalues, eigenvectors
      if (opt0 /= 1) then
        if (er(1) > er(2) .and. opt0 == 0  .or.
     .      er(1) < er(2) .and. opt0 == 2 .or.
     .    cdabs(zl(1,2)) > cdabs(zl(1,1)) .and. opt0 == 4) then
          xx = er(1); er(1) = er(2); er(2) = xx
          zx = zl(:,1) ; zl(:,1) = zl(:,2) ; zl(:,2) = zx
        endif
      endif

C      if (cdabs(zl(1,1)-1) > 1d-10) then
C        print *, 'hi'
C        stop
C      endif

C     Return zl as a real matrix
      if (opt == 10) then
        call dcopy(4,zl,2,z,1)
C     Return zl as a complex matrix
      else
        call dcopy(8,zl,1,z,1)
      endif

      end
C     Test
C a.out | head -6 | tee a ; a.out | tail -12 | head -6 | tee e; a.out | tail -6 | tee z
C mc -f2f15.10   z -i a -x z -x  e -v2dia --
C for real symmetric matrix, write a, e, z this way
C a.out | head -3 | tee a ;  a.out | tail -7 | head -3 | tee e ; a.out | tail -3 | tee z
C      subroutine fmain
C
C      integer opt
C      double precision er(2), areal(2,2), zreal(2,2)
C      double complex a(2,2), z(2,2), ec(2)
C
C      opt = 1
CC Hermitian matrix with a12 = 0 (uncomment next line)
CC     opt = 0
CC real symmetrix matrix (uncomment next line)
CC     opt = 10
C
CC General matrix
C      a(1,1) = (1d0,2d0)
C      a(1,2) = (2d0,4d0)
C      a(2,1) = (4d0,2d0)
C      a(2,2) = (3d0,5d0)
C
CC General matrix with a12 = 0
C      a(1,1) = (1d0,2d0)
C      a(1,2) = (2d0,4d0)*0
C      a(2,1) = (4d0,2d0)
C      a(2,2) = (3d0,5d0)
C
C      if (opt == 0 .or. opt == 10) a = a + transpose(dconjg(a))
C
C
C      if (opt == 10) then
C        areal(1,1) = 3
C        areal(1,2) = 4
C        areal(2,1) = 4
C        areal(2,2) = 9
C
C        write(*,"('% rows 2 cols 2  symm')")
C        print 333, areal(1,1:2)
C        print 333, areal(2,1:2)
C
C        call zdia22(opt,areal,er,ec,zreal)
C
C        print *, 'evals, evecs :'
C        write(*,"('% rows 2 cols 1')")
C        write(*,332) er
C        print *
C
C        write(*,"('% rows 2 cols 2')")
C        print 333, zreal(1,1:2)
C        print 333, zreal(2,1:2)
C
C        return
C      endif
C
C
C
C      write(*,"('% complex rows 2 cols 2  a')")
C      print 333, dble(a(1,1:2))
C      print 333, dble(a(2,1:2))
C      print *
C      print 333, dimag(a(1,1:2))
C      print 333, dimag(a(2,1:2))
C  333 format(2f15.10)
C
C      call zdia22(opt,a,er,ec,z)
C
C      print *, 'evals:'
C      if (opt == 0) then
C      write(*,"('% rows 2 cols 1')")
C      write(*,332) er
C      print *
C      print *
C      print *
C      else
C      write(*,"('% complex rows 2 cols 1')")
C      write(*,332) dble(ec)
C      print *
C      write(*,332) dimag(ec)
C  332 format(1f15.10)
C      endif
C
C      write(*,"('% complex rows 2 cols 2')")
C      print 333, dble(z(1,1:2))
C      print 333, dble(z(2,1:2))
C      print *
C      print 333, dimag(z(1,1:2))
C      print 333, dimag(z(2,1:2))
C
C
C      end
