C A collection of subroutines common to tetwtz and tetwtq
      subroutine midk(kk,ee,xx,i,j,kout,xout)
C- Evaluate xx,kk on the line kk(i)-kk(j), at the Fermi energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   kk    :k-points: only points i and j are used
Ci   ee    :vector of energies: only ee(i) and ee(j) are used
Ci   xx    :vector of some some function: only xx(i) and xx(j) are used
Ci   i     :index to point i
Ci   j     :index to point j
Co Outputs
Co   kout  :k at which band crosses Fermi level
Co   xout  :interpolated value of xx at Fermi level
Cl Local variables
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer(4):: i,j
      real(8) ::  kk(3,1:4),xx(1:4),ee(1:4),kout(3),xout
C ... Local parameters
      real(8) ::  ratio

      ratio     = ee(i)/(ee(i)-ee(j))
      xout      = xx(i)     + ratio * (xx(j)-xx(i))
      kout(1:3) = kk(1:3,i) + ratio * (kk(1:3,j)-kk(1:3,i))
      end

      subroutine midk3(kk,ee,xx,yy,i,j,kout,xout,yout)
C- Evaluate xx,yy,kk on the line kk(i)-kk(j), at the Fermi energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   kk    :k-points: only points i and j are used
Ci   ee    :vector of energies: only ee(i) and ee(j) are used
Ci   xx    :vector of some some function: only xx(i) and xx(j) are used
Ci   yy    :vector of some some function: only yy(i) and yy(j) are used
Ci   i     :index to point i
Ci   j     :index to point j
Co Outputs
Co   kout  :k at which band crosses Fermi level
Co   xout  :interpolated value of xx at Fermi level
Co   yout  :interpolated value of yy at Fermi level
Cl Local variables
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer(4):: i,j
      real(8) ::  kk(3,1:4),xx(1:4),yy(1:4),ee(1:4),kout(3),xout,yout
C ... Local parameters
      real(8) ::  ratio

      ratio     = ee(i)/(ee(i)-ee(j))
      xout      = xx(i)     + ratio * (xx(j)-xx(i))
      yout      = yy(i)     + ratio * (yy(j)-yy(i))
      kout(1:3) = kk(1:3,i) + ratio * (kk(1:3,j)-kk(1:3,i))
      end

      subroutine addsciss(delta, ef, nnn, eig)
C- Scissors operator
C ----------------------------------------------------------------------
Ci Inputs
Ci   delta :shift
Ci   ef    :Fermi level
Ci   nnn   :number of eigenvalues to shift
Co Outputs
Co   eig   :eigenvalues eig(i)>ef are shifted by delta
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nnn
      real(8):: eig(nnn),ef,delta
C ... Local parameters
      integer i

C     write(6,*)' asssciss delta=', delta
      do i=1,nnn
        if (eig(i)>ef) eig(i)= eig(i)+delta
      enddo
      end

      subroutine chkdgn(ene,ndat,nrank,ixini,ixend,iof,ipr)
C- Group energy levels into families of degenerate levels
C ----------------------------------------------------------------------
Ci Inputs
Ci   ene   :vector of energies, sorted
Ci   ndat  :number of energies
Ci   iof   :offset, which is added to ixini and ixend
Ci   ipr   :If T, print check data
Co Outputs
Co   ixini :indices to levels with energy distinct from predecessor
Co   ixend :indices to levels with energy same as predecessor
Co         :ixend(i)=ixini(i) => no degeneracy for that block
Co   nrank :number of distinct energies
Cl Local variables
Cl  epsx   :Two energies are considered equal when they differ by less than epsx
Cr Remarks
Cr   This routine returns indices to starting and ending points
Cr   in ene that have a distinct energy
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical ipr
      integer ndat,nrank,ixini(ndat),ixend(ndat),iof
      real(8) ene(ndat)
C ... Local parameters
C     integer stdo,nglob
      integer(4) :: i,ix
      real(8),parameter:: epsx=1d-4

C     stdo = nglob('stdo')
C     Handle special cases 0 or 1 point
      if (ipr) print *, 'chgdgn: ndat=',ndat
      if (ndat<1) then
        nrank = 0
        return
      endif
      ixini(1) = 1
      if (ndat==1) then
        ixend(1) = 1
        nrank = 1
        return
      endif

      i = 1
      do ix = 2, ndat
        if ( abs(ene(ix)-ene(ix-1)) >epsx ) then
          ixend(i) = ix-1
          i = i + 1
          ixini(i) = ixend(i-1)+1
          if (ix==ndat) then
            ixend(i) = ix
          endif
        elseif (ix==ndat) then
          ixend(i) = ndat
        endif
      enddo
c
      nrank = i
*poption noparallel
      do  i = 1, nrank
        ixini(i) = ixini(i)+iof
        ixend(i) = ixend(i)+iof
      enddo
c-check write
      if (ipr) then
        print *,' nrank=',nrank
        do  i = 1, ndat
          print "(' i ',i3,' ene=',d15.7)", i,ene(i)
        enddo
        print *
        do  i = 1, nrank
          print "(' i ',2i3,' e=',d15.7)",
     .      ixini(i),ixend(i),ene(ixini(i)-iof)
        enddo
      endif
      end

