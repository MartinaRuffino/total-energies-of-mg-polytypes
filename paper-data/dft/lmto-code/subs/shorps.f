      subroutine shorps(nbas,plat,mode,pin,pout)
C- Shift each basis vector by multiples of plat according to mode
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  : size of basis
Ci   nbas  : nbas<0: sign used as a switch.
Ci         : shorps returns in pout the change in pin, in units of plat
Ci   plat  : primitive lattice vectors, in units of alat
Ci   mode  : vector of length 3 governing shifts along selected axes.
Ci         : 0 suppresses shifts along plat(j)
Ci         : 1 shifts to unit cell at origin (pos in 1st quadrant)
Ci         : 2 shifts to minimize length of pos
Ci         : 10s digit of mode(1) (1st digit mode=2 only)
Ci         : If nonzero, do not shorten lattice vector unless
Ci         : change in length of vector is > 10^(-i), where i is
Ci         : 10s digit.
Ci   pin:  position (basis) vectors
Co Outputs
Co   pout  (may point to the same address space as pin).
Cr Remarks
Cr   pos = f . plat, with f = integer + fraction for each plat.
Cr   Integer part according to mode.
Cu Updates
Cu   23 Feb 09 New 10s mode in order to suppress change in basis
Ci             vectors when they are not shortened.
Cu   09 Jan 09 (A. Lozovoi) Do not change initial values of basis
Cu             vectors if they fall on the cell boundary
Cu   10 Apr 02 Patch to handle mixed boundary conditions
C ----------------------------------------------------------------
      implicit none
      integer nbas,mode(3)
      double precision plat(3,3),pin(3,abs(nbas)),pout(3,abs(nbas))
C Local
      double precision qlat(3,3),xx,x0,x(3),a2,ap,p0(3),xmin(3),amin
      integer ib,i,m,j1,j2,j3,nbas1,j1max,j2max,j3max,lmod(3)
      double precision tol,toll,dpin,dpout
      parameter (tol = 1d-12)

C     call tcn('shorps')
      nbas1 = iabs(nbas)
      lmod(1) = mod(mode(1),10)
      lmod(2) = mode(2)
      lmod(3) = mode(3)
      if (mode(1) >= 10) then
        i = mod(mode(1)/10,10)
        toll = 10d0**(-i)
      endif
C     call prmx('starting pos',pin,3,3,nbas1)

C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
      call mkqlat(plat,qlat,xx)

      do  ib = 1, nbas1

        call dpcopy(pin(1,ib),p0,1,3,1d0)

C   --- Reduce to unit cell centered at or near origin ---
        do  i = 1, 3
C   ... x0 is projection of pin along plat(i)
        x0 = pin(1,ib)*qlat(1,i)+pin(2,ib)*qlat(2,i)+pin(3,ib)*qlat(3,i)
        if (lmod(i) <= 0) then
          x(i) = x0
        else
C   ... leave component x0 intact if -.5<x0<.5 (within the 1st cell)
          if (dabs(x0) <= 0.5d0+tol) then
            xx = 0d0
          else
            xx = idnint(x0)
          endif
C     ... first octant for lmod=1
          if (lmod(i) == 1 .and. x0-xx < -tol) xx = xx-1
          x(i) = x0-xx
        endif
        enddo
        forall (m = 1:3) pout(m,ib) = x(1)*plat(m,1)+x(2)*plat(m,2)+x(3)*plat(m,3)

C   --- Try shortening by adding +/- lattice vectors ---
        j1max = 1
        if (lmod(1) <= 1) j1max = 0
        j2max = 1
        if (lmod(2) <= 1) j2max = 0
        j3max = 1
        if (lmod(3) <= 1) j3max = 0
    5   continue
        amin = 0
        do  j1 = -j1max, j1max
        do  j2 = -j2max, j2max
        do  j3 = -j3max, j3max

C     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
          do  i = 1, 3
            x(i) = plat(i,1)*j1+plat(i,2)*j2+plat(i,3)*j3
          enddo
          a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          ap = pout(1,ib)*x(1) + pout(2,ib)*x(2) + pout(3,ib)*x(3)
          if (a2+2*ap < amin) then
            xmin(1) = x(1)
            xmin(2) = x(2)
            xmin(3) = x(3)
            amin = a2+2*ap
          endif
        enddo
        enddo
        enddo
c       if (amin < 0) then
        if (amin < -tol) then
          pout(1,ib) = pout(1,ib) + xmin(1)
          pout(2,ib) = pout(2,ib) + xmin(2)
          pout(3,ib) = pout(3,ib) + xmin(3)
C         In cases w/ mixed boundary conditions, (-1,0,1) may not be enough
C         Patched 10 Apr 02
          if (lmod(1) == 0 .or. lmod(2) == 0 .or. lmod(3) == 0) goto 5
        endif

C   --- Check to see if p was sufficiently shortened ---
        if (mode(1) >= 10 .and. lmod(1) == 2) then
          dpin  = dsqrt(p0(1)**2 + p0(2)**2 + p0(3)**2)
          dpout = dsqrt(pout(1,ib)**2 + pout(2,ib)**2 + pout(3,ib)**2)
          if (abs(dpout-dpin) < toll) then
            pout(1,ib) = p0(1)
            pout(2,ib) = p0(2)
            pout(3,ib) = p0(3)
          endif
        endif

C   --- pout <- pout - pin, units of plat ---
        if (nbas < 0) then
          call dpadd(p0,pout(1,ib),1,3,-1d0)
          do  i = 1, 3
            xx = -p0(1)*qlat(1,i) - p0(2)*qlat(2,i) - p0(3)*qlat(3,i)
            if (dabs(xx-nint(xx)) > 1d-10) call rx('bug in shorps')
            pout(i,ib) = xx
          enddo
        endif
      enddo

C     call prmx('plat',plat,3,3,3)
C     call prmx('qlat',qlat,3,3,3)
C     call prmx('shortened pos',pout,3,3,nbas1)
C     mc out.dat -t qlat -x gives pos as multiples of plat

C     call tcx('shorps')
      end

C#ifdefC TEST
C      subroutine fmain
CC to see that the change in position vectors are multiples of
CC the lattice vector, copy input,output pos to 'pos','posf'; invoke
CC mc posf pos -- -t plat -t -i -x
C      implicit none
C      integer ix(3),i
C      double precision pos(3,48),pos2(3,48),plat(3,3),qlat(9),xx
C
C      data plat /
C     .  0.5d0,          .5d0, 0d0,
C     .  0.0d0,          0.d0, 1d0,
C     .  2.570990255d0, -2.570990255d0, 0d0/
C      data pos /
C     .  -0.697107d0,  1.197107d0,  0.250000d0,
C     .  -0.697107d0,  1.197107d0,  0.750000d0,
C     .  -0.770330d0,  0.770330d0,  0.000000d0,
C     .  -0.770330d0,  0.770330d0,  0.500000d0,
C     .  -0.343553d0,  0.843553d0,  0.250000d0,
C     .  -0.343553d0,  0.843553d0,  0.750000d0,
C     .  -0.416777d0,  0.416777d0,  0.000000d0,
C     .  -0.416777d0,  0.416777d0,  0.500000d0,
C     .   0.010000d0,  0.490000d0,  0.250000d0,
C     .   0.010000d0,  0.490000d0,  0.750000d0,
C     .   0.250000d0,  0.250000d0,  0.500000d0,
C     .   0.500000d0,  0.500000d0,  0.750000d0,
C     .   0.750000d0,  0.750000d0,  1.000000d0,
C     .   1.000000d0,  1.000000d0,  1.250000d0,
C     .   0.250000d0, -0.250000d0,  0.000000d0,
C     .   0.500000d0,  0.000000d0,  0.250000d0,
C     .   0.750000d0,  0.250000d0,  0.500000d0,
C     .   1.000000d0,  0.500000d0,  0.750000d0,
C     .   0.750000d0, -0.250000d0,  0.500000d0,
C     .   1.000000d0,  0.000000d0,  0.750000d0,
C     .   1.250000d0,  0.250000d0,  1.000000d0,
C     .   1.500000d0,  0.500000d0,  1.250000d0,
C     .   0.740000d0, -0.740000d0,  0.000000d0,
C     .   0.740000d0, -0.740000d0,  0.500000d0,
C     .   1.166777d0, -0.666777d0,  0.250000d0,
C     .   1.166777d0, -0.666777d0,  0.750000d0,
C     .   1.093553d0, -1.093553d0,  0.000000d0,
C     .   1.093553d0, -1.093553d0,  0.500000d0,
C     .   1.520330d0, -1.020330d0,  0.250000d0,
C     .   1.520330d0, -1.020330d0,  0.750000d0,
C     .   1.447107d0, -1.447107d0,  0.000000d0,
C     .   1.447107d0, -1.447107d0,  0.500000d0,
C     .  -1.050660d0,  1.550660d0,  0.250000d0,
C     .  -1.050660d0,  1.550660d0,  0.750000d0,
C     .  -1.123883d0,  1.123883d0,  0.000000d0,
C     .  -1.123883d0,  1.123883d0,  0.500000d0,
C     .   1.873883d0, -1.373883d0,  0.250000d0,
C     .   1.873883d0, -1.373883d0,  0.750000d0,
C     .   1.800660d0, -1.800660d0,  0.000000d0,
C     .   1.800660d0, -1.800660d0,  0.500000d0,
C     .  -1.404214d0,  1.904214d0,  0.250000d0,
C     .  -1.404214d0,  1.904214d0,  0.750000d0,
C     .  -1.477437d0,  1.477437d0,  0.000000d0,
C     .  -1.477437d0,  1.477437d0,  0.500000d0,
C     .   2.227437d0, -1.727437d0,  0.250000d0,
C     .   2.227437d0, -1.727437d0,  0.750000d0,
C     .   2.154214d0, -2.154214d0,  0.000000d0,
C     .   2.154214d0, -2.154214d0,  0.500000d0/
C
C
C      call prmx('plat',plat,3,3,3)
C      call prmx('starting pos',pos,3,3,48)
C      ix(1) = 2
C      ix(2) = 2
C      ix(3) = 2
C      call shorps(48,plat,ix,pos,pos2)
C      call prmx('final pos',pos2,3,3,48)
C
C      call mkqlat(plat,qlat,xx)
C      do  10  i = 1, 48
C   10 call shorbz(pos(1,i),pos2(1,i),plat,qlat)
C
C      call prmx('from shorbz',pos2,3,3,48)
C      end
C#endif
C#ifdefC TEST2
C      subroutine fmain
CC Check special case in which a bug fix, mixed boundary conditions
C      implicit none
C      integer ix(3),i
C      double precision pos(3,1),pos2(3,1),plat(3,3),qlat(9),xx
C      double precision dd1,dd2
C
C      call wkinit(10000)
C
C      data plat /-0.5d0,0.5d0,0.0d0,
C     .            0.0d0,0.0d0,1.0d0,
C     .            7.0d0,7.0d0,4.0d0/
C      data ix /2,2,0/
C
C      pos(1,1) = 2.5d0
C      pos(2,1) = 3.0d0
C      pos(3,1) = 0.0d0
C
C      dd1 = dsqrt(pos(1,1)**2 + pos(2,1)**2 + pos(3,1)**2)
C      call shorps(1,plat,ix,pos,pos2)
C      dd2 = dsqrt(pos2(1,1)**2 + pos2(2,1)**2 + pos2(3,1)**2)
C      print *, dd1, dd2
C
C      call mkqlat(plat,qlat,xx)
C      do  10  i = 1, 1
C   10 call shorbz(pos(1,i),pos2(1,i),plat,qlat)
C      dd2 = dsqrt(pos2(1,1)**2 + pos2(2,1)**2 + pos2(3,1)**2)
C      print *, dd1, dd2
C
C      end
C#endif
C#ifdefC TEST3
C      subroutine fmain
CC- shorps for a set of vectors
C      implicit none
C      integer ifi,i,fopng,rdm,nr,nc,ix(3),nsw
C      double precision xx,plat(3,3*2),dd,dd2
C      real(8),allocatable:: pos(:,:),pos2(:,:)
C      logical cmdstr
C      character posf*256,platf*256,sw*128
C      data ix /2,2,2/
C
C      nsw = 0
C
C      if (.not. cmdstr(1,sw)) goto 999
C      if (sw == '-quad1') then
C        ix = 1
C        nsw = nsw + 1
C      endif
C
C      if (.not. cmdstr(2+nsw,platf)) goto 999
C      ifi = fopng(platf,-1,0)
C      nr = 0
C      nc = 0
C      i = rdm(ifi,0,0,' ',xx,nr,nc)
C      if (nc /= 3) then
C        call info2(0,0,0,
C     .    ' shorps expected plat file with 3 columns, found %i',nc,0)
C        goto 999
C      endif
C      if (nr /= 3) then
C        call info2(0,0,0,
C     .    ' shorps expected plat file with 3 rows, found %i',nr,0)
C        goto 999
C      endif
C      rewind ifi
C      i = rdm(ifi,10000,nr*3,' ',plat,nr,nc)
C      if (i < 0) call rx('shorps failed to read plat file')
C      call fclose(ifi)
C
C      write(*,321) plat(:,1), plat(:,2), plat(:,3)
C  321 format(' P1=',3f12.6/' P2=',3f12.6/' P3=',3f12.6)
C
C      if (.not. cmdstr(1+nsw,posf)) goto 999
C      ifi = fopng(posf,-1,0)
C      nr = 0
C      nc = 0
C      i = rdm(ifi,0,0,' ',xx,nr,nc)
C      if (nc /= 3) then
C        call info2(0,0,0,
C     .    ' shorps expected pos file with 3 columns, found %i',nc,0)
C        goto 999
C      endif
C      call info2(0,0,0,' pos file has %i vectors',nr,0)
C      allocate(pos(3,nr*2))
C      allocate(pos2(3,nr))
C      rewind ifi
C      i = rdm(ifi,10000,nr*6,' ',pos,nr,nc)
C      if (i < 0) call rx('shorps failed to read pos file')
C      call fclose(ifi)
C
C      ix(1) = ix(1) + 60
C      call shorps(nr,plat,ix,pos,pos2)
C      call info0(0,0,0,'%18ppos%44pd%66pshortened%95pd')
C      do  i = 1, nr
C        dd  = dsqrt(pos(1,i)**2 + pos(2,i)**2 + pos(3,i)**2)
C        dd2 = dsqrt(pos2(1,i)**2 + pos2(2,i)**2 + pos2(3,i)**2)
C        print 543, pos(:,i),dd,pos2(:,i),dd2
C  543   format(3f12.6,2x,f11.5,2x,3f12.6,2x,f11.5)
C      enddo
C
C      return
C  999 continue
C      print *, 'Usage:'
C      print *, 'shorps pos plat'
C      print *, 'plat as row vectors (same as lm output)'
C
C      end
C#endif
