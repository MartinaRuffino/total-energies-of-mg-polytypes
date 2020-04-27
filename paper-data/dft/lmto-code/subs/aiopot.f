      logical function aiopot(nr,nsp,a,rmax,bhat,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci Inputs/Outputs
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   rmax  :augmentation radius, in a.u.
Ci         :If any of these are zero, aiopot sets those with value 0
Ci         :to file value and returns without reading potential
Ci   the next two are read for ifi>0 and written for ifi<0
Ci   bhat: :orientation of magnetic field (optional)
Ci         :If file write AND bhat(1) <> NULLR, write bhat in header
Ci   v     :potential
Ci   v(1)  :(read only) If input v(1)=NULLR, do not rewind file
Co Outputs
Co   v   :potential read from idst if ifi>0 (file read)
Co   bhat:orientation of magnetic field (optional)
Co        If file read AND bhat(1)<>NULLR, read bhat from POT: line
Cr Remarks
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, aiopot returns true only if the category is found,
Cr    the file's value of a and nr match input and rmax is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
Cu Updates
Cu   25 Sep 17 aiopot can suppress rewinding file
Cu   07 Aug 16 aiopot can interpolate radial meshes
Cu    4 Apr 04 Optional I/O of bfield orientation
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,bhat(3)
      real(8),target :: v(nr,nsp)
C ... Dynamically allocated arrays
      real(8),allocatable :: rofi(:),rofi2(:)
      real(8),pointer :: vloc(:,:)
C ... Local parameters
      integer i,j,isp,nr2,nsp2,ix2(3),jx,npoly
      double precision a2,rmax2,dvmx,dv
      real(8),parameter :: NULLR=-99999
      logical ltmp
      character strn*72
      procedure(integer) a2vec
      procedure(logical) scat

      aiopot = .false.
      if (ifi > 0) then
        ltmp = .true.; if (v(1,1) == NULLR) ltmp = .false.
        if (.not. scat(ifi,'POT:',':',ltmp)) return
        if (.not. (bhat(1) == -99d0 .or. bhat(1) == NULLR)) then
          bhat(1) = 0
          bhat(2) = 0
          bhat(3) = 0
          backspace ifi
          strn = ' '
          read(ifi,'(a72)') strn
          call words(strn,i)
          if (i >= 5) then
            call word(strn,3,i,j)
            i = i-1
            j = a2vec(strn,len(strn),i,4,' ',1,-2,-3,ix2,bhat)
            if (j /= 3) call rxs(
     .        'aiopot: failed to parse Bfield in line: ',strn)
C           print *, bhat
          endif
        endif
        read(ifi,*,err=15) nr2,nsp2,a2,rmax2
        ltmp = nr*nsp*a*rmax == 0
        if (nr == 0) nr=nr2
        if (nsp == 0) nsp=nsp2
        if (a == 0) a=a2
        if (rmax == 0) rmax=rmax2
        if (ltmp) then ! dimensioning parameters read; v cannot be read
          backspace ifi; backspace ifi
          aiopot = .true.
          return
        endif
C       Set up for interpolation
        vloc => v
        if (.not. (a2 == a .and. nr == nr2 .and. dabs(rmax2-rmax) <= 1d-5)) allocate(rofi(nr),rofi2(nr2),vloc(nr2,2))
C       Read the potential
        do  isp = 1, min0(nsp2, nsp)
          read(ifi,101) (vloc(i,isp), i = 1, nr2)
          do  i = 1, nr2
            vloc(i,nsp) = vloc(i,isp)
          enddo
        enddo
        aiopot = .true.
C       Do the interpolation
        if (.not. (a2 == a .and. nr == nr2 .and. dabs(rmax2-rmax) <= 1d-5)) then
C         call prrmsh('V(init) ',rofi2,vloc,nr2,nr2,nsp)

C         General polynomial interpolation.  Works, but slow.
          call radmsh(rmax,a,nr,rofi); call radmsh(rmax2,a2,nr2,rofi2)
          dvmx = 0
          jx = 0
          do  i = 1, nr
            npoly = 7 ; if (rofi(i) > rofi2(nr2)) npoly = 5
            do  isp = 1, nsp
              call polint(rofi2,vloc(1,isp),nr2,npoly,rofi(i),0d0,0,jx,v(i,isp),dv)
              dvmx = max(dvmx,abs(dv))
            enddo
          enddo
C         call prrmsh('V(final) ',rofi,v,nr,nr,nsp)
          deallocate(rofi,rofi2,vloc)
          call info8(20,0,0,' interpolate V(a=%d,nr=%i,rmax=%,6;6d) to '//
     .      ' v(a=%d,nr=%i,rmax=%,6;6d)  errmx=%,3;3g',
     .      a2,nr2,rmax2,a,nr,rmax,dvmx,8)
        endif
      else
        if (bhat(1) == -99d0 .or. bhat(1) == NULLR) then
          write(-ifi,'(''POT:'')')
        else
          write(-ifi,'(''POT:   bhat='',3f12.7)') bhat
        endif
        write(-ifi,102) nr,nsp,a,rmax
        do isp = 1, nsp
          write(-ifi,101) (v(i,isp),i = 1,nr)
        enddo
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f15.10)
   15 continue
      end
