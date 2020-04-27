      subroutine dlmmxy(mode,nl,nspc,nclass,nangl,dclabl,mxy,qnu,th,bxc)
C- Print off-diagonal components of local moments
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode   :0 do not update bxc
Ci          :1 update bxc
Ci   ic1,ic2:CPA class range to be processed
Ci   nclass :number of main classes
Ci   nangl  :number of CPA classes
Ci   nl     :number of l channels
Ci   nspc   :number of coupled spins
Ci   dclabl :class name, packed as a real number
Ci   mxy    :x,y components of local moments
Ci   qnu    :charge moments (to find mz)
Ci   th     :angles for all classes
Cio  bxc    :directions of xc fields
Cr Remarks
Cb Bugs
Cu Updates
Cu   10 Apr 12 (Belashchenko) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,nspc,nclass,nangl
      real(8) dclabl(nclass+nangl),mxy(nl,2,*),qnu(3,nl,2,*),bxc(3,*),
     .  th(nangl,nspc)
C ... Local parameters
      integer stdo,ipr,lgunit,l,ic,ic1,ic2
      character*8 clabl
      real(8) mx,my,mz,mp,dsum,bxold,phi,proj
      integer mpipid,procid,master

      procid = mpipid(1)
      master = 0

      stdo = lgunit(1)
      call getpr(ipr)

      if (procid == master .and. ipr >= 30) then
        write(stdo,300)
 300    format(/' Off-diagonal local moments and constraining fields')
        if (nl == 3) then
          write(stdo,301)
 301      format(12x,'Mx-s',4x,'Mx-p',4x,'Mx-d',5x,'Mx',6x,'My',6x,
     .      'Mz',4x,'Old Bx  New Bx')
        else
          write(stdo,302)
 302      format(12x,'Mx-s',4x,'Mx-p',4x,'Mx-d',4x,'Mx-f',5x,'Mx',6x,
     .      'My',6x,'Mz',4x,'Old Bx  New Bx')
        endif
      endif

      ic1 = nclass + 1
      ic2 = nclass + nangl

      do  ic = ic1, ic2
        call r8tos8(dclabl(ic),clabl)
        mx = dsum(nl,mxy(1,1,ic),1)
        my = dsum(nl,mxy(1,2,ic),1)
        mz = dsum(nl,qnu(1,1,1,ic),3) - dsum(nl,qnu(1,1,2,ic),3)
        bxold = bxc(1,ic)
C   ... Update the constraining field
        if (mode == 1) then
          if (nspc == 1) then
            bxc(1,ic) = dtan(datan(bxc(1,ic))-datan(mx/mz))
          else
            mp = sqrt(mx*mx+my*my)
            phi = th(ic-nclass,2)
            proj = dcos(phi) * mx + dsin(phi) * my
            if (proj < 0) mp = - mp
            bxc(1,ic) = dtan(datan(bxc(1,ic))-datan(mp/mz))
          endif
        endif
C   ... Printout
        if (procid == master .and. ipr >= 30) then
          if (nl == 3) then
C            write(stdo,311)clabl,(mxy(l,1,ic),l=1,3),
C     .        mx,my,mz,bxold,bxc(1,ic)
            call info8(30,0,0,' '//clabl//'%3;8,4D%;8,4D%;8,4D%;8,4D%;8,4D%;8,4D',
     .        mxy(1,1,ic),mx,my,mz,bxold,bxc(1,ic),7,8)
          else
C            write(stdo,312)clabl,(mxy(l,1,ic),l=1,4),
C     .        mx,my,mz,bxold,bxc(1,ic)
            call info8(30,0,0,' '//clabl//'%4;8,4D%;8,4D%;8,4D%;8,4D%;8,4D%;8,4D',
     .        mxy(1,1,ic),mx,my,mz,bxold,bxc(1,ic),7,8)
          endif
c         write(stdo,313)qnu(1,1,1,ic)-qnu(1,1,2,ic),
c    .      qnu(1,2,1,ic)-qnu(1,2,2,ic),
c    .      qnu(1,3,1,ic)-qnu(1,3,2,ic),qnu(1,4,1,ic)-qnu(1,4,2,ic)
        endif
      enddo

 311  format(1x,a,4F8.4,2F8.4,2F8.4)
 312  format(1x,a,5F8.4,2F8.4,2F8.4)
 313  format('  mz:    '4F8.4)

      end
