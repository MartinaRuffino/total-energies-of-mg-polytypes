      subroutine shorho(opt,strn,nr,nlm,nsp,ri,rwgt,rl1,rl2)
C- Prints out quantities related to site densities rho1, rho2
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1 info for rho1
Ci         :2 info for rho2
Ci         :3 combination 1+2
Ci   strn  :printout string
Ci   a     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nlm   :number of L channels
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rl1   :First density * r**2, YL repsn'
Ci   rl2   :Second density * r**2, YL repsn'
Co Outputs
Co   printout only
Cu Updates
Cu   23 Oct 13 (DMT) 'rli' vec2->scalar not to give ifort v14 hassle
Cu   24 Jan 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character strn*(*)
      integer opt,nr,nlm,nsp
      double precision ri(nr),rwgt(nr),rl1(nr,nlm,nsp),rl2(nr,nlm,nsp)
C ... Local parameters
      integer ir,ip,i,is,opt0,opt1,stdo,nglob,imaxmn(2,2,2),ineg(2)
      double precision ddot,pi,srfpi,q1p,q1m,q2p,q2m,q1,q2,a1,a2
      double precision r1p,r1m,r2p,r2m,rli,rmaxmn(2,2,2)
      integer nnn,lmax,ll,nxl(0:7),nth,nph,np,iprint
      parameter(nnn=122)
      double precision p(3,nnn),wp(nnn),p2(nnn*3),r2(nnn)
      real(8), allocatable :: yl(:),rps(:,:,:)
      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/

      if (iprint() == 0) return
      call info0(0,1,0,trim(strn))

      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4*pi)
      opt0 = mod(opt,2)
      opt1 = mod(opt/2,2)
      stdo = nglob('stdo')
      lmax = ll(nlm)
      nth = nxl(min(lmax,7))
      nph = 0
      call fpiint(nth,nph,np,p,wp)
      call dmcpy(p,1,3,p2,np,1,np,3)
      allocate (yl(np*nlm),rps(nr,np,nsp))
      call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax,np,yl,r2)

C     nsp = 1

      if (nsp == 1) then
        write(stdo,332)
  332   format(10x,'q',12x,'rho(MT)')
      else
        print 333
  333   format(10x,'q+',11x,'q-',11x,'q',12x,'mom',10x,'rho+(MT)',
     .    5x,'rho-(MT)')
      endif
      if (opt0 == 1) then
        q1p = srfpi*ddot(nr,rwgt,1,rl1,1)
        q1m = srfpi*ddot(nr,rwgt,1,rl1(1,1,nsp),1)
        q1  = (q1p+q1m)/(3-nsp)
        a1  = q1p-q1m
        r1p = rl1(nr,1,1)/ri(nr)**2/srfpi
        r1m = rl1(nr,1,nsp)/ri(nr)**2/srfpi
        if (nsp == 1) then
          write(stdo,334) 'rho1',q1p,r1p
        else
          write(stdo,334) 'rho1',q1p,q1m,q1,a1,r1p,r1m
  334     format(1x,a,6f13.7)
        endif
      endif

      if (opt1 == 1) then
        q2p = srfpi*ddot(nr,rwgt,1,rl2,1)
        q2m = srfpi*ddot(nr,rwgt,1,rl2(1,1,nsp),1)
        q2  = (q2p+q2m)/(3-nsp)
        a2  = q2p-q2m
        r2p = rl2(nr,1,1)/ri(nr)**2/srfpi
        r2m = rl2(nr,1,nsp)/ri(nr)**2/srfpi
        if (nsp == 1) then
          write(stdo,334) 'rho2',q2p,r2p
        else
          write(stdo,334) 'rho2',q2p,q2m,q2,a2,r2p,r2m
        endif
      endif

      if (opt0 == 1 .and. opt1 == 1) then
        if (nsp == 1) then
          write(stdo,334) 'diff',q1p-q2p,r1p-r2p
        else
          write(stdo,334) 'sum ',q1p+q2p,q1m+q2m,q1+q2,a1+a2,
     .      r1p+r2p,r1m+r2m
          write(stdo,334) 'diff',q1p-q2p,q1m-q2m,q1-q2,a1-a2,
     .      r1p-r2p,r1m-r2m
        endif
      endif

C ... Scale rl to true density
      if (opt0 == 1) call vxcns3(nr,nlm,nsp,ri,rl1,1)
      if (opt1 == 1) call vxcns3(nr,nlm,nsp,ri,rl2,1)

C ... Find minimum, maximum in l=0 density
      if (nsp == 1) then
        write(stdo,331)
  331   format('  l=0',
     .    6x,'rhmin ',6x,'r(rhmin)',6x,'rhmax ',6x,'r(rhmax)')
      else
        write(stdo,335)
  335   format('  l=0',
     .    6x,'rhmin+',6x,'r(rhmin)',6x,'rhmin-',6x,'r(rhmin)',
     .    6x,'rhmax+',6x,'r(rhmax)',6x,'rhmax-',6x,'r(rhmax)')
      endif
      call dvset(rmaxmn(1,1,1),1,4,9d9)
      call dvset(rmaxmn(1,1,2),1,4,-9d9)
      do  i = 1, 2
        if (opt0 == 0 .and. i == 1) cycle
        if (opt1 == 0 .and. i == 2) cycle
        do  ir = 1, nr
          do  is = 1, nsp
            if (i == 1) rli = rl1(ir,1,is)/srfpi
            if (i == 2) rli = rl2(ir,1,is)/srfpi
            if (rmaxmn(is,i,1) > rli) then
              imaxmn(is,i,1) = ir
              rmaxmn(is,i,1) = rli
            endif
            if (rmaxmn(is,i,2)*ri(ir)**2 < rli*ri(ir)**2) then
              imaxmn(is,i,2) = ir
              rmaxmn(is,i,2) = rli
            endif
          enddo
        enddo

        if (nsp == 1) then
          write(stdo,336) i,
     .      rmaxmn(1,i,1),ri(imaxmn(1,i,1)),
     .      rmaxmn(1,i,2),ri(imaxmn(1,i,2))
        else
          write(stdo,336) i,
     .      rmaxmn(1,i,1),ri(imaxmn(1,i,1)),
     .      rmaxmn(2,i,1),ri(imaxmn(2,i,1)),
     .      rmaxmn(1,i,2),ri(imaxmn(1,i,2)),
     .      rmaxmn(2,i,2),ri(imaxmn(2,i,2))
  336     format(' rho',i1,4(f13.6,f13.7))
        endif

      enddo

      call dvset(rmaxmn(1,1,1),1,4,9d9)
      call dvset(rmaxmn(1,1,2),1,4,-9d9)
      write(stdo,338)
  338 format(' all l')

      do  i = 1, 2
        if (opt0 == 0 .and. i == 1) cycle
        if (opt1 == 0 .and. i == 2) cycle
        do  is = 1, nsp
          if (i == 1)
     .      call dgemm('N','T',nr,np,nlm,1d0,rl1(1,1,is),nr,
     .      yl,np,0d0,rps(1,1,is),nr)
          if (i == 2)
     .      call dgemm('N','T',nr,np,nlm,1d0,rl2(1,1,is),nr,
     .      yl,np,0d0,rps(1,1,is),nr)

          do  ip = 1, np
          do  ir = 1, nr
            if (rps(ir,ip,is) < 0d0) then
              ineg(is) = ineg(is) + 1
            endif
            if (rmaxmn(is,i,1) > rps(ir,ip,is)) then
              imaxmn(is,i,1) = ir
              rmaxmn(is,i,1) = rps(ir,ip,is)
            endif
            if (rmaxmn(is,i,2)*ri(ir)**2 < rps(ir,ip,is)*ri(ir)**2)
     .        then
              imaxmn(is,i,2) = ir
              rmaxmn(is,i,2) = rps(ir,ip,is)
            endif
          enddo
          enddo
        enddo

        if (nsp == 1) then
          write(stdo,336) i,
     .      rmaxmn(1,i,1),ri(imaxmn(1,i,1)),
     .      rmaxmn(1,i,2),ri(imaxmn(1,i,2))
        else
          write(stdo,336) i,
     .      rmaxmn(1,i,1),ri(imaxmn(1,i,1)),
     .      rmaxmn(2,i,1),ri(imaxmn(2,i,1)),
     .      rmaxmn(1,i,2),ri(imaxmn(1,i,2)),
     .      rmaxmn(2,i,2),ri(imaxmn(2,i,2))
        endif

      enddo

C ... Undo the rl scaling
      if (opt0 == 1) call vxcns3(nr,nlm,nsp,ri,rl1,0)
      if (opt1 == 1) call vxcns3(nr,nlm,nsp,ri,rl2,0)

      end
