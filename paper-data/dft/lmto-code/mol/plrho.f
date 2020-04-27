      subroutine plrho(n0,lsw,vals,rhoi,nbas,ips,ioff,lmxl,a,nr,rsm,rmt,
     .  lcore,nxi,lxi,exi,alat,pos,orho,cy,
     .  nphi,ephi,lphi,nel,nhs,evl,evc)
C- PLots charge density or wave function
C        prm i1,2,3 i0 x1 x2  y1 y2  nx ny   dist rfac
C  vals   123 1 2 3 1  -5  5  -5 5   31 31   0.0  0.7  2  41
C  Notes: prm: permutation index permuting x,y,z
C         i1,i2,i3: 3 atoms defining plane of contour
C         i0,dist:  atom defining origin, and offset normal to plane
C         x1,x2,y1,y2,n1,n2: plotting range in plane and no. points
C         last three not used now.
C   lsw: 0 for charge density, 1 for wave function
      implicit none
      integer n0,lsw,nbas,ips(1),orho(1),lmxl(1),nr(1),nxi(1),lxi(n0,1),
     .  ioff(1),nphi(1),lphi(n0,1),lcore,nel,nhs,lsp
      double precision rhoi(1),pos(3,1),rsm(1),rmt(1),cy(1),
     .  a(1),exi(n0,1),ephi(n0,1),evl(1),evc(nhs,1),alat,atrans(3,3)
      double precision vals(15),vx(3),vy(3),dist,rfac,x1,x2,y1,y2
      integer ox,oy,oz,ou,i1mach,np,nx,ny,m,i0,i1,i2,i3,
     .  iperm,ist1,ib0,lsqr
      character*(256) strn,fn
      logical cmdopt
      real w(1)
      common /w/ w

c ---------- check if spin-pol ----------------
      if(lsp() /= 0)
     .   call rx('PLRHO: the program is not yet ready for spins')


C ... default values for contour plots
      call dpzero(vals,15)
      vals(1) = 123
      vals(2) = 1
      vals(3) = 2
      vals(4) = 3
      vals(5) = 1
      vals(6) = -10
      vals(7) = 10
      vals(8) = -10
      vals(9) = 10
      vals(10) = 51
      vals(11) = 51
      vals(12) = 0
      fn = 'rho'

      if (cmdopt('--plot',6,0,strn)) then
        call parvls(strn(7:),vals,fn)
      endif

c --------- input specifications for plot -----
      iperm=idnint(vals(1))
      i1=idnint(vals(2))
      i2=idnint(vals(3))
      i3=idnint(vals(4))
      i0=idnint(vals(5))
      x1=vals(6)
      x2=vals(7)
      y1=vals(8)
      y2=vals(9)
      nx=idnint(vals(10))
      ny=idnint(vals(11))
      dist=vals(12)
      rfac=vals(13)

      call awrit5(' plrho: i0: %d i1,i2,i3:%3:1d  x:%2:1d  y:%2:1d'//
     .  '  nx,ny:%2:1d',
     .  ' ',80,i1mach(2),vals(5),vals(2),vals(6),vals(8),vals(10))

c --------- define mesh of points ------
      do 33 m=1,3
      vx(m)=pos(m,i2)-pos(m,i1)
   33 vy(m)=pos(m,i3)-pos(m,i1)
      call defrr(ou,    nx*ny)
      call defrr(ox,    nx*ny)
      call defrr(oy,    nx*ny)
      call defrr(oz,    nx*ny)
      call plch1(iperm,vx,vy,x1,x2,y1,y2,nx,ny,pos(1,i0),dist,
     .   np,w(ox),w(oy),w(oz),atrans)
      call trpos(nbas,pos,i1,atrans)
      call dpzero(w(ou),  np)

C --- Plot charge density or wave function ---
   50 continue
      lsqr = 0
      if (lsw == 0) then
C   ... Interstitial density
        call plch2(nbas,ips,pos,rsm,nxi,lxi,exi,n0,ioff,cy,
     .   rhoi,np,w(ox),w(oy),w(oz),w(ou))
c   ... Overwrite points inside spheres with full density
        call plch3(rmt,a,nr,lmxl,nbas,pos,ips,cy,orho,np,
     .    w(ox),w(oy),w(oz),w(ou))
      elseif (lsw == 1) then
        call dpzero(w(ou),  np)
        ib0=0
        write(*,*) 'enter evl # to plot (0 to quit, -# for psi^2)'
        read (*,*) ist1
        if (ist1 < 0) lsqr = 1
        ist1 = iabs(ist1)
        if (ist1 == 0) call rx('plrho')
        write(6,566) ist1, evl(ist1)
  566   format(/'state',i4,'    eigenvalue=',f10.5)
        call plwv2(nbas,ib0,ips,pos,nel,rsm,nphi,lphi,ephi,n0,cy,
     .    ist1,nhs,evc,np,w(ox),w(oy),w(oz),w(ou))
      else
        call rx('plrho: bad lsw')
      endif

C --- write to file 81 ---
      call mopen(81,fn,'f')
      rewind 81
      call awrit2('%% rows %i cols %i',' ',80,81,nx,ny)
      if (lsqr == 0) call writ81(w(ou),nx*ny)
      if (lsqr == 1) call writ81a(w(ou),nx*ny)
      call mclose(81)
CL      write(71,722) vx,x1,x2,nx,vy,y1,y2,ny,i1,i2,i3,dist,rfac
  722 format(' ro  ----------- chden plot -------------------'/
     .  ' ro  x: (',3f8.4,' ) ',2f9.5,i5/
     .  ' ro  y: (',3f8.4,' ) ',2f9.5,i5/
     .  ' ro  i1,i2,i3',3i5,'   dist',f9.4,'   rfac',f8.3)

      if (lsw == 1) goto 50
      end
      subroutine writ81(u,n)
      real*8 u(n)
      write(81,810) u
  810 format(1p,5e14.6)
      end
      subroutine writ81a(u,n)
      real*8 u(n)
      do  10  i = 1, n
   10 u(i) = u(i)**2
      write(81,810) u
  810 format(1p,5e14.6)
      end
