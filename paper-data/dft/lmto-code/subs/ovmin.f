      module ovstat
! Not the best way to deal with the old common block alignment but probably the quickest.
        implicit none

        integer, parameter :: mxlst = 256
        real(8) :: alato,plato(9)
        integer :: nbaso,nbaspo,nlst,modeo(3),ilst(mxlst)
        integer, pointer :: p_ntab(:),p_iax(:),ipsl(:)
        real(8), pointer :: posl(:),zl(:),rmaxl(:)
      end module ovstat

      subroutine ovmin(sovmin,nbas,nbasp,alat,plat,rmax,hcr,clabl,
     .  ips,mode,z,ntab,iax,pos,iprtbl)
C- Check volume and sphere overlaps, optionally minimizing them
C ----------------------------------------------------------------
Ci Inputs
Ci   sovmin: a set of modifiers, with the syntax
Ci         : --mino[:dxmx=#][:xtol=#][:style=#]:site-list
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   rmax  :potential radius, in a.u.
Ci   hcr   :augmentation radius, in a.u.
Ci   clabl :species or class labels
Ci   ips   :species (or class) table: site ib belongs to species ips(ib)
Ci   mode  :vector of length 3 governing how pos shortened (see shorps)
Ci   z     :nuclear charge
Ci   ntab  :ntab(ib)=# pairs in iax table preceding ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   pos   :basis vectors
Ci   iprtbl: nonzero if to call ovlchk and print table of overlaps
Co Outputs
Co   Sphere overlaps are printed out
Cr Remarks
Cr   hcr   :hard core radius (printout only).  hcr(1)=0 => not used.
Cu Updates
Cu   08 Nov 17 prints out ES max overlaps separately
Cu   25 Sep 17 New nkill
Cu   24 Jul 15 Replace dclabl with clabl
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   14 Sep 11 Restore use of hcr
Cu   22 Oct 02 weight ES-ES and atom-ES overlaps differently when
Cu             minimizing sphere overlap positions
Cu    9 Dec 98 replace call to frpmin with call to gradzr.
Cu    8 Sep 98 small patches in minimizing algorithm
Cu   24 Nov 97 changed ovmin to call fovlp for fast execution
C ----------------------------------------------------------------
      use ovstat
      implicit none
C ... Passed parameters
      integer nbas,nbasp,iprtbl
      integer niax
      parameter (niax=10)
      integer,target :: ntab(nbas+1),iax(niax)
      double precision plat(3,3),pos(3,nbasp),rmax(*),hcr(*),z(*),alat
      integer ips(nbas),mode(3)
      character sovmin*(*)
      character*8 clabl(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: list(:)
      real(8), allocatable :: p(:)
C ... Local parameters
      double precision fovl,xtol,gtol,dxmn,dxmx,ovcall,fovmx(3),xx
      double precision wk(0:27)
      integer isw,ir,i,j,j1,j2,ls,m,lstyle,iter,
     .  iv,nlstc,mxint,nclass,ib,ic,maxit,ipr,n,nkill
      integer, parameter :: NULLI=-99999
      character dc*1
      external ovcall,mxint
C     static for ovcall
      integer novl
      procedure(integer) :: i1mach,parg,iclbsj,wordsw

C --- Print out positions and overlaps ---
      call getpr(ipr)
C      if (iprtbl > 0)
      call ovlchk(nbas,nbasp,pos,alat,rmax,hcr,clabl,ips,mode,plat,fovmx,xx)
      call fovlp(1,nbas,ntab,iax,plat,pos,ips,alat,rmax,z,6d0,
     .  1d0,.75d0,.5d0,fovmx,fovl,novl)
      i = 3; if (fovmx(3) == NULLI) i = 2; if (fovmx(2) == NULLI) i = 1
      if (novl == 0) novl = 1
      if (ipr >= 10 .or. iprtbl > 0)
     . call info5(2,1,0,' OVMIN, %i pairs:  '//
     .  'fovl = %;6g   <ovlp> = %;1d%%   max ovlp = %s,%n;1d%%',
     .  ntab(nbas+1),fovl/novl,(fovl/novl)**(1/6d0)*100,i,fovmx*100)

C --- Minimize overlaps wrt positions in list ---
      if (sovmin /= ' ') then
C   ... Default values for gradzr call
        xtol = 2d-4
        gtol = 1d-5
        dxmn = 1d-6
        dxmx = .10d0
        maxit = 20
        isw = 10051
        nkill = 0

        ls = len(sovmin)
        j1 = 1
        dc = sovmin(j1:j1)
        j1 = j1+1
        lstyle = 0

C   ... Return here to resume parsing for arguments
   40   continue
        call nwordg(sovmin,0,dc//' ',1,j1,j2)

C   ... Parse special arguments
        if (j2 < len(sovmin))  then
        if (sovmin(j2+1:j2+1) /= ' ')  then
          m = j1-1
          i = parg('dxmx=',4,sovmin,m,ls,dc,1,1,iv,dxmx)
          m = j1-1
          i = parg('xtol=',4,sovmin,m,ls,dc,1,1,iv,xtol)
          m = j1-1
          i = parg('style=',2,sovmin,m,ls,dc,1,1,iv,lstyle)
          m = j1-1
          i = parg('maxit=',2,sovmin,m,ls,dc,1,1,iv,maxit)
          m = j1-1
          i = parg('nkill=',2,sovmin,m,ls,dc,1,1,iv,nkill)

          j1 = j2+2
          goto 40
        endif
        endif

C   ... List of all sites to move
        if (lstyle > 0) then
          nclass = mxint(nbas,ips)
          allocate(list(nclass))
          call clist(lstyle,sovmin(j1:j2+1),clabl,z,nclass,nlstc,list)
          nlst = 0
          do  i = 1, nlstc
            ic = list(i)
            do  j = 1, nbas
              ib = iclbsj(ic,ips,-nbas,j)
              if (ib < 0) exit
              nlst = nlst+1
              ilst(nlst) = ib
            enddo
          enddo
          deallocate(list)
        elseif (wordsw(sovmin,dc,'z',dc//' ',j)+wordsw(sovmin,dc,'Z',dc//' ',j) > 0) then
          nlst = 0
          do  ib = 1, nbasp
            ic = ips(ib)
            if (z(ic) == 0) then
              nlst = nlst+1
              ilst(nlst) = ib
            endif
          enddo
        else
          call mkilst(sovmin(j1:),nlst,ilst)
        endif
C        call awrit2(' min wrt:  %n:1i',' ',80,i1mach(2),nlst,ilst)
C        call awrit3(' setup:     xtol = %,2g   dxmx = %,2g   maxit = %i'
C     .    ,' ',80,i1mach(2),xtol,dxmx,maxit)
        call info2(10,0,0,' min wrt:  %n:1i',nlst,ilst)
        call info5(10,0,0,' setup:     xtol = %,2g   dxmx = %,2g   maxit = %i',xtol,dxmx,maxit,4,5)
        if (nlst <= 0) then
          call info0(10,0,0,' no sites in list ... no minimization')
          return
        endif

C  ...  set up static block for ovcall
        alato = alat
        nbaso = nbas
        nbaspo = nbasp
        p_ntab => ntab
        p_iax  => iax
        allocate(posl(3*nbasp))
        call dpcopy(pos,posl,1,3*nbasp,1d0)
        nclass = mxint(nbas,ips)
        allocate(zl(nclass))
        call dpcopy(z,zl,1,nclass,1d0)
        allocate(rmaxl(nbasp))
        call dpcopy(rmax,rmaxl,1,nbasp,1d0)
        allocate(ipsl(nbasp))
        call icopy(nbasp,ips,1,ipsl,1)
        call icopy(3,mode,1,modeo,1)
        call dpcopy(plat,plato,1,9,1d0)

C  ...  Initialization for gradzr call
        n = 3*nlst
        allocate(p(10*n)); call dpzero(p,10*n)
        ir = 0
        do  i = 1, nlst
          j = ilst(i)
          call dpscop(posl,p,3,3*j-2,3*i-2,1d0)
        enddo
        xx = ovcall(n,0d0,p,ir)
        iter = 0
        call pshpr(ipr-5)
   22   continue
        if (nkill > 0) then
          if (mod(iter,nkill) == 0) ir=0
        endif
        call gradzr(n,p,xx,dxmn,dxmx,xtol,gtol,1.0d0,wk,isw,ir)
        iter = iter+1
        xx = ovcall(n,0d0,p,ir)
        if (ir < 0 .and. iter < maxit) goto 22
        call poppr

C   ... Update positions
        do  i = 1, nlst
          j = ilst(i)
          call dpscop(p,pos,3,3*i-2,3*j-2,1d0)
        enddo
        deallocate(posl,zl,rmaxl,ipsl,p)

C   --- Print out updated positions and overlaps ---
        call info2(10,1,0,' OVMIN:  %?#n#incomplete ##convergence in %i iteration(s).'//
     .    '%N%9fUpdated site positions:',ir,iter)
        if (iprtbl > 0) call ovlchk(nbas,nbasp,pos,alat,rmax,0d0,
     .    clabl,ips,mode,plat,fovmx,xx)
        call fovlp(1,nbas,ntab,iax,plat,pos,ips,alat,rmax,z,6d0,
     .    1d0,.75d0,.5d0,fovmx,fovl,novl)
        if (novl == 0) novl = 1
        i = 3; if (fovmx(3) == NULLI) i = 2; if (fovmx(2) == NULLI) i = 1
        call info5(2,1,0,' minimized: fovl = %;6g   <ovlp> = %;1d%%'//
     .    '   max ovlp = %s,%n;1d%%',
     .    fovl/novl,(fovl/novl)**(1/6d0)*100,i,fovmx*100,5)
      endif

      end
      double precision function ovcall(n,x,p,ir)
C- Generic function call for projection grad fovl in a spec'd direction
Ci x,p,ir see gradzr
Cu   08 Jul 13 Replace f77 pointers with f90 ones
      use ovstat
      implicit none
C ... Passed parameters
      integer ir,n
      double precision x,p(3*n)
C ... Dynamically allocated local arrays
      real(8), allocatable :: posb(:)
C static:
      integer novl,novlp,novlm
C ... Local parameters
      logical cmdopt
      integer j,i,ix,ipr,lgunit,novl0
      double precision fovl,dx,val,fovp,fovm,pos(3),xx(3),fov0
      procedure(real(8)) :: ddot
      character*120 outs
      parameter (dx=1d-4)

C ... Save pos, other initialization
      call getpr(ipr)
      allocate(posb(3*nbaspo))
      call dpcopy(posl,posb,1,3*nbaspo,1d0)
      call pshpr(1)

      do  i = 1, nlst
        j = ilst(i)
        call dpscop(p,posl,3,3*i-2,3*j-2,1d0)
      enddo

      call ovlchk(nbaso,nbaspo,posl,alato,rmaxl,0d0,' ',
     .  ipsl,modeo,plato,fovl,xx)
      call fovlp(1,nbaso,p_ntab,p_iax,plato,posl,ipsl,alato,
     .  rmaxl,zl,6d0,1d0,.75d0,.5d0,xx,fovl,novl)

      if (fovl == 0) then
        print *, 'ovmin: no spheres overlap:'
        call poppr
C        call fovlp(1,nbaso,p_ntab,p_iax,plato,posl,ipsl,
C     .    alato,rmaxl,zl,6d0,1d0,.75d0,.5d0,xx,fovl,novl)
        call ovlchk(nbaso,nbaspo,posl,alato,rmaxl,0d0,' ',
     .    ipsl,modeo,plato,fovp,xx)
        if (cmdopt('--wpos=',7,0,outs))
     .    call iopos(.true.,0,outs(8:),nbaso,posl,xx)
        call rx('ovmin: no spheres overlap')
      endif

C ... Gradient of fovl wrt pos
      do  i = 1, nlst
      j = ilst(i)
      call fovlp(j,j,p_ntab,p_iax,plato,posl,ipsl,
     .    alato,rmaxl,zl,6d0,1d0,.75d0,.5d0,xx,fov0,novl0)
      do  ix = 1, 3
        val = p(3*i-3+ix)
        call dvset(posl,3*j-3+ix,3*j-3+ix,val+dx)
C        call ovlchk(nbaso,nbaspo,posl,alato,rmaxl,0d0,' ',
C     .    ipsl,modeo,plato,fovp,xx)
        call fovlp(j,j,p_ntab,p_iax,plato,posl,ipsl,
     .    alato,rmaxl,zl,6d0,1d0,.75d0,.5d0,xx,fovp,novlp)
        call dvset(posl,3*j-3+ix,3*j-3+ix,val-dx)
C        call ovlchk(nbaso,nbaspo,posl,alato,rmaxl,0d0,' ',
C     .    ipsl,modeo,plato,fovm,xx)
        call fovlp(j,j,p_ntab,p_iax,plato,posl,ipsl,
     .    alato,rmaxl,zl,6d0,1d0,.75d0,.5d0,xx,fovm,novlm)
        call dvset(posl,3*j-3+ix,3*j-3+ix,val)
        fovp = fovl + 2*(fovp-fov0)
        fovm = fovl + 2*(fovm-fov0)
        p(n+3*i-3+ix) = dlog(fovp/fovm)/2/dx
*       print *, '... i,j,ix=',i,j,ix,fovp,fovm,p(n+3*i-3+ix)
        enddo
      enddo
      ovcall = ddot(n,p(n+1),1,p(2*n+1),1)
      if (ipr >= 50) then
        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
     .    'glf.x %;4g',' ',80,lgunit(1),x,fovl/novl,dlog(fovl/novl),
     .    dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
     .    'glf.x %;4g',' ',80,lgunit(2),x,fovl/novl,dlog(fovl/novl),
     .    dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
        do  i = 1, nbaspo
          call dpscop(posl,pos,3,3*i-2,1,1d0)
          write (lgunit(2),1) pos
    1     format(3F12.6)
        enddo
        call query('continue',-1,0)
      endif

C      call prmx('grad fovl',p(1+n),n,n,1)
C      call prmx('pos now',posl,3,3,nbaspo)

C ... restore pos
      call dpcopy(posb,posl,1,3*nbaspo,1d0)
      deallocate(posb)
      call poppr

      end
