      subroutine mkppar(dclabl,nl,nsp,nc,lmx,z,rmax,avw,amsh,nrmsh,
     .  pnu,idmod,ves,qnu,pp,a)
C- Make pp's for all atoms for which potential is available.
C ----------------------------------------------------------------------
Ci Inputs
Ci   dclabl:class name, packed as a real number
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nc    :number of classes
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.,
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   amsh  :radial mesh points are given by rofi(i) = b [e^(amsh(i-1)) -1]
Ci   nrmsh :number of radial points
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
Ci   ves   :es potential at MT boundary
Ci   qnu   :energy-weighted moments of the sphere charges
Co Outputs
Co   pp    :potential parameters, with b.c. ves(ic)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   04 May 11 modify call to potpar
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nsp,lmx(*),nrmsh(*),nc,idmod(nl,nc)
      double precision z(1),rmax(*),avw,amsh(*),pnu(nl,nsp,1),
     .  qnu(3,nl,nsp,nc),ves(nc),pp(6,nl,nsp,nc),dclabl(1)
C ... Dynamically allocated local arrays
      real(8), allocatable :: rofi(:)
      real(8), allocatable :: v(:)
      real(8), allocatable :: g(:)
      real(8), allocatable :: gp(:)
C ... Local parameters
      logical aiopot,sw,aiopar
      integer nrmx,nr,ic,ir,ifi,fopna,iprint,i1mach,n0
      integer idu(5)
C     integer lrel,lgrad
      parameter (nrmx=5001,n0=10)
      double precision a,b,ea,thrpv,exc(2),thrpvl(10),rpb,rmx,t(n0,2)
      double precision ekap,sop,xx,GFr(1)
      logical lso
      character*8 clabl

      call iinit(idu,5)
      lso = .false.
      ekap = 0d0

      do  ic = 1, nc
        call r8tos8(dclabl(ic),clabl)
        rmx = rmax(ic)
        nr = nrmsh(ic)
        a  = amsh(ic)
        call rx('lrel,lgrad not set')
        allocate(rofi(nr),v(nr*nsp))
        ifi = fopna(clabl,30,0)
        if (aiopot(nr,nsp,a,rmx,-99d0,v,ifi)) then
          ea = dexp(a)
          b = rmx/(dexp(a*(nr-1)) - 1d0)
          rpb = b
          do  ir = 1, nr
          call dvset(rofi,ir,ir,rpb-b)
            rpb = rpb*ea
          enddo
          allocate(g(nr*2),gp(nr*2*4))
          if (iprint() >= 60) then
            print *, 'mkppar: potential parms before calling potpar'
            sw = aiopar(clabl,0,pp(1,1,1,ic),xx,ves(ic),nl,lmx(ic),nsp,
     .        -i1mach(2))
          endif
          call rx('mkppar needs mpolp')
          call potpar(nl,nsp,lmx(ic),z(ic),rmx,avw,ekap,lso,.false.,
     .      .false.,a,nr,rofi,v,pnu(1,1,ic),idmod(1,ic),
     .      exc,qnu(1,1,1,ic),idu,xx,xx,1d0,thrpv,thrpvl,g,gp,
     .      pp(1,1,1,ic),xx,sop,xx,xx,t,GFr)

C     ... Shift enu and c by ves
          call daxpy(nl*nsp,1d0,ves(ic),0,pp(1,1,1,ic),6)
          call daxpy(nl*nsp,1d0,ves(ic),0,pp(2,1,1,ic),6)

          if (iprint() >= 60) then
            print *, 'mkppar: potential parms after calling potpar'
            sw = aiopar(clabl,0,pp(1,1,1,ic),xx,ves(ic),nl,lmx(ic),nsp,
     .        -i1mach(2))
          endif
          deallocate(g,gp)
        else
          if (iprint() > 20) print *,
     .      'mkppar: missing potential for class ',clabl
        endif
        call fclose(ifi)
        deallocate(rofi,v)
      enddo
      end
