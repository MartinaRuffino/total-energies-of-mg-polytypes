      subroutine hansz(rsml,nrx,lmax,nxi,lxi,exi,nr,r2,yl,shan,sbes)
C- Solid smoothed Hankels for a vector of points, any energy.
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsml  : vector of l-dependent smoothing radii of Hankels
Ci         : EITHER must be specified for lmin..lmax
Ci         : OR     rsmh(0) = const < 0. Implies rsmh(l) = -const for all l
Ci   nrx   :leading dimension of yl, xi, and h
Ci   nxi   :number of trial energies
Ci exi,lxi :trial energies and corresponding angular momenta
Ci   lmax  :maximum angular momentum, must be >= {lxi(ie),ie=1,nxi}
Ci   nr    :actual number of points, n must be =< nrx
Ci   r2    :r^2 = x^2 + y^2 + z^2 for the mesh points
Ci   yl    :Ylm(i,ilm), the (real) spherical harmonic polynomials
Ci  x,y,z  :Cartesian coordinates of points at which solid functions are
Ci         :to be evaluated
Co Outputs
Co   shan  :solid smoothed Hankel (e<0) or Neumann (e>0) functions, H=xi*Ylm
Ci   sbes  :if e>0, solid unsmoothed Bessel functions * e^{l+1/2)
Ci         :if e<0, not referenced
Cl Local variables
Co   xi    :e>0 radial part of smoothed Neumann function / r^l
Co         :e<0 radial part of smoothed Hankel function / r^l
Co   fi    :e>0 radial part of unsmoothed Bessel function * e^{l+1/2}/r^l
Cu Updates
Cu   13 Mar 09 (S. Lozovoi) sph. harmonics are passed from outside
Cu   04 May 07 (S. Lozovoi) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax,nrx,nxi,nr
      integer lxi(nxi)
c     double precision x(nrx),y(nrx),z(nrx)
      double precision exi(nxi),rsml(0:lmax)
      double precision r2(nrx),yl(nrx,(lmax+1)**2)
      double precision shan(nrx,(lmax+1)**2,nxi)
      double precision sbes(nrx,(lmax+1)**2,nxi)
C ... Local parameters
      integer il,im,ilm,job,ie,lmx
      logical lpos
      double precision ee
c ... work auto-arrays
c     double precision xi(nrx,0:lmax),fi(nrx,0:lmax)
c ... Allocatable work arrays
      double precision, allocatable :: xi(:,:),fi(:,:)


      call tcn('hansz')
      if (nrx < nr)
     . call rxi('hansz: increase nrx, needed at least ',nr)

c ... Allocate xi and fi (the latter is a true array if any energy > 0)
      lpos = .false.
      do  ie = 1, nxi
        if (exi(ie) > 0d0) then
          lpos = .true.
          exit
        endif
      enddo

      allocate (xi(nrx,0:lmax))
      if (lpos) then
        allocate (fi(nrx,0:lmax))
      else
        allocate (fi(1,1))
      endif

      job = 00
c Begin the cycle over energies
      do  ie = 1, nxi
        lmx = lxi(ie)
        if (lmx > lmax)
     .    call rxi('hansz: increase lmax, need at least ',lmx)
        ee = exi(ie)
        lpos = (ee > 0d0)

c Make radial parts of smooth Hankels
c         call tcn('rad. sm. hankels')
        if (lpos) then
          call hansrz(rsml,0,lmx,ee,r2,nrx,nr,job,xi,fi)
        else
          call hansrz(rsml,0,lmx,ee,r2,nrx,nr,job,xi,xi)
        endif
c         call tcx('rad. sm. hankels')

c Make smooth solid Hankels
c         call tcn('hs*Y')
        ilm = 0
        do  il = 0, lmx
           do  im = -il, il
             ilm = ilm+1
             shan(1:nr,ilm,ie) = xi(1:nr,il)*yl(1:nr,ilm)
             if (lpos)
     .         sbes(1:nr,ilm,ie) = fi(1:nr,il)*yl(1:nr,ilm)
           enddo
        enddo
c         call tcx('hs*Y')

c End the cycle over energies
      enddo

      deallocate (xi,fi)
      call tcx('hansz')
      end
