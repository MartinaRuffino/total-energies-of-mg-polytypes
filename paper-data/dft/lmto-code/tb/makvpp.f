      subroutine makvpp(alat,d,v0,lsec,lscale,ppmode,poly,cutmod,cutpp,
     .   erep,derep,dderep,prep)
C- Makes pair repulsive interaction
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat  : lattice constant (a.u.)
Ci   d     : distance at which PP is sought (in a.u.)
Ci   v0    : pair potential parameters (see vppder.f)
Ci   lsec  : return second derivative of pair interation
Ci   lscale: scale (lscale=.true.)/do not scale (lscale=.false.)
Ci           cutoff distances with alat
Ci   ppmode: par potential mode (see vppder.f)
Ci   poly  : the order of the cutoff polynomial (currently 4 or 5)
Ci   cutmod: cutoff mode
Ci           0 no cutoff
Ci           1 augmentative cutoff: v0 is augmented with the cutoff polynomial
Ci           2 multiplicative cutoff: v0 is multiplied by the cutoff polynomial
Ci           (see poly45.f for more comments)
Ci   cutpp : cutoff distances r1 = cutpp(1) and rc = cutpp(2), 0 =< r1 =< rc
Ci           r1 and rc are in units of alat if lscale = .true.  or in a.u.
Ci           if lscale = .false.
Co Outputs
Co   erep,derep,prep: repulsive interaction, its derivative and pressure
Co   dderep, second derivative of erep
Cr Remarks
Cr   * makvpp calls vppder (which actually evaluates pair potential) and
Cr     then applies the specified cutoff. This is similar to how makvme works.
Cr
Cu Updates
Cu   19 Apr 11 (SL) first created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical, intent(in) ::  lsec,lscale
      integer, intent(in) ::  ppmode,poly,cutmod
      double precision, intent(in)  :: alat,d,v0(3,3),cutpp(2)
      double precision, intent(out) :: erep,derep,dderep,prep
C Local paramters
      double precision r1,r2,val,slo,curv,e,de,dde

!       call tcn('makvpp')

      if (cutmod /= 0) then
C ... Set cuttof distances
        r1 = cutpp(1)
        r2 = cutpp(2)
        if (lscale) then
          r1 = r1 * alat
          r2 = r2 * alat
        endif
      endif

      if (cutmod == 0) then
C ... case no cutoff
        call vppder(lsec,ppmode,d,v0,erep,derep,dderep)
      elseif (d <= r1) then
C ... case cutoff /= 0 but d < r1 => same as above
        call vppder(lsec,ppmode,d,v0,erep,derep,dderep)
      elseif (d >= r2) then
C ... case cutoff /= 0 and d > r2 => set output to zero
        erep = 0d0
        derep = 0d0
        dderep = 0d0
      else
C ... case cutoff /= 0 and r1 < d < r2
        if (cutmod == 1) then
C ...     augmentative cutoff
          call vppder(.true.,ppmode,r1,v0,val,slo,curv)
          call pcut45(poly,d,r1,r2,val,slo,curv,erep,derep,dderep)
        elseif (cutmod == 2) then
C ...     multiplicative cutoff
          call vppder(lsec,ppmode,d,v0,val,slo,curv)
          call pcut45(poly,d,r1,r2,1d0,0d0,0d0,e,de,dde)
          erep = val * e
          derep = slo * e + val * de
          if (lsec)
     .    dderep = val * dde + 2d0 * slo * de + curv * e
        else
          call rxi(' makvpp: cutoff mode cutmod = %i not implemented',
     .      cutmod)
        endif
      endif
      prep = derep*d

!       call tcx('makvpp')
      end
