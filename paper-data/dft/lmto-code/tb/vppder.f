      subroutine vppder(lsec,ppmode,d,v0,erep,derep,dderep)
C- Evaluate pair potential and its derivatives at given distance
C ----------------------------------------------------------------------
Ci Inputs
Ci   lsec   : if lsec = .true. make the second derivative
Ci            otherwise make only the value and the first derivative
Ci   ppmode : the switch to select between particular forms of pair potential
Ci            (see Remarks)
Ci   d      : distance at which the pair potential is evaluated
Ci   v0     : pair potential parameters (see Remarks)
Co Outputs
Co   erep,derep,dderep: repulsive interaction and its first and second
Co                      derivatives
Co                      if lsec = .false. dderep is not referenced
Cr Remarks
Cr   * v0 is a vector length 9 allowing the pairwise energy to be
Cr     parameterised as specified by ppmode. v0 is dimensioned like this:
Cr         ! (1,1) (2,1) (3,1)  (1,2) (2,2) (3,2)  (1,3) (2,3) (3,3)
Cr
Cr   * The following set of ppmode switches is currently implemented:
Cr     ppmode =  0   no pair  potential
Cr       V0(r) = 0
Cr
Cr     ppmode = 10  sum of power law times exponentials, up to 3 terms
Cr       V0(r) = sum_i A_i r^B_i exp(-C_i r)
Cr       input:
Cr         !  A_1   B_1   C_1    A_2   B_2   C_2    A_3   B_3   C_3
Cr       restrictions:
Cr         (1) A_1 >= 0
Cr         (2) B_i <= 0, i = 1, 3
Cr         (3) C_i >= 0, i = 1, 3
Cr         (4) A_2 = 0 signals that the i = 3 term is zero as well
Cr            this is done for protection against parameter-switches
Cr            used previosly (see makv0.f)
Cr
Cr     ppmode = 20   Quadratic "Jim Chadi" potential
Cr       V0(r) = A_1 eps + A_2 eps^2, eps = (r - r0)/r0
Cr       input (dash means ignored positions):
Cr         !  A_1    -    r0     A_2    -     -      -     -     -
Cr       restrictions:
Cr         (1) requires cutoff (see makvme.f)
Cr
Cr     ppmode = 21  sum of a polynomial defined in 0 =< r =< rc] and two power law functions
Cr       V0(r) = A_1 r^B_1 + A_2 r^B_2 + P_n(r/rc - 1)*(r/rc - 1)^3*\theta(r < rc), n =< 3
Cr               where P_n(x) = p0 * (1 + p1*x + p2*x^2 + p3*x^3)
Cr               and \theta(r < rc) is the step function: \theta(x < x0) = 1, if x < x0
Cr                                                                       = 0, if x > x0
Cr     (the job of the polynomial is to adjust V0(r) in the range [0, rc], whereas at rc
Cr      its value and first two derivatives turn to zero)
Cr
Cr       input:
Cr         !  A_1   B_1   A_2    B_2   r0    p0     p1    p2    p3
Cr       restrictions:
Cr         (1) at least one of A_1 or A_2 >= 0
Cr         (2) B_i <= 0, i = 1, 2
Cr         (3) rc  >  0
Cr
Cr     ppmode = 30   Goodwin-Skinner-Pettifor potential (GSP)
Cr       V0(r) = A (r0/r)^m exp[m (-{r/rc}^mc + {r0/rc}^mc)]
Cr       input (dash means ignored positions):
Cr         !   A     -     -      m    mc    r0     rc     -     -
Cr       restrictions:
Cr         (1) A, m, mc, r0  >= 0
Cr         (2) rc  > 0
Cr         (3) rc  >= r0
Cr
Cr     ppmode = 31   GSP + exponential
Cr       V0(r) = A (r0/r)^m exp[m (-{r/rc}^mc + {r0/rc}^mc)] + A_2 exp(-C_2 r)
Cr       input (dash means ignored positions):
Cr         !   A     -     -      m    mc    r0     rc    A_2   C_2
Cr       restrictions:
Cr         (1)-(3) as in ppmode = 30
Cr         (4) C_2 >= 0
Cr
Cr     ppmode = 32   GSP + power law decay
Cr       V0(r) = A (r0/r)^m exp[m (-{r/rc}^mc + {r0/rc}^mc)] + A_2 r^B_2
Cr       input (dash means ignored positions):
Cr         !   A     -     -      m    mc    r0     rc    A_2   B_2
Cr       restrictions:
Cr         (1)-(3) as in ppmode = 30
Cr         (4) B_2 <= 0
Cr
Cr     ppmode = 33   GSP + (a1*exp(-2r) + a2/r^6)*left_cutoff
Cr       V0(r) = A (r0/r)^m exp[m (-{r/rc}^mc + {r0/rc}^mc)] + P5(rm1,rm2)*(A1*exp(-2r) + A2/r^6)
Cr       input (dash means ignored positions):
Cr         !   A     rm1   rm2    m    mc    r0     rc    A_1   A_2
Cr       restrictions:
Cr         (1)-(3) as in ppmode = 30
Cr         (4) rm2 < rm1
Cr         (5) A2 is usually negative but this is not compulsary
Cr       Comment: ppmode = 33 serves to imitate a Force Field separation of the PP into the
Cr                intra molecular (short range) part, represented by GSP in this case, and
Cr                inter molecular (long range) part, represented by the term in brackets,
Cr                l(r) = A1*exp(-2r) + A2/r^6.
Cr                The LR part is then smoothly cut off to zero on the left so that it doesn't
Cr                interfere with the SR part. This is done using the multiplicative cutoff
Cr                polynomial P5 (see p45.f)
Cr
Cr   * compatibility with makv0:
Cr       (1) positions of the parameters are preserved for compatibility
Cr           with existing ctrl files
Cr       (2) parameter-switches are removed. The type of pair potential,
Cr           type of cutoff, and cutoff distances should be specified in
Cr           ppmode, cutmod, and cutpp, respectively (see makvpp.f)
Cr       (3) a few types of PP were abondoned, such as shifted exponentials
Cr           or those with augmented hard core
Cr
Cr   * Note that r0 and rc are input in a.u., not in units of alat
Cr
Cu Updates
Cu   19 Apr 11 (SL)  first created from makv0.f
Cu   30 Jan 12 (SL)  option 21 added
Cu   24 Jul 13 (SL)  option 33 added
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical, intent(in) ::  lsec
      integer, intent(in) ::  ppmode
      double precision, intent(in)  :: d,v0(3,3)
      double precision, intent(out) :: erep,derep,dderep
C Local paramters
      integer iexp,nexp
      double precision, parameter :: tol=1d-8
      double precision e,de,dde,dle,eps,m,mc,r0,rc
      double precision ep,dep,ddep
      double precision val,slo,curv,rm1,rm2,v1,v2
      double precision frac,d1
      double precision a1,a2,b1,b2,rc1,p0,p1,p2,p3,
     .                 pol,dpol,ddpol,e1,e2,de1,de2,x,xx,xxx,iv
c----------
      integer mode0,mode1

C --- exponentials, GSP or quadratic?
      mode1 = mod(ppmode/10,10)
      mode0 = mod(ppmode,10)

      if (mode1 == 3) then
C --- GSP potential ---
          m  = v0(1,2)
          mc = v0(2,2)
          r0 = v0(3,2)
          rc = v0(1,3)
          e = v0(1,1)*exp(m*(log(r0)+(r0/rc)**mc))
          frac = (d/rc)**mc
          d1 = 1d0/d
          erep = e*exp(m*(log(d1)-frac))
          dle = -m*d1*(1d0 + mc*frac)      ! dle = erep'/erep
          derep = dle * erep
          if (lsec) dderep =
c    .      (dle + d1*(1d0 + mc*mc*frac/(1d0 + frac)))*derep
     .      erep * (dle*dle + m*d1*d1*(1d0 - mc*(mc-1d0)*frac))
C ...     add a power (mode0 == 1) or exponential (mode0 == 2) to GSP ---
          if (mode0 /= 0) then
C           *** mode 31 ***
            if (mode0 == 1) then
              e = v0(2,3)*exp(-v0(3,3)*d)
              de = -e*v0(3,3)
              if (lsec) dde = -de*v0(3,3)
C           *** mode 32 ***
            elseif (mode0 == 2) then
              e = v0(2,3)*d**v0(3,3)
              de = e*v0(3,3)*d1
              if (lsec) dde = de*(v0(3,3)-1d0)*d1
C           *** mode 33 ***
            elseif (mode0 == 3) then
              rm1 = v0(2,1)
              rm2 = v0(3,1)
              A1 = v0(2,3)
              A2 = v0(3,3)
C             *** Case d =< rm2 ***
              if (d <= rm2) then
                e = 0d0
                de = 0d0
                if (lsec) dde = 0d0
              else
                v1 = A1*exp(-2d0*d)
C               v2 = A2/d**6
                v2 = A2*d1**6
                val = v1 + v2
                slo = -2d0*v1 - 6d0*v2*d1
                if (lsec) curv = 4d0*v1 + 42d0*v2*d1*d1
C               *** Case rm2 < d =< rm1 ***
                if (d <= rm1) then
                  call pcut45(5,d,rm1,rm2,1d0,0d0,0d0,ep,dep,ddep)
                  e = val*ep
                  de = slo*ep + val*dep
                  if (lsec) dde = val*ddep + 2d0*slo*dep + curv*ep
C               *** Case d > rm1 ***
                else
                  e = val
                  de = slo
                  if (lsec) dde = curv
                endif
              endif
            else
              call rxi(' makvpp: ppmode = %i not implemented',ppmode)
            endif
            erep = erep + e
            derep = derep + de
            if (lsec) dderep = dderep + dde
          endif

      elseif (mode1 == 2) then
        if (mode0 == 0) then
C ---   Quadratic "Jim Chadi" potential ----
          iv = 1.0d0/v0(3,1)
          eps   = d*iv - 1d0
          erep  = (v0(1,1) + v0(1,2)*eps)*eps
          derep = (v0(1,1) + 2d0*v0(1,2)*eps)*iv
          if (lsec) dderep = 2d0*v0(1,2)*iv*iv
        elseif (mode0 == 1) then
C ---   Polynomial + two power law functions ----
          a1 = v0(1,1)
          b1 = v0(2,1)
          a2 = v0(3,1)
          b2 = v0(1,2)

          rc = v0(2,2)
          p0 = v0(3,2)
          p1 = v0(1,3)
          p2 = v0(2,3)
          p3 = v0(3,3)

C ...     Contribution from the power law functions
          d1 = 1d0/d
          e1 = a1*d1**(-b1)
          e2 = a2*d1**(-b2)
          de1 = e1*b1*d1
          de2 = e2*b2*d1
          erep = e1 + e2
          derep = de1 + de2
          if (lsec) dderep = dderep + (de1*(b1-1.d0) + de2*(b2-1.d0))*d1

C ...     Contribution from the polynomial
          if (d < rc .and. abs(p0) > tol) then
            rc1 = 1.d0/rc
            x = d*rc1-1d0
            xx = x*x
            xxx = xx*x
C ...       f(r) = P(x) x^3, x = r/rc - 1
C           rc (d/dr) f(r) = P'(x) x^3 + 3 P(x) x^2
C           rc^2 (d^2/dr^2) f(r) = P"(x) x^3 + 6 P'(x) x^2 + 6 P(x) x
            pol  = p0 * (1d0 + x * (p1 + x * (p2 + x*p3)))
            dpol = p0 * (p1 + x * (2d0*p2 + x*3d0*p3))

            erep = erep + pol*xxx
            derep = derep + (dpol*x + 3d0*pol)*xx*rc1
            if (lsec) then
              ddpol = 2d0*p0 * (p2 + x*3d0*p3)
              dderep = dderep +
     .                (ddpol*xx + 6d0*(dpol*x + pol))*x*rc1*rc1
            endif
          endif
        else
          call rxi(' makvpp: ppmode = %i not implemented',ppmode)
        endif
      elseif (mode1 == 1) then
C --- Sum of power law times exponentials ----
        erep = 0d0
        derep = 0d0
        if (lsec) dderep = 0d0

C ...   how many terms are there?
        if (abs(v0(1,2)) <= tol) then
          nexp = 1
        elseif (abs(v0(1,3)) <= tol
     .    .or. abs(v0(1,3)+1d0) <= tol) then
C ...   need to protect against parameters-switches
          nexp = 2
        else
          nexp = 3
        endif

        d1 = 1d0/d
        do   iexp = 1, nexp
          e = v0(1,iexp)*exp(v0(2,iexp)*log(d)-v0(3,iexp)*d)
          erep = erep + e
          dle = v0(2,iexp)*d1 - v0(3,iexp)
          de = dle * e
          derep = derep + de
          if (lsec) dderep = dderep - v0(2,iexp)*e*d1*d1 + dle*de
        enddo
      elseif (mode1 == 0) then
C --- No pair interactions ----
        erep = 0d0
        derep = 0d0
        if (lsec) dderep = 0d0
C --- end of if-loop over potential modes ---
      else
        call rxi(' makvpp: ppmode = %i not implemented',ppmode)
      endif

      end
