      subroutine makv0(alat,d,v0,erep,derep,p)
C- Makes pair repulsive interaction
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat,v0,d (au)
Co Outputs
Co   erep: repulsive interaction and its derivative and pressure
Cr Remarks
Cr   V0 is a vector length 9 allowing the pairwise energy to be a sum of
Cr   three terms of the form A r^B exp(-Cr).
Cr
Cr   If the second power exponent is positive it is interpreted as a
Cr   shift in origin to -rs of the second exponential so the first two
Cr   functions take the form A1 r^B1 exp(-C1r) + A2 exp(-C2(r+rs))
Cr   the third function is either A3 r^B3 exp(-C3r) or the last three
Cr   numbers indicate a smooth cut off: see below.
Cr
Cr   If the first power exponent is positive, then the pair potential
Cr   is of Chadi's form, namely A1 eps + A2 eps^2: the third number
Cr   in the 1st set is the equilibrium bond length.
Cr
Cr   If the first power exponent is positive and the third number
Cr   of the set (C) is negative then the Goodwin-Skinner-Pettifor
Cr   form is used: A (r0/d)^m exp[m (-{d/rc}^mc + {r0/rc}^mc)].
Cr   For example
Cr
Cr       ! A 1 -1 m mc r0 rc 0 0
Cr
Cr   In the GSP there is a remaining two parameters. These are used
Cr   to add an additional power or exponential depending whether
Cr   V0(3,3) is positive or negative. To enable this v0(2,3) must be
Cr   positive. If it's negative then the GSP is augmented by a
Cr   hard core (see below)
Cr
Cr   In the case of the quadratic potential, this may be augmented with
Cr   a polynomial to cut it off smoothly to zero between distances
Cr   r1 and rc. These are passed as the third members of the second
Cr   and third sets. For example
Cr
Cr       !  A1 1 r0 A2 0 r1 rc nh rcore
Cr
Cr   The entry rcore allows a hard core radius so be set for quadratic
Cr   and GSP potentials. For r < rcore, the pair potential becomes
Cr   Aexp(qr)/r^nh and A and q are determined by the matching of value
Cr   and slope at rcore. To augment the GSP with a hard core "head"
Cr   use, for example,
Cr
Cr       ! A 1 -1 m mc r0 rc nh rcore
Cr
Cr   ensuring that nh < 0. In this case no power or exponential is
Cr   added to the GSP potential.
Cr
Cr   The pair potential may also be of the form A r^B exp(-Cr) cut off
Cr   smoothly between r1 and rc. In this case the first six elements
Cr   are A>0, B<0 and C>0, the 7th element is -1 and the next two
Cr   are r1 and rc. Thus
Cr
Cr       ! A1 B2 C3 A2 B2 C2 -1 r1 rc
Cr
Cr   Note that r1 and rc are input in units of alat and multiplied here
Cr
Cr   v0 is dimensioned like this
Cr       ! (1,1) (2,1) (3,1)  (1,2) (2,2) (3,2)  (1,3) (2,3) (3,3)
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer i
      double precision alat,d,v0(3,3),erep,derep,p
C Local paramters
      double precision e,de,eps,m,mc,r0,rc,c(0:5),pl,ppl,vsc(2,3),
     .                 val,slo,curv,r1,rcore,qh,Ah,nh,epsc,erepc,derepc
      logical head,tail
C     integer i1mach

C 5th order polynomial and its derivatives
      pl(d) = c(0)+c(1)*d+c(2)*d**2+c(3)*d**3+c(4)*d**4+c(5)*d**5
      ppl(d) = c(1)+2d0*c(2)*d+3d0*c(3)*d**2+4d0*c(4)*d**3+5d0*c(5)*d**4

      head = .false.
      tail = .false.
      erep = 0
      derep = 0d0

      if (v0(2,1) > 0) then
C --- GSP or quadratic ---
        if (v0(3,1) < 0) then
C --- GSP potential ---
          m  = v0(1,2)
          mc = v0(2,2)
          r0 = v0(3,2)
          rc = v0(1,3)
C --- augment with head if r < r_core ---
          if (v0(2,3) < -1d-4 .and. d < v0(3,3)) then
            nh = v0(2,3)
            rcore = v0(3,3)
            e = v0(1,1)*r0**m*exp(m*(r0/rc)**mc)
            erepc = e*exp(-m*(rcore/rc)**mc)/rcore**m
            derepc = -m*erep*(1d0 + mc*(rcore/rc)**mc)/rcore
            val = erepc
            slo = derepc
            qh = - (slo/val + nh/rcore)
            Ah = val*rcore**nh*exp(qh*rcore)
            erep = Ah*exp(-qh*d)/(d**nh)
            derep = - erep * (qh + nh/d)
          else
            e = v0(1,1)*r0**m*exp(m*(r0/rc)**mc)
            erep = e*exp(-m*(d/rc)**mc)/d**m
            derep = -m*erep*(1d0 + mc*(d/rc)**mc)/d
          endif
C --- add a power or exponential to GSP ---
          if (v0(2,3) > 1d-4) then
            if (v0(3,3) > 0) then
              e = v0(2,3)*dexp(-v0(3,3)*d)
              de = -e*v0(3,3)
            else
              e = v0(2,3)*d**v0(3,3)
              de = e*v0(3,3)/d
            endif
            erep = erep + e
            derep = derep + de
          endif
        else
C --- quadratic "Jim Chadi" potential ---
          if (v0(3,3) > 1d-4) then
            head = .true.
            rcore = v0(3,3)
            nh = v0(2,3)
          endif
          if (v0(3,2) > 1d-4) then
            tail = .true.
            r1 = v0(3,2)
            rc = v0(1,3)
          endif
          eps   = (d - v0(3,1))/v0(3,1)
          erep  = v0(1,1)*eps + v0(1,2)*eps**2
          derep = (v0(1,1) + 2d0*v0(1,2)*eps)/v0(3,1)
C --- augment with heads and tails --
          if (tail .and. d > r1) then
            if (d > rc) then
              erep = 0d0
              derep = 0d0
            else
C --- tail ---
              epsc   = (r1 - v0(3,1))/v0(3,1)
              erepc  = v0(1,1)*epsc + v0(1,2)*epsc**2
              derepc = (v0(1,1) + 2d0*v0(1,2)*epsc)/v0(3,1)
              val = erepc
              slo = derepc
              curv = 2*v0(1,2)/v0(3,1)**2
              call rcut(val,slo,curv,r1,rc,c)
              erep = pl(d)
              derep = ppl(d)
            endif
          endif
C --- head ---
          if (head .and. d < rcore) then
            epsc   = (rcore - v0(3,1))/v0(3,1)
            erepc  = v0(1,1)*epsc + v0(1,2)*epsc**2
            derepc = (v0(1,1) + 2d0*v0(1,2)*epsc)/v0(3,1)
            val = erepc
            slo = derepc
            qh = - (slo/val + nh/rcore)
            Ah = val*rcore**nh*exp(qh*rcore)
            erep = Ah*exp(-qh*d)/(d**nh)
            derep = - erep * (qh + nh/d)
          endif
        endif
      else
C --- sum of power laws times exponentials ---
        if (abs(v0(1,3)+1) < 1d-8) then
          r1 = v0(2,3)*alat
          rc = v0(3,3)*alat
          if (d < r1) then
            if (v0(2,2) <= 0d0) then
              do   i = 1, 2
                e = v0(1,i)*d**v0(2,i)*dexp(-v0(3,i)*d)
                erep = erep + e
                derep = derep + e*(v0(2,i)/d - v0(3,i))
              enddo
            else
              e = v0(1,1)*d**v0(2,1)*dexp(-v0(3,1)*d)
     .          + v0(1,2)*dexp(-v0(3,2)*(d+v0(2,2)))
              erep = e
              derep = v0(1,1)*d**v0(2,1)*dexp(-v0(3,1)*d)
     .              * (v0(2,1)/d - v0(3,1))
     .              - v0(3,2)*v0(1,2)*dexp(-v0(3,2)*(d+v0(2,2)))
            endif
          elseif (d > rc) then
            erep = 0d0
            derep = 0d0
          else
            val  = 0d0
            slo  = 0d0
            curv = 0d0
            if (v0(2,2) <= 0d0) then
              do   i = 1, 2
                vsc(i,1) = v0(1,i)*r1**v0(2,i)*dexp(-v0(3,i)*r1)
                vsc(i,2) = vsc(i,1)*(v0(2,i)/r1 - v0(3,i))
                vsc(i,3) = vsc(i,2)*(v0(2,i)/r1 - v0(3,i))
     .                   - vsc(i,1)*v0(2,i)/r1**2
              enddo
            else
                vsc(1,1) = v0(1,1)*r1**v0(2,1)*dexp(-v0(3,1)*r1)
                vsc(1,2) = vsc(1,1)*(v0(2,1)/r1 - v0(3,1))
                vsc(1,3) = vsc(1,2)*(v0(2,1)/r1 - v0(3,1))
     .                   - vsc(1,1)*v0(2,1)/r1**2
                vsc(2,1) = v0(1,2) * exp(-v0(3,2)*(r1 + v0(2,2)))
                vsc(2,2) = -v0(3,2) * vsc(2,1)
                vsc(2,3) = v0(3,2)**2 * vsc(2,1)
            endif
            val  = vsc(1,1) + vsc(2,1)
            slo  = vsc(1,2) + vsc(2,2)
            curv = vsc(1,3) + vsc(2,3)
            call rcut(val,slo,curv,r1,rc,c)
            erep = pl(d)
            derep = ppl(d)
          endif
        else
          if (v0(2,2) <= 0d0) then
            do   i = 1, 3
              e = v0(1,i)*d**v0(2,i)*dexp(-v0(3,i)*d)
              erep = erep + e
              derep = derep + e*(v0(2,i)/d - v0(3,i))
            enddo
          else
            do   i = 1, 3, 2
              e = v0(1,i)*d**v0(2,i)*dexp(-v0(3,i)*d)
              erep = erep + e
              derep = derep + e*(v0(2,i)/d - v0(3,i))
            enddo
            e = v0(1,2) * exp(-v0(3,2)*(d + v0(2,2)))
            erep = erep + e
            derep = derep - e*v0(3,2)
          endif
        endif
      endif
      p = derep*d
      end
