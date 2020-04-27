      subroutine pgnei(mode,e1,e2,g1,g2,v,enu,i,f)
C- f=\int_{e1-v}^{e1} dE (E-enu)**i * max(g(E),0), where
C- g(E)=g1 for mode=0 and g(E)=[e2*g1-e1*g2+E*(g2-g1)]/(e2-e1) for mode=1
C- also, f=0 if g1=0
Ci Inputs
Ci mode,e1,e2,g1,g2,v,enu,i
Co Outputs
Co f
Cu Updates
Cu    28 Jan 04 (S.Faleev) First created
      implicit none
C Passed parameters
      integer mode,i
      double precision e1,e2,g1,g2,v,enu,f
C Local variables
      double precision g0,a,b

      if(mode == 0 .or. (mode == 1 .and. abs(g1-g2) < 1d-8)) then
         f = max(g1,0d0)*((e1-enu)**(i+1) - (e1-v-enu)**(i+1))/(i+1)
      elseif (mode == 1) then
         if (abs(e2-e1) < 1d-8) call rx('pgnei: illegal energy mesh')
         a = (e2*g1-e1*g2)/(e2-e1)
         b = (g2-g1)/(e2-e1)
         g0 = a + b*(e1-v)
C        g1 < 0
         if (g1 <= 0d0) then
            f = 0d0
C        g1 > 0
         else
C           g0 > 0
            if (g0 > 0d0) then
               f=( (e1-enu)**(i+2)-(e1-v-enu)**(i+2) )/(i+2)*b +
     .           ( (e1-enu)**(i+1)-(e1-v-enu)**(i+1) )/(i+1)*(a+b*enu)
C           g0 < 0
            else
               if ( (e1-(-a/b))*(-a/b-(e1-v)) < 0d0)
     .           call rx('pgnei: wrong node')
               f=( (e1-enu)**(i+2)-(-a/b-enu)**(i+2) )/(i+2)*b +
     .           ( (e1-enu)**(i+1)-(-a/b-enu)**(i+1) )/(i+1)*(a+b*enu)
            endif
         endif
      else
         call rx('pgnei: wrong mode')
      endif

      end
