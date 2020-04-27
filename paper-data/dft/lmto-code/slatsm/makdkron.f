      subroutine makdkron(opt,m,n,p,q,a,b,kron)
C- Make Kroneker product of a pair of matrices
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :0 a and b are both given
Ci         :1 a is the unit matrix
Ci         :2 b is the unit matrix
Ci         :10s digit
Ci         :0 Copy A \otimes B into kron
Ci         :1 Add  A \otimes B into kron
Ci   m     :first dimension of a
Ci   n     :2nd dimension of a
Ci   p     :first dimension of b
Ci   q     :2nd dimension of b
Ci   a     :left matrix entering into the Kroneker product
Ci   b     :right matrix entering into the Kroneker product
Co Outputs
Co   kron  :Kroneker product
Cr Remarks
Cr   kron = (A\otimes B)_{p(r-1)+v, q(s-1)+w} = a_{rs} b_{vw}
Cu Updates
Cu   30 Mar 16  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,m,n,p,q
      double precision a(m,n), b(p,q), kron(m*p,n*q)
C ... Local parameters
      integer r,s,v,w
      logical unita, unitb
      double precision ars, bvw

      unita = mod(mod(opt,10),2) /= 0
      unitb = mod(opt,10) >= 2
      if (mod(opt/10,10) == 0) call dpzero(kron,size(kron))

      do  r = 1, m
      do  s = 1, n

        if (.not. unita) then
          ars = a(r,s)
        else
          ars = 0 ; if (r == s) ars = 1
        endif

        do  v = 1, p
        do  w = 1, q

          if (.not. unitb) then
            bvw = b(v,w)
          else
            bvw = 0; if (v == w) bvw = 1
          endif

          kron(p*(r-1)+v,q*(s-1)+w) = kron(p*(r-1)+v,q*(s-1)+w) + ars * bvw
C         kron(p*(r-1)+v,q*(s-1)+w) = 100*r+10*s+v+dble(w)/10
C         print 333, r,s,v,w, p*(r-1)+v,q*(s-1)+w, kron(p*(r-1)+v,q*(s-1)+w)
  333     format(4i3,2x,2i3,f12.6)

        enddo
        enddo

      enddo
      enddo

      end

