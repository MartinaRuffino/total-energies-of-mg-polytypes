      subroutine wronkj(e1,e2,r,lmax,fkk,fkj,fjk,fjj)
C- Wronskians for spherical hankels and bessels on one sphere.
C ----------------------------------------------------------------------
Ci Inputs
Ci   e1    :Wronskian for k or j (e1) with k or j (e2)
Ci         :where k = Hankel for e<0 and Neumann for e>0
Ci   e2    :Wronskian for k or j (e1) with k or j (e2)
Ci   r     :radius
Ci   lmax  :maximum l for a given site
Co Outputs
Co   fkk   :W(k,k)
Co   fkj   :W(k,j)
Co   fjk   :W(j,k)
Co   fjj   :W(j,j)
Cr Remarks
Cr  fxy is continuous as e2->e1.
Cr  fkk,fjj are symmetric in e1,e2.  fkj(e1,e2)=fjk(e2,e1).
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax
      double precision e1,e2,r
      real(8) :: fkk(*),fkj(*),fjk(*),fjj(*),ak1(100),aj1(100),
     .   ak2(100),aj2(100),dk2(100),dj2(100),dk1(100),dj1(100)
C ... Local parameters
      integer l,lp1
      real(8) :: r3,rk,rj,efac
      real(8),parameter :: tol=1d-8

C --- Case e1=e2=0 ---
      if (dabs(e1) <= tol .and. dabs(e2) <= tol) then
        r3 = r*r*r
        rk = -1d0
        rj = 1d0/r
        do  l = 0, lmax
          lp1 = l+1
          rk = rk*(2*l-1)/r
          rj = rj*r/(2*l+1)
          fkk(lp1) = rk*rk*r3/(2*l-1)
          fjj(lp1) = -rj*rj*r3/(2*l+3)
          fkj(lp1) = -0.5d0*rj*rk*r3
          fjk(lp1) = fkj(lp1)
        enddo

C --- Case e1 /= e2 ---
      elseif (dabs(e1-e2) > tol) then
        efac = 1d0/(e2-e1)
        call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,0)
        call radkj(e2,r,lmax,ak2,aj2,dk2,dj2,0)
        do  lp1 = 1, lmax+1
          fkk(lp1) = efac*r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
          fjj(lp1) = efac*r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
          fkj(lp1) = efac*r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))-efac
          fjk(lp1) = efac*r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))+efac
        enddo

C --- Case e1 == e2 but not zero ---
      else
        call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,0)
        call radkj(e1,r,lmax,ak2,aj2,dk2,dj2,1)
        do  lp1 = 1, lmax+1
          fkk(lp1) = r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
          fjj(lp1) = r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
          fkj(lp1) = r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))
          fjk(lp1) = r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))
        enddo
      endif

      end

