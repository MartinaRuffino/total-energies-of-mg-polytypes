      subroutine vsl0(job,rmt,eh,lmaxg,lmaxh,ndim,sf,
     . fval,fslo,flap,aj,dj)
C- lm-resolved value, slope and Laplacian of a function F_L derived from
C  its 1c expansion into Bessel functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :if 1s digit = 0, take aj and dj on input
Ci         :if 1s digit = 1, calculate Bessel functions and their radial
Ci            derivatives and return them as aj and dj
Ci         :if 1s digit = 2, calculate Hankel functions and their radial
Ci            derivatives and return them as aj and dj
Ci         :if 10s digit = 0, all radial functions (aj,dj,fval,fslo, and flap)
Ci         :are divided by r^l, otherwise not.
Ci         :if 1s digit = 0 then 10s digit is ignored; the meaning of radial
Ci         :functions is defined by input arrays aj and dj
Ci   rmt   :radius at which value, slope, and Laplacian of F_L are evaluated
Ci   eh    :energy of F_L and Bessel functions implied in the expansion
Ci   lmaxg :l-cutoff for Bessel functions
Ci   lmaxh :l-cutoff for F_L
Ci   ndim  :leading dimension of sf and output arrays fval, fslo, and flap
Ci   sf    :coefficients sf(M,L) of decomposition of head function F_L
Ci         :into Bessel functions J_M, M=1..(lmaxg+1)^2, L=1..(lmaxh+1)^2
Ci         :(if F_L := H_L then sf are unscreened structure constants)
Co Outputs
Co   fval  : lm-resolved value of F_L       |
Co   fslo  : lm-resolved slope of F_L       | at rmt (see Remarks)
Co   flap  : lm-resolved Laplacian of F_L   |
Co   aj    : radial part of Bessel functions j_l / r^l (if job's 10s digit = 0)
Co   dj    : radial part of [d/dr (j_l)] / r^l         (     -    "    -      )
Cl Local variables
Cl   ak,dk : radial Hankel (eh<0) or Neumann (e>0) functions and their radial derivative
Cl         : ak and dk are not used here
Cr Remarks
Cr   fval, fslo, and flap are defined so that:
Cr             F_L(r-R) = \sum_M fval(M,L) * r^m * Y_M(r)  | if job's 10s digit = 0
Cr      (d/dr) F_L(r-R) = \sum_M fslo(M,L) * r^m * Y_M(r)  | otherwise without r^m
Cr        \lap F_L(r-R) = \sum_M flap(M,L) * r^m * Y_M(r)  |
Cr   are approximately fulfilled at |r| = rmt
Cb Bugs
Cu Updates
Cu   11 Dec 08 (S. Lozovoi) job = 2 added
Cu   25 Jul 07 (S. Lozovoi) created from vsl.f
C ----------------------------------------------------------------------
      implicit none
C Input parameters
      integer job,lmaxg,lmaxh,ndim
      double precision rmt,eh
      double precision sf(ndim,(lmaxh+1)**2)
      double precision aj(0:lmaxg),dj(0:lmaxg)
C Output parameters
      double precision fval(ndim,(lmaxh+1)**2),
     .   fslo(ndim,(lmaxh+1)**2),flap(ndim,(lmaxh+1)**2)
C Local variables
      integer n0
      parameter (n0 = 10)
      integer ilm,jlm,il,im
      double precision ak(0:n0),dk(0:n0),rfac,r1

      if (lmaxg > n0) call rxi('vsl0: lmaxg gt n0, lmaxg = ',lmaxg)
      if (ndim > (n0+1)**2)
     .   call rxi('vsl0: ndim gt (n0+1)**2, ndim = ',ndim)

c --- If job=1, make Bessel functions aj and their radial derivatives dj
      if (mod(job,10) > 0) then
        call radkj(eh,rmt,lmaxg,ak,aj,dk,dj,0)
        if (mod(job,10) == 2) then
          call dcopy(lmaxg+1,ak(0),1,aj(0),1)
          call dcopy(lmaxg+1,dk(0),1,dj(0),1)
        endif
        if (mod(job/10,10) == 0 .and. lmaxg > 0) then
          rfac = 1d0
          r1 = 1d0/rmt
          do  il = 1, lmaxg
            rfac = rfac * r1
            aj(il) = aj(il) * rfac
            dj(il) = dj(il) * rfac
          enddo
        endif
      endif

c --- Make value and slope
      do  jlm = 1, (lmaxh+1)**2
        ilm = 0
        do  il = 0, lmaxg
          do  im = -il, il
            ilm = ilm + 1
            fval(ilm,jlm) = sf(ilm,jlm)*aj(il)
            fslo(ilm,jlm) = sf(ilm,jlm)*dj(il)
          enddo
        enddo
      enddo

c --- Make Laplacians using \lap J_L = - e J_L
      call dcopy(ndim*(lmaxh+1)**2,fval,1,flap,1)
      call dscal(ndim*(lmaxh+1)**2,-eh,flap,1)

      end
