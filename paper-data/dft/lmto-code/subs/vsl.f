      subroutine vsl(job,rmtg,rsmg,kmin,kmax,ikh,eh,lmaxg,nf,k0,
     .  ndim,p0,sf,fval,fslo,flap,pkl,gpkl)
C- lm-resolved value, slope and Laplacian of function F_L at |r| = rmt
C  derived from its P_kM expansion
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :1s digit
Ci         : 1, generates pkl and gpkl; otherwise uses whatever is input
Ci         :10s digit
Ci         : 0, all radial functions (pkl,gpk,fval,fslo, and flap)
Ci         :    are divided by r^l (ignored if 1s digit = 0)
Ci         :>0, Standard radial functions
Ci   rmtg  :l-dependent radii at which value, slope, and Laplacian of F_L are
Ci         :evaluated. Must either be specified for 0..lmaxg, or rmtg(0) could
Ci         :be set to some negative number.
Ci         :The latter case implies that rmtg(l) = - rmtg(0) for all l.
Ci   rsmg  :smoothing radius entering into polynomials P_kL
Ci   kmax  :polynomial cutoff: Coefficients to polynomial expansion of f
Ci         :are contained in sf(0:kmax,...)
Ci   ikh   :If <0, coefficients to Bessel expansion of f
Ci         :are contained in sf(ikh,...); see Remarks.
Ci         :Note: it is an error for ikh to be less than kmin
Ci   eh    :Energy of Bessel expansion; see Remarks.
Ci   lmaxg :l-cutoff for P_kM expansion
Ci   nf    :f_j is returned for j=1..nf
Ci   kmin,k0,p0:dimensions coefficients array sf
Ci         :kmin is normally 0, but kmin must be <0 if ikh<0.
Ci         :p0 holds dimensions for powers of the Laplacian.
Ci         :p0 must be at least 1.  Only sf(:,:,p=0:1,:) are referenced.
Ci   ndim  :dimensions sf,fval,fslo,flap; must be at least (lmaxg+1)**2.
Ci   sf    :coefficients sf(k,M,0,L) contain contribution to function f_L
Ci         :in channel M from polynomial P_kM; see Remarks
Ci         :sf(k,M,1,L) hold corresponding coefficients for lap F_L.
Ci         :Here M=1..(lmaxg+1)^2, L=1..nf, k=0..kmax
Ci         :s(:,:,2..,:) are not used.
Ci         :If ikh<0, sf(ikh,M,1,L) contain Bessel contribution to f_L;
Ci         :see Remarks.
Co Outputs
Co   fval   : lm-resolved value of F_L       |
Co   fslo   : lm-resolved slope of F_L       | at rmt (see Remarks)
Co   flap   : lm-resolved Laplacian of F_L   |
Co   pkl    : radial part of polynomials P_kl / r^l (if job's 10s digit = 0)
Co   gpkl   : radial part of [d/dr (P_kl)] / r^l    (  -    "    -         )
Co
Cr Remarks
Cr   This routine returns f(M,j): the one-center expansion of function f_j,
Cr   i.e. the projection of f_j at radius rmt in spherical onto
Cr   spherical harmonic Y_L, and its derivatives.  Thus
Cr      f_j(r) = sum(M) f_M,j * Y_M(rhat) where |r|=rmt.
Cr
Cr   The radial part of f is expanded in polynomials P_kL, Bessel functions J_L,
Cr   or a sum of them.  sf contains coefficients to both parts of the expansion.
Cr
Cr   Cases:
Cr     ikh >=0  f(M,j) is expanded strictly as
Cr              sum(k) s(k,M,0,j) * P_kM(r), with k=0,...kmax
Cr              kmax<0 => there is no contribution from this term.
Cr     ikh < 0  The following term is added to f(M,j):
Cr              s(ikh,M,0,j) * B_M(eh,rmt)
Cr
Cr   This routine also returns slopes and laplacians of f at rmt.
Cr   fval, fslo, and flap are defined so that:
Cr             F_L(r-R) = \sum_M fval(M,L) * r^m * Y_M(r)  | if job's 10s digit = 0
Cr      (d/dr) F_L(r-R) = \sum_M fslo(M,L) * r^m * Y_M(r)  | otherwise without r^m
Cr         lap F_L(r-R) = \sum_M flap(M,L) * r^m * Y_M(r)  |
Cr   are (approximately) fulfilled at |r| = rmt
Cr
Cr   in other words, fval(M,L), fslo(M,L), and flap(M,L) are M-resolved
Cr   value, slope, and laplacian at the expansion site of F_L at the head site
Cr
Cr   fval, fslo, and flap correspond to E^(0)_LL', E^(1)_LL', and E^(2)_LL'
Cr   in MvS's notes on interstitial fits
Cr
Cb Bugs
Cu Updates
Cu   26 Aug 11 Handles possible contributions from Bessel functions
Cu   15 May 08 rmtg made l-dependent (S. Lozovoi)
Cu   26 May 07 first written (S. Lozovoi)
C ----------------------------------------------------------------------
      implicit none
C Input parameters
      integer job,kmin,kmax,lmaxg,nf,k0,ndim,p0,ikh
      double precision rsmg,rmtg(0:lmaxg),eh
      double precision sf(kmin:k0,ndim,0:p0,nf)
C Output parameters
      double precision fval(ndim,nf),fslo(ndim,nf),flap(ndim,nf)
      double precision pkl(0:k0,0:lmaxg),gpkl(0:k0,0:lmaxg)
C Local variables
      integer kmax0,n0,jvec
      parameter (n0 = 10, kmax0 = 20)
      integer ilm,jlm,ik,l,m
      double precision sval,sslo,slap
c     double precision rsmgl(0:n0),rsmg0,rsx
      double precision jl(0:lmaxg),djl(0:lmaxg)
      double precision ak(0:lmaxg),dk(0:lmaxg)
      double precision rmtgl(0:n0),rmtg0,rmx
      double precision wkl(0:kmax0,0:n0+1),gwkl(0:kmax0,0:n0),wk(1)
      double precision tol
      data tol/1.d-15/

      if (lmaxg > n0) call rx('vsl: lmaxg gt n0')
      if (kmax > kmax0) call rx('vsl: kmax gt kmax0')

C --- Handle negative rmt
c     if (rsmg(0) < 0d0) then
c       call dvset(rsmgl(0),1,lmaxg+1,-rsmg(0))
c     else
c       call dcopy(lmaxg+1,rsmg(0),1,rsmgl(0),1)
c     endif
      if (rmtg(0) < 0d0) then
        call dvset(rmtgl(0),1,lmaxg+1,-rmtg(0))
      else
        call dcopy(lmaxg+1,rmtg(0),1,rmtgl(0),1)
      endif

c --- Case job=1: make polynomials pkl and their radial derivatives gpkl
      if (mod(job,10) == 1) then
        if (mod(job/10,10) == 0) then
          jvec = 10
        else
          jvec = 11
        endif
        rmx = -1d2
        do  l = lmaxg, 0, -1
          rmtg0 = rmtgl(l)
          if (dabs(rmtg0-rmx) > tol) then

c           lmax+1 needed for gradients
            call vecpkl(rmtg0,rsmg,1,kmax,l+1,1,kmax0,wk,jvec,wkl,gwkl)

            if (ikh < 0) then
              call radkj(eh,rmtg0,l,ak,jl,dk,djl,0)
              if (mod(job/10,10) == 0) then
                jl(l)  =  jl(l)/rmtg0**l
                djl(l) =  djl(l)/rmtg0**l
              endif
            endif

            rmx = rmtg0
          endif
c         do  ik = 0, kmax
c           pkl(ik,l) = wkl(ik,l)
c           gpkl(ik,l) = gwkl(ik,l)
c         enddo
          call dcopy(kmax+1,wkl(0,l),1,pkl(0,l),1)
          call dcopy(kmax+1,gwkl(0,l),1,gpkl(0,l),1)
        enddo
      endif

C --- Sum up polynomials to make value, slope, and Laplacian
      if (kmax >= 0) then
      do  jlm = 1, nf
        ilm = 0
        do  l = 0, lmaxg
          do  m = -l, l
            ilm = ilm + 1
            sval = 0d0
            sslo = 0d0
            slap = 0d0
            do  ik = 0, kmax
              sval = sval + sf(ik,ilm,0,jlm)*pkl(ik,l)
              sslo = sslo + sf(ik,ilm,0,jlm)*gpkl(ik,l)
              slap = slap + sf(ik,ilm,1,jlm)*pkl(ik,l)
            enddo
            fval(ilm,jlm) = sval
            fslo(ilm,jlm) = sslo
            flap(ilm,jlm) = slap
          enddo
        enddo
      enddo
      endif

C --- Add Bessel contribution
      if (ikh < 0) then
      if (ikh < kmin) call rxi('vsl: illegal value, ikh-',ikh)
      do  jlm = 1, nf
        ilm = 0
        do  l = 0, lmaxg
          do  m = -l, l
            ilm = ilm + 1
            sval = sf(ikh,ilm,0,jlm)*jl(l)
            fval(ilm,jlm) = fval(ilm,jlm) + sval
            fslo(ilm,jlm) = fslo(ilm,jlm) + sf(ikh,ilm,0,jlm)*djl(l)
            flap(ilm,jlm) = flap(ilm,jlm) - eh*sval
          enddo
        enddo
      enddo
      endif

      end
