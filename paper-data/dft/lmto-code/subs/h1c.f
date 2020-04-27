      subroutine h1c(ph,nlmh,eh,pg,nlmg,rhc,
     .  cy,cg,indxcg,jcg,ndim,hvl)
C- L-resolved value and Laplacian of unsmoothed Hankels centered
C  at ph over a sphere centered at pg using 1c theorem, e < 0
C  (for testing)
C ----------------------------------------------------------------
Ci Inputs
Ci   ph    :coordinates of the head site
Ci   nlmh  :L-cutoff for Hankels at the head site
Ci   eh    :l-dependent Hankel energies
Ci   pg    :coordinates of the expansion site
Ci   nlmg  :L-cutoff for Bessels at the expansion site
Ci   nds   :leading dimensions of sg
Ci   rhc   :l-dependent augmentation radii at the expansion site
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordan coefficients
Ci   ndim  :Leading dimensions of hvl
Co Outputs
Co   hvl   :L-decomposed value and Laplacian of unsmoothed Hankel
Cu Updates
Cu   18 Dec 08 (S. Lozovoi) First written
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nlmh,nlmg,ndim
      double precision ph(3),pg(3)
      double precision eh(0:*),rhc(0:*)
      double precision hvl(ndim,ndim,0:*)
      integer indxcg(*),jcg(*)
      double precision cy(*),cg(*)
C Local parameters
      integer n0
      parameter (n0=10)
      integer lmxg,lmxh
      integer il,ilm,jlm
      integer ir,iprint
      double precision ehx
      double precision sr(ndim,ndim)
      double precision srdot(ndim,ndim)
      double precision sss
      integer ll
      double precision dh0(3)
c for value-slope-laplacian program
      integer jvsl,jrad
      double precision aj(0:n0,0:n0),dj(0:n0,0:n0)
      double precision fval(ndim,ndim),fslo(ndim,ndim),flap(ndim,ndim)
      double precision tol

      data tol /1d-15/

C ------ Initialization ----------
      call tcn('h1c')
C     loka = 0
c     iprint = 40
      iprint = 30
c Use true (j_l) or reduced (j_l/r^l) radial parts in the
c value-slope-laplacian test: jrad = 0 if the former
      jrad = 0
c     jrad = 1
c           dist0 = .7d0
c... jvsl = 1 if val-slope-laplacian radial functions are divided by r^l,
c... jvsl = 11 if not
      if (jrad == 0) then
        jvsl = 10
      else
        jvsl = 0
      endif

      lmxh = ll(nlmh)
      lmxg = ll(nlmg)

c Set up energies
      do  il = 0, lmxh
c ... l-dependent energies not yet implemented
c       eh(il) = -0.7d0
        eh(il) = eh(0)
      enddo


      if (iprint > 30) then
c Print out input parameters
        write(6,810) ' h1c input:   '
        write(6,810) ' Verbosity:   ',iprint
        write(6,810) ' ndim:        ',ndim
        write(6,*) '       Function to expand:'
        write(6,800) ' centre (ph): ',ph
        write(6,810) ' lmax:        ',lmxh
        write(6,800) ' eh:          ',(eh(il), il = 0, lmxh)
        write(6,800) ' rhc:         ',(rhc(il), il = 0, lmxh)
        write(6,*) '       Expansion site:'
        write(6,800) ' centre (pg): ',pg
        write(6,810) ' lmax:        ',lmxg
        write(6,*) ' '
      endif

  800 format(1x,a,50f8.4)
  810 format(1x,a,50i4)


c --- Find coefficients

      ehx = eh(0)
      sss = 0.d0
      do  ir = 1, 3
        dh0(ir) = - (pg(ir)-ph(ir))
        sss = sss + dh0(ir)**2
      enddo

c...deb
c     write(6,*) 'h1c: Distance between the origin and expansion sites:'
c     write(6,800) ' h1c: |pg-ph| =     ',dsqrt(sss)
c...deb

      if (sss > tol) then
        call mstrx2(ehx,dh0,nlmg,nlmh,ndim,cg,indxcg,jcg,cy,1,sr,srdot)

        if (iprint >= 50) then
          print *,' '
          print *,
     .     " Expansion coefficients L, s(L',L), L' = 1,nlmg :"
          print *,' '
          do  jlm = 1, nlmh
              write (6,841) jlm,(sr(ilm,jlm),ilm=1,min0(9,nlmg))
          enddo
        endif

         call vsl0(jvsl+1,rhc(0),eh,lmxg,lmxh,ndim,sr,
     .     fval,fslo,flap,aj,dj)

  841    format(i3,1x,16g10.2)
      else
c if sites coinside sr(i,j) = delta_i,j
        call dpzero(sr,ndim*ndim)
        do  jlm = 1, max(nlmg,nlmh)
          sr(jlm,jlm) = 1d0
        enddo
         call vsl0(jvsl+2,rhc(0),eh,lmxg,lmxh,ndim,sr,
     .     fval,fslo,flap,aj,dj)
      endif

      call dpzero(hvl,2*ndim*ndim)
      do  jlm = 1, nlmh
        do  ilm = 1, nlmg
          hvl(ilm,jlm,0) = fval(ilm,jlm)
          hvl(ilm,jlm,1) = flap(ilm,jlm)
        enddo
      enddo

      call tcx('h1c')
      end
