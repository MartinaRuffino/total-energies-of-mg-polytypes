      subroutine maadot(s_spec,nl,nbas,mode,ips,avw,alpha,adot,
     .  tral,trad)
C- Makes vector of alpha, adot and tral matrix 2nd gen kappa=0 strux
C ----------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxb rmt alpha
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci   nl    :(global maximum l) + 1
Ci   nsp   :number of spins (used only for 1s digit mode=1)
Ci   nbas  :size of basis
Ci   rw    :rmax (used to adot, 2nd generation modes 0 and 1 only)
Ci   mode  :1s digit 0: alpha=alp0, (which is input)
Ci         :         1: alpha=gamma, calculated from pot pars.
Ci                      NO LONGER SUPPORTED
Ci         :         2: alpha calculated for hard spheres.
Ci                      NO LONGER SUPPORTED; see mktral
Ci         10s digit 0: adot=Andersen's localized
Ci         :         1: adot so Kanpur notes Eq 3.87 is zero
Ci         :         2: adot calculated for hard spheres.
Ci                      NO LONGER SUPPORTED; see mktral
Co Outputs
Co   alpha :tight-binding screening parameters
Co   adot  :(kappa*avw)^2-derivative of alpha
Co   tral  :tral matrix, for compatibility with NMTO; see mktral
Co   trad  :trad matrix, for compatibility with NMTO; see mktral
Cr Remarks
Cr   This routine produces screening parameters for second-generation
Cr   kappa=0 structure constants.
Cu Updates
Cu   06 Sep 11 Started migration to f90 structures
Cu   31 Aug 00 Eliminated extensions to 2nd generation data
C ----------------------------------------------------------------
      use structures
      implicit none
C Passed parameters
      integer nl,nbas,mode,ips(nbas)
      double precision avw,alpha(nl*nl,nbas),adot(nl*nl,nbas),
     .  tral(4,nl**2,nbas),trad(4,nl**2,nbas)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C Local parameters
      integer l,m,lm,i,l2,i1mach,iprint,moda,modad,getdig,n0,is,lmxb
      parameter(n0=10)
      double precision t,alp0(n0),rw

      moda  = getdig(mode,0,10)
      modad = getdig(mode,1,10)

C --- For each atom in the basis, do ---
      do  i = 1, nbas
        is = ips(i)
        lmxb = s_spec(is)%lmxb
        rw = s_spec(is)%rmt
        alp0(1:1+lmxb) = s_spec(is)%alpha(1:1+lmxb)
C  Patch for now, floating orbitals ...
        if (rw == 0) then
          rw = 0.01d0
        endif
        rw = rw/avw
        lm = 0
        do  l = 0, nl-1
        do  m = -l, l
          l2 = l+l
          lm = lm+1
          if (l > lmxb) then
            alpha(lm,i) = 1d-6

          else

C       ... alpha
            if (moda == 0) then
              alpha(lm,i) = alp0(1+l)
            else
              call rxi(' MAADOT: mode not supported:',mode)
            endif

C       ... adot
            if (modad == 0) then
              t = 12*(l2+1)*(l2+3) / dble((l2-1)*(l2+5)*(l2+7))
              adot(lm,i) = -alpha(lm,i)*
     .          (rw**2/(l2+5) + alpha(lm,i)*t/rw**(l2-1))
            elseif (modad == 1) then
              t = -rw**(l2+3)/(2d0*(l2+1)**2*(l2+3)) +
     .          alpha(lm,i)*(rw**2/(l2+1) +
     .          alpha(lm,i)*2/rw**(l2-1)/(l2-1))
              adot(lm,i) = -t
            else
              call rxi(' MAADOT: mode not supported:',mode)
            endif
          endif

C     ... tral matrix, for compatibility with 3rd generation
          tral(1,lm,i) = 1d0
          tral(2,lm,i) = 0d0
          tral(3,lm,i) = -alpha(lm,i)
          tral(4,lm,i) = 1d0
          trad(1,lm,i) = 0d0
          trad(2,lm,i) = 0d0
          trad(3,lm,i) = -adot(lm,i)
          trad(4,lm,i) = 0d0

          enddo
        enddo

C   --- Output ---
        if (iprint() >= 50) then
          call awrit2(' MAADOT: alpha for atom %i (mode=%i)',
     .      ' ',80,i1mach(2),i,mod(mode,10))
          write(*,885) (alpha(lm,i), lm=1,nl*nl)
  885     format(5f12.6)
          call awrit2(' MAADOT: a-dot for atom %i (mode=%i)',
     .      ' ',80,i1mach(2),i,mode/10)
          write(*,885) (adot(lm,i), lm=1,nl*nl)
        endif

      enddo
    1 format(5F12.6)

      end
