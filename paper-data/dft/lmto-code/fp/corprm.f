      subroutine corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,
     .  lfoc,rfoc,z)
C- Returns parameters for smooth core+nucleus representation
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lfoca rfoca qc z ctail etail stc lmxb p pz rmt rg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   is    :species index
Co Outputs
Co   cofg  :coefficient to Gaussian part of pseudocore density
Co         :assigned so that pseudocore charge = true core charge
Co   cofh  :coefficient to Hankel part of pseudocore density
Co         :Hankel contribution is determined by inputs
Co         :(qcorh,ceh,rfoc) and should accurately represent the
Co         :true core density for r>rmt
Co   qcorg :charge in the gaussian part; see Remarks
Co   qcorh :charge in the Hankel part; see Remarks
Co   qsc   :number of electrons in semicore treated by local orbitals
Co   lfoc  :switch specifying treatment of core density.
Co          0 => val,slo = 0 at sphere boundary
Co          1 => core tails included explicitly with valence
Co          2 => tails included perturbatively
Co
Co   rfoc :smoothing radius for hankel head fitted to core tail
Co   z     :nuclear charge
Cr Remarks
Cr   qcorg and qcorh are the charges in the gaussian and hankel parts.
Cr   The hankel part is used when the core is allowed to spill out of
Cr   the augmentation sphere.
Cr
Cr   cofg and cofh are the coefficients in front of the standard
Cr   gaussian and smoothed hankel functions for l=0.
Cr   That is: the pseudocore density is
Cr      cofg*g0(rg;r)*Y0 + cofh*h0(rfoca;r)*Y0        (1)
Cr   ceh and rfoc are the sm- hankel energy and radius for the hankel term
Cr   cofg is set so that qc = integral of eq. 1 above.
Cr
Cr   For lfoc=0 there is no Hankel part; qc carried entirely by Gausian
Cr   For lfoc>0 there is a  Hankel part; Gaussian carries difference
Cr              between qc and charge in Hankel part.
Cr
Cr   To add to the radial density 4*pi*r**2*rho_true, multiply
Cr   cofg,cofh by srfpi.
Cl Local variables
Cl    ccof :coefficient for core tail, for a smoothed Hankel.
Cl          ccof is differs from spec->ctail because ctail is
Cl          constructed for an unsmoothed Hankel.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Sep 01 Generates qsc.  Argument list changed.
Cu   24 Apr 00 Adapted from nfp corpars.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer is
      double precision qcorg,qcorh,qsc,cofg,cofh,ceh,rfoc,z
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer n0,lfoc,lmxb,l
      parameter (n0=10)
      double precision ccof,fpi,q0,q1,qc,rmt,rsm,srfpi,stc,y0,x0(0:n0),
     .  xi(0:n0),pnu(n0,2),pz(n0,2)

      fpi = 16d0*datan(1d0)
      srfpi = dsqrt(fpi)
      y0 = 1d0/srfpi

C ... Extract parameters including qc and qsc
      lfoc = s_spec(is)%lfoca
      rfoc = s_spec(is)%rfoca
      qc = s_spec(is)%qc
      z = s_spec(is)%z
      ccof = s_spec(is)%ctail
      ceh = s_spec(is)%etail
      stc = s_spec(is)%stc
      lmxb = s_spec(is)%lmxb
      pnu = s_spec(is)%p
      pz = s_spec(is)%pz
      rmt = s_spec(is)%rmt
      if (rfoc <= 1d-5) rfoc = s_spec(is)%rg
      qsc = 0
      do  l = 0, lmxb
        if (int(pz(l+1,1)) /= 0) then
          if (int(mod(pz(l+1,1),10d0)) < int(pnu(l+1,1)))
     .      qsc = qsc + 4*l+2
        endif
      enddo

C ... Scale smoothed hankel coeff for exact core spillout charge
C     q1 = spillout charge in sm. Hankel fitting core density
C     q0 = spillout charge in Hankel with rsm = 0.05
      if (ccof /= 0) then
      call hansmr(rmt,0d0,1/rfoc,x0,1)
      call hansmr(rmt,ceh,1/rfoc,xi,1)
      q1 = srfpi/ceh*(-dexp(rfoc**2/4*ceh)
     .   - rmt**3*(xi(1)-dexp(rfoc**2/4*ceh)*x0(1)))
      rsm = 0.05d0
      call hansmr(rmt,0d0,1/rsm,x0,1)
      call hansmr(rmt,ceh,1/rsm,xi,1)
      q0 = srfpi/ceh*(-dexp(rsm**2/4*ceh)
     .   - rmt**3*(xi(1)-dexp(rsm**2/4*ceh)*x0(1)))
      q0 = q0*y0
      q1 = q1*y0
      ccof = ccof*q0/q1
      endif

C ... Set gaussian and hankel charges
      qcorg = qc
      qcorh = 0d0
      if (lfoc > 0) then
        qcorh = -ccof*dexp(ceh*rfoc*rfoc/4d0)/ceh
        qcorg = qc-qcorh
      endif

C ... Coeffients to the the gaussian and hankel terms
      cofh = -y0*qcorh*ceh*dexp(-ceh*rfoc*rfoc/4d0)
      cofg = y0*qcorg

c      write (6,352) is,qcorg,qcorh,cofg,cofh
c  352 format(' spec',i3,'  qcorg,qcorh=',2f10.6,'  cofg,cofh=',2f12.4)

      end