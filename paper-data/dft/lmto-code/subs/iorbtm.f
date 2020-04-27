      subroutine iorbtm(s_spec,ltso,ls,ics,nl,nlo,nclass,nsp,orbtm,sopertp,tso)
C- Printout of orbital moments and possibly site-resolved SO energy
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: name
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   ltso  :if 1, tso contains decomposition of SO coupling by site
Ci   ls    :if 0, data supplied by class
Ci         :   1, data supplied by site
Ci   ics   :species table:
Ci         :class ic belongs to species ics(ic) (ls=0)
Ci         :site ic belongs to species ics(ic) (ls=1)
Ci   nl    :(global maximum l) + 1
Ci   nlo   :number of l or lm
Ci   nclass:number of inequivalent classes
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   orbtm :orbital moments, ordered m=l..-l
Ci  sopertp:parameters used to estimate Fermi surface contribution to change in band sum.
Ci         :sopertp is used only if ltso=1
Ci          (1) sampling Fermi level
Ci          (2) sampling sumev
Ci          (3) sampling DOS(Ef)
Ci          (4) sum [delta w] sum of change in weights, 1st order perturbation
Ci          (5) Same as (4), but delta w from 2nd order pert theory
Ci          (6) sum [delta (w e)] change in weight-eval product, sampling
Ci          (7) sum [w (delta e)] weighted change in evals by sampling
Ci          (8) sum [w (delta e)] (weighted change in evals by tetrahedron
Ci          (9) Same as (6), but delta e from 2nd order pert theory
Ci         (10) Same as (7), but delta e from 2nd order pert theory
Ci         (11) Same as (8), but delta e from 2nd order pert theory
Ci   tso   :If ltso=1,
Ci         :site-resolved contribution to band energy from SO coupling
Co Outputs
Cr Remarks
Cr   Should be a term related to Fermi level shift:
Cr   delta(wi Ei) = wi(new)*Ei(new) - wi(old)*Ei(old)
Cr   wi = wi(eF) + dw/dEi dEi + dw/dEf dEf
Cr   dEf = dEf/dq dq = 1/D(Ef) dq
Cr   2nd term in dw then dw = dw/dEf 1/D(Ef) dq
Cr   Contributes a term:   sum_i dwi/dEf 1/D(Ef) dq ei
Cu Updates
Cu    9 Aug 13 Print SO coupling energies, if available
Cu   10 Nov 11 Begin migration to f90 structures
Cu   09 Aug 04 (A. Chantis) Correct sign of orbl
Cu   08 Dec 00 First implementation
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ltso,ls,nlo,nl,nsp,nclass,ics(*)
      double precision orbtm(nlo,nsp,nclass),tso(2,2,2,0:nclass)
      double precision sopertp(12)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ic,isp,l,im,lm,m,l1,ipr,stdo,is,nglob
      double precision amom,orbl(10),efermis,xx
      character*8 slabl
      character*5 strni(2)
      data strni /'class','site'/
C ... External calls
      external dpzero,getpr,spacks

      call getpr(ipr)
      if (ipr < 20) return
      stdo = nglob('stdo')
      write(stdo,332) strni(ls+1)
  332 format(/' IORBTM:  orbital moments :'/1x,a5,
     .  ' spec        spin   Moment decomposed by l ...')

      do  ic = 1, nclass
        is = ics(ic)
        if (s_spec(is)%lmxa < 0) cycle
        slabl = s_spec(is)%name
        amom = 0
        do  isp = 1, nsp
          call dpzero(orbl,nl)
          lm = 0
          do  l = 0, nl-1
            l1 = l+1
            im = l
            if (nl == nlo) im = 0
            do  m = -im, im
              lm = lm+1
              orbl(l1) = orbl(l1) + orbtm(lm,isp,ic)
              amom = amom + orbtm(lm,isp,ic)
            enddo
          enddo
C          write(stdo,333) ic,slabl,isp,(orbl(l1),l1=1,nl)
C  333     format(i5,4x,a8,i6,8f12.6)
          call info5(10,0,0,'%,5i    '//slabl//'%,6i%n;12,6D',ic,isp,nl,orbl,0)
        enddo
        if (ltso /= 0) then
          call info5(20,0,0,'%,5i  LS++ %;13,8D  +- %;13,8D  '//
     .      '-+ %;13,8D  -- %;13,8D',ic,
     .      tso(1,1,1,ic),tso(1,1,2,ic),tso(1,2,1,ic),tso(1,2,2,ic))
          call info5(20,0,0,'%,5i  L+ - L- %;10,6D'//
     .      '   <L.S> =%;13,8D  2nd order pert =%;13,8D',
     .      ic,amom,sum(tso(1,1:2,1:2,ic)),tso(2,1,1,ic),0)
        else
          call info5(20,0,0,'%,5i  L+ - L- %;10,6D',ic,amom,ltso,0,0)
        endif
      enddo
      if (ltso /= 0) then
        call info5(20,0,0,' Total SO coupling: <L.S>%;13,8D  '//
     .    '2nd order pert%;13,8D   sum-of-sites%;13,8D',
     .    sum(tso(1,1:2,1:2,0)),tso(2,1,1,0),
     .    sum(tso(2,1,1,1:nclass)),0,0)
        if (sopertp(8) /= 0) then
        efermis = sopertp(1) ! Fermi energy from sampling integration

        xx = sopertp(6)-sopertp(4)*efermis*0-sopertp(7)
        call info5(20,0,0,' FS: '//
     .  '<tet>%;13,8D  <sam>%;13,8D    wt corr%;13,8D'//
     .  '  <tet>+wt corr%;13,8D',
     .  sopertp(8),sopertp(7),xx,sopertp(8)+xx,0d0)

        xx = sopertp(9)-sopertp(5)*efermis*0-sopertp(10)
        call info5(20,0,0,' FS2 '//
     .  '<tet>%;13,8D  <sam>%;13,8D    wt corr%;13,8D'//
     .  '  <tet>+wt corr%;13,8D',
     .  sopertp(11),sopertp(10),xx,sopertp(11)+xx,0d0)
      endif

      endif

      end
