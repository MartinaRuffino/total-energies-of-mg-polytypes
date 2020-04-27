      subroutine pp2enu(s_spec,nclass,ics,nrclas,nsp,nl,rmax,avw,
     .  amom,pnu,qnu,sumevm,initc,dclabl,enu,pp)
C- Interactively shift pot pars to specified enu
C ----------------------------------------------------------------
Ci Inputs
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: idmod
Co     Stored:    idmod
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci   nclass:number of inequivalent classes
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   nrclas:number of atoms in the ith class
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :(global maximum l) + 1
Ci   rmax  :augmentation radius, in a.u.,
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   amom  :input moments qnu
Ci   initc : record of what parameters are available.
Ci           1 P,Q   2 pp   4 sop   8 vintra  16 pmpol  32 gradrm
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :final moments qnu
Ci   enu: default linearization energy to take.
Co Outputs
Co   pp   : enu set for each channel user elects to change
Co   idmod: set to 2 for each change user elects to change .
Co  sumevm:
Cu Updates
Cu   05 Oct 01 Adapted to work with lm v6.11
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nclass,nsp,nl,nrclas(nclass),initc(nclass)
      integer ics(nclass)
      integer n0
      parameter (n0=10)
      integer idmod(n0)
      double precision amom(3,nl,nsp,nclass),enu,dclabl(nclass)
      double precision pnu(nl,nsp,nclass),qnu(3,nl,nsp,nclass)
      double precision pp(6,nl,nsp,nclass),rmax(nclass),avw,sumevm
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      double precision xx,eb
      integer ic,isp,l,ln,lgunit,lc,awrite,j1,j2,idmodx,i,parg,is
      external lgunit
      character outs*256,strn*72

      call awrit1(' pp2enu: for each channel, enter one of the'//
     .  ' following:%N'//
     .  '%9f<ret> to set enu to %,1d%N'//
     .  '%7fk <ret> to keep existing enu%N'//
     .  '%5fval <ret> to shift enu to position ''val''',
     .  outs,len(outs),lgunit(1),enu)
      lc = 8
      call awrit1('%N  spec%nfl    Old enu   Default  Set to',
     .  outs,len(outs),lgunit(1),lc)
      sumevm = 0d0
      do  ic = 1, nclass
      is = ics(ic)
      idmod = s_spec(is)%idmod

      strn = ' '
      call r8tos8(dclabl(ic),strn)
      if (mod(initc(ic)/2,2) == 1) then
        do  isp = 1, nsp
        do  l = 1, nl
          idmodx = idmod(l)
          xx = pp(1,l,isp,ic)
          outs = ' '
          ln = awrite('%,3i '//strn(1:lc)//'  %i  %;9,5D%;10,5D%3f?',
     .      outs,80,0,ic,l-1,xx,enu,xx,xx,xx,xx)
          call cwrite(outs,0,ln,0)
          read(*,'(a80)') outs
          call word(outs,1,j1,j2)
          if (j2 < j1) then
            xx = enu
            idmodx = 2 + (idmod(l)-mod(idmod(l),10))
          elseif (outs(j1:j2) /= 'k') then
            i = parg(' ',4,outs,j1-1,len(outs),' ',1,1,j1,xx)
            if (i == -1)
     .        call rxs('PP2ENU: failed to parse input',outs)
            idmodx = 2 + (idmod(l)-mod(idmod(l),10))
          endif
          if (idmodx /= idmod(l) .or. xx /= pp(1,l,isp,ic)) then
            call awrit3('%nf... set  idmod = %i  enu = %,1d',
     .        ' ',80,lgunit(1),lc+27,idmodx,xx)
            enu = xx
            eb = enu - pp(1,l,isp,ic)
            call enutod(.false.,nl,l,ic,rmax(ic),avw,pp(1,l,isp,ic),amom(1,l,isp,ic),[xx],
     .        0d0,idmodx,0d0,pnu(l,isp,ic),qnu(1,l,isp,ic),[xx],eb)
            idmod(l) = idmodx
            pp(1,l,isp,ic) = enu
          endif
          sumevm = sumevm + (qnu(2,l,isp,ic) +
     .                     qnu(1,l,isp,ic)*pp(1,l,isp,ic))*nrclas(ic)
          s_spec(is)%idmod(1:nl) = idmod(1:nl)
        enddo
        enddo
      else
        call awrit1('%,3i '//strn(1:lc)//'  is missing ppars',' ',80,
     .    lgunit(1),ic)
      endif
      enddo
      end
